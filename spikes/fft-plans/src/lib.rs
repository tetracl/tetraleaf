//! Spike S9 (Helmholtz CH0.9): per-box FFT plans of the φ step, their bits and their time.
//!
//! The φ step of S → F and F → R (proposal §5.4, §6.7) maps, on each of the 2P − 1 rings of a
//! box's grid (P = p + 1), the 2P − 1 Fourier modes μ = −p..p to N_φ azimuthal samples and
//! back. This crate does it three ways on one layout:
//!
//! - [`FftwPhi`]: FFTW through RLST 0.12's `rlst::fftw` (`fftw_system`), the code path of
//!   CH4.6: one plan per (N_φ, rings, direction), `Rigor::Estimate`, executed box by box on
//!   the box's own slice of a packed, 64-byte aligned buffer; interleaved or split storage.
//! - [`RustPhi`]: rustfft, one plan per (N_φ, direction), executed ring by ring of a box.
//! - [`naive_synthesis`], [`naive_analysis`]: the DFT written as a sum, the trusted path of
//!   the accuracy checks.
//!
//! The product as shipped (`tetraleaf_tables::helmholtz::diagonal::Conversions`) fuses the
//! φ step with the Legendre step; [`ProductBoxes`] times it whole, and its φ share comes from
//! `Conversions::multiply_adds` (CH4.1a's method).
//!
//! # Layout
//!
//! A box's buffer holds `rings` rings of `azimuths` interleaved complex samples, ring after
//! ring, φ contiguous: sample j of ring i at `i · azimuths + j`, the order of a signature
//! (CONVENTIONS §3.15). [`Shape::rings`] is 2P − 1, or 2P with a last ring of zeros (the even
//! batch §6.7 leaves to S9). A box's modes are its 2P − 1 live rings of 2P − 1 modes,
//! ring-major, mode k at μ = k − p, which the φ step puts in slot μ mod N_φ. Boxes are packed
//! back to back in one allocation aligned to `rlst::FFTW_ALIGNED` = 64 bytes; a box takes
//! 16 · azimuths · rings bytes, a multiple of 64 since N_φ is a multiple of 4, so every box
//! starts in alignment class 0 ([`Packed::box_class`]).
//!
//! # Directions
//!
//! Synthesis (S → F's φ step) is F(φⱼ) = Σ_μ c_μ e^{+i μ φⱼ}, φⱼ = 2π j / N_φ: FFTW's and
//! rustfft's backward transform of the scattered modes. Analysis (F → R's) is
//! c_μ = Σⱼ F(φⱼ) e^{−i μ φⱼ}: the forward transform, gathered. Neither is scaled.

use std::sync::Arc;

use rlst::dense::array::AlignedDynArray;
use rlst::fftw::{
    AlignmentClass, C2c, Dim, Direction, FftwError, Layout, Placement, Plan, Planner, Preserving,
    SplitC2c,
};
use rlst::{FFTW_ALIGNED, c64};
use rustfft::{Fft, FftPlanner};
use tetraleaf_ref::helmholtz::diagonal::azimuth_count;
use tetraleaf_tables::helmholtz::diagonal::{ConversionScratch, Conversions};

/// The degrees of the timing runs: P = p + 1 from 34 to 262, the high-frequency levels'
/// degrees of CH4.0e–CH4.5's runs.
pub const DEGREES: [usize; 12] = [33, 41, 50, 66, 82, 100, 123, 135, 144, 200, 245, 261];

/// The azimuthal grid of a level.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Grid {
    /// The product's grid: 4P − 3 rounded up to a multiple of 4 (`azimuth_count`), the
    /// grid the product runs on today; its N_φ / 4 may have a large prime factor.
    Product,
    /// The smallest multiple of 4 at or above 4P − 3 whose other prime factors are at most 7
    /// (proposal §6.7, "Grid sizes").
    Smooth,
    /// The smallest power of two at or above 4P − 3.
    PowerOfTwo,
}

impl Grid {
    /// Every grid, in the order of the tables.
    pub const ALL: [Grid; 3] = [Grid::Product, Grid::Smooth, Grid::PowerOfTwo];

    /// N_φ at degree p.
    pub fn azimuths(self, p: usize) -> usize {
        let least = 4 * p + 1; // 4P − 3
        match self {
            Grid::Product => azimuth_count(p),
            Grid::Smooth => (least.div_ceil(4)..)
                .find(|&q| is_smooth(q))
                .map(|q| 4 * q)
                .unwrap(),
            Grid::PowerOfTwo => least.next_power_of_two(),
        }
    }

    /// The name in the tables and on the command line.
    pub fn name(self) -> &'static str {
        match self {
            Grid::Product => "product",
            Grid::Smooth => "smooth",
            Grid::PowerOfTwo => "pow2",
        }
    }

    /// The grid of a name.
    pub fn parse(name: &str) -> Option<Self> {
        Self::ALL.into_iter().find(|g| g.name() == name)
    }
}

/// Whether n has no prime factor above 7.
pub fn is_smooth(mut n: usize) -> bool {
    for f in [2, 3, 5, 7] {
        while n.is_multiple_of(f) {
            n /= f;
        }
    }
    n == 1
}

/// The shape of a box's φ step.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Shape {
    /// The degree p.
    pub p: usize,
    /// N_φ.
    pub azimuths: usize,
    /// The rings of a box's buffer: 2P − 1, or 2P with a zero ring.
    pub rings: usize,
}

impl Shape {
    /// The shape at degree p on a grid, with 2P − 1 rings or (`even`) 2P.
    pub fn new(p: usize, grid: Grid, even: bool) -> Self {
        Self {
            p,
            azimuths: grid.azimuths(p),
            rings: 2 * p + 1 + usize::from(even),
        }
    }

    /// The live rings and the modes per ring: 2P − 1.
    pub fn modes(&self) -> usize {
        2 * self.p + 1
    }

    /// A box's modes: (2P − 1)².
    pub fn mode_len(&self) -> usize {
        self.modes() * self.modes()
    }

    /// A box's samples: N_φ · rings.
    pub fn box_len(&self) -> usize {
        self.azimuths * self.rings
    }

    /// The slot of mode k (μ = k − p) on a ring.
    pub fn slot(&self, k: usize) -> usize {
        (k + self.azimuths - self.p) % self.azimuths
    }

    /// Writes a box's modes into its slots of `buffer` (interleaved), leaving every other
    /// element as it is.
    pub fn scatter(&self, modes: &[c64], buffer: &mut [c64]) {
        let m = self.modes();
        for (ring, row) in modes.chunks_exact(m).enumerate() {
            let line = &mut buffer[ring * self.azimuths..(ring + 1) * self.azimuths];
            // Slots 0..=p take μ = 0..=p, slots N_φ − p.. take μ = −p..−1.
            line[..=self.p].copy_from_slice(&row[self.p..]);
            line[self.azimuths - self.p..].copy_from_slice(&row[..self.p]);
        }
    }

    /// Reads a box's modes from its slots of `buffer` (interleaved).
    pub fn gather(&self, buffer: &[c64], modes: &mut [c64]) {
        let m = self.modes();
        for (ring, row) in modes.chunks_exact_mut(m).enumerate() {
            let line = &buffer[ring * self.azimuths..(ring + 1) * self.azimuths];
            row[self.p..].copy_from_slice(&line[..=self.p]);
            row[..self.p].copy_from_slice(&line[self.azimuths - self.p..]);
        }
    }

    /// As [`scatter`](Self::scatter), into split planes: the real parts in `buffer[..box_len]`,
    /// the imaginary parts in `buffer[box_len..]`.
    pub fn scatter_split(&self, modes: &[c64], buffer: &mut [f64]) {
        let (re, im) = buffer.split_at_mut(self.box_len());
        let m = self.modes();
        for (ring, row) in modes.chunks_exact(m).enumerate() {
            let base = ring * self.azimuths;
            for (k, z) in row.iter().enumerate() {
                let at = base + self.slot(k);
                re[at] = z.re;
                im[at] = z.im;
            }
        }
    }

    /// As [`gather`](Self::gather), from split planes.
    pub fn gather_split(&self, buffer: &[f64], modes: &mut [c64]) {
        let (re, im) = buffer.split_at(self.box_len());
        let m = self.modes();
        for (ring, row) in modes.chunks_exact_mut(m).enumerate() {
            let base = ring * self.azimuths;
            for (k, z) in row.iter_mut().enumerate() {
                let at = base + self.slot(k);
                *z = c64::new(re[at], im[at]);
            }
        }
    }
}

/// Boxes of one shape packed back to back in one allocation aligned to 64 bytes (`rlst`'s
/// `AlignedDynArray` at `FFTW_ALIGNED`), each box `stride` elements from the last.
pub struct Packed<E> {
    data: AlignedDynArray<E, 1, FFTW_ALIGNED>,
    stride: usize,
    boxes: usize,
}

impl<E: Copy + Default> Packed<E> {
    /// `boxes` boxes of `stride` elements, zeroed.
    pub fn new(stride: usize, boxes: usize) -> Self {
        Self {
            data: AlignedDynArray::<E, 1, FFTW_ALIGNED>::from_shape([stride * boxes.max(1)]),
            stride,
            boxes,
        }
    }

    /// The number of boxes.
    pub fn boxes(&self) -> usize {
        self.boxes
    }

    /// Box b.
    pub fn get(&self, b: usize) -> &[E] {
        &self.data.data().unwrap()[b * self.stride..(b + 1) * self.stride]
    }

    /// Box b, mutable.
    pub fn get_mut(&mut self, b: usize) -> &mut [E] {
        &mut self.data.data_mut().unwrap()[b * self.stride..(b + 1) * self.stride]
    }

    /// Every element.
    pub fn all(&self) -> &[E] {
        self.data.data().unwrap()
    }

    /// Every element, mutable.
    pub fn all_mut(&mut self) -> &mut [E] {
        self.data.data_mut().unwrap()
    }

    /// The boxes, mutable, one slice each (for threads).
    pub fn boxes_mut(&mut self) -> std::slice::ChunksExactMut<'_, E> {
        let stride = self.stride;
        self.all_mut().chunks_exact_mut(stride)
    }

    /// The alignment class of box b's first element.
    pub fn box_class(&self, b: usize) -> AlignmentClass {
        AlignmentClass::of_ptr(self.get(b).as_ptr())
    }
}

/// The layout of one box: `rings` transforms of length N_φ, stride 1, ring distance `distance`
/// (N_φ for the packed layout), out of place.
pub fn box_layout(shape: &Shape, distance: usize) -> Result<Layout, FftwError> {
    let n = shape.azimuths as isize;
    let d = distance as isize;
    Layout::new(
        &[Dim { n, is: 1, os: 1 }],
        &[Dim {
            n: shape.rings as isize,
            is: d,
            os: d,
        }],
        Placement::OutOfPlace,
    )
}

/// The layout of `boxes` boxes in one plan (howmany = boxes · rings): for the record only.
pub fn batch_layout(shape: &Shape, boxes: usize) -> Result<Layout, FftwError> {
    let n = shape.azimuths as isize;
    Layout::new(
        &[Dim { n, is: 1, os: 1 }],
        &[Dim {
            n: (shape.rings * boxes) as isize,
            is: n,
            os: n,
        }],
        Placement::OutOfPlace,
    )
}

/// FFTW's φ step of one shape: a synthesis and an analysis plan for one box, interleaved.
pub struct FftwPhi {
    /// The shape.
    pub shape: Shape,
    /// Backward: scattered modes to samples.
    pub synthesis: Plan<f64, C2c, Preserving>,
    /// Forward: samples to slots.
    pub analysis: Plan<f64, C2c, Preserving>,
}

impl FftwPhi {
    /// Plans for a box of `shape`, alignment class 0, `Rigor::Estimate`.
    pub fn new(shape: Shape) -> Result<Self, FftwError> {
        Self::with_planner(shape, &Planner::<f64>::new(), shape.azimuths)
    }

    /// Plans with a given planner (its alignment classes, rigor) and ring distance.
    pub fn with_planner(
        shape: Shape,
        planner: &Planner<f64>,
        distance: usize,
    ) -> Result<Self, FftwError> {
        let layout = box_layout(&shape, distance)?;
        Ok(Self {
            shape,
            synthesis: planner.c2c(&layout, Direction::Backward)?,
            analysis: planner.c2c(&layout, Direction::Forward)?,
        })
    }

    /// A per-thread scratch of two boxes: box 0 the synthesis input, zeroed once, so that
    /// only its mode slots are ever written; box 1 the analysis output. (One box for both
    /// leaves the analysis output in the slots the synthesis takes as zeros.)
    pub fn scratch(&self) -> Packed<c64> {
        Packed::<c64>::new(self.shape.box_len(), 2)
    }

    /// S → F's φ step of one box: scatter `modes` into `scratch`, transform into `samples`.
    pub fn synthesize(&self, modes: &[c64], scratch: &mut Packed<c64>, samples: &mut [c64]) {
        self.shape.scatter(modes, scratch.get_mut(0));
        self.synthesis.execute(scratch.get(0), samples).unwrap();
    }

    /// F → R's φ step of one box: transform `samples` into `scratch`, gather `modes`.
    pub fn analyze(&self, samples: &[c64], scratch: &mut Packed<c64>, modes: &mut [c64]) {
        self.analysis.execute(samples, scratch.get_mut(1)).unwrap();
        self.shape.gather(scratch.get(1), modes);
    }

    /// The fingerprints of the two plans (synthesis, analysis).
    pub fn fingerprints(&self) -> [u64; 2] {
        [self.synthesis.fingerprint(), self.analysis.fingerprint()]
    }
}

/// FFTW's φ step on split storage: a box is a plane of real parts and one of imaginary parts,
/// each N_φ · rings reals.
pub struct FftwSplitPhi {
    /// The shape.
    pub shape: Shape,
    /// Backward.
    pub synthesis: Plan<f64, SplitC2c, Preserving>,
    /// Forward.
    pub analysis: Plan<f64, SplitC2c, Preserving>,
}

impl FftwSplitPhi {
    /// Plans for a split box of `shape`, alignment class 0, `Rigor::Estimate`.
    pub fn new(shape: Shape) -> Result<Self, FftwError> {
        let n = shape.azimuths as isize;
        let gap = shape.box_len() as isize;
        let layout = Layout::new_split(
            &[Dim { n, is: 1, os: 1 }],
            &[Dim {
                n: shape.rings as isize,
                is: n,
                os: n,
            }],
            (gap, gap),
            Placement::OutOfPlace,
        )?;
        let planner = Planner::<f64>::new();
        Ok(Self {
            shape,
            synthesis: planner.split_c2c(&layout, Direction::Backward)?,
            analysis: planner.split_c2c(&layout, Direction::Forward)?,
        })
    }

    /// The reals of a split box: 2 · N_φ · rings.
    pub fn box_reals(&self) -> usize {
        2 * self.shape.box_len()
    }

    /// A per-thread scratch of two split boxes, as [`FftwPhi::scratch`].
    pub fn scratch(&self) -> Packed<f64> {
        Packed::<f64>::new(self.box_reals(), 2)
    }

    /// S → F's φ step of one split box.
    pub fn synthesize(&self, modes: &[c64], scratch: &mut Packed<f64>, samples: &mut [f64]) {
        self.shape.scatter_split(modes, scratch.get_mut(0));
        self.synthesis.execute(scratch.get(0), samples).unwrap();
    }

    /// F → R's φ step of one split box.
    pub fn analyze(&self, samples: &[f64], scratch: &mut Packed<f64>, modes: &mut [c64]) {
        self.analysis.execute(samples, scratch.get_mut(1)).unwrap();
        self.shape.gather_split(scratch.get(1), modes);
    }

    /// The fingerprints of the two plans (synthesis, analysis).
    pub fn fingerprints(&self) -> [u64; 2] {
        [self.synthesis.fingerprint(), self.analysis.fingerprint()]
    }
}

/// rustfft's φ step of one shape: one plan per direction of length N_φ, run over a box's live
/// rings in place (`rings` beyond 2P − 1 are not transformed).
pub struct RustPhi {
    /// The shape.
    pub shape: Shape,
    /// Backward (inverse) of length N_φ.
    pub synthesis: Arc<dyn Fft<f64>>,
    /// Forward of length N_φ.
    pub analysis: Arc<dyn Fft<f64>>,
}

/// rustfft's scratch for one thread: the plans' own scratch and a box for F → R's copy.
pub struct RustScratch {
    fft: Vec<c64>,
    work: Vec<c64>,
}

impl RustPhi {
    /// Plans from a fresh `FftPlanner` (it chooses the algorithm by length and CPU features).
    pub fn new(shape: Shape) -> Self {
        let mut planner = FftPlanner::<f64>::new();
        Self {
            shape,
            synthesis: planner.plan_fft_inverse(shape.azimuths),
            analysis: planner.plan_fft_forward(shape.azimuths),
        }
    }

    /// The live part of a box: 2P − 1 rings.
    fn live(&self) -> usize {
        self.shape.modes() * self.shape.azimuths
    }

    /// A per-thread scratch.
    pub fn scratch(&self) -> RustScratch {
        let len = self
            .synthesis
            .get_inplace_scratch_len()
            .max(self.analysis.get_inplace_scratch_len());
        RustScratch {
            fft: vec![c64::new(0.0, 0.0); len],
            work: vec![c64::new(0.0, 0.0); self.live()],
        }
    }

    /// S → F's φ step of one box, in place in `samples`: zero, scatter, transform.
    pub fn synthesize(&self, modes: &[c64], scratch: &mut RustScratch, samples: &mut [c64]) {
        let live = &mut samples[..self.live()];
        live.fill(c64::new(0.0, 0.0));
        self.shape.scatter(modes, live);
        self.synthesis.process_with_scratch(live, &mut scratch.fft);
    }

    /// F → R's φ step of one box: copy (the samples are kept), transform in place, gather.
    pub fn analyze(&self, samples: &[c64], scratch: &mut RustScratch, modes: &mut [c64]) {
        scratch.work.copy_from_slice(&samples[..self.live()]);
        self.analysis
            .process_with_scratch(&mut scratch.work, &mut scratch.fft);
        self.shape.gather(&scratch.work, modes);
    }
}

/// The table e^{2πi t/N}, t = 0..N, each entry from f64's `sin_cos` of the reduced angle.
pub fn twiddles(n: usize) -> Vec<c64> {
    (0..n)
        .map(|t| {
            let (s, c) = (2.0 * std::f64::consts::PI * t as f64 / n as f64).sin_cos();
            c64::new(c, s)
        })
        .collect()
}

/// The synthesis as a sum: samples of the live rings (the rest untouched).
pub fn naive_synthesis(shape: &Shape, modes: &[c64], samples: &mut [c64]) {
    let n = shape.azimuths;
    let w = twiddles(n);
    let m = shape.modes();
    for (ring, row) in modes.chunks_exact(m).enumerate() {
        for j in 0..n {
            let mut sum = c64::new(0.0, 0.0);
            for (k, &c) in row.iter().enumerate() {
                let t = (shape.slot(k) * j) % n;
                sum += c * w[t];
            }
            samples[ring * n + j] = sum;
        }
    }
}

/// The analysis as a sum: modes of the live rings.
pub fn naive_analysis(shape: &Shape, samples: &[c64], modes: &mut [c64]) {
    let n = shape.azimuths;
    let w = twiddles(n);
    let m = shape.modes();
    for (ring, row) in modes.chunks_exact_mut(m).enumerate() {
        for (k, out) in row.iter_mut().enumerate() {
            let mut sum = c64::new(0.0, 0.0);
            for j in 0..n {
                let t = (n - (shape.slot(k) * j) % n) % n;
                sum += samples[ring * n + j] * w[t];
            }
            *out = sum;
        }
    }
}

/// The relative ℓ2 distance of `a` from `b`.
pub fn relative_error(a: &[c64], b: &[c64]) -> f64 {
    let (mut num, mut den) = (0.0_f64, 0.0_f64);
    for (x, y) in a.iter().zip(b) {
        num += (x - y).norm_sqr();
        den += y.norm_sqr();
    }
    (num / den).sqrt()
}

/// A deterministic pseudo-random complex vector (splitmix64), entries in [−1/2, 1/2)².
pub fn random(len: usize, seed: u64) -> Vec<c64> {
    let mut state = seed.wrapping_mul(0x9e37_79b9_7f4a_7c15) ^ 0xd1b5_4a32_d192_ed03;
    let mut next = || {
        state = state.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = state;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^= z >> 31;
        (z >> 11) as f64 / (1u64 << 53) as f64 - 0.5
    };
    (0..len).map(|_| c64::new(next(), next())).collect()
}

/// FNV-1a over the bits of complex values.
pub fn hash(values: &[c64]) -> u64 {
    let mut h = 0xcbf2_9ce4_8422_2325_u64;
    for z in values {
        for x in [z.re, z.im] {
            for b in x.to_bits().to_le_bytes() {
                h ^= u64::from(b);
                h = h.wrapping_mul(0x0000_0100_0000_01b3);
            }
        }
    }
    h
}

/// Whether two slices have the same bits.
pub fn same_bits(a: &[c64], b: &[c64]) -> bool {
    a.len() == b.len()
        && a.iter()
            .zip(b)
            .all(|(x, y)| x.re.to_bits() == y.re.to_bits() && x.im.to_bits() == y.im.to_bits())
}

/// The product's S → F and F → R as shipped, on a batch of boxes at degree p (scale 1).
pub struct ProductBoxes {
    /// The conversions' tables.
    pub conversions: Conversions<f64>,
    /// The boxes' multipoles (also their locals).
    pub expansions: Vec<f64>,
    /// The boxes' signatures.
    pub signatures: Vec<c64>,
    /// One thread's scratch.
    pub scratch: ConversionScratch<f64>,
}

impl ProductBoxes {
    /// `boxes` boxes of random multipoles at degree p.
    pub fn new(p: usize, boxes: usize) -> Self {
        let conversions = Conversions::<f64>::new(p, 1.0);
        let scratch = ConversionScratch::new(&conversions);
        let len = conversions.expansion_len() * boxes;
        let expansions = random(len.div_ceil(2), 7)
            .iter()
            .flat_map(|z| [z.re, z.im])
            .take(len)
            .collect();
        let signatures = vec![c64::new(0.0, 0.0); conversions.samples() * boxes];
        Self {
            conversions,
            expansions,
            signatures,
            scratch,
        }
    }

    /// S → F of every box.
    pub fn s2f(&mut self) {
        self.conversions
            .s2f(&self.expansions, &mut self.signatures, &mut self.scratch);
    }

    /// F → R of every box (adds into the locals).
    pub fn f2r(&mut self) {
        self.conversions
            .f2r(&self.signatures, &mut self.expansions, &mut self.scratch);
    }

    /// The φ step's share of a conversion's multiply-adds.
    pub fn phi_share(&self) -> f64 {
        let (legendre, phi) = self.conversions.multiply_adds();
        phi as f64 / (legendre + phi) as f64
    }
}

/// Seconds per call of `f` on one thread: one warm-up call, then samples of at least
/// `min_seconds` each (at least 2 calls), the median of `samples`.
pub fn time_per_call(samples: usize, min_seconds: f64, mut f: impl FnMut()) -> f64 {
    f();
    let mut times: Vec<f64> = (0..samples)
        .map(|_| {
            let start = std::time::Instant::now();
            let mut calls = 0usize;
            while calls < 2 || start.elapsed().as_secs_f64() < min_seconds {
                f();
                calls += 1;
            }
            start.elapsed().as_secs_f64() / calls as f64
        })
        .collect();
    times.sort_by(f64::total_cmp);
    times[samples / 2]
}

/// The boxes of a timing batch: at least 2, at most what keeps a batch's samples near
/// 64 MiB.
pub fn timing_boxes(box_len: usize) -> usize {
    ((64 << 20) / (16 * box_len)).clamp(2, 16)
}

/// The machine line of a results file: the architecture, the SIMD features rustfft and FFTW
/// can use, `rlst::fftw::version()`.
pub fn machine_line() -> String {
    let host = std::env::var("HOSTNAME")
        .ok()
        .or_else(|| {
            std::process::Command::new("hostname")
                .output()
                .ok()
                .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string())
        })
        .unwrap_or_else(|| "unknown".into());
    format!(
        "host `{host}`, {}, features {}, FFTW `{}`, rustfft 6.4 (Cargo.lock)",
        std::env::consts::ARCH,
        simd_features(),
        rlst::fftw::version()
    )
}

/// The SIMD features detected at run time that rustfft dispatches on.
pub fn simd_features() -> String {
    let mut found: Vec<&str> = Vec::new();
    #[cfg(target_arch = "x86_64")]
    {
        if std::arch::is_x86_feature_detected!("sse4.1") {
            found.push("sse4.1");
        }
        if std::arch::is_x86_feature_detected!("avx") {
            found.push("avx");
        }
        if std::arch::is_x86_feature_detected!("fma") {
            found.push("fma");
        }
        if std::arch::is_x86_feature_detected!("avx512f") {
            found.push("avx512f");
        }
    }
    #[cfg(target_arch = "aarch64")]
    {
        if std::arch::is_aarch64_feature_detected!("neon") {
            found.push("neon");
        }
        if std::arch::is_aarch64_feature_detected!("sve") {
            found.push("sve");
        }
    }
    if found.is_empty() {
        "none".into()
    } else {
        found.join(", ")
    }
}
