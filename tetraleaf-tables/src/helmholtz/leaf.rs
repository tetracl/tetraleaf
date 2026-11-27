//! P2M of a whole leaf by a batched kernel (Helmholtz CH4.9a2b; CONVENTIONS §3.15, "Radial
//! functions and scaling" and "Storage").
//!
//! [`tetraleaf_ref::helmholtz::leaf::p2m`] forms a box's stored multipole point by point:
//! per point the radial functions ĵₙ (`bessel_j`), the harmonics (`harmonics`), the factor
//! e^{Im(k)ρ − β}, and then the sum over the (p + 1)² terms,
//!
//! ```text
//! M̂ₙᵐ += Σⱼ qⱼ ĵₙ(kρⱼ) e^{Im(k)ρⱼ − β} Yₙᵐ(ρ̂ⱼ)
//! ```
//!
//! (its module documentation). [`LeafP2m`] forms the same sum for all the points of a leaf,
//! [`LANES`] at a time as lanes, on the leaf routines of `tetraleaf-math` (Helmholtz
//! CH4.9a2a): the harmonics of the batch column by column from tabulated coefficients
//! ([`LeafHarmonics::columns`]), and the radial functions in real arithmetic for real k
//! ([`bessel_j_real`]; complex k keeps [`bessel_j`], lane by lane). Per batch:
//!
//! 1. ρⱼ = √(x² + y² + z²) per lane (the reference takes `hypot`; the two agree to rounding);
//! 2. ĵₙ(kρⱼ) for n ≤ p, and eⱼ = e^{Im(k)ρⱼ − β} (for real k one value for the leaf);
//! 3. per charge vector and degree the weights wₙⱼ = (qⱼ eⱼ) ĵₙ(kρⱼ), complex in general,
//!    real for real k with real charges;
//! 4. per order m the column Q̄ₙᵐ of the batch (n = m … p) and (cos mφ, sin mφ), and per term
//!    the harmonics yⱼ = Q̄ₙᵐ cos mφ (and Q̄ₙᵐ sin mφ for −m) and the sums Σⱼ wₙⱼ yⱼ of each
//!    plane, the lanes added pairwise in a fixed order (eight pairs, then four, two, one), each
//!    sum added (+=) into the coefficient.
//!
//! The products in steps 3 and 4 are those of the reference in its order of operations
//! (q e, then times ĵₙ, then times Yₙᵐ = Q̄ₙᵐ cos mφ); only ρ, the routines' values and the
//! order of the sum over points differ. The last batch of a leaf is padded with the centre
//! and zero charges, whose terms are zeros. Every step is a loop over the lanes of fixed-size
//! arrays without a branch on the data, which the compiler vectorises without `unsafe`; the
//! start orders of the radial recursion, [`bessel_j`] of complex k and the factor's `exp`
//! run per lane.
//!
//! # Values and bits
//!
//! A leaf's multipole depends on its points and charges in their order alone, so it is the
//! same at every thread count and on any number of ranks that hold the leaf's points in the
//! same order. Each charge vector gets the arithmetic of the one-vector call (the vectors
//! share the routines' values, never a sum), so vector v of m is the one-vector kernel's bits.
//! For real k the real plane of complex charges takes the operations of real charges, so
//! complex charges with Im q = 0 give real charges' real plane bit for bit.
//!
//! Against the reference the values agree to a contract per term (`tests/helmholtz/leaf.rs`):
//! with N points and xⱼ = kρⱼ,
//!
//! ```text
//! |M̂ₙᵐ − M̂ₙᵐ_ref| ≤ C u √((2n + 1)/(4π)) Σⱼ |qⱼ eⱼ| (((n + 1)² + N) |ĵₙ(xⱼ)| + |xⱼ| |ĵₙ′(xⱼ)|)
//! ```
//!
//! the (n + 1)² from the harmonics' recurrence (CH4.9a2a), N from the two sums over points, and
//! the last term ĵₙ's condition at xⱼ: ρ by `sqrt` and the reference's by `hypot` differ in the
//! last bit, which near a zero of ĵₙ moves it by far more than u |ĵₙ|. The tests hold C = 10;
//! the measured worst is about 1 in f64 and f32.
//! In f32 the routines flush under-range values (the harmonics' columns below 2⁻¹⁰⁰, ĵₙ below
//! 2⁻¹⁰⁰ beyond the turning point; [`tetraleaf_math::orthonormal::leaf`], "Under-range
//! values"), which moves a term by at most about 2⁻⁵⁷ |wₙⱼ|, far below u.
//!
//! # Example
//!
//! ```
//! use tetraleaf_math::{Complex, Layout};
//! use tetraleaf_ref::helmholtz::{Level, Workspace, leaf};
//! use tetraleaf_ref::Frame;
//! use tetraleaf_tables::helmholtz::leaf::{LeafCharges, LeafP2m, P2mScratch};
//!
//! let p = 10;
//! let level = Level::new(Complex::new(3.0, 0.0), 1.0);
//! let points = [[0.3, -0.2, 0.5], [-0.7, 0.1, 0.9], [0.0, 0.0, 0.0]];
//! // Complex charges, point-major: Re q, Im q per point.
//! let values = [1.0, 0.5, -0.25, 0.0, 2.0, -1.0];
//! let charges = LeafCharges { values: &values, vectors: 1, complex: true };
//! let kernel = LeafP2m::<f64>::new(p);
//! let mut scratch = P2mScratch::new(p, 1);
//! let len = Layout::new(p).len();
//! let mut multipole = vec![0.0; 2 * len];
//! kernel.p2m(p, &level, &points, charges, &mut scratch, &mut multipole);
//!
//! let q: Vec<Complex<f64>> = values.chunks(2).map(|c| Complex::new(c[0], c[1])).collect();
//! let mut reference = vec![0.0; 2 * len];
//! let frame = Frame::new([0.0; 3], 1.0);
//! leaf::p2m(p, &level, &frame, &points, &q, &mut Workspace::new(p), &mut reference);
//! for (a, b) in multipole.iter().zip(&reference) {
//!     assert!((a - b).abs() < 1e-13);
//! }
//! ```

use core::array;

use tetraleaf_math::orthonormal::{LeafHarmonics, LeafScratch};
use tetraleaf_math::radial::{bessel_j, bessel_j_real};
use tetraleaf_math::{Complex, Layout, RealScalar};
use tetraleaf_ref::helmholtz::Level;

/// The points per batch: 16, the fastest of 4, 8 and 16 lanes for the routines on Kathleen
/// and locust, or within 5 % (Helmholtz CH4.9a2a, its report §5).
pub const LANES: usize = 16;

const B: usize = LANES;

/// The charges of a leaf's points for [`LeafP2m::p2m`]: point-major, point j's `vectors`
/// vectors one after the other, each one real value or two (Re q, Im q) if `complex`, as the
/// source chunks of CONVENTIONS §3.13 and `tetraleaf-exec`'s Helmholtz densities hold them:
/// vector v's real part of point j at j c m + c v, its imaginary part after it, c = 1 or 2.
#[derive(Clone, Copy, Debug)]
pub struct LeafCharges<'a, T> {
    /// The values, c m per point.
    pub values: &'a [T],
    /// The number m ≥ 1 of charge vectors.
    pub vectors: usize,
    /// Whether each charge is complex (two values) or real (one).
    pub complex: bool,
}

impl<T> LeafCharges<'_, T> {
    /// The values per point, c m.
    pub fn columns(&self) -> usize {
        if self.complex {
            2 * self.vectors
        } else {
            self.vectors
        }
    }
}

/// The batched P2M of a leaf up to a degree (the [module documentation](self)): the leaf
/// harmonics' table. It allocates once, in [`LeafP2m::new`], and is read-only after, so one can
/// serve every thread.
#[derive(Clone, Debug)]
pub struct LeafP2m<T> {
    harmonics: LeafHarmonics<T>,
}

impl<T: RealScalar> LeafP2m<T> {
    /// The kernel for degrees up to `p`.
    pub fn new(p: usize) -> Self {
        Self {
            harmonics: LeafHarmonics::new(p),
        }
    }

    /// The largest degree.
    pub fn p(&self) -> usize {
        self.harmonics.p()
    }

    /// The planes of a multipole of `level` with these charges: 1 for real k with real
    /// charges (a real multipole, CONVENTIONS §3.15's storage), 2 otherwise.
    pub fn planes(level: &Level<T>, complex: bool) -> usize {
        if level.k().im == T::zero() && !complex {
            1
        } else {
            2
        }
    }

    /// P2M of a leaf: adds the stored multipoles of degree `p` of the points `points` (offsets
    /// from the leaf's centre in the units of `level`, whose half-width is the leaf's) with the
    /// charges `charges` to `multipoles`, for every charge vector:
    /// M̂ₙᵐ += Σⱼ qⱼ ĵₙ(kρⱼ) e^{Im(k)ρⱼ − β} Yₙᵐ(ρ̂ⱼ), [`planes`](Self::planes) (p + 1)²
    /// values per vector in the layout of §3.6, the real plane then the imaginary plane, the
    /// vectors one after the other. The sum of
    /// [`tetraleaf_ref::helmholtz::leaf::p2m`] over the points, to the contract of the [module
    /// documentation](self).
    ///
    /// # Panics
    ///
    /// If `p` exceeds the kernel's or the scratch's degree, the charges have more vectors than
    /// the scratch, the level is Yukawa's (k = iκ: H2's path), `charges.values` does not hold
    /// [`columns`](LeafCharges::columns) values per point, or `multipoles` does not have
    /// m [`planes`](Self::planes) (p + 1)² values.
    pub fn p2m(
        &self,
        p: usize,
        level: &Level<T>,
        points: &[[T; 3]],
        charges: LeafCharges<'_, T>,
        scratch: &mut P2mScratch<T>,
        multipoles: &mut [T],
    ) {
        assert!(p <= self.p(), "degree {p} above the kernel's {}", self.p());
        assert!(
            p <= scratch.p,
            "degree {p} above the scratch's {}",
            scratch.p
        );
        assert!(
            charges.vectors >= 1 && charges.vectors <= scratch.vectors,
            "{} charge vectors, the scratch holds 1 to {}",
            charges.vectors,
            scratch.vectors
        );
        assert!(
            !level.is_yukawa(),
            "the batched P2M takes Helmholtz's complex path, not Yukawa's"
        );
        let columns = charges.columns();
        assert_eq!(
            charges.values.len(),
            columns * points.len(),
            "the charges must hold {columns} values per point"
        );
        let len = Layout::new(p).len();
        let planes = Self::planes(level, charges.complex);
        assert_eq!(
            multipoles.len(),
            charges.vectors * planes * len,
            "the multipoles must hold {} values",
            charges.vectors * planes * len
        );
        let k = level.k();
        let real_k = k.im == T::zero();
        let (s, beta) = (level.scale(), level.damping());
        // For real k the factor is e^{0·ρ − β}, the same for every point.
        let e_real = (T::zero() - beta).exp();
        let P2mScratch {
            harmonics: harmonic_scratch,
            radial,
            radial_im,
            one,
            weights,
            ..
        } = scratch;
        let degrees = p + 1;
        for start in (0..points.len()).step_by(B) {
            let count = B.min(points.len() - start);
            let lanes: [[T; 3]; B] = array::from_fn(|j| {
                if j < count {
                    points[start + j]
                } else {
                    [T::zero(); 3]
                }
            });
            let rho: [T; B] = array::from_fn(|j| {
                (lanes[j][0] * lanes[j][0] + lanes[j][1] * lanes[j][1] + lanes[j][2] * lanes[j][2])
                    .sqrt()
            });
            // The radial functions and the factor.
            let e: [T; B] = if real_k {
                let x: [T; B] = array::from_fn(|j| k.re * rho[j]);
                bessel_j_real(&x, s, &mut radial[..degrees]);
                [e_real; B]
            } else {
                for j in 0..B {
                    bessel_j(k * rho[j], s, &mut one[..degrees], None);
                    for (n, value) in one[..degrees].iter().enumerate() {
                        radial[n][j] = value.re;
                        radial_im[n][j] = value.im;
                    }
                }
                array::from_fn(|j| (k.im * rho[j] - beta).exp())
            };
            // The weights per vector and degree: plane 0 Re, plane 1 Im.
            for v in 0..charges.vectors {
                let charge = |j: usize, part: usize| {
                    if j < count {
                        charges.values
                            [(start + j) * columns + v * (columns / charges.vectors) + part]
                    } else {
                        T::zero()
                    }
                };
                let w_re: [T; B] = array::from_fn(|j| charge(j, 0) * e[j]);
                let w_im: [T; B] = if charges.complex {
                    array::from_fn(|j| charge(j, 1) * e[j])
                } else {
                    [T::zero(); B]
                };
                let out = &mut weights[2 * v * degrees..2 * (v + 1) * degrees];
                let (re, im) = out.split_at_mut(degrees);
                for n in 0..degrees {
                    let j_re = &radial[n];
                    if real_k {
                        re[n] = array::from_fn(|j| w_re[j] * j_re[j]);
                        if planes == 2 {
                            im[n] = array::from_fn(|j| w_im[j] * j_re[j]);
                        }
                    } else {
                        // (q e) ĵₙ as the reference's complex product.
                        let j_im = &radial_im[n];
                        re[n] = array::from_fn(|j| w_re[j] * j_re[j] - w_im[j] * j_im[j]);
                        im[n] = array::from_fn(|j| w_re[j] * j_im[j] + w_im[j] * j_re[j]);
                    }
                }
            }
            let weights = &*weights;
            self.harmonics
                .columns(p, &lanes, harmonic_scratch, |m, column, cos, sin| {
                    for (n, q) in (m..=p).zip(column) {
                        let centre = n * n + n;
                        let y_cos: [T; B] = if m == 0 {
                            *q
                        } else {
                            array::from_fn(|j| q[j] * cos[j])
                        };
                        let y_sin: [T; B] = array::from_fn(|j| q[j] * sin[j]);
                        for v in 0..charges.vectors {
                            let w = &weights[2 * v * degrees..2 * (v + 1) * degrees];
                            let multipole =
                                &mut multipoles[v * planes * len..(v + 1) * planes * len];
                            for plane in 0..planes {
                                let wn = &w[plane * degrees + n];
                                let coefficients = &mut multipole[plane * len..(plane + 1) * len];
                                let i = centre + m;
                                coefficients[i] = coefficients[i] + dot(wn, &y_cos);
                                if m > 0 {
                                    let i = centre - m;
                                    coefficients[i] = coefficients[i] + dot(wn, &y_sin);
                                }
                            }
                        }
                    }
                });
        }
    }
}

/// Σⱼ aⱼ bⱼ over the lanes, the products added pairwise in a fixed order: lane j with j + 8,
/// then j with j + 4, j + 2, j + 1.
#[inline(always)]
fn dot<T: RealScalar>(a: &[T; B], b: &[T; B]) -> T {
    let eight: [T; 8] = array::from_fn(|j| a[j] * b[j] + a[j + 8] * b[j + 8]);
    let four: [T; 4] = array::from_fn(|j| eight[j] + eight[j + 4]);
    let two: [T; 2] = array::from_fn(|j| four[j] + four[j + 2]);
    two[0] + two[1]
}

/// The work space of [`LeafP2m::p2m`] up to a degree and a number of charge vectors: the
/// harmonics' column, the radial functions and the weights of one batch. It allocates once,
/// in [`P2mScratch::new`]; its size grows with p and m, not with the leaf's points.
#[derive(Clone, Debug)]
pub struct P2mScratch<T> {
    p: usize,
    vectors: usize,
    harmonics: LeafScratch<T, B>,
    /// ĵₙ per lane (its real part for complex k), n ≤ p.
    radial: Vec<[T; B]>,
    /// Im ĵₙ per lane for complex k.
    radial_im: Vec<[T; B]>,
    /// One lane's ĵₙ for complex k.
    one: Vec<Complex<T>>,
    /// Per vector the weights' real then imaginary planes, p + 1 rows each.
    weights: Vec<[T; B]>,
}

impl<T: RealScalar> P2mScratch<T> {
    /// Work space for degrees up to `p` and up to `vectors` charge vectors.
    pub fn new(p: usize, vectors: usize) -> Self {
        let degrees = p + 1;
        Self {
            p,
            vectors,
            harmonics: LeafScratch::new(p),
            radial: vec![[T::zero(); B]; degrees],
            radial_im: vec![[T::zero(); B]; degrees],
            one: vec![Complex::new(T::zero(), T::zero()); degrees],
            weights: vec![[T::zero(); B]; 2 * vectors * degrees],
        }
    }

    /// The largest degree.
    pub fn p(&self) -> usize {
        self.p
    }

    /// The most charge vectors.
    pub fn vectors(&self) -> usize {
        self.vectors
    }
}
