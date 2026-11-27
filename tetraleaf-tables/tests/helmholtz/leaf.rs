//! The batched P2M of a leaf (Helmholtz CH4.9a2b, `tetraleaf_tables::helmholtz::leaf`) against
//! `tetraleaf_ref::helmholtz::leaf::p2m`, its trusted path.
//!
//! Error measure: per coefficient, |M̂ₙᵐ − M̂ₙᵐ_ref| over the contract of the module
//! documentation, C u √((2n + 1)/(4π)) Σⱼ |qⱼ eⱼ| (((n + 1)² + N) |ĵₙ(xⱼ)| + |xⱼ| |ĵₙ′(xⱼ)|),
//! xⱼ = kρⱼ, formed in f64 from the reference's radial functions, C = [`C`]; in f32 against the f64 reference on
//! the f32 inputs, u = 2⁻²⁴, plus the flushed values' bound. The worst ratio is printed.
//! Bit identities (vectors, real charges, accumulation) by `to_bits`.

use tetraleaf_math::radial::bessel_j;
use tetraleaf_math::{Complex, Layout, RealScalar};
use tetraleaf_ref::Frame;
use tetraleaf_ref::helmholtz::{Level, Workspace, leaf};
use tetraleaf_tables::helmholtz::leaf::{LANES, LeafCharges, LeafP2m, P2mScratch};

/// The contract's constant, about ten times the measured worst (root `CLAUDE.md`'s glibc
/// lesson).
const C: f64 = 10.0;

/// A deterministic generator of uniform values in [−1, 1).
struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((self.0 >> 11) as f64) / ((1u64 << 52) as f64) - 1.0
    }
}

/// A leaf of `n` points: the special points of the leaf frame first (the centre, the poles,
/// near-pole points, the axes, a corner), then uniform random points of [−1, 1)³.
fn points(n: usize, seed: u64) -> Vec<[f64; 3]> {
    let special = [
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 0.0, -0.7],
        [1e-9, 0.0, 0.8],
        [0.0, 1e-12, -0.5],
        [0.6, 0.0, 0.0],
        [0.0, -0.9, 0.0],
        [1.0, 1.0, 1.0],
        [-1.0, 0.5, 0.0],
    ];
    let mut rng = Lcg(seed);
    (0..n)
        .map(|j| {
            if j < special.len() && seed % 2 == 0 {
                special[j]
            } else {
                [rng.next(), rng.next(), rng.next()]
            }
        })
        .collect()
}

/// m vectors of charges, point-major, complex or real.
fn charges(n: usize, vectors: usize, complex: bool, seed: u64) -> Vec<f64> {
    let mut rng = Lcg(seed ^ 0x9e37_79b9);
    let c = if complex { 2 } else { 1 };
    (0..n * c * vectors).map(|_| rng.next()).collect()
}

/// The reference's multipoles: per vector one `leaf::p2m` call of every point (2 planes).
fn reference(
    p: usize,
    k: Complex<f64>,
    points: &[[f64; 3]],
    values: &[f64],
    vectors: usize,
    complex: bool,
) -> Vec<Vec<f64>> {
    let level = Level::new(k, 1.0);
    let frame = Frame::new([0.0; 3], 1.0);
    let mut ws = Workspace::new(p);
    let c = if complex { 2 } else { 1 };
    (0..vectors)
        .map(|v| {
            let q: Vec<Complex<f64>> = (0..points.len())
                .map(|j| {
                    let at = j * c * vectors + c * v;
                    Complex::new(values[at], if complex { values[at + 1] } else { 0.0 })
                })
                .collect();
            let mut m = vec![0.0; 2 * Layout::new(p).len()];
            leaf::p2m(p, &level, &frame, points, &q, &mut ws, &mut m);
            m
        })
        .collect()
}

/// Σⱼ |qⱼ eⱼ ĵₙ(kρⱼ)| and Σⱼ |qⱼ eⱼ| |kρⱼ| |ĵₙ′(kρⱼ)| (ĵₙ's condition at kρⱼ: the kernel's ρ and
/// the reference's differ in the last bit) per vector and degree, and Σⱼ |qⱼ eⱼ| per vector, in
/// f64.
fn weights(
    p: usize,
    k: Complex<f64>,
    points: &[[f64; 3]],
    values: &[f64],
    vectors: usize,
    complex: bool,
) -> (Vec<Vec<f64>>, Vec<Vec<f64>>, Vec<f64>) {
    let level = Level::new(k, 1.0);
    let c = if complex { 2 } else { 1 };
    let mut radial = vec![Complex::new(0.0, 0.0); p + 1];
    let mut derivative = vec![Complex::new(0.0, 0.0); p + 1];
    let mut sums = vec![vec![0.0; p + 1]; vectors];
    let mut conditioning = vec![vec![0.0; p + 1]; vectors];
    let mut charge = vec![0.0; vectors];
    for (j, x) in points.iter().enumerate() {
        let r = x[0].hypot(x[1]).hypot(x[2]);
        bessel_j(k * r, level.scale(), &mut radial, Some(&mut derivative));
        let e = (k.im * r - level.damping()).exp();
        for v in 0..vectors {
            let at = j * c * vectors + c * v;
            let q = Complex::new(values[at], if complex { values[at + 1] } else { 0.0 });
            let qe = q.norm() * e;
            charge[v] += qe;
            for n in 0..=p {
                sums[v][n] += qe * radial[n].norm();
                conditioning[v][n] += qe * (k * r).norm() * derivative[n].norm();
            }
        }
    }
    (sums, conditioning, charge)
}

/// The kernel at precision T on the inputs rounded to T, its multipoles in f64.
fn kernel<T: RealScalar>(
    p: usize,
    k: Complex<f64>,
    points: &[[f64; 3]],
    values: &[f64],
    vectors: usize,
    complex: bool,
) -> (Vec<f64>, usize) {
    let level = Level::new(Complex::new(T::from_f64(k.re), T::from_f64(k.im)), T::one());
    let pts: Vec<[T; 3]> = points.iter().map(|x| x.map(T::from_f64)).collect();
    let q: Vec<T> = values.iter().map(|&v| T::from_f64(v)).collect();
    let planes = LeafP2m::planes(&level, complex);
    let len = Layout::new(p).len();
    let mut out = vec![T::zero(); vectors * planes * len];
    let charges = LeafCharges {
        values: &q,
        vectors,
        complex,
    };
    LeafP2m::new(p).p2m(
        p,
        &level,
        &pts,
        charges,
        &mut P2mScratch::new(p, vectors),
        &mut out,
    );
    (out.iter().map(|&x| x.to_f64()).collect(), planes)
}

/// Rounds a leaf's inputs to f32 and back, so that the f64 reference sees the f32 kernel's.
fn to_f32(
    points: &[[f64; 3]],
    values: &[f64],
    k: Complex<f64>,
) -> (Vec<[f64; 3]>, Vec<f64>, Complex<f64>) {
    let round = |x: f64| x as f32 as f64;
    (
        points.iter().map(|x| x.map(round)).collect(),
        values.iter().map(|&v| round(v)).collect(),
        Complex::new(round(k.re), round(k.im)),
    )
}

/// The worst ratio of the error to the contract over every coefficient of one leaf.
#[allow(clippy::too_many_arguments)]
fn worst<T: RealScalar>(
    p: usize,
    k: Complex<f64>,
    points: &[[f64; 3]],
    values: &[f64],
    vectors: usize,
    complex: bool,
    u: f64,
    floor: (f64, f64),
) -> f64 {
    let (got, planes) = kernel::<T>(p, k, points, values, vectors, complex);
    let want = reference(p, k, points, values, vectors, complex);
    let (sums, conditioning, charge) = weights(p, k, points, values, vectors, complex);
    let layout = Layout::new(p);
    let len = layout.len();
    let mut ratio: f64 = 0.0;
    for v in 0..vectors {
        if planes == 1 {
            assert!(
                want[v][len..].iter().all(|&x| x == 0.0),
                "the reference's imaginary plane of real k and real charges is not zero"
            );
        }
        for plane in 0..planes {
            for i in 0..len {
                let (n, _) = layout.nm(i);
                let a = got[(v * planes + plane) * len + i];
                let b = want[v][plane * len + i];
                assert!(
                    a.is_finite(),
                    "p = {p}, k = {k}: coefficient {i} not finite"
                );
                let norm = ((2 * n + 1) as f64 / (4.0 * core::f64::consts::PI)).sqrt();
                let size = ((n + 1) * (n + 1) + points.len()) as f64;
                let bound = C * u * norm * (size * sums[v][n] + conditioning[v][n])
                    + floor.0 * sums[v][n]
                    + floor.1 * charge[v]
                    + f64::MIN_POSITIVE;
                let r = (a - b).abs() / bound * C;
                assert!(
                    (a - b).abs() <= bound,
                    "p = {p}, k = {k}, N = {}, vector {v}, plane {plane}, (n, i) = ({n}, {i}): \
                     {a:e} against {b:e}, bound {bound:e}",
                    points.len()
                );
                ratio = ratio.max(r);
            }
        }
    }
    ratio
}

const KS: [(f64, f64); 6] = [
    (0.3, 0.0),
    (4.0, 0.0),
    (16.0, 0.0),
    (4.0, 0.5),
    (16.0, 2.0),
    (1e-6, 0.0),
];

#[test]
fn p2m_within_the_contract_in_f64() {
    let mut most: f64 = 0.0;
    for (seed, &p) in [0usize, 1, 2, 5, 8, 17, 26, 45].iter().enumerate() {
        for &(re, im) in &KS {
            let k = Complex::new(re, im);
            for &n in &[1usize, 5, 16, 17, 40, 200] {
                if p >= 26 && n == 200 {
                    continue;
                }
                for (complex, vectors) in [(true, 1), (false, 1), (true, 2), (false, 3)] {
                    let s = (seed * 1000 + n) as u64;
                    let pts = points(n, s);
                    let q = charges(n, vectors, complex, s + 7);
                    let r = worst::<f64>(
                        p,
                        k,
                        &pts,
                        &q,
                        vectors,
                        complex,
                        f64::EPSILON / 2.0,
                        (0.0, 0.0),
                    );
                    most = most.max(r);
                }
            }
        }
    }
    println!("f64: the worst error over the contract without C: {most:.3}");
}

#[test]
fn p2m_within_the_contract_in_f32() {
    let mut most: f64 = 0.0;
    for (seed, &p) in [0usize, 2, 8, 17, 32].iter().enumerate() {
        for &(re, im) in &KS {
            for &n in &[1usize, 16, 40] {
                for (complex, vectors) in [(true, 1), (false, 2)] {
                    let s = (seed * 1000 + n + 1) as u64;
                    let (pts, q, k) = to_f32(
                        &points(n, s),
                        &charges(n, vectors, complex, s),
                        Complex::new(re, im),
                    );
                    // The flushed values: harmonics below 2⁻⁵⁷ of |w|, ĵₙ below 2⁻¹⁰⁰ times
                    // |Y| ≤ 1 of |q e| (module documentation).
                    let floor = (2f64.powi(-57), 2f64.powi(-100));
                    let r = worst::<f32>(
                        p,
                        k,
                        &pts,
                        &q,
                        vectors,
                        complex,
                        f32::EPSILON as f64 / 2.0,
                        floor,
                    );
                    most = most.max(r);
                    let (got, _) = kernel::<f32>(p, k, &pts, &q, vectors, complex);
                    assert!(
                        got.iter()
                            .all(|&x| x == 0.0 || x.abs() >= f32::MIN_POSITIVE as f64),
                        "p = {p}, k = {k}: a subnormal coefficient"
                    );
                }
            }
        }
    }
    println!("f32: the worst error over the contract without C: {most:.3}");
}

#[test]
fn each_vector_is_the_one_vector_kernels_bits() {
    let p = 12;
    for &(re, im) in &KS[..5] {
        let k = Complex::new(re, im);
        for complex in [true, false] {
            let n = 37;
            let vectors = 3;
            let pts = points(n, 3);
            let q = charges(n, vectors, complex, 5);
            let (all, planes) = kernel::<f64>(p, k, &pts, &q, vectors, complex);
            let c = if complex { 2 } else { 1 };
            let size = planes * Layout::new(p).len();
            for v in 0..vectors {
                let one: Vec<f64> = (0..n)
                    .flat_map(|j| q[j * c * vectors + c * v..][..c].to_vec())
                    .collect();
                let (alone, _) = kernel::<f64>(p, k, &pts, &one, 1, complex);
                let same = all[v * size..(v + 1) * size]
                    .iter()
                    .zip(&alone)
                    .all(|(a, b)| a.to_bits() == b.to_bits());
                assert!(
                    same,
                    "k = {k}, complex {complex}: vector {v} differs from its one-vector call"
                );
            }
        }
    }
}

#[test]
fn real_charges_are_complex_ones_real_plane_for_real_k() {
    let p = 15;
    for re in [0.3, 4.0, 16.0] {
        let k = Complex::new(re, 0.0);
        let n = 41;
        let pts = points(n, 11);
        let real = charges(n, 1, false, 13);
        let complex: Vec<f64> = real.iter().flat_map(|&q| [q, 0.0]).collect();
        let (a, planes_a) = kernel::<f64>(p, k, &pts, &real, 1, false);
        let (b, planes_b) = kernel::<f64>(p, k, &pts, &complex, 1, true);
        assert_eq!((planes_a, planes_b), (1, 2));
        let len = Layout::new(p).len();
        assert!(
            a.iter()
                .zip(&b[..len])
                .all(|(x, y)| x.to_bits() == y.to_bits()),
            "k = {k}"
        );
        assert!(
            b[len..].iter().all(|&x| x == 0.0),
            "k = {k}: an imaginary plane not zero"
        );
    }
}

#[test]
fn p2m_accumulates() {
    let p = 9;
    let level = Level::new(Complex::new(5.0, 0.5), 1.0);
    let len = Layout::new(p).len();
    let kernel = LeafP2m::<f64>::new(p);
    let mut scratch = P2mScratch::new(p, 1);
    for n in [LANES, 3 * LANES + 5] {
        let pts = points(n, 17);
        let q = charges(n, 1, true, 19);
        let charges = LeafCharges {
            values: &q,
            vectors: 1,
            complex: true,
        };
        let mut once = vec![0.0; 2 * len];
        kernel.p2m(p, &level, &pts, charges, &mut scratch, &mut once);
        let start: Vec<f64> = (0..2 * len).map(|i| (i as f64).sin()).collect();
        let mut added = start.clone();
        kernel.p2m(p, &level, &pts, charges, &mut scratch, &mut added);
        for i in 0..2 * len {
            let want = start[i] + once[i];
            assert!(
                (added[i] - want).abs() <= 1e-15 * (start[i].abs() + once[i].abs()) * n as f64,
                "N = {n}, {i}: {} against {want}",
                added[i]
            );
        }
        if n == LANES {
            // One batch: each coefficient one sum added, so twice is 2× exactly.
            let mut twice = once.clone();
            kernel.p2m(p, &level, &pts, charges, &mut scratch, &mut twice);
            assert!(
                twice
                    .iter()
                    .zip(&once)
                    .all(|(t, o)| t.to_bits() == (2.0 * o).to_bits())
            );
        }
    }
    // No points: nothing changes.
    let mut zero = vec![1.0; 2 * len];
    let none = LeafCharges {
        values: &[],
        vectors: 1,
        complex: true,
    };
    kernel.p2m(p, &level, &[], none, &mut scratch, &mut zero);
    assert!(zero.iter().all(|&x| x == 1.0));
}

#[test]
fn a_leafs_bits_do_not_depend_on_the_scratch() {
    // A scratch built for a larger degree and more vectors, and reused, gives the same bits.
    let p = 11;
    let level = Level::new(Complex::new(7.0, 0.0), 1.0);
    let pts = points(29, 23);
    let q = charges(29, 1, true, 29);
    let charges = LeafCharges {
        values: &q,
        vectors: 1,
        complex: true,
    };
    let len = Layout::new(p).len();
    let mut a = vec![0.0; 2 * len];
    LeafP2m::new(p).p2m(p, &level, &pts, charges, &mut P2mScratch::new(p, 1), &mut a);
    let mut big = P2mScratch::new(30, 4);
    let wide = LeafP2m::new(30);
    let other = points(50, 31);
    let other_q = self::charges(50, 1, true, 37);
    let mut junk = vec![0.0; 2 * Layout::new(30).len()];
    let other_charges = LeafCharges {
        values: &other_q,
        vectors: 1,
        complex: true,
    };
    wide.p2m(30, &level, &other, other_charges, &mut big, &mut junk);
    let mut b = vec![0.0; 2 * len];
    wide.p2m(p, &level, &pts, charges, &mut big, &mut b);
    assert!(a.iter().zip(&b).all(|(x, y)| x.to_bits() == y.to_bits()));
}

#[test]
#[should_panic(expected = "Yukawa")]
fn yukawa_is_refused() {
    let level = Level::new(Complex::new(0.0, 2.0), 1.0);
    let mut out = vec![0.0; 2 * 16];
    let charges = LeafCharges {
        values: &[1.0, 0.0],
        vectors: 1,
        complex: true,
    };
    LeafP2m::new(3).p2m(
        3,
        &level,
        &[[0.1, 0.2, 0.3]],
        charges,
        &mut P2mScratch::new(3, 1),
        &mut out,
    );
}

#[test]
#[should_panic(expected = "the multipoles must hold")]
fn a_short_multipole_is_refused() {
    let level = Level::new(Complex::new(2.0, 0.0), 1.0);
    let mut out = vec![0.0; 16];
    let charges = LeafCharges {
        values: &[1.0, 0.0],
        vectors: 1,
        complex: true,
    };
    LeafP2m::new(3).p2m(
        3,
        &level,
        &[[0.1, 0.2, 0.3]],
        charges,
        &mut P2mScratch::new(3, 1),
        &mut out,
    );
}

#[test]
#[should_panic(expected = "values per point")]
fn short_charges_are_refused() {
    let level = Level::new(Complex::new(2.0, 0.0), 1.0);
    let mut out = vec![0.0; 16];
    let charges = LeafCharges {
        values: &[1.0],
        vectors: 1,
        complex: false,
    };
    LeafP2m::new(3).p2m(
        3,
        &level,
        &[[0.1, 0.2, 0.3], [0.0; 3]],
        charges,
        &mut P2mScratch::new(3, 1),
        &mut out,
    );
}
