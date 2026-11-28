//! The batched P2M of a leaf on Yukawa's real path (Helmholtz CH2.1c,
//! `tetraleaf_tables::helmholtz::leaf`, "Yukawa") against
//! `tetraleaf_ref::helmholtz::leaf::p2m_yukawa`, its trusted path.
//!
//! Error measure: per coefficient, |M̂ᵞₙᵐ − M̂ᵞₙᵐ_ref| over the contract of the module documentation,
//! C u √((2n + 1)/(4π)) Σⱼ |qⱼ eⱼ| (((n + 1)² + N) îₙ(xⱼ) + xⱼ (|îₙ′(xⱼ)| + îₙ(xⱼ))), xⱼ = κρⱼ and
//! eⱼ = e^{κρⱼ − β}, formed in f64 from the reference's radial functions, C = [`C`]; in f32 against
//! the f64 reference on the f32 inputs, u = 2⁻²⁴, plus the flushed values' bound and the f32 level's
//! own rounding (s and β formed in f32: [`LEVEL_C`] (β + n) u √((2n + 1)/(4π)) Σⱼ |qⱼ eⱼ| îₙ, beside
//! the contract and not multiplied by C). The worst ratio
//! is printed. Bit identities (vectors, complex charges as two real columns, accumulation) by
//! `to_bits`.

use tetraleaf_math::radial::modified_i;
use tetraleaf_math::{Complex, Layout, RealScalar};
use tetraleaf_ref::Frame;
use tetraleaf_ref::helmholtz::{Level, Workspace, leaf};
use tetraleaf_tables::helmholtz::leaf::{LANES, LeafCharges, LeafP2m, P2mScratch};

use super::leaf::{C, Lcg, points};

/// κ̂ of the leaf frame (half-width 1): from the Laplace limit to strong damping (κ̂ a ≈ 104), and
/// κ̂ = 250 (κ̂ a ≈ 433, a shallow leaf at κD = 10³) in f64.
const KAPPAS: [f64; 5] = [1e-6, 0.3, 4.0, 16.0, 60.0];
const KAPPA_LARGE: f64 = 250.0;

/// The f32 level's own rounding (its s and β formed in f32 against the f64 reference's level):
/// e^{−β} and 1/sⁿ relative (β + n) u, held apart from the contract's C (module documentation).
const LEVEL_C: f64 = 2.0;

/// m vectors of charges, point-major, complex or real.
fn charges(n: usize, vectors: usize, complex: bool, seed: u64) -> Vec<f64> {
    let mut rng = Lcg(seed ^ 0x7f4a_7c15);
    let c = if complex { 2 } else { 1 };
    (0..n * c * vectors).map(|_| rng.next()).collect()
}

/// The real charge column (vector v, part 0 = Re or 1 = Im) of point-major charges.
fn column(
    values: &[f64],
    n: usize,
    vectors: usize,
    complex: bool,
    v: usize,
    part: usize,
) -> Vec<f64> {
    let c = if complex { 2 } else { 1 };
    (0..n)
        .map(|j| values[j * c * vectors + c * v + part])
        .collect()
}

/// The reference's multipoles: per vector and real column one `leaf::p2m_yukawa` call of every
/// point, the columns' planes one after the other (Re q then Im q for complex charges).
fn reference(
    p: usize,
    kappa: f64,
    points: &[[f64; 3]],
    values: &[f64],
    vectors: usize,
    complex: bool,
) -> Vec<f64> {
    let level = Level::new(Complex::new(0.0, kappa), 1.0);
    let frame = Frame::new([0.0; 3], 1.0);
    let mut ws = Workspace::new(p);
    let len = Layout::new(p).len();
    let parts = if complex { 2 } else { 1 };
    let mut out = vec![0.0; vectors * parts * len];
    for v in 0..vectors {
        for part in 0..parts {
            let q = column(values, points.len(), vectors, complex, v, part);
            let at = (v * parts + part) * len;
            leaf::p2m_yukawa(
                p,
                &level,
                &frame,
                points,
                &q,
                &mut ws,
                &mut out[at..at + len],
            );
        }
    }
    out
}

/// Per real column (vector, part) and degree Σⱼ |qⱼ eⱼ| îₙ(xⱼ) and Σⱼ |qⱼ eⱼ| xⱼ (|îₙ′(xⱼ)| + îₙ(xⱼ)),
/// and per column Σⱼ |qⱼ eⱼ| and Σⱼ |qⱼ|, in f64.
#[allow(clippy::type_complexity)]
fn weights(
    p: usize,
    kappa: f64,
    points: &[[f64; 3]],
    values: &[f64],
    vectors: usize,
    complex: bool,
) -> (Vec<Vec<f64>>, Vec<Vec<f64>>, Vec<f64>, Vec<f64>) {
    let level = Level::new(Complex::new(0.0, kappa), 1.0);
    let parts = if complex { 2 } else { 1 };
    let columns = vectors * parts;
    let mut radial = vec![0.0; p + 1];
    let mut derivative = vec![0.0; p + 1];
    let mut sums = vec![vec![0.0; p + 1]; columns];
    let mut conditioning = vec![vec![0.0; p + 1]; columns];
    let mut charge = vec![0.0; columns];
    let mut plain = vec![0.0; columns];
    for (j, x) in points.iter().enumerate() {
        let r = x[0].hypot(x[1]).hypot(x[2]);
        let xj = kappa * r;
        modified_i(xj, level.scale(), &mut radial, Some(&mut derivative));
        let e = (xj - level.damping()).exp();
        for v in 0..vectors {
            for part in 0..parts {
                let o = v * parts + part;
                let c = if complex { 2 } else { 1 };
                let q = values[j * c * vectors + c * v + part].abs();
                let qe = q * e;
                charge[o] += qe;
                plain[o] += q;
                for n in 0..=p {
                    sums[o][n] += qe * radial[n].abs();
                    conditioning[o][n] += qe * xj * (derivative[n].abs() + radial[n].abs());
                }
            }
        }
    }
    (sums, conditioning, charge, plain)
}

/// The kernel at precision T on the inputs rounded to T, its multipoles in f64.
fn kernel<T: RealScalar>(
    p: usize,
    kappa: f64,
    points: &[[f64; 3]],
    values: &[f64],
    vectors: usize,
    complex: bool,
) -> (Vec<f64>, usize) {
    let level = Level::new(Complex::new(T::zero(), T::from_f64(kappa)), T::one());
    assert!(level.is_yukawa());
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

/// The worst ratio of the error to the contract over every coefficient of one leaf.
#[allow(clippy::too_many_arguments)]
fn worst<T: RealScalar>(
    p: usize,
    kappa: f64,
    points: &[[f64; 3]],
    values: &[f64],
    vectors: usize,
    complex: bool,
    u: f64,
    floor: (f64, f64),
    level_rounding: f64,
) -> f64 {
    let (got, planes) = kernel::<T>(p, kappa, points, values, vectors, complex);
    assert_eq!(planes, if complex { 2 } else { 1 });
    let want = reference(p, kappa, points, values, vectors, complex);
    assert_eq!(got.len(), want.len());
    let (sums, conditioning, charge, plain) = weights(p, kappa, points, values, vectors, complex);
    let layout = Layout::new(p);
    let len = layout.len();
    let beta = Level::new(Complex::new(0.0, kappa), 1.0).damping();
    let mut ratio: f64 = 0.0;
    for o in 0..vectors * planes {
        for i in 0..len {
            let (n, _) = layout.nm(i);
            let a = got[o * len + i];
            let b = want[o * len + i];
            assert!(
                a.is_finite(),
                "p = {p}, κ = {kappa}: coefficient {i} not finite"
            );
            let norm = ((2 * n + 1) as f64 / (4.0 * core::f64::consts::PI)).sqrt();
            let size = ((n + 1) * (n + 1) + points.len()) as f64;
            // The level's s and β formed in T against the f64 reference's: e^{−β} and 1/sⁿ
            // relative (β + n) u_T (f32 only; in f64 the two share the level).
            let level = LEVEL_C * level_rounding * u * norm * (beta + n as f64) * sums[o][n];
            let bound = C * u * norm * (size * sums[o][n] + conditioning[o][n])
                + level
                + floor.0 * sums[o][n]
                + floor.1 * norm * (charge[o] + plain[o])
                + f64::MIN_POSITIVE;
            assert!(
                (a - b).abs() <= bound,
                "p = {p}, κ = {kappa}, N = {}, column {o}, (n, i) = ({n}, {i}): {a:e} against \
                 {b:e}, bound {bound:e}",
                points.len()
            );
            ratio = ratio.max((a - b).abs() / bound * C);
        }
    }
    ratio
}

#[test]
fn yukawa_p2m_within_the_contract_in_f64() {
    let mut most: f64 = 0.0;
    for (seed, &p) in [0usize, 1, 2, 5, 8, 17, 26, 45].iter().enumerate() {
        for &kappa in KAPPAS.iter().chain([&KAPPA_LARGE]) {
            for &n in &[1usize, 5, 16, 17, 40, 200] {
                if p >= 26 && n == 200 {
                    continue;
                }
                for (complex, vectors) in [(false, 1), (true, 1), (false, 3), (true, 3)] {
                    let s = (seed * 1000 + n) as u64;
                    let pts = points(n, s);
                    let q = charges(n, vectors, complex, s + 7);
                    let r = worst::<f64>(
                        p,
                        kappa,
                        &pts,
                        &q,
                        vectors,
                        complex,
                        f64::EPSILON / 2.0,
                        (0.0, 0.0),
                        0.0,
                    );
                    most = most.max(r);
                }
            }
        }
    }
    println!("Yukawa, f64: the worst error over the contract without C: {most:.3}");
}

#[test]
fn yukawa_p2m_within_the_contract_in_f32() {
    let round = |x: f64| x as f32 as f64;
    let mut most: f64 = 0.0;
    for (seed, &p) in [0usize, 2, 8, 17, 32].iter().enumerate() {
        for &kappa in &KAPPAS {
            let kappa = round(kappa);
            for &n in &[1usize, 16, 40] {
                for (complex, vectors) in [(false, 1), (true, 2)] {
                    let s = (seed * 1000 + n + 1) as u64;
                    let pts: Vec<[f64; 3]> = points(n, s).iter().map(|x| x.map(round)).collect();
                    let q: Vec<f64> = charges(n, vectors, complex, s)
                        .iter()
                        .map(|&v| round(v))
                        .collect();
                    // The flushed values: harmonics below 2⁻⁵⁷ of |w|; îₙ below 2⁻¹⁰⁰ times
                    // |Y| ≤ √((2n + 1)/(4π)) of |q e|; and e below 2⁻¹⁰⁰ (strong damping) times
                    // îₙ ≤ 1 and |Y| of |q| (module documentation, "Yukawa").
                    let floor = (2f64.powi(-57), 2f64.powi(-100));
                    let r = worst::<f32>(
                        p,
                        kappa,
                        &pts,
                        &q,
                        vectors,
                        complex,
                        f32::EPSILON as f64 / 2.0,
                        floor,
                        1.0,
                    );
                    most = most.max(r);
                    let (got, _) = kernel::<f32>(p, kappa, &pts, &q, vectors, complex);
                    assert!(
                        got.iter()
                            .all(|&x| x == 0.0 || x.abs() >= f32::MIN_POSITIVE as f64),
                        "p = {p}, κ = {kappa}: a subnormal coefficient"
                    );
                }
            }
        }
    }
    println!("Yukawa, f32: the worst error over the contract without C: {most:.3}");
}

#[test]
fn yukawa_each_vector_is_the_one_vector_kernels_bits() {
    let p = 12;
    for &kappa in &KAPPAS {
        for complex in [false, true] {
            let n = 37;
            let vectors = 3;
            let pts = points(n, 3);
            let q = charges(n, vectors, complex, 5);
            let (all, planes) = kernel::<f64>(p, kappa, &pts, &q, vectors, complex);
            let c = if complex { 2 } else { 1 };
            let size = planes * Layout::new(p).len();
            for v in 0..vectors {
                let one: Vec<f64> = (0..n)
                    .flat_map(|j| q[j * c * vectors + c * v..][..c].to_vec())
                    .collect();
                let (alone, _) = kernel::<f64>(p, kappa, &pts, &one, 1, complex);
                assert!(
                    all[v * size..(v + 1) * size]
                        .iter()
                        .zip(&alone)
                        .all(|(a, b)| a.to_bits() == b.to_bits()),
                    "κ = {kappa}, complex {complex}: vector {v} differs from its one-vector call"
                );
            }
        }
    }
}

#[test]
fn yukawa_complex_charges_are_two_real_columns() {
    // A complex vector's planes are the real-charge kernel on Re q and on Im q, bit for bit
    // (`tetraleaf-exec`'s densities: a complex vector on Yukawa is two real vectors).
    let p = 15;
    for &kappa in &KAPPAS {
        let n = 41;
        let pts = points(n, 11);
        let q = charges(n, 2, true, 13);
        let (both, planes) = kernel::<f64>(p, kappa, &pts, &q, 2, true);
        assert_eq!(planes, 2);
        let len = Layout::new(p).len();
        for v in 0..2 {
            for part in 0..2 {
                let real = column(&q, n, 2, true, v, part);
                let (alone, planes) = kernel::<f64>(p, kappa, &pts, &real, 1, false);
                assert_eq!(planes, 1);
                let at = (2 * v + part) * len;
                assert!(
                    both[at..at + len]
                        .iter()
                        .zip(&alone)
                        .all(|(a, b)| a.to_bits() == b.to_bits()),
                    "κ = {kappa}: vector {v}, part {part}"
                );
            }
        }
    }
}

#[test]
fn yukawa_p2m_accumulates() {
    let p = 9;
    let level = Level::new(Complex::new(0.0, 3.0), 1.0);
    let len = Layout::new(p).len();
    let kernel = LeafP2m::<f64>::new(p);
    let mut scratch = P2mScratch::new(p, 1);
    let pts = points(LANES, 17);
    let q = charges(LANES, 1, false, 19);
    let charges = LeafCharges {
        values: &q,
        vectors: 1,
        complex: false,
    };
    let mut once = vec![0.0; len];
    kernel.p2m(p, &level, &pts, charges, &mut scratch, &mut once);
    // One batch: each coefficient one sum added, so twice is 2× exactly.
    let mut twice = once.clone();
    kernel.p2m(p, &level, &pts, charges, &mut scratch, &mut twice);
    assert!(
        twice
            .iter()
            .zip(&once)
            .all(|(t, o)| t.to_bits() == (2.0 * o).to_bits())
    );
    assert!(once.iter().any(|&x| x != 0.0));
}

#[test]
#[should_panic(expected = "values per point")]
fn yukawa_short_charges_are_refused() {
    let level = Level::new(Complex::new(0.0, 2.0), 1.0);
    let mut out = vec![0.0; 2 * 16];
    let charges = LeafCharges {
        values: &[1.0, 0.5, 2.0],
        vectors: 1,
        complex: true,
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

#[test]
#[should_panic(expected = "the multipoles must hold")]
fn yukawa_two_planes_for_real_charges_are_refused() {
    let level = Level::new(Complex::new(0.0, 2.0), 1.0);
    let mut out = vec![0.0; 2 * 16];
    let charges = LeafCharges {
        values: &[1.0],
        vectors: 1,
        complex: false,
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
