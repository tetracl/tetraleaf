//! The batched L2P of a leaf on Yukawa's real path (Helmholtz CH2.1e,
//! `tetraleaf_tables::helmholtz::leaf::LeafL2p`, "L2P") against
//! `tetraleaf_ref::helmholtz::leaf::l2p_yukawa`, its trusted path.
//!
//! Error measure: per target, real local and output (φ̃, each gradient component), the absolute
//! difference over the contract of the module documentation ("L2P", Yukawa), formed in f64 from the
//! reference's radial functions at the reference's r (`hypot`): with wₙ = √((2n + 1)/(4π)),
//! Sₙ = Σₘ |L̂ᵞₙᵐ|, x = κr, tₙ = îₙ/r (κ/(3s) for n = 1 at r = 0) and e = e^{κr − β},
//! C u e Σₙ wₙ Sₙ (((n + 1)² + 2n + 1) îₙ + x (|îₙ′| + îₙ)) for φ̃ and
//! C u e Σₙ wₙ Sₙ (κ (((n + 1)² + 2n + 3) |îₙ′| + x îₙ + x |îₙ′|) + ((n + 1)³ + n(n + 1) + x) tₙ) for
//! each gradient component, C = [`C`]. In f32 against the f64 reference on the f32 inputs,
//! u = 2⁻²⁴, beside the contract (not multiplied by C): the f32 level's own rounding of s and β
//! ([`LEVEL_C`] (β + n + 1) u e Σₙ wₙ Sₙ times îₙ, or κ |îₙ′| + tₙ), the routines' flushed harmonics
//! (2⁻⁵⁰ per degree's sum) and the flushed îₙ, îₙ′ and e (τ Σₙ wₙ Sₙ (1 + κ + 1/r)). The worst ratio
//! is printed. Bit identities (vectors, complex charges' two locals, positions, potentials with and
//! without gradients, accumulation) by `to_bits`.

use core::f64::consts::PI;

use tetraleaf_math::radial::modified_i;
use tetraleaf_math::{Complex, Layout, RealScalar};
use tetraleaf_ref::Frame;
use tetraleaf_ref::helmholtz::{Level, Workspace, leaf};
use tetraleaf_tables::helmholtz::leaf::{L2pScratch, LANES, LeafL2p, LeafOutput};

use super::leaf::{C, Lcg, points};

/// κ̂ of the leaf frame (half-width 1), CH2.1c's: from the Laplace limit to strong damping
/// (κ̂ a ≈ 104), and κ̂ = 250 (κ̂ a ≈ 433, a shallow leaf at κD = 10³) in f64.
const KAPPAS: [f64; 5] = [1e-6, 0.3, 4.0, 16.0, 60.0];
const KAPPA_LARGE: f64 = 250.0;

/// The f32 level's own rounding (its s and β formed in f32 against the f64 reference's level):
/// e^{−β} and 1/sⁿ relative (β + n) u, held apart from the contract's C (CH2.1c's constant).
const LEVEL_C: f64 = 2.0;

/// The under-range threshold below which the kernel flushes îₙ, îₙ′ and e.
fn tau<T: RealScalar>() -> f64 {
    if T::max_value().to_f64() > 1e300 {
        2f64.powi(-900)
    } else {
        2f64.powi(-100)
    }
}

/// c m real locals of degree p (c = `planes`), uniform in [−1, 1).
fn locals(p: usize, vectors: usize, planes: usize, seed: u64) -> Vec<f64> {
    let mut rng = Lcg(seed ^ 0x2c1b_3c6d);
    (0..planes * vectors * Layout::new(p).len())
        .map(|_| rng.next())
        .collect()
}

/// The targets of a leaf of n: `points`' (the centre, the poles, near-pole points, the axes and a
/// corner for an even seed), with three points next to the z-axis in place of the last three random
/// ones when n ≥ 12, as `l2p.rs`' targets.
fn targets(n: usize, seed: u64) -> Vec<[f64; 3]> {
    let mut x = points(n, seed);
    if n >= 12 {
        x[n - 3] = [1e-160, 0.0, 1.0];
        x[n - 2] = [-1e-160, 1e-160, -0.5];
        x[n - 1] = [-1e-40, 1e-40, -0.5];
    }
    x
}

/// One leaf's L2P by the kernel at precision T on the inputs rounded to T: the output chunk in f64.
fn kernel<T: RealScalar>(
    p: usize,
    kappa: f64,
    targets: &[[f64; 3]],
    locals: &[f64],
    vectors: usize,
    planes: usize,
    gradients: bool,
) -> Vec<f64> {
    let level = Level::new(Complex::new(T::zero(), T::from_f64(kappa)), T::one());
    assert!(level.is_yukawa());
    let x: Vec<[T; 3]> = targets.iter().map(|x| x.map(T::from_f64)).collect();
    let l: Vec<T> = locals.iter().map(|&v| T::from_f64(v)).collect();
    let per = planes * vectors * if gradients { 4 } else { 1 };
    let mut out = vec![T::zero(); per * targets.len()];
    LeafL2p::new(p).l2p(
        p,
        &level,
        &x,
        &l,
        &mut L2pScratch::new(p, vectors),
        LeafOutput {
            values: &mut out,
            vectors,
            planes,
            gradients,
        },
    );
    out.iter().map(|&v| v.to_f64()).collect()
}

/// A target's φ̃ and gradient.
type Evaluation = (f64, [f64; 3]);

/// The reference's outputs per real local (column) and target: `l2p_yukawa` with the gradient.
fn reference(p: usize, kappa: f64, targets: &[[f64; 3]], locals: &[f64]) -> Vec<Vec<Evaluation>> {
    let level = Level::new(Complex::new(0.0, kappa), 1.0);
    let frame = Frame::new([0.0; 3], 1.0);
    let mut ws = Workspace::new(p);
    let len = Layout::new(p).len();
    locals
        .chunks_exact(len)
        .map(|local| {
            targets
                .iter()
                .map(|&x| {
                    let (mut phi, mut g) = ([0.0], [[0.0; 3]]);
                    leaf::l2p_yukawa(
                        p,
                        &level,
                        &frame,
                        local,
                        &[x],
                        &mut ws,
                        &mut phi,
                        Some(&mut g),
                    );
                    (phi[0], g[0])
                })
                .collect()
        })
        .collect()
}

/// The bounds of one target and column without the factor C u, in f64 at the reference's r.
#[derive(Clone, Copy)]
struct Bound {
    /// φ̃'s and the gradient's contract sums.
    phi: f64,
    grad: f64,
    /// The f32 level's rounding: φ̃'s and the gradient's (β + n + 1) sums.
    level_phi: f64,
    level_grad: f64,
    /// The routines' flushed harmonics: Σₙ Sₙ (n + 1)(îₙ + κ |îₙ′| + tₙ) e.
    floor: f64,
    /// The flushed îₙ, îₙ′ and e: Σₙ wₙ Sₙ (1 + κ + 1/r).
    flush: f64,
}

/// The contract's sums per real local (column) and target.
fn bounds(p: usize, kappa: f64, targets: &[[f64; 3]], locals: &[f64]) -> Vec<Vec<Bound>> {
    let level = Level::new(Complex::new(0.0, kappa), 1.0);
    let (s, beta) = (level.scale(), level.damping());
    let len = Layout::new(p).len();
    let mut radial = vec![0.0; p + 1];
    let mut derivative = vec![0.0; p + 1];
    locals
        .chunks_exact(len)
        .map(|local| {
            let size: Vec<f64> = (0..=p)
                .map(|n| (n * n..(n + 1) * (n + 1)).map(|i| local[i].abs()).sum())
                .collect();
            targets
                .iter()
                .map(|x| {
                    let r = x[0].hypot(x[1]).hypot(x[2]);
                    let xr = kappa * r;
                    modified_i(xr, s, &mut radial, Some(&mut derivative));
                    let e = (xr - beta).exp();
                    let mut b = Bound {
                        phi: 0.0,
                        grad: 0.0,
                        level_phi: 0.0,
                        level_grad: 0.0,
                        floor: 0.0,
                        flush: 0.0,
                    };
                    for n in 0..=p {
                        let w = ((2 * n + 1) as f64 / (4.0 * PI)).sqrt() * size[n];
                        let (i, di) = (radial[n].abs(), derivative[n].abs());
                        let t = if r > 0.0 {
                            i / r
                        } else if n == 1 {
                            kappa / (3.0 * s)
                        } else {
                            0.0
                        };
                        let n1 = (n + 1) as f64;
                        b.phi += w * ((n1 * n1 + (2 * n + 1) as f64) * i + xr * (di + i));
                        b.grad += w
                            * (kappa * ((n1 * n1 + (2 * n + 3) as f64) * di + xr * i + xr * di)
                                + (n1 * n1 * n1 + (n * (n + 1)) as f64 + xr) * t);
                        let shift = beta + n1;
                        b.level_phi += w * shift * i;
                        b.level_grad += w * shift * (kappa * di + t);
                        b.floor += size[n] * n1 * (i + kappa * di + t);
                        b.flush += w * (1.0 + kappa + if r > 0.0 { 1.0 / r } else { 0.0 });
                    }
                    b.phi *= e;
                    b.grad *= e;
                    b.level_phi *= e;
                    b.level_grad *= e;
                    b.floor *= e;
                    b
                })
                .collect()
        })
        .collect()
}

/// The precision's terms beside the contract: u, the level's rounding (0 in f64, where kernel and
/// reference share the level's s and β), the floor's factor and τ.
#[derive(Clone, Copy)]
struct Precision {
    u: f64,
    level: f64,
    floor: f64,
    tau: f64,
}

const F64: Precision = Precision {
    u: f64::EPSILON / 2.0,
    level: 0.0,
    floor: 0.0,
    tau: 0.0,
};

fn f32_precision() -> Precision {
    Precision {
        u: f32::EPSILON as f64 / 2.0,
        level: 1.0,
        // The routines' flushed harmonics move a term by at most about 2⁻⁵⁷ of its size (CH4.9b2).
        floor: 2f64.powi(-50),
        tau: tau::<f32>(),
    }
}

/// The worst ratio of the error to the contract (times C) over the outputs of one leaf; panics on a
/// breach.
#[allow(clippy::too_many_arguments)]
fn worst<T: RealScalar>(
    p: usize,
    kappa: f64,
    targets: &[[f64; 3]],
    locals: &[f64],
    vectors: usize,
    planes: usize,
    gradients: bool,
    precision: Precision,
) -> f64 {
    let got = kernel::<T>(p, kappa, targets, locals, vectors, planes, gradients);
    let want = reference(p, kappa, targets, locals);
    let bound = bounds(p, kappa, targets, locals);
    let n = targets.len();
    let columns = planes * vectors;
    assert_eq!(want.len(), columns);
    let Precision {
        u,
        level,
        floor,
        tau,
    } = precision;
    let mut ratio: f64 = 0.0;
    let mut check = |a: f64, b: f64, sum: f64, level_sum: f64, b_: &Bound, what: &str| {
        assert!(a.is_finite(), "{what}: not finite");
        let beside = LEVEL_C * level * u * level_sum + floor * b_.floor + tau * b_.flush;
        let limit = C * u * sum + beside + f64::MIN_POSITIVE;
        let d = (a - b).abs();
        assert!(
            d <= limit,
            "p = {p}, κ = {kappa}, N = {n}, {what}: {a:e} against {b:e}, bound {limit:e}"
        );
        // The ratio in units of the contract, beside the terms held apart.
        let over = (d - beside).max(0.0);
        ratio = ratio.max(over / (u * sum + f64::MIN_POSITIVE));
    };
    for o in 0..columns {
        for i in 0..n {
            let (phi, g) = want[o][i];
            let b = bound[o][i];
            check(
                got[o * n + i],
                phi,
                b.phi,
                b.level_phi,
                &b,
                &format!("column {o}, target {i} {:?}, φ", targets[i]),
            );
            if gradients {
                for (axis, &ga) in g.iter().enumerate() {
                    check(
                        got[columns * n + 3 * (o * n + i) + axis],
                        ga,
                        b.grad,
                        b.level_grad,
                        &b,
                        &format!("column {o}, target {i} {:?}, ∇φ[{axis}]", targets[i]),
                    );
                }
            }
        }
    }
    ratio
}

#[test]
fn yukawa_l2p_within_the_contract_in_f64() {
    let mut most: f64 = 0.0;
    for (case, p) in [0, 1, 2, 5, 8, 17, 26, 45].into_iter().enumerate() {
        let sizes: &[usize] = if p >= 26 {
            &[1, 17, 40]
        } else {
            &[1, 5, 16, 17, 40]
        };
        let kappas = KAPPAS.iter().copied().chain([KAPPA_LARGE]);
        for (c, kappa) in kappas.enumerate() {
            for (s, &n) in sizes.iter().enumerate() {
                let seed = (case * 31 + c * 7 + s) as u64 * 2;
                let x = targets(n, seed);
                let vectors = 1 + (case + s) % 3;
                let planes = 1 + (case + c + s) % 2;
                let l = locals(p, vectors, planes, seed);
                for gradients in [false, true] {
                    most = most.max(worst::<f64>(
                        p, kappa, &x, &l, vectors, planes, gradients, F64,
                    ));
                }
            }
        }
    }
    println!("Yukawa L2P against the reference, f64: worst {most:.3} of the contract (C = {C})");
    assert!(most <= C);
}

#[test]
fn yukawa_l2p_within_the_contract_in_f32() {
    let round = |x: f64| x as f32 as f64;
    let mut most: f64 = 0.0;
    for (case, p) in [0, 1, 2, 5, 8, 17, 32].into_iter().enumerate() {
        for (c, &kappa) in KAPPAS.iter().enumerate() {
            for (s, n) in [1, 16, 40].into_iter().enumerate() {
                let seed = (case * 17 + c * 5 + s) as u64 * 2;
                let kappa = round(kappa);
                let x: Vec<[f64; 3]> = targets(n, seed).iter().map(|x| x.map(round)).collect();
                let planes = 1 + (case + s) % 2;
                let l: Vec<f64> = locals(p, 1, planes, seed).into_iter().map(round).collect();
                for gradients in [false, true] {
                    most = most.max(worst::<f32>(
                        p,
                        kappa,
                        &x,
                        &l,
                        1,
                        planes,
                        gradients,
                        f32_precision(),
                    ));
                }
            }
        }
    }
    println!("Yukawa L2P against the reference, f32: worst {most:.3} of the contract (C = {C})");
}

/// Targets whose r² underflows to zero although they are not the centre (the kernel takes the
/// centre's path, the tangential limit κ/(3s) and î at 0, where the reference's r > 0; the two agree
/// to rounding), and the centre alone and among others.
#[test]
fn yukawa_targets_whose_squares_underflow() {
    let x = [
        [1e-170, 0.0, 0.0],
        [0.0, -3e-170, 2e-170],
        [4e-200, 1e-190, -2e-180],
        [0.0; 3],
        [0.3, -0.2, 0.5],
    ];
    let mut most: f64 = 0.0;
    for p in [0, 1, 2, 8, 17] {
        for &kappa in &KAPPAS {
            let l = locals(p, 2, 1, p as u64);
            for gradients in [false, true] {
                most = most.max(worst::<f64>(p, kappa, &x, &l, 2, 1, gradients, F64));
            }
            let one = locals(p, 1, 2, 7);
            most = most.max(worst::<f64>(p, kappa, &x[3..4], &one, 1, 2, true, F64));
        }
    }
    println!("Yukawa L2P where r² underflows: worst {most:.3} of the contract (C = {C})");
}

/// The bits: each vector the one-vector call's; a complex vector's two columns the one-plane
/// kernel on each of its locals; each target's outputs the same alone, in another position and in a
/// leaf of 1, 16 and 40; the potentials with gradients those without.
#[test]
fn a_yukawa_targets_bits_do_not_depend_on_the_leaf_or_the_vectors() {
    let p = 12;
    let bits = |values: &[f64]| values.iter().map(|v| v.to_bits()).collect::<Vec<_>>();
    // Column o of target i with c columns in all: potential and the three gradient components.
    let outputs = |out: &[f64], columns: usize, n: usize, o: usize, i: usize| {
        let mut v = vec![out[o * n + i]];
        if out.len() == 4 * columns * n {
            v.extend((0..3).map(|a| out[columns * n + 3 * (o * n + i) + a]));
        }
        v
    };
    let len = Layout::new(p).len();
    for &kappa in &KAPPAS {
        let x = targets(40, 4);
        let n = x.len();
        for planes in [1, 2] {
            let l = locals(p, 3, planes, 4);
            let all = kernel::<f64>(p, kappa, &x, &l, 3, planes, true);
            // Each column the one-plane, one-vector kernel on its local.
            for o in 0..3 * planes {
                let one = kernel::<f64>(p, kappa, &x, &l[o * len..(o + 1) * len], 1, 1, true);
                for i in 0..n {
                    assert_eq!(
                        bits(&outputs(&all, 3 * planes, n, o, i)),
                        bits(&outputs(&one, 1, n, 0, i)),
                        "κ = {kappa}, {planes} planes, column {o}, target {i}"
                    );
                }
            }
        }
        let l1 = locals(p, 1, 1, 9);
        let full = kernel::<f64>(p, kappa, &x, &l1, 1, 1, true);
        let potentials = kernel::<f64>(p, kappa, &x, &l1, 1, 1, false);
        assert_eq!(bits(&full[..n]), bits(&potentials), "κ = {kappa}");
        // The centre, two random targets and the three next to the z-axis.
        for i in [0, 5, 17, 37, 38, 39] {
            for size in [1, LANES, 40] {
                let mut leaf: Vec<[f64; 3]> = x
                    .iter()
                    .copied()
                    .cycle()
                    .skip(i + 1)
                    .take(size - 1)
                    .collect();
                leaf.push(x[i]);
                let out = kernel::<f64>(p, kappa, &leaf, &l1, 1, 1, true);
                assert_eq!(
                    bits(&outputs(&out, 1, size, 0, size - 1)),
                    bits(&outputs(&full, 1, n, 0, i)),
                    "κ = {kappa}, target {i} in a leaf of {size}"
                );
            }
        }
    }
}

#[test]
fn yukawa_l2p_accumulates() {
    let p = 8;
    let kappa = 4.0;
    let x = targets(21, 6);
    let l = locals(p, 2, 2, 6);
    let level = Level::new(Complex::new(0.0, kappa), 1.0);
    let fresh = kernel::<f64>(p, kappa, &x, &l, 2, 2, true);
    let start: Vec<f64> = (0..fresh.len()).map(|i| i as f64 * 0.25 - 3.0).collect();
    let mut out = start.clone();
    LeafL2p::<f64>::new(p).l2p(
        p,
        &level,
        &x,
        &l,
        &mut L2pScratch::new(p, 2),
        LeafOutput {
            values: &mut out,
            vectors: 2,
            planes: 2,
            gradients: true,
        },
    );
    for ((o, s), f) in out.iter().zip(&start).zip(&fresh) {
        assert_eq!(o.to_bits(), (s + f).to_bits());
    }
}

#[test]
#[should_panic(expected = "1 output columns per vector")]
fn one_column_off_yukawa_is_refused() {
    let p = 4;
    let level = Level::new(Complex::new(2.0, 0.5), 1.0);
    let l = vec![0.0; Layout::new(p).len()];
    let mut out = vec![0.0; 1];
    LeafL2p::<f64>::new(p).l2p(
        p,
        &level,
        &[[0.1, 0.2, 0.3]],
        &l,
        &mut L2pScratch::new(p, 1),
        LeafOutput {
            values: &mut out,
            vectors: 1,
            planes: 1,
            gradients: false,
        },
    );
}

#[test]
#[should_panic(expected = "the locals must hold")]
fn short_yukawa_locals_are_refused() {
    let p = 4;
    let level = Level::new(Complex::new(0.0, 2.0), 1.0);
    // Two vectors of one plane need 2 (p + 1)² values.
    let l = vec![0.0; Layout::new(p).len()];
    let mut out = vec![0.0; 2];
    LeafL2p::<f64>::new(p).l2p(
        p,
        &level,
        &[[0.1, 0.2, 0.3]],
        &l,
        &mut L2pScratch::new(p, 2),
        LeafOutput {
            values: &mut out,
            vectors: 2,
            planes: 1,
            gradients: false,
        },
    );
}

#[test]
#[should_panic(expected = "values per target")]
fn a_short_yukawa_output_is_refused() {
    let p = 4;
    let level = Level::new(Complex::new(0.0, 2.0), 1.0);
    let l = vec![0.0; Layout::new(p).len()];
    let mut out = vec![0.0; 2];
    LeafL2p::<f64>::new(p).l2p(
        p,
        &level,
        &[[0.1, 0.2, 0.3]],
        &l,
        &mut L2pScratch::new(p, 1),
        LeafOutput {
            values: &mut out,
            vectors: 1,
            planes: 1,
            gradients: true,
        },
    );
}

#[test]
fn the_planes() {
    let yukawa = Level::new(Complex::new(0.0, 2.0), 1.0);
    let real = Level::new(Complex::new(2.0, 0.0), 1.0);
    let complex = Level::new(Complex::new(2.0, 0.5), 1.0);
    assert_eq!(LeafL2p::<f64>::planes(&yukawa, false), 1);
    assert_eq!(LeafL2p::<f64>::planes(&yukawa, true), 2);
    for level in [real, complex] {
        for c in [false, true] {
            assert_eq!(LeafL2p::<f64>::planes(&level, c), 2);
        }
    }
}
