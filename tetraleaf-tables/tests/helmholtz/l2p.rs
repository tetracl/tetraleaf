//! The batched L2P of a leaf (Helmholtz CH4.9b2, `tetraleaf_tables::helmholtz::leaf::LeafL2p`)
//! against `tetraleaf_ref::helmholtz::leaf::l2p`, its trusted path.
//!
//! Error measure: per target and output (φ̃, each gradient component), the complex modulus of
//! the difference over the contract of the module documentation ("L2P"), formed in f64 from the
//! reference's radial functions at the reference's r (`hypot`), C = [`C`]; in f32 against the f64
//! reference on the f32 inputs, u = 2⁻²⁴, plus a floor for the routines' flushed values. The worst
//! ratio is printed. Bit identities (vectors, positions, potentials with and without gradients,
//! accumulation) by `to_bits`.

use core::f64::consts::PI;

use tetraleaf_math::radial::bessel_j;
use tetraleaf_math::{Complex, Layout, RealScalar};
use tetraleaf_ref::Frame;
use tetraleaf_ref::helmholtz::{Level, Workspace, leaf};
use tetraleaf_tables::helmholtz::leaf::{L2pScratch, LANES, LeafL2p, LeafOutput};

use super::leaf::{C, KS, Lcg, points};

/// m complex locals of degree p, the vectors one after the other, uniform in [−1, 1).
fn locals(p: usize, vectors: usize, seed: u64) -> Vec<f64> {
    let mut rng = Lcg(seed ^ 0x51ed_270b);
    (0..2 * vectors * Layout::new(p).len())
        .map(|_| rng.next())
        .collect()
}

/// The targets of a leaf of n: `points`' (the centre, the poles, near-pole points, the axes and
/// a corner for an even seed), with two points next to the z-axis (squares below the normal
/// range) in place of the last two random ones when n ≥ 12.
fn targets(n: usize, seed: u64) -> Vec<[f64; 3]> {
    let mut x = points(n, seed);
    if n >= 12 {
        x[n - 2] = [1e-160, 0.0, 1.0];
        x[n - 1] = [-1e-40, 1e-40, -0.5];
    }
    x
}

/// One leaf's L2P by the kernel at precision T on the inputs rounded to T: the output chunk in f64.
fn kernel<T: RealScalar>(
    p: usize,
    k: Complex<f64>,
    targets: &[[f64; 3]],
    locals: &[f64],
    vectors: usize,
    gradients: bool,
) -> Vec<f64> {
    let level = Level::new(Complex::new(T::from_f64(k.re), T::from_f64(k.im)), T::one());
    let x: Vec<[T; 3]> = targets.iter().map(|x| x.map(T::from_f64)).collect();
    let l: Vec<T> = locals.iter().map(|&v| T::from_f64(v)).collect();
    let per = 2 * vectors * if gradients { 4 } else { 1 };
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
            gradients,
        },
    );
    out.iter().map(|&v| v.to_f64()).collect()
}

/// The reference's outputs per vector and target: φ̃ and the gradient.
fn reference(
    p: usize,
    k: Complex<f64>,
    targets: &[[f64; 3]],
    locals: &[f64],
    vectors: usize,
) -> Vec<Vec<(Complex<f64>, [Complex<f64>; 3])>> {
    let level = Level::new(k, 1.0);
    let frame = Frame::new([0.0; 3], 1.0);
    let mut ws = Workspace::new(p);
    let len = 2 * Layout::new(p).len();
    let zero = Complex::new(0.0, 0.0);
    (0..vectors)
        .map(|v| {
            let local = &locals[v * len..(v + 1) * len];
            targets
                .iter()
                .map(|&x| {
                    let (mut phi, mut g) = ([zero], [[zero; 3]]);
                    leaf::l2p(
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

/// The contract's sums per vector and target without the factor C u: (φ̃'s, the gradient's), and
/// the floor's (Σₙ Sₙ (n + 1)(|ĵₙ| + |k| |ĵₙ′| + tₙ) e), all in f64 at the reference's r.
fn bounds(
    p: usize,
    k: Complex<f64>,
    targets: &[[f64; 3]],
    locals: &[f64],
    vectors: usize,
) -> Vec<Vec<(f64, f64, f64)>> {
    let level = Level::new(k, 1.0);
    let s = level.scale();
    let layout = Layout::new(p);
    let len = layout.len();
    let zero = Complex::new(0.0, 0.0);
    let mut radial = vec![zero; p + 1];
    let mut derivative = vec![zero; p + 1];
    (0..vectors)
        .map(|v| {
            let local = &locals[2 * v * len..2 * (v + 1) * len];
            let size: Vec<f64> = (0..=p)
                .map(|n| {
                    (n * n..(n + 1) * (n + 1))
                        .map(|i| Complex::new(local[i], local[len + i]).norm())
                        .sum()
                })
                .collect();
            targets
                .iter()
                .map(|x| {
                    let r = x[0].hypot(x[1]).hypot(x[2]);
                    let z = k * r;
                    bessel_j(z, s, &mut radial, Some(&mut derivative));
                    let e = (k.im * r - level.damping()).exp();
                    let (mut phi, mut grad, mut floor) = (0.0, 0.0, 0.0);
                    for n in 0..=p {
                        let w = ((2 * n + 1) as f64 / (4.0 * PI)).sqrt() * size[n];
                        let (j, dj) = (radial[n].norm(), derivative[n].norm());
                        let t = if r > 0.0 {
                            j / r
                        } else if n == 1 {
                            k.norm() / (3.0 * s)
                        } else {
                            0.0
                        };
                        let n1 = (n + 1) as f64;
                        phi += w * ((n1 * n1 + (2 * n + 1) as f64) * j + z.norm() * dj);
                        grad += w
                            * (k.norm() * ((n1 * n1 + (2 * n + 3) as f64) * dj + z.norm() * j)
                                + (n1 * n1 * n1 + (n * (n + 1)) as f64) * t);
                        floor += size[n] * n1 * (j + k.norm() * dj + t);
                    }
                    (phi * e, grad * e, floor * e)
                })
                .collect()
        })
        .collect()
}

/// The worst ratio of the error to the contract (times C) over the outputs of one leaf; panics on a
/// breach. `floor` multiplies the floor's sum (the flushed values in f32).
#[allow(clippy::too_many_arguments)]
fn worst<T: RealScalar>(
    p: usize,
    k: Complex<f64>,
    targets: &[[f64; 3]],
    locals: &[f64],
    vectors: usize,
    gradients: bool,
    u: f64,
    floor: f64,
) -> f64 {
    let got = kernel::<T>(p, k, targets, locals, vectors, gradients);
    let want = reference(p, k, targets, locals, vectors);
    let bound = bounds(p, k, targets, locals, vectors);
    let n = targets.len();
    let mut ratio: f64 = 0.0;
    let mut check = |a: Complex<f64>, b: Complex<f64>, sum: f64, fl: f64, what: &str| {
        assert!(a.re.is_finite() && a.im.is_finite(), "{what}: not finite");
        let limit = C * u * sum + floor * fl + f64::MIN_POSITIVE;
        let d = (a - b).norm();
        assert!(
            d <= limit,
            "p = {p}, k = {k}, N = {n}, {what}: {a} against {b}, bound {limit:e}"
        );
        ratio = ratio.max(d / limit * C);
    };
    for v in 0..vectors {
        for i in 0..n {
            let (phi, g) = want[v][i];
            let (bp, bg, fl) = bound[v][i];
            let at = |column: usize| got[column * n + i];
            let a = Complex::new(at(2 * v), at(2 * v + 1));
            check(
                a,
                phi,
                bp,
                fl,
                &format!("vector {v}, target {i} {:?}, φ", targets[i]),
            );
            if gradients {
                for (axis, &ga) in g.iter().enumerate() {
                    let at = |column: usize| got[2 * vectors * n + 3 * (column * n + i) + axis];
                    let a = Complex::new(at(2 * v), at(2 * v + 1));
                    check(
                        a,
                        ga,
                        bg,
                        fl,
                        &format!("vector {v}, target {i} {:?}, ∇φ[{axis}]", targets[i]),
                    );
                }
            }
        }
    }
    ratio
}

#[test]
fn l2p_within_the_contract_in_f64() {
    let u = f64::EPSILON / 2.0;
    let mut most: f64 = 0.0;
    for (case, p) in [0, 1, 2, 5, 8, 17, 26, 45].into_iter().enumerate() {
        let sizes: &[usize] = if p >= 26 {
            &[1, 17, 40]
        } else {
            &[1, 5, 16, 17, 40]
        };
        for &(re, im) in &KS {
            let k = Complex::new(re, im);
            for (s, &n) in sizes.iter().enumerate() {
                let seed = (case * 31 + s) as u64 * 2;
                let x = targets(n, seed);
                let vectors = 1 + (case + s) % 3;
                let l = locals(p, vectors, seed);
                for gradients in [false, true] {
                    most = most.max(worst::<f64>(p, k, &x, &l, vectors, gradients, u, 0.0));
                }
            }
        }
    }
    println!("L2P against the reference, f64: worst {most:.3} of the contract's u units (C = {C})");
}

#[test]
fn l2p_within_the_contract_in_f32() {
    let u = f32::EPSILON as f64 / 2.0;
    let mut most: f64 = 0.0;
    for (case, p) in [0, 1, 2, 5, 8, 17, 32].into_iter().enumerate() {
        for &(re, im) in &KS {
            for (s, n) in [1, 16, 40].into_iter().enumerate() {
                let seed = (case * 17 + s) as u64 * 2;
                let round = |x: f64| x as f32 as f64;
                let k = Complex::new(round(re), round(im));
                let x: Vec<[f64; 3]> = targets(n, seed).iter().map(|x| x.map(round)).collect();
                let l: Vec<f64> = locals(p, 1, seed).into_iter().map(round).collect();
                // The routines' flushed values (CH4.9a2a, CH4.9b1) move a term by at most about
                // 2⁻⁵⁷ of its size; 2⁻⁵⁰ per degree's sum is far below u.
                let floor = 2f64.powi(-50);
                for gradients in [false, true] {
                    most = most.max(worst::<f32>(p, k, &x, &l, 1, gradients, u, floor));
                }
            }
        }
    }
    println!("L2P against the reference, f32: worst {most:.3} of the contract's u units (C = {C})");
}

/// The bits: each vector the one-vector call's; each target's outputs the same alone, in another
/// position and in a leaf of 1, 16 and 40; the potentials with gradients those without.
#[test]
fn a_targets_bits_do_not_depend_on_the_leaf_or_the_vectors() {
    let p = 12;
    for &(re, im) in &KS {
        let k = Complex::new(re, im);
        let x = targets(40, 4);
        let l = locals(p, 3, 4);
        let len = 2 * Layout::new(p).len();
        let all = kernel::<f64>(p, k, &x, &l, 3, true);
        let n = x.len();
        let bits = |values: &[f64]| values.iter().map(|v| v.to_bits()).collect::<Vec<_>>();
        // Output column o of target i, potential and the three gradient components.
        let outputs = |out: &[f64], vectors: usize, n: usize, o: usize, i: usize| {
            let mut v = vec![out[o * n + i]];
            if out.len() == 8 * vectors * n {
                v.extend((0..3).map(|a| out[2 * vectors * n + 3 * (o * n + i) + a]));
            }
            v
        };
        for v in 0..3 {
            let one = kernel::<f64>(p, k, &x, &l[v * len..(v + 1) * len], 1, true);
            for i in 0..n {
                for part in 0..2 {
                    assert_eq!(
                        bits(&outputs(&all, 3, n, 2 * v + part, i)),
                        bits(&outputs(&one, 1, n, part, i)),
                        "k = {k}, vector {v}, target {i}"
                    );
                }
            }
        }
        let l1 = &l[..len];
        let full = kernel::<f64>(p, k, &x, l1, 1, true);
        let potentials = kernel::<f64>(p, k, &x, l1, 1, false);
        assert_eq!(bits(&full[..2 * n]), bits(&potentials), "k = {k}");
        for i in [0, 5, 17, 39] {
            for size in [1, LANES, 40] {
                // Target i at position `size − 1`, the others shifted.
                let mut leaf: Vec<[f64; 3]> = x
                    .iter()
                    .copied()
                    .cycle()
                    .skip(i + 1)
                    .take(size - 1)
                    .collect();
                leaf.push(x[i]);
                let out = kernel::<f64>(p, k, &leaf, l1, 1, true);
                for part in 0..2 {
                    assert_eq!(
                        bits(&outputs(&out, 1, size, part, size - 1)),
                        bits(&outputs(&full, 1, n, part, i)),
                        "k = {k}, target {i} in a leaf of {size}"
                    );
                }
            }
        }
    }
}

#[test]
fn l2p_accumulates() {
    let p = 8;
    let k = Complex::new(4.0, 0.5);
    let x: Vec<[f64; 3]> = targets(21, 6);
    let l = locals(p, 2, 6);
    let level = Level::new(k, 1.0);
    let kernel = LeafL2p::<f64>::new(p);
    let mut scratch = L2pScratch::new(p, 2);
    let fresh = self::kernel::<f64>(p, k, &x, &l, 2, true);
    let start: Vec<f64> = (0..fresh.len()).map(|i| i as f64 * 0.25 - 3.0).collect();
    let mut out = start.clone();
    kernel.l2p(
        p,
        &level,
        &x,
        &l,
        &mut scratch,
        LeafOutput {
            values: &mut out,
            vectors: 2,
            gradients: true,
        },
    );
    for ((o, s), f) in out.iter().zip(&start).zip(&fresh) {
        assert_eq!(o.to_bits(), (s + f).to_bits());
    }
    // A larger scratch and kernel give the same bits.
    let mut again = vec![0.0; fresh.len()];
    LeafL2p::<f64>::new(20).l2p(
        p,
        &level,
        &x,
        &l,
        &mut L2pScratch::new(30, 3),
        LeafOutput {
            values: &mut again,
            vectors: 2,
            gradients: true,
        },
    );
    assert!(
        again
            .iter()
            .zip(&fresh)
            .all(|(a, b)| a.to_bits() == b.to_bits())
    );
}

#[test]
#[should_panic(expected = "not Yukawa's")]
fn yukawa_is_refused() {
    let p = 4;
    let level = Level::new(Complex::new(0.0, 2.0), 1.0);
    let l = vec![0.0; 2 * Layout::new(p).len()];
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
            gradients: false,
        },
    );
}

#[test]
#[should_panic(expected = "the locals must hold")]
fn short_locals_are_refused() {
    let p = 4;
    let level = Level::new(Complex::new(2.0, 0.0), 1.0);
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
            gradients: false,
        },
    );
}

#[test]
#[should_panic(expected = "values per target")]
fn a_short_output_is_refused() {
    let p = 4;
    let level = Level::new(Complex::new(2.0, 0.0), 1.0);
    let l = vec![0.0; 2 * Layout::new(p).len()];
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
            gradients: true,
        },
    );
}

#[test]
#[should_panic(expected = "the scratch holds 1 to 1")]
fn more_vectors_than_the_scratch_are_refused() {
    let p = 4;
    let level = Level::new(Complex::new(2.0, 0.0), 1.0);
    let l = vec![0.0; 4 * Layout::new(p).len()];
    let mut out = vec![0.0; 4];
    LeafL2p::<f64>::new(p).l2p(
        p,
        &level,
        &[[0.1, 0.2, 0.3]],
        &l,
        &mut L2pScratch::new(p, 1),
        LeafOutput {
            values: &mut out,
            vectors: 2,
            gradients: false,
        },
    );
}

#[test]
#[should_panic(expected = "above the scratch's")]
fn a_small_scratch_is_refused() {
    let p = 4;
    let level = Level::new(Complex::new(2.0, 0.0), 1.0);
    let l = vec![0.0; 2 * Layout::new(p).len()];
    let mut out = vec![0.0; 2];
    LeafL2p::<f64>::new(p).l2p(
        p,
        &level,
        &[[0.1, 0.2, 0.3]],
        &l,
        &mut L2pScratch::new(p - 1, 1),
        LeafOutput {
            values: &mut out,
            vectors: 1,
            gradients: false,
        },
    );
}
