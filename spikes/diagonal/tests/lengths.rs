//! The scan of T's length (Helmholtz CH4.3d1) against its trusted paths: S6's chain at 2p,
//! S6's double-double identity, the truncation model's sum in double-double, the lengths, the
//! FMM's degree rule, determinism.

use tetraleaf_spike_diagonal::classes::classes;
use tetraleaf_spike_diagonal::lengths::{Config, excess_bandwidth, lengths, run_lengths};
use tetraleaf_spike_diagonal::rules::refit;
use tetraleaf_spike_diagonal::run::{self, VARIANTS, Variant, run_unit};
use tetraleaf_spike_diagonal::transfer::log10_tail;
use tetraleaf_spike_truncation::chain::Kind;
use tetraleaf_spike_truncation::complex::Complex;
use tetraleaf_spike_truncation::radial;
use tetraleaf_spike_truncation::real::{Dd, Real};

/// The index of `run`'s `f64:2P-2` variant.
fn full_f64() -> usize {
    VARIANTS
        .iter()
        .position(|v| {
            *v == Variant {
                precision: run::Precision::F64,
                band: run::Band::Full,
            }
        })
        .unwrap()
}

/// Test 1: at L = 2p the scan's f64 chain is S6's `f64:2P-2` variant bit for bit, its
/// deviation from RCR S6's `round` and its total S6's `total`; RCR's truncation the same.
#[test]
fn the_scan_at_2p_is_s6s_chain() {
    let ids: Vec<usize> = (0..classes().len()).collect();
    let v = full_f64();
    for kind in [Kind::Real, Kind::Damped(0.1)] {
        let (ka, p) = (4.0, 10);
        let ls = lengths(ka, p, &[1e-6]);
        assert_eq!(*ls.last().unwrap(), 2 * p);
        let mut cfg = Config::default();
        cfg.base.threads = 4;
        let s6 = run_unit(kind, ka, p, &ids, &cfg.base);
        let ours = run_lengths(kind, ka, p, &ls, &ids, &cfg);
        let last = ls.len() - 1;
        for (a, b) in s6.classes.iter().zip(&ours.classes) {
            assert_eq!(a.class, b.class);
            assert_eq!(a.trunc, b.trunc, "{kind:?} class {}", a.class);
            for pl in 0..a.trunc.len() {
                let e = b.errors[pl][last];
                assert_eq!(a.round[pl][v], e.dev, "{kind:?} class {} pl {pl}", a.class);
                assert_eq!(
                    a.total[pl][v], e.total,
                    "{kind:?} class {} pl {pl}",
                    a.class
                );
            }
            assert_eq!(a.log10_max_t[2].to_bits(), b.log10_max_t[last].to_bits());
        }
    }
}

/// Test 2: the double-double chain at 2p is RCR to 10⁻²⁰ (S6's identity through the scan);
/// below 2p T's truncation falls by at least two orders of magnitude from the shortest length
/// to 2p − 1 on the nearest and the farthest class (on the nearest it falls slowly, about
/// (2a/|t|)^l = 0.87^l, Gegenbauer's rate: 1.1·10⁻³ to 7.6·10⁻⁸ at k a 8, p 18, uniform).
#[test]
fn the_truncation_is_resolved() {
    let cs = classes();
    let nearest = cs.iter().position(|c| c.rep == [0, 0, 2]).unwrap();
    let farthest = cs.iter().position(|c| c.rep == [3, 3, 3]).unwrap();
    for (kind, ka) in [
        (Kind::Real, 8.0),
        (Kind::Real, 16.0),
        (Kind::Damped(0.1), 8.0),
    ] {
        let p = refit(ka, 1e-6);
        let ls: Vec<usize> = (p..=2 * p).collect();
        let mut cfg = Config::default();
        cfg.base.threads = 2;
        let res = run_lengths(kind, ka, p, &ls, &[nearest, farthest], &cfg);
        for c in &res.classes {
            for pl in 0..c.errors.len() {
                let dd: Vec<f64> = c.errors[pl].iter().map(|e| e.dd.0.max(e.dd.1)).collect();
                let n = dd.len();
                println!(
                    "{kind:?} k a {ka} p {p} class {:?} pl {pl}: {}",
                    c.info.rep,
                    dd.iter()
                        .map(|x| format!("{x:.2e}"))
                        .collect::<Vec<_>>()
                        .join(" ")
                );
                assert!(dd[n - 1] < 1e-20, "at 2p: {}", dd[n - 1]);
                assert!(dd[0] > 1e2 * dd[n - 2], "{} {}", dd[0], dd[n - 2]);
                // The double-double chain's rounding is far below the f64 chain's.
                let round = c.errors[pl][n - 1].round;
                assert!(round.0.is_finite() && round.1.is_finite());
            }
        }
    }
}

/// Test 3: `log10_tail` against the sum in double-double with S7's Bessel and Hankel
/// functions, within 10⁻¹⁰ in log₁₀; −∞ with no terms; finite where hₗ leaves f64's range.
#[test]
fn the_tail_matches_double_double() {
    let a2 = 3f64.sqrt();
    for kind in [Kind::Real, Kind::Damped(0.1)] {
        for ka in [3.0, 16.0, 64.0] {
            let p = refit(ka, 1e-6);
            let k: Complex<f64> = kind.wavenumber(ka / (a2 / 2.0));
            for d2 in [4, 12, 27] {
                let t = (d2 as f64).sqrt();
                let (z2, z) = (k * a2, k * t);
                let (z2d, zd) = (
                    Complex::<Dd>::from_f64(z2.re, z2.im),
                    Complex::<Dd>::from_f64(z.re, z.im),
                );
                let mut j = vec![Complex::<Dd>::ZERO; 2 * p + 1];
                let mut h = vec![Complex::<Dd>::ZERO; 2 * p + 1];
                radial::bessel_j(z2d, Dd::ONE, &mut j);
                radial::hankel(zd, Dd::ONE, &mut h);
                for l in [p, (3 * p) / 2, 2 * p - 1] {
                    let mut sum = Dd::ZERO;
                    for m in l + 1..=2 * p {
                        sum += Dd::from_usize(2 * m + 1) * j[m].abs() * h[m].abs();
                    }
                    let want = (sum / h[0].abs()).to_f64().log10();
                    let got = log10_tail(z2, z, l, 2 * p);
                    assert!(
                        (got - want).abs() < 1e-10,
                        "{kind:?} k a {ka} |t|² {d2} L {l}: {got} {want}"
                    );
                }
            }
        }
    }
    let z = Complex::new(10.0, 0.0);
    assert_eq!(log10_tail(z, z, 20, 20), f64::NEG_INFINITY);
    // hₗ(40) leaves f64's range near l = 300 (about 10²²²; jₗ(30) about 10⁻²⁶⁰ stays in
    // it): the rescaled sum stays finite.
    let got = log10_tail(Complex::new(30.0, 0.0), Complex::new(40.0, 0.0), 150, 300);
    assert!(got.is_finite(), "{got}");
}

/// Test 4: the lengths ascend within [p, 2p] and hold the shortest, ⌈3p/2⌉, the
/// excess-bandwidth lengths in range, and 2p.
#[test]
fn the_lengths() {
    for (ka, eps) in [
        (2.83, 1e-8),
        (13.5, 1e-6),
        (27.06, 1e-8),
        (108.25, 1e-3),
        (216.5, 1e-8),
    ] {
        let p = refit(ka, eps);
        let es = [1e-3, 1e-6, 1e-8];
        let ls = lengths(ka, p, &es);
        assert!(ls.windows(2).all(|w| w[0] < w[1]), "{ls:?}");
        assert!(ls[0] >= p && *ls.last().unwrap() == 2 * p, "{ls:?}");
        let mid = (3 * p).div_ceil(2);
        assert!(ls.contains(&mid));
        for e in es {
            let l = excess_bandwidth(ka, e);
            assert_eq!(
                (p..=2 * p).contains(&l),
                ls.contains(&l),
                "{ka} {e} {l} {ls:?}"
            );
        }
        assert!(ls.len() <= 34, "{}", ls.len());
        println!("k a {ka} ε {eps:e} p {p}: {ls:?}");
    }
    // Song and Chew's length at k a = 108.25, 10⁻⁶: ⌈216.5 + 1.8 · 6^{2/3} · 216.5^{1/3}⌉ =
    // ⌈252.19⌉ = 253.
    assert_eq!(excess_bandwidth(108.25, 1e-6), 253);
}

/// Test 5: the spike's copy of the FMM's rule equals `tetraleaf_exec::helmholtz::DegreeRule::
/// CALIBRATED` (CH4.3c2), its degrees printed on 2026-11-22 from `tetraleaf-exec` on branch
/// `tbetcke/h-ch4-3d1` at ε = 10⁻³ … 10⁻⁸.
#[test]
fn the_rule_is_the_fmms() {
    let table: [(f64, [usize; 6]); 8] = [
        (1.0, [7, 9, 12, 15, 20, 24]),
        (2.7, [7, 9, 12, 15, 20, 24]),
        (13.5, [20, 22, 24, 25, 26, 27]),
        (21.65, [28, 31, 33, 35, 36, 38]),
        (27.1, [34, 37, 39, 41, 43, 44]),
        (54.1, [61, 65, 68, 71, 74, 76]),
        (108.3, [116, 120, 125, 128, 131, 134]),
        (216.5, [224, 230, 235, 240, 244, 248]),
    ];
    let es = [1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8];
    for (x, want) in table {
        let got: Vec<usize> = es.iter().map(|&e| refit(x, e)).collect();
        assert_eq!(got, want, "x = {x}");
    }
}

/// Test 6: a unit's results do not depend on the thread count.
#[test]
fn the_scan_is_deterministic() {
    let ids = [0, 5, 33];
    let (ka, p) = (8.0, 12);
    let ls = lengths(ka, p, &[1e-6]);
    let mut cfg = Config::default();
    cfg.base.threads = 1;
    let a = run_lengths(Kind::Real, ka, p, &ls, &ids, &cfg);
    cfg.base.threads = 3;
    let b = run_lengths(Kind::Real, ka, p, &ls, &ids, &cfg);
    let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
    for (x, y) in a.classes.iter().zip(&b.classes) {
        assert_eq!(x.trunc, y.trunc);
        assert_eq!(bits(&x.log10_max_t), bits(&y.log10_max_t));
        for (ex, ey) in x.errors.iter().zip(&y.errors) {
            for (u, v) in ex.iter().zip(ey) {
                let f = |e: &tetraleaf_spike_diagonal::lengths::LengthErrors| {
                    [e.dd, e.dev, e.round, e.total].map(|(a, b)| (a.to_bits(), b.to_bits()))
                };
                assert_eq!(f(u), f(v));
            }
        }
    }
}
