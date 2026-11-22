//! The transfer function's length (Helmholtz CH4.3d1): one unit of the scan of lengths L at
//! one (kind, k a, p), every class and placement.
//!
//! S6's chain (`crate::run`) truncates T at 2P − 2 = 2p, where in exact arithmetic it is the
//! truncated M2L that RCR applies. At a shorter length L it drops the (S|R) terms of degree
//! n + n′ > L: a truncation of its own, traded against the rounding of T's largest terms,
//! which grow exponentially once l passes k|t| (proposal §3.6). Per class and placement this
//! module measures, at each length of [`lengths`]:
//!
//! - `trunc`: S7's RCR against the direct sum (the degree's own truncation, as in `run`);
//! - `dd`: the chain in double-double at length L against double-double RCR (T's truncation
//!   alone; computed where p ≤ [`Config::dd_max_p`], NaN elsewhere);
//! - `dev`: the f64 chain (T rounded once from double-double, as the FMM builds it) at length
//!   L against double-double RCR (truncation and rounding);
//! - `round`: the f64 chain against the double-double chain at the same length (the rounding
//!   alone; NaN where the double-double chain is not run);
//! - `total`: the f64 chain against the direct sum;
//!
//! with log₁₀ of E's term at L ([`log10_term`]) and of the truncation model's sum
//! ([`log10_tail`]). The f64 chain at L = 2p is `run`'s `f64:2P-2` variant bit for bit: the
//! same samples (T's partial sums do not depend on which snapshots are taken), the same
//! signature and the same operations in the same order (`tests/lengths.rs`).

use std::time::Instant;

use tetraleaf_spike_truncation::calibrate::parallel_map;
use tetraleaf_spike_truncation::chain::{self, Field, Kind, Setup, Work, e2p, p2e, rcr};
use tetraleaf_spike_truncation::complex::Complex;
use tetraleaf_spike_truncation::geometry::{Placement, dist2};
use tetraleaf_spike_truncation::harmonics::{Harmonics, len};
use tetraleaf_spike_truncation::real::{Dd, Real};
use tetraleaf_spike_truncation::rotation::{Direction, Rotation};

use crate::classes::{Class, classes};
use crate::grid::Grid;
use crate::run::{self, cloud, l2p_many, narrow, widen};
use crate::signature::{f2r_ring, finish, s2f_ring};
use crate::transfer::{Transfer, log10_tail, log10_term};

/// The excess-bandwidth length ⌈2ka + 1.8 (log₁₀ 1/ε)^{2/3} (2ka)^{1/3}⌉ of the
/// high-frequency FMM literature (review 2026-11-12, M8), with Song and Chew's constant 1.8
/// and κa = Re(k) a.
pub fn excess_bandwidth(ka: f64, eps: f64) -> usize {
    let x = 2.0 * ka;
    (x + 1.8 * (1.0 / eps).log10().powf(2.0 / 3.0) * x.cbrt() - 1e-9).ceil() as usize
}

/// The lengths of a unit at degree p: from L_lo = max(p, min(⌈3p/2⌉, ⌊2ka⌋ − ⌈(2ka)^{1/3}⌉))
/// to 2p in steps of max(1, ⌈(2p − L_lo)/24⌉), plus ⌈3p/2⌉ and the excess-bandwidth lengths
/// of `epsilons` ([`excess_bandwidth`]) where they lie in [p, 2p]; ascending and distinct.
pub fn lengths(ka: f64, p: usize, epsilons: &[f64]) -> Vec<usize> {
    let x = 2.0 * ka;
    let below = (x.floor() - x.cbrt().ceil()).max(0.0) as usize;
    let lo = p.max((3 * p).div_ceil(2).min(below)).min(2 * p);
    let step = (2 * p - lo).div_ceil(24).max(1);
    let mut out: Vec<usize> = (lo..=2 * p).step_by(step).collect();
    out.push(2 * p);
    out.push((3 * p).div_ceil(2));
    out.extend(epsilons.iter().map(|&e| excess_bandwidth(ka, e)));
    out.retain(|&l| (p..=2 * p).contains(&l));
    out.sort_unstable();
    out.dedup();
    out
}

/// Points, seeds and limits.
#[derive(Clone, Copy, Debug)]
pub struct Config {
    /// S6's configuration (points, seeds, threads, the term limit).
    pub base: run::Config,
    /// The largest p at which the double-double chain runs at every length.
    pub dd_max_p: usize,
}

impl Default for Config {
    fn default() -> Self {
        Config {
            base: run::Config::default(),
            dd_max_p: 100,
        }
    }
}

/// The errors (potential, gradient) of one length on one placement.
#[derive(Clone, Copy, Debug)]
pub struct LengthErrors {
    /// The double-double chain against double-double RCR (NaN where not run).
    pub dd: (f64, f64),
    /// The f64 chain against double-double RCR.
    pub dev: (f64, f64),
    /// The f64 chain against the double-double chain (NaN where not run).
    pub round: (f64, f64),
    /// The f64 chain against the direct sum.
    pub total: (f64, f64),
}

/// The results of one class.
#[derive(Clone, Debug)]
pub struct ClassLengths {
    /// Index into [`classes`].
    pub class: usize,
    /// The class.
    pub info: Class,
    /// Per length: log₁₀ |k||t| maxⱼ≤L (2j + 1)|hⱼ(k|t|)|.
    pub log10_term: Vec<f64>,
    /// Per length: [`log10_tail`] (−∞ at L = 2p).
    pub log10_tail: Vec<f64>,
    /// Per length: log₁₀ max |T_L| over the samples (NaN where not run).
    pub log10_max_t: Vec<f64>,
    /// Per placement: RCR against the direct sum.
    pub trunc: Vec<(f64, f64)>,
    /// Per placement and length.
    pub errors: Vec<Vec<LengthErrors>>,
    /// Wall time of the class, seconds.
    pub seconds: f64,
}

/// The results of a unit.
#[derive(Clone, Debug)]
pub struct UnitLengths {
    /// The kind of k.
    pub kind: Kind,
    /// κ a at the M2L level.
    pub ka: f64,
    /// The degree.
    pub p: usize,
    /// The lengths scanned, ascending.
    pub lengths: Vec<usize>,
    /// Per class run.
    pub classes: Vec<ClassLengths>,
    /// Wall time, seconds.
    pub seconds: f64,
}

/// Points and their charges.
type Cloud = (Vec<[Dd; 3]>, Vec<Complex<Dd>>);

/// Per placement: the sources, the double-double multipole and its f64 signature (ring
/// after ring), as `run` forms them.
struct Shared {
    sources: Vec<Cloud>,
    multipoles: Vec<Vec<Complex<Dd>>>,
    sig64: Vec<Vec<Complex<f64>>>,
}

/// Runs the classes `class_ids` (indices into [`classes`]) at (kind, k a, p) and the lengths
/// `ls` (ascending, within [p, 2p]).
pub fn run_lengths(
    kind: Kind,
    ka: f64,
    p: usize,
    ls: &[usize],
    class_ids: &[usize],
    cfg: &Config,
) -> UnitLengths {
    assert!(p >= 2, "p ≥ 2");
    assert!(
        !ls.is_empty()
            && ls.windows(2).all(|w| w[0] < w[1])
            && ls[0] >= p
            && *ls.last().unwrap() <= 2 * p,
        "lengths ascending within [p, 2p]"
    );
    let start = Instant::now();
    let base = &cfg.base;
    let mut setup = Setup::<Dd>::new(kind, ka);
    setup.s = [Dd::ONE; 3];
    let grid = Grid::new(p);
    let nphi = grid.nphi;
    let sources: Vec<_> = Placement::ALL
        .iter()
        .map(|&pl| cloud(base, pl, 1, [0; 3], base.n_src))
        .collect();
    let multipoles: Vec<Vec<Complex<Dd>>> = sources
        .iter()
        .map(|(pts, q)| p2e(true, setup.k, Dd::ONE, [Dd::ZERO; 3], pts, q, p))
        .collect();
    let tab64 = grid.tables::<f64>();
    let m64: Vec<Vec<Complex<f64>>> = multipoles
        .iter()
        .map(|m| m.iter().map(|&z| narrow(z)).collect())
        .collect();
    let rings: Vec<usize> = (0..grid.rings()).collect();
    let per_ring: Vec<Vec<Vec<Complex<f64>>>> = parallel_map(&rings, base.threads, |&i| {
        let mut harm = Harmonics::<Dd>::new(p);
        grid.eval_ring(i, &mut harm);
        let mut th64 = Vec::new();
        grid.theta_of::<f64>(&harm, &mut th64);
        let mut ab64 = vec![Complex::ZERO; 2 * (p + 1)];
        (0..Placement::ALL.len())
            .map(|pl| {
                let mut f = vec![Complex::ZERO; nphi];
                s2f_ring(p, &th64, &m64[pl], &tab64, &mut ab64, &mut f);
                f
            })
            .collect()
    });
    let mut sig64 = vec![Vec::with_capacity(grid.samples()); Placement::ALL.len()];
    for r in per_ring {
        for (pl, f) in r.into_iter().enumerate() {
            sig64[pl].extend_from_slice(&f);
        }
    }
    let shared = Shared {
        sources,
        multipoles,
        sig64,
    };
    let all = classes();
    let results = parallel_map(class_ids, base.threads, |&ci| {
        run_class(&setup, &grid, &shared, ls, ci, all[ci], cfg)
    });
    UnitLengths {
        kind,
        ka,
        p,
        lengths: ls.to_vec(),
        classes: results,
        seconds: start.elapsed().as_secs_f64(),
    }
}

fn run_class(
    setup: &Setup<Dd>,
    grid: &Grid,
    shared: &Shared,
    ls: &[usize],
    ci: usize,
    info: Class,
    cfg: &Config,
) -> ClassLengths {
    let start = Instant::now();
    let base = &cfg.base;
    let p = grid.p;
    let nphi = grid.nphi;
    let k = setup.k;
    let rep = info.rep;
    let c_t: [Dd; 3] = rep.map(|x| Dd::from_f64(x as f64));
    let np = shared.sources.len();
    // The references, as `run`: S7's RCR and the direct sum, in double-double.
    let rot = Rotation::new(Direction::new(c_t), p);
    let table = setup.m2l_table(dist2(rep), p);
    let mut work = Work::new(p);
    let mut targets = Vec::with_capacity(np);
    let mut exact = Vec::with_capacity(np);
    let mut reference: Vec<Vec<Field<Dd>>> = Vec::with_capacity(np);
    let mut trunc = Vec::with_capacity(np);
    for (pl, &placement) in Placement::ALL.iter().enumerate() {
        let mut local = vec![Complex::ZERO; len(p)];
        let (tgt, _) = cloud(base, placement, 4, rep, base.n_tgt);
        let (src, q) = &shared.sources[pl];
        let ex = chain::direct(k, src, q, &tgt);
        rcr(
            &rot,
            &table,
            p,
            p,
            None,
            &shared.multipoles[pl],
            &mut local,
            &mut work,
        );
        let mut fields = Vec::new();
        e2p(
            true,
            k,
            Dd::ONE,
            c_t,
            &local,
            p,
            &tgt,
            &mut work,
            &mut fields,
        );
        trunc.push(chain::errors(&fields, &ex));
        targets.push(tgt);
        exact.push(ex);
        reference.push(fields);
    }
    let kc = Complex::new(k.re.to_f64(), k.im.to_f64());
    let t_abs = (dist2(rep) as f64).sqrt();
    let kt = kc * t_abs;
    let two_ka = kc * (2.0 * chain::unit_radius::<f64>());
    let log10_term: Vec<f64> = ls.iter().map(|&l| log10_term(kt, l)).collect();
    let log10_tail: Vec<f64> = ls
        .iter()
        .map(|&l| log10_tail(two_ka, kt, l, 2 * p))
        .collect();
    // The lengths whose T stays in the double-double range, as `run`'s bands.
    let snaps: Vec<usize> = ls
        .iter()
        .zip(&log10_term)
        .filter(|(_, t)| **t <= base.max_log10_term)
        .map(|(l, _)| *l)
        .collect();
    let nb = snaps.len();
    let nan = (f64::NAN, f64::NAN);
    let inf = (f64::INFINITY, f64::INFINITY);
    let mut errors = vec![
        vec![
            LengthErrors {
                dd: nan,
                dev: inf,
                round: nan,
                total: inf,
            };
            ls.len()
        ];
        np
    ];
    let mut log10_max_t = vec![f64::NAN; ls.len()];
    if nb > 0 {
        let with_dd = p <= cfg.dd_max_p;
        let tr = Transfer::new(k, rep, &snaps);
        let tab64 = grid.tables::<f64>();
        let tabdd = grid.tables::<Dd>();
        let mdd = &shared.multipoles;
        let mut a64 = vec![vec![vec![Complex::<f64>::ZERO; len(p)]; nb]; np];
        let mut add = if with_dd {
            vec![vec![vec![Complex::<Dd>::ZERO; len(p)]; nb]; np]
        } else {
            Vec::new()
        };
        let mut harm = Harmonics::<Dd>::new(p);
        let (mut th64, mut thdd) = (Vec::new(), Vec::new());
        let mut tv = vec![Complex::ZERO; nb];
        let mut tdd = vec![vec![Complex::<Dd>::ZERO; nphi]; nb];
        let mut t64 = vec![vec![Complex::<f64>::ZERO; nphi]; nb];
        let mut g64 = vec![Complex::<f64>::ZERO; nphi];
        let mut gdd = vec![Complex::<Dd>::ZERO; nphi];
        let mut fdd = vec![Complex::<Dd>::ZERO; nphi];
        let mut cs64 = vec![Complex::<f64>::ZERO; 2 * (p + 1)];
        let mut csdd = vec![Complex::<Dd>::ZERO; 2 * (p + 1)];
        let mut max_t = vec![0.0f64; nb];
        for i in 0..grid.rings() {
            grid.eval_ring(i, &mut harm);
            grid.theta_of::<f64>(&harm, &mut th64);
            if with_dd {
                grid.theta_of::<Dd>(&harm, &mut thdd);
            }
            for j in 0..nphi {
                let mu = tr.mu(grid.direction(i, j));
                tr.eval(mu, &mut tv);
                for b in 0..nb {
                    t64[b][j] = narrow(tv[b]);
                    tdd[b][j] = tv[b];
                    max_t[b] = max_t[b].max(t64[b][j].abs());
                }
            }
            let range = i * nphi..(i + 1) * nphi;
            for pl in 0..np {
                let f64s = &shared.sig64[pl][range.clone()];
                for b in 0..nb {
                    for j in 0..nphi {
                        g64[j] = t64[b][j] * f64s[j];
                    }
                    f2r_ring(
                        p,
                        &th64,
                        tab64.w[i],
                        &g64,
                        &tab64,
                        &mut cs64,
                        &mut a64[pl][b],
                    );
                }
                if with_dd {
                    s2f_ring(p, &thdd, &mdd[pl], &tabdd, &mut csdd, &mut fdd);
                    for b in 0..nb {
                        for j in 0..nphi {
                            gdd[j] = tdd[b][j] * fdd[j];
                        }
                        f2r_ring(
                            p,
                            &thdd,
                            tabdd.w[i],
                            &gdd,
                            &tabdd,
                            &mut csdd,
                            &mut add[pl][b],
                        );
                    }
                }
            }
        }
        // Snapshot b → index into `ls`.
        let at: Vec<usize> = snaps
            .iter()
            .map(|l| ls.iter().position(|x| x == l).unwrap())
            .collect();
        for (b, &li) in at.iter().enumerate() {
            log10_max_t[li] = max_t[b].log10();
        }
        for pl in 0..np {
            let mut locals: Vec<Vec<Complex<Dd>>> = Vec::with_capacity(2 * nb);
            for acc in a64[pl].iter_mut() {
                finish(p, acc);
                locals.push(acc.iter().map(|&z| widen(z)).collect());
            }
            if with_dd {
                for acc in add[pl].iter_mut() {
                    finish(p, acc);
                    locals.push(acc.clone());
                }
            }
            let views: Vec<&[Complex<Dd>]> = locals.iter().map(|l| &l[..]).collect();
            let fields = l2p_many(k, c_t, &views, p, &targets[pl]);
            for (b, &li) in at.iter().enumerate() {
                let e = &mut errors[pl][li];
                e.dev = chain::errors(&fields[b], &reference[pl]);
                e.total = chain::errors(&fields[b], &exact[pl]);
                if with_dd {
                    e.dd = chain::errors(&fields[nb + b], &reference[pl]);
                    e.round = chain::errors(&fields[b], &fields[nb + b]);
                }
            }
        }
    }
    ClassLengths {
        class: ci,
        info,
        log10_term,
        log10_tail,
        log10_max_t,
        trunc,
        errors,
        seconds: start.elapsed().as_secs_f64(),
    }
}
