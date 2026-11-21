//! The calibration of the Helmholtz FMM's per-level degree p_l over (k a_l, ε) in the full
//! FMM (Helmholtz CH3.3; docs/design/helmholtz-fmm-proposal.md §3.5 and §10, CH3.3;
//! docs/helmholtz/tasks/CH3.3.md). Free of MPI and of the FMM: the `helmholtz_calibrate`
//! example measures the error curves with [`helmholtz_accuracy`](crate::helmholtz_accuracy),
//! and this module finds the degrees in them and fits the rule.
//!
//! - **Curves.** A [`Curve`] holds the relative ℓ2 errors [φ, ∇φ] of the FMM at consecutive
//!   degrees p = `first`, `first` + 1, …, with every other setting fixed: one common
//!   degree on every level ≥ 2 (a plateau scan), or one level's degree with the others held
//!   high (a level scan).
//! - **The smallest degree.** [`Curve::smallest_degree`] is the smallest p at which the
//!   larger of the two errors meets the target at p, p + 1 and p + 2 ([`RUN`]), as S7's
//!   search (`spikes/truncation/SPIKE_REPORT.md` §1.4): the errors are not monotone in p to
//!   the last digit, so one lucky degree below the threshold is not taken.
//! - **Points.** A [`Point`] is one measured degree: the level's x = Re(k) a_l (a_l the
//!   circumradius), the tolerance ε and the degree that meets the target (ε/2, the plan's
//!   margin).
//! - **The rule.** [`DegreeRule`] has S7's shape with fitted constants and a quadratic
//!   plateau (the full FMM's plateau grows faster than linearly in L = log₁₀(1/ε)),
//!
//!   ```text
//!   p(x, ε) = max(⌈γ L² + α L + β⌉, ⌈x + C (ln 1/ε)^{2/3} x^{1/3}⌉),
//!   ```
//!
//!   the low-frequency plateau and the high-frequency branch. [`fit_plateau`] fits γ, α
//!   and β to the plateau points, [`fit_high`] fits C to the points the plateau does not cover,
//!   and [`covers`] checks that the rule is at least every point's degree.
//! - **The branch refitted** (Helmholtz CH4.3c2). [`fit_branch`] fits the branch
//!   ⌈x + (C (ln 1/ε)^{2/3} + a L + b) x^{1/3} + d⌉ in one of the [`BranchForm`]s beside a
//!   given plateau, optionally at most a cap's rule on [`cap_grid`]; [`CH3_3_LEVELS`] are
//!   CH3.3's level points as data.

use tetraleaf_exec::helmholtz::DegreeRule;

/// The consecutive degrees that must each meet a target (module documentation).
pub const RUN: usize = 3;

/// The errors of the FMM at consecutive degrees (module documentation).
#[derive(Clone, Debug, PartialEq)]
pub struct Curve {
    /// The degree of `errors[0]`.
    pub first: usize,
    /// The relative ℓ2 errors [φ, ∇φ] at `first`, `first` + 1, ….
    pub errors: Vec<[f64; 2]>,
}

impl Curve {
    /// An empty curve starting at degree `first`.
    pub fn new(first: usize) -> Self {
        Self {
            first,
            errors: Vec::new(),
        }
    }

    /// The degree the next measurement is for.
    pub fn next_degree(&self) -> usize {
        self.first + self.errors.len()
    }

    /// Appends the errors at [`Curve::next_degree`].
    pub fn push(&mut self, errors: [f64; 2]) {
        self.errors.push(errors);
    }

    /// The smallest degree p of the curve whose larger error meets `target` (≤) at p and at
    /// the next [`RUN`] − 1 degrees, all of them measured; `None` if there is none. A NaN
    /// error never meets a target.
    pub fn smallest_degree(&self, target: f64) -> Option<usize> {
        let meets: Vec<bool> = self
            .errors
            .iter()
            .map(|e| e[0].max(e[1]) <= target)
            .collect();
        meets
            .windows(RUN)
            .position(|w| w.iter().all(|&m| m))
            .map(|i| self.first + i)
    }

    /// Whether the last [`RUN`] degrees all meet `target`: the scan for every target at
    /// least as large is complete.
    pub fn reached(&self, target: f64) -> bool {
        self.errors.len() >= RUN
            && self.errors[self.errors.len() - RUN..]
                .iter()
                .all(|e| e[0].max(e[1]) <= target)
    }
}

/// One measured degree (module documentation).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Point {
    /// x = Re(k) a_l of the level, a_l its circumradius; for a plateau scan the largest x
    /// of the scanned levels.
    pub x: f64,
    /// The tolerance ε.
    pub epsilon: f64,
    /// The smallest degree that meets the target at ε.
    pub degree: usize,
}

/// The step of the grids of α, β and C: the constants are multiples of it, rounded up.
pub const STEP: f64 = 0.05;

/// Whether `rule` is at least every point's degree.
pub fn covers(rule: &DegreeRule, points: &[Point]) -> bool {
    points
        .iter()
        .all(|p| rule.degree(p.x, p.epsilon) >= p.degree)
}

/// The smallest multiple of [`STEP`] above `value` + 10⁻⁹: a constant at which a
/// [`DegreeRule`] ceiling of a quantity that must exceed p − 1 reaches p.
fn step_above(value: f64) -> f64 {
    (((value + 1e-9) / STEP).floor() + 1.0) * STEP
}

/// The step of the grid of γ.
pub const GAMMA_STEP: f64 = 0.01;

/// The plateau of the points: over γ on the grid of [`GAMMA_STEP`] in [0, 0.5] and α on the
/// grid of [`STEP`] in [0, 10] (so the plateau never falls as ε does), β the smallest
/// multiple of [`STEP`] with ⌈γ L² + α L + β⌉ ≥ p at every point (L = log₁₀(1/ε)), the
/// constants with the least total excess Σ (⌈γ L² + α L + β⌉ − p) over the points, ties to
/// the least Σ (γ L² + α L + β). The x of the points is not used; C is 0. `None` without
/// points.
pub fn fit_plateau(points: &[Point]) -> Option<DegreeRule> {
    if points.is_empty() {
        return None;
    }
    let decades: Vec<f64> = points.iter().map(|p| (1.0 / p.epsilon).log10()).collect();
    let mut best: Option<(usize, f64, DegreeRule)> = None;
    for g in 0..=(0.5 / GAMMA_STEP).round() as usize {
        let gamma = g as f64 * GAMMA_STEP;
        for i in 0..=(10.0 / STEP).round() as usize {
            let alpha = i as f64 * STEP;
            let curve = |l: f64| (gamma * l + alpha) * l;
            let beta = step_above(
                points
                    .iter()
                    .zip(&decades)
                    .map(|(p, &l)| p.degree as f64 - 1.0 - curve(l))
                    .fold(f64::NEG_INFINITY, f64::max),
            );
            let rule = DegreeRule::gd(gamma, alpha, beta, 0.0);
            let excess = points
                .iter()
                .map(|p| rule.plateau(p.epsilon).saturating_sub(p.degree))
                .sum::<usize>();
            let total = decades.iter().map(|&l| curve(l) + beta).sum::<f64>();
            if best.is_none_or(|(e, t, _)| excess < e || (excess == e && total < t)) {
                best = Some((excess, total, rule));
            }
        }
    }
    best.map(|(_, _, rule)| rule)
}

/// C of the high-frequency branch beside `plateau`'s plateau: the smallest multiple of
/// [`STEP`] with ⌈x + C (ln 1/ε)^{2/3} x^{1/3}⌉ ≥ p at every point whose degree the plateau
/// does not reach; 0 if the plateau reaches every one. Points with x = 0 that the plateau
/// misses cannot be covered: `None`.
pub fn fit_high(points: &[Point], plateau: &DegreeRule) -> Option<f64> {
    let mut c: Option<f64> = None;
    for p in points
        .iter()
        .filter(|p| plateau.plateau(p.epsilon) < p.degree)
    {
        if p.x <= 0.0 {
            return None;
        }
        let needed =
            (p.degree as f64 - 1.0 - p.x) / ((1.0 / p.epsilon).ln().powf(2.0 / 3.0) * p.x.cbrt());
        c = Some(c.map_or(needed, |c: f64| c.max(needed)));
    }
    Some(c.map_or(0.0, |c| step_above(c).max(0.0)))
}

/// The rule fitted to the plateau points (all levels low-frequency) and the other points
/// (the level scans, and plateau scans of trees with a high-frequency level): γ, α and β
/// from the first, C from both. `None` if either fit fails.
pub fn fit(plateau: &[Point], others: &[Point]) -> Option<DegreeRule> {
    let plateau_rule = fit_plateau(plateau)?;
    let all: Vec<Point> = plateau.iter().chain(others).copied().collect();
    let c = fit_high(&all, &plateau_rule)?;
    let rule = DegreeRule { c, ..plateau_rule };
    covers(&rule, &all).then_some(rule)
}

/// The ε of CH3.3's level scans, [`CH3_3_LEVELS`]' columns: 10⁻³ … 10⁻⁸.
pub const CH3_3_EPSILONS: [f64; 6] = [1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8];

/// CH3.3's level scans (`tetraleaf-validate/results/helmholtz-ch3-3.md` §3; Helmholtz CH4.3c2 keeps
/// them as data, since the scans' outputs were not archived): (workload, kD, level, x, the
/// smallest degree meeting ε/2 at each of [`CH3_3_EPSILONS`]). x = Re(k) a_l to four places
/// (`tetraleaf-exec`'s `tests/helmholtz/degree.rs`, `MEASURED`).
pub const CH3_3_LEVELS: [(&str, f64, usize, f64, [usize; 6]); 16] = [
    ("w3:108849", 30.0, 4, 1.6238, [6, 8, 11, 14, 18, 22]),
    ("w1", 100.0, 5, 2.7063, [7, 8, 11, 14, 17, 20]),
    ("w1", 30.0, 3, 3.2474, [7, 9, 11, 14, 17, 20]),
    ("w3:108849", 30.0, 3, 3.2474, [7, 9, 11, 14, 17, 22]),
    ("w1", 100.0, 4, 5.4127, [9, 11, 13, 15, 18, 21]),
    ("w4-cube:100000", 100.0, 4, 5.4127, [10, 11, 13, 15, 18, 21]),
    ("w1", 30.0, 2, 6.4948, [10, 12, 14, 16, 19, 22]),
    ("w3:108849", 30.0, 2, 6.4948, [10, 12, 14, 16, 19, 22]),
    (
        "w4-clusters:100000",
        100.0,
        3,
        6.7862,
        [8, 10, 11, 13, 14, 16],
    ),
    ("w1", 100.0, 3, 10.8253, [15, 17, 19, 20, 22, 24]),
    (
        "w4-cube:100000",
        100.0,
        3,
        10.8253,
        [15, 17, 19, 21, 23, 24],
    ),
    (
        "w4-clusters:100000",
        100.0,
        2,
        13.5725,
        [13, 16, 18, 20, 22, 24],
    ),
    ("w1", 100.0, 2, 21.6507, [25, 28, 31, 33, 35, 37]),
    (
        "w4-cube:100000",
        100.0,
        2,
        21.6507,
        [26, 28, 31, 33, 35, 37],
    ),
    ("w1", 200.0, 2, 43.3018, [47, 51, 54, 57, 60, 62]),
    ("w1", 300.0, 2, 64.9528, [69, 73, 77, 80, 83, 86]),
];

/// [`CH3_3_LEVELS`] as points, one per row and ε.
pub fn ch3_3_level_points() -> Vec<Point> {
    CH3_3_LEVELS
        .iter()
        .flat_map(|&(_, _, _, x, degrees)| {
            CH3_3_EPSILONS
                .iter()
                .zip(degrees)
                .map(move |(&epsilon, degree)| Point { x, epsilon, degree })
        })
        .collect()
}

/// The step of b in [`fit_branch`]: b is the smallest multiple of it that covers the points.
pub const BRANCH_STEP: f64 = 0.01;

/// The forms of the high-frequency branch that [`fit_branch`] fits (Helmholtz CH4.3c2), each
/// ⌈x + (C (ln 1/ε)^{2/3} + a L + b) x^{1/3} + d⌉ with some constants fixed at 0, and the
/// grids of the others (b is never on a grid: [`BRANCH_STEP`]).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BranchForm {
    /// ⌈x + (a L + b) x^{1/3}⌉, the excess linear in L (CH3.3's report §6): a in [0, 2] by
    /// 0.01.
    Linear,
    /// ⌈x + (a L + b) x^{1/3} + d⌉: a in [0, 2] by 0.01, d in [0, 10] by 0.25.
    LinearOffset,
    /// ⌈x + (C (ln 1/ε)^{2/3} + a L + b) x^{1/3} + d⌉: C in [−1, 3] by 0.1, a in [−1, 2] by
    /// 0.05, d in [0, 4] by 0.5; the excess may then be concave in L, as the measured one is
    /// at small ε.
    Full,
}

impl BranchForm {
    /// The grids of (C, a, d).
    fn grids(self) -> [Vec<f64>; 3] {
        let grid = |lo: f64, hi: f64, step: f64| -> Vec<f64> {
            let n = ((hi - lo) / step).round() as usize;
            (0..=n).map(|i| lo + i as f64 * step).collect()
        };
        match self {
            Self::Linear => [vec![0.0], grid(0.0, 2.0, 0.01), vec![0.0]],
            Self::LinearOffset => [vec![0.0], grid(0.0, 2.0, 0.01), grid(0.0, 10.0, 0.25)],
            Self::Full => [
                grid(-1.0, 3.0, 0.1),
                grid(-1.0, 2.0, 0.05),
                grid(0.0, 4.0, 0.5),
            ],
        }
    }
}

/// The (x, ε) on which [`fit_branch`] holds a cap: x = 4 · 2^{k/16} for k = 0 … 112 (4 to
/// 512) and ε = 10⁻³ … 10⁻⁸ (D1's range), at the decades. Between them a rule fitted under
/// a cap may exceed it: CH4.3c2's fit by one degree at 2 of the 5,763 points of the same x and
/// ε = 10^{−j/10}, j = 30 … 80; holding the cap there too raised the fit's degrees at
/// x = 216 by 2–4. The default rule takes CH3.3's branch as its `DegreeRule::cap` instead, at
/// most CH3.3's everywhere (Helmholtz CH4.3c2's report).
pub fn cap_grid() -> Vec<(f64, f64)> {
    (0..=112)
        .flat_map(|k| {
            let x = 4.0 * 2f64.powf(f64::from(k) / 16.0);
            (3..=8).map(move |d| (x, 10f64.powi(-d)))
        })
        .collect()
}

/// The high-frequency branch beside `plateau`'s plateau (its γ, α and β kept) in `form`: over
/// the form's grids of C, a and d, b the smallest multiple of [`BRANCH_STEP`] with the rule
/// at least every point's degree, the constants with the least total excess Σ (rule − p)
/// over the points, ties to the smaller Σ of the branch's values. With a `cap`, only rules
/// at most the cap's degree on [`cap_grid`] count (Helmholtz CH4.3c2: no degree above CH3.3's
/// rule in D1). Points the plateau reaches constrain nothing but count in the excess. Points
/// with x = 0 above the plateau cannot be covered: `None`, as without points or without a
/// rule under the cap.
pub fn fit_branch(
    points: &[Point],
    plateau: &DegreeRule,
    form: BranchForm,
    cap: Option<&DegreeRule>,
) -> Option<DegreeRule> {
    if points.is_empty() {
        return None;
    }
    let base = DegreeRule {
        c: 0.0,
        a: 0.0,
        b: 0.0,
        d: 0.0,
        cap: None,
        ..*plateau
    };
    let above: Vec<&Point> = points
        .iter()
        .filter(|p| base.plateau(p.epsilon) < p.degree)
        .collect();
    if above.iter().any(|p| p.x <= 0.0) {
        return None;
    }
    let caps: Vec<(f64, f64, usize)> = cap.map_or_else(Vec::new, |cap| {
        cap_grid()
            .into_iter()
            .map(|(x, e)| (x, e, cap.degree(x, e)))
            .collect()
    });
    let [cs, alphas, offsets] = form.grids();
    let mut best: Option<(usize, f64, DegreeRule)> = None;
    for &c in &cs {
        for &a in &alphas {
            for &d in &offsets {
                // ⌈x + (C G + a L + b) x^{1/3} + d − 10⁻⁹⌉ ≥ p ⟺
                // b > (p − 1 − x − d + 10⁻⁹) / x^{1/3} − C G − a L.
                let partial = DegreeRule { c, a, d, ..base };
                let needed = above
                    .iter()
                    .map(|p| {
                        (p.degree as f64 - 1.0 - p.x - d + 1e-9) / p.x.cbrt()
                            - partial.excess(p.epsilon)
                    })
                    .fold(f64::NEG_INFINITY, f64::max);
                let b = if needed.is_finite() {
                    ((needed / BRANCH_STEP).floor() + 1.0) * BRANCH_STEP
                } else {
                    -10.0
                };
                let rule = DegreeRule { b, ..partial };
                if caps.iter().any(|&(x, e, p)| rule.degree(x, e) > p) {
                    continue;
                }
                let excess = points
                    .iter()
                    .map(|p| rule.degree(p.x, p.epsilon).saturating_sub(p.degree))
                    .sum::<usize>();
                let total = points
                    .iter()
                    .map(|p| p.x + rule.excess(p.epsilon) * p.x.cbrt() + d)
                    .sum::<f64>();
                if best.is_none_or(|(e, t, _)| excess < e || (excess == e && total < t)) {
                    best = Some((excess, total, rule));
                }
            }
        }
    }
    best.map(|(_, _, rule)| rule).filter(|r| covers(r, points))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A curve whose larger error is 10^(−p/3), with optional overrides.
    fn curve(first: usize, last: usize, overrides: &[(usize, f64)]) -> Curve {
        let mut c = Curve::new(first);
        for p in first..=last {
            let e = overrides
                .iter()
                .find(|(q, _)| *q == p)
                .map_or(10f64.powf(-(p as f64) / 3.0), |&(_, e)| e);
            c.push([e / 10.0, e]);
        }
        c
    }

    #[test]
    fn smallest_degree_of_a_monotone_curve() {
        let c = curve(1, 30, &[]);
        // 10^(−p/3) ≤ 10⁻³ from p = 9 on.
        assert_eq!(c.smallest_degree(1e-3), Some(9));
        assert_eq!(c.smallest_degree(1.0), Some(1));
        assert_eq!(c.next_degree(), 31);
        // The last three are 28, 29 and 30; 10^(−28/3) ≈ 4.6e-10.
        assert_eq!(c.smallest_degree(5e-10), Some(28));
        assert!(c.reached(5e-10));
        assert!(!c.reached(1e-10));
        assert_eq!(c.smallest_degree(1e-10), None);
    }

    #[test]
    fn a_lucky_degree_below_the_threshold_is_not_taken() {
        // p = 6 dips below 10⁻³ but p = 7 does not.
        let c = curve(1, 20, &[(6, 1e-4), (7, 2e-3)]);
        assert_eq!(c.smallest_degree(1e-3), Some(9));
        // A bump after the threshold moves it past the bump.
        let c = curve(1, 20, &[(10, 5e-3)]);
        assert_eq!(c.smallest_degree(1e-3), Some(11));
    }

    #[test]
    fn nan_never_meets_and_short_curves_have_none() {
        let c = curve(1, 20, &[(10, f64::NAN)]);
        assert_eq!(c.smallest_degree(1e-3), Some(11));
        assert_eq!(curve(5, 6, &[]).smallest_degree(1.0), None);
        assert!(!curve(5, 6, &[]).reached(1.0));
        assert_eq!(Curve::new(3).smallest_degree(1.0), None);
    }

    #[test]
    fn the_larger_error_decides() {
        let mut c = Curve::new(1);
        for _ in 0..3 {
            c.push([1e-9, 2e-3]);
        }
        assert_eq!(c.smallest_degree(1e-3), None);
        assert_eq!(c.smallest_degree(2e-3), Some(1));
    }

    /// Points made by a known rule at the decades 3..=10 and a few x.
    fn points_of(rule: &DegreeRule, xs: &[f64]) -> Vec<Point> {
        (3..=10)
            .flat_map(|d| {
                let epsilon = 10f64.powi(-d);
                xs.iter().map(move |&x| Point {
                    x,
                    epsilon,
                    degree: rule.degree(x, epsilon),
                })
            })
            .collect()
    }

    #[test]
    fn the_plateau_fit_recovers_a_rule_and_covers_its_points() {
        let truth = DegreeRule::gd(0.3, 0.2, 2.8, 0.0);
        let points = points_of(&truth, &[0.5]);
        let fitted = fit_plateau(&points).unwrap();
        assert!(covers(&fitted, &points));
        // The ceilings of the data are matched exactly: no excess.
        for p in &points {
            assert_eq!(fitted.plateau(p.epsilon), p.degree, "{p:?} {fitted:?}");
        }
        assert_eq!(fit_plateau(&[]), None);
    }

    #[test]
    fn the_plateau_fit_is_never_below_a_point() {
        let points =
            [(1e-3, 6), (1e-4, 9), (1e-5, 13), (1e-6, 15), (1e-8, 23)].map(|(epsilon, degree)| {
                Point {
                    x: 1.0,
                    epsilon,
                    degree,
                }
            });
        let rule = fit_plateau(&points).unwrap();
        assert!(covers(&rule, &points), "{rule:?}");
        assert!(rule.gamma >= 0.0 && rule.alpha >= 0.0, "{rule:?}");
    }

    #[test]
    fn the_full_fit_covers_plateau_and_level_points() {
        let truth = DegreeRule::gd(0.0, 3.0, -2.0, 0.8);
        let plateau = points_of(&truth, &[0.1, 1.0]);
        let levels = points_of(&truth, &[8.0, 21.7, 43.3, 64.9]);
        let rule = fit(&plateau, &levels).unwrap();
        assert!(
            covers(&rule, &plateau) && covers(&rule, &levels),
            "{rule:?}"
        );
        assert!(rule.c <= truth.c + 1e-12, "{rule:?}");
        // The plateau's constants come from the plateau alone.
        let alone = fit_plateau(&plateau).unwrap();
        assert_eq!(
            (alone.gamma, alone.alpha, alone.beta),
            (rule.gamma, rule.alpha, rule.beta)
        );
    }

    #[test]
    fn a_plateau_point_at_x_zero_above_the_plateau_cannot_be_fitted() {
        let points = [Point {
            x: 0.0,
            epsilon: 1e-3,
            degree: 10,
        }];
        let low = DegreeRule::gd(0.0, 0.0, 0.0, 0.0);
        assert_eq!(fit_high(&points, &low), None);
        let high = DegreeRule { beta: 10.0, ..low };
        assert_eq!(fit_high(&points, &high), Some(0.0));
    }

    /// Level points of a rule at the decades 3..=8 and the given x.
    fn level_points(rule: &DegreeRule, xs: &[f64]) -> Vec<Point> {
        (3..=8)
            .flat_map(|d| {
                let epsilon = 10f64.powi(-d);
                xs.iter().map(move |&x| Point {
                    x,
                    epsilon,
                    degree: rule.degree(x, epsilon),
                })
            })
            .collect()
    }

    #[test]
    fn the_branch_fit_covers_and_matches_a_rule_of_its_form() {
        let xs = [5.4, 10.8, 21.7, 43.3, 64.9, 108.3, 216.5];
        let plateau = DegreeRule::CH3_3;
        for (form, truth) in [
            (
                BranchForm::Linear,
                DegreeRule {
                    c: 0.0,
                    a: 0.8,
                    b: -0.3,
                    ..plateau
                },
            ),
            (
                BranchForm::LinearOffset,
                DegreeRule {
                    c: 0.0,
                    a: 0.8,
                    b: -2.2,
                    d: 3.5,
                    ..plateau
                },
            ),
        ] {
            let points = level_points(&truth, &xs);
            let fitted = fit_branch(&points, &plateau, form, None).unwrap();
            assert!(covers(&fitted, &points), "{form:?}: {fitted:?}");
            // The truth is on the grid, so the fit's excess is at most its (none).
            for p in &points {
                assert_eq!(fitted.degree(p.x, p.epsilon), p.degree, "{form:?}: {p:?}");
            }
            assert_eq!(
                (fitted.gamma, fitted.alpha, fitted.beta, fitted.c),
                (plateau.gamma, plateau.alpha, plateau.beta, 0.0)
            );
            if form == BranchForm::Linear {
                assert_eq!(fitted.d, 0.0);
            }
        }
    }

    #[test]
    fn the_branch_fit_on_ch3_3s_points() {
        let points = ch3_3_level_points();
        assert_eq!(points.len(), 96);
        // CH3.3's rule covers its own points.
        assert!(covers(&DegreeRule::CH3_3, &points));
        let excess = |rule: &DegreeRule| {
            points
                .iter()
                .map(|p| rule.degree(p.x, p.epsilon) - p.degree)
                .sum::<usize>()
        };
        let a = fit_branch(&points, &DegreeRule::CH3_3, BranchForm::Linear, None).unwrap();
        let b = fit_branch(&points, &DegreeRule::CH3_3, BranchForm::LinearOffset, None).unwrap();
        // The offset form contains the linear one (d = 0), so it fits at least as closely; the
        // linear one fits CH3.3's points more closely than CH3.3's own branch (its report §6).
        assert!(excess(&b) <= excess(&a), "{a:?} {b:?}");
        assert!(excess(&a) <= excess(&DegreeRule::CH3_3), "{a:?}");
        assert_eq!(
            fit_branch(&[], &DegreeRule::CH3_3, BranchForm::Linear, None),
            None
        );
    }

    /// The full form recovers a rule of its form on the grid (no excess at the points), and a
    /// fit under a cap stays at or below the cap on [`cap_grid`] while it covers the points.
    #[test]
    fn the_full_form_and_the_cap() {
        let xs = [5.4, 10.8, 21.7, 43.3, 64.9, 108.3, 216.5];
        let plateau = DegreeRule::CH3_3;
        let truth = DegreeRule {
            c: 1.8,
            a: -0.4,
            b: -5.18,
            d: 3.5,
            ..plateau
        };
        let points = level_points(&truth, &xs);
        let full = fit_branch(&points, &plateau, BranchForm::Full, None).unwrap();
        assert!(covers(&full, &points), "{full:?}");
        for p in &points {
            assert_eq!(full.degree(p.x, p.epsilon), p.degree, "{full:?}: {p:?}");
        }
        let points = ch3_3_level_points();
        let capped =
            fit_branch(&points, &plateau, BranchForm::LinearOffset, Some(&plateau)).unwrap();
        assert!(covers(&capped, &points), "{capped:?}");
        for (x, epsilon) in cap_grid() {
            assert!(
                capped.degree(x, epsilon) <= plateau.degree(x, epsilon),
                "{capped:?} at x = {x}, ε = {epsilon:e}"
            );
        }
    }

    #[test]
    fn steps_lie_above() {
        assert!((step_above(0.0) - 0.05).abs() < 1e-12);
        assert!((step_above(0.051) - 0.10).abs() < 1e-12);
        assert!((step_above(0.05) - 0.10).abs() < 1e-12);
        assert!((step_above(-0.07) - -0.05).abs() < 1e-12);
    }
}
