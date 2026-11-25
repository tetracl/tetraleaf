//! The degree rules: S7's closed form (Helmholtz CH3.1a; spike S7,
//! `spikes/truncation/results/summary.md`, "Per-level degree p*") at its grid points, never
//! below p*; the calibrated default (Helmholtz CH3.3, its high-frequency branch refitted in
//! CH4.3c2 and on the worst density in CH4.3c3, one degree added after the cap in CH4.8b;
//! `tetraleaf-validate/results/helmholtz-ch3-3.md`, `helmholtz-ch4-3c2.md`, `helmholtz-ch4-3c3.md`,
//! `helmholtz-ch4-8a.md`, `helmholtz-ch4-8b.md`) never above S7's, monotone, at its fitted points
//! of both densities, CH4.3c3's rule plus its margin and CH4.8a's degrees as run; CH3.3's,
//! CH4.3c2's and CH4.3c3's rules kept (CH4.3c3's never above its cap's); and `Degrees::resolve`.
//! Exact integers.

use tetraleaf_exec::geometry::{Domain, radius};
use tetraleaf_exec::helmholtz::{DegreeError, DegreeRule, Degrees, Wavenumber, degree, s7_degree};
use tetraleaf_octree::PhysicalBox;

const EPSILONS: [f64; 8] = [1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10];

/// Rows of S7's real-k table at k a = 2^j: (k a, [(p*, closed form) per ε]).
const ROWS: [(f64, [(usize, usize); 8]); 7] = [
    (
        1.0,
        [
            (12, 16),
            (17, 22),
            (23, 28),
            (29, 34),
            (35, 40),
            (42, 46),
            (49, 52),
            (54, 58),
        ],
    ),
    (
        4.0,
        [
            (11, 16),
            (15, 22),
            (20, 28),
            (26, 34),
            (31, 40),
            (38, 46),
            (44, 52),
            (50, 58),
        ],
    ),
    (
        16.0,
        [
            (24, 26),
            (26, 28),
            (28, 29),
            (30, 34),
            (33, 40),
            (38, 46),
            (44, 52),
            (50, 58),
        ],
    ),
    (
        32.0,
        [
            (41, 44),
            (44, 46),
            (46, 49),
            (49, 51),
            (51, 53),
            (53, 55),
            (55, 56),
            (57, 58),
        ],
    ),
    (
        64.0,
        [
            (75, 79),
            (79, 82),
            (82, 85),
            (85, 88),
            (88, 90),
            (90, 92),
            (93, 95),
            (95, 97),
        ],
    ),
    (
        128.0,
        [
            (141, 147),
            (146, 151),
            (150, 154),
            (153, 158),
            (157, 161),
            (161, 164),
            (164, 167),
            (167, 169),
        ],
    ),
    (
        256.0,
        [
            (272, 280),
            (278, 284),
            (284, 289),
            (288, 293),
            (293, 297),
            (297, 301),
            (300, 304),
            (304, 308),
        ],
    ),
];

#[test]
fn the_closed_form_at_s7_grid_points() {
    for (ka, row) in ROWS {
        let k = Wavenumber::new(ka, 0.0).unwrap();
        for (epsilon, (p_star, closed)) in EPSILONS.into_iter().zip(row) {
            let p = s7_degree(&k, 1.0, epsilon);
            assert_eq!(p, closed, "k a = {ka}, ε = {epsilon:e}");
            assert!(
                p >= p_star,
                "k a = {ka}, ε = {epsilon:e}: {p} < p* = {p_star}"
            );
        }
    }
}

/// The largest measured degree at ε/2 per ε (10⁻³ … 10⁻⁸) at the levels' x of CH3.3's
/// scans (`tetraleaf-validate/results/helmholtz-ch3-3.md` §2 and §3) and CH4.3c2's
/// (`tetraleaf-validate/results/helmholtz-ch4-3c2.md`): (x, [p per ε]); the plateau scans'
/// envelope at their largest x, then the level scans. CH4.3c2's rows take the largest degree
/// over the workloads at one x and the smallest of their x (the cube at N = 10⁴ has a
/// slightly smaller box).
const MEASURED: [(f64, [usize; 6]); 16] = [
    (2.1651, [6, 9, 12, 15, 20, 24]),
    (1.6238, [6, 8, 11, 14, 18, 22]),
    (2.7063, [7, 8, 11, 14, 17, 20]),
    (3.2474, [7, 9, 11, 14, 17, 22]),
    (5.4127, [10, 11, 13, 15, 18, 21]),
    (6.4948, [10, 12, 14, 16, 19, 22]),
    (6.7862, [8, 10, 11, 13, 14, 16]),
    (10.8253, [15, 17, 19, 21, 23, 24]),
    (13.5725, [13, 16, 18, 20, 22, 24]),
    (21.6507, [26, 28, 31, 33, 35, 37]),
    (43.3018, [47, 51, 54, 57, 60, 62]),
    (64.9528, [69, 73, 77, 80, 83, 86]),
    (27.0637, [32, 34, 37, 39, 41, 44]),
    (54.1273, [58, 62, 65, 68, 71, 74]),
    (108.2478, [112, 117, 121, 126, 129, 132]),
    (216.4956, [218, 225, 230, 235, 240, 244]),
];

/// The largest measured degree at ε/2 per ε (10⁻³ … 10⁻⁸) with W1's constant density
/// (Helmholtz CH4.3c3, `tetraleaf-validate/results/helmholtz-ch4-3c3.md`: level scans at kD = 250,
/// 500 and 1000), at the smallest x of the scans of one level width.
const MEASURED_CONSTANT: [(f64, [usize; 6]); 6] = [
    (6.7659, [11, 13, 15, 16, 19, 22]),
    (13.5318, [18, 21, 23, 25, 26, 28]),
    (27.0636, [32, 35, 37, 39, 42, 44]),
    (54.1273, [60, 63, 67, 70, 72, 75]),
    (108.2547, [114, 119, 123, 127, 131, 134]),
    (216.5096, [225, 230, 236, 240, 245, 248]),
];

#[test]
fn the_calibrated_rule_covers_the_measured_degrees() {
    for (x, row) in MEASURED.into_iter().chain(MEASURED_CONSTANT) {
        for (epsilon, p) in EPSILONS.into_iter().zip(row) {
            let rule = DegreeRule::CALIBRATED.degree(x, epsilon);
            assert!(rule >= p, "x = {x}, ε = {epsilon:e}: {rule} < {p}");
        }
    }
}

/// CH4.3c2's rule covered the random charges but missed four of the constant density's degrees
/// by one (Helmholtz CH4.3c3): x = 216.5 at 10⁻³ and 10⁻⁷, x = 13.5 at 10⁻⁸ (W1 at kD = 250 and
/// 500; the third tree, kD = 1000, needed 27).
#[test]
fn ch4_3c2s_rule_misses_the_constant_density() {
    let mut misses = Vec::new();
    for (x, row) in MEASURED.into_iter().chain(MEASURED_CONSTANT) {
        for (epsilon, p) in EPSILONS.into_iter().zip(row) {
            let rule = DegreeRule::CH4_3C2.degree(x, epsilon);
            if rule < p {
                misses.push((x, epsilon, p - rule));
            }
        }
    }
    assert_eq!(
        misses,
        [(13.5318, 1e-8, 1), (216.5096, 1e-3, 1), (216.5096, 1e-7, 1)]
    );
}

#[test]
fn the_calibrated_rule_is_at_most_s7s_and_monotone() {
    let xs: Vec<f64> = (-30..=33).map(|j| 2f64.powf(j as f64 / 4.0)).collect();
    let epsilons: Vec<f64> = (30..=100).map(|j| 10f64.powf(-j as f64 / 10.0)).collect();
    for &x in &xs {
        let k = Wavenumber::new(x, 0.0).unwrap();
        for (i, &epsilon) in epsilons.iter().enumerate() {
            let p = degree(&k, 1.0, epsilon);
            assert_eq!(p, DegreeRule::CALIBRATED.degree(x, epsilon));
            assert!(
                p <= s7_degree(&k, 1.0, epsilon),
                "x = {x}, ε = {epsilon:e}: {p} above S7's"
            );
            if i > 0 {
                assert!(p >= degree(&k, 1.0, epsilons[i - 1]), "not monotone in ε");
            }
            assert!(p >= degree(&k, 0.9, epsilon), "not monotone in x");
        }
    }
}

/// The refit's cap (Helmholtz CH4.3c3): CH4.3c3's rule, the default before CH4.8b's margin, at or
/// below its cap's rule,
/// max(⌈0.30 L² + 0.25 L + 2.60⌉, ⌈x + 0.82 (ln 1/ε)^{2/3} x^{1/3}⌉), everywhere on
/// 4 ≤ x ≤ 512 and 10⁻⁸ ≤ ε ≤ 10⁻³ by tenths of a decade; above CH3.3's rule by one at 39 of
/// those 5,763 points (x = 8–21 and ε ≤ 10^{−6.7}, where the constant density needs it), and
/// below it at x = 54–216 by 6–9 degrees at 10⁻³ and 1–3 at 10⁻⁸. CH4.3c2's rule at most CH3.3's
/// everywhere (its construction), equal to it at x ≤ 7.
#[test]
fn ch4_3c3s_rule_is_at_most_its_cap() {
    let cap = DegreeRule::gd(0.30, 0.25, 2.60, 0.82);
    let mut above_ch3_3 = 0;
    for k in 0..=112 {
        let x = 4.0 * 2f64.powf(f64::from(k) / 16.0);
        for j in 30..=80 {
            let epsilon = 10f64.powf(-f64::from(j) / 10.0);
            let (p, ch3_3) = (
                DegreeRule::CH4_3C3.degree(x, epsilon),
                DegreeRule::CH3_3.degree(x, epsilon),
            );
            assert!(
                p <= cap.degree(x, epsilon),
                "x = {x}, ε = {epsilon:e}: {p} above the cap's"
            );
            if p > ch3_3 {
                assert_eq!(p, ch3_3 + 1, "x = {x}, ε = {epsilon:e}");
                assert!(
                    (8.0..=21.0).contains(&x) && epsilon <= 2e-7,
                    "x = {x}, ε = {epsilon:e}"
                );
                above_ch3_3 += 1;
            }
            let before = DegreeRule::CH4_3C2.degree(x, epsilon);
            assert!(
                before <= ch3_3,
                "x = {x}, ε = {epsilon:e}: CH4.3c2's {before} above CH3.3's {ch3_3}"
            );
            if x <= 7.0 {
                assert_eq!(before, ch3_3, "x = {x}, ε = {epsilon:e}");
            }
        }
    }
    assert_eq!(above_ch3_3, 39);
    let below = |x: f64, epsilon: f64| {
        DegreeRule::CH3_3.degree(x, epsilon) - DegreeRule::CH4_3C3.degree(x, epsilon)
    };
    assert_eq!([54.1, 108.3, 216.5].map(|x| below(x, 1e-3)), [6, 8, 9]);
    assert_eq!([54.1, 108.3, 216.5].map(|x| below(x, 1e-8)), [1, 1, 3]);
}

/// The excess is taken at least 0 (the refit's e(ε) < 0 above ε ≈ 2 · 10⁻³), so at any ε the
/// high-frequency branch is at least ⌈x⌉; the cap is CH3.3's branch with coefficient 0.82; the
/// default's branch is CH4.3c3's plus one.
#[test]
fn the_high_branch_is_never_below_x() {
    let rule = DegreeRule::CH4_3C3;
    let cap = DegreeRule::gd(0.30, 0.25, 2.60, 0.82);
    assert!(rule.excess(1e-2) < 0.0 && rule.excess(1e-3) > 0.0);
    for x in [10.0, 100.0, 216.5, 1000.0] {
        for epsilon in [0.5, 0.1, 1e-2, 2e-3, 1e-3] {
            let high = rule.high(x, epsilon);
            assert!(high >= x.ceil() as usize, "x = {x}, ε = {epsilon:e}");
            assert!(high <= cap.high(x, epsilon), "x = {x}, ε = {epsilon:e}");
            assert_eq!(DegreeRule::CALIBRATED.high(x, epsilon), high + 1);
        }
    }
    // At ε = 10⁻² and x = 216.5: ⌈216.5 + 0 + 2⌉ = 219 (unclamped 216.5 − 0.22 · 6.0 + 2, 218).
    assert_eq!(rule.high(216.5, 1e-2), 219);
    let uncapped = DegreeRule { cap: None, ..rule };
    // The cap binds at small x: at x = 5 and 10⁻⁸, ⌈5 + 4.89 · 1.71 + 2⌉ = 16 against the cap's
    // ⌈5 + 0.82 · 6.98 · 1.71⌉ = 15; at x = 13.5 it does not, 28 both (CH3.3's 27).
    assert_eq!((uncapped.high(5.0, 1e-8), rule.high(5.0, 1e-8)), (16, 15));
    assert_eq!((uncapped.high(13.5, 1e-8), rule.high(13.5, 1e-8)), (28, 28));
    assert_eq!(DegreeRule::CH3_3.high(13.5, 1e-8), 27);
}

/// CH3.3's rule is kept as [`DegreeRule::CH3_3`]: the degrees of the default before CH4.3c2.
#[test]
fn ch3_3s_rule_is_kept() {
    let rule = DegreeRule::CH3_3;
    assert_eq!(rule, DegreeRule::gd(0.30, 0.25, 2.60, 0.80));
    // CH3.3's formula as it stood before CH4.3c2, restated: max(⌈0.30 L² + 0.25 L + 2.60⌉,
    // ⌈x + 0.80 (ln 1/ε)^{2/3} x^{1/3}⌉), the ceilings forgiving 10⁻⁹; the same degrees on the
    // tests' grid.
    let before = |x: f64, epsilon: f64| {
        let l = (1.0 / epsilon).log10();
        let up = |v: f64| (v - 1e-9).ceil().max(0.0) as usize;
        up((0.30 * l + 0.25) * l + 2.60).max(up(x + 0.80
            * (1.0 / epsilon).ln().powf(2.0 / 3.0)
            * x.cbrt()))
    };
    for j in -30..=36 {
        let x = 2f64.powf(f64::from(j) / 4.0);
        for i in 30..=100 {
            let epsilon = 10f64.powf(-f64::from(i) / 10.0);
            assert_eq!(
                rule.degree(x, epsilon),
                before(x, epsilon),
                "x = {x}, ε = {epsilon:e}"
            );
        }
    }
    // kD = 100 on the unit cube at 10⁻⁶ (the module example): 35 and 22, CH4.3c2's too (its
    // cap binds there); CH4.3c3's 34 and 21, the default's 35 and 22.
    let x = 100.0 * 3f64.sqrt() / 8.0;
    assert_eq!(rule.degree(x, 1e-6), 35);
    assert_eq!(rule.degree(x / 2.0, 1e-6), 22);
    assert_eq!(DegreeRule::CH4_3C2.degree(x, 1e-6), 35);
    assert_eq!(DegreeRule::CH4_3C2.degree(x / 2.0, 1e-6), 22);
    assert_eq!(DegreeRule::CH4_3C3.degree(x, 1e-6), 34);
    assert_eq!(DegreeRule::CH4_3C3.degree(x / 2.0, 1e-6), 21);
    assert_eq!(DegreeRule::CALIBRATED.degree(x, 1e-6), 35);
    assert_eq!(DegreeRule::CALIBRATED.degree(x / 2.0, 1e-6), 22);
    // x = 216.5 at 10⁻³ and 10⁻⁸: CH3.3's 234 and 251, CH4.3c2's 224 and 248, CH4.3c3's 225 and
    // 248, the default's 226 and 249.
    let at = |r: &DegreeRule| [1e-3, 1e-8].map(|e| r.degree(216.5, e));
    assert_eq!(at(&rule), [234, 251]);
    assert_eq!(at(&DegreeRule::CH4_3C2), [224, 248]);
    assert_eq!(at(&DegreeRule::CH4_3C3), [225, 248]);
    assert_eq!(at(&DegreeRule::CALIBRATED), [226, 249]);
}

/// CH4.3c3's rule is kept as [`DegreeRule::CH4_3C3`] (Helmholtz CH4.8b): the default before, its
/// constants those of CH4.3c3's report; the default is it with one degree added to the
/// high-frequency branch after the cap, the plateau unchanged (CH4.8a's question 2, signed off on
/// #260), so at every x and ε it is max(plateau, CH4.3c3's branch + 1).
#[test]
fn the_default_is_ch4_3c3s_rule_with_one_degree_after_the_cap() {
    let rule = DegreeRule::CH4_3C3;
    assert_eq!(
        rule,
        DegreeRule {
            c: 2.00,
            a: -0.55,
            b: -4.66,
            d: 2.0,
            cap: Some(0.82),
            ..DegreeRule::CH3_3
        }
    );
    assert_eq!(rule.margin, 0);
    assert_eq!(DegreeRule::CALIBRATED, DegreeRule { margin: 1, ..rule });
    for k in 0..=160 {
        let x = 2f64.powf(f64::from(k) / 16.0 - 1.0);
        for j in 30..=80 {
            let epsilon = 10f64.powf(-f64::from(j) / 10.0);
            let p = DegreeRule::CALIBRATED.degree(x, epsilon);
            assert_eq!(
                DegreeRule::CALIBRATED.plateau(epsilon),
                rule.plateau(epsilon)
            );
            assert_eq!(
                p,
                rule.plateau(epsilon).max(rule.high(x, epsilon) + 1),
                "x = {x}, ε = {epsilon:e}"
            );
            let before = rule.degree(x, epsilon);
            assert!(p == before || p == before + 1, "x = {x}, ε = {epsilon:e}");
        }
    }
}

/// The degrees CH4.8a ran (`tetraleaf-validate/results/helmholtz-ch4-8a.md` §5 and §9, W1 and W2
/// on the unit cube at kD = 250 and 500): `h1` the default's, and the default before CH4.8b
/// CH4.3c3's rule's, on every level from 2 whose degree the waves did not raise to their resolving
/// degree p_w (28 at 10⁻⁷, 32 at 10⁻⁸, levels of x ≤ 13.5 there).
#[test]
fn the_degrees_ch4_8a_ran() {
    let domain = Domain::new(&PhysicalBox::new([0.0, 0.0, 0.0, 1.0, 1.0, 1.0])).unwrap();
    // W2 at kD = 500, levels 2–8 (x = 108.25, 54.13, …, 1.69); `None` where p_w ran.
    let h1: [(f64, [Option<usize>; 7]); 6] = [
        (1e-3, [116, 61, 33, 19, 12, 8, 7].map(Some)),
        (1e-4, [121, 65, 36, 22, 14, 10, 9].map(Some)),
        (1e-5, [125, 68, 39, 24, 16, 12, 12].map(Some)),
        (1e-6, [129, 71, 41, 26, 17, 15, 15].map(Some)),
        (
            1e-7,
            [Some(132), Some(74), Some(43), None, None, None, None],
        ),
        (
            1e-8,
            [Some(135), Some(76), Some(45), None, None, None, None],
        ),
    ];
    let before: [(f64, [Option<usize>; 7]); 6] = [
        (1e-3, [115, 60, 32, 18, 11, 7, 7].map(Some)),
        (1e-4, [120, 64, 35, 21, 13, 9, 9].map(Some)),
        (1e-5, [124, 67, 38, 23, 15, 12, 12].map(Some)),
        (1e-6, [128, 70, 40, 25, 16, 15, 15].map(Some)),
        (
            1e-7,
            [Some(131), Some(73), Some(42), None, None, None, None],
        ),
        (
            1e-8,
            [Some(134), Some(75), Some(44), None, None, None, None],
        ),
    ];
    for (kd, first) in [(500.0, 0), (250.0, 1)] {
        let k = Wavenumber::new(kd, 0.0).unwrap();
        for (table, rule) in [
            (&h1, DegreeRule::CALIBRATED),
            (&before, DegreeRule::CH4_3C3),
        ] {
            for (epsilon, row) in table {
                let degrees = Degrees::Tolerance(*epsilon)
                    .resolve_with(&rule, &k, &domain, 9 - first)
                    .unwrap();
                for (l, ran) in (2..9 - first).zip(&row[first..]) {
                    if let Some(p) = ran {
                        assert_eq!(degrees[l], *p, "kD = {kd}, ε = {epsilon:e}, level {l}");
                    }
                }
            }
        }
        assert_eq!(
            Degrees::Tolerance(1e-6)
                .resolve(&k, &domain, 9 - first)
                .unwrap(),
            Degrees::Tolerance(1e-6)
                .resolve_with(&DegreeRule::CALIBRATED, &k, &domain, 9 - first)
                .unwrap()
        );
    }
}

/// `h1` as CH4.8a ran it (the default's constants with d = 3 and no cap, no margin) against the
/// default (the margin after the cap): the same degree on 4 ≤ x ≤ 512 and 10⁻⁸ ≤ ε ≤ 10⁻³ but where
/// CH4.3c3's cap binds, at x = 7.1–13.5 and ε ≤ 2 · 10⁻⁷, where `h1` is one degree more (CH4.8a's
/// report §1); there every level of CH4.8a's trees ran p_w above both.
#[test]
fn h1_as_run_against_the_default() {
    let h1 = DegreeRule {
        d: 3.0,
        cap: None,
        ..DegreeRule::CH4_3C3
    };
    let mut more = 0;
    for k in 0..=112 {
        let x = 4.0 * 2f64.powf(f64::from(k) / 16.0);
        for j in 30..=80 {
            let epsilon = 10f64.powf(-f64::from(j) / 10.0);
            let (ran, p) = (
                h1.degree(x, epsilon),
                DegreeRule::CALIBRATED.degree(x, epsilon),
            );
            if ran != p {
                assert_eq!(ran, p + 1, "x = {x}, ε = {epsilon:e}");
                assert!(
                    (7.1..=13.5).contains(&x) && epsilon <= 2e-7,
                    "x = {x}, ε = {epsilon:e}"
                );
                assert!(p <= 28, "x = {x}, ε = {epsilon:e}: {p} above p_w");
                more += 1;
            }
        }
    }
    assert!(more > 0);
}

/// CH4.3c2's rule is kept as [`DegreeRule::CH4_3C2`] (Helmholtz CH4.3c3): the default before.
#[test]
fn ch4_3c2s_rule_is_kept() {
    let rule = DegreeRule::CH4_3C2;
    assert_eq!(
        rule,
        DegreeRule {
            c: 2.30,
            a: -0.75,
            b: -5.73,
            d: 5.0,
            cap: Some(0.80),
            ..DegreeRule::CH3_3
        }
    );
    // Its degrees at CH4.3c2's measured x (its report §2, "the rule").
    let at = |x: f64| [1e-3, 1e-8].map(|e| rule.degree(x, e));
    assert_eq!(at(27.0637), [34, 44]);
    assert_eq!(at(54.1273), [61, 76]);
    assert_eq!(at(108.2478), [115, 134]);
    assert_eq!(at(216.4956), [224, 248]);
}

#[test]
fn degree_rule_shapes() {
    // S7's closed form is γ = 0, α = 6, β = −2, C = 1 (up to the forgiving ceilings).
    let s7 = DegreeRule::gd(0.0, 6.0, -2.0, 1.0);
    for (ka, row) in ROWS {
        for (epsilon, (_, closed)) in EPSILONS.into_iter().zip(row) {
            assert_eq!(
                s7.degree(ka, epsilon),
                closed,
                "k a = {ka}, ε = {epsilon:e}"
            );
        }
    }
    // The ceilings forgive 10⁻⁹: 2.5 · 6 + 1 = 16 exactly at ε = 10⁻⁶.
    let rule = DegreeRule::gd(0.0, 2.5, 1.0, 0.0);
    assert_eq!(rule.plateau(1e-6), 16);
    assert_eq!(rule.high(0.0, 1e-6), 0);
    assert_eq!(rule.degree(20.0, 1e-6), 20);
}

#[test]
fn complex_k_takes_its_real_part() {
    let real = Wavenumber::new(30.0, 0.0).unwrap();
    let complex = Wavenumber::new(30.0, 3.0).unwrap();
    for epsilon in EPSILONS {
        assert_eq!(degree(&complex, 0.7, epsilon), degree(&real, 0.7, epsilon));
    }
}

#[test]
fn resolve() {
    let domain = Domain::new(&PhysicalBox::new([-1.0, -1.0, -1.0, 3.0, 3.0, 3.0])).unwrap();
    let k = Wavenumber::new(50.0, 0.0).unwrap();
    // Side 4: r₂ = 0.5, x₂ = 50 · 0.866 = 43.3. CH4.3c3's branch, e(10⁻⁸) = 2.00 · 18.42^{2/3}
    // − 4.4 − 4.66 = 4.89, gives 43.3 + 4.89 · 43.3^{1/3} + 2 = 62.48, below the cap
    // 43.3 + 0.82 · 6.975 · 3.512 = 63.39: 63, the default's 64; level 3 38 and 39; level 4
    // (x = 10.8) the cap 10.83 + 0.82 · 6.975 · 2.213 = 23.49: 24 and the default's 25 above the
    // plateau ⌈0.30 · 64 + 0.25 · 8 + 2.60⌉ = 24; level 5 the plateau.
    let degrees = Degrees::Tolerance(1e-8).resolve(&k, &domain, 6).unwrap();
    assert_eq!(degrees, [0, 0, 64, 39, 25, 24]);
    // CH3.3's rule, named: CH4.3c3's degrees here; at 10⁻³ the refit is lower on level 2.
    let ch3_3 = Degrees::Tolerance(1e-8)
        .resolve_with(&DegreeRule::CH3_3, &k, &domain, 6)
        .unwrap();
    assert_eq!(ch3_3, [0, 0, 63, 38, 24, 24]);
    assert_eq!(
        Degrees::Tolerance(1e-8)
            .resolve_with(&DegreeRule::CH4_3C3, &k, &domain, 6)
            .unwrap(),
        ch3_3
    );
    assert_eq!(
        Degrees::Tolerance(1e-8)
            .resolve_with(&DegreeRule::CALIBRATED, &k, &domain, 6)
            .unwrap(),
        degrees
    );
    let at_1e3 = |rule: &DegreeRule| {
        Degrees::Tolerance(1e-3)
            .resolve_with(rule, &k, &domain, 6)
            .unwrap()
    };
    assert_eq!(at_1e3(&DegreeRule::CALIBRATED), [0, 0, 50, 28, 16, 11]);
    assert_eq!(at_1e3(&DegreeRule::CH4_3C3), [0, 0, 49, 27, 15, 10]);
    assert_eq!(at_1e3(&DegreeRule::CH4_3C2), [0, 0, 50, 28, 17, 11]);
    assert_eq!(at_1e3(&DegreeRule::CH3_3), [0, 0, 54, 30, 18, 11]);
    // S7's closed form: 68 on level 2, its floor ⌈6 · 8⌉ − 2 = 46 below.
    let s7: Vec<usize> = (2..6)
        .map(|l| s7_degree(&k, 3f64.sqrt() * radius(l, &domain), 1e-8))
        .collect();
    assert_eq!(s7, [68, 46, 46, 46]);

    let named = Degrees::PerLevel(vec![1, 2, 3, 4, 5]);
    assert_eq!(named.resolve(&k, &domain, 3).unwrap(), [1, 2, 3]);
    assert_eq!(
        named.resolve(&k, &domain, 6),
        Err(DegreeError::TooFewLevels {
            given: 5,
            levels: 6
        })
    );
    for epsilon in [0.0, 1.0, -1e-6, f64::NAN] {
        assert!(matches!(
            Degrees::Tolerance(epsilon).resolve(&k, &domain, 3),
            Err(DegreeError::Tolerance(_))
        ));
    }
}
