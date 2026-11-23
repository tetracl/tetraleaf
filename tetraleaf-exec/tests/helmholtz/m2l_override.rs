//! The operator's M2L overrides (Helmholtz CH4.4b; `HelmholtzOperator::with_m2l_overrides`), a
//! measurement setting: a level's degree or M2L in place of the setting's, every other level as
//! the setting gives it. Under Auto with waves from every degree and the time rule at free form
//! costs, k = 40 on the unit cube, at 10⁻⁴, 10⁻⁷ and 10⁻⁸ (levels with waves, levels the form
//! took from them, raised levels): RCR on each level (no waves, no form, no batched tables, the
//! degree as run), its classes off (all, or some), a degree with the setting's M2L (the operator
//! given that degree), and a degree with RCR; the overrides refused. The FMM with overrides is
//! tested in `tests/helmholtz_fmm.rs` (`check_m2l_override`).

use tetraleaf_exec::geometry::{Domain, radius};
use tetraleaf_exec::helmholtz::{
    ChargeKind, Densities, DiagonalCosts, HelmholtzDiagonal, HelmholtzM2l, HelmholtzOperator,
    M2lOverride, OverrideError, OverrideM2l, Wavenumber, degree,
};
use tetraleaf_octree::PhysicalBox;

const LEVELS: usize = 5;

/// Costs under which the form takes every level the bound reaches (as `diagonal.rs`'s).
const FREE: DiagonalCosts = DiagonalCosts {
    s2f: [0.0, 0.0],
    f2r: [0.0, 0.0],
    pair: [0.0, 0.0],
    rotation: [1.0, 1.0],
    plane_wave: [1e12, 1e12],
};

fn domain() -> Domain {
    Domain::new(&PhysicalBox::new([0.0, 0.0, 0.0, 1.0, 1.0, 1.0])).unwrap()
}

/// The default rule's degrees at ε, k = 40.
fn degrees(k: &Wavenumber, epsilon: f64) -> Vec<usize> {
    (0..LEVELS)
        .map(|l| match l {
            0 | 1 => 0,
            _ => degree(k, 3f64.sqrt() * radius(l, &domain()), epsilon),
        })
        .collect()
}

/// The operator of `p` under Auto (waves from degree 0) and the time rule at free form costs,
/// with `overrides`.
fn operator(
    k: Wavenumber,
    p: &[usize],
    epsilon: f64,
    overrides: Vec<M2lOverride>,
) -> HelmholtzOperator<f64> {
    let densities = Densities::new(k.kind(), ChargeKind::Complex, 1);
    HelmholtzOperator::<f64>::new(k, densities, &domain(), p, false, 64)
        .unwrap()
        .with_plane_wave_tolerance(Some(epsilon))
        .with_m2l_plane_wave_min_degree(Some(0))
        .with_m2l(HelmholtzM2l::Auto, 1 << 30)
        .with_diagonal_costs(Some(FREE))
        .with_diagonal(HelmholtzDiagonal::Auto)
        .with_m2l_overrides(overrides)
}

/// A level's M2L, classes, degree and translation degree.
fn level_of(
    op: &HelmholtzOperator<f64>,
    l: usize,
) -> (HelmholtzM2l, Option<Vec<usize>>, usize, usize) {
    (
        op.m2l_on(l),
        op.diagonal_classes(l),
        op.degree(l),
        op.translation_degree(l),
    )
}

/// Every override's level as named, every other level as the setting gives it.
#[test]
fn the_overrides_per_level() {
    let k = Wavenumber::new(40.0, 0.0).unwrap();
    let (mut waves, mut form, mut raised) = (0, 0, 0);
    for epsilon in [1e-4, 1e-7, 1e-8] {
        let p = degrees(&k, epsilon);
        let base = operator(k, &p, epsilon, Vec::new());
        let same_elsewhere = |what: &str, op: &HelmholtzOperator<f64>, level: usize| {
            for l in (2..LEVELS).filter(|&l| l != level) {
                assert_eq!(level_of(op, l), level_of(&base, l), "{what}: level {l}");
            }
        };
        for l in 2..LEVELS {
            let (m2l, classes, p_run, translation) = level_of(&base, l);
            println!(
                "ε {epsilon:e}, level {l}: {m2l:?}, classes {classes:?}, p {p_run} ({})",
                p[l]
            );
            waves += usize::from(m2l == HelmholtzM2l::PlaneWave);
            form += usize::from(classes.is_some());
            raised += usize::from(p_run > p[l]);
            let what = format!("ε {epsilon:e}, level {l}");
            // RCR: no waves, no form, the degree and the translations' degree as run.
            let rcr = operator(
                k,
                &p,
                epsilon,
                vec![M2lOverride {
                    level: l,
                    degree: None,
                    m2l: OverrideM2l::Rcr,
                }],
            );
            assert_eq!(
                level_of(&rcr, l),
                (HelmholtzM2l::Rotation, None, p_run, translation),
                "{what}: RCR"
            );
            assert!(rcr.plane_wave_precision(l).is_none(), "{what}: RCR");
            same_elsewhere(&format!("{what}: RCR"), &rcr, l);
            // Every class off: RCR's tables; some classes off: the rest.
            let all = operator(
                k,
                &p,
                epsilon,
                vec![M2lOverride {
                    level: l,
                    degree: None,
                    m2l: OverrideM2l::RcrClasses((0..34).collect()),
                }],
            );
            if let Some(classes) = &classes {
                assert_eq!(
                    level_of(&all, l),
                    level_of(&rcr, l),
                    "{what}: every class off"
                );
                let off: Vec<usize> = classes.iter().step_by(3).copied().collect();
                let some = operator(
                    k,
                    &p,
                    epsilon,
                    vec![M2lOverride {
                        level: l,
                        degree: None,
                        m2l: OverrideM2l::RcrClasses(off.clone()),
                    }],
                );
                let rest: Vec<usize> = classes
                    .iter()
                    .copied()
                    .filter(|c| !off.contains(c))
                    .collect();
                assert_eq!(
                    level_of(&some, l),
                    (m2l, Some(rest), p_run, translation),
                    "{what}: classes {off:?} off"
                );
                same_elsewhere(&format!("{what}: some classes off"), &some, l);
            } else {
                // No form there: the classes change nothing.
                assert_eq!(level_of(&all, l), level_of(&base, l), "{what}: no form");
            }
            // A degree with the setting's M2L: the operator given that degree.
            let q = p[l] + 3;
            let mut given = p.clone();
            given[l] = q;
            let named = operator(k, &given, epsilon, Vec::new());
            let setting = operator(
                k,
                &p,
                epsilon,
                vec![M2lOverride {
                    level: l,
                    degree: Some(q),
                    m2l: OverrideM2l::Setting,
                }],
            );
            for j in 2..LEVELS {
                assert_eq!(
                    level_of(&setting, j),
                    level_of(&named, j),
                    "{what}: degree {q}, the setting, level {j}"
                );
            }
            // A degree with RCR: RCR at that degree, the translations too.
            let rcr_q = operator(
                k,
                &p,
                epsilon,
                vec![M2lOverride {
                    level: l,
                    degree: Some(q),
                    m2l: OverrideM2l::Rcr,
                }],
            );
            assert_eq!(
                level_of(&rcr_q, l),
                (HelmholtzM2l::Rotation, None, q, q),
                "{what}: RCR at degree {q}"
            );
            same_elsewhere(&format!("{what}: RCR at degree {q}"), &rcr_q, l);
        }
        // No override, and the setter before the others: the setting's tables.
        let early = {
            let densities = Densities::new(k.kind(), ChargeKind::Complex, 1);
            HelmholtzOperator::<f64>::new(k, densities, &domain(), &p, false, 64)
                .unwrap()
                .with_m2l_overrides(vec![M2lOverride {
                    level: 2,
                    degree: None,
                    m2l: OverrideM2l::Rcr,
                }])
                .with_plane_wave_tolerance(Some(epsilon))
                .with_m2l_plane_wave_min_degree(Some(0))
                .with_m2l(HelmholtzM2l::Auto, 1 << 30)
                .with_diagonal_costs(Some(FREE))
                .with_diagonal(HelmholtzDiagonal::Auto)
        };
        let late = operator(
            k,
            &p,
            epsilon,
            vec![M2lOverride {
                level: 2,
                degree: None,
                m2l: OverrideM2l::Rcr,
            }],
        );
        for l in 0..LEVELS {
            assert_eq!(
                level_of(&early, l),
                level_of(&late, l),
                "the order, level {l}"
            );
        }
        // A level beyond the tree: ignored.
        let beyond = operator(
            k,
            &p,
            epsilon,
            vec![M2lOverride {
                level: LEVELS + 2,
                degree: Some(3),
                m2l: OverrideM2l::Rcr,
            }],
        );
        for l in 0..LEVELS {
            assert_eq!(
                level_of(&beyond, l),
                level_of(&base, l),
                "beyond, level {l}"
            );
        }
    }
    assert!(
        waves > 0 && form > 0 && raised > 0,
        "the cases cover waves ({waves}), the form ({form}) and raised levels ({raised})"
    );
}

/// The refused lists.
#[test]
fn the_overrides_refused() {
    let o = |level, m2l| M2lOverride {
        level,
        degree: None,
        m2l,
    };
    assert_eq!(M2lOverride::check(&[]), Ok(()));
    assert_eq!(
        M2lOverride::check(&[o(1, OverrideM2l::Rcr)]),
        Err(OverrideError::Level(1))
    );
    assert_eq!(
        M2lOverride::check(&[o(2, OverrideM2l::Rcr), o(2, OverrideM2l::Setting)]),
        Err(OverrideError::Repeated(2))
    );
    assert_eq!(
        M2lOverride::check(&[o(3, OverrideM2l::RcrClasses(vec![0, 34]))]),
        Err(OverrideError::Class {
            class: 34,
            count: 34
        })
    );
    let refused = std::panic::catch_unwind(|| {
        let k = Wavenumber::new(40.0, 0.0).unwrap();
        operator(k, &degrees(&k, 1e-4), 1e-4, vec![o(0, OverrideM2l::Rcr)])
    });
    assert!(refused.is_err(), "the operator panics on a refused list");
}
