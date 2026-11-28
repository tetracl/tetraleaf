//! The RCR tables' parts in a store (Helmholtz CH2.2a, `tetraleaf_tables::helmholtz::store`):
//! tables built with a store, cold and warm, on any map, are the tables built without one bit
//! for bit through every application (every octant's M2M and L2L, every offset's M2L, every
//! level; one plane for Yukawa, two otherwise; compared by `to_bits`); an octave shift of k̂₀
//! reuses every part of the finer levels; the blocks of one degree are built once; the cap
//! evicts the least recently used parts, the same on every run; builds on two threads share a
//! store.

use std::sync::Arc;

use tetraleaf_math::{Complex, Layout, RealScalar};
use tetraleaf_tables::geometry::{M2L_OFFSET_COUNT, OCTANT_COUNT};
use tetraleaf_tables::helmholtz::store::TableStore;
use tetraleaf_tables::helmholtz::{DEFAULT_BLOCK_BUDGET, HelmholtzTables, TableChoice};
use tetraleaf_tables::parts::{Parts, Serial};

use super::input;
use super::parts::Threads;

/// k̂₀ as (Re, Im), its degrees and 2k̂₀'s, one level deeper.
type Octave = ((f64, f64), [usize; 5], [usize; 6]);

/// k̂₀ as (Re, Im) and the degrees.
type Case = ((f64, f64), &'static [usize]);

/// Complex and real k above the projection threshold on level 2, complex k below it, and
/// Yukawa (Re k̂₀ = 0).
const CASES: [Case; 4] = [
    ((120.0, 6.0), &[0, 0, 12, 10, 8]),
    ((120.0, 0.0), &[0, 0, 12, 10, 8]),
    ((6.0, 0.5), &[0, 0, 9, 8, 7]),
    ((0.0, 6.0), &[0, 0, 9, 9, 9]),
];

fn k(case: Case) -> Complex<f64> {
    Complex::new(case.0.0, case.0.1)
}

/// Every application of `tables` on fixed inputs (`parts.rs`' `applications` in T, with the
/// tables' own planes), as the bits of each value widened to f64 (exact).
fn applications<T: RealScalar>(tables: &HelmholtzTables<T>) -> Vec<u64> {
    let planes = tables.local_planes();
    let mut scratch = tables.scratch();
    let mut out = Vec::new();
    let len = |l: usize| planes * Layout::new(tables.m2l_degree(l)).len();
    let typed = |p: usize, seed: u64| -> Vec<T> {
        input(p, planes, seed)
            .into_iter()
            .map(T::from_f64)
            .collect()
    };
    let mut push =
        |values: &[T]| out.extend(values.iter().map(|v| RealScalar::to_f64(*v).to_bits()));
    for level in 0..tables.nlevels() {
        if level + 1 < tables.nlevels() && tables.has_translations() {
            for o in 0..OCTANT_COUNT {
                let child = typed(tables.m2l_degree(level + 1), o as u64);
                let mut parent = vec![T::zero(); len(level)];
                tables.m2m(level, o, &child, &mut parent, planes, 1, &mut scratch);
                push(&parent);
                let parent = typed(tables.m2l_degree(level), 100 + o as u64);
                let mut child = vec![T::zero(); len(level + 1)];
                tables.l2l(level + 1, o, &parent, &mut child, planes, 1, &mut scratch);
                push(&child);
            }
        }
    }
    for level in (2..tables.nlevels()).filter(|&l| tables.has_m2l(l)) {
        let q = tables.m2l_degree(level);
        for d in 0..M2L_OFFSET_COUNT {
            let x = typed(q, 1000 + d as u64);
            let mut y = vec![T::zero(); len(level)];
            tables.m2l(level, d, &x, planes, &mut y, 1, &mut scratch);
            push(&y);
        }
    }
    out
}

/// The tables of `degrees` built without a store.
fn plain<T: RealScalar>(k_root: Complex<f64>, degrees: &[usize]) -> HelmholtzTables<T> {
    HelmholtzTables::build_by(
        k_root,
        degrees,
        &TableChoice::new(DEFAULT_BLOCK_BUDGET),
        &Serial,
    )
}

/// The tables of `degrees` built with `store`, on `parts`.
fn stored<T: RealScalar>(
    k_root: Complex<f64>,
    degrees: &[usize],
    parts: &impl Parts,
    store: &TableStore<T>,
) -> HelmholtzTables<T> {
    HelmholtzTables::build_in(
        k_root,
        degrees,
        &TableChoice::new(DEFAULT_BLOCK_BUDGET),
        parts,
        store,
    )
}

#[test]
fn a_build_in_a_store_is_the_build_without_one() {
    for case in CASES {
        let expected = plain::<f64>(k(case), case.1);
        let bits = applications(&expected);
        let parts = expected.build_report().built();
        assert_eq!(expected.build_report().reused(), 0, "{case:?}: no store");
        assert_eq!(expected.build_report().store_bytes, None);
        let store = TableStore::<f64>::unbounded();
        let cold = stored(k(case), case.1, &Threads, &store);
        assert_eq!(applications(&cold), bits, "{case:?}: cold");
        assert_eq!(cold.bytes(), expected.bytes(), "{case:?}: cold bytes");
        let report = cold.build_report();
        // Cold: the distinct parts computed once (blocks of one degree and angle shared by
        // levels), the rest reused within the build.
        assert_eq!(report.built() + report.reused(), parts, "{case:?}: cold");
        assert!(report.built() <= parts && report.built() > 0, "{case:?}");
        assert_eq!(report.store_bytes, Some(store.bytes()));
        assert_eq!(
            store.entries(),
            report.built(),
            "{case:?}: one entry per part"
        );
        // Warm: every part from the store, on another map.
        let warm = stored(k(case), case.1, &Serial, &store);
        assert_eq!(applications(&warm), bits, "{case:?}: warm");
        assert_eq!(
            warm.build_report().built(),
            0,
            "{case:?}: warm builds nothing"
        );
        assert_eq!(warm.build_report().reused(), parts, "{case:?}: warm");
        assert_eq!(warm.build_report().seconds(), Default::default());
        let levels: Vec<usize> = (0..warm.nlevels()).map(|l| warm.level_bytes(l)).collect();
        let reported: Vec<usize> = warm.build_report().levels.iter().map(|l| l.bytes).collect();
        assert_eq!(reported, levels, "{case:?}: the bytes per level");
        assert_eq!(levels.iter().sum::<usize>(), warm.bytes());
    }
}

#[test]
fn the_store_holds_its_parts_once() {
    // Two levels at one degree: their M2L blocks are one set, built once and shared.
    let k_root = Complex::new(0.0, 6.0);
    let degrees = [0, 0, 9, 9, 9];
    let store = TableStore::<f64>::unbounded();
    let tables = stored(k_root, &degrees, &Serial, &store);
    let report = tables.build_report();
    let plain = plain::<f64>(k_root, &degrees);
    for l in 3..5 {
        assert!(tables.m2l_blocks_stored(l));
        // The blocks (49 angles) and the octant blocks of degree 9 are level 2's.
        assert!(
            report.levels[l].reused >= 49,
            "level {l}: {:?}",
            report.levels[l]
        );
    }
    // The store counts each part once; the tables count each level whole.
    assert!(store.bytes() < tables.bytes());
    assert_eq!(tables.bytes(), plain.bytes());
}

#[test]
fn an_octave_reuses_every_finer_level() {
    // k̂₀ and 2k̂₀ on one domain: level l + 1 at 2k̂₀ has level l's k̂ at k̂₀ exactly, so with
    // the degrees moved down one level the second build computes only its new coarse level's
    // parts (level 2's M2L, the translations and blocks of the levels above level 3).
    let cases: [Octave; 3] = [
        ((60.0, 3.0), [0, 0, 12, 10, 8], [0, 0, 14, 12, 10, 8]),
        ((60.0, 0.0), [0, 0, 12, 10, 8], [0, 0, 14, 12, 10, 8]),
        ((0.0, 3.0), [0, 0, 9, 9, 9], [0, 0, 9, 9, 9, 9]),
    ];
    for ((re, im), coarse, fine) in cases {
        let k_root = Complex::new(re, im);
        let k_double = Complex::new(2.0 * re, 2.0 * im);
        let store = TableStore::<f64>::unbounded();
        let first = stored(k_root, &coarse, &Threads, &store);
        let second = stored(k_double, &fine, &Threads, &store);
        let expected = plain::<f64>(k_double, &fine);
        assert_eq!(
            applications(&second),
            applications(&expected),
            "{k_double:?}"
        );
        let report = second.build_report();
        for l in 3..fine.len() {
            assert_eq!(
                report.levels[l].built, 0,
                "{k_double:?}: level {l} reuses: {:?}",
                report.levels
            );
            assert_eq!(
                report.levels[l].reused,
                first.build_report().levels[l - 1].built
                    + first.build_report().levels[l - 1].reused,
                "{k_double:?}: level {l} holds level {}'s parts",
                l - 1
            );
        }
        // Level 2's M2L tables (15 distances) are new at 2k̂₀ unless Yukawa's plateau degree
        // matches level 1 at k̂₀, which has no M2L: new in every case.
        assert!(
            report.levels[2].built >= 15,
            "{k_double:?}: {:?}",
            report.levels
        );
    }
}

#[test]
fn set_m2l_blocks_in_a_store_is_set_m2l_blocks() {
    for case in CASES {
        let n = case.1.len();
        let none = vec![false; n];
        let choice = TableChoice {
            block_levels: Some(&none),
            ..TableChoice::new(DEFAULT_BLOCK_BUDGET)
        };
        let all: Vec<bool> = (0..n).map(|l| l >= 2).collect();
        let mut expected = HelmholtzTables::<f64>::build_by(k(case), case.1, &choice, &Serial);
        expected.set_m2l_blocks(&all, DEFAULT_BLOCK_BUDGET, &Serial);
        let store = TableStore::<f64>::unbounded();
        let mut tables = HelmholtzTables::build_in(k(case), case.1, &choice, &Threads, &store);
        let before = tables.build_report().built();
        tables.set_m2l_blocks_in(&all, DEFAULT_BLOCK_BUDGET, &Threads, &store);
        assert_eq!(applications(&tables), applications(&expected), "{case:?}");
        assert_eq!(tables.bytes(), expected.bytes());
        assert!(
            tables.build_report().built() > before,
            "{case:?}: blocks added"
        );
        assert_eq!(tables.build_report().store_bytes, Some(store.bytes()));
        // A second operator's blocks come from the store.
        let mut again = HelmholtzTables::build_in(k(case), case.1, &choice, &Serial, &store);
        again.set_m2l_blocks_in(&all, DEFAULT_BLOCK_BUDGET, &Serial, &store);
        assert_eq!(again.build_report().built(), 0, "{case:?}");
        assert_eq!(applications(&again), applications(&expected), "{case:?}");
    }
}

#[test]
fn the_cap_keeps_the_most_recent_parts() {
    let a = CASES[0];
    let b = ((240.0, 12.0), &[0, 0, 0, 12, 10, 8][..]);
    // A store of cap 0 keeps nothing and changes no bit.
    let none = TableStore::<f64>::new(0);
    let tables = stored(k(a), a.1, &Serial, &none);
    assert_eq!((none.bytes(), none.entries()), (0, 0));
    assert_eq!(
        applications(&tables),
        applications(&plain::<f64>(k(a), a.1))
    );
    let again = stored(k(a), a.1, &Serial, &none);
    assert_eq!(
        again.build_report().reused(),
        tables.build_report().reused()
    );
    // B alone, to size the cap.
    let alone = TableStore::<f64>::unbounded();
    stored(k(b), b.1, &Serial, &alone);
    // A, then B, under B's bytes: the parts B used are kept (those it shares with A
    // included), A's others evicted, the same on every run.
    let mut runs = Vec::new();
    for _ in 0..2 {
        let store = TableStore::<f64>::new(alone.bytes());
        stored(k(a), a.1, &Serial, &store);
        let tables = stored(k(b), b.1, &Serial, &store);
        assert!(tables.build_report().reused() > 0, "B shares parts with A");
        assert_eq!(store.bytes(), alone.bytes());
        assert_eq!(store.entries(), alone.entries());
        runs.push((store.bytes(), store.entries()));
        // Everything B needs is there.
        let warm = stored(k(b), b.1, &Serial, &store);
        assert_eq!(warm.build_report().built(), 0);
        assert_eq!(applications(&warm), applications(&plain::<f64>(k(b), b.1)));
    }
    assert_eq!(runs[0], runs[1]);
    // `clear` empties the store; the tables keep their parts.
    let store = TableStore::<f64>::unbounded();
    let tables = stored(k(a), a.1, &Serial, &store);
    store.clear();
    assert_eq!((store.bytes(), store.entries()), (0, 0));
    assert_eq!(
        applications(&tables),
        applications(&plain::<f64>(k(a), a.1))
    );
}

#[test]
fn builds_on_two_threads_share_a_store() {
    let store = Arc::new(TableStore::<f64>::unbounded());
    let (a, b) = (CASES[0], CASES[2]);
    let built: Vec<Vec<u64>> = std::thread::scope(|scope| {
        let handles: Vec<_> = [a, b, a, b]
            .into_iter()
            .map(|case| {
                let store = Arc::clone(&store);
                scope.spawn(move || applications(&stored(k(case), case.1, &Serial, &*store)))
            })
            .collect();
        handles
            .into_iter()
            .map(|h| h.join().expect("a build"))
            .collect()
    });
    assert_eq!(built[0], applications(&plain::<f64>(k(a), a.1)));
    assert_eq!(built[1], applications(&plain::<f64>(k(b), b.1)));
    assert_eq!(built[2], built[0]);
    assert_eq!(built[3], built[1]);
}

#[test]
fn a_store_of_f32_tables() {
    for case in [CASES[0], CASES[3]] {
        let expected = applications(&plain::<f32>(k(case), case.1));
        let store = TableStore::<f32>::unbounded();
        let cold = stored(k(case), case.1, &Threads, &store);
        assert_eq!(applications(&cold), expected, "{case:?}: cold");
        let warm = stored(k(case), case.1, &Serial, &store);
        assert_eq!(warm.build_report().built(), 0);
        assert_eq!(applications(&warm), expected, "{case:?}: warm");
    }
}

#[test]
fn another_k_at_the_same_degrees_builds_its_own() {
    // Complex and real k, and two real k, at the same degrees in one store: only the y-blocks
    // (no k̂ in their key) are shared; every coaxial table is the second k's own.
    for (a, b) in [(CASES[0], CASES[1]), (CASES[1], ((121.0, 0.0), CASES[1].1))] {
        let store = TableStore::<f64>::unbounded();
        let first = stored(k(a), a.1, &Serial, &store);
        let second = stored(k(b), b.1, &Serial, &store);
        assert_eq!(
            applications(&second),
            applications(&plain::<f64>(k(b), b.1)),
            "{b:?} after {a:?}"
        );
        assert_eq!(applications(&first), applications(&plain::<f64>(k(a), a.1)));
        for l in 2..b.1.len() {
            // 15 M2L tables, and the translations below every level but the finest.
            let coaxial = 15 + if l + 1 < b.1.len() { 2 } else { 0 };
            assert!(
                second.build_report().levels[l].built >= coaxial,
                "{b:?}: level {l}: {:?}",
                second.build_report().levels[l]
            );
        }
    }
}

#[test]
fn another_degree_at_the_same_k_builds_its_own() {
    // One k̂, level 2's degree changed, then its M2L degree raised above its translations'
    // (`TableChoice::m2l_degrees`, Helmholtz CH6.1h): only the parts of the changed degrees
    // are built again, and each build is `build_by`'s bit for bit.
    let k_root = Complex::new(120.0, 6.0);
    let store = TableStore::<f64>::unbounded();
    stored(k_root, &[0, 0, 12, 10, 8], &Serial, &store);
    let degrees = [0, 0, 13, 10, 8];
    let second = stored(k_root, &degrees, &Serial, &store);
    assert_eq!(
        applications(&second),
        applications(&plain::<f64>(k_root, &degrees))
    );
    let report = second.build_report();
    assert!(
        report.levels[2].built >= 15,
        "level 2's M2L at 13: {report:?}"
    );
    assert_eq!(report.levels[3].built, 0, "level 3 unchanged: {report:?}");
    assert_eq!(report.levels[4].built, 0, "level 4 unchanged: {report:?}");
    let raised = [0, 0, 14, 11, 8];
    let choice = TableChoice {
        m2l_degrees: Some(&raised),
        ..TableChoice::new(DEFAULT_BLOCK_BUDGET)
    };
    let expected = HelmholtzTables::<f64>::build_by(k_root, &degrees, &choice, &Serial);
    let third = HelmholtzTables::build_in(k_root, &degrees, &choice, &Threads, &store);
    assert_eq!(applications(&third), applications(&expected));
    let report = third.build_report();
    assert!(
        report.levels[2].built >= 15 && report.levels[3].built >= 15,
        "{report:?}"
    );
    assert_eq!(report.levels[4].built, 0, "level 4 unchanged: {report:?}");
}

#[test]
fn a_cap_below_one_build_keeps_the_same_parts_on_every_run() {
    // A cap of half a build's parts: the eviction within one build (all its parts of one age,
    // dropped in key order) keeps the same parts on every run, and changes no bit.
    let case = CASES[0];
    let full = TableStore::<f64>::unbounded();
    stored(k(case), case.1, &Serial, &full);
    let kept: Vec<(usize, usize)> = (0..2)
        .map(|_| {
            let store = TableStore::<f64>::new(full.bytes() / 2);
            let tables = stored(k(case), case.1, &Threads, &store);
            assert_eq!(
                applications(&tables),
                applications(&plain::<f64>(k(case), case.1))
            );
            assert!(store.bytes() <= full.bytes() / 2 && store.entries() > 0);
            // The parts kept are found by the next build.
            let warm = stored(k(case), case.1, &Serial, &store);
            assert!(warm.build_report().reused() > tables.build_report().reused());
            (store.bytes(), store.entries())
        })
        .collect();
    assert_eq!(kept[0], kept[1]);
}
