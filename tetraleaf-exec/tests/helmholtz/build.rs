//! The operator's tables built once (Helmholtz CH4.0h): an operator from
//! `HelmholtzOperator::unbuilt` with the settings, then a pool and `build_tables`, has the
//! tables of `with_block_budget`'s operator with the same settings applied after it (the
//! chain of setters that rebuilds its tables, on a pool of the same threads, so that the
//! scratch sets compare), level for level and bit for bit through the
//! pair bodies; its RCR tables store M2L blocks only on the levels whose M2L is RCR's, and
//! the unbuilt operator reports the degrees given and no table. And (Helmholtz CH2.2a) the
//! operator's RCR tables in a store: the operator without one bit for bit, cold, warm and at
//! k doubled, and the parts reused per level.

use std::sync::Arc;

use tetraleaf_exec::geometry::Domain;
use tetraleaf_exec::helmholtz::{
    ChargeKind, Densities, HelmholtzM2l, HelmholtzOperator, Wavenumber,
};
use tetraleaf_math::Layout;
use tetraleaf_octree::{MortonKey, PhysicalBox, morton};
use tetraleaf_plan::operator::{FmmSizes, PairOperator};
use tetraleaf_tables::geometry::{M2L_OFFSET_COUNT, m2l_offsets};
use tetraleaf_tables::helmholtz::DEFAULT_BLOCK_BUDGET;
use tetraleaf_tables::helmholtz::store::TableStore;

use super::common::Rng;

fn domain() -> Domain {
    Domain::new(&PhysicalBox::new([0.0, 0.0, 0.0, 1.0, 1.0, 1.0])).unwrap()
}

/// A source and a target on `level` whose offset is the V-list offset d.
fn pair(level: usize, d: usize) -> (MortonKey, MortonKey) {
    let offset = m2l_offsets()[d];
    let target = offset.map(|x| x.max(0) as usize);
    let source: [usize; 3] = core::array::from_fn(|a| (target[a] as i64 - offset[a]) as usize);
    let key = |i: [usize; 3]| morton::from_index_and_level(i, level);
    (key(source), key(target))
}

/// The dense budget of the cases: Batched's, Dense's and Classes' tables on the two finest
/// levels (p = 6 and 4), so that the debug run stays short.
const DENSE_BUDGET: usize = 64 << 20;

/// The settings of a case, applied in the builder's order.
fn settings(
    op: HelmholtzOperator<f64>,
    m2l: HelmholtzM2l,
    tolerance: Option<f64>,
) -> HelmholtzOperator<f64> {
    op.with_m2l_max_degree(Some(10))
        .with_m2l_plane_wave_min_degree(Some(11))
        .with_plane_wave_tolerance(tolerance)
        .with_m2l(m2l, DENSE_BUDGET)
}

/// Every M2L pair body on every level, every seventh offset, on random multipoles; then
/// every M2M and L2L pair body of one child per level.
fn pair_outputs(op: &mut HelmholtzOperator<f64>) -> Vec<f64> {
    let mut rng = Rng(11);
    let mut out = Vec::new();
    for level in 2..op.nlevels() {
        let multipole: Vec<f64> = (0..op.multipole_size(level)).map(|_| rng.next()).collect();
        for d in (0..M2L_OFFSET_COUNT).step_by(7) {
            let (source, target) = pair(level, d);
            let mut local = vec![0.0; op.local_size(level)];
            op.m2l(source, target, d, &multipole, &mut local);
            out.extend(local);
        }
    }
    for level in 1..op.nlevels() {
        let child = morton::from_index_and_level([1, 0, 1], level);
        let parent = morton::parent(child).expect("level ≥ 1");
        let multipole: Vec<f64> = (0..op.multipole_size(level)).map(|_| rng.next()).collect();
        let mut up = vec![0.0; op.multipole_size(level - 1)];
        op.m2m(
            child,
            parent,
            morton::child_index(child),
            &multipole,
            &mut up,
        );
        out.extend(up);
        let local: Vec<f64> = (0..op.local_size(level - 1)).map(|_| rng.next()).collect();
        let mut down = vec![0.0; op.local_size(level)];
        op.l2l(parent, child, morton::child_index(child), &local, &mut down);
        out.extend(down);
    }
    out
}

#[test]
fn the_unbuilt_operator_builds_the_settings_tables() {
    // k̂₀ = 15; below 10⁻⁷ Auto raises its levels with waves to 32.
    let k = Wavenumber::new(30.0, 0.0).unwrap();
    let degrees = [0, 0, 14, 9, 6, 4];
    let pool = Arc::new(
        rayon::ThreadPoolBuilder::new()
            .num_threads(3)
            .build()
            .unwrap(),
    );
    let cases = [
        (HelmholtzM2l::Rotation, None),
        (HelmholtzM2l::Auto, Some(1e-6)),
        (HelmholtzM2l::Auto, Some(5e-8)),
        (HelmholtzM2l::PlaneWave, Some(1e-6)),
        (HelmholtzM2l::Batched, None),
        (HelmholtzM2l::Dense, None),
        (HelmholtzM2l::Classes, None),
    ];
    for (m2l, tolerance) in cases {
        let what = format!("{m2l:?} at {tolerance:?}");
        let densities = Densities::new(k.kind(), ChargeKind::Complex, 1);
        let mut chain = settings(
            HelmholtzOperator::with_block_budget(
                k,
                densities,
                &domain(),
                &degrees,
                false,
                8,
                DEFAULT_BLOCK_BUDGET,
            )
            .unwrap()
            .with_pool(Arc::clone(&pool)),
            m2l,
            tolerance,
        );
        let unbuilt = settings(
            HelmholtzOperator::unbuilt(
                k,
                densities,
                &domain(),
                &degrees,
                false,
                8,
                DEFAULT_BLOCK_BUDGET,
            )
            .unwrap(),
            m2l,
            tolerance,
        );
        // Unbuilt: the degrees given, no M2L table, nothing timed.
        for (l, &p) in degrees.iter().enumerate() {
            assert_eq!(unbuilt.degree(l), p, "{what}");
            assert!(!unbuilt.tables().has_m2l(l), "{what}");
        }
        assert!(!unbuilt.tables().has_translations(), "{what}");
        assert_eq!(unbuilt.table_timings().rcr.as_nanos(), 0, "{what}");
        let mut once = unbuilt.with_pool(Arc::clone(&pool)).build_tables();
        assert_eq!(once.nlevels(), chain.nlevels());
        for l in 0..degrees.len() {
            assert_eq!(once.degree(l), chain.degree(l), "{what}, level {l}");
            assert_eq!(
                once.translation_degree(l),
                chain.translation_degree(l),
                "{what}, level {l}"
            );
            assert_eq!(once.m2l_on(l), chain.m2l_on(l), "{what}, level {l}");
            let stored = once.tables().m2l_blocks_stored(l);
            assert_eq!(
                stored,
                chain.tables().m2l_blocks_stored(l),
                "{what}, level {l}"
            );
            // Blocks only where the M2L is RCR's, and there within the budget (these fit).
            assert_eq!(
                stored,
                l >= 2 && once.m2l_on(l) == HelmholtzM2l::Rotation,
                "{what}, level {l}"
            );
        }
        // Each M2L of its own runs on some level (so its tables were built).
        if m2l != HelmholtzM2l::Auto {
            assert!((2..degrees.len()).any(|l| once.m2l_on(l) == m2l), "{what}");
        }
        // And RCR on the levels Batched's and Dense's budget leaves out (Auto and PlaneWave
        // may give every level waves or Batched, and the class tables fit every level).
        if matches!(m2l, HelmholtzM2l::Batched | HelmholtzM2l::Dense) {
            assert!(
                (2..degrees.len()).any(|l| once.m2l_on(l) == HelmholtzM2l::Rotation),
                "{what}: some level on RCR"
            );
        }
        assert_eq!(once.tables().bytes(), chain.tables().bytes(), "{what}");
        assert_eq!(once.dense_bytes(), chain.dense_bytes(), "{what}");
        assert_eq!(once.class_bytes(), chain.class_bytes(), "{what}");
        assert_eq!(once.batched_bytes(), chain.batched_bytes(), "{what}");
        assert_eq!(once.plane_wave_bytes(), chain.plane_wave_bytes(), "{what}");
        let (a, b) = (pair_outputs(&mut once), pair_outputs(&mut chain));
        assert_eq!(a.len(), b.len());
        for (i, (x, y)) in a.iter().zip(&b).enumerate() {
            assert!(x == y, "{what}: value {i}: {x:e} against {y:e}");
        }
        let n = Layout::new(degrees[2]).len();
        assert!(a.len() > n, "{what}");
    }
}

/// An operator from `with_block_budget` is built: `build_tables` changes nothing.
#[test]
fn build_tables_on_a_built_operator_does_nothing() {
    let k = Wavenumber::new(4.0, 0.5).unwrap();
    let densities = Densities::new(k.kind(), ChargeKind::Complex, 1);
    let op =
        HelmholtzOperator::<f64>::new(k, densities, &domain(), &[0, 0, 6, 5], false, 8).unwrap();
    let before = op.table_timings();
    let op = op.build_tables();
    assert_eq!(op.table_timings(), before);
    assert!(op.tables().has_translations());
}

/// An operator of `k`'s kind and `degrees`, under `m2l` at `tolerance`, its tables built on
/// `pool` in `store` (none: without one).
fn built(
    k: Wavenumber,
    degrees: &[usize],
    m2l: HelmholtzM2l,
    tolerance: Option<f64>,
    pool: &Arc<rayon::ThreadPool>,
    store: Option<&Arc<TableStore<f64>>>,
) -> HelmholtzOperator<f64> {
    let densities = Densities::new(k.kind(), ChargeKind::Complex, 1);
    let op = HelmholtzOperator::unbuilt(
        k,
        densities,
        &domain(),
        degrees,
        false,
        8,
        DEFAULT_BLOCK_BUDGET,
    )
    .unwrap()
    .with_table_store_option(store.cloned());
    settings(op, m2l, tolerance)
        .with_pool(Arc::clone(pool))
        .build_tables()
}

/// The operator's RCR tables in a store (Helmholtz CH2.2a): cold and warm, every pair body,
/// level and byte count the operator's without a store, bit for bit; warm, every part from
/// the store; and at k doubled on the same domain, with the degrees moved down one level,
/// every RCR part of the levels below the new coarse one taken from the store.
#[test]
fn the_operator_builds_its_tables_in_a_store() {
    let pool = Arc::new(
        rayon::ThreadPoolBuilder::new()
            .num_threads(3)
            .build()
            .unwrap(),
    );
    let cases = [
        (30.0, 2.0, HelmholtzM2l::Rotation, None),
        (30.0, 0.0, HelmholtzM2l::Auto, Some(1e-6)),
        (0.0, 30.0, HelmholtzM2l::Rotation, None),
    ];
    let degrees = [0, 0, 14, 9, 6, 4];
    for (re, im, m2l, tolerance) in cases {
        let what = format!("k = {re} + {im}i, {m2l:?} at {tolerance:?}");
        let k = Wavenumber::new(re, im).unwrap();
        let mut plain = built(k, &degrees, m2l, tolerance, &pool, None);
        assert_eq!(plain.table_build().reused(), 0, "{what}");
        assert!(plain.table_store().is_none());
        let expected = pair_outputs(&mut plain);
        let store = Arc::new(TableStore::<f64>::unbounded());
        let mut cold = built(k, &degrees, m2l, tolerance, &pool, Some(&store));
        let mut warm = built(k, &degrees, m2l, tolerance, &pool, Some(&store));
        for op in [&mut cold, &mut warm] {
            assert_eq!(op.nlevels(), plain.nlevels());
            for l in 0..degrees.len() {
                assert_eq!(op.degree(l), plain.degree(l), "{what}");
                assert_eq!(op.m2l_on(l), plain.m2l_on(l), "{what}");
                assert_eq!(
                    op.tables().m2l_blocks_stored(l),
                    plain.tables().m2l_blocks_stored(l),
                    "{what}"
                );
            }
            assert_eq!(op.tables().bytes(), plain.tables().bytes(), "{what}");
            let got = pair_outputs(op);
            assert_eq!(got.len(), expected.len());
            assert!(
                got.iter()
                    .zip(&expected)
                    .all(|(a, b)| a.to_bits() == b.to_bits()),
                "{what}: the pair bodies with a store"
            );
        }
        let parts = plain.table_build().built();
        assert!(cold.table_build().built() > 0, "{what}");
        assert_eq!(
            cold.table_build().built() + cold.table_build().reused(),
            parts,
            "{what}"
        );
        assert_eq!(warm.table_build().built(), 0, "{what}: warm");
        assert_eq!(warm.table_build().reused(), parts, "{what}: warm");
        assert_eq!(warm.table_build().store_bytes, Some(store.bytes()));
        // The octave: k doubled, the degrees one level down.
        if m2l == HelmholtzM2l::Rotation {
            let double = Wavenumber::new(2.0 * re, 2.0 * im).unwrap();
            let finer = [0, 0, 20, 14, 9, 6, 4];
            let mut octave = built(double, &finer, m2l, tolerance, &pool, Some(&store));
            let mut fresh = built(double, &finer, m2l, tolerance, &pool, None);
            let report = octave.table_build();
            for l in 3..finer.len() {
                assert_eq!(report.levels[l].built, 0, "{what}: level {l}: {report:?}");
                assert!(report.levels[l].reused > 0, "{what}: level {l}");
            }
            assert!(report.levels[2].built > 0, "{what}: the new coarse level");
            let (a, b) = (pair_outputs(&mut octave), pair_outputs(&mut fresh));
            assert!(
                a.iter().zip(&b).all(|(x, y)| x.to_bits() == y.to_bits()),
                "{what}: the octave's pair bodies"
            );
        }
    }
}

/// `TableStores` (Helmholtz CH2.2a): one store per precision, each with the cap, shared by the
/// handle's clones and equal only to them.
#[test]
fn table_stores_hold_a_store_per_precision() {
    use tetraleaf_exec::helmholtz::TableStores;
    let stores = TableStores::new(1 << 20);
    assert_eq!(stores.of::<f32>().cap(), 1 << 20);
    assert_eq!(stores.of::<f64>().cap(), 1 << 20);
    let clone = stores.clone();
    assert!(Arc::ptr_eq(stores.of::<f32>(), clone.of::<f32>()));
    assert!(Arc::ptr_eq(stores.of::<f64>(), clone.of::<f64>()));
    assert_eq!(stores, clone);
    assert_ne!(stores, TableStores::new(1 << 20));
    // An f32 operator builds its tables in the f32 store, the f64 store untouched.
    let k = Wavenumber::new(30.0, 2.0).unwrap();
    let densities = Densities::new(k.kind(), ChargeKind::Complex, 1);
    let op = HelmholtzOperator::<f32>::unbuilt(
        k,
        densities,
        &domain(),
        &[0, 0, 8, 6],
        false,
        8,
        DEFAULT_BLOCK_BUDGET,
    )
    .unwrap()
    .with_table_store(Arc::clone(stores.of::<f32>()))
    .build_tables();
    assert!(op.table_build().built() > 0);
    assert!(stores.of::<f32>().entries() > 0);
    assert_eq!(stores.of::<f64>().entries(), 0);
}
