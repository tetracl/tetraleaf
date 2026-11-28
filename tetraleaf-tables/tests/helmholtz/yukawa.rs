//! Yukawa's real tables (Helmholtz CH2.1a, `tetraleaf_tables::helmholtz` at k̂₀ = iκ̂₀): M2M
//! and L2L on the 8 octants and M2L on the 316 offsets against
//! `tetraleaf_ref::helmholtz::rcr`'s `_yukawa` forms, rectangular and raised degrees, from
//! κ̂ a ≈ 10⁻³ to strong damping; Dense against Rotation; vectors, stored blocks,
//! contributions and the build in parts bit for bit; the refusals.
//!
//! Error measure: the FMM-weighted error of `main.rs` (`compare`, `weights`), each one-plane
//! expansion widened by a zero imaginary plane; |M̂^Y_n| = |M̂_n| and |L̂^Y_n| = |L̂_n| (CONVENTIONS
//! §3.15, "Yukawa's real form"), so the weights of the complex expansions at k = iκ apply.
//! Within [`BOUND`] times [`conditioning`], the bound of the real and complex tables, also
//! under strong damping (`under_strong_damping`), where an M2L's output is e^{−κ(|t| − 2a)}
//! times its input's size and the coaxial sum cancels as `rcr`'s does (its error relative to
//! the input printed beside). Bit for bit: `==` per value.

use tetraleaf_math::{Complex, Layout};
use tetraleaf_ref::Frame;
use tetraleaf_ref::helmholtz::{Level, Workspace, rcr};
use tetraleaf_tables::geometry::{M2L_OFFSET_COUNT, OCTANT_COUNT, m2l_offsets, octant_direction};
use tetraleaf_tables::helmholtz::batched::BatchedM2l;
use tetraleaf_tables::helmholtz::classes::ClassM2l;
use tetraleaf_tables::helmholtz::dense::DenseM2l;
use tetraleaf_tables::helmholtz::{
    DEFAULT_BLOCK_BUDGET, HelmholtzTables, TableChoice, add_m2l_contribution,
};
use tetraleaf_tables::parts::Serial;

use super::parts::Threads;
use super::{
    BOUND, assert_same, compare, conditioning, embed, input, prefix, weighted_norm, weights,
    widened,
};

/// κ̂₀ of the tests, k̂₀ = iκ̂₀: κ̂ a from about 10⁻³ on level 4 to about 17 on level 2.
const KAPPA_ROOTS: [f64; 3] = [0.01, 2.0, 40.0];

fn yukawa(kappa: f64) -> Complex<f64> {
    Complex::new(0.0, kappa)
}

/// Yukawa M2L on `offsets` of `level` against `rcr::m2l_yukawa`: the largest FMM-weighted
/// error relative to the output, and relative to the input (the multipole's weighted norm),
/// the locals started at `initial` (zero where the output is far below it).
fn check_m2l(
    tables: &HelmholtzTables<f64>,
    level: usize,
    offsets: impl Iterator<Item = usize>,
    initial: f64,
) -> (f64, f64) {
    let p = tables.m2l_degree(level);
    let len = Layout::new(p).len();
    let k_level = Level::new(tables.k_hat(level), 1.0);
    let (w_local, w_multipole) = (weights(&k_level, false, p), weights(&k_level, true, p));
    let from = Frame::new([0.0; 3], 1.0);
    let mut ws = Workspace::new(p);
    let mut scratch = tables.scratch();
    let (mut output, mut input_relative) = (0.0f64, 0.0f64);
    for d in offsets {
        let to = Frame::new(m2l_offsets()[d].map(|x| 2.0 * x as f64), 1.0);
        let table = rcr::m2l_table(p, p, &k_level, &from, &to);
        let x = input(p, 1, d as u64);
        let mut got = vec![initial; len];
        tables.m2l(level, d, &x, 1, &mut got, 1, &mut scratch);
        let mut expected = vec![initial; len];
        rcr::m2l_yukawa(
            p,
            p,
            &k_level,
            &from,
            &to,
            &table,
            &mut ws,
            &x,
            &mut expected,
        );
        let (got, expected) = (widened(&got), widened(&expected));
        output = output.max(compare(&got, &expected, initial, &w_local));
        let diff: Vec<f64> = got.iter().zip(&expected).map(|(a, b)| a - b).collect();
        input_relative = input_relative
            .max(weighted_norm(&w_local, &diff) / weighted_norm(&w_multipole, &widened(&x)));
    }
    (output, input_relative)
}

/// Yukawa M2M and L2L of every octant between `level` and `level` + 1 against `rcr`'s
/// `_yukawa` forms, as `main.rs`'s `check_shifts`: the child's degrees p_{l+1} within its
/// M2L degree's expansion (NaN above them for M2M, kept values for L2L), L2L against ½
/// `rcr`'s. The largest FMM-weighted error.
fn check_shifts(tables: &HelmholtzTables<f64>, level: usize) -> f64 {
    let (p_parent, p_child) = (tables.m2l_degree(level), tables.degree(level + 1));
    let q_child = tables.m2l_degree(level + 1);
    let (len_parent, len_child) = (Layout::new(p_parent).len(), Layout::new(p_child).len());
    let k = tables.k_hat(level + 1);
    let child_level = Level::new(k, 1.0);
    let parent_level = Level::new(k, 2.0);
    let parent_frame = Frame::new([0.0; 3], 2.0);
    let m2m_table = rcr::m2m_table(p_child, p_parent, &child_level);
    let l2l_table = rcr::l2l_table(p_parent, p_child, &child_level);
    let mut ws = Workspace::new(p_parent.max(p_child));
    let mut scratch = tables.scratch();
    let mut worst = 0.0f64;
    for o in 0..OCTANT_COUNT {
        let child_frame = Frame::new(octant_direction(o).map(|s| s as f64), 1.0);
        let x = input(p_child, 1, o as u64);
        let mut got = vec![0.25; len_parent];
        let block = embed(&x, 1, p_child, q_child, f64::NAN);
        tables.m2m(level, o, &block, &mut got, 1, 1, &mut scratch);
        let mut expected = vec![0.25; len_parent];
        rcr::m2m_yukawa(
            p_child,
            p_parent,
            &child_level,
            &parent_level,
            &child_frame,
            &parent_frame,
            &m2m_table,
            &mut ws,
            &x,
            &mut expected,
        );
        worst = worst.max(compare(
            &widened(&got),
            &widened(&expected),
            0.25,
            &weights(&parent_level, true, p_parent),
        ));

        let x = input(p_parent, 1, 8 + o as u64);
        let mut block = embed(&vec![0.0; len_child], 1, p_child, q_child, 7.0);
        tables.l2l(level + 1, o, &x, &mut block, 1, 1, &mut scratch);
        let got = prefix(&block, 1, p_child, q_child);
        assert_same(
            &block,
            &embed(&got, 1, p_child, q_child, 7.0),
            "Yukawa L2L beyond the child's degree",
        );
        let mut expected = vec![0.0; len_child];
        rcr::l2l_yukawa(
            p_parent,
            p_child,
            &parent_level,
            &child_level,
            &parent_frame,
            &child_frame,
            &l2l_table,
            &mut ws,
            &x,
            &mut expected,
        );
        let half: Vec<f64> = expected.iter().map(|v| 0.5 * v).collect();
        worst = worst.max(compare(
            &widened(&got),
            &widened(&half),
            0.0,
            &weights(&child_level, false, p_child),
        ));
    }
    worst
}

#[test]
fn yukawa_tables_are_real_with_one_local_plane() {
    let tables = HelmholtzTables::<f64>::build(yukawa(2.0), &[0, 0, 4], DEFAULT_BLOCK_BUDGET);
    assert!(tables.is_yukawa() && tables.real_shifts());
    assert_eq!(tables.local_planes(), 1);
    assert_eq!(tables.m2l_contribution_len(2, 3), 3 * 25);
    for k in [Complex::new(2.0, 0.0), Complex::new(2.0, 0.5)] {
        let tables = HelmholtzTables::<f64>::build(k, &[0, 0, 4], DEFAULT_BLOCK_BUDGET);
        assert!(!tables.is_yukawa());
        assert_eq!(tables.local_planes(), 2);
        assert_eq!(tables.m2l_contribution_len(2, 3), 3 * 2 * 25);
    }
}

#[test]
fn yukawa_m2l_on_every_offset_matches_rcr() {
    for kappa in KAPPA_ROOTS {
        let tables =
            HelmholtzTables::<f64>::build(yukawa(kappa), &[0, 0, 7, 6, 5], DEFAULT_BLOCK_BUDGET);
        let worst = (2..5)
            .map(|level| {
                check_m2l(&tables, level, 0..M2L_OFFSET_COUNT, 0.5).0 / conditioning(&tables, level)
            })
            .fold(0.0, f64::max);
        println!(
            "κ̂₀ = {kappa}: Yukawa M2L against rcr, worst {worst:.2e} per unit of conditioning"
        );
        assert!(worst <= BOUND, "κ̂₀ = {kappa}: {worst:e}");
    }
}

/// Rectangular degrees between every pair of levels, and the M2L degrees raised above the
/// translations' (Helmholtz CH6.1h): M2M and L2L on the first degrees, M2L bit for bit the
/// tables built at the raised degrees.
#[test]
fn yukawa_m2m_and_l2l_on_every_octant_match_rcr() {
    let (degrees, raised) = ([2, 6, 9, 4, 7], [2, 8, 12, 4, 10]);
    for kappa in KAPPA_ROOTS {
        let k = yukawa(kappa);
        for q in [degrees, raised] {
            let tables = HelmholtzTables::<f64>::build_with_m2l_degrees(
                k,
                &degrees,
                &q,
                DEFAULT_BLOCK_BUDGET,
            );
            let worst = (0..4)
                .map(|level| check_shifts(&tables, level) / conditioning(&tables, level))
                .fold(0.0, f64::max);
            println!(
                "κ̂₀ = {kappa}, M2L degrees {q:?}: Yukawa M2M and L2L against rcr, worst \
                 {worst:.2e} per unit of conditioning"
            );
            assert!(worst <= BOUND, "κ̂₀ = {kappa}, {q:?}: {worst:e}");
        }
        let tables = HelmholtzTables::<f64>::build_with_m2l_degrees(k, &degrees, &raised, 0);
        let at_raised = HelmholtzTables::<f64>::build(k, &raised, 0);
        let (mut a, mut b) = (tables.scratch(), at_raised.scratch());
        for (level, &q) in raised.iter().enumerate().skip(2) {
            for d in (0..M2L_OFFSET_COUNT).step_by(5) {
                let x = input(q, 1, d as u64);
                let mut got = vec![0.5; Layout::new(q).len()];
                let mut expected = got.clone();
                tables.m2l(level, d, &x, 1, &mut got, 1, &mut a);
                at_raised.m2l(level, d, &x, 1, &mut expected, 1, &mut b);
                assert_same(&got, &expected, &format!("raised M2L, level {level}, {d}"));
            }
        }
    }
}

/// Strong damping: κ̂ a ≈ 69 and 173 on level 2, where an M2L's output is about
/// e^{−0.54 κ̂ a} of its input (10⁻¹⁶ and 10⁻⁴¹) and the coaxial sum cancels. The tables
/// cancel as `rcr` does: relative to the output within the bound (relative to the input
/// printed); M2M and L2L of the pair (2, 3) within the bound. The locals start at zero, which
/// a start of 0.5 would swamp.
#[test]
fn under_strong_damping() {
    for kappa in [160.0, 400.0] {
        let tables =
            HelmholtzTables::<f64>::build(yukawa(kappa), &[0, 0, 10, 8], DEFAULT_BLOCK_BUDGET);
        let c = conditioning(&tables, 2);
        let (output, input) = check_m2l(&tables, 2, 0..M2L_OFFSET_COUNT, 0.0);
        let (output, input) = (output / c, input / c);
        let shifts = check_shifts(&tables, 2) / c;
        println!(
            "κ̂₀ = {kappa}: Yukawa M2L {output:.2e} relative to the output, {input:.2e} to the \
             input; M2M and L2L {shifts:.2e} (per unit of conditioning)"
        );
        assert!(output <= BOUND, "κ̂₀ = {kappa}: M2L {output:e}");
        assert!(shifts <= BOUND, "κ̂₀ = {kappa}: M2M and L2L {shifts:e}");
    }
}

/// Several vectors against one, stored blocks against per-call, the contribution then its
/// add against M2L, and the build on a reversed three-thread map against the serial build:
/// bit for bit, through every application.
#[test]
fn yukawa_applications_bit_for_bit() {
    let k = yukawa(6.0);
    let (degrees, raised) = ([0, 3, 9, 8, 7], [0, 3, 10, 8, 9]);
    let stored = HelmholtzTables::<f64>::build_with_m2l_degrees(k, &degrees, &raised, usize::MAX);
    let per_call = HelmholtzTables::<f64>::build_with_m2l_degrees(k, &degrees, &raised, 0);
    let choice = TableChoice {
        m2l_degrees: Some(&raised),
        ..TableChoice::new(usize::MAX)
    };
    let threads = HelmholtzTables::<f64>::build_by(k, &degrees, &choice, &Threads);
    let serial = HelmholtzTables::<f64>::build_by(k, &degrees, &choice, &Serial);
    assert!(stored.m2l_blocks_stored(2) && !per_call.m2l_blocks_stored(2));
    assert_eq!(threads.bytes(), serial.bytes());
    let all = applications(&stored);
    assert_same(&all, &applications(&per_call), "per-call blocks");
    assert_same(&all, &applications(&threads), "built in parts");
    assert_same(&all, &applications(&serial), "build_by on Serial");

    let mut scratch = stored.scratch();
    for (level, &q) in raised.iter().enumerate().skip(2) {
        let len = Layout::new(q).len();
        for d in (0..M2L_OFFSET_COUNT).step_by(7) {
            let x = [input(q, 1, d as u64), input(q, 1, 1000 + d as u64)].concat();
            let mut two = vec![0.5; 2 * len];
            stored.m2l(level, d, &x, 1, &mut two, 2, &mut scratch);
            for v in 0..2 {
                let mut one = vec![0.5; len];
                stored.m2l(level, d, &x[v * len..][..len], 1, &mut one, 1, &mut scratch);
                assert_same(&two[v * len..][..len], &one, "vector of two");
            }
            let mut contribution = vec![f64::NAN; stored.m2l_contribution_len(level, 2)];
            stored.m2l_contribution(level, d, &x, 1, &mut contribution, 2, &mut scratch);
            let mut added = vec![0.5; 2 * len];
            add_m2l_contribution(&contribution, &mut added);
            assert_same(&added, &two, "contribution then its add");
        }
    }
}

/// Every application of Yukawa tables on fixed inputs, one plane each: per level with a child
/// level the M2M and L2L of every octant, then per level from 2 the M2L of every offset.
fn applications(tables: &HelmholtzTables<f64>) -> Vec<f64> {
    let mut scratch = tables.scratch();
    let mut out = Vec::new();
    let len = |l: usize| Layout::new(tables.m2l_degree(l)).len();
    for level in 0..tables.nlevels() {
        if level + 1 < tables.nlevels() {
            for o in 0..OCTANT_COUNT {
                let child = input(tables.m2l_degree(level + 1), 1, o as u64);
                let mut parent = vec![0.0; len(level)];
                tables.m2m(level, o, &child, &mut parent, 1, 1, &mut scratch);
                out.extend(parent);
                let parent = input(tables.m2l_degree(level), 1, 100 + o as u64);
                let mut child = vec![0.0; len(level + 1)];
                tables.l2l(level + 1, o, &parent, &mut child, 1, 1, &mut scratch);
                out.extend(child);
            }
        }
        if level >= 2 {
            for d in 0..M2L_OFFSET_COUNT {
                let multipole = input(tables.m2l_degree(level), 1, 200 + d as u64);
                let mut local = vec![0.0; len(level)];
                tables.m2l(level, d, &multipole, 1, &mut local, 1, &mut scratch);
                out.extend(local);
            }
        }
    }
    out
}

/// Dense against Rotation on every offset of levels 2–4 (each pair to rounding, the
/// FMM-weighted measure over the conditioning), against `rcr` on a sample, columns Rotation's
/// M2L of e_j bit for bit, two vectors each the one-vector call, f32 the f64 matrices rounded,
/// the bytes and the selection with one local plane.
#[test]
fn yukawa_dense_matches_rotation() {
    let degrees = [0, 0, 7, 6, 5];
    let levels = [false, false, true, true, true];
    for kappa in KAPPA_ROOTS {
        let k = yukawa(kappa);
        let tables = HelmholtzTables::<f64>::build(k, &degrees, DEFAULT_BLOCK_BUDGET);
        let dense = DenseM2l::<f64>::build(k, &degrees, &levels);
        let dense32 = DenseM2l::<f32>::build(k, &degrees, &levels);
        assert_eq!(
            dense.bytes(),
            (2..5)
                .map(|l| DenseM2l::<f64>::level_bytes_for(degrees[l], 1))
                .sum::<usize>()
        );
        let mut scratch = tables.scratch();
        let mut worst = 0.0f64;
        for (level, &p) in degrees.iter().enumerate().skip(2) {
            let n = Layout::new(p).len();
            let w = weights(&Level::new(tables.k_hat(level), 1.0), false, p);
            let mut identity = vec![0.0; n * n];
            for j in 0..n {
                identity[j * n + j] = 1.0;
            }
            for d in 0..M2L_OFFSET_COUNT {
                let a = dense.matrix(level, d);
                assert_eq!(a.len(), n * n);
                let rounded: Vec<f32> = a.iter().map(|&v| v as f32).collect();
                assert!(
                    rounded == dense32.matrix(level, d),
                    "f32, level {level}, {d}"
                );
                if d % 29 == 0 {
                    let mut columns = vec![0.0; n * n];
                    tables.m2l(level, d, &identity, 1, &mut columns, n, &mut scratch);
                    assert_same(a, &columns, "columns");
                }
                let x = [input(p, 1, 7 * d as u64 + 1), input(p, 1, 7 * d as u64 + 2)].concat();
                let mut got = vec![0.5; 2 * n];
                dense.m2l(level, d, &x, 1, &mut got, 2);
                for v in 0..2 {
                    let mut one = vec![0.5; n];
                    dense.m2l(level, d, &x[v * n..][..n], 1, &mut one, 1);
                    assert_same(&got[v * n..][..n], &one, "Dense, vector of two");
                    let mut expected = vec![0.5; n];
                    tables.m2l(
                        level,
                        d,
                        &x[v * n..][..n],
                        1,
                        &mut expected,
                        1,
                        &mut scratch,
                    );
                    let error = compare(&widened(&one), &widened(&expected), 0.5, &w);
                    worst = worst.max(error / conditioning(&tables, level));
                }
            }
        }
        println!("κ̂₀ = {kappa}: Yukawa Dense against Rotation, worst {worst:.2e}");
        assert!(worst <= BOUND, "κ̂₀ = {kappa}: {worst:e}");
    }

    // The bytes and the selection: one local plane halves a level's tables.
    for p in [0, 5, 12] {
        assert_eq!(
            2 * DenseM2l::<f64>::level_bytes_for(p, 1),
            DenseM2l::<f64>::level_bytes(p)
        );
        assert_eq!(
            DenseM2l::<f64>::level_bytes_for(p, 2),
            DenseM2l::<f64>::level_bytes(p)
        );
    }
    // Two one-plane levels fit where one two-plane level does.
    let degrees = [0, 0, 4, 4, 4];
    let budget = DenseM2l::<f64>::level_bytes_for(4, 2);
    assert_eq!(
        DenseM2l::<f64>::select_for(&degrees, budget, 1),
        [false, false, false, true, true]
    );
    assert_eq!(
        DenseM2l::<f64>::select_for(&degrees, budget, 2),
        [false, false, false, false, true]
    );
    assert_eq!(
        DenseM2l::<f64>::select_for(&degrees, budget, 2),
        DenseM2l::<f64>::select(&degrees, budget)
    );
    assert_eq!(DenseM2l::<f64>::select_for(&degrees, 0, 1), [false; 5]);
}

/// Dense against `rcr::m2l_yukawa` on a sample of offsets.
#[test]
fn yukawa_dense_sample_matches_rcr() {
    let k = yukawa(8.0);
    let degrees = [0, 0, 0, 9];
    let dense = DenseM2l::<f64>::build(k, &degrees, &[false, false, false, true]);
    let tables = HelmholtzTables::<f64>::build(k, &degrees, DEFAULT_BLOCK_BUDGET);
    let (level, p) = (3, 9);
    let len = Layout::new(p).len();
    let k_level = Level::new(tables.k_hat(level), 1.0);
    let w = weights(&k_level, false, p);
    let from = Frame::new([0.0; 3], 1.0);
    let mut ws = Workspace::new(p);
    let mut worst = 0.0f64;
    for d in (0..M2L_OFFSET_COUNT).step_by(13) {
        let to = Frame::new(m2l_offsets()[d].map(|x| 2.0 * x as f64), 1.0);
        let table = rcr::m2l_table(p, p, &k_level, &from, &to);
        let x = input(p, 1, d as u64);
        let mut got = vec![0.0; len];
        dense.m2l(level, d, &x, 1, &mut got, 1);
        let mut expected = vec![0.0; len];
        rcr::m2l_yukawa(
            p,
            p,
            &k_level,
            &from,
            &to,
            &table,
            &mut ws,
            &x,
            &mut expected,
        );
        worst = worst.max(compare(&widened(&got), &widened(&expected), 0.0, &w));
    }
    let worst = worst / conditioning(&tables, level);
    println!("Yukawa Dense against rcr, worst {worst:.2e}");
    assert!(worst <= BOUND, "{worst:e}");
}

#[test]
#[should_panic(expected = "one multipole plane")]
fn yukawa_m2l_of_two_planes_is_refused() {
    let tables = HelmholtzTables::<f64>::build(yukawa(2.0), &[0, 0, 3], DEFAULT_BLOCK_BUDGET);
    let mut scratch = tables.scratch();
    let n = Layout::new(3).len();
    let mut local = vec![0.0; 2 * n];
    tables.m2l(2, 0, &vec![0.0; 2 * n], 2, &mut local, 1, &mut scratch);
}

#[test]
#[should_panic(expected = "one multipole plane")]
fn yukawa_dense_m2l_of_two_planes_is_refused() {
    let dense = DenseM2l::<f64>::build(yukawa(2.0), &[0, 0, 3], &[false, false, true]);
    let n = Layout::new(3).len();
    let mut local = vec![0.0; 2 * n];
    dense.m2l(2, 0, &vec![0.0; 2 * n], 2, &mut local, 1);
}

#[test]
#[should_panic(expected = "Yukawa")]
fn batched_refuses_yukawa() {
    let _ = BatchedM2l::<f64>::build(yukawa(2.0), &[0, 0, 3], &[false, false, true]);
}

#[test]
#[should_panic(expected = "Yukawa")]
fn classes_refuse_yukawa() {
    let _ = ClassM2l::<f64>::build(yukawa(2.0), &[0, 0, 3], &[false, false, true]);
}

#[test]
fn other_wavenumbers_are_still_refused() {
    for (re, im) in [
        (-1.0, 0.0),
        (0.0, -1.0),
        (1.0, -0.5),
        (0.0, 0.0),
        (f64::NAN, 1.0),
        (0.0, f64::INFINITY),
    ] {
        let k = Complex::new(re, im);
        let refused = std::panic::catch_unwind(|| {
            HelmholtzTables::<f64>::build(k, &[0, 0, 2], DEFAULT_BLOCK_BUDGET)
        })
        .is_err();
        assert!(refused, "k̂₀ = {re} + {im}i was accepted");
    }
}
