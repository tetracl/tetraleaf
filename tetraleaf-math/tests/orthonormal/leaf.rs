//! The leaf harmonics (Helmholtz CH4.9a2a, docs/helmholtz/tasks/CH4.9a2a.md, "Tests that prove
//! it") against `harmonics` per point and term, over the leaf frame's directions.

use tetraleaf_math::orthonormal::{LeafHarmonics, LeafScratch, harmonics};
use tetraleaf_math::{Layout, RealScalar};

use crate::common::SplitMix64;

/// The contract per term: |Y − Y_ref| ≤ C (n + 1)² u √((2n + 1)/(4π)) + τ_abs, with u the type's
/// unit roundoff: the plain recurrence's n² u near the poles (module docs). The worst measured
/// on the M6 is 1.02 in f64 (p ≤ 100) and 0.87 in f32 (p ≤ 64), both at the lowest degrees (at
/// p = 100 the degree 100 terms reach about 0.6), printed by each test; C keeps a margin of
/// about ten (root `CLAUDE.md`: libm differs between machines).
const C: f64 = 10.0;

/// The directions of the plan: random points of [−1, 1)³, the poles, near-pole points, the
/// axes, the equator, the plane x = 0, the centre and the cube's corners.
fn directions() -> Vec<[f64; 3]> {
    let mut rng = SplitMix64::new(0x4c48_4152);
    let mut points: Vec<[f64; 3]> = (0..1000)
        .map(|_| {
            [
                rng.range(-1.0, 1.0),
                rng.range(-1.0, 1.0),
                rng.range(-1.0, 1.0),
            ]
        })
        .collect();
    points.extend([
        [0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0],
        [0.0, 0.0, 0.3],
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, -0.7, 0.0],
        [0.6, -0.8, 0.0],
        [0.0, 0.6, 0.8],
        [0.0, -0.6, -0.8],
    ]);
    for e in [1, 2, 3, 5, 8, 12, 20, 40, 100, 160, 300] {
        let t = 10f64.powi(-e);
        points.push([t, 0.0, 1.0]);
        points.push([-t, t, -0.5]);
        points.push([0.0, t, 0.9]);
    }
    for corner in 0..8 {
        let s = |bit: usize| if corner & bit == 0 { -1.0 } else { 1.0 };
        points.push([s(1), s(2), s(4)]);
    }
    points
}

/// The largest |Y − Y_ref| / ((n + 1)² u √((2n + 1)/(4π))) over the points and terms of degree
/// at most p, with the values of `T` taken from B lanes at a time, checked against the contract
/// (plus `floor` absolute); and that no value is subnormal.
fn worst<T: RealScalar, const B: usize>(p: usize, unit: f64, floor: f64) -> f64 {
    let layout = Layout::new(p);
    let table = LeafHarmonics::<T>::new(p);
    let mut scratch = LeafScratch::<T, B>::new(p);
    let mut lanes = vec![[T::zero(); B]; layout.len()];
    let mut reference = vec![0.0; layout.len()];
    let points = directions();
    let mut worst = 0.0f64;
    for batch in points.chunks(B) {
        let mut block = [[T::zero(); 3]; B];
        for (b, x) in block.iter_mut().zip(batch) {
            *b = x.map(T::from_f64);
        }
        table.values(p, &block, &mut scratch, &mut lanes);
        for (j, x) in batch.iter().enumerate() {
            // The reference at the point as `T` holds it, in f64.
            harmonics(p, block[j].map(|c| c.to_f64()), &mut reference);
            for (n, range) in layout.degrees() {
                let scale =
                    ((n + 1) * (n + 1)) as f64 * unit * ((2 * n + 1) as f64 / (4.0 * PI)).sqrt();
                for i in range {
                    let y = lanes[i][j].to_f64();
                    assert!(
                        y == 0.0 || y.abs() >= T::min_positive_value().to_f64(),
                        "subnormal Y at {x:?}, slot {i}, p = {p}: {y:e}"
                    );
                    let err = (y - reference[i]).abs();
                    assert!(
                        err <= C * scale + floor,
                        "p = {p}, point {x:?}, slot {i}: {y:e} against {:e}",
                        reference[i]
                    );
                    worst = worst.max(err / scale);
                }
            }
        }
    }
    worst
}

use std::f64::consts::PI;

#[test]
fn the_leaf_harmonics_within_the_contract_in_f64() {
    let u = f64::EPSILON / 2.0;
    for p in [0, 1, 2, 5, 8, 17, 32, 45, 64, 100] {
        let w = worst::<f64, 8>(p, u, 0.0);
        println!("f64 p = {p}: worst {w:.3} (n + 1)² u √((2n + 1)/(4π))");
    }
}

#[test]
fn the_leaf_harmonics_within_the_contract_in_f32() {
    let u = f64::from(f32::EPSILON) / 2.0;
    // f32's under-range threshold 2⁻¹⁰⁰ leaves at most about 2⁻⁵⁷ (module docs).
    let floor = 2f64.powi(-57);
    for p in [0, 1, 2, 5, 8, 17, 32, 45, 64] {
        let w = worst::<f32, 16>(p, u, floor);
        println!("f32 p = {p}: worst {w:.3} (n + 1)² u √((2n + 1)/(4π))");
    }
}

/// A point's values are the same bits at every lane count and position.
#[test]
fn the_lanes_do_not_depend_on_each_other() {
    let p = 23;
    let len = Layout::new(p).len();
    let points = directions();
    let table = LeafHarmonics::<f64>::new(p);
    let one: Vec<Vec<f64>> = points
        .iter()
        .map(|&x| {
            let mut s = LeafScratch::<f64, 1>::new(p);
            let mut out = vec![[0.0; 1]; len];
            table.values(p, &[x], &mut s, &mut out);
            out.iter().map(|v| v[0]).collect()
        })
        .collect();
    fn check<const B: usize>(
        p: usize,
        table: &LeafHarmonics<f64>,
        points: &[[f64; 3]],
        one: &[Vec<f64>],
    ) {
        let len = Layout::new(p).len();
        let mut s = LeafScratch::<f64, B>::new(p);
        let mut out = vec![[0.0; B]; len];
        // Shifted by one so that every point also sits in another lane.
        for shift in [0, 1] {
            for (start, batch) in points[shift..].chunks(B).enumerate() {
                let mut block = [[0.5, 0.5, 0.5]; B];
                block[..batch.len()].copy_from_slice(batch);
                table.values(p, &block, &mut s, &mut out);
                for j in 0..batch.len() {
                    let k = shift + start * B + j;
                    for (row, &want) in out.iter().zip(&one[k]) {
                        assert_eq!(row[j].to_bits(), want.to_bits(), "B = {B}, point {k}");
                    }
                }
            }
        }
    }
    check::<4>(p, &table, &points, &one);
    check::<8>(p, &table, &points, &one);
    check::<16>(p, &table, &points, &one);
}

/// `columns` hands out the values `values` stores; a table of a larger degree gives the same
/// bits at a smaller p; a degree above the table panics.
#[test]
fn columns_values_and_degrees() {
    let points = [
        [0.3, -0.4, 0.2],
        [0.0, 0.0, -0.9],
        [-0.7, 0.1, 0.0],
        [0.0; 3],
    ];
    let p = 12;
    let layout = Layout::new(p);
    let small = LeafHarmonics::<f64>::new(p);
    let large = LeafHarmonics::<f64>::new(40);
    let mut s = LeafScratch::<f64, 4>::new(40);
    let (mut a, mut b) = (vec![[0.0; 4]; layout.len()], vec![[0.0; 4]; layout.len()]);
    small.values(p, &points, &mut s, &mut a);
    large.values(p, &points, &mut s, &mut b);
    assert_eq!(a, b);
    let mut seen = 0;
    large.columns(p, &points, &mut s, |m, q, cos, sin| {
        assert_eq!(q.len(), p + 1 - m);
        for (n, v) in (m..=p).zip(q) {
            for j in 0..4 {
                let (c, sn) = if m == 0 { (1.0, 0.0) } else { (cos[j], sin[j]) };
                let mi = m as isize;
                assert_eq!(
                    a[layout.idx(n, mi)][j],
                    if m == 0 { v[j] } else { v[j] * c }
                );
                if m > 0 {
                    assert_eq!(a[layout.idx(n, -mi)][j], v[j] * sn);
                }
            }
            seen += 1;
        }
    });
    assert_eq!(seen, (p + 1) * (p + 2) / 2);
    let refused = std::panic::catch_unwind(|| {
        let mut s = LeafScratch::<f64, 4>::new(13);
        let mut out = vec![[0.0; 4]; Layout::new(13).len()];
        small.values(13, &points, &mut s, &mut out);
    });
    assert!(refused.is_err(), "a degree above the table's is refused");
}
