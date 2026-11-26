//! S → F and F → R around the φ step (Helmholtz CH4.6a1, `Conversions::{s2f_modes, f2r_modes}`):
//! the modes with the φ step written as a plain DFT against the reference
//! `tetraleaf_ref::helmholtz::diagonal::{s2f, f2r}` and against the product
//! (`Conversions::{s2f, f2r}`), the layout of a ring's run, and a box's bits in every batch.
//!
//! Error measure: as `diagonal.rs`, the largest absolute difference of a value over the largest
//! absolute value of the reference's result, within 256 u. The DFT here is a sum in f64 over the
//! grid's own cos and sin tables (`Grid::azimuth_tables`, at (m J) mod N_φ), the trusted path of
//! the FFT that tetraleaf-exec puts in its place. A box in a batch against the box alone: bit for
//! bit (`to_bits`).

use tetraleaf_math::Complex;
use tetraleaf_ref::helmholtz::Level;
use tetraleaf_ref::helmholtz::diagonal::{self, Grid, SignatureWorkspace};
use tetraleaf_tables::helmholtz::diagonal::{ConversionScratch, Conversions};

/// Deterministic values in [−1, 1).
fn values(count: usize, seed: u64) -> Vec<f64> {
    let mut state = seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1;
    (0..count)
        .map(|_| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state >> 11) as f64 / (1u64 << 52) as f64 - 1.0
        })
        .collect()
}

fn complex_values(count: usize, seed: u64) -> Vec<Complex<f64>> {
    let v = values(2 * count, seed);
    v.as_chunks::<2>()
        .0
        .iter()
        .map(|&[re, im]| Complex::new(re, im))
        .collect()
}

/// A level with power scale 1 (k a ≥ 1) or below it.
fn level(below_one: bool) -> Level<f64> {
    if below_one {
        Level::new(Complex::new(0.7, 0.0), 0.5)
    } else {
        Level::new(Complex::new(20.0, 0.0), 0.5)
    }
}

fn largest(values: impl Iterator<Item = f64>) -> f64 {
    values.fold(0.0, f64::max)
}

/// e^{±i m φ_J} from the grid's tables, exactly the reference's angles.
fn exponential(grid: &Grid<f64>, m: i64, j: usize, sign: f64) -> Complex<f64> {
    let (cos, sin) = grid.azimuth_tables();
    let n = cos.len() as i64;
    let at = ((m * j as i64).rem_euclid(n)) as usize;
    Complex::new(cos[at], sign * sin[at])
}

/// The inverse DFT of each ring's run of `ring_len` modes into the ring's N_φ samples.
fn synthesis(
    grid: &Grid<f64>,
    p: usize,
    modes: &[Complex<f64>],
    ring_len: usize,
) -> Vec<Complex<f64>> {
    let azimuths = grid.azimuths();
    let mut out = vec![Complex::new(0.0, 0.0); (2 * p + 1) * azimuths];
    for ring in 0..2 * p + 1 {
        let run = &modes[ring * ring_len..(ring + 1) * ring_len];
        for j in 0..azimuths {
            let mut sum = Complex::new(0.0, 0.0);
            for m in -(p as i64)..=p as i64 {
                let slot = if m >= 0 {
                    m as usize
                } else {
                    ring_len - (-m) as usize
                };
                sum += run[slot] * exponential(grid, m, j, 1.0);
            }
            out[ring * azimuths + j] = sum;
        }
    }
    out
}

/// The forward DFT of each ring's N_φ samples into a run of `ring_len` (only the slots of
/// |m| ≤ p; the others NaN, so that a read of them shows).
fn analysis(
    grid: &Grid<f64>,
    p: usize,
    samples: &[Complex<f64>],
    ring_len: usize,
) -> Vec<Complex<f64>> {
    let azimuths = grid.azimuths();
    let mut out = vec![Complex::new(f64::NAN, f64::NAN); (2 * p + 1) * ring_len];
    for ring in 0..2 * p + 1 {
        let g = &samples[ring * azimuths..(ring + 1) * azimuths];
        for m in -(p as i64)..=p as i64 {
            let mut sum = Complex::new(0.0, 0.0);
            for (j, &value) in g.iter().enumerate() {
                sum += value * exponential(grid, m, j, -1.0);
            }
            let slot = if m >= 0 {
                m as usize
            } else {
                ring_len - (-m) as usize
            };
            out[ring * ring_len + slot] = sum;
        }
    }
    out
}

/// S → F by modes and the DFT against the reference and the product, F → R from the forward DFT
/// against the reference and the product, box by box, at ring length N_φ and at a longer one.
fn check(p: usize, below_one: bool, boxes: usize) {
    let level = level(below_one);
    let grid = Grid::<f64>::new(p);
    let conversions = Conversions::<f64>::new(p, level.scale());
    let mut scratch = ConversionScratch::new(&conversions);
    let (len, samples, rings) = (
        conversions.expansion_len(),
        conversions.samples(),
        conversions.rings(),
    );
    let multipoles = values(boxes * len, p as u64 + 3);
    let mut product = vec![Complex::new(0.0, 0.0); boxes * samples];
    conversions.s2f(&multipoles, &mut product, &mut scratch);
    let signatures = complex_values(boxes * samples, p as u64 + 5);
    let start = values(boxes * len, p as u64 + 7);
    let mut product_locals = start.clone();
    conversions.f2r(&signatures, &mut product_locals, &mut scratch);

    let mut ws = SignatureWorkspace::new(p);
    for ring_len in [
        conversions.azimuths(),
        2 * p + 1,
        conversions.azimuths() + 5,
    ] {
        let stride = rings * ring_len + 3;
        let mut modes = vec![Complex::new(f64::NAN, f64::NAN); boxes * stride];
        conversions.s2f_modes(&multipoles, &mut modes, ring_len, stride, &mut scratch);
        let mut locals = start.clone();
        let mut forward = vec![Complex::new(f64::NAN, f64::NAN); boxes * stride];
        for b in 0..boxes {
            // The padding untouched, every slot of a run written.
            assert!(
                modes[b * stride + rings * ring_len..(b + 1) * stride]
                    .iter()
                    .all(|z| z.re.is_nan())
            );
            assert!(
                modes[b * stride..b * stride + rings * ring_len]
                    .iter()
                    .all(|z| z.re.is_finite() && z.im.is_finite())
            );
            for ring in 0..rings {
                let run = &modes[b * stride + ring * ring_len..][..ring_len];
                assert!(
                    run[p + 1..ring_len - p]
                        .iter()
                        .all(|z| z.re == 0.0 && z.im == 0.0)
                );
            }
            let signature = synthesis(&grid, p, &modes[b * stride..], ring_len);
            let mut expected = vec![Complex::new(0.0, 0.0); samples];
            diagonal::s2f(
                &grid,
                &level,
                &multipoles[b * len..(b + 1) * len],
                &mut ws,
                &mut expected,
            );
            let scale = largest(expected.iter().map(|z| z.norm()));
            for (against, name) in [
                (&expected[..], "reference"),
                (&product[b * samples..(b + 1) * samples], "product"),
            ] {
                let error = largest(signature.iter().zip(against).map(|(a, e)| (a - e).norm()));
                assert!(
                    error <= 256.0 * f64::EPSILON * scale,
                    "S → F p {p} ring length {ring_len} box {b} against the {name}: {:.2} u",
                    error / (f64::EPSILON * scale)
                );
            }
            let x = analysis(
                &grid,
                p,
                &signatures[b * samples..(b + 1) * samples],
                ring_len,
            );
            forward[b * stride..b * stride + rings * ring_len].copy_from_slice(&x);
        }
        conversions.f2r_modes(&forward, ring_len, stride, &mut locals, &mut scratch);
        for b in 0..boxes {
            let mut expected = start[b * len..(b + 1) * len].to_vec();
            diagonal::f2r(
                &grid,
                &level,
                &signatures[b * samples..(b + 1) * samples],
                &mut ws,
                &mut expected,
            );
            let scale = largest(expected.iter().map(|v| v.abs()));
            let got = &locals[b * len..(b + 1) * len];
            for (against, name) in [
                (&expected[..], "reference"),
                (&product_locals[b * len..(b + 1) * len], "product"),
            ] {
                let error = largest(got.iter().zip(against).map(|(a, e)| (a - e).abs()));
                assert!(
                    error <= 256.0 * f64::EPSILON * scale,
                    "F → R p {p} ring length {ring_len} box {b} against the {name}: {:.2} u",
                    error / (f64::EPSILON * scale)
                );
            }
        }
    }
}

#[test]
fn modes_against_the_reference_and_the_product() {
    for p in [0, 1, 2, 3, 6, 12] {
        check(p, false, 5);
    }
    check(5, true, 3);
    check(17, false, 2);
}

#[test]
#[ignore = "large degrees: run in release"]
fn modes_at_large_degrees() {
    for p in [24, 41, 72] {
        check(p, false, 2);
    }
}

/// A box's modes and local in a batch bit for bit the box alone, at every position, and the
/// Legendre step's bits not depending on the ring length.
#[test]
fn a_box_has_the_same_bits_in_every_batch() {
    let p = 9;
    let conversions = Conversions::<f64>::new(p, 1.0);
    let mut scratch = ConversionScratch::new(&conversions);
    let (len, rings, azimuths) = (
        conversions.expansion_len(),
        conversions.rings(),
        conversions.azimuths(),
    );
    let run = rings * azimuths;
    for boxes in [1, 2, 5, 7, 9] {
        let multipoles = values(boxes * len, 11);
        let mut modes = vec![Complex::new(0.0, 0.0); boxes * run];
        conversions.s2f_modes(&multipoles, &mut modes, azimuths, run, &mut scratch);
        let x = complex_values(boxes * run, 13);
        let start = values(boxes * len, 17);
        let mut locals = start.clone();
        conversions.f2r_modes(&x, azimuths, run, &mut locals, &mut scratch);
        for b in 0..boxes {
            let mut alone = vec![Complex::new(0.0, 0.0); run];
            conversions.s2f_modes(
                &multipoles[b * len..(b + 1) * len],
                &mut alone,
                azimuths,
                run,
                &mut scratch,
            );
            assert!(
                alone
                    .iter()
                    .zip(&modes[b * run..(b + 1) * run])
                    .all(|(a, z)| {
                        a.re.to_bits() == z.re.to_bits() && a.im.to_bits() == z.im.to_bits()
                    }),
                "S → F: box {b} of {boxes}"
            );
            // At the compact ring length the same modes.
            let compact = 2 * p + 1;
            let mut short = vec![Complex::new(0.0, 0.0); rings * compact];
            conversions.s2f_modes(
                &multipoles[b * len..(b + 1) * len],
                &mut short,
                compact,
                rings * compact,
                &mut scratch,
            );
            for ring in 0..rings {
                for m in 0..=p {
                    let (a, z) = (short[ring * compact + m], alone[ring * azimuths + m]);
                    assert_eq!(
                        (a.re.to_bits(), a.im.to_bits()),
                        (z.re.to_bits(), z.im.to_bits())
                    );
                    if m > 0 {
                        let (a, z) = (
                            short[ring * compact + compact - m],
                            alone[ring * azimuths + azimuths - m],
                        );
                        assert_eq!(
                            (a.re.to_bits(), a.im.to_bits()),
                            (z.re.to_bits(), z.im.to_bits())
                        );
                    }
                }
            }
            let mut local = start[b * len..(b + 1) * len].to_vec();
            conversions.f2r_modes(
                &x[b * run..(b + 1) * run],
                azimuths,
                run,
                &mut local,
                &mut scratch,
            );
            assert!(
                local
                    .iter()
                    .zip(&locals[b * len..(b + 1) * len])
                    .all(|(a, z)| a.to_bits() == z.to_bits()),
                "F → R: box {b} of {boxes}"
            );
        }
    }
}

/// f32 against the f64 product, within 256 u of f32.
#[test]
fn f32_against_the_product() {
    let p = 8;
    let level = level(false);
    let c64 = Conversions::<f64>::new(p, level.scale());
    let c32 = Conversions::<f32>::new(p, level.scale());
    let (mut s64, mut s32) = (ConversionScratch::new(&c64), ConversionScratch::new(&c32));
    let grid = Grid::<f64>::new(p);
    let len = c64.expansion_len();
    let multipoles = values(2 * len, 23);
    let m32: Vec<f32> = multipoles.iter().map(|&v| v as f32).collect();
    let azimuths = c32.azimuths();
    let run = c32.rings() * azimuths;
    let mut modes = vec![Complex::new(0.0f32, 0.0); 2 * run];
    c32.s2f_modes(&m32, &mut modes, azimuths, run, &mut s32);
    let mut product = vec![Complex::new(0.0, 0.0); 2 * c64.samples()];
    c64.s2f(&multipoles, &mut product, &mut s64);
    for b in 0..2 {
        let wide: Vec<Complex<f64>> = modes[b * run..(b + 1) * run]
            .iter()
            .map(|z| Complex::new(z.re as f64, z.im as f64))
            .collect();
        let signature = synthesis(&grid, p, &wide, azimuths);
        let expected = &product[b * c64.samples()..(b + 1) * c64.samples()];
        let scale = largest(expected.iter().map(|z| z.norm()));
        let error = largest(signature.iter().zip(expected).map(|(a, e)| (a - e).norm()));
        assert!(error <= 256.0 * f32::EPSILON as f64 * scale, "f32 box {b}");
    }
}

#[test]
#[should_panic(expected = "cannot hold")]
fn refuses_a_short_ring() {
    let conversions = Conversions::<f64>::new(4, 1.0);
    let mut scratch = ConversionScratch::new(&conversions);
    let multipoles = vec![0.0; conversions.expansion_len()];
    let mut modes = vec![Complex::new(0.0, 0.0); 100];
    conversions.s2f_modes(&multipoles, &mut modes, 8, 81, &mut scratch);
}

#[test]
#[should_panic(expected = "need")]
fn refuses_too_few_modes() {
    let conversions = Conversions::<f64>::new(4, 1.0);
    let mut scratch = ConversionScratch::new(&conversions);
    let multipoles = vec![0.0; 2 * conversions.expansion_len()];
    let run = conversions.rings() * 9;
    let mut modes = vec![Complex::new(0.0, 0.0); 2 * run - 1];
    conversions.s2f_modes(&multipoles, &mut modes, 9, run, &mut scratch);
}
