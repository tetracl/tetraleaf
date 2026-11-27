//! `bessel_j_real` against `bessel_j` at real arguments (Helmholtz CH4.9a2a,
//! docs/helmholtz/tasks/CH4.9a2a.md, "Tests that prove it"): bit for bit, in f32 but for the
//! values it flushes beyond the turning point.

use tetraleaf_math::radial::{bessel_j, bessel_j_real};
use tetraleaf_math::{Complex, RealScalar};

/// The arguments of the plan, mixed within each batch.
const XS: [f64; 13] = [
    0.0, 1e-12, 1e-6, 0.3, 0.499_999, 0.5, 1.0, 2.5, 4.0, 17.0, 33.3, 60.0, 300.0,
];

/// Every argument of `XS` in B lanes (rotated so that each sits in several positions), at the
/// scales and lengths of the plan: each lane's values equal `bessel_j`'s real parts bit for bit,
/// or, for τ > 0, are zero where the reference is below τ at an order n ≥ x, and then no value is
/// subnormal. τ = 0 asks for every bit. Values the
/// reference does not hold finite (f32 at s < 1 and x of hundreds, an overflow of ĵₙ/sⁿ; a level
/// with s < 1 has |k| a < 1, so x < 1) are skipped.
fn check<T: Real, const B: usize>(tau: f64) -> usize {
    let mut flushed = 0;
    for s in [1.0, 0.3, 1e-3] {
        for len in [1, 2, 9, 46, 101] {
            for rotation in 0..XS.len() {
                let x: [T; B] =
                    core::array::from_fn(|j| T::from_f64(XS[(rotation + j) % XS.len()]));
                let s = T::from_f64(s);
                let mut lanes = vec![[T::zero(); B]; len];
                bessel_j_real(&x, s, &mut lanes);
                for j in 0..B {
                    let mut one = vec![Complex::new(T::zero(), T::zero()); len];
                    bessel_j(Complex::new(x[j], T::zero()), s, &mut one, None);
                    for n in 0..len {
                        let (got, want) = (lanes[n][j], one[n].re);
                        if !want.is_finite() {
                            // Outside T's range (f32 at s < 1 and large x, where no level has
                            // s < 1): nothing to compare.
                            continue;
                        }
                        let small = want.abs().to_f64() < tau && n as f64 >= x[j].to_f64();
                        if got.to_bits_u64() != want.to_bits_u64() {
                            assert!(
                                small && got == T::zero(),
                                "x = {}, s = {}, len {len}, n = {n}: {:e} against {:e}",
                                x[j].to_f64(),
                                s.to_f64(),
                                got.to_f64(),
                                want.to_f64()
                            );
                            flushed += 1;
                        }
                        let g = got.to_f64().abs();
                        // Where the routine flushes (τ > 0), nothing is subnormal; f64 keeps the
                        // reference's bits, subnormal ones included.
                        assert!(
                            tau == 0.0 || g == 0.0 || g >= T::min_positive_value().to_f64(),
                            "subnormal value at x = {}, n = {n}",
                            x[j].to_f64()
                        );
                    }
                }
            }
        }
    }
    flushed
}

/// The bits of a value of either precision, for the comparison.
trait Bits {
    fn to_bits_u64(self) -> u64;
}

impl Bits for f64 {
    fn to_bits_u64(self) -> u64 {
        self.to_bits()
    }
}

impl Bits for f32 {
    fn to_bits_u64(self) -> u64 {
        u64::from(self.to_bits())
    }
}

trait Real: RealScalar + Bits {}
impl<T: RealScalar + Bits> Real for T {}

#[test]
fn the_real_radial_functions_are_bessel_j_bit_for_bit_in_f64() {
    // f64 flushes nothing: every value is the reference's bits.
    let flushed = check::<f64, 1>(0.0) + check::<f64, 4>(0.0) + check::<f64, 8>(0.0);
    assert_eq!(flushed, 0, "f64 keeps every bit");
}

#[test]
fn the_real_radial_functions_are_bessel_j_bit_for_bit_in_f32_but_the_flushed() {
    let tau = 2f64.powi(-100);
    let flushed = check::<f32, 1>(tau) + check::<f32, 8>(tau) + check::<f32, 16>(tau);
    println!("f32: {flushed} values flushed below 2^-100 beyond the turning point");
    assert!(flushed > 0, "the f32 cases reach the flush");
}
