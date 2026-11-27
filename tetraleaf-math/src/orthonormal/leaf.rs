//! The real orthonormal harmonics of a batch of points at the leaf degrees (Helmholtz
//! CH4.9a2a; CONVENTIONS §3.15, "Angular basis").
//!
//! [`harmonics`](super::harmonics) forms each point's harmonics by the telescoped recurrence of
//! [`super::legendre`], with a binary exponent per column: robust at the degrees in the hundreds
//! that the high-frequency translations need, and two square roots and four divisions per term
//! (Helmholtz CH4.9a1, `tetraleaf-validate/results/helmholtz-ch4-9a1.md` §3). The leaf
//! operators need the harmonics of every point of a leaf at the leaf's degree, at most a few
//! dozen. [`LeafHarmonics`] forms the same values from coefficients tabulated once, for B points
//! at a time as lanes, so that each step is a few multiply-adds over B lanes.
//!
//! # The recurrence
//!
//! With Q̄ₙᵐ = √((2n + 1)/(4π)) fₘ P̄ₙᵐ(cos θ), f₀ = 1 and fₘ = √2 otherwise (so that
//! Yₙᵐ = Q̄ₙᵐ cos mφ and Yₙ⁻ᵐ = Q̄ₙᵐ sin mφ, the values and layout of
//! [`harmonics`](super::harmonics)), the plain three-term recurrence of the normalised Legendre
//! functions becomes
//!
//! ```text
//! Q̄₀⁰ = 1/√(4π),   Q̄ₘᵐ = dₘ sin θ Q̄ₘ₋₁ᵐ⁻¹,   Q̄ₘ₊₁ᵐ = √(2m + 3) cos θ Q̄ₘᵐ,
//! Q̄ₙᵐ = aₙₘ cos θ Q̄ₙ₋₁ᵐ − bₙₘ Q̄ₙ₋₂ᵐ   (n ≥ m + 2),
//!
//! dₘ = √((2m + 1)/(2m)) (√3 for m = 1),   aₙₘ = √((4n² − 1)/(n² − m²)),
//! bₙₘ = √((2n + 1)((n − 1)² − m²)/((2n − 3)(n² − m²))),
//! ```
//!
//! each coefficient formed in f64 and rounded to `T` once, in [`LeafHarmonics::new`]. The
//! recurrence's rounding grows like n² u near the poles (spike S1, which is why
//! [`super::legendre`] telescopes it); at the leaf degrees that is far below the leaf's
//! tolerances, and `tests/orthonormal/leaf.rs` holds it to a contract per term against
//! [`harmonics`](super::harmonics) (f64 to p = 100, f32 to p = 64). cos θ, sin θ, cos φ and
//! sin φ come from the coordinates, (cos mφ, sin mφ) by the recursion of (x + iy)/ρ, with no
//! trigonometric call (§3.5); at the centre θ = 0 and φ = 0, on the z-axis φ = 0, as for
//! [`harmonics`](super::harmonics).
//!
//! # Under-range values
//!
//! The diagonal Q̄ₘᵐ falls like sinᵐ θ. Where it falls below τ (2⁻⁹⁰⁰ in f64, 2⁻¹⁰⁰ in f32) the
//! lane's column is set to zero, and so are the later columns, whose diagonals are smaller: no
//! operand of the recurrence is then subnormal on that account (subnormal operands are slow on
//! x86; CH4.9a1 measured f32 at p = 41–45 up to 1.75× f64's time there). Near the pole
//! P̄ₙᵐ / P̄ₘᵐ → √C(n + m, 2m), at most about 2⁴¹ for n ≤ 64, so with the normalisation a
//! zeroed column held values below about τ 2⁴³: 2⁻⁸⁵⁷ in f64 and 2⁻⁵⁷ in f32, far below the unit
//! roundoff of either. (A tiny cos θ can still make a single product subnormal; points that
//! close to the equator are rare and cost time only.)
//!
//! # Lanes
//!
//! The B points of a call are lanes: every per-point quantity is a `[T; B]`, and every step a
//! loop over the B lanes without a branch on the data (the special cases are selects), which
//! the compiler vectorises without `unsafe`. A point's values do not depend on the other lanes,
//! so they are the same bits at every B and in every position.
//!
//! ```
//! use tetraleaf_math::Layout;
//! use tetraleaf_math::orthonormal::{LeafHarmonics, LeafScratch, harmonics};
//!
//! let p = 6;
//! let table = LeafHarmonics::<f64>::new(p);
//! let mut scratch = LeafScratch::<f64, 2>::new(p);
//! let points = [[0.3, -0.2, 0.5], [0.0, 0.0, -1.0]];
//! let mut lanes = vec![[0.0; 2]; Layout::new(p).len()];
//! table.values(p, &points, &mut scratch, &mut lanes);
//!
//! let mut one = vec![0.0; Layout::new(p).len()];
//! harmonics(p, points[0], &mut one);
//! for (lane, y) in lanes.iter().zip(&one) {
//!     assert!((lane[0] - y).abs() < 1e-14);
//! }
//! ```

use crate::{Layout, RealScalar};

/// The under-range threshold τ of the [module documentation](self): 2⁻⁹⁰⁰ for a type with
/// f64's range, 2⁻¹⁰⁰ for f32's.
#[inline]
pub(crate) fn under_range<T: RealScalar>() -> T {
    if T::max_value().to_f64() > 1e300 {
        T::from_f64(2f64.powi(-900))
    } else {
        T::from_f64(2f64.powi(-100))
    }
}

/// The coefficients of the leaf harmonics' recurrence up to a degree ([module
/// documentation](self)), formed in f64 and rounded to `T`. It allocates once, in
/// [`LeafHarmonics::new`]; the values of degree n do not depend on the table's degree.
#[derive(Clone, Debug)]
pub struct LeafHarmonics<T> {
    p: usize,
    /// Q̄₀⁰ = 1/√(4π).
    y00: T,
    /// dₘ at index m (index 0 unused).
    diagonal: Vec<T>,
    /// √(2m + 3) at index m.
    first: Vec<T>,
    /// (aₙₘ, bₙₘ) of column m at `offsets[m] + n − m − 2`, n = m + 2 … p.
    steps: Vec<(T, T)>,
    offsets: Vec<usize>,
}

impl<T: RealScalar> LeafHarmonics<T> {
    /// The coefficients for degrees up to `p`.
    pub fn new(p: usize) -> Self {
        let diagonal = (0..=p)
            .map(|m| {
                let d = match m {
                    0 => 1.0,
                    1 => 3f64.sqrt(),
                    _ => ((2 * m + 1) as f64 / (2 * m) as f64).sqrt(),
                };
                T::from_f64(d)
            })
            .collect();
        let first = (0..=p)
            .map(|m| T::from_f64(((2 * m + 3) as f64).sqrt()))
            .collect();
        let mut offsets = Vec::with_capacity(p + 1);
        let mut steps = Vec::new();
        for m in 0..=p {
            offsets.push(steps.len());
            for n in m + 2..=p {
                let (nf, mf) = (n as f64, m as f64);
                let nm = nf * nf - mf * mf;
                let a = ((4.0 * nf * nf - 1.0) / nm).sqrt();
                let b = ((2.0 * nf + 1.0) * ((nf - 1.0) * (nf - 1.0) - mf * mf)
                    / ((2.0 * nf - 3.0) * nm))
                    .sqrt();
                steps.push((T::from_f64(a), T::from_f64(b)));
            }
        }
        Self {
            p,
            y00: T::from_f64(1.0 / (4.0 * core::f64::consts::PI).sqrt()),
            diagonal,
            first,
            steps,
            offsets,
        }
    }

    /// The largest degree of the table.
    pub fn p(&self) -> usize {
        self.p
    }

    /// The harmonics of degree at most `p` of the B points `points` (offsets from the
    /// expansion's centre), column by column: for each order m = 0 … p, calls
    /// `column(m, q, cos, sin)` with q[n − m] = Q̄ₙᵐ per lane for n = m … p, and the lanes'
    /// (cos mφ, sin mφ), so that Yₙᵐ = Q̄ₙᵐ cos mφ and Yₙ⁻ᵐ = Q̄ₙᵐ sin mφ (for m = 0, cos 0φ = 1
    /// and sin 0φ = 0). The columns are those of the [module documentation](self), under-range
    /// lanes zeroed.
    ///
    /// # Panics
    ///
    /// If `p` exceeds the table's degree or `scratch`'s.
    pub fn columns<const B: usize>(
        &self,
        p: usize,
        points: &[[T; 3]; B],
        scratch: &mut LeafScratch<T, B>,
        mut column: impl FnMut(usize, &[[T; B]], &[T; B], &[T; B]),
    ) {
        assert!(p <= self.p, "degree {p} above the table's {}", self.p);
        assert!(
            p <= scratch.p,
            "degree {p} above the scratch's {}",
            scratch.p
        );
        let lanes = Lanes::of(points);
        let tau = under_range::<T>();
        let q = &mut scratch.column[..=p];
        let mut diag = [self.y00; B];
        let (mut cos_m, mut sin_m) = ([T::one(); B], [T::zero(); B]);
        for m in 0..=p {
            if m > 0 {
                let d = self.diagonal[m];
                for j in 0..B {
                    let next = d * lanes.sin_theta[j] * diag[j];
                    diag[j] = if next.abs() < tau { T::zero() } else { next };
                    let (c, s) = (cos_m[j], sin_m[j]);
                    cos_m[j] = c * lanes.cos_phi[j] - s * lanes.sin_phi[j];
                    sin_m[j] = s * lanes.cos_phi[j] + c * lanes.sin_phi[j];
                }
            }
            q[m] = diag;
            if m < p {
                let f = self.first[m];
                q[m + 1] = core::array::from_fn(|j| f * lanes.cos_theta[j] * diag[j]);
            }
            let steps = &self.steps[self.offsets[m]..];
            for (n, &(a, b)) in (m + 2..=p).zip(steps) {
                let (one, two) = (q[n - 1], q[n - 2]);
                q[n] = core::array::from_fn(|j| a * lanes.cos_theta[j] * one[j] - b * two[j]);
            }
            column(m, &q[m..=p], &cos_m, &sin_m);
        }
    }

    /// The harmonics of degree at most `p` of the B points `points` in the layout of §3.6, per
    /// lane: `out[idx(n, m)][j]` = Yₙᵐ of point j, the values of [`harmonics`](super::harmonics)
    /// to the contract of the [module documentation](self). From [`columns`](Self::columns).
    ///
    /// # Panics
    ///
    /// If `p` exceeds the table's or the scratch's degree, or `out.len() != (p + 1)²`.
    pub fn values<const B: usize>(
        &self,
        p: usize,
        points: &[[T; 3]; B],
        scratch: &mut LeafScratch<T, B>,
        out: &mut [[T; B]],
    ) {
        let layout = Layout::new(p);
        assert_eq!(
            out.len(),
            layout.len(),
            "`out` must have length (p + 1)² = {}",
            layout.len()
        );
        self.columns(p, points, scratch, |m, q, cos, sin| {
            for (n, values) in (m..=p).zip(q) {
                if m == 0 {
                    out[layout.idx(n, 0)] = *values;
                } else {
                    let mi = m as isize;
                    out[layout.idx(n, mi)] = core::array::from_fn(|j| values[j] * cos[j]);
                    out[layout.idx(n, -mi)] = core::array::from_fn(|j| values[j] * sin[j]);
                }
            }
        });
    }
}

/// The work space of [`LeafHarmonics`] for B lanes up to a degree: one column, p + 1 rows of B
/// values. It allocates once, in [`LeafScratch::new`].
#[derive(Clone, Debug)]
pub struct LeafScratch<T, const B: usize> {
    p: usize,
    column: Vec<[T; B]>,
}

impl<T: RealScalar, const B: usize> LeafScratch<T, B> {
    /// Work space for degrees up to `p`.
    pub fn new(p: usize) -> Self {
        Self {
            p,
            column: vec![[T::zero(); B]; p + 1],
        }
    }
}

/// cos θ, sin θ, cos φ and sin φ of each lane, from its coordinates, with
/// [`harmonics`](super::harmonics)' conventions at the centre and on the z-axis.
struct Lanes<T, const B: usize> {
    cos_theta: [T; B],
    sin_theta: [T; B],
    cos_phi: [T; B],
    sin_phi: [T; B],
}

impl<T: RealScalar, const B: usize> Lanes<T, B> {
    #[inline]
    fn of(points: &[[T; 3]; B]) -> Self {
        let mut lanes = Self {
            cos_theta: [T::one(); B],
            sin_theta: [T::zero(); B],
            cos_phi: [T::one(); B],
            sin_phi: [T::zero(); B],
        };
        for (j, &[x, y, z]) in points.iter().enumerate() {
            let rho2 = x * x + y * y;
            let rho = rho2.sqrt();
            let r = (rho2 + z * z).sqrt();
            // Selects: a zero radius keeps θ = 0, a zero ρ keeps φ = 0.
            let (r_safe, rho_safe) = (
                if r > T::zero() { r } else { T::one() },
                if rho > T::zero() { rho } else { T::one() },
            );
            lanes.cos_theta[j] = if r > T::zero() { z / r_safe } else { T::one() };
            lanes.sin_theta[j] = if r > T::zero() {
                rho / r_safe
            } else {
                T::zero()
            };
            lanes.cos_phi[j] = if rho > T::zero() {
                x / rho_safe
            } else {
                T::one()
            };
            lanes.sin_phi[j] = if rho > T::zero() {
                y / rho_safe
            } else {
                T::zero()
            };
        }
        lanes
    }
}
