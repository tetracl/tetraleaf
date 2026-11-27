//! The real orthonormal harmonics of a batch of points at the leaf degrees, and their angular
//! derivatives (Helmholtz CH4.9a2a and CH4.9b1; CONVENTIONS §3.15, "Angular basis").
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
//! [`harmonics`](super::harmonics). ρ and r come from x² + y² and x² + y² + z², not `hypot`, so
//! the coordinates' squares must not overflow: |x|, |y|, |z| up to about 10¹⁵⁰ in f64 (10¹⁸ in f32;
//! a leaf's offsets lie within √3). Where ρ² (or r²) falls below the normal range, x and y (and z)
//! are first multiplied by 2⁶⁰⁰ in f64 (2¹⁰⁰ in f32), exactly, so that cos φ and sin φ (and cos θ)
//! keep their precision next to the z-axis (and the centre): the values there carry sinᵐ θ and
//! would not show it, but the derivatives of order 1 do (Helmholtz CH4.9b1). The scale is a select,
//! 1 elsewhere, so every other point's values are the same bits as without it; only x = y = 0 is
//! on the z-axis, and only the centre itself has θ = 0.
//!
//! # Under-range values
//!
//! The diagonal Q̄ₘᵐ falls like sinᵐ θ. Where it falls below τ (2⁻⁹⁰⁰ in f64, 2⁻¹⁰⁰ in f32) the
//! lane's column is set to zero, and so are the later columns, whose diagonals are smaller: no
//! operand of the recurrence is then subnormal on that account (subnormal operands are slow on
//! x86; CH4.9a1 measured f32 at p = 41–45 up to 1.75× f64's time there). Near the pole
//! P̄ₙᵐ / P̄ₘᵐ → √C(n + m, 2m), at most about 2⁴¹ for n ≤ 64, so with the normalisation a
//! zeroed column held values below about τ 2⁴³: 2⁻⁸⁵⁷ in f64 and 2⁻⁵⁷ in f32, far below the unit
//! roundoff of either. In f32 the ratio, and with it the bound, grows to about 2⁶⁸ at n = 100
//! (2⁻³² left), so f32 is held to its contract to p = 64 only (CONVENTIONS §3.15 reserves f32 for
//! low degrees); f64's margin covers every degree the leaves take. (A tiny cos θ can still make a
//! single product subnormal; points that close to the equator are rare and cost time only.)
//!
//! # Derivatives
//!
//! The gradient of L2P needs ∂θ Yₙᵐ and (1/sin θ) ∂φ Yₙᵐ, the outputs of
//! [`harmonics_derivatives`](super::harmonics_derivatives) (Helmholtz CH4.9b1).
//! [`LeafHarmonics::derivative_columns`] forms them per order from a second column,
//! Uₙᵐ = Q̄ₙᵐ / sin θ (m ≥ 1): the recurrence is linear, so U obeys it too, from
//! Uₘᵐ = dₘ Q̄ₘ₋₁ᵐ⁻¹ (the diagonal before its factor sin θ), finite at the poles, with no division
//! by sin θ. From sin θ ∂θ Pₙᵐ = n cos θ Pₙᵐ − (n + m) Pₙ₋₁ᵐ in the columns' normalisation,
//!
//! ```text
//! ∂θ Q̄ₙᵐ = n cos θ Uₙᵐ − cₙₘ Uₙ₋₁ᵐ,   cₙₘ = √((2n + 1)(n² − m²)/(2n − 1))   (m ≥ 1, cₘₘ = 0),
//! ∂θ Q̄ₙ⁰ = −√(n(n + 1)/2) Q̄ₙ¹,
//! ```
//!
//! the second the ladder of `harmonics_derivatives` with the √2 of the columns, which is why
//! order 1 is formed before order 0 is handed out. Then ∂θ Yₙ^{±m} = ∂θ Q̄ₙᵐ {cos mφ, sin mφ} and
//! (1/sin θ) ∂φ Yₙ^{±m} = m Uₙᵐ {−sin mφ, cos mφ}, in the frame φ = 0 on the z-axis (θ = 0 at the
//! centre). U's diagonal is flushed below τ as Q̄'s is, but on its own value: at a pole Q̄₁¹ = 0
//! while U₁¹ = d₁ Y₀⁰ carries ∂θ Y₁¹. U carries one power of sin θ less than Q̄, so a flushed U
//! column held values below about τ 2⁴³ as Q̄'s does, and ∂θ at most n times that, far below u in
//! both precisions. Both terms of ∂θ are of the size n |U|, so its rounding is about n times U's,
//! the recurrence's n² u: `tests/orthonormal/leaf.rs` holds ∂θ and (1/sin θ) ∂φ per term to
//! 10 (n + 1)³ u √((2n + 1)/(4π)) against `harmonics_derivatives` (0.87 and 0.59 of it measured on
//! the M6 in f64 to p = 100, 0.76 and 0.26 in f32 to p = 64), the values to the contract above.
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

/// The power of two by which [`Lanes::of`] scales coordinates whose ρ² or r² falls below the
/// normal range: 2⁶⁰⁰ for a type with f64's range, 2¹⁰⁰ for f32's, so that the squares of the
/// smallest subnormal coordinate are normal and those of the largest scaled one finite.
#[inline]
fn rescale<T: RealScalar>() -> T {
    if T::max_value().to_f64() > 1e300 {
        T::from_f64(2f64.powi(600))
    } else {
        T::from_f64(2f64.powi(100))
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
    /// n at index n, for ∂θ (Helmholtz CH4.9b1).
    degrees: Vec<T>,
    /// √(n(n + 1)/2) at index n: ∂θ Q̄ₙ⁰ = −√(n(n + 1)/2) Q̄ₙ¹.
    ladder: Vec<T>,
    /// cₙₘ = √((2n + 1)(n² − m²)/(2n − 1)) of column m ≥ 1 at `over_sin_offsets[m] + n − m − 1`,
    /// n = m + 1 … p: ∂θ Q̄ₙᵐ = n cos θ Uₙᵐ − cₙₘ Uₙ₋₁ᵐ.
    over_sin: Vec<T>,
    over_sin_offsets: Vec<usize>,
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
        let degrees = (0..=p).map(|n| T::from_f64(n as f64)).collect();
        let ladder = (0..=p)
            .map(|n| T::from_f64(((n * (n + 1)) as f64 / 2.0).sqrt()))
            .collect();
        let mut over_sin_offsets = Vec::with_capacity(p + 1);
        let mut over_sin = Vec::new();
        for m in 0..=p {
            over_sin_offsets.push(over_sin.len());
            if m == 0 {
                continue;
            }
            for n in m + 1..=p {
                let (nf, mf) = (n as f64, m as f64);
                let c = ((2.0 * nf + 1.0) * (nf * nf - mf * mf) / (2.0 * nf - 1.0)).sqrt();
                over_sin.push(T::from_f64(c));
            }
        }
        Self {
            p,
            y00: T::from_f64(1.0 / (4.0 * core::f64::consts::PI).sqrt()),
            diagonal,
            first,
            steps,
            offsets,
            degrees,
            ladder,
            over_sin,
            over_sin_offsets,
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
            self.recur(m, p, diag, &lanes.cos_theta, q);
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

impl<T: RealScalar> LeafHarmonics<T> {
    /// The column of order m from its diagonal `start` (rows m … p of `q`, indexed by n): the
    /// first step and the three-term steps of the [module documentation](self). The recurrence is
    /// linear, so it also forms Uₙᵐ = Q̄ₙᵐ / sin θ from Uₘᵐ ("Derivatives").
    #[inline(always)]
    fn recur<const B: usize>(
        &self,
        m: usize,
        p: usize,
        start: [T; B],
        cos_theta: &[T; B],
        q: &mut [[T; B]],
    ) {
        q[m] = start;
        if m < p {
            let f = self.first[m];
            q[m + 1] = core::array::from_fn(|j| f * cos_theta[j] * start[j]);
        }
        let steps = &self.steps[self.offsets[m]..];
        for (n, &(a, b)) in (m + 2..=p).zip(steps) {
            let (one, two) = (q[n - 1], q[n - 2]);
            q[n] = core::array::from_fn(|j| a * cos_theta[j] * one[j] - b * two[j]);
        }
    }

    /// The harmonics of degree at most `p` of the B points `points` and their angular derivatives,
    /// column by column (Helmholtz CH4.9b1; the [module documentation](self), "Derivatives"): for
    /// each order m = 0 … p, calls `column(m, q, d, u, cos, sin)` with, per lane and for
    /// n = m … p at index n − m, q = Q̄ₙᵐ ([`columns`](Self::columns)' values, bit for bit),
    /// d = ∂θ Q̄ₙᵐ and u = Uₙᵐ = Q̄ₙᵐ / sin θ (zero for m = 0, whose terms have no φ derivative),
    /// and the lanes' (cos mφ, sin mφ). Then
    ///
    /// ```text
    /// ∂θ Yₙᵐ = d cos mφ,   ∂θ Yₙ⁻ᵐ = d sin mφ,
    /// (1/sin θ) ∂φ Yₙᵐ = −m u sin mφ,   (1/sin θ) ∂φ Yₙ⁻ᵐ = m u cos mφ,
    /// ```
    ///
    /// `harmonics_derivatives`' values, in the frame φ = 0 on the z-axis and θ = 0 at the centre.
    ///
    /// # Panics
    ///
    /// If `p` exceeds the table's degree or `scratch`'s.
    pub fn derivative_columns<const B: usize>(
        &self,
        p: usize,
        points: &[[T; 3]; B],
        scratch: &mut LeafScratch<T, B>,
        mut column: impl FnMut(usize, &[[T; B]], &[[T; B]], &[[T; B]], &[T; B], &[T; B]),
    ) {
        assert!(p <= self.p, "degree {p} above the table's {}", self.p);
        assert!(
            p <= scratch.p,
            "degree {p} above the scratch's {}",
            scratch.p
        );
        let lanes = Lanes::of(points);
        let cos_theta = &lanes.cos_theta;
        let tau = under_range::<T>();
        let LeafScratch {
            column: q,
            next,
            over_sin: u,
            over_sin_next: u_next,
            d_theta: d,
            zeros,
            ..
        } = scratch;
        let (q, next, u, u_next, d) = (
            &mut q[..=p],
            &mut next[..=p],
            &mut u[..=p],
            &mut u_next[..=p],
            &mut d[..=p],
        );
        let zeros = &zeros[..=p];
        let mut diag = [self.y00; B];
        let (mut cos_m, mut sin_m) = ([T::one(); B], [T::zero(); B]);
        // Moves the diagonal and the phase to order m ≥ 1: returns Uₘᵐ = dₘ Q̄ₘ₋₁ᵐ⁻¹ (the old
        // diagonal), and sets the diagonal to Q̄ₘᵐ = dₘ sin θ Q̄ₘ₋₁ᵐ⁻¹, each flushed below τ, the
        // operations of `columns`.
        let advance = |m: usize, diag: &mut [T; B], cos_m: &mut [T; B], sin_m: &mut [T; B]| {
            let d = self.diagonal[m];
            let mut start = [T::zero(); B];
            for j in 0..B {
                let over = d * diag[j];
                start[j] = if over.abs() < tau { T::zero() } else { over };
                let next = d * lanes.sin_theta[j] * diag[j];
                diag[j] = if next.abs() < tau { T::zero() } else { next };
                let (c, s) = (cos_m[j], sin_m[j]);
                cos_m[j] = c * lanes.cos_phi[j] - s * lanes.sin_phi[j];
                sin_m[j] = s * lanes.cos_phi[j] + c * lanes.sin_phi[j];
            }
            start
        };
        // ∂θ Q̄ₙᵐ = n cos θ Uₙᵐ − cₙₘ Uₙ₋₁ᵐ for m ≥ 1 (cₘₘ = 0).
        let d_theta = |m: usize, u: &[[T; B]], d: &mut [[T; B]]| {
            let dm = self.degrees[m];
            d[m] = core::array::from_fn(|j| dm * cos_theta[j] * u[m][j]);
            let cs = &self.over_sin[self.over_sin_offsets[m]..];
            for (n, &c) in (m + 1..=p).zip(cs) {
                let dn = self.degrees[n];
                let (one, two) = (u[n], u[n - 1]);
                d[n] = core::array::from_fn(|j| dn * cos_theta[j] * one[j] - c * two[j]);
            }
        };
        // Order 0, with order 1 formed first: ∂θ Q̄ₙ⁰ = −√(n(n + 1)/2) Q̄ₙ¹.
        self.recur(0, p, diag, cos_theta, q);
        d[0] = [T::zero(); B];
        if p >= 1 {
            let (mut cos_one, mut sin_one) = (cos_m, sin_m);
            let start = advance(1, &mut diag, &mut cos_one, &mut sin_one);
            self.recur(1, p, diag, cos_theta, next);
            self.recur(1, p, start, cos_theta, u_next);
            for n in 1..=p {
                let l = self.ladder[n];
                d[n] = core::array::from_fn(|j| -(l * next[n][j]));
            }
            column(0, q, d, zeros, &cos_m, &sin_m);
            (cos_m, sin_m) = (cos_one, sin_one);
            d_theta(1, u_next, d);
            column(1, &next[1..], &d[1..], &u_next[1..], &cos_m, &sin_m);
        } else {
            column(0, q, d, zeros, &cos_m, &sin_m);
        }
        for m in 2..=p {
            let start = advance(m, &mut diag, &mut cos_m, &mut sin_m);
            self.recur(m, p, diag, cos_theta, q);
            self.recur(m, p, start, cos_theta, u);
            d_theta(m, u, d);
            column(m, &q[m..], &d[m..], &u[m..], &cos_m, &sin_m);
        }
    }

    /// The harmonics and their angular derivatives of degree at most `p` of the B points `points`
    /// in the layout of §3.6, per lane: `value`, `d_theta` and `d_phi` hold Yₙᵐ, ∂θ Yₙᵐ and
    /// (1/sin θ) ∂φ Yₙᵐ at `[idx(n, m)][j]`, the outputs of
    /// [`harmonics_derivatives`](super::harmonics_derivatives) to the contract of the [module
    /// documentation](self). From [`derivative_columns`](Self::derivative_columns).
    ///
    /// # Panics
    ///
    /// If `p` exceeds the table's or the scratch's degree, or an output's length is not (p + 1)².
    pub fn derivatives<const B: usize>(
        &self,
        p: usize,
        points: &[[T; 3]; B],
        scratch: &mut LeafScratch<T, B>,
        value: &mut [[T; B]],
        d_theta: &mut [[T; B]],
        d_phi: &mut [[T; B]],
    ) {
        let layout = Layout::new(p);
        for (out, name) in [
            (&*value, "value"),
            (&*d_theta, "d_theta"),
            (&*d_phi, "d_phi"),
        ] {
            assert_eq!(
                out.len(),
                layout.len(),
                "`{name}` must have length (p + 1)² = {}",
                layout.len()
            );
        }
        self.derivative_columns(p, points, scratch, |m, q, d, u, cos, sin| {
            let mi = m as isize;
            let mt = T::from_f64(m as f64);
            for (row, n) in (m..=p).enumerate() {
                let pos = layout.idx(n, mi);
                if m == 0 {
                    value[pos] = q[row];
                    d_theta[pos] = d[row];
                    d_phi[pos] = [T::zero(); B];
                    continue;
                }
                let neg = layout.idx(n, -mi);
                let w: [T; B] = core::array::from_fn(|j| mt * u[row][j]);
                value[pos] = core::array::from_fn(|j| q[row][j] * cos[j]);
                value[neg] = core::array::from_fn(|j| q[row][j] * sin[j]);
                d_theta[pos] = core::array::from_fn(|j| d[row][j] * cos[j]);
                d_theta[neg] = core::array::from_fn(|j| d[row][j] * sin[j]);
                d_phi[pos] = core::array::from_fn(|j| -(w[j] * sin[j]));
                d_phi[neg] = core::array::from_fn(|j| w[j] * cos[j]);
            }
        });
    }
}

/// The work space of [`LeafHarmonics`] for B lanes up to a degree: the columns of
/// [`columns`](LeafHarmonics::columns) and [`derivative_columns`](LeafHarmonics::derivative_columns),
/// six of p + 1 rows of B values. It allocates once, in [`LeafScratch::new`].
#[derive(Clone, Debug)]
pub struct LeafScratch<T, const B: usize> {
    p: usize,
    column: Vec<[T; B]>,
    /// Order 1's column, formed before order 0 is handed out (`derivative_columns`).
    next: Vec<[T; B]>,
    /// Uₙᵐ of the current order, and of order 1.
    over_sin: Vec<[T; B]>,
    over_sin_next: Vec<[T; B]>,
    /// ∂θ Q̄ₙᵐ of the current order.
    d_theta: Vec<[T; B]>,
    /// Order 0's u, zero.
    zeros: Vec<[T; B]>,
}

impl<T: RealScalar, const B: usize> LeafScratch<T, B> {
    /// Work space for degrees up to `p`.
    pub fn new(p: usize) -> Self {
        let rows = vec![[T::zero(); B]; p + 1];
        Self {
            p,
            column: rows.clone(),
            next: rows.clone(),
            over_sin: rows.clone(),
            over_sin_next: rows.clone(),
            d_theta: rows.clone(),
            zeros: rows,
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
        let (one, min) = (T::one(), T::min_positive_value());
        let big = rescale::<T>();
        for (j, &[x, y, z]) in points.iter().enumerate() {
            // Where ρ² or r² falls below the normal range, the coordinates are scaled by a power of
            // two first, so that the direction keeps its precision (module docs); elsewhere the
            // scale is 1 and every operation the same bits as without it.
            let rho2 = x * x + y * y;
            let r2 = rho2 + z * z;
            let s_phi = if rho2 < min { big } else { one };
            let (xs, ys) = (x * s_phi, y * s_phi);
            let rho_s = (xs * xs + ys * ys).sqrt();
            let s_theta = if r2 < min { big } else { one };
            let (xt, yt, zt) = (x * s_theta, y * s_theta, z * s_theta);
            let r = (xt * xt + yt * yt + zt * zt).sqrt();
            let rho = rho_s / s_phi * s_theta;
            // Selects: a zero radius keeps θ = 0, a zero ρ keeps φ = 0.
            let (r_safe, rho_safe) = (
                if r > T::zero() { r } else { T::one() },
                if rho_s > T::zero() { rho_s } else { T::one() },
            );
            lanes.cos_theta[j] = if r > T::zero() { zt / r_safe } else { T::one() };
            lanes.sin_theta[j] = if r > T::zero() {
                rho / r_safe
            } else {
                T::zero()
            };
            lanes.cos_phi[j] = if rho_s > T::zero() {
                xs / rho_safe
            } else {
                T::one()
            };
            lanes.sin_phi[j] = if rho_s > T::zero() {
                ys / rho_safe
            } else {
                T::zero()
            };
        }
        lanes
    }
}
