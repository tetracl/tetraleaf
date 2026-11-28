//! P2M and L2P of a whole leaf by batched kernels (Helmholtz CH4.9a2b and CH4.9b2; CONVENTIONS
//! §3.15, "Radial functions and scaling" and "Storage").
//!
//! [`tetraleaf_ref::helmholtz::leaf::p2m`] forms a box's stored multipole point by point:
//! per point the radial functions ĵₙ (`bessel_j`), the harmonics (`harmonics`), the factor
//! e^{Im(k)ρ − β}, and then the sum over the (p + 1)² terms,
//!
//! ```text
//! M̂ₙᵐ += Σⱼ qⱼ ĵₙ(kρⱼ) e^{Im(k)ρⱼ − β} Yₙᵐ(ρ̂ⱼ)
//! ```
//!
//! (its module documentation). [`LeafP2m`] forms the same sum for all the points of a leaf,
//! [`LANES`] at a time as lanes, on the leaf routines of `tetraleaf-math` (Helmholtz
//! CH4.9a2a): the harmonics of the batch column by column from tabulated coefficients
//! ([`LeafHarmonics::columns`]), and the radial functions in real arithmetic for real k
//! ([`bessel_j_real`]; complex k keeps [`bessel_j`], lane by lane). Per batch:
//!
//! 1. ρⱼ = √(x² + y² + z²) per lane (the reference takes `hypot`; the two agree to rounding);
//! 2. ĵₙ(kρⱼ) for n ≤ p, and eⱼ = e^{Im(k)ρⱼ − β} (for real k one value for the leaf);
//! 3. per charge vector and degree the weights wₙⱼ = (qⱼ eⱼ) ĵₙ(kρⱼ), complex in general,
//!    real for real k with real charges;
//! 4. per order m the column Q̄ₙᵐ of the batch (n = m … p) and (cos mφ, sin mφ), and per term
//!    the harmonics yⱼ = Q̄ₙᵐ cos mφ (and Q̄ₙᵐ sin mφ for −m) and the sums Σⱼ wₙⱼ yⱼ of each
//!    plane, the lanes added pairwise in a fixed order (eight pairs, then four, two, one), each
//!    sum added (+=) into the coefficient.
//!
//! The products in steps 3 and 4 are those of the reference in its order of operations
//! (q e, then times ĵₙ, then times Yₙᵐ = Q̄ₙᵐ cos mφ); only ρ, the routines' values and the
//! order of the sum over points differ. The last batch of a leaf is padded with the centre
//! and zero charges, whose terms are zeros. Every step is a loop over the lanes of fixed-size
//! arrays without a branch on the data, which the compiler vectorises without `unsafe`; the
//! start orders of the radial recursion, [`bessel_j`] of complex k and the factor's `exp`
//! run per lane.
//!
//! # Values and bits
//!
//! A leaf's multipole depends on its points and charges in their order alone, so it is the
//! same at every thread count and on any number of ranks that hold the leaf's points in the
//! same order. Each charge vector gets the arithmetic of the one-vector call (the vectors
//! share the routines' values, never a sum), so vector v of m is the one-vector kernel's bits.
//! For real k the real plane of complex charges takes the operations of real charges, so
//! complex charges with Im q = 0 give real charges' real plane bit for bit.
//!
//! Against the reference the values agree to a contract per term (`tests/helmholtz/leaf.rs`):
//! with N points and xⱼ = kρⱼ,
//!
//! ```text
//! |M̂ₙᵐ − M̂ₙᵐ_ref| ≤ C u √((2n + 1)/(4π)) Σⱼ |qⱼ eⱼ| (((n + 1)² + N) |ĵₙ(xⱼ)| + |xⱼ| |ĵₙ′(xⱼ)|)
//! ```
//!
//! the (n + 1)² from the harmonics' recurrence (CH4.9a2a), N from the two sums over points, and
//! the last term ĵₙ's condition at xⱼ: ρ by `sqrt` and the reference's by `hypot` differ in the
//! last bit, which near a zero of ĵₙ moves it by far more than u |ĵₙ|. The tests hold C = 10;
//! the measured worst is about 1 in f64 and f32.
//! In f32 the routines flush under-range values (the harmonics' columns below 2⁻¹⁰⁰, ĵₙ below
//! 2⁻¹⁰⁰ beyond the turning point; [`tetraleaf_math::orthonormal::leaf`], "Under-range
//! values"), which moves a term by at most about 2⁻⁵⁷ |wₙⱼ|, far below u.
//!
//! # Yukawa
//!
//! On Yukawa's real path (k = iκ, Helmholtz CH2.1c) the kernel restates
//! [`tetraleaf_ref::helmholtz::leaf::p2m_yukawa`], the real multipole of real charges,
//!
//! ```text
//! M̂ᵞₙᵐ += Σⱼ qⱼ îₙ(κρⱼ) e^{κρⱼ − β} Yₙᵐ(ρ̂ⱼ)
//! ```
//!
//! with the same steps: îₙ by [`modified_i`] per lane (the reference's routine at the reference's
//! length p + 1, so a lane's values are the reference's at its ρ), eⱼ = e^{κρⱼ − β} per lane, the
//! real weights wₙⱼ = (qⱼ eⱼ) îₙ(κρⱼ) as for real k, then the columns and the lane sums. Each real
//! charge column is one real multipole: [`LeafP2m::planes`] is 1 for real charges, and for complex
//! charges plane 0 is Re q's multipole and plane 1 Im q's (`tetraleaf-exec`'s densities: a complex
//! vector on Yukawa is two real vectors), each the real-charge kernel's bits on that column. Values
//! of îₙ and of eⱼ below the under-range threshold τ (2⁻¹⁰⁰ in f32, 2⁻⁹⁰⁰ in f64) are set to zero,
//! as [`bessel_j_real`] flushes ĵₙ: in f32 îₙ underflows at the leaf degrees near the centre and eⱼ
//! under strong damping (κ(a − ρ) ≳ 70), and no operand of the sum is then subnormal. A flushed
//! term is below τ |qⱼ| √((2n + 1)/(4π)) (îₙ ≤ 1).
//!
//! The contract is the one above with the factor's condition added: eⱼ now depends on ρⱼ, whose
//! `sqrt` and `hypot` differ in the last bit, so
//!
//! ```text
//! |M̂ᵞₙᵐ − M̂ᵞₙᵐ_ref| ≤ C u √((2n + 1)/(4π)) Σⱼ |qⱼ eⱼ| (((n + 1)² + N) îₙ(xⱼ) + xⱼ (|îₙ′(xⱼ)| + îₙ(xⱼ)))
//! ```
//!
//! with xⱼ = κρⱼ (`tests/helmholtz/leaf_yukawa.rs`: C = 10, about 0.8 measured in f64 to p = 45 and
//! κ̂ a to 104, 1.1 in f32 against the f64 reference with the f32 level's own rounding of s and β).
//!
//! # L2P
//!
//! [`tetraleaf_ref::helmholtz::leaf::l2p`] evaluates a box's stored local point by point:
//! φ̃ += Σ L̂ₙᵐ ĵₙ(k|v|) e^{Im(k)|v| − β} Yₙᵐ(v̂), and with gradients ∇φ̃ (its module documentation:
//! per degree the sums Σₘ L̂ₙᵐ Yₙᵐ, Σₘ L̂ₙᵐ ∂θ Yₙᵐ and Σₘ L̂ₙᵐ (1/sin θ) ∂φ Yₙᵐ, then the radial
//! functions, their derivatives and the tangential factor tₙ = ĵₙ/|v|, whose limit at the centre is
//! k/(3s) for n = 1 and 0 otherwise). [`LeafL2p`] evaluates all the targets of a leaf, [`LANES`] at
//! a time as lanes (Helmholtz CH4.9b2), on the same routines: the harmonics by
//! [`LeafHarmonics::columns`], or with gradients by [`LeafHarmonics::derivative_columns`] (Q̄ₙᵐ,
//! ∂θ Q̄ₙᵐ and Uₙᵐ = Q̄ₙᵐ / sin θ, CH4.9b1), the radial functions by [`bessel_j_real`] (with
//! derivatives [`bessel_j_real_derivatives`]) for real k and [`bessel_j`] per lane for complex k. Per
//! batch:
//!
//! 1. r = √(x² + y² + z²) per lane, ĵₙ(kr) (and ĵₙ′), e = e^{Im(k)r − β};
//! 2. per order m, degree n, charge vector and plane, from the coefficients c⁺ = L̂ₙᵐ and
//!    c⁻ = L̂ₙ⁻ᵐ per lane A = c⁺ cos mφ + c⁻ sin mφ and B = c⁻ cos mφ − c⁺ sin mφ, and the per-degree
//!    sums valueₙ += Q̄ₙᵐ A, θₙ += ∂θ Q̄ₙᵐ A, φₙ += m Uₙᵐ B (for m = 0: Q̄ₙ⁰ c and ∂θ Q̄ₙ⁰ c), which
//!    are the reference's three sums of degree n with the terms of ±m taken together;
//! 3. per lane the reference's degree loop in its order: φ̃ = Σₙ ĵₙ valueₙ, the radial part
//!    Σₙ ĵₙ′ valueₙ, the polar and azimuthal parts Σₙ tₙ θₙ and Σₙ tₙ φₙ (tₙ = ĵₙ (1/r), or the
//!    centre's limit by a select), the gradient k (radial) r̂ + (polar) θ̂ + (azimuthal) φ̂ in
//!    Cartesian components from the harmonics' own directions ([`LeafDirections`]), all times e,
//!    added (+=) into the target's outputs.
//!
//! The sum of step 2 is the kernel's inner loop: with gradients, per pair of terms ±m and plane
//! seven products and five sums per lane (A, B and the three sums), six operations per term, over
//! the lanes of fixed-size arrays with no reduction across lanes, which the compiler vectorises
//! without `unsafe`. The last batch is padded with the centre,
//! whose outputs are not written. No lane reads another: a target's outputs depend on its offset and
//! the locals alone, so they are the same bits in any leaf, position and batch, at every thread
//! count and on any number of ranks. Each charge vector gets the one-vector call's arithmetic (the
//! vectors share the routines' values), so vector v of m is the one-vector kernel's bits; the
//! potentials with gradients are the potentials without them bit for bit (the same columns, whose
//! values are [`LeafHarmonics::columns`]' bits).
//!
//! Against the reference the values agree to a contract per target (`tests/helmholtz/l2p.rs`): with
//! wₙ = √((2n + 1)/(4π)), Sₙ = Σₘ |L̂ₙᵐ|, x = kr and tₙ = |ĵₙ(x)|/r (its limit at r = 0),
//!
//! ```text
//! |φ̃ − φ̃_ref|  ≤ C u e Σₙ wₙ Sₙ (((n + 1)² + 2n + 1) |ĵₙ| + |x| |ĵₙ′|)
//! |gₐ − gₐ_ref| ≤ C u e Σₙ wₙ Sₙ (|k| (((n + 1)² + 2n + 3) |ĵₙ′| + |x| |ĵₙ|) + ((n + 1)³ + n(n + 1)) tₙ)
//! ```
//!
//! per component a of the gradient: the harmonics' (n + 1)² u and their derivatives' (n + 1)³ u
//! (CH4.9a2a, CH4.9b1), 2n + 1 terms per degree summed in another order, and r by `sqrt` against the
//! reference's `hypot` through ĵₙ′ and ĵₙ″ (bounded by the radial equation). The tests hold C = 10.
//!
//! # Example
//!
//! ```
//! use tetraleaf_math::{Complex, Layout};
//! use tetraleaf_ref::helmholtz::{Level, Workspace, leaf};
//! use tetraleaf_ref::Frame;
//! use tetraleaf_tables::helmholtz::leaf::{LeafCharges, LeafP2m, P2mScratch};
//!
//! let p = 10;
//! let level = Level::new(Complex::new(3.0, 0.0), 1.0);
//! let points = [[0.3, -0.2, 0.5], [-0.7, 0.1, 0.9], [0.0, 0.0, 0.0]];
//! // Complex charges, point-major: Re q, Im q per point.
//! let values = [1.0, 0.5, -0.25, 0.0, 2.0, -1.0];
//! let charges = LeafCharges { values: &values, vectors: 1, complex: true };
//! let kernel = LeafP2m::<f64>::new(p);
//! let mut scratch = P2mScratch::new(p, 1);
//! let len = Layout::new(p).len();
//! let mut multipole = vec![0.0; 2 * len];
//! kernel.p2m(p, &level, &points, charges, &mut scratch, &mut multipole);
//!
//! let q: Vec<Complex<f64>> = values.chunks(2).map(|c| Complex::new(c[0], c[1])).collect();
//! let mut reference = vec![0.0; 2 * len];
//! let frame = Frame::new([0.0; 3], 1.0);
//! leaf::p2m(p, &level, &frame, &points, &q, &mut Workspace::new(p), &mut reference);
//! for (a, b) in multipole.iter().zip(&reference) {
//!     assert!((a - b).abs() < 1e-13);
//! }
//! ```

use core::array;

use tetraleaf_math::orthonormal::{LeafDirections, LeafHarmonics, LeafScratch};
use tetraleaf_math::radial::{bessel_j, bessel_j_real, bessel_j_real_derivatives, modified_i};
use tetraleaf_math::{Complex, Layout, RealScalar};
use tetraleaf_ref::helmholtz::Level;

/// The points per batch: 16, the fastest of 4, 8 and 16 lanes for the routines on Kathleen
/// and locust, or within 5 % (Helmholtz CH4.9a2a, its report §5).
pub const LANES: usize = 16;

const B: usize = LANES;

/// The charges of a leaf's points for [`LeafP2m::p2m`]: point-major, point j's `vectors`
/// vectors one after the other, each one real value or two (Re q, Im q) if `complex`, as the
/// source chunks of CONVENTIONS §3.13 and `tetraleaf-exec`'s Helmholtz densities hold them:
/// vector v's real part of point j at j c m + c v, its imaginary part after it, c = 1 or 2.
#[derive(Clone, Copy, Debug)]
pub struct LeafCharges<'a, T> {
    /// The values, c m per point.
    pub values: &'a [T],
    /// The number m ≥ 1 of charge vectors.
    pub vectors: usize,
    /// Whether each charge is complex (two values) or real (one).
    pub complex: bool,
}

impl<T> LeafCharges<'_, T> {
    /// The values per point, c m.
    pub fn columns(&self) -> usize {
        if self.complex {
            2 * self.vectors
        } else {
            self.vectors
        }
    }
}

/// The batched P2M of a leaf up to a degree (the [module documentation](self)): the leaf
/// harmonics' table. It allocates once, in [`LeafP2m::new`], and is read-only after, so one can
/// serve every thread.
#[derive(Clone, Debug)]
pub struct LeafP2m<T> {
    harmonics: LeafHarmonics<T>,
}

impl<T: RealScalar> LeafP2m<T> {
    /// The kernel for degrees up to `p`.
    pub fn new(p: usize) -> Self {
        Self {
            harmonics: LeafHarmonics::new(p),
        }
    }

    /// The largest degree.
    pub fn p(&self) -> usize {
        self.harmonics.p()
    }

    /// The planes of a multipole of `level` with these charges: 1 for real charges with real k
    /// or on Yukawa's real path (k = iκ; a real multipole, CONVENTIONS §3.15's storage), 2
    /// otherwise: a complex multipole, or on Yukawa the real multipoles of Re q and Im q (the
    /// [module documentation](self), "Yukawa").
    pub fn planes(level: &Level<T>, complex: bool) -> usize {
        if (level.k().im == T::zero() || level.is_yukawa()) && !complex {
            1
        } else {
            2
        }
    }

    /// P2M of a leaf: adds the stored multipoles of degree `p` of the points `points` (offsets
    /// from the leaf's centre in the units of `level`, whose half-width is the leaf's) with the
    /// charges `charges` to `multipoles`, for every charge vector:
    /// M̂ₙᵐ += Σⱼ qⱼ ĵₙ(kρⱼ) e^{Im(k)ρⱼ − β} Yₙᵐ(ρ̂ⱼ), [`planes`](Self::planes) (p + 1)²
    /// values per vector in the layout of §3.6, the real plane then the imaginary plane, the
    /// vectors one after the other. The sum of
    /// [`tetraleaf_ref::helmholtz::leaf::p2m`] over the points, to the contract of the [module
    /// documentation](self). On Yukawa's real path (k = iκ, Helmholtz CH2.1c) it is the sum of
    /// [`tetraleaf_ref::helmholtz::leaf::p2m_yukawa`], M̂ᵞₙᵐ += Σⱼ qⱼ îₙ(κρⱼ) e^{κρⱼ − β} Yₙᵐ(ρ̂ⱼ),
    /// per real charge column: one plane per vector for real charges, the planes of Re q and Im q
    /// for complex ones.
    ///
    /// # Panics
    ///
    /// If `p` exceeds the kernel's or the scratch's degree, the charges have more vectors than
    /// the scratch, `charges.values` does not hold [`columns`](LeafCharges::columns) values per
    /// point, or `multipoles` does not have m [`planes`](Self::planes) (p + 1)² values.
    pub fn p2m(
        &self,
        p: usize,
        level: &Level<T>,
        points: &[[T; 3]],
        charges: LeafCharges<'_, T>,
        scratch: &mut P2mScratch<T>,
        multipoles: &mut [T],
    ) {
        assert!(p <= self.p(), "degree {p} above the kernel's {}", self.p());
        assert!(
            p <= scratch.p,
            "degree {p} above the scratch's {}",
            scratch.p
        );
        assert!(
            charges.vectors >= 1 && charges.vectors <= scratch.vectors,
            "{} charge vectors, the scratch holds 1 to {}",
            charges.vectors,
            scratch.vectors
        );
        let columns = charges.columns();
        assert_eq!(
            charges.values.len(),
            columns * points.len(),
            "the charges must hold {columns} values per point"
        );
        let len = Layout::new(p).len();
        let planes = Self::planes(level, charges.complex);
        assert_eq!(
            multipoles.len(),
            charges.vectors * planes * len,
            "the multipoles must hold {} values",
            charges.vectors * planes * len
        );
        let k = level.k();
        let real_k = k.im == T::zero();
        let yukawa = level.is_yukawa();
        // Real weights, wₙ = (q e) ĵₙ or (q e) îₙ per real charge column.
        let real_radial = real_k || yukawa;
        let (s, beta) = (level.scale(), level.damping());
        let tau = under_range::<T>();
        // For real k the factor is e^{0·ρ − β}, the same for every point.
        let e_real = (T::zero() - beta).exp();
        let P2mScratch {
            harmonics: harmonic_scratch,
            radial,
            radial_im,
            one,
            one_real,
            weights,
            ..
        } = scratch;
        let degrees = p + 1;
        for start in (0..points.len()).step_by(B) {
            let count = B.min(points.len() - start);
            let lanes: [[T; 3]; B] = array::from_fn(|j| {
                if j < count {
                    points[start + j]
                } else {
                    [T::zero(); 3]
                }
            });
            let rho: [T; B] = array::from_fn(|j| {
                (lanes[j][0] * lanes[j][0] + lanes[j][1] * lanes[j][1] + lanes[j][2] * lanes[j][2])
                    .sqrt()
            });
            // The radial functions and the factor.
            let e: [T; B] = if real_k {
                let x: [T; B] = array::from_fn(|j| k.re * rho[j]);
                bessel_j_real(&x, s, &mut radial[..degrees]);
                [e_real; B]
            } else if yukawa {
                // îₙ(κρ) per lane by the reference's routine and length, and e = e^{κρ − β} per
                // lane, the values below τ flushed (module documentation, "Yukawa").
                let kappa = k.im;
                for j in 0..B {
                    modified_i(kappa * rho[j], s, &mut one_real[..degrees], None);
                    for (n, &value) in one_real[..degrees].iter().enumerate() {
                        radial[n][j] = if value < tau { T::zero() } else { value };
                    }
                }
                array::from_fn(|j| {
                    let e = (kappa * rho[j] - beta).exp();
                    if e < tau { T::zero() } else { e }
                })
            } else {
                for j in 0..B {
                    bessel_j(k * rho[j], s, &mut one[..degrees], None);
                    for (n, value) in one[..degrees].iter().enumerate() {
                        radial[n][j] = value.re;
                        radial_im[n][j] = value.im;
                    }
                }
                array::from_fn(|j| (k.im * rho[j] - beta).exp())
            };
            // The weights per vector and degree: plane 0 Re, plane 1 Im.
            for v in 0..charges.vectors {
                let charge = |j: usize, part: usize| {
                    if j < count {
                        charges.values
                            [(start + j) * columns + v * (columns / charges.vectors) + part]
                    } else {
                        T::zero()
                    }
                };
                let w_re: [T; B] = array::from_fn(|j| charge(j, 0) * e[j]);
                let w_im: [T; B] = if charges.complex {
                    array::from_fn(|j| charge(j, 1) * e[j])
                } else {
                    [T::zero(); B]
                };
                let out = &mut weights[2 * v * degrees..2 * (v + 1) * degrees];
                let (re, im) = out.split_at_mut(degrees);
                for n in 0..degrees {
                    let j_re = &radial[n];
                    if real_radial {
                        re[n] = array::from_fn(|j| w_re[j] * j_re[j]);
                        if planes == 2 {
                            im[n] = array::from_fn(|j| w_im[j] * j_re[j]);
                        }
                    } else {
                        // (q e) ĵₙ as the reference's complex product.
                        let j_im = &radial_im[n];
                        re[n] = array::from_fn(|j| w_re[j] * j_re[j] - w_im[j] * j_im[j]);
                        im[n] = array::from_fn(|j| w_re[j] * j_im[j] + w_im[j] * j_re[j]);
                    }
                }
            }
            let weights = &*weights;
            self.harmonics
                .columns(p, &lanes, harmonic_scratch, |m, column, cos, sin| {
                    for (n, q) in (m..=p).zip(column) {
                        let centre = n * n + n;
                        let y_cos: [T; B] = if m == 0 {
                            *q
                        } else {
                            array::from_fn(|j| q[j] * cos[j])
                        };
                        let y_sin: [T; B] = array::from_fn(|j| q[j] * sin[j]);
                        for v in 0..charges.vectors {
                            let w = &weights[2 * v * degrees..2 * (v + 1) * degrees];
                            let multipole =
                                &mut multipoles[v * planes * len..(v + 1) * planes * len];
                            for plane in 0..planes {
                                let wn = &w[plane * degrees + n];
                                let coefficients = &mut multipole[plane * len..(plane + 1) * len];
                                let i = centre + m;
                                coefficients[i] = coefficients[i] + dot(wn, &y_cos);
                                if m > 0 {
                                    let i = centre - m;
                                    coefficients[i] = coefficients[i] + dot(wn, &y_sin);
                                }
                            }
                        }
                    }
                });
        }
    }
}

/// The under-range threshold τ of `tetraleaf_math::orthonormal::leaf` ("Under-range values"),
/// restated (it is private there): 2⁻⁹⁰⁰ for a type with f64's range, 2⁻¹⁰⁰ for f32's. Yukawa's
/// îₙ below it are flushed to zero, as `bessel_j_real` flushes ĵₙ.
#[inline]
fn under_range<T: RealScalar>() -> T {
    if T::max_value().to_f64() > 1e300 {
        T::from_f64(2f64.powi(-900))
    } else {
        T::from_f64(2f64.powi(-100))
    }
}

/// Σⱼ aⱼ bⱼ over the lanes, the products added pairwise in a fixed order: lane j with j + 8,
/// then j with j + 4, j + 2, j + 1.
#[inline(always)]
fn dot<T: RealScalar>(a: &[T; B], b: &[T; B]) -> T {
    let eight: [T; 8] = array::from_fn(|j| a[j] * b[j] + a[j + 8] * b[j + 8]);
    let four: [T; 4] = array::from_fn(|j| eight[j] + eight[j + 4]);
    let two: [T; 2] = array::from_fn(|j| four[j] + four[j + 2]);
    two[0] + two[1]
}

/// The work space of [`LeafP2m::p2m`] up to a degree and a number of charge vectors: the
/// harmonics' column, the radial functions and the weights of one batch. It allocates once,
/// in [`P2mScratch::new`]; its size grows with p and m, not with the leaf's points.
#[derive(Clone, Debug)]
pub struct P2mScratch<T> {
    p: usize,
    vectors: usize,
    harmonics: LeafScratch<T, B>,
    /// ĵₙ per lane (its real part for complex k), n ≤ p.
    radial: Vec<[T; B]>,
    /// Im ĵₙ per lane for complex k.
    radial_im: Vec<[T; B]>,
    /// One lane's ĵₙ for complex k.
    one: Vec<Complex<T>>,
    /// One lane's îₙ on Yukawa's real path.
    one_real: Vec<T>,
    /// Per vector the weights' real then imaginary planes, p + 1 rows each.
    weights: Vec<[T; B]>,
}

impl<T: RealScalar> P2mScratch<T> {
    /// Work space for degrees up to `p` and up to `vectors` charge vectors.
    pub fn new(p: usize, vectors: usize) -> Self {
        let degrees = p + 1;
        Self {
            p,
            vectors,
            harmonics: LeafScratch::new(p),
            radial: vec![[T::zero(); B]; degrees],
            radial_im: vec![[T::zero(); B]; degrees],
            one: vec![Complex::new(T::zero(), T::zero()); degrees],
            one_real: vec![T::zero(); degrees],
            weights: vec![[T::zero(); B]; 2 * vectors * degrees],
        }
    }

    /// The largest degree.
    pub fn p(&self) -> usize {
        self.p
    }

    /// The most charge vectors.
    pub fn vectors(&self) -> usize {
        self.vectors
    }
}

/// A leaf's target output for [`LeafL2p::l2p`]: two columns per charge vector (Re, then Im), as
/// `tetraleaf-exec`'s Helmholtz densities lay out a target output chunk of n targets (CONVENTIONS
/// §3.13's chunk of 2m output vectors): column o's potential of target i at o n + i, and with
/// gradients its gradient's component a at 2m n + 3 (o n + i) + a, column o = 2v + part for vector
/// v. The kernel adds (+=) into it.
#[derive(Debug)]
pub struct LeafOutput<'a, T> {
    /// The values, 2m n, or 8m n with gradients.
    pub values: &'a mut [T],
    /// The number m ≥ 1 of charge vectors.
    pub vectors: usize,
    /// Whether the chunk holds gradients after the potentials.
    pub gradients: bool,
}

impl<T> LeafOutput<'_, T> {
    /// The values per target: 2m, or 8m with gradients.
    pub fn per_target(&self) -> usize {
        2 * self.vectors * if self.gradients { 4 } else { 1 }
    }
}

/// The batched L2P of a leaf up to a degree (the [module documentation](self), "L2P"): the leaf
/// harmonics' table. It allocates once, in [`LeafL2p::new`], and is read-only after, so one can
/// serve every thread.
#[derive(Clone, Debug)]
pub struct LeafL2p<T> {
    harmonics: LeafHarmonics<T>,
}

/// The per-degree sums of [`LeafL2p::l2p`], each `p + 1` rows of lanes per charge vector: the
/// value's real and imaginary planes, then ∂θ's and (1/sin θ) ∂φ's.
const SUMS: usize = 6;

impl<T: RealScalar> LeafL2p<T> {
    /// The kernel for degrees up to `p`.
    pub fn new(p: usize) -> Self {
        Self {
            harmonics: LeafHarmonics::new(p),
        }
    }

    /// The largest degree.
    pub fn p(&self) -> usize {
        self.harmonics.p()
    }

    /// L2P of a leaf: adds the potentials of the stored locals `locals` of degree `p` of `level`
    /// (whose half-width is the leaf's) at the targets `targets` (offsets from the leaf's centre in
    /// the level's units), and with `output.gradients` their gradients, into `output`, for every
    /// charge vector: φ̃ += Σ L̂ₙᵐ ĵₙ(k|v|) e^{Im(k)|v| − β} Yₙᵐ(v̂). `locals` holds per vector one
    /// complex local, 2 (p + 1)² values in the layout of §3.6, the real plane then the imaginary
    /// plane, the vectors one after the other. The sum of
    /// [`tetraleaf_ref::helmholtz::leaf::l2p`] per target, to the contract of the [module
    /// documentation](self), "L2P".
    ///
    /// # Panics
    ///
    /// If `p` exceeds the kernel's or the scratch's degree, the output has more vectors than the
    /// scratch or none, the level is Yukawa's (k = iκ: H2's path), `locals` does not hold
    /// 2m (p + 1)² values, or `output.values` does not hold
    /// [`per_target`](LeafOutput::per_target) values per target.
    pub fn l2p(
        &self,
        p: usize,
        level: &Level<T>,
        targets: &[[T; 3]],
        locals: &[T],
        scratch: &mut L2pScratch<T>,
        output: LeafOutput<'_, T>,
    ) {
        assert!(p <= self.p(), "degree {p} above the kernel's {}", self.p());
        assert!(
            p <= scratch.p,
            "degree {p} above the scratch's {}",
            scratch.p
        );
        let vectors = output.vectors;
        assert!(
            vectors >= 1 && vectors <= scratch.vectors,
            "{vectors} charge vectors, the scratch holds 1 to {}",
            scratch.vectors
        );
        assert!(
            !level.is_yukawa(),
            "the batched L2P takes Helmholtz's complex path, not Yukawa's"
        );
        let len = Layout::new(p).len();
        assert_eq!(
            locals.len(),
            2 * vectors * len,
            "the locals must hold {} values",
            2 * vectors * len
        );
        let n = targets.len();
        assert_eq!(
            output.values.len(),
            output.per_target() * n,
            "the output must hold {} values per target",
            output.per_target()
        );
        let gradients = output.gradients;
        let k = level.k();
        let real_k = k.im == T::zero();
        let (s, beta) = (level.scale(), level.damping());
        let e_real = (T::zero() - beta).exp();
        // The tangential factor's limit at the centre, k ĵ₁(z)/z → k/(3s).
        let at_centre = k * (T::one() / (T::from_f64(3.0) * s));
        let degrees = p + 1;
        let L2pScratch {
            harmonics: harmonic_scratch,
            radial,
            radial_im,
            derivative,
            derivative_im,
            one,
            one_derivative,
            sums,
            ..
        } = scratch;
        let (radial, radial_im) = (&mut radial[..degrees], &mut radial_im[..degrees]);
        let (derivative, derivative_im) =
            (&mut derivative[..degrees], &mut derivative_im[..degrees]);
        if real_k {
            radial_im.fill([T::zero(); B]);
            derivative_im.fill([T::zero(); B]);
        }
        let sums = &mut sums[..vectors * SUMS * degrees];
        let columns = 2 * vectors;
        for start in (0..n).step_by(B) {
            let count = B.min(n - start);
            let lanes: [[T; 3]; B] = array::from_fn(|j| {
                if j < count {
                    targets[start + j]
                } else {
                    [T::zero(); 3]
                }
            });
            let r: [T; B] = array::from_fn(|j| {
                (lanes[j][0] * lanes[j][0] + lanes[j][1] * lanes[j][1] + lanes[j][2] * lanes[j][2])
                    .sqrt()
            });
            // The radial functions (and derivatives) and the factor.
            let e: [T; B] = if real_k {
                let x: [T; B] = array::from_fn(|j| k.re * r[j]);
                if gradients {
                    bessel_j_real_derivatives(&x, s, radial, derivative);
                } else {
                    bessel_j_real(&x, s, radial);
                }
                [e_real; B]
            } else {
                for j in 0..B {
                    let d = gradients.then_some(&mut one_derivative[..degrees]);
                    bessel_j(k * r[j], s, &mut one[..degrees], d);
                    for n in 0..degrees {
                        radial[n][j] = one[n].re;
                        radial_im[n][j] = one[n].im;
                        if gradients {
                            derivative[n][j] = one_derivative[n].re;
                            derivative_im[n][j] = one_derivative[n].im;
                        }
                    }
                }
                array::from_fn(|j| (k.im * r[j] - beta).exp())
            };
            // The per-degree sums of every vector.
            sums.fill([T::zero(); B]);
            if gradients {
                self.harmonics.derivative_columns(
                    p,
                    &lanes,
                    harmonic_scratch,
                    |m, q, d, u, cos, sin| {
                        let mt = T::from_f64(m as f64);
                        for (row, n) in (m..=p).enumerate() {
                            let centre = n * n + n;
                            let (q, d) = (&q[row], &d[row]);
                            let mu: [T; B] = array::from_fn(|j| mt * u[row][j]);
                            for v in 0..vectors {
                                let local = &locals[2 * v * len..2 * (v + 1) * len];
                                let acc = &mut sums[v * SUMS * degrees..(v + 1) * SUMS * degrees];
                                for plane in 0..2 {
                                    let c = &local[plane * len..(plane + 1) * len];
                                    let (value, rest) = acc.split_at_mut(2 * degrees);
                                    let (theta, phi) = rest.split_at_mut(2 * degrees);
                                    let at = plane * degrees + n;
                                    if m == 0 {
                                        let c0 = c[centre];
                                        add_product(&mut value[at], q, c0);
                                        add_product(&mut theta[at], d, c0);
                                    } else {
                                        let (plus, minus) = (c[centre + m], c[centre - m]);
                                        let a: [T; B] =
                                            array::from_fn(|j| plus * cos[j] + minus * sin[j]);
                                        let b: [T; B] =
                                            array::from_fn(|j| minus * cos[j] - plus * sin[j]);
                                        add_lanes(&mut value[at], q, &a);
                                        add_lanes(&mut theta[at], d, &a);
                                        add_lanes(&mut phi[at], &mu, &b);
                                    }
                                }
                            }
                        }
                    },
                );
            } else {
                self.harmonics
                    .columns(p, &lanes, harmonic_scratch, |m, q, cos, sin| {
                        for (row, n) in (m..=p).enumerate() {
                            let centre = n * n + n;
                            let q = &q[row];
                            for v in 0..vectors {
                                let local = &locals[2 * v * len..2 * (v + 1) * len];
                                let acc = &mut sums[v * SUMS * degrees..(v + 1) * SUMS * degrees];
                                for plane in 0..2 {
                                    let c = &local[plane * len..(plane + 1) * len];
                                    let at = plane * degrees + n;
                                    if m == 0 {
                                        add_product(&mut acc[at], q, c[centre]);
                                    } else {
                                        let (plus, minus) = (c[centre + m], c[centre - m]);
                                        let a: [T; B] =
                                            array::from_fn(|j| plus * cos[j] + minus * sin[j]);
                                        add_lanes(&mut acc[at], q, &a);
                                    }
                                }
                            }
                        }
                    });
            }
            // Per lane the reference's degree loop, then the output.
            let directions = gradients.then(|| LeafDirections::of(&lanes));
            let inverse: [T; B] = array::from_fn(|j| {
                if r[j] > T::zero() {
                    T::one() / r[j]
                } else {
                    T::zero()
                }
            });
            for v in 0..vectors {
                let acc = &sums[v * SUMS * degrees..(v + 1) * SUMS * degrees];
                let row =
                    |which: usize, plane: usize, n: usize| &acc[(2 * which + plane) * degrees + n];
                let zero = [T::zero(); B];
                let (mut phi_re, mut phi_im) = (zero, zero);
                for n in 0..degrees {
                    let (a, b) = (row(0, 0, n), row(0, 1, n));
                    let (jr, ji) = (&radial[n], &radial_im[n]);
                    for j in 0..B {
                        phi_re[j] = phi_re[j] + (jr[j] * a[j] - ji[j] * b[j]);
                        phi_im[j] = phi_im[j] + (jr[j] * b[j] + ji[j] * a[j]);
                    }
                }
                let (re, im) = (2 * v, 2 * v + 1);
                for j in 0..count {
                    let i = start + j;
                    output.values[re * n + i] = output.values[re * n + i] + phi_re[j] * e[j];
                    output.values[im * n + i] = output.values[im * n + i] + phi_im[j] * e[j];
                }
                let Some(directions) = &directions else {
                    continue;
                };
                // The radial, polar and azimuthal parts, each complex per lane.
                let mut parts = [[zero; 2]; 3];
                for n in 0..degrees {
                    let (dr, di) = (&derivative[n], &derivative_im[n]);
                    let (jr, ji) = (&radial[n], &radial_im[n]);
                    let value = (row(0, 0, n), row(0, 1, n));
                    let theta = (row(1, 0, n), row(1, 1, n));
                    let phi = (row(2, 0, n), row(2, 1, n));
                    // tₙ = ĵₙ (1/r), or the centre's limit (n = 1) or zero.
                    let (limit_re, limit_im) = if n == 1 {
                        (at_centre.re, at_centre.im)
                    } else {
                        (T::zero(), T::zero())
                    };
                    let t_re: [T; B] = array::from_fn(|j| {
                        if r[j] > T::zero() {
                            jr[j] * inverse[j]
                        } else {
                            limit_re
                        }
                    });
                    let t_im: [T; B] = array::from_fn(|j| {
                        if r[j] > T::zero() {
                            ji[j] * inverse[j]
                        } else {
                            limit_im
                        }
                    });
                    let terms = [
                        ((dr, di), value),
                        ((&t_re, &t_im), theta),
                        ((&t_re, &t_im), phi),
                    ];
                    for (part, ((fr, fi), (a, b))) in parts.iter_mut().zip(terms) {
                        for j in 0..B {
                            part[0][j] = part[0][j] + (fr[j] * a[j] - fi[j] * b[j]);
                            part[1][j] = part[1][j] + (fr[j] * b[j] + fi[j] * a[j]);
                        }
                    }
                }
                let [radial_part, polar, azimuthal] = parts;
                // k times the radial part.
                let kr: [[T; B]; 2] = [
                    array::from_fn(|j| k.re * radial_part[0][j] - k.im * radial_part[1][j]),
                    array::from_fn(|j| k.re * radial_part[1][j] + k.im * radial_part[0][j]),
                ];
                let LeafDirections {
                    cos_theta: ct,
                    sin_theta: st,
                    cos_phi: cp,
                    sin_phi: sp,
                } = *directions;
                let base = columns * n;
                for (column, plane) in [(re, 0), (im, 1)] {
                    let (kr, pol, az) = (&kr[plane], &polar[plane], &azimuthal[plane]);
                    for j in 0..count {
                        let g = [
                            kr[j] * (st[j] * cp[j]) + pol[j] * (ct[j] * cp[j]) + az[j] * (-sp[j]),
                            kr[j] * (st[j] * sp[j]) + pol[j] * (ct[j] * sp[j]) + az[j] * cp[j],
                            kr[j] * ct[j] + pol[j] * (-st[j]),
                        ];
                        let at = base + 3 * (column * n + start + j);
                        for (axis, g) in g.into_iter().enumerate() {
                            output.values[at + axis] = output.values[at + axis] + g * e[j];
                        }
                    }
                }
            }
        }
    }
}

/// acc += a c per lane.
#[inline(always)]
fn add_product<T: RealScalar>(acc: &mut [T; B], a: &[T; B], c: T) {
    for j in 0..B {
        acc[j] = acc[j] + a[j] * c;
    }
}

/// acc += a b per lane.
#[inline(always)]
fn add_lanes<T: RealScalar>(acc: &mut [T; B], a: &[T; B], b: &[T; B]) {
    for j in 0..B {
        acc[j] = acc[j] + a[j] * b[j];
    }
}

/// The work space of [`LeafL2p::l2p`] up to a degree and a number of charge vectors: the harmonics'
/// columns, the radial functions and derivatives and the per-degree sums of one batch. It allocates
/// once, in [`L2pScratch::new`]; its size grows with p and m, not with the leaf's targets.
#[derive(Clone, Debug)]
pub struct L2pScratch<T> {
    p: usize,
    vectors: usize,
    harmonics: LeafScratch<T, B>,
    /// ĵₙ per lane, its real and imaginary parts (zero for real k), n ≤ p.
    radial: Vec<[T; B]>,
    radial_im: Vec<[T; B]>,
    /// ĵₙ′ per lane, its real and imaginary parts.
    derivative: Vec<[T; B]>,
    derivative_im: Vec<[T; B]>,
    /// One lane's ĵₙ and ĵₙ′ for complex k.
    one: Vec<Complex<T>>,
    one_derivative: Vec<Complex<T>>,
    /// Per vector the [`SUMS`] per-degree sums, p + 1 rows each.
    sums: Vec<[T; B]>,
}

impl<T: RealScalar> L2pScratch<T> {
    /// Work space for degrees up to `p` and up to `vectors` charge vectors.
    pub fn new(p: usize, vectors: usize) -> Self {
        let degrees = p + 1;
        let rows = vec![[T::zero(); B]; degrees];
        let zero = Complex::new(T::zero(), T::zero());
        Self {
            p,
            vectors,
            harmonics: LeafScratch::new(p),
            radial: rows.clone(),
            radial_im: rows.clone(),
            derivative: rows.clone(),
            derivative_im: rows,
            one: vec![zero; degrees],
            one_derivative: vec![zero; degrees],
            sums: vec![[T::zero(); B]; SUMS * vectors * degrees],
        }
    }

    /// The largest degree.
    pub fn p(&self) -> usize {
        self.p
    }

    /// The most charge vectors.
    pub fn vectors(&self) -> usize {
        self.vectors
    }
}
