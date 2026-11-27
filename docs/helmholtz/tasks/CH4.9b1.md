# CH4.9b1 The leaf harmonics' derivatives and the real radial derivatives, batched over points

Started 2026-11-27 on `tbetcke/h-ch4-9b1`. Roadmap: CH4.9b1, the first half of CH4.9b (item 14 of
the priority list; CH4.9 added on 2026-11-26, the user: "Yes. Add P2M and L2P to the priority list
as tasks"). No proposal §10 entry: the acceptance is the roadmap's *Done when* with constraint 1.
Proposal §2 (constraints 1, 6, 9); CONVENTIONS §3.15 ("Angular basis", "Radial functions and
scaling"); CH4.9a1's report `tetraleaf-validate/results/helmholtz-ch4-9a1.md` (L2P with gradients:
the harmonics with their derivatives 58–84 % of a point), CH4.9a2a's
`tetraleaf-validate/results/helmholtz-ch4-9a2a.md` and its record `docs/helmholtz/tasks/CH4.9a2a.md`.

## Size: CH4.9b split

CH4.9b as written holds what CH4.9a2 held for P2M: routines in `tetraleaf-math` (here the
derivatives that L2P's gradient needs: ∂θ and (1/sin θ) ∂φ of the leaf harmonics, ĵₙ′ of the real
radial functions), a kernel in `tetraleaf-tables`, its dispatch in `tetraleaf-exec` with the FMM's
checks, and times per point and per level on Kathleen and locust, besides the evaluation before and
after CH4.9a–b and the ISA question. CH4.9a2a and CH4.9a2b together were about 3,000 lines; this is
the skill's rule for a split, along the same line (CH4.9a1's question 2: the routines first, then
the kernel):

- **CH4.9b1** (this plan): the derivative routines in `tetraleaf-math`, tested against
  `harmonics_derivatives` and `bessel_j`'s derivatives, timed per point against the reference's L2P
  parts, with L2P with gradients projected per point.
- **CH4.9b2**: the L2P kernel in `tetraleaf-tables::helmholtz::leaf` on them, the operator's
  dispatch with `reference_leaves`, the FMM's checks, L2P per level and the evaluation before and
  after CH4.9a–b on both machines, the ISA question and the f32 sum.

Both lines, with *Done when* and *Needs*, are in the roadmap (committed with this plan); CH4.9b's
original *Done when* stays as the parent's text.

## Goal
Make the parts that dominate a leaf point's L2P with gradients (the harmonics with their angular
derivatives, the radial functions with theirs) cheap and vectorisable over a batch of points, and
measure how far they fall, before the kernel.

## Acceptance (the roadmap's *Done when*)
`LeafHarmonics` gives per order the column Q̄ₙᵐ (`columns`' bits), ∂θ Q̄ₙᵐ and Q̄ₙᵐ / sin θ of a
batch as lanes, without a division by sin θ, held per term to a contract against
`harmonics_derivatives` over the leaf frame's directions (the poles, the centre and the axes
included) in f64 to p = 100 and f32 to p = 64, no subnormal output, lanes and positions bit for bit;
`bessel_j_real`'s derivatives (`bessel_j`'s derivatives' real parts bit for bit in f64, in f32 but
for flushed values); `helmholtz_leaf --lanes` times them per point beside the reference's L2P parts
with the reference's sum on their values (the projected L2P with gradients) on one core of Kathleen
and locust; no `tetraleaf-simd` change (constraint 1).

## Read first
- `tetraleaf-math/src/orthonormal/leaf.rs` (`LeafHarmonics::columns`, the recurrence and the
  under-range rule) and `harmonics.rs` (`harmonics_derivatives`: ∂θ by the ladder over orders
  m ± 1, (1/sin θ) ∂φ from P̄ₙᵐ / sin θ, the frame φ = 0 on the z-axis).
- `tetraleaf-math/src/radial.rs` (`bessel_j`'s derivatives from the neighbouring orders and
  `r_top`; `bessel_j_real`).
- `tetraleaf-ref/src/helmholtz/leaf.rs` (`l2p`, `evaluate` with gradients, `Polar::cartesian`, the
  tangential limit at the centre): what CH4.9b2's kernel restates and what the example's projected
  sum restates now.
- `tetraleaf-validate/examples/helmholtz_leaf.rs` (CH4.9a1–CH4.9a2b).

## Changes
- `tetraleaf-math/src/orthonormal/leaf.rs`:
  - The derivatives. With Uₙᵐ = Q̄ₙᵐ / sin θ (m ≥ 1), the column recurrence is linear, so U obeys
    it too, from Uₘᵐ = dₘ Q̄ₘ₋₁ᵐ⁻¹ (the diagonal before its factor sin θ) and
    Uₘ₊₁ᵐ = √(2m + 3) cos θ Uₘᵐ: finite at the poles, no division by sin θ. From
    sin θ ∂θ Pₙᵐ = n cos θ Pₙᵐ − (n + m) Pₙ₋₁ᵐ in the normalisation of the columns,
    ∂θ Q̄ₙᵐ = n cos θ Uₙᵐ − cₙₘ Uₙ₋₁ᵐ with cₙₘ = √((2n + 1)(n² − m²)/(2n − 1)) (m ≥ 1; cₘₘ = 0), and
    for m = 0 the ladder ∂θ Q̄ₙ⁰ = −√(n(n + 1)/2) Q̄ₙ¹ (the √2 of fₘ). Then ∂θ Yₙ^{±m} = ∂θ Q̄ₙᵐ
    {cos mφ, sin mφ} and (1/sin θ) ∂φ Yₙ^{±m} = m Uₙᵐ {−sin mφ, cos mφ}, as
    `harmonics_derivatives` gives them (the frame φ = 0 on the z-axis, θ = 0 at the centre).
    `LeafHarmonics::new` also tabulates cₙₘ, n and √(n(n + 1)/2) in f64, rounded to T once.
  - `LeafHarmonics::derivative_columns::<B>(p, points, scratch, f)`: per order m the column
    q[n − m] = Q̄ₙᵐ (the same operations as `columns`, so its bits), d[n − m] = ∂θ Q̄ₙᵐ and
    u[n − m] = Uₙᵐ (all zero for m = 0, whose terms have no φ derivative), and (cos mφ, sin mφ),
    handed to `f(m, q, d, u, cos, sin)`. Order 1's column is formed before order 0 is handed out
    (order 0's ∂θ reads it). Under-range: U's diagonal is set to zero below τ as Q̄'s is (U is not
    zeroed because Q̄ is: at the pole Q̄₁¹ = 0 but U₁¹ = d₁ Y₀⁰ carries ∂θ Y₁¹); the bound on what
    that drops goes into the module docs beside Q̄'s.
  - `LeafHarmonics::derivatives::<B>(p, points, scratch, value, d_theta, d_phi)`: the layout of
    §3.6 per lane, `harmonics_derivatives`' three outputs, from `derivative_columns`; for tests and
    simple callers.
  - `LeafScratch` gains the buffers (U, ∂θ, order 1's columns), allocated in `new`.
- `tetraleaf-math/src/radial.rs`: `bessel_j_real_derivatives::<T, B>(x, s, values, derivatives)`:
  `bessel_j_real`'s values (its bits) and ĵₙ′ per lane by `bessel_j`'s formula from the neighbouring
  orders, ĵ_len from the ratio r_len that the backward recursion passes (as `bessel_j`'s `r_top`),
  so each lane's derivatives are `bessel_j`'s derivatives' real parts bit for bit in f64. In f32 a
  derivative beyond the turning point below τ is zeroed as the values are (no subnormal output);
  one formed from a flushed neighbour differs from the reference by at most about (n + 1) τ / s, far
  below u. `bessel_j_real` itself keeps its code path and bits (a shared private routine).
- `tetraleaf-math/CLAUDE.md`: a paragraph (CH4.9b1).
- `tetraleaf-validate/examples/helmholtz_leaf.rs`: per case also `derivative_columns` (every
  column summed into a sink, as `columns` is timed) and `derivatives` over the batch,
  `bessel_j_real_derivatives` (real k), each checked first (the derivatives' worst against
  `harmonics_derivatives` in the contract's units; the radial derivatives against `bessel_j`'s bits
  but for flushed values); "sum ∇, leaf values": the reference's L2P sum with gradients (`evaluate`'s
  loop and `Polar::cartesian`, restated as `sum` restates `form`) on the routines' values of the
  batch's first 16 points cycled over the targets; and the projected L2P with gradients:
  `derivative_columns` + the radial functions with derivatives + the factor's `exp` + that sum,
  over the reference's L2P with gradients. `tetraleaf-validate/CLAUDE.md`: a sentence.
- `tools/helmholtz/ch4-9b1-runs.sh` (`smoke`, `parts`, as CH4.9a2a's), `ch4-9b1-tables.py` on
  `common.py`.
- Report `tetraleaf-validate/results/helmholtz-ch4-9b1.md`.
- Conventions: none (the derivatives are those of §3.15's harmonics and `harmonics_derivatives`'
  outputs; only the recurrence that forms them is new). No dependency. No Laplace path:
  `tetraleaf-math`'s Laplace modules untouched; the new items are additions to the Helmholtz
  modules (`orthonormal::leaf`, `radial`).

## Tests that prove it
- `tetraleaf-math/tests/orthonormal/leaf.rs` (CH4.9a2a's directions: 1,000 random points of
  [−1, 1)³, the poles, near-pole points with sin θ from 10⁻¹ to 10⁻³⁰⁰, the axes, the equator, the
  plane x = 0, the centre, the cube's corners):
  - `derivatives` against `harmonics_derivatives` per point and term, f64, p ∈ {0, 1, 2, 5, 8, 17,
    32, 45, 64, 100}: the values as before; ∂θ and (1/sin θ) ∂φ within
    C (n + 1)³ u √((2n + 1)/(4π)) (a first guess, one power of n above the values' contract for
    the factor n of the derivative; the measured worst printed, C set to about ten times it as for
    the values; the first measurement decides, recorded under "As built").
  - f32 against the f64 reference, p ≤ 64, the same contract with u = 2⁻²⁴ plus the flush's floor;
    no subnormal in any output (asserted).
  - `derivative_columns`' q bit for bit `columns`' columns; its d and u give `derivatives`' layout
    bit for bit; lanes B = 1, 4, 8, 16 and shifted positions the same bits per point.
  - At the north and south poles and the centre: ∂θ Y₁¹ and (1/sin θ) ∂φ Y₁⁻¹ equal √(3/(4π))
    (signs by the pole, as `harmonics_derivatives`' doctest), the others of degree 1 zero.
- `tetraleaf-math/tests/radial/real.rs`:
  - `bessel_j_real_derivatives` against `bessel_j(Complex::new(x, 0), s, …, Some(…))` bit for bit
    in f64 (values and derivatives), lanes B = 1, 4, 8, CH4.9a2a's arguments (x = 0, 10⁻¹², 10⁻⁶,
    0.3, 0.5, 1, 2.5, 4, 17, 33.3, 60, 300), s ∈ {1, 0.3, 10⁻³}, lengths 1, 2, 9, 46, 101.
  - Its values bit for bit `bessel_j_real`'s.
  - f32: bit for bit but where a value is flushed or a derivative's neighbour was (there within
    (n + 1) τ / s of the reference, and zero beyond the turning point below τ); no subnormal
    output.

## Checks
`/verify` (the root checks, the strict workspace checks; `tetraleaf-math`'s and
`tetraleaf-validate`'s own); the example's smoke run. Not expected: the MPI list (no change to exec,
plan or octree), `cargo deny` (no manifest change).

## Machine runs
Kathleen (one core, `singlenode`, slow nodes excluded) and locust (one core, gated on a quiet
node): `tools/helmholtz/ch4-9b1-runs.sh parts`, p ∈ {8, 12, 15, 17, 19, 26, 32, 41, 45} × k̂ ∈
{0.5, 2, 4, 8, 16}, f64 and f32, one precision and degree per launch (B12). It decides how far the
derivative routines fall against the reference's L2P parts and how much of the projected L2P with
gradients is left to the sum, which CH4.9b2's kernel and its ISA question read. The M6 not timed
(ask first; not needed).

## Constraints and risks
- Constraint 1: no Laplace path; `tetraleaf-math` gains items in its Helmholtz modules only.
  Constraint 6: each routine against its reference (`harmonics_derivatives`, `bessel_j`).
  Constraint 9: nothing depends on data.
- Risk: ∂θ from U loses accuracy where the two terms cancel. n cos θ Uₙ and cₙₘ Uₙ₋₁ are both of
  the size of n |U|, so the error is about n times U's own (the recurrence's n² u): the contract
  above allows (n + 1)³ u; the test measures it near the poles and the equator.
- Risk: the flush of U's diagonal. U carries one power of sin θ less than Q̄, so its column under
  τ holds values below about τ 2⁴³ as Q̄'s (CH4.9a2a's bound), times at most n for ∂θ: still far
  below u in both precisions; the f32 test against the f64 reference checks it.
- Risk: the m = 0 ladder needs order 1 first; a slip there shows in the pole tests.

## Pending runs
(filled in step 6)

## Open questions
(filled when they arise)
