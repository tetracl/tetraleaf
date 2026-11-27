# CH4.9a2a The leaf harmonics and the real radial functions, batched over points

Started 2026-11-27 on `tbetcke/h-ch4-9a2a`. Roadmap: CH4.9a2a, the first half of CH4.9a2 (item 14
of the priority list; CH4.9 added on 2026-11-26, the user: "Yes. Add P2M and L2P to the priority
list as tasks"). No proposal §10 entry: the acceptance is the roadmap's *Done when* with constraint
1. Proposal §2 (constraints 1, 6, 9); CONVENTIONS §3.15 ("Angular basis", "Radial functions and
scaling"); CH4.9a1's report `tetraleaf-validate/results/helmholtz-ch4-9a1.md` §2–§4.

## Size: CH4.9a2 split

CH4.9a2 as written holds four pieces: a shared harmonics routine and a real radial recursion
(basis functions, `tetraleaf-math`), a P2M kernel on them (an operator formula, `tetraleaf-tables`,
whose kernel `tetraleaf-exec` may only reach through it), its dispatch under `P2pChoice` with the
FMM's checks (`tetraleaf-exec`), and times per point and per level on Kathleen and locust. That is
three library crates and a new operator plus measurements on two machines, about 1,800 lines: the
skill's rule for a split. The routines are the part CH4.9a1 found to cost 44–84 % of a point, and
CH4.9b reuses them, so they come first, measured per point:

- **CH4.9a2a** (this plan): `LeafHarmonics` and `bessel_j_real` in `tetraleaf-math`, tested
  against `harmonics` and `bessel_j`, timed per point against the reference's parts.
- **CH4.9a2b**: the P2M kernel in `tetraleaf-tables::helmholtz` on them, the operator's dispatch,
  the FMM's checks, P2M per level on both machines, the ISA question if needed.

Both lines, with *Done when* and *Needs*, are in the roadmap (committed with this plan). CH4.9b
needs CH4.9a2b.

## Goal
Make the two parts that dominate a leaf point's P2M and L2P (the harmonics, the radial functions)
cheap and vectorisable over a batch of points, and measure how far they fall, before the kernel.

## Acceptance (the roadmap's *Done when*)
`tetraleaf_math::orthonormal::LeafHarmonics` (the orthonormal harmonics by the plain three-term
recurrence from coefficients tabulated once in f64 up to a degree and rounded to T, for a batch of
points as lanes, column by column, written to vectorise without `unsafe`; values below an
under-range threshold flushed to zero so that no operand is subnormal) within a contract per term
against `harmonics` over the leaf frame's directions (random points of [−1, 1)³, the poles,
near-pole points, the axes, the equator, the centre), f64 to p = 100 and f32 to p = 64;
`tetraleaf_math::radial::bessel_j_real` (ĵₙ of real arguments in real arithmetic, lanes with
per-lane start orders) bit for bit `bessel_j` at real z in f64, and in f32 but for the values it
flushes beyond the turning point; the `helmholtz_leaf` example times both per point against the
reference's parts (f64 and f32) on one core of Kathleen and locust, with P2M per point projected
from CH4.9a1's parts; CH4.9a1's f32 slowdown checked on the new routines. No Laplace path, no
dependency.

## Read first
- `tetraleaf-math/src/orthonormal/harmonics.rs` (`harmonics`, `Direction`, `normalisation`,
  `order_factor`, `store`) and `legendre.rs` (`Columns`: the telescoped recurrence, the exponent
  per column).
- `tetraleaf-math/src/radial.rs` (`bessel_j`, `start_order_j`, `j0_j1`, `div`).
- `tetraleaf-ref/src/helmholtz/leaf.rs` (`p2m`, `form`): what CH4.9a2b's kernel restates.
- `tetraleaf-validate/examples/helmholtz_leaf.rs` (CH4.9a1).

## Changes
- `tetraleaf-math/src/orthonormal/leaf.rs` (new; re-exported from `orthonormal`):
  - `LeafHarmonics<T>`: `new(p)` tabulates in f64, rounded to T once: Y₀⁰ = 1/√(4π); per order m
    the diagonal factor (Q̄ₘᵐ = d_m sin θ Q̄ₘ₋₁ᵐ⁻¹, with √2 for m = 1 and the normalisation
    √((2n + 1)/(4π)) folded in), the first step √(2m + 3), and per (n, m), n ≥ m + 2, the
    three-term coefficients a_nm = √((4n² − 1)/(n² − m²)) and
    b_nm = √((2n + 1)((n − 1)² − m²)/((2n − 3)(n² − m²))) of Q̄ₙᵐ = a_nm cos θ Q̄ₙ₋₁ᵐ − b_nm Q̄ₙ₋₂ᵐ.
    Allocates once.
  - `LeafHarmonics::columns::<B>(p, points, scratch, f)`: for B points as lanes (`[[T; 3]; B]`
    offsets from the centre), cos θ, sin θ, cos φ, sin φ per lane (selects at the centre and on the
    z-axis, as `Direction`: θ = 0 at the centre, φ = 0 on the axis; no trigonometric call, §3.5),
    then per order m the column Q̄ₙᵐ, n = m … p, as `[T; B]` rows and the phase (cos mφ, sin mφ),
    handed to `f(m, column, cos, sin)`. Lanes whose diagonal falls below the threshold τ
    (f64 2^−900, f32 2^−100) are set to zero, so a column never holds a subnormal operand from
    the diagonal (CH4.9a1's f32 slowdown); the bound on what that drops is in the module docs.
  - `LeafHarmonics::values::<B>(p, points, scratch, out)`: the full layout of §3.6 per lane
    (`out[idx(n, m)][j]`), from `columns`; for tests and simple callers.
  - `LeafScratch<T, B>`: the column buffer, allocated once.
- `tetraleaf-math/src/radial.rs`: `bessel_j_real::<T, B>(x: &[T; B], s, values: &mut [[T; B]])`,
  ĵₙ(x) of real x ≥ 0 per lane in real arithmetic: per lane the start order of `start_order_j`
  restated in real f64 arithmetic (the same orders: Smith's division of reals is the real quotient),
  the backward ratio recursion run from the largest start order with each lane entering at its own
  (a lane's ratio stays 0 above its start, as the reference starts at 0), the zero-denominator guard
  per lane, ĵ₀ or ĵ₁ from the closed forms per lane (the series below 0.5), the upward pass. Each
  lane's operations are `bessel_j`'s on its real parts, so the values are its bits. In f32 a value
  below τ beyond the turning point (n ≥ x, where jₙ has no zero and decreases) is flushed to zero
  and stays zero. Derivatives are CH4.9b's.
- `tetraleaf-math/CLAUDE.md`: a paragraph (CH4.9a2a).
- `tetraleaf-validate/examples/helmholtz_leaf.rs`: `--lanes B` (default 8); per case also the new
  harmonics (`values` and `columns` with a sum of every column, to keep it honest) and
  `bessel_j_real` over the batch, µs per point, and the projected P2M per point (CH4.9a1's
  reference P2M less its harmonics and radial parts plus the new ones); a check per case that
  `bessel_j_real` is the reference's bits (f64) and the harmonics' largest difference.
  `tetraleaf-validate/CLAUDE.md`: a sentence.
- `tools/helmholtz/ch4-9a2a-runs.sh` (`smoke`, `parts`, as CH4.9a1's), `ch4-9a2a-tables.py` on
  `common.py`.
- Report `tetraleaf-validate/results/helmholtz-ch4-9a2a.md`.
- Conventions: none (the harmonics are §3.15's; only the recurrence that forms them is new).
  No dependency. No Laplace path: `tetraleaf-math`'s Laplace modules (`harmonics`, `rotation`,
  `legendre` of §3.3) untouched; the new items are additions to the Helmholtz modules.

## Tests that prove it
- `tetraleaf-math/tests/orthonormal/leaf.rs`:
  - `values` against `harmonics` per point and term, f64, p ∈ {0, 1, 2, 5, 8, 17, 32, 45, 64,
    100}, on 1,000 random points of [−1, 1)³, the poles (0, 0, ±1), near-pole points (sin θ from
    10⁻¹ to 10⁻³⁰⁰), the axes, the equator (z = 0), the plane x = 0, the centre, the cube's
    corners; the contract |Y − Y_ref| ≤ c (n + 1) u √((2n + 1)/(4π)) (c from the measured worst
    with a margin of about 10, root `CLAUDE.md`'s glibc lesson), printed.
  - f32 against the f64 reference, p ≤ 64, the same contract with u = 2⁻²⁴ plus τ; no subnormal in
    any output (asserted).
  - lanes: B = 1, 4, 8, 16 give the same bits per point (a point's values do not depend on its
    neighbours).
  - `columns` and `values` agree bit for bit; p below the table's degree gives the first degrees
    of the table's values bit for bit; p above it panics.
- `tetraleaf-math/tests/radial/real.rs`:
  - `bessel_j_real` against `bessel_j(Complex::new(x, 0), s, …)` bit for bit in f64, lanes
    B = 1, 4, 8 with mixed x in a batch (x = 0, 10⁻¹², 10⁻⁶, 0.3, 0.5 (the series' edge), 1, 4,
    17, 60, 300), s ∈ {1, 0.3, 10⁻³}, lengths 1, 2, 9, 46, 101.
  - f32: bit for bit `bessel_j` but where the reference is below τ beyond the turning point (there
    0), and no subnormal output.

## Checks
`/verify` (the root checks, the strict workspace checks; tetraleaf-math's and tetraleaf-validate's
own); `cargo test -p tetraleaf-math --release` for the large-degree cases if they are ignored in
debug; the example's smoke run.

## Machine runs
Kathleen (one core, `singlenode`, slow nodes excluded) and locust (one core, gated on a quiet
node): `tools/helmholtz/ch4-9a2a-runs.sh parts`, p ∈ {8, 12, 15, 17, 19, 26, 32, 41, 45} × k̂ ∈
{0.5, 2, 4, 8, 16}, f64 and f32, one precision and degree per launch (B12). It decides how far the
harmonics and the radial functions fall against the sum, which CH4.9a2b's ISA question reads. The
M6 not timed (ask first; not needed).

## Constraints and risks
- Constraint 1: no Laplace path; `tetraleaf-math` gains items in its Helmholtz modules only.
  Constraint 6: each routine is tested against its reference (`harmonics`, `bessel_j`).
  Constraint 9: nothing depends on data.
- Risk: the plain recurrence's rounding near the poles grows like n² u (S1); at p ≤ 100 that is
  about 10⁴ u at worst, far below the leaf's tolerances; the test measures it.
- Risk: the flush drops values that matter. Near the pole P̄ₙᵐ/P̄ₘᵐ → √C(n + m, 2m), at most
  about 2^41 for n ≤ 64 (m = 32), so a diagonal below τ leaves values below about τ 2^43 with the
  normalisation: 2^−857 in f64 (τ = 2^−900) and 2^−57 in f32 (τ = 2^−100), far below u. The test
  checks f32 at p = 64 against the f64 reference, and the bound is stated in the docs.
- Risk: vectorisation not achieved by the compiler. The timing shows it; the lanes are fixed-size
  arrays so that LLVM sees whole vectors.

## Pending runs
None yet.

## Open questions
None yet.
