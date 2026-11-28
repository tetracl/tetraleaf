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

## As built
- The contract of the derivatives: the first guess, C (n + 1)³ u √((2n + 1)/(4π)) with C = 10, held.
  The worst measured on the M6 is 0.87 (∂θ) and 0.59 ((1/sin θ) ∂φ) in those units in f64 to
  p = 100, and 0.76 and 0.26 in f32 to p = 64. The self-review's probe, per degree at p = 100, put ∂θ
  at 0.61–0.80 at every degree, so the third power is the growth, not a low-degree effect.
- Not in "Changes" above: `LeafHarmonics`' directions (`Lanes::of`, shared with `columns` and so with
  CH4.9a2b's P2M kernel) scale x and y, and z, by 2⁶⁰⁰ (f64) or 2¹⁰⁰ (f32) where ρ² or r² is below
  the normal range. The first test run failed at (10⁻¹⁶⁰, 0, 1) in f64 and (−10⁻⁴⁰, 10⁻⁴⁰, −0.5) in
  f32: cos φ and sin φ from a subnormal ρ² lost precision. The values never showed it, since they
  carry sinᵐ θ, but the derivatives of order 1 do not. The scale is a select, 1 elsewhere, so every
  other point keeps its bits. A point whose ρ² underflowed used to be taken on the z-axis; now only
  x = y = 0 is. In the self-review's fix, ρ is scaled back by an exact power of two, not by a
  division, so sin θ keeps its precision where both scales apply, and the test directions reach both
  branches.
- `LeafScratch` holds six columns (the values, order 1's, U and order 1's U, ∂θ, and a zero column
  that order 0 hands out as u).
- `bessel_j_real` and `bessel_j_real_derivatives` share a private routine (`real_j`); the values'
  path is unchanged, plus the ratio at the requested length, which the derivative of the last order
  reads.
- The example checks, before timing, that its restated L2P sum with gradients on the reference's own
  values is `leaf::l2p` with gradients bit for bit at 16 points (f64, f32, real and complex k). It also
  asserts the derivatives' contract (C = 10).

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
None.
- **Kathleen:** job 673246 (`run.sbatch`, `singlenode`, `--exclude=node-c11b-004,node-c11b-013`) ran on
  node-c11a-016, 2026-11-27 23:51 – 2026-11-28 00:05, after 1 h 55 min in the queue.
  - Step 1 was the release build of f3d24bc, synced while the job was pending.
  - Step 2 was `ch4-9b1-runs.sh smoke parts`; every launch exited 0.
- **locust:** on f3d24bc 22:28–22:41 (`parts`, every launch exit 0; the one-minute load at most 1.21
  at a launch's start, from the build), and the rescale's ABBA 22:47–22:58.
- **Superseded:** locust 21:57–22:10 on ae26e8f (`ch4-9b1-v1`, before the self-review's fix; within
  about 1 % per row).
- Outputs in `~/projects/fmm-archive/{kathleen,locust}/runs/ch4-9b1` (locust also `-v1` and `-ab`).

## Results (the report `tetraleaf-validate/results/helmholtz-ch4-9b1.md`)
- **Accuracy:** the derivatives are within 0.87 and 0.59 (n + 1)³ u √((2n + 1)/(4π)) in f64, 0.76 and
  0.26 in f32; the radial derivatives are `bessel_j`'s bits in f64.
- **L2P with gradients projected per point** (medians over k̂):
  - locust: 0.17–0.23× (f64), 0.18–0.22× (f32);
  - Kathleen: 0.24–0.27× (f64), 0.21–0.25× (f32).
- **The derivative columns:** 0.064–0.100× the reference's harmonics with derivatives (f64, both
  machines).
- **The sum's share of the projection:** 37–62 % on locust, 49–73 % on Kathleen. It costs 1.12–1.57
  ns per term on locust and 2.8–3.4 on Kathleen, in either precision there.
- **The rescale:** P2M's kernel within 1 % (ABBA, locust).

## Checks run (the M6)
`/verify` on f3d24bc:
- fmt.
- Clippy: default, examples, `--workspace --all-targets`.
- `cargo doc` with `-D warnings`.
- At opt-level 1: `cargo test` (1,429 passed) and `cargo test --workspace` (1,473 passed).
- The MPI test `helmholtz_fmm` at 2 and 4 ranks (54.3 s and 49.1 s), since the batched P2M reads
  `LeafHarmonics`' changed directions.
- The debug time of `tetraleaf-math`'s `orthonormal` at opt-level 0: 22.7 s, against `main`'s 22.3 s
  in a worktree of its own at the same load, about 6. A first reading of 29.2 s was taken under a
  load of 20–54.
- The example's checks on the M6, short smoke runs with their times not used: the restated gradient
  sum bit for bit `leaf::l2p` (f64, f32, real and complex k) and the contracts.

The same checks ran on ae26e8f before the self-review's fixes, with the same results.

Not run:
- the rest of the MPI list (no change to exec, plan or octree);
- `cargo deny` (no manifest change);
- Laplace's benchmarks (no Laplace path; `tetraleaf-simd` untouched).

A self-review by a subagent found nothing blocking. Fixed:
- `Lanes::of`'s ρ through a subnormal where both scales apply;
- the example's gradient sum checked against `leaf::l2p`;
- the derivatives' contract asserted in the example;
- the record of the rescale;
- the test nits: the signs of zero in the bit test, the f32 radial flush reached, degree 1's zeros.

## Open questions
None.
