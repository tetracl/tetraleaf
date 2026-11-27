# Helmholtz CH4.9b1: the leaf harmonics' derivatives and the real radial derivatives, batched over points

Task: `docs/helmholtz/tasks/CH4.9b1.md` (roadmap CH4.9b1, the first half of CH4.9b, split on
2026-11-27 by the skill's size rule; item 14 of the priority list). It builds the derivatives that a
Helmholtz leaf point's L2P with gradients needs, on CH4.9a2a's leaf routines, and measures them per
point against the reference's parts, before CH4.9b2 builds the L2P kernel on them.

## 0. Summary

1. **The routines** (§1): `LeafHarmonics::derivative_columns` forms, per order, the column Q̄ₙᵐ
   (`columns`' bits), ∂θ Q̄ₙᵐ and Uₙᵐ = Q̄ₙᵐ / sin θ for B points as lanes. U comes from a second
   column of the same recurrence, so nothing is divided by sin θ. `bessel_j_real_derivatives` adds
   ĵₙ′ to `bessel_j_real`. Safe Rust, no `tetraleaf-simd` change, no dependency.
2. **Accuracy** (§2):
   - The derivatives are within 0.87 (∂θ) and 0.59 ((1/sin θ) ∂φ) (n + 1)³ u √((2n + 1)/(4π)) per
     term of `harmonics_derivatives` in f64 to p = 100, and 0.76 and 0.26 in f32 to p = 64. The test
     holds them to 10 times that. The points cover the leaf frame, the poles, near-pole points, the
     axes, the equator, the centre, and coordinates whose squares are subnormal.
   - The radial derivatives are `bessel_j`'s bit for bit in f64. In f32 they differ only beside the
     values flushed below 2⁻¹⁰⁰ beyond the turning point.
3. **Time per point on locust** (§3; one core, 16 lanes):
   - The derivative columns take **0.064–0.100×** the reference's harmonics with derivatives in f64
     (0.056–0.077× in f32), and the real radial functions with derivatives 0.15–0.23× the
     reference's (medians over k̂; 0.13–0.25× per case).
   - With the reference's sum kept (restated, on the routines' values), L2P with gradients per point
     is projected at **0.17–0.24×** of today's in f64 (0.18–0.24× in f32).
   - The sum is then 25–63 % of the projection in f64 (28–70 % in f32), 37–62 % as medians over k̂,
     growing with p: the part CH4.9b2's kernel vectorises.
4. **Kathleen** (§3): pending (job 673246, queued behind a full partition).
5. **A fix to CH4.9a2a's directions** (§4): next to the z-axis, cos φ and sin φ came from a
   subnormal ρ² and lost precision. The values never showed it, but the derivatives of order 1 do.
   The coordinates are now scaled by a power of two there, through a select, so every other point
   keeps its bits. P2M's kernel time is unchanged within 1 % (an ABBA comparison on locust).

## 1. Setup

- **Code**: branch `tbetcke/h-ch4-9b1` at f3d24bc, release builds, BLAS threads 1.
- **The routines** (`tetraleaf-math/src/orthonormal/leaf.rs`, "Derivatives"; `radial.rs`):
  - U obeys the columns' recurrence, from Uₘᵐ = dₘ Q̄ₘ₋₁ᵐ⁻¹, the diagonal before its factor sin θ.
  - From sin θ ∂θ Pₙᵐ = n cos θ Pₙᵐ − (n + m) Pₙ₋₁ᵐ in the columns' normalisation:
    ∂θ Q̄ₙᵐ = n cos θ Uₙᵐ − cₙₘ Uₙ₋₁ᵐ, cₙₘ = √((2n + 1)(n² − m²)/(2n − 1)), for m ≥ 1.
  - For m = 0, ∂θ Q̄ₙ⁰ = −√(n(n + 1)/2) Q̄ₙ¹, so order 1 is formed before order 0 is handed out.
  - Then ∂θ Yₙ^{±m} = ∂θ Q̄ₙᵐ {cos mφ, sin mφ} and (1/sin θ) ∂φ Yₙ^{±m} = m Uₙᵐ {−sin mφ, cos mφ}:
    `harmonics_derivatives`' outputs, in its frame on the z-axis and at the centre.
  - U's diagonal is flushed below τ on its own value. At a pole Q̄₁¹ = 0, but U₁¹ carries ∂θ Y₁¹.
  - `bessel_j_real_derivatives` runs `bessel_j_real`'s code (one private routine, its values'
    bits) and `bessel_j`'s derivative formula from the neighbouring orders on the real parts. The
    order above the last comes from the ratio the backward recursion passes there.
- **The example** `helmholtz_leaf` (extended; its module documentation): per degree p and k̂ on one
  thread, 256 uniform points of the leaf frame, real k, f64 or f32 per launch.
  - The reference's L2P with gradients and its parts as in CH4.9a1: `bessel_j` with derivatives,
    and `harmonics_derivatives`.
  - The derivative columns over the same points, 16 at a time, every value summed into an
    accumulator.
  - The real radial functions with derivatives.
  - "sum ∇, leaf values": the reference's L2P sum with gradients (its loop over the terms and the
    gradient's Cartesian assembly, restated) on the routines' values of 16 points, cycled over the
    targets.
  - Checked before timing:
    - the restated sum on the reference's own values is `leaf::l2p` with gradients bit for bit;
    - the derivatives are within their contract (asserted);
    - the radial derivatives are bit for bit, or within (n + 1) τ / s beside a flushed value.
  - "projected L2P ∇" adds the derivative columns, the radial functions with derivatives, the
    factor's `exp` and that sum. It leaves out the call's overhead, as CH4.9a2a's projection did.
- **Runs** (`tools/helmholtz/ch4-9b1-runs.sh parts`, one precision and degree per launch, review
  2026-11-24 B12): p ∈ {8, 12, 15, 17, 19, 26, 32, 41, 45} (the leaf levels' degrees of CH4.6a2's
  cases) × k̂ ∈ {0.5, 2, 4, 8, 16}, f64 then f32, pinned to core 2.
  - locust: 2026-11-27 22:28–22:41. The one-minute load was at most 1.21 at a launch's start, from
    the build before it; no other process of this account was busy. Every launch exited 0.
  - A first locust run on ae26e8f (21:57–22:10, before the self-review's fix to `Lanes::of`) agrees
    within about 1 % per row. It is superseded and kept in `ch4-9b1-v1`.
  - Kathleen: pending.
  - Outputs: `~/projects/fmm-archive/{locust,kathleen}/runs/ch4-9b1/`.
  - Tables: `python3 tools/helmholtz/ch4-9b1-tables.py summary|cases <machine>=<dir>`.

## 2. Accuracy

`tetraleaf-math/tests/orthonormal/leaf.rs`, on the M6, against `harmonics_derivatives` per point and
term over CH4.9a2a's directions:

- The directions: 1,000 random points of [−1, 1)³, the poles, near-pole points with sin θ from 10⁻¹
  to 10⁻³⁰⁰, the axes, the equator, the plane x = 0, the centre and the cube's corners.
- This task adds (3, −4, 5)·10⁻²⁰⁰ and (3·10⁻¹⁷⁰, 4·10⁻¹⁷⁰, 10⁻²⁰⁰), whose ρ² and r² lie below the
  normal range.

The worst error in units of (n + 1)³ u √((2n + 1)/(4π)) (the values in (n + 1)² units, as CH4.9a2a):

| precision | p | values | ∂θ | (1/sin θ) ∂φ |
| --- | ---: | ---: | ---: | ---: |
| f64 | 1 | 1.023 | 0.640 | 0.512 |
| f64 | 8 | 1.023 | 0.713 | 0.587 |
| f64 | 17 | 1.023 | 0.811 | 0.587 |
| f64 | 45 | 1.023 | 0.846 | 0.587 |
| f64 | 100 | 1.023 | 0.872 | 0.587 |
| f32 | 8 | 0.872 | 0.760 | 0.259 |
| f32 | 64 | 0.872 | 0.760 | 0.259 |

- The power of n: the self-review's probe put ∂θ's worst at 0.61–0.80 of the unit at every degree
  from 1 to 100 at p = 100 (1.3 to 73 in (n + 1)² units). So (n + 1)³ is the growth: both terms of
  ∂θ are about n |U|, and U carries the recurrence's n² u.
- The test holds C = 10.
- In the example's runs the worst "derivatives check" over every case is 0.66 on locust. The example
  compares against `harmonics_derivatives` in the run's precision, so in f32 the measure includes
  the reference's own rounding.
- The radial derivatives (`tests/radial/real.rs`) are `bessel_j`'s derivatives bit for bit in f64
  (B = 1, 4, 8; x from 0 to 300 mixed in a batch; s = 1, 0.3, 10⁻³; lengths 1–101), and so are the
  values against `bessel_j_real`.
  - In f32, 7,975 derivatives of the test's grid differ, each beside a flushed value and within
    (n + 1) τ / s, and none is subnormal.
  - In the runs, f32 differs at 5–9,303 derivatives per degree from p = 15 (the column "radial ′
    differ"), every one so checked, and f64 at none.

## 3. Time per point

µs per point, one core, 16 lanes, the median over k̂ ∈ {0.5, 2, 4, 8, 16}. "harmonics ′" is the
reference's `harmonics_derivatives`, "leaf ′" the derivative columns. "radial ′" is `bessel_j` with
derivatives, "real ′" `bessel_j_real_derivatives`. "sum share" is the restated sum's share of the
projection.

### locust (Neoverse-V2, NEON)

| precision | p | L2P ∇ | harmonics ′ | leaf ′ | ÷ | radial ′ | real ′ | ÷ | sum ∇, leaf values | projected | ÷ L2P ∇ | sum share |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| f64 | 8 | 1.602 | 0.713 | 0.071 | 0.100 | 0.713 | 0.147 | 0.23 | 0.127 | 0.346 | 0.23 | 37 % |
| f64 | 12 | 2.519 | 1.470 | 0.129 | 0.088 | 0.768 | 0.148 | 0.21 | 0.233 | 0.511 | 0.21 | 46 % |
| f64 | 15 | 3.425 | 2.242 | 0.180 | 0.080 | 0.798 | 0.145 | 0.20 | 0.333 | 0.659 | 0.20 | 51 % |
| f64 | 17 | 4.139 | 2.849 | 0.221 | 0.078 | 0.825 | 0.150 | 0.20 | 0.410 | 0.782 | 0.19 | 52 % |
| f64 | 19 | 4.931 | 3.526 | 0.266 | 0.075 | 0.856 | 0.149 | 0.19 | 0.495 | 0.911 | 0.19 | 54 % |
| f64 | 26 | 8.388 | 6.516 | 0.455 | 0.070 | 0.959 | 0.168 | 0.17 | 0.859 | 1.484 | 0.18 | 58 % |
| f64 | 32 | 12.123 | 9.755 | 0.660 | 0.068 | 1.059 | 0.185 | 0.17 | 1.265 | 2.104 | 0.17 | 60 % |
| f64 | 41 | 19.169 | 15.892 | 1.032 | 0.065 | 1.221 | 0.203 | 0.15 | 1.997 | 3.228 | 0.17 | 62 % |
| f64 | 45 | 22.835 | 19.118 | 1.228 | 0.064 | 1.288 | 0.206 | 0.15 | 2.379 | 3.815 | 0.17 | 62 % |
| f32 | 8 | 1.390 | 0.533 | 0.040 | 0.075 | 0.683 | 0.126 | 0.20 | 0.125 | 0.293 | 0.22 | 43 % |
| f32 | 12 | 2.152 | 1.144 | 0.077 | 0.067 | 0.731 | 0.125 | 0.18 | 0.231 | 0.434 | 0.21 | 53 % |
| f32 | 15 | 2.888 | 1.748 | 0.113 | 0.065 | 0.758 | 0.121 | 0.18 | 0.330 | 0.567 | 0.20 | 58 % |
| f32 | 17 | 3.468 | 2.222 | 0.141 | 0.063 | 0.784 | 0.122 | 0.17 | 0.406 | 0.672 | 0.20 | 60 % |
| f32 | 19 | 4.107 | 2.743 | 0.173 | 0.063 | 0.817 | 0.126 | 0.16 | 0.489 | 0.791 | 0.19 | 62 % |
| f32 | 26 | 6.886 | 5.054 | 0.298 | 0.059 | 0.910 | 0.127 | 0.15 | 0.848 | 1.274 | 0.19 | 67 % |
| f32 | 32 | 9.914 | 7.593 | 0.449 | 0.059 | 1.005 | 0.132 | 0.14 | 1.234 | 1.819 | 0.18 | 68 % |
| f32 | 41 | 15.725 | 12.537 | 0.707 | 0.056 | 1.157 | 0.146 | 0.13 | 1.947 | 2.809 | 0.18 | 69 % |
| f32 | 45 | 18.725 | 15.093 | 0.841 | 0.056 | 1.221 | 0.155 | 0.13 | 2.317 | 3.319 | 0.18 | 70 % |

- Over every case the projection is 0.17–0.24× (f64) and 0.18–0.24× (f32) of the reference's L2P
  with gradients.
- The derivative columns cost about 2.2–2.4× CH4.9a2b's plain columns (0.513 → 1.228 µs at p = 45 in
  f64). That is three columns' work (Q̄, U, ∂θ) for one.
- The radial functions with derivatives cost 0.145–0.206 µs in f64 (medians over k̂), 1.05–1.25×
  `bessel_j_real` alone in the same launches (0.131–0.166 µs); they grow with k̂ through the start
  order.
- What is left is the sum: per target and term three complex products (the value and two
  derivatives) on each plane. It costs 1.12–1.57 ns per term in f64 (from p = 45 down to p = 8),
  against 0.27–0.40 ns per term for the batched P2M kernel's sum ("kernel sum") in the same launches.
  - The sum is restated scalar code. In CH4.9b2's kernel it runs per target over terms, with no
    lane reduction (the roadmap's ISA question).

### Kathleen (Cascade Lake, SSE2 baseline build)

Pending: job 673246 (`singlenode`, one core, slow nodes excluded; the same stages on f3d24bc).
Submitted 2026-11-27 21:55; at 22:50 still waiting on a full partition (189 nodes allocated).

## 4. The directions next to the z-axis, and P2M

The first test run failed at (10⁻¹⁶⁰, 0, 1) in f64 and at (−10⁻⁴⁰, 10⁻⁴⁰, −0.5) in f32:

- `Lanes::of` (CH4.9a2a) forms cos φ = x/ρ with ρ = √(x² + y²). There ρ² was subnormal, so cos φ
  and sin φ kept only a few digits (5·10⁻⁶ relative in f64; in f32 ρ² was zero, so φ = 0).
- The values Yₙ^{±m} carry sinᵐ θ, about 10⁻¹⁶⁰ there, so they never showed it.
- ∂θ Y₁^{±1} and (1/sin θ) ∂φ Y₁^{±1} carry cos φ and sin φ undamped.

The fix scales x and y by 2⁶⁰⁰ (f64) or 2¹⁰⁰ (f32) where ρ² is below the normal range, and x, y and z
where r² is:

- The scale is exact, a power of two, chosen by a select.
- ρ in r's units is ρ's own where the scales agree, else scaled back by the exact 2⁻⁶⁰⁰. That was the
  self-review's fix: the first version divided and passed through a subnormal at points within
  10⁻³⁰⁸ of the centre.
- Elsewhere the scale is 1, so every point with normal squares, which is every point a leaf holds
  away from 10⁻¹⁵⁴ of an axis, keeps CH4.9a2a's bits. P2M's kernel, which reads `columns`,
  therefore changes no bit of any FMM test (`helmholtz_fmm` at 2 and 4 ranks unchanged, and the
  P2M kernel's tests pass).

The time on locust, P2M's column pass ("leaf columns", µs per point) and the P2M kernel, medians
over k̂, against CH4.9a2b's locust launches:

| precision | p | columns, CH4.9a2b | columns, here | ÷ | P2M kernel, CH4.9a2b | here | ÷ |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| f64 | 17 | 0.091 | 0.095 | 1.044 | 0.322 | 0.333 | 1.034 |
| f64 | 45 | 0.513 | 0.542 | 1.057 | 1.284 | 1.299 | 1.012 |
| f32 | 17 | 0.052 | 0.056 | 1.077 | 0.245 | 0.253 | 1.033 |
| f32 | 45 | 0.342 | 0.353 | 1.032 | 0.969 | 0.996 | 1.028 |

The kernel ran 1.2–4.1 % slower than in CH4.9a2b's launches, and the column pass 3–14 %. To separate
the rescale from the binary's code placement (pure code moves shift host times by ±5–13 %), the same
commit was built twice on locust: once as it is and once with CH4.9a2b's `Lanes::of` restored, the
rescale alone removed.
- Each was launched four times in ABBA order (without, with, with, without) at p = 17 and 45, in f64
  and f32, k̂ ∈ {0.5, 2, 4, 8, 16}, one setting per launch, on core 2.
- 2026-11-27 22:47–22:58, the one-minute load at most 1.25 at a launch's start.
- Outputs in `~/projects/fmm-archive/locust/runs/ch4-9b1-ab/`.

The means of each pair of launches, µs per point, medians over k̂:

| precision | p | P2M kernel, without | with | ÷ | columns, without | with | ÷ | derivative columns ÷ |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| f64 | 17 | 0.331 | 0.331 | 1.000 | 0.091 | 0.095 | 1.044 | 0.998 |
| f64 | 45 | 1.301 | 1.300 | 1.000 | 0.515 | 0.542 | 1.053 | 1.000 |
| f32 | 17 | 0.251 | 0.253 | 1.008 | 0.055 | 0.056 | 1.018 | 1.004 |
| f32 | 45 | 0.984 | 0.992 | 1.008 | 0.352 | 0.353 | 1.003 | 1.001 |

- The rescale leaves P2M's kernel within 1 %: 1.000× in f64 and 1.008× in f32, repeated in both
  pairs.
- The kernel's 1–4 % against CH4.9a2b's launches is the binary: the build without the rescale ran
  1.301 µs at p = 45 against CH4.9a2b's 1.284.
- The example's standalone column pass is 4–5 % slower in f64 with the rescale, consistently in both
  pairs. The pass is O(p²) per point and the rescale a few operations per point, so this reads as
  code generation in the example's own loop, not the arithmetic; the kernel, which reads the same
  columns, shows none.

## 5. Questions

None for sign-off. CH4.9b2 builds the L2P kernel on these routines; its FMM errors are their FMM
evidence, as CH4.9a2b's were for CH4.9a2a's.
