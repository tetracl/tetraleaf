# Helmholtz CH4.9a1: where a leaf point's P2M and L2P time goes

Task: `docs/helmholtz/tasks/CH4.9a1.md` (roadmap CH4.9a1, the first half of CH4.9a, split on
2026-11-26: "Split a1/a2 (Recommended)"; item 14 of the priority list). It measures where the
Helmholtz leaf expansions spend their time before a kernel is designed, and puts the per-level
kind times into the harness.

## 0. Summary

1. **The harness** (§1): `helmholtz_accuracy --levels` prints P2M, M2M, L2L and L2P per level (the
   table "Translations and leaf expansions"); the `helmholtz_leaf` example times one point's P2M
   and L2P and their parts.
2. **The harmonics dominate** (§2): from p = 12 on, the orthonormal harmonics take 44–81 % of a
   point's P2M in f64 on both machines (40–79 % in f32) and, with their derivatives, 58–84 % of L2P
   with gradients; the radial functions take most of the rest (63 % at p = 8 on locust, 9–17 % at
   p = 32–45); the sum
   over the (p + 1)² terms, the part a SIMD kernel over points would speed up directly, is only
   4–10 %.
3. **Why** (§3): the harmonics routine (`tetraleaf_math::orthonormal::harmonics`) takes a square root
   and two divisions per term (n, m) and rescales every value through `powi`; the coefficients
   depend on (n, m) alone. It is the robust recurrence, telescoped near the poles, built for
   degrees in the hundreds; at the leaf degrees (p ≤ 45) they could be computed once per degree.
   The radial functions run complex Smith divisions per step of the Miller recursion, also for
   real k.
4. **f32 is barely faster** (§2): 0.93–1.08× f64's P2M time on Kathleen at p ≤ 32 and 0.77–0.90× on
   locust;
   the cost is not arithmetic throughput. On Kathleen f32 at p = 41–45 is slower than f64 (subnormal
   values at small k̂, §2).
5. **Recommendation (question 1)**: CH4.9a2 as a restructured kernel in safe Rust, not an ISA kernel
   first: leaf harmonics (and their derivatives) from coefficients precomputed per degree, the
   radial functions in real arithmetic for real k, the points of a leaf batched so that the sum
   can vectorise; `tetraleaf-simd` untouched (constraint 1). ISA kernels only if that falls short.
   Since the harmonics also dominate L2P, the shared routine is built once (question 2).

## 1. Setup

- **Code**: branch `tbetcke/h-ch4-9a1` at 5ed8c44, release builds, BLAS threads 1.
- **The example** `helmholtz_leaf` (module documentation): per degree p and k̂ = Re(k) r on one
  thread, 256 uniform random points of the leaf frame `Densities` uses (half-width 1, [−1, 1)³)
  with complex charges, real k: the reference `leaf::p2m` and `leaf::l2p` (without and with
  gradients) over the batch, and their parts timed apart by calling the same public functions in
  the same order; the parts reassemble the reference's multipole **bit for bit in all 180 cases**
  (90 per machine), checked before any timing. µs per point, the median over 15 batches of at
  least 20 ms (`bench::median_time_per_call`).
- **The grid**: p ∈ {8, 12, 15, 17, 19, 26, 32, 41, 45} (the leaf levels' degrees of CH4.6a2's
  cases under the defaults) × k̂ ∈ {0.5, 2, 4, 8, 16} (W1–W3's finest levels have k̂ 3.9–7.8),
  f64 and f32; one precision and degree per launch (`tools/helmholtz/ch4-9a1-runs.sh parts`).
- **Kathleen**: job 673085, node-c11a-142, one core (`taskset -c 2`), 2026-11-26 20:45–20:51.
- **locust**: 2026-11-26 20:45–20:52, one core (`taskset -c 2`), started by the quiet gate (no
  other process of this account above 20 % of a core, one-minute load below 2); load at most 2.7
  during the runs.
- **The M6**: not timed (the user's word for CH4.6a2, "Skip the M6"; the example's smoke ran there
  for correctness only).
- Raw outputs: `~/projects/fmm-archive/{kathleen,locust}/runs/ch4-9a1/`. Tables:
  `tools/helmholtz/ch4-9a1-tables.py {shares,summary}`.

## 2. The parts, f64 (median over k̂)

Shares of the reference call's time per point; "L2P ∇" is L2P with gradients, its parts the radial
functions with derivatives, the harmonics with derivatives, and the rest (the sum, the gradient's
assembly, the call).

| machine | p | P2M µs | ns per term | radial | harmonics | sum | L2P ∇ µs | ns per term | radial ′ | harmonics ′ | rest |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Kathleen | 8 | 1.363 | 16.8 | 54 % | 36 % | 6 % | 2.137 | 26.4 | 36 % | 50 % | 15 % |
| Kathleen | 12 | 1.960 | 11.6 | 40 % | 50 % | 7 % | 3.582 | 21.2 | 23 % | 61 % | 17 % |
| Kathleen | 17 | 2.975 | 9.2 | 28 % | 62 % | 8 % | 6.022 | 18.6 | 15 % | 68 % | 17 % |
| Kathleen | 26 | 5.435 | 7.5 | 18 % | 73 % | 9 % | 12.371 | 17.0 | 9 % | 73 % | 19 % |
| Kathleen | 32 | 7.650 | 7.0 | 14 % | 76 % | 9 % | 18.322 | 16.8 | 6 % | 73 % | 20 % |
| Kathleen | 45 | 13.763 | 6.5 | 9 % | 81 % | 10 % | 33.885 | 16.0 | 4 % | 77 % | 19 % |
| locust | 8 | 1.101 | 13.6 | 63 % | 30 % | 4 % | 1.604 | 19.8 | 45 % | 45 % | 11 % |
| locust | 12 | 1.533 | 9.1 | 48 % | 44 % | 5 % | 2.520 | 14.9 | 31 % | 58 % | 11 % |
| locust | 17 | 2.234 | 6.9 | 35 % | 57 % | 6 % | 4.144 | 12.8 | 20 % | 69 % | 11 % |
| locust | 26 | 4.062 | 5.6 | 22 % | 70 % | 6 % | 8.397 | 11.5 | 11 % | 78 % | 11 % |
| locust | 32 | 5.689 | 5.2 | 17 % | 75 % | 7 % | 12.160 | 11.2 | 9 % | 80 % | 11 % |
| locust | 45 | 10.297 | 4.9 | 12 % | 81 % | 7 % | 22.895 | 10.8 | 6 % | 84 % | 11 % |

(Every degree of the grid, both precisions and every k̂: `ch4-9a1-tables.py summary` and `shares`.)

- **The cost per point matches the FMM's.** Kathleen's P2M at p = 17, 2.98 µs per point on one core,
  is the in-FMM thread time per point measured on W1 at kD = 1000's finest level (about 3.3 µs,
  CH4.6a2's per-level table, job 673083); the leaf operators run at one core's per-point cost.
- **k̂ matters little**: at p = 17 P2M takes 2.79–3.36 µs over k̂ = 0.5–16 on Kathleen, the radial
  share rising with k̂ (25–37 %).
- **The parts add up**: their sum is 0.96–0.99 of the reference call on locust and 0.95–1.05 on
  Kathleen, but for Kathleen's f32 at p = 41–45 (0.89–1.52, below).
- **f32**: on Kathleen 0.93–1.08× f64's P2M time at p ≤ 32, on locust 0.77–0.90×; L2P ∇ 0.96–1.09×
  on Kathleen, 0.82–0.86× on locust. At p = 41 and 45 on Kathleen f32's P2M takes 1.37–1.40× f64's
  and its sum part up to 11.6 µs per point: at small k̂ the high-degree radial values fall below
  f32's normal range, and x86 runs subnormal arithmetic slowly (locust does not show it). Relevant
  to f32 under D3 (CH7.1), not to the f64 kernel.

## 3. Why the harmonics cost what they do

`tetraleaf_math::orthonormal::harmonics` (CH1.2a, spike S1) runs, per column m, the normalised
Legendre recurrence in its telescoped form (the differences E_n, stable near the poles at degrees
in the hundreds), with a mantissa and binary exponent per column. Per term (n, m) it takes
c = √((n − m)(n + m)), two divisions (m²/(n + c), m²/(n − 1 + c′)) and a third (diff/c), plus the
normalisation √((2n + 1)/4π) per degree and, wherever the exponent is not zero, two `powi` calls in
`scaled`. None of these depends on the point. At the leaf degrees (p ≤ 45, far from the range the
telescoped form was built for) the coefficients can be tabulated once per degree, leaving a few
multiply-adds per term.

`tetraleaf_math::radial::bessel_j` runs Miller's backward recursion from a start order found by a
forward recursion until a growth factor (`start_order_j`), with Smith's complex division per step,
also for real z. For real k (D1's range, every case here) the recursion can run in real arithmetic,
and the start order can come from a closed-form estimate.

## 4. Question for sign-off: CH4.9a2's design

1. **The kernel** (recommended: restructured safe Rust first).
   - Leaf harmonics (and, for L2P, their angular derivatives) from coefficients tabulated once per
     degree, a new routine beside `tetraleaf_math::orthonormal::harmonics` (or in
     `tetraleaf-tables::helmholtz`), tested against it per point within a contract over the leaf
     frame's directions, poles included.
   - The radial functions for real k in real arithmetic (complex k keeps the reference's path).
   - A leaf's points batched, the sum written so that it vectorises over points without `unsafe`.
   - No `tetraleaf-simd` change (constraint 1 untouched); an ISA kernel only if the measured time
     stays far from the sum's.
   - If the harmonics fell to about the sum's cost, P2M would take roughly 0.3–0.5× today's at
     p ≥ 17 (an estimate from §2's shares, to be measured).
   - Alternatives: ISA kernels in `tetraleaf-simd` now (constraint 1 and D9 apply, and the per-term
     setup the harmonics cost would remain unless rewritten too); or both at once.
2. **The order**: the leaf harmonics routine is shared by P2M and L2P (it is 58–84 % of L2P with
   gradients), so CH4.9a2 builds it with P2M, and CH4.9b reuses it for L2P (recommended, the
   roadmap's order); or one task for both.
