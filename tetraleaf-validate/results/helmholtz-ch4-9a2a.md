# Helmholtz CH4.9a2a: the leaf harmonics and the real radial functions, batched over points

Task: `docs/helmholtz/tasks/CH4.9a2a.md` (roadmap CH4.9a2a, the first half of CH4.9a2, split on
2026-11-27 by the skill's size rule; item 14 of the priority list). It builds the two routines that
CH4.9a1 found to take most of a Helmholtz leaf point's P2M and L2P, and measures them per point
against the reference's parts on Kathleen and locust, before CH4.9a2b builds the kernel on them.

## 0. Summary

1. **The routines** (§1): `tetraleaf_math::orthonormal::LeafHarmonics` forms the orthonormal
   harmonics of B points at once (lanes) by the plain three-term recurrence from coefficients
   tabulated once; `tetraleaf_math::radial::bessel_j_real` forms ĵₙ of B real arguments in real
   arithmetic. Both flush under-range values to zero, so that no operand is subnormal. Safe Rust,
   no `tetraleaf-simd` change, no dependency.
2. **Accuracy** (§2): the harmonics are within 1.02 (n + 1)² u √((2n + 1)/(4π)) per term of
   `harmonics` in f64 to p = 100 (0.87 in f32 to p = 64), over random points of the leaf frame, the
   poles, near-pole points (sin θ to 10⁻¹⁶⁰; below about 10⁻¹⁶² a point is taken on the axis), the
   axes, the equator and the centre; the test holds them to 10 times that. The radial functions are
   `bessel_j`'s real parts **bit for bit** in f64, and in f32 but for the values flushed beyond the
   turning point (below 2⁻¹⁰⁰).
3. **Time per point** (§3, one core, 16 lanes, medians over k̂): the harmonics take **0.06–0.09×**
   the reference's in f64 on both machines (0.03–0.06× in f32), the radial functions 0.13–0.24×.
   With the reference's sum loop kept, P2M per point is projected at **0.16–0.23× (Kathleen) and
   0.14–0.21× (locust)** of today's in f64 (0.13–0.25× over every case), and the sum is then
   20–55 % of it: the next part to go, CH4.9a2b's batched sum.
4. **f32** (§4): CH4.9a1's slowdown (f32 at p = 41–45 1.37–1.40× f64 on Kathleen) is gone from the
   routines and from the sum on their values: the projected f32 P2M there is 1.07–1.54 µs, about
   0.63× f64's at p = 45. The reference's sum on its own values stays slow in f32 (subnormal operands).
5. **Lanes** (§5): 16 lanes is the fastest of 4, 8 and 16 on both machines in both precisions but
   in Kathleen's f64 column pass, where 8 lanes are 4–5 % faster.

## 1. Setup

- **Code**: branch `tbetcke/h-ch4-9a2a` at f9bc99a, release builds, BLAS threads 1. After the runs
  (the self-review) f64's radial pass stopped flushing, to keep every bit: no f64 value of the runs
  was flushed (the column "flushed" is 0 for every f64 case), so its times stand; and the example's
  radial check gained the τ and n ≥ x conditions of the unit test.
- **The routines** (`tetraleaf-math/src/orthonormal/leaf.rs`, `radial.rs`, their module docs):
  - `LeafHarmonics::new(p)` tabulates in f64 and rounds to T once: Q̄₀⁰ = 1/√(4π), the diagonal
    factors dₘ (√3 for m = 1, √((2m + 1)/(2m)) after), √(2m + 3), and per (n, m) the
    coefficients aₙₘ = √((4n² − 1)/(n² − m²)) and bₙₘ = √((2n + 1)((n − 1)² − m²)/((2n − 3)(n² − m²)))
    of Q̄ₙᵐ = aₙₘ cos θ Q̄ₙ₋₁ᵐ − bₙₘ Q̄ₙ₋₂ᵐ, with the normalisation √((2n + 1)/(4π)) and √2 folded in.
  - `columns::<B>` forms each order's column for B points, a multiply-add chain over B lanes per
    term, and hands it with the lanes' (cos mφ, sin mφ) to a consumer; `values::<B>` stores the full
    layout of §3.6. A lane whose diagonal falls below τ is set to zero (the bound on what that
    drops: the module docs, 2⁻⁸⁵⁷ in f64, 2⁻⁵⁷ in f32).
  - `bessel_j_real::<T, B>` runs `bessel_j`'s operations on the lanes' real parts: the start order
    (restated in real f64, the same orders), the backward ratio recursion from the largest start
    with each lane entering at its own, the closed forms or the series for ĵ₀ and ĵ₁, the upward
    pass; the only departure, in f32, is the flush beyond the turning point.
- **The example** `helmholtz_leaf` (CH4.9a1's, extended; module documentation): per degree p and
  k̂ on one thread, 256 uniform points of the leaf frame with complex charges, real k; the
  reference's P2M and its parts as in CH4.9a1, then the leaf routines on the same points 16 at a
  time: the harmonics stored (`values`) and column by column (`columns`, each value added into an
  accumulator), the real radial functions, and the reference's sum loop on the routines' own values
  of 16 points cycled over the charges ("sum, leaf values"). Checked before timing: the radial
  values bit for bit (the differing ones counted, each a flushed zero), the harmonics' largest
  difference in the contract's units. "projected P2M" adds the new parts (column-wise harmonics,
  real radial functions, the factor's `exp`, the sum on the leaf values); the call's overhead is
  left out (the reference's parts add up to 0.96–1.00 of its call in f64 here).
- **The grid**: p ∈ {8, 12, 15, 17, 19, 26, 32, 41, 45} × k̂ ∈ {0.5, 2, 4, 8, 16}, f64 and f32,
  one precision and degree per launch (`tools/helmholtz/ch4-9a2a-runs.sh parts`, `LANES=16`); the
  lane count from `lanes` (p = 17 and 45 at k̂ = 4, 4, 8 and 16 lanes).
- **Kathleen**: job 673107 (`parts`), node-c11a-044, 08:06–08:16; job 673103 (`smoke`, `lanes`),
  node-c11a-043, 07:47–07:49; one core (`taskset -c 2`), 2026-11-27.
- **locust**: `parts` 08:06–08:15, `smoke` and `lanes` 07:45–07:47, one core (`taskset -c 2`),
  no other process of this account; one-minute load at most 1.08 during `parts`, and 1.6–2.9 at the
  start of `smoke` and the `lanes` launches, which followed the build (five-minute 0.67–0.74).
- **Superseded runs**, kept in the archive and not used: Kathleen job 673105 and a locust run at
  07:52–08:05 (`ch4-9a2a-v1`) timed the sum on one point's flushed radial values, whose f32 figures
  CH4.9a1's one-point sum distorted (§4); a locust run at 07:47–07:53 (`ch4-9a2a-aborted`) ran two
  copies on one core (a kill missed the script; none of its numbers is used).
- **The M6**: not timed (correctness smoke runs only).
- Raw outputs: `~/projects/fmm-archive/{kathleen,locust}/runs/ch4-9a2a/`. Tables:
  `tools/helmholtz/ch4-9a2a-tables.py {summary,cases,lanes}`.

## 2. Accuracy

`tetraleaf-math/tests/orthonormal/leaf.rs` and `tests/radial/real.rs` (in CI's `cargo test`):

| check | measure | result |
| --- | --- | --- |
| `LeafHarmonics` against `harmonics`, f64, p ∈ {0, 1, 2, 5, 8, 17, 32, 45, 64, 100}, 1,052 directions | max |Y − Y_ref| / ((n + 1)² u √((2n + 1)/(4π))) | 1.02 (bound 10) |
| the same in f32 against the f64 reference, p ≤ 64 | the same, plus 2⁻⁵⁷ | 0.87 (bound 10); no subnormal value |
| lanes B = 1, 4, 8, 16, each point in two positions | bits | the same bits |
| `columns` against `values`; a larger table at a smaller p | bits | the same bits |
| `bessel_j_real` against `bessel_j` at real z, f64, 13 arguments 0 to 300 mixed in lanes, s ∈ {1, 0.3, 10⁻³}, lengths 1–101, B = 1, 4, 8 | bits | bit for bit, every value |
| the same in f32, B = 1, 8, 16 | bits | bit for bit but 6,650 values flushed below 2⁻¹⁰⁰ at n ≥ x; no subnormal value |

The harmonics' error grows like (n + 1)² u, the plain recurrence's n² u near the poles (spike S1,
the reason `legendre` telescopes); the worst lies at the lowest degrees (at p = 100 the degree-100
terms reach about 0.6). In the timing runs every case's harmonics stayed within 1.28 (f64) and 1.77
(f32, against the f32 reference) of the same unit, and every radial value was the reference's bits
or a flushed zero.

## 3. Time per point, f64 (medians over k̂, one core, 16 lanes)

| machine | p | P2M µs | harmonics µs | leaf µs | ÷ | radial µs | real µs | ÷ | sum, leaf values µs | projected µs | ÷ P2M | sum share |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Kathleen | 8 | 1.337 | 0.476 | 0.044 | 0.092 | 0.732 | 0.159 | 0.24 | 0.081 | 0.288 | 0.23 | 28 % |
| Kathleen | 12 | 1.953 | 0.974 | 0.080 | 0.083 | 0.785 | 0.163 | 0.23 | 0.136 | 0.384 | 0.21 | 35 % |
| Kathleen | 17 | 2.908 | 1.759 | 0.140 | 0.079 | 0.844 | 0.155 | 0.20 | 0.225 | 0.524 | 0.19 | 42 % |
| Kathleen | 26 | 5.478 | 3.895 | 0.290 | 0.074 | 0.981 | 0.165 | 0.17 | 0.444 | 0.896 | 0.17 | 49 % |
| Kathleen | 32 | 7.649 | 5.787 | 0.422 | 0.073 | 1.082 | 0.178 | 0.16 | 0.647 | 1.267 | 0.16 | 52 % |
| Kathleen | 45 | 13.760 | 11.131 | 0.793 | 0.071 | 1.316 | 0.201 | 0.14 | 1.227 | 2.215 | 0.16 | 55 % |
| locust | 8 | 1.102 | 0.326 | 0.030 | 0.092 | 0.691 | 0.142 | 0.22 | 0.044 | 0.218 | 0.21 | 20 % |
| locust | 12 | 1.535 | 0.673 | 0.053 | 0.079 | 0.734 | 0.138 | 0.21 | 0.076 | 0.269 | 0.18 | 28 % |
| locust | 17 | 2.242 | 1.271 | 0.092 | 0.072 | 0.782 | 0.135 | 0.19 | 0.128 | 0.357 | 0.17 | 36 % |
| locust | 26 | 4.078 | 2.851 | 0.188 | 0.066 | 0.899 | 0.145 | 0.17 | 0.254 | 0.589 | 0.15 | 43 % |
| locust | 32 | 5.704 | 4.263 | 0.273 | 0.064 | 0.986 | 0.159 | 0.15 | 0.375 | 0.810 | 0.14 | 46 % |
| locust | 45 | 10.304 | 8.313 | 0.513 | 0.062 | 1.187 | 0.171 | 0.13 | 0.725 | 1.406 | 0.14 | 51 % |

(Every degree, both precisions: `ch4-9a2a-tables.py summary`; every case: `cases`.)

- **The harmonics** fall 11–16× in f64 (0.062–0.092× over every case on both machines): the
  column pass is now a few multiply-adds per term over 16 lanes, where the reference paid two square
  roots and four divisions per term per point (CH4.9a1 §3). At p = 45 they take 0.37 ns per term
  per point on Kathleen, 0.24 on locust. This is the column pass, timed with every value added into
  an accumulator as a streaming kernel would consume it; storing the full layout (`values`, the
  reference's output form) costs more on Kathleen, 1.06–1.16 against 0.79 µs at p = 45 (about
  0.10×), from the stores, and about the same on locust (0.48–0.54 against 0.51 µs).
- **The radial functions** fall 4–8× (0.12–0.27× over every case): real arithmetic in place of
  Smith's complex division, and the ratio recursion over 16 lanes. What remains is per lane: the
  start order (a scalar forward recursion in f64) and the closed forms' `sin` and `cos`; they now
  cost 0.13–0.20 µs per point, about flat in p.
- **The sum is next.** With the routines replaced, the reference's sum loop is 28–55 % of the
  projected P2M on Kathleen and 20–51 % on locust, rising with p; together the two routines now cost
  about the sum's time at p = 45 (0.99 against 1.23 µs on Kathleen, 0.68 against 0.73 on locust).
  CH4.9a2b's kernel batches the sum over points, so that the per-point charge weights multiply
  the column's lanes directly; the ISA question of CH4.9a2b reads this: the routines no longer stay
  far from the sum.
- k̂ matters little: at p = 45 the projected f64 P2M spans 2.20–2.27 µs over k̂ = 0.5–16 on
  Kathleen, the radial part rising with k̂ (larger start orders).

## 4. f32

| machine | p | P2M µs (reference) | ÷ f64's | leaf harmonics ÷ | sum µs (one point's reference values) | sum, leaf values µs | projected µs | ÷ P2M |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Kathleen | 17 | 2.881 | 0.99 | 0.037 | 0.217 | 0.201 | 0.411 | 0.15 |
| Kathleen | 32 | 8.043 | 1.05 | 0.034 | 0.559 | 0.431 | 0.793 | 0.10 |
| Kathleen | 41 | 16.641 | 1.42 | 0.033 | 6.642 | 0.620 | 1.110 | 0.07 |
| Kathleen | 45 | 19.289 | 1.40 | 0.033 | 6.418 | 0.797 | 1.390 | 0.07 |
| locust | 17 | 1.853 | 0.83 | 0.056 | 0.103 | 0.101 | 0.265 | 0.15 |
| locust | 45 | 7.803 | 0.76 | 0.056 | 0.421 | 0.427 | 0.933 | 0.12 |

- The reference repeats CH4.9a1's slowdown on Kathleen (f32 1.40–1.42× f64 at p = 41–45; 12.4–23.1 µs
  per point over k̂). Its sum part, timed on one point's reference values, takes up to 11.6 µs there:
  the reference's f32 harmonics and radial values hold subnormal numbers, and x86 runs subnormal
  operands slowly; locust does not show it.
- On the routines' own values the same sum loop takes 0.62–0.99 µs at p = 41–45 on Kathleen over k̂
  (0.80 µs median at p = 45, against 1.23 µs in f64), and the projected f32 P2M is 1.07–1.54 µs,
  0.63× f64's at p = 45 (medians): the slowdown is gone from the routines and from the sum on their values.
- f32's harmonics gain more than f64's (0.033–0.057×): twice the lanes per vector.
- The one-point sum part of CH4.9a1, kept in the table for comparison, is not a measure of f32 at
  these degrees (its cost depends on how many of that point's values are subnormal): a first version
  of the projection subtracted it and gave negative times (Kathleen job 673105, superseded).
- Not settled here: the products q e ĵₙ Yₙᵐ of two normal values may themselves be subnormal in f32
  (a radial value near 2⁻¹⁰⁰ times a small harmonic), and so may the multipole's high-degree
  coefficients at small k̂. Kathleen's f32 sum on the leaf values at k̂ = 8 (0.99 µs at p = 45,
  against 0.75–0.84 at the other k̂) may show it. CH4.9a2b's f32 kernel avoids or reports it (an Also
  on its line).

## 5. Lanes

| machine | precision | p | 4 lanes: columns µs | 8 | 16 | real radial µs at 4 / 8 / 16 |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| Kathleen | f64 | 17 | 0.187 | 0.133 | 0.139 | 0.215 / 0.171 / 0.156 |
| Kathleen | f64 | 45 | 1.266 | 0.782 | 0.814 | 0.294 / 0.203 / 0.184 |
| Kathleen | f32 | 45 | 1.040 | 0.624 | 0.397 | 0.262 / 0.176 / 0.129 |
| locust | f64 | 17 | 0.148 | 0.108 | 0.092 | 0.187 / 0.145 / 0.132 |
| locust | f64 | 45 | 1.047 | 0.680 | 0.514 | 0.254 / 0.172 / 0.146 |
| locust | f32 | 45 | 0.852 | 0.523 | 0.349 | 0.222 / 0.148 / 0.108 |

16 lanes is the fastest or within 5 % of it everywhere (Kathleen's f64 column pass at 8 lanes is
4–5 % faster; AVX-512 holds 8 f64 lanes, so 16 is two vectors per step), and on NEON (2 f64 lanes per
vector) the gain from 8 to 16 is 15–25 %: more independent chains hide the recurrence's latency.
CH4.9a2b takes 16 lanes in both precisions unless its own measurement says otherwise.

## 6. Question for sign-off

1. **The routines for CH4.9a2b** (recommended: adopt them). CH4.9a2b builds the P2M kernel on
   `LeafHarmonics` and `bessel_j_real`. The radial values are the reference's bits in f64; the harmonics
   differ from the telescoped reference by at most 1.02 (n + 1)² u √((2n + 1)/(4π)) per term in
   f64 (at p = 45 at most about 6·10⁻¹³ absolute at the highest degree, by that bound), the pair-level evidence. The FMM
   evidence follows in CH4.9a2b: its acceptance holds the FMM's errors within 1.01× the reference
   path's (the skill's rule: a numerical change on both measures). Alternative: the telescoped
   recurrence of `legendre` with its coefficients tabulated (no change in accuracy, slower: its
   exponent per column and the division on its loop-carried chain remain).
