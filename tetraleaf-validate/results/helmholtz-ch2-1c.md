# Helmholtz CH2.1c: Yukawa's batched P2M

Task: docs/helmholtz/tasks/CH2.1c.md (roadmap CH2.1c, item 15 of the priority list; the third part
of CH2.1, proposal §10 CH2.1, §5.5). The kernel is `tetraleaf_tables::helmholtz::leaf::LeafP2m` on
Yukawa's real path (its module docs, "Yukawa"), run by `Densities::p2m_kernel` in the FMM. Raw
outputs are in `~/projects/fmm-archive/kathleen/runs/ch2-1c/` (`acc`, `time`); the tables come from
`tools/helmholtz/ch2-1c-tables.py`, the launches from `tools/helmholtz/ch2-1c-runs.sh`.

## 0. Summary

- **The kernel.** `LeafP2m` now restates `p2m_yukawa` (k = iκ): îₙ by `modified_i` per lane at the
  reference's length, e^{κρ − β} per lane, real weights per real charge column, then CH4.9a2b's
  columns and lane sums. Every coefficient is within the contract
  C u √((2n + 1)/(4π)) Σⱼ |qⱼ eⱼ| (((n + 1)² + N) îₙ + κρⱼ (|îₙ′| + îₙ)), C = 10: the worst measured
  ratio is 0.80 of C in f64 (p ≤ 45, κ̂ ≤ 250, i.e. κ̂ a ≤ 433) and 1.54 of C in f32 (§2).
- **The FMM's errors within 1.01× the per-point path's, and within ε.** On CH2.1b's 114 cases (W1,
  W3 at 10⁵ and 10⁶, W4's four distributions; κD = 10⁻² … 10³; ε = 10⁻³, 10⁻⁶, 10⁻⁸) the batched
  P2M's errors are the reference P2M's to the harness's three printed digits (1.0000×) in every
  case, the largest 0.27 ε as in CH2.1b; also with complex charges and on the constant density (§3).
  The tests: 1.0000× and outputs within 4.1·10⁻¹⁶ of the reference leaves' on Yukawa's three
  scenarios.
- **Deterministic.** A leaf's multipole depends on its chunk alone. The FMM tests' threads, pair
  calls, vectors, complex charges as two real vectors and the one-rank FMM at 2 and 4 ranks hold bit
  for bit with the kernel.
- **Time.** P2M per run 0.29–0.34× the per-point path's on Kathleen (1 × 40), from 7–10 % of an
  evaluation to 2–3 %; the evaluation 0.93–0.95× (§4). Yukawa's L2P, still the per-point
  reference's, is now 15–25 % of an evaluation: the question of §5.

## 1. Setup

- **Machine.** Kathleen, one exclusive node (2 × 20 Cascade Lake cores), `singlenode`, slow nodes
  excluded, both jobs on node-c11a-011: job 673284 (build, `smoke`, `accuracy`; its `worst` and W1
  `levels` launches held no point, below) and job 673285 (`worst` and `levels` again, on W1n). One
  rank × 40 threads, BLAS threads 1, release builds of 4893d27 (the library as d08db08).
- **Settings.** The library's defaults named (`--m2l auto-pw --diagonal auto --diagonal-build
  interpolated`, the V and W/X rules of the cost model; on Yukawa Auto runs Rotation on every level),
  gradients, §9.1's samples, real charges unless named. One P2M per launch: `--leaves kernel` (the
  library's default) or `--leaves reference` (the per-point reference calls; P2P the AVX-512 kernel
  in both). Yukawa is `--kd 0 --kim κ` with D = 1.
- **Workloads.** W1 and W2 size their points by Re k, so on Yukawa their defaults hold no point (the
  first job's `worst` and W1 `levels` launches printed "no rank has a source or a target point").
  W1 runs as W1n (`w1:<n>`, CH2.1b) at a fixed N; W2 has no fixed-N form (its cap is refined by
  wavelength), so the constant density runs on W1 and W3 at 10⁶ points instead.

## 2. The contract per leaf (the M6)

`tetraleaf-tables/tests/helmholtz/leaf_yukawa.rs`, against `tetraleaf_ref::helmholtz::leaf::p2m_yukawa`
per coefficient, in the FMM's leaf frame (half-width 1):

- f64: p ∈ {0, 1, 2, 5, 8, 17, 26, 45}, κ̂ ∈ {10⁻⁶, 0.3, 4, 16, 60, 250} (κ̂ = 250, κ̂ a = 433, a
  shallow leaf at κD = 10³), leaves of N = 1, 5, 16, 17, 40, 200 points (for even seeds the centre,
  the poles, points within 10⁻⁹ of the axis, the axes and a corner first; then uniform points), real
  and complex charges, one to three vectors: the worst error is 0.80 of the contract without C. f64's
  flush (2⁻⁹⁰⁰) is not reached by these points.
- f32 against the f64 reference on the f32 inputs, p ≤ 32, κ̂ ≤ 60: 1.54. Two terms beside the
  contract:
  - the f32 level's own rounding, held apart with its own constant 2 (not multiplied by C). The
    level's s and β are formed in f32, so e^{−β} and 1/sⁿ differ from the f64 reference's by
    (β + n) u. At the centre (ρ = 0, no conditioning term) and κ̂ = 16 (β = 27.7) this gave 5.8 of
    the contract without the term;
  - the flushes. îₙ and e below τ = 2⁻¹⁰⁰ are set to zero, as `bessel_j_real` flushes ĵₙ. Without
    the flush of e, κ̂ = 60 (β = 104) made e^{κρ − β} ≈ 10⁻⁴⁶: f32 underflowed to zero against the
    reference's 1.1·10⁻⁴⁶. A flushed term is below τ |q| √((2n + 1)/(4π)).

  No f32 coefficient is subnormal.
- Bits: vector v of three the one-vector call's; a complex vector's two planes the real-charge kernel
  on Re q and on Im q; a second call adds exactly; a multipole of the wrong length refused.
- `tetraleaf-exec/tests/helmholtz/densities.rs`: `p2m_kernel` on Yukawa within 10⁻¹³ of the largest
  coefficient of `p2m`'s, every path, three vectors; complex vectors as two real ones under the
  kernel bit for bit.

The real- and complex-k paths' operations are unchanged: their tests hold as before, and the FMM's
real and complex checksums in the tests are the same.

## 3. The FMM against the direct sum (Kathleen)

The batched P2M's error over the reference P2M's, the larger departure from 1 of φ's and ∇φ's
(`ch2-1c-tables.py accuracy acc`):

| Workload | ε | κD = 0.01 | κD = 0.1 | κD = 1 | κD = 10 | κD = 100 | κD = 1000 |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| W1 (N = 10⁵) | 1e-3 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| W1 (N = 10⁵) | 1e-6 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| W1 (N = 10⁵) | 1e-8 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| W3 (N = 10⁵) | 1e-3 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| W3 (N = 10⁵) | 1e-6 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| W3 (N = 10⁵) | 1e-8 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| W4 cube | 1e-3 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| W4 cube | 1e-6 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| W4 cube | 1e-8 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| W4 sphere | 1e-3 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| W4 sphere | 1e-6 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| W4 sphere | 1e-8 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| W4 Plummer | 1e-3 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| W4 Plummer | 1e-6 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| W4 Plummer | 1e-8 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| W4 clusters | 1e-3 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| W4 clusters | 1e-6 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| W4 clusters | 1e-8 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| W3 (N = 10⁶) | 1e-6 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |

114 cases, none above 1.01; the harness prints errors to three digits, so 1.0000 means equal there.
Every error is within ε, the largest 0.267 ε (∇φ, W4's Plummer sphere at 10⁻⁸), CH2.1b's. In 13
cases the two outputs' checksums are equal: the cases at κD = 1000 whose far field changes no bit of
the output (W1, W3 at 10⁵ and 10⁶, the W4 cube and sphere; CH2.1b §2: their errors are the P2P's
round-off). Elsewhere they differ, as the contract allows.

Complex charges, W1 at 10⁵, 10⁻⁶: κD = 1 φ 4.91·10⁻⁸ under both P2Ms, κD = 100 2.34·10⁻⁸ under
both (1.0000×).

The constant density (the worst density the record names, CH4.3c3 and CH4.8; `worst`), 10⁻⁶:

| case | φ, kernel | φ, reference | ÷ | ∇φ, kernel | ∇φ, reference | ÷ |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| W1 (N = 10⁶), κD = 1 | 9.95e-10 | 9.95e-10 | 1.0000 | 1.51e-7 | 1.51e-7 | 1.0000 |
| W1 (N = 10⁶), κD = 100 | 8.32e-9 | 8.32e-9 | 1.0000 | 2.57e-6 | 2.57e-6 | 1.0000 |
| W3 (N = 10⁶), κD = 1 | 4.89e-10 | 4.89e-10 | 1.0000 | 3.10e-8 | 3.10e-8 | 1.0000 |
| W3 (N = 10⁶), κD = 100 | 2.35e-8 | 2.35e-8 | 1.0000 | 3.12e-8 | 3.12e-8 | 1.0000 |

W1's ∇φ at κD = 100 is 2.6 ε under both P2Ms: the constant density on Yukawa misses ε on the
gradient whatever the P2M. The plateau degree (15 at 10⁻⁶) is the rule's for Re k > 0; CH2.3
calibrates Yukawa's degrees and already names the constant density among its cases (its Also from the
H4 gate's second run). It is reported there as a finding (§5), not a defect of this task.

## 4. Time (Kathleen, 1 × 40)

Per case four launches, kernel, reference, reference, kernel; the means of each setting's two (the
spread the larger of the two settings' max ÷ min of the evaluation). P2M the sum over the levels of
the second evaluation's P2M (`ch2-1c-tables.py runs time`):

| case | κD | P2M s, kernel | P2M s, reference | P2M ÷ | P2M share, reference | evaluate s, kernel | evaluate s, reference | evaluate ÷ (spread) | φ ÷ | ∇φ ÷ |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | ---: | ---: |
| W1 (N = 2·10⁶) | 1 | 0.0406 | 0.1330 | 0.305 | 0.099 | 1.245 | 1.340 | 0.929 (1.058) | 1.0000 | 1.0000 |
| W1 (N = 2·10⁶) | 100 | 0.0441 | 0.1400 | 0.315 | 0.104 | 1.270 | 1.340 | 0.948 (1.099) | 1.0000 | 1.0000 |
| W3 (N = 10⁶) | 1 | 0.0214 | 0.0744 | 0.288 | 0.067 | 1.035 | 1.105 | 0.937 (1.028) | 1.0000 | 1.0000 |
| W3 (N = 10⁶) | 100 | 0.0255 | 0.0759 | 0.337 | 0.070 | 1.000 | 1.080 | 0.926 (1.038) | 1.0000 | 1.0000 |

Per level (`levels time`; levels whose reference P2M is below 1 ms left out):

| case | κD | level | p | P2M s, kernel | P2M s, reference | ÷ |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| W1 (N = 2·10⁶) | 1 | 5 | 15 | 0.0034 | 0.0105 | 0.319 |
| W1 (N = 2·10⁶) | 1 | 6 | 15 | 0.0368 | 0.1220 | 0.302 |
| W1 (N = 2·10⁶) | 100 | 5 | 15 | 0.0036 | 0.0111 | 0.326 |
| W1 (N = 2·10⁶) | 100 | 6 | 15 | 0.0401 | 0.1286 | 0.312 |
| W3 (N = 10⁶) | 1 | 4 | 15 | 0.0153 | 0.0561 | 0.274 |
| W3 (N = 10⁶) | 1 | 5 | 15 | 0.0058 | 0.0181 | 0.320 |
| W3 (N = 10⁶) | 100 | 4 | 15 | 0.0191 | 0.0587 | 0.326 |
| W3 (N = 10⁶) | 100 | 5 | 15 | 0.0063 | 0.0170 | 0.371 |

The ratio is close to real k's on Kathleen (CH4.9a2b: 0.25–0.31×), although îₙ runs per lane by the
scalar `modified_i` (real k has the batched `bessel_j_real`): at p = 15 the harmonics and the sum
dominate. A batched `modified_i` would be a change to `tetraleaf-math`, not measured here.

Yukawa's L2P with gradients, still the reference's per-target calls, takes 0.30–0.32 s of W1's
evaluations (24–25 %) and 0.16 s of W3's (15–16 %) in the kernel runs; real and complex k's L2P runs
CH4.9b2's kernel, 0.17–0.27× its reference.

## 5. Questions for sign-off

Signed off 2026-11-28, the user on #275: "Q1 and Q2 as recommended".

1. **Yukawa's batched L2P.** Add CH2.1e, `LeafL2p` on Yukawa's real path (restating `l2p_yukawa`,
   with gradients; îₙ and îₙ′ per lane; the contract as CH4.9b2's; the FMM within 1.01×), to the
   priority list in item 15 right after CH2.1c (recommended: L2P is 15–25 % of Yukawa's evaluation
   after this task, and CH4.9b2's kernel gained 0.17–0.27× on it for real k); or add it off the list;
   or leave Yukawa's L2P on the reference.
2. **The constant density's gradient on Yukawa.** W1 at 10⁶ points, κD = 100, 10⁻⁶: ∇φ 2.6 ε under
   both P2Ms (§3). Record it on CH2.3's line as a case its calibration must bring within ε
   (recommended: CH2.3 already runs the constant density and owns Yukawa's degrees; done on this
   branch, with a note that W2 needs a fixed-N form for Yukawa); or open a task of its own before
   CH2.3.
