# Helmholtz CH2.1e: Yukawa's batched L2P

Task: docs/helmholtz/tasks/CH2.1e.md (roadmap CH2.1e, item 15 of the priority list, after CH2.1c;
proposal §10 CH2.1, §5.5). The kernel is `tetraleaf_tables::helmholtz::leaf::LeafL2p` on Yukawa's
real path (its module docs, "L2P"), run by `Densities::l2p_kernel` in the FMM. Raw outputs are in
`~/projects/fmm-archive/kathleen/runs/ch2-1e/` (`acc`, `time`, the job's log); the tables come from
`tools/helmholtz/ch2-1e-tables.py`, the launches from `tools/helmholtz/ch2-1e-runs.sh`.

## 0. Summary

- **The kernel.** `LeafL2p` now restates `l2p_yukawa` (k = iκ): îₙ and îₙ′ by `modified_i` per lane
  at the reference's length, e^{κr − β} per lane, the three flushed below τ as CH2.1c's P2M flushes
  them, CH4.9b2's per-degree sums per real plane, then the reference's real degree loop with κ and
  the centre's limit κ/(3s). `LeafOutput::planes` gives the output columns per vector (1 on Yukawa's
  real path, 2 otherwise). Every output is within CH4.9b2's contract with ĵ → î, k → κ and e's
  condition added (C = 10): the worst measured is 1.94 of C in f64 (p ≤ 45, κ̂ ≤ 250) and 0.64 in f32
  (§2). Real and complex k's L2P is `main`'s bit for bit (a hash over 192 calls).
- **The FMM's errors within 1.01× the per-point path's, and within ε.** On CH2.1b's 114 cases the
  batched leaves' errors are the reference leaves' to the harness's printed digits (1.0000×) in every
  case, the largest 0.27 ε as in CH2.1b and CH2.1c; also with complex charges and on the constant
  density (§3). The tests: 1.0000×, outputs within 4.1·10⁻¹⁶ on Yukawa's three scenarios.
- **Deterministic.** A target's outputs depend on its offset and the locals alone (the bit tests:
  vectors, columns, positions, leaves of 1, 16 and 40). The FMM tests' threads, pair calls, vectors,
  complex charges as two real vectors and the one-rank FMM at 2 and 4 ranks hold bit for bit.
- **Time.** L2P per level 0.21–0.28× the per-point path's on Kathleen (1 × 40), per run 0.23–0.24×;
  Yukawa's L2P falls from 14–22 % of an evaluation to 4–6 %. The evaluation with both leaf kernels
  (CH2.1c and CH2.1e) is 0.77–0.89× the reference leaves' (§4).

## 1. Setup

- **Machine.** Kathleen, one exclusive node (2 × 20 Cascade Lake cores), `singlenode`, slow nodes
  excluded: job 673289 on node-c11b-046, 2026-11-28 14:56–15:02 (the release build of
  `helmholtz_accuracy`, then `smoke`, `accuracy`, `worst` and `levels`; all 114 launches exit 0).
  One rank × 40 threads, BLAS threads 1, the release build of f5656e8 (the library as 5e7850b; the
  later commits change docs and tests only).
- **Settings.** As CH2.1c's: the library's defaults named (`--m2l auto-pw --diagonal auto
  --diagonal-build interpolated`, the V and W/X rules of the cost model; on Yukawa Auto runs Rotation
  on every level), gradients, §9.1's samples, real charges unless named, Yukawa as `--kd 0 --kim κ`
  with D = 1. One leaf setting per launch: `--leaves kernel` (the batched P2M and L2P, the library's
  default) or `--leaves reference` (the per-point reference calls of both; P2P the AVX-512 kernel in
  both). The reference setting is therefore the state before CH2.1c and CH2.1e together; the per-level
  L2P times isolate CH2.1e.
- **Workloads.** W1 runs as W1n (`w1:<n>`, CH2.1b) at a fixed N, W2 is not run (no fixed-N form on
  Yukawa; CH2.3's Also), as at CH2.1c.

## 2. The contract per target (the M6)

`tetraleaf-tables/tests/helmholtz/l2p_yukawa.rs`, against `tetraleaf_ref::helmholtz::leaf::l2p_yukawa`
per target and real local, in the FMM's leaf frame (half-width 1). With wₙ = √((2n + 1)/(4π)),
Sₙ = Σₘ |L̂ᵞₙᵐ|, x = κr, tₙ = îₙ/r (κ/(3s) for n = 1 at r = 0), e = e^{κr − β}:

```text
|φ̃ − φ̃_ref|  ≤ C u e Σₙ wₙ Sₙ (((n + 1)² + 2n + 1) îₙ + x (|îₙ′| + îₙ))
|gₐ − gₐ_ref| ≤ C u e Σₙ wₙ Sₙ (κ (((n + 1)² + 2n + 3) |îₙ′| + x îₙ + x |îₙ′|) + ((n + 1)³ + n(n + 1) + x) tₙ)
```

CH4.9b2's terms with î for ĵ and κ for k, plus e's condition (x times each term: r by `sqrt` against
the reference's `hypot`), as CH2.1c added it to P2M's.

| Test | Range | Worst, of C = 10 |
| --- | --- | ---: |
| f64 | p ∈ {0, 1, 2, 5, 8, 17, 26, 45}, κ̂ ∈ {10⁻⁶, 0.3, 4, 16, 60, 250}, N = 1–40 (centre, poles, axes, next to the z-axis), one or two planes, 1–3 vectors, with and without gradients | 1.94 |
| f32 against the f64 reference on the f32 inputs | p ≤ 32, κ̂ ≤ 60, one or two planes | 0.64 |
| r² underflowing (f64, r to 10⁻²⁰⁰) and the centre | p ≤ 17, κ̂ ≤ 60 | 0.22 |
| next to the centre in f32 (r = 10⁻²⁰, 10⁻²⁵, 3·10⁻³⁰) | p ≤ 17, κ̂ ≤ 4 | 0.18 |

In f32 three terms are held apart from C (not multiplied by it), as at CH2.1c: the f32 level's own
rounding of s and β (2 (β + n + 1) u e Σₙ wₙ Sₙ times îₙ, or κ|îₙ′| + tₙ), the harmonics' flushed
columns (2⁻⁵⁰ per degree's sum, CH4.9b2), and the flushed îₙ, îₙ′ and e (τ Σₙ wₙ Sₙ (1 + κ + 1/r),
absolute: at κ̂ = 60 in f32 e is flushed for r below about 0.58, where the kernel's outputs are
exact zeros and the check is this term alone; the near-centre f32 case runs at κ̂ ≤ 4 to keep the
arithmetic in view). Two mutations of the kernel (κ dropped from the radial gradient; the centre's
limit set to zero) fail the f64, f32 and underflow tests.

The bits, by `to_bits`: each column of m vectors of one or two planes the one-plane, one-vector
call's on its local; a target alone, in another position and in leaves of 1, 16 and 40; the
potentials with gradients those without; accumulation. Real and complex k: a hash of `LeafL2p`'s
outputs over k̂ ∈ {0.3, 4, 16, 4 + 0.5i, 16 + 2i, 10⁻⁶}, p ∈ {0, 3, 12, 26}, 1 and 3 vectors, with and
without gradients, `ee41f6d8571cd78e` on the branch and on `main`'s kernel.

## 3. The FMM's errors (Kathleen, job 673289)

kernel ÷ reference of the error against the direct sum, the larger of φ's and ∇φ's ratio's distance
from 1, per case (`ch2-1e-tables.py accuracy`):

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

114 cases: kernel ÷ reference 1.0000 in every case (φ and ∇φ, the harness's three printed digits),
none above 1.01; the largest error over ε 0.267 (within ε in all 114: CH2.1's acceptance); the
checksums equal in 13 (the outputs differ by rounding elsewhere, as expected of a contract).

| Complex charges, W1 (N = 10⁵), 10⁻⁶ | φ, kernel | φ, reference | ÷ |
| --- | ---: | ---: | ---: |
| κD = 1 | 4.91e-8 | 4.91e-8 | 1.0000 |
| κD = 100 | 2.34e-8 | 2.34e-8 | 1.0000 |

The constant density (the worst the record names; CH4.3c3, CH4.8), N = 10⁶, 10⁻⁶:

| case | φ, kernel | φ, reference | ÷ | ∇φ, kernel | ∇φ, reference | ÷ |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| W1n, κD = 1 | 9.95e-10 | 9.95e-10 | 1.0000 | 1.51e-7 | 1.51e-7 | 1.0000 |
| W1n, κD = 100 | 8.32e-9 | 8.32e-9 | 1.0000 | 2.57e-6 | 2.57e-6 | 1.0000 |
| W3, κD = 1 | 4.89e-10 | 4.89e-10 | 1.0000 | 3.10e-8 | 3.10e-8 | 1.0000 |
| W3, κD = 100 | 2.35e-8 | 2.35e-8 | 1.0000 | 3.12e-8 | 3.12e-8 | 1.0000 |

W1n at κD = 100 misses ε on ∇φ (2.6 ε) under either setting, CH2.1c's finding, already an Also on
CH2.3's line; the leaf kernels change nothing in it.

## 4. Time (Kathleen, job 673289, 1 × 40)

`--levels`, two evaluations (the rows the second's), per case kernel, reference, reference, kernel;
a setting's time the mean of its two launches; 10⁻⁶, real charges, gradients. Per run
(`ch2-1e-tables.py runs`):

| case | κD | L2P s, kernel | L2P s, reference | L2P ÷ | L2P share, reference | P2M ÷ | evaluate s, kernel | evaluate s, reference | evaluate ÷ (spread) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| W1n (N = 2·10⁶) | 1 | 0.0671 | 0.2893 | 0.232 | 0.218 | 0.338 | 1.045 | 1.325 | 0.789 (1.023) |
| W1n (N = 2·10⁶) | 100 | 0.0698 | 0.3030 | 0.231 | 0.219 | 0.336 | 1.070 | 1.385 | 0.773 (1.022) |
| W3 (N = 10⁶) | 1 | 0.0352 | 0.1557 | 0.226 | 0.140 | 0.329 | 0.960 | 1.115 | 0.861 (1.046) |
| W3 (N = 10⁶) | 100 | 0.0372 | 0.1550 | 0.240 | 0.141 | 0.356 | 0.980 | 1.100 | 0.891 (1.018) |

Per level (`ch2-1e-tables.py levels`; levels with under 1 ms of reference L2P left out):

| case | κD | level | p | L2P s, kernel | L2P s, reference | ÷ |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| W1n (N = 2·10⁶) | 1 | 5 | 15 | 0.0054 | 0.0234 | 0.233 |
| W1n (N = 2·10⁶) | 1 | 6 | 15 | 0.0614 | 0.2658 | 0.231 |
| W1n (N = 2·10⁶) | 100 | 5 | 15 | 0.0058 | 0.0246 | 0.235 |
| W1n (N = 2·10⁶) | 100 | 6 | 15 | 0.0638 | 0.2782 | 0.229 |
| W3 (N = 10⁶) | 1 | 4 | 15 | 0.0257 | 0.1205 | 0.214 |
| W3 (N = 10⁶) | 1 | 5 | 15 | 0.0095 | 0.0353 | 0.270 |
| W3 (N = 10⁶) | 100 | 4 | 15 | 0.0274 | 0.1197 | 0.229 |
| W3 (N = 10⁶) | 100 | 5 | 15 | 0.0099 | 0.0353 | 0.279 |

- L2P 0.21–0.28× per level, 0.23–0.24× per run: at p = 15 the same gain CH4.9b2 measured for real
  and complex k on Kathleen (0.19–0.23×), with `modified_i` per lane (scalar, as CH2.1c's P2M).
- Yukawa's L2P falls from 14–22 % of an evaluation (the reference) to 3.7–6.5 % (the kernel:
  0.0352/0.960 to 0.0698/1.070); P2M is 0.33–0.36× (CH2.1c's 0.29–0.34×).
- The evaluation 0.77–0.89×, spread within 1.05 per setting: the two leaf kernels together.
- The M6 and locust are not timed (one machine class, the kernel safe Rust without ISA dispatch, as
  at CH2.1b and CH2.1c).

## 5. Questions

None: the task is the roadmap line as signed off at CH2.1c (#275's question 1), and nothing here
moves a rule, a bound or a default other than the one it names (Yukawa's L2P follows the P2P choice
as P2M does since CH2.1c).
