# Helmholtz CH2.1b: the Yukawa operator and FMM

Task: docs/helmholtz/tasks/CH2.1b.md (roadmap CH2.1b, item 15 of the priority list; proposal §10
CH2.1, §5.5, §6.3, §12 D4). The operator is `tetraleaf_exec::helmholtz::HelmholtzOperator` on
CH2.1a's one-plane tables (the module docs, "Yukawa"); the FMM is `HelmholtzFmmBuilder` with
k = iκ. Raw outputs are in `~/projects/fmm-archive/kathleen/runs/ch2-1b/` (`acc`, `time`, `time2`).
The tables come from `tools/helmholtz/ch2-1b-tables.py`; the launches from
`tools/helmholtz/ch2-1b-runs.sh`.

## 0. Summary

- **The acceptance's runs.** On W1 (its sphere at 10⁵ points), W3 (10⁵ and 10⁶ points) and W4's
  four distributions (10⁵ points), at κD = 10⁻², 10⁻¹, 1, 10, 10² and 10³ and ε = 10⁻³, 10⁻⁶ and
  10⁻⁸, the FMM under the library's defaults is within ε against the sampled direct sum in all 114
  cases. The largest error is 0.27 ε (∇φ, W4's Plummer sphere at 10⁻⁸). The degree rule gives its
  plateau at every κ (x = Re(k) a = 0): within target holds, and CH2.3's calibration has room (§2).
- **Deterministic.** The checksum of the whole output is the same at 1 × 1, 1 × 40 and 4 × 10 ranks ×
  threads on W1 and W4's Plummer sphere at κD = 1 and 100 (§3); the tests add 2 and 4 threads, pair
  by pair, RCR by pairs and the one-rank FMM at 2 and 4 ranks, bit for bit.
- **Real and complex k unchanged.** On six cases (real k with complex and with real charges, complex
  k), under the defaults and under Rotation, `main`'s and the branch's checksums agree in all 12 pairs
  (§4).
- **The complex path's time (D4): no cost measured.** Real and complex k on four cases, `main` against
  the branch, ABBA twice: the evaluation's median ratio is 0.991–1.040 (geometric mean 1.015) on the
  default builds and 0.947–1.035 (1.003) on builds with `-C llvm-args=-align-loops=64`; the medians'
  ranges overlap in every case and the sign of each case's difference changes between the two builds,
  which is run-to-run spread and code placement, not work (§5). The branch adds one test of the kind
  per level call or leaf body to the real and complex paths and nothing to any inner loop.

## 1. Setup

- **Machine.** Kathleen, one exclusive node (2 × 20 Cascade Lake cores), `singlenode`, slow nodes
  excluded (`node-c11b-004`, `node-c11b-013`), all on node-c11b-046: job 673280 (accuracy,
  determinism, Dense, complex charges, checksums; branch 503ad58), 673278 (the first timing; branch
  e0f648b) and 673279 (the timing with aligned builds; branch 864ebc0, the library as e0f648b), and
  `main` at 2413c5b (job 673276). Job 673277 ran the accuracy stages on e0f648b, before the
  self-review's fix of the V rule's prices on Yukawa (§2). One rank × 40 threads unless named, BLAS
  threads 1, release builds.
- **Settings.** The library's defaults named (`--m2l auto-pw --diagonal auto --diagonal-build
  interpolated`, the V and W/X rules of the cost model). On Yukawa, Auto runs Rotation on every level,
  no level takes the waves or the diagonal form, and the V rule prices a pair at Rotation's cost (the
  machine's constants, measured for the real and complex paths; its conversions are exact).
- **Problem.** Yukawa is `--kd 0 --kim κ` (D = 1). Real charges, §9.1's uniform parts, sources equal
  targets, gradients. The oracle is `tetraleaf_ref::helmholtz::p2p::direct_sum` at k = iκ (Yukawa's
  e^{−κr}/r, bit for bit `direct_sum_yukawa`), at 1,000 sampled targets. W1 at a fixed N
  (`Workload::W1n`, the example's `w1:<n>`): Yukawa has no wavelength to set W1's N. W3 and W4's
  cube are both uniform in the cube from the same seed, so their rows agree to the printed digit.

## 2. The acceptance: the FMM against the direct sum

max(φ, ∇φ) relative ℓ2 over the samples, per workload, ε and κD (`ch2-1b-tables.py accuracy acc2`):

| Workload | ε | κD = 0.01 | κD = 0.1 | κD = 1 | κD = 10 | κD = 100 | κD = 1000 | max / ε |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| W1 (N = 10⁵) | 1e-3 | 4.41e-05 | 4.41e-05 | 4.41e-05 | 4.14e-05 | 4.85e-06 | 4.78e-14 | 4.4e-02 |
| W1 (N = 10⁵) | 1e-6 | 8.04e-08 | 8.04e-08 | 8.04e-08 | 7.62e-08 | 2.31e-08 | 4.78e-14 | 8.0e-02 |
| W1 (N = 10⁵) | 1e-8 | 3.64e-10 | 3.64e-10 | 3.64e-10 | 3.57e-10 | 1.92e-10 | 4.78e-14 | 3.6e-02 |
| W3 (N = 10⁵) | 1e-3 | 4.73e-05 | 4.73e-05 | 4.72e-05 | 6.47e-05 | 3.55e-07 | 1.28e-13 | 6.5e-02 |
| W3 (N = 10⁵) | 1e-6 | 1.39e-07 | 1.39e-07 | 1.39e-07 | 1.31e-07 | 1.61e-08 | 1.28e-13 | 1.4e-01 |
| W3 (N = 10⁵) | 1e-8 | 6.43e-10 | 6.43e-10 | 6.43e-10 | 6.19e-10 | 1.39e-10 | 1.28e-13 | 6.4e-02 |
| W4 cube | 1e-3 | 4.73e-05 | 4.73e-05 | 4.72e-05 | 6.47e-05 | 3.55e-07 | 8.55e-14 | 6.5e-02 |
| W4 cube | 1e-6 | 1.39e-07 | 1.39e-07 | 1.39e-07 | 1.31e-07 | 1.61e-08 | 8.55e-14 | 1.4e-01 |
| W4 cube | 1e-8 | 6.43e-10 | 6.43e-10 | 6.43e-10 | 6.19e-10 | 1.39e-10 | 8.55e-14 | 6.4e-02 |
| W4 sphere | 1e-3 | 2.10e-05 | 2.11e-05 | 2.19e-05 | 1.91e-05 | 1.24e-06 | 2.19e-12 | 2.2e-02 |
| W4 sphere | 1e-6 | 3.62e-08 | 3.64e-08 | 3.78e-08 | 3.64e-08 | 5.88e-09 | 2.19e-12 | 3.8e-02 |
| W4 sphere | 1e-8 | 1.44e-10 | 1.45e-10 | 1.51e-10 | 1.55e-10 | 4.83e-11 | 2.19e-12 | 1.5e-02 |
| W4 Plummer | 1e-3 | 8.29e-05 | 8.29e-05 | 8.29e-05 | 8.27e-05 | 7.31e-05 | 4.45e-06 | 8.3e-02 |
| W4 Plummer | 1e-6 | 2.09e-07 | 2.09e-07 | 2.09e-07 | 2.09e-07 | 1.88e-07 | 2.79e-08 | 2.1e-01 |
| W4 Plummer | 1e-8 | 2.67e-09 | 2.67e-09 | 2.67e-09 | 2.67e-09 | 2.61e-09 | 6.42e-10 | 2.7e-01 |
| W4 clusters | 1e-3 | 5.90e-05 | 5.90e-05 | 5.90e-05 | 5.93e-05 | 7.30e-05 | 2.09e-05 | 7.3e-02 |
| W4 clusters | 1e-6 | 1.99e-07 | 1.99e-07 | 1.99e-07 | 1.99e-07 | 1.96e-07 | 9.02e-08 | 2.0e-01 |
| W4 clusters | 1e-8 | 5.20e-10 | 5.20e-10 | 5.20e-10 | 5.20e-10 | 5.10e-10 | 2.31e-10 | 5.2e-02 |

| W3 (N = 10⁶), ε = 1e-6 | κD = 0.01 | κD = 0.1 | κD = 1 | κD = 10 | κD = 100 | κD = 1000 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| max(φ, ∇φ) | 1.28e-07 | 1.28e-07 | 1.28e-07 | 1.71e-07 | 9.66e-08 | 1.14e-13 |

Readings:

- **The Laplace limit** (κD ≤ 1): the errors do not move with κ, as the kernel is Laplace's to
  O(κ r) on every level. 0.02–0.27 ε.
- **Strong damping** (κD = 1000 on the uniform and sphere trees): the far field is e^{−κ r} below
  the near field at the leaf distance, so the error is the P2P's round-off (10⁻¹⁴–10⁻¹²) at every
  ε. Where leaves are small (Plummer, clusters: leaves of a few 10⁻³), the far field still counts
  and the errors stay near the Laplace limit's or fall.
- **No miss** under the plateau's degrees (7, 15 and 24 at 10⁻³, 10⁻⁶ and 10⁻⁸ on every level from
  2, at every κ); the room (0.02–0.27 ε, and degrees that ignore the damping) is CH2.3's
  calibration and pruning (θ_D). The expansion's loss under strong damping (the roadmap's CH2.3, G3)
  shows in the operator's far chain at κD = 100 (a pair |t| ≈ 1.17 apart on level 2, about 3·10⁻⁷
  relative to a field e^{−119}; `tests/helmholtz/operator.rs`), not in these FMM errors: such pairs
  carry nothing measurable in the sum.
- **Before the self-review's fix** (job 673277: the V rule priced Yukawa's Auto levels of degree at
  most Batched's crossover at Batched's cost, though they run Rotation, so it converted fewer V pairs
  to exact P2P) every case was within ε as well, at most 0.27 ε; the fix lowers some errors by up to
  13 % (W4 clusters at 10⁻⁶: 2.29·10⁻⁷ to 1.99·10⁻⁷) and changes no case's standing.

## 3. Determinism

Checksums of the whole output (`ch2-1b-tables.py determinism acc2`), 10⁻⁶:

| case | 1 × 1 | 1 × 40 | 4 × 10 | the matrix's 1 × 40 | equal |
| --- | --- | --- | --- | --- | --- |
| w1, κD = 1, 1e-6 | f055d4a00dece261 | f055d4a00dece261 | f055d4a00dece261 | f055d4a00dece261 | yes |
| w4-plummer, κD = 1, 1e-6 | 5673a5887a692faf | 5673a5887a692faf | 5673a5887a692faf | 5673a5887a692faf | yes |
| w1, κD = 100, 1e-6 | f4ddc8fa0ca03df7 | f4ddc8fa0ca03df7 | f4ddc8fa0ca03df7 | f4ddc8fa0ca03df7 | yes |
| w4-plummer, κD = 100, 1e-6 | 85e3289d7ae50851 | 85e3289d7ae50851 | 85e3289d7ae50851 | 85e3289d7ae50851 | yes |

The tests (`tetraleaf-exec/tests/helmholtz_fmm.rs`: `check`, `check_vectors` and `check_yukawa` on
three Yukawa scenarios, `check_rcr_split` on one) add 2 and 4 threads, pair by pair, RCR by pairs and,
at 2 and 4 ranks, the one-rank FMM bit for bit.

Dense and complex charges (`ch2-1b-tables.py extra acc2`), 10⁻⁶:

| case | M2L | φ | ∇φ |
| --- | --- | ---: | ---: |
| W4 cube, κD = 1 | Rotation | 6.19e-8 | 1.39e-7 |
| W4 cube, κD = 1 | Dense (levels 2, 3, 331 MB) | 6.19e-8 | 1.39e-7 |
| W4 cube, κD = 30 | Rotation | 1.29e-7 | 8.82e-8 |
| W4 cube, κD = 30 | Dense (levels 2, 3, 331 MB) | 1.29e-7 | 8.82e-8 |
| W1 (N = 10⁵), κD = 1, complex charges | Rotation | 4.91e-8 | 8.37e-8 |
| W1 (N = 10⁵), κD = 100, complex charges | Rotation | 2.34e-8 | 3.01e-8 |

Dense differs from Rotation by rounding (the checksums differ; the tests hold the outputs within
`DENSE_TOLERANCE` 1e-14, 2.8·10⁻¹⁵ measured).

## 4. Real and complex k: `main` against the branch

Checksums at 1 × 40, 10⁻⁶, complex charges but `w1real` (`ch2-1b-tables.py checksums acc2`):

| case | M2L | main | branch | equal |
| --- | --- | --- | --- | --- |
| w1-100-0 | defaults | 3b2814f58c3e9de5 | 3b2814f58c3e9de5 | yes |
| w1-100-0 | rotation | 064b2e7a182b2216 | 064b2e7a182b2216 | yes |
| w1-30-0 | defaults | a106373b16a89858 | a106373b16a89858 | yes |
| w1-30-0 | rotation | 2a700e6360395bf5 | 2a700e6360395bf5 | yes |
| w1real-100-0 | defaults | 597b438b6807fbd0 | 597b438b6807fbd0 | yes |
| w1real-100-0 | rotation | d606340e0c6ea5a5 | d606340e0c6ea5a5 | yes |
| w3-30-0 | defaults | 490ac0e6b28fc825 | 490ac0e6b28fc825 | yes |
| w3-30-0 | rotation | e7ac9a8be62faebf | e7ac9a8be62faebf | yes |
| w4-cube-20-2 | defaults | 58673d70a8f5b47c | 58673d70a8f5b47c | yes |
| w4-cube-20-2 | rotation | 943472c4faf747fb | 943472c4faf747fb | yes |
| w4-plummer-10-0 | defaults | c5bf688e230268e0 | c5bf688e230268e0 | yes |
| w4-plummer-10-0 | rotation | 58f93ca6ade84f78 | 58f93ca6ade84f78 | yes |

The defaults' rows run the waves, Batched, the diagonal form and RCR by pairs; every existing test of
the workspace passes unchanged. Every timing launch of §5 also gave `main`'s checksum on the branch.

## 5. The complex path's time (D4)

`main` (2413c5b) against the branch, real and complex k (`ch2-1b-tables.py time <dir> <prefix>`): one
setting per launch, the order main, branch, branch, main twice (four launches each), 10⁻⁶, the
defaults, `--levels`, three evaluations (the last's time; the "Build" table's, to the millisecond),
medians and their ranges, the slowest rank's. W1 at kD = 100, 250 and 500 (real k; the defaults take
the waves, Batched and the diagonal form) and W3 at 10⁶ points at kD = 50 + 5i (complex k).

Default builds (job 673279, `time-*`):

| case | N | build s main | branch | ratio | evaluate s main | branch | ratio | evaluations | checksums equal |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| w1-100-0 | 79578 | 0.124 | 0.124 | 1.000 | 0.089 | 0.090 | 1.011 | 4 + 4 | yes |
|   evaluate s range | | | | | 0.082–0.098 | 0.087–0.096 | | | |
| w1-250-0 | 497360 | 0.639 | 0.622 | 0.973 | 0.551 | 0.546 | 0.991 | 4 + 4 | yes |
|   evaluate s range | | | | | 0.535–0.566 | 0.530–0.561 | | | |
| w1-500-0 | 1989437 | 2.742 | 2.730 | 0.995 | 2.468 | 2.566 | 1.040 | 4 + 4 | yes |
|   evaluate s range | | | | | 2.442–2.503 | 2.430–2.733 | | | |
| w3-50-5 | 1000000 | 0.780 | 0.784 | 1.006 | 1.312 | 1.336 | 1.018 | 4 + 4 | yes |
|   evaluate s range | | | | | 1.281–1.352 | 1.316–1.348 | | | |

Geometric mean of the ratios (branch / main): build 0.994, evaluate 1.015

Both trees built with `RUSTFLAGS="-C llvm-args=-align-loops=64"` (job 673279, `aligned-*`; the method
of CH2.0 for telling code placement from work):

| case | N | build s main | branch | ratio | evaluate s main | branch | ratio | evaluations | checksums equal |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| w1-100-0 | 79578 | 0.124 | 0.125 | 1.004 | 0.086 | 0.089 | 1.035 | 4 + 4 | yes |
|   evaluate s range | | | | | 0.081–0.094 | 0.087–0.091 | | | |
| w1-250-0 | 497360 | 0.607 | 0.627 | 1.031 | 0.523 | 0.542 | 1.035 | 4 + 4 | yes |
|   evaluate s range | | | | | 0.515–0.531 | 0.519–0.572 | | | |
| w1-500-0 | 1989437 | 2.736 | 2.757 | 1.008 | 2.465 | 2.458 | 0.998 | 4 + 4 | yes |
|   evaluate s range | | | | | 2.449–2.501 | 2.449–2.552 | | | |
| w3-50-5 | 1000000 | 0.782 | 0.783 | 1.001 | 1.320 | 1.250 | 0.947 | 4 + 4 | yes |
|   evaluate s range | | | | | 1.286–1.362 | 1.236–1.305 | | | |

Geometric mean of the ratios (branch / main): build 1.011, evaluate 1.003

The first timing (job 673278, two evaluations, three cases) gave 1.073, 1.051 and 0.986 for the
evaluation (geometric mean 1.036); the second, with heavier cases, does not repeat it, and in both
second sets every case's ranges overlap. W1 at kD = 250 is 0.991× on the default builds and 1.035× on
the aligned ones; W3's complex case 1.018× and 0.947×. Read with the record's run-to-run and
placement effects on Kathleen (±5–13 % from code moves, CH2.0; 1.09–1.20× for a second setting in a
launch, CH4.3b2c), no cost of the kind test to the complex path is measured, nor would one be expected
of a branch per leaf body or level call. Not measured: locust and the M6 (the change has no
ISA-dependent part; one machine class).

## 6. What the branch changes on the real and complex paths

- `Kernel::{multipole_columns, local_columns}` return the densities' planes and m as the calls passed
  before (one match on the kind per M2M, L2L and M2L call); P2L, M2P and P2P test the kind once per
  leaf body and return before the existing code on Yukawa.
- `select_tables` maps Auto to Rotation and skips the diagonal form on Yukawa only; `wave_levels`
  returns early for Re k = 0 only; Dense's selection takes `select_for(…, 2)`, `select`'s own body.
- The builder's V rule maps Auto to Rotation on Yukawa only.
