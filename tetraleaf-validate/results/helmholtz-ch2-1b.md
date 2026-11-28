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
  plateau at every κ (x = Re(k) a = 0): within target holds, but CH2.3's calibration has room (§2).
- **Deterministic.** The checksum of the whole output is the same at 1 × 1, 1 × 40 and 4 × 10 ranks ×
  threads on W1 and W4's Plummer sphere at κD = 1 and 100 (§3); the tests add 2 and 4 threads, pair
  by pair and the one-rank FMM at 2 and 4 ranks, bit for bit.
- **Real and complex k unchanged.** On five cases (real k: W1 at kD = 30 and 100, W3 at kD = 30,
  W4's Plummer sphere at kD = 10; complex k: W4's cube at kD = 20 + 2i), under the defaults and under
  Rotation, `main`'s and the branch's checksums agree in all 10 pairs (§4).
- **The complex path's time (D4).** §5.

## 1. Setup

- **Machine.** Kathleen, one exclusive node (2 × 20 Cascade Lake cores), `singlenode`, slow nodes
  excluded (`node-c11b-004`, `node-c11b-013`): jobs 673277 (accuracy, determinism, Dense, complex
  charges, checksums), 673278 (the first timing) and 673279 (the timing with aligned builds), all on
  node-c11b-046. One rank × 40 threads unless named, BLAS threads 1, release builds of branch
  commit e0f648b (673277, 673278) and 864ebc0 (673279; no library change since e0f648b) and of `main`
  at 2413c5b (job 673276).
- **Settings.** The library's defaults named (`--m2l auto-pw --diagonal auto --diagonal-build
  interpolated`, the V and W/X rules of the cost model). On Yukawa, Auto runs Rotation on every level,
  no level takes the waves or the diagonal form, and the near rules use the machine's cost constants
  as for real k (their conversions are exact).
- **Problem.** Yukawa is `--kd 0 --kim κ` (D = 1). Real charges, §9.1's uniform parts, sources equal
  targets, gradients. The oracle is `tetraleaf_ref::helmholtz::p2p::direct_sum` at k = iκ (Yukawa's
  e^{−κr}/r, bit for bit `direct_sum_yukawa`), at 1,000 sampled targets. W1 at a fixed N
  (`Workload::W1n`, the example's `w1:<n>`): Yukawa has no wavelength to set W1's N. W3 and W4's
  cube are both uniform in the cube from the same seed, so their rows agree to the printed digit.

## 2. The acceptance: the FMM against the direct sum

max(φ, ∇φ) relative ℓ2 over the samples, per workload, ε and κD (`ch2-1b-tables.py accuracy acc`):

| Workload | ε | κD = 0.01 | κD = 0.1 | κD = 1 | κD = 10 | κD = 100 | κD = 1000 | max / ε |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| W1 (N = 10⁵) | 1e-3 | 4.52e-05 | 4.52e-05 | 4.52e-05 | 4.25e-05 | 4.95e-06 | 4.78e-14 | 4.5e-02 |
| W1 (N = 10⁵) | 1e-6 | 8.28e-08 | 8.28e-08 | 8.27e-08 | 7.86e-08 | 2.42e-08 | 4.78e-14 | 8.3e-02 |
| W1 (N = 10⁵) | 1e-8 | 3.64e-10 | 3.64e-10 | 3.64e-10 | 3.57e-10 | 1.92e-10 | 4.78e-14 | 3.6e-02 |
| W3 (N = 10⁵) | 1e-3 | 4.73e-05 | 4.73e-05 | 4.72e-05 | 6.47e-05 | 3.55e-07 | 1.28e-13 | 6.5e-02 |
| W3 (N = 10⁵) | 1e-6 | 1.39e-07 | 1.39e-07 | 1.39e-07 | 1.31e-07 | 1.61e-08 | 1.28e-13 | 1.4e-01 |
| W3 (N = 10⁵) | 1e-8 | 6.43e-10 | 6.43e-10 | 6.43e-10 | 6.19e-10 | 1.39e-10 | 1.28e-13 | 6.4e-02 |
| W4 cube | 1e-3 | 4.73e-05 | 4.73e-05 | 4.72e-05 | 6.47e-05 | 3.55e-07 | 8.55e-14 | 6.5e-02 |
| W4 cube | 1e-6 | 1.39e-07 | 1.39e-07 | 1.39e-07 | 1.31e-07 | 1.61e-08 | 8.55e-14 | 1.4e-01 |
| W4 cube | 1e-8 | 6.43e-10 | 6.43e-10 | 6.43e-10 | 6.19e-10 | 1.39e-10 | 8.55e-14 | 6.4e-02 |
| W4 sphere | 1e-3 | 2.09e-05 | 2.10e-05 | 2.18e-05 | 1.90e-05 | 1.25e-06 | 2.19e-12 | 2.2e-02 |
| W4 sphere | 1e-6 | 3.64e-08 | 3.66e-08 | 3.80e-08 | 3.67e-08 | 5.92e-09 | 2.19e-12 | 3.8e-02 |
| W4 sphere | 1e-8 | 1.44e-10 | 1.45e-10 | 1.51e-10 | 1.55e-10 | 4.83e-11 | 2.19e-12 | 1.5e-02 |
| W4 Plummer | 1e-3 | 8.29e-05 | 8.29e-05 | 8.29e-05 | 8.27e-05 | 7.32e-05 | 4.45e-06 | 8.3e-02 |
| W4 Plummer | 1e-6 | 2.16e-07 | 2.16e-07 | 2.16e-07 | 2.16e-07 | 1.95e-07 | 3.22e-08 | 2.2e-01 |
| W4 Plummer | 1e-8 | 2.67e-09 | 2.67e-09 | 2.67e-09 | 2.67e-09 | 2.61e-09 | 6.42e-10 | 2.7e-01 |
| W4 clusters | 1e-3 | 5.91e-05 | 5.91e-05 | 5.91e-05 | 5.94e-05 | 7.31e-05 | 2.10e-05 | 7.3e-02 |
| W4 clusters | 1e-6 | 2.29e-07 | 2.29e-07 | 2.29e-07 | 2.29e-07 | 2.26e-07 | 1.16e-07 | 2.3e-01 |
| W4 clusters | 1e-8 | 5.20e-10 | 5.20e-10 | 5.20e-10 | 5.20e-10 | 5.10e-10 | 2.31e-10 | 5.2e-02 |

| W3 (N = 10⁶), ε = 1e-6 | κD = 0.01 | κD = 0.1 | κD = 1 | κD = 10 | κD = 100 | κD = 1000 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| max(φ, ∇φ) | 1.28e-07 | 1.28e-07 | 1.28e-07 | 1.71e-07 | 9.66e-08 | 1.14e-13 |

Degrees by level, the same at every κ (the rule's plateau at x = 0): 7, 15 and 24 at 10⁻³, 10⁻⁶ and
10⁻⁸ on every level from 2. Readings:

- **The Laplace limit** (κD ≤ 1): the errors do not move with κ, as the kernel is Laplace's to
  O(κ r) on every level. 0.02–0.27 ε.
- **Strong damping** (κD = 1000 on the uniform and sphere trees): the far field is e^{−κ r} below
  the near field at the leaf distance, so the error is the P2P's round-off (10⁻¹⁴–10⁻¹²) at every
  ε. Where leaves are small (Plummer, clusters: leaves of a few 10⁻³), the far field still counts
  and the errors stay near the Laplace limit's or fall (Plummer 3.2·10⁻⁸ at 10⁻⁶).
- **No miss** under the plateau's degrees; the room (0.02–0.27 ε, and degrees that ignore the
  damping) is CH2.3's calibration and pruning (θ_D). The RCR M2L's loss under strong damping (the
  roadmap's CH2.3, G3) shows in the operator's far chain at κD = 100 (a pair at |t| ≈ 1.17 on level
  2, 3·10⁻⁷ relative to a field e^{−117} of the undamped one; `tests/helmholtz/operator.rs`), not in
  these FMM errors: such pairs carry nothing measurable in the sum.

## 3. Determinism

Checksums of the whole output (`ch2-1b-tables.py determinism acc`), 10⁻⁶:

| case | 1 × 1 | 1 × 40 | 4 × 10 | the matrix's 1 × 40 | equal |
| --- | --- | --- | --- | --- | --- |
| W1 (N = 10⁵), κD = 1 | cff5767e8ff1ee43 | cff5767e8ff1ee43 | cff5767e8ff1ee43 | cff5767e8ff1ee43 | yes |
| W4 Plummer, κD = 1 | 3a1b3cbfca6da01e | 3a1b3cbfca6da01e | 3a1b3cbfca6da01e | 3a1b3cbfca6da01e | yes |
| W1 (N = 10⁵), κD = 100 | 2f011f2816e20a0b | 2f011f2816e20a0b | 2f011f2816e20a0b | 2f011f2816e20a0b | yes |
| W4 Plummer, κD = 100 | acb78b0151783060 | acb78b0151783060 | acb78b0151783060 | acb78b0151783060 | yes |

The tests (`tetraleaf-exec/tests/helmholtz_fmm.rs`, `check_yukawa` and `check` on three Yukawa
scenarios) add 2 and 4 threads, pair by pair and, at 2 and 4 ranks, the one-rank FMM bit for bit.

Dense and complex charges (`ch2-1b-tables.py extra acc`), 10⁻⁶:

| case | M2L | φ | ∇φ |
| --- | --- | ---: | ---: |
| W4 cube, κD = 1 | Rotation | 6.19e-8 | 1.39e-7 |
| W4 cube, κD = 1 | Dense (levels 2, 3, 331 MB) | 6.19e-8 | 1.39e-7 |
| W4 cube, κD = 30 | Rotation | 1.29e-7 | 8.82e-8 |
| W4 cube, κD = 30 | Dense (levels 2, 3, 331 MB) | 1.29e-7 | 8.82e-8 |
| W1 (N = 10⁵), κD = 1, complex charges | Rotation | 4.94e-8 | 8.61e-8 |
| W1 (N = 10⁵), κD = 100, complex charges | Rotation | 2.41e-8 | 3.09e-8 |

Dense differs from Rotation by rounding (the checksums differ; the tests hold the outputs within
`DENSE_TOLERANCE` 1e-14, 2.8·10⁻¹⁵ measured).

## 4. Real and complex k: `main` against the branch

Checksums at 1 × 40, 10⁻⁶ (`ch2-1b-tables.py checksums acc`):

| case | M2L | main | branch | equal |
| --- | --- | --- | --- | --- |
| W1, kD = 100 | defaults | 3b2814f58c3e9de5 | 3b2814f58c3e9de5 | yes |
| W1, kD = 100 | Rotation | 064b2e7a182b2216 | 064b2e7a182b2216 | yes |
| W1, kD = 30 | defaults | a106373b16a89858 | a106373b16a89858 | yes |
| W1, kD = 30 | Rotation | 2a700e6360395bf5 | 2a700e6360395bf5 | yes |
| W3 (N = 10⁵), kD = 30 | defaults | 490ac0e6b28fc825 | 490ac0e6b28fc825 | yes |
| W3 (N = 10⁵), kD = 30 | Rotation | e7ac9a8be62faebf | e7ac9a8be62faebf | yes |
| W4 cube, kD = 20 + 2i | defaults | 58673d70a8f5b47c | 58673d70a8f5b47c | yes |
| W4 cube, kD = 20 + 2i | Rotation | 943472c4faf747fb | 943472c4faf747fb | yes |
| W4 Plummer, kD = 10 | defaults | c5bf688e230268e0 | c5bf688e230268e0 | yes |
| W4 Plummer, kD = 10 | Rotation | 58f93ca6ade84f78 | 58f93ca6ade84f78 | yes |

The defaults' rows run the waves, Batched, the diagonal form and RCR by pairs; every existing test of
the workspace passes unchanged.

## 5. The complex path's time (D4)

TIMING
