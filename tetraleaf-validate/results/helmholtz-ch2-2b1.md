# Helmholtz CH2.2b1: RCR's M2L tables only where an RCR pair reads them

Task: docs/helmholtz/tasks/CH2.2b1.md (roadmap CH2.2b1, the first of CH2.2b's three parts, item 15
of the priority list; the Also of #277's question 2, from `helmholtz-ch2-2a.md` §5). Proposal
§6.3 ("Table cache"), §7 ("Build cost"). Raw outputs are in
`~/projects/fmm-archive/kathleen/runs/ch2-2b1/` (the launches, `runs.log`, the jobs' logs; the
failed W3 launch at its default N in `failed-w3-5e8/`); the tables come from
`tools/helmholtz/ch2-2b1-tables.py`, the launches from `tools/helmholtz/ch2-2b1-runs.sh`.

## 0. Summary

- **The change.** `HelmholtzOperator::select_tables` now chooses every level's M2L, the diagonal
  form's classes included, before it builds the RCR tables. It then builds RCR's M2L coaxial
  tables, and stores RCR's M2L y-blocks, only on the levels with an RCR pair: a level whose M2L
  is Rotation's and that the diagonal form does not take whole. Before, every level from 2 got
  RCR M2L tables, including the levels of the waves, Batched, Dense and Classes, and the
  diagonal form's levels with all 34 classes; none of them reads those tables.
- **Bit for bit.** Every output is the output with RCR M2L tables on every level:
  - in the tests, against the doc-hidden hook `with_m2l_tables_everywhere` (the tables before):
    the operator's pair bodies under every M2L setting, and the FMM at 1, 2 and 4 ranks, 2 threads
    and against the one-rank FMM;
  - on Kathleen, the checksums of `main` (63a7a4f) and of the branch agree in every launch (§2).
- **The build.** Under the library's defaults on W1's sphere at 10⁶ points every level from 2
  takes all 34 classes at kD ≥ 500, so no level reads RCR's M2L tables (§2):
  - kD = 1000: the RCR tables 0.455 s against 1.933 s (0.24×), the build 4.25 s against 5.61 s
    (0.76×), the RCR tables 363 MB against 1,842 MB;
  - kD = 500: the RCR tables 0.24×, the build 0.91×, 57 MB against 398 MB;
  - kD = 250: the RCR tables 0.25×, the build 0.98× (level 6 takes 21 classes and keeps its M2L
    tables).
  - W3 at 10⁶ points and kD = 500 (§3): the RCR tables 0.27×, the build 0.95×, 64 MB against
    395 MB.
- **What remains of the RCR tables** is the translations (M2M, L2L and the octant blocks), whose
  coarse pair dominates: at kD = 1000 the pair between levels 1 and 2 (p = 241) holds 151 MB and
  the most thread-seconds.

## 1. Setup

- **Machine.** Kathleen, one exclusive node, `singlenode`, node-c11b-004 and node-c11b-013 excluded.
  One rank × 40 threads, BLAS threads 1.

  | Job | Node | Time (2026-11-28) | Stages | Builds |
  | --- | --- | --- | --- | --- |
  | 673308 | node-c11a-058 | 19:13–19:54 | `smoke`, `w1` (W3 at its default N failed, then cancelled) | `main` 63a7a4f in `src/main`, the branch b20d0a6 (library as c423598) |
  | 673310 | node-c11b-037 | 19:55–19:56 | `w3` | the binaries of 673308 |

- **Launches.** Each case runs four launches, `main`, branch, branch, `main`, one setting per
  launch (review 2026-11-24, B12). The sides are the mean of their two launches.
- **Workloads.** W1's sphere at 10⁶ points (`w1:1000000`, CH2.2a's tree: 7 levels, 256 points per
  leaf), kD 250, 500, 1000; W3 at 10⁶ points (`w3:1000000`, CH4.3b2c's) at kD = 500. ε = 10⁻⁶,
  gradients, §9.1's samples, complex charges, the library's defaults named (CH2.2a's).
- **Measures.** "RCR tables s" is the build's part `TableTimings::rcr` (wall time, the slowest
  rank's); "diagonal tables s" `TableTimings::diagonal`. MB and thread-seconds are the table "RCR
  tables by level" summed over the levels, of each side's first launch.

## 2. W1 (job 673308)

| case | build s (main) | build s (branch) | ÷ | RCR tables s (main) | RCR tables s (branch) | ÷ | diagonal tables s (main) | diagonal tables s (branch) | RCR MB (main) | RCR MB (branch) | RCR thread-s (main) | RCR thread-s (branch) | φ | ∇φ | checksums |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| smoke | 0.071 | 0.043 | 0.61 | 0.001 | 0.001 | 0.50 | 0.012 | 0.012 | 3.6 | 1.6 | 0.014 | 0.005 | 6.07e-8 | 1.51e-7 | same |
| w1-250 | 1.354 | 1.323 | 0.98 | 0.045 | 0.011 | 0.25 | 0.145 | 0.150 | 174.0 | 11.6 | 0.895 | 0.037 | 4.65e-8 | 9.93e-8 | same |
| w1-500 | 1.891 | 1.717 | 0.91 | 0.234 | 0.061 | 0.26 | 0.472 | 0.463 | 397.6 | 57.4 | 4.191 | 0.176 | 3.64e-8 | 5.81e-8 | same |
| w1-1000 | 5.613 | 4.253 | 0.76 | 1.933 | 0.455 | 0.24 | 1.730 | 1.826 | 1842.2 | 363.2 | 33.538 | 1.192 | 4.73e-8 | 5.61e-8 | same |

`smoke` is W1 at 20,000 points, kD = 30 (milliseconds). "checksums same": the four launches of the
case agree.

Per level at kD = 1000 and kD = 250 (`main`'s first launch against the branch's):

| case | level | k̂_l | p | q | built (main) | built (branch) | thread-s (main) | thread-s (branch) | MB (main) | MB (branch) |
| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| w1-1000 | 0 | 500.0 | 0 | 0 | 4 | 4 | 0.0001 | 0.0002 | 0.00 | 0.00 |
| w1-1000 | 1 | 250.0 | 0 | 0 | 4 | 4 | 0.2375 | 0.1198 | 151.18 | 151.18 |
| w1-1000 | 2 | 125.0 | 241 | 241 | 19 | 4 | 27.9921 | 0.9144 | 1319.12 | 178.29 |
| w1-1000 | 3 | 62.5 | 129 | 129 | 19 | 4 | 4.2309 | 0.1312 | 205.70 | 27.91 |
| w1-1000 | 4 | 31.25 | 71 | 71 | 68 | 4 | 0.9348 | 0.0236 | 132.86 | 4.82 |
| w1-1000 | 5 | 15.63 | 41 | 41 | 68 | 4 | 0.1093 | 0.0027 | 26.50 | 0.99 |
| w1-1000 | 6 | 7.81 | 26 | 26 | 64 | 0 | 0.0329 | 0.0000 | 6.81 | 0.00 |
| w1-250 | 2 | 31.25 | 71 | 71 | 68 | 4 | 0.7876 | 0.0234 | 132.86 | 4.82 |
| w1-250 | 3 | 15.63 | 41 | 41 | 68 | 4 | 0.0679 | 0.0027 | 26.50 | 0.99 |
| w1-250 | 4 | 7.81 | 26 | 26 | 68 | 4 | 0.0177 | 0.0009 | 7.08 | 0.27 |
| w1-250 | 5 | 3.91 | 17 | 17 | 68 | 4 | 0.0110 | 0.0004 | 2.12 | 0.09 |
| w1-250 | 6 | 1.95 | 15 | 15 | 64 | 64 | 0.0038 | 0.0040 | 1.43 | 1.43 |

- **The choice per level.** Under the defaults the time rule puts all 34 classes of levels 2–6 on
  the diagonal form at kD = 500 and 1000 ("diagonal 2:34,…,6:34" in the M2L column), so the
  branch builds no M2L table on any level: the 4 parts per level left are the translations of the
  pair below it (level 6, the leaves, has none). At kD = 250 level 6 takes 21 classes; its RCR
  pairs read its M2L tables, which the branch builds (64 parts, as on `main`).
- **Time.** The RCR tables fall to 0.24–0.26× at kD 250–1000. At kD = 1000 the build loses 1.36 s
  (0.76×): 1.48 s fewer in the RCR tables, 0.10 s more in the diagonal tables (untouched code). The diagonal tables'
  time, which this task does not touch, is 1.73 s on `main` and 1.83 s on the branch, so the
  diagonal form's tables are now the largest table part of the build (CH2.2b2 keys them).
- **The evaluation** (no code of it changed): 8.05 and 8.12 s on `main`'s launches at kD = 1000,
  8.24 and 8.43 s on the branch's; at kD = 500 2.09 and 2.13 against 2.05 and 2.18. Not decided on.

## 3. W3 at 10⁶ points (job 673310)

| case | build s (main) | build s (branch) | ÷ | RCR tables s (main) | RCR tables s (branch) | ÷ | diagonal tables s (main) | diagonal tables s (branch) | RCR MB (main) | RCR MB (branch) | RCR thread-s (main) | RCR thread-s (branch) | φ | ∇φ | checksums |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| w3-500 | 3.099 | 2.934 | 0.95 | 0.251 | 0.068 | 0.27 | 0.449 | 0.457 | 395.3 | 64.0 | 4.180 | 0.226 | 7.24e-8 | 9.94e-8 | same |

| level | k̂_l | p | q | built (main) | built (branch) | thread-s (main) | thread-s (branch) | MB (main) | MB (branch) |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 2 | 62.5 | 129 | 129 | 19 | 4 | 3.1749 | 0.1538 | 205.70 | 27.91 |
| 3 | 31.25 | 71 | 71 | 68 | 4 | 0.8792 | 0.0314 | 132.86 | 4.82 |
| 4 | 15.63 | 41 | 41 | 68 | 4 | 0.0633 | 0.0029 | 26.50 | 0.99 |
| 5 | 7.81 | 26 | 26 | 64 | 64 | 0.0165 | 0.0160 | 6.81 | 6.81 |

- Levels 2–4 take all 34 classes; level 5 (the leaves' level, p = 26) runs RCR and keeps its M2L
  tables, which its pairs read. The RCR tables fall to 0.27×, the build to 0.95× (0.17 s of
  3.10 s); the RCR tables hold 64 MB against 395 MB. The checksums of the four launches agree.
- The evaluation: 5.17 and 5.16 s on `main`, 5.38 and 5.07 s on the branch.
- W3 at its default N (`w3`, 5.0·10⁸ points at kD = 500) failed in job 673308: the redistribution
  overflows MPI's i32 counts, after a 1,010 s oracle (exit 0; `failed-w3-5e8/`). The plan named
  CH4.3b2c's 10⁶ points from then on.

## 4. What the tests hold

- `tetraleaf-exec/tests/helmholtz/build.rs`, `the_rcr_m2l_tables_are_built_where_read`: Rotation, Auto
  at 10⁻⁶ and 5·10⁻⁸, PlaneWave, Batched, Dense, Classes, and the diagonal form at 10⁻¹⁰ (levels 4 and
  5 taken whole, levels 2 and 3 with 32 classes): M2L tables and stored blocks exactly on the levels
  with an RCR pair (restated from `m2l_on` and the class count), every level from 2 with the hook,
  every pair body bit for bit the hook's, fewer bytes where a level lost its tables. With the old
  choice it fails on its first case (Auto at 10⁻⁶, level 2).
- `tetraleaf-exec/tests/helmholtz_fmm.rs`, `check_m2l_tables` (in `check_diagonal`, the cube at
  kD = 40 + 2i and the clusters at kD = 20): the diagonal rule under Rotation and under the
  library's Auto, and the form at 0.1 under Rotation, which takes levels whole (levels 2 and 3 on
  the cube, 2–4 on the clusters): the levels on every rank, the output bit for bit the hook's; the
  whole-level run at 2 threads and against the one-rank FMM. At 1, 2 and 4 ranks on the M6.

## 5. Questions

None: the task changes which tables are built, not a rule, a default's output or a convention.
