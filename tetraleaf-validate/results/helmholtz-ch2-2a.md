# Helmholtz CH2.2a: the RCR tables per level in a store

Task: docs/helmholtz/tasks/CH2.2a.md (roadmap CH2.2a, the first half of CH2.2, item 15 of the
priority list; proposal §10 CH2.2, §6.3 "Table cache", §7 "Build cost", §8 R10). The store is
`tetraleaf_tables::helmholtz::store` (its module docs), used by `HelmholtzTables::build_in` and by the
operator through `HelmholtzFmmBuilder::table_store`. Raw outputs are in
`~/projects/fmm-archive/kathleen/runs/ch2-2a/` (`real-*`, `yukawa-*`, `smoke-*`, `rotation/`,
`repeat/`, the jobs' logs); the tables come from `tools/helmholtz/ch2-2a-tables.py`, the launches
from `tools/helmholtz/ch2-2a-runs.sh`.

## 0. Summary

- **Per-level keys.** Every part of an RCR table build is keyed by its own inputs. A translation
  is keyed by the child's k̂ (as bits) and its two degrees. An M2L coaxial table is keyed by k̂_l,
  q_l and |d|². A y-block set is keyed by its degree and the polar angle's integer key, with no
  k̂. Each part is held in a `TableStore` that builds share, behind a byte cap with
  least-recently-used eviction. Tables built with a store are the tables built without one, bit
  for bit, cold and warm. This holds in the tests (every application, real and complex k,
  Yukawa, f32, any map) and in every Kathleen run, whose 16 checksums each equal the storeless
  launch's.
- **Reuse across an octave sweep.** On one domain k̂_{l+1}(2k) = k̂_l(k) exactly. A sweep of kD in
  octaves on W1's sphere at 10⁶ points therefore reuses every RCR part of levels 3–6 and builds
  only the new coarse level 2 and the pair above it (§2, §3):
  - under the defaults at kD 500 and 1000, 23 parts of 246–295 are built;
  - for Yukawa, 21 of 344.
- **The time it saves is small, because the coarse level dominates.** The new level 2 holds
  77–80 % of a cold build's RCR thread-seconds at kD ≥ 500 (p = 129 and 241; O(p³) per coaxial
  table), and more of its wall time: its 15 coaxial tables are the longest parts, which the other
  parts finish beside. The RCR tables of a warm octave take 0.88–0.94× the cold build's wall time
  at kD = 500 and 1000, under the defaults and under Rotation alike (§2). The whole build is
  0.96–1.00×.
  Yukawa's RCR tables cost about 1 ms at ε = 10⁻⁶ (p = 15 on every level), so the store saves
  nothing measurable there (§3).
- **The same k again is the store's real gain.** Building kD = 1000 a second time in one process
  takes every part from the store. The build takes 4.17 s against 6.28 s cold (0.66×), and the RCR
  tables take 0 s against 2.28 s (§4).
- **Memory.** The tables of W1 at kD = 1000 under the defaults hold 1,842 MB, 1,319 MB of them on
  level 2. The store after the whole sweep holds 1,691 MB: the parts shared within a build are
  counted once. Yukawa's tables hold 6.6 MB per FMM and its store 3.3 MB after seven κ (§2, §3).
- **A finding outside this task (question 2).** Under the defaults at kD ≥ 500 every level takes
  all 34 classes on the diagonal form. Those levels' RCR M2L coaxial tables are still built, and
  they are never read. They are most of the RCR tables' 2.0 s in W1's 6.0 s build at kD = 1000.

## 1. Setup

- **Machine.** Kathleen, one exclusive node per job (2 × 20 Cascade Lake cores), `singlenode`,
  with node-c11b-004 and node-c11b-013 excluded. Each job ran the release build of
  `helmholtz_accuracy`, then the stages, and every launch exited 0:

  | Job | Node | Time (2026-11-28) | Stages | Release build |
  | --- | --- | --- | --- | --- |
  | 673293 | node-c11b-046 | 16:55–16:59 | `smoke`, `real`, `yukawa` | 8cd5834 (library as 655e8db) |
  | 673294 | node-c11a-012 | 16:55–16:59 | `rotation` | 8cd5834 |
  | 673295 | node-c11b-046 | 17:03–17:04 | `repeat` | 4feef81 (library as 655e8db; later commits change tests and docs only) |

  One rank × 40 threads, BLAS threads 1.
- **Workload.** W1's sphere at 10⁶ points (`w1:1000000`), the same points at every kD. The FMM's
  domain is therefore the same, k̂₀ = k r₀ with r₀ = 0.50000727, and doubling kD doubles k̂₀
  exactly. ε = 10⁻⁶, gradients, §9.1's samples, complex charges.
- **Settings.**
  - `real`: the library's defaults named (`--m2l auto-pw --diagonal auto --diagonal-build
    interpolated`, the V and W/X rules of the cost model).
  - `rotation`: RCR's M2L on every level.
  - `yukawa`: `--kd 0 --kim κ` under the defaults, which run Rotation on Yukawa.
- **Launches.** Each sweep ran twice, in two launches. The first used `--table-store 65536`: one
  store for the launch, so every kD after the first is warm. The second ran without a store, every
  kD cold. Cold and warm are compared across the two launches, never within one (review
  2026-11-24, B12).
- **The tree.** At 10⁶ points the tree has 7 levels at every kD (256 points per leaf). In an octave,
  level l + 1 at 2k is level l at k, for levels 3–6. Level 2 at 2k matches level 1 at k, which has
  no M2L: it is the new coarse level.
- **Measures.**
  - "RCR tables s" is the build's part `TableTimings::rcr`, wall time.
  - "Thread-s" sums each built part's own clock from `TableBuild`. It is an inflated measure on 40
    threads (contention): level 3 at kD = 1000 took 6.25 thread-s cold, while the same parts as
    level 2 at kD = 500 took 3.83.
  - MB is a level's or the store's bytes, each level counted whole.

## 2. Real k: the sweep under the defaults and under Rotation (jobs 673293, 673294)

Per kD, the launch with a store (warm from 125 on) against the launch without; `÷` is the RCR tables'
wall time, store ÷ cold.

The defaults:

| kD | build s (store) | build s (cold) | RCR tables s (store) | RCR tables s (cold) | ÷ | parts built | reused | thread-s (store) | thread-s (cold) | tables MB | store MB | checksums |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 62.5 | 1.138 | 1.141 | 0.003 | 0.002 | 1.50 | 236 | 10 | 0.068 | 0.063 | 11.7 | 11.13 | same |
| 125 | 1.153 | 1.183 | 0.002 | 0.003 | 0.67 | 70 | 225 | 0.045 | 0.072 | 38.4 | 36.84 | same |
| 250 | 1.332 | 1.413 | 0.037 | 0.059 | 0.63 | 70 | 274 | 0.699 | 1.031 | 174.0 | 165.72 | same |
| 500 | 2.049 | 2.065 | 0.259 | 0.275 | 0.94 | 23 | 272 | 3.478 | 5.005 | 397.6 | 371.42 | same |
| 1000 | 6.011 | 6.112 | 2.005 | 2.124 | 0.94 | 23 | 223 | 30.049 | 38.992 | 1842.2 | 1690.55 | same |

Rotation:

| kD | build s (store) | build s (cold) | RCR tables s (store) | RCR tables s (cold) | ÷ | parts built | reused | thread-s (store) | thread-s (cold) | tables MB | store MB | checksums |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 62.5 | 1.051 | 1.024 | 0.003 | 0.003 | 1.00 | 236 | 108 | 0.063 | 0.065 | 13.8 | 11.13 | same |
| 125 | 1.044 | 1.038 | 0.003 | 0.004 | 0.75 | 70 | 274 | 0.079 | 0.099 | 39.4 | 36.84 | same |
| 250 | 1.099 | 1.087 | 0.057 | 0.058 | 0.98 | 70 | 274 | 0.856 | 0.963 | 174.0 | 165.72 | same |
| 500 | 1.304 | 1.318 | 0.253 | 0.287 | 0.88 | 23 | 272 | 3.569 | 5.017 | 397.6 | 371.42 | same |
| 1000 | 3.660 | 3.826 | 2.058 | 2.230 | 0.92 | 23 | 223 | 31.105 | 42.165 | 1842.2 | 1690.55 | same |

The first kD is cold in both launches. Its "reused" counts parts shared within the build: y-blocks of
one degree and angle on two levels, and under Rotation the M2L blocks that the octants' angles
share.

Per level at kD = 1000 (the defaults; Rotation's agree in parts and bytes):

| level | k̂_l | p | q | built | reused | thread-s (store) | thread-s (cold) | MB |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0 | 500.0 | 0 | 0 | 2 | 2 | 0.0001 | 0.0001 | 0.00 |
| 1 | 250.0 | 0 | 0 | 4 | 0 | 0.2512 | 0.2829 | 151.18 |
| 2 | 125.0 | 241 | 241 | 17 | 2 | 29.7977 | 31.2078 | 1319.12 |
| 3 | 62.5 | 129 | 129 | 0 | 19 | 0 | 6.2546 | 205.70 |
| 4 | 31.25 | 71 | 71 | 0 | 68 | 0 | 1.1681 | 132.86 |
| 5 | 15.63 | 41 | 41 | 0 | 68 | 0 | 0.0629 | 26.50 |
| 6 | 7.81 | 26 | 26 | 0 | 64 | 0 | 0.0155 | 6.81 |

- **Which parts are built.** A warm octave builds:
  - level 2's 15 M2L coaxial tables at its new degree;
  - the translations of the pairs above level 3;
  - the y-blocks of level 2's new degree.

  Every part of levels 3–6 comes from the store. Level 3 holds 19 parts and not 68 because its M2L
  blocks do not fit the block budget (512 MiB, spent from the finest level up), so its pairs form
  them per call.
- **Time.** Level 2 alone holds 31.2 of the cold build's 39.0 thread-s at kD = 1000 (80 %), and
  3.8 of 5.0 at 500 (77 %). The proposal's "a sweep in octaves rebuilds only the top level per frequency" (§6.3, §7)
  holds part for part. But the top level is most of the cost: the wall time of the RCR tables falls
  by 6–12 % at kD = 500 and 1000. The 0.63–0.75 at kD = 125–250 are milliseconds.
- **Accuracy.** Every row is within ε:
  - the defaults: φ 3.6–5.5·10⁻⁸, ∇φ 5.6–12·10⁻⁸;
  - Rotation: φ 3.9–5.6·10⁻⁸.

  The checksums equal the storeless launch's at every kD.

## 3. Yukawa: the κ sweep (job 673293)

| κD | build s (store) | build s (cold) | RCR tables s (store) | RCR tables s (cold) | parts built | reused | thread-s (store) | thread-s (cold) | tables MB | store MB | checksums |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 15.625 | 1.078 | 1.075 | 0.001 | 0.002 | 138 | 206 | 0.012 | 0.037 | 6.6 | 2.07 | same |
| 31.25 | 1.070 | 1.066 | 0.001 | 0.001 | 21 | 323 | 0.002 | 0.022 | 6.6 | 2.27 | same |
| 62.5 | 1.067 | 1.068 | 0.001 | 0.001 | 21 | 323 | 0.002 | 0.025 | 6.6 | 2.47 | same |
| 125 | 1.047 | 1.081 | 0.001 | 0.001 | 21 | 323 | 0.002 | 0.018 | 6.6 | 2.68 | same |
| 250 | 1.069 | 1.060 | 0.001 | 0.001 | 21 | 323 | 0.002 | 0.017 | 6.6 | 2.88 | same |
| 500 | 1.055 | 1.064 | 0.001 | 0.001 | 21 | 323 | 0.002 | 0.019 | 6.6 | 3.08 | same |
| 1000 | 1.065 | 1.065 | 0.001 | 0.001 | 21 | 323 | 0.002 | 0.019 | 6.6 | 3.29 | same |

- **Degrees.** Yukawa's degree rule gives the plateau, p = 15, on every level at 10⁻⁶. Each level
  differs only in k̂_l, so a warm κ builds 21 parts: level 2's 15 M2L tables and the translations
  above level 3. Already the first, cold κ builds only 138 of 344 parts, since the y-blocks of p = 15
  serve every level.
- **Cost.** The RCR tables take about 1 ms of a 1.07 s build. The rest of the build is the octree,
  the plan and the leaf operators, which no store touches. Under Yukawa's defaults the store saves
  nothing measurable at this ε.
- **Accuracy.** Within ε at every κ: φ 1.9·10⁻⁹–7.0·10⁻⁸. Checksums equal.

## 4. The same k again (job 673295)

W1 at 10⁶ points, kD = 1000 twice in one launch under the defaults: with a store (the second build
fully warm), and without (both cold).

| run | build s (store) | build s (cold) | RCR tables s (store) | RCR tables s (cold) | parts built | reused | thread-s (store) | thread-s (cold) | store MB | checksums |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| first | 5.701 | 5.745 | 2.040 | 1.993 | 240 | 6 | 35.218 | 34.162 | 1686.23 | same |
| second | 4.168 | 6.278 | 0.000 | 2.283 | 0 | 246 | 0.000 | 41.260 | 1686.23 | same |

- **Saving.** The second build takes every part from the store: 2.1 s of a 6.3 s build, 0.66×.
- **What a store across processes would give.** This is the bound on what an on-disk store
  (proposal §6.3) would give a second process at the same k, before the cost of reading 1.7 GB.
- **What remains.** The rest of the build is the diagonal form's tables (`TableTimings::diagonal`,
  1.7–2.2 s at kD = 1000 in these launches), the plan (1.05 s), the octree (0.22 s) and the
  scratch. CH2.2b keys the
  diagonal tables.

## 5. Questions

Both signed off 2026-11-28 on #277 (the user: "Q1 and Q2 as recommended"): no on-disk store, the
proposal's §6.3 and §7 amended; the unread RCR M2L tables an Also on CH2.2b.

1. **The on-disk store.**
   - **The proposal.** §6.3 names an on-disk table cache, keyed per level, with a byte cap and LRU
     eviction.
   - **What CH2.2a built.** The in-memory form of it, one store per process (per rank), sharing
     parts by `Arc`. It meets §10's acceptance: the build cost and memory are reported, and the
     reuse across an octave sweep is shown.
   - **What a disk store would gain.** Measured here: an octave saves 6–12 % of the RCR tables'
     wall time, and the same k again saves 2.1 s of 6.3 s at kD = 1000 before any read cost.
   - **Recommended.** No on-disk store for now. Amend §6.3 to the in-memory store, with the
     measured gains, and leave a disk store to a later task if a user needs repeated processes at
     one k. The alternative is a CH2.2c, a disk store after CH2.2b.
2. **The RCR M2L tables the diagonal form never reads.**
   - **The finding.** Under the defaults at kD ≥ 500 every level of W1 takes all 34 classes on the
     diagonal form, yet `select_tables` builds the RCR M2L coaxial tables of every RCR-eligible
     level. It builds them before `diagonal::select` decides the classes, which reads k̂ from them.
   - **The cost.** Level 2's 15 tables at p = 241 are most of the RCR tables' 2.0 s in a 6.0 s build
     at kD = 1000.
   - **Recommended.** An "Also" on CH2.2b, which keys the diagonal tables in the same
     `select_tables`: decide the diagonal classes first, and build RCR M2L tables only for the
     levels with an RCR pair. This is the same output bit for bit, since no read changes. The
     alternative is a task of its own in H4, off the list.
