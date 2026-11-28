# CH2.2a The RCR tables per level in a store

Started 2026-11-28 on `tbetcke/h-ch2-2a`. Roadmap: CH2.2 (item 15 of the priority list, after
CH2.1e): "Per-level tables with per-level cache keys. *Done when:* build cost and memory reported;
levels reused across an octave sweep. *Needs:* CH2.1, D4." Split on 2026-11-28 by `/next-task`'s
size rule into CH2.2a (this task: the RCR tables of every path, Yukawa's whole M2L) and CH2.2b (the
M2L families: Dense, Classes, Batched, the exponential M2L, the diagonal form). Proposal: §10 CH2.2
("Per-level tables in `tetraleaf-tables` with per-level cache keys | build cost and memory reported;
reuse of levels across an octave sweep shown | CH2.1"); §6.3 ("Table cache"); §7 ("Build cost");
§8 R10; §2 (constraints 1, 2, 5, 7, 9).

## Goal
Key every part of a `HelmholtzTables` build by what it depends on (k̂ of its level, its degrees, its
geometric class), keep the parts in a store that operators share, so that a second FMM on the same
domain, at the same k or at k doubled, builds only the parts no earlier build made; report the
build per level (built or reused, time, bytes).

## Acceptance (the roadmap's *Done when* for CH2.2a; proposal §10 CH2.2 for its share)
- Tables built with a store, cold (empty store) and warm (every part in it), are the tables built
  without one bit for bit, through every application (M2M and L2L of every octant, M2L of every
  offset, every level), on real and complex k and Yukawa, f64 and f32, on any `Parts` map.
- The FMM with a store gives the output without one bit for bit (cold and warm; threads; ranks).
- Across an octave sweep on one domain (k, 2k, 4k, …: k̂_{l+1}(2k) = k̂_l(k) exactly), every level's
  parts but the new coarse level's are reused, shown per level in the build's report: in the tests
  (tables and FMM) and on Kathleen for real k (W1n at fixed N, kD 62.5–1000, under the library's
  defaults and under Rotation) and for Yukawa (κD 15.625–1000), with each level's build time and
  bytes, cold against warm.
- Build cost and memory reported: per level the parts built and reused, their thread-seconds, the
  level's bytes, and the store's bytes (distinct parts).

## Read first
- `tetraleaf-tables/src/helmholtz.rs`: `HelmholtzTables`, `LevelTables`, `Shift`, `Task`, `Part`,
  `build_by`, `set_m2l_blocks`, `bytes`, `block_selection`.
- `tetraleaf-tables/src/parts.rs`; `tetraleaf-tables/tests/helmholtz/parts.rs` (`Threads`,
  `applications`, `CASES`).
- `tetraleaf-exec/src/helmholtz/operator.rs`: `unbuilt`, `select_tables`, `TableTimings`,
  `table_timings`; `fmm.rs`: `HelmholtzFmmBuilder` (`block_budget`), `HelmholtzFmm::build_timings`;
  `tests/helmholtz/build.rs`; `tests/helmholtz_fmm.rs` (`check_build_timings`, `yukawa_scenarios`).
- `tetraleaf-validate/src/helmholtz_accuracy.rs` (`Case`, `Row::build_parts`) and
  `examples/helmholtz_accuracy.rs` (`--kd`, `--kim`, the table "Build").
- Reports: `helmholtz-ch4-4a.md` §7 (build per k; RCR tables 9.2–12.9 s of W1's build at
  kD = 1000), `helmholtz-ch4-0h.md` (the build by part).

## Changes
- `tetraleaf-tables/src/helmholtz/store.rs` (new, `tetraleaf_tables::helmholtz::store`):
  - `TableStore<T>`: a map from `PartKey` to the part (`Arc`), behind a `Mutex`, with a byte cap
    (`TableStore::new(cap)`, `TableStore::unbounded()`) and least-recently-used eviction by a build
    counter (ties by key order, so deterministic): after each build the oldest entries go until the
    store holds at most `cap` bytes. Evicting drops only the store's reference: tables in use keep
    theirs. `bytes`, `entries`, `clear`.
  - `PartKey` (crate-visible variants): `Up { k_hat, p_in, p_out }`, `Down { k_hat, p_in, p_out }`
    (k̂ of the child level as the bits of its real and imaginary parts), `M2l { k_hat, q, distance }`
    (|d|², an integer), `Blocks { degree, angle }` (the polar angle's integer key; octant and M2L
    blocks alike, no k̂: y-blocks depend on the degree and the angle alone). T is the store's type
    parameter. No convention version or format: the store is in memory, one process.
  - `LevelBuild`: per level the parts built and reused, the built parts' thread-seconds (each part's
    clock, in the part's closure), the level's bytes. `TableBuild` (per level, and the store's bytes
    after the build).
- `tetraleaf-tables/src/helmholtz.rs`:
  - `LevelTables` and `Shift` hold their parts as `Arc` (the M2M, L2L, M2L shifts and the block
    sets), so tables from a store share memory; the applications unchanged (they read through the
    `Arc`).
  - `HelmholtzTables::build_in(k̂₀, degrees, &TableChoice, &impl Parts, &TableStore<T>) -> Self`
    (as built: the report through `build_report() -> &TableBuild`, kept with the tables, for builds
    with or without a store; `level_bytes`): the parts the store holds are taken from it, the
    others (each distinct key once) computed by `parts` and put in; `build_by` keeps its bits and
    computes every part (no sharing without a store). `set_m2l_blocks_in` likewise. `bytes` counts
    each level's values as before (shared parts counted per level).
  - The module docs: "The store".
- `tetraleaf-exec`:
  - `HelmholtzOperator::with_table_store(Arc<TableStore<T>>)` (recorded by `unbuilt`, read by
    `select_tables`'s RCR build and `set_m2l_blocks`), `table_build()` (the last RCR build's
    `TableBuild`), `TableTimings` unchanged.
  - `HelmholtzFmmBuilder::table_store(TableStores)` (as built: the builder is not generic over T, so
    `TableStores` holds one store per precision and `build::<T>` takes T's; default none: today's
    builds, bit for bit and part for part), `HelmholtzFmm::table_build()`. The store is per process: each rank keeps its
    own (or shares one between the FMMs it builds).
  - Module docs ("The table store"), crate `CLAUDE.md`.
- `tetraleaf-validate`: `helmholtz_accuracy::Case::table_store` (an `Option<Arc<TableStore<f64>>>`,
  `None` in `Case::new` and `Case::defaults`) and `Row::table_build` (per level, the slowest rank's
  seconds, parts summed over the ranks); the example's `--table-store MB` (one store for the whole
  launch, shared by every case in order), `--kim` as a list (a Yukawa sweep in one launch), and with
  `--levels` a table "RCR tables by level" (k̂_l, p_l, q_l, parts built and reused, seconds, MB).
- `tools/helmholtz/ch2-2a-runs.sh` (stages `smoke`, `real`, `rotation`, `yukawa`; each with the
  store and without, in separate launches) and `ch2-2a-tables.py` (on `common.py`).
- Report `tetraleaf-validate/results/helmholtz-ch2-2a.md`.
- No convention change; no dependency; `tetraleaf-plan` untouched.

## Tests that prove it
- `tetraleaf-tables/tests/helmholtz/store.rs` (trusted path: `build_by` without a store, through
  `parts.rs`' `applications`, bit for bit):
  - cold and warm builds on `CASES` (real, complex) and a Yukawa case, on `Serial` and `Threads`:
    the applications bit for bit `build_by`'s, `bytes` equal; cold: every part built; warm: every
    part reused, none built;
  - the octave shift: k̂₀ then 2k̂₀ with the degrees moved down one level (the coarse level new): the
    second build reuses every part of levels ≥ 3 and of the pairs below level 2, builds level 2's M2L
    and the top pair's, and its tables are a cold build's bit for bit (real, complex, Yukawa);
  - blocks keyed by degree: two levels at one degree build their M2L blocks once;
  - `set_m2l_blocks_in` against `set_m2l_blocks`;
  - the cap: 0 keeps nothing and changes no bit; a cap between one level's parts and two keeps the
    last build's most recent parts, the same entries on every run; f32 against f32 without a store;
  - two builds on two threads sharing one store: both bit for bit a plain build.
- `tetraleaf-exec/tests/helmholtz/build.rs`: the operator with a store (cold, warm) against without:
  degrees, M2L per level, bytes and every pair body bit for bit, under Rotation, Auto and Yukawa;
  `table_build` per level.
- `tetraleaf-exec/tests/helmholtz_fmm.rs`, `check_table_store`: one real and one Yukawa scenario, the
  FMM with a store cold and warm and at k then 2k on the same points, each output bit for bit the
  storeless FMM's of the same setting (same ranks), at 1 and 2 threads; the 2k build's report reuses
  every RCR part of levels ≥ 3; every rank.
- `tetraleaf-validate/tests/helmholtz_accuracy.rs`: a two-case sweep with `Case::table_store`, the
  rows' `table_build` (the second reuses), checksums equal to the storeless rows'.

The H6 and H4 gate rules: no bound, degree, tolerance or rule changes, so there is no wider scope
and no pair-or-FMM evidence to pair; no bound is gated.

## Checks
`/verify` (fmt; Clippy default, all-targets, workspace; `cargo test`, `--workspace` at opt-level 1;
docs; the MPI list at 2 and 4 ranks; tables' and exec's crate checks); `cargo deny` not needed (no
manifest change unless one appears).

## Machine runs
Kathleen, `singlenode`, one node, 1 × 40 (the store per process; build times on 40 threads):
`tools/helmholtz/ch2-2a-runs.sh real|rotation|yukawa`, each a launch with `--table-store` over the
whole sweep (W1n at N = 10⁶, kD 62.5, 125, 250, 500, 1000 at ε = 10⁻⁶; Yukawa `--kd 0 --kim
15.625,…,1000`) and a launch without (the cold builds), checksums equal. Decides nothing: it shows
the reuse and reports the build per level. Times in one launch are not compared across settings
(CH4.3f's drift): cold and warm are separate launches.

## Constraints and risks
- Constraint 1: no Laplace path (the Helmholtz `helmholtz` module of `tetraleaf-tables`, the Helmholtz
  operator and builder, the harness). Laplace's `TableCache` and its format untouched.
- Constraint 2: tables bit for bit with and without the store, so the output is; the store reads no
  rank's data.
- Constraint 5: the parts are built as before (coaxial tables in f64 rounded to T, y-blocks in T).
- Constraint 7: a store keeps tables beyond one FMM's life; its cap bounds it, and the node budget
  counts the FMM's tables as before (shared parts counted per FMM).
- Constraint 9: keys from k̂_l, the degrees and integer classes only.
- Risk: k̂ keyed by bits; a sweep reuses only where 2k is exact and the domain the same (both hold
  for doublings of kD on one point set). Different domains share only the y-blocks.
- Risk: in an octave sweep the new coarse level is the most expensive (p³ per coaxial table), so the
  time saved may be a small part of the RCR build; reported as measured.

## Pending runs
None. Kathleen, `run.sbatch`, `singlenode`, `--exclude=node-c11b-004,node-c11b-013`, each the release
build of `helmholtz_accuracy` then the stages; every launch exit 0:
- job 673293 on node-c11b-046, 2026-11-28 16:55–16:59 (build of 8cd5834): `smoke`, `real`, `yukawa`;
- job 673294 on node-c11a-012, 16:55–16:59 (8cd5834): `rotation`;
- job 673295 on node-c11b-046, 17:03–17:04 (4feef81): `repeat` (added after the sweeps: the same k
  twice, the bound of an on-disk store's gain).
Outputs in `~/projects/fmm-archive/kathleen/runs/ch2-2a/` (`rotation/`, `repeat/`, the jobs' logs).

## Results (the report `tetraleaf-validate/results/helmholtz-ch2-2a.md`)
- Bit for bit: the tests (tables through every application, cold, warm, octave, f32, threads; the
  operator's pair bodies; the FMM at 1 rank and in CI's list) and all 16 Kathleen checksums against
  the storeless launches.
- An octave on W1 at 10⁶ points reuses every RCR part of levels 3–6 (23 of 246–295 parts built at
  kD 500 and 1000 under the defaults; Yukawa 21 of 344): the RCR tables' wall time 0.88–0.94× cold at
  kD 500–1000 (the new level 2 holds 77–80 % of the thread-seconds), the build 0.96–1.00×.
- The same k again: the build 0.66× (4.17 s against 6.28 s at kD = 1000), the RCR tables 0 s.
- Memory: 1,842 MB of RCR tables at kD = 1000 (1,319 MB on level 2); the store 1.69 GB after the
  sweep; Yukawa 6.6 MB per FMM.
- Finding: under the defaults at kD ≥ 500 the RCR M2L tables of levels that run every class on the
  diagonal form are built and never read (most of 2.0 s of a 6.0 s build at kD = 1000).

## Open questions
For sign-off (the PR's questions; the report §5):
1. The on-disk store of §6.3: recommended none for now, §6.3 amended to the in-memory store with
   the measured gains; or a CH2.2c after CH2.2b.
2. The RCR M2L tables the diagonal form never reads: recommended an Also on CH2.2b (written on the
   branch); or a task of its own in H4, off the list.
