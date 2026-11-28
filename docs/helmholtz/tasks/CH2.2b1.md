# CH2.2b1 The RCR M2L tables only where an RCR pair reads them

Started 2026-11-28 on `tbetcke/h-ch2-2b1`. Roadmap: CH2.2b (item 15 of the priority list, after
CH2.2a), split on 2026-11-28 by `/next-task`'s size rule into CH2.2b1 (this task: the Also of
#277's question 2), CH2.2b2 (the diagonal form's tables in the store) and CH2.2b3 (Dense, Classes,
Batched and the exponential M2L in the store, and the sweep under the defaults). Proposal: §10
CH2.2 ("Per-level tables in `tetraleaf-tables` with per-level cache keys | build cost and memory
reported; reuse of levels across an octave sweep shown | CH2.1"), for which this task is
preparation, not a share of the criterion; §6.3 ("Table cache"), §7 ("Build cost"); §2
(constraints 1, 2, 7, 9).

## Goal
Build the RCR M2L coaxial tables, and store their y-blocks, only on the levels whose M2L reads
them (an RCR pair under the setting), deciding every level's M2L, the diagonal form's classes
included, before the RCR tables, with every output unchanged bit for bit.

## Acceptance (the roadmap's *Done when*)
- The operator's and the FMM's outputs bit for bit those with RCR M2L tables on every level from 2
  (the build before this task): in the tests (the pair bodies of every setting; the FMM at 2
  threads and on P ranks against the one-rank FMM); on Kathleen, the checksums `main`'s.
- On every level with an RCR pair the tables hold the RCR M2L, and on no other (the levels of the
  waves, Batched, Dense, Classes, and the diagonal form's levels that take all 34 classes).
- The build of W1 at 10⁶ points, kD 250, 500, 1000, and of W3 at 10⁶ points at kD = 500, under the library's
  defaults against `main`'s on Kathleen: the RCR tables' time, the build, the tables' bytes.

## Read first
- `tetraleaf-exec/src/helmholtz/operator.rs`: `select_tables` (the order: waves, Dense, Classes,
  Batched, the RCR tables, then `diagonal::select`), `Kernel::m2l` (the pair dispatch: waves, a
  diagonal offset, Dense, Classes, Batched, else RCR), `size_rcr_buffer`.
- `tetraleaf-exec/src/helmholtz/operator_diagonal.rs`: `select`, `DiagonalLevel`, the level call's
  RCR pairs (`kernel.tables.m2l` for offsets without `has_offset`).
- `tetraleaf-tables/src/helmholtz.rs`: `TableChoice::{m2l_levels, block_levels}`, `has_m2l`,
  `block_selection` (blocks only on levels with M2L tables).
- Report `tetraleaf-validate/results/helmholtz-ch2-2a.md` §0, §2, §5 (the finding).

## Changes
- `tetraleaf-exec/src/helmholtz/operator.rs`, `select_tables`: the diagonal form's levels chosen
  before the RCR tables (it reads only k̂_l, a function of k̂₀ and the level, as `HelmholtzTables`
  forms it); then `read` per level: l ≥ 2, the level's M2L RCR's (no waves, Dense, Classes or
  Batched tables) and, on a diagonal level, some V offset outside its classes. The RCR tables
  are built with `TableChoice::m2l_levels` and `block_levels` both `read`, and built again when a
  level's M2L tables differ from `read` (today: when a level from 2 has none). The timings keep
  their parts (`TableTimings::{rcr, diagonal}`), in the new order.
- A doc-hidden test hook, `HelmholtzOperator::with_m2l_tables_everywhere(bool)` and
  `HelmholtzFmmBuilder::m2l_tables_everywhere(bool)` (default `false`): RCR M2L tables on every
  level from 2 and blocks as before this task (the trusted path of the tests).
- Module docs ("M2L", "The table store"), the crate `CLAUDE.md` paragraph.
- `tools/helmholtz/ch2-2b1-runs.sh` (stages `smoke`, `w1`, `w3`; `main`'s binary from `MAIN`),
  `ch2-2b1-tables.py` (on `common.py`). Report `tetraleaf-validate/results/helmholtz-ch2-2b1.md`.
- No harness code (the build parts, "RCR tables by level" and the tables' MB are CH2.2a's); no
  convention; no dependency; `tetraleaf-plan` and `tetraleaf-tables` untouched.

## Tests that prove it
- `tetraleaf-exec/tests/helmholtz/build.rs`, `the_rcr_m2l_tables_are_built_where_read`: under
  Rotation, Auto at 10⁻⁶ and 5·10⁻⁸ (waves, Batched), PlaneWave, Batched, Dense, Classes, and the
  diagonal form (`Tolerance`, every class on some level; and a level with some classes), at k̂₀ =
  15 and a real k̂ where the form takes all 34 classes: `has_m2l(l)` exactly on the levels with an
  RCR pair, blocks only there; every pair body bit for bit the operator with
  `with_m2l_tables_everywhere(true)`; that operator's `has_m2l` on every level from 2; fewer table
  bytes where some level drops. Plus the existing `the_unbuilt_operator_builds_the_settings_tables`
  (the unbuilt operator's one build against the chain of setters, now under the new rule).
- `tetraleaf-exec/tests/helmholtz_fmm.rs`, `check_m2l_tables` (in `check_diagonal`, on its cube at
  kD = 40 + 2i and clusters at kD = 20: `HelmholtzDiagonal::Rule` under Rotation and under the
  library's M2L Auto, whose levels with waves or Batched lose their tables, and
  `Tolerance(0.1)` under Rotation, which takes some level whole): the levels' M2L tables as above
  on every rank, the output bit for bit the FMM with `m2l_tables_everywhere(true)`; 2 threads and
  the one-rank FMM bit for bit (the Rule runs by `check_diagonal`'s own checks, the whole-level run
  by its own).
- Every existing operator and FMM check runs on the new tables (a read of a missing table panics).

H6/H4 gate rules: no bound, degree, tolerance or rule changes (only which tables are built), so no
wider scope, no pair or FMM evidence to pair and no bound to gate.

## Checks
`/verify` (fmt; Clippy default, all-targets, workspace; `cargo test`, `--workspace` at opt-level 1;
docs; the MPI list at 2 and 4 ranks; exec's crate checks); `cargo deny` not needed (no manifest
change).

## Machine runs
Kathleen, `singlenode`, one node, 1 × 40, `--exclude=node-c11b-004,node-c11b-013`: `main`'s
`helmholtz_accuracy` (63a7a4f, built in `src/main`) and the branch's, one setting per launch, in the
order main, branch, branch, main per case: W1 at 10⁶ points (`w1:1000000`) at kD 250, 500, 1000 and
W3 at 10⁶ points (`w3:1000000`, CH4.3b2c's; its default N at kD = 500 is 5·10⁸) at kD = 500, ε = 10⁻⁶, the library's defaults named (CH2.2a's), `--levels`. Decides nothing: it
shows the RCR tables' time and bytes saved, and the checksums `main`'s.

## Constraints and risks
- Constraint 1: no Laplace path (the Helmholtz operator's table choice, its builder hook, tests).
- Constraint 2: the RCR M2L of a level is the same function of the same inputs wherever built, so
  the output is bit for bit; the choice per level reads the rank-agreed settings only.
- Constraint 7: fewer bytes; the block budget may now store blocks on coarser RCR levels (no bit
  changes: stored and per-call blocks agree, CH4.0e's `check_block_budget`).
- Constraint 9: the choice is the setting's per-level M2L, from k̂_l, the degrees and the costs.
- Risk: a reader of RCR M2L tables outside `Kernel::m2l`'s dispatch would panic on a level without
  them; the tests run every M2L path, and the hook's operator compares bit for bit.

## Pending runs
(filled in step 6)

## Open questions
None yet.
