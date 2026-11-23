# CH4.5 Host performance against FMM3D at matched accuracy

Started 2026-11-23 on `tbetcke/h-ch4-5`. Roadmap: "CH4.5 Host performance against FMM3D at matched
accuracy" (H4; priority list item 10, the last of H4's tasks; its *Needs* CH4.4 ticked on `main`
with CH4.4a and CH4.4b). Proposal: §10 CH4.5; §9.2 ("What is compared": time and memory at matched
accuracy, as C6.8, setup apart; "Limits": FMM3D's eps chosen by its error, amended at H6.G); §9.1
(the workloads and norms); §2 (constraints 1, 2, 6, 9). C6.8: docs/phase6/README.md, its method in
`tools/fmm3d/README.md` and `tetraleaf-validate/results/phase6-fmm3d.md`. Builds on CH4.0e (FMM3D
under the defaults at matched accuracy, CH3.6's matrix, both machines) and CH4.4a (the sweep, FMM3D
beside it at the same eps).

**CH4.0k reconsidered before it** (the roadmap's Also, review 2026-11-19, A6): asked on 2026-11-23,
the user: "Keep off; CH4.5 shows it (Recommended)". CH4.0k stays off the list; this task reports our
build by part at one rank per node beside the best split, the shared steps (octree, routes, plan)
named, so that the H4 gate sees their weight.

## Goal
Time our host FMM under the library's defaults against FMM3D's `hfmm3d` on CH4.4's sweep, each at
the accuracy of the other, on one node of Kathleen and on locust, with the setup apart from the
evaluation and the memory of both.

## Acceptance (proposal §10, verbatim)
"as C6.8: for each of our tolerances the `eps` of FMM3D with the closest error; time and memory,
setup apart, one node, on the sweep of CH4.4".

Roadmap Alsos: (review 2026-11-19, A6, B6, B8) ours at 1 × N (one rank, N threads, as a threaded
caller runs it) beside the best split; FMM3D's OpenMP binding stated; CH4.0k reconsidered before it
(above); the H4 gate states one-rank-per-node figures (this report supplies them). (CH4.3e's report
§4) c_leaf 4, measured at one rank per node only, at the best split and at N × 1 (W3 at N = 10⁶,
kD = 500). (CH4.4a's questions 3 and 4) at kD = 1000 FMM3D `fastker` alone, its failure and peaks
stated; our memory one ε per launch; timed jobs exclude node-c11b-004 and node-c11b-013.

How it is read (to be confirmed in the PR):
- **The sweep**: CH4.4a's 26 cases (W1 at kD = 10–1000 and W2 at 10–500 with §9.1's random charges;
  W3 resolved at kD = 10, 30, 100 and at N = 10⁶ for kD = 200, 500; the plane wave's trace on W1 at
  kD = 500, 1000; the constant density on W1 and W2 at kD = 10, 30, 100 and W3 at 10, 30), at our
  ε = 10⁻³, 10⁻⁶, 10⁻⁸. One complex charge vector, potentials and gradients, §9.1's 1,000 samples.
- **Ours**: the library's defaults (`helmholtz_fmm3d compare --settings defaults`, `Case::defaults`,
  held to `HelmholtzFmmBuilder::new` bit for bit by `check_defaults`), one ε per launch (memory per
  launch, CH4.4a's question 4), one build and two evaluations: "setup" is the build, "evaluation"
  the second evaluation (the first reported too).
- **FMM3D's eps "with the closest error"**: over a grid of eps (3·10⁻³, 10⁻³, 3·10⁻⁴, …, 10⁻⁹, two
  per decade), per case, ε and variant, the eps minimising d = max over φ and ∇φ of |log₁₀(FMM3D's
  error / ours)|, the faster eps on a tie (C6.8 took φ alone; with two fields the worse of the two).
  Beside it CH4.0e's one-sided match (the fastest eps whose φ and ∇φ errors are both at most ours),
  which never gives FMM3D a looser result than ours. Where ours is round-off (the direct sums of the
  smallest trees) neither exists within the grid; the most accurate eps is shown and marked.
- **Time, setup apart** (C6.8's measure): ours the build and the evaluation apart, the slowest rank's;
  FMM3D keeps nothing between calls, so its call is timed whole and its setup (the tree:
  `hndiv`, `pts_tree_mem`, `pts_tree_build`, `pts_tree_sort` as `hfmm3d` calls them) timed apart by
  the driver, as T9 timed Laplace's. Ratios: ours' build plus one evaluation over FMM3D's call (one
  solve), and ours' evaluation over FMM3D's call less its tree (repeated solves).
- **One node, two splits of ours**: per case and ε the best split (the smallest build plus one
  evaluation over 1 × 40, 2 × 20, 4 × 10, 8 × 5 on Kathleen and 1 × 72, 2 × 36, 4 × 18, 8 × 9 on
  locust, each in two launches, the orders reversed: CH4.0e's rule) and 1 × N (one rank, every
  core), both against FMM3D on every core.
- **FMM3D's binding** (B8): OpenMP at its defaults (no `OMP_PROC_BIND`, `OMP_PLACES`: libgomp
  leaves the threads unbound), stated in every log by the driver; bound (`OMP_PLACES=cores
  OMP_PROC_BIND=close`) beside it on five large cases (stage `bind`), so the report says whether
  binding would change the comparison.
- **Memory**: ours the FMM's resident growth over its build and evaluations (`Row::memory`, the sum
  over the node's ranks and the largest rank's), one ε per launch; FMM3D's the driver's peak RSS per
  eps (`tools/fmm3d/maxrss.py`, one launch per eps, its input and output arrays included).
- **kD = 1000**: FMM3D `fastker` alone on Kathleen (`default` exceeds the node's 180 GB, CH4.4a). On
  locust (572 GB, shared) `default` too where the node has 250 GB available at the launch, else
  stated; locust's `fastker` is the slow variant there (CH4.0e §2.2), so its rows flatter ours and
  FMM3D's faster variant per case is the comparison's.
- **c_leaf** (CH4.3e's Also): W3 at N = 10⁶, kD = 500, c_leaf 2 and 4 at every split of the best-split
  set and at N × 1 (40 × 1, 72 × 1), both orders: c_leaf 4 over 2 per split.

## Read first
- Proposal §9.2, §10 CH4.5; docs/phase6/README.md C6.8; `tools/fmm3d/README.md` ("The method").
- `tetraleaf-validate/results/helmholtz-ch4-0e.md` §1–§2 (the matched rule and the splits),
  `helmholtz-ch4-4a.md` §1, §4–§6, `helmholtz-ch4-3e.md` §4, `helmholtz-ch4-0h.md` §3 (the shared
  steps at one rank per node), `phase6-fmm3d.md` §3–§4.
- `tetraleaf-validate/examples/helmholtz_fmm3d.rs` (`compare`, `errors`), `src/helmholtz_fmm3d.rs`
  (`seconds`, `maxrss`), `tools/fmm3d/driver.f90` (`helmholtz` mode), `tools/fmm3d/run.sh`
  (`helmholtz`), `tools/helmholtz/ch4-0e-runs.sh`, `ch4-0e-tables.py`, `ch4-4a-runs.sh`.

## Changes
`tetraleaf-validate` and `tools/` only (development tooling; no library crate, table, layout,
convention, dependency or manifest change):
1. `tools/fmm3d/driver.f90`, `helmholtz` mode: FMM3D's setup timed once per eps token, apart from
   the calls (`hndiv` of the linked variant with `hfmm3d_s_c_g`'s flags, then `pts_tree_mem`,
   `pts_tree_build`, `pts_tree_sort` as `hfmm3d` calls them); a field `tree_s=` at the end of the
   `helmholtz:` line (after `out=`, so `seconds`' parser and earlier logs are unaffected); and the
   OpenMP binding printed (`omp_get_proc_bind`, `OMP_PLACES`) on the `helmholtz: problem` line.
2. `src/helmholtz_fmm3d.rs`: `tree_seconds(log, token)` (the `tree_s=` field; `None` on older logs),
   unit-tested beside `seconds`.
3. `examples/helmholtz_fmm3d.rs`: `compare --evaluations n` (default 1, CH3.6's runs repeat; the row
   the last evaluation's, the first printed in the log) and `--c-leaf c` (`Case::c_leaf`, the
   library's when absent); `errors` prints FMM3D's tree seconds in a column.
4. `tools/helmholtz/ch4-5-runs.sh` (stages below; ours through `compare --settings defaults` with no
   `--fmm3d`, FMM3D through the driver one eps per launch under `maxrss.py`, `errors` over the grid)
   and `tools/helmholtz/ch4-5-tables.py` (the grid's errors, the matches by both rules, the eps list
   of the timed stage, and the report's tables).
5. `tools/fmm3d/README.md` (the driver's `tree_s` and binding line), `tools/kathleen/README.md` if
   a job needs a note, `tetraleaf-validate/CLAUDE.md`; the report
   `tetraleaf-validate/results/helmholtz-ch4-5.md`.

## Tests that prove it
- The acceptance is the runs; every error against the sampled f64 direct sum
  (`helmholtz_accuracy::Oracle`), ours and FMM3D's against the same oracle on byte-identical problem
  files (`compare --problem`, `errors --problem`).
- `tree_seconds` and `seconds` on a log line with and without `tree_s=` (unit tests in
  `src/helmholtz_fmm3d.rs`).
- `compare --evaluations 2` and `--c-leaf`: `tests/helmholtz_accuracy.rs` already holds
  `Case::evaluations` (two evaluations, the second's output the first's) and `Case::c_leaf`; the
  example passes them to `Case` only. A smoke run of `compare` on W1 at kD = 30 (1 and 2 ranks, the
  M6) with `--evaluations 2 --c-leaf 2` and without: the checksum with two evaluations equals one's,
  and `--c-leaf` changes the tree as printed.
- The driver's `tree_s` against FMM3D's own: a smoke run on the M6 (W1 at kD = 30) under both
  variants; `tree_s` below the call's seconds, the outputs unchanged bit for bit against a launch of
  the driver before the change (`cmp` of the `.out` files).
- `ch4-5-tables.py`'s matching on CH4.0e's three hand-checked rows (its one-sided rule must give
  CH4.0e's eps there from CH4.0e's archived grid), and `bash -n` on the script.

## Checks
`/verify` (format, Clippy default, examples and workspace all-targets, the default and workspace
tests at opt-level 1, `cargo doc`), `cargo clippy -p tetraleaf-validate --all-targets -- -D warnings`,
`cargo test -p tetraleaf-validate`; `bash -n` and `python3 -m py_compile` on the scripts; the driver
built by `tools/fmm3d/build.sh` on the M6, Kathleen and locust. No MPI crate's code beyond the
example; no Laplace file, so no D9.

## Machine runs
Release, BLAS threads 1. Kathleen: `singlenode` jobs (≤ 6 h), `--exclude=node-c11b-004,node-c11b-013`,
`run.sbatch 'tools/helmholtz/ch4-5-runs.sh <logs>/ch4-5 <stage…>'`; ours bound by socket
(`ppr:<r/2>:socket:PE=<t>`, 1 × 40 by `slot:PE=40`, 40 × 1 by `ppr:20:socket:PE=1`), FMM3D on 40 OpenMP
threads. locust: `tools/gh200/remote.sh`, the load and `free -g` recorded before and after every
launch, ours `slot:PE=<t>`, FMM3D on 72 threads. Stages:
1. `grid` (Kathleen; decides the matches): FMM3D over the grid at REPEAT 1, both variants
   (`fastker` alone at kD = 1000; `default` at 3·10⁻³ only where it fits), one eps per launch under
   `maxrss.py`, then `errors`. Estimated 7–9 h of node time (kD = 1000: 14 eps × about 350 s per
   density), in three or four jobs.
2. `ours` (Kathleen and locust): every case at each ε alone, every split in order (launch a) and
   reversed (b), two evaluations. About 2–3 h per machine.
3. `cleaf` (Kathleen and locust): W3 at N = 10⁶, kD = 500, c_leaf 2 and 4, every split and N × 1, a
   and b.
4. `time` (Kathleen and locust; after `grid` and both `ours`): FMM3D at each case's eps from
   `ch4-5-tables.py eps` (both matches on either machine, and ε itself), each variant, OpenMP at
   its defaults, REPEAT 2 below N = 10⁶ and 1 from there (calls of seconds to minutes, the setup
   in each), launches a and b in opposite orders, one eps per launch under `maxrss.py`. On locust
   the matches are made with locust's own errors of ours and the grid's of FMM3D; locust's
   FMM3D errors at those eps are measured and any change of match reported (CH4.0e found
   FMM3D's errors equal on the two machines to the digits shown). On locust also `default` at
   kD = 1000 at `fastker`'s eps, where 250 GB are available.
5. `bind` (Kathleen and locust; B8): FMM3D bound (`OMP_PLACES=cores OMP_PROC_BIND=close`)
   beside unbound on W1 at kD = 250, 500, 1000, W3 at kD = 100 and at N = 10⁶, kD = 500, at the
   closest eps of 10⁻⁶, each variant, launches a and b. Every core-count comparison in the
   report is at FMM3D's defaults (unbound), as CH4.0e and T9 ran it; `bind` says whether binding
   changes that. (Planned first: every timed eps under both bindings with REPEAT 2; about a day
   of node time at kD = 1000, so trimmed to this before the stage ran.)
Raw outputs to `~/projects/fmm-archive/{kathleen,locust}/runs/ch4-5/`. No run on the M6 (smoke
runs of the tools only, small, not timed).

## Constraints and risks
- Constraint 1 (Laplace): untouched; no library code changes. Constraints 2 and 9 (determinism over
  ranks and threads): each case's checksum is compared across splits and launches.
- Risk: node time. kD = 1000's grid dominates; if a job runs out, the grid is trimmed there to the
  eps around our errors (3·10⁻³–10⁻³, 3·10⁻⁷–10⁻⁸) and the trim stated.
- Risk: locust is shared; a launch under load is repeated, and the load is in the report.
- Risk: FMM3D `default` at loose eps on large trees fails to allocate (CH4.0e: W1 at kD = 500 at
  10⁻²); a failure is a row, not a stop.
- The matched rule moves in FMM3D's steps (its iprec and order tables): several eps give the same
  error, the match as close as the grid allows (C6.8, CH4.0e).

## Pending runs

## Open questions
