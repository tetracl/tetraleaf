# CH4.4a The kD sweep beside FMM3D

Started 2026-11-22 on `tbetcke/h-ch4-4a`. Roadmap: "CH4.4 kD sweep, with the FMM3D gate" (H4;
priority list item 10, after CH4.3c3; its *Needs* CH4.3a–CH4.3e and CH1.4 are ticked on `main`),
split on 2026-11-22 by the size rule into CH4.4a (the sweep: errors, crosscheck, time, memory,
build) and CH4.4b (the crosscheck mapped to levels and classes). Proposal: §10 CH4.4; §9.1 (the
workloads, norms, the box measure and its 100 ε bound, amended at CH4.0b); §9.2 (FMM3D, the H4 row);
§7 ("Build cost": measured and reported per frequency in CH4.4); §2 (constraints 1, 2, 6, 9).
Builds on CH3.6 (the comparison's tooling, `tetraleaf-validate/results/helmholtz-ch3-6.md`), CH4.0b
(the densities), CH4.0e (FMM3D under the defaults) and CH4.3c3 (the refit; the constant density's
plateau levels left to CH4.4).

## Goal
Run the Helmholtz FMM under the library's defaults (the H4 configuration: the diagonal form's time
rule, the refit degrees, c_leaf 4) and FMM3D's `hfmm3d` on the kD sweep of §10 CH4.4, against the
same sampled direct sum, with the full-target crosscheck, and report time, memory and build cost per
k.

## Acceptance (proposal §10, verbatim)
"FMM against the direct sum (sampled targets) on W1 to kD = 1000, W2 to 500 and W3 to 200, beside
`hfmm3d` against the same sum; the full-target crosscheck mapped to levels and classes, with no error
concentrated at a class boundary; time, memory and build cost per k". Roadmap Also: the worst box at
most 100 ε (§9.1, amended at CH4.0b); the plane wave's trace at kD = 500 and 1000 beside §9.1's
charges (CH4.0b's sign-off 2); the constant density on the sweep's low-frequency trees (CH4.3c3's
report §4).

CH4.4a takes everything but "mapped to levels and classes, with no error concentrated at a class
boundary", which is CH4.4b's. How it is read (to be confirmed in the PR):
- **Our configuration**: the library's defaults, every one (`HelmholtzFmmBuilder::new`'s), through a
  new `Case::defaults`; the harness's `Case::new` keeps its named settings (Rotation, the form off)
  so that the earlier reports' runs repeat. One complex charge vector (CH3.6's reading of §9.1's
  eight, kept), gradients on, 1,000 samples.
- **The sweep** (ε = 10⁻³, 10⁻⁶, 10⁻⁸, D1's range at CH3.6's three points):

  | Workload | kD | N |
  | --- | --- | --- |
  | W1, random charges | 10, 30, 100, 250, 500, 1000 | 796 … 7.96 · 10⁶ |
  | W2, random charges | 10, 30, 100, 250, 500 | about 1.6 N(W1) |
  | W3 resolved, random charges | 10, 30, 100 | 4,032 … 4.03 · 10⁶ |
  | W3 at N = 10⁶ (under-resolved), random charges | 200, 500 | 10⁶ |
  | W1, the plane wave's trace | 500, 1000 | |
  | W1, W2, the constant density | 10, 30, 100 | |
  | W3 resolved, the constant density | 10, 30 | |

  "W3 to 200" is read as §9.1 states W3: resolved where feasible (to kD = 100, 4 · 10⁶ points; at
  kD = 200 it would be 3.2 · 10⁷) and under-resolved at N = 10⁶ for kD = 200 and 500 (500 added: it
  costs little and §9.1 names it). The low-frequency trees of the constant density are kD ≤ 100,
  every tree whose levels are mostly plateau (x < 6.8); every tree's deepest levels are plateau
  levels, so the random-charge rows cover them too.
- **Accuracy**: ours within ε in §9.1's relative ℓ2 of φ and of ∇φ at the samples against the f64
  direct sum (`tetraleaf_ref::helmholtz::p2p::direct_sum` × 1/(4π)), every case and ε; FMM3D
  (`default` and `fastker`, its eps = ε) against the same sum, reported, not gated (§9.2: compared
  at matched accuracy in CH4.5).
- **The full-target crosscheck**: our output and FMM3D's at every point against FMM3D `default` at
  eps = 10⁻¹² (its own error at the samples printed): the relative ℓ2 over all targets and the box
  measure of §9.1 (per box of the deepest level, the largest error at a target over the reference's
  RMS); ours gated at the worst box ≤ 100 ε, FMM3D's reported.
- **Time per k**: one launch each on one Kathleen node; ours build and evaluation (the slowest
  rank), FMM3D's whole call (it keeps nothing between calls), as CH3.6.
- **Build cost per k**: our build by part (the harness's `BUILD_PARTS`: tree, plan, the tables by
  kind, …), the slowest rank's.
- **Memory per k**: ours the resident growth of the FMM over its build and evaluation (Linux: the
  peak `VmHWM` after resetting it through `/proc/self/clear_refs` right before the build, less
  `VmRSS` then), the largest rank's and the sum over the node's ranks; FMM3D's the driver's peak
  resident size per eps (`tools/fmm3d/maxrss.py`: getrusage's peak of the child, one driver launch per eps; the driver's input
  and output arrays included, 13 doubles a point, stated beside it).

## Read first
- Proposal §9.1, §9.2, §7 "Build cost"; `results/helmholtz-ch3-6.md` (method), `helmholtz-ch4-0e.md`
  (FMM3D under the defaults; FMM3D misses its eps at 10⁻³ on kD = 100), `helmholtz-ch4-3c3.md` §4.
- `tetraleaf-validate/examples/helmholtz_fmm3d.rs`, `src/helmholtz_fmm3d.rs`,
  `src/helmholtz_accuracy.rs` (`Case`, `run_output`, `BUILD_PARTS`), `src/scaling.rs`
  (`resident_memory`).
- `tools/fmm3d/run.sh` (the `helmholtz` step), `tools/kathleen/jobs/fmm3d.sbatch`.

## Changes
`tetraleaf-validate` and `tools/` only (development tooling; no library crate, table, layout,
convention, dependency or manifest change):
1. `helmholtz_accuracy`: `Case::defaults` (every setting the library's default: `HelmholtzM2l::Auto`,
   `HelmholtzDiagonal::Auto` with `DEFAULT_DIAGONAL_BUILD`, the V, W and X rules of the cost model,
   the rest `None`, which is already the library's); `Row::memory` (the FMM's resident growth: per
   rank the bytes before the build and the peak after it, reduced to the largest growth and the sums
   over the ranks; `None` where no peak can be read or reset, macOS) and a function that reads it.
2. `examples/helmholtz_fmm3d.rs`: `--settings harness|defaults` (default `harness`, CH3.6's runs
   repeat) and `--density random|planewave|constant` in `problem` and `compare` (the problem file
   carries the charges); `compare` prints our memory and FMM3D's peak RSS (from the driver's log,
   a line `maxrss: token=<eps> kb=<n>`), a table "Build" (a row per ε: the build, its parts, the
   evaluation) and in the log each level's degree, k a, M2L and the form's classes.
3. `tools/fmm3d/run.sh` (`helmholtz`): the driver launched once per eps under a new
   `tools/fmm3d/maxrss.py` (getrusage's peak of the child; Kathleen's compute nodes have no GNU
   `time`), appending the `maxrss:` line to the log; a case
   `<workload>@<kD>[/<density>]`; `H_SETTINGS` (passed as `--settings`); the README's
   `helmholtz` section.
4. `tools/helmholtz/ch4-4a-runs.sh` (the stages below, each a call of `run.sh helmholtz` with its
   cases and FMM3D's variants, `VARIANTS`) and `ch4-4a-tables.py` (the report's tables from the
   cases' Markdown).
5. `tests/helmholtz_accuracy.rs`: `check_defaults` (below).
6. `tetraleaf-validate/CLAUDE.md`; the report `tetraleaf-validate/results/helmholtz-ch4-4a.md`.

## Tests that prove it
- The acceptance is the runs, each against the sampled f64 direct sum (and the crosscheck against
  FMM3D at 10⁻¹², whose own error at the samples is printed).
- `check_defaults` (the MPI-owning test, one rank in CI, 2 and 4 ranks by hand; W1 at kD = 30,
  10⁻⁶): the harness under `Case::defaults` against an FMM built with `HelmholtzFmmBuilder::new`
  and nothing set but charges, gradients and threads: the checksum bit for bit (the trusted path
  for "the library's defaults"); `Row::memory` present on Linux with peak ≥ before, absent on macOS.
- Unit tests: the `maxrss:` line parsed (and a log without one); the densities' problem files
  differ in the charges only.
- The FMM3D side of `compare` is CH3.6's, covered by `check_fmm3d_measures`.

## Checks
`/verify` (format, Clippy default, examples and workspace all-targets, the default and workspace
tests at opt-level 1, `cargo doc`, the MPI list at 2 and 4 ranks), with
`cargo clippy -p tetraleaf-validate --all-targets -- -D warnings` and `cargo test -p
tetraleaf-validate`; `bash -n` on the scripts (`shellcheck` is not installed on the M6). No Laplace file changes, so no D9.

## Machine runs
Kathleen, `singlenode` jobs (1 × 40, ≤ 6 h each), release, BLAS threads 1, each
`run.sbatch 'tools/helmholtz/ch4-4a-runs.sh <logs>/ch4-4a <stage…>'` (`H_SETTINGS=defaults`); ours at 4 × 10 ranks × threads bound by socket
(CH3.6's split, T9's), FMM3D on 40 OpenMP threads, `default` and `fastker` at ε, `default` also at
10⁻¹²; the density cases with `default` only. Stages of `ch4-4a-runs.sh`:
- `smoke`: W1 at kD = 10 and W2 at kD = 10 with the constant density (the script, the memory lines);
- `low`: W1 and W2 at kD = 10, 30, 100, 250; W3 at kD = 10, 30;
- `const`: the constant density's cases;
- `mid`: W1 and W2 at kD = 500, W1 at kD = 500 with the plane wave's trace;
- `w1-1000`, `w1-1000-pw`: W1 at kD = 1000, random charges and the plane wave's trace;
- `w3`: W3 resolved at kD = 100, W3 at N = 10⁶ for kD = 200, 500.
Decides the acceptance. Raw outputs to `~/projects/fmm-archive/kathleen/runs/ch4-4a/`.

## Constraints and risks
- Constraint 1: no Laplace path touched (the harness, scripts, docs).
- Constraint 2: unchanged library code; the checksums printed.
- Constraint 6: every accuracy number against the f64 direct sum; FMM3D at 10⁻¹² only as the
  crosscheck's reference.
- Risks: FMM3D at 10⁻¹² on W1 at kD = 1000 (8 · 10⁶ points) may take an hour or exceed a job:
  W1 at kD = 1000 has jobs of its own; fallback the reference at 10⁻¹⁰. The full outputs (512 MB
  each at kD = 1000) are read on every rank: about 10 GB per case on 4 ranks, within 188 GB.
  A worst box above 100 ε is a finding, reported as such, not a test change.

## Pending runs
(filled in step 6)

## Open questions
(filled when they arise)
