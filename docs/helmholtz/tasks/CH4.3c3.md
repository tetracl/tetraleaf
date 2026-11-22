# CH4.3c3 The constant density in the calibration

Started 2026-11-22 on `tbetcke/h-ch4-3c3`. Roadmap: CH4.3c3 "The constant density in the
calibration" (added at CH4.3c2, its question 3: "q2 and Q3 as you recommend", the user on #249,
2026-11-21). Proposal: §10 CH4.3 (the per-(level, class) rule, calibrated); §3.5 ("Calibrated",
amended at CH4.3c2: "raised by 2 degrees because W1's constant density … missed ε"); §9.1
(workloads, norms; ε/2 in the scans, ≤ ε in the verification); §2 (constraints 1, 9).

## Goal
Measure, by level scans with W1's constant density, the margin that CH4.3c2 chose by hand (+2
degrees over the random-charge fit) and refit the high-frequency branch on the worst density if
the margin is wrong either way.

## Acceptance (the roadmap's CH4.3c3)
"(CH4.3c2's question 3) `helmholtz_calibrate` scans a chosen density (`--density`); W1's level
scans with the constant density at x = 13.5–216 (kD = 250, 500, 1000) measure the margin CH4.3c2
chose (+2 degrees over the random-charge fit); the rule refitted on the worst density if the
margin is wrong either way, verified as CH4.3c2, for sign-off."

Checked against the tree on 2026-11-22: `DegreeRule::CALIBRATED` (`tetraleaf-exec/src/helmholtz/degree.rs`)
is CH4.3c2's (C = 2.30, a = −0.75, b = −5.73, d = 5, cap 0.80); `helmholtz_calibrate` builds its
`Case` with `Density::Random` (no option); `helmholtz_accuracy::Case::density` exists (CH4.0b).

What the items mean here:
- **`--density`.** `helmholtz_calibrate --density random|planewave|constant` sets `Case::density`
  of the scans (complex charges for the plane wave, as `helmholtz_accuracy` refuses otherwise);
  the output's header names it, and a point of a density other than the random charges carries
  it in its workload field (`W1/constant`), so `fit` and `refit` read old and new files alike.
- **The scans.** CH3.3's and CH4.3c2's method (`--scan levels`, one level's degree from
  ⌊x⌋ − 2 upward, every other level held at `--hold 0.29,0.45,1.65,1 --hold-eps 1e-10`, Rotation,
  the smallest degree whose larger error of φ and ∇φ meets ε/2 at three degrees in a row, ε =
  10⁻³ … 10⁻⁸) on W1 with the constant density at every level of x from 13.5 to 216: kD = 250
  levels 2–4 (x = 54.1, 27.1, 13.5), kD = 500 levels 2–5 (108.3, 54.1, 27.1, 13.5), kD = 1000
  levels 2–6 (216.5, 108.3, 54.1, 27.1, 13.5): 12 scans, the largest degree at one x the
  measured one. (The plan first split kD = 1000's levels 2 and 3 into two launches by ε; a
  degree there costs about 110 s on a Kathleen node, not CH4.3c2's estimate of 10 min, so every
  scan is one launch.)
  The trees are the library's defaults of today (c_leaf 4 since CH4.3e; CH4.3c2 scanned at 2):
  the rule serves the default, and the scans' levels at x ≥ 13.5 hold the same boxes' far fields.
- **The margin.** Per scanned (x, ε): the constant density's degree p_c against CH4.3c2's fit
  (d = 3, capped: `DegreeRule::CALIBRATED` with d = 3) and against the rule (d = 5): the margin
  needed p_c − fit and the rule's spare rule − p_c. Beside them CH4.3c2's random-charge degrees
  at the same x (x = 27, 54, 108, 216), the difference the density makes.
- **Wrong either way.** Too small: some p_c above the rule. Too large: the refit on the worst
  density (below) lower than the rule by 2 degrees or more at some measured (x, ε) (one degree is
  the curve's resolution: the three-in-a-row criterion and the ceiling move it). Otherwise the
  margin stands and the rule does not change.
- **Refitted on the worst density.** `helmholtz_calibrate refit` over CH3.3's points, CH4.3c2's
  random-charge scans (its archived outputs) and this task's constant-density points, form C
  capped at CH3.3's rule on `cap_grid` (CH4.3c2's procedure), the excess clamped at 0, then CH3.3's
  branch as its `cap` (the default's construction); no hand-chosen offset, since the worst density
  is now among the points.
- **Verified as CH4.3c2.** For a changed rule: pair evidence (`ch4-3c2-pairs.py`, S7's worst pair
  at the refit's degree within the largest ratio at a degree the FMM confirmed, k a = 45 … 256)
  and FMM evidence (CH3.1c's matrix and CH3.3's held-out cases, Rotation and the defaults, and
  W1 at kD = 250, 500, 1000 with the three densities at 10⁻⁸, 10⁻⁷, 10⁻⁶, 10⁻³, all within ε).
  For the rule unchanged: W1 at kD = 250, 500, 1000 with the constant density at the same four ε
  under today's defaults (c_leaf 4, the time rule), Rotation and `helmholtz_boxes`' default,
  within ε (CH4.3c2 verified at c_leaf 2); its pair evidence stays CH4.3c2's.
- **For sign-off.** The finding (the margin measured), and a changed rule if any, as the PR's
  questions.

## Read first
- `tetraleaf-validate/examples/helmholtz_calibrate.rs` (scans, `point:` lines, `fit`, `refit`);
  `tetraleaf-validate/src/helmholtz_calibration.rs` (`fit_branch`, `cap_grid`, `CH3_3_LEVELS`).
- `tetraleaf-validate/results/helmholtz-ch4-3c2.md` (§2 the random-charge degrees, §3 the fit,
  §4 the constant density, §6 the verification); `docs/helmholtz/tasks/CH4.3c2.md`.
- `tetraleaf-exec/src/helmholtz/degree.rs`; `tetraleaf-exec/tests/helmholtz/degree.rs` (`MEASURED`).
- `tools/helmholtz/ch4-3c2-runs.sh`, `ch4-3c2-pairs.py`, `ch4-3c2-tables.py`.

## Changes
- `tetraleaf-validate`:
  - `helmholtz_calibration`: `point_workload(name, density)` (the workload field of a point line,
    `W1` or `W1/constant`) and `PointLine` with `parse` (the nine fields of a `point:` line, the
    density from the workload field), shared by the example's `fit` and `refit`; unit tests.
  - `examples/helmholtz_calibrate.rs`: `--density`; the header's "Density:" line and the scan
    heading's density; `refit` also prints CH4.3c2's fit (d = 3) and the rule (`CALIBRATED`) as
    rules, the refit with CH3.3's branch as its cap ("C ≤ CH3.3, capped"), and a table of the
    margin per point of a density other than the random charges (p, the fit's, the rule's, the
    refit's).
- `tetraleaf-exec`: only if the rule changes (`DegreeRule::CALIBRATED` and its docs, for sign-off).
  Either way `tests/helmholtz/degree.rs` gains the constant density's measured degrees
  (`MEASURED_CONSTANT`), which the default must cover (CH4.3c2's `MEASURED` is the random charges').
- `tools/helmholtz/ch4-3c3-runs.sh` (the stages below), `tools/helmholtz/ch4-3c3-tables.py` (the
  margin table from the scans' outputs and CH4.3c2's random-charge degrees); report
  `tetraleaf-validate/results/helmholtz-ch4-3c3.md`.
- Conventions: none. Dependencies: none. `tetraleaf-plan`: untouched. Laplace: untouched.

## Tests that prove it
- `helmholtz_calibration`'s unit tests: a point line written with each density parses back to
  its fields and density; an old line (no density) parses as the random charges; malformed lines
  are refused.
- `tetraleaf-exec/tests/helmholtz/degree.rs`: `the_calibrated_rule_covers_the_measured_degrees`
  also over `MEASURED_CONSTANT` (exact integers, from the scans); the existing monotonicity, S7
  and cap tests unchanged (and, for a changed rule, rerun on it).
- The density reaches the scans: the example's output header and the scans' errors (the
  constant density's curve differs from the random charges' on the same tree); `Case::density`
  itself is CH4.0b's, tested in `tests/helmholtz_accuracy.rs`.
- Trusted paths: the direct sum (`tetraleaf_ref::helmholtz::p2p::direct_sum`, the harness's
  `Oracle`) for every FMM error; S7's extended-precision chain for the pair measure.

## Checks
`/verify` (format, Clippy default and workspace all-targets, rustdoc `-D warnings`, the default
and workspace tests at opt-level 1, the MPI list at 2 and 4 ranks where an MPI crate's tests
change); `cargo test -p tetraleaf-validate --lib helmholtz_calibration`; `cargo test -p
tetraleaf-exec --test helmholtz degree`.

## Machine runs
Kathleen, one rank × 40 threads per launch, BLAS threads 1, release; `tools/kathleen/submit.sh`
with `run.sbatch` (`singlenode`, 6 h, two jobs at a time) or a 4-node `small` job running one
launch per node (`srun -N1 -n1` per launch, in parallel); raw outputs to
`~/projects/fmm-archive/kathleen/runs/ch4-3c3/`. Accuracy runs; times read, never asserted.

1. **Scans** (`ch4-3c3-runs.sh scans`, ONLY= per launch): the 12 scans above. Cost measured by
   the first jobs: about 20 s per degree at kD = 500 and 110 s at kD = 1000. Decides p_c.
2. **Margin and refit** (on the Mac): `ch4-3c3-tables.py margin`, `helmholtz_calibrate refit`
   over CH4.3c2's archived scans and step 1's. Decides whether the margin stands.
3. **Verification**: the rule unchanged, `ch4-3c3-runs.sh verify-w1` (W1, constant density,
   kD = 250, 500, 1000, today's defaults); changed, CH4.3c2's `verify-matrix` and `verify-w1` on
   the refit, and the pair evidence. Decides the acceptance ≤ ε.

## Constraints and risks
- Constraint 1 (Laplace): no Laplace path or default changes; the rule is the Helmholtz FMM's
  alone. Constraint 9: the degrees still depend on k, the domain and ε only; a density is a
  calibration input, never read at run time.
- A changed rule moves every downstream choice that reads the degrees (Auto's M2L by degree, p_w,
  the form's classes, the V rule): the verification runs the defaults too.
- Risk: kD = 1000's scans are long (the hold at 10⁻¹⁰ on every other level of an 8-million-point
  tree); hence the split launches and the measured cost before the jobs are packed.
- Risk: the constant density's curves may be less monotone than the random charges' (a coherent
  sum can cancel at one degree); the three-in-a-row criterion guards it, and the curves are in the
  report.

## Results
- The `--density` option, the point lines with the density and `refit`'s margin table, the cap
  search and CH4.3c2's rule beside the default (`tetraleaf-validate`); tested.
- Scans (jobs 663686–663689; 663696 added x = 6.8, below the plan's range, once the refit proved to
  lower the degrees there): the constant density needs 0–3 degrees over CH4.3c2's fit; CH4.3c2's
  rule one short at x = 216 (10⁻³, 10⁻⁷) and x = 13.5, 10⁻⁸ (kD = 250 and 500), 1–2 too many at
  x ≤ 54: wrong both ways.
- CH4.3c2's rule on today's trees with the constant density (job 663690): within ε, at most 0.94 ε.
- The refit: the plan's capped refit (CH3.3's branch) has no covering rule (x = 13.5, 10⁻⁸ needs 28,
  CH3.3's 27); `refit` raises the cap's coefficient to the least that admits one, 0.82: C = 2.00,
  a = −0.55, b = −4.66, d = 2 (a departure from the plan's construction, put to the user in the
  report's question 2 with the alternatives). Made the default on the branch; CH4.3c2's rule kept
  as `DegreeRule::CH4_3C2`.
- Pair evidence: at most 4.8 ε at k a = 45 … 256 (the guard 19–42 ε).
- FMM evidence (jobs 663699, 663700): the matrix's 174 runs within ε (at most 0.44 ε), W1 with the
  three densities within ε (at most 0.65 ε). Times (job 663701): 0.86–0.93× at 10⁻³, 0.98–1.02×
  at 10⁻⁶ and 10⁻⁸.

## Pending runs
- None.

## Open questions
- For sign-off (the report §8): the acceptance as read; the refit as the default (or CH4.3c2's
  rule kept, or d = 6 with the cap 0.82); proposal §3.5 amended.
