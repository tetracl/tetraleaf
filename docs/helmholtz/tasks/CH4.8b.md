# CH4.8b The constant density on W2: the remedy as the default

Started 2026-11-25 on `tbetcke/h-ch4-8b`. Roadmap: "CH4.8b The constant density on W2: the remedy
as the default" (H4; priority list item 11; its *Needs*, CH4.8a, merged as #260). Proposal: §10
CH4.4 (the row whose bound CH4.8 restores); §3.5 (the degree rule, amended here); §9.1 (ε the
relative ℓ2 tolerance; the worst box at most 100 ε within D1); §2 (constraints 1, 2, 6, 9). Builds
on CH4.8a (`tetraleaf-validate/results/helmholtz-ch4-8a.md`, its §8 signed off on #260: "All
decisions as recommended."), CH4.4a (the random-charge sweep, `results/helmholtz-ch4-4a.md`) and
CH4.3c3 (the rule it amends, `results/helmholtz-ch4-3c3.md`).

## Goal
Make CH4.8a's signed-off remedy `h1` (one degree added to the high-frequency branch after its cap,
the plateau unchanged) the library's default rule, and verify it on CH4.8's whole set with its
cost, for sign-off.

## Acceptance (proposal §10, verbatim)
§10 CH4.4's row, whose bound CH4.8 restores: "FMM against the direct sum (sampled targets) on W1 to
kD = 1000, W2 to 500 and W3 to 200, beside `hfmm3d` against the same sum; the full-target
crosscheck mapped to levels and classes, with no error concentrated at a class boundary; time,
memory and build cost per k". With §9.1: ε the relative ℓ2 tolerance, the worst box at most 100 ε
within D1 (kD ≤ 500).

The roadmap's *Done when* of CH4.8b: CH4.8a's signed-off remedy `h1` (one degree added after the
cap of the high-frequency branch, as a margin of `DegreeRule`; the plateau unchanged) the library's
default (proposal §3.5 amended); under the defaults W1–W3 with the constant density at kD = 10–500
within ε in ℓ2 at 10⁻³–10⁻⁸ and every worst box ≤ 100 ε, the random-charge cases of CH4.4a
unchanged within ε; its cost (time, memory, build per k, as §10's CH4.4) on Kathleen and locust;
for sign-off.

## Read first
- `tetraleaf-validate/results/helmholtz-ch4-8a.md` §1 (`h1` as run: d = 3 without the cap), §5,
  §7, §8 (question 2: "+1 after the cap", the cap kept).
- `tetraleaf-exec/src/helmholtz/degree.rs` (`DegreeRule`, `CALIBRATED`, `high`),
  `tetraleaf-exec/tests/helmholtz/degree.rs` (`MEASURED`, the cap test, `resolve`).
- `tetraleaf-exec/CLAUDE.md`, "The degree rule": a change of the constants is a recalibration
  (rerun the verification, update `MEASURED`'s test).
- `tools/helmholtz/ch4-8a-runs.sh`, `ch4-4a-runs.sh`, `ch4-5-runs.sh` (the launches reused).

## Changes
- `tetraleaf-exec` (Helmholtz only, `src/helmholtz/degree.rs`): `DegreeRule` gains a public field
  `margin: usize`, degrees added to the high-frequency branch after its cap and ceiling,
  ⌈min(branch, cap)⌉ + margin (0 in `gd`, `CH3_3`, `CH4_3C2`); `DegreeRule::CH4_3C3`, the default
  until now (margin 0), kept for the runs before and after; `DegreeRule::CALIBRATED` = `CH4_3C3`
  with margin 1, the builder's default through `Degrees::resolve` and `degree`. The module docs,
  the formula and the doctest's degrees ([0, 0, 34, 21] → [0, 0, 35, 22] at kD = 100, 10⁻⁶).
  `tetraleaf-exec/CLAUDE.md`'s paragraph.
- `tetraleaf-validate`: `--rule ch4-3c3` in `helmholtz_accuracy`, `helmholtz_boxes` and (new there)
  `helmholtz_fmm3d compare` (`--rule calibrated|ch4-3c3`, `Case::degree_rule`), for the before and
  after; `helmholtz_boxes`' `--rule γ,…` constants keep margin 0, so CH4.8a's named rules repeat.
  The crate's `CLAUDE.md` line.
- Tests whose asserted degrees, levels or plans follow from the rule are updated where the new
  degrees change them, each with its reason stated; none weakened.
- `tools/helmholtz/ch4-8b-runs.sh` (stages below) and `tools/helmholtz/ch4-8b-tables.py` (the
  report's tables; `tools/helmholtz/ch4-4b-tables.py` and `ch4-8a-tables.py` reused where they
  read the outputs). Report `tetraleaf-validate/results/helmholtz-ch4-8b.md`.
- Proposal §3.5 amended ("amended at CH4.8b"), its Appendix A line; the roadmap.
- No convention change (the degrees are no convention of `docs/CONVENTIONS.md`); no dependency; no
  `tetraleaf-plan` interface change; no Laplace path.

## Tests that prove it
- `DegreeRule::CALIBRATED` = max(plateau, `CH4_3C3`'s capped branch + 1) on a grid of
  4 ≤ x ≤ 512 and 10⁻⁸ ≤ ε ≤ 10⁻³ (by sixteenths of an octave and tenths of a decade); the plateau
  unchanged; `CH4_3C3` at most its cap and its earlier properties (the cap test moved to it).
- The trusted path, CH4.8a's runs: the default's degrees per level equal `h1`'s as run (the report's
  tables §5, W1 and W2 at kD = 250 and 500, all six ε for W2, the levels the waves do not raise to
  p_w), and `CH4_3C3`'s the default's as run there.
- `h1` as run (d = 3, no cap) against the new default: equal on that grid but where the cap binds
  (x = 7.1–13.5, ε ≤ 2·10⁻⁷), where `h1` is one degree more: the points listed, as CH4.8a §1 states.
- `MEASURED` and `MEASURED_CONSTANT` still covered; S7's closed form still at or above the default
  (or the points where it is not, stated), monotone in x and ε.
- The FMM tests (`tests/helmholtz_fmm.rs`, the harness's `tests/helmholtz_accuracy.rs`) pass under
  the new default degrees against the direct sum, at 2 threads and on P ranks against the one-rank
  FMM bit for bit, as they stand.

## Checks
`/verify` (fmt; Clippy default, examples, workspace all-targets; rustdoc `-D warnings`; the default
and workspace tests at opt-level 1; the MPI list at 2 and 4 ranks for `tetraleaf-exec`'s
`helmholtz_fmm`), `bash -n`, ShellCheck and dry runs of the runs script, `py_compile` of the tables
script, `tools/roadmap/sync_issues.py --out`.

## Machine runs
Kathleen `singlenode`, one rank × 40 threads, release, BLAS threads 1, `--exclude=node-c11b-004,
node-c11b-013`, raw outputs to `~/projects/fmm-archive/kathleen/runs/ch4-8b/`; locust one rank × 72
on a quiet node (`nvidia-smi`, `uptime` first), outputs to `~/projects/fmm-archive/locust/runs/ch4-8b/`.
Stages of `tools/helmholtz/ch4-8b-runs.sh`:
1. `const` (the acceptance's constant density; the third rule: the worst density over the whole
   range): `helmholtz_boxes --m2l defaults --study single --refine 8 --refine-faces --worst 8` under
   the new default at all six ε of 10⁻³–10⁻⁸ (CH4.8a's lesson: the gate's three ε missed 10⁻⁴ and
   10⁻⁷), on W1 and W2 at kD = 10, 30, 100, 175, 250, 500, W2 also at 375 (between the records'
   points: CH4.8a ran no W2 between 100 and 250), W3 resolved at kD = 10, 30, 100 and W3 at N = 10⁶
   for kD = 200 and 500; the sample cached (CH4.8a's and the H4 gate's caches reused). Decides the
   constant-density acceptance (ℓ2 ≤ ε, worst box ≤ 100 ε).
2. `const-before`: the same matrix under `--rule ch4-3c3` (the default until now), so that every row
   has its before; the cached samples make the rows the records ran equal to them (a check).
3. `random`: CH4.4a's random-charge and plane-wave cases (W1 at kD = 10–1000, W2 at 10–500, W3 at
   10, 30, 100 and at N = 10⁶ for 200 and 500, W1 at 500 and 1000 with the plane wave's trace) under
   the new default, `helmholtz_fmm3d compare --settings defaults` without FMM3D, all six ε in one
   launch per case (accuracy may share a launch). Decides "the random-charge cases unchanged within
   ε".
4. `cost` (Kathleen and locust): `compare --settings defaults --rule calibrated` (B) against
   `--rule ch4-3c3` (A), one ε per launch, A B B A, on W1 at kD = 100, 250, 500, 1000, W2 at 500 and
   W3 at N = 10⁶, kD = 500, at 10⁻³, 10⁻⁶, 10⁻⁸: the build, the evaluation, the memory and the build
   per k, as medians with their spread (one setting per launch, CH4.3b2c §6).
Estimated: `const` and `const-before` about 3–4 node-hours (W3's references the dearest), `random`
about 1–2, `cost` about 2 per machine. Jobs of at most 6 h.

## Constraints and risks
- Constraint 1: no Laplace path; `tetraleaf-exec`'s Helmholtz degree rule and the harness only.
- Constraint 2 (the degree rule changes on both measures): CH4.8a's pair and FMM evidence chose
  `h1`; this task's FMM runs over the whole set confirm it, and the pair table of CH4.8a §6 is `h1`'s.
- Constraint 9: the rule reads x and ε only.
- The verification set widens (the H6 gate's rule): every ε of D1 on W1–W3 with the constant
  density and the random charges, W2 between the records' kD, W1 to kD = 1000 for the random charges.
- Risk: the degrees change every Helmholtz test's FMM; asserted plans or per-level choices may move
  (a level's M2L by the time rule, a crossover). Each such change is stated with its reason.
- Risk: the cost at small kD, where the high-frequency levels are few, is about zero; at kD = 500
  CH4.8a measured 1.017× (geometric mean). Memory grows with (p + 1)² on the high-frequency levels.

## Pending runs
(filled when submitted)

## Open questions
None yet.
