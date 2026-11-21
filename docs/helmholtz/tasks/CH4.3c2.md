# CH4.3c2 The high-frequency degrees refit

Started 2026-11-21 on `tbetcke/h-ch4-3c2`. Roadmap: CH4.3c2 "The high-frequency degrees refit"
(the second of CH4.3c1–CH4.3c2, split on 2026-11-21; a *Need* of CH4.3b2c since CH4.3c1's
question 2). Proposal: §10 CH4.3 (the per-(level, class) rule, calibrated); §3.5 (truncation, the
calibrated rule, "CH4.3 recalibrates the high-frequency levels", "the excess p − x … grows about
linearly in L … which is CH4.3's to refit"); §9.1 (workloads, norms; ε/2, verify ≤ ε); §2
(constraints 1, 2, 9). Reviews: 2026-11-19 A8 ("levels 2–4 carry 0.6–1.0 % of Σ|e|² … CH4.3's
refit searches downward at x ≥ 54"); 2026-11-12 M7 (the high-frequency branch has the wrong
ε-dependence: one C over-predicts by 8 degrees at 10⁻³, leaves 2 at 10⁻⁸ at x = 65).

## Goal
Refit the calibrated rule's high-frequency branch in the full FMM, with the excess p − x linear
in L = log₁₀(1/ε) and the degrees searched downward at x ≥ 54, so that the high-frequency levels
are not over-resolved, and put the refit to the user as the new default.

## Acceptance (the roadmap's CH4.3c2, from CH4.3c's *Done when*)
"(review 2026-11-19, A8) the refit of the high-frequency degrees searches downward at x ≥ 54
(levels 2–4 carry 0.6–1.0 % of the error on W1 at kD = 1000, 10⁻⁸), the excess p − x refitted in
L (proposal §3.5; CH3.3's report §6 names the branch x + (a L + b) x^{1/3}); pair evidence and
FMM evidence, both; the rule for sign-off."

Checked against the tree on 2026-11-21: `DegreeRule::CALIBRATED` (`tetraleaf-exec/src/helmholtz/degree.rs`)
is still CH3.3's, max(⌈0.30 L² + 0.25 L + 2.60⌉, ⌈x + 0.80 (ln 1/ε)^{2/3} x^{1/3}⌉); CH4.0c's
§2.2 (W1 at kD = 1000, 10⁻⁸) is still the only share-by-level measurement beyond x = 65.

What the items mean here:
- **Searches downward at x ≥ 54.** Level scans (CH3.3's method: one level's degree, every other
  level held far below ε) at x = 54, 108 and 216 that start below the current rule's degree,
  so the measured degree is the smallest that meets ε/2, not "the rule's degree suffices".
- **Refitted in L.** A high-frequency branch whose excess is linear in L, fitted to cover every
  measured degree (CH3.3's 176 points with this task's); the plateau unchanged.
- **Pair evidence.** S7's per-pair worst case (`spikes/truncation/results/c1_curves.csv`, the
  one-level chain in extended precision, the worst of 316 offsets and three placements, φ and
  ∇φ, to k a = 256): the per-pair error at the FMM's measured degrees and at the refit's. The
  FMM is known to need less than the per-pair worst case (CH3.3: at x = 65 the FMM's degree
  leaves the worst pair at 2–5 ε), so the per-pair measure is used as a guard on the
  extrapolation, not as the criterion: at every grid k a from 45 to 256 and every ε of D1, the
  refit's per-pair worst error over ε must not exceed the largest such ratio at a degree the FMM
  confirmed at x ≥ 43 (CH3.3's and this task's points). A refit that passes the FMM but fails
  this is reported, and not proposed beyond the x the FMM measured.
- **FMM evidence.** The level scans above, and the refit as the default verified against ε
  (not ε/2) on CH3.1c's matrix and CH3.3's held-out cases, on W1 at kD = 500 and 1000 with the
  three densities (CH4.0b's sign-off: CH4.3 verifies on the smooth densities too), under Rotation
  and the library's default.
- **For sign-off.** The refit is made the default on this branch (`DegreeRule::CALIBRATED`),
  CH3.3's rule kept as a named constant for reproduction; the PR asks whether to adopt it.

## Read first
- `tetraleaf-exec/src/helmholtz/degree.rs` (the rule, its documentation of the measured range);
  `tetraleaf-exec/tests/helmholtz/degree.rs` (`MEASURED`, the S7 bounds, monotonicity).
- `tetraleaf-validate/src/helmholtz_calibration.rs` (`Curve`, `Point`, `fit_plateau`, `fit_high`,
  `fit`, `covers`); `tetraleaf-validate/examples/helmholtz_calibrate.rs` (scans, `--hold`, `fit`).
- `tetraleaf-validate/results/helmholtz-ch3-3.md` (§1 method, §3 level scans, §4 fit, §6);
  `helmholtz-ch4-0c.md` (§2.2 levels, §2.3 times at kD = 1000).
- `spikes/truncation/SPIKE_REPORT.md` §3 and `results/c1.csv`, `results/c1_curves.csv`.
- The two restatements of the rule in `tetraleaf-tables`' tests
  (`tests/helmholtz/plane_wave_pair.rs::calibrated_degree`,
  `tests/helmholtz/diagonal_transfer.rs::calibrated_degree`) and the module docs of
  `tetraleaf_tables::helmholtz::plane_wave::merge` that quote it.

## Changes
- `tetraleaf-exec` (`helmholtz::degree`): `DegreeRule` gains the linear-in-L branch. Shape (public
  API of `tetraleaf-exec`, not of `tetraleaf-plan`): fields `gamma, alpha, beta` (plateau, unchanged)
  and the high-frequency branch's coefficients, CH3.3's C (ln 1/ε)^{2/3} kept expressible so that
  CH3.3's rule stays a constant (`DegreeRule::CH3_3`) and S7's closed form a `DegreeRule`;
  `CALIBRATED` the refit. The exact form is the fit's outcome (below); the candidates are
  (A) x + (a L + b) x^{1/3}, CH3.3's §6, and, if (A) covering the small-x points leaves the
  x ≥ 54 points over-resolved by more than 2 degrees, (B) one more constant term. Module docs:
  the measured range and the refit's record. Plateau and Laplace untouched.
- `tetraleaf-validate`: `helmholtz_calibration`'s fit of the new branch (grids over the branch's
  constants, least total excess over the points above the plateau, covering every point);
  `helmholtz_calibrate fit` printing CH3.3's rule and the refit against each point;
  `helmholtz_accuracy --rule calibrated|ch3-3|s7` (CH3.3's rule named, for the times and errors
  before). The measured points stay constants in the tests, as `MEASURED` is (no fixture file).
- `tetraleaf-tables`' two test restatements of the rule follow `CALIBRATED`.
- `tools/helmholtz/ch4-3c2-runs.sh` (the stages below), `tools/helmholtz/ch4-3c2-pairs.py` (the
  pair evidence from S7's CSVs), report `tetraleaf-validate/results/helmholtz-ch4-3c2.md`.
- Conventions: none. Dependencies: none. `tetraleaf-plan`: untouched.

## Tests that prove it
- `tetraleaf-exec/tests/helmholtz/degree.rs`: `MEASURED` extended with this task's points; the refit
  covers every measured degree (exact integers); at most S7's closed form and monotone in x and
  ε over the existing grids; `CH3_3` reproduces CH3.3's rule's degrees on the old grid (exact);
  the shapes (S7's closed form, the ceilings' 10⁻⁹) still expressible; `resolve` and the doctest's
  degrees at the refit.
- `fmm-validate`'s `helmholtz_calibration` unit tests: the new fit recovers a known rule from
  synthetic points, covers the points, and its least excess (as the existing `fit` tests do).
- The FMM tests that use the default rule (`tests/helmholtz_fmm.rs`, `helmholtz_accuracy`,
  `tetraleaf-tables`' helmholtz tests) pass on the refit unchanged in their bounds: they compare
  against the direct sum and Rotation at the same degrees, so a degree change moves no tolerance.
- Trusted paths: the direct sum (`tetraleaf_ref::helmholtz::p2p::direct_sum`) for every FMM error;
  S7's extended-precision chain for the pair measure.

## Checks
`/verify` (format, Clippy default and workspace all-targets, rustdoc `-D warnings`, the default
and workspace tests at opt-level 1, the MPI list at 2 and 4 ranks); `cargo test -p tetraleaf-exec
--test helmholtz degree`; `cargo test -p tetraleaf-validate`; `cargo test -p tetraleaf-tables --test
helmholtz` (the restatements).

## Machine runs
Kathleen, `singlenode`, one rank × 40 threads, BLAS threads 1, release; jobs through
`tools/kathleen/submit.sh run.sbatch`; raw outputs to
`~/projects/fmm-archive/kathleen/runs/ch4-3c2*/`. Accuracy runs; times reported, never asserted.

1. **Scans** (`helmholtz_calibrate --scan levels`, ε = 10⁻³ … 10⁻⁸, the others held at
   `--hold 0.29,0.45,1.65,1 --hold-eps 1e-10` as CH3.3, Rotation), each starting below the
   current rule's degree at 10⁻³ (the default start ⌊x⌋ − 2):
   - W4 cube (N = 10⁵) at kD = 250 (level 2, x = 54.1), 500 (levels 2–4: 108.3, 54.1, 27.1) and
     1000 (levels 2–4: 216.5, 108.3, 54.1);
   - W1 at kD = 250 (level 2, x = 54.1) and 500 (levels 2 and 3, x = 108.3 and 54.1);
   - W1 at kD = 1000 level 2 is not scanned (about 10 min per degree, CH4.0c §2.3): the cube's
     scan there and the verification below stand for it.
   Decides the measured degrees at x ≥ 54.
2. **Fit** (on the Mac, no MPI): `helmholtz_calibrate fit` over CH3.3's points and step 1's;
   the pair evidence from S7's curves (`ch4-3c2-pairs.py`). Decides the refit.
3. **Verification** of the refit as the default, against ε: CH3.1c's matrix and CH3.3's
   held-out cases (`helmholtz_accuracy`, Rotation and the default), with `--rule ch3-3` on the
   same cases for the errors and times before; W1 at kD = 500 and 1000 at 10⁻⁸, 10⁻⁷, 10⁻⁶, 10⁻³
   with the random charges, the plane wave's trace and the constant density (`helmholtz_boxes`
   as CH4.0c, the samples only), Rotation and the default. Decides the acceptance ≤ ε.

## Constraints and risks
- Constraint 1 (Laplace): no Laplace path or default changes; the rule is the Helmholtz FMM's
  alone. Constraint 9: the degrees still depend on k, the domain and ε only.
- Every Helmholtz default downstream reads the degrees (Auto's M2L by degree, the waves' p_w =
  max(p_l, 28/32), the diagonal form's classes, the V rule's crossovers by p): changes on the
  high-frequency levels move these choices. The plateau is unchanged, so the low-frequency
  levels with waves keep their degrees; levels at x ≈ 5–30, where the high branch rises above
  the plateau, can move. The verification runs the default too.
- Risk: the scans at small x (5–14) bind a single (a L + b) x^{1/3} at 10⁻³ (x = 5.4 needs a
  coefficient 2.05 at L = 3, x = 65 only 0.76), so form (A) may gain little at large x at 10⁻³;
  hence candidate (B). Risk: W1 at kD = 1000 measured only by verification.
- Kathleen: single-node jobs within `singlenode`'s 6 h, as CH3.3 and CH4.0c ran.

## Pending runs
(filled in step 6)

## Open questions
(filled when they arise)
