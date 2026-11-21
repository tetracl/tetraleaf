# CH4.3c2 The high-frequency degrees refit

Started 2026-11-21 on `tbetcke/h-ch4-3c2`; moved the same day onto the tetraleaf rename (the
commits replayed with the new crate names on `tbetcke/h-ch4-3c2-tetraleaf`; the pre-rename
branch is bundled in `~/projects/fmm-archive/ch4-3c2-pre-rename/`). Roadmap: CH4.3c2 "The high-frequency degrees refit"
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
  As built: the branch ⌈x + max(0, C (ln 1/ε)^{2/3} + a L + b) x^{1/3} + d⌉ (fields `c, a, b,
  d`; `DegreeRule::gd` for a = b = d = 0); the fit compares (A), (B) = (A) + d and (C) = (B) with
  C free, each also capped at CH3.3's rule on 4 ≤ x ≤ 512 at the decades of D1, since (B) and (C)
  uncapped raise the degree by one at x = 10–27, 10⁻⁸ (the costly lower levels). The fit is
  (C) capped: C = 2.30, a = −0.75, b = −5.73, d = 3. The excess is clamped at 0: the fitted e(ε)
  is negative above ε ≈ 2 · 10⁻³, outside D1. After the constant density's failure (below) the
  default is the fit raised by 2 (d = 5) and capped by CH3.3's branch everywhere
  (`DegreeRule::cap`, the high branch min(own, x + 0.80 (ln 1/ε)^{2/3} x^{1/3})): the user's
  suggestion on 2026-11-21 (raise level 2's degree slightly, at a small cost in time) applied to
  every level above the plateau, since the rule knows x, not levels.
- `tetraleaf-validate`: `helmholtz_calibration`'s fit of the new branch (grids over the branch's
  constants, least total excess over the points above the plateau, covering every point);
  `helmholtz_calibrate fit` printing CH3.3's rule and the refit against each point;
  `helmholtz_accuracy --rule calibrated|ch3-3|s7` (CH3.3's rule named, for the times and errors
  before). The measured points stay constants in the tests, as `MEASURED` is (no fixture file).
- `tetraleaf-tables`' two test restatements of the rule: used only by ignored studies that
  reproduce CH6.1f's, CH4.3a's and CH4.3c1's recorded tables, so they stay CH3.3's rule, their
  docs saying so (as built; the plan first had them follow `CALIBRATED`).
- `tools/helmholtz/ch4-3c2-runs.sh` (the stages below), `tools/helmholtz/ch4-3c2-pairs.py` (the
  pair evidence from S7's CSVs), report `tetraleaf-validate/results/helmholtz-ch4-3c2.md`.
- Conventions: none. Dependencies: none. `tetraleaf-plan`: untouched.

## Tests that prove it
- `tetraleaf-exec/tests/helmholtz/degree.rs`: `MEASURED` extended with this task's points; the refit
  covers every measured degree (exact integers); at most S7's closed form and monotone in x and
  ε over the existing grids; `CH3_3` reproduces CH3.3's rule's degrees on the old grid (exact);
  the shapes (S7's closed form, the ceilings' 10⁻⁹) still expressible; `resolve` and the doctest's
  degrees at the refit.
- `tetraleaf-validate`'s `helmholtz_calibration` unit tests: the new fit recovers a known rule from
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

## Results so far
- Scans (Kathleen jobs 653796 `scans-cube`, cancelled after its kD = 500 launch; 653797
  `scans-w1`; 653806 `scans-cube-1000`; 653817 `scans-cube-small`; the cube at N = 10⁵ costs 760 s per
  degree at kD = 1000, so x = 216 and kD = 1000's x = 108 were scanned on the cube at N = 10⁴,
  `scans-cube-small`; 653806, the N = 10⁵ cube at kD = 1000 level 3, ends at its time limit with
  part of its ε, and 653817's last launch, kD = 1000 level 3 on the N = 10⁴ cube, a cross-check
  of x = 108 that three workloads measured, was cancelled for the verification's slot). Measured
  at ε/2, 10⁻³ … 10⁻⁸: x = 54.1: 58, 62, 65, 68, 71, 74 (largest of four workloads); x = 108.3:
  112, 117, 121, 126, 129, 132 (three); x = 216.5: 218, 225, 230, 235, 240, 244; x = 27.1: 32, 34,
  37, 39, 41, 44. CH3.3's rule: 66–76, 123–135, 234–251, 36–44.
- Fit over 150 points (CH3.3's 96 and 54): total excess Σ (rule − p) CH3.3 556, (A) 486, (B)
  185, (C) 169, (C) capped 188 (the default; at x = 216 222–246, at x = 108 113–132). Capping at
  tenths of a decade as well costs 328 (at x = 216 226–248); at the decades the default stays at
  most CH3.3's, between them it is one degree above at 2 of 5,763 grid points (x = 9.1 at
  10^{−7.2}, x = 10.8 at 10^{−7.8}), tested.
- Pair evidence (`ch4-3c2-pairs.py`): the worst pair at the FMM's measured degrees at x ≥ 43 is
  up to 19–42 ε per ε; the default's at k a = 45 … 256 at most 18.5 ε (k a = 128, 10⁻⁸), within
  the FMM's largest at every grid k a and ε.

- Verification (jobs 663560, matrix: 174 runs within ε, at most 0.74 ε; 663562 and 663563, W1):
  the random charges and the plane wave within ε at kD = 500 and 1000; **the constant density
  misses ε**: kD = 1000 φ 1.22, 1.03, 2.16, 2.60 ε at 10⁻³, 10⁻⁶, 10⁻⁷, 10⁻⁸ (Rotation and the
  default alike), kD = 500 ∇φ 1.12 ε at 10⁻⁶ and 1.58 ε at 10⁻⁸. CH3.3's rule on the same cases
  and code (job 663571): at most 0.46 ε at kD = 1000, 0.91 ε (∇φ, 10⁻⁸) at kD = 500.
- Per level (job 663572, `helmholtz_boxes --study levels`, kD = 1000, constant, Rotation): at
  10⁻⁸ level 2 (x = 216, p 246) carries 0.44 of Σ|e|² in φ, level 3 0.14, level 6 0.09 of ∇φ's;
  level 2 raised to 255 alone halves the error (2.60 → 1.18 ε). At 10⁻³ the error spreads over
  every level (2: 0.31, 3: 0.22, 4: 0.13, 6–7 0.16–0.26 of ∇φ's). The level scans used random
  charges only; the constant density on W1's sphere adds the boxes' errors coherently and needs
  more at every level. The refit as fitted is not acceptable.

- The rule (d = 5, capped): the offset study (job 663585: W1, constant density, kD = 250, 500,
  1000, d = 5, 7, 9 uncapped) found d = 5 enough (at most 0.60 ε); the capped rule verified (jobs
  663586, 663587): the matrix's 174 runs within ε (at most 0.64 ε), W1 at kD = 250, 500 and 1000
  with the three densities within ε (at most 0.94 ε, ∇φ of the constant density at kD = 500,
  10⁻⁸; CH3.3's rule 0.91 ε there). Pair evidence: at most 4.6 ε at k a = 45 … 256. Over the 150
  points the rule's total excess is 367 (CH3.3's 556, the fit's 188). The report is
  `tetraleaf-validate/results/helmholtz-ch4-3c2.md`.

## Pending runs
- None: every run of the plan finished; the same-node times (§7 of the report) are in the report.

## Open questions
- The refit fails the constant density (above): rescan with the worst density and refit, a
  margin, or keep CH3.3's rule? Asked 2026-11-21; the user suggested raising level 2's degree
  slightly and accepting a small cost in time. Answered by the raise of 2 on every level above
  the plateau, capped by CH3.3's rule (above).
- For sign-off (the report §8): adopt the rule; amend the proposal's §3.5; scan the constant
  density in a later calibration.
- The proposal's §3.5 states CH3.3's rule as the calibrated one: amended only on the sign-off
  of the refit (the PR's question).
