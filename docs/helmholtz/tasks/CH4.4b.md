# CH4.4b The crosscheck mapped to levels and classes

Started 2026-11-23 on `tbetcke/h-ch4-4b`. Roadmap: "CH4.4b The crosscheck mapped to levels and
classes" (H4; priority list item 10, after CH4.4a, its *Needs*, merged as #255 on 2026-11-23).
Proposal: §10 CH4.4 (the part "mapped to levels and classes, with no error concentrated at a class
boundary"); §9.1 (the box measure, its 100 ε bound, amended at CH4.0b, "CH4.4's 'no error
concentrated at a class boundary' is read with it"); §3.6 (the rule per (level, class), the bound E,
θ = 1/2, the regimes); §2 (constraints 1, 2, 6, 9). Builds on CH4.0b (the box study and the levels'
parts, `tetraleaf-validate/results/helmholtz-ch4-0b.md`) and CH4.4a (the cases, the 114.7 ε box,
`tetraleaf-validate/results/helmholtz-ch4-4a.md`).

## Goal
Map the full-target crosscheck's box errors under the library's defaults to the levels and, on the
diagonal form's levels, to the form's classes against RCR; show whether any error concentrates at a
class boundary or where the M2L changes; and find the cause and a remedy of W1's 114.7 ε box at
kD = 1000 (plane wave's trace, 10⁻⁸).

## Acceptance (proposal §10, verbatim)
CH4.4's row: "FMM against the direct sum (sampled targets) on W1 to kD = 1000, W2 to 500 and W3 to
200, beside `hfmm3d` against the same sum; the full-target crosscheck mapped to levels and classes,
with no error concentrated at a class boundary; time, memory and build cost per k". CH4.4a took all
but the mapping. The roadmap's *Done when* of CH4.4b: on CH4.4a's cases whose trees run the
diagonal form or RCR at k a ≥ 2.8, the full-target crosscheck's box errors split by level (CH4.0b's
level parts, one level's degree raised at a time) and, on the form's levels, by the form's classes
against RCR (the level run all-RCR), with no error concentrated at a class boundary: neither in the
classes beside the bound's threshold (§3.6's choice per (level, class)) nor at the levels where the
M2L changes (the regimes, the waves, the form); the worst boxes' levels and classes named. Also
(CH4.4a's question 2, signed off): the worst box of W1 at kD = 1000 with the plane wave's trace at
10⁻⁸, 114.7 ε (∇φ), its levels and classes, and a remedy (a degree on a level, or the bound) on pair
and FMM evidence.

How it is read (to be confirmed in the PR, question 1):
- **The cases.** CH4.4a's cases whose trees hold a level of k a ≥ 2.8 with the diagonal form or RCR
  under the defaults (CH4.4a's logs, "by level"): W1 at kD = 30 (10⁻⁶ only: Batched at 10⁻³, waves
  at 10⁻⁸), 100, 250, 500, 1000, with the plane wave's trace at 500 and 1000 and the constant density
  at 100 (and 30, 10⁻⁶); W2 at kD = 30 (10⁻⁶), 100, 250, 500, the constant density at 100 (and 30,
  10⁻⁶); W3 at N = 108,849, kD = 30 (10⁻⁶), resolved at kD = 100, at N = 10⁶ for kD = 200 and 500.
  ε = 10⁻³, 10⁻⁶, 10⁻⁸, CH4.4a's. One complex vector, gradients on.
- **The reference.** The f64 direct sum (`tetraleaf_ref::helmholtz::p2p::direct_sum` × 1/(4π)) at
  every point up to N = 2 · 10⁵ and above at §9.1's samples, two targets per box of the deepest
  level, and every point of the base run's worst boxes (new: `--refine n`, the n worst of φ and of
  ∇φ by the sampled measure, so the worst boxes' maxima are exact, not lower bounds). CH4.4a's
  reference was FMM3D at 10⁻¹²; the worst box is checked against CH4.4a's figure (the same output,
  the same box: equal to the digits FMM3D's own error allows, ≈ 4 · 10⁻¹³ relative).
- **Levels.** Each level's part of the error as CH4.0b defined it (the run minus the run with that
  level raised; the parts add up, the rest printed), but under the defaults: the level raised to the
  rule's degree at ε/1000 (at least one more than as run) and run RCR, every other level as the
  defaults choose. CH4.0b ran its parts on all-Rotation, since named degrees change the waves'
  precision and take the form off; a new measurement setting (below) raises one level while the
  tolerance, the rules and the plan stay.
- **The M2L.** On each level whose M2L under the defaults is not RCR (the form, the waves,
  Batched), the run with that level on RCR at its degree as run (its resolving degree p_w where
  raised): the M2L's own part, its share and its effect at the worst boxes. "No error concentrated
  where the M2L changes": the default's worst box and ℓ2 at most 1.1× that run's
  (`DIAGONAL_RATIO` and `PLANE_WAVE_RATIO`, the tests' factor for an M2L against RCR). One-sided
  (amended 2026-11-23 from "within 1.1× either way", after the stage `low` showed it, for
  question 1: the waves
  approximate the kernel, not RCR's truncation at p_w, so RCR at p_w may err more than they do,
  1.19× on W2 at kD = 30, which is no error the waves add).
- **The classes.** On each form level the classes by the bound's E(l, c)/(θ ε), in decades: the
  passing classes within a decade of the threshold ("beside the threshold") run RCR as one run, then
  the next decade; on the worst case (W1, kD = 1000, plane wave, 10⁻⁸) also each class of the top
  decade alone. "No error concentrated at a class boundary": the threshold decade's part within the
  same 1.1×, and no worst box whose error comes from it. The classes that fail the bound run RCR in
  the default already (exact), so the boundary's other side is RCR's truncation, which the level
  parts cover.
- **The worst boxes.** For each case and ε the worst boxes of φ and ∇φ with their place (cell,
  centre, the levels' |v|/r) and the levels' and the M2L runs' parts at their worst targets (the
  class decades' runs give the classes; a box's V offsets span nearly every class, so they are
  not listed).
- **The 114.7 ε box.** The level (and class) that carries it; the remedy proposed from (a) the FMM
  with that level's degree raised by 1 and 2 under the defaults' M2L (the form at the new degree, the
  new setting), its worst box and ℓ2, (b) the per-pair evidence (spike S7's curves,
  `spikes/truncation/results/c1_curves.csv`: the worst pair over the 316 offsets at the level's k a
  and the degrees, over ε), both measures as the H6 gate's rule asks; or §9.1's bound, if no degree
  is justified. A remedy that changes the degree rule is a recalibration: its sign-off decides it and
  a follow-up task implements it (CH4.3c's scans and verification), not this PR.

## Read first
- `tetraleaf-validate/results/helmholtz-ch4-0b.md` §1–§2 (the method), `helmholtz-ch4-4a.md` §2–§3.
- `tetraleaf-validate/src/helmholtz_boxes.rs`, `examples/helmholtz_boxes.rs` (the box study).
- `tetraleaf-exec/src/helmholtz/operator.rs` `select_tables` (each level's M2L, degree and waves);
  `operator_diagonal.rs` `select`, `diagonal_rule_classes`; `fmm.rs` (the builder's settings).
- `tetraleaf_tables::helmholtz::diagonal::class_bounds` (E per class).

## Changes
- `tetraleaf-exec` (Helmholtz only; no Laplace path, no `tetraleaf-plan` change):
  - `helmholtz::M2lOverride { level, degree: Option<usize>, m2l: OverrideM2l }`,
    `OverrideM2l::{Setting, Rcr, RcrClasses(Vec<usize>)}`: a measurement setting.
    `degree` replaces the level's base degree before the rules read it (the waves' resolving
    degree, the time rule and the form then work from it); `Setting` keeps the setting's choice on
    the level; `Rcr` runs the level's M2L as RCR at its degree as run (the base degree if named,
    else p_w where the waves would have raised it), no waves, form, Batched, Dense or Classes there;
    `RcrClasses(c)` takes the classes c off the form on that level. Its translations follow the
    degree as `PlaneWaveRaise` gives them; the plan is not changed (the builder's V and W/X rules
    read the setting before the override, so the base run and the override's run share the plan).
  - `HelmholtzOperator::with_m2l_overrides`, `HelmholtzFmmBuilder::m2l_overrides`, reported by
    `HelmholtzOperator::m2l_overrides` (through `HelmholtzFmm::operator`); a level below 2, a
    level named twice or a class ≥ 34 is `HelmholtzFmmError::M2lOverride` (before any collective,
    the same on every rank); a level beyond the tree is ignored (documented). Every public item
    documented; the crate's CLAUDE.md.
  - As built, the override's degree is applied by the builder to the resolved degrees (so the
    plan, the near rules and the waves' levels read it, and a degree under `Setting` cannot add
    waves on a level the plan did not group for them); its M2L acts in the operator alone.
- `tetraleaf-validate`:
  - `helmholtz_accuracy::Case::m2l_overrides` (empty in `Case::new` and `Case::defaults`), passed to
    the builder.
  - `helmholtz_boxes`: `Reference::extend` (more targets, the values at the old targets kept, bit
    for bit what `Reference::at` gives at the union); per-level class bounds for a run
    (`class_ratios`: E(l, c)/(θ ε) of each class at the level's degree, k̂ and the build's constants).
  - The `helmholtz_boxes` example: `--m2l defaults` (`Case::defaults`), `--refine n`,
    `--reference-cache dir` (the sampled reference shared by a case's ε: about 12 min a launch
    at kD = 1000), `--study
    levels` under the defaults by overrides (under `rotation` CH4.0b's named degrees, so its runs
    repeat), `--study m2l` (each non-RCR level on RCR, the class decades on the form's levels,
    `--class-runs` for each class of the top decade), `--raise-by 1,2` (the degree-raised runs
    under the setting's M2L, for the remedy), and tables: the levels' M2L as run (k a, p, M2L,
    classes, the threshold decades), the M2L's parts (share of Σ|e|², worst box and ℓ2 over the
    default's), the worst boxes' parts by level and by M2L.
  - The crate's CLAUDE.md entry.
- `tools/helmholtz/ch4-4b-runs.sh` (the jobs' stages) and `ch4-4b-tables.py` (the report's tables
  from the outputs); report `tetraleaf-validate/results/helmholtz-ch4-4b.md`.
- Conventions: none.

## Tests that prove it
- `tetraleaf-exec`, `tests/helmholtz/` (no MPI): the override in the operator: a level's degree,
  its M2L as named (`plane_wave_level` none, no diagonal level, no batched tables), the classes left
  out, the other levels' choices unchanged, an empty list the setting's tables; the errors refused.
- `tests/helmholtz_fmm.rs`, `check_m2l_override` (1, 2, 4 ranks in CI's MPI list), each against a
  trusted path:
  - under Rotation with the form off, `{level l, degree q, Rcr}` bit for bit the FMM with
    `Degrees::PerLevel` naming q on level l (the named-degree path; the leaf size named);
  - under the form at free costs (`check_diagonal_auto`'s case), `Rcr` on every form level bit for
    bit the form off, `RcrClasses(all)` bit for bit `Rcr`, `RcrClasses(some)` the rule's classes less
    those, errors within `DIAGONAL_RATIO` of the form off;
  - under the defaults on a tree with waves, `Rcr` on a wave level: no waves there, its degree p_w,
    errors within `PLANE_WAVE_RATIO` of the default; `{degree p + 1, Setting}` the level's degree;
  - an empty list the default's bits; 2 threads and the one-rank FMM bit for bit for each.
- `tetraleaf-validate`, `tests/helmholtz_accuracy.rs`: `Case::m2l_overrides` reaches the builder
  (a row's degrees and M2L by level); unit tests of `Reference::extend` (bit for bit `Reference::at`
  of the union) and `class_ratios` (against `class_bounds`).
- The level parts' rest (CH4.0b's check that they add up) printed for every run, ≤ 1 % of Σ|e|²
  expected; a larger rest is reported, not hidden.

## Checks
`/verify` (fmt; Clippy default, workspace all-targets; rustdoc `-D warnings`; the default and
workspace tests at opt-level 1; the MPI list at 2 and 4 ranks, `helmholtz_fmm` changed);
`cargo clippy -p tetraleaf-validate --all-targets -- -D warnings`, `cargo test -p tetraleaf-validate`;
`tests/helmholtz_fmm.rs` at one rank, opt-level 0 (its debug cost recorded). No Laplace path: no
D9 benchmarks.

## Machine runs
- The M6 (accuracy only, no timing): the example's smoke runs, W1 at kD = 100 at the three ε
  (every-point reference), each study.
- Kathleen, `singlenode`, one rank × 40 threads (the box study is one rank), release, BLAS threads
  1, node-c11b-004 and node-c11b-013 excluded: every case above at its ε, `--m2l defaults --study
  levels,m2l --refine 8`; on W1 at kD = 1000 and 500 with the plane wave's trace at 10⁻⁸ also
  `--class-runs` and `--raise-by 1,2` on the level that carries the worst box. Estimated 5–6
  node-hours (the kD = 1000 cases about 2.5 h: two densities × 3 ε × about 16 runs of about 1 min,
  the reference about 10 min each), in several jobs of at most 6 h. The worst box of each base run
  compared with CH4.4a's.

## Constraints and risks
- Constraint 1: no Laplace path changes (Helmholtz operator and builder, the harness); Laplace's
  outputs unchanged by construction, its benchmarks not run.
- Constraint 9: the override is a named setting, not a rule; the default (no override) unchanged
  bit for bit (`check_m2l_override`'s empty list; every other check unchanged).
- Risk: raising a level under the time rule could change that level's M2L choice (the form at a
  higher degree may fail more classes); the level-raised runs run RCR there by construction, so
  the part is the level's whole error (truncation and M2L), as CH4.0b's.
- Risk: memory at kD = 1000 on one rank (CH4.4a: 29 GB over four ranks, plus the reference and
  the extra outputs: one `Errors` per run, 64 B per reference target): within the node's 180 GB.
- Risk: the reference's cost at kD = 1000 (about 1.5 · 10⁵ targets × 8 · 10⁶ sources).

## Pending runs
Kathleen, `singlenode`, 1 × 40, `--exclude=node-c11b-004,node-c11b-013`, outputs in
`/scratch/scratch/ucahtbe/fmm/logs/ch4-4b/` (from 51c914d's tree). Job 663774 (0d6749c) built and
ran the smoke case (W1 at kD = 30, 10⁻⁶: 3.34 / 8.26 ε, CH4.4a's figures). The first stage jobs
(663775 `low` completed, 663776–663779 cancelled) ran without `--refine-faces`: the sample of two
per box missed CH4.4a's worst box at kD = 500 (29.1 against 40.4 ε); their outputs are kept in
`superseded/`. Job 663780 rebuilt with it and reran W1 at kD = 500, 10⁻⁸: 27.59 / 40.47 ε against
CH4.4a's 27.38 / 40.42 ε (FMM3D's reference at 4 × 10). Pending:
- 663781 `tools/helmholtz/ch4-4b-runs.sh … low`
- 663782 `… mid`
- 663783 `… w3`
- 663784 `… w1-1000`
- 663785 `… pw1000 remedy` (the 114.7 ε case and its remedy runs)

Each decides the mapping of its cases; 663785 also the cause and the remedy (question 2).

## Open questions
(filled when they arise)
