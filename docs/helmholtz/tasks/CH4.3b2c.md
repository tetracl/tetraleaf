# CH4.3b2c The diagonal form by default

Started 2026-11-22 on `tbetcke/h-ch4-3b2c`. Roadmap: CH4.3b2c "The diagonal form by default"
(the last of CH4.3b2a–d; its *Needs* CH4.3b2a, CH4.3b2b, CH4.3b2d and CH4.3c2 are merged).
Proposal: §10 CH4.3 ("the time crossover per class; the fallback of a whole level to RCR when the
scratch budget binds"); §3.6 (the rule per (level, class)); §5.3 and §6.4 (D8, Auto's exponential
M2L); §5.4 (the level call); §2 (constraints 1, 2, 7, 9).

## Goal
Make the diagonal form's time rule (`HelmholtzDiagonal::Auto`, as CH4.3b2a and CH4.3b2b leave it)
the builder's default, on the evidence of W1 and W3 at 10⁻³–10⁻⁸ on Kathleen and locust against
the defaults (with CH4.3b2d's `VRaisedPrice::PerBox` and CH4.3c2's degrees), for sign-off.

## Acceptance (the roadmap's CH4.3b2c; proposal §10 CH4.3)
Proposal §10 CH4.3: "C_T and θ from S6 confirmed in the full FMM on W1 and W3; the time crossover
per class; the fallback of a whole level to RCR when the scratch budget binds" (C_T and θ:
CH4.3a; the crossover: CH4.3b1, CH4.3b2a; this task the default and the fallback measured).

The roadmap's line: "W1 and W3 at 10⁻³–10⁻⁸ on Kathleen and locust, the defaults with
`HelmholtzDiagonal::Auto` (as CH4.3b2a and CH4.3b2b leave it) against the defaults (with
CH4.3b2d's `VRaisedPrice::PerBox`, the default since its sign-off): the errors within 1.1×, the
evaluation and the build per level; the order between the waves, the form and RCR per level (D8
amended here if it moves); the whole-level fallback at 40 × 1 (W1 at kD = 1000, levels 2–3);
`Auto` the default (`HelmholtzDiagonal`, `Off` until here) for sign-off. Also (CH4.3b1's report):
the model's misses on W3's small dense levels (the waves kept where the form was 1.1–1.9× faster;
in CH4.3b2a's final runs 1.07–2.10×, and two raised W1 levels); the constants at many ranks per
node (fitted at one rank per node). Also (CH4.3b2a's report §5.2, its question 4): W3 at kD = 200
on Kathleen ran 1.03–1.12× the defaults' evaluation with every choice within 1.05× (the same waves
level slower in the time rule's run): repeated, and spread separated from a placement effect of
the form's buffers (CH4.3b2b's runs saw it too: same-plan runs of one launch up to 30 % apart,
later runs slower; its report §1)."

How it is read:
- **"The errors within 1.1×"**: relative ℓ2 of φ and ∇φ at §9.1's samples, the defaults with
  `--diagonal auto` over the defaults with `--diagonal off`, every case, both machines.
- **"The evaluation and the build per level"**: the evaluation (the second's, as CH4.3b2a) and
  per level the M2L's seconds, the diagonal tables' build per level and the build in all, both
  settings, from the same launch.
- **"The order per level"**: per level with V pairs at k a ≥ 2.8 which M2L the default runs, and
  the options timed at the level's degree in the same launch (the time rule forced to the form and
  to RCR, CH4.3b2a's `FORM_COSTS` and `RCR_COSTS`), scored within 1.05× of the fastest as CH4.3b2a.
  Under the default a level with waves goes to the form where the model has the form below the
  waves at p_w: that moves D8's choice, so D8's amendment is a question here.
- **"The whole-level fallback at 40 × 1"**: W1 at kD = 1000, 10⁻⁶, 40 ranks of one thread on one
  Kathleen node: levels 2–3 (whose scratch does not fit a rank's share of
  `DEFAULT_SIGNATURE_NODE_BUDGET`, CH4.1b) run RCR whole, the same on every rank, and the errors
  with `--diagonal auto` within 1.1× of `--diagonal off`'s in the same launch. (The output differs
  from 1 × 40's where the fallback changes a level's M2L; the bit identity across ranks is the
  tests', at equal choices.)
- **"The constants at many ranks per node"**: the time rule's choices at 40 × 1 (Kathleen) and
  72 × 1 (locust) on W1 at kD = 500 and W3 at kD = 500, 10⁻⁶ and 10⁻⁸, scored against the forced
  options in the same launch, as at one rank. A miss rate well above one rank's is a follow-up
  (constants per ranks × threads), not a change here.
- **"The model's misses"**: listed per level as CH4.3b2a did, on the current defaults; a refit of
  the waves' per-box constant on W3's dense levels is proposed as a follow-up task if the misses
  cost more than a few per cent of the levels' time, not made here (a cost-model change wants its
  own pair and FMM evidence).
- **"Repeated, and spread separated from placement"**: W3 at kD = 200 at the four ε, each launch
  `--diagonal off,auto,off,auto` (one build each, two evaluations each): the waves level's time in
  the two `off` runs gives the spread, `off` after `auto` against `off` before gives the placement
  effect; and the w1/w3 stages again with the order reversed (`--diagonal auto,off`).

## Read first
- `tetraleaf-exec/src/helmholtz/fmm.rs` (`HelmholtzFmmBuilder::new`, `diagonal`),
  `operator_diagonal.rs` (`HelmholtzDiagonal`, `DiagonalSettings`, `wave_level_choice`),
  `operator.rs` module docs "M2L".
- CH4.3b2a's record and report (§4, §5), CH4.3b2b's and CH4.3b2d's reports (§1),
  `tools/helmholtz/ch4-3b2a-{runs.sh,tables.py}`.
- `tetraleaf-exec/tests/helmholtz_fmm.rs` (`default_builder`, `builder`, `check_defaults`,
  `check_diagonal*`), `tetraleaf-validate/src/helmholtz_accuracy.rs` (`Case::new`).

## Changes
- `tetraleaf-exec` (Helmholtz modules only):
  - `HelmholtzFmmBuilder::new`: `diagonal: HelmholtzDiagonal::Auto` (for sign-off); the enum's
    `Default` `Auto` too, as `HelmholtzM2l`'s since CH3.2d; the operator's own default before
    `with_diagonal` stays `Off` (`DiagonalSettings` names it), so the operator tests keep theirs.
  - Docs: `new`, `diagonal`, `HelmholtzDiagonal`, the operator module's "M2L";
    `tetraleaf-exec/CLAUDE.md`.
- `tetraleaf-exec/tests/helmholtz_fmm.rs`: `default_builder` names `HelmholtzDiagonal::Off` (the
  settings before, as it names the tree and near rules), so every check keeps its reference;
  `check_defaults` asserts the default `Auto` (`operator().diagonal()`), and its named builder
  names it, so the default is checked against the direct sum, named bit for bit, 2 threads and
  the one-rank FMM.
- `tetraleaf-validate`: `Case::new` keeps `Off` named (the earlier reports' runs repeat; the
  example's `--diagonal` default `off`); its docs and `CLAUDE.md` say the library's default.
- `tools/helmholtz/ch4-3b2c-runs.sh` (stages `w1`, `w3`, `reverse`, `repeat`, `ranks`, `smoke`)
  and `ch4-3b2c-tables.py` (CH4.3b2a's `errors`, `choices`, `summary` imported, plus `build`,
  `repeat` and `order`); the report `tetraleaf-validate/results/helmholtz-ch4-3b2c.md`.
- On sign-off: D8 amended in the roadmap's decision table and the proposal's §5.3 and §6.4, with
  a line in Appendix A.
- No convention, no dependency, no `tetraleaf-plan` interface, no Laplace file.

## Tests that prove it
- `check_defaults` (the sphere at 10⁻⁶ with a tolerance, the clusters with named degrees): the
  default is `Auto`; its errors within the scenario's bound against the direct sum; the settings
  named bit for bit the default's; 2 threads and the one-rank FMM bit for bit (on 2 and 4 ranks
  in CI's MPI list).
- `check_diagonal_auto`, `check_time_rule_raised` (CH4.3b1, CH4.3b2a) keep testing `Auto` against
  the setting off within `DIAGONAL_RATIO` (1.1), threads, chunks and one rank, unchanged.
- The existing checks unchanged in their references (`default_builder` names `Off`).
- `tetraleaf-validate`'s `tests/helmholtz_accuracy.rs` unchanged (`Case::new` names `Off`).

## Checks
`/verify` (fmt, Clippy default and workspace all-targets, rustdoc `-D warnings`, the default and
workspace tests at opt-level 1, the MPI list at 2 and 4 ranks: `helmholtz_fmm` changes); the
one-rank `helmholtz_fmm` at opt-level 0 against its budget.

## Machine runs
- **Kathleen** (`singlenode`, one node per job): `ch4-3b2c-runs.sh` stages `w1 w3` at 1 × 40
  (W1 at kD = 100, 200, 500 at 10⁻³, 10⁻⁶, 10⁻⁷, 10⁻⁸ and kD = 1000 at 10⁻⁶; W3 at N = 10⁶,
  kD = 200 and 500 at the same ε; `--diagonal off,auto --diagonal-costs machine,form,rcr`),
  `reverse` (the same with `--diagonal auto,off`, the machine's costs), `repeat` (W3 at kD = 200,
  `--diagonal off,auto,off,auto`), `ranks` (40 × 1: W1 at kD = 1000, 10⁻⁶, the fallback; W1 and
  W3 at kD = 500, 10⁻⁶ and 10⁻⁸ with the three costs). Decides: the acceptance, the fallback, the
  constants at 40 ranks per node, the W3 Also.
- **locust** (1 × 72, the load sampled every 30 s): `w1 w3 reverse repeat`, and `ranks` at
  72 × 1 without the kD = 1000 fallback case. Decides: the acceptance on NEON's constants.
- No M6 timing.

## Constraints and risks
- Constraint 1: no Laplace path (Helmholtz modules, the harness, scripts, docs).
- Constraints 2 and 9: the default's choices are CH4.3b1–CH4.3b2a's, functions of agreed counts,
  (p, k̂, ε, T), agreed constants and the machine class; tested bit for bit on threads and ranks.
- Constraint 7: the form's scratch within `DEFAULT_SIGNATURE_NODE_BUDGET` per node, the fallback
  measured at 40 × 1.
- Risk: the default makes the library's output depend on the machine class in one more place
  (as Batched's crossover and the V costs already do); documented.
- Risk: the debug test time of `check_defaults` grows by the diagonal tables of the sphere.

## Pending runs
None. The runs are done (raw outputs in `~/projects/fmm-archive/{kathleen,locust}/runs/ch4-3b2c*/`):
- Kathleen job 663643 (node-c11b-048, `singlenode`, 2026-11-22 11:44–12:14), `w1 w3` at 1 × 40,
  21 launches, every one exit 0; job 663644 (node-c11a-058, 11:45–12:21), `repeat ranks reverse`,
  30 launches, every one exit 0;
- locust 11:45–12:08, `w1 w3 repeat reverse ranks`, 50 launches, every one exit 0, the load sampled
  every 30 s (1-minute load within the runs' own processes: at most 53 under 72 threads, 73 under
  72 ranks). A first start at 11:45 was stopped after 8 launches by mistake and restarted from
  scratch; its outputs were deleted.

## Results
The report: tetraleaf-validate/results/helmholtz-ch4-3b2c.md. The default errs 0.991–1.034× the
defaults' on every case (both machines, 1 × N and N × 1); the evaluation 0.37–0.99× on Kathleen and
0.44–0.99× on locust (the geometric mean of both run orders); the choices within 1.05× on 75/78 and
82/90 levels; the fallback at 40 × 1 runs W1 at kD = 1000's levels 2–3 on RCR whole, and at one rank
per core it binds on more levels (W3 at kD = 500 at 40 × 1 and 72 × 1, W1 at kD = 500 at 72 × 1). At
40 × 1 the waves on W3 at kD = 500's level 4 (10⁻⁸) cost 8× their one-rank time (18 s lost, the
defaults too) and every rank builds the diagonal tables (about 10 s). W3 at kD = 200 on Kathleen is a
run-order drift (the same computation 9–20 % slower the second time in a launch), not the form.

The self-review (a subagent, with the plan, proposal §2 and §10 and the branch) found no blocking
code issue; it corrected the report's build-plus-evaluation range (both run orders), added the
fallback at one rank per core on more levels (locust's 72 × 1), and asked for the rustdoc of
`HelmholtzFmmBuilder::diagonal` to say that the levels depend on the ranks per node where the budget
binds; all applied.

Debug test time (opt-level 0, one rank, the M6): `helmholtz_fmm` 700.5 s at load 2–4 (CH4.3b2d
703 s); at opt-level 1 39.5 s.

## Open questions
The report's §9, questions 1–5, for sign-off.
