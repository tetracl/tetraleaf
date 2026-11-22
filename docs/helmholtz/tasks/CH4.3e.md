# CH4.3e c_leaf re-decided with the rule

Started 2026-11-22 on `tbetcke/h-ch4-3e`. Roadmap: "CH4.3e c_leaf re-decided with the rule" (H4;
priority list item 10, after CH4.3d1; CH4.3c2's sign-off put CH4.3c3 after it; its *Need* CH4.3b2c
is merged). Proposal: §10 CH4.3 (the per-(level, class) rule, calibrated; CH4.3e is the roadmap's
sub-task of it, no row of its own); §3.7 (leaves and tree depth, D11); §5.1 (the leaf operators'
arguments); §5.7 (the P2P's phase range); §9.1 (W1–W3, the error norms); §2 (constraints 1, 2, 9).
Builds on CH4.0d3 (c_leaf 2, `tetraleaf-validate/results/helmholtz-ch4-0d3.md` §4) and CH4.3b2c (the
time rule the default, `results/helmholtz-ch4-3b2c.md`).

## Goal
Decide the builder's c_leaf (`DEFAULT_C_LEAF`, D11) again now that the diagonal form's time rule is
the default: measure where it puts leaves on the levels whose M2L is high frequency (the form, or
RCR at k a ≥ 2.8), check that such trees stay within ε and the tests' bounds, and time c_leaf 1, 2,
4 and 8 on W1–W3, for sign-off.

## Acceptance (the roadmap's *Done when*; no proposal §10 row of its own)
"(since CH4.0d3's c_leaf 2, signed off on #219) c_leaf re-decided with CH4.3b's default (at c_leaf
2 the census's first pass put leaves on high-frequency levels on 9 trees at 10⁻⁶; proposal §3.7
keeps W and X lists off them; CH4.0d3's report §4), on W1–W3's trees, for sign-off (D11)."

How it is read:
- **"Leaves on high-frequency levels"**, as built: a level is high frequency where the default's
  M2L is the diagonal form (`diagonal_classes` > 0) or RCR at k a ≥ 2.8 (S6's range, the time rule's
  `DIAGONAL_AUTO_MIN_KA`); the census's first-pass classifier is replaced by what the FMM runs. Per
  tree the occupied leaves on such levels and the W and X entries that touch them (M2P into them,
  P2L out of them), from the harness's `--levels` tables.
- **"§3.7 keeps W and X lists off them"**: §3.7's premise is that the leaf operators and their
  arguments stay low frequency. The leaf operators (P2M, L2P, P2L, M2P) do not depend on the M2L's
  regime: they run at the level's degree from the degree rule, calibrated to x = 216 (CH4.0c,
  CH4.3c2), and their arguments stay below k a_leaf ≤ 5.4 c_leaf (§5.1). So the check is
  empirical and structural: a test with leaves, W and X entries on a level that runs the form
  (below), and every sweep run's errors within ε. If both hold, §3.7's consequence is amended for
  sign-off (leaves may lie on the form's levels), not enforced.
- **"Re-decided"**: c_leaf 1, 2, 4 and 8 under the defaults (the time rule, the leaf-size rule, the
  V, W and X rules, the interpolated build, `auto-pw`), each run's build plus one evaluation and the
  evaluation alone over c_leaf 2's (the default), geometric means per ε over the cases, both run
  orders (Kathleen's run-order drift, CH4.3b2c's report §6). As CH4.0d3 (T11's rule): another value
  only where its geometric mean over a machine class's cases is at least 5 % lower, at errors within
  ε on every case; else c_leaf 2 stays. The losers per case listed.
- **The wider scope** (H6 gate §7, "a wider scope widens the verification set"): c_leaf above 2
  widens the P2P's phase range (Re(k) r ≤ 44 c_leaf with V pairs near, 176 at 4, 352 at 8, within
  the kernel's contract domain 2¹⁰, `P2pKernel::evaluate_helmholtz`, per-term bound (32 + 12 |k| r) u
  ≈ 9.4·10⁻¹³ at 352) and the leaves' x (≤ 21.8 at 4, 43.5 at 8, inside the degree rule's calibrated
  range). Each run reports its phase bound (`Row::phase_bound`) and its errors at 10⁻⁸, the tightest
  ε; a value above 2 is proposed only with those.

## Read first
- Proposal §3.7, §5.1, §5.7; roadmap D11; CH4.0d3's report §4 and record; CH4.3b2c's report §4, §6.
- `spikes/census/SPIKE_REPORT.md` §3–§4 (the first pass, 43 trees).
- `tetraleaf-exec/src/helmholtz/mod.rs` (`min_leaf_level`, `C_LEAF`, `DEFAULT_C_LEAF`),
  `fmm.rs` (`c_leaf`, `leaves_by_level`, `far_interactions`), `operator_diagonal.rs`
  (`DIAGONAL_AUTO_MIN_KA`).
- `tetraleaf-exec/tests/helmholtz_fmm.rs` (`check_leaf_sizes`, `check_diagonal_auto`,
  `diagonal_scenarios`, `FREE_DIAGONAL`).
- `tetraleaf-validate/examples/helmholtz_accuracy.rs` (`--c-leaf`, the tables "Time rule", "Leaves",
  "Leaf operators"), `tools/helmholtz/ch4-0d3-*` and `ch4-3b2c-*`.

## Changes
- `tetraleaf-exec/tests/helmholtz_fmm.rs`: `check_form_leaves` (below), on a clusters scenario.
- `tools/helmholtz/ch4-3e-runs.sh` (stages `w1`, `w2`, `w3`, `large`, `smoke`; each case in both
  c_leaf orders) and `tools/helmholtz/ch4-3e-tables.py` (`census`, `summary`, `cases`, `errors`).
- The report `tetraleaf-validate/results/helmholtz-ch4-3e.md`.
- `DEFAULT_C_LEAF` changed only if the measurements propose it, built in this PR for sign-off and
  reverted if not signed off (as CH4.0d3).
- On sign-off: D11 and proposal §3.7 (and §5.1 where it says every leaf is low frequency) amended,
  with a line in Appendix A.
- No harness code (the rows hold the leaves per level, the classes per level, the M2L that ran, the
  W and X entries per level and the phase bound); no convention, no dependency, no `tetraleaf-plan`
  interface, no Laplace file.

## Tests that prove it
- `check_form_leaves` (`tests/helmholtz_fmm.rs`; 1, 2, 4 ranks in CI): Gaussian clusters at
  kD = 40, ε = 10⁻⁴, c_leaf 2, leaves of up to 64 points, the W and X rules off (so W and X pairs stay
  far), the time rule at free form costs (the form wherever the bound passes a class, whatever the
  machine class): asserts that some level running the form holds occupied leaves, that M2P (W
  entries) and P2L (X entries) touch those leaves;
  the errors against the direct sum within the scenario's bound and within `DIAGONAL_RATIO` (1.1) of
  the same tree with the form off; 2 threads and the one-rank FMM bit for bit. The same at the
  default costs printed (no assertion on which levels the machine's costs take).
- Every existing check unchanged (`default_builder` names `C_LEAF`; `check_defaults` checks the
  library's c_leaf, whatever it is).

## Checks
`/verify` (fmt; Clippy default and workspace all-targets; rustdoc `-D warnings`; the default and
workspace tests at opt-level 1; the MPI list at 2 and 4 ranks: `helmholtz_fmm` changes; the
one-rank `helmholtz_fmm` at opt-level 0 against its budget). Shell and Python scripts checked with
`bash -n`, `shellcheck` if present, and `python3 -m py_compile`.

## Machine runs
All release, BLAS threads 1, `helmholtz_accuracy --levels --evaluations 2 --m2l auto-pw --diagonal
auto --diagonal-build interpolated --near-v cost --near-wx cost` (the library's defaults named; the
leaf size the library's rule), each case once with `--c-leaf 1,2,4,8` and once with `8,4,2,1`.
- **Kathleen** (`singlenode`, 1 × 40, one job per stage group): W1 at kD = 100, 200, 500 at 10⁻³,
  10⁻⁶, 10⁻⁷, 10⁻⁸; W2 at kD = 100, 200, 500 at the same ε; W3 at N = 10⁶, kD = 200 and 500, at the
  same ε; `large`: W1 at kD = 1000 at 10⁻³ and 10⁻⁶, W3 resolved at kD = 100 (N = 4·10⁶) at 10⁻³,
  10⁻⁶ and 10⁻⁸. Decides: the census, the errors, the timing rule.
- **locust** (1 × 72, the load checked before and sampled during): the W1 and W3 stages (NEON's
  leaf sizes and costs). Decides: the same on the other machine class.
- No M6 timing.

## Constraints and risks
- Constraint 1: no Laplace path (a Helmholtz test, scripts, docs; `DEFAULT_C_LEAF` is the Helmholtz
  builder's).
- Constraints 2 and 9: c_leaf changes the tree, not the guarantees; the new test checks 2 threads
  and the one-rank FMM bit for bit at c_leaf 2 with leaves on the form's level; the leaf rule reads
  k, the domain and c_leaf only.
- Risk: c_leaf 4 and 8 change nothing where the points per leaf bind (W1's surface at the rule's
  256 points); they matter on under-resolved trees (W3 at N = 10⁶, kD = 500), where the leaves of
  c_leaf 2 hold a few points. A rule by tree would read the points (constraint 9); only one value
  per machine class (or per ε) is proposed.
- Risk: Kathleen's run-order drift (9–20 %, CH4.3b2c §6); both orders, geometric mean.
- Risk: the debug cost of the new check (a clusters FMM at kD = 40 with the diagonal tables).

## Pending runs
None. The runs are done (raw outputs in `~/projects/fmm-archive/{kathleen,locust}/runs/ch4-3e/`):
- Kathleen, `singlenode`, one rank on 40 cores, at fc3b38d: job 663653 (node-c11b-012, 2026-11-22
  14:19–14:53, `w1` and W1 at kD = 1000), 663654 (node-c11a-021, 14:20–14:34, `w2`), 663655
  (node-c11a-021, 14:34–15:18, `w3` and W3 resolved at kD = 100); 74 launches, every one exit 0;
- locust 14:21–14:40, `w1 w3` at 1 × 72, 40 launches, every one exit 0, the 1-minute load 1.3
  before and at most 44.7 during (the runs' own threads).

## Results
The report: tetraleaf-validate/results/helmholtz-ch4-3e.md. Leaves lie on levels of k a ≥ 2.8 (the
form's range) at every c_leaf from 1 to 8 (a leaf of c_leaf λ has k a ≤ 5.4 c_leaf), whatever their
M2L (the form, RCR, the plane waves or none), on W2 with W and X lists; every run within ε (at most
0.27 ε in φ, 0.44 ε in ∇φ). c_leaf 4 puts build plus one evaluation at 0.874× (Kathleen)
and 0.790× (locust) c_leaf 2's, from the under-resolved W3 at kD = 500 (0.26–0.36×), the rest the
same tree or 0.93–0.98×; c_leaf 8 only 1–3 % more, from the same case; c_leaf 1 1.42–1.52×. Built:
`DEFAULT_C_LEAF` = 4, for sign-off.

## Open questions
Signed off 2026-11-22 by the user on #253: "Q1-3 as recommended" (D11's row, proposal §3.7, §5.1,
§5.7 and Appendix A amended). The report's §6: 1, the acceptance as read; 2, c_leaf 4 the default (or 8, or 2),
with an Also on CH4.5 (c_leaf 4 at many ranks per node); 3, proposal §3.7, §5.1 and §5.7 amended to
the finding (leaves may lie on levels of k a ≥ 2.8; under-resolved inputs get leaves up to c_leaf
λ), D11's row amended.

The self-review (a subagent with the plan, proposal §2 and §10 and the branch) found nothing
blocking; fixed: the census counted only the form's and RCR's levels (now every level of k a ≥ 2.8,
the M2L named), the test counted X entries of empty leaves (now the point interactions), the
report's wording (§3.7's third bullet in question 3, the tie-break between 4 and 8, the losers,
one ∇φ range), the doc comment's claims; not changed: the locust load sampled without its
processes (stated in the report).
