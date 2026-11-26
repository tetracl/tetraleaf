# CH4.6a2 The time per level against the product, both grids priced

Started 2026-11-26 on `tbetcke/h-ch4-6a2`. Roadmap: CH4.6a2 (item 13 of the priority list; CH4.6a
split at CH4.6a1, the user: "Split a1/a2/a3 (Recommended)"). Proposal: §10 CH4.6a (its last two
clauses); §6.7 ("Grid sizes", "Host back ends"), §12 D7; §2 constraints 1, 2, 9. Spike S9:
`spikes/fft-plans/SPIKE_REPORT.md` §5–§7. Record of CH4.6a1: `docs/helmholtz/tasks/CH4.6a1.md`.

## Goal
Measure what rustfft's φ step gains per diagonal level against the product on the product's
grid, price the smooth grid from S9's per-box FFT times and the level's own pair cost, and ask
whether CH4.6a3 (the smooth grid) is built.

## Acceptance (roadmap *Done when*; proposal §10 CH4.6a, its last clauses)
Roadmap: "each diagonal level's conversions and level call, the product against rustfft on the
product's grid, on locust, Kathleen and the M6 (the M6 with the user's go-ahead; S9 did not time
it), on W1 at kD = 250–1000 and W3 at kD = 500; the smooth grid priced from S9's per-box FFT times
and the level's measured pair cost at its larger N_φ; a question for sign-off: build the smooth
grid (CH4.6a3) or keep the product's grid and close CH4.6a3 unrun. Also (from CH4.6a1): the
harness (`tetraleaf-validate`'s `helmholtz_accuracy`) takes the φ step as a setting and records the
back end and rustfft's CPU features that ran (proposal §6.7: results files record the back end);
and whether the time rule's `DiagonalCosts` need an FFT variant".

Proposal §10 CH4.6a (verbatim, the clauses left to CH4.6a2 and CH4.6a3): "a level's grid (the
product's or the smooth one) chosen by its measured time; the time per level against the product
on the M6, locust and Kathleen".

**The M6:** asked on 2026-11-26 before any run ("May I run the timing launches on the M6"); the
user: "Skip the M6". No M6 times in CH4.6a2; the PR asks to sign off the acceptance without them
(question 3), as S9 did.

## Read first
- `tetraleaf-validate/src/helmholtz_accuracy.rs` (`Case`, `Row`, `run`'s builder),
  `examples/helmholtz_accuracy.rs` (`--diagonal-build`'s loop, `diagonal_suffix`, the tables
  "Diagonal M2L" and "Time rule").
- `tetraleaf-exec/src/helmholtz/operator_phi.rs` (`HelmholtzPhiStep`, `rustfft_cpu_features`,
  `RUSTFFT_VERSION`), `operator.rs` (`phi_step_on`), `fmm.rs` (`HelmholtzFmmBuilder::phi_step`),
  `operator_diagonal.rs` (`DiagonalCosts`, the model, `DIAGONAL_COSTS_{X86,NEON}`).
- `spikes/fft-plans/examples/times.rs` (`times <variant> <grid> [p,…]`), its results per machine.
- `tools/helmholtz/ch4-3b2c-runs.sh` (the launch helpers), `tools/helmholtz/common.py`.

## Changes
- **tetraleaf-validate** (a measurement crate; no library, no Laplace path):
  - `Cargo.toml`: a feature `rustfft = ["tetraleaf-exec/rustfft"]`, in no default set (no new
    dependency: it forwards tetraleaf-exec's optional one; `cargo tree -e normal -i rustfft` of
    the default members stays empty).
  - `Case::phi_step` (`HelmholtzPhiStep`, `Product` in `Case::new` and `Case::defaults`, the
    library's default, named always) into `HelmholtzFmmBuilder::phi_step`.
  - `Row::phi_step` (the setting as built, `HelmholtzOperator::phi_step`) and `Row::phi_steps`
    (per level the φ step that ran, `phi_step_on`: `None` without the form, the product where the
    ranks' agreement put the level on it).
  - The example: `--phi-step product|rustfft,…` (each `--diagonal` run other than `off` with each
    in turn, inside `--diagonal-build`'s loop; `rustfft` refused before MPI starts in a build
    without the feature), the column "M2L" ending `/rustfft`, a line "φ step:" under the
    machine's (rustfft's version and `rustfft_cpu_features()` on this rank; "product" alone
    without the feature), and a column "φ" in "Diagonal M2L" (the level's step as run).
  - `tetraleaf-validate/CLAUDE.md`: a paragraph for CH4.6a2.
- **Report scripts**: `tools/helmholtz/ch4-6a2-runs.sh` (the stages below) and
  `tools/helmholtz/ch4-6a2-tables.py` on `common.py` (the level tables, the smooth-grid price,
  the time rule's question).
- **CI**: one line in `run-tests-rustfft`: `cargo clippy -p tetraleaf-validate --all-targets
  --features rustfft -- -D warnings`, so that the harness with the feature keeps building; root
  `CLAUDE.md`'s description of the job.
- No library code, no convention, no dependency or pin, no `tetraleaf-plan` interface change.

## Tests that prove it
- `tests/helmholtz_accuracy.rs` (`check_diagonal`): the existing interpolated-build or `auto` row
  reports `phi_step == Product` and `phi_steps` `Some(Product)` exactly on the levels with
  diagonal classes, `None` elsewhere. With the feature (`#[cfg(feature = "rustfft")]`): the same
  case under `RustFft`: every diagonal level `Some(RustFft)`, the errors within 1.01× of the
  product row's (CH4.6a1's FMM bound), the classes and pairs the product row's (the φ step changes
  no choice: the time rule's costs are the product's). Run at one rank in debug, and with the
  feature at 1 and 2 ranks by hand.
- The example: `--phi-step rustfft` refused without the feature (by hand), the smoke run of the
  runs script (`smoke` stage) at W1 kD = 30 checks the tables' columns.

## Checks
`/verify` (default, strict workspace, the MPI list at 2 and 4 ranks, `cargo deny`), plus
`cargo clippy -p tetraleaf-validate --all-targets --features rustfft -- -D warnings`,
`RUST_MIN_STACK=8388608 cargo test -p tetraleaf-validate --features rustfft --test
helmholtz_accuracy` (1 rank; and 2 ranks by hand), `cargo tree --workspace -e normal -i rustfft`
(none), `python3 tools/helmholtz/common.py selftest`-style self-check of the tables script on the
smoke output.

## Machine runs
Release builds with `--features rustfft`, BLAS threads 1, one setting per launch (review
2026-11-24, B12), §9.1's complex charges and samples, gradients, two evaluations (the rows the
second's), the library's defaults named (`--m2l auto-pw --diagonal auto --diagonal-build
interpolated --near-v cost --near-wx cost`), `--levels`.

- Stage `levels` (1 × cores): W1 at kD = 250, 500, 1000 at 10⁻⁶, W1 at kD = 500 at 10⁻³ and
  10⁻⁸ (D1's edges), W3 at N = 10⁶, kD = 500 at 10⁻⁶; per case six launches in the order
  product, rustfft, rustfft, product (ABBA), then `--diagonal-costs form` with the product and
  with rustfft (every level the bound reaches on the form, for the time rule's question).
- Stage `waves` (added after Kathleen's `levels`, where every level in the rule's range already
  took the form, so the time rule's question had no level to decide): the cases where CH4.3b2c
  found the waves kept against the form, W1 at kD = 100, 200 at 10⁻⁷ and 10⁻⁸, W3 at N = 10⁶,
  kD = 200 at 10⁻⁶ and 10⁻⁸ and kD = 500 at 10⁻⁸; per case the defaults with each φ step, then
  `--diagonal-costs form` with each.
- Stage `ffts`: S9's `times rustfft product <ps>`, `times rustfft smooth <ps>` and `times product
  product <ps>` at the diagonal levels' degrees from `levels` (one launch each, one thread).
- **locust** (1 × 72, `tools/gh200/remote.sh`, after `nvidia-smi` and `uptime`), and **Kathleen**
  (1 × 40, `singlenode`, `--exclude=node-c11b-004,node-c11b-013`, ≤ 6 h).
- Decides: (1) the gain per level and per run of rustfft on the product's grid; (2) the smooth
  grid's price per level, thread-seconds:
  S → F + B_s (t_s^smooth − t_s^product) + pairs × N_φ^smooth / N_φ + F → R + B_t (t_a^smooth −
  t_a^product) + RCR pairs, with B_s, B_t the level's conversions and t the per-box FFT times
  from `ffts`; (3) the time rule: on each level the default did not put on the form, the forced
  form with rustfft against the option it ran, and the model's form with S → F and F → R scaled
  by the measured rustfft ÷ product ratio against its RCR and waves.
- Raw outputs to `~/projects/fmm-archive/<machine>/runs/ch4-6a2/`; the report
  `tetraleaf-validate/results/helmholtz-ch4-6a2.md`.

## Constraints and risks
- Constraint 1: no Laplace path is touched (the harness and its scripts; one CI line).
- Constraint 9: nothing measured at run time; the time rule's question is answered from these
  runs, and any FFT variant of `DiagonalCosts` would be a follow-up task's.
- Risk: Kathleen's run-order drift (CH4.3b2c §6): one setting per launch and ABBA; the two
  product launches give the spread.
- Risk: memory at kD = 1000 on locust (W1, 1 × 72) — CH4.5 ran it.

## Pending runs
None. locust's `levels`, `waves` and `ffts` ran on 2026-11-26 19:35–19:44, started by a gate on the
M6 once the node was quiet (no other process of this account above 20 % of a core, a one-minute
load below 2; a vLLM server of another user, 0.3–0.6 core, ran throughout); the first gate's
condition was wrong and started one launch at 18:18 on a busy node, which was stopped and its
outputs set aside (`/data/ucahtbe/fmm/runs/ch4-6a2-discarded-not-quiet`). The M6: not run (the user,
2026-11-26: "Skip the M6"; question 3).

## Results (the report `tetraleaf-validate/results/helmholtz-ch4-6a2.md`)
Jobs 673075 (`smoke`, `levels`), 673077 (`ffts`), 673078 (`waves`), node-c11b-009, 1 × 40. rustfft's
conversions 0.23–0.60× the product's per conversion; the diagonal levels' M2L 0.63–0.93×; the
evaluation 0.89–0.95×, build plus one evaluation 0.94–0.97×; errors 1.0000×, checksums agree per
setting. The smooth grid: 0.98–1.03× per level, 1.000–1.004× per evaluation (0.991–1.001× at the
FFT's upper bound). The time rule with rustfft's conversions flips one level of 13 cases, 1.6 % of
that level, within noise. locust (1 × 72): the conversions 0.31–0.71×, the levels 0.60–0.95×, the
evaluation 0.91–0.99×, errors 0.993–1.000×; the smooth grid 1.000–1.005× per evaluation; the time
rule moves two levels (a 3 % loss and a gain where the form beats the waves with the product too).

## Checks run (the M6)
`/verify`: fmt; Clippy default, examples, `--workspace --all-targets`, and tetraleaf-validate with
`--features rustfft`; `cargo doc` with `-D warnings`; `cargo check --workspace --all-targets`;
`cargo deny --locked check licenses`; `cargo tree --workspace -e normal -i rustfft` (none); the
one-rank debug run of `tests/helmholtz_accuracy.rs` at opt-level 0, 400 s at load 11 (437.8 s
recorded at load 19–24); at opt-level 1 `cargo test` (1,402 passed) and `cargo test --workspace`
(1,446 passed); `cargo test -p tetraleaf-validate --features rustfft` (92 passed); the harness test
with the feature at 2 ranks by hand (passed, 17 s). Not run: the MPI list (no change to
tetraleaf-exec, -plan or -octree), Laplace benchmarks (no Laplace path changed).

## Self-review
A subagent (the plan, §2, §10, the branch): no code defect. Fixed: the report's counts (four levels
kept on the waves, 13 cases), the non-smooth degrees (ten of fifteen, not p ≥ 116 alone), the
smooth-grid caveat (the level's FFT runs several times S9's per box: the price given at both ends),
the flip's gain against the same setting's waves (1.6 %, within noise), W3's gain within its
spread, two statements of §2, the plan's pending runs, the misplaced doc comment of
`check_diagonal`, the GF/s of rustfft rows called product-equivalent, the stage list in
`tetraleaf-validate/CLAUDE.md`, a guard in the tables script's choice, and the quoted output line.

## Open questions
Asked before any run, 2026-11-26: the M6's timing launches; the user: "Skip the M6".

For sign-off on #266, signed off 2026-11-26, the user: "Q1-3 as recommended":
1. Keep the product's grid for rustfft and close CH4.6a3 unrun.
2. No FFT variant of `DiagonalCosts`; the product's costs stay the rule's under either φ step.
3. CH4.6a's acceptance without the M6's times.
Proposal §6.7, §10 (CH4.6a) and §12 (D7) amended inline ("amended at CH4.6a2"), Appendix A.
