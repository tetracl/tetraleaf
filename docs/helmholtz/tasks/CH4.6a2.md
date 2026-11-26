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
(filled in step 6)

## Open questions
Asked before any run, 2026-11-26: the M6's timing launches; the user: "Skip the M6".
