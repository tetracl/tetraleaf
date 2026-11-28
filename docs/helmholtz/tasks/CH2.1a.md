# CH2.1a Yukawa's RCR and dense M2L tables

Started 2026-11-28 on `tbetcke/h-ch2-1a`. Roadmap: CH2.1a, the first part of CH2.1 (item 15 of the
priority list, Yukawa on the host; added on 2026-11-24 at the H4 gate). Proposal: §10 CH2.1; §2
(constraints 1, 5, 6, 9); §5.5 (Yukawa: real coefficients, real tables, (p + 1)² reals per
expansion); §6.2 (the Yukawa row of the sharing table); §6.3 (one operator per kind, D4).
CONVENTIONS §3.15 ("Yukawa's real form", "Coaxial translations", "Level units").

## Size: CH2.1 split

CH2.1 as written holds three pieces of work in three crates: the tables of Yukawa's real path in
`tetraleaf-tables` (refused there today: `HelmholtzTables::build_by` asserts Re k̂₀ > 0), the
operator and the FMM in `tetraleaf-exec` (refused there: `OperatorError::Yukawa`) with the runs on
W1, W3 and W4 over κD from 10⁻² to 10³ and the D4 question, and the Also from CH4.9a2b, the batched
P2M for Yukawa. That is the skill's rule for a split (more than two crates, an operator plus runs),
along the layers:

- **CH2.1a** (this plan): the tables, tested against `tetraleaf_ref::helmholtz::rcr`'s `_yukawa` forms.
- **CH2.1b**: the operator and the FMM on them, the acceptance's runs, D4 for sign-off.
- **CH2.1c**: the batched P2M.

The three lines, with *Done when* and *Needs*, are in the roadmap (committed with this plan); D4's
"Prepared by" names CH2.1b.

## Goal
Build the per-level RCR and dense M2L tables of Yukawa's real path (k = iκ, real expansions) in
`tetraleaf_tables::helmholtz`, so that CH2.1b's operator runs Rotation and Dense on them.

## Acceptance (proposal §10 CH2.1, verbatim, for the parent)
"Real `HelmholtzOperator` path on the drivers of CH2.0 (or on copies); Rotation and Dense M2L" —
"FMM against the direct sum on W1, W3 and W4 for κ D from 10⁻² to 10³; deterministic".

CH2.1a's part (the roadmap's *Done when*): M2M and L2L on the 8 octants and M2L on the 316 offsets
against `rcr`'s `_yukawa` forms within the real and complex tables' bound (2e-13 max(1, |k̂| a) in
the FMM-weighted measure) from κ̂ a ≈ 10⁻³ to the strongly damped levels, rectangular and raised
degrees included; Dense against Rotation pair by pair; vectors, stored blocks, contributions and the
build in parts bit for bit; real and complex k's tables bit for bit as before (by hash).

## Read first
- `tetraleaf-ref/src/helmholtz/rcr.rs`: `m2m_yukawa`, `l2l_yukawa`, `m2l_yukawa` (the real tables
  i^{n−n′} Ê, i^{n′−n} Ê, −i^{n+n′} Ê with 4πκ e^{2β − κ|t|}; `real_part_times_i_power`,
  `power_ratio`), and the module docs' "Stored form" and "Tested range" (strong damping: an M2L's
  output relative error grows, round-off relative to the input).
- `tetraleaf-tables/src/helmholtz.rs`: `build_by` (the parts `Up`, `Down`, `M2l`), `Shift::new`,
  `translate` (the real step on each plane), `m2l`, `m2l_contribution(_len)`.
- `tetraleaf-tables/src/helmholtz/dense.rs` (the 2n × n layout, `apply`), `batched.rs`,
  `classes.rs` (both read Dense's layout).
- `tetraleaf-tables/tests/helmholtz/{main,dense,contribution,parts}.rs` (the measure, the
  comparisons with `rcr`, `print_the_hashes`).
- `tetraleaf-ref/tests/helmholtz/translation_yukawa.rs` (the two ways at κ a up to 300).

## Changes
`tetraleaf-tables` only, in its Helmholtz module (no Laplace table, layout or cache touched):

- `helmholtz.rs`:
  - `build_by` accepts k̂₀ = iκ̂₀ (Re = 0, Im > 0; still refuses Re < 0, Im < 0, zero, non-finite).
    Yukawa's tables are real (`real_shifts` true): `Up` the entries
    `real_part_times_i_power(n − n′, Ê)` × (s_c/s_P)ⁿ s_P^{n−n′}, `Down`
    `real_part_times_i_power(n′ − n, Ê)` × (s_c/s_P)^{n′} s_P^{n′−n} × ½, `M2l`
    −`real_part_times_i_power(n + n′, Ê)` × 4πκ e^{2β − κ|t|}, all formed in f64 as `rcr`'s
    `_yukawa` forms form them and rounded to T once. A private `real_part_times_i_power` restates
    `rcr`'s (a stored-form factor, the module's exception of `tetraleaf-tables/CLAUDE.md`).
  - `HelmholtzTables::is_yukawa()`; `local_planes()` (1 for Yukawa, 2 otherwise): the planes of a
    local. `m2l` and `m2l_contribution` write `multipole_planes` planes for Yukawa (a real step,
    `translate`'s assertion), two otherwise; `m2l_contribution_len` counts the local's planes.
  - Module docs: Yukawa's stored form in the table of "The stored form, folded in"; the panics.
- `helmholtz/dense.rs`: for Yukawa one real n × n matrix per offset (column j Rotation's M2L of e_j,
  one plane), `matrix` of n × n; `m2l` takes and writes one plane per vector; `apply`'s one-plane
  chain over the local's rows (2n complex, n Yukawa: the same loop, the complex path's bits).
  `level_bytes_for(p, local_planes)` and `select_for(degrees, budget, local_planes)`, the
  existing `level_bytes` and `select` their two-plane case (unchanged signatures and values).
- `helmholtz/batched.rs`, `helmholtz/classes.rs`: refuse Yukawa with a panic naming CH2.1 (their
  layout is Dense's complex one; Yukawa's Batched and Classes are left to a later task, named in
  the roadmap at CH2.1b's line if CH2.1b's Auto wants them).
- `tetraleaf-tables/CLAUDE.md`: the Yukawa tables in the `helmholtz` and `dense` rules.

No public API removed or changed in meaning; no dependency; no convention change (CONVENTIONS §3.15
already defines Yukawa's real form); `tetraleaf-plan` untouched; no report script.

## Tests that prove it
New `tetraleaf-tables/tests/helmholtz/yukawa.rs` (in `main.rs`'s module list):

1. **Against `rcr`'s `_yukawa` forms**, the FMM-weighted measure of `main.rs` (`weights`,
   `compare`), within `BOUND` × max(1, κ̂ a) (BOUND = 2e-13, the real and complex tables'):
   M2M and L2L on all 8 octants between every level pair, M2L on all 316 offsets, one plane, at
   k̂₀ = i κ̂₀ with κ̂₀ ∈ {0.01, 2, 40} over levels 0–4 (κ̂ a from about 10⁻³ on level 4 to about 17
   on level 2), degrees 0–12 as `main.rs` takes them; rectangular degrees (p_{l+1} ≠ p_l) and
   raised M2L degrees (`build_with_m2l_degrees`, the first degrees read and filled, the rest NaN
   or untouched as `main.rs` checks). Strong damping (κ̂ a ≈ 70 and 170 on level 2: M2L outputs
   e^{−0.54 κ̂ a} of their input): measured relative to the output and to the input; if the output
   measure exceeds the bound there (`rcr`'s "Tested range" expects round-off relative to the input
   only), the test holds the input-relative error to the bound and prints the output's, and the
   report says so (CH2.3 prunes such pairs, proposal §5.5).
2. **Bit for bit**: several vectors against one; stored M2L blocks against per-call; the
   contribution then `add_m2l_contribution` against `m2l`; the build on a reversed three-thread
   map against the serial build (`parts.rs`' map), through every application.
3. **Dense against Rotation**: every offset of levels 2–4, each column Rotation's M2L of e_j to
   rounding (dense.rs's measure and bound), a pair's application against Rotation's on a random
   multipole, two vectors each the one-vector call, f32 the f64 matrices rounded, `level_bytes_for`
   and `select_for` with one local plane, a zero budget nothing.
4. **Refusals**: Batched and Classes panic on Yukawa; `build` still panics for Re < 0, Im < 0, zero
   and non-finite k̂₀; `m2l` with two multipole planes into one local plane panics.
5. **Real and complex k unchanged**: `parts.rs`' ignored `print_the_hashes` on `main` and on the
   branch, the same hashes; every existing test of the crate unchanged and passing.

## Checks
`/verify` (fmt; Clippy default and workspace; `cargo test -p tetraleaf-tables`, the ignored tests
of this change in release; docs; the root and workspace tests). No MPI-crate code changes, so the
MPI list is not needed for this change, but `/verify` decides by what changed.

## Machine runs
None: tables against the reference on the M6, no timing.

## Constraints and risks
- Constraint 1 (Laplace): no Laplace code, table, layout or cache touched; `tetraleaf-tables`'
  Laplace modules unchanged.
- Constraint 5: tables built in f64 and rounded once, as the real and complex tables.
- Constraint 6: every fast path against `rcr`'s `_yukawa` forms, and Dense against Rotation.
- Constraint 9: Yukawa's tables are a function of k̂₀ and the degrees.
- Risk: strong damping's cancellation (above); low κ: the power factors overflow where
  p log₁₀(1/s_P) ≳ 300 (`rcr`'s limit, `build`'s documented panic; CH2.3 refuses it in the
  builder, review 2026-11-19 M13). The tests stay inside it.
- The real and complex tables must keep their bits: the Yukawa branches are separate entries of
  the same closures; checked by hash.

## Pending runs
None.

## Open questions
None yet.
