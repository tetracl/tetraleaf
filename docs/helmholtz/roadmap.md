# Helmholtz FMM roadmap

The tasks of `docs/design/helmholtz-fmm-proposal.md` (draft 0.8), Section 10, in the order
`/next-task` takes them. The proposal is the full record: each task's *Done when* condenses
its acceptance criterion in the §10 table, and where the two differ the proposal wins.
Progress is logged in `docs/helmholtz/progress.md`; each task's plan is
`docs/helmholtz/tasks/<ID>.md`; each milestone's gate report is
`docs/helmholtz/gates/<milestone>.md`.

Status (2026-11-27): **H1, H3, H4 and H6 closed** (`docs/helmholtz/gates/H1.md`, 2026-11-19;
`gates/H3.md`, 2026-11-02; `gates/H4.md`, 2026-11-25; `gates/H6.md`, 2026-11-12). Item 14, the leaf
expansions: CH4.9a1 and CH4.9a2a merged; CH4.9a2b (P2M by the batched kernel under every ISA, 0.25–0.31×
the reference's per run on Kathleen and 0.15–0.27× on locust, the FMM's errors 0.993–1.000×) done in
its PR, its ISA question signed off as recommended (carried to CH4.9b). Next CH4.9b (L2P), then H2 and H5 (items 15–16), CH2.1 first. Eight H4 tasks stay open off the list, and
CH4.9c. The project was renamed tetraleaf on 2026-11-21. HR, H0, H2, H5, H7, H8 and H9 are open; of
H0, CH0.4, CH0.6 and CH0.9 are done, and of H7, CH7.5a–c. The history of this line until H6.G is in
`progress.md` and in this file's history (`git show 4b5b37f:docs/helmholtz/roadmap.md`); until H1.G
in `git show 8475009:docs/helmholtz/roadmap.md`; until the H4 gate in `git show
4d28d41:docs/helmholtz/roadmap.md`; until its second run in `git show 03bc2bd:docs/helmholtz/roadmap.md`;
until CHA.10a in `git show 6bad3c2:docs/helmholtz/roadmap.md`; until CH0.9 in `git show
dc69cd4:docs/helmholtz/roadmap.md`; until CH4.6a1 in `git show 34e3caf:docs/helmholtz/roadmap.md`;
until CH4.6a2 in `git show 7c8177d:docs/helmholtz/roadmap.md`; until CH4.9a1 in `git show
9758028:docs/helmholtz/roadmap.md`; until CH4.9a2a in `git show a62c5b9:docs/helmholtz/roadmap.md`;
until CH4.9a2b in `git show 9401720:docs/helmholtz/roadmap.md`.

## How `/next-task` reads this file

- `- [ ]` is open, `- [x]` done. The tick is part of the task's own pull request, so a task
  counts as done once that PR is merged into `main`.
- Work in flight is not marked here. A task with a branch `tbetcke/h-<id>` or an open PR is
  in progress: a **draft** PR waits for machine runs, a sign-off or an open question (its
  plan file says which); a **ready** PR waits for the user's merge.
- A task is **runnable** when it is open, has no branch or open PR, and everything in its
  *Needs* is done: the tasks ticked on `main`, the decisions signed off in the table below.
  Every task also needs **D0**.
- `/next-task` takes the first runnable task of the **priority list** below, in its order.
  If the first open task of the list waits only for decisions, it asks the user for them
  first. If no task of the list is runnable (the rest wait for a merge, a machine run or a
  sign-off), it says what the list waits for and stops; it takes a task off the list only
  when the user names one. With the list empty or all ticked, it takes the first runnable
  task in the order of this file, as before.
- A task too big for one reviewable PR is split here into lettered sub-tasks (CH3.1a,
  CH3.1b), each with its own *Done when* and *Needs*, before it starts.
- Each milestone is mirrored as a GitHub issue (labels `enhancement`, `helmholtz`), ticked
  and closed from this file on every push to `main` (tools/roadmap/README.md). Change this
  file, never the issues; a task's PR ends with `Part of #<n>`, its milestone's issue.
- *Machines* names the machines its runs need besides the development Mac (the M6 since
  2026-10-31, the M3 Max before). *Sign-off* marks a task whose PR ends with questions for
  the user; the tasks that need it wait for the answer.
- "Review <date>, <finding>" cites a finding of a project review in `docs/reviews/`
  (`docs/reviews/2026-11-12-H6.md`: M5, A2, L1, …). The architecture tasks that come out of a
  review are the section "HA" at the end of this file, with IDs CHA.<n>; it is not a
  milestone, has no gate and no issue, so a CHA task's PR has no `Part of` line.

## Priority list: the path to H4, then H2 and H5

Set by the user on 2026-11-12 at the H6 gate ("As proposed (Recommended)";
`docs/helmholtz/gates/H6.md` §5 and question 4), D6's order: with the low-frequency M2L done,
the far field at kD = 500 is mostly the high-frequency levels' RCR. It replaces the list of
2026-11-02, the path to H6, which finished with CH6.1h and H6.G. H4's follow-ups first, then
the signatures' convention and trusted path, then H4's tasks, then the FFTW track (corrected at
the H4 gate: CH4.6 needs H4, H4.G does not need CH4.6; `gates/H4.md` §4, R1). Amended on 2026-11-12 by the triage of the review after the H6 gate
(`docs/reviews/2026-11-12-H6.md` §3; questions 1 and 3 of #209, the user, 2026-11-12, on #209: "I agree", to the
recommendations): CH4.0e and CHA.1 to CHA.3 before CH1.3d. Amended on 2026-11-16 by CH4.0e's
sign-off (#220's question 1, the user: "Both questions as recommended."): H4 goes ahead, for kD
beyond D1 and the high-frequency classes CH4.3 calibrates, after CH4.0h and CH4.0i; the reason
for H4 above is replaced by CH4.0e's report §3 (inside D1 the high-frequency RCR levels are
14–31 % of build plus one evaluation, on W1 at kD = 500 alone). Amended on 2026-11-16 by
CH4.0h's sign-off (#221's question 2, the user: "Please sign off questions as recommended."): CH4.0j beside CH4.0i in
item 6, CH4.0k off the list. Amended on 2026-11-19 at the H1 gate (`docs/helmholtz/gates/H1.md`
§7, question 5; the user on #230: "CHA.11 + measure first (Recommended)"): CHA.11 as item 9,
and CH4.1's first step a measurement of the product conversions with a stop to ask (review
2026-11-19, A5: with the product φ step a high-frequency level gains about 2.3× at most).
Confirmed on 2026-11-19 by CH4.1a's sign-off (#234's question 1, the user: "Q1-3 as recommended. Merge when green"):
the measured product leaves an FFT about 1–2 % of the run, so items 10 and 11 keep their order.
Amended on 2026-11-24 at the H4 gate (`docs/helmholtz/gates/H4.md` §7, question 5; the user on
#258: "FFT, then H2 and H5 (Recommended)"): after the FFT track the Helmholtz family on the host,
H2 then H5; CHA.7a before H8; CH4.7, CH4.3f and CH4.0k stay off the list. Amended on the same day
by the gate's run (`gates/H4.md` §9; the user: "Keep H4 open, CH4.8 (Recommended)"): CH4.8 and
H4.G as item 11, before the FFT track. Amended on 2026-11-25 at the H4 gate's second run
(`gates/H4.md` §0.8, question 4; the user on #262: "CHA.10a, then FFT, H2, H5 (Recommended)"):
CHA.10a as item 12, before the FFT track. Amended on 2026-11-26 after CH4.6a2 (the user: "Yes. Add P2M
and L2P to the priority list as tasks"): CH4.9a and CH4.9b as item 14, before H2 (P2M and L2P were
13–31 % of an evaluation under the defaults on Kathleen, job 673083, against 8–17 % for M2M and L2L
on the diagonal levels).

1. CH4.0 The M2L rule's sparse-level term
2. CH4.0b The worst-box concentration
3. CH4.0c The degree rule beyond its measured range
4. CH4.0d1, CH4.0d1b, CH4.0d2, CH4.0d3 The Helmholtz leaf size (added on 2026-11-12, the user:
   "After CH4.0c (Recommended)"; it changes the trees that CH4.0's rule and CH4.4's and
   CH4.5's timings run on; split on 2026-11-13): the empty boxes, their build cost (added at
   CH4.0d1, next after it: #216's question 2, signed off 2026-11-13), the V pairs of few
   points, the sweep
5. CH4.0e The comparison under the default, and H4's expected gain (review A1 and B2; its
   sign-off is H4's go-ahead and order)
6. CH4.0h The Helmholtz build by part, and parallel; CH4.0i The coarse RCR levels (added on
   2026-11-16 by CH4.0e's sign-off); CH4.0j The exponential M2L's tables within a level in
   parts (added on 2026-11-16 by CH4.0h's sign-off)
7. CHA.1, CHA.2, CHA.3 The spikes out of the workspace, `unsafe` enforced, the superseded
   configuration removed
8. CH1.3d1, CH1.3d2 The diagonal form in extended precision (CONVENTIONS §3.15 "Signatures",
   CH4.2's trusted path; split on 2026-11-19)
9. CHA.11 The Helmholtz wave levels in one place (added on 2026-11-19 at the H1 gate, the
   user on #230: "CHA.11 + measure first (Recommended)"; review 2026-11-19, L1)
10. CH4.1 (CH4.1a, CH4.1b; split on 2026-11-19), CH4.2 (CH4.2a, CH4.2b1, CH4.2b2; split on 2026-11-19 and 2026-11-20), CH4.2c (CH4.2c1, CH4.2c2; added at CH4.2b2, its question 2, split on 2026-11-20), CH4.3 (CH4.3a–CH4.3e, CH4.3b as CH4.3b1 and CH4.3b2a–d; split on 2026-11-20; CH4.3b2d added at CH4.3b2b, its question 3; CH4.3c as CH4.3c1 and CH4.3c2, split on 2026-11-21; CH4.3c3 added at CH4.3c2, its question 3; CH4.3d as CH4.3d1 and CH4.3d2, split on 2026-11-22), CH4.4 (CH4.4a, CH4.4b; split on 2026-11-22), CH4.5 (CH4.1a first
    measures the product conversions against RCR and stops to ask before CH4.2 whether item 11
    moves ahead: review 2026-11-19, A5)
11. CH4.8 The constant density on W2 (CH4.8a, CH4.8b; split on 2026-11-24), then H4.G again
    (added on 2026-11-24 at the H4 gate; H4.G closed on 2026-11-25)
12. CHA.10a The report scripts' shared module (added on 2026-11-25 at the H4 gate's second run)
13. CH0.9 S9, CH4.6a the optional rustfft back end (CH4.6a1–a3; split on 2026-11-26). D7, signed off at CH0.9 (the user on #264,
    2026-11-26: "Q1-3 as recommended"), moved CHR.3 and CH4.6b, the FFTW route, off the list
14. CH4.9a (CH4.9a1, CH4.9a2a, CH4.9a2b; split on 2026-11-26 and 2026-11-27), CH4.9b (CH4.9b1, CH4.9b2; split on 2026-11-27) The Helmholtz leaf expansions P2M and
    L2P by SIMD kernels (added on 2026-11-26 after CH4.6a2, the user: "Yes. Add P2M and L2P to the
    priority list as tasks")
15. CH2.1, CH2.2, CH2.3, CH2.5 and H2.G Yukawa on the host (added on 2026-11-24 at the H4 gate;
    CH2.1 and CH2.2 wait for D4, which CH2.1 prepares)
16. CH5.1, CH5.2, CH5.3 and H5.G Complex wavenumber (added on 2026-11-24 at the H4 gate)

Off the list, in file order: CH0.8b, CH0.5 and H0.G; CHR.3 and CH4.6b (moved off at CH0.9 by D7, the user on #264, 2026-11-26: "Q1-3 as recommended"); CH4.0f (the user, 2026-11-12, on #211);
CH4.7 (added at CH4.5, its question 4; the user on #257, 2026-11-24: "Accept all 4 decisions as recommended."); CH4.3f (added at CH4.3b2c, off the list until CH4.5, kept off at CH4.5, its question 3, by the user on #257, 2026-11-24: "Accept all 4 decisions as recommended."; the user on #251, 2026-11-22: "Q1-5 as recommended"); CH4.2d (added at CH4.2b2, off the list until CH4.5, kept off at CH4.5, its question 2, by the user on #257, 2026-11-24: "Accept all 4 decisions as recommended."; the user on #239, 2026-11-20: "All 4 questions as recommended.");
CH4.0k (added at CH4.0h; the user, 2026-11-16, on #221: off the list; kept off at CH4.0j, on
#223); CH4.0l (added at CH4.0i;
the user, 2026-11-16, on #222: off the list);
CH4.0g (after CH4.0d1, through its *Needs*); HR.G; H7 (CH7.5d now has D8 and waits for CH7.1; CH7.5f before CH7.5e), H8, H9 (CH9.1,
after every other milestone; CH9.0, the dipoles, added at the H1 gate); CHA.4 to CHA.10, CHA.12,
CHA.13 and CHA.14 (CHA.7 before CH8.1, through its *Needs*, its first letter CHA.7a before H8 by the H4
gate; CHA.6 recommended before CH7.5d).

## Decisions

The proposal's Section 12, plus D0. A decision is recorded here when the user takes it:
`signed off <YYYY-MM-DD>: <the answer as given>`.

| ID | Decision | Proposal's recommendation | Prepared by | Needed by | Status |
| --- | --- | --- | --- | --- | --- |
| D0 | Draft 0.8 of the proposal as the plan of record | as written; later changes by sign-off | review of draft 0.7 | every task | signed off 2026-10-26: as written |
| D1 | Target kD range and tolerances | for example kD ≤ 500 at 10⁻³ to 10⁻⁸ (proposal §12) | – | CH0.8b, CH3.0 | signed off 2026-10-31: kD ≤ 500 at 10⁻³ to 10⁻⁸, as recommended |
| D2 | Placement: modules in the layer crates (A) or separate crates (B) | A (§6.1) | – | CH1.2 | signed off 2026-10-29: A, modules in the layer crates, as recommended |
| D3 | f32 scope | f64 only until H7; then f32 for low-frequency configurations that pass an accuracy contract; high-frequency levels in f32 only if S6 finds stable classes | CH0.6 | CH7.1 | signed off 2026-11-06 (CH0.6, question 1; the user: "i agree with the recommendations"): as recommended, high-frequency levels may run the diagonal form in f32 at ε ≥ 10⁻⁴ for real k and weak damping, each class by the bound in u = 2⁻²⁴ with C_T = 0.18 and θ = 1/2, failing classes on RCR; f64 below 10⁻⁴ and for η ≈ 1; CH7.1 settles the f32 RCR and the low-frequency levels |
| D4 | Yukawa: its own operator type, or the wavenumber kind as a type parameter | one `HelmholtzOperator<T, K>` if CH2.1 shows no cost to the complex path | CH2.1 | CH2.2 | open |
| D5 | Source types of the first release | real and complex charges, potentials and gradients; dipoles once C6.4 exists for Laplace (note: Laplace dipoles were merged in Phase 6 T8 and T8d, so the condition is met) | – | CH1.1 | signed off 2026-10-29: charges now, dipoles later; §3.15 fixes real and complex charges with potentials and gradients, and dipoles get a later, separate task; that task is CH9.0, added at the H1 gate (the user, 2026-11-19, on #230: "Add CH9.0 (Recommended)") |
| D6 | Order of H4 and H6 | H6 first | H3.G | CH6.1, CH4.1 | signed off 2026-11-02 (H3 gate, question 2): "i agree with the recommendations", H6 first, as recommended; on the gate report's reason (§3.3): the calibrated degrees cut the low-frequency share of the far field at kD = 500 to about half (an estimate), and the FMM is slowest against FMM3D on low-frequency trees at tight ε, which H6 targets |
| D7 | Host back ends of the φ step | the product as default and reference; FFTW behind the optional `fftw` feature; rustfft if S9 finds it deterministic per box and faster, decided per machine | CH0.9 | CH4.6a, CHR.3, CH4.6b | signed off 2026-11-26: "Q1-3 as recommended" (#264; `spikes/fft-plans/SPIKE_REPORT.md` §6): the product the default and reference on every machine; rustfft the one optional FFT back end (CH4.6a), with a check across ranks and the M6's times; FFTW (CHR.3, CH4.6b) off the priority list, kept in the roadmap. Amended 2026-11-26 (CH4.6a2, questions 1–3; the user on #266: "Q1-3 as recommended"): rustfft on the product's grid (CH4.6a3 closed unrun), the time rule's costs the product's under either φ step, the M6's times not taken (the user: "Skip the M6"); proposal §6.7, §10 and §12 amended |
| D8 | Exponential M2L as a default on low-frequency levels | after CH6.2 | CH6.2h | CH7.5d | signed off 2026-11-08 (CH6.2h, questions 1–5; the user: "1–5 as recommended"): `Auto` takes the exponential M2L by default where p_l ≥ d_pw (Batched's crossover + 1; 12 on AVX-512 with two planes), CH6.1e's rule gives a quadrature (ε ≥ 10⁻⁶) and Im κ_l < 5π, f64 only; no boxes term (Batched's few-box levels to CH4.0); proposal §5.3 and §6.4 amended. Amended 2026-11-11 (CH6.1g, question 2; the user: "Lift for 10⁻⁸–10⁻⁷ only (Recommended)"): `Auto` also takes it at 10⁻⁸ ≤ ε < 10⁻⁷, where the rule runs iprec 3 at the raised degree max(p_l, 32), and not at 10⁻⁷ ≤ ε < 10⁻⁶ until CH6.1h; proposal §5.3 amended. Amended 2026-11-11 (CH6.1h, question 2; the user: "Lift it (Recommended)"): `Auto` also takes it at 10⁻⁷ ≤ ε < 10⁻⁶, with the raised degree on the M2L alone there (`PlaneWaveRaise::Rule`); proposal §5.3 amended. Amended 2026-11-22 (CH4.3b2c, question 3; the user on #251, 2026-11-22: "Q1-5 as recommended"): with the diagonal form's time rule the builder's default, a level of D8's rule in its range (k a ≥ 2.8) keeps the waves unless the form beats them at the level's resolving degree, then runs the cheaper of the form and RCR; proposal §5.3, §6.4 and §12 amended |
| D9 | What "within noise" means for Laplace's benchmarks (constraint 1) | the run-to-run spread of the Laplace baseline, repeated on each machine before the change, per benchmark; plus the C3S.4 throughput gate | – | CHR.2, CH2.0 | signed off 2026-10-26: as recommended, the measured spread of the baseline per benchmark and machine, plus the C3S.4 gate. Clarified 2026-10-26 (CHR.2): one-sided, slower than `main` by more than the spread fails, faster passes when the outputs are bit for bit identical; from CH2.0 on five runs per side in the order A B B A A B B A A B, the spread max − min of `main`'s five run medians, never pooled with the branch's. Amended 2026-10-31 (CH2.0, the user: "as recommended"): a failing row may be repeated once and the repeat counts, with both results reported. Amended 2026-10-31 (the machine move, the user: "Amend: M6 replaces M3 Max"): the M6 (`tools/m6/machine.md`) replaces the M3 Max as a machine of record, its C3S.4 reference `spikes/p2p-simd/results-m6.md` |
| D10 | cubek-fft on the device | after S5: adopt, contribute upstream first, or keep the transform as a GEMM | CH0.5 | CH7.4 | open |
| D11 | Maximum leaf width c_leaf | from S8: the largest c_leaf with P2P and leaf costs within 10 % of their best on W1 and W3 (expected 0.5 to 1) | CH0.8a | CH3.0 | signed off 2026-10-28: as recommended (CH0.8a report §9, question 1), c_leaf = 1 by the revised rule, the largest c_leaf with no high-frequency leaf on any workload at 10⁻⁶ and 10⁻⁸, replacing the rule of this row, which is degenerate; a first-pass value that CH0.8b may revise with S6's constants and S7's p; amended at CH4.0d3 (the user on #219, 2026-11-14: "I agree with the recommendation"): the builder's default c_leaf = 2 (`DEFAULT_C_LEAF`), `C_LEAF` = 1 stays the leaf rule's constant, CH4.3 re-decides; amended 2026-11-22 at CH4.3e (the user on #253: "Q1-3 as recommended"): the builder's default c_leaf = 4, leaves on levels of k a ≥ 2.8 accepted (the leaf operators independent of the regime; proposal §3.7, §5.1, §5.7) |
| D12 | The reference workloads | the set of §9.1 | CH0.8a | CH2.1 | signed off 2026-10-28: as recommended (CH0.8a report §9, question 2), the workloads of §9.1 unchanged, with W3 at 10⁶ points and kD = 500 named as the node-budget and leaf-rule-cost case |
| D13 | RLST's `fftw_mkl` back end | – | – | – | withdrawn (RLST ≥ 0.11 has none) |
| D14 | Ownership and schedule of Phase HR | – | – | – | closed (RLST 0.11.0 and 0.12.0) |

## HR: RLST's FFTW interface and its adoption

**Gate:** tetraleaf on RLST 0.12 with Laplace unchanged; the FFTW of every machine recorded;
the optional FFTW wrapper of the φ step behind its feature.

- [x] **CHR.1 RLST's FFTW rework.** Done in RLST 0.11.0 (2026-10-25) and 0.12.0
  (2026-10-26).
- [x] **CHR.2 Move to RLST 0.12.** Report: `tetraleaf-validate/results/helmholtz-chr-2.md`. The workspace pin and `tetraleaf-octree/Cargo.toml`'s own pin to
  0.12.0 together; green-kernels' `[patch]` left on 0.10.2 (it uses `rlst::simd`, removed in
  0.12); `CLAUDE.md`, `tetraleaf-plan/CLAUDE.md` and `docs/design/workspace-structure.md` updated.
  *Done when:* the default members build, lint, test and document; the MPI tests and the
  weekly examples pass; `cargo check --workspace --all-targets` builds the spikes;
  `cargo tree -e normal` lists no fftw-sys; a cargo-deny licence check of the default set is
  in CI; CI's and locust's `cmake` installs checked; `cargo upgrades` no longer reports
  `rlst`; Laplace bit for bit and its benchmarks within noise (D9). Fallback: stay on 0.10.2.
  *Needs:* D9. *Machines:* locust, Kathleen.
- [x] **CHR.4 FFTW on the machines.** Report: `tetraleaf-validate/results/helmholtz-chr-4.md`;
  the probe `tools/fftw/probe.sh`. *Done when:* `pkg-config --modversion fftw3 fftw3f`
  and `rlst::fftw::version()` recorded for the M3 Max, locust (spack's FFTW, not the RPM's
  3.3.8), Kathleen and CI, in `tools/*/machine.md`; RLST 0.12's own
  `cargo test --features fftw_system` passes on each, `fftw_determinism` included.
  Fallback: a machine without FFTW in both precisions uses the product. *Needs:* –.
  *Machines:* locust, Kathleen.
- [ ] **CHR.3 FFTW wrapper of the φ step.** In `tetraleaf-exec` behind
  `fftw = ["rlst/fftw_system"]`: one plan per (level, N_φ), interleaved per-box scratch,
  `Rigor::Estimate`, fingerprints all-gathered at build. *Done when:* against the product to
  round-off; a box's bits the same across batch positions, batch sizes and thread counts;
  ranks with differing fingerprints both use the product; the feature off by default; a CI
  job with the feature (`pkg-config`, the FFTW version report, cargo-deny's entry for
  fftw-sys). *Needs:* CHR.2, CH0.9, D7. Off the priority list since CH0.9 (D7;
  the user on #264, 2026-11-26: "Q1-3 as recommended").
- [ ] **HR.G Gate.** `/milestone-gate HR`. *Needs:* every task above.

## H0: Spikes

**Gate:** S1–S9 answered, each with a `SPIKE_REPORT.md`; the decisions they prepare (D3,
D7, D10, D11, D12) put to the user.

- [x] **CH0.1 S1: rotation blocks at large degree.** Report:
  `spikes/rotation-blocks/SPIKE_REPORT.md`. The Gumerov–Duraiswami (2014) recursion
  in the orthonormal real basis. *Done when:* orthogonal and in agreement with an
  arbitrary-precision reference to p = 1000; in agreement with N-conjugated
  `tetraleaf_math::rotation::blocks` at p ≤ 30; the degree up to which conjugation stays
  accurate; the relation of the β and π − β blocks; time per block set. *Needs:* –.
- [x] **CH0.2 S2: coaxial coefficients.** Report: `spikes/coaxial/SPIKE_REPORT.md`.
  (S|R), (S|S), (R|R) by recursion. *Done when:*
  relative errors against an arbitrary-precision reference for p and k a up to 300, real,
  complex and imaginary k; the imaginary parts of (S|S) and (R|R) zero to round-off for real
  k; the index symmetry that halves the tables. *Needs:* –.
- [x] **CH0.3 S3: radial functions.** Report: `spikes/radial/SPIKE_REPORT.md`.
  j_n, h_n, i_n, k_n for real and complex arguments, with
  the scaling of §3.4 (exponential for Yukawa above κ r ≈ 1). *Done when:* accurate against
  arbitrary precision on every operator's argument and degree ranges; the scaling chosen.
  *Needs:* –.
- [x] **CH0.8a S8: workload census, first pass.** Report:
  `spikes/census/SPIKE_REPORT.md`. Trees of W1–W4 over the kD sweep from
  `tetraleaf-octree` with a prototype maximum leaf width (c_leaf 0.25 to 2), p from the formula; no
  FMM is run. *Done when:* per level the occupied boxes, U, V (by distance class), W and X
  sizes, p_l, regime and per-class choice, the memory model of §7 and the operation counts;
  c_leaf proposed (D11); no high-frequency leaf; the workloads proposed (D12); batch sizes
  for S5 and S9. *Needs:* –. *Sign-off:* D11, D12.
- [x] **CH0.7 S7: truncation under realistic geometry.** Report:
  `spikes/truncation/SPIKE_REPORT.md`, tables in `spikes/truncation/results/summary.md`.
  *Done when:* the P2M → M2M → M2L →
  L2L → L2P chain over two levels and M2M and L2L over three, on all 316 offsets, with points
  filling the boxes (uniform, corners, W1 patches), real, complex and imaginary k, in
  extended precision with the S1 and S2 prototypes; p as a function of (k a, ε, operator)
  meeting ε with a stated margin on every placement. *Needs:* CH0.1, CH0.2, CH0.3.
- [x] **CH0.6 S6: stability of the diagonal form.** Seeded by Appendix B.1. *Done when:* the
  S → F, diagonal, F → R chain in f64 and f32 on all 34 classes against extended-precision
  RCR and the direct sum, with realistic placements, k a 3 to 400, ε 10⁻³ to 10⁻¹⁰, P from the
  formula and from S7, Im k > 0, |t| = 3w as a control; C_T and θ; the classes on RCR per
  (k a, ε); the f32 verdict (D3); whether interpolating T keeps the bound; prototype times.
  Accepted when the bound never under-predicts by more than 2×. *Needs:* CH0.3, CH0.7.
  *Sign-off:* D3. Report `spikes/diagonal/SPIKE_REPORT.md` (results `spikes/diagonal/results/`;
  the raw rows archived in `~/projects/fmm-archive/locust/runs/ch0-6/`). Signed off
  2026-11-06 (questions 1–7 as recommended): D3, the proposal's §3.5 and §3.6 amended.
- [ ] **CH0.8b S8: workload census, final pass.** *Done when:* the first pass repeated with
  the calibrated p of CH3.3 (`DegreeRule::CALIBRATED`, not S7's p, on the low-frequency
  levels; amended at H3.G, proposal §7) and S6's constants, counting the exponential M2L's
  waves buffer and the raised stores of its levels below 10⁻⁶ (amended at H6.G, proposal §7);
  its numbers replace §7's. *Needs:* CH0.8a, CH0.6, CH0.7, D1.
- [x] **CH0.9 S9: FFT plans, bits and host timing.** Report: `spikes/fft-plans/SPIKE_REPORT.md`
  (raw `spikes/fft-plans/results/`). Per-box FFTW and rustfft plans keep a box's bits on the M6,
  locust and Kathleen; fingerprints agree across two Kathleen nodes (job 672796); §6.7's layout
  confirmed, the batch 2P − 1; FFTW 0.04–0.25× and rustfft 0.05–0.32× the product's φ step on the
  smooth grid; M6 times skipped (the user, 2026-11-26: "Skip the M6"). Signed off on #264 (the user,
  2026-11-26: "Q1-3 as recommended"): D7 (CH4.6a rustfft, the FFTW route off the list), §6.7 amended.
- [x] **CH0.4 S4: plane-wave quadratures per k w_l.** Report:
  `spikes/pw-quadrature/SPIKE_REPORT.md`. *Done when:* how FMM3D obtains them;
  a generator or table for this library with its accuracy over the k w_l of low-frequency
  levels. *Needs:* –.
- [ ] **CH0.5 S5: cubek-fft for the φ step.** *Done when:* f64 accuracy against the product
  and an extended-precision DFT; time against the transform as a GEMM for N_φ = 256 to 2048
  at S8's batch sizes, on the GH200 and the CPU runtime; the cost of per-butterfly twiddles
  and power-of-two padding; D10 proposed. *Needs:* CH0.8a. *Machines:* locust.
  *Sign-off:* D10.
- [ ] **H0.G Gate.** `/milestone-gate H0`. *Needs:* every task above.

## H1: Conventions, math, reference

**Gate:** CONVENTIONS §3.15 signed off; `tetraleaf-math` and `tetraleaf-ref` carry the Helmholtz
family, checked against FMM3D's fixtures.

Its ticked lines were condensed at H1.G (2026-11-19, as at H6.G); each task's full line is in
its record or report and its progress entry (`docs/helmholtz/progress.md`), and in this file's
history (`git show 8475009:docs/helmholtz/roadmap.md`).

- [x] **CH1.4 FMM3D convention check and fixtures.** Report:
  `tetraleaf-validate/results/helmholtz-ch1-4.md`; fixtures `tetraleaf-ref/fixtures/helmholtz_fmm3d.json`.
  FMM3D computes the draft's q e^{ikr}/(4πr) for real, complex and imaginary k (one of 16
  candidates, to 2.6e-16); R12: `hfmm3d` reliable for imaginary zk at κD ≤ 30 and eps ≥ 10⁻⁶.
  Signed off 2026-10-28, "1–3 as recommended" (CH2.5 limited to where FMM3D met its eps).
- [x] **CH1.1 CONVENTIONS §3.15.** Record: `docs/helmholtz/tasks/CH1.1.md`. §3.15 and
  `CONVENTION_VERSION = 3`, `tools/fixtures/check_3_15.py` (101 checks); 4πik on M2L, M2P and
  P2L, the stored coefficients with e^{∓β_l}. Signed off 2026-10-29 with D2 and D5.
- CH1.2 `tetraleaf-math`: harmonics, radial functions, rotation blocks. Split on 2026-10-29 into
  CH1.2a–c, one per spike (S1, S3, S2).
- [x] **CH1.2a Orthonormal harmonics and rotation blocks.** Record: `docs/helmholtz/tasks/CH1.2a.md`.
  `tetraleaf_math::orthonormal`; against mpmath to n = 1000 (5.7e-14), blocks orthogonal to
  4.3e-14.
- [x] **CH1.2b Radial functions and scaling.** Record: `docs/helmholtz/tasks/CH1.2b.md`.
  `tetraleaf_math::radial` and `Complex` (`num-complex`, signed off 2026-10-30); against mpmath on
  every operator's range (7.1e-14), Yukawa's relations bit for bit.
- [x] **CH1.2c Coaxial coefficients.** Record: `docs/helmholtz/tasks/CH1.2c.md`; report
  `tetraleaf-validate/results/helmholtz-ch1-2c.md`. `tetraleaf_math::coaxial`, recursion and projection by
  §3.15's rule; one shift table for M2M and L2L (signed off 2026-10-30, "1 and 2 as
  recommended").
- CH1.3 `tetraleaf-ref`: oracle, reference operators, translations, the diagonal form. Split on
  2026-10-30 into CH1.3a–d; CH1.3d deferred after CH0.6 and split on 2026-11-19 into CH1.3d1–d2.
- [x] **CH1.3a Direct-sum oracle.** Record: `docs/helmholtz/tasks/CH1.3a.md`.
  `tetraleaf_ref::helmholtz::p2p`; FMM3D's fixtures to 9.9e-16, k → 0 within §9's tolerances,
  Yukawa two ways bit for bit.
- [x] **CH1.3b Leaf operators.** Record: `docs/helmholtz/tasks/CH1.3b.md`.
  `tetraleaf_ref::helmholtz::{leaf, Level, Workspace}` in §3.15's stored form; series at
  round-off, k → 0 within 0.54 of the bound, finite at Im(k) a = 1000.
- [x] **CH1.3c Translations.** Record: `docs/helmholtz/tasks/CH1.3c.md`.
  `tetraleaf_ref::helmholtz::{rcr, direct}`; RCR against direct 8.4e-15 over the 316 offsets,
  Yukawa two ways 9.5e-15; the damping finding went to CH2.3 at H1.G.
- [x] **CH1.3d1 The diagonal form in tetraleaf-ref.** Record: `docs/helmholtz/tasks/CH1.3d1.md`.
  `tetraleaf_ref::helmholtz::diagonal`, against RCR within 2E on the 34 classes (k a 2.83–54,
  p 7–82). "≤ 2E" and the split signed off 2026-11-19 on #228.
- [x] **CH1.3d2 CONVENTIONS §3.15 revised.** Record: `docs/helmholtz/tasks/CH1.3d2.md`.
  `CONVENTION_VERSION` 4: "Signatures", "Level units", "Exponential M2L" exempt, Fₖ, §3.9's M2L
  range; proposal §3.1, §3.4, §4 amended. Signed off 2026-11-19 on #229 ("All signed off").
- [x] **H1.G Gate.** `/milestone-gate H1`. Report: `docs/helmholtz/gates/H1.md`; review
  `docs/reviews/2026-11-19-H1.md`. Closed 2026-11-19 (the user on #230, every question as
  recommended): the proposal amended (C1–C3, its Appendix A; C4 in the CH4.0e line); CH9.0
  (D5's dipoles) added; the review's triage applied (CHA.11, CHA.12 and the "Also" lines);
  CHA.11 on the priority list before CH4.1. *Needs:* every task above.

## H2: Yukawa on the host

**Gate:** CH2.5, Yukawa against FMM3D (or skipped with its reason).

- [x] **CH2.0 Laplace's real host drivers made generic over the kernel.** Report:
  `tetraleaf-validate/results/helmholtz-ch2-0.md`; the drivers `tetraleaf-exec/src/driver.rs`. Dense, Classes,
  Rotation, Batched (§6.2). *Done when:* constraint 1 of §2: no new dispatch below the level
  call, layouts and tables unchanged, Laplace bit for bit, its benchmarks within noise on the
  three machines; the PR states how p reaches Laplace's inner loops. Else dropped, and Yukawa
  gets copies. *Needs:* CH1.3a–c, CHR.2, D9. *Machines:* locust, Kathleen.
- CH2.4 Yukawa P2P in `tetraleaf-simd`. On a loop skeleton shared with Laplace if
  constraint 1 holds. *Done when:* an accuracy contract like simd-p2p.md's; cost per pair
  against Laplace; Laplace's P2P bit for bit and at least 90 % of the spike (C3S.4).
  *Needs:* CH1.2. *Machines:* locust, Kathleen. Split on 2026-10-31 into two sub-tasks,
  as Phase 3S split the inverse square root (T4) from the kernel (T5) (the user, "Split
  a/b"); CH2.4 is done when both are:
- [x] **CH2.4a The vectorised exponential.** Report: `tetraleaf-validate/results/helmholtz-ch2-4a.md`;
  design `docs/design/simd-p2p.md` §4.9. `tetraleaf_simd::exp` (`exp_slice`), the vector
  layer's `max` and `scale_pow2`, and `arch::exp` (Cody–Waite reduction and a polynomial), for
  x ≤ 0, f32 and f64, every ISA. *Done when:* within 4 u_T relative on [x_min, 0] (normal
  results), f32 exhaustively and f64 on 10⁷ samples against double-double, on every ISA a
  machine offers; below x_min the value at x_min; position and length invariance;
  no call in the inner loop; the existing entry points' inner loops unchanged instruction
  for instruction. *Needs:* CH1.2. *Machines:* Kathleen (AVX-512, the contract).
- [x] **CH2.4b Yukawa P2P on a shared skeleton.** Report: `tetraleaf-validate/results/helmholtz-ch2-4b.md`;
  design `docs/design/simd-p2p.md` §4.10. `P2pKernel::evaluate_yukawa`; the P2P
  loop skeleton generic over the pair arithmetic, Laplace's the first. *Done when:*
  CH2.4's: an accuracy contract like simd-p2p.md's (§4.9 onwards) on every ISA; cost per
  pair against Laplace on the three machines, beside green-kernels' modified Helmholtz;
  Laplace's P2P bit for bit, its benchmarks within noise (D9), the plane-wave M2L steps
  (`plane_wave_parts`) on locust included (CH2.4a's question 2, signed off 2026-10-31), and
  at least 90 % of the spike (C3S.4). *Needs:* CH2.4a. *Machines:* locust, Kathleen.
- [ ] **CH2.1 The real `HelmholtzOperator` path.** On CH2.0's drivers (or copies), Rotation
  and Dense M2L. *Done when:* the FMM against the direct sum on W1, W3 and W4 for κD from
  10⁻² to 10³; deterministic. Also (from CH4.9a2b): Yukawa's leaf P2M by a batched kernel,
  `p2m_yukawa` restated as `tetraleaf_tables::helmholtz::leaf::LeafP2m` restates `p2m` (`LeafP2m`
  refuses Yukawa; `Densities::p2m_kernel` keeps the per-point `p2m_yukawa` there), with its contract
  and the FMM's errors within 1.01× the per-point path's. *Needs:* CH2.0, D12. *Sign-off:* D4.
- [ ] **CH2.2 Per-level tables with per-level cache keys.** *Done when:* build cost and
  memory reported; levels reused across an octave sweep. *Needs:* CH2.1, D4.
- [ ] **CH2.3 Pruning and the calibration of p_l.** *Done when:* the relative error within
  target for κD from 10⁻² to 10³; on pruned levels the V, W and X lists empty and the buffers
  zeroed; the same levels pruned on both ranks of a two-rank run; θ_D set. Also (CH1.3c, the
  H1 gate's G3): under strong damping RCR's M2L loses relative accuracy (2.1·10⁻¹⁰ at κa = 30,
  |t| = 3w), round-off relative to its input; θ_D prunes such pairs. Also (review 2026-11-19,
  M13): the two f64 limits, RCR's M2M and L2L overflowing once p log₁₀(1/s) ≳ 300 and the
  (S|R) start column near p_src + p_tgt ≈ 175 at s < 1, refused by the builder (deep trees at
  κD = 10⁻²) and stated in §3.15's "Precision and range" at its next revision. Also (the H4 gate's second run, the user on #262, 2026-11-25: "Accept all (Recommended)"; review 2026-11-25, A12 and M18): the FMM checked
  with the constant density beside the random charges, on W2 beside W1, W3 and W4, at all six ε of D1, with the
  box measure beside ℓ2; and the constant density on a tree with more plateau levels than W1–W4 have
  (a cap or surface refined 100×, or 30–40 points per wavelength), where the plateau's levels may
  stack as W2's cap stacked on the high-frequency branch (CH4.8a). *Needs:* CH2.1,
  CH0.7.
- [ ] **CH2.5 Gate: Yukawa against FMM3D.** *Done when:* both codes against the same sampled
  direct sum over κD on W1 and W3, and the full-target crosscheck. CH1.4's R12 found FMM3D
  reliable for imaginary zk only in part. Signed off 2026-10-28 (CH1.4 report §5,
  question 2): the comparison with FMM3D is limited to where it met its eps, κD ≤ 30 at
  every eps and eps ≥ 10⁻⁶ (up to 2×) beyond; the sampled direct sum is the measure, and
  CH2.5 is not skipped. *Needs:* CH2.3, CH1.4. *Machines:* Kathleen.
- [ ] **H2.G Gate.** `/milestone-gate H2`. *Needs:* every task above.

## H3: Helmholtz on the host, RCR on every level

**Gate:** CH3.6, RCR on every level against FMM3D. The first Helmholtz deliverable, correct
at any kD.

Closed 2026-11-02 (`docs/helmholtz/gates/H3.md`). Its ticked lines were condensed at H6.G
(2026-11-12, the user: "Condense the roadmap"); each task's full line is in its record or
report and its progress entry (`docs/helmholtz/progress.md`), and in this file's history
(`git show ca938a5:docs/helmholtz/roadmap.md`).

- CH3.0 Builder: the leaf rule and complex densities. Split on 2026-10-31 into CH3.0a–c (the
  user, "Split a/b/c, proceed").
- [x] **CH3.0a tetraleaf-octree's maximum leaf width.** Record: `docs/helmholtz/tasks/CH3.0a.md`.
  `OctreeOptions::with_min_leaf_level`, default none; Laplace's outputs byte for byte (the
  user: "Bytes only, no D9").
- [x] **CH3.0b The Helmholtz builder's leaf rule.** Record: `docs/helmholtz/tasks/CH3.0b.md`.
  No leaf wider than c_leaf λ on W1–W4 over D1's range; W1–W4 in `tetraleaf-validate`.
- [x] **CH3.0c Complex densities.** Record: `docs/helmholtz/tasks/CH3.0c.md`. Complex charges
  in, split potentials and gradients out, one complex expansion per vector.
- [x] **CH3.5 The shared drivers extended to complex expansions.** Record:
  `docs/helmholtz/tasks/CH3.5.md`. Closed as a record on 2026-11-01 (the user, "Close as a
  record"): no driver change needed; the batched two-plane case moved to CH3.2.
- CH3.1 Complex planes and RCR on every level. Split on 2026-11-01 into CH3.1a–c (the user,
  "a/b/c").
- [x] **CH3.1a The Helmholtz operator.** Record: `docs/helmholtz/tasks/CH3.1a.md`.
  `HelmholtzOperator` on the shared drivers, real and complex k and charges, every body
  against `tetraleaf_ref::helmholtz`; Yukawa refused.
- [x] **CH3.1b The Helmholtz FMM.** Record: `docs/helmholtz/tasks/CH3.1b.md`. The point
  pipeline shared with `Fmm` (`src/pipeline.rs`); bit for bit across threads and ranks;
  Laplace bit for bit, minimal D9 on the M6 (the user).
- [x] **CH3.1c The acceptance runs.** Report: `tetraleaf-validate/results/helmholtz-ch3-1c.md`.
  `helmholtz_accuracy`; 61 cases within ε, deterministic.
- CH3.2 Dense, Classes and Batched where the node budget allows. Split on 2026-11-02 into
  CH3.2a–c; CH3.2d added at CH3.2c's sign-off.
- [x] **CH3.2a Dense.** Report: `tetraleaf-validate/results/helmholtz-ch3-2a.md`. Rotation to
  rounding, per pair faster only to p ≈ 6; opt-in.
- [x] **CH3.2b Classes.** Report: `tetraleaf-validate/results/helmholtz-ch3-2b.md`. Dense to
  rounding, tables about 18× smaller, faster than Rotation per pair only to p = 4; opt-in.
- [x] **CH3.2c Batched, and the measurements per p.** Report:
  `tetraleaf-validate/results/helmholtz-ch3-2c.md`. Its own block (the user: "Own Helmholtz block",
  "Split planes, 2n²"); fused, ahead of Rotation per pair to p ≈ 10 (M6, locust), 15 (Kathleen).
- [x] **CH3.2d The M2L rule per level.** Report: `tetraleaf-validate/results/helmholtz-ch3-2d.md`.
  `HelmholtzM2l::Auto` by degree and machine class, the builder's default (signed off
  2026-11-02); best on 34 of 36 rows, one miss at 1.106 (CH4.0).
- [x] **CH3.3 Calibration of p_l over (k a_l, ε).** Report:
  `tetraleaf-validate/results/helmholtz-ch3-3.md`. `DegreeRule::CALIBRATED` the default (the user:
  "Calibrate, new default", "Both branches", "ε/2, verify ≤ ε"); S7 an upper bound.
- CH3.4 Helmholtz P2P, real and complex k. Split on 2026-11-01 into CH3.4a–c.
- [x] **CH3.4a The vectorised sine and cosine.** Record: `docs/helmholtz/tasks/CH3.4a.md`.
  `tetraleaf_simd::sincos`, within 4 u_T for |θ| ≤ 2¹⁰ on every ISA.
- [x] **CH3.4b Helmholtz P2P on the skeleton.** Report:
  `tetraleaf-validate/results/helmholtz-ch3-4b.md`. `evaluate_helmholtz`, contract (32 + 12|k|r) u_T;
  Laplace instruction for instruction, D9 on Kathleen within noise.
- [x] **CH3.4c The Helmholtz operator's P2P.** Record: `docs/helmholtz/tasks/CH3.4c.md`.
  The operator's P2P on the vectorised kernel, the reference kept as the trusted path.
- [x] **CH3.6 Gate: RCR on every level against FMM3D.** Report:
  `tetraleaf-validate/results/helmholtz-ch3-6.md`. Ours within ε on all 57 case × ε points (largest
  0.63 ε), Kathleen jobs 617649–617651.
- [x] **H3.G Gate.** `/milestone-gate H3`. Report: `docs/helmholtz/gates/H3.md`. Closed
  2026-11-02 (the user, "i agree with the recommendations"): D6 signed off (H6 first); the
  proposal amended (C1–C6, its Appendix A); CH4.0b and CH4.0c added; the priority list for
  H6. *Needs:* every task above.

## H6: Exponential M2L on low-frequency levels (before H4 by default, D6)

**Gate:** the exponential M2L adopted where it wins, or kept opt-in (D8).

Its ticked lines were condensed at H6.G (2026-11-12, the user: "Condense the roadmap"); each
task's full line is in its record or report and its progress entry
(`docs/helmholtz/progress.md`), and in this file's history
(`git show 4b5b37f:docs/helmholtz/roadmap.md`).

- CH6.1 The §3.14 machinery generalised: six directions, per-level quadratures, merge and
  shift; errors within about 1.1× Rotation's. Split on 2026-11-03 into CH6.1a–e, and CH6.1f–h
  added (CH6.1e's question 3; CH6.1f split on 2026-11-09; CH6.1h by CH6.1g's sign-off 3).
- [x] **CH6.1a The quadratures.** Record: `docs/helmholtz/tasks/CH6.1a.md`. FMM3D's tables for
  iprec 1–3 with NIST's notice, our α counts by a Bessel bound to the corner; the kernel within
  1.50× and the chain within 1.18× the stated error on every domain.
- [x] **CH6.1b The generator.** Record: `docs/helmholtz/tasks/CH6.1b.md`.
  `tetraleaf_ref::helmholtz::plane_wave::generate`, the trusted path, bit for bit S4's.
- [x] **CH6.1c The pair chain.** Record: `docs/helmholtz/tasks/CH6.1c.md`. CONVENTIONS §3.15
  "Exponential M2L" (signed off 2026-11-03, as written); within 0.76× the stated error of RCR
  at resolving degrees.
- [x] **CH6.1d Merge and shift in the FMM.** Record: `docs/helmholtz/tasks/CH6.1d.md`.
  `HelmholtzM2l::PlaneWave` (opt-in); its pair path to 2.0e-16; deterministic, rank-invariant.
- [x] **CH6.1e The quadrature per level, measured.** Report:
  `tetraleaf-validate/results/helmholtz-ch6-1e.md`. iprec 2 at ε ≥ 10⁻⁶, within 1.070× Rotation's
  on 87 cases; signed off 2026-11-03 ("1–4 as recommended").
- [x] **CH6.1f The degree-matched quadrature on the pair chain.** Report:
  `tetraleaf-validate/results/helmholtz-ch6-1f.md`. The resolving degree recommended (candidate 4);
  signed off 2026-11-09 ("1–3 as recommended").
- [x] **CH6.1g The quadrature below 10⁻⁶ in the FMM.** Report:
  `tetraleaf-validate/results/helmholtz-ch6-1g.md`. iprec 3 at max(p_l, 28) and max(p_l, 32); errors
  0.002–0.80× Rotation's on 26 cases; signed off 2026-11-11 (D8 amended; CH6.1h added).
- [x] **CH6.1h The resolving degree on the M2L alone.** Report:
  `tetraleaf-validate/results/helmholtz-ch6-1h.md`. `PlaneWaveRaise::Rule` (`M2l` at 10⁻⁷, `Leaves`
  at 10⁻⁸); evaluation 0.43–0.83× and 0.31–0.71× Rotation's; D8 lifted at 10⁻⁷ (2026-11-11).
- CH6.2 Measured against Rotation, the table strategies and `hfmm3d`'s exponential M2L; D8
  proposed. Split from 2026-11-03 to 2026-11-08 into CH6.2a–h; CH6.2i added by CH6.2f's
  question 2.
- [x] **CH6.2a The waves buffer bounded.** Report: `tetraleaf-validate/results/helmholtz-ch6-2a.md`.
  Chunks of whole parents under a node budget, 8 GiB by default (signed off 2026-11-04);
  bit for bit at every bound.
- [x] **CH6.2b The M2L timed per level.** Report: `tetraleaf-validate/results/helmholtz-ch6-2b.md`.
  The M6 and locust: at 10⁻⁶ the exponential M2L the fastest on almost every level below the
  coarsest.
- [x] **CH6.2c CH6.2b's Kathleen runs.** Report: `tetraleaf-validate/results/helmholtz-ch6-2b.md`
  (§2, §4, §5). 0.24–0.61× Rotation's M2L at 10⁻⁶; level 2 slow at 1 × 40 (CH6.2f).
- [x] **CH6.2d Beside FMM3D's exponential M2L.** Report:
  `tetraleaf-validate/results/helmholtz-ch6-2d.md`. Per level FMM3D's over ours 0.44/0.89/1.11 at
  10⁻⁶ on Kathleen; W1 at kD = 500 22.4 s, our evaluation with the build apart
  (not reported at 10⁻⁶; 10.6–14.8 s on this case at 10⁻⁷ and 10⁻⁸, H6 §8), against 61.6 and
  43.2 s for FMM3D's whole call (review 2026-11-12, B1).
- [x] **CH6.2e The counts' corner.** Report: `tetraleaf-validate/results/helmholtz-ch6-2e.md`. The FMM
  held at reach 4.5; its adoption signed off 2026-11-07 and withdrawn by CH6.2g.
- [x] **CH6.2g The counts' reach adopted.** Report: `tetraleaf-validate/results/helmholtz-ch6-2g.md`.
  **Not adopted**: the worst pair at reach 4.5 errs up to 2·10⁴× the stated error; the corner
  kept (the user: "Keep the corner (Recommended)").
- [x] **CH6.2f The coarse levels' merge and shift.** Report:
  `tetraleaf-validate/results/helmholtz-ch6-2f.md`. `PlaneWaveSplit::Auto`, by child where parents
  leave threads idle; levels 2–3 0.38–0.88× Rotation's; signed off 2026-11-08.
- [x] **CH6.2h D8: the exponential M2L's level rule.** Report:
  `tetraleaf-validate/results/helmholtz-ch6-2h.md`. D8 signed off 2026-11-08 ("1–5 as
  recommended"): `Auto`'s default; the M2L 0.03–0.70× the previous `Auto`'s at 10⁻⁶.
- [x] **CH6.2i Laplace's plane-wave M2L on the coarse levels.** Report:
  `tetraleaf-validate/results/helmholtz-ch6-2i.md`. The split by child in Laplace's level call, bit
  for bit, D9 on three machines; signed off 2026-11-09 ("ok, 1–3 as recommended").
- [x] **H6.G Gate.** `/milestone-gate H6`. Report: `docs/helmholtz/gates/H6.md`. Closed
  2026-11-12 (the user: "Run first, then close (Recommended)", then all as recommended): the
  large cases at 10⁻⁷ and 10⁻⁸ verified (Kathleen job 633000: errors 0.45–0.84× Rotation's,
  the evaluation 0.27–0.65×); the proposal amended (C1–C8, its Appendix A); the priority list
  for H4; the roadmap condensed. *Needs:* every task above.

## H4: High-frequency levels

**Gate:** CH4.4 and CH4.5, the kD sweep and host performance against FMM3D.

The done lines are condensed at the H4 gate (2026-11-24); in full in `git show 4d28d41:docs/helmholtz/roadmap.md`;
CH4.8a, CH4.8b and H4.G condensed at its second run (2026-11-25), in full in `git show 03bc2bd:docs/helmholtz/roadmap.md`.

- [x] **CH4.0 The M2L rule's sparse-level term.** `Auto` batches a level by m̂ from the global
  tree; dense budget 8 GiB per node; CH3.2d's sweeps within 1.049 of the best. Report `tetraleaf-validate/results/helmholtz-ch4-0.md`.
  Signed off on #211 (2026-11-12): CH8.0 added, the residual CH4.0f's; §6.4 amended.
- [x] **CH4.0b The worst-box concentration.** The cause is placement (box corners, W1's tangent
  points); smooth densities within ε in ℓ2. Report `tetraleaf-validate/results/helmholtz-ch4-0b.md`. Signed off on #212: ε stays ℓ2,
  the box measure gated at 100 ε (§9.1 amended).
- [x] **CH4.0c The degree rule beyond its measured range.** W1 at kD = 1000 and the
  under-resolved W3 within ε, no correction. Report `tetraleaf-validate/results/helmholtz-ch4-0c.md`. Signed off on #214 (§3.5 amended
  in #215).
- CH4.0d The Helmholtz leaf size. Added on 2026-11-12 (the user: "Can you add leaf size tuning to
  the roadmap?"); split on 2026-11-13 into CH4.0d1–d3 (the user: "d1/d2/d3 (Recommended)").
- [x] **CH4.0d1 Empty boxes left out of the lists.** `PlanRules::empty_boxes`, on by default;
  evaluation 0.20–0.56 on sparse trees. Report `tetraleaf-validate/results/helmholtz-ch4-0d1.md`. Signed off on #216 (2026-11-13).
- [x] **CH4.0d1b The empty-box rule's build cost.** The plan built once on the octree's index
  (`Plan::with_rules_on`); build 0.80–1.03× the rule off's. Report `tetraleaf-validate/results/helmholtz-ch4-0d1b.md`.
- [x] **CH4.0d2 V pairs between leaves of few points as near pairs.** `VNear` per (parent,
  source), opt-in, its default with CH4.0d3. Report `tetraleaf-validate/results/helmholtz-ch4-0d2.md`. Signed off on #218.
- [x] **CH4.0d3 The leaf size and c_leaf.** W/X pairs as near pairs; leaf size, V/W/X rules and
  c_leaf 2 the defaults (0.41–0.83× build plus evaluation). Report `tetraleaf-validate/results/helmholtz-ch4-0d3.md`. Signed off on
  #219 (2026-11-14); D11 amended.
- [x] **CH4.0e The comparison under the default, and H4's expected gain.** CH3.6's matrix beside
  `hfmm3d`; H4 inside D1 bounded at 14–31 % of W1 at kD = 500. Report `tetraleaf-validate/results/helmholtz-ch4-0e.md`. Signed off on
  #220 (2026-11-16): H4 after CH4.0h and CH4.0i; block budget 512 MiB.
- [ ] **CH4.0f Batched's split across the matrix's rows.** Added on 2026-11-12 at CH4.0
  (report §4, §8): with CH4.0's split by blocks, Batched on levels of 64–5,000 boxes at 40–72
  threads is still 1.1–5× Rotation's level time (0.2–2 ms per level, up to 1.28× of the
  per-level best on ms-scale evaluations), because narrower blocks read each table for fewer
  columns. Split a block's column-kernel calls across the threads by rows of the n × n table
  (each thread a slice of the local's coefficients for every column), so that a level of few
  blocks keeps its columns per read; every element keeps its chain. *Done when:* the level call
  the pair chain bit for bit at every split and thread count; on CH6.2h's matrix at 1 × 40 and
  1 × 72 the levels of 64–5,000 boxes within 1.10 of Rotation's level time or Rotation's, the
  evaluation within 1.10 of the per-level best. It needs a row-slice entry point of
  `tetraleaf_simd::MatrixKernel` (Laplace's crate: constraint 1, bit for bit and D9 on the three
  machines). *Needs:* CH4.0. *Machines:* Kathleen, locust (the M6 with the user's go-ahead).
- [ ] **CH4.0g Empty boxes in Laplace's trees.** Added on 2026-11-13 at CH4.0c (the user, on
  #214: "Yes. Please."). Laplace's leaves come from the points alone, so an empty leaf is the
  empty sibling of a split box or a split of the 2:1 balance: probably rare on the uniform
  cube, probably common on the Plummer sphere, the Gaussian clusters and points on a surface
  (not counted). First count, without timing (the M6 needs no quiet gate): per tree and level
  the leaves and boxes without source or target points and the V, W and M2M pairs CH4.0d1's
  rule leaves out, on the cube, the Plummer sphere, the clusters and a sphere's surface at
  Phase 6 T11's default leaf sizes. Then, where the pairs drop noticeably, the rule under
  Laplace's `FmmBuilder` timed against the rule off. The V conversion is not for Laplace: its
  leaf size is tuned against P2P and the near rule already converts W and X. *Done when:*
  the counts reported; and either the rule proposed as Laplace's default with its gain per case
  (the output bit for bit, −0 and +0 equal, on host, device and P ranks; D9's benchmarks within
  noise or better on the M6, locust and Kathleen; the work cut's counts, the device tuner's
  pair buckets and Batched's columns per read checked), or the report shows it does not pay.
  *Needs:* CH4.0d1 (the rule). *Machines:* Kathleen, locust (the M6 with the user's go-ahead
  for timings). *Sign-off:* Laplace's default.
- [x] **CH4.0h The Helmholtz build by part, and parallel.** Report `tetraleaf-validate/results/helmholtz-ch4-0h.md`. Signed off on
  #221 (2026-11-16): CH4.0j on the list, CH4.0k off.
- [x] **CH4.0i The coarse RCR levels.** The Rotation M2L by pairs (`RcrSplit::Auto` by degree,
  the default). Report `tetraleaf-validate/results/helmholtz-ch4-0i.md`. Signed off on #222 (2026-11-16): CH4.0l off the list.
- [ ] **CH4.0l The RCR pair: a class's pairs at once, and vectorised chains.** Added on
  2026-11-16 at CH4.0i (report §3, §4): by pairs, level 2 of W1 at kD = 500 runs at 0.40–0.54
  of its own one core, each pair at p = 131 reading a 12 MB coaxial table and 12 MB of y-blocks
  from memory on every thread; and the pair itself is 1.8–2.4× S6's prototype on the same
  core, in latency-bound sum chains. Apply each coaxial row and each angle's blocks to all
  pairs of a class at once (each pair's chain unchanged, so bit for bit), then vectorise the
  chains (a new order of operations: `rcr` and `HelmholtzTables` in step, the Helmholtz
  outputs change in their last bits). *Done when:* level 2 at ≥ 70 % of linear against its
  one core on the node, or the report says why not; the pair against S6's on the same core;
  the FMM's errors within 1.1× of before and every bit identity between the RCR paths kept.
  Helmholtz only. *Needs:* CH4.0i. *Machines:* Kathleen, locust. Off the priority list
  (CH4.0i's question 3, the user, 2026-11-16, on #222): H4's diagonal form takes these levels'
  passing classes.
- [x] **CH4.0j The exponential M2L's tables within a level in parts.** Bit for bit; build part
  0.05–0.29 of CH4.0h's. Report `tetraleaf-validate/results/helmholtz-ch4-0j.md`. Signed off on #223 (2026-11-16): CH4.0k stays off.
- [ ] **CH4.0k The shared build steps on the node.** Added on 2026-11-16 at CH4.0h (report §3,
  §6): the plan with its rules (`Plan::with_rules_on`), the routes (`Points::route`) and the
  octree run serially in each rank: 2.0–2.4 s of W1 at kD = 500's build at one rank per node,
  2.8–3.8 s of W3 at kD = 100's (4·10⁶ points); several ranks split them by points. Time them
  by part, then thread what pays. Laplace's paths (`tetraleaf-plan`, `pipeline.rs`, `tetraleaf-octree`):
  constraint 1, its outputs bit for bit, D9's benchmarks on the M6, locust and Kathleen.
  Since CH4.0j also most of the build of CH3.6's cases under 1 s at one rank per node (the plan
  5–66 ms, the routes 13–24 ms, the octree 5–19 ms; CH4.0j's report §4).
  *Done when:* W1 at kD = 500's build at 1 × 40 and 1 × 72 at most half its one-core time, the
  output bit for bit; Laplace's builds no slower (D9). *Needs:* CH4.0h. *Machines:* Kathleen,
  locust (the M6 with the user's go-ahead). *Sign-off:* the Laplace change.
- CH4.1 Split on 2026-11-19 into CH4.1a (measure first, the H1 gate's question 5) and CH4.1b.
- [x] **CH4.1a The product conversions, measured.** An FFT about 1–2 % of a run. Report
  `tetraleaf-validate/results/helmholtz-ch4-1a.md`. Signed off on #234 (2026-11-19): H4 goes on with the product, the FFT stays item 11.
- [x] **CH4.1b The signature buffers and the level's scratch.** `diagonal::level`: buffers, the
  level driver, `scratch_shape` under the node budget. Report `tetraleaf-validate/results/helmholtz-ch4-1b.md`. Signed off on #236.
- CH4.2 Diagonal M2L per (level, class). Split on 2026-11-19 into CH4.2a and CH4.2b.
- [x] **CH4.2a The transfer functions per (level, class) and the pair.** `DiagonalM2l` in
  tetraleaf-tables. Record `docs/helmholtz/tasks/CH4.2a.md`. Signed off on #237 (2026-11-20).
- CH4.2b The diagonal M2L in the operator. Split on 2026-11-20 into CH4.2b1 and CH4.2b2.
- [x] **CH4.2b1 The diagonal M2L in the operator: the level call and the rule.**
  `HelmholtzDiagonal`, opt-in; the whole-level fallback. Record `docs/helmholtz/tasks/CH4.2b1.md`.
  Signed off on #238.
- [x] **CH4.2b2 The diagonal M2L on W1, measured.** Errors 1.000× Rotation's. Report
  `tetraleaf-validate/results/helmholtz-ch4-2b2.md`. Signed off on #239 (2026-11-20): CH4.2c before CH4.3, CH4.2d off the list.
- CH4.2c The transfer functions interpolated, and built in parts. Added at CH4.2b2; split on
  2026-11-20 into CH4.2c1 and CH4.2c2.
- [x] **CH4.2c1 The transfer functions interpolated, and built by ring.** Record
  `docs/helmholtz/tasks/CH4.2c1.md`. Signed off on #240.
- [x] **CH4.2c2 The interpolated build in the operator, and W1's tables timed.** Tables 4–6×
  faster. Report `tetraleaf-validate/results/helmholtz-ch4-2c2.md`. Signed off on #241: `Interpolated` the default build.
- [ ] **CH4.2d Dispatched kernels for the diagonal pair and the conversions.** (Added at
  CH4.2b2, its question 3, signed off 2026-11-20; off the list, kept off at CH4.5, its question 2, signed off 2026-11-24.) *Done when:*
  the pair and the conversions' column kernel on `tetraleaf-simd`'s ISAs, bit for bit or within
  rounding of the scalar path as the crate's rules say (constraint 1); the diagonal levels'
  call timed against CH4.2b2's (pairs 50–80 %, conversions 24–49 % of its thread time).
  *Needs:* CH4.2b2. *Machines:* Kathleen, locust.
- CH4.3 The per-(level, class) rule, calibrated. Split on 2026-11-20 into CH4.3a–CH4.3e.
- [x] **CH4.3a The bound's constants in the full FMM.** C_T (0.15 direct, 0.23 interpolated) and
  θ = 1/2 hold, margin 16 at 10⁻⁸. Report `tetraleaf-validate/results/helmholtz-ch4-3a.md`. Signed off on #242 (2026-11-20).
- CH4.3b The time rule, and the diagonal form by default. Split on 2026-11-20 into CH4.3b1 and
  CH4.3b2 (CH4.3b2a–d).
- [x] **CH4.3b1 The time rule per level, built and calibrated.** `HelmholtzDiagonal::Auto`
  (opt-in), `DiagonalCosts` per machine class; none below k a 2.8. Report `tetraleaf-validate/results/helmholtz-ch4-3b1.md`. Signed
  off on #243.
- CH4.3b2 The diagonal form by default. Split on 2026-11-20 into CH4.3b2a–d.
- [x] **CH4.3b2a The time rule at matched accuracy.** A level with waves keeps p_w and the waves
  unless the form beats them. Report `tetraleaf-validate/results/helmholtz-ch4-3b2a.md`. Signed off on #244 (2026-11-21).
- [x] **CH4.3b2b The V rule on the form's levels, and the defaults' errors at 10⁻⁶.**
  `VAfterPlan` opt-in. Report `tetraleaf-validate/results/helmholtz-ch4-3b2b.md`. Signed off on #245: CH4.3b2d added.
- [x] **CH4.3b2d The waves' V crossover on raised levels, in one plan.** `VRaisedPrice::PerBox`
  the default. Report `tetraleaf-validate/results/helmholtz-ch4-3b2d.md`. Signed off on #247 (2026-11-21).
- [x] **CH4.3b2c The diagonal form by default.** `Auto` the builder's default; errors
  0.991–1.034×, evaluation 0.37–0.99×. Report `tetraleaf-validate/results/helmholtz-ch4-3b2c.md`. Signed off on #251 (2026-11-22): D8
  amended (§5.3, §6.4), CH4.3f added off the list, CH8.0's Also.
- CH4.3c The high-frequency degrees refit, and C₀ at large degrees. Split on 2026-11-21 into
  CH4.3c1 and CH4.3c2 (CH4.3c3 added at CH4.3c2).
- [x] **CH4.3c1 C₀ at large degrees, against a double-double reference.** C₀ = 200 confirmed to
  p = 251. Report `tetraleaf-validate/results/helmholtz-ch4-3c1.md`. Signed off on #246 (2026-11-21).
- [x] **CH4.3c2 The high-frequency degrees refit.** 4–10 degrees lower at 10⁻³. Report
  `tetraleaf-validate/results/helmholtz-ch4-3c2.md`. Signed off on #249: the refit the default, §3.5 amended, CH4.3c3 added.
- [x] **CH4.3c3 The constant density in the calibration.** Refit under a cap of 0.82, every case
  within 0.65 ε. Report `tetraleaf-validate/results/helmholtz-ch4-3c3.md`. Signed off on #254 (2026-11-22): the default, §3.5 amended.
- CH4.3d The transfer function's length. Split on 2026-11-22 into CH4.3d1 and CH4.3d2.
- [x] **CH4.3d1 The transfer function's length, per pair.** No length rule reaches the 0–16 %.
  Report `tetraleaf-validate/results/helmholtz-ch4-3d1.md`. Signed off on #252: 2P − 2 kept, §3.6 amended, review M8 closed.
- [x] **CH4.3d2 The shorter transfer function in the FMM.** Closed unrun by CH4.3d1's sign-off
  (#252, 2026-11-22).
- [x] **CH4.3e c_leaf re-decided with the rule.** c_leaf 4 the default (0.874× and 0.790× c_leaf
  2's, from the under-resolved W3). Report `tetraleaf-validate/results/helmholtz-ch4-3e.md`. Signed off on #253: D11, §3.7, §5.1, §5.7
  amended.
- [ ] **CH4.3f The time rule's costs at many ranks per node.** (Added at CH4.3b2c, its question
  4, signed off 2026-11-22; off the list, kept off at CH4.5, its question 3, signed off 2026-11-24.) *Done when:* the waves' cost per box at one rank per
  core on W3's dense levels explained (Kathleen 40 × 1, W3 at kD = 500, 10⁻⁸, level 4: 20.5 s
  against 2.5 s at 1 × 40, the form 2.6 s; CH4.3b2c's report §8) and the time rule's costs fitted per
  ranks per node (or read from the rank's share), the dense-level misses at one rank per node (3.5 %
  and 2.2 % of the levels' time, its §5) refitted or closed, verified as CH4.3b2a on Kathleen and
  locust at 1 × N and N × 1. Also (the H4 gate's triage, the user on #258, 2026-11-24: "Accept all (Recommended)"; review
  2026-11-24, B12): first Kathleen's drift within a launch diagnosed (CH4.3b2c §6: the same
  setting 1.09–1.20× slower the second time; an `off, off` launch, `numastat`, 2 × 20 against
  1 × 40), before constants are fitted again. *Needs:* CH4.3b2c. *Machines:* Kathleen, locust.
- CH4.4 kD sweep, with the FMM3D gate. Split on 2026-11-22 into CH4.4a and CH4.4b.
- [x] **CH4.4a The kD sweep beside FMM3D.** Within ε on W1 to kD = 1000, W2 to 500, W3 to 200
  (at most 0.76 ε); one box at 114.7 ε beyond D1. Report `tetraleaf-validate/results/helmholtz-ch4-4a.md`. Signed off on #255.
- [x] **CH4.4b The crosscheck mapped to levels and classes.** No error at a class boundary; the
  100 ε bound gated within D1. Report `tetraleaf-validate/results/helmholtz-ch4-4b.md`. Signed off on #256: §9.1 amended.
- [x] **CH4.5 Host performance against FMM3D at matched accuracy.** 0.13–0.18× FMM3D's call for
  N ≥ 10⁶ at the best split, at matched accuracy under CH4.3c3's rule (the default 1.01–1.02× its
  cost, CH4.8b §3; R10 of the H4 gate's second run), its faster variant (`default` alone on the plane wave at kD = 500 and
  the eight constant-density cases, CH4.4a's runs; review 2026-11-24, B11). Report `tetraleaf-validate/results/helmholtz-ch4-5.md`.
  Signed off on #257 (2026-11-24): CH4.2d, CH4.3f off the list, CH4.7 added.
- [ ] **CH4.7 The form choice prices the tables' build for small trees.** (Added at CH4.5, its
  question 4, signed off 2026-11-24; off the list.) *Done when:* the time rule (or the leaf rule) weighs the build of
  Batched's matrices and the diagonal tables against the evaluation where one solve is the use, so
  that the six Kathleen rows they explain (CH4.5 §3.4: W1 and W2 at kD = 30 and W3 at kD = 10, both
  densities, 10⁻³; Batched at p = 11, 15–101 ms of the build) are not slower than FMM3D's call,
  verified on CH4.5's small cases on Kathleen and locust with the large cases unchanged. Also (the
  H4 gate's triage, the user on #258, 2026-11-24: "Accept all (Recommended)"; review 2026-11-24, A9, A10): a density-aware leaf rule
  weighed beside c_leaf (the width cap only on boxes above a point count, the P2P phase contract
  the hard cap), since c_leaf 4's gain is the under-resolved W3's alone (CH4.3e: 0.995 and 1.000
  without it); the time rule's choices scored per level on cases outside its fit (W2; the M6
  with the user's go-ahead). Also (the H4 gate's second run, the user on #262, 2026-11-25: "Accept all (Recommended)"; review 2026-11-25, A11): 128 against 256 points per
  leaf on NEON after H4 (CH4.0d3 chose 128 at ε ≥ 10⁻⁷ for the evaluation alone; at 1 × 72 W1 at
  kD = 1000 builds in 13.4–14.0 s at 128 against 7.6 s at 256, CH4.5 §5): W1 and W2 at kD = 250–1000
  on locust at 1 × 72 and 4 × 18, the build and the evaluation apart. Also (from CH4.6a2, its report
  §5–§6): the levels the time rule keeps on the waves where the form runs faster under either φ
  step, Kathleen's W3 at kD = 500, 10⁻⁸, level 4 (0.87× with the product) and locust's W1 at
  kD = 200, 10⁻⁷, level 4 (0.61×), W1 at kD = 500, 10⁻⁸, level 5 (0.91×) and W3 at kD = 200, 10⁻⁸,
  level 3 (0.88×), in that scoring. *Needs:* CH4.5. *Machines:* Kathleen, locust.
- CH4.6a The optional rustfft back end of the φ step. Added at CH0.9 (D7; the user on #264,
  2026-11-26: "Q1-3 as recommended"; proposal §6.7, §10). Split on 2026-11-26 into CH4.6a1–a3,
  measuring before the second grid is built (the user: "Split a1/a2/a3 (Recommended)"); rustfft
  pinned `=6.4.1` (the user, 2026-11-26: "Exact =6.4.1 (Recommended)").
- [x] **CH4.6a1 The trait and rustfft on the product's grid.** `HelmholtzPhiStep::RustFft` behind
  `tetraleaf-exec`'s feature `rustfft`: ≤ 5.1 u from the product per conversion (p ≤ 245); in the FMM
  the form's own 1 u sensitivity, errors 1.0000×; fingerprints agreed over the ranks; CI
  `run-tests-rustfft`. Record `docs/helmholtz/tasks/CH4.6a1.md`. That reading of §10's
  "round-off" signed off 2026-11-26 (the user on #265: "I accept PR Question 1").
- [x] **CH4.6a2 The time per level against the product, both grids priced.** rustfft's conversions
  0.23–0.71× the product's, the diagonal levels 0.60–0.95×, the evaluation 0.89–0.95× (Kathleen,
  1 × 40) and 0.91–0.99× (locust, 1 × 72), errors 0.993–1.000×; the smooth grid 1.000–1.005× per
  evaluation; the time rule's costs need no FFT variant; the harness's `--phi-step`. Report
  `tetraleaf-validate/results/helmholtz-ch4-6a2.md`, record `docs/helmholtz/tasks/CH4.6a2.md`. The M6
  skipped (the user, 2026-11-26: "Skip the M6"). Signed off on #266 (2026-11-26, the user: "Q1-3 as
  recommended"): the product's grid kept and CH4.6a3 closed unrun; no FFT costs; the acceptance
  without the M6; proposal §6.7, §10 and D7 amended.
- [x] **CH4.6a3 The smooth grid, chosen per level by its time.** Closed unrun by CH4.6a2's sign-off
  (#266, 2026-11-26: "Q1-3 as recommended"): the smooth grid priced at 1.000–1.005× of an
  evaluation on Kathleen and locust; rustfft stays on the product's grid.
- [ ] **CH4.6b The optional FFTW back end of the φ step.** CH4.6 until CH0.9; off the priority list
  since then (D7). *Done when:* against the product to round-off; per-box plans with bits the same
  across batch positions, sizes and threads, and fingerprints all-gathered; off by default, a default
  build linking no FFTW, cargo-deny, a CI job with the feature; the GPL warning in the crate
  documentation, the README, the builder option and `Cargo.toml`. *Needs:* CH4.1b, CHR.3, CH0.9,
  D7, CH4.6a1 (the trait).
- [x] **CH4.8a The constant density on W2: the error mapped and a remedy chosen.** The degree rule's
  truncation on W2's deep tree against a face; `h1` (one degree after the cap) chosen. Report
  `tetraleaf-validate/results/helmholtz-ch4-8a.md`. Added at the H4 gate (M15). Signed off on #260 (2026-11-24).
- [x] **CH4.8b The constant density on W2: the remedy as the default.** `DegreeRule::margin` the default
  (§3.5 amended); W1–W3's constant density within 0.67 ε and 76.3 ε at all six ε, 1.01–1.02× the time.
  Report `tetraleaf-validate/results/helmholtz-ch4-8b.md`. Signed off on #261 (2026-11-25).
- CH4.9 The Helmholtz leaf expansions by SIMD kernels. Added on 2026-11-26 after CH4.6a2 (the user:
  "Yes. Add P2M and L2P to the priority list as tasks"; item 14 of the priority list; no proposal §10
  entry, so the acceptance is the *Done when* below with constraint 1). Measured on Kathleen (job
  673083, 1 × 40, the defaults, `helmholtz_accuracy --levels` with P2M, M2M, L2L and L2P per level,
  outputs `~/projects/fmm-archive/kathleen/runs/explore-kinds/`): P2M + L2P 13–31 % of an
  evaluation on W1 at kD = 250–1000 and W3 at kD = 500 (largest at 10⁻⁸ and at kD = 250; L2P the
  larger, 1.58 s of W1 at kD = 1000's 15.2 s on the finest level), where they run one scalar
  `tetraleaf_ref::helmholtz::leaf` call per point (`Densities::{p2m, l2p}`).
- CH4.9a P2M by a SIMD kernel. Split on 2026-11-26 into CH4.9a1 and CH4.9a2, measuring before the
  kernel is built (the user: "Split a1/a2 (Recommended)"): on Kathleen the reference took about
  10 ns per term for P2M and 24 ns for L2P with gradients (W1 at kD = 1000, p = 17), so most of a
  point's time is setup (the radial recursion, the harmonics, the exponential), and the kernel's
  shape is not yet known. The original *Done when*, with the table first: `helmholtz_accuracy
  --levels` prints P2M, M2M, L2L and L2P per level; a P2M kernel vectorised over a leaf's points on
  the P2P's ISAs, for real and complex k, real and complex charges and m vectors, run by the
  operator's `P2pChoice` (`Reference` the per-point reference calls), within a contract per leaf
  against `tetraleaf_ref::helmholtz::leaf::p2m` in f64 and f32, a leaf's bits the same at every
  thread count and on P ranks the one-rank FMM's, the FMM's errors within 1.01× the reference
  path's, P2M's time per level on Kathleen and locust (the M6 with the user's go-ahead); constraint
  1; Yukawa's `p2m_yukawa` an Also for H2.
- [x] **CH4.9a1 Where a point's P2M and L2P time goes, and the per-level table.** The harmonics take
  44–81 % of a leaf point's P2M (f64, p ≥ 12, medians over k̂) and 58–84 % of L2P with gradients,
  the sum 4–10 %, on Kathleen and locust; `helmholtz_accuracy --levels`' per-level kinds and the `helmholtz_leaf`
  example. Report `tetraleaf-validate/results/helmholtz-ch4-9a1.md`, record
  `docs/helmholtz/tasks/CH4.9a1.md`. Questions 1–2 signed off 2026-11-27 (the user on #267: "Q1-2
  as recommended"): CH4.9a2's kernel restructured safe Rust first, the leaf harmonics built with
  P2M in CH4.9a2 and reused by CH4.9b.
- CH4.9a2 P2M by the kernel chosen. Split on 2026-11-27 into CH4.9a2a and CH4.9a2b (the skill's size
  rule: a new routine in `tetraleaf-math`, a kernel in `tetraleaf-tables`, its dispatch in
  `tetraleaf-exec` and runs on two machines, more than one reviewable PR). The *Done when* before the
  split: CH4.9a's *Done when* less the table, the kernel as signed off at CH4.9a1 (its questions 1–2;
  the user on #267, 2026-11-27: "Q1-2 as recommended"): a restructured kernel in safe Rust, not an
  ISA kernel first, with (1) the leaf harmonics from coefficients tabulated once per degree, a
  routine shared with L2P beside `tetraleaf_math::orthonormal::harmonics` (or in
  `tetraleaf-tables::helmholtz`), tested against it per point within a contract over the leaf
  frame's directions, the poles included; (2) the radial functions in real arithmetic for real k
  (complex k keeps the reference's path); (3) a leaf's points batched, the sum written to vectorise
  over points without `unsafe`; no `tetraleaf-simd` change (constraint 1), and an ISA kernel only if
  the measured time stays far from the sum's (then a question in the PR). Also (from CH4.9a1): in
  f32 at p = 41–45 the reference ran 1.37–1.40× f64's time on Kathleen as the median over k̂, up to
  1.75× at k̂ = 2–4 (likely subnormal values on x86, untested); the f32 kernel avoids or reports it.
- [x] **CH4.9a2a The leaf harmonics and the real radial functions, batched over points.**
  `LeafHarmonics` within 1.02 (n + 1)² u √((2n + 1)/(4π)) per term of `harmonics`, `bessel_j_real` its
  bits in f64 (f32 flushes under-range values); on one core the harmonics 0.06–0.09× and the radial functions 0.13–0.24× the reference's, P2M
  projected at 0.14–0.23× (f64), Kathleen's f32 slowdown gone. Report
  `tetraleaf-validate/results/helmholtz-ch4-9a2a.md`, record `docs/helmholtz/tasks/CH4.9a2a.md`.
  Question 1 signed off 2026-11-27 (the user on #268: "Q1 as recommended"): CH4.9a2b's kernel on
  these routines, the FMM's errors in CH4.9a2b confirming them.
- [x] **CH4.9a2b P2M by the kernel.** `tetraleaf_tables::helmholtz::leaf::LeafP2m` (16 points per
  batch on CH4.9a2a's routines, safe Rust) under every ISA, the reference under `Reference` or
  `reference_leaves`; within C = 10 of its contract per coefficient (about 1 measured, f64 and f32), the
  FMM's errors 0.993–1.000× the reference path's (W1's and W2's constant density included); P2M per run
  0.25–0.31× on Kathleen and 0.15–0.27× on locust, now 0.6–2.4 % of an evaluation. Report
  `tetraleaf-validate/results/helmholtz-ch4-9a2b.md`, record `docs/helmholtz/tasks/CH4.9a2b.md`.
  Question 1 signed off 2026-11-27 (the user on #269: "Q1 as recommended"): no ISA kernel for P2M
  now, the ISA question carried to CH4.9b.
- CH4.9b L2P with gradients by a SIMD kernel. Split on 2026-11-27 into CH4.9b1 and CH4.9b2 (the
  skill's size rule, as CH4.9a2: routines in `tetraleaf-math`, a kernel in `tetraleaf-tables`, its
  dispatch in `tetraleaf-exec` and runs on two machines). The *Done when* before the split: as CH4.9a
  for L2P, potentials and gradients, against `tetraleaf_ref::helmholtz::leaf::l2p`, reusing CH4.9a2's
  leaf harmonics and radial functions and adding their derivatives (CH4.9a1's question 2); the
  evaluation's time under the defaults before and after CH4.9a–b on CH4.6a2's cases. Also (from
  CH4.9a2b, its question 1, signed off 2026-11-27, the user on #269: "Q1 as recommended"): the ISA
  question carried here. The workspace builds for each architecture's baseline (SSE2, NEON: two f64
  lanes), on which P2M's batched sum, with a lane reduction per term, was no faster than the scalar
  sum on Kathleen (0.96–1.23×, report §5); L2P sums over terms per target, with no reduction across
  lanes, so measure its sum against its floor on both machines and, if it stays far from it on x86,
  decide one ISA kernel for both sums. Also (from CH4.9a2b): the f32 kernel's sum only 0.88× its f64
  sum on Kathleen at p = 45 (subnormal products possible, not diagnosed); reuse `reference_leaves`
  for L2P's reference path.
- [ ] **CH4.9b1 The leaf harmonics' derivatives and the real radial derivatives, batched over
  points.** *Done when:* `LeafHarmonics` gives per order the column Q̄ₙᵐ (`columns`' bits), ∂θ Q̄ₙᵐ
  and Q̄ₙᵐ / sin θ of a batch as lanes, without a division by sin θ, held per term to a contract
  against `harmonics_derivatives` over the leaf frame's directions (the poles, the centre and the
  axes included) in f64 to p = 100 and f32 to p = 64, no subnormal output, lanes and positions bit
  for bit; `bessel_j_real`'s derivatives (`bessel_j`'s derivatives' real parts bit for bit in f64,
  in f32 but for flushed values); `helmholtz_leaf --lanes` times them per point beside the
  reference's L2P parts with the reference's sum on their values (the projected L2P with gradients)
  on one core of Kathleen and locust; no `tetraleaf-simd` change (constraint 1). *Needs:* CH4.9a2b.
  *Machines:* Kathleen, locust.
- [ ] **CH4.9b2 L2P with gradients by the kernel.** *Done when:* CH4.9a's *Done when* for L2P,
  potentials and gradients, on CH4.9b1's routines: the kernel in `tetraleaf-tables::helmholtz::leaf`
  against `tetraleaf_ref::helmholtz::leaf::l2p` within a contract per target in f64 and f32 (real
  and complex k, the centre's tangential limit), run by the operator's `P2pChoice` with
  `reference_leaves` (`Reference` the per-point calls), m vectors, a leaf's bits the same at every
  thread count and on P ranks the one-rank FMM's, the FMM's errors within 1.01× the reference
  path's (W1's and W2's constant density included), L2P per level and the evaluation under the
  defaults before and after CH4.9a–b on CH4.6a2's cases on Kathleen and locust; the ISA question
  (the sum against its floor on both machines) and the f32 sum as above; constraint 1. *Needs:*
  CH4.9b1. *Machines:* Kathleen, locust.
- [ ] **CH4.9c Laplace's P2M and L2P by SIMD kernels.** Added on 2026-11-26 (the user: "Please add it
  to the roadmap"); off the priority list. A Laplace task: constraint 1 binds in full. Measured on
  Kathleen (job 673084, the cube at N = 10⁷, f64, the library's defaults, gradients, `scaling --kinds`
  at 1 × 40 and 40 × 1; outputs `~/projects/fmm-archive/kathleen/runs/laplace-kinds/`): P2M + L2P
  19 % of the compute at p = 3, 34–35 % at p = 8, 39–42 % at p = 12 and 36–40 % at p = 18, 1.4–7.1×
  the P2P from p = 8 on, L2P 1.5–2.3× P2M; they run `tetraleaf_ref::leaf::{p2m_many, l2p_many}`
  per point (`LaplaceOperator::{p2m, l2p}`). *Done when:* P2M (charges and dipoles) and L2P (with
  gradients) kernels in `tetraleaf-simd` on the P2P's ISAs, vectorised over a leaf's points and
  sharing CH4.9a–b's design where it fits, run under `P2pChoice` (`Reference` the per-point calls,
  the trusted path); within a contract per leaf against `tetraleaf_ref::leaf` in f64 and f32; the
  FMM's errors within the C3.2 and C3.3 gates; every bit-identity check of `mpi_exec`,
  `invariance` and `multi_rank` (threads, overlap, ranks, several vectors) holding for the new
  kernel; Laplace's outputs on `Reference` unchanged bit for bit; on the default kernel new bits
  only with the user's sign-off (constraint 1), and the Laplace fixtures kept; the device path
  unchanged; D9's benchmarks before and after on the M6 (with the user's go-ahead), locust and
  Kathleen. *Needs:* CH4.9b (the design shared). *Machines:* Kathleen, locust, the M6.
- [x] **H4.G Gate.** `/milestone-gate H4`. Report `docs/helmholtz/gates/H4.md`, reviews
  `docs/reviews/2026-11-24-H4.md` and `docs/reviews/2026-11-25-H4.md`. The first run (#258, 2026-11-24)
  kept H4 open for CH4.8; the second (#262, 2026-11-25) found the condition met on both machine
  classes' trees and on W1–W4 with the constant density (locust and Kathleen job 672715). Signed off on
  #262 (the user, 2026-11-25: "Close H4 now (Recommended)", corrections R7–R11 and the triage "Accept
  all (Recommended)", the list "CHA.10a, then FFT, H2, H5 (Recommended)"). Closed with CH4.0f, CH4.0g,
  CH4.0k, CH4.0l, CH4.2d, CH4.3f, CH4.6 and CH4.7 open (the first run's R1, the user, 2026-11-24).

## H5: Complex wavenumber

**Gate:** CH5.3, complex k against FMM3D.

- [ ] **CH5.1 Classifier and pruning for Im k > 0.** *Done when:* regimes as §3.6 with S6's
  damped results; the error within target over a grid of (Re kD, Im kD) on W1. Also (CH4.3a's
  report §5, its question 2, signed off 2026-11-20): a damping factor of the bound fitted
  in the input-relative measure (CH1.3c's weighting; with e^{Im(k)(|t| + 2a)} E over-predicts by
  up to 10¹⁴ at η = 0.1 and 10²⁵ at η = 1, without it E under-predicts by up to 3.4·10⁶),
  confirmed in the FMM on damped W1. *Needs:*
  CH2.3, CH4.3.
- [ ] **CH5.2 Damping-dependent truncation and bound.** *Done when:* the error within target
  for weakly damped cases. Also (the H4 gate's second run, the user on #262, 2026-11-25: "Accept all (Recommended)"; review 2026-11-25, A12): checked with the
  constant density beside the random charges, on W2, at all six ε of D1, with the box measure beside
  ℓ2. *Needs:* CH5.1, CH0.6.
- [ ] **CH5.3 Gate: complex k against FMM3D.** *Done when:* both codes against the same
  sampled direct sum over the (Re kD, Im kD) grid, and the full-target crosscheck. The
  comparison with FMM3D is limited for Im kD ≳ 50 as in CH2.5 (CH1.4 report §5,
  question 2, signed off 2026-10-28). *Needs:* CH5.2, CH1.4. *Machines:* Kathleen.
- [ ] **H5.G Gate.** `/milestone-gate H5`. *Needs:* every task above.

## H7: Device path

**Gate:** CH7.6, the device path against FMM3D at matched accuracy.

- [ ] **CH7.1 Rotation and translation kernels on complex planes.** Real coaxial steps of
  M2M and L2L, the complex (S|R) step of the M2L, the stacked GEMM where Dense applies.
  *Done when:* device output within the host's error. Also (review 2026-11-12, M9): in
  f32 the stored low-frequency quantities carry 1/(2n + 1)!!, about 10⁻⁴¹ at n = 30, below
  f32's normal range, and the raised degrees are 28 and 32; the degrees f32 supports
  stated and refused beyond. Also (the H4 gate's triage, the user on #258, 2026-11-24: "Accept all (Recommended)"; review
  2026-11-24, L9): new kernels in the single-body `#[comptime]` form of
  `tetraleaf_kernels::plane_wave`, not per-layout pairs. *Needs:* CH3.1, D3. *Machines:* locust.
- [ ] **CH7.2 Device P2P for Helmholtz and Yukawa.** *Done when:* an accuracy contract as the
  host's; cost per pair against the host and Laplace's device P2P. Also (review 2026-11-19,
  A7; corrected at the H4 gate for c_leaf 4, review 2026-11-24, M17): in f32 the phase bound at
  the default c_leaf = 4 with V pairs as near pairs (Re(k) r ≤ 176, about 1.3·10⁻⁴ per term,
  above D3's ε ≥ 10⁻⁴): the f32 contract decided with c_leaf or the V pairs as near pairs.
  *Needs:* CH3.4.
  *Machines:* locust.
- [ ] **CH7.3 Device leaf operators.** P2M, L2P, M2P, P2L with radial recursions. *Done
  when:* against the host operators within error. *Needs:* CH3.1. *Machines:* locust.
- [ ] **CH7.4 Diagonal M2L and the transforms on the device.** cubek-fft where D10 adopted it,
  else the transform as a GEMM. *Done when:* as CH4.2 on the device; device locals against
  the host's within error. *Needs:* CH4.2b, CH7.1, CH0.5, D10. *Machines:* locust.
- CH7.5 Plane-wave M2L on the device: Laplace's first. Added on 2026-11-05 at the user's
  request, splitting CH7.5: the device kernels of merge and shift (the complex products over
  a box's waves, the merge lists per parent, the shifts to the children), built first on
  Laplace, whose host strategy exists (Phase 6 T13, `M2lStrategy::PlaneWave`) and which on
  a device ran as `Dense` until CH7.5c ("PlaneWave, run as dense on the device",
  `tetraleaf-exec/src/device.rs`); CH7.5d reuses them. The task T13's brief left for afterwards
  (docs/phase6/T13-plane-wave-m2l.md, "Afterwards"): T10 modelled it at 2.4–20× T7's dense
  M2L stage on the H100, larger at high p (`spikes/plane-wave/SPIKE_REPORT.md` §3), with the
  irregular merge lists the open device question. Split on 2026-11-05 (the user, "a/b/c,
  Helmholtz → d"; constraint 1 agreed: new kernels, the device's `PlaneWave` opt-in until
  CH7.5c's sign-off) into:
- [x] **CH7.5a The plane-wave kernels.** Record: `docs/helmholtz/tasks/CH7.5a.md`; tests
  `tetraleaf-kernels/tests/kernels/translate/plane_wave.rs`. `tetraleaf-kernels` only (with read-only
  accessors of `tetraleaf_tables::PlaneWaveTables`), module `plane_wave`: the
  steps of `tetraleaf_tables::plane_wave` (CONVENTIONS §3.14) as kernels, one launch per step
  and level call (the outgoing waves per box and frame, the merge per parent and
  accumulator over a merge list built on the host, the gather with the child shifts, the
  local with the rotation back), the tables uploaded once from plain arrays, explicit `fma`
  chains, no atomics. *Done when:* each kernel bit for bit a host `mul_add` replica of its
  order and within a stated bound of the table function it restates, f32 and f64 on the
  CubeCL CPU runtime (CI), f32 on Metal (the M6, by hand), both on CUDA (locust, by hand),
  every layout; a level call over a full and an adaptive level against
  `PlaneWaveTables::m2l` per pair to rounding and against `tetraleaf_ref::direct::m2l` on
  levels 2, 9 and 16 within T13's quadrature bound (p ≤ 20 f64 in an ignored sweep); two
  launches bit for bit. No `tetraleaf-exec` change. *Needs:* Phase 6 T13 (merged).
  *Machines:* locust, the M6 (Metal).
- [x] **CH7.5b The device's `PlaneWave`.** Record: `docs/helmholtz/tasks/CH7.5b.md`;
  `fmm::DevicePlaneWave`, `FmmBuilder::{device_plane_wave, device_plane_wave_budget}`,
  `device::PlaneWaveReport`, `PlaneWaveScratch::covering`; tests
  `device_common::check_plane_waves`, `tests/mpi_exec.rs` (`check_plane_wave`), the gate in
  `tests/device_fmm.rs`. "The CubeCL CPU runtime in CI" met by hand for `tetraleaf-exec` (no
  CI job builds it with `cpu`; the kernels' tests run in CI): signed off 2026-11-05 as
  recommended. `tetraleaf-exec`'s device path runs M2L under
  `PlaneWave` on CH7.5a's kernels when a builder setting names it (otherwise "PlaneWave,
  run as dense on the device", as before): the merge lists from the V views at build, the
  waves buffer bounded by chunks of parents within the device's scratch budget, the
  launches in the formula of `tests/device_common`; the lists built in the order of CH7.5a's
  test (`Lists::new` in `tetraleaf-kernels/tests/kernels/translate/plane_wave.rs`, the order of
  `src/operator_plane_wave.rs`), several vectors as columns, each chunk's buffers within
  the kernels' 32-bit index range (CH7.5a). *Done when:* `PlaneWave` on the device
  (CUDA f32 and f64, Metal f32, the CubeCL CPU runtime in CI) within 1.1× `Dense`'s errors
  against the direct sum at the C3.2 cube and the C3.3 Plummer sphere (T13's requirement 1,
  f64 p = 8, 12, 18 and f32 p = 8), device locals against the host `PlaneWave`'s within the
  stated bound, bit for bit from run to run and on P ranks against the one-rank device run;
  constraint 1: without the setting every strategy, the device's `PlaneWave` and Laplace's
  defaults bit for bit, the device `Auto` unchanged. *Needs:* CH7.5a. *Machines:* locust,
  the M6 (Metal).
- [x] **CH7.5c Measured, and adopted or dropped.** Record: `docs/helmholtz/tasks/CH7.5c.md`;
  report `tetraleaf-validate/results/helmholtz-ch7-5c.md`; `fmm::DevicePlaneWaveLayout`,
  `FmmBuilder::device_plane_wave_layout`, `device::{default_plane_wave_budget,
  CUDA_PLANE_WAVE_WIDE_UNITS, CUDA_PLANE_WAVE_WIDE_FROM, METAL_PLANE_WAVE_BUDGET}`,
  `layout_sweep --part planewave`. Adopted (signed off 2026-11-06: "sign off all three
  questions as recommended"): `DevicePlaneWave::Kernels` the default of `PlaneWave` on every
  device (`Dense` the named opt-out), cubes of 128 units on Metal and on CUDA to p = 12 and
  of 256 on CUDA from p = 13, Metal's plane-wave budget 1 GB, the kernels kept; the device
  `Auto` unchanged (CH7.5e); D9 by the device executables and gates on the CPU runtime,
  Metal and CUDA and one ABBA run of `tetraleaf-bench --backend cuda` against `main`.
  The device M2L stage and the evaluation
  under CH7.5b's setting against T7's dense M2L at f64 p = 8, 12, 18 on locust's H100 (f32
  p = 8 on the M6's Metal), with the kernels' layouts measured (CH7.5a's `GPU_UNITS` = 128
  is a first choice; the merge's irregular lists are the spike's open question) and the
  chunk budget (CH7.5b: `device_plane_wave_budget`, default the translation scratch
  budget; CH7.5b's gate printed, unmeasured, the cube at f64 p = 8 on the H100 slower than
  `Dense` and the Plummer sphere faster); Laplace's
  device gates and benchmarks within noise (D9). *Done when:* the report gives both; adopted as the device's `PlaneWave` where it wins at
  p ≥ 12 (T13's gate; the setting becomes the default of `PlaneWave` on the device), else
  `PlaneWave` stays `Dense` on the device and the kernels are dropped; the device `Auto`
  unchanged (a change to it is a sign-off). *Needs:* CH7.5b. *Machines:* locust, the M6.
  *Sign-off:* adoption on the device.
- [ ] **CH7.5d Exponential M2L on the device**, if D8 adopted it on the host. On CH7.5a's
  merge and shift kernels, with the Helmholtz conversions (CH6.1c's pair chain) and the
  quadrature per level (CH6.1e). *Done when:* as CH6.1 on the device. *Needs:* CH7.5c,
  CH6.2, CH7.1, D8. *Machines:* locust.
- [ ] **CH7.5f CubeCL 0.11.3 and cubek 0.3.1.** Added on 2026-11-06 (the user: "Please
  schedule right before CH7.5e"; off the priority list), so that CH7.5e's device baselines
  are taken on the new releases. The workspace pins `cubecl` `=0.11.0`, `cubek-matmul` and
  `cubek-std` `=0.3.0` move together to `=0.11.3` and `=0.3.1` (root `CLAUDE.md`, "Working
  rules"; `run-dependency-checks` lists them as exempt, in #193). The releases fix none of
  `docs/design/cubecl-upstream.md`'s entries; they bring LLVM for the native CPU on the CPU
  runtime (cubecl #1833, which may move the CPU runtime's results in the last bits and
  changes the `tracel-llvm` bundle, 23.1.0-3 to 23.1.0-4), a parked idle device server
  (#1806, launch latency), a CUDA PTX placeholder fix (#1691), loop conditions in their own
  scope (#1788) and the plane matmul's planes per cube floored at one (cubek #784).
  *Done when:* the pins and `Cargo.lock` moved, with the `tetraleaf-kernels` checks of root
  `CLAUDE.md` passing on the CPU runtime and Metal (the M6) and CUDA (locust); every
  device test and gate passing, any changed bit identity explained; Laplace's device
  benchmarks and CH7.5c's merge and shift times within noise of `main` (D9) on the M6's
  Metal and locust's H100; `cargo deny` clean; `cubecl-upstream.md`'s "0.11.0: still
  present" entries re-checked against 0.11.3, and the pin's sentences in root `CLAUDE.md`
  updated. *Needs:* CH7.5c. *Machines:* locust, the M6 (Metal).
- [ ] **CH7.5e The device `Auto` on `PlaneWave`.** Added on 2026-11-06 (CH7.5c's question
  2, signed off as recommended; off the priority list): the device's `M2lStrategy::Auto`
  resolves to `Dense` or `Rotation` by the static rule and never to `PlaneWave`, although
  `PlaneWave` on merge and shift took 0.19–0.21× `Dense`'s M2L time on the H100 at f64
  p = 12 and 18 (CH7.5c). At N = 10⁵ the device gate's evaluation times (printed, not
  asserted) put merge and shift behind `Dense` on the cube (11.8 against 3.7 ms at f32
  p = 8 and 13.4 against 10.2 ms at f64 p = 12 on the H100; 16.5 against 8.2 ms on Metal)
  and ahead on the Plummer sphere, so the rule needs N as well. Measure the device's
  strategies against each other by backend, precision, p and N (the static rule's crossovers and the near rule's cost models on the
  device, as Phase 6 T11 did for the host's `Auto`) and propose a rule. *Done when:* the
  report gives the device `Auto`'s rule by backend, precision and p with `PlaneWave` in it
  where it wins, adopted on sign-off; Laplace's device gates and benchmarks within noise
  (D9). *Needs:* CH7.5c, CH7.5f. *Machines:* locust, the M6 (Metal). *Sign-off:* the device
  `Auto`.
- [ ] **CH7.6 Device performance against FMM3D at matched accuracy.** As CH4.5, the device
  path beside the host's. *Needs:* CH7.2, CH7.3, CH7.4, CH4.5. *Machines:* locust.
- [ ] **H7.G Gate.** `/milestone-gate H7`. *Needs:* every task above.

## H8: Distributed path

**Gate:** Kathleen scaling at kD of a few hundred, rank-invariant.

- [ ] **CH8.0 One copy of the M2L tables per node.** Added on 2026-11-12 at CH4.0 (question 1
  of #211; the user: "With regards to a shared table per node can you add a task to the roadmap
  to do this?"). Every rank builds and holds its own copy of the tables that are the same on
  every rank: Helmholtz's dense, class and Batched tables per level, and Laplace's dense
  tables. Since CH4.0 Helmholtz's share a node budget (8 GiB), so at 40 ranks per node a rank
  holds 205 MB and the levels beyond it run Rotation; threads per rank are today's remedy
  (node-m2l.md §7). Laplace has no shared copy either: Phase 5N's spike measured one in an
  `MPI_Win_allocate_shared` window (node-m2l.md §2.6: Dense's per-pair cost 2.4–5.8× below
  private copies at one rank per core, about what 4 × 10 threads gives), and its decision 6
  left it out ("no shared window in this phase", to revisit in Phase 5S, which did not).
  - One rank per node (`split_shared`) builds the tables into a node window, the others read
    them after a barrier; the node budget then counts one copy per node, not one per rank.
  - The FFI through `mpi-sys` (rsmpi 0.8.2 has no window API; the user's permission of
    2026-10-10): the crate and module that hold it named, every public function safe, a
    `// SAFETY:` comment on every block, root CLAUDE.md's `unsafe` rule amended in the same PR.
  - Laplace's tables in the same window behind an option, its default and outputs unchanged
    (constraint 1), or left out with the reason.
  - Also (added at CH4.3b2c, its question 4, signed off 2026-11-22): the diagonal form's transfer
    functions and conversions, built on every rank (about 10 s at 40 × 1 on Kathleen), and its
    signatures' scratch budget per rank, which at one rank per core keeps the form off W1 at
    kD = 1000's levels 2–3 and W3 at kD = 500's level 5 (40 × 1) and levels 2–5 (72 × 1 on
    locust; CH4.3b2c's report §7, §8).
  *Done when:* the outputs bit for bit the private tables' at 1, 2 and 4 ranks (CI) and at
  40 ranks per node on Kathleen; the resident memory per node at 40 × 1, 10 × 4 and 1 × 40
  against private copies; the evaluation at 40 × 1 against 4 × 10 and 1 × 40 on W1 at
  kD = 100 and 500 at 10⁻⁶, where the node budget now sends levels to Rotation; a default
  proposed (the window on or off, the node budget with it). Also (CH4.1b's report §3, question
  1, signed off 2026-11-19 on #236): a coarse level's source signatures of the diagonal form formed once per node and shared
  (every rank's V rows on levels 2–3 need nearly every source: 413–433 MiB per rank at
  kD = 1000), so that the diagonal form keeps those levels at many ranks per node. *Needs:*
  CH4.0. *Machines:* Kathleen, locust. *Sign-off:* the `unsafe` exception and the default.
- [ ] **CH8.1 Byte cap for the replicated top; exchange of the high-frequency levels.**
  *Done when:* bit for bit the one-rank output on P ranks; memory per rank and per node
  bounded. *Needs:* CH4.4, CHA.7. *Machines:* Kathleen.
- [ ] **CH8.2 Kathleen scaling at kD of a few hundred.** *Done when:* efficiency and memory
  reported as in Phase 5S, at most 4 nodes per job. *Needs:* CH8.1. *Machines:* Kathleen.
- [ ] **H8.G Gate.** `/milestone-gate H8`. *Needs:* every task above.

## H9: The final comparison with FMM3D

**Gate:** CH9.1, our FMM with every optimisation and default beside FMM3D, in depth.

- [ ] **CH9.0 Helmholtz dipole sources.** D5's task (signed off 2026-10-29: "dipoles get a
  later, separate task"), added at the H1 gate (G2; the user, 2026-11-19, on #230: "Add CH9.0"). Dipoles
  (Burton–Miller, double-layer operators) through the source-type interface of Laplace's C6.4
  (Phase 6 T8, T8d). *Done when:* CONVENTIONS §3.15's dipole convention (a version bump,
  proposed in the PR first); the dipole direct sum and leaf operators in `tetraleaf-ref` against
  the direct sum and k → 0 against Laplace's dipoles; the host FMM's dipole source type within
  ε on W1 and W3 for real, complex and imaginary k; FMM3D's dipoles as a check where it meets
  its eps. Split into lettered sub-tasks before it starts. *Needs:* H3.G, CH1.3c.
- [ ] **CH9.1 The final comparison with FMM3D.** Added on 2026-11-07 (the user: "can you
  schedule an in depth comparison with fmm3d when all optimizations are applied?"), after the
  last optimisation: H4 to H8 closed and every default signed off (D7, D8, the device
  `Auto`). Our FMM with its defaults beside `hfmm3d` (the `default` and `fastker` builds,
  and CH6.2d's patched `timed` builds for the breakdown) on one node of the M6, locust and
  Kathleen, as Laplace's T9 (`tetraleaf-validate/results/phase6-fmm3d.md`) and CH4.5, in depth:
  - accuracy: for each of our tolerances (10⁻³ to 10⁻⁸, D1) FMM3D at the same requested ε
    and at the eps of the closest error (matched accuracy); both against the same sampled
    direct sum, with the full-target crosscheck and the worst box;
  - the reference workloads W1–W4 (D12) over D1's range (kD ≤ 500), complex k and Yukawa
    where FMM3D meets its eps (CH2.5, CH5.3), real and complex charges, potentials and
    gradients, several vectors (FMM3D's vectorised interface), dipoles (CH9.0);
  - time: the build (FMM3D's setup) apart from the evaluation, and the evaluations after
    which ours is ahead; per stage and per level (FMM3D's from CH6.2d's patch: the M2L per
    level, lists 3 and 4, P2P), at one thread and on the whole node (FMM3D's OpenMP on every
    core, ours at the node's best ranks × threads split and at 1 × N, review 2026-11-19, B6),
    with FMM3D's OpenMP binding stated (B8), each machine quiet for the whole
    run (no other jobs; locust's load sampled during the runs);
  - memory: the peak resident size per node;
  - the device path (CH7.6) beside the host's, and our multi-node runs (H8) labelled apart:
    FMM3D runs on one node.
  *Done when:* the report has these tables on the three machines; every case where FMM3D is
  more than 1.1× faster at matched accuracy is traced to its stage and level, with its cause
  or a follow-up task; CH4.5's and CH7.6's results cited where unchanged. Also (the H4 gate's
  triage, the user on #258, 2026-11-24: "Accept all (Recommended)"; review 2026-11-24): both FMM3D variants on every case, the density cases
  included (B11); the match also as the fastest eps within 5 % of the closest error, beside the
  closest (B13); W1 with the constant density at kD = 350 and 700, three ε, off the degree
  rule's calibration octaves (M16); Laplace at f64 p = 8 and 12 under its defaults (PlaneWave)
  beside `lfmm3d`, with a sphere's surface and a constant density beside the random charges
  (B14, B15). Also (the H4 gate's second run, the user on #262, 2026-11-25: "Accept all (Recommended)"; review 2026-11-25): W2 with the constant density at kD = 750 and
  1000 with the box measure, where the worst box's growth with kD (29–76 ε from kD = 175 to 500)
  may pass 100 ε (M21); W1 with the constant density at kD = 100.53 and 502.65, near the single
  layer's zeros kD = 2πn (M22); CH4.5's comparison and CH4.4b's mapping by level and class under
  the default, which CH4.8b raised by one degree after the cap (B16; the gate's §0.4). *Needs:* H4.G, H5.G,
  H6.G, H7.G, H8.G, CH9.0. *Machines:* locust, Kathleen, the M6.
- [ ] **H9.G Gate.** `/milestone-gate H9`. *Needs:* every task above.

## Architecture and cleanup (HA)

From the architecture pass of the project reviews (`docs/reviews/`; the first,
`docs/reviews/2026-11-12-H6.md` §2, after the H6 gate), by its triage of 2026-11-12
(questions 3 and 6 of #209, the user, 2026-11-12, on #209: "I agree": the section and its IDs as written,
CHA.1 to CHA.3 on the priority list, CHA.7 a *Need* of CH8.1, CHA.6 recommended before
CH7.5d but not a *Need* of it, no issue). Not a milestone: no gate, no issue, and
`/milestone-gate` does not close it; each later review adds to it.

Every task here keeps constraint 1 of the proposal's §2 in its own *Done when*: Laplace's
outputs bit for bit, its benchmarks within D9's spread wherever a host or device path's code
moves, `/verify`'s Laplace protection run. No behaviour or default changes, and a public
item is removed only where the task names it.

- [x] **CHA.1 The spikes out of the workspace.** Review K3, U1. Record: `spikes/README.md`,
  `docs/helmholtz/tasks/CHA.1.md`. All twelve kept, a workspace of their own (`spikes/Cargo.toml`,
  one `spikes/Cargo.lock`, the CubeCL pins repeated, green-kernels, pulp and the `[patch]` moved
  there); the root lock lost 30 packages, no version changed. Signed off 2026-11-19 (the user:
  "Keep all, none archived (Recommended)", "One spikes workspace (Recommended)").
- [x] **CHA.2 `unsafe` enforced.** Review C1. Record: `docs/helmholtz/tasks/CHA.2.md`.
  `[lints.rust] unsafe_code` in each `Cargo.toml` (every target): `forbid` in seven crates,
  tetraleaf-octree included, `deny` in tetraleaf-plan, whose one exception (the example's counting
  allocator) is allowed in its file, named in root CLAUDE.md. `points_to_morton` made safe (1–2 ms
  per 10⁷ points). Signed off 2026-11-19 (the user: at the plan "[lints.rust] per crate
  (Recommended)", "M6, now (Recommended)", "Skip D9 runs (Recommended)"; on #225 "Please go ahead
  with forbid.").
- [x] **CHA.3 The superseded M2L configuration removed.** Review P1. Record:
  `docs/helmholtz/tasks/CHA.3.md`, `docs/design/node-m2l.md` §11. `NodeContext` keeps p and the
  value's bytes; `FmmBuilder::dense_node_budget`, `AUTO_DENSE_*`, `resolve(p)` and `scaling
  --budget` gone; `Auto` the same on all 42 (precision, p) cells. Signed off 2026-11-19 at the
  plan (the user, all four as recommended: "Go ahead", keep `one_rank`, "Drop budget_mb only",
  "Skip D9 runs").
- [ ] **CHA.4 A test-support crate.** Review T1, C3. A small crate without workspace
  dependencies (the generator, the point sets, the error norms, the machine facts) that
  every crate may dev-depend on and tetraleaf-bench depends on in place of tetraleaf-validate
  where that is all it needs. *Done when:* the copies of the generator (41 at the first review,
  52 files at the second, review 2026-11-19) and the copies of
  `relative_l2` (eleven), `plummer` and `cube` replaced where the draws stay the same (checked per
  file by the tests' own fixed values); `workspace-structure.md` has the crate. Also (the H4 gate's second run, the user on #262, 2026-11-25: "Accept all (Recommended)"; review 2026-11-25,
  K9, T1): the device examples' shared helpers (`time` and `warm_up`, byte-identical in three and four
  of `m2l_kernels.rs`, `rotation_kernels.rs`, `gemm_variants.rs`, `translation_kernels.rs`; the spike
  reference constants; `relative_l2`) in one `gpu`-gated module of tetraleaf-validate; the MPI tests'
  `union`, `share` and `union_points` (`mpi_exec.rs`, `helmholtz_fmm.rs`, `multi_rank.rs`,
  `device_ranks.rs`) in one place.
- [ ] **CHA.5 The Helmholtz builder's measurement settings apart.** Review P2, P5. *Done
  when:* `HelmholtzFmmBuilder` shows the user's settings; the measurement settings
  (`plane_wave_split`, `auto_plane_wave`, `m2l_plane_wave_min_degree`, `m2l_max_degree`,
  `block_budget`, `PlaneWavePrecision::PerLevel`, `pair_calls`) sit behind one hidden
  tuning struct; `plane_wave_reach` and `PlaneWaveRaise::Level` removed if their guard
  tests (CH6.2g, CH6.1h) can keep their checks without them, else kept there; the harness
  and the reports' reproduction scripts still run. Also (review 2026-11-19, P2, P6, P7): the
  settings added since (`empty_boxes`, `near_v`, `near_wx`, `batched_term`, `m2l_block_targets`,
  `rcr_split`, `rcr_budget`, `plane_wave_budget`, `dense_node_budget`, `c_leaf`) sorted the same
  way, the on/off switches of signed-off rules behind the tuning struct; `Plan`'s
  `with_near_rule` and `from_key_types_with_near_rule` folded into the `with_rules` forms;
  `HelmholtzOperator`'s unbuilt state a type that `build_tables` consumes, not a run-time
  `bool`. Also (the H4 gate's triage, the user on #258, 2026-11-24: "Accept all (Recommended)"; review 2026-11-24, P2, P9,
  P10, K6): H4's settings sorted the same way (`diagonal`, `diagonal_build`, `diagonal_costs`,
  `diagonal_budget`, `diagonal_chunk`, `degree_rule`, `v_after_plan`, `v_raised_price`,
  `v_costs`, `m2l_overrides`); `VAfterPlan` removed if CH4.3b2d's guard test keeps its check
  without it, `VRaisedPrice::PerPair` removed, `TransferBuild::Direct` kept as the interpolated
  build's oracle and documented so; `VNear`'s and `WxNear`'s `#[default]` moved to `CostModel`,
  the builder's; one Helmholtz settings parser in `tetraleaf-validate/src` for
  `helmholtz_accuracy`, `helmholtz_boxes` and `helmholtz_fmm3d`. Also (the H4 gate's second run, the user on #262, 2026-11-25: "Accept all (Recommended)"; review 2026-11-25, P12): the
  superseded degree rules (`DegreeRule::CH3_3`, `CH4_3C2`, `CH4_3C3`), used by no library path of
  tetraleaf-exec, and one rule-name parse (now copied in three examples) moved to
  `tetraleaf_validate::helmholtz_calibration`; exec keeps `CALIBRATED` and `gd`. *Sign-off:* the settings
  removed.
- [ ] **CHA.6 One plane-wave M2L traversal.** Review L2. The slot table, the parent runs
  and the split by child of `operator_plane_wave.rs` and
  `helmholtz/operator_plane_wave.rs` in `driver.rs`, generic over a small kernel trait
  (outgoing, merge for the parent, convert for the child), as `HostKernel` is. *Done when:*
  both operators bit for bit; Laplace within D9 on the three machines; the two split
  thresholds (a quarter, a half) kept as each operator's own constant, or unified by a
  measurement. Recommended before CH7.5d, which restates the traversal on the device
  either way. Also (the H4 gate's triage, the user on #258, 2026-11-24: "Accept all (Recommended)"; review 2026-11-24, L2): one
  slot map (the level's distinct sources to slots) and one chunk budget (bound, peak) for the
  three chunked Helmholtz M2L buffers (`PlaneBuffer`, `RcrBuffer`, `DiagonalBuffer`).
  *Machines:* locust, Kathleen.
- [ ] **CHA.7 One distributed build for both facades.** Review L1, L3. The steps of
  `FmmBuilder::build` (the settings' agreement, the node split, the plan, the near rule and
  the second cut, the routing by node, the evaluator, the timings) in `pipeline`, behind a
  crate-private trait the two builders implement; `ExecOperator`'s kind and stage timings
  generic over the operator. *Done when:* Laplace bit for bit and within D9; `HelmholtzFmm`
  routes by node (`Evaluator::with_routing` with the node of its `split_shared`), with the
  work-weighted cut and the timings, rank-invariant as before; what CH8.1 still has to add
  listed. Also (review 2026-11-19, L1, P8): the Helmholtz builder's second cost model
  (`VCosts`, `WxCosts`) beside `cost.rs`, and its copies of `agree`, `Calls` and `delegate!`;
  `kind_times`' literal indices replaced by an operator kind. Split into lettered sub-tasks
  before it starts. Also (review 2026-11-24, L1): since H4 `HelmholtzFmmBuilder::build` is 680
  lines with 17 collectives, builds its plan once on a shared index with the empty-box, V and
  W/X rules, where Laplace builds it up to three times without them. *Machines:* locust, Kathleen.
- [ ] **CHA.7a The Helmholtz build's steps named.** Added at the H4 gate (the user on #258, 2026-11-24: "Accept all (Recommended)";
  review 2026-11-24, L1), CHA.7's first letter, before H8. The agreements and steps of
  `HelmholtzFmmBuilder::build` as named crate-private steps of `pipeline.rs` (the settings'
  agreement, the node split, the tree, the plan with its rules, the operator's tables, the
  evaluator), the facades not yet merged; the copies of `agree`, `Calls` and `delegate!` one.
  *Done when:* every Helmholtz output bit for bit (the MPI list at 2 and 4 ranks, `invariance`,
  a `helmholtz_accuracy` matrix by checksum); Laplace untouched; `build` under 200 lines.
  *Machines:* Kathleen (the checksums at 40 × 1 and 1 × 40).
- [ ] **CHA.8 The large files of tetraleaf-exec split.** Review C2, K4, T2, U2. `device.rs`
  into modules and `DeviceOperator::new` into staged constructors; `fmm.rs`'s option enums
  and report types into submodules; `tests/mpi_exec.rs` into `tests/mpi_exec/main.rs` with
  one module per scenario, and the same for `tetraleaf-exec/tests/helmholtz_fmm.rs` and
  `tetraleaf-plan/tests/mpi_regressions.rs` (review 2026-11-19, T2), `helmholtz_fmm.rs` first
  (review 2026-11-24, T2: 4,826 lines, one `#[test]`). Moves only. *Done when:* bit for bit; within D9 with the
  placement checked (`-align-loops=64` where a time moves); one PR per file.
  *Machines:* locust, Kathleen.
- [ ] **CHA.9 Documentation of what is, not of how it came.** Review K1, L4. *Done when:*
  the rustdoc of public items says what the item does and why its default is what it is,
  the history in `progress.md` and the phase documents; each crate CLAUDE.md cut to its
  rules and checks (the time rule of `tetraleaf-exec/CLAUDE.md` restated to what holds);
  `docs/design/workspace-structure.md` has the Helmholtz modules; the crate descriptions
  and `tetraleaf-exec/src/lib.rs` name Helmholtz; `laplace-fmm-plan.md` §2.5's 0.76ᵖ named as
  the M2L's ratio (M10). Also (review 2026-11-19, C4): `tetraleaf-ref/CLAUDE.md` and
  `tetraleaf-exec/CLAUDE.md` record as a decision that the host FMM's leaf operators call
  `tetraleaf-ref`'s per point (no second trusted path there; a change to the oracle moves the
  FMM's times). Also (the H4 gate's triage, the user on #258, 2026-11-24: "Accept all (Recommended)"; review 2026-11-24, K1, C9,
  M17): tetraleaf-exec first (its CLAUDE.md of 1,538 lines, `lib.rs`'s "the Laplace FMM");
  `workspace-structure.md`'s graph with the exec → ref and exec → math edges; the device path's
  downcast (`device.rs`, `Box<dyn Any>`) recorded as the one place; the proposal's §3.6 table
  with §3.5's current degrees (128, 70, 40). Comments and documents only, one crate per PR.
- [ ] **CHA.10 Housekeeping.** Review K2, U3. *Done when:* the one-off scripts of
  `tools/helmholtz/` and the task-named cases of `tools/scaling/run.sh` are apart from the
  living tools, each report's reproduction line still valid; `tetraleaf-octree/Cargo.toml` takes its
  versions from `[workspace.dependencies]` (the migration the root CLAUDE.md leaves open;
  *Sign-off*) and points at the GitHub repository; its unused features removed if no
  example or test names them (`battleship` is empty). Also (review 2026-11-19, K2, K5): a
  shared `tools/helmholtz/common.py` for new report scripts (the `*-tables.py` define `cells`
  11 times), the old ones left; `approx` dropped from `[workspace.dependencies]`, which no
  member uses, and from the two CLAUDE.md files that allow it. Also (the H4 gate's triage,
  the user on #258, 2026-11-24: "Accept all (Recommended)"; review 2026-11-24, K2, K7, U6): `common.py` first, before the other
  housekeeping (69 scripts and 12,232 lines, fourteen `*-tables.py` defining `cells`); `mpi`'s
  `complex` feature in `tetraleaf-plan/Cargo.toml` dropped if nothing needs it; the crate
  descriptions naming Helmholtz; `tetraleaf-validate/Cargo.toml`'s `[[example]]` list complete
  (ten examples are auto-discovered). Also (the H4 gate's second run, the user on #262, 2026-11-25: "Accept all (Recommended)"; review 2026-11-25, C11, P13, U7, K7): `common.py` split off
  as CHA.10a; in tetraleaf-octree the empty `types` module, the `battleship` feature and the
  unreferenced `_test_sphere.vtk` removed, `tools` (a test key generator, the library's only use of
  `rand` and `rand_chacha`) behind `cfg(test)` or a feature; `Plan`'s constructors cut to `new`,
  `with_rules`, `with_rules_on` and `from_key_types_with_rules` (beside CHA.5's fold);
  `run-dependency-checks.yml` given a `timeout-minutes`, `submodules: "true"` dropped from it and
  `run-examples.yml` (no `.gitmodules`), `run-examples.yml` installing through `apt-install.sh`.
- [x] **CHA.10a The report scripts' shared module.** Review 2026-11-25, K2. Record:
  `docs/helmholtz/tasks/CHA.10a.md`. `tools/helmholtz/common.py` (`cells`, `tables`, the parsers of
  `helmholtz_accuracy`, `helmholtz_boxes` and `helmholtz_fmm3d compare`, `gmean`; a `selftest`),
  Python 3.9 standard library; `ch4-8b-tables.py` on it, its 34 outputs byte-identical; `/next-task`
  sends new report scripts to it. The old scripts left as they are.
- [x] **CHA.11 The Helmholtz wave levels in one place.** Review 2026-11-19, L1. Record:
  `docs/helmholtz/tasks/CHA.11.md`. `operator::wave_levels` chooses the levels with waves for
  `select_tables` and `HelmholtzFmmBuilder::build`; the restatement, its assert and the second
  f64 gate gone. Bit for bit on 348 `helmholtz_accuracy` rows and the MPI list. It also makes
  named precisions with the V rule run where a level lies outside the quadrature's tables, which
  panicked on `main` (#233's question 1; the user, 2026-11-19: "Accept as recommended"). H4's
  diagonal-form choice per level goes in the same place.
- [ ] **CHA.12 One entry-point macro set in tetraleaf-simd.** Review 2026-11-19, L5. The eleven
  macros that `arch/avx2.rs`, `avx512.rs` and `neon.rs` each define after their `Simd` impls
  (about 670 lines each, about 140 differing) once in `arch/mod.rs`, taking the token, the
  `target_feature` literal and the block constants. *Done when:* every kernel bit for bit on
  both architectures (CI's x86_64 and arm64 jobs); D9's protocol on the three machines
  (code placement). *Machines:* locust, Kathleen.
- [ ] **CHA.13 The machine-class constants in one place.** Added at the H4 gate (the user on #258, 2026-11-24: "Accept all (Recommended)"; review 2026-11-24, C8). The constants keyed on `Isa` as a proxy for the machine class
  (`V_COSTS_*`, `WX_COSTS_*` in `helmholtz/fmm.rs`, `DIAGONAL_COSTS_*` in
  `helmholtz/operator_diagonal.rs`, the Auto thresholds of `helmholtz/operator.rs`,
  `default_max_points_per_leaf`'s choice) in one crate-private `machine` module of
  `helmholtz/` with one `MachineClass::of(isa)`, moved, not changed; Laplace's `CALIBRATIONS`
  (`cost.rs`) named beside it. *Done when:* every output bit for bit on both architectures (CI's
  x86_64 and arm64 jobs and the MPI list); Laplace untouched.
- [ ] **CHA.14 The octree's neighbour map built on first use.** Added at the H4 gate's second run
  (the user on #262, 2026-11-25: "Accept all (Recommended)"; review 2026-11-25, C10).
  `Octree::new` and `repartition` (`tetraleaf-octree/src/octree.rs:629`, `:797`) build
  `compute_neighbours`' `HashMap<MortonKey, Vec<MortonKey>>` (up to 26 keys per key), which no
  product crate reads (tetraleaf-plan builds its lists from `all_keys()`). First measure its time and
  memory at N = 10⁶ and 10⁸ (Kathleen, the octree's build by part); if it pays, build it on the first
  `neighbour_map()` call (`OnceLock`) or behind an option the FMM leaves off. *Done when:* the
  measurement reported; if changed, Laplace's and Helmholtz's outputs bit for bit on host and P ranks
  (constraint 1) and the octree's examples and tests unchanged in what they check. *Needs:* none.
  *Machines:* Kathleen.
