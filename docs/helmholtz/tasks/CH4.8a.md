# CH4.8a The constant density on W2: the error mapped and a remedy chosen

Started 2026-11-24 on `tbetcke/h-ch4-8a`. Roadmap: "CH4.8 The constant density on W2" (H4;
priority list item 11, added at the H4 gate; its *Needs*, CH4.4b, merged as #256), split here
into CH4.8a (this task: the mapping, the cause, the remedy chosen on FMM and pair evidence, for
sign-off) and CH4.8b (the remedy as the default and the verification matrix under it, for
sign-off), as CH4.4b left a degree-rule remedy to "a follow-up task". Proposal: §10 CH4.4 (the
row whose bound CH4.8 restores); §3.5 (the degree rule, amended at CH4.3c2 and CH4.3c3); §9.1
(the box measure, its 100 ε bound within D1, amended at CH4.0b and CH4.4b); §2 (constraints 1, 2,
6, 9). Builds on the H4 gate's run (`docs/helmholtz/gates/H4.md` §9, Kathleen job 672173; review
`docs/reviews/2026-11-24-H4.md`, M15), CH4.4b (the box study under the defaults,
`tetraleaf-validate/results/helmholtz-ch4-4b.md`) and CH4.3c3 (the refit on W1's constant density,
`results/helmholtz-ch4-3c3.md`).

## Goal
Find where and why W2 with the constant density misses ε (∇φ 1.28 ε in ℓ2 at kD = 500, 10⁻³) and
the 100 ε box bound (up to 143.4 ε) under the defaults, complete the list of misses on the
verification set's cases not yet run, and choose a remedy on FMM and pair evidence with its cost,
for sign-off; CH4.8b makes it the default.

## Acceptance (proposal §10, verbatim)
§10 CH4.4's row, whose bound this task restores: "FMM against the direct sum (sampled targets) on
W1 to kD = 1000, W2 to 500 and W3 to 200, beside `hfmm3d` against the same sum; the full-target
crosscheck mapped to levels and classes, with no error concentrated at a class boundary; time,
memory and build cost per k". With §9.1: ε the relative ℓ2 tolerance, the worst box at most 100 ε
within D1 (kD ≤ 500).

The roadmap's *Done when* of CH4.8, split:
- **CH4.8a (this task):** the error mapped to levels and classes (`helmholtz_boxes --study
  levels,m2l`, as CH4.4b) and its cause stated (the degree rule on the dense cap's levels, the
  tangent-point placement, the form or Batched at p = 7); the verification set's cases not yet run
  under the defaults run, so that every miss is known before a remedy is chosen; a remedy (the
  degree rule's high- or low-frequency branch refitted on W2's constant density, or a level's
  degree raised by rule) chosen on both measures (the H6 gate's rule): the FMM's errors on every
  case that misses and on W1's constant density, and the worst pair per level (spike S7's curves),
  with its cost (build and evaluation, one setting per launch on Kathleen); for sign-off.
- **CH4.8b:** the signed-off remedy the library's default; under the defaults W1–W3 with the
  constant density at kD = 10–500 within ε in ℓ2 at 10⁻³–10⁻⁸ and every worst box ≤ 100 ε, the
  random-charge cases of CH4.4a unchanged within ε; its cost on Kathleen and locust; for sign-off.

## Read first
- `docs/helmholtz/gates/H4.md` §9 (the miss), `docs/reviews/2026-11-24-H4.md` M15.
- `tetraleaf-validate/results/helmholtz-ch4-4b.md` (the method; the 114.7 ε box's remedy, every
  level +1), `helmholtz-ch4-3c3.md` §2–§4 (the constant density's degrees; the plateau unscanned on
  it), `helmholtz-ch3-3.md` (the plateau's calibration: W1, W3, W4, never W2).
- `tetraleaf-exec/src/helmholtz/degree.rs` (`DegreeRule`, `CALIBRATED`).
- `tetraleaf-validate/examples/helmholtz_boxes.rs` (`--study`, `--raise-*`, `--rule`).
- `~/projects/fmm-archive/kathleen/runs/h4-gate/` (the gate's twelve launches).

## What the record already says (before any run)
- The gate's run (`--study single`): at kD = 500, 10⁻³ the worst boxes of φ and ∇φ are the
  tangent point on z = 0, in the dense cap, every level's target near a box corner (|v|/r
  1.2–1.73 on levels 2–8); the error grows with the levels on whose box face a target lies (RMS of
  ∇φ 0.38 ε with none, 1.73 ε with six, 4.28 ε with seven). The dense cap (30°, ten times the
  density, about 42 % of the points and so of §9.1's samples) lies within z ≤ 0.067 D, beside the
  z = 0 face on levels 2–4; it adds levels 7–8 to the tree, on the plateau (p = 7 at 10⁻³, Batched;
  15 at 10⁻⁶; 32 under the waves at 10⁻⁸).
- The worst pair per level at the degrees as run (S7's curves, `ch4-4b-tables.py pairs`), W2 at
  kD = 500: at 10⁻³ 2.6–6.3 ε on the high-frequency levels 2–6 and **16 and 10 ε** on the plateau
  levels 7–8 (p = 7); at 10⁻⁶ 2.3–3.7 ε on 2–5 and 32–126 ε on 6–8 (the plateau's p = 15, 16, half
  of S7's p*, calibrated on the FMM's ℓ2, not per pair); at 10⁻⁸ 2.7–5.5 ε on 2–7, 13.9 ε on 8.
  The plateau's pairs are far from ε by design; whether they carry W2's miss is what the mapping
  decides.

## Changes
- No library change in CH4.8a (the remedy is CH4.8b's). No convention change.
- `tools/helmholtz/ch4-8a-runs.sh`: the machine runs' stages (below), after
  `tools/helmholtz/ch4-4b-runs.sh`; the report's tables from `tools/helmholtz/ch4-4b-tables.py`
  (`summary`, `levels`, `m2l`, `worst`, `pairs`) on the outputs, with a script of this task's own
  only for what those lack (`tools/helmholtz/ch4-8a-tables.py`, if needed).
- If a candidate remedy cannot be stated as `--rule γ,α,β,C,a,b,d[,cap]` or by `--raise-levels`
  (for example a degree raised by the tree's depth), a measurement setting in the harness to run
  it (`tetraleaf-validate` only), with a test; stated here before it is written.
- Report `tetraleaf-validate/results/helmholtz-ch4-8a.md`.
- The roadmap: CH4.8 split into CH4.8a and CH4.8b (this commit), the priority list's item 11.

## Tests that prove it
- No library code: the measurements are the evidence. The harness's settings used are tested
  (CH4.4b's `check_m2l_override`, `Reference::extend`, `class_ratios`; CH4.3c2's `--rule`).
- Each base run's worst box and §9.1 ℓ2 checked against the gate's run on the same case and ε (the
  same build, the same reference sample: equal to the digits).
- The levels' parts add up (the rest printed, ≤ 0.04 of Σ|e|² as in CH4.4b; a larger rest
  reported).
- Pair evidence: `ch4-4b-tables.py pairs` (S7's worst pair over the 316 offsets) at each level's
  x and degree, as run and raised.

## Checks
`/verify` for what the branch changes (scripts and documents: `bash -n`, a dry run of each stage,
`py_compile` of any new Python, the roadmap through `tools/roadmap/sync_issues.py --out`); cargo's
checks if the harness changes (fmt, Clippy, `cargo test -p tetraleaf-validate`, rustdoc).

## Machine runs
Kathleen, `singlenode`, one rank × 40 threads, release, BLAS threads 1, node-c11b-004 and
node-c11b-013 excluded, every launch under the library's defaults (`--m2l defaults`) with the
reference refined (`--refine 8 --refine-faces --worst 8`) and the sample cached (the gate's cache
reused for W2 at kD = 250 and 500):
1. `map`: W2 with the constant density at kD = 250 and 500, ε = 10⁻³, 10⁻⁶, 10⁻⁸: `--study
   levels,m2l --all-levels --raise-by 1,2` (each level's part, each non-RCR level's M2L on RCR, the
   form's classes by decade, the level of the worst box's largest part raised by 1 and 2).
2. `all`: the same cases and W1 with the constant density at kD = 500: every level raised at once
   by 1 and 2 (`--study single --raise-levels v --raise-together --raise-by 1,2`).
3. `survey` (the third rule: the bound on the worst density over its whole range): W3 with the
   constant density, resolved at kD = 100 and at N = 10⁶ for kD = 200 and 500, the three ε, `--study
   levels,m2l --all-levels` (CH4.4a ran it at kD = 10 and 30 only). W1 and W2 with it at kD = 10–100
   are CH4.4a's and CH4.4b's (within ε, worst box ≤ 37.6 ε), at 250 and 500 the gate's.
4. `remedy` (after 1–3, its candidates from them): each candidate as a rule (`--rule …`, or the
   setting above), `--study single`, on every case that misses and on W1's constant density at
   kD = 250 and 500 (the three ε); then the cost: one launch per (candidate, case), the build and
   the evaluation of `--study single` (one setting per launch, CH4.3b2c §6), on W2 and W1 at kD = 500
   and W3 at N = 10⁶, kD = 500, beside the default.
Estimated: stages 1–3 about 2 node-hours (the W2 kD = 500 reference is cached; W3's about 10 min
each), stage 4 about 2–3. Jobs of at most 6 h. Raw outputs to
`~/projects/fmm-archive/kathleen/runs/ch4-8a/`, never the repository.

## Constraints and risks
- Constraint 1: no Laplace path (scripts and documents; at most the Helmholtz harness).
- Constraint 2: the degree rule changes only in CH4.8b, on both measures.
- The verification set widens (the H6 gate's first rule): CH4.8b's matrix is the whole of D1's
  kD range with the constant density on W1–W3 and the random charges of CH4.4a; this task runs
  the cases no record covered (W3's constant density beyond kD = 30), so that the remedy is
  chosen on every known miss.
- Risk: the remedy raises the plateau, which every low-frequency tree reads (W3, W4, the low kD):
  its cost then reaches every case; the cost runs measure it on W3 at N = 10⁶ too.
- Risk: no single degree fixes a tangent-point box whose error spreads over levels (CH4.4b: every
  level +1 took 114.8 ε to 48.4 ε); the remedy may be a raise on every level of a deep tree.

## Pending runs
None. Kathleen, `singlenode`, 1 × 40, `--exclude=node-c11b-004,node-c11b-013`, every one of 208
launches exit 0: 672424 (`smoke map all survey`, node-c11a-079), 672432 (`remedy cost`,
node-c11a-021), 672477 (`repeat`, node-c11a-021), 672488 (`mid-eps`, node-c11b-012; added at the
self-review, which found W2's ε between the gate's three unrun). The `remedy`, `cost` and `repeat` stages were added
to the script after `map` and `all` had shown the remedy's shape (commits cf85e11, f5114e1). Outputs
in `~/projects/fmm-archive/kathleen/runs/ch4-8a/`.

## Results
Report `tetraleaf-validate/results/helmholtz-ch4-8a.md`. The cause: the degree rule's truncation on
W2's deep tree, whose dense cap lies against the domain's z = 0 face, with the constant density adding
the levels' errors; at 10⁻³ the plateau's edge (k a = 3.38, Batched at p = 7) carries half of ∇φ's ℓ2,
at 10⁻⁶ the level at k a = 6.77, at 10⁻⁸ levels 2–4 the worst box; the M2L adds nothing (49 runs,
1.000 in the worst box). W3 with the constant density holds under the defaults (ℓ2 ≤ 0.51 ε, worst
box ≤ 17.3 ε). The remedy `h1` (the high-frequency branch one degree higher, the plateau unchanged)
puts W2 within ℓ2 0.48 ε and worst box 76.3 ε at all six ε of 10⁻³–10⁻⁸ (the default misses at
10⁻⁴ and 10⁻⁷ too), W1 within 0.28 ε and 17.6 ε, the worst pair on the
levels of k a ≥ 13.5 from 1.4–6.4 ε to 0.7–2.3 ε (on the levels that carry the ℓ2 misses 16 → 5.4 ε
and 32 → 20 ε: only the FMM shows those closed), at 0.989–1.047× the default's build plus one evaluation (A B B A,
geometric mean 1.017). Questions 1–3 of the report's §8 signed off as recommended (the user on #260, 2026-11-24: "All decisions as recommended.").

## Checks run
`/verify` on the M6: `cargo fmt --all -- --check`, `cargo clippy -- -D warnings`, `cargo clippy
--examples -- -D warnings`, `cargo doc --no-deps`, `cargo clippy --workspace --all-targets -- -D
warnings`, all clean; `cargo test` (1,392 passed) and `cargo test --workspace` (1,436 passed) at
opt-level 1, none failed. `bash -n` and ShellCheck on the runs script, `py_compile` on the tables
script, a dry run of every stage (208 launches) locally and on Kathleen, the roadmap through
`tools/roadmap/sync_issues.py --out`. Not run: the MPI list (no MPI crate's code changed), `cargo
deny` (no manifest change), Laplace's benchmarks (no Laplace path).

## Open questions
None. The report's §8 signed off as recommended (the user on #260, 2026-11-24: "All decisions as recommended."):
the acceptance as read; `h1` the remedy (CH4.8b: a margin of `DegreeRule` after the cap, §3.5
amended there); the plateau unchanged.
