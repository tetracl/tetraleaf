# CH2.1c Yukawa's batched P2M

Started 2026-11-28 on `tbetcke/h-ch2-1c`. Roadmap: CH2.1c, the third part of CH2.1 (item 15 of the
priority list; split on 2026-11-28): "CH2.1's Also: `LeafP2m` for Yukawa as it restates `p2m`,
`p2m_yukawa` restated, `Densities::p2m_kernel` on it. *Done when:* its contract as CH4.9a2b's, the
FMM's errors within 1.01× the per-point path's. *Needs:* CH2.1b." Proposal: §10 CH2.1 (the FMM
against the direct sum on W1, W3 and W4 for κD from 10⁻² to 10³; deterministic); §2 (constraints 1,
6, 9); §5.5 (Yukawa, real arithmetic); CONVENTIONS §3.15 ("Radial functions and scaling",
"Storage": M̂ᵞ = i⁻ⁿ M̂, (p + 1)² reals).

## Goal
Run the Yukawa FMM's P2M as CH4.9a2b's batched kernel per leaf (`LeafP2m`), restating
`tetraleaf_ref::helmholtz::leaf::p2m_yukawa`, with the per-point reference calls kept as the trusted
path.

## Acceptance (the roadmap's *Done when*; proposal §10 CH2.1)
- The contract as CH4.9a2b's: every coefficient of a leaf's Yukawa multipole within
  C u √((2n + 1)/(4π)) Σⱼ |qⱼ eⱼ| (((n + 1)² + N) |îₙ(xⱼ)| + |xⱼ| (|îₙ′(xⱼ)| + |îₙ(xⱼ)|)) of
  `p2m_yukawa`'s, xⱼ = κρⱼ, C = 10, in f64 and in f32 (against the f64 reference on the f32 inputs).
  The extra |xⱼ| |îₙ| is the factor e^{κρ − β}'s condition (ρ by `sqrt` against the reference's
  `hypot`), which for real k is one value per leaf and for complex k sits inside ĵ's scaling.
- The FMM's errors within 1.01× the per-point path's (`reference_leaves`), on the tests' Yukawa
  scenarios and on CH2.1b's cases on Kathleen; the FMM still within ε there (CH2.1's acceptance).
- Deterministic: a leaf's bits a function of its points and charges in order; the FMM bit for bit
  at every thread count and on P ranks the one-rank FMM's (the existing checks, now on the kernel).

## Read first
- `tetraleaf-ref/src/helmholtz/leaf.rs`: `p2m_yukawa`, `form`, `offset`, `Polar::of` (the trusted
  path: per point `modified_i(κr, s, …)` of length p + 1, `harmonics`, e = e^{κr − β}, then
  (q e) îₙ times Yₙᵐ).
- `tetraleaf-tables/src/helmholtz/leaf.rs` (`LeafP2m::p2m`, `planes`, `P2mScratch`) and its tests
  `tetraleaf-tables/tests/helmholtz/leaf.rs`.
- `tetraleaf-math/src/radial.rs` (`modified_i`).
- `tetraleaf-exec/src/helmholtz/densities.rs` (`Densities::p2m`, `p2m_kernel`'s Yukawa fall-back) and
  `operator.rs` (`leaves_batched`, `Kernel::p2m`).
- `tetraleaf-exec/tests/helmholtz_fmm.rs` (`check_leaf_kernel`, `yukawa_scenarios`, `check_kernels`).
- `tetraleaf-validate/results/helmholtz-ch4-9a2b.md`, `helmholtz-ch2-1b.md`.

## Changes
- `tetraleaf-tables/src/helmholtz/leaf.rs`:
  - `LeafP2m::planes`: 1 for real charges with real k or Yukawa (k = iκ), 2 otherwise; for Yukawa
    with complex charges plane 0 is Re q's real multipole and plane 1 Im q's (the layout of
    `Densities`: a complex vector on Yukawa is two real vectors).
  - `LeafP2m::p2m` on Yukawa: per lane îₙ(κρ) by `modified_i` (the reference's routine, of length
    p + 1, so a lane's values are the reference's at its ρ) into a real buffer, e = e^{κρ − β} per
    lane, the weights wₙ = (q e) îₙ per plane as the real-k branch forms them, then the same
    columns and lane sums. The refusal goes.
  - `P2mScratch`: one real buffer of p + 1 for a lane's îₙ.
  - The module docs: Yukawa's formula, its contract term, the bits (a complex vector's planes the
    real-charge kernel on Re q and on Im q bit for bit).
- `tetraleaf-exec/src/helmholtz/densities.rs`: `p2m_kernel` runs the kernel on Yukawa (the fall-back
  goes); the module docs and `p2m_kernel`'s docs. `operator.rs`: the module docs' P2M row and
  `leaves_batched`'s doc (Yukawa's L2P stays the reference's).
- `tetraleaf-tables/CLAUDE.md`, `tetraleaf-exec/CLAUDE.md`: the lines for the change.
- `tools/helmholtz/ch2-1c-runs.sh`, `ch2-1c-tables.py` (on `common.py`).
- Report `tetraleaf-validate/results/helmholtz-ch2-1c.md`.
- Conventions: none. No dependency. No Laplace path: `tetraleaf-tables`' Helmholtz `leaf` module and
  `tetraleaf-exec`'s Helmholtz densities only.

## Tests that prove it
- `tetraleaf-tables/tests/helmholtz/leaf.rs` (or a sibling `leaf_yukawa.rs`):
  - the contract per coefficient against `p2m_yukawa` in f64: p ∈ {0, 1, 2, 5, 8, 17, 26, 45},
    κ̂ ∈ {10⁻⁶, 0.3, 4, 16, 60} (κ̂ a to about 100: strong damping, where îₙ e is large and the
    multipole's values span many orders), leaves of N = 1, 5, 16, 17, 40, 200, real and complex
    charges, m = 1, 2, 3; the worst ratio printed;
  - f32 against the f64 reference on the f32 inputs (p ≤ 32), the subnormal check as CH4.9a2b's (or
    the flush it needs, decided by the test);
  - bits: vector v the one-vector call's; a complex vector's two planes the real-charge kernel on
    Re q and on Im q; accumulation; the scratch's size;
  - the refusal test replaced by these; the wrong-length refusals on Yukawa.
- `tetraleaf-exec/tests/helmholtz/densities.rs`: `p2m_kernel` against `p2m` within the contract on
  Yukawa, real and complex charges, m = 2.
- `tetraleaf-exec/tests/helmholtz_fmm.rs`: `check_leaf_kernel` on the Yukawa scenarios (errors within
  `LEAF_RATIO` = 1.01× of `reference_leaves(true)`, outputs within `LEAF_TOLERANCE`, the reference
  leaves' one-rank FMM bit for bit on P ranks); the existing Yukawa checks (threads, pairs, vectors,
  one rank) now run the kernel.

## Checks
`/verify` (root, strict workspace, `tetraleaf-tables`', `tetraleaf-exec`'s and
`tetraleaf-validate`'s own; the MPI list at 2 and 4 ranks: `helmholtz_fmm` changes).

## Machine runs
Kathleen, `singlenode`, 1 × 40, slow nodes excluded, release, BLAS threads 1:
- `accuracy` (the 1.01× and ε on the FMM): CH2.1b's cases (W1 at 10⁵, W3 at 10⁵ and 10⁶, W4's four
  distributions at 10⁵; κD = 10⁻² … 10³; ε = 10⁻³, 10⁻⁶, 10⁻⁸), `--leaves kernel` and
  `--leaves reference`, one launch each.
- `worst` (the worst density the record names: the constant density, CH4.3c3, CH4.8; CH2.3's Also):
  W1 and W3 at 10⁶ points, constant density, κD = 1 and 100, 10⁻⁶, kernel and reference (errors).
  As built: W1's and W2's defaults size their points by Re k and hold none on Yukawa; W1 runs as
  W1n, and W2, which has no fixed-N form, is replaced by W3.
- `levels` (P2M's time, reported, not a criterion): W1 at 2·10⁶ points and W3 at 10⁶, κD = 1 and
  100, 10⁻⁶, `--levels`, ABBA, one setting per launch (B12).
The M6 is not timed.

## Constraints and risks
- Constraint 1: no Laplace path touched.
- Constraint 6: the kernel against `p2m_yukawa` per leaf; the FMM against the reference leaves and the
  direct sum.
- Constraint 9: the choice is the operator's existing one (`P2pChoice`, `reference_leaves`).
- Risk: `modified_i` per lane keeps the radial functions scalar; Yukawa's gain may be smaller than
  real k's (CH4.9a2b's real k used `bessel_j_real`). A batched `modified_i_real` in `tetraleaf-math`
  would be a follow-up, not this task.
- Risk: f32 îₙ at small κρ and large n is subnormal (iₙ ~ xⁿ/(2n + 1)!!), so the f32 sum may slow
  and the subnormal check may fail; the kernel then flushes îₙ below 2⁻¹⁰⁰ as `bessel_j_real` does
  (a change below u, stated in the contract). f32 is not an FMM path for Yukawa before H7 (D3).
- Yukawa's L2P stays the reference's (`l2p_kernel` falls back): a follow-up line in the roadmap.

## Pending runs
None. Kathleen, `singlenode`, slow nodes excluded, both on node-c11a-011, 2026-11-28: job 673284
(build, `smoke`, `accuracy`, and `worst` and `levels` on the default W1 and W2, which held no point
on Yukawa; 4 min 47 s) and job 673285 (`worst` and `levels` on W1n and W3, 2 min 9 s, after the
script's fix 4893d27). Every launch exited 0. Outputs in `~/projects/fmm-archive/kathleen/runs/ch2-1c/`
(`acc`, `time`, and the jobs' logs).

## Results (the report `tetraleaf-validate/results/helmholtz-ch2-1c.md`)
- The contract per leaf: 0.80 of C = 10 in f64 (p ≤ 45, κ̂ ≤ 250), 1.54 in f32 against the f64
  reference beside the f32 level's own rounding (2 (β + n) u, held apart) and the flushes' floor.
- The FMM: 114 of 114 of CH2.1b's cases within ε (at most 0.27 ε), the batched P2M's errors 1.0000×
  the reference P2M's (three printed digits) in every case, with complex charges and on the constant
  density; the tests 1.0000×, outputs within 4.1·10⁻¹⁶.
- P2M per run 0.29–0.34× (Kathleen 1 × 40), from 7–10 % of the evaluation to 2–3 %; the evaluation
  0.93–0.95×. Yukawa's reference L2P is 15–25 % of an evaluation after it.
- A finding: the constant density on W1 at 10⁶, κD = 100, 10⁻⁶ misses ε on ∇φ (2.6 ε) under either
  P2M: an Also on CH2.3's line.

## Self-review
A subagent (plan, proposal §2 and §10, the branch) found no blocking defect in the code. Fixed after
it: the f32 level's rounding held apart from C with its own constant 2 (inside the C bracket it could
hide a kernel error at low degree), the flushes' floor with the harmonics' norm, κ̂ = 250 and three
complex vectors in the f64 contract, a Yukawa charge-length refusal, a kind check in `p2m_kernel`
(beside the planes' check), the module docs' claim on subnormal operands, the operator's and
`tetraleaf-exec/CLAUDE.md`'s wording, and the plan's runs as built. Left as noted: f64's flush
(2⁻⁹⁰⁰) is not reached by the test points.

## Open questions
For sign-off on the PR (the report §5), signed off 2026-11-28, the user on #275: "Q1 and Q2 as
recommended":
1. CH2.1e (Yukawa's batched L2P) on the priority list in item 15 right after CH2.1c (recommended), off
   the list, or not at all.
2. The constant density's gradient on Yukawa recorded on CH2.3's line (recommended; done on the
   branch), or a task of its own before CH2.3.
