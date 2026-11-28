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
  W1 and W2, constant density, κD = 1 and 100, 10⁻⁶, kernel and reference (errors).
- `levels` (P2M's time, reported, not a criterion): W1 at 10⁶ points and W3 at 10⁶, κD = 1 and 100,
  10⁻⁶, `--levels`, ABBA, one setting per launch (B12).
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
(filled in step 6)

## Open questions
(filled when they arise)
