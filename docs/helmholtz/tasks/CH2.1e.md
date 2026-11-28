# CH2.1e Yukawa's batched L2P

Started 2026-11-28 on `tbetcke/h-ch2-1e`. Roadmap: CH2.1e (item 15 of the priority list, after
CH2.1c; added at CH2.1c, #275's question 1, the user, 2026-11-28: "Q1 and Q2 as recommended"):
"`LeafL2p` on Yukawa's real path as CH2.1c's `LeafP2m`: `l2p_yukawa` restated with gradients (îₙ and
îₙ′ per lane by `modified_i`, e^{κr − β} per lane, the tangential factor's limit at the centre),
`Densities::l2p_kernel` on it, the real output columns of Yukawa's real path. *Done when:* its
contract as CH4.9b2's, a target's bits independent of its leaf, batch and position, the FMM's errors
within 1.01× the per-point path's on CH2.1b's cases, L2P's time per level on Kathleen. *Needs:*
CH2.1c." Proposal: §10 CH2.1 (the FMM against the direct sum on W1, W3 and W4 for κD from 10⁻² to
10³; deterministic); §2 (constraints 1, 6, 9); §5.5 (Yukawa, real arithmetic); CONVENTIONS §3.15
("Radial functions and scaling", "Storage": L̂ᵞ = iⁿ L̂, (p + 1)² reals).

## Goal
Run the Yukawa FMM's L2P, potentials and gradients, as CH4.9b2's batched kernel per leaf
(`LeafL2p`), restating `tetraleaf_ref::helmholtz::leaf::l2p_yukawa`, with the per-target reference
calls kept as the trusted path.

## Acceptance (the roadmap's *Done when*; proposal §10 CH2.1)
- The contract as CH4.9b2's, per target in f64 and in f32 (against the f64 reference on the f32
  inputs), with ĵ → î, k → κ and the factor's condition added (e^{κr − β} per lane, r by `sqrt`
  against the reference's `hypot`, as CH2.1c's P2M contract added it): with wₙ = √((2n + 1)/(4π)),
  Sₙ = Σₘ |L̂ᵞₙᵐ|, x = κr, tₙ = îₙ(x)/r (κ/(3s) for n = 1 at r = 0, else 0), e = e^{κr − β},
  |Δφ| ≤ C u e Σₙ wₙ Sₙ (((n + 1)² + 2n + 1) îₙ + x (|îₙ′| + îₙ)),
  |Δgₐ| ≤ C u e Σₙ wₙ Sₙ (κ (((n + 1)² + 2n + 3) |îₙ′| + x îₙ + x |îₙ′|) + ((n + 1)³ + n(n + 1) + x) tₙ),
  C = 10 (the extra x terms are e's condition on each output term).
- A target's bits independent of its leaf, batch and position; vector v the one-vector call's; a
  complex vector's two columns the real-charge kernel on Re q's and on Im q's local bit for bit;
  potentials with gradients those without; accumulation.
- The FMM's errors within 1.01× the per-point path's (`reference_leaves`) on the tests' Yukawa
  scenarios and on CH2.1b's 114 cases on Kathleen; the FMM still within ε there (CH2.1's
  acceptance); deterministic (threads, ranks, vectors: the existing checks, now on the kernel).
- L2P's time per level on Kathleen, kernel against the reference (reported, not a criterion).

## Read first
- `tetraleaf-ref/src/helmholtz/leaf.rs`: `l2p_yukawa`, `evaluate` (F = T), `Polar`.
- `tetraleaf-tables/src/helmholtz/leaf.rs`: `LeafL2p::l2p`, `LeafOutput`, `L2pScratch`, and CH2.1c's
  Yukawa branch of `LeafP2m::p2m` (`modified_i` per lane, the flush below τ).
- `tetraleaf-math/src/radial.rs` (`modified_i` with derivatives: î′₁(0) = 1/(3s)).
- `tetraleaf-tables/tests/helmholtz/l2p.rs` (CH4.9b2's contract and bit tests),
  `leaf_yukawa.rs` (CH2.1c's κ̂ set, `LEVEL_C`).
- `tetraleaf-exec/src/helmholtz/densities.rs` (`l2p`, `l2p_kernel`'s Yukawa fall-back),
  `operator.rs` (`leaves_batched`'s docs), `tests/helmholtz/densities.rs`
  (`the_batched_l2p_is_the_reference_to_rounding`), `tests/helmholtz_fmm.rs` (`check_leaf_kernel`
  on `yukawa_scenarios`).
- `tools/helmholtz/ch2-1c-runs.sh`, `ch2-1c-tables.py`; report `helmholtz-ch2-1c.md`.

## Changes
- `tetraleaf-tables/src/helmholtz/leaf.rs`:
  - `LeafOutput` gains `planes`: the output columns per vector, 2 (Re, Im; every path but Yukawa's
    real one) or 1 (Yukawa with real charges). `per_target` is planes m (× 4 with gradients); the
    gradients start at planes m n. Its callers name `planes: 2` (exec, the tests, `helmholtz_leaf`).
  - `LeafL2p::planes(level, complex)`: 1 on Yukawa with real charges, else 2, as `LeafP2m::planes`.
  - `LeafL2p::l2p` on Yukawa (k = iκ): the locals `planes` real (p + 1)² planes per vector (1 for
    real charges; Re q's and Im q's locals for complex ones, `Densities`' layout), each into its own
    output column. Per batch: îₙ (and îₙ′ with gradients) by `modified_i` per lane at the reference's
    length p + 1, values below τ flushed as CH2.1c's P2M flushes them; e = e^{κr − β} per lane,
    flushed below τ; the same per-degree sums per plane; then per lane and plane the reference's
    real degree loop (φ = Σ îₙ valueₙ, radial Σ îₙ′ valueₙ, tₙ = îₙ (1/r) or κ/(3s) at the centre for
    n = 1), the gradient κ radial r̂ + polar θ̂ + azimuthal φ̂, all times e. Real and complex k keep
    their arithmetic unchanged (their bits as `main`'s). The refusal goes; a non-Yukawa level with
    `planes` 1 is refused.
  - The module docs: "L2P" gains Yukawa's formula, contract and bits.
- `tetraleaf-exec/src/helmholtz/densities.rs`: `l2p_kernel` runs the kernel on Yukawa (the fall-back
  goes, a kind check as `p2m_kernel`'s), `LeafOutput { planes: self.local_planes() }`; the module and
  method docs. `operator.rs`: the module docs' L2P row and `leaves_batched`'s doc.
- `tetraleaf-validate/examples/helmholtz_leaf.rs`: `planes: 2` in its `LeafOutput` (no other change).
- `tetraleaf-tables/CLAUDE.md`, `tetraleaf-exec/CLAUDE.md`: a line each.
- `tools/helmholtz/ch2-1e-runs.sh` (CH2.1c's stages: `smoke`, `accuracy`, `worst`, `levels`) and
  `ch2-1e-tables.py` (on `common.py`; L2P's time per level).
- Report `tetraleaf-validate/results/helmholtz-ch2-1e.md`.
- Conventions: none. No dependency. No Laplace path: `tetraleaf-tables`' Helmholtz `leaf` module,
  `tetraleaf-exec`'s Helmholtz densities and docs, one literal in a validate example.

## Tests that prove it
- `tetraleaf-tables/tests/helmholtz/l2p_yukawa.rs`, against `leaf::l2p_yukawa` per target:
  - the contract (f64): p ∈ {0, 1, 2, 5, 8, 17, 26, 45}; κ̂ ∈ {10⁻⁶, 0.3, 4, 16, 60} and 250 (CH2.1c's
    set: κ̂ a to 433); targets N = 1, 5, 16, 17, 40 (the centre, poles, axes, next to the z-axis);
    real and complex charges' locals (1 and 2 planes), m = 1, 2, 3; with and without gradients; the
    worst ratio printed;
  - f32 against the f64 reference on the f32 inputs (p ≤ 32, κ̂ ≤ 60), with CH2.1c's level term
    (`LEVEL_C` (β + n) u, held apart) and the flushes' floor;
  - targets whose r² underflows, and the centre (the tangential limit κ/(3s));
  - bits: vector v the one-vector call's; a complex vector's two columns the one-plane kernel on each
    plane; a target alone, in any position, in leaves of 1, 16 and 40; potentials with gradients
    those without; accumulation;
  - refusals: a non-Yukawa level with one plane; Yukawa's wrong lengths. `l2p.rs`'s
    `yukawa_is_refused` replaced.
- `tetraleaf-tables/tests/helmholtz/l2p.rs`: unchanged in substance (`planes: 2`); real and complex k
  bit for bit `main`'s (the existing bit tests; the arithmetic of those paths is untouched).
- `tetraleaf-exec/tests/helmholtz/densities.rs`: `the_batched_l2p_is_the_reference_to_rounding` on
  Yukawa within 1e-13 like the other paths (the "Yukawa keeps the reference" asserts go); Yukawa's
  complex vectors two real ones under the kernel (L2P too).
- `tetraleaf-exec/tests/helmholtz_fmm.rs`: `check_leaf_kernel` on Yukawa's three scenarios now covers
  L2P (`LEAF_TOLERANCE`, `LEAF_RATIO` 1.01, the one-rank FMM); the existing Yukawa checks (threads,
  pairs, vectors, one rank) run the kernel.

## Checks
`/verify` (root, strict workspace, `tetraleaf-tables`', `tetraleaf-exec`'s and
`tetraleaf-validate`'s own; the MPI list at 2 and 4 ranks: `helmholtz_fmm` changes).

## Machine runs
Kathleen, `singlenode`, 1 × 40, slow nodes excluded, release, BLAS threads 1 (CH2.1c's stages):
- `accuracy` (the 1.01× and ε): CH2.1b's 114 cases (W1n at 10⁵, W3 at 10⁵ and 10⁶, W4's four
  distributions at 10⁵; κD = 10⁻² … 10³; ε = 10⁻³, 10⁻⁶, 10⁻⁸), `--leaves kernel` and
  `--leaves reference`, one launch each; W1n at κD = 1 and 100 with complex charges.
- `worst` (the worst density the record names: the constant density, CH4.3c3, CH4.8; CH2.3's Also):
  W1n and W3 at 10⁶ points, κD = 1 and 100, 10⁻⁶, kernel and reference.
- `levels` (L2P's time per level, reported): W1n at 2·10⁶ and W3 at 10⁶, κD = 1 and 100, 10⁻⁶,
  `--levels`, ABBA, one setting per launch (B12).
The M6 is not timed.

## Constraints and risks
- Constraint 1: no Laplace path touched; real and complex k's L2P arithmetic unchanged.
- Constraint 6: the kernel against `l2p_yukawa` per target; the FMM against the reference leaves and
  the direct sum. Constraint 9: the choice is the operator's existing one.
- Risk: `modified_i` per lane stays scalar (as CH2.1c's P2M); the gain may be below real k's.
- Risk: the gradient's contract form is new for Yukawa; it may need a term, as P2M's did.
- Risk: f32 îₙ′ subnormal near the centre at large n; flushed below τ as îₙ.

## Pending runs
None. Kathleen job 673289 (`run.sbatch`, `singlenode`, `--exclude=node-c11b-004,node-c11b-013`) on
node-c11b-046, 2026-11-28 14:56–15:02, after 1 s in the queue: the release build of f5656e8, then
`tools/helmholtz/ch2-1e-runs.sh … smoke accuracy worst` and `… levels`; all 114 launches exit 0.
Outputs in `~/projects/fmm-archive/kathleen/runs/ch2-1e/` (`acc`, `time`, the job's log).

## Results (the report `tetraleaf-validate/results/helmholtz-ch2-1e.md`)
- The contract per target: 1.94 of C = 10 in f64 (p ≤ 45, κ̂ ≤ 250), 0.64 in f32 against the f64
  reference beside the terms held apart; 0.22 where r² underflows, 0.18 next to the centre in f32.
- The FMM: 114 of 114 of CH2.1b's cases within ε (at most 0.267 ε), the batched leaves' errors
  1.0000× the reference leaves' in every case, with complex charges and on the constant density;
  the tests 1.0000×, outputs within 4.1·10⁻¹⁶.
- L2P per level 0.21–0.28× (Kathleen 1 × 40), from 14–22 % of an evaluation to 3.7–6.5 %; the
  evaluation with both leaf kernels 0.77–0.89×.

## Self-review
A subagent (plan, proposal §2 and §10, the branch) found no blocking defect. Fixed after it: the
stale doc comments of `SUMS` and `L2pScratch`, the flush's absolute bound in the module docs (it
holds for the FMM's O(1) stored locals), an f32 case next to the centre (r² subnormal or zero), the
hash comparison and the runs recorded here. Left as noted: the position and leaf bit tests run one
plane and one vector (each column is the one-plane call's bits, checked apart); f32 runs one vector.

## Bits of real and complex k
Real and complex k's L2P is `main`'s bit for bit: a scratch test (not committed) hashed `LeafL2p`'s
outputs over k̂ ∈ {0.3, 4, 16, 4 + 0.5i, 16 + 2i, 10⁻⁶}, p ∈ {0, 3, 12, 26}, 1 and 3 vectors, with
and without gradients, 37 targets, on the branch and on `main`'s `leaf.rs`: `ee41f6d8571cd78e` both
(the M6, 2026-11-28).

## Open questions
None yet.
