# CH4.9b2 L2P with gradients by the kernel

Started 2026-11-28 on `tbetcke/h-ch4-9b2`. Roadmap: CH4.9b2, the second half of CH4.9b (item 14 of
the priority list; CH4.9 added on 2026-11-26, the user: "Yes. Add P2M and L2P to the priority list as
tasks"). No proposal §10 entry: the acceptance is the roadmap's *Done when* with constraint 1.
Proposal §2 (constraints 1, 6, 9); CONVENTIONS §3.15 ("Radial functions and scaling", "Storage");
CH4.9a2b's report `tetraleaf-validate/results/helmholtz-ch4-9a2b.md` (§5 the ISA question, §6 f32,
§7 question 1) and CH4.9b1's `tetraleaf-validate/results/helmholtz-ch4-9b1.md` (§3 the restated sum).

## Goal
Run the Helmholtz FMM's L2P, potentials and gradients, as one batched kernel per leaf on CH4.9a2a's
and CH4.9b1's routines, the per-point reference calls kept as the trusted path; measure L2P per level
and the evaluation before and after CH4.9a–b on Kathleen and locust, and answer the ISA question.

## Acceptance (the roadmap's *Done when*)
CH4.9a's *Done when* for L2P, potentials and gradients, on CH4.9b1's routines: the kernel in
`tetraleaf-tables::helmholtz::leaf` against `tetraleaf_ref::helmholtz::leaf::l2p` within a contract per
target in f64 and f32 (real and complex k, the centre's tangential limit), run by the operator's
`P2pChoice` with `reference_leaves` (`Reference` the per-point calls), m vectors, a leaf's bits the
same at every thread count and on P ranks the one-rank FMM's, the FMM's errors within 1.01× the
reference path's (W1's and W2's constant density included), L2P per level and the evaluation under
the defaults before and after CH4.9a–b on CH4.6a2's cases on Kathleen and locust; the ISA question
(the sum against its floor on both machines) and the f32 sum as above; constraint 1. Also (from
CH4.9b1, its report §3): the restated L2P sum with gradients is 37–62 % of a point on locust and
49–73 % on Kathleen, 1.12–1.57 ns per term in f64 on locust and 2.8–3.4 on Kathleen; the kernel's sum
is the part to vectorise, per target over terms, its three products on each plane. From CH4.9b's
parent line: the ISA question carried from CH4.9a2b (one ISA kernel for both sums if L2P's stays far
from its floor on x86), the f32 sum (CH4.9a2b §6), `reference_leaves` reused for L2P.

## Read first
- `tetraleaf-ref/src/helmholtz/leaf.rs` (`l2p`, `evaluate`, `Polar::{of, cartesian}`): the trusted
  path, its order of operations and the centre's limit k/(3s).
- `tetraleaf-math/src/orthonormal/leaf.rs` (`derivative_columns`, `columns`, `Lanes::of`),
  `radial.rs` (`bessel_j_real`, `bessel_j_real_derivatives`, `bessel_j`).
- `tetraleaf-tables/src/helmholtz/leaf.rs` (`LeafP2m`, `LeafCharges`): CH4.9a2b's kernel and layout.
- `tetraleaf-exec/src/helmholtz/densities.rs` (`l2p`, `p2m_kernel`, `LeafScratch`, the output chunk)
  and `operator.rs` (`Kernel::{p2m, l2p}`, `p2m_batched`, `reference_leaves`).
- `tetraleaf-exec/tests/helmholtz_fmm.rs` (`check_leaf_kernel`, `check_kernels`), `tests/helmholtz/
  densities.rs`; `tetraleaf-validate/examples/helmholtz_leaf.rs` (CH4.9b1's restated sum),
  `helmholtz_accuracy.rs` (`--leaves`).

## Changes
- `tetraleaf-math/src/orthonormal/leaf.rs`: `Lanes` becomes the public `LeafDirections<T, B>`
  (`of(points)`, fields cos θ, sin θ, cos φ, sin φ), the directions `columns` and
  `derivative_columns` use, so that the kernel's Cartesian gradient takes the harmonics' own frame
  (the rescale next to the z-axis included). No value changes (a rename and `pub`).
- `tetraleaf-tables/src/helmholtz/leaf.rs`:
  - `LeafOutput<'a, T> { values, vectors, gradients }`: a leaf's target output chunk as
    `tetraleaf-exec`'s Helmholtz densities hold it (two columns per vector, Re then Im; column o's φ̂
    of target i at o n + i, its ĝ at 2m n + 3 (o n + i) + axis), accumulated (+=).
  - `LeafL2p<T>`: `new(p)` holds a `LeafHarmonics` up to p; `l2p(p, level, targets, locals, scratch,
    output)` adds Σ L̂ₙᵐ ĵₙ(k|v|) e^{Im(k)|v| − β} Yₙᵐ(v̂) and, with gradients, its gradient for
    every vector, the locals complex (two planes of (p + 1)² per vector, the vectors one after the
    other). Per batch of 16 targets (the last padded with the centre, its outputs not written):
    r = √(x² + y² + z²); ĵₙ (and ĵₙ′ with gradients) by `bessel_j_real(_derivatives)` for real k,
    `bessel_j` per lane for complex k; e per lane (once for real k); the harmonics by `columns`
    (potentials only) or `derivative_columns`, and per order m and degree n the per-degree sums
    per lane and plane, from the coefficients of ±m as A = c⁺ cos mφ + c⁻ sin mφ and
    B = c⁻ cos mφ − c⁺ sin mφ: value += Q̄ A, ∂θ += (∂θ Q̄) A, (1/sin θ) ∂φ += m U B (m = 0: Q̄ c, ∂θ
    Q̄ c). Then per lane the reference's degree loop: φ = Σ ĵₙ valueₙ, radial = Σ ĵₙ′ valueₙ, and
    with tₙ = ĵₙ (1/r) (at r = 0: k/(3s) for n = 1, else 0) the polar and azimuthal sums, the
    Cartesian gradient by `LeafDirections` (k radial r̂ + polar θ̂ + azimuthal φ̂), all times e, added
    into the output. Lanes never meet: a target's values depend on it and the locals alone. No
    `unsafe`, no `mul_add`, no allocation (`L2pScratch::new(p, vectors)`). Yukawa refused.
  - The module docs: the L2P section (arithmetic, contract, bits).
- `tetraleaf-exec/src/helmholtz/densities.rs`: `Densities::l2p_kernel` (Yukawa keeps `l2p`);
  `LeafScratch` holds `LeafL2p` and its scratch.
- `tetraleaf-exec/src/helmholtz/operator.rs`: `Kernel::l2p` runs the kernel when the leaves are
  batched; `p2m_batched` becomes `leaves_batched` (P2M and L2P together; `reference_leaves` gives
  both references); the module docs' table rows and "P2M kernel" → "Leaf kernels".
- `tetraleaf-validate`: `Row::p2m_batched` → `Row::leaves_batched`, the example's `/ref-p2m` →
  `/ref-leaves`; `helmholtz_leaf`: the L2P kernel per point with and without gradients over the
  batch as one leaf (checked against the reference first, the contract's worst ratio printed), its
  ratio to the reference, its sum (the kernel less the routines' parts), ns per term and the sum's
  arithmetic rate (GFLOP/s by the kernel's operation count per term).
- `tools/helmholtz/ch4-9b2-runs.sh` (`smoke`, `levels`, `worst`, `parts`, `native`) and
  `ch4-9b2-tables.py` on `common.py`.
- Crate `CLAUDE.md`s (math, tables, exec, validate): a paragraph each.
- Report `tetraleaf-validate/results/helmholtz-ch4-9b2.md`.
- Conventions: none. No dependency. No Laplace path (`P2pChoice` read, not changed; `tetraleaf-simd`
  untouched; `LeafDirections` is Helmholtz's leaf module).

## Tests that prove it
- `tetraleaf-math/tests/orthonormal/leaf.rs`: `LeafDirections::of` on the poles, the centre and a
  rescaled direction (the values the harmonics use; existing contracts unchanged).
- `tetraleaf-tables/tests/helmholtz/l2p.rs`, against `leaf::l2p` per target:
  - the contract (f64), with w_n = √((2n + 1)/(4π)), Sₙ = Σₘ |L̂ₙᵐ|, x = k r, tₙ = |ĵₙ|/r (the limit
    at r = 0):
    |Δφ| ≤ C u e Σₙ wₙ Sₙ (((n + 1)² + 2n + 1) |ĵₙ| + |x| |ĵₙ′|),
    |Δgₐ| ≤ C u e Σₙ wₙ Sₙ (|k| (((n + 1)² + 2n + 3) |ĵₙ′| + |x| |ĵₙ|) + ((n + 1)³ + n(n + 1)) tₙ),
    C = 10, the worst ratio printed; p ∈ {0, 1, 2, 5, 8, 17, 26, 45}; targets N = 1, 5, 16, 17, 40
    (random in [−1, 1)³, the centre, the poles, the axes, two next to the z-axis); k̂ real (0.3, 4,
    16), complex (4 + 0.5i, 16 + 2i), tiny (10⁻⁶); with and without gradients; m = 1, 2, 3.
  - f32 against the f64 reference on the f32 inputs, u = 2⁻²⁴, p ≤ 32.
  - the centre: the tangential limit (the gradient at (0, 0, 0) within the contract).
  - bits: vector v of m the one-vector call's; a target's outputs the same alone, in any position
    and in a leaf of 1, 16 and 40; a second call adds exactly; the potentials with gradients the
    potentials without them bit for bit.
  - Yukawa refused; wrong lengths refused; the scratch's degree and vectors.
- `tetraleaf-exec/tests/helmholtz/densities.rs`: `the_batched_l2p_is_the_reference_to_rounding`
  (every path, with and without gradients, three vectors, within 1e-13 of the largest output,
  vector v alone bit for bit, accumulation, Yukawa the reference's bits).
- `tetraleaf-exec/tests/helmholtz_fmm.rs`: `check_leaf_kernel` now covers both kernels (the default
  against `reference_leaves`: within `LEAF_TOLERANCE` 1e-13, errors within `LEAF_RATIO` 1.01, the
  one-rank FMM); `check_kernels`' `Isa::Scalar` with `reference_leaves` bit for bit `Reference`;
  threads, ranks and vectors through every existing check on the default.

## Checks
`/verify` (root, strict workspace; `tetraleaf-math`'s, `tetraleaf-tables`', `tetraleaf-exec`'s and
`tetraleaf-validate`'s own; the MPI list at 2 and 4 ranks: `helmholtz_fmm` changes); the
`helmholtz_leaf` and `helmholtz_accuracy` smoke runs.

## Machine runs
- Kathleen (one node, `singlenode`, slow nodes excluded) and locust (gated on a quiet node):
  - `levels`: `helmholtz_accuracy --levels` at 1 × 40 and 1 × 72 on CH4.6a2's cases (W1 at kD = 250,
    500, 1000 at 10⁻⁶; W1 at kD = 500 at 10⁻³ and 10⁻⁸; W3 at N = 10⁶, kD = 500 at 10⁻⁶; complex
    charges, gradients), `--leaves kernel` (after CH4.9a–b) and `reference` (before), one setting
    per launch (B12), ABBA: L2P and P2M per level, the evaluation, the errors (the 1.01×).
  - `worst`: W1's and W2's constant density at kD = 500, 10⁻⁶ (CH4.3c3, CH4.8), errors.
  - `parts`: `helmholtz_leaf` on one core, p ∈ {8, 12, 17, 26, 32, 45} × k̂ ∈ {0.5, 2, 4, 8, 16},
    f64 and f32: the kernel per point against the reference, its sum's ns per term and GFLOP/s.
  - `native`: `parts` in f64 from a build with `-C target-cpu=native` (its own target directory):
    what wider vectors give the same safe Rust (AVX-512 on Kathleen, SVE2/NEON on locust), the
    ISA question's second measure beside the floor (the baseline's nominal peak: SSE2 two f64 lanes
    on two ports at 2.5 GHz on Kathleen; NEON on locust).
- The M6 not timed (ask first; not needed).

## Constraints and risks
- Constraint 1: no Laplace path (above). Constraint 6: the kernel against `leaf::l2p` per target,
  the FMM against the reference leaves and the direct sum. Constraint 9: the choice depends on the
  P2P choice and the setting alone.
- Risk: the per-degree accumulators (6 (p + 1) × 16 values per vector, 35 kB at p = 45) and the
  harmonics' three columns pass through L1/L2 per order; measured.
- Risk: the contract's form is new (the gradient's terms); it may need a term, as P2M's did.

## Pending runs
(filled at step 6)

## Open questions
(filled when they arise)
