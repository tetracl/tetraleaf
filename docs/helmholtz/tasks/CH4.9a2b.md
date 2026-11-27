# CH4.9a2b P2M by the kernel

Started 2026-11-27 on `tbetcke/h-ch4-9a2b`. Roadmap: CH4.9a2b, the second half of CH4.9a2 (item 14
of the priority list; CH4.9 added on 2026-11-26, the user: "Yes. Add P2M and L2P to the priority
list as tasks"). No proposal §10 entry: the acceptance is the roadmap's *Done when* with constraint
1. Proposal §2 (constraints 1, 6, 9); CONVENTIONS §3.15 ("Radial functions and scaling",
"Storage"); CH4.9a1's report `tetraleaf-validate/results/helmholtz-ch4-9a1.md` and CH4.9a2a's
`tetraleaf-validate/results/helmholtz-ch4-9a2a.md` (§4 f32, §5 lanes, §6 question 1).

## Goal
Run the Helmholtz FMM's P2M as one batched kernel per leaf on CH4.9a2a's routines, the per-point
reference calls kept as the trusted path, and measure P2M per level on Kathleen and locust.

## Acceptance (the roadmap's *Done when*)
CH4.9a's *Done when* less the table, on CH4.9a2a's routines (its question 1, the user on #268,
2026-11-27: "Q1 as recommended"; the FMM's errors within 1.01× here are their FMM evidence): the
P2M kernel (in `tetraleaf-tables::helmholtz`, against `tetraleaf_ref::helmholtz::leaf::p2m`, a
leaf's points batched, the sum vectorised over points without `unsafe`, real k on
`bessel_j_real`, complex k on `bessel_j`) run by the operator's `P2pChoice` (`Reference` the
per-point reference calls), real and complex k, real and complex charges, m vectors, within a
contract per leaf in f64 and f32, a leaf's bits the same at every thread count and on P ranks the
one-rank FMM's, the FMM's errors within 1.01× the reference path's, P2M's time per level on
Kathleen and locust (the M6 with the user's go-ahead); the ISA question in the PR if the time
stays far from the sum's; no `tetraleaf-simd` change; Yukawa's `p2m_yukawa` an Also for H2. Also:
the f32 slowdown, avoided or reported (subnormal products in the f32 sum, CH4.9a2a §4); 16 lanes;
the batched sum decides the ISA question.

## Read first
- `tetraleaf-ref/src/helmholtz/leaf.rs` (`p2m`, `form`, `offset`, `Polar::of`): the trusted path.
- `tetraleaf-math/src/orthonormal/leaf.rs` (`LeafHarmonics::columns`), `radial.rs`
  (`bessel_j_real`, `bessel_j`).
- `tetraleaf-exec/src/helmholtz/densities.rs` (`Densities::p2m`, `LeafScratch`, the chunk layout)
  and `operator.rs` (`Kernel::p2m`, `p2p: Option<P2pKernel<T>>`, `with_p2p`).
- `tetraleaf-exec/tests/helmholtz_fmm.rs` (`check_kernels`: `Isa::Scalar` bit for bit `Reference`).
- `tetraleaf-validate/examples/helmholtz_leaf.rs`, `helmholtz_accuracy.rs` (`--levels`, `--p2p`).

## Changes
As built (the self-review): the kernel sums each term over the 16 lanes per batch (a pairwise
reduction) and adds it to the coefficient, rather than per-lane accumulators summed once per leaf,
so that the scratch does not grow with the coarse levels' degrees (the operator sizes it by the
largest degree of every level); `P2mScratch::new(p, vectors)`. The accumulators per leaf are the
report's third alternative (§7). The plan as written follows.

- `tetraleaf-tables/src/helmholtz/leaf.rs` (new, `pub mod leaf` of `helmholtz`):
  - `LANES` = 16 (CH4.9a2a §5).
  - `LeafCharges<'a, T> { values, vectors, complex }`: a leaf's strengths point-major, as the
    chunks of §3.13 / CH3.0c hold them (point j's vector v at j c m + c v, Re then Im; c = 2 if
    complex).
  - `LeafP2m<T>`: `new(p)` holds a `LeafHarmonics` up to p; `p2m(p, level, points, charges,
    scratch, multipoles)` adds the stored multipoles of the leaf's points (offsets from the
    centre in the level's units) for every vector: `planes` (1 for real k with real charges, else
    2) times (p + 1)² values per vector, §3.15's storage, accumulated (+=). Per batch of 16
    points (the last padded with the centre and zero charges): r by `sqrt` of the squares, ĵₙ by
    `bessel_j_real` (real k) or `bessel_j` per lane (complex k), e = e^{Im(k) r − β} (once per
    leaf for real k), the weights wₙ = q e ĵₙ per vector and degree, then the harmonics column
    by column (`LeafHarmonics::columns`) with every term added into a per-lane accumulator
    (`[T; 16]` per coefficient and plane), the lanes summed once per leaf in a fixed pairwise
    order and added to the multipole. Every step a loop over the lanes without a branch on the
    data, no `unsafe`, no allocation (`P2mScratch::new(p, vectors, complex)` allocates once).
    Yukawa is refused (a panic: H2's path; `p2m_yukawa` an Also for H2).
  - The module docs: the arithmetic, the contract, the bits (a leaf's values depend on its
    points in order alone; vector v the one-vector kernel's bits).
- `tetraleaf-exec/src/helmholtz/densities.rs`: `Densities::p2m_kernel` (the kernel in the leaf's
  frame; the reference stays `Densities::p2m`); `LeafScratch` holds the kernel's table and
  scratch (built for the path in `LeafScratch::for_path`, or lazily sized; decided while
  building).
- `tetraleaf-exec/src/helmholtz/operator.rs`: `Kernel::p2m` runs the kernel unless the P2P choice
  is `Reference` or the operator's `reference_leaves` is set; `HelmholtzOperator::
  with_reference_leaves(bool)`, `reference_leaves()`; `HelmholtzFmmBuilder::reference_leaves`
  (a measurement and test setting: the leaf expansions by the per-point reference calls whatever
  the P2P kernel; default off). The module docs' table row for P2M. Yukawa never reaches it (the
  operator refuses Yukawa).
- `tetraleaf-exec/CLAUDE.md`, `tetraleaf-tables/CLAUDE.md`: a paragraph each.
- `tetraleaf-validate/examples/helmholtz_accuracy.rs`: `--leaves kernel|reference` (the builder's
  setting), shown in the run's settings.
- `tetraleaf-validate/examples/helmholtz_leaf.rs`: the kernel per point ("P2M kernel", µs per
  point over the batch as one leaf of `--points`), its check against the reference (the contract),
  and its ratio to the reference call and to the sum's part.
- `tools/helmholtz/ch4-9a2b-runs.sh` (`smoke`, `parts`, `levels`), `ch4-9a2b-tables.py` on
  `common.py`.
- Report `tetraleaf-validate/results/helmholtz-ch4-9a2b.md`.
- Conventions: none. No dependency (`tetraleaf-tables` already depends on `tetraleaf-math` and
  `tetraleaf-ref`). No Laplace path: `P2pChoice` is read, not changed; `tetraleaf-tables` gains a
  Helmholtz module; `tetraleaf-simd` untouched.

## Tests that prove it
- `tetraleaf-tables/tests/helmholtz/leaf.rs`:
  - the contract per term against `leaf::p2m` (f64): |M − M_ref| ≤ C u ((n + 1)² + N)
    √((2n + 1)/(4π)) Σⱼ |qⱼ eⱼ ĵₙ(k̂ρⱼ)| (+ the f32 flush term), C = 10, the worst ratio printed;
    p ∈ {0, 1, 2, 5, 8, 17, 26, 45}; leaves of N = 1, 5, 16, 17, 40, 200 points (random in
    [−1, 1)³, the centre, the poles, the axes); k̂ real (0.3, 4, 16), complex (4 + 0.5i, 16 + 2i),
    tiny (10⁻⁶); real and complex charges; m = 1, 2, 3.
  - f32 against the f64 reference on the f32 inputs with u = 2⁻²⁴ (p ≤ 32), and no subnormal
    coefficient.
  - vector v of m vectors bit for bit the one-vector kernel's; a second call adds (+=) exactly;
    real charges with real k the real plane of complex charges with Im q = 0 bit for bit.
  - Yukawa refused; wrong lengths refused.
- `tetraleaf-exec/tests/helmholtz/densities.rs`: `Densities::p2m_kernel` against `Densities::p2m`
  within the contract on every path (real and complex k and charges, m = 2).
- `tetraleaf-exec/tests/helmholtz_fmm.rs`:
  - every scenario now runs the kernel by default (`Auto`): 2 and 4 threads, pair by pair, two
    vectors and the one-rank FMM bit for bit hold for it (the existing checks);
  - `check_kernels`: `Isa::Scalar` with `reference_leaves` bit for bit `Reference` (the P2P
    kernel's identity, as before); every ISA within `KERNEL_TOLERANCE` of `Reference` (now with
    the P2M kernel too);
  - `check_leaf_kernel` (new; scenarios() and the cube at kD = 20 + 2i with real charges): the
    default against `reference_leaves(true)` on the same ISA: errors within `LEAF_RATIO` = 1.01×
    and the outputs within 1e-13 relative L2, the one-rank FMM with the reference leaves bit for
    bit on P ranks.

## Checks
`/verify` (root, strict workspace; `tetraleaf-tables`', `tetraleaf-exec`'s and
`tetraleaf-validate`'s own; the MPI list at 2 and 4 ranks: `helmholtz_fmm` changes); the
`helmholtz_leaf` and `helmholtz_accuracy` smoke runs.

## Machine runs
- Kathleen (one node, `singlenode`, slow nodes excluded) and locust (gated on a quiet node):
  - `levels`: `helmholtz_accuracy --levels` at 1 × 40 (Kathleen) and 1 × 72 (locust) on CH4.9's
    cases (W1 at kD = 250, 500, 1000 at 10⁻⁶; W1 at kD = 500 at 10⁻³ and 10⁻⁸; W3 at kD = 500 at
    10⁻⁶; N about 2·10⁶, complex charges), `--leaves kernel` and `reference`, one setting per
    launch (B12), ABBA: P2M per level, the evaluation, the errors (the 1.01× of the acceptance).
    The worst case the record names: W1's and W2's constant density at kD = 500, 10⁻⁶ (CH4.3c3,
    CH4.8), errors only.
  - `parts`: `helmholtz_leaf` on one core, p ∈ {8, 12, 17, 26, 32, 45} × k̂ ∈ {0.5, 2, 4, 8, 16},
    f64 and f32: the kernel per point against the reference and the sum's part (the ISA question,
    the f32 slowdown).
- The M6 not timed (ask first; not needed).

## Constraints and risks
- Constraint 1: no Laplace path (above). Constraint 6: the kernel against `leaf::p2m` per leaf,
  the FMM against the reference leaves and the direct sum. Constraint 9: the choice depends on the
  P2P choice and the setting alone.
- Risk: the lane accumulators ((p + 1)² × 16 × planes per vector, 540 kB at p = 45 with two
  planes) stream through L2; if the sum stays far from its floor, a per-batch reduction is the
  alternative (measured, then the ISA question).
- Risk: `Isa::Scalar` no longer `Reference`'s FMM bits; kept by `reference_leaves`.

## Pending runs
None. Kathleen job 673212 (node-c11a-080, 1 × 40; `smoke`, `levels`, `worst`, `parts`; 17 min 52 s,
2026-11-27 16:45–17:01) and locust (1 × 72, 16:45–16:57, after a gate: one-minute load below 2, no
other process of this account busy; a root `find` took 0.3–1 core throughout); every launch exit 0.
Outputs in `~/projects/fmm-archive/{kathleen,locust}/runs/ch4-9a2b/`. The M6 not timed (a short
smoke run of `helmholtz_leaf`, its numbers not used).

## Results (the report `tetraleaf-validate/results/helmholtz-ch4-9a2b.md`)
- Per leaf: within the contract C u √((2n + 1)/(4π)) Σⱼ |qⱼ eⱼ| (((n + 1)² + N) |ĵₙ| + |kρⱼ| |ĵₙ′|),
  C about 1.0 measured against the tests' 10 (as built: the first contract, without the conditioning
  term, failed at p = 2, k̂ = 16 near a zero of ĵ₀, where ρ by `sqrt` and by `hypot` differ in the
  last bit).
- FMM: errors 0.993–1.000× the reference P2M's (three printed digits) (W1's and W2's constant density at kD = 500 included);
  tests 1.0000×, outputs within 3.0·10⁻¹⁶.
- P2M per run 0.25–0.31× (Kathleen), 0.15–0.27× (locust); 0.6–2.4 % of an evaluation after; the
  evaluation 0.92–0.99× and 0.94–0.97×.
- Per point: the kernel 0.20–0.24× and 0.12–0.20× the reference (f64); its batched sum no faster than
  the scalar sum on Kathleen (baseline SSE2, two lanes, a lane reduction per term).
- f32: the kernel 0.72× its f64 time at p = 45 on Kathleen, but its batched sum shows the slowdown's
  pattern (0.88× the f64 sum, 1.45× the scalar sum); not diagnosed.

## Checks run (the M6)
`/verify`: fmt; Clippy default, examples, `--workspace --all-targets` (one `is_multiple_of` lint in the
new test, fixed); `cargo doc` with `-D warnings`; at opt-level 1 `cargo test` (1,421 passed) and
`cargo test --workspace` (1,465 passed); the MPI list at 2 and 4 ranks (every executable passed;
`helmholtz_fmm` 53.7 s and 49.4 s); `helmholtz_fmm` at one rank, opt-level 1, 53 s; the one-rank
debug runs at opt-level 0, load 3–6: `helmholtz_fmm` 951 s (`main` 934 s in a worktree of its own,
+1.8 %), `tetraleaf-exec`'s `helmholtz` 129 s, `tetraleaf-tables`' `helmholtz` 44 s. Not run:
`cargo deny` (no manifest change), Laplace's benchmarks (no Laplace path).

## Open questions
For sign-off on the PR (the report §7):
1. The ISA question: no ISA kernel for P2M now, the question carried to CH4.9b (recommended); or an
   ISA kernel of P2M's sum in `tetraleaf-simd` now; or lane accumulators per leaf in safe Rust.
