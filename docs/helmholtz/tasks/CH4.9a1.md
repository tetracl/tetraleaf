# CH4.9a1 Where a point's P2M and L2P time goes, and the per-level table

Started 2026-11-26 on `tbetcke/h-ch4-9a1`. Roadmap: CH4.9a1, the first half of CH4.9a (item 14 of the priority list, added on
2026-11-26 after CH4.6a2, the user: "Yes. Add P2M and L2P to the priority list as tasks"). No
proposal §10 entry: the acceptance is the roadmap's *Done when* with constraint 1. Proposal §2
(constraints 1, 2, 9), §6.2 (the drivers), §3.15 via CONVENTIONS (the leaf operators).

## Size: CH4.9a split

The *Done when* holds three pieces of work: the per-level kind table in the harness; a new kernel
on four ISAs in `tetraleaf-simd` with the operator's dispatch and its tests; and its times on two
machines. That is a new operator plus measurements on several machines, the skill's rule for a
split. And the evidence so far says the reference's time is mostly per-point setup, not the
O(p²) sum that SIMD across points would speed up most directly: on Kathleen (job 673083, 1 × 40)
the finest level of W1 at kD = 1000 (p = 17, 324 terms) took about 3.3 µs of thread time per point
for P2M (about 10 ns per term) and 7.9 µs for L2P with gradients (about 24 ns per term), where a
plain sum costs about 1 ns per term. Per point the reference runs a complex Miller recursion
(`tetraleaf_math::radial::bessel_j`: a start order that depends on the point, Smith's division
per step), the orthonormal harmonics (an exponent per column) and an `exp`. Which of these
dominates decides the kernel's shape, and whether it needs `tetraleaf-simd` (a crate on
Laplace's path, constraint 1, with `unsafe`) at all. So, measuring first (as CH4.1a and CH4.6a2
did):

- **CH4.9a1 Where a point's P2M and L2P time goes, and the per-level table.** *Done when:*
  `helmholtz_accuracy --levels` prints P2M, M2M, L2L and L2P per level (the table of job 673083,
  with its documentation and a test); an MPI-free example times one point's P2M and L2P (with and
  without gradients) by part (the radial functions, the harmonics and their derivatives, the
  sum, the rest) on one core, at the leaf levels' degrees and k ρ of W1–W3 under the defaults
  (p = 8–45, k ρ to about 22; f64 and f32), against the reference call, on Kathleen and locust;
  a question for sign-off: the kernel's design (a restructured kernel in safe Rust in
  `tetraleaf-exec` or `tetraleaf-tables`, an ISA kernel in `tetraleaf-simd`, or both) and the
  order of CH4.9a2 and CH4.9b. *Needs:* CH3.4c. *Machines:* Kathleen, locust.
- **CH4.9a2 P2M by the kernel chosen.** *Done when:* the roadmap's CH4.9a *Done when* less the
  table: the kernel as signed off at CH4.9a1, run by the operator's `P2pChoice` (`Reference` the
  per-point calls), within a contract per leaf against `tetraleaf_ref::helmholtz::leaf::p2m` in
  f64 and f32, a leaf's bits the same at every thread count and on P ranks the one-rank FMM's,
  the FMM's errors within 1.01× the reference path's, P2M's time per level on Kathleen and
  locust. *Needs:* CH4.9a1 and its sign-off.

CH4.9b (L2P) and CH4.9c (Laplace) keep their lines; CH4.9b needs CH4.9a2.

Signed off 2026-11-26, the user: "Split a1/a2 (Recommended)".

The rest of this plan is **CH4.9a1**.

## Goal (CH4.9a1)
Measure where a Helmholtz leaf point's P2M and L2P time goes, and put the per-level kind times in
the harness, so that the kernel's design rests on numbers.

## Acceptance (CH4.9a1, as split above)
See the split.

## Read first
- `tetraleaf-ref/src/helmholtz/leaf.rs` (`p2m`, `l2p`, `form`, `for_each_target`), its
  `Workspace`; `tetraleaf-math/src/radial.rs` (`bessel_j`, `start_order_j`, `j0_j1`, `div`) and
  `orthonormal/harmonics.rs` (`harmonics`, `harmonics_derivatives`).
- `tetraleaf-exec/src/helmholtz/densities.rs` (`Densities::{p2m, l2p}`, `LeafScratch`): the calls
  per point, vector and plane.
- `tetraleaf-validate/examples/helmholtz_conversions.rs` (CH4.1a: an MPI-free timing example of
  this kind) and `tetraleaf-validate/src/bench.rs` (its timing helpers).

## Changes
- `tetraleaf-validate/examples/helmholtz_accuracy.rs`: the table "Translations and leaf
  expansions" with `--levels` (from `Row::kind_seconds`, already collected), its module
  documentation; `tests/helmholtz_accuracy.rs` checks its source, `Row::kind_seconds`, on a small
  case: P2M and L2P only on levels with leaves, M2M and L2L zero on levels 0–1's inputs as
  recorded.
- A new MPI-free example `tetraleaf-validate/examples/helmholtz_leaf.rs` (not registered; it
  times): per (p, k ρ) on one thread, a batch of points in a leaf frame, the reference P2M and L2P
  per point, then their parts timed apart by calling the same public functions the reference
  calls (`bessel_j`, `harmonics`, `harmonics_derivatives`, the sum as a loop), each result checked
  against the reference's first; f64 and f32. A report script `tools/helmholtz/ch4-9a1-tables.py`
  on `common.py`; a runs script `tools/helmholtz/ch4-9a1-runs.sh`.
- `tetraleaf-validate/CLAUDE.md`: a paragraph.
- No library code, no Laplace path, no dependency.

## Tests that prove it
- The harness table: `tests/helmholtz_accuracy.rs` (the kinds per level present where the plan has
  that work).
- The example checks each part's result against the reference call before timing it (relative
  difference at most 10⁻¹² in f64, bit for bit where it calls the same function).

## Checks
`/verify`; the example's smoke run.

## Machine runs
Kathleen (one core, `singlenode`, slow nodes excluded) and locust (one core, gated on a quiet
node): the example at p ∈ {8, 12, 15, 17, 19, 26, 32, 41, 45} and k ρ ∈ {0.5, 2, 5, 10, 22}, f64
and f32; one setting per launch where a time decides. The M6 not run (the user, 2026-11-26, for
CH4.6a2's timing: "Skip the M6"; asked again for this task only if needed).

## Constraints and risks
- Constraint 1: no Laplace path touched in CH4.9a1. CH4.9a2's design may touch `tetraleaf-simd`;
  that is its sign-off.
- Risk: the parts timed apart may not add up to the whole (call overhead, cache); the example
  prints the sum of the parts against the whole.

## Pending runs
None. Kathleen job 673085 (node-c11a-142, one core, 20:45–20:51) and locust (one core, 20:45–20:52,
started by the quiet gate) ran `smoke` and `parts`; every launch exited 0, every case's parts
reassembled the reference bit for bit (90 per machine). Outputs in
`~/projects/fmm-archive/{kathleen,locust}/runs/ch4-9a1/`.

## Results (the report `tetraleaf-validate/results/helmholtz-ch4-9a1.md`)
From p = 12 on the harmonics take 44–81 % of a point's P2M (f64) and 58–84 % of L2P with gradients;
the radial functions most of the rest; the sum 4–10 %. f32 barely faster (subnormal slowdowns on
x86 at p = 41–45). The per-point cost on one core equals the in-FMM thread time per point.
Recommendation: CH4.9a2 as restructured safe Rust (harmonics from coefficients tabulated per degree,
real radial recursion for real k, points batched), no `tetraleaf-simd` change.

## Checks run (the M6)
`/verify`: fmt; Clippy default, examples, `--workspace --all-targets`; `cargo doc` with `-D
warnings`; at opt-level 1 `cargo test` (1,402 passed) and `cargo test --workspace` (1,446 passed);
`tests/helmholtz_accuracy.rs` at opt-level 1 (passed, 17.8 s; the strengthened `check_leaves`); the
tables script's `selftest`; the example's smoke run (bit for bit). Not run: the MPI list (no change
to tetraleaf-exec, -plan or -octree); `cargo deny` (no manifest change).

## Open questions
Asked before implementing, answered 2026-11-26: the split, "Split a1/a2 (Recommended)".

For sign-off on the PR (the report §4):
1. CH4.9a2's kernel: restructured safe Rust first (recommended), or ISA kernels in `tetraleaf-simd`,
   or both.
2. The order: the leaf harmonics routine built with P2M in CH4.9a2 and reused by CH4.9b
   (recommended), or one task for both.
