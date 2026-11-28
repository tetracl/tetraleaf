# CH2.1b The Yukawa operator and FMM

Started 2026-11-28 on `tbetcke/h-ch2-1b`. Roadmap: CH2.1b, the second part of CH2.1 (item 15 of the
priority list, Yukawa on the host; added on 2026-11-24 at the H4 gate). Proposal: §10 CH2.1; §2
(constraints 1, 2, 6, 9); §5.5 (Yukawa: real coefficients, real tables, no high-frequency regime);
§6.2 (the Yukawa row); §6.3 (the operator, D4); §9.1 (W1, W3, W4 and the norms); §12 D4.
CONVENTIONS §3.15 ("Yukawa's real form", "Inputs and outputs"). Builds on CH2.1a
(`docs/helmholtz/tasks/CH2.1a.md`: the one-plane RCR and Dense tables).

## Goal
`HelmholtzOperator` and `HelmholtzFmmBuilder` accept Yukawa (k = iκ) on CH2.1a's tables, the FMM is
measured against the direct sum on W1, W3 and W4 for κD from 10⁻² to 10³, and D4 goes to the user
with the complex path's time before and after.

## Acceptance (proposal §10 CH2.1, verbatim, for the parent)
"Real `HelmholtzOperator` path on the drivers of CH2.0 (or on copies); Rotation and Dense M2L" —
"FMM against the direct sum on W1, W3 and W4 for κ D from 10⁻² to 10³; deterministic".

CH2.1b's part (the roadmap's *Done when*): the FMM against the direct sum on W1, W3 and W4 for κD
from 10⁻² to 10³ (the errors reported; within target is CH2.3's); deterministic (threads, chunks,
the one-rank FMM bit for bit); real and complex k's output bit for bit as before; D4 put for
sign-off with the complex path's time before and after. Also (from CH2.1a): `BatchedM2l` and
`ClassM2l` refuse Yukawa, so `Auto` on Yukawa chooses among Rotation and Dense unless this task
adds their one-plane tables (say which).

## Read first
- `tetraleaf-exec/src/helmholtz/operator.rs`: `unbuilt` (the refusal), `select_tables`, the bodies
  `p2m` … `p2p`, `each_target`; `operator_rcr.rs` (the contribution's planes and vectors);
  `wave_levels`; `operator_diagonal.rs::select`.
- `tetraleaf-exec/src/helmholtz/densities.rs`: Yukawa's real path, and its complex charges as 2m
  real vectors value for value (the same chunks, expansions and output columns).
- `tetraleaf-exec/src/helmholtz/fmm.rs`: `check_settings` (the refusal), `build`.
- `tetraleaf-tables/src/helmholtz.rs` and `dense.rs` after CH2.1a: Yukawa's tables take one plane;
  `DenseM2l::select_for(…, local_planes)`.
- `tetraleaf-ref/src/helmholtz/leaf.rs` (`p2l_yukawa`, `m2p_yukawa`) and `p2p.rs` (`p2p_yukawa`);
  `tetraleaf-simd`'s `P2pKernel::evaluate_yukawa`.
- `tetraleaf-validate/src/helmholtz_accuracy.rs` (`Problem`, `errors`: the imaginary plane read as 0
  where the output has none), `helmholtz_workloads.rs` (W1's N is a function of kD).

## Changes
**`tetraleaf-exec`** (Helmholtz code only; no Laplace path, `driver.rs` and `pipeline.rs` unchanged):

- Operator, Yukawa accepted (`OperatorError::Yukawa` removed: no input gives it). On Yukawa every
  real charge column is one real vector of one plane: the tables (M2M, L2L, RCR's M2L and its
  contributions by pairs, Dense) take `planes` 1 and `columns()` vectors (m for real charges, 2m for
  complex: Re q⁰, Im q⁰, …, the same storage as Densities' two planes per vector). Two private
  helpers, `table_planes` (multipole and local) and `table_vectors`, give `(densities' planes, m)`
  otherwise, so the real and complex paths call the tables with the arguments they pass today.
- Bodies on Yukawa: P2M and L2P `Densities`' (the reference's `_yukawa` forms; `p2m_kernel` and
  `l2p_kernel` already keep them); P2L `leaf::p2l_yukawa` per column into its plane (the charges
  times r̂(t|s)); M2P `leaf::m2p_yukawa` per column into output column j (the column's potentials and
  gradient triples are contiguous: CONVENTIONS §3.13's layout with the output columns as vectors);
  P2P `P2pKernel::evaluate_yukawa(κ̂, …)` per column, the reference `p2p::p2p_yukawa` under
  `P2pChoice::Reference`. `PlaneWaveRaise::M2l`'s first-degree path as for the others.
- M2L on Yukawa: `Rotation` and `Dense` (`DenseM2l::select_for` with one local plane). `Auto` runs
  Rotation on every level: no Batched or Classes tables (CH2.1a refuses them), no exponential M2L
  (its tables are for Re k > 0; `wave_levels` gives none for Yukawa) and no diagonal form on any
  level whatever `HelmholtzDiagonal` says (proposal §5.5: the high-frequency regime does not occur).
  `with_m2l` panics on Yukawa under `Batched`, `Classes` and `PlaneWave` (documented).
  Their one-plane tables are left to a follow-up (below), so Auto's choice on Yukawa is Rotation's.
- Builder: Yukawa accepted; `Batched`, `Classes` and `PlaneWave` named on Yukawa refused before any
  collective (`HelmholtzFmmError::YukawaM2l`); the rest of `build` unchanged (min_leaf_level 0,
  the degree rule at x = Re(k) a = 0, i.e. its plateau, which CH2.3 calibrates).
- Module docs and `tetraleaf-exec/CLAUDE.md`: the Yukawa path.

**`tetraleaf-validate`**: `helmholtz_workloads::Workload::W1n { n }` (W1's sphere at n points: W1's N
follows the wavelength, which Yukawa has not), the example's `w1:<n>`; the harness runs Yukawa
through `--kd 0 --kim κ`; `tools/helmholtz/ch2-1b-runs.sh` (stages below) and `ch2-1b-tables.py` on
`common.py`.

No dependency, no convention, no `tetraleaf-plan` or Laplace change.

## Tests that prove it
1. `tests/helmholtz/operator.rs` (no MPI): the chains (far, v4, x, w, near) on Yukawa at κD = 0.01,
   1, 20 and 100, real and complex charges, ε = 10⁻⁶ and 10⁻¹⁰ at S7's degrees, against the direct
   sum within ε (near within 1e-13); vector v of three each the one-vector operator bit for bit;
   complex charges (q_r, q_i) bit for bit the real path with the two real vectors q_r and q_i; Dense's
   chain within 1e-14 of Rotation's; the refusals of `with_m2l` under Batched, Classes, PlaneWave;
   the `low_degree_chains_are_the_reference_chains` comparison with `tetraleaf_ref`'s `_yukawa`
   chains (leaf operators and `rcr`) within its 2e-13.
2. `tests/helmholtz/p2p.rs`: every ISA against `Reference` on Yukawa within CH2.4b's contract summed.
3. `tests/helmholtz_fmm.rs` (`check_yukawa`, MPI): Yukawa on the cube and the clusters at κD = 0.01,
   1 and 30 with a tolerance: against `direct_sum_yukawa` within ε; 2 threads, pair calls and on P
   ranks the one-rank FMM bit for bit; two vectors each the one-vector FMM; complex charges equal to
   the real FMM of the 2m real vectors bit for bit; Dense within `DENSE_TOLERANCE` of Rotation; the
   default builder (Auto) Rotation's bits; every ISA within `KERNEL_TOLERANCE` of `Reference`;
   `check_errors` the refused M2L.
4. Real and complex k unchanged: every existing test of `tetraleaf-exec` and `tetraleaf-validate`
   unchanged and passing; the harness's checksums on CH3.1c-style cases (W1 at kD = 30 and 100, W3
   at kD = 30, W4 at kD = 10, real and complex k, Auto and Rotation) equal on `main` and the branch
   (the Kathleen job below).
5. `tetraleaf-validate`: `W1n` point for point W1's construction (`w1_count`'s N gives W1's points),
   a Yukawa case in `tests/helmholtz_accuracy.rs` (checksum at 1 and 2 threads, errors below ε).

## Checks
`/verify` (fmt; Clippy default, all-targets, workspace; `cargo test`, `--workspace`; docs; the MPI
list at 2 and 4 ranks; `cargo deny` only if a manifest changes).

## Machine runs
Kathleen, `singlenode`, 1 × 40, `--exclude=node-c11b-004,node-c11b-013`, release, BLAS threads 1:
- `accuracy` (decides the acceptance's report): W1 (`w1:100000`), W3 (N = 10⁵ and 10⁶), W4's four
  distributions (N = 10⁵) at κD = 10⁻², 10⁻¹, 1, 10, 10², 10³, ε = 10⁻³, 10⁻⁶ and 10⁻⁸ (W3 at
  10⁶ at 10⁻⁶ only), real charges, the library's defaults (Auto: Rotation), gradients; at 1 × 40
  and, for determinism, 1 × 1 and 4 × 10 on a subset (W1 and W4 plummer at κD = 1 and 100, 10⁻⁶);
  Dense against Rotation on W4 cube at κD = 1, 10⁻⁶; complex charges on W1 at κD = 1.
- `checksums` (decides "bit for bit as before"): item 4's cases, `main` and the branch.
- `time` (decides D4's evidence): the complex path, `main` against the branch, one setting per
  launch, ABBA twice: W1 at kD = 100 and 250 (real k), W3 at N = 10⁶ at kD = 50 + 5i (complex k), at
  10⁻⁶, the defaults, two evaluations; build and evaluation.
`main`'s binary from `src/main` on Kathleen (its own target directory).

## Constraints and risks
- Constraint 1: no Laplace code touched (`driver.rs`, `pipeline.rs`, `operator.rs` of Laplace,
  `tetraleaf-simd`, `tetraleaf-plan`, `tetraleaf-octree` unchanged).
- Constraint 2: Yukawa's level calls run on the same drivers and splits, so threads, chunks and
  ranks keep their bit identities; tested.
- Constraint 9: Auto on Yukawa reads only the kind.
- The complex path: the Yukawa branches are per leaf body or per level call (a test of the kind),
  never in an inner loop; its time is measured (D4), its output checked by checksum.
- Low κ on deep trees: RCR's power factors overflow at p log₁₀(1/s) ≳ 300 (review 2026-11-19, M13;
  CH2.3 refuses it). At κD = 10⁻² and N = 10⁶ the deepest level's s is about 10⁻⁵ and p ≤ 30, so
  p log₁₀(1/s) ≲ 150; a failure there is reported, not hidden.
- The degree rule's plateau is calibrated for Re k > 0; Yukawa's errors may sit far below ε at
  large κ (strong damping) or near it at small κ (the Laplace limit). Reported; CH2.3 calibrates.

## Pending runs
None. Kathleen, `singlenode`, `--exclude=node-c11b-004,node-c11b-013`, all on node-c11b-046,
2026-11-28: job 673276 (`main`'s release build at 2413c5b, `src/main`), 673277 (`smoke accuracy
determinism dense complex checksums` on e0f648b; superseded by 673280 after the self-review's fix of
the V rule's prices), 673278 (`time` on e0f648b), 673279 (the `-align-loops=64` builds of both trees,
then `time time-aligned` with `LEVELS=1 EVALS=3`, four cases), 673280 (`smoke accuracy determinism
dense complex checksums` on 503ad58, with the real-charge checksum case). Every launch exited 0.
Outputs in `~/projects/fmm-archive/kathleen/runs/ch2-1b/` (`acc`, `acc2`, `time`, `time2`).

## Results (the report `tetraleaf-validate/results/helmholtz-ch2-1b.md`)
- **The acceptance:** within ε in all 114 cases (W1 at 10⁵ points, W3 at 10⁵ and 10⁶, W4's four
  distributions at 10⁵; κD = 10⁻² to 10³; ε = 10⁻³, 10⁻⁶, 10⁻⁸), at most 0.27 ε (W4 Plummer, 10⁻⁸).
  The degrees are the rule's plateau at every κ (7, 15, 24); at κD = 1000 on the uniform and sphere
  trees the errors are the P2P's round-off.
- **Deterministic:** checksums equal at 1 × 1, 1 × 40 and 4 × 10 on four cases; the tests' threads,
  pair calls, RCR by pairs and the one-rank FMM at 2 and 4 ranks.
- **Real and complex k:** `main`'s and the branch's checksums equal in 12 of 12 pairs (six cases,
  real charges included, defaults and Rotation), and in every timing launch.
- **D4's time:** evaluation branch/main 0.991–1.040 (geometric mean 1.015) on default builds,
  0.947–1.035 (1.003) on aligned builds, the ranges overlapping in every case and the signs changing
  between the builds: no cost measured (job 673278's first, lighter set gave 1.036).
- **Tests (the M6):** the Yukawa chains within ε at κD = 0.01–100 but the far chain at κD = 100
  (3.0e-7 at every ε: the expansion's cancellation e^{20} u, held to 3e-6), our chains `rcr`'s
  `_yukawa` chains to 3e-14 even at κ = 100, Dense 1e-16–1.6e-15 from Rotation, complex charges two
  real vectors bit for bit, P2P on NEON within 0.009 of CH3.4b's bound; the FMM's Yukawa scenarios
  1e-5–1e-7 against bounds 1e-2–1e-5.

## Self-review
A subagent (plan, proposal §2 and §10, the branch) found no blocking defect. Fixed after it: the V
rule priced Yukawa's Auto levels at Batched's cost (now Rotation's; the accuracy stages rerun as job
673280, errors lower by up to 13 %); the far chain's exception explained by the expansion's own
cancellation and held to 3e-6 (was 1e-5); `check_rcr_split` on a Yukawa scenario; a real-charge real-k
checksum case; doc fixes. Left as noted: no f32 Yukawa FMM test (the plan names none; the P2P test
covers f32); Yukawa's P2P forms each pair's distance once per real column (c m calls, where the
complex path shares it over vectors), a performance note for CH2.1d; the `PlaneWaveRaise::M2l` branches
of Yukawa's P2L and M2P are unreachable today (no waves on Yukawa) and kept so the bodies follow the
others.

## Open questions
1. **D4.** The proposal recommends one operator `HelmholtzOperator<T, K>`, K a type parameter, "if
   CH2.1 shows no cost to the complex path". Built: one operator `HelmholtzOperator<T>` with the kind
   a runtime field resolved at build (as real and complex k since CH3.1a), read per level call or leaf
   body, never in an inner loop; no cost measured (report §5). Recommended: sign off one operator with
   the kind resolved at build as a field, and amend §5.5, §6.3 and §12 to say so. Alternatives: the
   type parameter K as the proposal sketches it (a refactor of the operator, no measured gain to
   expect); Yukawa its own operator type (copies of the operator's bodies and settings).
2. **CH2.1d's priority.** Recommended: on the priority list in item 15 after CH2.3 (its rules read the
   calibrated degrees), before CH2.5. Alternative: off the list until H2.G or later (Yukawa's FMM then
   runs Rotation on every level, correct and slower at low degrees).
