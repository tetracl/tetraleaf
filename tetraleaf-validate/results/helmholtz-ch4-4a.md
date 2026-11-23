# Helmholtz CH4.4a: the kD sweep beside FMM3D

Task plan and record: `docs/helmholtz/tasks/CH4.4a.md`. Proposal: §10 CH4.4 (CH4.4a takes all of
it but "mapped to levels and classes, with no error concentrated at a class boundary", which is
CH4.4b's), §9.1, §9.2, §7 "Build cost". Machine: Kathleen, one node per job (2 × 20 Cascade Lake
cores, Xeon Gold 6248, 188 GB), `singlenode`, release, BLAS threads 1; ours at 4 ranks × 10 threads
bound by socket (`--map-by ppr:2:socket:PE=10`, CH3.6's split), FMM3D v2.1.0 on 40 OpenMP threads.
Raw outputs: `~/projects/fmm-archive/kathleen/runs/ch4-4a/` (the cases' Markdown and logs, the jobs'
output, and `slow-nodes/`, §6).

## Summary

- **Accuracy.** Under the library's defaults (`Case::defaults`: the diagonal form's time rule, the
  refit degrees of CH4.3c3, c_leaf 4, the V, W and X rules) our FMM is within ε at §9.1's 1,000
  samples on every case and ε of the sweep: W1 at kD = 10–1000, W2 at kD = 10–500, W3 resolved at
  kD = 10–100 and at N = 10⁶ for kD = 200 and 500, the plane wave's trace on W1 at kD = 500 and
  1000, the constant density on W1 and W2 at kD = 10–100 and W3 at kD = 10, 30; 26 cases × 3 ε, the
  largest error 0.76 ε (∇φ, W2 with the constant density at kD = 100, 10⁻⁶).
- **The full-target crosscheck** against FMM3D at 10⁻¹² (its own error at the samples at most
  3.9 · 10⁻¹³): our all-target errors within ε everywhere; **our worst box at most 100 ε on every
  case but one**: W1 at kD = 1000 with the plane wave's trace at 10⁻⁸, **114.7 ε in ∇φ** (69.5 ε in
  φ), the box at the sphere's tangent point on the cube's +x face (CH4.0b's placement). Next largest:
  73.9 ε (the same case at 10⁻⁶), 67.8 ε (W1, plane wave, kD = 500, 10⁻⁸), 43.3 ε (W1, random
  charges, kD = 1000, 10⁻⁸). A finding, for sign-off (question 2); CH4.4b maps it to levels and
  classes.
- **Time per k.** Ours (build plus one evaluation, the slowest rank) is 0.09–0.95× FMM3D's faster
  variant's whole call on the 42 case × ε points with N ≥ 10⁵ (geometric mean 0.21), 0.08–0.17× on
  W1 at kD = 500 and 1000; our evaluation alone 0.03–0.62×. Below N = 10⁵ (geometric mean 0.63) ours
  is slower on six points, all at 10⁻³, the first ε of a launch, whose build carries about 0.1 s of
  one-time process costs (§4). These are times at the same requested ε, not at matched accuracy (CH4.5).
- **Memory per k.** Our FMM's resident growth on the node: 1.0–1.3 GB on W1 at kD = 250, 5.6–7.0 GB
  at 500, 21–29 GB at kD = 1000 (8 · 10⁶ points). FMM3D's driver peaks at 9.7–27 GB (`fastker`)
  and 44–108 GB (`default`) on W1 at kD = 500; at kD = 1000 `default` exceeds the job's 180 GB at
  10⁻³, 10⁻⁶ and 10⁻¹² (killed) and peaks at 191 GB at 10⁻⁸, `fastker` 45–118 GB. Our figures after a
  launch's first ε are lower bounds (§5).
- **Build per k.** W1: 0.09–0.13 s at kD = 100, 0.46–0.64 at 250, 2.6–2.9 at 500, 16–19 s at 1000,
  where the RCR tables are 9.2–12.9 s of it (67 % of build plus one evaluation at 10⁻³).
- **Two slow Kathleen nodes.** node-c11b-004 and node-c11b-013 ran both codes about 3× slower on
  problems of 10⁶ points and more, with the same outputs; their cases were run again elsewhere and
  only those runs are reported (§6).

## 1. Method

- **Cases.** `helmholtz_fmm3d problem` writes each case (D = 1, real k, one complex charge vector,
  §9.1's 1,000 sampled targets, the density named) to FMM3D's problem file; `compare` refuses a file
  that is not its case's byte for byte. kD = 10–1000 on W1 (N ≈ 7.96 (kD)²), W2 (W1 with the dense
  cap), W3 resolved (⌈1000 (kD/2π)³⌉ points) and W3 at N = 10⁶. "W3 to 200" is read as §9.1 states
  W3: resolved to kD = 100 (4.0 · 10⁶ points; 3.2 · 10⁷ at 200) and under-resolved at N = 10⁶ for
  kD = 200 and 500. The densities: §9.1's random charges on every case; the plane wave's trace
  (`PLANE_WAVE_DIRECTION`) on W1 at kD = 500 and 1000 (CH4.0b's sign-off 2); the constant density
  on the low-frequency trees, kD ≤ 100 (CH4.3c3's report §4).
- **Ours.** `helmholtz_fmm3d compare --settings defaults`: `Case::defaults`, every setting the
  library's (`HelmholtzFmmBuilder::new`; `tests/helmholtz_accuracy.rs`'s `check_defaults` holds it
  to the builder bit for bit), at ε = 10⁻³, 10⁻⁶, 10⁻⁸ in one launch per case, each ε built and
  evaluated once; the AVX-512 P2P (`Auto`).
- **FMM3D.** `hfmm3d_s_c_g` through `tools/fmm3d/driver.f90`, one launch per eps (under
  `tools/fmm3d/maxrss.py` for its peak memory, CH3.6 ran one launch for every eps), the variants
  `default` and `fastker` at ε (one variant on the density cases), FULL = 1. The reference at
  10⁻¹²: `default`, except at kD = 1000, where `default` runs out of the node's memory and `fastker`
  is the reference (plan, "Open questions").
- **Errors.** §9.1's relative ℓ2 of φ and ∇φ at the samples against the f64 direct sum
  (`tetraleaf_ref::helmholtz::p2p::direct_sum` × 1/(4π)), both codes. The crosscheck: every output
  against the reference at every point, the all-target relative ℓ2 and §9.1's box measure (per box
  of our tree's deepest level, the largest error at a target over the reference's RMS), gated at
  100 ε for ours.
- **Times.** Single launches; ours the build and the evaluation, each the slowest rank's; FMM3D's
  the driver's call, its setup included (it keeps nothing between calls).
- **Memory.** Ours (`Row::memory`): per rank the peak resident size `VmHWM` over the build and the
  evaluation after resetting it (`5` to `/proc/self/clear_refs`) less the resident size then,
  summed over the 4 ranks (and the largest rank's). FMM3D's: the driver process's peak resident set
  (getrusage), its input and output arrays (13 doubles a point, 0.83 GB at kD = 1000) included.
- **Jobs.** 663706 (W1, W2 to kD = 250, W3 to 30, the constant density; node-c11b-043), 663732
  (W2 at kD = 500; node-c11b-043), 663733 (W1 at kD = 1000; node-c11a-035), 663734 (the plane wave
  at kD = 1000; node-c11a-143), 663736 (W1 at kD = 500, both densities, and the three larger W3
  cases; node-c11a-119), 663737 (memory and ε order, §5; node-c11a-070), 663735 (the slow nodes'
  check, §6; node-c11a-143); 663705 the smoke run. Superseded: 663707, 663710 (node-c11b-004) and
  663708 (node-c11b-013), slow; 663708 and 663709 also lost FMM3D `default` at kD = 1000 to memory.

## 2. Accuracy

"ours ÷ ε" is the larger of φ and ∇φ at the samples over ε; the worst box is over ε, φ / ∇φ;
FMM3D's at its own eps = ε, by variant ("–": the run left no output, §4). On the smallest trees
our FMM sums directly (W1 and W2 at kD = 10: two levels, no far field, with c_leaf 4; W1 at kD = 30,
W3 at kD = 10 and 30 at 10⁻⁸), so its errors there are round-off.

| Case | ε | ours φ, ∇φ | ours ÷ ε | ours all φ, ∇φ | ours worst box ÷ ε | FMM3D φ, ∇φ (by variant) | FMM3D worst box ÷ ε |
| --- | --- | --- | --- | --- | --- | --- | --- |
| W1, kD = 10 | 1e-3 | 1.44e-15, 1.74e-15 | 0.00, 0.00 | 1.79e-15, 2.10e-15 | 0.00 / 0.00 | default 1.54e-5, 3.95e-5; fastker 9.03e-16, 1.06e-15 | default 0.09 / 0.20; fastker 0.00 / 0.00 |
|  | 1e-6 | 1.44e-15, 1.74e-15 | 0.00, 0.00 | 1.79e-15, 2.10e-15 | 0.00 / 0.00 | default 1.21e-15, 1.30e-15; fastker 1.29e-15, 1.39e-15 | default 0.00 / 0.00; fastker 0.00 / 0.00 |
|  | 1e-8 | 1.44e-15, 1.74e-15 | 0.00, 0.00 | 1.79e-15, 2.10e-15 | 0.00 / 0.00 | default 1.21e-15, 1.30e-15; fastker 1.29e-15, 1.39e-15 | default 0.00 / 0.00; fastker 0.00 / 0.00 |
| W1, kD = 30 | 1e-3 | 5.07e-5, 8.05e-5 | 0.05, 0.08 | 5.27e-5, 9.64e-5 | 0.58 / 1.60 | default 7.50e-5, 7.42e-5; fastker 7.48e-5, 7.37e-5 | default 0.35 / 0.36; fastker 0.35 / 0.35 |
|  | 1e-6 | 1.08e-7, 2.50e-7 | 0.11, 0.25 | 1.39e-7, 3.26e-7 | 3.34 / 8.26 | default 1.05e-8, 1.19e-8; fastker 2.29e-15, 3.06e-15 | default 0.08 / 0.16; fastker 0.00 / 0.00 |
|  | 1e-8 | 8.83e-13, 2.13e-12 | 0.00, 0.00 | 8.30e-13, 1.65e-12 | 0.00 / 0.00 | default 5.67e-12, 7.75e-12; fastker 2.29e-15, 3.06e-15 | default 0.00 / 0.01; fastker 0.00 / 0.00 |
| W1, kD = 100 | 1e-3 | 1.54e-4, 2.55e-4 | 0.15, 0.26 | 1.68e-4, 2.42e-4 | 5.43 / 7.20 | default 1.23e-3, 1.20e-3; fastker 1.23e-3, 1.20e-3 | default 6.73 / 5.60; fastker 6.70 / 5.57 |
|  | 1e-6 | 1.39e-7, 2.99e-7 | 0.14, 0.30 | 1.27e-7, 2.92e-7 | 4.51 / 14.72 | default 3.25e-7, 4.40e-7; fastker 3.23e-7, 4.36e-7 | default 4.31 / 6.08; fastker 4.30 / 6.16 |
|  | 1e-8 | 6.28e-10, 8.74e-10 | 0.06, 0.09 | 7.64e-10, 1.35e-9 | 5.04 / 8.67 | default 1.33e-10, 2.09e-10; fastker 1.33e-10, 2.07e-10 | default 0.21 / 0.41; fastker 0.21 / 0.41 |
| W1, kD = 250 | 1e-3 | 2.28e-4, 3.12e-4 | 0.23, 0.31 | 2.35e-4, 3.22e-4 | 13.23 / 13.57 | default 2.10e-4, 1.96e-4; fastker 2.09e-4, 1.95e-4 | default 1.39 / 1.26; fastker 1.42 / 1.24 |
|  | 1e-6 | 1.50e-7, 2.32e-7 | 0.15, 0.23 | 1.35e-7, 2.57e-7 | 8.85 / 16.98 | default 6.82e-8, 7.81e-8; fastker 6.60e-8, 7.33e-8 | default 1.07 / 1.37; fastker 1.13 / 1.37 |
|  | 1e-8 | 1.31e-9, 1.52e-9 | 0.13, 0.15 | 1.41e-9, 1.94e-9 | 16.12 / 21.55 | default 5.64e-11, 6.29e-11; fastker 5.63e-11, 6.26e-11 | default 0.13 / 0.16; fastker 0.13 / 0.16 |
| W1, kD = 500 | 1e-3 | 2.85e-4, 3.97e-4 | 0.28, 0.40 | 2.62e-4, 3.68e-4 | 21.28 / 22.99 | default 1.87e-4, 1.75e-4; fastker 1.86e-4, 1.75e-4 | default 1.72 / 1.55; fastker 1.72 / 1.53 |
|  | 1e-6 | 1.51e-7, 2.96e-7 | 0.15, 0.30 | 1.38e-7, 2.78e-7 | 14.22 / 20.01 | default 6.21e-8, 7.19e-8; fastker 6.09e-8, 6.95e-8 | default 1.47 / 1.45; fastker 1.42 / 1.35 |
|  | 1e-8 | 1.16e-9, 1.48e-9 | 0.12, 0.15 | 1.34e-9, 1.75e-9 | 27.38 / 40.42 | default 5.19e-11, 6.33e-11; fastker 5.17e-11, 6.31e-11 | default 0.17 / 0.21; fastker 0.17 / 0.21 |
| W1, kD = 1000 | 1e-3 | 2.68e-4, 3.93e-4 | 0.27, 0.39 | 2.87e-4, 4.06e-4 | 19.96 / 20.26 | fastker 1.75e-4, 1.67e-4; default –, – | fastker 1.81 / 1.32 |
|  | 1e-6 | 1.38e-7, 2.81e-7 | 0.14, 0.28 | 1.47e-7, 2.78e-7 | 25.94 / 28.00 | fastker 5.69e-8, 6.31e-8; default –, – | fastker 1.37 / 1.38 |
|  | 1e-8 | 2.28e-9, 1.91e-9 | 0.23, 0.19 | 1.46e-9, 1.93e-9 | 41.84 / 43.31 | fastker 4.21e-11, 5.03e-11; default 4.22e-11, 5.05e-11 | fastker 0.17 / 0.22; default 0.17 / 0.22 |
| W1, constant, kD = 10 | 1e-3 | 9.45e-16, 1.42e-15 | 0.00, 0.00 | 1.32e-15, 1.80e-15 | 0.00 / 0.00 | default 8.95e-6, 3.39e-5 | default 0.03 / 0.15 |
|  | 1e-6 | 9.45e-16, 1.42e-15 | 0.00, 0.00 | 1.32e-15, 1.80e-15 | 0.00 / 0.00 | default 9.11e-16, 1.22e-15 | default 0.00 / 0.00 |
|  | 1e-8 | 9.45e-16, 1.42e-15 | 0.00, 0.00 | 1.32e-15, 1.80e-15 | 0.00 / 0.00 | default 9.11e-16, 1.22e-15 | default 0.00 / 0.00 |
| W1, constant, kD = 30 | 1e-3 | 9.07e-5, 3.05e-4 | 0.09, 0.30 | 1.07e-4, 3.43e-4 | 1.54 / 5.55 | default 1.13e-4, 2.21e-4 | default 0.63 / 1.38 |
|  | 1e-6 | 4.99e-8, 2.03e-7 | 0.05, 0.20 | 5.41e-8, 2.55e-7 | 0.76 / 5.45 | default 1.43e-8, 1.94e-8 | default 0.04 / 0.12 |
|  | 1e-8 | 7.88e-13, 2.01e-12 | 0.00, 0.00 | 7.07e-13, 1.63e-12 | 0.00 / 0.00 | default 9.56e-12, 1.42e-11 | default 0.00 / 0.00 |
| W1, constant, kD = 100 | 1e-3 | 2.64e-4, 4.41e-4 | 0.26, 0.44 | 2.78e-4, 4.38e-4 | 2.92 / 8.24 | default 2.62e-3, 2.55e-3 | default 24.24 / 15.61 |
|  | 1e-6 | 2.01e-7, 4.41e-7 | 0.20, 0.44 | 2.11e-7, 4.34e-7 | 5.55 / 15.02 | default 6.33e-7, 9.91e-7 | default 7.28 / 11.60 |
|  | 1e-8 | 3.22e-9, 5.48e-9 | 0.32, 0.55 | 3.99e-9, 6.82e-9 | 22.48 / 34.85 | default 2.81e-10, 5.16e-10 | default 0.45 / 0.78 |
| W1, planewave, kD = 500 | 1e-3 | 2.11e-4, 4.19e-4 | 0.21, 0.42 | 1.90e-4, 3.95e-4 | 27.92 / 39.87 | default 9.21e-5, 1.10e-4 | default 3.04 / 3.20 |
|  | 1e-6 | 6.59e-8, 2.18e-7 | 0.07, 0.22 | 8.86e-8, 2.26e-7 | 26.02 / 42.19 | default 3.66e-8, 5.91e-8 | default 2.22 / 3.19 |
|  | 1e-8 | 6.82e-10, 1.19e-9 | 0.07, 0.12 | 9.93e-10, 1.88e-9 | 35.43 / 67.77 | default 2.64e-11, 4.32e-11 | default 0.43 / 0.67 |
| W1, planewave, kD = 1000 | 1e-3 | 2.53e-4, 4.59e-4 | 0.25, 0.46 | 2.14e-4, 4.38e-4 | 35.59 / 43.35 | fastker 8.43e-5, 9.65e-5 | fastker 2.37 / 2.77 |
|  | 1e-6 | 5.25e-8, 2.69e-7 | 0.05, 0.27 | 9.23e-8, 2.27e-7 | 46.68 / 73.86 | fastker 2.45e-8, 3.60e-8 | fastker 2.21 / 3.36 |
|  | 1e-8 | 4.50e-10, 9.45e-10 | 0.04, 0.09 | 9.96e-10, 1.96e-9 | 69.49 / 114.74 | fastker 1.69e-11, 2.61e-11 | fastker 0.36 / 0.56 |
| W2, kD = 10 | 1e-3 | 2.42e-15, 7.13e-15 | 0.00, 0.00 | 3.05e-15, 8.68e-15 | 0.00 / 0.00 | default 4.04e-5, 5.36e-5; fastker 1.45e-15, 4.73e-15 | default 0.26 / 0.37; fastker 0.00 / 0.00 |
|  | 1e-6 | 2.42e-15, 7.13e-15 | 0.00, 0.00 | 3.05e-15, 8.68e-15 | 0.00 / 0.00 | default 4.48e-9, 7.99e-9; fastker 1.45e-15, 4.73e-15 | default 0.03 / 0.07; fastker 0.00 / 0.00 |
|  | 1e-8 | 2.42e-15, 7.13e-15 | 0.00, 0.00 | 3.05e-15, 8.68e-15 | 0.00 / 0.00 | default 1.72e-15, 4.81e-15; fastker 1.45e-15, 4.73e-15 | default 0.00 / 0.00; fastker 0.00 / 0.00 |
| W2, kD = 30 | 1e-3 | 9.09e-5, 1.23e-4 | 0.09, 0.12 | 8.86e-5, 1.22e-4 | 1.44 / 2.39 | default 5.07e-5, 9.18e-5; fastker 3.74e-5, 2.34e-5 | default 0.34 / 0.68; fastker 0.22 / 0.19 |
|  | 1e-6 | 8.97e-8, 1.57e-7 | 0.09, 0.16 | 1.00e-7, 2.08e-7 | 2.34 / 7.99 | default 7.70e-9, 1.77e-8; fastker 5.10e-9, 3.57e-9 | default 0.05 / 0.17; fastker 0.06 / 0.07 |
|  | 1e-8 | 2.39e-12, 1.90e-12 | 0.00, 0.00 | 2.51e-12, 2.95e-12 | 0.00 / 0.02 | default 2.90e-12, 3.44e-12; fastker 2.77e-12, 3.35e-12 | default 0.00 / 0.01; fastker 0.00 / 0.01 |
| W2, kD = 100 | 1e-3 | 1.46e-4, 1.27e-4 | 0.15, 0.13 | 1.32e-4, 1.23e-4 | 3.18 / 2.47 | default 7.97e-4, 4.32e-4; fastker 7.94e-4, 4.18e-4 | default 5.99 / 2.81; fastker 6.03 / 2.82 |
|  | 1e-6 | 1.24e-7, 2.32e-7 | 0.12, 0.23 | 1.47e-7, 2.17e-7 | 10.59 / 11.31 | default 1.89e-7, 1.54e-7; fastker 1.89e-7, 1.51e-7 | default 4.87 / 3.66; fastker 4.85 / 3.65 |
|  | 1e-8 | 4.66e-10, 4.49e-10 | 0.05, 0.04 | 5.45e-10, 4.70e-10 | 3.32 / 2.88 | default 9.74e-11, 9.16e-11; fastker 9.71e-11, 9.13e-11 | default 0.24 / 0.29; fastker 0.24 / 0.29 |
| W2, kD = 250 | 1e-3 | 1.96e-4, 1.97e-4 | 0.20, 0.20 | 2.20e-4, 2.32e-4 | 10.11 / 5.36 | default 1.44e-4, 1.20e-4; fastker 1.39e-4, 7.10e-5 | default 1.19 / 1.21; fastker 1.09 / 0.51 |
|  | 1e-6 | 1.18e-7, 1.51e-7 | 0.12, 0.15 | 1.40e-7, 2.27e-7 | 11.93 / 15.11 | default 4.43e-8, 3.43e-8; fastker 4.40e-8, 2.48e-8 | default 0.79 / 0.59; fastker 0.79 / 0.54 |
|  | 1e-8 | 7.24e-10, 4.80e-10 | 0.07, 0.05 | 9.12e-10, 7.54e-10 | 10.84 / 8.47 | default 3.08e-11, 2.53e-11; fastker 3.07e-11, 1.95e-11 | default 0.08 / 0.07; fastker 0.08 / 0.07 |
| W2, kD = 500 | 1e-3 | 2.43e-4, 2.64e-4 | 0.24, 0.26 | 2.54e-4, 2.68e-4 | 13.44 / 10.78 | default –, –; fastker 1.50e-4, 8.45e-5 | fastker 1.66 / 0.85 |
|  | 1e-6 | 1.38e-7, 2.35e-7 | 0.14, 0.24 | 1.45e-7, 2.40e-7 | 11.57 / 15.72 | default 5.15e-8, 3.72e-8; fastker 5.09e-8, 3.00e-8 | default 1.10 / 0.75; fastker 1.10 / 0.73 |
|  | 1e-8 | 8.34e-10, 4.83e-10 | 0.08, 0.05 | 1.03e-9, 7.70e-10 | 13.81 / 13.74 | default 4.55e-11, 3.22e-11; fastker 4.54e-11, 2.88e-11 | default 0.16 / 0.10; fastker 0.16 / 0.10 |
| W2, constant, kD = 10 | 1e-3 | 1.14e-15, 4.58e-15 | 0.00, 0.00 | 1.45e-15, 5.65e-15 | 0.00 / 0.00 | default 1.97e-5, 1.51e-4 | default 0.14 / 0.81 |
|  | 1e-6 | 1.14e-15, 4.58e-15 | 0.00, 0.00 | 1.45e-15, 5.65e-15 | 0.00 / 0.00 | default 2.74e-9, 6.34e-9 | default 0.01 / 0.05 |
|  | 1e-8 | 1.14e-15, 4.58e-15 | 0.00, 0.00 | 1.45e-15, 5.65e-15 | 0.00 / 0.00 | default 9.68e-16, 3.33e-15 | default 0.00 / 0.00 |
| W2, constant, kD = 30 | 1e-3 | 8.04e-5, 6.37e-4 | 0.08, 0.64 | 7.86e-5, 7.01e-4 | 0.90 / 11.69 | default 6.72e-5, 3.90e-4 | default 1.03 / 2.48 |
|  | 1e-6 | 2.94e-8, 1.74e-7 | 0.03, 0.17 | 2.90e-8, 2.09e-7 | 0.44 / 7.53 | default 1.26e-8, 2.70e-8 | default 0.06 / 0.70 |
|  | 1e-8 | 8.51e-12, 1.59e-11 | 0.00, 0.00 | 8.41e-12, 1.60e-11 | 0.00 / 0.00 | default 8.70e-12, 1.70e-11 | default 0.00 / 0.01 |
| W2, constant, kD = 100 | 1e-3 | 1.28e-4, 6.35e-4 | 0.13, 0.64 | 1.24e-4, 6.39e-4 | 2.75 / 24.24 | default 1.35e-3, 3.59e-3 | default 15.92 / 33.92 |
|  | 1e-6 | 1.11e-7, 7.56e-7 | 0.11, 0.76 | 1.05e-7, 7.66e-7 | 3.18 / 34.68 | default 2.64e-7, 1.11e-6 | default 11.10 / 40.49 |
|  | 1e-8 | 8.83e-10, 3.69e-9 | 0.09, 0.37 | 8.05e-10, 3.76e-9 | 6.56 / 37.55 | default 9.57e-11, 4.24e-10 | default 0.49 / 2.12 |
| W3 (N = 4032), kD = 10 | 1e-3 | 3.88e-5, 6.58e-5 | 0.04, 0.07 | 3.82e-5, 6.31e-5 | 0.40 / 0.55 | default 1.73e-5, 2.27e-5; fastker 1.73e-5, 2.27e-5 | default 0.11 / 0.13; fastker 0.11 / 0.13 |
|  | 1e-6 | 4.01e-8, 1.24e-7 | 0.04, 0.12 | 4.46e-8, 1.38e-7 | 0.83 / 3.32 | default 9.66e-9, 2.09e-8; fastker 1.56e-15, 1.29e-14 | default 0.07 / 0.16; fastker 0.00 / 0.00 |
|  | 1e-8 | 2.52e-15, 1.96e-14 | 0.00, 0.00 | 3.38e-15, 1.71e-14 | 0.00 / 0.00 | default 7.52e-12, 1.85e-11; fastker 1.56e-15, 1.29e-14 | default 0.01 / 0.02; fastker 0.00 / 0.00 |
| W3 (N = 108849), kD = 30 | 1e-3 | 1.60e-4, 3.26e-4 | 0.16, 0.33 | 1.62e-4, 1.17e-4 | 2.30 / 1.65 | default 1.01e-4, 8.82e-5; fastker 1.00e-4, 8.17e-5 | default 0.87 / 0.22; fastker 0.86 / 0.22 |
|  | 1e-6 | 1.44e-7, 3.60e-7 | 0.14, 0.36 | 1.28e-7, 1.16e-7 | 4.10 / 3.88 | default 1.58e-8, 4.34e-8; fastker 1.43e-8, 1.78e-8 | default 0.16 / 0.15; fastker 0.16 / 0.09 |
|  | 1e-8 | 5.04e-12, 6.99e-12 | 0.00, 0.00 | 4.85e-12, 2.01e-12 | 0.01 / 0.01 | default 5.81e-12, 1.46e-11; fastker 5.81e-12, 1.46e-11 | default 0.00 / 0.01; fastker 0.00 / 0.01 |
| W3 (N = 4031442), kD = 100 | 1e-3 | 3.19e-4, 3.32e-4 | 0.32, 0.33 | 2.82e-4, 1.88e-4 | 11.05 / 6.34 | default 2.03e-3, 1.31e-3; fastker 2.03e-3, 1.30e-3 | default 19.47 / 6.96; fastker 19.47 / 6.95 |
|  | 1e-6 | 2.60e-7, 3.63e-7 | 0.26, 0.36 | 1.64e-7, 1.73e-7 | 13.54 / 9.64 | default 5.45e-7, 5.17e-7; fastker 5.45e-7, 5.17e-7 | default 12.06 / 6.53; fastker 12.04 / 6.53 |
|  | 1e-8 | 1.18e-9, 1.82e-9 | 0.12, 0.18 | 9.18e-10, 7.82e-10 | 9.72 / 7.97 | default 1.89e-10, 2.16e-10; fastker 1.88e-10, 2.16e-10 | default 0.57 / 0.58; fastker 0.57 / 0.58 |
| W3 (N = 1000000), kD = 200 | 1e-3 | 3.25e-4, 4.89e-4 | 0.32, 0.49 | 3.29e-4, 2.50e-4 | 10.82 / 7.07 | default 1.68e-3, 1.62e-3; fastker 1.68e-3, 1.62e-3 | default 13.04 / 6.20; fastker 13.03 / 6.19 |
|  | 1e-6 | 1.70e-7, 2.73e-7 | 0.17, 0.27 | 1.69e-7, 1.54e-7 | 10.21 / 10.58 | default 4.64e-7, 6.17e-7; fastker 4.62e-7, 6.14e-7 | default 7.40 / 5.64; fastker 7.46 / 5.67 |
|  | 1e-8 | 1.43e-9, 2.64e-9 | 0.14, 0.26 | 1.41e-9, 1.42e-9 | 9.42 / 7.67 | default 2.13e-10, 3.44e-10; fastker 2.13e-10, 3.44e-10 | default 0.49 / 0.53; fastker 0.49 / 0.53 |
| W3 (N = 1000000), kD = 500 | 1e-3 | 3.50e-4, 4.50e-4 | 0.35, 0.45 | 3.61e-4, 4.14e-4 | 16.43 / 15.86 | default 2.50e-4, 2.54e-4; fastker 2.31e-4, 2.27e-4 | default 2.42 / 1.47; fastker 2.61 / 1.35 |
|  | 1e-6 | 1.91e-7, 2.53e-7 | 0.19, 0.25 | 1.58e-7, 1.99e-7 | 15.40 / 15.41 | default 9.44e-8, 1.15e-7; fastker 8.75e-8, 1.01e-7 | default 1.86 / 1.38; fastker 1.84 / 1.30 |
|  | 1e-8 | 2.04e-9, 3.42e-9 | 0.20, 0.34 | 1.99e-9, 2.70e-9 | 24.59 / 24.76 | default 8.68e-11, 1.00e-10; fastker 8.68e-11, 1.00e-10 | default 0.22 / 0.19; fastker 0.22 / 0.19 |
| W3 (N = 4032), constant, kD = 10 | 1e-3 | 8.02e-6, 4.31e-5 | 0.01, 0.04 | 8.55e-6, 4.45e-5 | 0.10 / 0.48 | default 1.87e-5, 3.01e-5 | default 0.05 / 0.13 |
|  | 1e-6 | 1.39e-8, 6.40e-8 | 0.01, 0.06 | 1.45e-8, 8.10e-8 | 0.11 / 1.64 | default 1.23e-8, 1.69e-8 | default 0.03 / 0.10 |
|  | 1e-8 | 2.37e-15, 1.27e-14 | 0.00, 0.00 | 3.34e-15, 1.15e-14 | 0.00 / 0.00 | default 1.90e-11, 2.07e-11 | default 0.00 / 0.02 |
| W3 (N = 108849), constant, kD = 30 | 1e-3 | 1.28e-4, 4.00e-4 | 0.13, 0.40 | 1.05e-4, 2.16e-4 | 1.95 / 2.66 | default 1.92e-4, 2.72e-4 | default 1.46 / 0.86 |
|  | 1e-6 | 4.77e-8, 1.91e-7 | 0.05, 0.19 | 4.77e-8, 1.12e-7 | 1.43 / 3.15 | default 3.03e-8, 4.50e-8 | default 0.18 / 0.24 |
|  | 1e-8 | 2.38e-11, 2.92e-11 | 0.00, 0.00 | 2.44e-11, 1.68e-11 | 0.01 / 0.01 | default 2.38e-11, 3.04e-11 | default 0.01 / 0.01 |

Largest of ours over ε: 0.76; largest worst box over ε: 114.7; misses (over ε, or a worst box above 100 ε): W1, planewave, kD = 1000 at 1e-8.

## 3. The worst box

Our worst box over ε grows with kD and with tighter ε on W1 (random charges, 10⁻⁸: 8.7 ε at
kD = 100, 21.6 at 250, 40.4 at 500, 43.3 at 1000); the plane wave's trace concentrates more (67.8 ε
at kD = 500, 114.7 ε at kD = 1000, 10⁻⁸, both ∇φ). Over our own all-target ℓ2 error the worst box
stands 224–698× out on W1 at kD = 1000 at 10⁻⁸ (the log's five worst boxes), against FMM3D's
33–146×. The worst boxes are those where the sphere touches its bounding cube: (127,63,64) and
(63,63,127) on level 7 at kD = 1000, (31,63,32) on level 6 at kD = 500, the centres of the cube's
faces; every listed box is one of the six tangent points' neighbours. This is CH4.0b's finding (a
level's error at targets near a box corner is 19–45× its RMS, and W1's six tangent points are box
vertices on every level), now on the high-frequency trees: the ℓ2 error meets ε with a margin
(0.09 ε in ∇φ at the samples, 0.20 ε over all targets, on the case above 100 ε), the box measure
does not.

One case × ε of 78 exceeds §9.1's 100 ε, by 15 %. It is reported as a finding, not hidden by a
test change: question 2 asks whether CH4.4a closes with it and CH4.4b (which maps the box errors to
levels and offset classes) decides its cause and remedy. FMM3D's worst box is at most 3.4 ε at
kD ≥ 500, and up to 40 ε at kD = 100 (W2 with the constant density, 10⁻⁶), where its ℓ2 error can
miss its eps (2.6 · 10⁻³ at 10⁻³ on W1 with the constant density, as CH6.2d found on random charges).

## 4. Time and memory

At kD ≥ 250 and on every case of N ≥ 10⁶ ours is 0.08–0.51× FMM3D's faster variant's call (W1 at
kD = 1000: 28.0, 32.3 and 40.3 s against `fastker`'s 287, 362 and 416 s); its evaluation alone
0.03–0.42×. On the small trees ours is slower at 10⁻³ only, the first ε of each launch, where the
build includes the process's one-time work (0.11 s against 0.02 s at 10⁻⁶ on W1 at kD = 30): the
first ε pays for what later ε reuse, in time as in memory (§5). These ratios are at the same
requested ε; FMM3D's errors are often below ours there (W1 at kD = 1000, 10⁻⁸: 4.2 · 10⁻¹¹ against
2.3 · 10⁻⁹) and at 10⁻³ on kD = 100 above its eps. Matched accuracy is CH4.5's.

FMM3D `default` needs much more memory than `fastker` at 10⁻³ (W1 at kD = 500: 108 against 27 GB;
W3 at kD = 100: 120 against 16 GB) and did not fit the node at kD = 1000: killed by the job's
180 GB limit at 10⁻³ and 10⁻⁶ (getrusage's last peak 35 and 193 GB) and, in job 663708, at 10⁻¹²
(188 GB); 191 GB at 10⁻⁸; and on W2 at kD = 500 at 10⁻³ (193 GB). Those rows read "failed (no output)" with the peak; `compare`
reports them since this task instead of stopping. CH4.5's comparison at kD = 1000 therefore has
`fastker` alone on one node (question 3).

Our memory grows with N and the degrees: on W1 per point 2.4–6.3 kB at kD = 100, 2.7–3.6 kB at
kD = 500 and 1000 (8 · 10⁶ points: 21–29 GB on the node, at most 7.3 GB on one rank). At the
tightest ε ours stays below FMM3D `fastker`'s peak from kD = 250 on W1 (W1 at kD = 1000,
10⁻⁸: 28.7 against 57.2 GB), and above it on W3 at N = 10⁶ (10.5 against 6.4 GB at kD = 500), whose
dense levels hold many boxes per point.

## 5. The memory figures and the order of ε

One launch measures the three ε in turn. Freed memory that the allocator keeps is reused by the
next ε without raising its peak, so a later ε's growth is a lower bound. Job 663737 ran W1 at
kD = 500 alone (no FMM3D outputs read) with ε in both orders: 10⁻⁸ grows 8.59 GB first and 6.82 GB
third, 10⁻³ 6.47 GB first and 4.38 GB third; 10⁻⁶ 4.99–5.22 GB; the checksums the same in both
orders. So the sweep's 10⁻³ figures (its first ε) are a fresh process's, and its 10⁻⁶ and 10⁻⁸
figures may be up to about 20–30 % low. The first ε also carries the process's one-time
allocations: on the smallest trees with a far field 0.2–1.0 GB at 10⁻³ against 0.03–0.18 GB later
(W1 at kD = 30: 0.87 against 0.03 GB). FMM3D's figures are per launch, one eps each. Question 4 proposes one ε per
launch for CH4.5's memory.

## 6. The slow nodes

The first jobs ran W1 at kD = 500 (both densities) and the three larger W3 cases on node-c11b-004
(jobs 663707, 663710) and W1 at kD = 1000 on node-c11b-013 (663708). There FMM3D took about 3× the
time CH3.6 and CH4.0e measured with the same binaries (`fastker`, W1 at kD = 500, 10⁻³: 103 s
against 33.4 s), and ours too (W1 at kD = 500, 10⁻³: build 8.6 s against 2.6 s); small problems ran
at the usual speed on the same nodes, and the outputs were the same bit for bit (our checksums,
FMM3D's errors to every digit shown). Job 663735 ran `fastker` on that case four ways on
node-c11a-143: plain and under `maxrss.py`, FULL = 0 and 1: 32.9 s each, so neither the wrapper nor
the full output is the cause. Those cases were run again with both nodes excluded (663736) and only
those runs are reported; the slow runs are kept in `slow-nodes/` of the raw outputs.
`tools/kathleen/README.md` now names both nodes for exclusion from timed jobs. The FMM3D times
reported here agree with CH3.6 and CH4.0e where they overlap (W1 at kD = 500: `default` 49.3, 62.8,
70.9 s; CH3.6 49.3, 62.9, 71.7 s).

## 7. Build per k

The build follows the top level's degree, not N (W3 at N = 10⁶ and kD = 500 builds in 2.3 s, W1 at
kD = 500 with twice the points in 2.6–2.9 s). From kD = 500 to 1000 it grows 6–7× (2.6–2.9 s to
16–19 s) for 4× the points: the RCR tables of level 2 (p = 225–248) are 9.2–12.9 s of it, the
diagonal tables 3.4–4.6 s; the octree 0.45 s and the plan 0.8 s. At kD = 1000 the build is
40–67 % of build plus one evaluation (CH4.0e found 39–43 % on W1 at kD = 500 before H4's work).
The proposal's per-level cache keys (§6.3; a sweep in octaves rebuilding only the top level) are
not used by the harness, which builds each case afresh.

## 8. Questions for sign-off

1. **The acceptance as read** (plan, "Acceptance"): CH4.4a without the mapping to levels and
   classes (CH4.4b); our side under `Case::defaults`; one complex charge vector (CH3.6's reading of
   §9.1's eight); "W3 to 200" as resolved to kD = 100 and under-resolved at N = 10⁶ for kD = 200 and
   500; the constant density on kD ≤ 100; FMM3D `fastker` as the crosscheck's reference at
   kD = 1000, where `default` does not fit the node (its own error at the samples 3.9 · 10⁻¹³ and
   2.4 · 10⁻¹³). Recommended: as read.
2. **The worst box of W1 at kD = 1000 with the plane wave's trace, 10⁻⁸: 114.7 ε** (§3), above
   §9.1's 100 ε, the ℓ2 errors within ε. Recommended: CH4.4a closes with it recorded as a finding;
   CH4.4b maps it (its *Done when* names the worst boxes' levels and classes) and proposes the
   remedy (a degree on the responsible level, or §9.1's bound), confirmed by pair and FMM evidence
   as the H6 gate's rule asks. Alternative: CH4.4a stays open until a fix lands.
3. **FMM3D `default` at kD = 1000 does not fit one Kathleen node** (180 GB; §4). Recommended: an
   Also on CH4.5: at kD = 1000 `fastker` alone, `default`'s failure and peaks stated.
4. **Memory per ε.** Recommended: an Also on CH4.5: our memory one ε per launch (a later ε's
   figure in one launch is a lower bound, §5).


## 9. Tables

### Time

Seconds, single launches; ours the slowest rank's build and evaluation; FMM3D's call, its setup
included; ratios ours (build + evaluation, then evaluation alone) over each variant's call ("–": no
output).

| Case | N | ε | ours build s | ours run s | ours total s | FMM3D default s | FMM3D fastker s | ours ÷ default (total, run) | ours ÷ fastker (total, run) |
| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| W1, kD = 10 | 796 | 1e-3 | 0.01 | 0.00 | 0.01 | 0.01 | – | 1.00, 0.00 | – |
|  | 796 | 1e-6 | 0.00 | 0.00 | 0.00 | 0.01 | 0.01 | 0.00, 0.00 | 0.00, 0.00 |
|  | 796 | 1e-8 | 0.00 | 0.00 | 0.00 | 0.01 | 0.02 | 0.00, 0.00 | 0.00, 0.00 |
| W1, kD = 30 | 7,162 | 1e-3 | 0.11 | 0.01 | 0.12 | 0.06 | 0.03 | 2.00, 0.17 | 4.00, 0.33 |
|  | 7,162 | 1e-6 | 0.02 | 0.01 | 0.03 | 0.08 | 0.06 | 0.38, 0.12 | 0.50, 0.17 |
|  | 7,162 | 1e-8 | 0.02 | 0.03 | 0.05 | 0.11 | 0.06 | 0.45, 0.27 | 0.83, 0.50 |
| W1, kD = 100 | 79,578 | 1e-3 | 0.09 | 0.07 | 0.16 | 0.92 | 0.43 | 0.17, 0.08 | 0.37, 0.16 |
|  | 79,578 | 1e-6 | 0.13 | 0.09 | 0.22 | 0.95 | 0.51 | 0.23, 0.09 | 0.43, 0.18 |
|  | 79,578 | 1e-8 | 0.13 | 0.23 | 0.36 | 1.56 | 0.67 | 0.23, 0.15 | 0.54, 0.34 |
| W1, kD = 250 | 497,360 | 1e-3 | 0.46 | 0.39 | 0.85 | 8.12 | 4.46 | 0.10, 0.05 | 0.19, 0.09 |
|  | 497,360 | 1e-6 | 0.60 | 0.57 | 1.17 | 9.79 | 5.54 | 0.12, 0.06 | 0.21, 0.10 |
|  | 497,360 | 1e-8 | 0.64 | 1.24 | 1.88 | 11.86 | 6.91 | 0.16, 0.10 | 0.27, 0.18 |
| W1, kD = 500 | 1,989,437 | 1e-3 | 2.56 | 1.88 | 4.44 | 49.30 | 33.14 | 0.09, 0.04 | 0.13, 0.06 |
|  | 1,989,437 | 1e-6 | 2.88 | 2.71 | 5.59 | 62.82 | 43.25 | 0.09, 0.04 | 0.13, 0.06 |
|  | 1,989,437 | 1e-8 | 2.78 | 5.50 | 8.28 | 70.88 | 48.74 | 0.12, 0.08 | 0.17, 0.11 |
| W1, kD = 1000 | 7,957,748 | 1e-3 | 18.74 | 9.25 | 27.99 | – | 287.14 | – | 0.10, 0.03 |
|  | 7,957,748 | 1e-6 | 19.46 | 12.85 | 32.31 | – | 361.73 | – | 0.09, 0.04 |
|  | 7,957,748 | 1e-8 | 16.27 | 24.05 | 40.32 | 508.27 | 416.48 | 0.08, 0.05 | 0.10, 0.06 |
| W1, constant, kD = 10 | 796 | 1e-3 | 0.00 | 0.00 | 0.00 | 0.01 | – | 0.00, 0.00 | – |
|  | 796 | 1e-6 | 0.00 | 0.00 | 0.00 | 0.01 | – | 0.00, 0.00 | – |
|  | 796 | 1e-8 | 0.00 | 0.00 | 0.00 | 0.01 | – | 0.00, 0.00 | – |
| W1, constant, kD = 30 | 7,162 | 1e-3 | 0.11 | 0.01 | 0.12 | 0.06 | – | 2.00, 0.17 | – |
|  | 7,162 | 1e-6 | 0.02 | 0.01 | 0.03 | 0.09 | – | 0.33, 0.11 | – |
|  | 7,162 | 1e-8 | 0.02 | 0.02 | 0.04 | 0.12 | – | 0.33, 0.17 | – |
| W1, constant, kD = 100 | 79,578 | 1e-3 | 0.08 | 0.06 | 0.14 | 0.93 | – | 0.15, 0.06 | – |
|  | 79,578 | 1e-6 | 0.12 | 0.10 | 0.22 | 0.94 | – | 0.23, 0.11 | – |
|  | 79,578 | 1e-8 | 0.13 | 0.23 | 0.36 | 1.58 | – | 0.23, 0.15 | – |
| W1, planewave, kD = 500 | 1,989,437 | 1e-3 | 2.58 | 1.88 | 4.46 | 49.34 | – | 0.09, 0.04 | – |
|  | 1,989,437 | 1e-6 | 2.84 | 2.73 | 5.57 | 62.24 | – | 0.09, 0.04 | – |
|  | 1,989,437 | 1e-8 | 2.75 | 5.52 | 8.27 | 71.60 | – | 0.12, 0.08 | – |
| W1, planewave, kD = 1000 | 7,957,748 | 1e-3 | 18.42 | 9.08 | 27.50 | – | 275.79 | – | 0.10, 0.03 |
|  | 7,957,748 | 1e-6 | 19.34 | 12.49 | 31.83 | – | 354.76 | – | 0.09, 0.04 |
|  | 7,957,748 | 1e-8 | 15.96 | 23.46 | 39.42 | – | 408.50 | – | 0.10, 0.06 |
| W2, kD = 10 | 1,276 | 1e-3 | 0.00 | 0.00 | 0.00 | 0.02 | 0.01 | 0.00, 0.00 | 0.00, 0.00 |
|  | 1,276 | 1e-6 | 0.00 | 0.00 | 0.00 | 0.04 | 0.01 | 0.00, 0.00 | 0.00, 0.00 |
|  | 1,276 | 1e-8 | 0.00 | 0.00 | 0.00 | 0.02 | 0.01 | 0.00, 0.00 | 0.00, 0.00 |
| W2, kD = 30 | 11,480 | 1e-3 | 0.16 | 0.02 | 0.18 | 0.12 | 0.08 | 1.50, 0.17 | 2.25, 0.25 |
|  | 11,480 | 1e-6 | 0.03 | 0.03 | 0.06 | 0.22 | 0.16 | 0.27, 0.14 | 0.38, 0.19 |
|  | 11,480 | 1e-8 | 0.04 | 0.09 | 0.13 | 0.27 | 0.25 | 0.48, 0.33 | 0.52, 0.36 |
| W2, kD = 100 | 127,554 | 1e-3 | 0.14 | 0.09 | 0.23 | 1.34 | 0.59 | 0.17, 0.07 | 0.39, 0.15 |
|  | 127,554 | 1e-6 | 0.15 | 0.16 | 0.31 | 1.98 | 0.90 | 0.16, 0.08 | 0.34, 0.18 |
|  | 127,554 | 1e-8 | 0.15 | 0.43 | 0.58 | 2.63 | 1.34 | 0.22, 0.16 | 0.43, 0.32 |
| W2, kD = 250 | 797,211 | 1e-3 | 0.61 | 0.50 | 1.11 | 13.35 | 6.05 | 0.08, 0.04 | 0.18, 0.08 |
|  | 797,211 | 1e-6 | 0.85 | 0.87 | 1.72 | 14.52 | 7.86 | 0.12, 0.06 | 0.22, 0.11 |
|  | 797,211 | 1e-8 | 0.82 | 2.04 | 2.86 | 21.73 | 10.36 | 0.13, 0.09 | 0.28, 0.20 |
| W2, kD = 500 | 3,188,840 | 1e-3 | 3.32 | 2.22 | 5.54 | – | 40.51 | – | 0.14, 0.05 |
|  | 3,188,840 | 1e-6 | 4.01 | 3.72 | 7.73 | 79.77 | 51.76 | 0.10, 0.05 | 0.15, 0.07 |
|  | 3,188,840 | 1e-8 | 3.44 | 7.94 | 11.38 | 110.75 | 62.19 | 0.10, 0.07 | 0.18, 0.13 |
| W2, constant, kD = 10 | 1,276 | 1e-3 | 0.00 | 0.00 | 0.00 | 0.02 | – | 0.00, 0.00 | – |
|  | 1,276 | 1e-6 | 0.00 | 0.00 | 0.00 | 0.04 | – | 0.00, 0.00 | – |
|  | 1,276 | 1e-8 | 0.00 | 0.00 | 0.00 | 0.02 | – | 0.00, 0.00 | – |
| W2, constant, kD = 30 | 11,480 | 1e-3 | 0.14 | 0.02 | 0.16 | 0.13 | – | 1.23, 0.15 | – |
|  | 11,480 | 1e-6 | 0.03 | 0.03 | 0.06 | 0.24 | – | 0.25, 0.12 | – |
|  | 11,480 | 1e-8 | 0.04 | 0.10 | 0.14 | 0.28 | – | 0.50, 0.36 | – |
| W2, constant, kD = 100 | 127,554 | 1e-3 | 0.13 | 0.09 | 0.22 | 1.34 | – | 0.16, 0.07 | – |
|  | 127,554 | 1e-6 | 0.15 | 0.15 | 0.30 | 2.01 | – | 0.15, 0.07 | – |
|  | 127,554 | 1e-8 | 0.15 | 0.44 | 0.59 | 2.69 | – | 0.22, 0.16 | – |
| W3 (N = 4032), kD = 10 | 4,032 | 1e-3 | 0.03 | 0.00 | 0.03 | 0.02 | 0.01 | 1.50, 0.00 | 3.00, 0.00 |
|  | 4,032 | 1e-6 | 0.01 | 0.01 | 0.02 | 0.05 | 0.02 | 0.40, 0.20 | 1.00, 0.50 |
|  | 4,032 | 1e-8 | 0.01 | 0.01 | 0.02 | 0.07 | 0.02 | 0.29, 0.14 | 1.00, 0.50 |
| W3 (N = 108849), kD = 30 | 108,849 | 1e-3 | 0.13 | 0.08 | 0.21 | 0.36 | 0.22 | 0.58, 0.22 | 0.95, 0.36 |
|  | 108,849 | 1e-6 | 0.05 | 0.12 | 0.17 | 0.97 | 0.32 | 0.18, 0.12 | 0.53, 0.38 |
|  | 108,849 | 1e-8 | 0.06 | 0.28 | 0.34 | 1.26 | 0.45 | 0.27, 0.22 | 0.76, 0.62 |
| W3 (N = 4031442), kD = 100 | 4,031,442 | 1e-3 | 1.25 | 2.34 | 3.59 | 22.90 | 7.40 | 0.16, 0.10 | 0.49, 0.32 |
|  | 4,031,442 | 1e-6 | 1.18 | 4.09 | 5.27 | 31.95 | 22.96 | 0.16, 0.13 | 0.23, 0.18 |
|  | 4,031,442 | 1e-8 | 1.19 | 10.38 | 11.57 | 38.51 | 28.17 | 0.30, 0.27 | 0.41, 0.37 |
| W3 (N = 1000000), kD = 200 | 1,000,000 | 1e-3 | 0.48 | 1.41 | 1.89 | 8.83 | 3.93 | 0.21, 0.16 | 0.48, 0.36 |
|  | 1,000,000 | 1e-6 | 0.51 | 1.83 | 2.34 | 14.58 | 5.17 | 0.16, 0.13 | 0.45, 0.35 |
|  | 1,000,000 | 1e-8 | 0.56 | 2.84 | 3.40 | 16.91 | 6.72 | 0.20, 0.17 | 0.51, 0.42 |
| W3 (N = 1000000), kD = 500 | 1,000,000 | 1e-3 | 2.26 | 3.38 | 5.64 | 61.22 | 53.38 | 0.09, 0.06 | 0.11, 0.06 |
|  | 1,000,000 | 1e-6 | 2.35 | 4.29 | 6.64 | 86.72 | 73.23 | 0.08, 0.05 | 0.09, 0.06 |
|  | 1,000,000 | 1e-8 | 2.39 | 7.63 | 10.02 | 93.07 | 83.42 | 0.11, 0.08 | 0.12, 0.09 |
| W3 (N = 4032), constant, kD = 10 | 4,032 | 1e-3 | 0.03 | 0.00 | 0.03 | 0.02 | – | 1.50, 0.00 | – |
|  | 4,032 | 1e-6 | 0.01 | 0.01 | 0.02 | 0.04 | – | 0.50, 0.25 | – |
|  | 4,032 | 1e-8 | 0.01 | 0.01 | 0.02 | 0.08 | – | 0.25, 0.12 | – |
| W3 (N = 108849), constant, kD = 30 | 108,849 | 1e-3 | 0.13 | 0.08 | 0.21 | 0.36 | – | 0.58, 0.22 | – |
|  | 108,849 | 1e-6 | 0.05 | 0.13 | 0.18 | 0.95 | – | 0.19, 0.14 | – |
|  | 108,849 | 1e-8 | 0.05 | 0.29 | 0.34 | 1.30 | – | 0.26, 0.22 | – |


### Memory

Ours: the FMM's resident growth summed over the 4 ranks (the largest rank's), §1 and §5; FMM3D: the
driver's peak resident set at that eps, its arrays included (a failed run's last peak).

| Case | N | ε | ours GB (largest rank) | FMM3D default GB | FMM3D fastker GB |
| --- | ---: | --- | ---: | ---: | ---: |
| W1, kD = 10 | 796 | 1e-3 | 0.01 (rank ≤ 0.00) | 0.01 | 0.00 |
|  | 796 | 1e-6 | 0.01 (rank ≤ 0.00) | 0.00 | 0.00 |
|  | 796 | 1e-8 | 0.01 (rank ≤ 0.00) | 0.00 | 0.00 |
| W1, kD = 30 | 7,162 | 1e-3 | 0.87 (rank ≤ 0.22) | 0.06 | 0.02 |
|  | 7,162 | 1e-6 | 0.03 (rank ≤ 0.01) | 0.03 | 0.00 |
|  | 7,162 | 1e-8 | 0.07 (rank ≤ 0.02) | 0.05 | 0.01 |
| W1, kD = 100 | 79,578 | 1e-3 | 0.19 (rank ≤ 0.05) | 3.51 | 0.69 |
|  | 79,578 | 1e-6 | 0.23 (rank ≤ 0.06) | 1.39 | 0.36 |
|  | 79,578 | 1e-8 | 0.50 (rank ≤ 0.13) | 1.34 | 0.47 |
| W1, kD = 250 | 497,360 | 1e-3 | 1.04 (rank ≤ 0.26) | 25.08 | 5.95 |
|  | 497,360 | 1e-6 | 1.14 (rank ≤ 0.29) | 10.90 | 2.19 |
|  | 497,360 | 1e-8 | 1.33 (rank ≤ 0.35) | 9.65 | 2.79 |
| W1, kD = 500 | 1,989,437 | 1e-3 | 6.45 (rank ≤ 1.63) | 108.15 | 26.96 |
|  | 1,989,437 | 1e-6 | 5.63 (rank ≤ 1.45) | 49.43 | 9.72 |
|  | 1,989,437 | 1e-8 | 7.02 (rank ≤ 1.85) | 43.90 | 12.30 |
| W1, kD = 1000 | 7,957,748 | 1e-3 | 21.18 (rank ≤ 5.31) | 35.06 | 118.02 |
|  | 7,957,748 | 1e-6 | 22.66 (rank ≤ 5.74) | 192.87 | 45.19 |
|  | 7,957,748 | 1e-8 | 28.68 (rank ≤ 7.29) | 191.00 | 57.22 |
| W1, constant, kD = 10 | 796 | 1e-3 | 0.01 (rank ≤ 0.00) | 0.01 | – |
|  | 796 | 1e-6 | 0.01 (rank ≤ 0.00) | 0.01 | – |
|  | 796 | 1e-8 | 0.01 (rank ≤ 0.00) | 0.00 | – |
| W1, constant, kD = 30 | 7,162 | 1e-3 | 0.86 (rank ≤ 0.22) | 0.06 | – |
|  | 7,162 | 1e-6 | 0.03 (rank ≤ 0.01) | 0.03 | – |
|  | 7,162 | 1e-8 | 0.07 (rank ≤ 0.02) | 0.05 | – |
| W1, constant, kD = 100 | 79,578 | 1e-3 | 0.19 (rank ≤ 0.05) | 3.52 | – |
|  | 79,578 | 1e-6 | 0.23 (rank ≤ 0.06) | 1.39 | – |
|  | 79,578 | 1e-8 | 0.50 (rank ≤ 0.13) | 1.34 | – |
| W1, planewave, kD = 500 | 1,989,437 | 1e-3 | 6.49 (rank ≤ 1.63) | 108.14 | – |
|  | 1,989,437 | 1e-6 | 5.64 (rank ≤ 1.44) | 49.43 | – |
|  | 1,989,437 | 1e-8 | 7.04 (rank ≤ 1.86) | 43.89 | – |
| W1, planewave, kD = 1000 | 7,957,748 | 1e-3 | 21.16 (rank ≤ 5.32) | – | 118.02 |
|  | 7,957,748 | 1e-6 | 22.58 (rank ≤ 5.66) | – | 45.24 |
|  | 7,957,748 | 1e-8 | 28.81 (rank ≤ 7.27) | – | 57.22 |
| W2, kD = 10 | 1,276 | 1e-3 | 0.02 (rank ≤ 0.01) | 0.01 | 0.01 |
|  | 1,276 | 1e-6 | 0.01 (rank ≤ 0.01) | 0.02 | 0.01 |
|  | 1,276 | 1e-8 | 0.01 (rank ≤ 0.00) | 0.00 | 0.00 |
| W2, kD = 30 | 11,480 | 1e-3 | 0.87 (rank ≤ 0.22) | 0.12 | 0.04 |
|  | 11,480 | 1e-6 | 0.06 (rank ≤ 0.02) | 0.10 | 0.08 |
|  | 11,480 | 1e-8 | 0.18 (rank ≤ 0.05) | 0.12 | 0.12 |
| W2, kD = 100 | 127,554 | 1e-3 | 0.37 (rank ≤ 0.09) | 4.59 | 0.91 |
|  | 127,554 | 1e-6 | 0.27 (rank ≤ 0.07) | 2.83 | 0.60 |
|  | 127,554 | 1e-8 | 0.72 (rank ≤ 0.18) | 1.78 | 0.79 |
| W2, kD = 250 | 797,211 | 1e-3 | 1.62 (rank ≤ 0.41) | 48.72 | 11.56 |
|  | 797,211 | 1e-6 | 1.77 (rank ≤ 0.45) | 20.48 | 4.46 |
|  | 797,211 | 1e-8 | 2.51 (rank ≤ 0.64) | 17.34 | 5.74 |
| W2, kD = 500 | 3,188,840 | 1e-3 | 7.43 (rank ≤ 1.90) | 192.87 | 51.49 |
|  | 3,188,840 | 1e-6 | 8.26 (rank ≤ 2.08) | 87.65 | 18.08 |
|  | 3,188,840 | 1e-8 | 11.15 (rank ≤ 2.93) | 75.60 | 23.20 |
| W2, constant, kD = 10 | 1,276 | 1e-3 | 0.01 (rank ≤ 0.01) | 0.01 | – |
|  | 1,276 | 1e-6 | 0.01 (rank ≤ 0.00) | 0.02 | – |
|  | 1,276 | 1e-8 | 0.01 (rank ≤ 0.00) | 0.01 | – |
| W2, constant, kD = 30 | 11,480 | 1e-3 | 0.87 (rank ≤ 0.22) | 0.12 | – |
|  | 11,480 | 1e-6 | 0.05 (rank ≤ 0.01) | 0.10 | – |
|  | 11,480 | 1e-8 | 0.18 (rank ≤ 0.05) | 0.12 | – |
| W2, constant, kD = 100 | 127,554 | 1e-3 | 0.37 (rank ≤ 0.10) | 4.59 | – |
|  | 127,554 | 1e-6 | 0.27 (rank ≤ 0.07) | 2.83 | – |
|  | 127,554 | 1e-8 | 0.72 (rank ≤ 0.18) | 1.78 | – |
| W3 (N = 4032), kD = 10 | 4,032 | 1e-3 | 0.18 (rank ≤ 0.05) | 0.01 | 0.01 |
|  | 4,032 | 1e-6 | 0.07 (rank ≤ 0.02) | 0.02 | 0.00 |
|  | 4,032 | 1e-8 | 0.06 (rank ≤ 0.01) | 0.04 | 0.00 |
| W3 (N = 108849), kD = 30 | 108,849 | 1e-3 | 0.97 (rank ≤ 0.24) | 0.53 | 0.09 |
|  | 108,849 | 1e-6 | 0.12 (rank ≤ 0.03) | 0.89 | 0.17 |
|  | 108,849 | 1e-8 | 0.14 (rank ≤ 0.04) | 0.27 | 0.27 |
| W3 (N = 4031442), kD = 100 | 4,031,442 | 1e-3 | 2.79 (rank ≤ 0.70) | 120.00 | 15.67 |
|  | 4,031,442 | 1e-6 | 4.41 (rank ≤ 1.12) | 22.78 | 9.16 |
|  | 4,031,442 | 1e-8 | 7.98 (rank ≤ 2.03) | 30.55 | 12.16 |
| W3 (N = 1000000), kD = 200 | 1,000,000 | 1e-3 | 2.35 (rank ≤ 0.59) | 49.37 | 6.53 |
|  | 1,000,000 | 1e-6 | 2.00 (rank ≤ 0.51) | 62.25 | 8.19 |
|  | 1,000,000 | 1e-8 | 2.87 (rank ≤ 0.73) | 9.83 | 9.82 |
| W3 (N = 1000000), kD = 500 | 1,000,000 | 1e-3 | 7.33 (rank ≤ 1.88) | 23.34 | 3.96 |
|  | 1,000,000 | 1e-6 | 7.26 (rank ≤ 1.84) | 32.15 | 5.14 |
|  | 1,000,000 | 1e-8 | 10.51 (rank ≤ 2.73) | 6.39 | 6.38 |
| W3 (N = 4032), constant, kD = 10 | 4,032 | 1e-3 | 0.18 (rank ≤ 0.05) | 0.01 | – |
|  | 4,032 | 1e-6 | 0.07 (rank ≤ 0.02) | 0.02 | – |
|  | 4,032 | 1e-8 | 0.05 (rank ≤ 0.01) | 0.04 | – |
| W3 (N = 108849), constant, kD = 30 | 108,849 | 1e-3 | 0.95 (rank ≤ 0.24) | 0.53 | – |
|  | 108,849 | 1e-6 | 0.11 (rank ≤ 0.03) | 0.89 | – |
|  | 108,849 | 1e-8 | 0.14 (rank ≤ 0.04) | 0.27 | – |


### Build by part

Ours, seconds, each part the slowest rank's (`BUILD_PARTS`; "RCR tables" to "scratch" lie within
"operator"); W1's sweep and the cases of N ≥ 10⁶.

| Case | ε | build s | octree | plan | operator | RCR tables | matrices | plane-wave tables | diagonal tables | scratch | replan | evaluator | evaluation s |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| W1, kD = 10 | 1e-3 | 0.01 | 0.001 | 0.000 | 0.001 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.001 | 0.00 |
|  | 1e-6 | 0.00 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.00 |
|  | 1e-8 | 0.00 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.00 |
| W1, kD = 30 | 1e-3 | 0.11 | 0.001 | 0.000 | 0.101 | 0.000 | 0.100 | 0.000 | 0.000 | 0.000 | 0.000 | 0.004 | 0.01 |
|  | 1e-6 | 0.02 | 0.001 | 0.000 | 0.019 | 0.000 | 0.000 | 0.000 | 0.017 | 0.001 | 0.000 | 0.004 | 0.01 |
|  | 1e-8 | 0.02 | 0.001 | 0.000 | 0.018 | 0.001 | 0.000 | 0.016 | 0.000 | 0.001 | 0.000 | 0.002 | 0.03 |
| W1, kD = 100 | 1e-3 | 0.09 | 0.005 | 0.007 | 0.066 | 0.003 | 0.000 | 0.000 | 0.060 | 0.001 | 0.000 | 0.007 | 0.07 |
|  | 1e-6 | 0.13 | 0.005 | 0.007 | 0.108 | 0.003 | 0.000 | 0.000 | 0.100 | 0.001 | 0.000 | 0.013 | 0.09 |
|  | 1e-8 | 0.13 | 0.005 | 0.008 | 0.102 | 0.005 | 0.000 | 0.025 | 0.068 | 0.001 | 0.000 | 0.007 | 0.23 |
| W1, kD = 250 | 1e-3 | 0.46 | 0.027 | 0.033 | 0.351 | 0.100 | 0.000 | 0.000 | 0.241 | 0.003 | 0.000 | 0.018 | 0.39 |
|  | 1e-6 | 0.60 | 0.025 | 0.034 | 0.499 | 0.150 | 0.000 | 0.000 | 0.336 | 0.005 | 0.000 | 0.012 | 0.57 |
|  | 1e-8 | 0.64 | 0.025 | 0.039 | 0.523 | 0.105 | 0.000 | 0.020 | 0.381 | 0.008 | 0.000 | 0.013 | 1.24 |
| W1, kD = 500 | 1e-3 | 2.56 | 0.108 | 0.176 | 2.028 | 1.146 | 0.000 | 0.000 | 0.850 | 0.012 | 0.000 | 0.148 | 1.88 |
|  | 1e-6 | 2.88 | 0.104 | 0.178 | 2.353 | 1.198 | 0.000 | 0.000 | 1.102 | 0.024 | 0.000 | 0.162 | 2.71 |
|  | 1e-8 | 2.78 | 0.102 | 0.193 | 2.242 | 0.874 | 0.000 | 0.022 | 1.299 | 0.025 | 0.000 | 0.239 | 5.50 |
| W1, kD = 1000 | 1e-3 | 18.74 | 0.459 | 0.782 | 16.224 | 12.734 | 0.000 | 0.000 | 3.413 | 0.016 | 0.000 | 0.735 | 9.25 |
|  | 1e-6 | 19.46 | 0.443 | 0.813 | 17.013 | 12.947 | 0.000 | 0.000 | 4.246 | 0.019 | 0.000 | 0.537 | 12.85 |
|  | 1e-8 | 16.27 | 0.446 | 0.873 | 13.775 | 9.208 | 0.000 | 0.020 | 4.584 | 0.031 | 0.000 | 0.685 | 24.05 |
| W2, kD = 500 | 1e-3 | 3.32 | 0.194 | 0.334 | 2.248 | 1.348 | 0.055 | 0.000 | 0.816 | 0.015 | 0.000 | 0.069 | 2.22 |
|  | 1e-6 | 4.01 | 0.185 | 0.346 | 3.043 | 1.896 | 0.000 | 0.014 | 1.087 | 0.025 | 0.000 | 0.095 | 3.72 |
|  | 1e-8 | 3.44 | 0.190 | 0.372 | 2.445 | 1.168 | 0.000 | 0.027 | 1.194 | 0.027 | 0.000 | 0.058 | 7.94 |
| W3 (N = 4031442), kD = 100 | 1e-3 | 1.25 | 0.207 | 0.268 | 0.148 | 0.003 | 0.095 | 0.000 | 0.047 | 0.001 | 0.000 | 0.057 | 2.34 |
|  | 1e-6 | 1.18 | 0.197 | 0.273 | 0.104 | 0.004 | 0.000 | 0.016 | 0.080 | 0.002 | 0.000 | 0.048 | 4.09 |
|  | 1e-8 | 1.19 | 0.190 | 0.272 | 0.111 | 0.007 | 0.000 | 0.033 | 0.067 | 0.002 | 0.000 | 0.052 | 10.38 |
| W3 (N = 1000000), kD = 200 | 1e-3 | 0.48 | 0.051 | 0.049 | 0.269 | 0.029 | 0.073 | 0.015 | 0.147 | 0.002 | 0.000 | 0.030 | 1.41 |
|  | 1e-6 | 0.51 | 0.049 | 0.054 | 0.298 | 0.055 | 0.000 | 0.009 | 0.218 | 0.003 | 0.000 | 0.026 | 1.83 |
|  | 1e-8 | 0.56 | 0.047 | 0.059 | 0.335 | 0.057 | 0.000 | 0.022 | 0.241 | 0.006 | 0.000 | 0.020 | 2.84 |
| W3 (N = 1000000), kD = 500 | 1e-3 | 2.26 | 0.059 | 0.400 | 1.665 | 0.760 | 0.000 | 0.000 | 0.828 | 0.013 | 0.000 | 0.126 | 3.38 |
|  | 1e-6 | 2.35 | 0.058 | 0.402 | 1.764 | 0.610 | 0.000 | 0.000 | 1.063 | 0.021 | 0.000 | 0.161 | 4.29 |
|  | 1e-8 | 2.39 | 0.058 | 0.392 | 1.784 | 0.535 | 0.000 | 0.055 | 1.141 | 0.025 | 0.000 | 0.155 | 7.63 |
