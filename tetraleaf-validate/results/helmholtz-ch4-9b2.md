# Helmholtz CH4.9b2: L2P with gradients by the batched kernel

Task: docs/helmholtz/tasks/CH4.9b2.md (roadmap CH4.9b2, item 14 of the priority list). The kernel is
`tetraleaf_tables::helmholtz::leaf::LeafL2p`; the operator runs it through `Densities::l2p_kernel`
under every ISA, and the per-target reference under `P2pChoice::Reference` or `reference_leaves`.
Raw outputs are in `~/projects/fmm-archive/{kathleen,locust}/runs/ch4-9b2/`. The tables come from
`tools/helmholtz/ch4-9b2-tables.py`.

## 0. Summary

- **Accuracy per target.** The kernel stays within its contract against
  `tetraleaf_ref::helmholtz::leaf::l2p` (§2) in f64 up to p = 45 and in f32 up to p = 32. The worst
  measured ratio is 1.9 (f64) and 1.3 (f32) against the tests' C = 10. The cases cover real, complex
  and tiny k, the centre's tangential limit, the z-axis rescale and targets whose r² underflows.
- **Accuracy in the FMM.** Errors are 0.993–1.000× the reference leaves' on CH4.6a2's cases and on
  W1's and W2's constant density at kD = 500, on both machines (the acceptance asks for at most
  1.01×). In the tests the outputs agree within 2.3·10⁻¹⁶ (1.0000×).
- **L2P per run.** L2P takes 0.19–0.23× the reference's time on Kathleen and 0.17–0.27× on locust.
- **The evaluation before and after CH4.9a–b** (both leaf kernels against both references):
  0.77–0.94× on Kathleen and 0.82–0.94× on locust. P2M + L2P fall from 11–30 % of an evaluation to
  2.9–7.7 % (Kathleen) and from 9–22 % to 2.0–4.4 % (locust).
- **Per point, one core, f64.** The kernel with gradients takes 0.15–0.24× the reference's L2P with
  gradients (Kathleen) and 0.14–0.22× (locust). That is 0.68–0.87× and 0.82–0.94× CH4.9b1's
  projection: the vectorised sum beats the restated scalar sum, at 0.56–0.74× and 0.70–0.84× its
  time.
- **The ISA question** (§6).
  - On the workspace's baseline build (SSE2 and NEON, two f64 lanes, no fma), the sum runs at
    4.9–8.0 GFLOP/s on Kathleen and 9.1–15.5 on locust. That is about 0.5–0.8 and 0.35–0.6 of each
    core's nominal baseline peak, so it is not far from its floor on either machine.
  - The same source built with `-C target-cpu=native` halves the sum on Kathleen (AVX-512:
    0.47–0.52×; the L2P kernel 0.56–0.67×, the P2M kernel 0.73–0.77×). On locust it gains 5–10 %.
  - An ISA kernel could therefore save 0.9–2.9 % of an evaluation on Kathleen and 0.1–0.4 % on
    locust. Recommended: no ISA kernel (question 1).
- **f32.** The L2P kernel is 0.55–0.78× its f64 time on both machines, with no slowdown; on Kathleen
  the f32 sum costs 0.54–0.59× the f64 sum per term. P2M's f32 kernel stays at 0.69–0.92× its f64
  time, as CH4.9a2b reported.

## 1. Setup

- **Kernel** (the module docs of `tetraleaf_tables::helmholtz::leaf`, "L2P").
  - Batches of 16 targets as lanes.
  - r is √(x² + y² + z²). The radial functions are `bessel_j_real(_derivatives)` for real k and
    `bessel_j` per lane for complex k.
  - The harmonics come from `LeafHarmonics::columns`, or `derivative_columns` with gradients.
  - Per order m and degree n, the per-degree sums take the terms of ±m together:
    A = c⁺ cos mφ + c⁻ sin mφ and B = c⁻ cos mφ − c⁺ sin mφ, with value += Q̄ A, ∂θ += (∂θ Q̄) A and
    (1/sin θ) ∂φ += m U B.
  - Then the reference's degree loop runs per lane: the tangential factor ĵₙ/r, with the centre's
    limit k/(3s) as a select. The Cartesian gradient uses `LeafDirections`, the harmonics' own frame
    (now public in `tetraleaf-math`, its values unchanged). Everything is multiplied by e and added
    into the target output chunk.
  - Safe Rust with no lane reduction; no `tetraleaf-simd` change.
  - The sum performs (p + 1)(12.5 p + 8) operations per target (`l2p_sum_flops`), about 12.5 per
    term.
- **Kathleen.**
  - Job 673257 on node-c11b-046 (not a slow node), `singlenode`, 2026-11-28 00:47–01:13.
  - Steps:
    1. The release build of cddcd8a.
    2. The same example built with `-C target-cpu=native` into a target directory of its own.
    3. `ch4-9b2-runs.sh smoke levels worst parts native`; every launch exited 0.
  - `helmholtz_accuracy` ran at 1 × 40 (`--map-by slot:PE=40 --bind-to core`).
  - `helmholtz_leaf` was pinned to core 2. Each one-core launch started within seconds of the last
    harness launch, so its one-minute load figure is the 40-thread average decaying (28.7 at the
    first f64 `parts` launch, 1.07 by the native ones).
- **locust.**
  - The same commit, 2026-11-28 00:49–01:08: the baseline build, the native build, then the same
    stages, every launch exiting 0.
  - A gate held the start until the one-minute load was below 2.
  - `helmholtz_accuracy` ran at 1 × 72; `helmholtz_leaf` was pinned to core 2.
  - The first f64 `parts` launches again start under the decaying average of the 72-thread launches
    (52 at p = 8, 2.1 by p = 45, 1.0 for the native ones). The P2M kernel's times there are the
    native launches' within 3 % and CH4.9a2b's within 3 %, so the core was quiet.
- **Harness settings.**
  - The defaults named (`--m2l auto-pw --diagonal auto --diagonal-build interpolated`, the V and W/X
    rules of the cost model), §9.1's complex charges and samples, gradients, and two evaluations
    with the rows taken from the second.
  - One leaf setting per launch, in the order kernel, reference, reference, kernel (ABBA; B12).
  - `--leaves kernel` is the library's default: P2M and L2P by the batched kernels, after CH4.9a–b.
    `--leaves reference` runs both by the per-point reference calls, before CH4.9a–b. P2P is the
    vectorised kernel in both.
- **Code state.** The self-review's fixes (cb02854) came after these launches. They change no
  kernel arithmetic: tests, docs, and the `∇ kernel sum` column of `helmholtz_leaf`, which
  subtracted the per-point `exp` time for real k although the kernel forms that factor once per
  call. `ch4-9b2-tables.py` recomputes the sum from the raw columns, so §5's and §6's sums are the
  corrected ones.

## 2. Accuracy per target

`tetraleaf-tables/tests/helmholtz/l2p.rs` holds each output per target to the contract below.
Notation: wₙ = √((2n + 1)/(4π)), Sₙ = Σₘ |L̂ₙᵐ|, x = kr, and tₙ = |ĵₙ|/r (or its limit at r = 0).

```text
|φ̃ − φ̃_ref|  ≤ C u e Σₙ wₙ Sₙ (((n + 1)² + 2n + 1) |ĵₙ| + |x| |ĵₙ′|)
|gₐ − gₐ_ref| ≤ C u e Σₙ wₙ Sₙ (|k| (((n + 1)² + 2n + 3) |ĵₙ′| + |x| |ĵₙ|) + ((n + 1)³ + n(n + 1)) tₙ)
```

The terms come from:

- the harmonics' (n + 1)² u and the derivatives' (n + 1)³ u (CH4.9a2a, CH4.9b1);
- the 2n + 1 terms of a degree summed in another order;
- r by `sqrt` against the reference's `hypot`, through ĵₙ′ and ĵₙ″ (bounded by the radial equation).

| test | cases | worst, in units of the contract (C = 10) |
| --- | --- | ---: |
| f64 | p ∈ {0, 1, 2, 5, 8, 17, 26, 45}; N = 1–40; k̂ 0.3, 4, 16, 4 + 0.5i, 16 + 2i, 10⁻⁶; with and without gradients; 1–3 vectors | 1.898 |
| f32, against the f64 reference on the f32 inputs | p ≤ 32, the same k̂ | 1.335 |
| r² underflowing (10⁻¹⁷⁰ to 10⁻²⁰⁰ from the centre) and the centre alone | p ≤ 17, every k̂ | 0.173 |

The targets include the centre, the poles, the axes, a corner, and three points next to the z-axis
(x² + y² below the normal range: two in f64, one in f32).

The bit identities are checked by `to_bits`:

- each target's outputs in leaves of 1, 16 and 40 and in any position (the centre, random targets,
  the three rescaled ones);
- vector v of three is the one-vector call's;
- the potentials with gradients are the potentials without them;
- a second call adds exactly.

`helmholtz_leaf`'s check before timing, the largest difference over the largest value, is at most
4.6·10⁻¹⁵ (f64) and 2.3·10⁻⁶ (f32) on both machines.

## 3. Accuracy in the FMM

In the tests, `check_leaf_kernel` runs every scenario and the cube at kD = 20 + 2i with real
charges. Against `reference_leaves` (now the reference's P2M and L2P), the outputs agree within
2.3·10⁻¹⁶ (φ) and 1.4·10⁻¹⁶ (∇φ), the errors are 1.0000× the reference's, and the one-rank FMM
holds on P ranks.

On the machines, the errors of the kernels' run over the reference's (each setting's first launch):

| case | Kathleen φ ÷ | Kathleen ∇φ ÷ | locust φ ÷ | locust ∇φ ÷ |
| --- | ---: | ---: | ---: | ---: |
| W1, kD = 250, 10⁻⁶ | 1.0000 | 1.0000 | 0.9982 | 0.9932 |
| W1, kD = 500, 10⁻⁶ | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| W1, kD = 1000, 10⁻⁶ | 1.0000 | 1.0000 | 0.9985 | 1.0000 |
| W1, kD = 500, 10⁻³ | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| W1, kD = 500, 10⁻⁸ | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| W3 (10⁶), kD = 500, 10⁻⁶ | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| W1, kD = 500, 10⁻⁶, constant density | 1.0000 | 1.0000 | 0.9978 | 1.0000 |
| W2, kD = 500, 10⁻⁶, constant density | 0.9978 | 1.0000 | 1.0000 | 1.0000 |

The harness prints three digits, so a ratio below 1 is a change in the third digit. The two
launches of each setting gave one checksum.

## 4. L2P and P2M per run, and the evaluation before and after CH4.9a–b

"Kernel" is the library's default; "reference" runs both leaf operators by the per-point calls. Each
time is the mean of the setting's two launches (ABBA), from the second evaluation. "Share" is
P2M + L2P over the evaluation. The spread is the larger ratio of a setting's two launches.

**Kathleen** (1 × 40)

| case | ε | L2P s, kernel | L2P s, reference | L2P ÷ | P2M ÷ | share, reference | share, kernel | evaluate s, kernel | evaluate s, reference | evaluate ÷ (spread) |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| W1, kD = 250 | 1e-6 | 0.0263 | 0.1303 | 0.202 | 0.236 | 0.263 | 0.072 | 0.565 | 0.730 | 0.774 (1.056) |
| W1, kD = 500 | 1e-6 | 0.0892 | 0.4327 | 0.206 | 0.239 | 0.205 | 0.053 | 2.560 | 3.080 | 0.831 (1.013) |
| W1, kD = 1000 | 1e-6 | 0.3574 | 1.7188 | 0.208 | 0.232 | 0.168 | 0.040 | 13.305 | 14.770 | 0.901 (1.040) |
| W1, kD = 500 | 1e-3 | 0.0596 | 0.2666 | 0.224 | 0.256 | 0.187 | 0.051 | 1.880 | 2.165 | 0.868 (1.042) |
| W1, kD = 500 | 1e-8 | 0.2315 | 1.2274 | 0.189 | 0.223 | 0.298 | 0.077 | 4.490 | 5.815 | 0.772 (1.009) |
| W3 (10⁶), kD = 500 | 1e-6 | 0.0986 | 0.4209 | 0.234 | 0.271 | 0.109 | 0.029 | 5.090 | 5.445 | 0.935 (1.013) |

**locust** (1 × 72)

| case | ε | L2P s, kernel | L2P s, reference | L2P ÷ | P2M ÷ | share, reference | share, kernel | evaluate s, kernel | evaluate s, reference | evaluate ÷ (spread) |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| W1, kD = 250 | 1e-6 | 0.0062 | 0.0261 | 0.236 | 0.233 | 0.151 | 0.040 | 0.240 | 0.270 | 0.889 (1.000) |
| W1, kD = 500 | 1e-6 | 0.0228 | 0.0988 | 0.231 | 0.199 | 0.127 | 0.031 | 1.085 | 1.210 | 0.897 (1.009) |
| W1, kD = 1000 | 1e-6 | 0.0887 | 0.3919 | 0.226 | 0.189 | 0.106 | 0.024 | 5.375 | 5.745 | 0.936 (1.028) |
| W1, kD = 500 | 1e-3 | 0.0135 | 0.0490 | 0.274 | 0.235 | 0.109 | 0.031 | 0.705 | 0.765 | 0.922 (1.014) |
| W1, kD = 500 | 1e-8 | 0.0561 | 0.3407 | 0.165 | 0.150 | 0.223 | 0.044 | 1.850 | 2.270 | 0.815 (1.000) |
| W3 (10⁶), kD = 500 | 1e-6 | 0.0241 | 0.1194 | 0.201 | 0.183 | 0.092 | 0.020 | 1.800 | 1.960 | 0.918 (1.011) |

**Per level** (levels with occupied leaves; L2P and P2M ÷ the reference's)

| case | ε | Kathleen level (p): L2P ÷, P2M ÷ | locust level (p): L2P ÷, P2M ÷ |
| --- | --- | --- | --- |
| W1, kD = 250 | 1e-6 | 4 (26): 0.213, 0.262; 5 (17): 0.199, 0.224 | 5 (17): 0.208, 0.192; 6 (15): 0.246, 0.221 |
| W1, kD = 500 | 1e-6 | 5 (26): 0.199, 0.242; 6 (17): 0.207, 0.235 | 6 (17): 0.214, 0.183; 7 (15): 0.235, 0.192 |
| W1, kD = 1000 | 1e-6 | 6 (26): 0.200, 0.258; 7 (17): 0.209, 0.227 | 7 (17): 0.205, 0.174; 8 (15): 0.233, 0.186 |
| W1, kD = 500 | 1e-3 | 5 (19): 0.229, 0.229; 6 (12): 0.222, 0.255 | 6 (12): 0.241, 0.211; 7 (8): 0.280, 0.218 |
| W1, kD = 500 | 1e-8 | 5 (32): 0.192, 0.220; 6 (32): 0.188, 0.222 | 5 (32): 0.168, 0.147; 6 (32): 0.164, 0.147 |
| W3 (10⁶), kD = 500 | 1e-6 | 5 (26): 0.234, 0.268 | 5 (26): 0.201, 0.175 |

The leaf levels differ between the machines because the default leaf size is the ISA's (256 points
on x86, 128 on NEON at ε ≥ 10⁻⁷; CH4.0d3).

## 5. Time per point (one core, 16 lanes; medians over k̂ ∈ {0.5, 2, 4, 8, 16})

Column notes:

- "∇ kernel" is `LeafL2p` with gradients over 256 targets as one leaf.
- "÷ projected" compares against CH4.9b1's projection: the derivative columns, the radial functions
  with derivatives, the factor's `exp` and the restated scalar sum.
- "sum" is the kernel less the derivative columns and the real radial functions with derivatives.
- "scalar sum" is the reference's sum restated on the routines' values (CH4.9b1's "sum ∇, leaf
  values").

**Kathleen** (the SSE2 baseline build)

| precision | p | reference L2P ∇ µs | ∇ kernel µs | ÷ reference | ÷ projected | sum µs | sum share | sum ns per term | sum GFLOP/s | sum ÷ scalar sum | L2P (no ∇) ÷ reference |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| f64 | 8 | 2.108 | 0.481 | 0.239 | 0.87 | 0.197 | 0.41 | 2.432 | 4.93 | 0.74 | 0.233 |
| f64 | 12 | 3.511 | 0.712 | 0.208 | 0.80 | 0.346 | 0.49 | 2.047 | 5.94 | 0.68 | 0.201 |
| f64 | 17 | 6.109 | 1.117 | 0.186 | 0.75 | 0.590 | 0.53 | 1.821 | 6.73 | 0.62 | 0.182 |
| f64 | 26 | 12.507 | 2.107 | 0.169 | 0.70 | 1.211 | 0.57 | 1.661 | 7.42 | 0.58 | 0.164 |
| f64 | 32 | 18.788 | 2.986 | 0.160 | 0.69 | 1.748 | 0.59 | 1.605 | 7.70 | 0.57 | 0.159 |
| f64 | 45 | 35.593 | 5.427 | 0.152 | 0.68 | 3.264 | 0.60 | 1.543 | 8.04 | 0.56 | 0.153 |
| f32 | 8 | 2.062 | 0.375 | 0.170 | 0.71 | 0.117 | 0.31 | 1.444 | 8.31 | 0.43 | 0.190 |
| f32 | 17 | 5.974 | 0.707 | 0.117 | 0.52 | 0.325 | 0.46 | 1.003 | 12.21 | 0.34 | 0.123 |
| f32 | 32 | 18.559 | 1.730 | 0.093 | 0.44 | 0.951 | 0.55 | 0.873 | 14.16 | 0.31 | 0.094 |
| f32 | 45 | 37.033 | 3.096 | 0.087 | 0.43 | 1.763 | 0.57 | 0.833 | 14.89 | 0.30 | 0.085 |

**locust** (Neoverse-V2, NEON baseline)

| precision | p | reference L2P ∇ µs | ∇ kernel µs | ÷ reference | ÷ projected | sum µs | sum share | sum ns per term | sum GFLOP/s | sum ÷ scalar sum | L2P (no ∇) ÷ reference |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| f64 | 8 | 1.515 | 0.320 | 0.219 | 0.94 | 0.107 | 0.33 | 1.321 | 9.08 | 0.84 | 0.199 |
| f64 | 12 | 2.430 | 0.447 | 0.192 | 0.89 | 0.176 | 0.39 | 1.041 | 11.67 | 0.76 | 0.170 |
| f64 | 17 | 4.038 | 0.670 | 0.169 | 0.86 | 0.302 | 0.45 | 0.932 | 13.14 | 0.74 | 0.151 |
| f64 | 26 | 8.275 | 1.230 | 0.149 | 0.83 | 0.610 | 0.50 | 0.837 | 14.74 | 0.71 | 0.126 |
| f64 | 32 | 12.020 | 1.727 | 0.142 | 0.82 | 0.886 | 0.51 | 0.814 | 15.20 | 0.70 | 0.117 |
| f64 | 45 | 23.133 | 3.114 | 0.135 | 0.82 | 1.694 | 0.54 | 0.801 | 15.49 | 0.71 | 0.108 |
| f32 | 8 | 1.376 | 0.228 | 0.171 | 0.78 | 0.063 | 0.28 | 0.778 | 15.43 | 0.50 | 0.161 |
| f32 | 17 | 3.452 | 0.405 | 0.119 | 0.61 | 0.145 | 0.36 | 0.448 | 27.37 | 0.36 | 0.116 |
| f32 | 32 | 9.975 | 0.968 | 0.098 | 0.54 | 0.396 | 0.41 | 0.364 | 34.00 | 0.33 | 0.091 |
| f32 | 45 | 18.884 | 1.723 | 0.092 | 0.53 | 0.729 | 0.42 | 0.345 | 36.00 | 0.32 | 0.084 |

The f32 rows at p = 12 and 26 are in the archive (`ch4-9b2-tables.py parts`); they lie between their
neighbours. Only real k is timed per point; the complex-k path (`bessel_j` and `exp` per lane) is
tested for accuracy (§2, §3) and runs in none of the FMM's cases.

After the kernel the sum is still 41–60 % of a target's L2P with gradients on Kathleen and 33–54 %
on locust (CH4.9b1's restated sum was 49–73 % and 37–62 % of the projection).

## 6. The ISA question

The question, carried from CH4.9a2b (its question 1, signed off 2026-11-27): measure L2P's sum
against its floor on both machines, and if it stays far from that floor on x86, decide on one ISA
kernel for both sums.

- **The floor.** The workspace builds for each architecture's baseline, with no `target-cpu`:
  - Kathleen (SSE2): two FP ports × two f64 lanes = 4 operations per cycle without fma, about
    10 GFLOP/s at the nominal 2.5 GHz and 15.6 at the single-core turbo of 3.9 GHz;
  - locust (NEON): four 128-bit FP pipes = 8 per cycle, about 25–27 GFLOP/s at 3.1–3.4 GHz.

  The sum reaches 4.9–8.0 GFLOP/s on Kathleen (0.5–0.8 of the nominal figure, 0.3–0.5 of the turbo
  one) and 9.1–15.5 on locust (0.35–0.6). It stays within about 2–3× of the baseline's peak on both
  machines, so it is not far from its floor on either.
- **Wider vectors.** The same source built with `-C target-cpu=native` (stage `native`; f64 medians
  over k̂):

| p | Kathleen ∇ kernel ÷ baseline | sum ÷ | sum GFLOP/s, baseline → native | P2M kernel ÷ | locust ∇ kernel ÷ | sum ÷ | P2M kernel ÷ |
| ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: |
| 8 | 0.67 | 0.52 | 4.93 → 9.44 | 0.77 | 0.99 | 0.95 | 1.03 |
| 12 | 0.64 | 0.51 | 5.94 → 11.67 | 0.74 | 0.98 | 0.96 | 1.02 |
| 17 | 0.61 | 0.48 | 6.73 → 14.12 | 0.73 | 0.97 | 0.93 | 1.04 |
| 26 | 0.58 | 0.47 | 7.42 → 15.75 | 0.74 | 0.95 | 0.89 | 1.04 |
| 32 | 0.58 | 0.51 | 7.70 → 15.13 | 0.74 | 0.95 | 0.90 | 1.02 |
| 45 | 0.56 | 0.48 | 8.04 → 16.68 | 0.74 | 0.95 | 0.90 | 1.02 |

  - On Kathleen, AVX-512's eight lanes halve the safe-Rust sum without any change to the code. An
    ISA kernel in `tetraleaf-simd` (runtime dispatch, as P2P) would reach about that, plus fma.
  - On locust the native build adds SVE2 at the same 128-bit width, and gains 5–10 %.
- **What it would save.**
  - Taking the native ratios as what an ISA kernel would reach, it would remove 33–44 % of L2P and
    23–27 % of P2M on Kathleen, where the two together are now 2.9–7.7 % of an evaluation (§4).
  - That is 0.9–2.9 % of an evaluation on CH4.6a2's cases: W1 at 10⁻⁸ and at kD = 250 the most, W3
    the least.
  - On locust (5–10 % of 2.0–4.4 %) it is 0.1–0.4 %.
- **Recommendation: no ISA kernel** (question 1). Both kernels' sums stay in safe Rust. A user on
  x86 who builds for the machine (`-C target-cpu=native`) gets the halving at no cost in code. The
  ISA question closes here, as CH4.9a2b's sign-off carried it. The leaf expansions are no longer
  where an evaluation's time goes: the M2L and P2P take the rest.

## 7. f32

| p | Kathleen L2P ∇ kernel f32 ÷ f64 | L2P f32 ÷ f64 | P2M f32 ÷ f64 | locust L2P ∇ kernel f32 ÷ f64 | L2P f32 ÷ f64 | P2M f32 ÷ f64 |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 8 | 0.78 | 0.92 | 0.92 | 0.71 | 0.78 | 0.82 |
| 17 | 0.63 | 0.68 | 0.72 | 0.60 | 0.67 | 0.78 |
| 32 | 0.58 | 0.58 | 0.69 | 0.56 | 0.60 | 0.74 |
| 45 | 0.57 | 0.56 | 0.71 | 0.55 | 0.60 | 0.75 |

- **L2P.** The kernel's f32 sum has no lane reduction and gains from SSE2's four f32 lanes: on
  Kathleen 0.54–0.59× the f64 sum per term at the same operation count. CH4.9a1's f32 slowdown does
  not reach it.
- **P2M.** The P2M kernel's f32 figure is CH4.9a2b's again: 0.71× at p = 45 on Kathleen. Its sum
  pays a lane reduction per term, and its subnormal products are still not diagnosed.
- f32 is reserved for low degrees (CONVENTIONS §3.15), and the FMM runs f64 only. Nothing more is
  proposed for it here.

## 8. Question for sign-off

Signed off 2026-11-28 by the user on #272: "Q1 as recommended" (no ISA kernel; the question closes).

1. **The ISA question** (§6).
   - **Recommended: no ISA kernel for the leaf sums.** On the baseline build both sums stay within
     about 2–3× of the baseline's peak. An AVX-512 kernel would save 0.9–2.9 % of an evaluation on
     Kathleen and 0.1–0.4 % on locust. The question closes; the report records that
     `-C target-cpu=native` halves the L2P sum on x86.
   - Alternative: a follow-up task off the priority list (CH4.9d), the two leaf sums as
     `tetraleaf-simd` kernels with AVX-512 and AVX2 dispatch.
     - Constraint 1: a new kernel in Laplace's crate; its tests and benchmarks must stay unchanged.
     - It needs `unsafe` in `tetraleaf-simd`'s `arch` modules, where the root rules allow it.
   - Alternative: that kernel now, in this task.

## 9. Reproduction

```sh
cargo build --release -p tetraleaf-validate --example helmholtz_accuracy --example helmholtz_leaf
RUSTFLAGS='-C target-cpu=native' CARGO_TARGET_DIR=target/native \
    cargo build --release -p tetraleaf-validate --example helmholtz_leaf
NATIVE_DIR=target/native tools/helmholtz/ch4-9b2-runs.sh <out> smoke levels worst parts native
python3 tools/helmholtz/ch4-9b2-tables.py runs|levels|worst|parts|native <out>
```

On Kathleen it ran as `tools/kathleen/submit.sh --exclude=node-c11b-004,node-c11b-013
--time=02:30:00 -- run.sbatch '<build>' '<native build>' '<runs>'` (job 673257). On locust it ran
through `tools/gh200/remote.sh` after `tools/gh200/sync.sh`.
