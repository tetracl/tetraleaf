# Helmholtz CH4.9a2b: P2M by the batched kernel

Task: `docs/helmholtz/tasks/CH4.9a2b.md` (roadmap CH4.9a2b, the second half of CH4.9a2; item 14 of
the priority list). It builds the Helmholtz FMM's P2M as one batched kernel per leaf on CH4.9a2a's
routines, runs it under the operator's `P2pChoice` with the per-point reference calls kept as the
trusted path, and measures it per point and per level on Kathleen and locust.

## 0. Summary

1. **The kernel** (§1): `tetraleaf_tables::helmholtz::leaf::LeafP2m` forms a leaf's multipoles for
   every charge vector in one call, 16 points at a time, on `LeafHarmonics` and `bessel_j_real`
   (real k; complex k on `bessel_j` per lane), each term summed over the lanes in a fixed pairwise
   order. Safe Rust, no `tetraleaf-simd` change, no dependency. The FMM runs it under every ISA and
   the reference's per-point calls under `P2pChoice::Reference` or `reference_leaves`.
2. **Accuracy, per leaf** (§2): every coefficient within its contract against
   `tetraleaf_ref::helmholtz::leaf::p2m`, C u √((2n + 1)/(4π)) Σⱼ |qⱼ eⱼ| (((n + 1)² + N) |ĵₙ| +
   |kρⱼ| |ĵₙ′|), with the measured C about **1.0** in f64 (p ≤ 45) and in f32 (p ≤ 32) against the
   test's 10; the kernel's multipole within 1.3–2.2·10⁻¹⁵ (f64) and 0.8–1.3·10⁻⁶ (f32) of the
   reference's largest value on the runs' leaves.
3. **Accuracy, FMM** (§3): the errors against the direct sum **0.993–1.000×** the reference P2M's on
   every case on both machines (the harness prints three digits, so a printed 1.000 means at most
   about 1.007) (the acceptance: within 1.01×), W1's and W2's constant density at
   kD = 500 included; in the tests 1.0000× on five scenarios, the outputs within 3.0·10⁻¹⁶.
4. **Time per run** (§4, one rank, every core): P2M **0.25–0.31× (Kathleen) and 0.15–0.27×
   (locust)** the reference's, from 3.3–8.5 % (Kathleen) and 3.2–7.2 % (locust) of an evaluation to
   1.0–2.4 % and 0.6–1.4 %; the evaluation 0.92–0.99× (Kathleen) and 0.94–0.97× (locust), several
   of those within the runs' spread (§4); the P2M ratios are the robust result.
5. **Time per point** (§5, one core): the kernel 0.20–0.24× (Kathleen) and 0.12–0.20× (locust) the
   reference's P2M in f64. Its batched sum is 27–56 % of it on Kathleen, where it is **no faster than
   the reference's scalar sum** (0.96–1.23×), and 16–47 % on locust (0.70–0.83×): the build targets
   baseline x86-64 and NEON, so the safe-Rust sum runs on 128-bit vectors and pays a lane reduction
   per term. That is the ISA question (§7).
6. **f32** (§6): the kernel's f32 time is 0.72× its f64 time at p = 45 on Kathleen (0.76× on
   locust), where the reference's is 1.40×; but the slowdown's pattern moved into the batched sum,
   which in f32 grows from 0.71× to 1.45× the scalar sum between p = 8 and 45 on Kathleen (not
   diagnosed).

## 1. Setup

- **Code**: branch `tbetcke/h-ch4-9a2b` at 1e9c996, release builds, BLAS threads 1.
- **The kernel** (`tetraleaf-tables/src/helmholtz/leaf.rs`, its module docs): per batch of 16 points
  (the last padded with the centre and zero charges) ρ = √(x² + y² + z²), ĵₙ(kρ) for n ≤ p, the factor
  e^{Im(k)ρ − β} (one value for real k), per vector and degree the weights (q e) ĵₙ, then per order m
  the batch's column Q̄ₙᵐ and (cos mφ, sin mφ) from `LeafHarmonics::columns`, and per term the sum
  over the lanes of wₙ Q̄ₙᵐ cos mφ (and sin mφ) per plane, lane j with j + 8, then j + 4, j + 2,
  j + 1, added to the coefficient. The products are the reference's in its order; ρ (`sqrt`, the
  reference `hypot`), the routines' values and the order of the sum over points differ.
- **The operator** (`tetraleaf-exec`, `helmholtz::operator` module docs, "P2M kernel"): the choice of
  P2P chooses P2M: `Densities::p2m_kernel` under every ISA (`Auto` included), `Densities::p2m` (one
  reference call per point) under `Reference`; `HelmholtzFmmBuilder::reference_leaves` (a measurement
  and test setting) gives the reference's P2M with the vectorised P2P. Yukawa keeps the reference.
- **Runs** (`tools/helmholtz/ch4-9a2b-runs.sh`; tables by `tools/helmholtz/ch4-9a2b-tables.py`):
  - `levels`: `helmholtz_accuracy --levels` at one rank with every core under the library's defaults
    named, two evaluations (the second's rows), §9.1's complex charges, on CH4.9's cases (W1 at
    kD = 250, 500, 1000 at 10⁻⁶, W1 at kD = 500 at 10⁻³ and 10⁻⁸, W3 at N = 10⁶ and kD = 500 at
    10⁻⁶), one P2M per launch (review 2026-11-24, B12), ABBA: kernel, reference, reference, kernel.
  - `worst`: W1's and W2's constant density at kD = 500, 10⁻⁶, kernel then reference (the density
    the record calls worst, CH4.3c3 and CH4.8).
  - `parts`: `helmholtz_leaf` on one core, 256 points as one leaf, p ∈ {8, 12, 17, 26, 32, 45} ×
    k̂ ∈ {0.5, 2, 4, 8, 16}, f64 and f32, one precision and degree per launch.
- **Machines**: Kathleen job 673212 (node-c11a-080, 1 × 40, `--map-by slot:PE=40 --bind-to core`;
  `helmholtz_leaf` pinned to core 2), 2026-11-27 16:45–17:01, every launch exit 0. locust (1 × 72,
  `helmholtz_leaf` pinned to core 2), 16:45–16:57 the same day, started by a gate once the one-minute
  load was below 2 with no other process of this account busy; a root `find` took 0.3–1 core
  throughout; every launch exit 0. Outputs in `~/projects/fmm-archive/{kathleen,locust}/runs/ch4-9a2b/`.
  The M6 not timed (a smoke run of the example only, its numbers not used).

## 2. Accuracy per leaf

`tetraleaf-tables/tests/helmholtz/leaf.rs` holds every coefficient of every vector and plane to

```text
|M̂ₙᵐ − M̂ₙᵐ_ref| ≤ C u √((2n + 1)/(4π)) Σⱼ |qⱼ eⱼ| (((n + 1)² + N) |ĵₙ(kρⱼ)| + |kρⱼ| |ĵₙ′(kρⱼ)|),  C = 10
```

on p ∈ {0, 1, 2, 5, 8, 17, 26, 45}, leaves of N = 1, 5, 16, 17, 40 and 200 points (the centre, the
poles, near-pole points, the axes and a corner among them), k̂ ∈ {0.3, 4, 16, 4 + 0.5i, 16 + 2i,
10⁻⁶}, real and complex charges, 1–3 vectors, in f64, and in f32 against the f64 reference on the
f32 inputs (p ≤ 32, no subnormal coefficient). The measured C is **0.999 in f64 and 1.005 in f32**.
The first contract, ((n + 1)² + N) u √((2n + 1)/(4π)) Σⱼ |wₙⱼ|, failed at p = 2, k̂ = 16 on one point
near a zero of ĵ₀: ρ by `sqrt` and by `hypot` differ in the last bit, which moves ĵₙ there by
|kρ| |ĵₙ′| u, far more than u |ĵₙ|; the conditioning term is that function's own sensitivity to ρ,
not the kernel's error. The tests also check vector v of three the one-vector call's bits, real
charges with real k the real plane of complex charges with Im q = 0 bit for bit, accumulation, the
scratch's degree and vectors, and the refusals (Yukawa, lengths).

On the runs' leaves (§5's column "worst check", the largest difference over the reference's largest
value, over every k̂): 1.3–2.2·10⁻¹⁵ in f64 and 7.8·10⁻⁷–1.3·10⁻⁶ in f32 on both machines.

## 3. Accuracy in the FMM

`tests/helmholtz_fmm.rs` (`check_leaf_kernel`, on the M6 at one rank): the kernel against the
reference's P2M on the same P2P kernel, outputs within 1.6–3.0·10⁻¹⁶ relative L2 and the errors
1.0000× on the cube at kD = 20 + 2i (ε = 10⁻³ and p = 8), the cube at kD = 1, the clusters at
kD = 20 (real charges) and the cube at kD = 20 + 2i with real charges at 10⁻⁷; every other check of
the executable (threads, pair by pair, vectors, one rank) now runs the kernel.

The runs, the errors' ratio kernel ÷ reference (each setting's first launch; the harness prints
three digits, so 0.998 and 0.993 are one unit of the last digit: 4.63 against 4.64, 1.45 against
1.46):

| case | Kathleen φ | Kathleen ∇φ | locust φ | locust ∇φ |
| --- | ---: | ---: | ---: | ---: |
| W1, kD = 250, 10⁻⁶ | 1.0000 | 1.0000 | 0.9982 | 0.9932 |
| W1, kD = 500, 10⁻⁶ | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| W1, kD = 1000, 10⁻⁶ | 1.0000 | 1.0000 | 0.9985 | 1.0000 |
| W1, kD = 500, 10⁻³ | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| W1, kD = 500, 10⁻⁸ | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| W3 (10⁶), kD = 500, 10⁻⁶ | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| W1, kD = 500, 10⁻⁶, constant | 1.0000 | 1.0000 | 0.9978 | 1.0000 |
| W2, kD = 500, 10⁻⁶, constant | 0.9978 | 1.0000 | 1.0000 | 1.0000 |

The two launches of each setting give one checksum on both machines (the kernel's bits do not
depend on the launch); the kernel's and the reference's checksums differ, as they must.

## 4. P2M per run and per level

P2M is the sum over the levels of "Translations and leaf expansions"' P2M (the slowest rank's level
calls, wall seconds), each setting the mean of its two launches; "spread" the larger ratio of a
setting's two evaluations.

Kathleen (1 × 40):

| case | P2M s, kernel | P2M s, reference | P2M ÷ | P2M of the evaluation, reference → kernel | evaluate s, kernel | evaluate s, reference | evaluate ÷ (spread) |
| --- | ---: | ---: | ---: | --- | ---: | ---: | --- |
| W1, kD = 250, 10⁻⁶ | 0.0145 | 0.0562 | 0.257 | 8.5 % → 2.4 % | 0.610 | 0.660 | 0.924 (1.033) |
| W1, kD = 500, 10⁻⁶ | 0.0527 | 0.1980 | 0.266 | 6.7 % → 1.8 % | 2.870 | 2.975 | 0.965 (1.031) |
| W1, kD = 1000, 10⁻⁶ | 0.1872 | 0.7442 | 0.252 | 5.0 % → 1.3 % | 14.640 | 14.760 | 0.992 (1.065) |
| W1, kD = 500, 10⁻³ | 0.0356 | 0.1336 | 0.266 | 6.4 % → 1.8 % | 1.950 | 2.080 | 0.938 (1.031) |
| W1, kD = 500, 10⁻⁸ | 0.1222 | 0.4735 | 0.258 | 8.5 % → 2.3 % | 5.275 | 5.555 | 0.950 (1.025) |
| W3 (10⁶), kD = 500, 10⁻⁶ | 0.0536 | 0.1727 | 0.310 | 3.3 % → 1.0 % | 5.200 | 5.285 | 0.984 (1.043) |

locust (1 × 72):

| case | P2M s, kernel | P2M s, reference | P2M ÷ | P2M of the evaluation, reference → kernel | evaluate s, kernel | evaluate s, reference | evaluate ÷ (spread) |
| --- | ---: | ---: | ---: | --- | ---: | ---: | --- |
| W1, kD = 250, 10⁻⁶ | 0.0036 | 0.0149 | 0.242 | 5.4 % → 1.4 % | 0.265 | 0.275 | 0.964 (1.038) |
| W1, kD = 500, 10⁻⁶ | 0.0124 | 0.0568 | 0.219 | 4.6 % → 1.0 % | 1.190 | 1.245 | 0.956 (1.008) |
| W1, kD = 1000, 10⁻⁶ | 0.0436 | 0.2200 | 0.198 | 3.8 % → 0.8 % | 5.675 | 5.830 | 0.973 (1.010) |
| W1, kD = 500, 10⁻³ | 0.0095 | 0.0352 | 0.269 | 4.5 % → 1.3 % | 0.760 | 0.780 | 0.974 (1.000) |
| W1, kD = 500, 10⁻⁸ | 0.0249 | 0.1668 | 0.149 | 7.2 % → 1.1 % | 2.170 | 2.310 | 0.939 (1.009) |
| W3 (10⁶), kD = 500, 10⁻⁶ | 0.0123 | 0.0626 | 0.197 | 3.2 % → 0.6 % | 1.920 | 1.985 | 0.967 (1.010) |

On Kathleen the evaluation's ratios at kD = 1000 (0.992, spread 1.065) and on W3 (0.984, 1.043)
lie within the run-to-run spread; the P2M ratios do not.

Per level with occupied leaves (levels whose reference P2M takes less than 1 ms, the level call's
overhead alone, left out), Kathleen: 0.245–0.274× (p = 12–32) but W3's level 5, 0.308×; locust: 0.14–0.26× (W1 at 10⁻⁸, p = 32: 0.14–0.16×; p = 8–26: 0.18–0.26×). The full
tables: `python3 tools/helmholtz/ch4-9a2b-tables.py levels <out-dir>`.

## 5. Time per point (one core, 16 lanes, f64; medians over k̂)

"Kernel sum" is the kernel less its routines' parts (the columns, the radial functions, the factor's
`exp`, each timed apart as in CH4.9a2a), "scalar sum" the reference's sum loop on the routines' own
values (CH4.9a2a's "sum, leaf values").

| p | Kathleen reference µs | kernel µs | ÷ reference | ÷ projected | kernel sum share | kernel sum ÷ scalar sum | locust reference µs | kernel µs | ÷ reference | ÷ projected | kernel sum share | kernel sum ÷ scalar sum |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 8 | 1.310 | 0.289 | 0.24 | 0.99 | 0.27 | 0.96 | 1.080 | 0.196 | 0.20 | 0.94 | 0.16 | 0.70 |
| 12 | 1.861 | 0.390 | 0.22 | 1.02 | 0.37 | 1.05 | 1.488 | 0.242 | 0.17 | 0.93 | 0.24 | 0.75 |
| 17 | 2.809 | 0.569 | 0.21 | 1.05 | 0.44 | 1.14 | 2.214 | 0.322 | 0.15 | 0.92 | 0.31 | 0.78 |
| 26 | 5.240 | 1.036 | 0.20 | 1.09 | 0.52 | 1.18 | 4.098 | 0.541 | 0.13 | 0.92 | 0.40 | 0.83 |
| 32 | 7.377 | 1.459 | 0.20 | 1.10 | 0.54 | 1.22 | 5.768 | 0.738 | 0.13 | 0.92 | 0.42 | 0.83 |
| 45 | 13.437 | 2.660 | 0.20 | 1.12 | 0.56 | 1.23 | 10.342 | 1.284 | 0.12 | 0.91 | 0.47 | 0.83 |

The per-point runs and the FMM's cases take real k only; the complex-k path (`bessel_j` and the
factor's `exp` per lane) is tested for accuracy (§2, §3) but not timed.

On locust the kernel is a little below CH4.9a2a's projection (0.91–0.94×); on Kathleen it reaches it
at p = 8–12 and lies above it from p = 17 (1.05–1.12×), because its batched sum is slower than the
scalar sum it replaces. The workspace builds for each architecture's baseline (no `target-cpu`):
SSE2 on x86-64 and NEON on aarch64, two f64 lanes each, without fma. The batched sum's per-lane work
is the scalar sum's plus a pairwise reduction of 16 lanes per term and plane; on two lanes that
reduction costs about what the vectorisation saves. The harmonics' columns and the real radial
recursion vectorise well on the same build (CH4.9a2a: 0.06–0.09× and 0.13–0.24× their references).

## 6. f32

| p | Kathleen reference f32 ÷ f64 | Kathleen kernel f32 µs | kernel f32 ÷ f64 | locust kernel f32 µs | kernel f32 ÷ f64 |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 8 | 0.94 | 0.272 | 0.94 | 0.162 | 0.83 |
| 17 | 0.98 | 0.419 | 0.74 | 0.245 | 0.76 |
| 32 | 1.06 | 1.009 | 0.69 | 0.546 | 0.74 |
| 45 | 1.40 | 1.904 | 0.72 | 0.969 | 0.75 |

CH4.9a1's f32 slowdown of the reference on Kathleen (1.40× f64 at p = 45 here) does not reach the
kernel. The f32 kernel's sum is still only 0.88× its f64 sum at p = 45 on Kathleen (1.31 against
1.50 µs) and 1.45× the scalar sum on the routines' f32 values: four f32 lanes in SSE2 should give
more. Subnormal products q e ĵₙ Q̄ₙᵐ (CH4.9a2a §4's remaining case) may contribute; not diagnosed.
Between p = 8 and 45 on Kathleen the f32 kernel over its projection rises from 0.91 to 1.24 and its
sum's share to 0.69: the subnormal pattern of CH4.9a1 now in the batched sum rather than gone. f32
is reserved for low degrees (CONVENTIONS §3.15) and the FMM runs f64 only.

## 7. Question for sign-off

Signed off 2026-11-27, the user on #269: "Q1 as recommended" (no ISA kernel for P2M now; the question
carried to CH4.9b).

1. **The ISA question** (the roadmap's: an ISA kernel only if the measured time stays far from the
   sum's). On Kathleen the batched sum is 27–56 % of the kernel and no faster than the scalar sum;
   the kernel is 1.05–1.12× the projection from p = 17. An ISA kernel of the sum (AVX-512 and AVX2 on
   x86 by runtime dispatch, as P2P, in `tetraleaf-simd`) would see 4–8 f64 lanes and fma there.
   But P2M is now 1.0–2.4 % of an evaluation on Kathleen and 0.6–1.4 % on locust, so the most it could
   save is about half of that, under 1.3 % of an evaluation.
   - **Recommended: no ISA kernel for P2M now.** Carry the question to CH4.9b, whose L2P is the larger
     kind (L2P 1.58 s of W1 at kD = 1000's 15.2 s on Kathleen, CH4.9) and sums over terms per target,
     with no reduction across lanes. If CH4.9b's safe-Rust L2P also stays far from its floor on x86,
     decide one ISA kernel for both kernels' sums then, as one follow-up.
   - Alternative: an ISA kernel of P2M's sum in `tetraleaf-simd` now, a new kernel in Laplace's crate
     (constraint 1: no Laplace path changes, its tests and benchmarks unchanged).
   - Alternative: lane accumulators per leaf in safe Rust (no reduction per term; (p + 1)² × 16 values
     per plane, vector and thread at the leaf degrees), still on two lanes on x86.

## 8. Reproduction

```sh
cargo build --release -p tetraleaf-validate --example helmholtz_accuracy --example helmholtz_leaf
tools/helmholtz/ch4-9a2b-runs.sh <out-dir> smoke levels worst parts
python3 tools/helmholtz/ch4-9a2b-tables.py runs|levels|worst|parts <out-dir>
```
