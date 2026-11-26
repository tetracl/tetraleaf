# Helmholtz CH4.6a2: the φ step's time per level, the product against rustfft

Task: `docs/helmholtz/tasks/CH4.6a2.md` (roadmap CH4.6a2, item 13 of the priority list; proposal
§10 CH4.6a, its last two clauses; §6.7 "Grid sizes" and "Host back ends"; §12 D7). CH4.6a1 put the
φ step behind `HelmholtzPhiStep` with rustfft behind the feature `rustfft`, on the product's grid.
This task measures what rustfft gains per diagonal level, prices the smooth grid, and asks whether
CH4.6a3 (the smooth grid) is built.

## 0. Summary

1. **The harness** (§1): `helmholtz_accuracy --phi-step product|rustfft` (`Case::phi_step`), the
   back end in the output ("φ step: rustfft; rustfft 6.4.1, x86_64 features sse4.1, avx, fma"), the step
   each level ran (`Row::phi_steps`, a column "φ" of "Diagonal M2L"), and a feature `rustfft` of
   tetraleaf-validate that passes tetraleaf-exec's through.
2. **The conversions** (§2, §6): rustfft's S → F takes 0.27–0.46× the product's thread time per
   conversion on Kathleen (1 × 40) and 0.31–0.62× on locust (1 × 72), its F → R 0.23–0.60× and
   0.33–0.71×; over each case's conversions 0.30–0.39× and 0.41–0.47×.
3. **The levels and the run** (§2, §3, §6): the diagonal levels' M2L 0.63–0.93× (Kathleen) and
   0.60–0.95× (locust) the product's per level, 0.72–0.90× and 0.77–0.91× summed per case; the
   evaluation 0.89–0.95× and 0.91–0.99× (the product's own spread between its two launches
   1.000–1.068 and 1.000–1.036: Kathleen's W3, 0.936, is within its spread of 1.068, every other
   case's gain beyond its spread); build plus one evaluation 0.94–0.97× on Kathleen. The errors are
   0.993–1.000× the product's, the checksums of each setting's two launches agree, the build is
   unchanged. CH4.1a had estimated the FFT's gain at about 1–2 % of build plus one evaluation (W1,
   kD = 500); measured on Kathleen it is 3–6 %.
4. **The smooth grid** (§4, §6): priced from S9's per-box rustfft times at the levels' degrees on
   the node and the level's measured pair cost at the larger N_φ, it changes a level by 0.98–1.03×
   and the evaluation by 1.000–1.004× on Kathleen; with S9's FFT difference scaled up to the level's
   own rustfft conversion cost (an upper bound on the FFT's cost in the level), by 0.90–1.03× and
   0.991–1.001×. On locust 1.000–1.072× per level and 1.000–1.005× per evaluation: on NEON rustfft
   is often slower on the smooth grid. At 5 of the 15 degrees timed the product's grid is already
   7-smooth (N_φ = 72, 80, 108, 168, 288; also 36 and 64 at locust's p = 8 and 15), and where it is
   not the FFT saved is mostly less than the pairs' extra samples cost; on Kathleen the smooth grid
   gains at p = 61 (0.979×, up to 0.90× at the bound) and at p = 45 and 241 at the bound. Either way
   about 1 % of an evaluation or less, against a second set of tables per level. **Recommendation
   (question 1): keep the product's grid and close CH4.6a3 unrun.**
5. **The time rule** (§5, §6): on the 6 ABBA cases every level in the rule's range that runs the
   form keeps it. On the levels kept on the waves (7 more cases where CH4.3b2c found the waves kept,
   W1 at kD = 100, 200, W3 at kD = 200, 500, 10⁻⁶–10⁻⁸), the model with the form's conversions at
   rustfft's measured ratio moves one level on Kathleen (W3, kD = 200, 10⁻⁶, level 4: the forced
   form with rustfft 0.98× the same setting's waves, within noise) and two on locust (one a 3 %
   loss; one a 0.024 s gain where the form beats the waves with the product too). Where the form
   beats the waves by more (Kathleen 1 level, locust 3), it does so under either step: misses of the
   waves' model, CH4.3b2c's kind, not the φ step's (a CH4.7 *Also*). **Recommendation (question 2):
   no FFT variant of `DiagonalCosts`** (rustfft is opt-in, and the moves it would make are within
   noise or a mix of gain and loss).
6. **The M6** (§7): skipped at the user's word before any run (2026-11-26: "Skip the M6"; question
   3).

## 1. Setup

- **Code**: branch `tbetcke/h-ch4-6a2` at 68827ef (the harness; the later commits change only
  scripts and documents), `helmholtz_accuracy` built in release with `--features rustfft`, BLAS
  threads 1.
- **Runs** (`tools/helmholtz/ch4-6a2-runs.sh`; tables `tools/helmholtz/ch4-6a2-tables.py`): one
  rank with every core, §9.1's complex charges and 1,000 samples, gradients, two evaluations (the
  rows the second's), the library's defaults named (`--m2l auto-pw --diagonal auto
  --diagonal-build interpolated --near-v cost --near-wx cost`), `--levels`, **one φ step per
  launch** (review 2026-11-24, B12).
  - `levels`: W1 at kD = 250, 500, 1000 at 10⁻⁶, W1 at kD = 500 at 10⁻³ and 10⁻⁸, W3 at N = 10⁶,
    kD = 500 at 10⁻⁶; per case product, rustfft, rustfft, product (ABBA), then
    `--diagonal-costs form` with each step.
  - `waves`: W1 at kD = 100, 200 at 10⁻⁷, 10⁻⁸, W3 at N = 10⁶, kD = 200 at 10⁻⁶, 10⁻⁸ and kD = 500
    at 10⁻⁸; per case the defaults with each step, then `--diagonal-costs form` with each (one
    launch per setting).
  - `ffts`: spike S9's `times` (`spikes/fft-plans`), one thread, `rustfft product`, `rustfft
    smooth` and `product product` at the diagonal levels' degrees p = 12–241.
- **locust**: §6 (the same stages, 1 × 72, 2026-11-26 19:35–19:44, the harness built at 68827ef in
  the tree synced at ccbb9b7).
- **Kathleen**: node-c11b-009 (Xeon Gold 6248, AVX-512; rustfft dispatched on sse4.1, avx, fma),
  `singlenode`, slow nodes excluded; job 673075 (`smoke`, `levels`, 18:14–18:27), 673077 (`ffts`),
  673078 (`waves`), 2026-11-26. The node is exclusive; its load at a launch's start is the decay of
  the job's own build and launches.
- **Raw outputs**: `~/projects/fmm-archive/kathleen/runs/ch4-6a2/` (one `.md` and `.err` per launch,
  `runs.log`).

## 2. The conversions and the levels (Kathleen, 1 × 40)

Thread-µs per conversion (S → F per source signature, F → R per accumulator converted back) and
the level's M2L seconds, each the mean of the setting's two launches (`ch4-6a2-tables.py levels`).
P the product, R rustfft; every level listed ran rustfft under R (the ranks' agreement put none on
the product).

| case | level | p | conversions S → F, F → R | S → F µs P → R | ÷ | F → R µs P → R | ÷ | M2L s P → R | ÷ |
| --- | ---: | ---: | --- | --- | ---: | --- | ---: | --- | ---: |
| w1-250-1e-6 | 2 | 71 | 56, 56 | 2549.3 → 810.4 | 0.318 | 3190.7 → 978.3 | 0.307 | 0.0458 → 0.0290 | 0.634 |
| w1-250-1e-6 | 3 | 41 | 272, 272 | 1278.7 → 355.6 | 0.278 | 1282.2 → 305.8 | 0.238 | 0.0468 → 0.0323 | 0.689 |
| w1-250-1e-6 | 4 | 26 | 1160, 1160 | 419.3 → 112.5 | 0.268 | 394.3 → 123.2 | 0.312 | 0.0612 → 0.0445 | 0.728 |
| w1-250-1e-6 | 5 | 17 | 3724, 4081 | 154.4 → 49.5 | 0.321 | 140.6 → 51.2 | 0.364 | 0.0598 → 0.0471 | 0.788 |
| w1-500-1e-3 | 2 | 116 | 56, 56 | 8798.9 → 2985.5 | 0.339 | 11842.5 → 5691.5 | 0.481 | 0.1864 → 0.1540 | 0.826 |
| w1-500-1e-3 | 3 | 61 | 272, 272 | 2396.4 → 1096.3 | 0.457 | 2243.0 → 1108.2 | 0.494 | 0.1252 → 0.1071 | 0.855 |
| w1-500-1e-3 | 4 | 33 | 1160, 1160 | 586.4 → 213.3 | 0.364 | 549.2 → 223.0 | 0.406 | 0.1038 → 0.0891 | 0.859 |
| w1-500-1e-3 | 5 | 19 | 4760, 4760 | 171.7 → 62.6 | 0.364 | 175.1 → 59.3 | 0.339 | 0.1157 → 0.0952 | 0.822 |
| w1-500-1e-3 | 6 | 12 | 15177, 16333 | 75.4 → 25.7 | 0.340 | 60.6 → 20.8 | 0.344 | 0.1067 → 0.0723 | 0.677 |
| w1-500-1e-6 | 2 | 129 | 56, 56 | 9683.8 → 3752.6 | 0.388 | 15600.3 → 6575.9 | 0.422 | 0.2484 → 0.1612 | 0.649 |
| w1-500-1e-6 | 3 | 71 | 272, 272 | 3366.7 → 1218.7 | 0.362 | 3511.3 → 1241.0 | 0.353 | 0.1954 → 0.1597 | 0.817 |
| w1-500-1e-6 | 4 | 41 | 1160, 1160 | 957.2 → 353.0 | 0.369 | 913.0 → 286.9 | 0.314 | 0.1818 → 0.1577 | 0.867 |
| w1-500-1e-6 | 5 | 26 | 4760, 4760 | 325.6 → 117.8 | 0.362 | 300.2 → 100.0 | 0.333 | 0.2589 → 0.2252 | 0.870 |
| w1-500-1e-6 | 6 | 17 | 14825, 16333 | 130.6 → 51.8 | 0.396 | 122.2 → 41.3 | 0.338 | 0.2247 → 0.1696 | 0.755 |
| w1-500-1e-8 | 2 | 135 | 56, 56 | 12165.3 → 3810.7 | 0.313 | 20030.0 → 9502.0 | 0.474 | 0.3270 → 0.2156 | 0.659 |
| w1-500-1e-8 | 3 | 76 | 272, 272 | 3871.5 → 1425.9 | 0.368 | 4307.0 → 1851.7 | 0.430 | 0.2460 → 0.2067 | 0.840 |
| w1-500-1e-8 | 4 | 45 | 1160, 1160 | 1097.1 → 482.2 | 0.439 | 1118.3 → 458.0 | 0.410 | 0.2396 → 0.2091 | 0.873 |
| w1-500-1e-8 | 5 | 32 | 4760, 4760 | 520.8 → 183.7 | 0.353 | 493.3 → 143.5 | 0.291 | 0.7306 → 0.6476 | 0.886 |
| w1-1000-1e-6 | 2 | 241 | 56, 56 | 55369.4 → 16282.1 | 0.294 | 94786.2 → 56390.5 | 0.595 | 1.0316 → 0.8168 | 0.792 |
| w1-1000-1e-6 | 3 | 129 | 272, 272 | 13764.7 → 4820.1 | 0.350 | 23178.1 → 13185.2 | 0.569 | 1.1621 → 0.9385 | 0.808 |
| w1-1000-1e-6 | 4 | 71 | 1160, 1160 | 3289.2 → 1323.4 | 0.402 | 3854.4 → 1448.8 | 0.376 | 0.7741 → 0.7198 | 0.930 |
| w1-1000-1e-6 | 5 | 41 | 4760, 4760 | 893.1 → 325.4 | 0.364 | 933.6 → 286.1 | 0.306 | 0.7681 → 0.6300 | 0.820 |
| w1-1000-1e-6 | 6 | 26 | 19138, 19138 | 319.2 → 96.7 | 0.303 | 310.4 → 96.1 | 0.310 | 1.0841 → 0.8412 | 0.776 |
| w1-1000-1e-6 | 7 | 17 | 59526, 64756 | 136.5 → 41.2 | 0.302 | 122.0 → 42.3 | 0.347 | 0.9224 → 0.6533 | 0.708 |
| w3-500-1e-6 | 2 | 129 | 64, 64 | 9999.3 → 3494.5 | 0.349 | 13348.6 → 5844.1 | 0.438 | 0.2688 → 0.2183 | 0.812 |
| w3-500-1e-6 | 3 | 71 | 512, 512 | 3243.9 → 1217.6 | 0.375 | 3915.5 → 1050.6 | 0.268 | 0.6740 → 0.5933 | 0.880 |
| w3-500-1e-6 | 4 | 41 | 4096, 4096 | 881.7 → 308.9 | 0.350 | 1112.0 → 260.5 | 0.234 | 1.5235 → 1.4053 | 0.922 |

Read per row: S → F gains more than F → R on the coarsest levels (level 2: 56–64 boxes on 40
threads, p ≥ 116; and W1 at kD = 1000's level 3, p = 129), where F → R's thread time per
conversion stays 1.7–3.5× S → F's under rustfft; at the finer levels both are near a third of the
product's. The pairs' thread time per pair (not shown) moves 0.78–1.20× between the settings
without a trend (the largest moves on W1 at kD = 250's level 5, 1.20×, and kD = 1000's level 4,
1.15×); the pair does not depend on the φ step. The levels gain least where the conversions are a
small part of them (W1 at kD = 1000's level 4, 0.93×; W3's dense level 4, 0.92×).

## 3. The runs (Kathleen)

The evaluation (the second of each launch) and the build, per launch in the order run (P1, R1,
R2, P2); the ratios are of the means (`ch4-6a2-tables.py runs`).

| case | build s P1, R1, R2, P2 | evaluate s P1, R1, R2, P2 | evaluate R ÷ P | P spread | diagonal M2L s P, R | ÷ | errors R ÷ P (φ, ∇φ) | checksums |
| --- | --- | --- | ---: | ---: | --- | ---: | --- | --- |
| w1-250-1e-6 | 0.62, 0.62, 0.63, 0.63 | 0.700, 0.610, 0.630, 0.690 | 0.892 | 1.014 | 0.214, 0.153 | 0.716 | 1.0000, 1.0000 | agree |
| w1-500-1e-3 | 2.47, 2.54, 2.47, 2.53 | 2.080, 1.990, 1.980, 2.140 | 0.941 | 1.029 | 0.638, 0.518 | 0.812 | 1.0000, 1.0000 | agree |
| w1-500-1e-6 | 2.75, 2.75, 2.70, 2.65 | 3.040, 2.850, 2.880, 3.050 | 0.941 | 1.003 | 1.109, 0.873 | 0.787 | 1.0000, 1.0000 | agree |
| w1-500-1e-8 | 2.91, 2.87, 2.90, 2.91 | 5.810, 5.510, 5.470, 5.810 | 0.945 | 1.000 | 1.543, 1.279 | 0.829 | 1.0000, 1.0000 | agree |
| w1-1000-1e-6 | 13.45, 13.58, 13.53, 13.51 | 14.950, 13.750, 14.240, 15.110 | 0.931 | 1.011 | 5.742, 4.599 | 0.801 | 1.0000, 1.0000 | agree |
| w3-500-1e-6 | 3.11, 3.14, 3.15, 3.13 | 5.480, 4.620, 5.310, 5.130 | 0.936 | 1.068 | 2.466, 2.217 | 0.899 | 1.0000, 1.0000 | agree |

- Each rustfft launch is faster than both product launches of its case, except W3 at kD = 500,
  whose spread (1.068) covers one rustfft launch (5.31 s against 5.13 s); its mean ratio, 0.936,
  is within that spread.
- The build is the same (rustfft's plans are made in the conversions' build, a few ms); build plus
  one evaluation 0.94× (kD = 250) and 0.96–0.97× (kD ≥ 500).
- The errors differ from the product's in the fifth significant digit at most (CH4.6a1: the form's
  own conditioning of a few u per conversion).
- **Against CH4.1a's estimate.** CH4.1a modelled the FFT's whole gain (the φ step free) at about
  1–2 % of build plus one evaluation on W1 at kD = 500, from one-core times per box scaled to
  CH4.0e's level times. Measured with rustfft, which leaves 0.23–0.60× of the conversions, the gain
  is 3–6 % of build plus one evaluation and 5–11 % of the evaluation. The product's conversion in
  the level call costs about its one-core time (p = 129: 25.3 ms of thread time for S → F + F → R,
  against 20.1–24.7 ms per box in S9's one-core times at p = 123–135 on Kathleen), so the
  difference lies in the model's share of the conversions in the level, which this task did not
  analyse further.

## 4. The smooth grid, priced (Kathleen)

Per level, rustfft's thread-seconds as measured on the product's grid, with S → F and F → R moved
by the conversions times S9's per-box rustfft difference at the level's degree (smooth less product
grid, `ffts-rustfft-*.md`, one thread on the same node) and the pairs scaled by N_φ smooth / N_φ;
the level's M2L seconds scaled by the thread-seconds' ratio (`ch4-6a2-tables.py smooth`).

| p | N_φ product, smooth | rustfft per box S → F + F → R µs, product → smooth | levels (cases) | smooth ÷ measured |
| ---: | --- | --- | --- | --- |
| 12 | 52, 56 | 10.5 → 6.1 | w1-500-1e-3 L6 | 1.031 |
| 17 | 72, 72 | 12.2 → 12.5 | 3 levels | 1.001 |
| 19 | 80, 80 | 16.8 → 16.6 | w1-500-1e-3 L5 | 1.000 |
| 26 | 108, 108 | 36.4 → 35.3 | 3 levels | 0.999 |
| 32 | 132, 140 | 56.7 → 76.2 | w1-500-1e-8 L5 | 1.013 |
| 33 | 136, 140 | 110.9 → 78.4 | w1-500-1e-3 L4 | 1.014 |
| 41 | 168, 168 | 89.9 → 88.7 | 4 levels | 1.000 |
| 45 | 184, 192 | 254.5 → 108.2 | w1-500-1e-8 L4 | 1.006 |
| 61 | 248, 252 | 660.1 → 201.5 | w1-500-1e-3 L3 | 0.979 |
| 71 | 288, 288 | 284.1 → 296.4 | 4 levels | 1.000–1.001 |
| 76 | 308, 320 | 420.3 → 336.1 | w1-500-1e-8 L3 | 1.031 |
| 116 | 468, 480 | 1461.3 → 1106.9 | w1-500-1e-3 L2 | 1.017 |
| 129 | 520, 540 | 1919.9 → 1409.7 | 3 levels | 1.026–1.029 |
| 135 | 544, 560 | 2408.1 → 1577.9 | w1-500-1e-8 L2 | 1.017 |
| 241 | 968, 972 | 5945.5 → 5488.2 | w1-1000-1e-6 L2 | 1.002 |

Per case, the smooth grid adds 0.000–0.029 s to the evaluation: 1.000–1.004× rustfft's on the
product's grid. Where the product's grid has a large prime factor (p = 61: 248 = 8 · 31; p = 45:
184 = 8 · 23) rustfft gains up to 3.3× per box on the smooth grid, but the FFT is a third or less
of the level's thread time and the pairs, the larger part, grow by N_φ's ratio (up to 7.7 % at
p = 12, 3.8 % at p = 129). The price leaves out the transfer functions' and signatures' larger
tables (more memory, a longer build), which would weigh further against the smooth grid.

**At an upper bound of the FFT's cost in the level.** S9's per-box times are of one box on one
thread. In the level call rustfft's whole conversion (the Legendre step, the scatter and the FFT)
takes 3.3–12× S9's FFT alone per box (largest on the coarse levels: 12× at p = 241, 9.4× at
p = 129). Scaling S9's smooth-less-product difference by that factor, which charges the whole
conversion to the FFT, the levels move by 0.90–1.03× (the smooth grid gaining at p = 61, 0.90×;
45, 0.95×; 12, 0.95×; 135, 0.98×; 241, 0.99×) and the evaluations by 0.991–1.001×. The smooth grid
then gains at most 0.9 % of one evaluation (W1 at kD = 500, 10⁻³), and loses up to 0.1 % elsewhere.

## 5. The time rule (Kathleen)

The question (roadmap, from CH4.6a1): CH4.6a1 kept the product's `DiagonalCosts`, so rustfft never
changes which levels take the form. Does the rule need an FFT variant? Per level in the rule's
range (k a ≥ 2.8), the model with the form's S → F and F → R units (`DIAGONAL_COSTS_X86`, once per
box converted back) scaled by the case's measured rustfft ÷ product ratio of the conversions'
thread time (0.23–0.59), against the forced form's measured M2L (`--diagonal-costs form`)
(`ch4-6a2-tables.py rule`).

- **The ABBA cases**: every one of their 27 levels in the range already ran the form under the
  product's costs; a cheaper form keeps them. W3 at kD = 500, 10⁻⁶, level 5 (k a 13.5, a level of
  1.3 ms) runs RCR; the model keeps RCR, and the forced form with rustfft measures 5× RCR there.
- **The `waves` cases** (one launch per setting), the levels the defaults keep on the waves:

| case | level | p | k a | waves s (default P, R) | form P s | form R s | form R ÷ waves | model diag, RCR, waves s | model diag with R | choice |
| --- | ---: | ---: | ---: | --- | ---: | ---: | ---: | --- | ---: | --- |
| W1 kD 100, 10⁻⁸ | 3 | 32 | 10.83 | 0.0543, 0.0524 | 0.1024 | 0.0928 | 1.709 | 0.0758, 0.0922, 0.0589 | 0.0708 | waves |
| W1 kD 200, 10⁻⁸ | 4 | 32 | 10.83 | 0.2444, 0.2414 | 0.2984 | 0.3172 | 1.298 | 0.3311, 0.3994, 0.2513 | 0.3078 | waves |
| W3 kD 200, 10⁻⁶ | 4 | 22 | 10.83 | 0.3786, 0.3724 | 0.3897 | 0.3664 | 0.968 | 0.3623, 1.4925, 0.3611 | 0.3265 | **waves → form** |
| W3 kD 200, 10⁻⁸ | 4 | 32 | 10.83 | 0.7026, 0.7051 | 2.1504 | 2.0741 | 2.952 | 3.1994, 4.4083, 0.8873 | 3.1227 | waves |
| W3 kD 500, 10⁻⁸ | 4 | 45 | 27.06 | 2.4680, 2.4677 | 2.1383 | 2.0062 | 0.813 | 2.2582, 11.9401, 2.0875 | 2.1186 | waves |

  The one flip would gain 0.006 s against the same setting's waves (0.3664 s against 0.3724 s,
  1.6 % of that level; 0.012 s against the product launch's), within these single launches' noise.
  W3 at kD = 500, 10⁻⁸, level 4 is CH4.3b2c's known miss (the form beats the waves with the
  product too, 0.87×); rustfft does not change it. Every other level of these cases already ran the
  form; W1 at kD = 100 and 200 at 10⁻⁷ keep no level on the waves under the current defaults.
- **So**: on Kathleen an FFT variant of `DiagonalCosts` would change one level of the 13 cases, for
  about 1.6 % of that level, within noise. Recommendation: no variant now; the product's costs stay the rule's
  under either step (question 2).

## 6. locust (1 × 72)

The same stages on locust (Neoverse-V2, rustfft on NEON), 2026-11-26 19:35–19:44, started by the
gate of the plan once the node was quiet (no other process of this account above 20 % of a core, a
one-minute load below 2; a vLLM server of another user, about 0.3–0.6 core, ran throughout). The
load printed per launch is the decay of the launches' own 72 threads. A first `levels` launch at
18:18, on a busy node, was stopped and its outputs set aside
(`/data/ucahtbe/fmm/runs/ch4-6a2-discarded-not-quiet`). Raw outputs:
`~/projects/fmm-archive/locust/runs/ch4-6a2/`. NEON's time-rule costs put more levels on the form
than x86's (levels of p = 8 and 15 too).

**The levels** (as §2):

| case | level | p | S → F µs P → R | ÷ | F → R µs P → R | ÷ | M2L s P → R | ÷ |
| --- | ---: | ---: | --- | ---: | --- | ---: | --- | ---: |
| w1-1000-1e-6 | 2 | 241 | 26055.6 → 8380.0 | 0.322 | 43230.2 → 28504.3 | 0.659 | 0.3637 → 0.2816 | 0.774 |
| w1-1000-1e-6 | 3 | 129 | 6177.7 → 3211.7 | 0.520 | 11735.5 → 6186.3 | 0.527 | 0.3475 → 0.3130 | 0.901 |
| w1-1000-1e-6 | 4 | 71 | 1193.6 → 613.6 | 0.514 | 1971.9 → 905.8 | 0.459 | 0.1842 → 0.1722 | 0.935 |
| w1-1000-1e-6 | 5 | 41 | 298.5 → 141.7 | 0.475 | 480.1 → 170.7 | 0.356 | 0.1958 → 0.1790 | 0.914 |
| w1-1000-1e-6 | 6 | 26 | 106.1 → 64.6 | 0.609 | 139.2 → 53.6 | 0.385 | 0.2840 → 0.2663 | 0.938 |
| w1-1000-1e-6 | 7 | 17 | 43.0 → 26.8 | 0.622 | 41.6 → 19.5 | 0.469 | 0.3855 → 0.3662 | 0.950 |
| w1-1000-1e-6 | 8 | 15 | 34.5 → 12.1 | 0.351 | 25.2 → 11.1 | 0.441 | 0.7588 → 0.6637 | 0.875 |
| w1-250-1e-6 | 2 | 71 | 1157.3 → 397.3 | 0.343 | 1339.2 → 684.9 | 0.511 | 0.0095 → 0.0057 | 0.599 |
| w1-250-1e-6 | 3 | 41 | 304.2 → 112.4 | 0.369 | 370.9 → 127.0 | 0.342 | 0.0097 → 0.0074 | 0.765 |
| w1-250-1e-6 | 4 | 26 | 100.9 → 46.5 | 0.461 | 123.1 → 49.2 | 0.400 | 0.0113 → 0.0089 | 0.790 |
| w1-250-1e-6 | 5 | 17 | 41.5 → 20.2 | 0.486 | 40.6 → 19.2 | 0.472 | 0.0180 → 0.0147 | 0.819 |
| w1-250-1e-6 | 6 | 15 | 33.8 → 10.3 | 0.306 | 25.0 → 11.0 | 0.441 | 0.0469 → 0.0370 | 0.788 |
| w1-500-1e-3 | 2 | 116 | 3805.7 → 1408.7 | 0.370 | 4881.7 → 2996.3 | 0.614 | 0.0450 → 0.0347 | 0.772 |
| w1-500-1e-3 | 3 | 61 | 795.3 → 372.3 | 0.468 | 991.6 → 413.3 | 0.417 | 0.0249 → 0.0222 | 0.892 |
| w1-500-1e-3 | 4 | 33 | 181.5 → 84.2 | 0.464 | 201.5 → 83.7 | 0.415 | 0.0206 → 0.0174 | 0.845 |
| w1-500-1e-3 | 5 | 19 | 54.1 → 18.6 | 0.343 | 56.7 → 25.8 | 0.455 | 0.0258 → 0.0200 | 0.776 |
| w1-500-1e-3 | 6 | 12 | 21.7 → 8.3 | 0.383 | 21.5 → 8.8 | 0.410 | 0.0353 → 0.0291 | 0.826 |
| w1-500-1e-3 | 7 | 8 | 10.4 → 4.5 | 0.432 | 9.5 → 3.8 | 0.401 | 0.0312 → 0.0255 | 0.816 |
| w1-500-1e-6 | 2 | 129 | 4818.0 → 1833.3 | 0.381 | 6694.8 → 4350.3 | 0.650 | 0.0660 → 0.0527 | 0.799 |
| w1-500-1e-6 | 3 | 71 | 1218.5 → 550.2 | 0.452 | 1605.9 → 672.5 | 0.419 | 0.0446 → 0.0342 | 0.765 |
| w1-500-1e-6 | 4 | 41 | 296.7 → 133.7 | 0.451 | 432.5 → 143.9 | 0.333 | 0.0399 → 0.0316 | 0.792 |
| w1-500-1e-6 | 5 | 26 | 103.3 → 58.9 | 0.570 | 133.3 → 52.0 | 0.390 | 0.0662 → 0.0588 | 0.888 |
| w1-500-1e-6 | 6 | 17 | 42.2 → 25.2 | 0.597 | 41.9 → 19.4 | 0.463 | 0.0854 → 0.0779 | 0.911 |
| w1-500-1e-6 | 7 | 15 | 34.3 → 11.7 | 0.341 | 25.2 → 11.1 | 0.441 | 0.1883 → 0.1643 | 0.872 |
| w1-500-1e-8 | 2 | 135 | 5698.5 → 1847.6 | 0.324 | 8156.9 → 5091.0 | 0.624 | 0.0792 → 0.0607 | 0.766 |
| w1-500-1e-8 | 3 | 76 | 1472.0 → 668.6 | 0.454 | 1985.8 → 810.2 | 0.408 | 0.0547 → 0.0490 | 0.896 |
| w1-500-1e-8 | 4 | 45 | 371.0 → 177.1 | 0.477 | 495.2 → 186.0 | 0.376 | 0.0569 → 0.0492 | 0.864 |
| w3-500-1e-6 | 2 | 129 | 5021.7 → 1909.4 | 0.380 | 7842.3 → 5578.4 | 0.711 | 0.0808 → 0.0666 | 0.824 |
| w3-500-1e-6 | 3 | 71 | 1219.0 → 561.3 | 0.460 | 2091.0 → 740.7 | 0.354 | 0.1527 → 0.1457 | 0.954 |

**The runs** (as §3):

| case | build s P1, R1, R2, P2 | evaluate s P1, R1, R2, P2 | evaluate R ÷ P | P spread | diagonal M2L s P, R | ÷ | errors R ÷ P (φ, ∇φ) | checksums |
| --- | --- | --- | ---: | ---: | --- | ---: | --- | --- |
| w1-1000-1e-6 | 14.61, 14.39, 14.44, 14.18 | 5.830, 5.530, 5.590, 5.830 | 0.954 | 1.000 | 2.520, 2.242 | 0.890 | 0.9985, 1.0000 | agree |
| w1-250-1e-6 | 0.62, 0.61, 0.62, 0.62 | 0.290, 0.260, 0.260, 0.280 | 0.912 | 1.036 | 0.095, 0.074 | 0.773 | 1.0000, 0.9932 | agree |
| w1-500-1e-3 | 2.79, 2.80, 2.78, 2.77 | 0.780, 0.750, 0.740, 0.780 | 0.955 | 1.000 | 0.183, 0.149 | 0.815 | 1.0000, 1.0000 | agree |
| w1-500-1e-6 | 2.94, 2.88, 2.86, 2.88 | 1.230, 1.170, 1.160, 1.240 | 0.943 | 1.008 | 0.491, 0.419 | 0.855 | 1.0000, 1.0000 | agree |
| w1-500-1e-8 | 1.58, 1.58, 1.59, 1.56 | 2.280, 2.260, 2.260, 2.300 | 0.987 | 1.009 | 0.191, 0.159 | 0.833 | 1.0000, 1.0000 | agree |
| w3-500-1e-6 | 1.52, 1.51, 1.51, 1.51 | 1.980, 1.940, 1.940, 1.980 | 0.980 | 1.000 | 0.233, 0.212 | 0.909 | 1.0000, 1.0000 | agree |

- rustfft's S → F 0.31–0.62× and F → R 0.33–0.71× the product's per conversion (over each case's
  conversions 0.41–0.47×; Kathleen 0.30–0.39×); the diagonal levels 0.60–0.95×, summed per case
  0.77–0.91×; the evaluation 0.91–0.99× (each rustfft launch faster than both product launches of
  its case, the product's spread 1.000–1.036). The run gains less than on Kathleen: rustfft gains
  less against NEON's product (the conversions 0.41–0.47× against 0.30–0.39×), and at 10⁻⁸ and on
  W3 fewer of locust's levels run the form (the diagonal levels 8 % and 12 % of the evaluation,
  against 27 % and 46 % on Kathleen).
- The errors: 0.9932–1.0000× the product's (two cases move in the fourth digit, both down).
- **The smooth grid**: 1.000–1.072× per level, 1.000–1.005× per evaluation. On NEON rustfft is
  often slower on the smooth grid than on the product's (p = 241: 6,514 → 8,672 µs per box; p = 129:
  1,688 → 1,837 µs; p = 33: 68 → 88 µs), so the smooth grid loses at every degree where it differs.
  At p = 8 and 15 (not in the `ffts` list) the product's grid is the smooth one (N_φ = 36, 64).
- **The time rule** (the `rule` table; the case's rustfft ÷ product ratio of the conversions
  0.36–0.50): every level of the ABBA cases in the rule's range that runs the form keeps it. The
  18 levels kept on the waves (13 cases): the model with rustfft's conversions moves two, W1 at kD = 100,
  10⁻⁷, level 3 (the forced form with rustfft measures 1.03× the waves, a 3 % loss) and W1 at
  kD = 200, 10⁻⁷, level 4 (0.54× the waves, 0.024 s gained; but the form with the product measures
  0.61× the waves there too, a miss of the waves' model under either step). The form also beats
  the waves, under either step, on W1 at kD = 500, 10⁻⁸, level 5 (0.82× with rustfft, 0.91× with
  the product), W3 at kD = 200, 10⁻⁶, level 4 (0.86×, 1.02×) and 10⁻⁸, level 3 (0.66×, 0.88×),
  which the model keeps on the waves with either costs. These misses are the waves' model on
  locust's NEON constants, CH4.3b2c's kind ("every miss keeps the waves where the form was faster"),
  not the φ step's; a CH4.7 *Also* (roadmap) takes them to the time rule's scoring outside its fit.
  With the form forced on every level the rule reaches, W1 at kD = 200, 10⁻⁷ runs 0.36 s against the
  default's 0.33 s, under either step.

## 7. The M6

Not timed. Asked on 2026-11-26 before any run ("May I run the timing launches on the M6"); the user:
"Skip the M6". The roadmap's *Done when* names the M6 with the user's go-ahead, and §10 CH4.6a "the
time per level against the product on the M6, locust and Kathleen"; question 3 asks to sign off the
acceptance without it, as S9 did.

## 8. Caveats

- Two machine classes measured (x86 on Kathleen, NEON on locust); the M6 (NEON, another core) not.
- One rank per node: at many ranks per node each rank converts fewer boxes on fewer threads, and the
  ratio per conversion may differ (S9's one-core ratios on Kathleen were 0.05–0.11× for the φ step
  alone).
- The smooth grid's price takes S9's one-thread per-box difference as the level's; the level's FFT
  runs at some multiple of S9's per box that this task did not separate from the Legendre step (§4
  gives the price at both ends: S9's difference and the whole conversion's scale).
- The model in §5 counts the form's conversions by the level's measured F → R count, the operator
  by the boxes with a V row (`DiagonalShape::targets`); the two agree where all 34 classes pass, as
  on the flipped level, and the count is lower where fewer pass.
- The `waves` cases are single launches per setting: their small levels (10–150 ms) vary by up to
  ±10 % between settings that compute the same thing.
