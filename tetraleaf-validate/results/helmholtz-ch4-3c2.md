# Helmholtz CH4.3c2: the high-frequency degrees refitted

Helmholtz task CH4.3c2 ([docs/helmholtz/tasks/CH4.3c2.md](../../docs/helmholtz/tasks/CH4.3c2.md));
proposal [§3.5](../../docs/design/helmholtz-fmm-proposal.md) (the calibrated rule, "CH4.3
recalibrates the high-frequency levels") and §10 CH4.3; reviews 2026-11-19 A8 and 2026-11-12 M7;
CH4.0b's sign-off (CH4.3 verifies on the smooth densities). Measured on 2026-11-21 on Kathleen
(one node per job, 2 × 20 Xeon Gold 6248, one rank × 40 threads, AVX-512 P2P, BLAS threads 1,
release builds; the jobs and nodes in §8) from branch `tbetcke/h-ch4-3c2` (before and after the
tetraleaf rename); the fit and the pair evidence on the M6. Accuracy runs; the times are read,
never asserted.

## Summary

1. **The rule.** The high-frequency branch of `DegreeRule::CALIBRATED` becomes

   ```text
   p(x, ε) = max(⌈0.30 L² + 0.25 L + 2.60⌉,
                 ⌈min(x + max(0, e(ε)) x^{1/3} + 5,  x + 0.80 (ln 1/ε)^{2/3} x^{1/3})⌉),
   e(ε) = 2.30 (ln 1/ε)^{2/3} − 0.75 L − 5.73,      L = log₁₀(1/ε), x = Re(k) a:
   ```

   CH3.3's plateau; a branch fitted to the level scans with its excess over x linear in L
   beside CH3.3's (ln 1/ε)^{2/3}, raised by 2 degrees for the constant density (item 3), and
   capped by CH3.3's branch, so that no degree exceeds CH3.3's rule (kept as
   `DegreeRule::CH3_3`). It equals CH3.3's rule at x ≤ 7 and lies below it at x ≥ 43 by 4–10
   degrees at 10⁻³ and 0–3 at 10⁻⁸ (x = 216: 224 … 248 against 234 … 251).
2. **Downward scans at x ≥ 54** (random charges, from ⌊x⌋ − 2): at x = 54, 108 and 216 the FMM
   needs p − x = 1.5–28 over 10⁻³ … 10⁻⁸ where CH3.3's rule gives 12–35, the excess about
   linear in L (§2). The fit covering all 150 measured points (CH3.3's 96, this task's 54) over-
   resolves them by 188 degrees in all, against CH3.3's rule's 556 (§3).
3. **The constant density needs more.** That fit (offset d = 3) missed ε on W1 with the
   constant density: 2.6 ε at kD = 1000, 10⁻⁸, 1.2 ε at 10⁻³, 1.6 ε (∇φ) at kD = 500. The
   random charges' errors add incoherently, the constant density's on W1's sphere coherently;
   at kD = 1000, 10⁻⁸ level 2 (x = 216) carried 44 % of the error, at 10⁻³ every level a share
   (§4). Two more degrees (d = 5) bring it within 0.60 ε at kD = 250, 500 and 1000; with the
   cap (§3) the rule over-resolves the 150 points by 367 and keeps it within 0.94 ε (§6).
4. **Pair evidence.** S7's worst pair at the rule's degree is at most 4.6 ε at k a = 45 … 256,
   within the 5.6–42 ε the worst pair reaches at the degrees the FMM confirmed (§5).
5. **FMM evidence.** Every case is within ε: CH3.1c's matrix and CH3.3's held-out cases (174
   runs, Rotation and the defaults; at most 0.64 ε, CH3.3's rule 0.64 ε), and W1 at kD = 250, 500
   and 1000 with the random charges, the plane wave's trace and the constant density at 10⁻⁸,
   10⁻⁷, 10⁻⁶ and 10⁻³ (at most 0.94 ε, ∇φ of the constant density at kD = 500, 10⁻⁸, where
   CH3.3's rule gives 0.91 ε; §6).
6. **Times.** On one node, W1's evaluation under the rule against CH3.3's (two runs each,
   alternated): 11–15 % faster at 10⁻³ (kD = 1000 under Rotation 31.3 against 36.6 s), about 4 %
   at 10⁻⁶ (one of the four cases 2 % slower), no change at 10⁻⁸, where the degrees differ by at
   most 3 (§7). The build does not change.

## 1. Method

The method is CH3.3's ([helmholtz-ch3-3.md](helmholtz-ch3-3.md) §1): `helmholtz_calibrate --scan
levels` scans one level's degree upward from ⌊x⌋ − 2 while every other level is held at
`--hold 0.29,0.45,1.65,1 --hold-eps 1e-10` (far below ε), under Rotation, with §9.1's random
charges, and reports per ε the smallest degree whose larger error of φ and ∇φ (relative ℓ2 over
1,000 sampled targets, against the direct sum) meets ε/2 at three degrees in a row. Starting
below both rules makes the measured degree the smallest that suffices (A8's point).

The scans (`tools/helmholtz/ch4-3c2-runs.sh`): W1 at kD = 250 level 2 and kD = 500 levels 2 and 3
(job 653797); W4's cube at N = 10⁵ at kD = 250 level 2 and kD = 500 levels 2–4 (job 653796). At
N = 10⁵ the cube costs about 760 s per degree at kD = 1000 (job 653806, level 3, x = 108,
cancelled after 4 h with 10⁻³ … 10⁻⁵ settled: 111, 116 and 120, within the other workloads';
not in the fit), so x = 216 and a second x = 108 were scanned on the cube at N = 10⁴ (job
653817; its last launch, a third x = 108, cancelled). W1 at kD = 1000 is not scanned (about
10 min per degree, CH4.0c §2.3): the verification of §6 stands for it.

## 2. The measured degrees

The smallest degree meeting ε/2, the largest over the workloads at one x:

| x | workloads (kD, level) | 10⁻³ | 10⁻⁴ | 10⁻⁵ | 10⁻⁶ | 10⁻⁷ | 10⁻⁸ | CH3.3's rule | the rule |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- |
| 27.1 | cube 10⁵ (500, 4) | 32 | 34 | 37 | 39 | 41 | 44 | 36–44 | 34–44 |
| 54.1 | W1 (250, 2; 500, 3), cube 10⁵ (250, 2; 500, 3) | 58 | 62 | 65 | 68 | 71 | 74 | 66–76 | 61–76 |
| 108.3 | W1 (500, 2), cube 10⁵ and 10⁴ (500, 2) | 112 | 117 | 121 | 126 | 129 | 132 | 123–135 | 115–134 |
| 216.5 | cube 10⁴ (1000, 2) | 218 | 225 | 230 | 235 | 240 | 244 | 234–251 | 224–248 |

At one x the workloads agree within 1–2 degrees (x = 54: W1 57–74, the cube 58–74; x = 108: W1
111–132, the cube at N = 10⁵ 110–131, at N = 10⁴ 112–132). The excess p − x at x = 216 is 1.5–27.5
over the decades, about 5 per decade of L; CH3.3's 0.80 (ln 1/ε)^{2/3} x^{1/3} gives 17–35.

## 3. The fit

`helmholtz_calibrate refit` (`helmholtz_calibration::fit_branch`) keeps CH3.3's plateau and fits
the branch ⌈x + (C (ln 1/ε)^{2/3} + a L + b) x^{1/3} + d⌉ over grids of C, a and d, with b the
smallest multiple of 0.01 that covers every point above the plateau, minimising the total excess
Σ (rule − p) over the 150 points:

| form | C | a | b | d | Σ (rule − p) |
| --- | ---: | ---: | ---: | ---: | ---: |
| CH3.3's rule | 0.80 | 0 | 0 | 0 | 556 |
| A: (a L + b) x^{1/3} (CH3.3 §6) | 0 | 0.80 | −0.24 | 0 | 486 |
| B: A + d | 0 | 0.87 | −2.38 | 3.25 | 185 |
| C: B with C free | 1.80 | −0.40 | −5.18 | 3.50 | 169 |
| B, at most CH3.3's rule at the decades | 0 | 0.57 | 0.50 | 1.25 | 440 |
| C, at most CH3.3's rule at the decades (the fit) | 2.30 | −0.75 | −5.73 | 3 | 188 |
| **the fit raised by 2, capped by CH3.3's branch (the rule)** | **2.30** | **−0.75** | **−5.73** | **5** | **367** |

A alone gains little: covering x = 5–14 at 10⁻³ forces a coefficient there that over-resolves
the large x. B and C uncapped raise the degree by one at x = 10–27 near 10⁻⁸ (28 against 27 at
x = 13.5), on W1 at kD = 1000 the costly lower levels; hence a cap at CH3.3's rule. The fit held
it at the decades 10⁻³ … 10⁻⁸ on 4 ≤ x ≤ 512 (one degree above between them at 2 of 5,763
points); the rule holds it everywhere, as min(branch, CH3.3's branch) (`DegreeRule::cap`). The
raise by 2 (§4) with the cap leaves the rule CH3.3's at x ≤ 7 and within one of it to x ≈ 14.

e(ε) < 0 above ε ≈ 2·10⁻³, outside the measured range, where the fitted terms alone would put the
degree below x; the excess is clamped at 0 (`DegreeRule::high`), which changes no degree in D1
and none of CH3.3's rule. The rule is at most S7's closed form and monotone in x and ε over the
tests' grids (10⁻¹⁰ ≤ ε ≤ 10⁻³, x ≤ 304).

## 4. The constant density

The fit (d = 3), verified as the plan says on W1 at kD = 500 and 1000 with the three densities
(jobs 663562, 663563), kept the random charges and the plane wave's trace within ε but not the
constant density (Rotation; the library's default the same to two digits):

| W1, constant density | 10⁻³ | 10⁻⁶ | 10⁻⁷ | 10⁻⁸ |
| --- | --- | --- | --- | --- |
| kD = 1000, the fit: φ, ∇φ / ε | **1.22**, 0.78 | **1.03**, 0.76 | **2.16**, 0.92 | **2.60**, **1.48** |
| kD = 1000, CH3.3's rule (job 663571) | 0.02, 0.03 | 0.25, 0.43 | 0.19, 0.20 | 0.37, 0.46 |
| kD = 500, the fit | 0.32, 0.79 | 0.28, **1.12** | 0.46, 0.92 | 0.36, **1.58** |
| kD = 500, CH3.3's rule | 0.01, 0.04 | 0.07, 0.67 | 0.10, 0.46 | 0.14, 0.91 |

The level scans use random charges, whose per-box errors add with random phases; the constant
density on W1's sphere adds them coherently. `helmholtz_boxes --study levels` (job 663572;
each level raised in turn to the rule's degree at ε/1,000, its share of Σ|e|² over the samples)
at kD = 1000: at 10⁻⁸ level 2 (x = 216, p = 246) carries 0.44 of φ's, level 3 0.14, level 6 0.09
of ∇φ's, and level 2 raised to 255 alone halves the error (2.60 → 1.18 ε); at 10⁻³ every level
carries a share (levels 2–4: 0.31, 0.22, 0.13 of φ's; levels 6–7: 0.17 and 0.26 of ∇φ's).

The fit's offset raised by Δ, uncapped (job 663585), Rotation, φ, ∇φ / ε at the worst ε:

| Δ (d) | kD = 250 | kD = 500 | kD = 1000 |
| --- | --- | --- | --- |
| +2 (5) | 0.12, 0.19 | 0.11, 0.36 | 0.60, 0.54 |
| +4 (7) | 0.03, 0.06 | 0.03, 0.25 | 0.27, 0.16 |
| +6 (9) | 0.01, 0.04 | 0.01, 0.24 | 0.13, 0.06 |

Δ = 2 suffices at every ε and kD; the rule is the fit raised by 2 and capped (§3).

## 5. Pair evidence

`tools/helmholtz/ch4-3c2-pairs.py` reads S7's per-pair worst case
(`spikes/truncation/results/c1_curves.csv`: the one-level chain in extended precision, the worst of
316 offsets and three placements, φ and ∇φ, per-pair relative ℓ2) at the FMM's measured degrees
(x ≥ 43, interpolated between S7's half-octave grid) and at the rule's degree on the grid:

| ε | worst pair / ε at the FMM's degrees, largest | the rule's, largest over k a = 45 … 256 | the fit's (d = 3) |
| --- | ---: | ---: | ---: |
| 10⁻³ | 19.1 | 3.0 | 6.4 |
| 10⁻⁴ | 25.9 | 2.4 | 6.4 |
| 10⁻⁵ | 25.6 | 2.4 | 8.5 |
| 10⁻⁶ | 27.2 | 2.2 | 9.3 |
| 10⁻⁷ | 42.3 | 2.8 | 10.0 |
| 10⁻⁸ | 33.9 | 4.6 | 18.5 |

The FMM needs less than the worst pair: at its measured degrees (random charges) the worst pair is
5.6–42 ε, with no trend in x. The rule's stays within that at every grid k a and ε, as the fit's
did; the fit's FMM failure on the constant density (§4) was not visible to this measure, which
guards the extrapolation to k a = 256 and is coarse (one degree moves the per-pair error about
2×). Under CH3.3's rule the worst pair is 0.03–1.8 ε.

## 6. The FMM's errors under the rule

**The matrix** (job 663586 for the rule, 663560 for CH3.3's rule): `helmholtz_accuracy` on CH3.1c's
matrix (W1 at kD = 10, 30, 100, 500; W3 at kD = 10, 30, 100; W4's four distributions at kD = 1, 10,
100; W1 with real charges and with complex k = 100 + 10i) and CH3.3's held-out cases (W2 at kD = 30,
100; W1 at kD = 50, 150; W4's cube and Plummer at kD = 3, 30; W1 at 10⁻⁴ and 10⁻⁷), at 10⁻³, 10⁻⁶
and 10⁻⁸, under Rotation (`--m2l rotation`) and the library's defaults (`--m2l auto-pw --near-v cost
--near-wx cost`). All 174 runs are within ε; the largest error over ε is 0.638 under Rotation (W1
at kD = 500, 10⁻⁸; CH3.3's rule 0.636) and 0.353 under the defaults (CH3.3's 0.353). The runs whose
degrees differ (the larger of φ's and ∇φ's error over ε; the evaluation seconds ran on different
nodes, so they are not compared):

| M2L | case | N | kD | ε | p, default | p, CH3.3 | error / ε, default | CH3.3 | evaluate s, default / CH3.3 |
| --- | --- | ---: | --- | --- | --- | --- | ---: | ---: | --- |
| defaults | w1 | 1989437 | 500 | 1e-08 | 0×2, 134, 76, 44, 32×2 | 0×2, 135, 76, 44, 32×2 | 0.327 | 0.321 | 9.94 / 9.34 |
| defaults | w1 | 1989437 | 500 | 1e-06 | 0×2, 128, 71, 41, 25, 16 | 0×2, 131, 72, 41, 25, 16 | 0.325 | 0.318 | 6.87 / 6.32 |
| defaults | w1 | 1989437 | 500 | 1e-03 | 0×2, 115, 61, 34, 20, 13 | 0×2, 123, 66, 36, 21, 13 | 0.105 | 0.021 | 5.25 / 5.07 |
| defaults | w1 | 79578 | 100 | 1e-03 | 0×2, 28, 17, 11 | 0×2, 30, 18, 11 | 0.046 | 0.026 | 0.12 / 0.11 |
| defaults | w1-complex-k | 79578 | 100 + 10i | 1e-03 | 0×2, 28, 17, 11 | 0×2, 30, 18, 11 | 0.015 | 0.015 | 0.12 / 0.12 |
| defaults | w1-held | 179050 | 150 | 1e-03 | 0×2, 39, 23, 14, 9 | 0×2, 42, 24, 14, 9 | 0.065 | 0.038 | 0.29 / 0.28 |
| defaults | w1-held | 19895 | 50 | 1e-03 | 0×2, 17, 11 | 0×2, 18, 11 | 0.026 | 0.020 | 0.03 / 0.03 |
| defaults | w1-mid | 79578 | 100 | 1e-04 | 0×2, 31, 19, 12 | 0×2, 32, 19, 12 | 0.067 | 0.065 | 0.15 / 0.15 |
| defaults | w1-real | 79578 | 100 | 1e-03 | 0×2, 28, 17, 11 | 0×2, 30, 18, 11 | 0.044 | 0.027 | 0.11 / 0.10 |
| defaults | w2 | 127554 | 100 | 1e-03 | 0×2, 28, 17, 11, 7×2 | 0×2, 30, 18, 11, 7×2 | 0.076 | 0.074 | 0.15 / 0.14 |
| defaults | w3 | 4031442 | 100 | 1e-03 | 0×2, 28, 17, 11, 7 | 0×2, 30, 18, 11, 7 | 0.103 | 0.097 | 3.81 / 2.63 |
| defaults | w4-clusters | 100000 | 100 | 1e-03 | 0×2, 20, 13, 8, 7×3 | 0×2, 21, 13, 8, 7×3 | 0.064 | 0.064 | 0.10 / 0.10 |
| defaults | w4-cube | 100000 | 100 | 1e-03 | 0×2, 28, 17 | 0×2, 30, 18 | 0.051 | 0.011 | 0.16 / 0.16 |
| defaults | w4-plummer | 100000 | 100 | 1e-03 | 0×2, 28, 17, 11, 7×3 | 0×2, 30, 18, 11, 7×3 | 0.101 | 0.101 | 0.19 / 0.18 |
| defaults | w4-sphere | 100000 | 100 | 1e-03 | 0×2, 28, 17, 11 | 0×2, 30, 18, 11 | 0.022 | 0.007 | 0.12 / 0.12 |
| rotation | w1 | 1989437 | 500 | 1e-08 | 0×2, 134, 76, 44, 27, 24 | 0×2, 135, 76, 44, 27, 24 | 0.638 | 0.636 | 10.80 / 9.95 |
| rotation | w1 | 1989437 | 500 | 1e-06 | 0×2, 128, 71, 41, 25, 16 | 0×2, 131, 72, 41, 25, 16 | 0.362 | 0.356 | 7.30 / 6.78 |
| rotation | w1 | 1989437 | 500 | 1e-03 | 0×2, 115, 61, 34, 20, 13 | 0×2, 123, 66, 36, 21, 13 | 0.105 | 0.021 | 5.18 / 4.95 |
| rotation | w1 | 79578 | 100 | 1e-03 | 0×2, 28, 17, 11 | 0×2, 30, 18, 11 | 0.047 | 0.027 | 0.16 / 0.16 |
| rotation | w1-complex-k | 79578 | 100 + 10i | 1e-03 | 0×2, 28, 17, 11 | 0×2, 30, 18, 11 | 0.017 | 0.016 | 0.17 / 0.16 |
| rotation | w1-held | 179050 | 150 | 1e-03 | 0×2, 39, 23, 14, 9 | 0×2, 42, 24, 14, 9 | 0.067 | 0.042 | 0.37 / 0.35 |
| rotation | w1-held | 19895 | 50 | 1e-03 | 0×2, 17, 11 | 0×2, 18, 11 | 0.027 | 0.021 | 0.03 / 0.03 |
| rotation | w1-mid | 79578 | 100 | 1e-04 | 0×2, 31, 19, 12 | 0×2, 32, 19, 12 | 0.069 | 0.066 | 0.19 / 0.17 |
| rotation | w1-real | 79578 | 100 | 1e-03 | 0×2, 28, 17, 11 | 0×2, 30, 18, 11 | 0.045 | 0.028 | 0.15 / 0.15 |
| rotation | w2 | 127554 | 100 | 1e-03 | 0×2, 28, 17, 11, 7×2 | 0×2, 30, 18, 11, 7×2 | 0.079 | 0.077 | 0.21 / 0.20 |
| rotation | w3 | 4031442 | 100 | 1e-03 | 0×2, 28, 17, 11, 7 | 0×2, 30, 18, 11, 7 | 0.103 | 0.097 | 4.82 / 3.24 |
| rotation | w4-clusters | 100000 | 100 | 1e-03 | 0×2, 20, 13, 8, 7×3 | 0×2, 21, 13, 8, 7×3 | 0.069 | 0.069 | 0.18 / 0.17 |
| rotation | w4-cube | 100000 | 100 | 1e-03 | 0×2, 28, 17 | 0×2, 30, 18 | 0.051 | 0.011 | 0.18 / 0.19 |
| rotation | w4-plummer | 100000 | 100 | 1e-03 | 0×2, 28, 17, 11, 7×3 | 0×2, 30, 18, 11, 7×3 | 0.103 | 0.103 | 0.32 / 0.32 |
| rotation | w4-sphere | 100000 | 100 | 1e-03 | 0×2, 28, 17, 11 | 0×2, 30, 18, 11 | 0.022 | 0.007 | 0.14 / 0.14 |

**W1 at kD = 250, 500 and 1000** (jobs 663586 and 663587; `helmholtz_boxes`, §9.1's samples),
Rotation; under the library's default (Auto, the near rules off, as CH4.0c) every row is within
ε too, at most 0.74 ε (∇φ of the constant density at kD = 500, 10⁻⁸). CH3.3's rule from jobs 663562 (random) and 663571 (constant), the fit's from 663562
and 663563:

| kD | density | ε | p by level from 2, default | φ, ∇φ / ε | CH3.3: p | φ, ∇φ / ε | first fit (d = 3): φ, ∇φ / ε |
| ---: | --- | --- | --- | --- | --- | --- | --- |
| 250 | random | 1e-03 | 61, 34, 20, 13 | 0.07, 0.07 | – | – | – |
| 250 | random | 1e-06 | 71, 41, 25, 16 | 0.11, 0.22 | – | – | – |
| 250 | random | 1e-07 | 74, 43, 26, 20 | 0.06, 0.11 | – | – | – |
| 250 | random | 1e-08 | 76, 44, 27, 24 | 0.15, 0.24 | – | – | – |
| 250 | plane wave | 1e-03 | 61, 34, 20, 13 | 0.05, 0.08 | – | – | – |
| 250 | plane wave | 1e-06 | 71, 41, 25, 16 | 0.06, 0.19 | – | – | – |
| 250 | plane wave | 1e-07 | 74, 43, 26, 20 | 0.04, 0.11 | – | – | – |
| 250 | plane wave | 1e-08 | 76, 44, 27, 24 | 0.08, 0.23 | – | – | – |
| 250 | constant | 1e-03 | 61, 34, 20, 13 | 0.12, 0.15 | – | – | – |
| 250 | constant | 1e-06 | 71, 41, 25, 16 | 0.10, 0.49 | – | – | – |
| 250 | constant | 1e-07 | 74, 43, 26, 20 | 0.06, 0.27 | – | – | – |
| 250 | constant | 1e-08 | 76, 44, 27, 24 | 0.15, 0.59 | – | – | – |
| 500 | random | 1e-03 | 115, 61, 34, 20, 13 | 0.10, 0.10 | 123, 66, 36, 21, 13 | 0.01, 0.02 | 0.37, 0.44 |
| 500 | random | 1e-06 | 128, 71, 41, 25, 16 | 0.16, 0.36 | 131, 72, 41, 25, 16 | 0.15, 0.36 | 0.41, 0.59 |
| 500 | random | 1e-07 | 131, 74, 43, 26, 20 | 0.12, 0.28 | 133, 74, 43, 26, 20 | 0.10, 0.27 | 0.29, 0.41 |
| 500 | random | 1e-08 | 134, 76, 44, 27, 24 | 0.27, 0.64 | 135, 76, 44, 27, 24 | 0.27, 0.64 | 0.45, 0.74 |
| 500 | plane wave | 1e-03 | 115, 61, 34, 20, 13 | 0.09, 0.10 | – | – | 0.26, 0.47 |
| 500 | plane wave | 1e-06 | 128, 71, 41, 25, 16 | 0.07, 0.25 | – | – | 0.20, 0.39 |
| 500 | plane wave | 1e-07 | 131, 74, 43, 26, 20 | 0.06, 0.17 | – | – | 0.15, 0.26 |
| 500 | plane wave | 1e-08 | 134, 76, 44, 27, 24 | 0.14, 0.43 | – | – | 0.28, 0.49 |
| 500 | constant | 1e-03 | 115, 61, 34, 20, 13 | 0.10, 0.21 | 123, 66, 36, 21, 13 | 0.01, 0.04 | 0.32, 0.79 |
| 500 | constant | 1e-06 | 128, 71, 41, 25, 16 | 0.09, 0.71 | 131, 72, 41, 25, 16 | 0.07, 0.67 | 0.28, 1.12 |
| 500 | constant | 1e-07 | 131, 74, 43, 26, 20 | 0.14, 0.50 | 133, 74, 43, 26, 20 | 0.10, 0.46 | 0.46, 0.92 |
| 500 | constant | 1e-08 | 134, 76, 44, 27, 24 | 0.17, 0.94 | 135, 76, 44, 27, 24 | 0.14, 0.91 | 0.36, 1.58 |
| 1000 | random | 1e-03 | 224, 115, 61, 34, 20, 13 | 0.09, 0.10 | 234, 123, 66, 36, 21, 13 | 0.01, 0.02 | 0.43, 0.60 |
| 1000 | random | 1e-06 | 240, 128, 71, 41, 25, 16 | 0.13, 0.30 | 245, 131, 72, 41, 25, 16 | 0.12, 0.28 | 0.28, 0.47 |
| 1000 | random | 1e-07 | 244, 131, 74, 43, 26, 20 | 0.18, 0.17 | 248, 133, 74, 43, 26, 20 | 0.09, 0.14 | 0.25, 0.31 |
| 1000 | random | 1e-08 | 248, 134, 76, 44, 27, 24 | 0.30, 0.33 | 251, 135, 76, 44, 27, 24 | 0.20, 0.29 | 0.60, 0.58 |
| 1000 | plane wave | 1e-03 | 224, 115, 61, 34, 20, 13 | 0.04, 0.07 | – | – | 0.27, 0.49 |
| 1000 | plane wave | 1e-06 | 240, 128, 71, 41, 25, 16 | 0.06, 0.30 | – | – | 0.13, 0.40 |
| 1000 | plane wave | 1e-07 | 244, 131, 74, 43, 26, 20 | 0.06, 0.15 | – | – | 0.12, 0.21 |
| 1000 | plane wave | 1e-08 | 248, 134, 76, 44, 27, 24 | 0.11, 0.32 | – | – | 0.16, 0.35 |
| 1000 | constant | 1e-03 | 224, 115, 61, 34, 20, 13 | 0.53, 0.30 | 234, 123, 66, 36, 21, 13 | 0.02, 0.03 | 1.22, 0.78 |
| 1000 | constant | 1e-06 | 240, 128, 71, 41, 25, 16 | 0.45, 0.46 | 245, 131, 72, 41, 25, 16 | 0.25, 0.43 | 1.03, 0.76 |
| 1000 | constant | 1e-07 | 244, 131, 74, 43, 26, 20 | 0.61, 0.35 | 248, 133, 74, 43, 26, 20 | 0.19, 0.20 | 2.16, 0.92 |
| 1000 | constant | 1e-08 | 248, 134, 76, 44, 27, 24 | 0.69, 0.70 | 251, 135, 76, 44, 27, 24 | 0.37, 0.46 | 2.60, 1.48 |

## 7. Times

Job 663589 (node-c11a-144): W1 at kD = 1000 and 500, §9.1's random charges, `helmholtz_boxes` at
the samples, each case twice under each rule, the rules alternated (rule, CH3.3's, rule, CH3.3's),
Rotation and the library's default. The evaluation in seconds, the mean of the two runs (their
spread in brackets where above 2 %):

| kD | ε | M2L | p, levels 2–4, rule / CH3.3 | rule | CH3.3's | rule / CH3.3 |
| ---: | --- | --- | --- | ---: | ---: | ---: |
| 1000 | 10⁻³ | Rotation | 224, 115, 61 / 234, 123, 66 | 31.3 | 36.6 (8 %) | 0.85 |
| 1000 | 10⁻³ | default | | 32.3 | 36.7 | 0.88 |
| 1000 | 10⁻⁶ | Rotation | 240, 128, 71 / 245, 131, 72 | 42.6 | 44.3 | 0.96 |
| 1000 | 10⁻⁶ | default | | 41.5 | 43.0 (3 %) | 0.96 |
| 1000 | 10⁻⁸ | Rotation | 248, 134, 76 / 251, 135, 76 | 58.7 (7 %) | 57.8 | 1.02 |
| 1000 | 10⁻⁸ | default | | 59.8 | 60.1 (2 %) | 0.99 |
| 500 | 10⁻³ | Rotation | 115, 61, 34 / 123, 66, 36 | 4.65 (6 %) | 5.21 | 0.89 |
| 500 | 10⁻³ | default | | 4.94 | 5.67 | 0.87 |
| 500 | 10⁻⁶ | Rotation | 128, 71, 41 / 131, 72, 41 | 6.71 | 6.98 (6 %) | 0.96 |
| 500 | 10⁻⁶ | default | | 6.70 | 6.54 | 1.02 |
| 500 | 10⁻⁸ | Rotation | 134, 76, 44 / 135, 76, 44 | 9.89 | 9.86 | 1.00 |
| 500 | 10⁻⁸ | default | | 10.71 | 10.83 | 0.99 |

The levels whose degrees fall are the top ones (V pairs: level 2 2,504, level 3 12,216 at
kD = 1000), so the saving is a modest part of the evaluation; at 10⁻⁸ the 1–3 degrees are within
the runs' spread. The build (25.6–27.5 s at kD = 1000 on this node, the tables and the plan) is
the same under both rules. This node ran slower than job 663587's (the evaluation at 10⁻³ 31.3
against 27.8 s, the build 26 against 12–13 s); the ratios hold, the rules alternated on one
node. The earlier runs of §6 ran the two rules on different nodes and are not compared.

## 8. Questions for sign-off

1. Adopt the rule as the default (`DegreeRule::CALIBRATED`, §3): the random-charge fit raised by 2
   and capped by CH3.3's branch; CH3.3's rule kept as `DegreeRule::CH3_3`?
2. Amend the proposal's §3.5 ("Calibrated") with the rule, its measured range (x to 216 by level
   scans, the constant density's 2 degrees) and this report, as "(amended at CH4.3c2)" with an
   Appendix A line?
3. The calibration scans random charges only, and the constant density needed 2 more degrees at
   x ≥ 27: should a later calibration scan the constant density too (a `--density` for
   `helmholtz_calibrate`), so that the margin is measured rather than chosen?

## 9. Jobs, reproduce

Kathleen jobs (one node each, `singlenode`): scans 653796 (node-c11b-002, cancelled after its
kD = 500 launch), 653797 (c11b-004), 653806 (c11b-004, cancelled), 653817 (c11a-114, last launch
cancelled); the fit's verification 663560 (c11a-025), 663562 (c11a-027), 663563 (c11a-025); the
constant density's diagnosis 663571 (c11a-027), 663572 (c11a-044), 663585 (c11a-144); the rule's
verification 663586 (c11a-144), 663587 (c11b-043); times 663589 (c11a-144).

```sh
cargo build --release -p tetraleaf-validate --example helmholtz_calibrate \
    --example helmholtz_accuracy --example helmholtz_boxes
tools/helmholtz/ch4-3c2-runs.sh scans-w1 <out>           # and scans-cube, scans-cube-small
cargo run --release -p tetraleaf-validate --example helmholtz_calibrate -- refit <out>/*.md
tools/helmholtz/ch4-3c2-pairs.py <points> --rule rule=0.30,0.25,2.60,2.30,-0.75,-5.73,5,0.80
tools/helmholtz/ch4-3c2-runs.sh verify-matrix <out>      # and verify-w1 (ONLY= to split it)
tools/helmholtz/ch4-3c2-tables.py matrix <out>           # §6's tables; `w1 <files>` for W1
```

The raw outputs are in `~/projects/fmm-archive/kathleen/runs/ch4-3c2/`.
