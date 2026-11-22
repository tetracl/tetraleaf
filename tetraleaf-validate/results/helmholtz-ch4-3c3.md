# Helmholtz CH4.3c3: the constant density in the calibration

Helmholtz task CH4.3c3 ([docs/helmholtz/tasks/CH4.3c3.md](../../docs/helmholtz/tasks/CH4.3c3.md));
roadmap H4, item 10 (added at CH4.3c2, its question 3); proposal
[§3.5](../../docs/design/helmholtz-fmm-proposal.md) ("Calibrated", amended at CH4.3c2) and §10
CH4.3; builds on [helmholtz-ch4-3c2.md](helmholtz-ch4-3c2.md) (the refit, its +2 for the constant
density chosen by hand). Measured on 2026-11-22 on Kathleen (one node per job, 2 × 20 Xeon Gold
6248, one rank × 40 threads, AVX-512 P2P, BLAS threads 1, release builds; the jobs and nodes in
§9) from branch `tbetcke/h-ch4-3c3`; the fit and the pair evidence on the M6. Accuracy runs; the
times are read, never asserted.

## Summary

1. **The margin, measured.** Level scans of W1 with the constant density (q = 1 at every point)
   at every level of x = 6.8–216 at kD = 250, 500 and 1000 (14 scans, §2) need 0–3 degrees over
   CH4.3c2's random-charge fit (§3). The +2 that CH4.3c2 chose is wrong both ways: **one short**
   at x = 216 (10⁻³ and 10⁻⁷) and at x = 13.5, 10⁻⁸ (two of the three trees), exactly enough at
   x = 108 and 10⁻⁷–10⁻⁸, and **up to 2 too many** at x ≤ 54, where the constant density needs at
   most one degree more than the random charges at x ≤ 27 and 1–2 at x = 54.
2. **The refit on the worst density.** CH4.3c2's procedure over CH3.3's points, CH4.3c2's
   random-charge scans and this task's constant-density scans (234 points) cannot keep CH4.3c2's
   cap (CH3.3's branch): x = 13.5 at 10⁻⁸ needs 28, CH3.3's rule gives 27. The least raised cap
   that admits a covering fit is x + 0.82 (ln 1/ε)^{2/3} x^{1/3}, and the fit under it is

   ```text
   p(x, ε) = max(⌈0.30 L² + 0.25 L + 2.60⌉,
                 ⌈min(x + max(0, e(ε)) x^{1/3} + 2,  x + 0.82 (ln 1/ε)^{2/3} x^{1/3})⌉),
   e(ε) = 2.00 (ln 1/ε)^{2/3} − 0.55 L − 4.66,
   ```

   no hand-chosen offset. It covers every measured degree with a total excess of 349 (CH4.3c2's
   rule covers 230 of the 234 with 472); against CH4.3c2's rule it is 0–2 degrees lower at
   x ≤ 150 but one higher at x = 8–21 and ε ≤ 2 · 10⁻⁷ (where it exceeds CH3.3's rule, x = 13.5 at
   10⁻⁸ among them), and −1 to +2 beyond x = 150 (§4). It is the
   default (`DegreeRule::CALIBRATED`, signed off, §8; CH4.3c2's rule kept as
   `DegreeRule::CH4_3C2`).
3. **Pair evidence.** S7's worst pair at the refit's degree is at most 4.8 ε at k a = 45 … 256
   (CH4.3c2's rule 4.6 ε), within the 19–42 ε the worst pair reaches at the degrees the FMM
   confirmed at x ≥ 43 (§5).
4. **FMM evidence.** Every case within ε under the refit: CH3.1c's matrix and CH3.3's held-out
   cases (174 runs, Rotation and the defaults; at most 0.44 ε, CH4.3c2's rule 0.64 ε), and W1 at
   kD = 250, 500 and 1000 with the random charges, the plane wave's trace and the constant density
   at 10⁻⁸, 10⁻⁷, 10⁻⁶ and 10⁻³ (at most 0.65 ε; CH4.3c2's rule 0.94 ε) (§6). The lower degrees
   raise the errors at 10⁻³ (W1 at kD = 100: 0.05 → 0.26 ε; the constant density on W1 at
   kD = 1000: 0.53 → 0.65 ε) and the raised x = 13.5 level lowers them at 10⁻⁸ (W1 at kD = 500,
   constant: 0.94 → 0.49 ε).
5. **Times.** W1's evaluation on one node, the rules alternated: 7–14 % faster at 10⁻³
   (kD = 500 and 1000), unchanged within 2 % at 10⁻⁶ and 10⁻⁸ (§7).

## 1. Method

CH3.3's and CH4.3c2's ([helmholtz-ch4-3c2.md](helmholtz-ch4-3c2.md) §1): `helmholtz_calibrate
--scan levels` scans one level's degree upward from ⌊x⌋ − 2 while every other level is held at
`--hold 0.29,0.45,1.65,1 --hold-eps 1e-10` (far below ε), under Rotation, and reports per ε the
smallest degree whose larger error of φ and ∇φ (relative ℓ2 over 1,000 sampled targets, against
the direct sum) meets ε/2 at three degrees in a row. New here is `--density` (`Case::density`,
CH4.0b's): the scans run W1 with q_j = 1 (complex charges with Im q = 0), whose boxes' errors
add coherently on the sphere, where §9.1's random charges add them with random phases. A point
of a density other than the random charges carries it in its workload field (`W1/constant`),
which `fit` and `refit` read (`helmholtz_calibration::PointLine`).

The scans (`tools/helmholtz/ch4-3c3-runs.sh scans`, `scans-low`): kD = 250 levels 2–5 (x = 54.1,
27.1, 13.5, 6.8), kD = 500 levels 2–6 (108.3 … 6.8), kD = 1000 levels 2–6 (216.5 … 13.5), on the
library's trees of today (c_leaf 4 since CH4.3e; CH4.3c2 scanned at c_leaf 2). The plan's range
was x = 13.5–216; x = 6.8 was added once the refit (§4) proved to lower the degrees there. A
degree costs about 20 s at kD = 500 and 110 s at kD = 1000 on a node (CH4.3c2's estimate of 10 min
at kD = 1000 was for its cube at N = 10⁵), so every scan ran as one launch: about 4.8 node-hours in all.

## 2. The constant density's degrees

The smallest degree meeting ε/2 per scan:

| kD | level | x | 10⁻³ | 10⁻⁴ | 10⁻⁵ | 10⁻⁶ | 10⁻⁷ | 10⁻⁸ |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 250 | 2 | 54.13 | 60 | 63 | 67 | 70 | 72 | 75 |
| 250 | 3 | 27.06 | 32 | 34 | 37 | 39 | 42 | 44 |
| 250 | 4 | 13.53 | 18 | 20 | 22 | 24 | 26 | 28 |
| 250 | 5 | 6.77 | 11 | 13 | 14 | 16 | 19 | 22 |
| 500 | 2 | 108.26 | 113 | 119 | 123 | 127 | 131 | 134 |
| 500 | 3 | 54.13 | 58 | 62 | 66 | 69 | 71 | 74 |
| 500 | 4 | 27.06 | 32 | 35 | 37 | 39 | 42 | 44 |
| 500 | 5 | 13.53 | 18 | 21 | 23 | 25 | 26 | 28 |
| 500 | 6 | 6.77 | 11 | 13 | 15 | 16 | 18 | 22 |
| 1000 | 2 | 216.51 | 225 | 230 | 236 | 240 | 245 | 248 |
| 1000 | 3 | 108.26 | 114 | 119 | 123 | 127 | 130 | 133 |
| 1000 | 4 | 54.13 | 59 | 63 | 66 | 69 | 71 | 74 |
| 1000 | 5 | 27.06 | 32 | 35 | 37 | 39 | 41 | 43 |
| 1000 | 6 | 13.53 | 18 | 20 | 22 | 24 | 26 | 27 |

At one x the trees agree within 1–2 degrees, the larger tree not always the larger degree
(x = 54: kD = 250 needs 60–75, kD = 1000 59–74). The two 28s at x = 13.5 and 10⁻⁸ are narrow:
at p = 27 the gradient's error is 5.1·10⁻⁹ (kD = 250) and 8.5·10⁻⁹ (kD = 500) against the target
5·10⁻⁹, 0.51 ε and 0.85 ε; kD = 1000 meets it at 27 (4.4·10⁻⁹). The degree meeting ε itself is
at most CH4.3c2's rule everywhere.

## 3. The margin

The largest degree over the trees at one x (`ch4-3c3-tables.py margin`) against the random
charges' at the same x (CH4.3c2's §2; CH3.3's clusters at x = 13.57), CH4.3c2's fit (d = 3) and
CH4.3c2's rule (d = 5, the fit raised by 2):

| x | ε | constant | random | fit | CH4.3c2's rule | constant − fit | rule − constant |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 6.77 | 10⁻³ … 10⁻⁸ | 11, 13, 15, 16, 19, 22 | – | 11, 13, 15, 16, 20, 24 | 13, 14, 15, 16, 20, 24 | −2 … 0 | 0 … +2 |
| 13.53 | 10⁻³ … 10⁻⁸ | 18, 21, 23, 25, 26, 28 | 13, 16, 18, 20, 22, 24 | 18, 20, 22, 24, 26, 27 | 20, 22, 24, 25, 26, 27 | 0 … +1 | **−1** … +2 |
| 27.06 | 10⁻³ … 10⁻⁸ | 32, 35, 37, 39, 42, 44 | 32, 34, 37, 39, 41, 44 | 32, 35, 37, 40, 42, 44 | 34, 37, 39, 41, 43, 44 | −1 … 0 | 0 … +2 |
| 54.13 | 10⁻³ … 10⁻⁸ | 60, 63, 67, 70, 72, 75 | 58, 62, 65, 68, 71, 74 | 59, 63, 66, 69, 72, 74 | 61, 65, 68, 71, 74, 76 | 0 … +1 | +1 … +2 |
| 108.25 | 10⁻³ … 10⁻⁸ | 114, 119, 123, 127, 131, 134 | 112, 117, 121, 126, 129, 132 | 113, 118, 122, 126, 129, 132 | 115, 120, 124, 128, 131, 134 | +1 … +2 | 0 … +1 |
| 216.51 | 10⁻³ … 10⁻⁸ | 225, 230, 236, 240, 245, 248 | 218, 225, 230, 235, 240, 244 | 222, 228, 234, 238, 242, 246 | 224, 230, 236, 240, 244, 248 | +2 … +3 | **−1** … 0 |

The constant density's excess over the random charges grows with x: at most 1 at x ≤ 27 (at
x = 13.5 the random row is another workload's), 1–2 at x = 54–108, 4–7 at x = 216. CH4.3c2's
+2, one offset for every x, is too small at x = 216 (10⁻³ and 10⁻⁷) and x = 13.5 (10⁻⁸) and too
large at x ≤ 54 (spare up to 2). By the plan's criterion (too small: a degree above
the rule; too large: a refit on the worst density 2 or more below the rule at a measured point)
it is wrong both ways, so the rule is refitted.

## 4. The refit

`helmholtz_calibrate refit` over CH3.3's 96 level points, CH4.3c2's 54 random-charge points and
this task's 84 constant-density points (each x > 4), forms as CH4.3c2 (b the least multiple of
0.01 that covers, the least total excess Σ (rule − p)):

| rule | C | a | b | d | cap | covers | Σ (rule − p) |
| --- | ---: | ---: | ---: | ---: | ---: | --- | ---: |
| CH3.3's rule | 0.80 | 0 | 0 | 0 | – | no (x = 13.5, 10⁻⁸) | 793 |
| CH4.3c2's fit | 2.30 | −0.75 | −5.73 | 3 | 0.80 | no | 170 |
| CH4.3c2's rule | 2.30 | −0.75 | −5.73 | 5 | 0.80 | no (4 points) | 472 |
| A: (a L + b) x^{1/3} | 0 | 0.73 | 0.18 | 0 | – | yes | 707 |
| B: A + d | 0 | 0.80 | −1.23 | 2 | – | yes | 411 |
| C: B with C free | 1.80 | −0.40 | −4.41 | 2 | – | yes | 346 |
| B, C at most CH3.3's rule | | | | | 0.80 | no rule fits | |
| **C under the cap 0.82 (the refit)** | **2.00** | **−0.55** | **−4.66** | **2** | **0.82** | **yes** | **349** |

No rule at or below CH3.3's covers x = 13.5 at 10⁻⁸ (28 against 27). `refit` raises the cap's
coefficient from 0.80 by 0.01 to the first that admits a covering fit, 0.82, holds the fit under
it on `cap_grid` and then takes the cap everywhere (`DegreeRule::cap`), as CH4.3c2 built its
default. The cap costs 3 degrees of total excess against C uncapped (349 against 346) and keeps
the low levels at most x + 0.82 (ln 1/ε)^{2/3} x^{1/3}: C uncapped is above that at x ≈ 5 and
10⁻⁸.

The refit against CH4.3c2's rule (the difference refit − CH4.3c2's per ε, 10⁻³ … 10⁻⁸):

| x | CH4.3c2's rule | refit − CH4.3c2 | refit − CH3.3 |
| ---: | --- | --- | --- |
| 4 | 9, 10, 12, 15, 20, 24 | −1, 0, 0, 0, 0, 0 | −1, 0, 0, 0, 0, 0 |
| 6.8 | 13, 14, 15, 16, 20, 24 | −2, −1, 0, 0, 0, 0 | −2, −1, 0, 0, 0, 0 |
| 10.8 | 17, 19, 20, 21, 23, 24 | −2, −1, −1, 0, 0, 0 | −3, −1, −1, 0, 0, 0 |
| 13.5 | 20, 22, 24, 25, 26, 27 | −2, −1, −1, −1, 0, **+1** | −3, −1, −1, −1, 0, **+1** |
| 27.1 | 34, 37, 39, 41, 43, 44 | −2, −2, −1, −1, −1, 0 | −4, −3, −2, −1, −1, 0 |
| 54.1 | 61, 65, 68, 71, 74, 76 | −1, −1, −1, −1, −1, −1 | −6, −4, −3, −2, −1, −1 |
| 108.3 | 116, 120, 125, 128, 131, 134 | −1, 0, −1, 0, 0, 0 | −8, −6, −4, −3, −2, −1 |
| 216.5 | 224, 230, 235, 240, 244, 248 | **+1, +1, +1**, 0, **+1**, 0 | −9, −7, −5, −5, −3, −3 |
| 433 | 441, 449, 455, 461, 466, 471 | +2, +1, +2, +1, +2, +1 | −12, −10, −7, −6, −4, −4 |

On the grid of `cap_grid` by tenths of a decade (4 ≤ x ≤ 512, 10⁻⁸ … 10⁻³; 5,763 points) the
refit lies above CH3.3's rule at 39 points, by one degree each, at x = 8–21 and ε ≤ 10^{−6.7}
(tested: `the_calibrated_rule_is_at_most_its_cap`). The constant density's x = 6.8 scans, which
the plan did not have, bound the refit there exactly: 11 and 13 at 10⁻³ and 10⁻⁴, as measured.
The refit is at most S7's closed form and monotone in x and ε on the tests' grids. Below x = 6.8
the plateau decides most degrees and no scan used the constant density (CH3.3's plateau scans are
the random charges'); the refit is 1–2 degrees below CH4.3c2's at x ≈ 4–7 and 10⁻³, which CH4.4's
sweep checks with the constant density (an Also on CH4.4).

## 5. Pair evidence

`tools/helmholtz/ch4-3c2-pairs.py` (spike S7's per-pair worst case, the one-level chain in
extended precision over 316 offsets and three placements, φ and ∇φ) on the 150 FMM points at
x ≥ 43 (CH3.3's, CH4.3c2's and this task's): at the FMM's measured degrees the worst pair reaches
19.1, 25.9, 25.6, 27.2, 42.3 and 33.9 ε at 10⁻³ … 10⁻⁸. At the refit's degree, on S7's grid
k a = 45 … 256:

| k a | 10⁻³ | 10⁻⁴ | 10⁻⁵ | 10⁻⁶ | 10⁻⁷ | 10⁻⁸ |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 45.3 | 51: 2.96 | 55: 1.91 | 58: 2.12 | 60: 3.62 | 63: 3.27 | 65: 4.80 |
| 64.0 | 70: 3.13 | 74: 3.34 | 78: 1.95 | 81: 1.81 | 84: 1.23 | 86: 2.24 |
| 90.5 | 97: 2.27 | 102: 1.54 | 106: 2.40 | 109: 4.12 | 112: 4.10 | 115: 4.08 |
| 128.0 | 135: 3.01 | 140: 2.38 | 145: 2.04 | 148: 2.35 | 152: 2.76 | 155: 4.56 |
| 181.0 | 189: 2.32 | 194: 2.02 | 199: 1.88 | 204: 1.67 | 208: 1.56 | 211: 4.14 |
| 256.0 | 264: 1.54 | 271: 1.44 | 276: 2.20 | 281: 1.97 | 285: 2.55 | 290: 2.47 |

(each cell the degree, then the worst pair's error over ε). At most 4.8 ε (CH4.3c2's rule 4.6 ε,
k a = 128, 10⁻⁸): within the FMM's largest at every grid k a and ε, the plan's guard.

## 6. The FMM's errors

### 6.1 CH4.3c2's rule on today's trees

First, the rule as it stood, with W1's constant density under today's defaults (c_leaf 4 and the
diagonal form's time rule since CH4.3b2c and CH4.3e; CH4.3c2 verified at c_leaf 2), job 663690,
`helmholtz_boxes` at the samples, φ and ∇φ over ε:

| kD | ε | Rotation | default |
| ---: | --- | --- | --- |
| 250 | 10⁻⁸, 10⁻⁷, 10⁻⁶, 10⁻³ | 0.15, 0.59; 0.06, 0.27; 0.10, 0.47; 0.12, 0.15 | 0.11, 0.46; 0.05, 0.27; 0.10, 0.47; 0.12, 0.15 |
| 500 | 10⁻⁸, 10⁻⁷, 10⁻⁶, 10⁻³ | 0.17, **0.94**; 0.14, 0.48; 0.08, 0.48; 0.10, 0.21 | 0.16, 0.73; 0.13, 0.43; 0.08, 0.47; 0.10, 0.21 |
| 1000 | 10⁻⁸, 10⁻⁷, 10⁻⁶, 10⁻³ | 0.69, 0.70; 0.61, 0.35; 0.41, 0.35; 0.53, 0.30 | 0.66, 0.67; 0.61, 0.35; 0.41, 0.35; 0.53, 0.30 |

Within ε everywhere, at most 0.94 ε (∇φ at kD = 500, 10⁻⁸, Rotation; CH4.3c2's 0.94 ε at
c_leaf 2): the case whose level 5 (x = 13.5, p = 27) missed ε/2 in its scan (0.85 ε alone, §2).
The rule's misses of §3 are misses of the scans' ε/2, not of ε.

### 6.2 The refit

CH3.1c's matrix and CH3.3's held-out cases (job 663699; `ch4-3c3-tables.py matrix`): 174 runs
under the refit and the same 174 under CH4.3c2's rule, Rotation (near rules off) and the
library's defaults, against ε. Every run within ε; the largest error over ε, the larger of φ's and
∇φ's:

| M2L | refit | CH4.3c2's rule |
| --- | ---: | ---: |
| Rotation | 0.443 (W4's cube, kD = 100, 10⁻⁸) | 0.637 |
| defaults | 0.397 (W1, kD = 500, 10⁻³) | 0.353 |

The runs whose degrees changed are those with a level at x = 4–108 (kD ≤ 150 and W1 at
kD = 500): at 10⁻³ their degrees fall by 1–2 and their errors rise from 0.02–0.10 ε to 0.05–0.40 ε
(W1 at kD = 100: 0.05 → 0.26 ε; W4's cube at kD = 100: 0.05 → 0.33 ε), at 10⁻⁶ by 0.00–0.17 ε
(W4's cube: 0.06 → 0.23 ε). At 10⁻⁸ W1 at kD = 500 runs level 5 at 28 (CH4.3c2's 27): 0.64 →
0.42 ε under Rotation.

W1 at kD = 250, 500 and 1000 with the three densities (job 663700), φ and ∇φ over ε:

The larger of φ's and ∇φ's error over ε, the worst over Rotation and `helmholtz_boxes`' default,
at 10⁻⁸, 10⁻⁷, 10⁻⁶ and 10⁻³, under the refit and under CH4.3c2's rule (the same job, node and
trees):

| kD | density | refit: 10⁻⁸, 10⁻⁷, 10⁻⁶, 10⁻³ | CH4.3c2's rule |
| ---: | --- | --- | --- |
| 250 | random | 0.17, 0.15, 0.24, 0.32 | 0.24, 0.11, 0.21, 0.07 |
| 250 | plane wave | 0.19, 0.17, 0.24, 0.37 | 0.23, 0.11, 0.19, 0.08 |
| 250 | constant | 0.46, 0.39, 0.56, 0.54 | 0.59, 0.27, 0.47, 0.15 |
| 500 | random | 0.42, 0.30, 0.33, 0.40 | 0.64, 0.27, 0.31, 0.10 |
| 500 | plane wave | 0.15, 0.19, 0.22, 0.42 | 0.42, 0.17, 0.20, 0.10 |
| 500 | constant | 0.49, 0.49, 0.49, 0.63 | **0.94**, 0.48, 0.48, 0.21 |
| 1000 | random | 0.23, 0.18, 0.29, 0.40 | 0.33, 0.18, 0.28, 0.09 |
| 1000 | plane wave | 0.10, 0.15, 0.28, 0.46 | 0.32, 0.15, 0.27, 0.07 |
| 1000 | constant | 0.62, 0.37, 0.45, **0.65** | 0.70, 0.61, 0.41, 0.53 |

Every run within ε under both rules; the refit's worst 0.65 ε (φ of the constant density at
kD = 1000, 10⁻³), CH4.3c2's 0.94 ε (∇φ of the constant density at kD = 500, 10⁻⁸, Rotation, as in
§6.1). The refit trades error between the tolerances as its degrees do: at 10⁻³, where it is 1–2
degrees lower on the middle levels, its errors rise to 0.32–0.65 ε (CH4.3c2's 0.07–0.53 ε); at
10⁻⁸, where it raises the x = 13.5 level, they fall (kD = 500, constant: 0.94 → 0.49 ε).

## 7. Times

Job 663701 (node-c11b-009): W1 at kD = 1000 and 500, §9.1's random charges, `helmholtz_boxes`
at the samples, each case twice under each rule, the rules alternated (refit, CH4.3c2's, refit,
CH4.3c2's), Rotation and `helmholtz_boxes`' default (Auto, the near rules off), as CH4.3c2's §7.
The evaluation in seconds, the mean of the two runs (each pair within 2.5 %):

| kD | ε | M2L | p, levels 2–6, refit / CH4.3c2 | refit | CH4.3c2's | refit / CH4.3c2 |
| ---: | --- | --- | --- | ---: | ---: | ---: |
| 1000 | 10⁻³ | Rotation | 225, 115, 60, 32, 18 / 224, 115, 61, 34, 20 | 28.24 | 30.38 | 0.93 |
| 1000 | 10⁻³ | default | | 28.71 | 31.40 | 0.91 |
| 1000 | 10⁻⁶ | Rotation | 240, 128, 70, 40, 25 / 240, 128, 71, 41, 25 | 43.46 | 43.29 | 1.00 |
| 1000 | 10⁻⁶ | default | | 42.14 | 41.39 | 1.02 |
| 1000 | 10⁻⁸ | Rotation | 248, 134, 75, 44, 28 / 248, 134, 76, 44, 27 | 60.42 | 59.64 | 1.01 |
| 1000 | 10⁻⁸ | default | | 69.49 | 69.35 | 1.00 |
| 500 | 10⁻³ | Rotation | 115, 60, 32, 18, 11 / 115, 61, 34, 20, 13 | 4.20 | 4.80 | 0.87 |
| 500 | 10⁻³ | default | | 4.39 | 5.13 | 0.86 |
| 500 | 10⁻⁶ | Rotation | 128, 70, 40, 25, 16 / 128, 71, 41, 25, 16 | 7.11 | 7.23 | 0.98 |
| 500 | 10⁻⁶ | default | | 6.98 | 6.96 | 1.00 |
| 500 | 10⁻⁸ | Rotation | 134, 75, 44, 28, 24 / 134, 76, 44, 27, 24 | 11.31 | 11.16 | 1.01 |
| 500 | 10⁻⁸ | default | | 13.87 | 13.82 | 1.00 |

At 10⁻³ the refit lowers the dense middle levels by 1–2 degrees and the evaluation by 7–14 %; at
10⁻⁶ and 10⁻⁸ the changes (one degree down on levels 4–5, one up on the x = 13.5 level at 10⁻⁸)
cancel within the runs' spread. The errors of the same runs: the refit's larger at 10⁻³ (0.40 ε
against 0.09–0.10 ε) and within ε/2 everywhere, smaller at 10⁻⁸ under Rotation (kD = 500: 0.42
against 0.64 ε).

## 8. Questions for sign-off

Signed off 2026-11-22 by the user on #254, "1-3 as recommended": all three as asked, the refit
the default (neither alternative), proposal §3.5 amended with an Appendix A line.

1. The acceptance as read: the margin measured by level scans with W1's constant density at
   x = 13.5–216 (and 6.8) at kD = 250, 500 and 1000 (§2, §3); found wrong both ways (one short at
   x = 216 and x = 13.5, 10⁻⁸; up to 2 too many at x ≤ 54), so refitted on the worst density and
   verified as CH4.3c2 (§4–§6)?
2. Adopt the refit as the default (`DegreeRule::CALIBRATED`: C = 2.00, a = −0.55, b = −4.66,
   d = 2, capped by x + 0.82 (ln 1/ε)^{2/3} x^{1/3}; CH4.3c2's rule kept as
   `DegreeRule::CH4_3C2`)? The alternatives: keep CH4.3c2's rule (the FMM within ε everywhere,
   §6.1, but four scan points above it by one degree at ε/2), or CH4.3c2's rule with d = 6 and the
   cap 0.82 (covers every point, one degree more than CH4.3c2's wherever its cap does not bind, total excess 602).
3. Amend the proposal's §3.5 ("Calibrated") with the refit, the constant density's scans and this
   report, as "(amended at CH4.3c3)" with an Appendix A line?

## 9. Jobs, reproduce

Kathleen jobs (one node each, `singlenode`; 6.8 node-hours in all): builds 663685 and 663698
(node-c11b-043); scans 663686 (c11b-043: kD = 250, kD = 1000 level 6; 40 min), 663687 (c11a-074:
kD = 500; 36 min), 663688 (c11a-011: kD = 1000 levels 3–5; 2 h 22 min), 663689 (c11a-021: kD = 1000
level 2; 66 min), 663696 (c11a-007: x = 6.8; 6 min); CH4.3c2's rule on today's trees 663690
(c11a-007); the refit's matrix 663699 (c11b-043), W1 663700 (c11b-002), times 663701 (c11b-009).

```sh
cargo build --release -p tetraleaf-validate --example helmholtz_calibrate \
    --example helmholtz_accuracy --example helmholtz_boxes
tools/helmholtz/ch4-3c3-runs.sh scans <out>               # and scans-low
tools/helmholtz/ch4-3c3-tables.py scans <out>/*.md        # §2
tools/helmholtz/ch4-3c3-tables.py margin <out>/*.md       # §3
cargo run --release -p tetraleaf-validate --example helmholtz_calibrate -- \
    refit <CH4.3c2's scans>/*.md <out>/*.md               # §4
tools/helmholtz/ch4-3c2-pairs.py <points> \
    --rule refit=0.30,0.25,2.60,2.00,-0.55,-4.66,2,0.82   # §5
tools/helmholtz/ch4-3c3-runs.sh verify-w1 <out>           # §6.1 (RULE=ch4-3c2 after the change)
tools/helmholtz/ch4-3c3-runs.sh verify-matrix <out>       # §6.2; verify-w1-all, times (§7)
tools/helmholtz/ch4-3c3-tables.py matrix <out>
tools/helmholtz/ch4-3c2-tables.py w1 <out>/*.md
```

The raw outputs are in `~/projects/fmm-archive/kathleen/runs/ch4-3c3/` (CH4.3c2's scans in
`…/ch4-3c2/`).
