# Helmholtz CH4.8a: the constant density on W2, the error mapped and a remedy chosen

Task plan and record: `docs/helmholtz/tasks/CH4.8a.md`. Proposal: §10 CH4.4 (the row whose bound
this task restores), §3.5 (the degree rule), §9.1 (ε the relative ℓ2 tolerance; the worst box at most
100 ε within D1). Builds on the H4 gate's run (`docs/helmholtz/gates/H4.md` §9, Kathleen job 672173;
review `docs/reviews/2026-11-24-H4.md`, M15) and CH4.4b's box study under the defaults
(`helmholtz-ch4-4b.md`). Machine: Kathleen, one node per job (2 × 20 Cascade Lake cores, Xeon Gold
6248), `singlenode`, one rank × 40 threads, release, BLAS threads 1, node-c11b-004 and node-c11b-013
excluded. Jobs 672424 (`smoke`, `map`, `all`, `survey`; node-c11a-079), 672432 (`remedy`, `cost`;
node-c11a-021), 672477 (`repeat`; node-c11a-021), 672488 (`mid-eps`, added at the self-review;
node-c11b-012). Every launch exit 0 (208).
Raw outputs: `~/projects/fmm-archive/kathleen/runs/ch4-8a/` (one Markdown file per launch, its log
beside it). Script `tools/helmholtz/ch4-8a-runs.sh`; tables `tools/helmholtz/ch4-8a-tables.py` and
CH4.4b's `tools/helmholtz/ch4-4b-tables.py`.

## Summary

- **The cause.** W2's dense polar cap (half-angle 30°, ten times W1's density, about 42 % of the points
  and so of §9.1's samples) lies against the cube's z = 0 face (z ≤ 0.067 D), so its targets lie near a
  face of their boxes on many levels at once, and the cap adds two levels to the tree (levels 7–8 at
  kD = 500). With the constant density every level's truncation adds coherently there: ∇φ's RMS error
  grows from 0.38 ε at targets near no face to 1.7 ε near six and 4.3 ε near seven (kD = 500, 10⁻³). The
  M2L adds nothing: 49 runs with each non-RCR level (the form, its classes beside the threshold, the
  waves, Batched) on RCR change no worst box (1.000) and the ℓ2 by 0.999–1.005. It is the degree rule's
  truncation, per tolerance at different levels (§3):
  - **10⁻³, the ℓ2 miss (∇φ 1.28 ε at kD = 500):** the level at k a = 3.38, Batched at p = 7, where the
    high-frequency branch meets the plateau, carries about half of ∇φ's error (0.47, 0.50 of Σ|e|²;
    its worst pair 16 ε per pair); that level alone raised to 15 gives 0.82 ε (0.62 ε at kD = 250).
  - **10⁻⁶ (∇φ 0.99 ε, worst box 101 ε at kD = 500):** the level at k a = 6.77, p = 16, carries 0.42 of
    the ℓ2 (alone raised 0.34 ε); the worst box is levels 2–3's.
  - **10⁻⁸ (worst boxes 129.2 and 143.4 ε at kD = 250 and 500):** levels 2–4 at the tangent point, 0.78
    of ∇φ's error on level 2 at kD = 500; the levels at p_w = 32 under the waves add nothing.
  One level raised at a time moves ∇φ's worst box both ways (level 3 alone at kD = 500, 10⁻³: 56.6 →
  75.4 ε), since the levels' errors partly cancel there; every level raised together lowers it in
  every case (φ's worst box rose in two, by at most 1.1 ε: W1 at kD = 500, 10⁻³, 5.5 → 6.6 ε; W2 at
  kD = 500, 10⁻³, 20.2 → 20.9 ε).
- **Every miss found is W2's, and there are more than the gate's.** W2 at kD = 500 also misses at
  10⁻⁴ (∇φ 1.03 ε) and 10⁻⁷ (worst box 165.9 ε), at the tolerances the gate did not run (the
  self-review's `mid-eps` stage, §5). The cases of the verification set no record covered under the
  defaults, W3 with the constant density resolved at kD = 100 and at N = 10⁶ for kD = 200 and 500,
  hold at 10⁻³, 10⁻⁶ and 10⁻⁸: ℓ2 at most 0.51 ε, worst box at most 17.3 ε (§4). W1 with it at kD = 250
  and 500 holds (the gate: ℓ2 ≤ 0.64 ε, worst box ≤ 46.5 ε); W1 and W2 at kD ≤ 100 are CH4.4a's and
  CH4.4b's (within ε, ≤ 37.6 ε). Not run here: W1 and W3 at 10⁻⁴, 10⁻⁵ and 10⁻⁷, and W2 between kD =
  100 and 250; CH4.8b's verification covers them.
- **The remedy: the high-frequency branch one degree higher (`h1`).** p = max(plateau, high(x, ε) + 1),
  every level above the plateau one degree up, the plateau unchanged. On W2 with the constant density
  at kD = 250 and 500, at all six ε of 10⁻³–10⁻⁸, it gives ℓ2 at most 0.48 ε (from 1.28 ε) and worst
  boxes at most 76.3 ε (from 165.9 ε); on W1 with it ℓ2 at most 0.28 ε and worst boxes at most 17.6 ε (§5). Per pair (§6) it
  takes the worst pair on the levels of k a ≥ 13.5 from 1.4–6.4 ε to 0.7–2.3 ε, but on the levels that
  carry the ℓ2 misses only from 16 to 5.4 ε (k a = 3.38, 10⁻³) and from 32 to 20 ε (k a = 6.77, 10⁻⁶):
  per pair those stay well above ε, and only the FMM shows the misses closed. Raising the plateau too
  (`h1p1`) leaves `h1`'s worst boxes and its ℓ2 at 10⁻³ and 10⁻⁸ as they are and lowers ∇φ's ℓ2 at 10⁻⁶
  by a fifth (W2: 0.37 against 0.45 ε at kD = 500); two degrees (`h2`) lower ∇φ's ℓ2 to 0.32–0.79× and
  its worst box to 0.35–0.90× `h1`'s (at 10⁻³, 10⁻⁶, 10⁻⁸).
- **Its cost.** Build plus one evaluation 1.004–1.058× the default's at kD = 500 in single launches,
  inside the launches' own spread of about ±5 % (§7); the A B B A repeat, the medians of four
  launches each, gives 0.989–1.047× (geometric mean 1.017).
- **For sign-off (§8):** `h1` as the rule CH4.8b makes the default and verifies; the plateau kept.

## 1. Method

- **The cases.** W2 with the constant density (q = 1 at every point) at kD = 250 (N = 797,211) and
  500 (N = 3,188,840), the gate's misses; W1 with it at kD = 250 and 500 (the remedies' second
  workload); W3 with it resolved at kD = 100 and at N = 10⁶ for kD = 200 and 500 (the survey). ε =
  10⁻³, 10⁻⁶, 10⁻⁸, and for W2 also 10⁻⁴, 10⁻⁵, 10⁻⁷ under the default, `h1` and `h1p1` (`mid-eps`,
  added at the self-review: the H4 gate had run three ε only); one complex charge vector, gradients
  on.
- **The FMM.** The library's defaults (`--m2l defaults`, `Case::defaults`: CH4.3c3's degree rule, the
  waves under Auto, the diagonal form's time rule, c_leaf 4, the V, W and X rules), and for the remedies
  the same with another `DegreeRule` (`--rule γ,α,β,C,a,b,d[,cap]`).
- **The reference.** `helmholtz_boxes` (CH4.0b, CH4.4b): the f64 direct sum at §9.1's 1,000 samples,
  two targets per box of the deepest level and every point of the 8 worst boxes of φ and of ∇φ and of
  the boxes at the centres of the domain's faces (`--refine 8 --refine-faces`); 1.3·10⁵ targets at
  kD = 250, 5.3·10⁵ at 500. The sample is the H4 gate's, read from its cache (W1, W2), so every base run
  equals the gate's to the digits printed (W2 at kD = 500: 1.28·10⁻³, 56.63, 101.02 and 143.40 ε).
- **The mapping** (`map`), as CH4.4b: each level raised to the rule's degree at ε/1000 and on RCR, the
  others the defaults' (the level's part of the error; the parts' rest at most 0.024 of Σ|e|² at 10⁻³,
  0.003 at 10⁻⁶, 0.001 at 10⁻⁸), each non-RCR level on RCR, the form's classes by decade below the
  bound's threshold on RCR, and the level of the worst ∇φ box's largest part raised by 1 and 2.
- **Every level raised** (`all`): the levels with V pairs one and two degrees up at once, under the
  setting's M2L (`--raise-levels v --raise-together`): the shape of a rule's remedy before a rule.
- **The candidate rules** (`remedy`), each as constants of `DegreeRule` with the default's plateau
  and branch (γ = 0.30, α = 0.25, β = 2.60, C = 2.00, a = −0.55, b = −4.66):
  - `h1`: d = 3, no cap; `h2`: d = 4, no cap (the branch one and two degrees up);
  - `h1p1`, `h2p1`: the same with β = 3.60 (the plateau one degree up).
  `h1` is the default's capped branch + 1 at every x ≤ 400 and 10⁻⁸ ≤ ε ≤ 10⁻³ but at x = 7.1–13.5
  and ε ≤ 2·10⁻⁷, where it is one degree more still at isolated points (635 of 8·10⁶ on a grid of 0.01
  in x and 200 ε per five decades: the cap binds there and d = 3 without it lies above the cap + 1). At
  every one of them the waves' resolving degree p_w (28 or 32, CH6.1g) lies above both, so on a level
  with a quadrature, as every level there in these runs, the degrees as run are the same. CH4.8b
  implements the rule as "high + 1" with the cap kept (§8, question 2).
- **The cost** (`cost`, `repeat`): one launch per setting (CH4.3b2c §6), `--study single`, the run's
  build and evaluation as the example times them (one evaluation), on W2 and W1 with the constant
  density and W3 at N = 10⁶, kD = 500; `cost` each rule once in turn per case and ε, `repeat` the
  default and `h1` in the order A B B A A B B A.

## 2. The mapping: levels, M2L and the worst boxes

The levels' parts on W2 (share of Σ|e|² over the reference's targets; the shares overlap, since the
parts are not orthogonal, and their sum is the error to the rest printed):

| kD | ε | the largest part of ∇φ (k a, p, M2L) | the others above 0.1 |
| --- | --- | --- | --- |
| 250 | 10⁻³ | level 6 (3.38, 7, Batched) 0.50 | levels 2–5 (54.1–6.8) 0.16–0.25 each |
| 250 | 10⁻⁶ | level 5 (6.77, 16, the form) 0.34 | levels 2, 3 0.19, 0.13 |
| 250 | 10⁻⁸ | level 2 (54.13, 75, the form) 0.46 | level 3 0.42 |
| 500 | 10⁻³ | level 7 (3.38, 7, Batched) 0.47 | levels 5, 6 0.14, 0.15 |
| 500 | 10⁻⁶ | level 6 (6.77, 16, the form) 0.42 | levels 2, 3 0.16, 0.10 |
| 500 | 10⁻⁸ | level 2 (108.25, 134, the form) 0.78 | levels 3, 4 0.35, 0.26 |

- **The M2L.** 49 runs (37 levels on RCR, 12 class decades on RCR): every worst box ÷ the default's
  1.000, ℓ2 0.999–1.005, the largest share of ∇φ 0.0009. No error is the form's, its classes', the
  waves' or Batched's.
- **One level raised** (to the rule's degree at ε/1000, on RCR): at kD = 500, 10⁻³ level 7 alone takes
  ∇φ's ℓ2 from 1.28 to 0.82 ε and leaves the worst box (56.6 ε); level 2 alone takes the box to 32.4 ε,
  level 3 alone to 75.4 ε. At 10⁻⁶ level 6 alone takes ℓ2 from 0.99 to 0.34 ε, level 2 alone the box
  from 101 to 53 ε; at 10⁻⁸ level 4 alone the box from 143 to 78 ε, level 3 alone to 156 ε. The
  `--raise-by` runs on the worst box's level (+1, +2 under the setting's M2L) are as mixed (§9).
- **The worst boxes** are the tangent point's at z = 0 in every W2 run, cells (127–128, 127–128, 0) of
  the deepest level, every level's target near a box corner (|v|/r 1.2–1.73), as the gate found. The
  error grows with the levels on whose face a target lies (kD = 500):

| levels on a face | targets | ∇φ RMS ÷ ε at 10⁻³ | at 10⁻⁶ | at 10⁻⁸ |
| ---: | ---: | ---: | ---: | ---: |
| 0 | 104,980 | 0.38 | 0.18 | 0.22 |
| 3 | 86,585 | 1.22 | 0.75 | 0.26 |
| 6 | 13,454 | 1.73 | 1.47 | 0.49 |
| 7 | 6,836 | 4.28 | 5.56 | 6.74 |

The cap's targets lie within 0.067 D of z = 0, so near the face of their box on the coarse levels; the
ℓ2 at §9.1's samples (1.28 ε) exceeds the figure over the reference's 5.3·10⁵ targets (1.07 ε, two
per box of the deepest level and the refined boxes) because the samples follow the points, 42 % of
them in the cap.

## 3. The cause

The degree rule was calibrated for ε/2 in ℓ2 on W1, W3 and W4 (CH3.3's plateau) and refitted on W1's
constant density (CH4.3c3's branch); W2's constant density was never in it. Its tree is the deepest of
D1's workloads, its points crowd one face of the domain, and the constant density makes the levels'
truncation errors add. No single level carries the error: which level dominates depends on ε (the
plateau's edge at 10⁻³, k a ≈ 6.8 at 10⁻⁶, the top levels at 10⁻⁸), as the rule's per-pair margins do
(the worst pair 16 ε at k a = 3.38, p = 7, 10⁻³; 32 ε at 6.77, p = 16, 10⁻⁶; 3–5 ε on the top levels
at 10⁻⁸; §6). Of the roadmap's suspects: the degree rule on the dense cap's levels and the tangent-point
placement, together; not the form, the waves or Batched at p = 7 as an M2L (its RCR at p = 7 is the
same), though Batched's level is the one whose degree matters at 10⁻³.

## 4. The survey: W3 with the constant density

Under the defaults (`--study levels,m2l --all-levels`):

| case | ε | §9.1 φ, ∇φ ÷ ε | worst box φ, ∇φ ÷ ε |
| --- | --- | --- | --- |
| w3-kd100-constant | 1e-3 | 0.29, 0.36 | 7.0, 3.0 |
| w3-kd100-constant | 1e-6 | 0.14, 0.28 | 4.7, 4.2 |
| w3-kd100-constant | 1e-8 | 0.06, 0.15 | 4.5, 3.6 |
| w3:1000000-kd200-constant | 1e-3 | 0.33, 0.46 | 6.6, 10.6 |
| w3:1000000-kd200-constant | 1e-6 | 0.16, 0.28 | 6.6, 12.3 |
| w3:1000000-kd200-constant | 1e-8 | 0.15, 0.29 | 4.9, 11.3 |
| w3:1000000-kd500-constant | 1e-3 | 0.41, 0.51 | 9.4, 17.3 |
| w3:1000000-kd500-constant | 1e-6 | 0.15, 0.20 | 10.0, 12.7 |
| w3:1000000-kd500-constant | 1e-8 | 0.19, 0.26 | 10.1, 17.3 |

No miss; the M2L on RCR moves no worst box by more than 1.3 % (the waves on RCR at p_w, one-sided as
CH4.4b reads it). With W1 at kD = 250 and 500 (the gate) and W1, W2 at kD ≤ 100 and W3 at 10, 30
(CH4.4a, CH4.4b), every workload and kD of CH4.8's verification set with the constant density is now
measured under the defaults at 10⁻³, 10⁻⁶ and 10⁻⁸, and W2 at kD = 250 and 500 at all six ε; the
misses are W2's at kD = 250 and 500 (at 10⁻⁸; and at kD = 500 at 10⁻³, 10⁻⁴, 10⁻⁶ and 10⁻⁷). Not run:
W1 and W3 at 10⁻⁴, 10⁻⁵, 10⁻⁷ and W2 at 100 < kD < 250, CH4.8b's.

## 5. The candidate rules on W1 and W2 with the constant density

Bold: §9.1's ℓ2 above ε, or a worst box above 100 ε. The default's rows are the `map` runs (W2) and
the H4 gate's (W1).

| case | ε | rule | p by level | §9.1 φ, ∇φ ÷ ε | worst box φ, ∇φ ÷ ε |
| --- | --- | --- | --- | --- | --- |
| w1-kd250-constant | 1e-3 | default | 0×2, 60, 32, 18, 11 | 0.21, 0.54 | 2.6, 10.6 |
| w1-kd250-constant | 1e-3 | h1 | 0×2, 61, 33, 19, 12 | 0.09, 0.19 | 1.1, 3.1 |
| w1-kd250-constant | 1e-3 | h1p1 | 0×2, 61, 33, 19, 12 | 0.09, 0.19 | 1.1, 3.1 |
| w1-kd250-constant | 1e-3 | h2 | 0×2, 62, 34, 20, 13 | 0.04, 0.09 | 1.3, 1.2 |
| w1-kd250-constant | 1e-3 | h2p1 | 0×2, 62, 34, 20, 13 | 0.04, 0.09 | 1.3, 1.2 |
| w1-kd250-constant | 1e-6 | default | 0×2, 70, 40, 25, 16 | 0.17, 0.51 | 12.2, 33.7 |
| w1-kd250-constant | 1e-6 | h1 | 0×2, 71, 41, 26, 17 | 0.06, 0.20 | 6.3, 15.9 |
| w1-kd250-constant | 1e-6 | h1p1 | 0×2, 71, 41, 26, 17 | 0.06, 0.20 | 6.3, 15.9 |
| w1-kd250-constant | 1e-6 | h2 | 0×2, 72, 42, 27, 18 | 0.03, 0.09 | 2.2, 5.8 |
| w1-kd250-constant | 1e-6 | h2p1 | 0×2, 72, 42, 27, 18 | 0.03, 0.09 | 2.2, 5.8 |
| w1-kd250-constant | 1e-8 | default | 0×2, 75, 44, 32×2 | 0.21, 0.44 | 11.7, 46.5 |
| w1-kd250-constant | 1e-8 | h1 | 0×2, 76, 45, 32×2 | 0.05, 0.19 | 8.6, 17.6 |
| w1-kd250-constant | 1e-8 | h1p1 | 0×2, 76, 45, 32×2 | 0.05, 0.19 | 8.6, 17.6 |
| w1-kd250-constant | 1e-8 | h2 | 0×2, 77, 46, 32×2 | 0.03, 0.06 | 1.5, 7.3 |
| w1-kd250-constant | 1e-8 | h2p1 | 0×2, 77, 46, 32×2 | 0.03, 0.06 | 1.5, 7.3 |
| w1-kd500-constant | 1e-3 | default | 0×2, 115, 60, 32, 18, 11 | 0.22, 0.64 | 5.5, 23.4 |
| w1-kd500-constant | 1e-3 | h1 | 0×2, 116, 61, 33, 19, 12 | 0.09, 0.28 | 6.6, 10.1 |
| w1-kd500-constant | 1e-3 | h1p1 | 0×2, 116, 61, 33, 19, 12 | 0.09, 0.28 | 6.6, 10.1 |
| w1-kd500-constant | 1e-3 | h2 | 0×2, 117, 62, 34, 20, 13 | 0.05, 0.13 | 3.4, 8.9 |
| w1-kd500-constant | 1e-3 | h2p1 | 0×2, 117, 62, 34, 20, 13 | 0.05, 0.13 | 3.4, 8.9 |
| w1-kd500-constant | 1e-6 | default | 0×2, 128, 70, 40, 25, 16 | 0.10, 0.44 | 11.2, 33.6 |
| w1-kd500-constant | 1e-6 | h1 | 0×2, 129, 71, 41, 26, 17 | 0.04, 0.15 | 4.3, 11.7 |
| w1-kd500-constant | 1e-6 | h1p1 | 0×2, 129, 71, 41, 26, 17 | 0.04, 0.15 | 4.3, 11.7 |
| w1-kd500-constant | 1e-6 | h2 | 0×2, 130, 72, 42, 27, 18 | 0.02, 0.07 | 2.2, 7.1 |
| w1-kd500-constant | 1e-6 | h2p1 | 0×2, 130, 72, 42, 27, 18 | 0.02, 0.07 | 2.2, 7.1 |
| w1-kd500-constant | 1e-8 | default | 0×2, 134, 75, 44, 32×2 | 0.10, 0.34 | 13.6, 33.7 |
| w1-kd500-constant | 1e-8 | h1 | 0×2, 135, 76, 45, 32×2 | 0.04, 0.12 | 3.3, 12.3 |
| w1-kd500-constant | 1e-8 | h1p1 | 0×2, 135, 76, 45, 32×2 | 0.04, 0.12 | 3.3, 12.3 |
| w1-kd500-constant | 1e-8 | h2 | 0×2, 136, 77, 46, 32×2 | 0.02, 0.05 | 2.2, 4.8 |
| w1-kd500-constant | 1e-8 | h2p1 | 0×2, 136, 77, 46, 32×2 | 0.02, 0.05 | 2.1, 4.8 |
| w2-kd250-constant | 1e-3 | default | 0×2, 60, 32, 18, 11, 7×2 | 0.18, 0.98 | 13.4, 44.7 |
| w2-kd250-constant | 1e-3 | h1 | 0×2, 61, 33, 19, 12, 8, 7 | 0.07, 0.40 | 8.0, 19.4 |
| w2-kd250-constant | 1e-3 | h1p1 | 0×2, 61, 33, 19, 12, 8×2 | 0.07, 0.39 | 8.0, 19.4 |
| w2-kd250-constant | 1e-3 | h2 | 0×2, 62, 34, 20, 13, 9, 7 | 0.03, 0.18 | 6.6, 9.9 |
| w2-kd250-constant | 1e-3 | h2p1 | 0×2, 62, 34, 20, 13, 9, 8 | 0.03, 0.14 | 6.6, 9.9 |
| w2-kd250-constant | 1e-4 | default | 0×2, 64, 35, 21, 13, 9×2 | 0.10, 0.87 | 28.3, 66.8 |
| w2-kd250-constant | 1e-4 | h1 | 0×2, 65, 36, 22, 14, 10, 9 | 0.04, 0.34 | 12.7, 32.7 |
| w2-kd250-constant | 1e-4 | h1p1 | 0×2, 65, 36, 22, 14, 10×2 | 0.04, 0.28 | 12.7, 32.7 |
| w2-kd250-constant | 1e-5 | default | 0×2, 67, 38, 23, 15, 12×2 | 0.08, 0.41 | 30.8, 69.2 |
| w2-kd250-constant | 1e-5 | h1 | 0×2, 68, 39, 24, 16, 12×2 | 0.04, 0.29 | 14.2, 30.7 |
| w2-kd250-constant | 1e-5 | h1p1 | 0×2, 68, 39, 24, 16, 13×2 | 0.04, 0.19 | 14.2, 30.7 |
| w2-kd250-constant | 1e-6 | default | 0×2, 70, 40, 25, 16, 15×2 | 0.09, 0.53 | 38.7, 89.6 |
| w2-kd250-constant | 1e-6 | h1 | 0×2, 71, 41, 26, 17, 15×2 | 0.04, 0.29 | 15.9, 40.7 |
| w2-kd250-constant | 1e-6 | h1p1 | 0×2, 71, 41, 26, 17, 16×2 | 0.04, 0.23 | 15.9, 40.7 |
| w2-kd250-constant | 1e-6 | h2 | 0×2, 72, 42, 27, 18, 15×2 | 0.02, 0.23 | 6.5, 28.6 |
| w2-kd250-constant | 1e-6 | h2p1 | 0×2, 72, 42, 27, 18, 16×2 | 0.02, 0.14 | 6.4, 16.7 |
| w2-kd250-constant | 1e-7 | default | 0×2, 73, 42, 28×4 | 0.07, 0.25 | 23.8, 97.4 |
| w2-kd250-constant | 1e-7 | h1 | 0×2, 74, 43, 28×4 | 0.02, 0.08 | 13.5, 33.7 |
| w2-kd250-constant | 1e-7 | h1p1 | 0×2, 74, 43, 28×4 | 0.02, 0.08 | 13.5, 33.6 |
| w2-kd250-constant | 1e-8 | default | 0×2, 75, 44, 32×4 | 0.08, 0.26 | 51.2, **129.2** |
| w2-kd250-constant | 1e-8 | h1 | 0×2, 76, 45, 32×4 | 0.03, 0.11 | 17.3, 49.9 |
| w2-kd250-constant | 1e-8 | h1p1 | 0×2, 76, 45, 32×4 | 0.03, 0.11 | 17.3, 49.9 |
| w2-kd250-constant | 1e-8 | h2 | 0×2, 77, 46, 32×4 | 0.01, 0.04 | 6.2, 17.6 |
| w2-kd250-constant | 1e-8 | h2p1 | 0×2, 77, 46, 32×4 | 0.01, 0.04 | 6.2, 17.6 |
| w2-kd500-constant | 1e-3 | default | 0×2, 115, 60, 32, 18, 11, 7×2 | 0.18, **1.28** | 20.2, 56.6 |
| w2-kd500-constant | 1e-3 | h1 | 0×2, 116, 61, 33, 19, 12, 8, 7 | 0.08, 0.48 | 20.9, 24.0 |
| w2-kd500-constant | 1e-3 | h1p1 | 0×2, 116, 61, 33, 19, 12, 8×2 | 0.08, 0.47 | 20.9, 24.0 |
| w2-kd500-constant | 1e-3 | h2 | 0×2, 117, 62, 34, 20, 13, 9, 7 | 0.04, 0.21 | 10.1, 21.7 |
| w2-kd500-constant | 1e-3 | h2p1 | 0×2, 117, 62, 34, 20, 13, 9, 8 | 0.04, 0.18 | 10.1, 21.7 |
| w2-kd500-constant | 1e-4 | default | 0×2, 120, 64, 35, 21, 13, 9×2 | 0.12, **1.03** | 18.1, 60.5 |
| w2-kd500-constant | 1e-4 | h1 | 0×2, 121, 65, 36, 22, 14, 10, 9 | 0.04, 0.39 | 18.3, 30.8 |
| w2-kd500-constant | 1e-4 | h1p1 | 0×2, 121, 65, 36, 22, 14, 10×2 | 0.04, 0.33 | 18.3, 30.8 |
| w2-kd500-constant | 1e-5 | default | 0×2, 124, 67, 38, 23, 15, 12×2 | 0.08, 0.56 | 42.2, 94.6 |
| w2-kd500-constant | 1e-5 | h1 | 0×2, 125, 68, 39, 24, 16, 12×2 | 0.04, 0.33 | 25.9, 50.7 |
| w2-kd500-constant | 1e-5 | h1p1 | 0×2, 125, 68, 39, 24, 16, 13×2 | 0.04, 0.22 | 25.9, 50.7 |
| w2-kd500-constant | 1e-6 | default | 0×2, 128, 70, 40, 25, 16, 15×2 | 0.11, 0.99 | 46.1, **101.0** |
| w2-kd500-constant | 1e-6 | h1 | 0×2, 129, 71, 41, 26, 17, 15×2 | 0.05, 0.45 | 24.4, 47.8 |
| w2-kd500-constant | 1e-6 | h1p1 | 0×2, 129, 71, 41, 26, 17, 16×2 | 0.04, 0.37 | 24.4, 47.8 |
| w2-kd500-constant | 1e-6 | h2 | 0×2, 130, 72, 42, 27, 18, 15×2 | 0.03, 0.32 | 9.1, 33.7 |
| w2-kd500-constant | 1e-6 | h2p1 | 0×2, 130, 72, 42, 27, 18, 16×2 | 0.02, 0.20 | 9.1, 27.3 |
| w2-kd500-constant | 1e-7 | default | 0×2, 131, 73, 42, 28×4 | 0.07, 0.29 | 68.6, **165.9** |
| w2-kd500-constant | 1e-7 | h1 | 0×2, 132, 74, 43, 28×4 | 0.03, 0.11 | 25.0, 76.3 |
| w2-kd500-constant | 1e-7 | h1p1 | 0×2, 132, 74, 43, 28×4 | 0.03, 0.11 | 25.0, 76.5 |
| w2-kd500-constant | 1e-8 | default | 0×2, 134, 75, 44, 32×4 | 0.08, 0.29 | 45.1, **143.4** |
| w2-kd500-constant | 1e-8 | h1 | 0×2, 135, 76, 45, 32×4 | 0.03, 0.12 | 23.4, 45.5 |
| w2-kd500-constant | 1e-8 | h1p1 | 0×2, 135, 76, 45, 32×4 | 0.03, 0.12 | 23.4, 45.5 |
| w2-kd500-constant | 1e-8 | h2 | 0×2, 136, 77, 46, 32×4 | 0.02, 0.05 | 9.5, 20.4 |
| w2-kd500-constant | 1e-8 | h2p1 | 0×2, 136, 77, 46, 32×4 | 0.02, 0.05 | 9.5, 20.4 |

- `h1` brings every case inside both bounds: ℓ2 at most 0.48 ε (W2 at kD = 500, 10⁻³, just inside the
  rule's calibration target of ε/2, §3.5), worst box at most 76.3 ε (W2 at kD = 500, 10⁻⁷; at most
  50.7 ε at the other ε). The default misses at 10⁻⁴ (∇φ 1.03 ε) and 10⁻⁷ (165.9 ε) too, which the H4
  gate's three ε did not show. Where `h1p1` raises the same levels as "every level
  with V pairs raised by one" (§9; W2 at kD = 500, 10⁻³ and 10⁻⁶) the two give the same figures.
- `h1p1` (the plateau too) gives `h1`'s worst boxes, and its ℓ2 at 10⁻³ and 10⁻⁸ to 0.01 ε; at 10⁻⁶
  it lowers ∇φ's ℓ2 on W2 (0.23 and 0.37 against 0.29 and 0.45 ε at kD = 250 and 500): under `h1` the
  plateau levels (k a = 3.38 and 1.69, p = 15) carry about a third of ∇φ's Σ|e|² there (1 − (0.37/0.45)²
  = 0.32 at kD = 500, 0.37 at kD = 250), more than the 0.02–0.09 per level of the default's runs.
- `h2` lowers ∇φ's ℓ2 to 0.32–0.79× and the worst box to 0.35–0.90× `h1`'s (ℓ2 at most 0.32 ε, worst
  box at most 33.7 ε).

## 6. Pair evidence

Spike S7's worst pair over the 316 offsets (`spikes/truncation/results/c1_curves.csv`, read as
CH4.3c2's pair evidence reads it, `tools/helmholtz/ch4-4b-tables.py pairs`) at the default's and
`h1`'s degrees, over ε, at the levels of D1's trees (x = Re(k) a_l):

| ε | x = 216.5 | 108.25 | 54.13 | 27.06 | 13.53 | 6.77 | 3.38 | 1.69 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 10⁻³ | 1.6 → 1.0 | 2.6 → 1.5 | 3.0 → 1.7 | 4.8 → 2.2 | 6.4 → 2.3 | 6.3 → 1.8 | 16 → 5.4 | 10 |
| 10⁻⁴ | 1.4 → 0.94 | 2.1 → 1.3 | 2.6 → 1.3 | 4.2 → 1.8 | 2.7 → 0.86 | 5.0 → 1.8 | 19 → 7.5 | 24 |
| 10⁻⁵ | 1.6 → 1.0 | 2.7 → 1.5 | 3.0 → 1.3 | 2.7 → 0.98 | 2.5 → 0.77 | 7.0 → 3.2 | 18 | 37 |
| 10⁻⁶ | 2.5 → 1.5 | 2.3 → 1.2 | 2.6 → 1.2 | 3.7 → 1.3 | 2.6 → 0.97 | 32 → 20 | 44 | 130 |
| 10⁻⁷ | 1.6 → 1.0 | 3.4 → 1.8 | 2.0 → 0.71 | 4.9 → 1.6 | 9.7 → 4.9 | 39 | 44 | 150 |
| 10⁻⁸ | 3.9 → 2.3 | 4.2 → 1.7 | 3.2 → 1.3 | 5.5 → 1.8 | 31 → 17 (as run 4.2) | 71 | 90 | 270 |

(One figure where `h1` leaves the level on the plateau.) On the levels of x ≥ 13.5 the rule's worst
pair falls from 1.4–6.4 ε to 0.7–2.3 ε. On the levels that carry the ℓ2 misses it falls less: 16 → 5.4
ε (x = 3.38, 10⁻³), 32 → 20 ε (x = 6.77, 10⁻⁶), and at 10⁻⁴, 10⁻⁵ and 10⁻⁷ 19 → 7.5, 7.0 → 3.2 and 9.7 →
4.9 ε; the plateau's pairs stay 10–270 ε, by its calibration on the FMM's ℓ2, not per pair (CH3.3;
half of S7's p*). Below 10⁻⁶ a level of x ≤ 13.5 with a quadrature keeps the waves' resolving degree
p_w = 28 or 32 (CH6.1g, CH4.3b2a) whatever M2L the time rule gives it, above both rules' degrees: at
kD = 500, 10⁻⁸ level 5 (x = 13.53) runs the form at p = 32 under both rules, its worst pair 4.2 ε, and
the levels below run the waves. The two measures do not say the same: per pair `h1` is never worse
and is better on every level it raises, but the levels that carry W2's ℓ2 misses stay 5–20 ε per
pair under it; only the FMM shows those misses closed (§5). That is how the plateau and the branch's
low end were calibrated (on the FMM's ℓ2, not per pair), and why the FMM's verification over CH4.8's
whole set is CH4.8b's acceptance.

## 7. Cost

One launch per setting, build plus one evaluation over the default's (the default's in seconds):

| case | ε | rule | p by level | build s | evaluate s | build + evaluate ÷ default |
| --- | --- | --- | --- | ---: | ---: | ---: |
| w1-kd500-constant | 1e-3 | default | 0×2, 115, 60, 32, 18, 11 | 2.60 | 2.60 | 1.000 |
| w1-kd500-constant | 1e-3 | h1 | 0×2, 116, 61, 33, 19, 12 | 2.66 | 2.71 | 1.033 |
| w1-kd500-constant | 1e-3 | h1p1 | 0×2, 116, 61, 33, 19, 12 | 2.67 | 2.78 | 1.048 |
| w1-kd500-constant | 1e-3 | h2 | 0×2, 117, 62, 34, 20, 13 | 2.66 | 2.83 | 1.056 |
| w1-kd500-constant | 1e-3 | h2p1 | 0×2, 117, 62, 34, 20, 13 | 2.65 | 2.75 | 1.038 |
| w1-kd500-constant | 1e-6 | default | 0×2, 128, 70, 40, 25, 16 | 2.84 | 3.69 | 1.000 |
| w1-kd500-constant | 1e-6 | h1 | 0×2, 129, 71, 41, 26, 17 | 2.82 | 3.78 | 1.011 |
| w1-kd500-constant | 1e-6 | h1p1 | 0×2, 129, 71, 41, 26, 17 | 2.84 | 3.99 | 1.046 |
| w1-kd500-constant | 1e-6 | h2 | 0×2, 130, 72, 42, 27, 18 | 2.88 | 3.88 | 1.035 |
| w1-kd500-constant | 1e-6 | h2p1 | 0×2, 130, 72, 42, 27, 18 | 2.90 | 3.94 | 1.047 |
| w1-kd500-constant | 1e-8 | default | 0×2, 134, 75, 44, 32×2 | 3.01 | 6.85 | 1.000 |
| w1-kd500-constant | 1e-8 | h1 | 0×2, 135, 76, 45, 32×2 | 3.05 | 6.89 | 1.008 |
| w1-kd500-constant | 1e-8 | h1p1 | 0×2, 135, 76, 45, 32×2 | 3.03 | 6.83 | 1.000 |
| w1-kd500-constant | 1e-8 | h2 | 0×2, 136, 77, 46, 32×2 | 3.05 | 7.07 | 1.026 |
| w1-kd500-constant | 1e-8 | h2p1 | 0×2, 136, 77, 46, 32×2 | 3.07 | 7.41 | 1.063 |
| w2-kd500-constant | 1e-3 | default | 0×2, 115, 60, 32, 18, 11, 7×2 | 4.19 | 3.30 | 1.000 |
| w2-kd500-constant | 1e-3 | h1 | 0×2, 116, 61, 33, 19, 12, 8, 7 | 4.25 | 3.39 | 1.020 |
| w2-kd500-constant | 1e-3 | h1p1 | 0×2, 116, 61, 33, 19, 12, 8×2 | 4.13 | 3.37 | 1.001 |
| w2-kd500-constant | 1e-3 | h2 | 0×2, 117, 62, 34, 20, 13, 9, 7 | 4.23 | 3.56 | 1.040 |
| w2-kd500-constant | 1e-3 | h2p1 | 0×2, 117, 62, 34, 20, 13, 9, 8 | 4.25 | 3.64 | 1.053 |
| w2-kd500-constant | 1e-6 | default | 0×2, 128, 70, 40, 25, 16, 15×2 | 4.46 | 5.07 | 1.000 |
| w2-kd500-constant | 1e-6 | h1 | 0×2, 129, 71, 41, 26, 17, 15×2 | 4.43 | 5.19 | 1.009 |
| w2-kd500-constant | 1e-6 | h1p1 | 0×2, 129, 71, 41, 26, 17, 16×2 | 4.44 | 5.39 | 1.031 |
| w2-kd500-constant | 1e-6 | h2 | 0×2, 130, 72, 42, 27, 18, 15×2 | 4.49 | 5.43 | 1.041 |
| w2-kd500-constant | 1e-6 | h2p1 | 0×2, 130, 72, 42, 27, 18, 16×2 | 4.48 | 5.39 | 1.036 |
| w2-kd500-constant | 1e-8 | default | 0×2, 134, 75, 44, 32×4 | 4.58 | 10.09 | 1.000 |
| w2-kd500-constant | 1e-8 | h1 | 0×2, 135, 76, 45, 32×4 | 4.70 | 10.51 | 1.037 |
| w2-kd500-constant | 1e-8 | h1p1 | 0×2, 135, 76, 45, 32×4 | 4.58 | 9.72 | 0.975 |
| w2-kd500-constant | 1e-8 | h2 | 0×2, 136, 77, 46, 32×4 | 4.58 | 9.92 | 0.988 |
| w2-kd500-constant | 1e-8 | h2p1 | 0×2, 136, 77, 46, 32×4 | 4.64 | 10.27 | 1.016 |
| w3:1000000-kd500-constant | 1e-3 | default | 0×2, 115, 60, 32, 18 | 3.00 | 3.73 | 1.000 |
| w3:1000000-kd500-constant | 1e-3 | h1 | 0×2, 116, 61, 33, 19 | 2.93 | 4.19 | 1.058 |
| w3:1000000-kd500-constant | 1e-3 | h1p1 | 0×2, 116, 61, 33, 19 | 2.89 | 3.58 | 0.961 |
| w3:1000000-kd500-constant | 1e-3 | h2 | 0×2, 117, 62, 34, 20 | 2.90 | 4.16 | 1.049 |
| w3:1000000-kd500-constant | 1e-3 | h2p1 | 0×2, 117, 62, 34, 20 | 2.81 | 3.65 | 0.960 |
| w3:1000000-kd500-constant | 1e-6 | default | 0×2, 128, 70, 40, 25 | 3.06 | 4.83 | 1.000 |
| w3:1000000-kd500-constant | 1e-6 | h1 | 0×2, 129, 71, 41, 26 | 3.05 | 4.87 | 1.004 |
| w3:1000000-kd500-constant | 1e-6 | h1p1 | 0×2, 129, 71, 41, 26 | 3.06 | 5.41 | 1.074 |
| w3:1000000-kd500-constant | 1e-6 | h2 | 0×2, 130, 72, 42, 27 | 3.10 | 5.55 | 1.096 |
| w3:1000000-kd500-constant | 1e-6 | h2p1 | 0×2, 130, 72, 42, 27 | 3.14 | 5.41 | 1.084 |
| w3:1000000-kd500-constant | 1e-8 | default | 0×2, 134, 75, 44, 32 | 3.06 | 6.76 | 1.000 |
| w3:1000000-kd500-constant | 1e-8 | h1 | 0×2, 135, 76, 45, 32 | 3.06 | 7.02 | 1.026 |
| w3:1000000-kd500-constant | 1e-8 | h1p1 | 0×2, 135, 76, 45, 32 | 3.08 | 6.72 | 0.998 |
| w3:1000000-kd500-constant | 1e-8 | h2 | 0×2, 136, 77, 46, 32 | 3.09 | 7.20 | 1.048 |
| w3:1000000-kd500-constant | 1e-8 | h2p1 | 0×2, 136, 77, 46, 32 | 3.13 | 7.20 | 1.052 |

The single launches' spread is about ±5 %: `h1` and `h1p1` run the same degrees on W1, W3 and at W2's
10⁻⁸, and differ by up to 0.097 (W3, 10⁻³: 1.058 and 0.961); `h1` alone 1.004–1.058×. The repeat, A B B A A B B A, the median of
four launches each with its spread:

| case | ε | default: median (spread) s | h1: median (spread) s | h1 ÷ default |
| --- | --- | ---: | ---: | ---: |
| w1-kd500-constant | 1e-3 | 5.15 (0.11, n = 4) | 5.39 (0.11, n = 4) | 1.047 |
| w1-kd500-constant | 1e-6 | 6.50 (0.24, n = 4) | 6.69 (0.24, n = 4) | 1.030 |
| w1-kd500-constant | 1e-8 | 10.00 (0.45, n = 4) | 10.02 (0.20, n = 4) | 1.001 |
| w2-kd500-constant | 1e-3 | 7.36 (0.06, n = 4) | 7.55 (0.19, n = 4) | 1.027 |
| w2-kd500-constant | 1e-6 | 9.52 (0.43, n = 4) | 9.70 (0.32, n = 4) | 1.019 |
| w2-kd500-constant | 1e-8 | 14.58 (0.18, n = 4) | 14.43 (0.55, n = 4) | 0.989 |
| w3:1000000-kd500-constant | 1e-3 | 6.69 (0.45, n = 4) | 6.68 (0.43, n = 4) | 0.999 |
| w3:1000000-kd500-constant | 1e-6 | 8.18 (0.76, n = 4) | 8.20 (0.67, n = 4) | 1.003 |
| w3:1000000-kd500-constant | 1e-8 | 9.77 (0.05, n = 4) | 10.12 (0.35, n = 4) | 1.036 |

`h1` costs 0–5 % of build plus one evaluation at kD = 500, 1.017× in the geometric mean over the
nine cases: the clearest W1 at 10⁻³ (1.047×, its difference 0.24 s against both sides' spread of
0.11 s), most others within the spread of four launches (up to 0.76 s, W3 at 10⁻⁶). The cost is the
evaluation's: `h1`'s build is within 2.6 % of the default's, either way, in every `cost` row.

A model check: `h1` raises each high-frequency level by one degree, (p + 2)²/(p + 1)² of its
coefficients, 1.6 % at p = 128 to 17 % at p = 11, and the RCR and diagonal M2L by up to the cube of
that; with those levels a part of the evaluation only, a few per cent is the expected order. Memory
is not measured here: the tables grow with the degree on the same levels; CH4.8b reports it (§10's
"time, memory and build cost per k").

## 8. Questions for sign-off

1. **The acceptance as read.** CH4.8 split into CH4.8a (this report: the mapping, the cause, the
   survey, the remedy chosen on FMM and pair evidence with its cost) and CH4.8b (the remedy as the
   default, CH4.8's verification under it on Kathleen and locust). *Recommended:* as read.
2. **The remedy.** *Recommended:* `h1`, the high-frequency branch one degree higher: p(x, ε) =
   max(plateau, ⌈min(x + max(0, e(ε)) x^{1/3} + 2, x + 0.82 (ln 1/ε)^{2/3} x^{1/3})⌉ + 1), CH4.3c3's
   branch and cap with one degree added after the cap (as measured, `h1` is d = 3 without the cap,
   the same degrees on every level run; it differs from "+1 after the cap" only at x = 8.8–10.7 near
   10⁻⁸, by one degree more there), the plateau unchanged; implemented in CH4.8b as a constant of
   `DegreeRule` (a margin added after the cap), verified there on CH4.8's whole set, with proposal
   §3.5 amended. A hand-chosen margin, as CH4.3c2's +2 was: the refit of CH4.3c3 had none, and W2's
   constant density shows the rule needs one. Alternatives: `h2` (two degrees: the errors halved
   again, ℓ2 ≤ 0.32 ε and worst box ≤ 33.7 ε; 0.988–1.096× the default in single launches, not
   repeated); a refit of the branch on
   W2's constant density by level scans (CH4.3c3's method; the degrees would follow from scans not
   run here, hours per case at kD = 500); the rule kept and W2's constant density named outside the
   bound (not recommended: it is inside D1 on a reference workload).
3. **The plateau.** *Recommended:* unchanged. `h1` alone keeps every measured case within ε and the
   box bound (§5, at all six ε on W2); `h1p1` changes no worst box and lowers only ∇φ's ℓ2 at 10⁻⁶ (W2: 0.37
   against 0.45 ε at kD = 500; the plateau carries about a third of ∇φ's Σ|e|² there under `h1`), while
   a higher plateau costs every low-frequency tree (W3, W4, kD ≤ 30). Alternative: β + 1 as well
   (`h1p1`), for more margin in ∇φ's ℓ2 at 10⁻⁴–10⁻⁶ (W2 at kD = 500: 0.33, 0.22, 0.37 against 0.39,
   0.33, 0.45 ε).

## 9. Tables

### Every level raised at once (`all`)

| case | ε | run | p by level | §9.1 φ, ∇φ ÷ ε | worst box φ, ∇φ ÷ ε |
| --- | --- | --- | --- | --- | --- |
| w1-kd500-constant | 1e-3 | default | 0×2, 115, 60, 32, 18, 11 | 0.22, 0.64 | 5.5, 23.4 |
| w1-kd500-constant | 1e-3 | levels [2, 3, 4, 5, 6] at +1 | 0×2, 116, 61, 33, 19, 12 | 0.09, 0.28 | 6.6, 10.1 |
| w1-kd500-constant | 1e-3 | levels [2, 3, 4, 5, 6] at +2 | 0×2, 117, 62, 34, 20, 13 | 0.05, 0.13 | 3.4, 8.9 |
| w1-kd500-constant | 1e-6 | default | 0×2, 128, 70, 40, 25, 16 | 0.10, 0.44 | 11.2, 33.6 |
| w1-kd500-constant | 1e-6 | levels [2, 3, 4, 5, 6] at +1 | 0×2, 129, 71, 41, 26, 17 | 0.04, 0.15 | 4.3, 11.7 |
| w1-kd500-constant | 1e-6 | levels [2, 3, 4, 5, 6] at +2 | 0×2, 130, 72, 42, 27, 18 | 0.02, 0.07 | 2.2, 7.1 |
| w1-kd500-constant | 1e-8 | default | 0×2, 134, 75, 44, 32×2 | 0.10, 0.34 | 13.6, 33.7 |
| w1-kd500-constant | 1e-8 | levels [2, 3, 4, 5] at +1 | 0×2, 135, 76, 45, 33, 32 | 0.04, 0.12 | 3.5, 14.2 |
| w1-kd500-constant | 1e-8 | levels [2, 3, 4, 5] at +2 | 0×2, 136, 77, 46, 34, 32 | 0.02, 0.05 | 2.4, 5.3 |
| w2-kd250-constant | 1e-3 | default | 0×2, 60, 32, 18, 11, 7×2 | 0.18, 0.98 | 13.4, 44.7 |
| w2-kd250-constant | 1e-3 | levels [2, 3, 4, 5, 6, 7] at +1 | 0×2, 61, 33, 19, 12, 8×2 | 0.07, 0.39 | 8.0, 19.5 |
| w2-kd250-constant | 1e-3 | levels [2, 3, 4, 5, 6, 7] at +2 | 0×2, 62, 34, 20, 13, 9×2 | 0.03, 0.13 | 6.6, 9.9 |
| w2-kd250-constant | 1e-6 | default | 0×2, 70, 40, 25, 16, 15×2 | 0.09, 0.53 | 38.7, 89.6 |
| w2-kd250-constant | 1e-6 | levels [2, 3, 4, 5, 6, 7] at +1 | 0×2, 71, 41, 26, 17, 16×2 | 0.04, 0.23 | 15.9, 40.8 |
| w2-kd250-constant | 1e-6 | levels [2, 3, 4, 5, 6, 7] at +2 | 0×2, 72, 42, 27, 18, 17×2 | 0.01, 0.08 | 6.4, 16.7 |
| w2-kd250-constant | 1e-8 | default | 0×2, 75, 44, 32×4 | 0.08, 0.26 | 51.2, **129.2** |
| w2-kd250-constant | 1e-8 | levels [2, 3, 4, 5, 6] at +1 | 0×2, 76, 45, 33×3, 32 | 0.03, 0.10 | 18.0, 53.8 |
| w2-kd250-constant | 1e-8 | levels [2, 3, 4, 5, 6] at +2 | 0×2, 77, 46, 34×3, 32 | 0.01, 0.04 | 6.5, 18.2 |
| w2-kd500-constant | 1e-3 | default | 0×2, 115, 60, 32, 18, 11, 7×2 | 0.18, **1.28** | 20.2, 56.6 |
| w2-kd500-constant | 1e-3 | levels [2, 3, 4, 5, 6, 7, 8] at +1 | 0×2, 116, 61, 33, 19, 12, 8×2 | 0.08, 0.47 | 20.9, 24.0 |
| w2-kd500-constant | 1e-3 | levels [2, 3, 4, 5, 6, 7, 8] at +2 | 0×2, 117, 62, 34, 20, 13, 9×2 | 0.04, 0.16 | 10.1, 21.7 |
| w2-kd500-constant | 1e-6 | default | 0×2, 128, 70, 40, 25, 16, 15×2 | 0.11, 0.99 | 46.1, **101.0** |
| w2-kd500-constant | 1e-6 | levels [2, 3, 4, 5, 6, 7, 8] at +1 | 0×2, 129, 71, 41, 26, 17, 16×2 | 0.04, 0.37 | 24.4, 47.8 |
| w2-kd500-constant | 1e-6 | levels [2, 3, 4, 5, 6, 7, 8] at +2 | 0×2, 130, 72, 42, 27, 18, 17×2 | 0.02, 0.15 | 9.1, 27.3 |
| w2-kd500-constant | 1e-8 | default | 0×2, 134, 75, 44, 32×4 | 0.08, 0.29 | 45.1, **143.4** |
| w2-kd500-constant | 1e-8 | levels [2, 3, 4, 5, 6, 7] at +1 | 0×2, 135, 76, 45, 33×3, 32 | 0.03, 0.12 | 24.0, 50.3 |
| w2-kd500-constant | 1e-8 | levels [2, 3, 4, 5, 6, 7] at +2 | 0×2, 136, 77, 46, 34×3, 32 | 0.02, 0.05 | 9.6, 21.7 |

### The mapping per launch (`map`, `survey`)

| launch | worst box ÷ ε φ, ∇φ | ℓ2 φ, ∇φ | level (share of Σ|Δ∇φ|²) | rest φ, ∇φ | largest M2L change: worst box ∇φ (run) | ℓ2 ∇φ (run) | classes in (0.1, 1] |
| --- | --- | --- | --- | --- | --- | --- | --- |
| map-w2-kd250-constant-1e-3 | 13.39, 44.72 | 3.35e-4, 9.74e-4 | 6 (Batched, 0.50) | 0.018, 0.013 | 1.000 (level 2 RCR (was diagonal 34/34)) | 1.000 (level 2 RCR (was diagonal 34/34)) | – |
| map-w2-kd250-constant-1e-6 | 38.68, 89.56 | 4.47e-7, 1.28e-6 | 5 (diagonal 34/34, 0.34) | 0.002, 0.002 | 1.000 (level 2 RCR (was diagonal 34/34)) | 1.005 (level 7 RCR (was waves)) | 6: [13, 14, 15, 16, 17] |
| map-w2-kd250-constant-1e-8 | 51.19, 129.24 | 4.39e-9, 1.12e-8 | 2 (diagonal 34/34, 0.46) | 0.001, 0.001 | 1.000 (level 2 RCR (was diagonal 34/34)) | 1.000 (level 2 RCR (was diagonal 34/34)) | 4: [7, 8] |
| map-w2-kd500-constant-1e-3 | 20.21, 56.63 | 3.35e-4, 1.07e-3 | 7 (Batched, 0.47) | 0.024, 0.017 | 1.000 (level 2 RCR (was diagonal 34/34)) | 1.000 (level 2 RCR (was diagonal 34/34)) | – |
| map-w2-kd500-constant-1e-6 | 46.09, 101.02 | 2.65e-7, 9.69e-7 | 6 (diagonal 34/34, 0.42) | 0.003, 0.002 | 1.000 (level 2 RCR (was diagonal 34/34)) | 1.003 (level 8 RCR (was waves)) | 7: [13, 14, 15, 16, 17] |
| map-w2-kd500-constant-1e-8 | 45.08, 143.40 | 2.84e-9, 8.12e-9 | 2 (diagonal 34/34, 0.78) | 0.001, 0.001 | 1.000 (level 2 RCR (was diagonal 34/34)) | 1.000 (level 2 RCR (was diagonal 34/34)) | 5: [7, 8] |
| survey-w3-kd100-constant-1e-3 | 7.02, 3.05 | 3.10e-4, 1.50e-4 | 3 (diagonal 34/34, 0.69) | 0.008, 0.013 | 1.000 (level 2 RCR (was diagonal 34/34)) | 1.000 (level 2 RCR (was diagonal 34/34)) | – |
| survey-w3-kd100-constant-1e-6 | 4.71, 4.15 | 1.44e-7, 1.15e-7 | 4 (waves, 0.42) | 0.001, 0.001 | 1.013 (level 5 RCR (was waves)) | 1.004 (level 4 RCR (was waves)) | – |
| survey-w3-kd100-constant-1e-8 | 4.54, 3.61 | 7.58e-10, 6.37e-10 | 3 (waves, 0.50) | 0.000, 0.001 | 0.996 (level 3 RCR (was waves)) | 0.999 (level 5 RCR (was waves)) | – |
| survey-w3:1000000-kd200-constant-1e-3 | 6.56, 10.57 | 3.31e-4, 4.36e-4 | 4 (waves, 0.78) | 0.008, 0.010 | 1.000 (level 2 RCR (was diagonal 34/34)) | 1.000 (level 2 RCR (was diagonal 34/34)) | – |
| survey-w3:1000000-kd200-constant-1e-6 | 6.61, 12.33 | 1.74e-7, 2.72e-7 | 4 (waves, 0.58) | 0.001, 0.002 | 0.999 (level 4 RCR (was waves)) | 1.000 (level 2 RCR (was diagonal 34/34)) | 5: [0, 1] |
| survey-w3:1000000-kd200-constant-1e-8 | 4.86, 11.26 | 1.43e-9, 2.49e-9 | 4 (waves, 0.63) | 0.000, 0.000 | 1.010 (level 4 RCR (was waves)) | 1.000 (level 2 RCR (was diagonal 34/34)) | 2: [0, 1] |
| survey-w3:1000000-kd500-constant-1e-3 | 9.41, 17.26 | 3.61e-4, 5.03e-4 | 5 (diagonal 34/34, 0.55) | 0.014, 0.016 | 1.000 (level 2 RCR (was diagonal 34/34)) | 1.000 (level 2 RCR (was diagonal 34/34)) | – |
| survey-w3:1000000-kd500-constant-1e-6 | 10.03, 12.74 | 1.62e-7, 2.44e-7 | 4 (diagonal 34/34, 0.57) | 0.001, 0.002 | 1.000 (level 2 RCR (was diagonal 34/34)) | 1.000 (level 2 RCR (was diagonal 34/34)) | – |
| survey-w3:1000000-kd500-constant-1e-8 | 10.11, 17.30 | 1.93e-9, 3.32e-9 | 4 (waves, 0.58) | 0.000, 0.001 | 1.000 (level 2 RCR (was diagonal 34/34)) | 1.000 (level 2 RCR (was diagonal 34/34)) | – |

### Levels (`map`)

| launch | level | k a | p | M2L | share φ | share ∇φ | classes (0.1, 1] | (0.01, 0.1] | failing |
| --- | ---: | ---: | ---: | --- | ---: | ---: | --- | ---: | ---: |
| map-w2-kd250-constant-1e-3 | 2 | 54.13 | 60 | diagonal 34/34 | 0.655 | 0.246 | [] | 0 | 0 |
| map-w2-kd250-constant-1e-3 | 3 | 27.06 | 32 | diagonal 34/34 | 0.548 | 0.228 | [] | 0 | 0 |
| map-w2-kd250-constant-1e-3 | 4 | 13.53 | 18 | diagonal 34/34 | 0.236 | 0.156 | [] | 0 | 0 |
| map-w2-kd250-constant-1e-3 | 5 | 6.77 | 11 | diagonal 34/34 | 0.122 | 0.172 | [] | 0 | 0 |
| map-w2-kd250-constant-1e-3 | 6 | 3.38 | 7 | Batched | 0.171 | 0.497 |  |  |  |
| map-w2-kd250-constant-1e-3 | 7 | 1.69 | 7 | Batched | 0.001 | 0.016 |  |  |  |
| map-w2-kd250-constant-1e-6 | 2 | 54.13 | 70 | diagonal 34/34 | 0.329 | 0.193 | [] | 0 | 0 |
| map-w2-kd250-constant-1e-6 | 3 | 27.06 | 40 | diagonal 34/34 | 0.109 | 0.130 | [] | 0 | 0 |
| map-w2-kd250-constant-1e-6 | 4 | 13.53 | 25 | diagonal 34/34 | 0.010 | 0.022 | [] | 0 | 0 |
| map-w2-kd250-constant-1e-6 | 5 | 6.77 | 16 | diagonal 34/34 | 0.154 | 0.345 | [] | 2 | 0 |
| map-w2-kd250-constant-1e-6 | 6 | 3.38 | 15 | diagonal 21/34 | 0.005 | 0.051 | [13, 14, 15, 16, 17] | 4 | 13 |
| map-w2-kd250-constant-1e-6 | 7 | 1.69 | 15 | waves | 0.001 | 0.019 |  |  |  |
| map-w2-kd250-constant-1e-8 | 2 | 54.13 | 75 | diagonal 34/34 | 0.415 | 0.459 | [] | 2 | 0 |
| map-w2-kd250-constant-1e-8 | 3 | 27.06 | 44 | diagonal 32/34 | 0.314 | 0.422 | [] | 0 | 2 |
| map-w2-kd250-constant-1e-8 | 4 | 13.53 | 32 | diagonal 27/34 | 0.010 | 0.034 | [7, 8] | 4 | 7 |
| map-w2-kd250-constant-1e-8 | 5 | 6.77 | 32 | waves | 0.000 | 0.000 |  |  |  |
| map-w2-kd250-constant-1e-8 | 6 | 3.38 | 32 | waves | 0.000 | 0.000 |  |  |  |
| map-w2-kd250-constant-1e-8 | 7 | 1.69 | 32 | – | 0.000 | 0.000 |  |  |  |
| map-w2-kd500-constant-1e-3 | 2 | 108.25 | 115 | diagonal 34/34 | 0.337 | 0.094 | [] | 0 | 0 |
| map-w2-kd500-constant-1e-3 | 3 | 54.13 | 60 | diagonal 34/34 | 0.252 | 0.065 | [] | 0 | 0 |
| map-w2-kd500-constant-1e-3 | 4 | 27.06 | 32 | diagonal 34/34 | 0.132 | 0.079 | [] | 0 | 0 |
| map-w2-kd500-constant-1e-3 | 5 | 13.53 | 18 | diagonal 34/34 | 0.195 | 0.141 | [] | 0 | 0 |
| map-w2-kd500-constant-1e-3 | 6 | 6.77 | 11 | diagonal 34/34 | 0.098 | 0.154 | [] | 0 | 0 |
| map-w2-kd500-constant-1e-3 | 7 | 3.38 | 7 | Batched | 0.151 | 0.473 |  |  |  |
| map-w2-kd500-constant-1e-3 | 8 | 1.69 | 7 | Batched | 0.001 | 0.015 |  |  |  |
| map-w2-kd500-constant-1e-6 | 2 | 108.25 | 128 | diagonal 34/34 | 0.536 | 0.160 | [] | 0 | 0 |
| map-w2-kd500-constant-1e-6 | 3 | 54.13 | 70 | diagonal 34/34 | 0.223 | 0.100 | [] | 0 | 0 |
| map-w2-kd500-constant-1e-6 | 4 | 27.06 | 40 | diagonal 34/34 | 0.086 | 0.066 | [] | 0 | 0 |
| map-w2-kd500-constant-1e-6 | 5 | 13.53 | 25 | diagonal 34/34 | 0.013 | 0.014 | [] | 0 | 0 |
| map-w2-kd500-constant-1e-6 | 6 | 6.77 | 16 | diagonal 34/34 | 0.177 | 0.423 | [] | 2 | 0 |
| map-w2-kd500-constant-1e-6 | 7 | 3.38 | 15 | diagonal 21/34 | 0.011 | 0.086 | [13, 14, 15, 16, 17] | 4 | 13 |
| map-w2-kd500-constant-1e-6 | 8 | 1.69 | 15 | waves | 0.001 | 0.025 |  |  |  |
| map-w2-kd500-constant-1e-8 | 2 | 108.25 | 134 | diagonal 34/34 | 1.045 | 0.777 | [] | 0 | 0 |
| map-w2-kd500-constant-1e-8 | 3 | 54.13 | 75 | diagonal 34/34 | 0.261 | 0.348 | [] | 2 | 0 |
| map-w2-kd500-constant-1e-8 | 4 | 27.06 | 44 | diagonal 32/34 | 0.152 | 0.261 | [] | 0 | 2 |
| map-w2-kd500-constant-1e-8 | 5 | 13.53 | 32 | diagonal 27/34 | 0.018 | 0.041 | [7, 8] | 4 | 7 |
| map-w2-kd500-constant-1e-8 | 6 | 6.77 | 32 | waves | 0.000 | 0.000 |  |  |  |
| map-w2-kd500-constant-1e-8 | 7 | 3.38 | 32 | waves | 0.000 | 0.000 |  |  |  |
| map-w2-kd500-constant-1e-8 | 8 | 1.69 | 32 | – | 0.000 | 0.000 |  |  |  |
