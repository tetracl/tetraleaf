# Helmholtz CH4.5: host performance against FMM3D at matched accuracy

Task: `docs/helmholtz/tasks/CH4.5.md`. Roadmap: H4, CH4.5 (priority list item 10, the last of H4's
tasks). Proposal §10 CH4.5, §9.2 (what is compared, setup apart; FMM3D's eps chosen by its error),
§9.1, §2. Method as C6.8 (`tetraleaf-validate/results/phase6-fmm3d.md`, `tools/fmm3d/README.md`), on
CH4.4a's sweep (`helmholtz-ch4-4a.md`), with CH4.0e's splits and one-sided rule
(`helmholtz-ch4-0e.md`). Alsos: review 2026-11-19 A6, B6, B8; CH4.3e §4; CH4.4a questions 3–4.

## 0. Summary

- **Ours under the library's defaults against FMM3D's `hfmm3d` at matched accuracy, one node**: our
  build plus one evaluation over FMM3D's whole call, FMM3D at the eps whose error is closest to
  ours (§1), its faster variant; geometric means over CH4.4a's 26 cases, and over the 14 with
  N ≥ 10⁵ and the 8 with N ≥ 10⁶:

  | | Kathleen 10⁻³ | 10⁻⁶ | 10⁻⁸ | locust 10⁻³ | 10⁻⁶ | 10⁻⁸ |
  | --- | --- | --- | --- | --- | --- | --- |
  | best split, all 26 | 0.41 | 0.31 | 0.34 | 0.25 | 0.26 | 0.26 |
  | best split, N ≥ 10⁵ | 0.20 | 0.19 | 0.23 | 0.17 | 0.16 | 0.19 |
  | best split, N ≥ 10⁶ | 0.13 | 0.14 | 0.18 | 0.13 | 0.14 | 0.16 |
  | 1 × N, all 26 | 0.52 | 0.39 | 0.40 | 0.43 | 0.41 | 0.32 |
  | 1 × N, N ≥ 10⁵ | 0.26 | 0.24 | 0.27 | 0.34 | 0.27 | 0.23 |
  | evaluation ÷ (call − tree), N ≥ 10⁵ | 0.10 | 0.12 | 0.18 | 0.09 | 0.10 | 0.15 |

  The one-sided match (FMM3D's fastest eps at least as accurate as ours, CH4.0e's) gives the same
  picture (§3.1). For N ≥ 10⁵ ours is faster on every case at every ε on both machines at the best
  split (0.06–0.97× FMM3D's call on Kathleen, 0.08–0.43× on locust; W1 and the plane wave at
  kD = 1000 0.06–0.11× and 0.08–0.10×), and at 1 × N but for one row each, W3 at kD = 30 at 10⁻³
  (1.16× on Kathleen; 1.09× with the constant density on locust).
- **Where FMM3D is faster** (§3.4): ten rows of N ≤ 1.1·10⁴ on Kathleen at the best split, seven
  of them at 10⁻³ (1.1–4.4×; at 10⁻⁶ and 10⁻⁸ W2 and W3 at kD = 10, 1.3–1.4×), calls of 5–104 ms. Our
  evaluation is faster there too (0.15–0.74× FMM3D's call less its tree); our build is not. On six
  of them (W1 and W2 at kD = 30 and W3 at kD = 10, both densities, 10⁻³) Kathleen's costs choose a
  shallow tree with Batched at p = 11 at every split, whose matrices take 15–101 ms, where locust's
  choose one level deeper (the diagonal form, then Batched at p = 7, 6 ms); the other four (W2 at
  kD = 10, W3 at kD = 10 at 10⁻⁸) are a few milliseconds of fixed overhead. Proposed follow-up CH4.7
  (question 4). On locust no case loses at the best split.
- **Setup apart** (§4): FMM3D's setup (its tree, timed once apart) is 0.2–23 % of its call for
  N ≥ 10⁵; our evaluation alone over FMM3D's call less its tree is 0.03–0.52× there.
- **One rank per node** (§5, the H4 gate's figures; CH4.0k): at 1 × N the steps every rank repeats
  (octree, routes, plan) are 32–88 % of our build on the cases of N ≥ 10⁵ on Kathleen and 58–91 % on
  locust, against 24–66 % and 27–83 % at the best split; the build at 1 × N is 1.1–7.7× the best
  split's, the evaluation 0.80–1.32×. W1 at kD = 1000 builds in 12.6–14.3 s at 1 × 40 against
  8.8–10.8 s at 2 × 20 (Kathleen), 7.6–14.0 s at 1 × 72 against 4.4–5.5 s at 4 × 18 (locust). With
  the build in, 1 × N stays below FMM3D's call on every case of N ≥ 10⁵ but the two rows above.
- **Memory** (§6): FMM3D's peak over our FMM's resident growth at the best split (the node's sum),
  N ≥ 10⁵: 2.5× (0.18–24×) on Kathleen, 2.2× (0.38–13×) on locust, geometric means; the measures
  differ (FMM3D's whole process, its arrays included; our FMM's growth). FMM3D `default` needs
  442–447 GB at kD = 1000 at 3·10⁻³ and 10⁻³ (locust; beyond Kathleen's node, CH4.4a), and fails on
  Kathleen on W2 at kD = 500 at 3·10⁻³ and 10⁻³ (exit 137, the node's 180 GB).
- **c_leaf** (§7, CH4.3e's Also): on W3 at N = 10⁶, kD = 500, c_leaf 4 takes 0.29–0.74× c_leaf 2's
  build plus one evaluation on Kathleen at every split and at 40 × 1, and 0.24–0.89× on locust at the
  four splits; at 72 × 1 0.80–1.11× (1.11 at 10⁻⁶).
- **FMM3D's OpenMP binding** (§8, B8): at its defaults (libgomp: no `OMP_PROC_BIND`, threads
  unbound), stated in every log; bound (`OMP_PLACES=cores OMP_PROC_BIND=close`) it is 1.01–1.14×
  slower on Kathleen and 0.99–1.04× on locust on five large cases: the comparison stands unbound.
- **Checksums** (§9): every case's checksum is the same at every split and in both launches on
  both machines (624 launches each), and so are c_leaf's on Kathleen including 40 × 1; at 72 × 1 on
  locust the scratch budget per rank (114 MiB) sends W3's levels 2–3 to RCR instead of the diagonal
  form (CH4.3b2c §8), a different plan with its own checksum, the same in both launches.
- **FMM3D's errors** are the same on both machines (508 points of locust's timed runs against
  Kathleen's grid: within 1.06×, at round-off), so the matches are the grid's on both.

## 1. Method

**The sweep** is CH4.4a's: W1 at kD = 10, 30, 100, 250, 500, 1000 and W2 at 10–500 with §9.1's random
charges; W3 resolved at kD = 10, 30, 100 and at N = 10⁶ for kD = 200, 500; the plane wave's trace on W1
at kD = 500, 1000; the constant density on W1 and W2 at kD = 10, 30, 100 and W3 at 10, 30. 26 cases at
ε = 10⁻³, 10⁻⁶, 10⁻⁸, potentials and gradients, errors against the f64 direct sum at §9.1's 1,000
sampled targets, ours and FMM3D's on byte-identical problem files.

**Ours**: the library's defaults (`helmholtz_fmm3d compare --settings defaults`), one ε per launch,
one build and two evaluations (`--evaluations 2`): the build is the setup, the second evaluation the
evaluation (the first printed beside it), every time the slowest rank's. Splits: 1 × 40, 2 × 20,
4 × 10, 8 × 5 on Kathleen (two Cascade Lake sockets, 40 cores; even rank counts bound per socket),
1 × 72, 2 × 36, 4 × 18, 8 × 9 on locust (72 Grace cores), each in launch a (that order) and launch b
(reversed). **Best split**: the smallest build plus one evaluation over the splits and both launches
(CH4.0e's rule); **1 × N**: the one-rank split, the faster of its two launches.

**FMM3D** (v2.1.0 as built by `tools/fmm3d/build.sh`; `default` and `fastker`): one launch per eps
under `tools/fmm3d/maxrss.py`, OpenMP on every core at its defaults. The grid (Kathleen, REPEAT 1):
14 eps, 3·10⁻³ to 10⁻⁹, two per decade; both variants but `fastker` alone at kD = 1000 (CH4.4a's
question 3) and `default` alone on the plane wave and the constant density (CH4.4a's runs).
Timed (stage `time`): each case's matched eps (both rules, on either machine, and ε itself),
REPEAT 2 (the best of two calls) below N = 10⁶ and 1 from there, launches a and b in opposite orders;
FMM3D's time the smaller of the two launches. On locust also `default` at kD = 1000 (472 GB node,
455 GB available at every launch).

**The matches** (`tools/helmholtz/ch4-5-tables.py`), per case, ε and variant over the grid:
- **closest** (the acceptance's): the eps minimising d = max over φ and ∇φ of |log₁₀(FMM3D's error /
  ours)|, the faster eps on a tie. C6.8 took φ alone; with two fields the worse of the two. Of the
  variants, the faster call among those whose d is within one decade of the closer one's: on the
  smallest trees `fastker` sums directly at every eps (10⁻¹⁵), 10¹⁰ away from our 4·10⁻⁵ where
  `default`'s match is within 2×.
- **one-sided** (CH4.0e's): the fastest eps whose φ and ∇φ errors are both at most ours; where none
  is (ours at round-off), the most accurate, marked `>`. It never gives FMM3D a looser result.
FMM3D's error moves in steps (its order tables): on W1 at kD = 1000 1.7·10⁻⁴ for 3·10⁻³–10⁻³,
5.7·10⁻⁸ for 3·10⁻⁴–10⁻⁶, 4.2·10⁻¹¹ below. Where ours falls between two steps the closest match is
up to 36× away (eleven rows, nine at 10⁻⁸: FMM3D's step looser than ours, the one-sided match the
step below it); both are tabled.

**Setup apart** (C6.8's measure): FMM3D keeps nothing between calls, so its call is timed whole and
its setup (`hndiv`, `pts_tree_mem`, `pts_tree_build`, `pts_tree_sort` as `hfmm3d` calls them, the
driver's `tree_s`) timed once apart per eps. Ratios: ours' build plus one evaluation over FMM3D's
call (one solve); ours' evaluation over FMM3D's call less its tree (repeated solves), left out where
the tree is not below the call (one case, a call of milliseconds).

**Memory**: ours `Row::memory`, the FMM's resident growth over the build and evaluations, summed
over the node's ranks (the largest rank's beside it), one ε per launch; FMM3D's the launch's peak
RSS (`maxrss.py`), its input and output arrays included.

**Binding** (stage `bind`): FMM3D bound (`OMP_PLACES=cores OMP_PROC_BIND=close`) beside its defaults
on W1 at kD = 250, 500, 1000, W3 at kD = 100 and at N = 10⁶, kD = 500, at the closest eps of 10⁻⁶,
each variant, both bindings in launches a and b.

**c_leaf** (stage `cleaf`): W3 at N = 10⁶, kD = 500 at c_leaf 2 and 4 (`compare --c-leaf`), the four
splits and N × 1 (40 × 1, 72 × 1), launches a and b.

## 2. The matches

Per case, ε and variant, ours' errors and FMM3D's at each rule's match, with Kathleen's errors of
ours (locust's differ in the last digits; its matches, in the archive's tables, differ on two rows
of small cases, §11).

<details><summary>Kathleen's matches (one row per case, ε and variant)</summary>

| Case | N | ε | ours φ, ∇φ | variant | closest eps: φ, ∇φ (factor) | one-sided eps: φ, ∇φ |
| --- | ---: | --- | --- | --- | --- | --- |
| w1-kd10 | 796 | 1e-3 | 1.44e-15, 1.74e-15 | default | 1e-4: 1.22e-15, 1.30e-15 (1.3) | 1e-4: 1.22e-15, 1.30e-15 |
| w1-kd10 | 796 | 1e-3 | 1.44e-15, 1.74e-15 | fastker | 3e-4: 1.28e-15, 1.39e-15 (1.3) | 1e-3: 9.03e-16, 1.06e-15 |
| w1-kd10 | 796 | 1e-6 | 1.44e-15, 1.74e-15 | default | 1e-4: 1.22e-15, 1.30e-15 (1.3) | 1e-4: 1.22e-15, 1.30e-15 |
| w1-kd10 | 796 | 1e-6 | 1.44e-15, 1.74e-15 | fastker | 3e-4: 1.28e-15, 1.39e-15 (1.3) | 1e-3: 9.03e-16, 1.06e-15 |
| w1-kd10 | 796 | 1e-8 | 1.44e-15, 1.74e-15 | default | 1e-4: 1.22e-15, 1.30e-15 (1.3) | 1e-4: 1.22e-15, 1.30e-15 |
| w1-kd10 | 796 | 1e-8 | 1.44e-15, 1.74e-15 | fastker | 3e-4: 1.28e-15, 1.39e-15 (1.3) | 1e-3: 9.03e-16, 1.06e-15 |
| w1-kd30 | 7162 | 1e-3 | 5.07e-05, 8.05e-05 | default | 1e-3: 7.50e-05, 7.42e-05 (1.5) | 3e-4: 1.45e-08, 2.76e-08 |
| w1-kd30 | 7162 | 1e-3 | 5.07e-05, 8.05e-05 | fastker | 3e-3: 7.48e-05, 7.37e-05 (1.5) | 1e-4: 2.29e-15, 3.06e-15 |
| w1-kd30 | 7162 | 1e-6 | 1.08e-07, 2.50e-07 | default | 3e-4: 1.45e-08, 2.76e-08 (9.1) | 3e-4: 1.45e-08, 2.76e-08 |
| w1-kd30 | 7162 | 1e-6 | 1.08e-07, 2.50e-07 | fastker | 3e-3: 7.48e-05, 7.37e-05 (692.7) | 1e-4: 2.29e-15, 3.06e-15 |
| w1-kd30 | 7162 | 1e-8 | 8.83e-13, 2.13e-12 | default | 1e-8: 5.67e-12, 7.75e-12 (6.4) | >1e-8: 5.67e-12, 7.75e-12 |
| w1-kd30 | 7162 | 1e-8 | 8.83e-13, 2.13e-12 | fastker | 1e-4: 2.29e-15, 3.06e-15 (695.2) | 1e-4: 2.29e-15, 3.06e-15 |
| w1-kd100 | 79578 | 1e-3 | 1.54e-04, 2.55e-04 | default | 3e-3: 1.23e-03, 1.20e-03 (8.0) | 3e-4: 3.25e-07, 4.41e-07 |
| w1-kd100 | 79578 | 1e-3 | 1.54e-04, 2.55e-04 | fastker | 3e-3: 1.23e-03, 1.20e-03 (8.0) | 3e-4: 3.23e-07, 4.37e-07 |
| w1-kd100 | 79578 | 1e-6 | 1.39e-07, 2.99e-07 | default | 1e-5: 3.25e-07, 4.40e-07 (2.3) | 3e-7: 1.33e-10, 2.09e-10 |
| w1-kd100 | 79578 | 1e-6 | 1.39e-07, 2.99e-07 | fastker | 1e-5: 3.23e-07, 4.36e-07 (2.3) | 3e-7: 1.33e-10, 2.07e-10 |
| w1-kd100 | 79578 | 1e-8 | 6.28e-10, 8.74e-10 | default | 3e-7: 1.33e-10, 2.09e-10 (4.7) | 3e-7: 1.33e-10, 2.09e-10 |
| w1-kd100 | 79578 | 1e-8 | 6.28e-10, 8.74e-10 | fastker | 3e-7: 1.33e-10, 2.07e-10 (4.7) | 3e-7: 1.33e-10, 2.07e-10 |
| w1-kd250 | 497360 | 1e-3 | 2.28e-04, 3.12e-04 | default | 3e-3: 2.09e-04, 1.96e-04 (1.6) | 3e-3: 2.09e-04, 1.96e-04 |
| w1-kd250 | 497360 | 1e-3 | 2.28e-04, 3.12e-04 | fastker | 3e-3: 2.09e-04, 1.95e-04 (1.6) | 3e-3: 2.09e-04, 1.95e-04 |
| w1-kd250 | 497360 | 1e-6 | 1.50e-07, 2.32e-07 | default | 1e-5: 6.82e-08, 7.81e-08 (3.0) | 3e-4: 6.84e-08, 7.77e-08 |
| w1-kd250 | 497360 | 1e-6 | 1.50e-07, 2.32e-07 | fastker | 1e-5: 6.60e-08, 7.34e-08 (3.2) | 3e-4: 6.62e-08, 7.29e-08 |
| w1-kd250 | 497360 | 1e-8 | 1.31e-09, 1.52e-09 | default | 3e-7: 5.65e-11, 6.30e-11 (24.1) | 3e-7: 5.65e-11, 6.30e-11 |
| w1-kd250 | 497360 | 1e-8 | 1.31e-09, 1.52e-09 | fastker | 3e-7: 5.65e-11, 6.28e-11 (24.2) | 3e-7: 5.65e-11, 6.28e-11 |
| w1-kd500 | 1989437 | 1e-3 | 2.85e-04, 3.97e-04 | default | 3e-3: 1.87e-04, 1.75e-04 (2.3) | 3e-3: 1.87e-04, 1.75e-04 |
| w1-kd500 | 1989437 | 1e-3 | 2.85e-04, 3.97e-04 | fastker | 3e-3: 1.86e-04, 1.75e-04 (2.3) | 3e-3: 1.86e-04, 1.75e-04 |
| w1-kd500 | 1989437 | 1e-6 | 1.51e-07, 2.96e-07 | default | 3e-4: 6.22e-08, 7.22e-08 (4.1) | 3e-4: 6.22e-08, 7.22e-08 |
| w1-kd500 | 1989437 | 1e-6 | 1.51e-07, 2.96e-07 | fastker | 3e-4: 6.10e-08, 6.97e-08 (4.2) | 3e-4: 6.10e-08, 6.97e-08 |
| w1-kd500 | 1989437 | 1e-8 | 1.16e-09, 1.48e-09 | default | 1e-8: 5.19e-11, 6.34e-11 (23.4) | 3e-7: 5.16e-11, 6.32e-11 |
| w1-kd500 | 1989437 | 1e-8 | 1.16e-09, 1.48e-09 | fastker | 1e-8: 5.17e-11, 6.31e-11 (23.5) | 3e-7: 5.15e-11, 6.29e-11 |
| w1-kd1000 | 7957748 | 1e-3 | 2.68e-04, 3.93e-04 | fastker | 3e-3: 1.76e-04, 1.67e-04 (2.3) | 3e-3: 1.76e-04, 1.67e-04 |
| w1-kd1000 | 7957748 | 1e-6 | 1.38e-07, 2.81e-07 | fastker | 3e-4: 5.72e-08, 6.32e-08 (4.4) | 1e-4: 5.70e-08, 6.31e-08 |
| w1-kd1000 | 7957748 | 1e-8 | 2.28e-09, 1.91e-09 | fastker | 3e-5: 5.69e-08, 6.31e-08 (33.0) | 3e-7: 4.21e-11, 5.04e-11 |
| w2-kd10 | 1276 | 1e-3 | 2.42e-15, 7.13e-15 | default | 1e-7: 1.72e-15, 4.81e-15 (1.5) | 1e-7: 1.72e-15, 4.81e-15 |
| w2-kd10 | 1276 | 1e-3 | 2.42e-15, 7.13e-15 | fastker | 1e-3: 1.45e-15, 4.73e-15 (1.7) | 1e-3: 1.45e-15, 4.73e-15 |
| w2-kd10 | 1276 | 1e-6 | 2.42e-15, 7.13e-15 | default | 1e-7: 1.72e-15, 4.81e-15 (1.5) | 1e-7: 1.72e-15, 4.81e-15 |
| w2-kd10 | 1276 | 1e-6 | 2.42e-15, 7.13e-15 | fastker | 1e-3: 1.45e-15, 4.73e-15 (1.7) | 1e-3: 1.45e-15, 4.73e-15 |
| w2-kd10 | 1276 | 1e-8 | 2.42e-15, 7.13e-15 | default | 1e-7: 1.72e-15, 4.81e-15 (1.5) | 1e-7: 1.72e-15, 4.81e-15 |
| w2-kd10 | 1276 | 1e-8 | 2.42e-15, 7.13e-15 | fastker | 1e-3: 1.45e-15, 4.73e-15 (1.7) | 1e-3: 1.45e-15, 4.73e-15 |
| w2-kd30 | 11480 | 1e-3 | 9.09e-05, 1.23e-04 | default | 3e-3: 5.07e-05, 9.18e-05 (1.8) | 3e-3: 5.07e-05, 9.18e-05 |
| w2-kd30 | 11480 | 1e-3 | 9.09e-05, 1.23e-04 | fastker | 1e-3: 3.74e-05, 2.34e-05 (5.3) | 3e-3: 3.75e-05, 2.34e-05 |
| w2-kd30 | 11480 | 1e-6 | 8.97e-08, 1.57e-07 | default | 3e-4: 1.58e-08, 5.91e-08 (5.7) | 3e-4: 1.58e-08, 5.91e-08 |
| w2-kd30 | 11480 | 1e-6 | 8.97e-08, 1.57e-07 | fastker | 3e-4: 8.70e-09, 1.18e-08 (13.3) | 3e-4: 8.70e-09, 1.18e-08 |
| w2-kd30 | 11480 | 1e-8 | 2.39e-12, 1.90e-12 | default | 1e-9: 2.89e-12, 3.42e-12 (1.8) | >1e-9: 2.89e-12, 3.42e-12 |
| w2-kd30 | 11480 | 1e-8 | 2.39e-12, 1.90e-12 | fastker | 1e-9: 2.77e-12, 3.33e-12 (1.8) | >1e-9: 2.77e-12, 3.33e-12 |
| w2-kd100 | 127554 | 1e-3 | 1.46e-04, 1.27e-04 | default | 3e-3: 7.97e-04, 4.32e-04 (5.5) | 3e-4: 1.89e-07, 1.55e-07 |
| w2-kd100 | 127554 | 1e-3 | 1.46e-04, 1.27e-04 | fastker | 3e-3: 7.94e-04, 4.18e-04 (5.4) | 3e-4: 1.89e-07, 1.51e-07 |
| w2-kd100 | 127554 | 1e-6 | 1.24e-07, 2.32e-07 | default | 1e-4: 1.89e-07, 1.54e-07 (1.5) | 3e-7: 9.74e-11, 9.20e-11 |
| w2-kd100 | 127554 | 1e-6 | 1.24e-07, 2.32e-07 | fastker | 3e-4: 1.89e-07, 1.51e-07 (1.5) | 3e-7: 9.71e-11, 9.17e-11 |
| w2-kd100 | 127554 | 1e-8 | 4.66e-10, 4.49e-10 | default | 3e-7: 9.74e-11, 9.20e-11 (4.9) | 3e-7: 9.74e-11, 9.20e-11 |
| w2-kd100 | 127554 | 1e-8 | 4.66e-10, 4.49e-10 | fastker | 3e-7: 9.71e-11, 9.17e-11 (4.9) | 3e-7: 9.71e-11, 9.17e-11 |
| w2-kd250 | 797211 | 1e-3 | 1.96e-04, 1.97e-04 | default | 3e-3: 1.44e-04, 1.20e-04 (1.6) | 3e-4: 4.46e-08, 3.72e-08 |
| w2-kd250 | 797211 | 1e-3 | 1.96e-04, 1.97e-04 | fastker | 3e-3: 1.39e-04, 7.10e-05 (2.8) | 3e-3: 1.39e-04, 7.10e-05 |
| w2-kd250 | 797211 | 1e-6 | 1.18e-07, 1.51e-07 | default | 3e-4: 4.46e-08, 3.72e-08 (4.1) | 3e-4: 4.46e-08, 3.72e-08 |
| w2-kd250 | 797211 | 1e-6 | 1.18e-07, 1.51e-07 | fastker | 3e-4: 4.43e-08, 2.50e-08 (6.0) | 3e-4: 4.43e-08, 2.50e-08 |
| w2-kd250 | 797211 | 1e-8 | 7.24e-10, 4.80e-10 | default | 3e-7: 3.09e-11, 2.55e-11 (23.4) | 3e-7: 3.09e-11, 2.55e-11 |
| w2-kd250 | 797211 | 1e-8 | 7.24e-10, 4.80e-10 | fastker | 1e-7: 3.07e-11, 1.95e-11 (24.6) | 3e-7: 3.08e-11, 1.95e-11 |
| w2-kd500 | 3188840 | 1e-3 | 2.43e-04, 2.64e-04 | default | 3e-4: 5.15e-08, 3.73e-08 (7070.2) | 3e-4: 5.15e-08, 3.73e-08 |
| w2-kd500 | 3188840 | 1e-3 | 2.43e-04, 2.64e-04 | fastker | 3e-3: 1.50e-04, 8.45e-05 (3.1) | 3e-3: 1.50e-04, 8.45e-05 |
| w2-kd500 | 3188840 | 1e-6 | 1.38e-07, 2.35e-07 | default | 3e-4: 5.15e-08, 3.73e-08 (6.3) | 3e-4: 5.15e-08, 3.73e-08 |
| w2-kd500 | 3188840 | 1e-6 | 1.38e-07, 2.35e-07 | fastker | 3e-4: 5.09e-08, 3.00e-08 (7.8) | 3e-4: 5.09e-08, 3.00e-08 |
| w2-kd500 | 3188840 | 1e-8 | 8.34e-10, 4.83e-10 | default | 3e-7: 4.55e-11, 3.22e-11 (18.3) | 3e-7: 4.55e-11, 3.22e-11 |
| w2-kd500 | 3188840 | 1e-8 | 8.34e-10, 4.83e-10 | fastker | 3e-7: 4.54e-11, 2.89e-11 (18.4) | 3e-7: 4.54e-11, 2.89e-11 |
| w3-kd10 | 4032 | 1e-3 | 3.88e-05, 6.58e-05 | default | 3e-3: 1.72e-05, 2.28e-05 (2.9) | 3e-3: 1.72e-05, 2.28e-05 |
| w3-kd10 | 4032 | 1e-3 | 3.88e-05, 6.58e-05 | fastker | 3e-3: 1.72e-05, 2.28e-05 (2.9) | 3e-3: 1.72e-05, 2.28e-05 |
| w3-kd10 | 4032 | 1e-6 | 4.01e-08, 1.24e-07 | default | 3e-4: 4.18e-08, 1.38e-07 (1.1) | 1e-4: 2.32e-08, 8.34e-08 |
| w3-kd10 | 4032 | 1e-6 | 4.01e-08, 1.24e-07 | fastker | 3e-3: 1.72e-05, 2.28e-05 (429.9) | 1e-4: 1.56e-15, 1.29e-14 |
| w3-kd10 | 4032 | 1e-8 | 2.52e-15, 1.96e-14 | default | 1e-8: 7.52e-12, 1.85e-11 (2983.7) | >1e-8: 7.52e-12, 1.85e-11 |
| w3-kd10 | 4032 | 1e-8 | 2.52e-15, 1.96e-14 | fastker | 1e-4: 1.56e-15, 1.29e-14 (1.6) | 1e-4: 1.56e-15, 1.29e-14 |
| w3-kd30 | 108849 | 1e-3 | 1.60e-04, 3.26e-04 | default | 1e-3: 1.01e-04, 8.82e-05 (3.7) | 3e-3: 1.01e-04, 8.81e-05 |
| w3-kd30 | 108849 | 1e-3 | 1.60e-04, 3.26e-04 | fastker | 1e-3: 1.00e-04, 8.17e-05 (4.0) | 3e-3: 1.00e-04, 8.17e-05 |
| w3-kd30 | 108849 | 1e-6 | 1.44e-07, 3.60e-07 | default | 3e-4: 2.14e-08, 1.21e-07 (6.7) | 3e-4: 2.14e-08, 1.21e-07 |
| w3-kd30 | 108849 | 1e-6 | 1.44e-07, 3.60e-07 | fastker | 3e-4: 1.66e-08, 3.70e-08 (9.7) | 3e-4: 1.66e-08, 3.70e-08 |
| w3-kd30 | 108849 | 1e-8 | 5.04e-12, 6.99e-12 | default | 1e-9: 5.79e-12, 1.45e-11 (2.1) | >1e-9: 5.79e-12, 1.45e-11 |
| w3-kd30 | 108849 | 1e-8 | 5.04e-12, 6.99e-12 | fastker | 1e-9: 5.79e-12, 1.45e-11 (2.1) | >1e-9: 5.79e-12, 1.45e-11 |
| w3-kd100 | 4031442 | 1e-3 | 3.19e-04, 3.32e-04 | default | 3e-3: 2.03e-03, 1.31e-03 (6.4) | 3e-4: 5.45e-07, 5.18e-07 |
| w3-kd100 | 4031442 | 1e-3 | 3.19e-04, 3.32e-04 | fastker | 3e-3: 2.03e-03, 1.30e-03 (6.4) | 3e-4: 5.45e-07, 5.17e-07 |
| w3-kd100 | 4031442 | 1e-6 | 2.60e-07, 3.63e-07 | default | 3e-4: 5.45e-07, 5.18e-07 (2.1) | 3e-7: 1.89e-10, 2.17e-10 |
| w3-kd100 | 4031442 | 1e-6 | 2.60e-07, 3.63e-07 | fastker | 3e-4: 5.45e-07, 5.17e-07 (2.1) | 3e-7: 1.88e-10, 2.16e-10 |
| w3-kd100 | 4031442 | 1e-8 | 1.18e-09, 1.82e-09 | default | 3e-7: 1.89e-10, 2.17e-10 (8.4) | 3e-7: 1.89e-10, 2.17e-10 |
| w3-kd100 | 4031442 | 1e-8 | 1.18e-09, 1.82e-09 | fastker | 3e-7: 1.88e-10, 2.16e-10 (8.4) | 3e-7: 1.88e-10, 2.16e-10 |
| w3-n1000000-kd200 | 1000000 | 1e-3 | 3.25e-04, 4.89e-04 | default | 3e-3: 1.68e-03, 1.62e-03 (5.2) | 3e-4: 4.64e-07, 6.17e-07 |
| w3-n1000000-kd200 | 1000000 | 1e-3 | 3.25e-04, 4.89e-04 | fastker | 3e-3: 1.68e-03, 1.61e-03 (5.2) | 3e-4: 4.62e-07, 6.14e-07 |
| w3-n1000000-kd200 | 1000000 | 1e-6 | 1.70e-07, 2.73e-07 | default | 1e-4: 4.64e-07, 6.17e-07 (2.7) | 3e-7: 2.13e-10, 3.44e-10 |
| w3-n1000000-kd200 | 1000000 | 1e-6 | 1.70e-07, 2.73e-07 | fastker | 1e-4: 4.62e-07, 6.14e-07 (2.7) | 3e-7: 2.13e-10, 3.44e-10 |
| w3-n1000000-kd200 | 1000000 | 1e-8 | 1.43e-09, 2.64e-09 | default | 3e-7: 2.13e-10, 3.44e-10 (7.7) | 3e-7: 2.13e-10, 3.44e-10 |
| w3-n1000000-kd200 | 1000000 | 1e-8 | 1.43e-09, 2.64e-09 | fastker | 3e-7: 2.13e-10, 3.44e-10 (7.7) | 3e-7: 2.13e-10, 3.44e-10 |
| w3-n1000000-kd500 | 1000000 | 1e-3 | 3.50e-04, 4.50e-04 | default | 3e-3: 2.50e-04, 2.54e-04 (1.8) | 3e-3: 2.50e-04, 2.54e-04 |
| w3-n1000000-kd500 | 1000000 | 1e-3 | 3.50e-04, 4.50e-04 | fastker | 3e-3: 2.31e-04, 2.27e-04 (2.0) | 3e-3: 2.31e-04, 2.27e-04 |
| w3-n1000000-kd500 | 1000000 | 1e-6 | 1.91e-07, 2.53e-07 | default | 3e-5: 9.44e-08, 1.15e-07 (2.2) | 3e-4: 9.43e-08, 1.14e-07 |
| w3-n1000000-kd500 | 1000000 | 1e-6 | 1.91e-07, 2.53e-07 | fastker | 3e-5: 8.75e-08, 1.01e-07 (2.5) | 3e-4: 8.73e-08, 1.01e-07 |
| w3-n1000000-kd500 | 1000000 | 1e-8 | 2.04e-09, 3.42e-09 | default | 3e-7: 8.68e-11, 1.00e-10 (34.1) | 3e-7: 8.68e-11, 1.00e-10 |
| w3-n1000000-kd500 | 1000000 | 1e-8 | 2.04e-09, 3.42e-09 | fastker | 3e-7: 8.68e-11, 1.00e-10 (34.1) | 3e-7: 8.68e-11, 1.00e-10 |
| w1-kd500-planewave | 1989437 | 1e-3 | 2.11e-04, 4.19e-04 | default | 3e-3: 9.21e-05, 1.10e-04 (3.8) | 3e-3: 9.21e-05, 1.10e-04 |
| w1-kd500-planewave | 1989437 | 1e-6 | 6.59e-08, 2.18e-07 | default | 3e-5: 3.67e-08, 5.92e-08 (3.7) | 3e-4: 3.67e-08, 5.89e-08 |
| w1-kd500-planewave | 1989437 | 1e-8 | 6.82e-10, 1.19e-09 | default | 1e-7: 2.64e-11, 4.32e-11 (27.5) | 3e-7: 2.64e-11, 4.32e-11 |
| w1-kd1000-planewave | 7957748 | 1e-3 | 2.53e-04, 4.59e-04 | fastker | 1e-3: 8.43e-05, 9.64e-05 (4.8) | 1e-3: 8.43e-05, 9.64e-05 |
| w1-kd1000-planewave | 7957748 | 1e-6 | 5.25e-08, 2.69e-07 | fastker | 3e-4: 2.45e-08, 3.61e-08 (7.5) | 1e-4: 2.46e-08, 3.60e-08 |
| w1-kd1000-planewave | 7957748 | 1e-8 | 4.50e-10, 9.45e-10 | fastker | 3e-7: 1.69e-11, 2.61e-11 (36.2) | 3e-7: 1.69e-11, 2.61e-11 |
| w1-kd10-constant | 796 | 1e-3 | 9.45e-16, 1.42e-15 | default | 1e-4: 9.10e-16, 1.22e-15 (1.2) | 1e-4: 9.10e-16, 1.22e-15 |
| w1-kd10-constant | 796 | 1e-6 | 9.45e-16, 1.42e-15 | default | 1e-4: 9.10e-16, 1.22e-15 (1.2) | 1e-4: 9.10e-16, 1.22e-15 |
| w1-kd10-constant | 796 | 1e-8 | 9.45e-16, 1.42e-15 | default | 1e-4: 9.10e-16, 1.22e-15 (1.2) | 1e-4: 9.10e-16, 1.22e-15 |
| w1-kd30-constant | 7162 | 1e-3 | 9.07e-05, 3.05e-04 | default | 3e-3: 1.13e-04, 2.21e-04 (1.4) | 3e-4: 1.46e-08, 2.34e-08 |
| w1-kd30-constant | 7162 | 1e-6 | 4.99e-08, 2.03e-07 | default | 3e-4: 1.46e-08, 2.34e-08 (8.7) | 3e-4: 1.46e-08, 2.34e-08 |
| w1-kd30-constant | 7162 | 1e-8 | 7.88e-13, 2.01e-12 | default | 3e-9: 9.56e-12, 1.42e-11 (12.1) | >3e-9: 9.56e-12, 1.42e-11 |
| w1-kd100-constant | 79578 | 1e-3 | 2.64e-04, 4.41e-04 | default | 3e-3: 2.62e-03, 2.55e-03 (9.9) | 3e-5: 6.33e-07, 9.90e-07 |
| w1-kd100-constant | 79578 | 1e-6 | 2.01e-07, 4.41e-07 | default | 3e-4: 6.32e-07, 9.96e-07 (3.1) | 3e-7: 2.81e-10, 5.17e-10 |
| w1-kd100-constant | 79578 | 1e-8 | 3.22e-09, 5.48e-09 | default | 3e-7: 2.81e-10, 5.17e-10 (11.5) | 3e-7: 2.81e-10, 5.17e-10 |
| w2-kd10-constant | 1276 | 1e-3 | 1.14e-15, 4.58e-15 | default | 1e-7: 9.69e-16, 3.33e-15 (1.4) | 1e-7: 9.69e-16, 3.33e-15 |
| w2-kd10-constant | 1276 | 1e-6 | 1.14e-15, 4.58e-15 | default | 1e-7: 9.69e-16, 3.33e-15 (1.4) | 1e-7: 9.69e-16, 3.33e-15 |
| w2-kd10-constant | 1276 | 1e-8 | 1.14e-15, 4.58e-15 | default | 1e-7: 9.69e-16, 3.33e-15 (1.4) | 1e-7: 9.69e-16, 3.33e-15 |
| w2-kd30-constant | 11480 | 1e-3 | 8.04e-05, 6.37e-04 | default | 3e-3: 6.72e-05, 3.90e-04 (1.6) | 3e-3: 6.72e-05, 3.90e-04 |
| w2-kd30-constant | 11480 | 1e-6 | 2.94e-08, 1.74e-07 | default | 3e-4: 1.28e-08, 4.64e-08 (3.8) | 3e-4: 1.28e-08, 4.64e-08 |
| w2-kd30-constant | 11480 | 1e-8 | 8.51e-12, 1.59e-11 | default | 3e-9: 8.70e-12, 1.70e-11 (1.1) | >1e-9: 8.70e-12, 1.70e-11 |
| w2-kd100-constant | 127554 | 1e-3 | 1.28e-04, 6.35e-04 | default | 3e-3: 1.35e-03, 3.59e-03 (10.5) | 3e-4: 2.64e-07, 1.11e-06 |
| w2-kd100-constant | 127554 | 1e-6 | 1.11e-07, 7.56e-07 | default | 1e-4: 2.64e-07, 1.11e-06 (2.4) | 3e-7: 9.57e-11, 4.24e-10 |
| w2-kd100-constant | 127554 | 1e-8 | 8.83e-10, 3.69e-09 | default | 3e-8: 9.57e-11, 4.24e-10 (9.2) | 3e-7: 9.57e-11, 4.24e-10 |
| w3-kd10-constant | 4032 | 1e-3 | 8.02e-06, 4.31e-05 | default | 3e-3: 1.86e-05, 3.02e-05 (2.3) | 3e-4: 1.43e-08, 6.69e-08 |
| w3-kd10-constant | 4032 | 1e-6 | 1.39e-08, 6.40e-08 | default | 3e-4: 1.43e-08, 6.69e-08 (1.0) | 3e-5: 1.24e-08, 2.28e-08 |
| w3-kd10-constant | 4032 | 1e-8 | 2.37e-15, 1.27e-14 | default | 1e-8: 1.90e-11, 2.07e-11 (8025.3) | >1e-8: 1.90e-11, 2.07e-11 |
| w3-kd30-constant | 108849 | 1e-3 | 1.28e-04, 4.00e-04 | default | 3e-3: 1.92e-04, 2.72e-04 (1.5) | 3e-4: 3.16e-08, 8.12e-08 |
| w3-kd30-constant | 108849 | 1e-6 | 4.77e-08, 1.91e-07 | default | 3e-4: 3.16e-08, 8.12e-08 (2.4) | 3e-4: 3.16e-08, 8.12e-08 |
| w3-kd30-constant | 108849 | 1e-8 | 2.38e-11, 2.92e-11 | default | 3e-8: 2.38e-11, 3.04e-11 (1.0) | >3e-8: 2.38e-11, 3.04e-11 |

</details>

## 3. The comparison

### 3.1 Summary

Kathleen (40 cores; ours' splits 1 × 40 … 8 × 5; FMM3D on 40 threads). "ours slower" lists the cases
whose best split loses.

| rule | cases | ε | n | best split ÷ FMM3D: gmean (range) | 1 × N ÷ FMM3D: gmean (range) | evaluate ÷ (call − tree): gmean (range) [cases, where fewer: the tree below the call] | ours slower (best split) |
| --- | --- | --- | ---: | --- | --- | --- | --- |
| closest | all | 1e-3 | 26 | 0.41 (0.06–4.37) | 0.52 (0.08–4.37) | 0.14 (0.03–0.74) | w1-kd30, w2-kd10, w2-kd30, w3-kd10, w1-kd30-constant, w2-kd30-constant, w3-kd10-constant |
| closest | all | 1e-6 | 26 | 0.31 (0.07–1.26) | 0.39 (0.09–1.62) | 0.15 (0.04–0.57) | w2-kd10 |
| closest | all | 1e-8 | 26 | 0.34 (0.09–1.41) | 0.40 (0.10–1.78) | 0.20 (0.06–0.55) | w2-kd10, w3-kd10 |
| closest | N ≥ 10⁵ | 1e-3 | 14 | 0.20 (0.06–0.97) | 0.26 (0.08–1.16) | 0.10 (0.03–0.41) | – |
| closest | N ≥ 10⁵ | 1e-6 | 14 | 0.19 (0.07–0.70) | 0.24 (0.09–0.90) | 0.12 (0.04–0.51) | – |
| closest | N ≥ 10⁵ | 1e-8 | 14 | 0.23 (0.09–0.64) | 0.27 (0.10–0.78) | 0.18 (0.06–0.52) | – |
| closest | N ≥ 10⁶ | 1e-3 | 8 | 0.13 (0.06–0.48) | 0.18 (0.08–0.93) | 0.08 (0.03–0.37) | – |
| closest | N ≥ 10⁶ | 1e-6 | 8 | 0.14 (0.07–0.50) | 0.18 (0.09–0.61) | 0.09 (0.04–0.40) | – |
| closest | N ≥ 10⁶ | 1e-8 | 8 | 0.18 (0.09–0.56) | 0.20 (0.10–0.64) | 0.13 (0.06–0.47) | – |
| one-sided | all | 1e-3 | 26 | 0.38 (0.06–2.68) | 0.47 (0.08–3.40) | 0.13 (0.03–1.10) | w1-kd10, w1-kd30, w2-kd10, w2-kd30, w3-kd10, w1-kd30-constant, w2-kd30-constant |
| one-sided | all | 1e-6 | 26 | 0.30 (0.07–2.55) | 0.38 (0.09–3.45) | 0.15 (0.04–1.10) | w1-kd10, w2-kd10, w3-kd10 |
| one-sided | all | 1e-8 | 26 | 0.37 (0.09–2.55) | 0.44 (0.10–3.10) | 0.22 (0.06–1.20) | w1-kd10, w2-kd10, w3-kd10 |
| one-sided | N ≥ 10⁵ | 1e-3 | 14 | 0.17 (0.06–1.00) | 0.22 (0.08–1.20) | 0.08 (0.03–0.43) | – |
| one-sided | N ≥ 10⁵ | 1e-6 | 14 | 0.17 (0.07–0.70) | 0.22 (0.09–0.90) | 0.11 (0.04–0.51) | – |
| one-sided | N ≥ 10⁵ | 1e-8 | 14 | 0.23 (0.09–0.64) | 0.27 (0.10–0.78) | 0.18 (0.06–0.52) | – |
| one-sided | N ≥ 10⁶ | 1e-3 | 8 | 0.11 (0.06–0.41) | 0.16 (0.08–0.53) | 0.06 (0.03–0.31) | – |
| one-sided | N ≥ 10⁶ | 1e-6 | 8 | 0.13 (0.07–0.38) | 0.17 (0.09–0.47) | 0.08 (0.04–0.30) | – |
| one-sided | N ≥ 10⁶ | 1e-8 | 8 | 0.17 (0.09–0.56) | 0.20 (0.10–0.64) | 0.13 (0.06–0.47) | – |

locust (72 cores):

| rule | cases | ε | n | best split ÷ FMM3D: gmean (range) | 1 × N ÷ FMM3D: gmean (range) | evaluate ÷ (call − tree): gmean (range) [cases, where fewer: the tree below the call] | ours slower (best split) |
| --- | --- | --- | ---: | --- | --- | --- | --- |
| closest | all | 1e-3 | 26 | 0.25 (0.08–0.89) | 0.43 (0.12–1.18) | 0.12 (0.03–1.60) [25] | – |
| closest | all | 1e-6 | 26 | 0.26 (0.08–0.75) | 0.41 (0.13–1.12) | 0.14 (0.04–0.72) [25] | – |
| closest | all | 1e-8 | 26 | 0.26 (0.09–0.64) | 0.32 (0.12–1.47) | 0.18 (0.06–0.35) [25] | – |
| closest | N ≥ 10⁵ | 1e-3 | 14 | 0.17 (0.08–0.43) | 0.34 (0.12–1.09) | 0.09 (0.03–0.32) | – |
| closest | N ≥ 10⁵ | 1e-6 | 14 | 0.16 (0.08–0.33) | 0.27 (0.13–0.66) | 0.10 (0.04–0.30) | – |
| closest | N ≥ 10⁵ | 1e-8 | 14 | 0.19 (0.09–0.37) | 0.23 (0.12–0.51) | 0.15 (0.06–0.35) | – |
| closest | N ≥ 10⁶ | 1e-3 | 8 | 0.13 (0.08–0.36) | 0.27 (0.12–0.82) | 0.07 (0.03–0.32) | – |
| closest | N ≥ 10⁶ | 1e-6 | 8 | 0.14 (0.08–0.33) | 0.24 (0.13–0.66) | 0.08 (0.04–0.30) | – |
| closest | N ≥ 10⁶ | 1e-8 | 8 | 0.16 (0.09–0.37) | 0.21 (0.12–0.51) | 0.13 (0.06–0.35) | – |
| one-sided | all | 1e-3 | 26 | 0.25 (0.08–0.88) | 0.42 (0.12–1.25) | 0.11 (0.03–0.40) [25] | – |
| one-sided | all | 1e-6 | 26 | 0.25 (0.08–1.16) | 0.39 (0.14–1.62) | 0.12 (0.04–0.45) [25] | w2-kd10 |
| one-sided | all | 1e-8 | 26 | 0.25 (0.09–0.80) | 0.31 (0.11–1.47) | 0.17 (0.05–0.40) [25] | – |
| one-sided | N ≥ 10⁵ | 1e-3 | 14 | 0.15 (0.08–0.41) | 0.30 (0.12–0.91) | 0.08 (0.03–0.23) | – |
| one-sided | N ≥ 10⁵ | 1e-6 | 14 | 0.15 (0.08–0.33) | 0.26 (0.14–0.57) | 0.09 (0.04–0.27) | – |
| one-sided | N ≥ 10⁵ | 1e-8 | 14 | 0.19 (0.09–0.37) | 0.23 (0.11–0.51) | 0.15 (0.06–0.35) | – |
| one-sided | N ≥ 10⁶ | 1e-3 | 8 | 0.13 (0.08–0.27) | 0.25 (0.12–0.60) | 0.07 (0.03–0.23) | – |
| one-sided | N ≥ 10⁶ | 1e-6 | 8 | 0.14 (0.08–0.28) | 0.24 (0.14–0.55) | 0.08 (0.04–0.25) | – |
| one-sided | N ≥ 10⁶ | 1e-8 | 8 | 0.16 (0.09–0.37) | 0.21 (0.11–0.51) | 0.13 (0.06–0.35) | – |

### 3.2 Kathleen, per case

Times in seconds. "best split": its build, evaluation (the second) and their sum; "1 × N": the same at
1 × 40. FMM3D (closest): variant, eps (the factor 10^d), the call and its tree; then ours ÷ FMM3D at
the best split and at 1 × N, and our evaluation over the call less its tree; then the one-sided
match and its ratios.

| Case | N | ε | best split | build s | evaluate s | ours s | 1 × N: build s, evaluate s, s | FMM3D (closest): variant, eps (factor), call s, tree s | ours ÷ FMM3D: best, 1 × N | evaluate ÷ (call − tree) | FMM3D (one-sided): variant, eps, call s | ÷: best, 1 × N |
| --- | ---: | --- | --- | ---: | ---: | ---: | --- | --- | --- | ---: | --- | --- |
| w1-kd10 | 796 | 1e-3 | 2x20 | 0.004 | 0.001 | 0.005 | 0.005, 0.002, 0.007 | default, 1e-4 (1.3), 0.007, 0.001 | 0.74, 0.97 | 0.18 | fastker, 1e-3, 0.002 | 2.60, 3.40 |
| w1-kd10 | 796 | 1e-6 | 2x20 | 0.004 | 0.001 | 0.005 | 0.005, 0.002, 0.007 | default, 1e-4 (1.3), 0.007, 0.001 | 0.73, 0.99 | 0.18 | fastker, 1e-3, 0.002 | 2.55, 3.45 |
| w1-kd10 | 796 | 1e-8 | 2x20 | 0.004 | 0.001 | 0.005 | 0.005, 0.002, 0.006 | default, 1e-4 (1.3), 0.007, 0.001 | 0.73, 0.89 | 0.20 | fastker, 1e-3, 0.002 | 2.55, 3.10 |
| w1-kd30 | 7162 | 1e-3 | 1x40 | 0.081 | 0.015 | 0.096 | 0.081, 0.015, 0.096 | fastker, 3e-3 (1.5), 0.022, 0.002 | 4.37, 4.37 | 0.74 | fastker, 1e-4, 0.053 | 1.81, 1.81 |
| w1-kd30 | 7162 | 1e-6 | 2x20 | 0.018 | 0.008 | 0.026 | 0.019, 0.011, 0.030 | default, 3e-4 (9.1), 0.058, 0.002 | 0.45, 0.52 | 0.15 | fastker, 1e-4, 0.053 | 0.49, 0.57 |
| w1-kd30 | 7162 | 1e-8 | 2x20 | 0.022 | 0.024 | 0.046 | 0.025, 0.029, 0.054 | default, 1e-8 (6.4), 0.109, 0.002 | 0.43, 0.50 | 0.22 | fastker, 1e-4, 0.053 | 0.88, 1.02 |
| w1-kd100 | 79578 | 1e-3 | 2x20 | 0.079 | 0.060 | 0.139 | 0.099, 0.078, 0.177 | fastker, 3e-3 (8.0), 0.391, 0.008 | 0.36, 0.45 | 0.16 | fastker, 3e-4, 0.411 | 0.34, 0.43 |
| w1-kd100 | 79578 | 1e-6 | 2x20 | 0.103 | 0.085 | 0.188 | 0.116, 0.108, 0.224 | fastker, 1e-5 (2.3), 0.463, 0.006 | 0.41, 0.48 | 0.19 | fastker, 3e-7, 0.580 | 0.32, 0.39 |
| w1-kd100 | 79578 | 1e-8 | 2x20 | 0.113 | 0.205 | 0.318 | 0.124, 0.241, 0.365 | fastker, 3e-7 (4.7), 0.580, 0.006 | 0.55, 0.63 | 0.36 | fastker, 3e-7, 0.580 | 0.55, 0.63 |
| w1-kd250 | 497360 | 1e-3 | 2x20 | 0.379 | 0.381 | 0.760 | 0.529, 0.468, 0.997 | fastker, 3e-3 (1.6), 4.256, 0.040 | 0.18, 0.23 | 0.09 | fastker, 3e-3, 4.256 | 0.18, 0.23 |
| w1-kd250 | 497360 | 1e-6 | 2x20 | 0.460 | 0.555 | 1.015 | 0.594, 0.642, 1.236 | fastker, 1e-5 (3.2), 4.891, 0.032 | 0.21, 0.25 | 0.11 | fastker, 3e-4, 4.374 | 0.23, 0.28 |
| w1-kd250 | 497360 | 1e-8 | 2x20 | 0.519 | 1.250 | 1.768 | 0.666, 1.445, 2.112 | fastker, 3e-7 (24.2), 5.938, 0.032 | 0.30, 0.36 | 0.21 | fastker, 3e-7, 5.938 | 0.30, 0.36 |
| w1-kd500 | 1989437 | 1e-3 | 2x20 | 1.773 | 1.762 | 3.535 | 2.447, 2.032, 4.479 | fastker, 3e-3 (2.3), 31.625, 0.275 | 0.11, 0.14 | 0.06 | fastker, 3e-3, 31.625 | 0.11, 0.14 |
| w1-kd500 | 1989437 | 1e-6 | 2x20 | 1.997 | 2.538 | 4.535 | 2.716, 2.901, 5.617 | fastker, 3e-4 (4.2), 33.563, 0.231 | 0.14, 0.17 | 0.08 | fastker, 3e-4, 33.563 | 0.14, 0.17 |
| w1-kd500 | 1989437 | 1e-8 | 2x20 | 2.184 | 5.320 | 7.504 | 2.850, 5.642, 8.492 | fastker, 1e-8 (23.5), 48.023, 0.237 | 0.16, 0.18 | 0.11 | fastker, 3e-7, 44.121 | 0.17, 0.19 |
| w1-kd1000 | 7957748 | 1e-3 | 2x20 | 8.787 | 8.994 | 17.781 | 12.692, 10.441, 23.133 | fastker, 3e-3 (2.3), 280.074, 1.621 | 0.06, 0.08 | 0.03 | fastker, 3e-3, 280.074 | 0.06, 0.08 |
| w1-kd1000 | 7957748 | 1e-6 | 2x20 | 9.897 | 12.608 | 22.506 | 13.724, 14.166, 27.890 | fastker, 3e-4 (4.4), 310.586, 1.406 | 0.07, 0.09 | 0.04 | fastker, 1e-4, 310.205 | 0.07, 0.09 |
| w1-kd1000 | 7957748 | 1e-8 | 2x20 | 10.567 | 23.900 | 34.468 | 14.267, 25.143, 39.410 | fastker, 3e-5 (33.0), 317.425, 1.401 | 0.11, 0.12 | 0.08 | fastker, 3e-7, 391.382 | 0.09, 0.10 |
| w2-kd10 | 1276 | 1e-3 | 4x10 | 0.005 | 0.002 | 0.007 | 0.005, 0.003, 0.008 | fastker, 1e-3 (1.7), 0.005, 0.001 | 1.34, 1.66 | 0.52 | fastker, 1e-3, 0.005 | 1.34, 1.66 |
| w2-kd10 | 1276 | 1e-6 | 2x20 | 0.004 | 0.002 | 0.006 | 0.005, 0.003, 0.008 | fastker, 1e-3 (1.7), 0.005, 0.001 | 1.26, 1.62 | 0.57 | fastker, 1e-3, 0.005 | 1.26, 1.62 |
| w2-kd10 | 1276 | 1e-8 | 2x20 | 0.004 | 0.002 | 0.007 | 0.005, 0.003, 0.008 | fastker, 1e-3 (1.7), 0.005, 0.001 | 1.32, 1.60 | 0.55 | fastker, 1e-3, 0.005 | 1.32, 1.60 |
| w2-kd30 | 11480 | 1e-3 | 2x20 | 0.111 | 0.017 | 0.127 | 0.112, 0.023, 0.135 | fastker, 1e-3 (5.3), 0.064, 0.003 | 1.98, 2.11 | 0.27 | fastker, 3e-3, 0.058 | 2.19, 2.33 |
| w2-kd30 | 11480 | 1e-6 | 2x20 | 0.032 | 0.024 | 0.056 | 0.034, 0.035, 0.069 | fastker, 3e-4 (13.3), 0.106, 0.002 | 0.53, 0.65 | 0.23 | fastker, 3e-4, 0.106 | 0.53, 0.65 |
| w2-kd30 | 11480 | 1e-8 | 2x20 | 0.037 | 0.072 | 0.109 | 0.039, 0.091, 0.130 | fastker, 1e-9 (1.8), 0.246, 0.002 | 0.44, 0.53 | 0.30 | fastker, >1e-9, 0.246 | 0.44, 0.53 |
| w2-kd100 | 127554 | 1e-3 | 2x20 | 0.129 | 0.083 | 0.212 | 0.163, 0.108, 0.272 | fastker, 3e-3 (5.4), 0.538, 0.012 | 0.39, 0.51 | 0.16 | fastker, 3e-4, 0.689 | 0.31, 0.39 |
| w2-kd100 | 127554 | 1e-6 | 2x20 | 0.138 | 0.138 | 0.276 | 0.169, 0.161, 0.329 | fastker, 3e-4 (1.5), 0.689, 0.011 | 0.40, 0.48 | 0.20 | fastker, 3e-7, 1.120 | 0.25, 0.29 |
| w2-kd100 | 127554 | 1e-8 | 2x20 | 0.141 | 0.399 | 0.539 | 0.173, 0.433, 0.606 | fastker, 3e-7 (4.9), 1.120, 0.011 | 0.48, 0.54 | 0.36 | fastker, 3e-7, 1.120 | 0.48, 0.54 |
| w2-kd250 | 797211 | 1e-3 | 2x20 | 0.577 | 0.502 | 1.079 | 0.849, 0.627, 1.477 | fastker, 3e-3 (2.8), 5.789, 0.089 | 0.19, 0.26 | 0.09 | fastker, 3e-3, 5.789 | 0.19, 0.26 |
| w2-kd250 | 797211 | 1e-6 | 2x20 | 0.642 | 0.838 | 1.480 | 0.924, 0.903, 1.827 | fastker, 3e-4 (6.0), 6.077, 0.071 | 0.24, 0.30 | 0.14 | fastker, 3e-4, 6.077 | 0.24, 0.30 |
| w2-kd250 | 797211 | 1e-8 | 2x20 | 0.705 | 2.003 | 2.708 | 0.975, 2.121, 3.096 | fastker, 1e-7 (24.6), 9.266, 0.072 | 0.29, 0.33 | 0.22 | fastker, 3e-7, 8.923 | 0.30, 0.35 |
| w2-kd500 | 3188840 | 1e-3 | 2x20 | 2.622 | 2.224 | 4.846 | 4.040, 2.507, 6.547 | fastker, 3e-3 (3.1), 38.086, 0.624 | 0.13, 0.17 | 0.06 | fastker, 3e-3, 38.086 | 0.13, 0.17 |
| w2-kd500 | 3188840 | 1e-6 | 2x20 | 2.838 | 3.631 | 6.469 | 4.257, 3.991, 8.248 | fastker, 3e-4 (7.8), 40.063, 0.530 | 0.16, 0.21 | 0.09 | fastker, 3e-4, 40.063 | 0.16, 0.21 |
| w2-kd500 | 3188840 | 1e-8 | 2x20 | 3.045 | 7.968 | 11.012 | 4.426, 8.304, 12.730 | fastker, 3e-7 (18.4), 55.417, 0.525 | 0.20, 0.23 | 0.15 | fastker, 3e-7, 55.417 | 0.20, 0.23 |
| w3-kd10 | 4032 | 1e-3 | 2x20 | 0.021 | 0.003 | 0.024 | 0.023, 0.005, 0.029 | fastker, 3e-3 (2.9), 0.009, 0.002 | 2.68, 3.19 | 0.43 | fastker, 3e-3, 0.009 | 2.68, 3.19 |
| w3-kd10 | 4032 | 1e-6 | 2x20 | 0.016 | 0.006 | 0.022 | 0.019, 0.009, 0.028 | default, 3e-4 (1.1), 0.026, 0.002 | 0.84, 1.07 | 0.26 | fastker, 1e-4, 0.020 | 1.09, 1.39 |
| w3-kd10 | 4032 | 1e-8 | 2x20 | 0.019 | 0.009 | 0.028 | 0.023, 0.012, 0.036 | fastker, 1e-4 (1.6), 0.020, 0.001 | 1.41, 1.78 | 0.47 | fastker, 1e-4, 0.020 | 1.41, 1.78 |
| w3-kd30 | 108849 | 1e-3 | 2x20 | 0.122 | 0.079 | 0.201 | 0.137, 0.104, 0.241 | fastker, 1e-3 (4.0), 0.207, 0.017 | 0.97, 1.16 | 0.41 | fastker, 3e-3, 0.201 | 1.00, 1.20 |
| w3-kd30 | 108849 | 1e-6 | 4x10 | 0.055 | 0.115 | 0.170 | 0.081, 0.136, 0.217 | fastker, 3e-4 (9.7), 0.242, 0.018 | 0.70, 0.90 | 0.51 | fastker, 3e-4, 0.242 | 0.70, 0.90 |
| w3-kd30 | 108849 | 1e-8 | 2x20 | 0.066 | 0.243 | 0.309 | 0.087, 0.292, 0.379 | fastker, 1e-9 (2.1), 0.485, 0.017 | 0.64, 0.78 | 0.52 | fastker, >1e-9, 0.485 | 0.64, 0.78 |
| w3-kd100 | 4031442 | 1e-3 | 8x5 | 0.803 | 2.177 | 2.979 | 3.995, 2.545, 6.540 | fastker, 3e-3 (6.4), 6.996, 0.993 | 0.43, 0.93 | 0.36 | fastker, 3e-4, 20.085 | 0.15, 0.33 |
| w3-kd100 | 4031442 | 1e-6 | 8x5 | 0.756 | 3.784 | 4.541 | 3.976, 4.027, 8.003 | fastker, 3e-4 (2.1), 20.085, 0.889 | 0.23, 0.40 | 0.20 | fastker, 3e-7, 24.943 | 0.18, 0.32 |
| w3-kd100 | 4031442 | 1e-8 | 8x5 | 0.815 | 9.738 | 10.553 | 3.891, 9.410, 13.301 | fastker, 3e-7 (8.4), 24.943, 0.896 | 0.42, 0.53 | 0.40 | fastker, 3e-7, 24.943 | 0.42, 0.53 |
| w3-n1000000-kd200 | 1000000 | 1e-3 | 4x10 | 0.472 | 1.304 | 1.776 | 0.861, 1.428, 2.289 | fastker, 3e-3 (5.2), 3.711, 0.162 | 0.48, 0.62 | 0.37 | fastker, 3e-4, 4.323 | 0.41, 0.53 |
| w3-n1000000-kd200 | 1000000 | 1e-6 | 4x10 | 0.515 | 1.704 | 2.219 | 0.902, 1.808, 2.711 | fastker, 1e-4 (2.7), 4.428, 0.163 | 0.50, 0.61 | 0.40 | fastker, 3e-7, 5.784 | 0.38, 0.47 |
| w3-n1000000-kd200 | 1000000 | 1e-8 | 2x20 | 0.620 | 2.617 | 3.236 | 0.928, 2.759, 3.687 | fastker, 3e-7 (7.7), 5.784, 0.164 | 0.56, 0.64 | 0.47 | fastker, 3e-7, 5.784 | 0.56, 0.64 |
| w3-n1000000-kd500 | 1000000 | 1e-3 | 2x20 | 1.983 | 3.034 | 5.017 | 2.949, 3.745, 6.694 | fastker, 3e-3 (2.0), 50.332, 0.163 | 0.10, 0.13 | 0.06 | fastker, 3e-3, 50.332 | 0.10, 0.13 |
| w3-n1000000-kd500 | 1000000 | 1e-6 | 2x20 | 2.122 | 4.015 | 6.136 | 3.071, 5.147, 8.218 | fastker, 3e-5 (2.5), 61.909, 0.168 | 0.10, 0.13 | 0.07 | fastker, 3e-4, 55.908 | 0.11, 0.15 |
| w3-n1000000-kd500 | 1000000 | 1e-8 | 2x20 | 2.173 | 6.313 | 8.485 | 3.086, 6.721, 9.807 | fastker, 3e-7 (34.1), 75.353, 0.165 | 0.11, 0.13 | 0.08 | fastker, 3e-7, 75.353 | 0.11, 0.13 |
| w1-kd500-planewave | 1989437 | 1e-3 | 2x20 | 1.732 | 1.770 | 3.502 | 2.492, 2.049, 4.541 | default, 3e-3 (3.8), 46.245, 0.360 | 0.08, 0.10 | 0.04 | default, 3e-3, 46.245 | 0.08, 0.10 |
| w1-kd500-planewave | 1989437 | 1e-6 | 2x20 | 1.979 | 2.517 | 4.497 | 2.681, 2.828, 5.509 | default, 3e-5 (3.7), 52.194, 0.292 | 0.09, 0.11 | 0.05 | default, 3e-4, 47.792 | 0.09, 0.12 |
| w1-kd500-planewave | 1989437 | 1e-8 | 2x20 | 2.171 | 5.306 | 7.477 | 2.855, 5.610, 8.464 | default, 1e-7 (27.5), 64.778, 0.274 | 0.12, 0.13 | 0.08 | default, 3e-7, 63.579 | 0.12, 0.13 |
| w1-kd1000-planewave | 7957748 | 1e-3 | 2x20 | 8.858 | 8.931 | 17.789 | 12.551, 10.157, 22.708 | fastker, 1e-3 (4.8), 279.358, 1.633 | 0.06, 0.08 | 0.03 | fastker, 1e-3, 279.358 | 0.06, 0.08 |
| w1-kd1000-planewave | 7957748 | 1e-6 | 2x20 | 10.005 | 12.472 | 22.478 | 13.702, 14.236, 27.937 | fastker, 3e-4 (7.5), 304.303, 1.393 | 0.07, 0.09 | 0.04 | fastker, 1e-4, 310.323 | 0.07, 0.09 |
| w1-kd1000-planewave | 7957748 | 1e-8 | 2x20 | 10.771 | 23.893 | 34.664 | 14.299, 25.888, 40.187 | fastker, 3e-7 (36.2), 386.673, 1.412 | 0.09, 0.10 | 0.06 | fastker, 3e-7, 386.673 | 0.09, 0.10 |
| w1-kd10-constant | 796 | 1e-3 | 2x20 | 0.004 | 0.001 | 0.005 | 0.005, 0.002, 0.007 | default, 1e-4 (1.2), 0.007, 0.001 | 0.77, 0.96 | 0.20 | default, 1e-4, 0.007 | 0.77, 0.96 |
| w1-kd10-constant | 796 | 1e-6 | 2x20 | 0.004 | 0.001 | 0.005 | 0.005, 0.002, 0.006 | default, 1e-4 (1.2), 0.007, 0.001 | 0.77, 0.91 | 0.20 | default, 1e-4, 0.007 | 0.77, 0.91 |
| w1-kd10-constant | 796 | 1e-8 | 2x20 | 0.004 | 0.001 | 0.005 | 0.005, 0.002, 0.006 | default, 1e-4 (1.2), 0.007, 0.001 | 0.74, 0.91 | 0.18 | default, 1e-4, 0.007 | 0.74, 0.91 |
| w1-kd30-constant | 7162 | 1e-3 | 2x20 | 0.082 | 0.010 | 0.091 | 0.083, 0.014, 0.096 | default, 3e-3 (1.4), 0.050, 0.003 | 1.83, 1.92 | 0.21 | default, 3e-4, 0.062 | 1.47, 1.55 |
| w1-kd30-constant | 7162 | 1e-6 | 2x20 | 0.018 | 0.009 | 0.026 | 0.019, 0.011, 0.030 | default, 3e-4 (8.7), 0.062, 0.002 | 0.43, 0.49 | 0.14 | default, 3e-4, 0.062 | 0.43, 0.49 |
| w1-kd30-constant | 7162 | 1e-8 | 2x20 | 0.022 | 0.023 | 0.045 | 0.028, 0.027, 0.055 | default, 3e-9 (12.1), 0.115, 0.002 | 0.39, 0.47 | 0.20 | default, >3e-9, 0.115 | 0.39, 0.47 |
| w1-kd100-constant | 79578 | 1e-3 | 2x20 | 0.081 | 0.059 | 0.140 | 0.099, 0.079, 0.178 | default, 3e-3 (9.9), 0.866, 0.013 | 0.16, 0.21 | 0.07 | default, 3e-5, 0.772 | 0.18, 0.23 |
| w1-kd100-constant | 79578 | 1e-6 | 2x20 | 0.102 | 0.085 | 0.187 | 0.118, 0.110, 0.228 | default, 3e-4 (3.1), 0.713, 0.008 | 0.26, 0.32 | 0.12 | default, 3e-7, 1.333 | 0.14, 0.17 |
| w1-kd100-constant | 79578 | 1e-8 | 2x20 | 0.109 | 0.207 | 0.315 | 0.128, 0.246, 0.373 | default, 3e-7 (11.5), 1.333, 0.008 | 0.24, 0.28 | 0.16 | default, 3e-7, 1.333 | 0.24, 0.28 |
| w2-kd10-constant | 1276 | 1e-3 | 4x10 | 0.005 | 0.002 | 0.007 | 0.005, 0.003, 0.008 | default, 1e-7 (1.4), 0.020, 0.001 | 0.33, 0.38 | 0.11 | default, 1e-7, 0.020 | 0.33, 0.38 |
| w2-kd10-constant | 1276 | 1e-6 | 2x20 | 0.004 | 0.002 | 0.007 | 0.006, 0.003, 0.009 | default, 1e-7 (1.4), 0.020, 0.001 | 0.33, 0.43 | 0.12 | default, 1e-7, 0.020 | 0.33, 0.43 |
| w2-kd10-constant | 1276 | 1e-8 | 2x20 | 0.004 | 0.002 | 0.006 | 0.005, 0.003, 0.008 | default, 1e-7 (1.4), 0.020, 0.001 | 0.32, 0.40 | 0.12 | default, 1e-7, 0.020 | 0.32, 0.40 |
| w2-kd30-constant | 11480 | 1e-3 | 2x20 | 0.111 | 0.017 | 0.128 | 0.110, 0.024, 0.134 | default, 3e-3 (1.6), 0.104, 0.004 | 1.23, 1.29 | 0.17 | default, 3e-3, 0.104 | 1.23, 1.29 |
| w2-kd30-constant | 11480 | 1e-6 | 2x20 | 0.031 | 0.025 | 0.056 | 0.034, 0.035, 0.068 | default, 3e-4 (3.8), 0.157, 0.003 | 0.36, 0.43 | 0.16 | default, 3e-4, 0.157 | 0.36, 0.43 |
| w2-kd30-constant | 11480 | 1e-8 | 2x20 | 0.036 | 0.073 | 0.110 | 0.040, 0.099, 0.139 | default, 3e-9 (1.1), 0.279, 0.002 | 0.39, 0.50 | 0.27 | default, >1e-9, 0.295 | 0.37, 0.47 |
| w2-kd100-constant | 127554 | 1e-3 | 2x20 | 0.130 | 0.088 | 0.217 | 0.166, 0.101, 0.266 | default, 3e-3 (10.5), 1.237, 0.020 | 0.18, 0.22 | 0.07 | default, 3e-4, 1.399 | 0.16, 0.19 |
| w2-kd100-constant | 127554 | 1e-6 | 2x20 | 0.136 | 0.139 | 0.274 | 0.168, 0.167, 0.335 | default, 1e-4 (2.4), 1.487, 0.016 | 0.18, 0.23 | 0.09 | default, 3e-7, 2.215 | 0.12, 0.15 |
| w2-kd100-constant | 127554 | 1e-8 | 2x20 | 0.142 | 0.386 | 0.527 | 0.175, 0.416, 0.591 | default, 3e-8 (9.2), 2.469, 0.013 | 0.21, 0.24 | 0.16 | default, 3e-7, 2.215 | 0.24, 0.27 |
| w3-kd10-constant | 4032 | 1e-3 | 2x20 | 0.022 | 0.003 | 0.025 | 0.022, 0.005, 0.027 | default, 3e-3 (2.3), 0.022, 0.002 | 1.14, 1.22 | 0.15 | default, 3e-4, 0.029 | 0.86, 0.92 |
| w3-kd10-constant | 4032 | 1e-6 | 2x20 | 0.016 | 0.006 | 0.022 | 0.018, 0.011, 0.030 | default, 3e-4 (1.0), 0.029, 0.002 | 0.75, 1.02 | 0.22 | default, 3e-5, 0.031 | 0.70, 0.95 |
| w3-kd10-constant | 4032 | 1e-8 | 2x20 | 0.020 | 0.009 | 0.028 | 0.024, 0.013, 0.038 | default, 1e-8 (8025.3), 0.069, 0.002 | 0.41, 0.55 | 0.13 | default, >1e-8, 0.069 | 0.41, 0.55 |
| w3-kd30-constant | 108849 | 1e-3 | 2x20 | 0.124 | 0.081 | 0.204 | 0.142, 0.104, 0.246 | default, 3e-3 (1.5), 0.326, 0.023 | 0.63, 0.75 | 0.27 | default, 3e-4, 0.694 | 0.29, 0.35 |
| w3-kd30-constant | 108849 | 1e-6 | 4x10 | 0.057 | 0.116 | 0.173 | 0.083, 0.142, 0.225 | default, 3e-4 (2.4), 0.694, 0.022 | 0.25, 0.32 | 0.17 | default, 3e-4, 0.694 | 0.25, 0.32 |
| w3-kd30-constant | 108849 | 1e-8 | 2x20 | 0.065 | 0.245 | 0.310 | 0.087, 0.275, 0.363 | default, 3e-8 (1.0), 1.218, 0.018 | 0.25, 0.30 | 0.20 | default, >3e-8, 1.218 | 0.25, 0.30 |

### 3.3 locust, per case

| Case | N | ε | best split | build s | evaluate s | ours s | 1 × N: build s, evaluate s, s | FMM3D (closest): variant, eps (factor), call s, tree s | ours ÷ FMM3D: best, 1 × N | evaluate ÷ (call − tree) | FMM3D (one-sided): variant, eps, call s | ÷: best, 1 × N |
| --- | ---: | --- | --- | ---: | ---: | ---: | --- | --- | --- | ---: | --- | --- |
| w1-kd10 | 796 | 1e-3 | 4x18 | 0.002 | 0.002 | 0.003 | 0.003, 0.002, 0.005 | default, 1e-4 (1.3), 0.005, 0.000 | 0.66, 1.00 | 0.32 | fastker, 1e-3, 0.004 | 0.82, 1.25 |
| w1-kd10 | 796 | 1e-6 | 4x18 | 0.002 | 0.002 | 0.003 | 0.003, 0.002, 0.005 | default, 1e-4 (1.3), 0.005, 0.000 | 0.68, 1.04 | 0.34 | fastker, 1e-3, 0.004 | 0.85, 1.30 |
| w1-kd10 | 796 | 1e-8 | 4x18 | 0.002 | 0.002 | 0.003 | 0.003, 0.002, 0.005 | default, 1e-4 (1.3), 0.005, 0.000 | 0.64, 1.02 | 0.32 | fastker, 1e-3, 0.004 | 0.80, 1.28 |
| w1-kd30 | 7162 | 1e-3 | 4x18 | 0.014 | 0.004 | 0.018 | 0.015, 0.004, 0.020 | fastker, 3e-3 (3.8), 0.020, 0.002 | 0.88, 0.99 | 0.20 | fastker, 3e-3, 0.020 | 0.88, 0.99 |
| w1-kd30 | 7162 | 1e-6 | 4x18 | 0.012 | 0.006 | 0.018 | 0.016, 0.008, 0.024 | default, 3e-4 (10.7), 0.030, 0.001 | 0.60, 0.79 | 0.21 | default, 3e-4, 0.030 | 0.60, 0.79 |
| w1-kd30 | 7162 | 1e-8 | 4x18 | 0.011 | 0.011 | 0.022 | 0.015, 0.011, 0.026 | default, 1e-9 (1.2), 0.069, 0.002 | 0.32, 0.37 | 0.16 | fastker, 1e-4, 0.211 | 0.10, 0.12 |
| w1-kd100 | 79578 | 1e-3 | 4x18 | 0.035 | 0.031 | 0.066 | 0.066, 0.034, 0.101 | fastker, 3e-3 (7.9), 0.282, 0.006 | 0.23, 0.36 | 0.11 | default, 3e-4, 0.327 | 0.20, 0.31 |
| w1-kd100 | 79578 | 1e-6 | 4x18 | 0.051 | 0.039 | 0.090 | 0.077, 0.043, 0.120 | default, 1e-5 (2.2), 0.378, 0.007 | 0.24, 0.32 | 0.11 | fastker, 3e-7, 0.591 | 0.15, 0.20 |
| w1-kd100 | 79578 | 1e-8 | 2x36 | 0.053 | 0.108 | 0.161 | 0.067, 0.109, 0.176 | fastker, 3e-7 (4.7), 0.591, 0.005 | 0.27, 0.30 | 0.18 | fastker, 3e-7, 0.591 | 0.27, 0.30 |
| w1-kd250 | 497360 | 1e-3 | 4x18 | 0.222 | 0.164 | 0.386 | 0.560, 0.176, 0.736 | fastker, 3e-3 (2.3), 2.210, 0.029 | 0.17, 0.33 | 0.08 | fastker, 3e-3, 2.210 | 0.17, 0.33 |
| w1-kd250 | 497360 | 1e-6 | 4x18 | 0.256 | 0.260 | 0.516 | 0.589, 0.263, 0.852 | default, 1e-5 (3.4), 3.674, 0.031 | 0.14, 0.23 | 0.07 | default, 3e-4, 3.260 | 0.16, 0.26 |
| w1-kd250 | 497360 | 1e-8 | 4x18 | 0.206 | 0.588 | 0.794 | 0.338, 0.573, 0.911 | fastker, 3e-7 (24.2), 4.125, 0.023 | 0.19, 0.22 | 0.14 | fastker, 3e-7, 4.125 | 0.19, 0.22 |
| w1-kd500 | 1989437 | 1e-3 | 4x18 | 0.986 | 0.708 | 1.694 | 2.714, 0.728, 3.442 | fastker, 3e-3 (2.9), 14.522, 0.125 | 0.12, 0.24 | 0.05 | fastker, 3e-3, 14.522 | 0.12, 0.24 |
| w1-kd500 | 1989437 | 1e-6 | 4x18 | 1.085 | 1.156 | 2.242 | 2.835, 1.165, 4.001 | default, 3e-4 (5.3), 19.085, 0.131 | 0.12, 0.21 | 0.06 | default, 3e-4, 19.085 | 0.12, 0.21 |
| w1-kd500 | 1989437 | 1e-8 | 4x18 | 0.846 | 2.385 | 3.231 | 1.539, 2.297, 3.836 | fastker, 1e-8 (23.3), 25.727, 0.107 | 0.13, 0.15 | 0.09 | fastker, 3e-7, 24.330 | 0.13, 0.16 |
| w1-kd1000 | 7957748 | 1e-3 | 4x18 | 5.047 | 3.667 | 8.715 | 13.446, 3.693, 17.140 | fastker, 3e-3 (3.2), 109.704, 0.835 | 0.08, 0.16 | 0.03 | fastker, 3e-3, 109.704 | 0.08, 0.16 |
| w1-kd1000 | 7957748 | 1e-6 | 4x18 | 5.528 | 5.678 | 11.206 | 14.011, 5.610, 19.621 | fastker, 3e-4 (5.1), 136.648, 0.732 | 0.08, 0.14 | 0.04 | fastker, 1e-4, 139.757 | 0.08, 0.14 |
| w1-kd1000 | 7957748 | 1e-8 | 4x18 | 4.425 | 10.124 | 14.549 | 7.601, 9.867, 17.469 | fastker, 3e-5 (33.0), 143.756, 0.765 | 0.10, 0.12 | 0.07 | fastker, 3e-7, 165.372 | 0.09, 0.11 |
| w2-kd10 | 1276 | 1e-3 | 4x18 | 0.002 | 0.003 | 0.005 | 0.004, 0.005, 0.008 | default, 1e-3 (1.8), 0.028, 0.002 | 0.19, 0.30 | 0.11 | fastker, 1e-3, 0.012 | 0.44, 0.69 |
| w2-kd10 | 1276 | 1e-6 | 4x18 | 0.009 | 0.005 | 0.014 | 0.012, 0.007, 0.019 | default, 3e-4 (1.7), 0.022, 0.005 | 0.63, 0.88 | 0.32 | fastker, 1e-3, 0.012 | 1.16, 1.62 |
| w2-kd10 | 1276 | 1e-8 | 4x18 | 0.002 | 0.004 | 0.006 | 0.003, 0.005, 0.008 | fastker, 1e-3 (1.7), 0.012, 0.000 | 0.47, 0.63 | 0.32 | fastker, 1e-3, 0.012 | 0.47, 0.63 |
| w2-kd30 | 11480 | 1e-3 | 4x18 | 0.014 | 0.006 | 0.020 | 0.017, 0.007, 0.025 | fastker, 1e-3 (5.7), 0.066, 0.014 | 0.31, 0.37 | 0.11 | fastker, 3e-3, 0.058 | 0.35, 0.42 |
| w2-kd30 | 11480 | 1e-6 | 4x18 | 0.017 | 0.009 | 0.027 | 0.024, 0.011, 0.035 | default, 3e-4 (5.8), 0.095, 0.007 | 0.28, 0.37 | 0.11 | default, 3e-4, 0.095 | 0.28, 0.37 |
| w2-kd30 | 11480 | 1e-8 | 2x36 | 0.018 | 0.041 | 0.059 | 0.021, 0.042, 0.063 | default, 1e-9 (1.8), 0.188, 0.006 | 0.31, 0.34 | 0.23 | default, >1e-9, 0.188 | 0.31, 0.34 |
| w2-kd100 | 127554 | 1e-3 | 4x18 | 0.045 | 0.044 | 0.089 | 0.100, 0.050, 0.151 | fastker, 3e-3 (5.4), 0.473, 0.009 | 0.19, 0.32 | 0.09 | default, 3e-4, 0.668 | 0.13, 0.23 |
| w2-kd100 | 127554 | 1e-6 | 4x18 | 0.061 | 0.064 | 0.124 | 0.113, 0.070, 0.183 | default, 3e-4 (1.7), 0.668, 0.010 | 0.19, 0.27 | 0.10 | fastker, 3e-7, 1.101 | 0.11, 0.17 |
| w2-kd100 | 127554 | 1e-8 | 2x36 | 0.070 | 0.173 | 0.243 | 0.096, 0.172, 0.268 | fastker, 3e-7 (4.9), 1.101, 0.008 | 0.22, 0.24 | 0.16 | fastker, 3e-7, 1.101 | 0.22, 0.24 |
| w2-kd250 | 797211 | 1e-3 | 4x18 | 0.283 | 0.253 | 0.536 | 0.779, 0.267, 1.047 | fastker, 3e-3 (3.0), 3.299, 0.098 | 0.16, 0.32 | 0.08 | fastker, 3e-3, 3.299 | 0.16, 0.32 |
| w2-kd250 | 797211 | 1e-6 | 4x18 | 0.325 | 0.395 | 0.720 | 0.819, 0.404, 1.223 | default, 3e-4 (5.1), 4.792, 0.107 | 0.15, 0.26 | 0.08 | default, 3e-4, 4.792 | 0.15, 0.26 |
| w2-kd250 | 797211 | 1e-8 | 4x18 | 0.282 | 0.985 | 1.268 | 0.544, 0.964, 1.507 | fastker, 1e-7 (24.6), 7.216, 0.078 | 0.18, 0.21 | 0.14 | fastker, 3e-7, 7.082 | 0.18, 0.21 |
| w2-kd500 | 3188840 | 1e-3 | 4x18 | 1.330 | 1.056 | 2.386 | 3.892, 1.099, 4.990 | fastker, 3e-3 (3.3), 18.248, 0.353 | 0.13, 0.27 | 0.06 | fastker, 3e-3, 18.248 | 0.13, 0.27 |
| w2-kd500 | 3188840 | 1e-6 | 4x18 | 1.426 | 1.675 | 3.101 | 3.941, 1.669, 5.609 | default, 3e-4 (6.5), 25.285, 0.334 | 0.12, 0.22 | 0.07 | default, 3e-4, 25.285 | 0.12, 0.22 |
| w2-kd500 | 3188840 | 1e-8 | 4x18 | 1.152 | 3.825 | 4.978 | 2.615, 3.712, 6.327 | fastker, 3e-7 (18.4), 33.144, 0.299 | 0.15, 0.19 | 0.12 | fastker, 3e-7, 33.144 | 0.15, 0.19 |
| w3-kd10 | 4032 | 1e-3 | 2x36 | 0.009 | 0.003 | 0.012 | 0.009, 0.004, 0.013 | default, 3e-3 (2.9), 0.037, 0.016 | 0.32, 0.36 | 0.15 | default, 3e-3, 0.037 | 0.32, 0.36 |
| w3-kd10 | 4032 | 1e-6 | 4x18 | 0.007 | 0.004 | 0.011 | 0.011, 0.004, 0.015 | default, 3e-4 (1.1), 0.028, 0.016 | 0.40, 0.54 | 0.31 | default, 1e-4, 0.041 | 0.27, 0.37 |
| w3-kd10 | 4032 | 1e-8 | 4x18 | 0.010 | 0.007 | 0.017 | 0.013, 0.007, 0.020 | default, 3e-9 (7.9), 0.070, 0.018 | 0.25, 0.29 | 0.14 | fastker, 1e-4, 0.080 | 0.21, 0.26 |
| w3-kd30 | 108849 | 1e-3 | 8x9 | 0.041 | 0.035 | 0.076 | 0.128, 0.040, 0.168 | default, 1e-3 (3.8), 0.193, 0.045 | 0.39, 0.87 | 0.24 | default, 3e-3, 0.185 | 0.41, 0.91 |
| w3-kd30 | 108849 | 1e-6 | 8x9 | 0.036 | 0.068 | 0.104 | 0.130, 0.062, 0.192 | default, 3e-4 (6.8), 0.381, 0.036 | 0.27, 0.50 | 0.20 | default, 3e-4, 0.381 | 0.27, 0.50 |
| w3-kd30 | 108849 | 1e-8 | 4x18 | 0.029 | 0.145 | 0.174 | 0.054, 0.142, 0.195 | fastker, 1e-9 (2.1), 0.511, 0.036 | 0.34, 0.38 | 0.31 | fastker, >1e-9, 0.511 | 0.34, 0.38 |
| w3-kd100 | 4031442 | 1e-3 | 8x9 | 0.559 | 2.403 | 2.962 | 4.321, 2.366, 6.687 | fastker, 3e-3 (6.4), 8.167, 0.653 | 0.36, 0.82 | 0.32 | default, 3e-4, 11.084 | 0.27, 0.60 |
| w3-kd100 | 4031442 | 1e-6 | 8x9 | 0.560 | 3.098 | 3.659 | 4.321, 2.970, 7.291 | default, 3e-4 (2.1), 11.084, 0.657 | 0.33, 0.66 | 0.30 | default, 3e-7, 13.204 | 0.28, 0.55 |
| w3-kd100 | 4031442 | 1e-8 | 8x9 | 0.445 | 4.432 | 4.876 | 2.600, 4.156, 6.756 | default, 3e-7 (8.4), 13.204, 0.660 | 0.37, 0.51 | 0.35 | default, 3e-7, 13.204 | 0.37, 0.51 |
| w3-n1000000-kd200 | 1000000 | 1e-3 | 8x9 | 0.241 | 0.529 | 0.771 | 1.277, 0.474, 1.751 | default, 3e-3 (5.1), 3.579, 0.170 | 0.22, 0.49 | 0.16 | fastker, 3e-4, 4.644 | 0.17, 0.38 |
| w3-n1000000-kd200 | 1000000 | 1e-6 | 8x9 | 0.287 | 0.772 | 1.060 | 1.277, 0.684, 1.961 | fastker, 1e-4 (2.6), 4.688, 0.136 | 0.23, 0.42 | 0.17 | fastker, 3e-7, 5.256 | 0.20, 0.37 |
| w3-n1000000-kd200 | 1000000 | 1e-8 | 4x18 | 0.247 | 1.557 | 1.803 | 0.571, 1.454, 2.025 | fastker, 3e-7 (7.7), 5.256, 0.144 | 0.34, 0.39 | 0.30 | fastker, 3e-7, 5.256 | 0.34, 0.39 |
| w3-n1000000-kd500 | 1000000 | 1e-3 | 4x18 | 0.637 | 1.661 | 2.298 | 1.433, 1.439, 2.872 | default, 3e-3 (1.8), 23.491, 0.165 | 0.10, 0.12 | 0.07 | default, 3e-3, 23.491 | 0.10, 0.12 |
| w3-n1000000-kd500 | 1000000 | 1e-6 | 4x18 | 1.007 | 2.449 | 3.456 | 1.496, 2.131, 3.628 | fastker, 3e-5 (2.5), 28.332, 0.135 | 0.12, 0.13 | 0.09 | fastker, 3e-4, 25.869 | 0.13, 0.14 |
| w3-n1000000-kd500 | 1000000 | 1e-8 | 4x18 | 1.102 | 4.177 | 5.280 | 2.747, 3.708, 6.454 | fastker, 3e-7 (34.1), 33.526, 0.136 | 0.16, 0.19 | 0.13 | fastker, 3e-7, 33.526 | 0.16, 0.19 |
| w1-kd500-planewave | 1989437 | 1e-3 | 4x18 | 1.533 | 0.898 | 2.431 | 4.373, 0.886, 5.260 | default, 3e-3 (5.0), 18.703, 0.237 | 0.13, 0.28 | 0.05 | default, 3e-3, 18.703 | 0.13, 0.28 |
| w1-kd500-planewave | 1989437 | 1e-6 | 4x18 | 1.621 | 1.440 | 3.061 | 4.664, 1.389, 6.053 | default, 3e-5 (4.5), 20.859, 0.181 | 0.15, 0.29 | 0.07 | default, 3e-4, 19.265 | 0.16, 0.31 |
| w1-kd500-planewave | 1989437 | 1e-8 | 4x18 | 1.302 | 2.871 | 4.174 | 2.532, 2.792, 5.324 | default, 1e-7 (27.5), 26.232, 0.161 | 0.16, 0.20 | 0.11 | default, 3e-7, 25.112 | 0.17, 0.21 |
| w1-kd1000-planewave | 7957748 | 1e-3 | 4x18 | 5.836 | 4.284 | 10.121 | 16.816, 4.082, 20.898 | fastker, 1e-3 (7.1), 114.879, 0.922 | 0.09, 0.18 | 0.04 | fastker, 1e-3, 114.879 | 0.09, 0.18 |
| w1-kd1000-planewave | 7957748 | 1e-6 | 4x18 | 6.712 | 6.462 | 13.173 | 16.973, 6.376, 23.349 | fastker, 3e-4 (8.5), 136.496, 0.787 | 0.10, 0.17 | 0.05 | fastker, 1e-4, 139.792 | 0.09, 0.17 |
| w1-kd1000-planewave | 7957748 | 1e-8 | 4x18 | 4.351 | 10.032 | 14.383 | 9.585, 11.396, 20.981 | fastker, 3e-7 (36.2), 165.868, 0.780 | 0.09, 0.13 | 0.06 | fastker, 3e-7, 165.868 | 0.09, 0.13 |
| w1-kd10-constant | 796 | 1e-3 | 4x18 | 0.002 | 0.002 | 0.003 | 0.004, 0.003, 0.007 | default, 1e-4 (1.2), 0.006, 0.007 | 0.57, 1.18 | – | default, 1e-4, 0.006 | 0.57, 1.18 |
| w1-kd10-constant | 796 | 1e-6 | 4x18 | 0.002 | 0.002 | 0.003 | 0.003, 0.003, 0.006 | default, 1e-4 (1.2), 0.006, 0.007 | 0.57, 1.00 | – | default, 1e-4, 0.006 | 0.57, 1.00 |
| w1-kd10-constant | 796 | 1e-8 | 4x18 | 0.002 | 0.002 | 0.004 | 0.004, 0.004, 0.009 | default, 1e-4 (1.2), 0.006, 0.007 | 0.60, 1.47 | – | default, 1e-4, 0.006 | 0.60, 1.47 |
| w1-kd30-constant | 7162 | 1e-3 | 4x18 | 0.014 | 0.004 | 0.017 | 0.016, 0.006, 0.023 | default, 3e-3 (2.5), 0.022, 0.002 | 0.79, 1.02 | 0.19 | default, 3e-4, 0.028 | 0.62, 0.80 |
| w1-kd30-constant | 7162 | 1e-6 | 2x36 | 0.013 | 0.007 | 0.019 | 0.017, 0.008, 0.025 | default, 3e-4 (9.8), 0.028, 0.001 | 0.69, 0.88 | 0.26 | default, 3e-4, 0.028 | 0.69, 0.88 |
| w1-kd30-constant | 7162 | 1e-8 | 2x36 | 0.011 | 0.011 | 0.022 | 0.014, 0.011, 0.025 | default, 3e-8 (1.0), 0.055, 0.001 | 0.40, 0.46 | 0.20 | default, >1e-9, 0.067 | 0.33, 0.38 |
| w1-kd100-constant | 79578 | 1e-3 | 4x18 | 0.034 | 0.040 | 0.074 | 0.065, 0.044, 0.109 | default, 3e-3 (10.0), 0.413, 0.047 | 0.18, 0.26 | 0.11 | default, 3e-5, 0.379 | 0.20, 0.29 |
| w1-kd100-constant | 79578 | 1e-6 | 4x18 | 0.049 | 0.040 | 0.089 | 0.092, 0.059, 0.152 | default, 3e-4 (3.1), 0.368, 0.034 | 0.24, 0.41 | 0.12 | default, 3e-7, 0.658 | 0.13, 0.23 |
| w1-kd100-constant | 79578 | 1e-8 | 2x36 | 0.050 | 0.113 | 0.163 | 0.065, 0.135, 0.199 | default, 3e-7 (11.5), 0.658, 0.014 | 0.25, 0.30 | 0.17 | default, 3e-7, 0.658 | 0.25, 0.30 |
| w2-kd10-constant | 1276 | 1e-3 | 8x9 | 0.005 | 0.005 | 0.009 | 0.006, 0.008, 0.014 | default, 3e-3 (3.2), 0.030, 0.003 | 0.31, 0.47 | 0.17 | default, 1e-7, 0.017 | 0.55, 0.82 |
| w2-kd10-constant | 1276 | 1e-6 | 4x18 | 0.009 | 0.006 | 0.015 | 0.013, 0.008, 0.020 | default, 3e-4 (1.7), 0.023, 0.015 | 0.64, 0.88 | 0.72 | default, 1e-7, 0.017 | 0.87, 1.19 |
| w2-kd10-constant | 1276 | 1e-8 | 8x9 | 0.002 | 0.004 | 0.006 | 0.003, 0.005, 0.008 | default, 1e-7 (1.4), 0.017, 0.000 | 0.35, 0.46 | 0.22 | default, 1e-7, 0.017 | 0.35, 0.46 |
| w2-kd30-constant | 11480 | 1e-3 | 2x36 | 0.015 | 0.006 | 0.022 | 0.025, 0.007, 0.033 | default, 3e-3 (1.8), 0.072, 0.009 | 0.31, 0.45 | 0.10 | default, 3e-3, 0.072 | 0.31, 0.45 |
| w2-kd30-constant | 11480 | 1e-6 | 4x18 | 0.018 | 0.009 | 0.027 | 0.024, 0.013, 0.036 | default, 3e-4 (4.0), 0.096, 0.004 | 0.28, 0.38 | 0.10 | default, 3e-4, 0.096 | 0.28, 0.38 |
| w2-kd30-constant | 11480 | 1e-8 | 4x18 | 0.019 | 0.043 | 0.062 | 0.021, 0.073, 0.094 | default, 3e-9 (1.0), 0.177, 0.004 | 0.35, 0.53 | 0.25 | default, >1e-9, 0.197 | 0.31, 0.48 |
| w2-kd100-constant | 127554 | 1e-3 | 8x9 | 0.048 | 0.043 | 0.090 | 0.139, 0.052, 0.191 | default, 3e-3 (10.5), 0.568, 0.038 | 0.16, 0.34 | 0.08 | default, 3e-4, 0.745 | 0.12, 0.26 |
| w2-kd100-constant | 127554 | 1e-6 | 4x18 | 0.061 | 0.067 | 0.128 | 0.111, 0.085, 0.196 | default, 1e-4 (2.3), 0.777, 0.038 | 0.16, 0.25 | 0.09 | default, 3e-7, 1.166 | 0.11, 0.17 |
| w2-kd100-constant | 127554 | 1e-8 | 4x18 | 0.061 | 0.215 | 0.276 | 0.118, 0.173, 0.292 | default, 3e-8 (9.2), 1.286, 0.027 | 0.21, 0.23 | 0.17 | default, 3e-7, 1.166 | 0.24, 0.25 |
| w3-kd10-constant | 4032 | 1e-3 | 2x36 | 0.013 | 0.003 | 0.016 | 0.015, 0.006, 0.021 | default, 3e-3 (2.3), 0.018, 0.016 | 0.89, 1.16 | 1.60 | default, 3e-4, 0.022 | 0.73, 0.95 |
| w3-kd10-constant | 4032 | 1e-6 | 2x36 | 0.011 | 0.005 | 0.017 | 0.018, 0.007, 0.025 | default, 3e-4 (1.0), 0.022, 0.002 | 0.75, 1.12 | 0.26 | default, 3e-5, 0.025 | 0.66, 0.99 |
| w3-kd10-constant | 4032 | 1e-8 | 1x72 | 0.013 | 0.009 | 0.023 | 0.013, 0.009, 0.023 | default, 3e-9 (2.4), 0.054, 0.017 | 0.42, 0.42 | 0.25 | default, >3e-9, 0.054 | 0.42, 0.42 |
| w3-kd30-constant | 108849 | 1e-3 | 8x9 | 0.042 | 0.037 | 0.079 | 0.154, 0.048, 0.201 | default, 3e-3 (1.5), 0.184, 0.038 | 0.43, 1.09 | 0.25 | default, 3e-4, 0.367 | 0.22, 0.55 |
| w3-kd30-constant | 108849 | 1e-6 | 8x9 | 0.035 | 0.087 | 0.123 | 0.134, 0.077, 0.210 | default, 3e-4 (2.8), 0.367, 0.041 | 0.33, 0.57 | 0.27 | default, 3e-4, 0.367 | 0.33, 0.57 |
| w3-kd30-constant | 108849 | 1e-8 | 2x36 | 0.037 | 0.166 | 0.204 | 0.053, 0.173, 0.226 | default, 3e-8 (1.0), 0.607, 0.026 | 0.34, 0.37 | 0.29 | default, >3e-8, 0.607 | 0.34, 0.37 |

### 3.4 Where FMM3D is faster

On locust ours' best split is faster on every row (the largest ratio 0.89, W3 at kD = 10 with the
constant density, 10⁻³); at 1 × 72 FMM3D is faster on eight rows of N ≤ 1.1·10⁵ (up to 1.47×, W1
at kD = 10 with the constant density, a call of 6 ms).

On Kathleen FMM3D is faster at the best split on ten rows, all of N ≤ 1.2·10⁴: at 10⁻³ W1 and W2
at kD = 30 (4.37× and 1.98×), W2 and W3 at kD = 10, W1 and W2 at kD = 30 with the constant density,
W3 at kD = 10 with it (1.14–2.68×); at 10⁻⁶ and 10⁻⁸ W2 at kD = 10 (1.26–1.32×) and W3 at kD = 10 at
10⁻⁸ (1.41×). FMM3D's calls are 5–104 ms. (At 1 × 40 also W3 at kD = 30 at 10⁻³, N = 1.1·10⁵,
1.16×.) Our evaluation is faster on every one of them (0.15–0.74× the call less its tree); our build
is what loses. On W1 at kD = 30, 10⁻³ (Kathleen, 1 × 40) the build is 81 ms, 73 ms of it "matrices"
(Batched's, at p = 11 on leaf level 2, 2,448 V pairs), against 20 ms at 10⁻⁶ (the diagonal form
there). Batched's matrices are 73–101 ms of the build on W1 and W2 at kD = 30 (both densities) and
15–16 ms on W3 at kD = 10 (both densities), all at 10⁻³; on W2 at kD = 10 (every ε) and W3 at kD = 10
at 10⁻⁸ no matrices are built and the build is 4–20 ms of fixed steps. On locust the same W1 case
builds in 15 ms: its costs (NEON, not AVX-512) choose one level deeper at every split (the diagonal
form at level 2, Batched at p = 7 on level 3, 6 ms of tables), where Kathleen's choose the shallow
tree at every split too. The time rule prices the evaluation, not the tables' build; for one solve
on a small tree the build is the cost. Question 4 proposes it as CH4.7 for the six rows it explains.

Under the one-sided match the same rows lose, and also W1 at kD = 10 on Kathleen (N = 796, our
error 10⁻¹⁵, so only FMM3D's direct sum is as accurate: 2.55–2.60×, calls of 2 ms) and W2 at kD = 10 at 10⁻⁶ on
locust (1.16×, FMM3D's direct sum against our 5·10⁻⁸).

## 4. Setup apart

FMM3D's tree, timed once apart, is 0.2–14 % of its call on the trees of N ≥ 10⁵ on Kathleen and
0.4–23 % on locust (W1 at kD = 1000: 0.73–0.92 s of 110–166 s on locust; W3 at kD = 100: 0.65 s of
8–13 s). Our evaluation alone over FMM3D's
call less its tree (the repeated-solve measure), N ≥ 10⁵, geometric means: Kathleen 0.10, 0.12,
0.18 and locust 0.09, 0.10, 0.15 at 10⁻³, 10⁻⁶, 10⁻⁸ (ranges in §3.1). Over all 26 cases the range
on locust reaches 1.60 on one row (W3 at kD = 10 with the constant density, 10⁻³), where the tree,
timed apart, took 16 ms of an 18 ms call against 2 ms at the case's other eps: one noisy timing of a
call of milliseconds. Our build is 8–63 % of
build plus evaluation at the best split on the trees of N ≥ 10⁶ (the most at 10⁻³, the least on
W3 at kD = 100 at 10⁻⁸); the first evaluation beside the second is in each archived output's Build
table.

## 5. One rank per node: the build by part (CH4.0k)

The cases of N ≥ 10⁵ at 1 × N and at the best split. "shared share": octree, routes and plan over
the build, the steps every rank repeats over the whole tree (CH4.0h §3). Asked at CH4.5's start, the
user kept CH4.0k off the list ("Keep off; CH4.5 shows it"); this is what it shows. At 1 × N the shared
steps are 32–88 % of the build on Kathleen and 58–91 % on locust, at the best split 24–66 % and
27–83 %; the operator's tables and the smaller steps are the rest. The build at 1 × N is 1.1–7.7×
the best split's (W3 at kD = 100: 4.0 s against 0.8 s on Kathleen, 4.3 s against 0.56 s on locust),
the evaluation 0.80–1.32× the best split's. The H4 gate's one-rank-per-node figures are the 1 × N columns of §3.2–§3.3 and this table.

Kathleen:

| Case | ε | split | build s | octree | routes | plan | shared share | operator | evaluate s |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| w1-kd250 | 1e-3 | 1 × N 1x40 | 0.529 | 0.075 | 0.120 | 0.133 | 0.62 | 0.156 | 0.468 |
| w1-kd250 | 1e-3 | best 2x20 | 0.379 | 0.046 | 0.070 | 0.066 | 0.48 | 0.171 | 0.381 |
| w1-kd250 | 1e-6 | 1 × N 1x40 | 0.594 | 0.074 | 0.121 | 0.134 | 0.56 | 0.220 | 0.642 |
| w1-kd250 | 1e-6 | best 2x20 | 0.460 | 0.048 | 0.070 | 0.066 | 0.40 | 0.250 | 0.555 |
| w1-kd250 | 1e-8 | 1 × N 1x40 | 0.666 | 0.075 | 0.120 | 0.150 | 0.52 | 0.278 | 1.445 |
| w1-kd250 | 1e-8 | best 2x20 | 0.519 | 0.048 | 0.069 | 0.074 | 0.37 | 0.300 | 1.250 |
| w1-kd500 | 1e-3 | 1 × N 1x40 | 2.447 | 0.341 | 0.682 | 0.693 | 0.70 | 0.558 | 2.032 |
| w1-kd500 | 1e-3 | best 2x20 | 1.773 | 0.194 | 0.392 | 0.339 | 0.52 | 0.758 | 1.762 |
| w1-kd500 | 1e-6 | 1 × N 1x40 | 2.716 | 0.342 | 0.693 | 0.697 | 0.64 | 0.810 | 2.901 |
| w1-kd500 | 1e-6 | best 2x20 | 1.997 | 0.195 | 0.399 | 0.346 | 0.47 | 0.960 | 2.538 |
| w1-kd500 | 1e-8 | 1 × N 1x40 | 2.850 | 0.340 | 0.685 | 0.760 | 0.63 | 0.898 | 5.642 |
| w1-kd500 | 1e-8 | best 2x20 | 2.184 | 0.194 | 0.396 | 0.374 | 0.44 | 1.130 | 5.320 |
| w1-kd1000 | 1e-3 | 1 × N 1x40 | 12.692 | 1.461 | 4.044 | 3.252 | 0.69 | 3.243 | 10.441 |
| w1-kd1000 | 1e-3 | best 2x20 | 8.787 | 0.882 | 2.099 | 1.569 | 0.52 | 3.905 | 8.994 |
| w1-kd1000 | 1e-6 | 1 × N 1x40 | 13.724 | 1.485 | 4.042 | 3.277 | 0.64 | 4.223 | 14.166 |
| w1-kd1000 | 1e-6 | best 2x20 | 9.897 | 0.886 | 2.099 | 1.587 | 0.46 | 4.992 | 12.608 |
| w1-kd1000 | 1e-8 | 1 × N 1x40 | 14.267 | 1.472 | 4.053 | 3.524 | 0.63 | 4.588 | 25.143 |
| w1-kd1000 | 1e-8 | best 2x20 | 10.567 | 0.881 | 2.109 | 1.692 | 0.44 | 5.577 | 23.900 |
| w2-kd100 | 1e-3 | 1 × N 1x40 | 0.163 | 0.022 | 0.025 | 0.041 | 0.54 | 0.059 | 0.108 |
| w2-kd100 | 1e-3 | best 2x20 | 0.129 | 0.016 | 0.016 | 0.021 | 0.41 | 0.065 | 0.083 |
| w2-kd100 | 1e-6 | 1 × N 1x40 | 0.169 | 0.023 | 0.025 | 0.043 | 0.54 | 0.062 | 0.161 |
| w2-kd100 | 1e-6 | best 2x20 | 0.138 | 0.016 | 0.016 | 0.022 | 0.39 | 0.072 | 0.138 |
| w2-kd100 | 1e-8 | 1 × N 1x40 | 0.173 | 0.023 | 0.025 | 0.047 | 0.55 | 0.061 | 0.433 |
| w2-kd100 | 1e-8 | best 2x20 | 0.141 | 0.016 | 0.016 | 0.023 | 0.39 | 0.073 | 0.399 |
| w2-kd250 | 1e-3 | 1 × N 1x40 | 0.849 | 0.128 | 0.221 | 0.262 | 0.72 | 0.167 | 0.627 |
| w2-kd250 | 1e-3 | best 2x20 | 0.577 | 0.081 | 0.126 | 0.129 | 0.58 | 0.201 | 0.502 |
| w2-kd250 | 1e-6 | 1 × N 1x40 | 0.924 | 0.128 | 0.221 | 0.265 | 0.66 | 0.237 | 0.903 |
| w2-kd250 | 1e-6 | best 2x20 | 0.642 | 0.080 | 0.124 | 0.132 | 0.52 | 0.267 | 0.838 |
| w2-kd250 | 1e-8 | 1 × N 1x40 | 0.975 | 0.128 | 0.219 | 0.284 | 0.65 | 0.274 | 2.121 |
| w2-kd250 | 1e-8 | best 2x20 | 0.705 | 0.081 | 0.124 | 0.142 | 0.49 | 0.317 | 2.003 |
| w2-kd500 | 1e-3 | 1 × N 1x40 | 4.040 | 0.586 | 1.245 | 1.335 | 0.78 | 0.598 | 2.507 |
| w2-kd500 | 1e-3 | best 2x20 | 2.622 | 0.361 | 0.703 | 0.655 | 0.66 | 0.762 | 2.224 |
| w2-kd500 | 1e-6 | 1 × N 1x40 | 4.257 | 0.586 | 1.211 | 1.351 | 0.74 | 0.834 | 3.991 |
| w2-kd500 | 1e-6 | best 2x20 | 2.838 | 0.361 | 0.708 | 0.663 | 0.61 | 0.959 | 3.631 |
| w2-kd500 | 1e-8 | 1 × N 1x40 | 4.426 | 0.585 | 1.217 | 1.447 | 0.73 | 0.911 | 8.304 |
| w2-kd500 | 1e-8 | best 2x20 | 3.045 | 0.361 | 0.714 | 0.708 | 0.59 | 1.123 | 7.968 |
| w3-kd30 | 1e-3 | 1 × N 1x40 | 0.137 | 0.016 | 0.019 | 0.009 | 0.32 | 0.080 | 0.104 |
| w3-kd30 | 1e-3 | best 2x20 | 0.122 | 0.011 | 0.013 | 0.005 | 0.24 | 0.083 | 0.079 |
| w3-kd30 | 1e-6 | 1 × N 1x40 | 0.081 | 0.017 | 0.020 | 0.009 | 0.56 | 0.022 | 0.136 |
| w3-kd30 | 1e-6 | best 4x10 | 0.055 | 0.007 | 0.007 | 0.003 | 0.31 | 0.032 | 0.115 |
| w3-kd30 | 1e-8 | 1 × N 1x40 | 0.087 | 0.017 | 0.020 | 0.009 | 0.52 | 0.025 | 0.292 |
| w3-kd30 | 1e-8 | best 2x20 | 0.066 | 0.011 | 0.014 | 0.005 | 0.46 | 0.024 | 0.243 |
| w3-kd100 | 1e-3 | 1 × N 1x40 | 3.995 | 0.598 | 1.825 | 1.046 | 0.87 | 0.086 | 2.545 |
| w3-kd100 | 1e-3 | best 8x5 | 0.803 | 0.113 | 0.247 | 0.134 | 0.62 | 0.244 | 2.177 |
| w3-kd100 | 1e-6 | 1 × N 1x40 | 3.976 | 0.599 | 1.796 | 1.091 | 0.88 | 0.059 | 4.027 |
| w3-kd100 | 1e-6 | best 8x5 | 0.756 | 0.113 | 0.245 | 0.137 | 0.65 | 0.196 | 3.784 |
| w3-kd100 | 1e-8 | 1 × N 1x40 | 3.891 | 0.597 | 1.740 | 1.034 | 0.87 | 0.055 | 9.410 |
| w3-kd100 | 1e-8 | best 8x5 | 0.815 | 0.112 | 0.247 | 0.136 | 0.61 | 0.238 | 9.738 |
| w3-n1000000-kd200 | 1e-3 | 1 × N 1x40 | 0.861 | 0.140 | 0.299 | 0.184 | 0.72 | 0.149 | 1.428 |
| w3-n1000000-kd200 | 1e-3 | best 4x10 | 0.472 | 0.049 | 0.081 | 0.050 | 0.38 | 0.265 | 1.304 |
| w3-n1000000-kd200 | 1e-6 | 1 × N 1x40 | 0.902 | 0.141 | 0.303 | 0.201 | 0.71 | 0.167 | 1.808 |
| w3-n1000000-kd200 | 1e-6 | best 4x10 | 0.515 | 0.050 | 0.080 | 0.054 | 0.36 | 0.304 | 1.704 |
| w3-n1000000-kd200 | 1e-8 | 1 × N 1x40 | 0.928 | 0.141 | 0.308 | 0.208 | 0.71 | 0.177 | 2.759 |
| w3-n1000000-kd200 | 1e-8 | best 2x20 | 0.620 | 0.086 | 0.165 | 0.105 | 0.57 | 0.209 | 2.617 |
| w3-n1000000-kd500 | 1e-3 | 1 × N 1x40 | 2.949 | 0.178 | 0.373 | 1.515 | 0.70 | 0.781 | 3.745 |
| w3-n1000000-kd500 | 1e-3 | best 2x20 | 1.983 | 0.111 | 0.204 | 0.765 | 0.54 | 0.840 | 3.034 |
| w3-n1000000-kd500 | 1e-6 | 1 × N 1x40 | 3.071 | 0.179 | 0.381 | 1.526 | 0.68 | 0.889 | 5.147 |
| w3-n1000000-kd500 | 1e-6 | best 2x20 | 2.122 | 0.104 | 0.196 | 0.756 | 0.50 | 1.004 | 4.015 |
| w3-n1000000-kd500 | 1e-8 | 1 × N 1x40 | 3.086 | 0.178 | 0.398 | 1.519 | 0.68 | 0.845 | 6.721 |
| w3-n1000000-kd500 | 1e-8 | best 2x20 | 2.173 | 0.103 | 0.198 | 0.765 | 0.49 | 1.020 | 6.313 |
| w1-kd500-planewave | 1e-3 | 1 × N 1x40 | 2.492 | 0.343 | 0.686 | 0.688 | 0.69 | 0.603 | 2.049 |
| w1-kd500-planewave | 1e-3 | best 2x20 | 1.732 | 0.193 | 0.390 | 0.339 | 0.53 | 0.720 | 1.770 |
| w1-kd500-planewave | 1e-6 | 1 × N 1x40 | 2.681 | 0.340 | 0.670 | 0.698 | 0.64 | 0.799 | 2.828 |
| w1-kd500-planewave | 1e-6 | best 2x20 | 1.979 | 0.194 | 0.392 | 0.344 | 0.47 | 0.954 | 2.517 |
| w1-kd500-planewave | 1e-8 | 1 × N 1x40 | 2.855 | 0.339 | 0.684 | 0.762 | 0.63 | 0.903 | 5.610 |
| w1-kd500-planewave | 1e-8 | best 2x20 | 2.171 | 0.194 | 0.395 | 0.374 | 0.44 | 1.117 | 5.306 |
| w1-kd1000-planewave | 1e-3 | 1 × N 1x40 | 12.551 | 1.470 | 3.939 | 3.256 | 0.69 | 3.202 | 10.157 |
| w1-kd1000-planewave | 1e-3 | best 2x20 | 8.858 | 0.904 | 2.082 | 1.581 | 0.52 | 3.952 | 8.931 |
| w1-kd1000-planewave | 1e-6 | 1 × N 1x40 | 13.702 | 1.461 | 4.043 | 3.300 | 0.64 | 4.244 | 14.236 |
| w1-kd1000-planewave | 1e-6 | best 2x20 | 10.005 | 0.879 | 2.084 | 1.596 | 0.46 | 5.116 | 12.472 |
| w1-kd1000-planewave | 1e-8 | 1 × N 1x40 | 14.299 | 1.458 | 4.009 | 3.524 | 0.63 | 4.664 | 25.888 |
| w1-kd1000-planewave | 1e-8 | best 2x20 | 10.771 | 0.883 | 2.111 | 1.697 | 0.44 | 5.774 | 23.893 |
| w2-kd100-constant | 1e-3 | 1 × N 1x40 | 0.166 | 0.025 | 0.025 | 0.041 | 0.55 | 0.059 | 0.101 |
| w2-kd100-constant | 1e-3 | best 2x20 | 0.130 | 0.016 | 0.017 | 0.021 | 0.42 | 0.064 | 0.088 |
| w2-kd100-constant | 1e-6 | 1 × N 1x40 | 0.168 | 0.024 | 0.025 | 0.042 | 0.54 | 0.062 | 0.167 |
| w2-kd100-constant | 1e-6 | best 2x20 | 0.136 | 0.016 | 0.017 | 0.022 | 0.40 | 0.070 | 0.139 |
| w2-kd100-constant | 1e-8 | 1 × N 1x40 | 0.175 | 0.022 | 0.025 | 0.046 | 0.54 | 0.062 | 0.416 |
| w2-kd100-constant | 1e-8 | best 2x20 | 0.142 | 0.017 | 0.016 | 0.023 | 0.39 | 0.073 | 0.386 |
| w3-kd30-constant | 1e-3 | 1 × N 1x40 | 0.142 | 0.017 | 0.021 | 0.009 | 0.33 | 0.082 | 0.104 |
| w3-kd30-constant | 1e-3 | best 2x20 | 0.124 | 0.011 | 0.015 | 0.005 | 0.25 | 0.084 | 0.081 |
| w3-kd30-constant | 1e-6 | 1 × N 1x40 | 0.083 | 0.017 | 0.021 | 0.009 | 0.56 | 0.023 | 0.142 |
| w3-kd30-constant | 1e-6 | best 4x10 | 0.057 | 0.007 | 0.007 | 0.003 | 0.30 | 0.034 | 0.116 |
| w3-kd30-constant | 1e-8 | 1 × N 1x40 | 0.087 | 0.017 | 0.021 | 0.009 | 0.53 | 0.025 | 0.275 |
| w3-kd30-constant | 1e-8 | best 2x20 | 0.065 | 0.011 | 0.014 | 0.005 | 0.47 | 0.024 | 0.245 |

locust:

| Case | ε | split | build s | octree | routes | plan | shared share | operator | evaluate s |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| w1-kd250 | 1e-3 | 1 × N 1x72 | 0.560 | 0.071 | 0.118 | 0.280 | 0.84 | 0.057 | 0.176 |
| w1-kd250 | 1e-3 | best 4x18 | 0.222 | 0.024 | 0.034 | 0.069 | 0.57 | 0.084 | 0.164 |
| w1-kd250 | 1e-6 | 1 × N 1x72 | 0.589 | 0.071 | 0.120 | 0.283 | 0.81 | 0.080 | 0.263 |
| w1-kd250 | 1e-6 | best 4x18 | 0.256 | 0.023 | 0.034 | 0.070 | 0.49 | 0.117 | 0.260 |
| w1-kd250 | 1e-8 | 1 × N 1x72 | 0.338 | 0.046 | 0.109 | 0.080 | 0.70 | 0.081 | 0.573 |
| w1-kd250 | 1e-8 | best 4x18 | 0.206 | 0.017 | 0.035 | 0.020 | 0.35 | 0.126 | 0.588 |
| w1-kd500 | 1e-3 | 1 × N 1x72 | 2.714 | 0.330 | 0.599 | 1.442 | 0.87 | 0.200 | 0.728 |
| w1-kd500 | 1e-3 | best 4x18 | 0.986 | 0.095 | 0.193 | 0.340 | 0.64 | 0.319 | 0.708 |
| w1-kd500 | 1e-6 | 1 × N 1x72 | 2.835 | 0.332 | 0.605 | 1.458 | 0.84 | 0.296 | 1.165 |
| w1-kd500 | 1e-6 | best 4x18 | 1.085 | 0.094 | 0.184 | 0.342 | 0.57 | 0.423 | 1.156 |
| w1-kd500 | 1e-8 | 1 × N 1x72 | 1.539 | 0.199 | 0.556 | 0.414 | 0.76 | 0.295 | 2.297 |
| w1-kd500 | 1e-8 | best 4x18 | 0.846 | 0.060 | 0.182 | 0.102 | 0.41 | 0.476 | 2.385 |
| w1-kd1000 | 1e-3 | 1 × N 1x72 | 13.446 | 1.695 | 2.957 | 6.855 | 0.86 | 1.299 | 3.693 |
| w1-kd1000 | 1e-3 | best 4x18 | 5.047 | 0.420 | 0.978 | 1.574 | 0.59 | 1.913 | 3.667 |
| w1-kd1000 | 1e-6 | 1 × N 1x72 | 14.011 | 1.719 | 2.975 | 6.938 | 0.83 | 1.750 | 5.610 |
| w1-kd1000 | 1e-6 | best 4x18 | 5.528 | 0.418 | 0.985 | 1.593 | 0.54 | 2.373 | 5.678 |
| w1-kd1000 | 1e-8 | 1 × N 1x72 | 7.601 | 0.900 | 2.648 | 1.979 | 0.73 | 1.760 | 9.867 |
| w1-kd1000 | 1e-8 | best 4x18 | 4.425 | 0.265 | 0.912 | 0.470 | 0.37 | 2.688 | 10.124 |
| w2-kd100 | 1e-3 | 1 × N 1x72 | 0.100 | 0.015 | 0.022 | 0.040 | 0.77 | 0.013 | 0.050 |
| w2-kd100 | 1e-3 | best 4x18 | 0.045 | 0.006 | 0.007 | 0.010 | 0.52 | 0.017 | 0.044 |
| w2-kd100 | 1e-6 | 1 × N 1x72 | 0.113 | 0.015 | 0.022 | 0.041 | 0.69 | 0.025 | 0.070 |
| w2-kd100 | 1e-6 | best 4x18 | 0.061 | 0.006 | 0.007 | 0.010 | 0.37 | 0.033 | 0.064 |
| w2-kd100 | 1e-8 | 1 × N 1x72 | 0.096 | 0.014 | 0.021 | 0.025 | 0.63 | 0.026 | 0.172 |
| w2-kd100 | 1e-8 | best 2x36 | 0.070 | 0.009 | 0.013 | 0.013 | 0.50 | 0.028 | 0.173 |
| w2-kd250 | 1e-3 | 1 × N 1x72 | 0.779 | 0.106 | 0.205 | 0.371 | 0.88 | 0.051 | 0.267 |
| w2-kd250 | 1e-3 | best 4x18 | 0.283 | 0.036 | 0.061 | 0.092 | 0.67 | 0.079 | 0.253 |
| w2-kd250 | 1e-6 | 1 × N 1x72 | 0.819 | 0.108 | 0.205 | 0.377 | 0.84 | 0.083 | 0.404 |
| w2-kd250 | 1e-6 | best 4x18 | 0.325 | 0.035 | 0.060 | 0.093 | 0.58 | 0.122 | 0.395 |
| w2-kd250 | 1e-8 | 1 × N 1x72 | 0.544 | 0.078 | 0.197 | 0.157 | 0.80 | 0.080 | 0.964 |
| w2-kd250 | 1e-8 | best 4x18 | 0.282 | 0.028 | 0.062 | 0.039 | 0.46 | 0.142 | 0.985 |
| w2-kd500 | 1e-3 | 1 × N 1x72 | 3.892 | 0.530 | 1.083 | 1.877 | 0.90 | 0.206 | 1.099 |
| w2-kd500 | 1e-3 | best 4x18 | 1.330 | 0.154 | 0.328 | 0.450 | 0.70 | 0.339 | 1.056 |
| w2-kd500 | 1e-6 | 1 × N 1x72 | 3.941 | 0.511 | 1.042 | 1.895 | 0.88 | 0.301 | 1.669 |
| w2-kd500 | 1e-6 | best 4x18 | 1.426 | 0.154 | 0.328 | 0.451 | 0.65 | 0.437 | 1.675 |
| w2-kd500 | 1e-8 | 1 × N 1x72 | 2.615 | 0.363 | 0.985 | 0.827 | 0.83 | 0.321 | 3.712 |
| w2-kd500 | 1e-8 | best 4x18 | 1.152 | 0.115 | 0.312 | 0.198 | 0.54 | 0.489 | 3.825 |
| w3-kd30 | 1e-3 | 1 × N 1x72 | 0.128 | 0.012 | 0.020 | 0.065 | 0.76 | 0.011 | 0.040 |
| w3-kd30 | 1e-3 | best 8x9 | 0.041 | 0.003 | 0.003 | 0.008 | 0.34 | 0.022 | 0.035 |
| w3-kd30 | 1e-6 | 1 × N 1x72 | 0.130 | 0.012 | 0.019 | 0.066 | 0.75 | 0.012 | 0.062 |
| w3-kd30 | 1e-6 | best 8x9 | 0.036 | 0.003 | 0.004 | 0.008 | 0.41 | 0.017 | 0.068 |
| w3-kd30 | 1e-8 | 1 × N 1x72 | 0.054 | 0.009 | 0.017 | 0.005 | 0.59 | 0.013 | 0.142 |
| w3-kd30 | 1e-8 | best 4x18 | 0.029 | 0.004 | 0.006 | 0.002 | 0.39 | 0.014 | 0.145 |
| w3-kd100 | 1e-3 | 1 × N 1x72 | 4.321 | 0.456 | 2.011 | 1.447 | 0.91 | 0.026 | 2.366 |
| w3-kd100 | 1e-3 | best 8x9 | 0.559 | 0.077 | 0.210 | 0.172 | 0.82 | 0.045 | 2.403 |
| w3-kd100 | 1e-6 | 1 × N 1x72 | 4.321 | 0.449 | 2.027 | 1.445 | 0.91 | 0.021 | 2.970 |
| w3-kd100 | 1e-6 | best 8x9 | 0.560 | 0.078 | 0.211 | 0.175 | 0.83 | 0.043 | 3.098 |
| w3-kd100 | 1e-8 | 1 × N 1x72 | 2.600 | 0.349 | 1.339 | 0.636 | 0.89 | 0.022 | 4.156 |
| w3-kd100 | 1e-8 | best 8x9 | 0.445 | 0.066 | 0.199 | 0.076 | 0.77 | 0.061 | 4.432 |
| w3-n1000000-kd200 | 1e-3 | 1 × N 1x72 | 1.277 | 0.111 | 0.316 | 0.623 | 0.82 | 0.050 | 0.474 |
| w3-n1000000-kd200 | 1e-3 | best 8x9 | 0.241 | 0.019 | 0.035 | 0.074 | 0.53 | 0.088 | 0.529 |
| w3-n1000000-kd200 | 1e-6 | 1 × N 1x72 | 1.277 | 0.111 | 0.300 | 0.640 | 0.82 | 0.048 | 0.684 |
| w3-n1000000-kd200 | 1e-6 | best 8x9 | 0.287 | 0.019 | 0.032 | 0.076 | 0.44 | 0.134 | 0.772 |
| w3-n1000000-kd200 | 1e-8 | 1 × N 1x72 | 0.571 | 0.081 | 0.269 | 0.123 | 0.83 | 0.051 | 1.454 |
| w3-n1000000-kd200 | 1e-8 | best 4x18 | 0.247 | 0.031 | 0.084 | 0.031 | 0.59 | 0.082 | 1.557 |
| w3-n1000000-kd500 | 1e-3 | 1 × N 1x72 | 1.433 | 0.110 | 0.303 | 0.644 | 0.74 | 0.197 | 1.439 |
| w3-n1000000-kd500 | 1e-3 | best 4x18 | 0.637 | 0.035 | 0.080 | 0.155 | 0.42 | 0.320 | 1.661 |
| w3-n1000000-kd500 | 1e-6 | 1 × N 1x72 | 1.496 | 0.112 | 0.310 | 0.644 | 0.71 | 0.253 | 2.131 |
| w3-n1000000-kd500 | 1e-6 | best 4x18 | 1.007 | 0.051 | 0.131 | 0.275 | 0.45 | 0.473 | 2.449 |
| w3-n1000000-kd500 | 1e-8 | 1 × N 1x72 | 2.747 | 0.186 | 0.489 | 1.486 | 0.79 | 0.457 | 3.708 |
| w3-n1000000-kd500 | 1e-8 | best 4x18 | 1.102 | 0.037 | 0.095 | 0.370 | 0.46 | 0.547 | 4.177 |
| w1-kd500-planewave | 1e-3 | 1 × N 1x72 | 4.373 | 0.555 | 0.936 | 2.420 | 0.89 | 0.252 | 0.886 |
| w1-kd500-planewave | 1e-3 | best 4x18 | 1.533 | 0.167 | 0.297 | 0.611 | 0.70 | 0.390 | 0.898 |
| w1-kd500-planewave | 1e-6 | 1 × N 1x72 | 4.664 | 0.487 | 0.938 | 2.591 | 0.86 | 0.424 | 1.389 |
| w1-kd500-planewave | 1e-6 | best 4x18 | 1.621 | 0.162 | 0.294 | 0.569 | 0.63 | 0.527 | 1.440 |
| w1-kd500-planewave | 1e-8 | 1 × N 1x72 | 2.532 | 0.346 | 0.897 | 0.736 | 0.78 | 0.436 | 2.792 |
| w1-kd500-planewave | 1e-8 | best 4x18 | 1.302 | 0.104 | 0.284 | 0.179 | 0.44 | 0.695 | 2.871 |
| w1-kd1000-planewave | 1e-3 | 1 × N 1x72 | 16.816 | 2.104 | 3.711 | 8.583 | 0.86 | 1.658 | 4.082 |
| w1-kd1000-planewave | 1e-3 | best 4x18 | 5.836 | 0.425 | 0.977 | 1.930 | 0.57 | 2.293 | 4.284 |
| w1-kd1000-planewave | 1e-6 | 1 × N 1x72 | 16.973 | 1.942 | 3.510 | 8.694 | 0.83 | 2.091 | 6.376 |
| w1-kd1000-planewave | 1e-6 | best 4x18 | 6.712 | 0.510 | 1.153 | 2.064 | 0.56 | 2.814 | 6.462 |
| w1-kd1000-planewave | 1e-8 | 1 × N 1x72 | 9.585 | 1.233 | 3.326 | 2.402 | 0.73 | 2.299 | 11.396 |
| w1-kd1000-planewave | 1e-8 | best 4x18 | 4.351 | 0.266 | 0.903 | 0.466 | 0.38 | 2.625 | 10.032 |
| w2-kd100-constant | 1e-3 | 1 × N 1x72 | 0.139 | 0.016 | 0.036 | 0.059 | 0.80 | 0.016 | 0.052 |
| w2-kd100-constant | 1e-3 | best 8x9 | 0.048 | 0.004 | 0.004 | 0.005 | 0.27 | 0.031 | 0.043 |
| w2-kd100-constant | 1e-6 | 1 × N 1x72 | 0.111 | 0.014 | 0.021 | 0.041 | 0.69 | 0.025 | 0.085 |
| w2-kd100-constant | 1e-6 | best 4x18 | 0.061 | 0.006 | 0.006 | 0.010 | 0.37 | 0.034 | 0.067 |
| w2-kd100-constant | 1e-8 | 1 × N 1x72 | 0.118 | 0.016 | 0.034 | 0.033 | 0.70 | 0.025 | 0.173 |
| w2-kd100-constant | 1e-8 | best 4x18 | 0.061 | 0.006 | 0.007 | 0.007 | 0.32 | 0.037 | 0.215 |
| w3-kd30-constant | 1e-3 | 1 × N 1x72 | 0.154 | 0.012 | 0.019 | 0.080 | 0.72 | 0.018 | 0.048 |
| w3-kd30-constant | 1e-3 | best 8x9 | 0.042 | 0.003 | 0.004 | 0.008 | 0.34 | 0.022 | 0.037 |
| w3-kd30-constant | 1e-6 | 1 × N 1x72 | 0.134 | 0.014 | 0.020 | 0.066 | 0.74 | 0.012 | 0.077 |
| w3-kd30-constant | 1e-6 | best 8x9 | 0.035 | 0.003 | 0.004 | 0.008 | 0.41 | 0.016 | 0.087 |
| w3-kd30-constant | 1e-8 | 1 × N 1x72 | 0.053 | 0.009 | 0.016 | 0.005 | 0.58 | 0.013 | 0.173 |
| w3-kd30-constant | 1e-8 | best 2x36 | 0.037 | 0.006 | 0.011 | 0.003 | 0.53 | 0.011 | 0.166 |

## 6. Memory

Ours at the best split (the node's sum, the largest rank's in brackets) and at 1 × N, beside
FMM3D's peak at each match. FMM3D's peak is its whole process (inputs, outputs, tree, expansions);
ours the FMM's growth over its build and evaluations, the problem's arrays not included: the ratio
favours ours by the arrays (W1 at kD = 1000: 8·10⁶ points, about 0.8 GB of coordinates, charges,
potentials and gradients).

Kathleen:

| Case | N | ε | ours best split: node GB (largest rank) | ours 1 × N GB | FMM3D closest: variant eps GB | FMM3D one-sided: variant eps GB |
| --- | ---: | --- | --- | --- | --- | --- |
| w1-kd250 | 497360 | 1e-3 | 0.96 (rank ≤ 0.48) | 0.75 (rank ≤ 0.75) | fastker 3e-3 5.97 | fastker 3e-3 5.97 |
| w1-kd250 | 497360 | 1e-6 | 1.36 (rank ≤ 0.68) | 1.02 (rank ≤ 1.02) | fastker 1e-5 2.23 | fastker 3e-4 2.16 |
| w1-kd250 | 497360 | 1e-8 | 1.85 (rank ≤ 0.92) | 1.37 (rank ≤ 1.37) | fastker 3e-7 2.80 | fastker 3e-7 2.80 |
| w1-kd500 | 1989437 | 1e-3 | 4.55 (rank ≤ 2.29) | 3.26 (rank ≤ 3.26) | fastker 3e-3 26.88 | fastker 3e-3 26.88 |
| w1-kd500 | 1989437 | 1e-6 | 4.99 (rank ≤ 2.51) | 3.78 (rank ≤ 3.78) | fastker 3e-4 9.50 | fastker 3e-4 9.50 |
| w1-kd500 | 1989437 | 1e-8 | 6.66 (rank ≤ 3.34) | 5.11 (rank ≤ 5.11) | fastker 1e-8 12.48 | fastker 3e-7 12.20 |
| w1-kd1000 | 7957748 | 1e-3 | 15.48 (rank ≤ 7.75) | 12.07 (rank ≤ 12.07) | fastker 3e-3 117.02 | fastker 3e-3 117.02 |
| w1-kd1000 | 7957748 | 1e-6 | 18.97 (rank ≤ 9.51) | 15.24 (rank ≤ 15.24) | fastker 3e-4 43.53 | fastker 1e-4 43.73 |
| w1-kd1000 | 7957748 | 1e-8 | 24.52 (rank ≤ 12.28) | 20.08 (rank ≤ 20.08) | fastker 3e-5 43.99 | fastker 3e-7 55.80 |
| w2-kd100 | 127554 | 1e-3 | 0.30 (rank ≤ 0.15) | 0.22 (rank ≤ 0.22) | fastker 3e-3 0.95 | fastker 3e-4 0.61 |
| w2-kd100 | 127554 | 1e-6 | 0.39 (rank ≤ 0.20) | 0.30 (rank ≤ 0.30) | fastker 3e-4 0.61 | fastker 3e-7 0.86 |
| w2-kd100 | 127554 | 1e-8 | 0.70 (rank ≤ 0.35) | 0.51 (rank ≤ 0.51) | fastker 3e-7 0.86 | fastker 3e-7 0.86 |
| w2-kd250 | 797211 | 1e-3 | 1.25 (rank ≤ 0.62) | 0.95 (rank ≤ 0.95) | fastker 3e-3 11.53 | fastker 3e-3 11.53 |
| w2-kd250 | 797211 | 1e-6 | 1.84 (rank ≤ 0.92) | 1.43 (rank ≤ 1.43) | fastker 3e-4 4.38 | fastker 3e-4 4.38 |
| w2-kd250 | 797211 | 1e-8 | 2.70 (rank ≤ 1.36) | 2.08 (rank ≤ 2.08) | fastker 1e-7 5.77 | fastker 3e-7 5.72 |
| w2-kd500 | 3188840 | 1e-3 | 5.43 (rank ≤ 2.73) | 4.00 (rank ≤ 4.00) | fastker 3e-3 51.06 | fastker 3e-3 51.06 |
| w2-kd500 | 3188840 | 1e-6 | 6.80 (rank ≤ 3.40) | 5.48 (rank ≤ 5.48) | fastker 3e-4 17.39 | fastker 3e-4 17.39 |
| w2-kd500 | 3188840 | 1e-8 | 9.64 (rank ≤ 4.83) | 7.88 (rank ≤ 7.88) | fastker 3e-7 22.66 | fastker 3e-7 22.66 |
| w3-kd30 | 108849 | 1e-3 | 0.55 (rank ≤ 0.28) | 0.31 (rank ≤ 0.31) | fastker 1e-3 0.10 | fastker 3e-3 0.10 |
| w3-kd30 | 108849 | 1e-6 | 0.29 (rank ≤ 0.07) | 0.17 (rank ≤ 0.17) | fastker 3e-4 0.17 | fastker 3e-4 0.17 |
| w3-kd30 | 108849 | 1e-8 | 0.38 (rank ≤ 0.19) | 0.28 (rank ≤ 0.28) | fastker 1e-9 0.34 | fastker 1e-9 0.34 |
| w3-kd100 | 4031442 | 1e-3 | 3.99 (rank ≤ 0.51) | 2.28 (rank ≤ 2.28) | fastker 3e-3 15.52 | fastker 3e-4 8.76 |
| w3-kd100 | 4031442 | 1e-6 | 6.09 (rank ≤ 0.77) | 4.17 (rank ≤ 4.17) | fastker 3e-4 8.76 | fastker 3e-7 11.77 |
| w3-kd100 | 4031442 | 1e-8 | 10.79 (rank ≤ 1.36) | 6.71 (rank ≤ 6.71) | fastker 3e-7 11.77 | fastker 3e-7 11.77 |
| w3-n1000000-kd200 | 1000000 | 1e-3 | 2.43 (rank ≤ 0.61) | 1.42 (rank ≤ 1.42) | fastker 3e-3 6.51 | fastker 3e-4 8.02 |
| w3-n1000000-kd200 | 1000000 | 1e-6 | 2.58 (rank ≤ 0.65) | 1.60 (rank ≤ 1.60) | fastker 1e-4 8.06 | fastker 3e-7 9.70 |
| w3-n1000000-kd200 | 1000000 | 1e-8 | 2.80 (rank ≤ 1.40) | 2.13 (rank ≤ 2.13) | fastker 3e-7 9.70 | fastker 3e-7 9.70 |
| w3-n1000000-kd500 | 1000000 | 1e-3 | 5.21 (rank ≤ 2.68) | 3.75 (rank ≤ 3.75) | fastker 3e-3 3.92 | fastker 3e-3 3.92 |
| w3-n1000000-kd500 | 1000000 | 1e-6 | 5.35 (rank ≤ 2.71) | 3.79 (rank ≤ 3.79) | fastker 3e-5 5.03 | fastker 3e-4 4.95 |
| w3-n1000000-kd500 | 1000000 | 1e-8 | 8.16 (rank ≤ 4.11) | 5.88 (rank ≤ 5.88) | fastker 3e-7 6.25 | fastker 3e-7 6.25 |
| w1-kd500-planewave | 1989437 | 1e-3 | 4.53 (rank ≤ 2.28) | 3.25 (rank ≤ 3.25) | default 3e-3 107.20 | default 3e-3 107.20 |
| w1-kd500-planewave | 1989437 | 1e-6 | 5.00 (rank ≤ 2.50) | 3.83 (rank ≤ 3.83) | default 3e-5 48.30 | default 3e-4 47.53 |
| w1-kd500-planewave | 1989437 | 1e-8 | 6.71 (rank ≤ 3.37) | 5.10 (rank ≤ 5.10) | default 1e-7 43.15 | default 3e-7 42.78 |
| w1-kd1000-planewave | 7957748 | 1e-3 | 15.49 (rank ≤ 7.75) | 12.09 (rank ≤ 12.09) | fastker 1e-3 118.02 | fastker 1e-3 118.02 |
| w1-kd1000-planewave | 7957748 | 1e-6 | 18.98 (rank ≤ 9.50) | 15.24 (rank ≤ 15.24) | fastker 3e-4 43.50 | fastker 1e-4 43.72 |
| w1-kd1000-planewave | 7957748 | 1e-8 | 24.58 (rank ≤ 12.32) | 20.09 (rank ≤ 20.09) | fastker 3e-7 55.77 | fastker 3e-7 55.77 |
| w2-kd100-constant | 127554 | 1e-3 | 0.30 (rank ≤ 0.15) | 0.22 (rank ≤ 0.22) | default 3e-3 4.59 | default 3e-4 2.75 |
| w2-kd100-constant | 127554 | 1e-6 | 0.39 (rank ≤ 0.19) | 0.30 (rank ≤ 0.30) | default 1e-4 2.78 | default 3e-7 1.83 |
| w2-kd100-constant | 127554 | 1e-8 | 0.70 (rank ≤ 0.35) | 0.52 (rank ≤ 0.52) | default 3e-8 1.88 | default 3e-7 1.83 |
| w3-kd30-constant | 108849 | 1e-3 | 0.55 (rank ≤ 0.27) | 0.32 (rank ≤ 0.32) | default 3e-3 0.54 | default 3e-4 0.84 |
| w3-kd30-constant | 108849 | 1e-6 | 0.29 (rank ≤ 0.08) | 0.17 (rank ≤ 0.17) | default 3e-4 0.84 | default 3e-4 0.84 |
| w3-kd30-constant | 108849 | 1e-8 | 0.38 (rank ≤ 0.19) | 0.27 (rank ≤ 0.27) | default 3e-8 0.28 | default 3e-8 0.28 |

locust:

| Case | N | ε | ours best split: node GB (largest rank) | ours 1 × N GB | FMM3D closest: variant eps GB | FMM3D one-sided: variant eps GB |
| --- | ---: | --- | --- | --- | --- | --- |
| w1-kd250 | 497360 | 1e-3 | 1.69 (rank ≤ 0.43) | 1.03 (rank ≤ 1.03) | fastker 3e-3 6.05 | fastker 3e-3 6.05 |
| w1-kd250 | 497360 | 1e-6 | 2.64 (rank ≤ 0.66) | 1.65 (rank ≤ 1.65) | default 1e-5 10.91 | default 3e-4 10.64 |
| w1-kd250 | 497360 | 1e-8 | 3.57 (rank ≤ 0.89) | 2.06 (rank ≤ 2.06) | fastker 3e-7 2.93 | fastker 3e-7 2.93 |
| w1-kd500 | 1989437 | 1e-3 | 7.98 (rank ≤ 2.02) | 4.18 (rank ≤ 4.18) | fastker 3e-3 26.90 | fastker 3e-3 26.90 |
| w1-kd500 | 1989437 | 1e-6 | 9.87 (rank ≤ 2.48) | 6.26 (rank ≤ 6.26) | default 3e-4 47.72 | default 3e-4 47.72 |
| w1-kd500 | 1989437 | 1e-8 | 11.73 (rank ≤ 2.96) | 7.07 (rank ≤ 7.07) | fastker 1e-8 12.57 | fastker 3e-7 12.27 |
| w1-kd1000 | 7957748 | 1e-3 | 25.70 (rank ≤ 6.45) | 15.71 (rank ≤ 15.71) | fastker 3e-3 117.50 | fastker 3e-3 117.50 |
| w1-kd1000 | 7957748 | 1e-6 | 36.06 (rank ≤ 9.03) | 25.14 (rank ≤ 25.14) | fastker 3e-4 44.04 | fastker 1e-4 44.24 |
| w1-kd1000 | 7957748 | 1e-8 | 40.43 (rank ≤ 10.14) | 27.67 (rank ≤ 27.67) | fastker 3e-5 44.51 | fastker 3e-7 56.36 |
| w2-kd100 | 127554 | 1e-3 | 0.33 (rank ≤ 0.08) | 0.21 (rank ≤ 0.21) | fastker 3e-3 0.98 | default 3e-4 2.77 |
| w2-kd100 | 127554 | 1e-6 | 0.68 (rank ≤ 0.17) | 0.43 (rank ≤ 0.43) | default 3e-4 2.77 | fastker 3e-7 0.85 |
| w2-kd100 | 127554 | 1e-8 | 0.88 (rank ≤ 0.44) | 0.68 (rank ≤ 0.68) | fastker 3e-7 0.85 | fastker 3e-7 0.85 |
| w2-kd250 | 797211 | 1e-3 | 1.78 (rank ≤ 0.45) | 1.14 (rank ≤ 1.14) | fastker 3e-3 11.63 | fastker 3e-3 11.63 |
| w2-kd250 | 797211 | 1e-6 | 3.23 (rank ≤ 0.81) | 2.13 (rank ≤ 2.13) | default 3e-4 19.89 | default 3e-4 19.89 |
| w2-kd250 | 797211 | 1e-8 | 4.33 (rank ≤ 1.09) | 2.58 (rank ≤ 2.58) | fastker 1e-7 5.90 | fastker 3e-7 5.84 |
| w2-kd500 | 3188840 | 1e-3 | 8.75 (rank ≤ 2.20) | 4.92 (rank ≤ 4.92) | fastker 3e-3 51.23 | fastker 3e-3 51.23 |
| w2-kd500 | 3188840 | 1e-6 | 11.85 (rank ≤ 2.99) | 8.10 (rank ≤ 8.10) | default 3e-4 84.46 | default 3e-4 84.46 |
| w2-kd500 | 3188840 | 1e-8 | 14.45 (rank ≤ 3.65) | 9.34 (rank ≤ 9.34) | fastker 3e-7 22.91 | fastker 3e-7 22.91 |
| w3-kd30 | 108849 | 1e-3 | 0.69 (rank ≤ 0.09) | 0.17 (rank ≤ 0.17) | default 1e-3 0.56 | default 3e-3 0.55 |
| w3-kd30 | 108849 | 1e-6 | 0.94 (rank ≤ 0.13) | 0.43 (rank ≤ 0.43) | default 3e-4 0.87 | default 3e-4 0.87 |
| w3-kd30 | 108849 | 1e-8 | 0.64 (rank ≤ 0.16) | 0.34 (rank ≤ 0.34) | fastker 1e-9 0.36 | fastker 1e-9 0.36 |
| w3-kd100 | 4031442 | 1e-3 | 4.08 (rank ≤ 0.51) | 2.75 (rank ≤ 2.75) | fastker 3e-3 15.54 | default 3e-4 21.76 |
| w3-kd100 | 4031442 | 1e-6 | 9.23 (rank ≤ 1.18) | 6.95 (rank ≤ 6.95) | default 3e-4 21.76 | default 3e-7 29.58 |
| w3-kd100 | 4031442 | 1e-8 | 10.95 (rank ≤ 1.37) | 6.83 (rank ≤ 6.83) | default 3e-7 29.58 | default 3e-7 29.58 |
| w3-n1000000-kd200 | 1000000 | 1e-3 | 3.95 (rank ≤ 0.50) | 2.18 (rank ≤ 2.18) | default 3e-3 49.31 | fastker 3e-4 8.07 |
| w3-n1000000-kd200 | 1000000 | 1e-6 | 6.65 (rank ≤ 0.83) | 3.80 (rank ≤ 3.80) | fastker 1e-4 8.11 | fastker 3e-7 9.75 |
| w3-n1000000-kd200 | 1000000 | 1e-8 | 3.95 (rank ≤ 0.99) | 2.36 (rank ≤ 2.36) | fastker 3e-7 9.75 | fastker 3e-7 9.75 |
| w3-n1000000-kd500 | 1000000 | 1e-3 | 13.09 (rank ≤ 3.29) | 7.72 (rank ≤ 7.72) | default 3e-3 23.41 | default 3e-3 23.41 |
| w3-n1000000-kd500 | 1000000 | 1e-6 | 13.51 (rank ≤ 3.39) | 8.11 (rank ≤ 8.11) | fastker 3e-5 5.29 | fastker 3e-4 5.19 |
| w3-n1000000-kd500 | 1000000 | 1e-8 | 17.14 (rank ≤ 4.32) | 11.12 (rank ≤ 11.12) | fastker 3e-7 6.56 | fastker 3e-7 6.56 |
| w1-kd500-planewave | 1989437 | 1e-3 | 7.96 (rank ≤ 2.01) | 4.19 (rank ≤ 4.19) | default 3e-3 107.37 | default 3e-3 107.37 |
| w1-kd500-planewave | 1989437 | 1e-6 | 9.77 (rank ≤ 2.47) | 6.26 (rank ≤ 6.26) | default 3e-5 48.50 | default 3e-4 47.72 |
| w1-kd500-planewave | 1989437 | 1e-8 | 11.77 (rank ≤ 2.97) | 7.10 (rank ≤ 7.10) | default 1e-7 43.38 | default 3e-7 43.01 |
| w1-kd1000-planewave | 7957748 | 1e-3 | 25.88 (rank ≤ 6.51) | 15.66 (rank ≤ 15.66) | fastker 1e-3 118.50 | fastker 1e-3 118.50 |
| w1-kd1000-planewave | 7957748 | 1e-6 | 36.11 (rank ≤ 9.05) | 25.15 (rank ≤ 25.15) | fastker 3e-4 44.06 | fastker 1e-4 44.24 |
| w1-kd1000-planewave | 7957748 | 1e-8 | 40.53 (rank ≤ 10.17) | 27.64 (rank ≤ 27.64) | fastker 3e-7 56.39 | fastker 3e-7 56.39 |
| w2-kd100-constant | 127554 | 1e-3 | 0.48 (rank ≤ 0.06) | 0.21 (rank ≤ 0.21) | default 3e-3 4.63 | default 3e-4 2.78 |
| w2-kd100-constant | 127554 | 1e-6 | 0.69 (rank ≤ 0.17) | 0.43 (rank ≤ 0.43) | default 1e-4 2.77 | default 3e-7 1.87 |
| w2-kd100-constant | 127554 | 1e-8 | 1.22 (rank ≤ 0.31) | 0.68 (rank ≤ 0.68) | default 3e-8 1.88 | default 3e-7 1.87 |
| w3-kd30-constant | 108849 | 1e-3 | 0.70 (rank ≤ 0.09) | 0.17 (rank ≤ 0.17) | default 3e-3 0.55 | default 3e-4 0.87 |
| w3-kd30-constant | 108849 | 1e-6 | 0.92 (rank ≤ 0.12) | 0.43 (rank ≤ 0.43) | default 3e-4 0.87 | default 3e-4 0.87 |
| w3-kd30-constant | 108849 | 1e-8 | 0.45 (rank ≤ 0.22) | 0.34 (rank ≤ 0.34) | default 3e-8 0.30 | default 3e-8 0.30 |

FMM3D `default` at kD = 1000 (locust; not run on Kathleen, CH4.4a's question 3): 135–191 s per call
against `fastker`'s 110–178 s, peaks 442–447 GB at 3·10⁻³ and 10⁻³ and 186–198 GB below, against
`fastker`'s 118 GB and 44–58 GB. It is timed (`time-w1-kd1000*.md` in locust's archive) but not in
the comparison's tables, whose variants are the grid's (`fastker` alone at kD = 1000): it is 1.2 %
and 0.2 % faster than `fastker` at 3·10⁻⁴ and 10⁻⁴ and slower at every other eps, so no ratio changes
by more than about 1 %. On Kathleen
`default` failed on W2 at kD = 500 at 3·10⁻³ and 10⁻³ (exit 137, Slurm's oom_kill, the node's 180 GB;
job 663814's peak 188.7 GB) and in both timed launches at 10⁻³; it is no row's match there.

## 7. c_leaf (CH4.3e's Also)

W3 at N = 10⁶, kD = 500; build, evaluation and their sum at c_leaf 2 and 4, and 4 ÷ 2.

Kathleen:

| ε | split | c_leaf 2: build, evaluate, s | c_leaf 4: build, evaluate, s | 4 ÷ 2 | φ, ∇φ (4) |
| --- | --- | --- | --- | ---: | --- |
| 1e-3 | 1x40 | 13.492, 6.081, 19.573 | 2.944, 3.426, 6.370 | 0.33 | 3.50e-04, 4.50e-04 |
| 1e-3 | 2x20 | 6.544, 5.393, 11.937 | 1.951, 3.036, 4.987 | 0.42 | 3.50e-04, 4.50e-04 |
| 1e-3 | 4x10 | 4.574, 5.665, 10.240 | 2.230, 3.199, 5.429 | 0.53 | 3.50e-04, 4.50e-04 |
| 1e-3 | 8x5 | 5.241, 6.342, 11.583 | 4.062, 3.680, 7.742 | 0.67 | 3.50e-04, 4.50e-04 |
| 1e-3 | 40x1 | 11.311, 10.135, 21.447 | 11.117, 4.798, 15.915 | 0.74 | 3.50e-04, 4.50e-04 |
| 1e-6 | 1x40 | 17.016, 9.239, 26.255 | 3.043, 5.213, 8.257 | 0.31 | 1.91e-07, 2.53e-07 |
| 1e-6 | 2x20 | 8.538, 8.100, 16.638 | 2.106, 4.032, 6.139 | 0.37 | 1.91e-07, 2.53e-07 |
| 1e-6 | 4x10 | 5.468, 8.522, 13.990 | 2.274, 4.296, 6.570 | 0.47 | 1.91e-07, 2.53e-07 |
| 1e-6 | 8x5 | 6.553, 9.417, 15.971 | 5.279, 5.034, 10.313 | 0.65 | 1.91e-07, 2.53e-07 |
| 1e-6 | 40x1 | 14.440, 14.927, 29.367 | 13.878, 6.467, 20.345 | 0.69 | 1.91e-07, 2.53e-07 |
| 1e-8 | 1x40 | 17.339, 16.277, 33.615 | 3.102, 6.755, 9.857 | 0.29 | 2.04e-09, 3.42e-09 |
| 1e-8 | 2x20 | 8.664, 15.387, 24.051 | 2.176, 6.327, 8.503 | 0.35 | 2.04e-09, 3.42e-09 |
| 1e-8 | 4x10 | 5.793, 16.492, 22.285 | 2.329, 7.072, 9.401 | 0.42 | 2.04e-09, 3.42e-09 |
| 1e-8 | 8x5 | 6.114, 19.677, 25.791 | 3.904, 9.400, 13.304 | 0.52 | 2.04e-09, 3.42e-09 |
| 1e-8 | 40x1 | 16.045, 44.922, 60.967 | 15.287, 22.810, 38.098 | 0.62 | 2.04e-09, 3.42e-09 |

locust:

| ε | split | c_leaf 2: build, evaluate, s | c_leaf 4: build, evaluate, s | 4 ÷ 2 | φ, ∇φ (4) |
| --- | --- | --- | --- | ---: | --- |
| 1e-3 | 1x72 | 12.075, 1.747, 13.822 | 1.717, 1.593, 3.311 | 0.24 | 3.65e-04, 4.65e-04 |
| 1e-3 | 2x36 | 5.979, 1.855, 7.834 | 1.092, 1.684, 2.776 | 0.35 | 3.65e-04, 4.65e-04 |
| 1e-3 | 4x18 | 3.403, 2.011, 5.414 | 0.752, 1.844, 2.597 | 0.48 | 3.65e-04, 4.65e-04 |
| 1e-3 | 8x9 | 2.080, 2.160, 4.240 | 0.897, 1.978, 2.875 | 0.68 | 3.65e-04, 4.65e-04 |
| 1e-3 | 72x1 | 3.977, 7.533, 11.510 | 3.194, 7.741, 10.935 | 0.95 | 3.65e-04, 4.65e-04 |
| 1e-6 | 1x72 | 9.813, 2.099, 11.912 | 1.508, 1.911, 3.420 | 0.29 | 1.92e-07, 2.54e-07 |
| 1e-6 | 2x36 | 4.887, 2.214, 7.101 | 1.195, 2.296, 3.491 | 0.49 | 1.92e-07, 2.54e-07 |
| 1e-6 | 4x18 | 2.685, 2.431, 5.117 | 0.931, 2.496, 3.427 | 0.67 | 1.92e-07, 2.54e-07 |
| 1e-6 | 8x9 | 1.804, 2.638, 4.442 | 1.173, 2.763, 3.935 | 0.89 | 1.92e-07, 2.54e-07 |
| 1e-6 | 72x1 | 2.241, 10.372, 12.614 | 3.099, 10.941, 14.040 | 1.11 | 1.92e-07, 2.54e-07 |
| 1e-8 | 1x72 | 9.916, 3.505, 13.421 | 1.634, 3.219, 4.853 | 0.36 | 2.04e-09, 3.42e-09 |
| 1e-8 | 2x36 | 4.934, 3.631, 8.566 | 1.082, 3.406, 4.487 | 0.52 | 2.04e-09, 3.42e-09 |
| 1e-8 | 4x18 | 2.755, 3.863, 6.618 | 0.838, 3.619, 4.456 | 0.67 | 2.04e-09, 3.42e-09 |
| 1e-8 | 8x9 | 1.883, 4.553, 6.436 | 0.979, 3.868, 4.847 | 0.75 | 2.04e-09, 3.42e-09 |
| 1e-8 | 72x1 | 2.720, 18.746, 21.466 | 2.499, 14.770, 17.269 | 0.80 | 2.04e-09, 3.42e-09 |

c_leaf 4 is faster at every split and at 40 × 1 on Kathleen (0.29–0.74×), and at every split on
locust (0.24–0.89×); at 72 × 1 0.80–1.11×. Its gain is the build (Kathleen 1 × 40: 13–17 s → 3 s);
the errors are c_leaf 4's within ε (3.5·10⁻⁴, 1.9·10⁻⁷, 2.0·10⁻⁹). The default stands.

## 8. FMM3D's OpenMP binding (B8)

The driver prints `omp_get_proc_bind()`, `OMP_PROC_BIND` and `OMP_PLACES` on every run: at the
defaults `proc_bind=0` (false), neither variable set, 40 or 72 threads; bound `proc_bind=3` (close). Bound beside unbound,
the smaller of two launches each:

Kathleen:

| Case | variant | eps | unbound s (launches) | bound s (launches) | bound ÷ unbound |
| --- | --- | --- | ---: | ---: | ---: |
| w1-kd250 | default | 1e-5 | 8.719 (2) | 9.175 (2) | 1.05 |
| w1-kd250 | fastker | 1e-5 | 4.919 (2) | 5.077 (2) | 1.03 |
| w1-kd500 | default | 3e-4 | 47.662 (2) | 50.852 (2) | 1.07 |
| w1-kd500 | fastker | 3e-4 | 33.368 (2) | 34.823 (2) | 1.04 |
| w1-kd1000 | fastker | 3e-4 | 305.105 (2) | 347.203 (2) | 1.14 |
| w3-kd100 | default | 3e-4 | 28.730 (2) | 28.950 (2) | 1.01 |
| w3-kd100 | fastker | 3e-4 | 20.137 (2) | 20.277 (2) | 1.01 |
| w3-n1000000-kd500 | default | 3e-5 | 73.408 (2) | 74.848 (2) | 1.02 |
| w3-n1000000-kd500 | fastker | 3e-5 | 61.755 (2) | 62.891 (2) | 1.02 |

locust:

| Case | variant | eps | unbound s (launches) | bound s (launches) | bound ÷ unbound |
| --- | --- | --- | ---: | ---: | ---: |
| w1-kd250 | default | 1e-5 | 3.767 (2) | 3.922 (2) | 1.04 |
| w1-kd250 | fastker | 1e-5 | 3.784 (2) | 3.884 (2) | 1.03 |
| w1-kd500 | default | 3e-4 | 19.191 (2) | 19.509 (2) | 1.02 |
| w1-kd500 | fastker | 3e-4 | 19.375 (2) | 19.680 (2) | 1.02 |
| w1-kd1000 | fastker | 3e-4 | 136.291 (2) | 136.556 (2) | 1.00 |
| w3-kd100 | default | 3e-4 | 10.829 (2) | 10.923 (2) | 1.01 |
| w3-kd100 | fastker | 3e-4 | 37.506 (2) | 37.315 (2) | 0.99 |
| w3-n1000000-kd500 | default | 3e-5 | 29.454 (2) | 29.556 (2) | 1.00 |
| w3-n1000000-kd500 | fastker | 3e-5 | 28.176 (2) | 28.138 (2) | 1.00 |

Binding never helps FMM3D here (Kathleen 1.01–1.14× slower, the most on W1 at kD = 1000; locust
0.99–1.04×), so every comparison is at its defaults, as CH4.0e and T9 ran it. On Kathleen the eight
smallest cases' first timed launches also ran bound (§11): their bound times are in the archive and
not used.

## 9. Checksums

Kathleen:

ours: 78 case-ε pairs, 624 launches; checksums differing: 0
cleaf2: 3 case-ε pairs, 30 launches; checksums differing: 0
cleaf4: 3 case-ε pairs, 30 launches; checksums differing: 0

locust:

ours: 78 case-ε pairs, 624 launches; checksums differing: 0
cleaf2: 3 case-ε pairs, 30 launches; checksums differing: 3
  ('w3-n1000000-kd500', '1e-3'): ['480d3794ce374250', '5c694c492d568e77']
  ('w3-n1000000-kd500', '1e-6'): ['8d17722e2be47947', 'c14ecf6d58d6d0fe']
  ('w3-n1000000-kd500', '1e-8'): ['83e07b9a4b74ad72', 'cb05a822bab3fd44']
cleaf4: 3 case-ε pairs, 30 launches; checksums differing: 3
  ('w3-n1000000-kd500', '1e-3'): ['63b8cbe4b21cc5c5', '7c79e694fa7b2e23']
  ('w3-n1000000-kd500', '1e-6'): ['3581e2114c15f0af', '844472e608c0b6b7']
  ('w3-n1000000-kd500', '1e-8'): ['57f34d85ca090362', '757353c1a55303a9']

Ours is rank- and thread-count invariant at the defaults on every case over the four splits of each
machine and, on Kathleen, at 40 × 1 (constraints 2 and 9). At 72 × 1 on locust the defaults build a
different plan (§0, §7): its own checksum, the same in both launches.

## 10. FMM3D's errors on both machines

The matches are made on Kathleen's grid. locust's FMM3D errors at its timed eps against the grid at
the same eps and variant: 508 values (φ and ∇φ), the largest ratio 1.064 (W3 at kD = 10, `fastker`,
1.46·10⁻¹⁵ against 1.56·10⁻¹⁵); no match changes.

## 11. Runs

Release builds, BLAS threads 1, one ε per launch. Raw outputs in
`~/projects/fmm-archive/kathleen/runs/ch4-5/` and `~/projects/fmm-archive/locust/runs/ch4-5/` (not in
the repository): every launch's `.md`, `.err` and `.log`, `runs.log` (each launch's exit, seconds and,
on locust, load and free memory), the eps lists and the job logs.

| machine | job / run | stage | nodes | elapsed | state |
| --- | --- | --- | --- | --- | --- |
| Kathleen | 663813 | grid, kD = 1000 | node-c11a-017 | 2:47 | completed |
| Kathleen | 663814 | grid, the other 24 | node-c11a-118 | 2:15 | `default` killed twice (W2, kD = 500), else completed |
| Kathleen | 663816 | ours, cleaf | node-c11a-118 | 2:13 | completed |
| Kathleen | 663952, 663953 | time | node-c11b-028, node-c11a-118 | 2:00 | cancelled (below) |
| Kathleen | 664298 | time, kD = 1000 | node-c11a-083 | 2:40 | completed |
| Kathleen | 664306 | time (16 cases), bind | node-c11b-028 | 2:21 | `default` at 10⁻³ killed (W2, kD = 500), else completed |
| locust | detached, 13:55–14:56 | ours, cleaf | – | 1:01 | 453–456 GB available |
| locust | detached, 17:28–21:14 | time, bind | – | 3:46 | |
| locust | 21:16 | time, W2 at kD = 100 and W1 constant at kD = 30 | – | 0:02 | Kathleen's errors add one eps each |

Every Kathleen job excludes node-c11b-004 and node-c11b-013 (CH4.4a). Jobs 663952 and 663953 ran a
copy of `ch4-5-runs.sh` synced before the `time` stage was trimmed to one binding (both bindings,
REPEAT 2 everywhere; the kD = 1000 job would have run about 10 h): cancelled after 2 h; the eight
cases 663953 had finished (W1 at kD = 10–500, W2 at 10–100, every N < 10⁶) are kept, their unbound
runs being the trimmed stage's runs; the rest reran as 664298 and 664306. locust's timed list came from
locust's errors of ours (the plan's "On locust the matches are made with locust's own errors");
Kathleen's adds two eps, run on locust after the stage.

Commands (`tools/helmholtz/ch4-5-runs.sh <out> <stage…>`, `CASES` to restrict; Kathleen through
`tools/kathleen/submit.sh --time=06:00:00 --exclude=node-c11b-004,node-c11b-013 -- run.sbatch …`;
locust detached under `setsid nohup`), the eps lists by `ch4-5-tables.py eps <grid> <ours…>`
(`--default-1000` on locust) and `bind`, the tables by `ch4-5-tables.py <table> <grid> <dir>`.

## 12. Questions for sign-off

Signed off as recommended (the user on #257, 2026-11-24: "Accept all 4 decisions as recommended.").

1. **The acceptance as read.** "For each of our tolerances the eps of FMM3D with the closest error":
   per case, ε and variant the eps minimising the worse of φ's and ∇φ's |log₁₀ ratio| over a grid of
   two eps per decade; FMM3D's faster variant among those whose match is within a decade of the
   closer one's; CH4.0e's one-sided match beside it; time as one solve (build plus one evaluation
   over the call) and setup apart (evaluation over the call less its tree); memory as §6 measures it;
   one node, best split and 1 × N. *Recommended:* accept; H4's tasks are then all done and the H4
   gate (`/milestone-gate H4`) is next.
2. **CH4.2d** (dispatched kernels for the diagonal pair and the conversions; off the list,
   reconsidered at CH4.5). Ours' evaluation alone is 0.03–0.52× FMM3D's call less its tree on every
   case of N ≥ 10⁵, and the evaluation is not where FMM3D wins (§3.4). *Recommended:* keep it off the
   list. Alternative: put it on the list after H4's gate.
3. **CH4.3f** (the time rule's costs at many ranks per node; off the list until CH4.5). N × 1 is not
   an operating point this comparison needs: the best split and 1 × N are (c_leaf 4's N × 1 rows
   are 3.2–4.5× the best split's on Kathleen and 3.9–4.2× on locust). *Recommended:* keep it off the list.
   Alternative: put it on the list after H4's gate.
4. **CH4.7, new: the form choice prices the tables' build for small trees** (§3.4). *Done when:* the
   time rule (or the leaf rule) weighs the build of Batched's matrices and the diagonal tables
   against the evaluation where one solve is the use, so that the six Kathleen rows the matrices
   explain (W1 and W2 at kD = 30 and W3 at kD = 10, both densities, 10⁻³; N ≤ 1.2·10⁴) are not
   slower than FMM3D's call, verified on CH4.5's small cases on Kathleen and locust with the large
   cases unchanged; the other four rows (a few milliseconds of fixed steps) are out of its scope.
   *Recommended:* add it, off the list (the ten rows lose 0.002–0.074 s). Alternative: on the list
   before the H4 gate.
