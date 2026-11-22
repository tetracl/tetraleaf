# Helmholtz CH4.3b2c: the diagonal form by default

Task: docs/helmholtz/tasks/CH4.3b2c.md (roadmap CH4.3b2c, the last of CH4.3b2a–d; proposal §10
CH4.3, "the time crossover per class; the fallback of a whole level to RCR when the scratch budget
binds"; §3.6, §5.3, §5.4). CH4.3b1 built the diagonal form's time rule (`HelmholtzDiagonal::Auto`,
opt-in), CH4.3b2a ran it at matched accuracy, CH4.3b2b and CH4.3b2d settled the V near rule on its
levels and CH4.3c2 refitted the degrees. This task measures the rule against the defaults as they
now are and proposes it as the default.

## 0. Summary

1. **The change** (§2): `HelmholtzFmmBuilder::new` takes `HelmholtzDiagonal::Auto`, and so does
   the enum's `Default`; the operator's own default before `with_diagonal` stays `Off`. For
   sign-off (question 2).
2. **Matched accuracy** (§3): the default errs 0.991–1.034× the current defaults' (φ and ∇φ) on
   every case of W1 (kD = 100–1000) and W3 (N = 10⁶, kD = 200, 500) at 10⁻³–10⁻⁸, on Kathleen
   (1 × 40) and locust (1 × 72), and at 40 × 1 and 72 × 1.
3. **The evaluation** (§4): over both run orders, 0.37–0.99× the defaults' on Kathleen and
   0.44–0.99× on locust; no case slower. Largest on W1 at kD = 1000 (0.37×, 0.44×) and kD = 500
   (0.47–0.74×) and W3 at kD = 500 (0.64–0.72×); least at 10⁻⁷–10⁻⁸ on W1 at kD ≤ 200 and on W3
   at kD = 200 (0.86–0.99×). The build grows by the diagonal tables (0.01–2.0 s at one rank per
   node); build plus one evaluation 0.59–1.11× at one rank per node (above 1.01 only W3 at
   kD = 200 on Kathleen in the off-first order, §6, and W1 at kD = 100, 10⁻⁷, 1.08× of 0.29 s).
4. **The order per level** (§5): the form takes 45 of Kathleen's and 40 of locust's levels with
   waves, and every RCR level of W1 at kD ≥ 500 and W3 at kD = 500 it reaches; RCR takes one level
   from the waves on Kathleen and two on locust. Within 1.05× of the fastest option on 75 of 78
   levels (Kathleen; 1.46 s lost of 41.7 s) and 82 of 90 (locust; 0.28 s of 12.5 s); every miss
   keeps the waves where the form was faster, on W3's dense levels and two raised W1 levels, as in
   CH4.3b2a. This moves D8's choice under `Auto` (question 3).
5. **The whole-level fallback** (§7): at 40 × 1 on W1 at kD = 1000, 10⁻⁶, levels 2 and 3 (whose
   scratch does not fit a rank's 205 MiB) run RCR whole on every rank although the model takes the
   form there; levels 4–7 take the form; errors 1.000×, the evaluation 53.0 s against 63.7 s.
6. **Many ranks per node** (§8): at 40 × 1 and 72 × 1 the choices are within 1.05× on 16 of 17 and
   19 of 20 levels, but the one miss on Kathleen costs 18 s: on W3 at kD = 500, 10⁻⁸, level 4 the
   waves take 20.5 s at 40 × 1 against 2.5 s at 1 × 40 (the waves pass 19.0 s against 0.83 s),
   while the form takes 2.6 s; the one-rank constant does not see it. And every rank builds every
   diagonal table (about 10 s at 40 × 1 on Kathleen, 1.0 s at 72 × 1 on locust): build plus one
   evaluation 0.91–1.22× at 40 × 1. Question 4: a follow-up task.
7. **W3 at kD = 200 on Kathleen** (§6, CH4.3b2a's question 4): a run-order effect of Kathleen's
   node, not the form. With the defaults first the time rule ran 1.06–1.09× their evaluation, with
   it first 0.68–0.90×; within one launch of `off, auto, off, auto` the second `off` ran 1.09–1.20×
   the first at every ε. Locust shows neither. Question 5.

## 1. Setup

- **Runs.** `helmholtz_accuracy` (release, BLAS threads 1), launched by
  `tools/helmholtz/ch4-3b2c-runs.sh` at 0272673 (the branch's code before the default changed;
  every launch names `--diagonal`): §9.1's complex charges and 1,000 samples,
  gradients, the library's defaults (`--m2l auto-pw`, the V, W and X rules of the cost model,
  CH4.3b2d's `PerBox`, CH4.3c2's degrees, the interpolated build), two evaluations (the rows the
  second's). Stages:
  - `run-*` (1 × cores): W1 at kD = 100, 200, 500 at 10⁻³, 10⁻⁶, 10⁻⁷, 10⁻⁸ and kD = 1000 at 10⁻⁶,
    W3 at N = 10⁶, kD = 200 and 500 at the same ε; `--diagonal off,auto --diagonal-costs
    machine,form,rcr`: the defaults, the time rule (the proposed default) and the time rule forced
    to the form and to RCR where a class passes (CH4.3b2a's `FORM_COSTS`, `RCR_COSTS`), each level
    at its degree as the defaults run it;
  - `reverse-*`: the same cases, `--diagonal auto,off`, the machine's costs;
  - `repeat-*`: W3 at kD = 200 at the four ε, `--diagonal off,auto,off,auto`;
  - `ranks-*` (one rank per core): W1 and W3 at kD = 500 at 10⁻⁶ and 10⁻⁸ with the three costs,
    and on Kathleen W1 at kD = 1000, 10⁻⁶ with `--diagonal off,auto`.
- **Machines.**
  - Kathleen job 663643 (node-c11b-048, `singlenode`, 2026-11-22 11:44–12:14, `run`; one rank bound
    to every core) and job 663644 (node-c11a-058, 11:45–12:21, `repeat`, `ranks`, `reverse`; the
    ranks 20 per socket, one core each); every launch exit 0 (21 and 30).
  - locust (GH200, 72 cores), 11:45–12:08, every stage, every launch exit 0 (50); the load sampled
    every 30 s: 1-minute load at most 53 under the 72-thread launches and 73 under 72 ranks, within
    the runs' own processes.
- **Raw outputs**: `~/projects/fmm-archive/kathleen/runs/ch4-3b2c/` (`run`, the jobs' logs) and
  `ch4-3b2c-b/` (`repeat`, `ranks`, `reverse`); `~/projects/fmm-archive/locust/runs/ch4-3b2c/`
  (every stage and `load.log`).
- **Tables**: `tools/helmholtz/ch4-3b2c-tables.py errors|choices|summary|build|repeat|order`
  (the first three CH4.3b2a's, with a level a forced run does not list left out of its options).

## 2. The change

`HelmholtzFmmBuilder::new` sets `diagonal: HelmholtzDiagonal::Auto`; `HelmholtzDiagonal`'s
`#[default]` moves to `Auto` (as `HelmholtzM2l`'s did at CH3.2d); `DiagonalSettings::default`
names `Off`, so `HelmholtzOperator`'s default before `with_diagonal` is unchanged and the operator
tests keep theirs. With named degrees the default takes no level (no tolerance). The rule is
CH4.3b1–CH4.3b2a's unchanged: on a level in range (k a ≥ 2.8) the bound's classes at the level's
degree; an RCR level takes the form where `DiagonalCosts` puts it below RCR; a level with Auto's
waves keeps its resolving degree p_w and keeps the waves unless the form beats them, then the
cheaper of the form and RCR; the form's scratch within `DEFAULT_SIGNATURE_NODE_BUDGET` per node,
else RCR whole. Its choices depend on the machine class of the costs (`diagonal_costs(isa)`), as
Batched's crossover and the near rules' do. The tests' `default_builder` names `Off` (the settings
before); `check_defaults` checks the new default against the direct sum, the same settings named
bit for bit, 2 threads and the one-rank FMM, and prints its levels (on the M6 the sphere at
kD = 30, 10⁻³: level 2, all 34 classes). `tetraleaf-validate`'s `Case::new` keeps `Off` named, so
the earlier reports' runs repeat.

## 3. Matched accuracy

The default's errors over the defaults' (relative ℓ2 at §9.1's samples, the same launch):

| machine | layout | cases | φ ÷ defaults | ∇φ ÷ defaults |
| --- | --- | ---: | --- | --- |
| Kathleen | 1 × 40 | 21 | 0.991–1.026 | 0.995–1.034 |
| locust | 1 × 72 | 21 | 1.000–1.019 | 0.997–1.016 |
| Kathleen | 40 × 1 | 5 | 1.000–1.026 | 0.997–1.034 |
| locust | 72 × 1 | 4 | 1.000–1.019 | 1.000–1.011 |

Above 1.01 only at 10⁻⁶ (W1 at kD = 100 and 500 on both machines), where the waves at iprec 2 and
the form both err near ε; within 1.1× everywhere.

## 4. The evaluation and the build

Per case, the default (`auto`) against the defaults (`off`): the errors, the evaluations of the
`run` launch, the ratio with the defaults first (`run`) and with the default first (`reverse`)
and their geometric mean, the build, and the levels whose M2L changed.

| machine | case | φ ÷ | ∇φ ÷ | evaluate off s | evaluate auto s | ÷, off first | ÷, auto first | mean | build s off → auto | M2L changed (level: off → auto) |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- |
| Kathleen | w1-100-1e-3 | 1.000 | 1.000 | 0.109 | 0.095 | 0.86 | 0.76 | 0.81 | 0.18 → 0.10 | 2:waves→form, 3:waves→form, 4:rcr→form |
| Kathleen | w1-100-1e-6 | 1.006 | 1.022 | 0.193 | 0.117 | 0.61 | 0.57 | 0.59 | 0.10 → 0.11 | 2:waves→form, 3:waves→form, 4:waves→form |
| Kathleen | w1-100-1e-7 | 1.000 | 1.000 | 0.189 | 0.191 | 1.01 | 0.98 | 0.99 | 0.10 → 0.12 | 2:waves→form, 3:waves→form |
| Kathleen | w1-100-1e-8 | 1.000 | 1.000 | 0.268 | 0.257 | 0.96 | 0.99 | 0.97 | 0.11 → 0.11 | 2:waves→form |
| Kathleen | w1-1000-1e-6 | 1.000 | 1.000 | 37.888 | 14.277 | 0.38 | 0.36 | 0.37 | 12.96 → 15.95 | 2–4:rcr→form, 5–7:waves→form |
| Kathleen | w1-200-1e-3 | 1.000 | 1.000 | 0.499 | 0.310 | 0.62 | 0.72 | 0.67 | 0.33 → 0.36 | 2–4:waves→form, 5:rcr→form |
| Kathleen | w1-200-1e-6 | 1.000 | 1.007 | 0.711 | 0.483 | 0.68 | 0.61 | 0.65 | 0.35 → 0.42 | 2–5:waves→form |
| Kathleen | w1-200-1e-7 | 1.000 | 1.000 | 0.927 | 0.652 | 0.70 | 0.78 | 0.74 | 0.37 → 0.44 | 2–4:waves→form |
| Kathleen | w1-200-1e-8 | 1.000 | 1.000 | 1.140 | 0.954 | 0.84 | 0.88 | 0.86 | 0.37 → 0.43 | 2–3:waves→form |
| Kathleen | w1-500-1e-3 | 1.000 | 1.000 | 4.521 | 2.316 | 0.51 | 0.43 | 0.47 | 2.49 → 3.04 | 2–3:rcr→form, 4–6:waves→form |
| Kathleen | w1-500-1e-6 | 1.026 | 1.034 | 5.828 | 3.100 | 0.53 | 0.46 | 0.49 | 2.49 → 3.25 | 2–3:rcr→form, 4–6:waves→form |
| Kathleen | w1-500-1e-7 | 1.000 | 1.000 | 7.755 | 4.495 | 0.58 | 0.52 | 0.55 | 2.61 → 3.32 | 2–3:rcr→form, 4–5:waves→form, 6:waves→rcr |
| Kathleen | w1-500-1e-8 | 1.000 | 0.997 | 8.879 | 6.058 | 0.68 | 0.64 | 0.66 | 2.57 → 3.42 | 2–3:rcr→form, 4–5:waves→form |
| Kathleen | w3-200-1e-3 | 1.000 | 1.000 | 1.805 | 1.941 | 1.08 | 0.75 | 0.90 | 0.86 → 0.91 | 2–3:waves→form |
| Kathleen | w3-200-1e-6 | 0.991 | 0.996 | 2.143 | 2.346 | 1.09 | 0.68 | 0.86 | 0.81 → 0.93 | 2–3, 5:waves→form |
| Kathleen | w3-200-1e-7 | 1.000 | 0.995 | 2.418 | 2.629 | 1.09 | 0.90 | 0.99 | 0.84 → 0.94 | 2–3:waves→form |
| Kathleen | w3-200-1e-8 | 0.993 | 0.996 | 2.980 | 3.156 | 1.06 | 0.90 | 0.97 | 0.84 → 0.94 | 2–3:waves→form |
| Kathleen | w3-500-1e-3 | 1.000 | 1.000 | 10.826 | 7.256 | 0.67 | 0.67 | 0.67 | 16.68 → 17.72 | 2–3:rcr→form, 4:waves→form |
| Kathleen | w3-500-1e-6 | 1.000 | 1.000 | 14.975 | 9.582 | 0.64 | 0.65 | 0.64 | 17.25 → 17.06 | 2–3:rcr→form, 4:waves→form |
| Kathleen | w3-500-1e-7 | 1.000 | 1.000 | 21.017 | 14.700 | 0.70 | 0.67 | 0.68 | 17.14 → 17.65 | 2–3:rcr→form, 4:waves→form |
| Kathleen | w3-500-1e-8 | 1.000 | 1.000 | 23.715 | 17.185 | 0.72 | 0.72 | 0.72 | 17.48 → 17.65 | 2–3:rcr→form |
| locust | w1-100-1e-3 | 1.000 | 1.000 | 0.047 | 0.035 | 0.74 | 0.75 | 0.75 | 0.07 → 0.06 | 2–4:waves→form |
| locust | w1-100-1e-6 | 1.007 | 1.016 | 0.055 | 0.042 | 0.77 | 0.80 | 0.78 | 0.07 → 0.07 | 2–4:waves→form |
| locust | w1-100-1e-7 | 1.000 | 1.000 | 0.092 | 0.088 | 0.96 | 0.97 | 0.96 | 0.07 → 0.07 | 2:waves→form |
| locust | w1-100-1e-8 | 1.000 | 1.000 | 0.111 | 0.108 | 0.98 | 0.96 | 0.97 | 0.06 → 0.06 | 2:waves→form |
| locust | w1-1000-1e-6 | 1.000 | 0.997 | 12.576 | 5.555 | 0.44 | 0.44 | 0.44 | 13.05 → 13.82 | 2–4:rcr→form, 5–8:waves→form |
| locust | w1-200-1e-3 | 1.000 | 1.000 | 0.181 | 0.128 | 0.71 | 0.71 | 0.71 | 0.26 → 0.27 | 2–5:waves→form |
| locust | w1-200-1e-6 | 1.009 | 1.007 | 0.214 | 0.156 | 0.73 | 0.73 | 0.73 | 0.26 → 0.28 | 2–5:waves→form |
| locust | w1-200-1e-7 | 1.000 | 1.000 | 0.347 | 0.326 | 0.94 | 0.94 | 0.94 | 0.28 → 0.29 | 2–3:waves→form |
| locust | w1-200-1e-8 | 1.000 | 1.000 | 0.439 | 0.423 | 0.96 | 0.95 | 0.96 | 0.22 → 0.24 | 2–3:waves→form |
| locust | w1-500-1e-3 | 1.000 | 1.000 | 1.256 | 0.779 | 0.62 | 0.62 | 0.62 | 2.64 → 2.73 | 2–3:rcr→form, 4–6:waves→form, 7:rcr→form |
| locust | w1-500-1e-6 | 1.019 | 1.011 | 1.965 | 1.207 | 0.61 | 0.59 | 0.60 | 2.68 → 2.85 | 2–3:rcr→form, 4–7:waves→form |
| locust | w1-500-1e-7 | 1.000 | 1.000 | 2.536 | 1.717 | 0.68 | 0.66 | 0.67 | 2.71 → 2.83 | 2–3:rcr→form, 4–5:waves→form |
| locust | w1-500-1e-8 | 1.000 | 1.000 | 2.937 | 2.180 | 0.74 | 0.73 | 0.74 | 1.59 → 1.78 | 2–3:rcr→form, 4:waves→form |
| locust | w3-200-1e-3 | 1.000 | 1.000 | 0.585 | 0.579 | 0.99 | 0.99 | 0.99 | 1.27 → 1.25 | 2:waves→form |
| locust | w3-200-1e-6 | 1.000 | 1.000 | 0.695 | 0.679 | 0.98 | 0.96 | 0.97 | 1.28 → 1.31 | 2–3:waves→form |
| locust | w3-200-1e-7 | 1.000 | 1.000 | 1.437 | 1.385 | 0.96 | 0.98 | 0.97 | 1.36 → 1.38 | 2–3:waves→form |
| locust | w3-200-1e-8 | 1.000 | 1.000 | 1.494 | 1.447 | 0.97 | 0.99 | 0.98 | 0.56 → 0.58 | 2:waves→form |
| locust | w3-500-1e-3 | 1.000 | 1.000 | 2.298 | 1.659 | 0.72 | 0.72 | 0.72 | 9.87 → 9.70 | 2–3:rcr→form, 6:waves→rcr |
| locust | w3-500-1e-6 | 1.000 | 1.000 | 3.265 | 2.088 | 0.64 | 0.64 | 0.64 | 9.99 → 9.81 | 2–3:rcr→form, 6:waves→rcr |
| locust | w3-500-1e-7 | 1.000 | 1.000 | 4.332 | 3.013 | 0.70 | 0.70 | 0.70 | 9.62 → 9.93 | 2–3:rcr→form |
| locust | w3-500-1e-8 | 1.000 | 1.000 | 4.948 | 3.528 | 0.71 | 0.71 | 0.71 | 9.71 → 9.97 | 2–3:rcr→form |

The build's diagonal tables (the build's part "diagonal tables", the interpolated build) cost
0.03–0.63 s on Kathleen at 1 × 40 (2.0 s on W1 at kD = 1000) and 0.01–0.15 s on locust (0.61 s);
the plane-wave tables of the levels taken are no longer built. Build plus one evaluation:
0.59–1.11× on Kathleen (above 1.0 W3 at kD = 200 in the off-first order, §6, and W1 at kD = 100,
10⁻⁷, 1.08× of 0.29 s) and 0.76–1.01× on locust.

## 5. The order per level

Per level with V pairs in the rule's range (k a ≥ 2.8) and not Batched's, the options timed in the
same `run` launch at the level's degree as run, scored by the option the default took as timed
there (CH4.3b2a §5.1):

| machine | ε | levels | within 1.05× | seconds lost | seconds at the best | taken from the waves: form, RCR |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| Kathleen | 10⁻³ | 19 | 18 | 0.13 | 5.65 | 11, 0 |
| Kathleen | 10⁻⁶ | 26 | 26 | 0.00 | 13.58 | 17, 0 |
| Kathleen | 10⁻⁷ | 17 | 16 | 0.47 | 10.34 | 10, 1 |
| Kathleen | 10⁻⁸ | 16 | 15 | 0.86 | 12.13 | 7, 0 |
| Kathleen | all | 78 | 75 | 1.46 | 41.71 | 45, 1 |
| locust | 10⁻³ | 22 | 18 | 0.15 | 1.43 | 11, 1 |
| locust | 10⁻⁶ | 29 | 29 | 0.00 | 4.45 | 17, 1 |
| locust | 10⁻⁷ | 20 | 18 | 0.07 | 3.52 | 7, 0 |
| locust | 10⁻⁸ | 19 | 17 | 0.05 | 3.09 | 5, 0 |
| locust | all | 90 | 82 | 0.28 | 12.48 | 40, 2 |

The misses, every one the waves kept where the form was faster at p_w:

| machine | case | level | p | k a | waves s | form s | RCR s | ÷ best |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Kathleen | W3, kD = 200, 10⁻³ | 4 | 17 | 10.83 | 0.338 | 0.211 | 0.752 | 1.61 |
| Kathleen | W3, kD = 500, 10⁻⁷ | 5 | 28 | 13.53 | 6.171 | 5.573 | 23.673 | 1.11 |
| Kathleen | W3, kD = 500, 10⁻⁸ | 4 | 44 | 27.06 | 2.516 | 1.655 | 9.641 | 1.52 |
| locust | W1, kD = 200, 10⁻⁷ | 4 | 28 | 10.83 | 0.051 | 0.033 | 0.050 | 1.54 |
| locust | W1, kD = 500, 10⁻⁸ | 5 | 32 | 13.53 | 0.186 | 0.167 | 0.295 | 1.12 |
| locust | W3, kD = 200, 10⁻³ | 3–5 | 28–11 | 21.65–5.41 | | | | 1.38–2.12 |
| locust | W3, kD = 200, 10⁻⁸ | 3 | 38 | 21.65 | 0.057 | 0.031 | 0.130 | 1.86 |
| locust | W3, kD = 500, 10⁻³, 10⁻⁷ | 4 | 34–43 | 27.06 | | | | 1.13–1.23 |

These are CH4.3b1's and CH4.3b2a's misses (the waves' constant per box overestimates the form's
competitor on W3's dense levels, and two raised W1 levels where the form passes 27 classes at p_w):
3.5 % of the levels' best time on Kathleen and 2.2 % on locust. RCR takes a level from the waves
only where the model has it below the form (CH4.3b2a's order): W1 at kD = 500, 10⁻⁷, level 6 on
Kathleen, W3 at kD = 500, level 6 at 10⁻³ and 10⁻⁶ on locust, all within 1.05×.

**D8.** Under the default the form replaces Auto's waves on 45 (Kathleen) and 40 (locust) of the
levels where D8's rule gives them, at the level's resolving degree, and RCR on 1 and 2 more. D8's
row says `Auto` takes the exponential M2L on every level of its rule; with this default it takes
it on those levels unless the time rule moves the level (question 3).

## 6. W3 at kD = 200 on Kathleen: run order, not the form

CH4.3b2a saw the time rule 1.03–1.12× the defaults' evaluation on W3 at kD = 200 on Kathleen with
every choice within 1.05×; the `run` launches here repeat it (1.06–1.09×). Two tests separate the
setting from its place in the launch:

- **Order reversed** (`reverse`, the time rule first): 0.68–0.90× the defaults' (§4's table). The
  geometric mean of the two orders is 0.86–0.99×. On locust the two orders agree within 0.03.
- **Same setting twice** (`repeat`, `off, auto, off, auto` in one launch, one build each):

| machine | ε | off (1st) s | auto (2nd) s | off (3rd) s | auto (4th) s | 3rd ÷ 1st |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Kathleen | 10⁻³ | 1.774 | 1.931 | 1.926 | 1.970 | 1.086 |
| Kathleen | 10⁻⁶ | 2.126 | 2.275 | 2.556 | 2.005 | 1.202 |
| Kathleen | 10⁻⁷ | 2.436 | 2.637 | 2.865 | 2.571 | 1.176 |
| Kathleen | 10⁻⁸ | 2.941 | 3.234 | 3.305 | 3.472 | 1.124 |
| locust | 10⁻³ | 0.588 | 0.577 | 0.583 | 0.579 | 0.991 |
| locust | 10⁻⁶ | 0.696 | 0.677 | 0.700 | 0.678 | 1.006 |
| locust | 10⁻⁷ | 1.417 | 1.399 | 1.418 | 1.389 | 1.000 |
| locust | 10⁻⁸ | 1.500 | 1.451 | 1.469 | 1.444 | 0.979 |

The same computation (`off`) runs 9–20 % slower the second time in a launch on Kathleen, and the
waves' level 4 within it 13–22 % slower (0.314 → 0.364 s at 10⁻³, 0.383 → 0.469 s at 10⁻⁶,
0.682 → 0.818 s at 10⁻⁷, 0.741 → 0.838 s at 10⁻⁸), while the form's levels 2 and 3 run alike in
both `auto` runs: the drift is the launch's on Kathleen's two-socket node, not a placement effect of
the form's buffers (on W1 the time rule's ratio is the lower one when it runs first in 9 of 13
cases, at most 0.10 apart). CH4.3b2b's "later runs slower"
is the same effect. Kathleen reports of evaluation ratios should alternate the order, as this one
and CH4.3b2d's did.

## 7. The whole-level fallback at 40 × 1

W1 at kD = 1000, 10⁻⁶, 40 ranks of one thread (`ranks-w1-1000-1e-6`, Kathleen job 663644): the
form's scratch per rank against the rank's budget (8 GiB per node over 40 ranks, 205 MiB), and
the M2L per level:

| level | p | k a | scratch MiB (most) | ranks over | model's choice | ran (off) | s | ran (auto) | s |
| ---: | ---: | ---: | ---: | ---: | --- | --- | ---: | --- | ---: |
| 2 | 240 | 216.5 | (none fits) | 40 | diagonal | RCR by pairs | 35.05 | RCR by pairs | 34.79 |
| 3 | 128 | 108.3 | 204 | 7 | diagonal | RCR by pairs | 10.45 | RCR by pairs | 10.42 |
| 4 | 71 | 54.1 | 137 | 0 | diagonal | RCR by pairs | 4.99 | diagonal 34/34 | 0.96 |
| 5 | 41 | 27.1 | | 0 | diagonal | waves | 5.38 | diagonal 34/34 | 0.84 |
| 6 | 25 | 13.5 | | 0 | diagonal | waves | 2.81 | diagonal 34/34 | 0.94 |
| 7 | 16 | 6.8 | | 0 | diagonal | waves | 3.29 | diagonal 34/34 | 0.95 |

Levels 2 and 3 run RCR whole on every rank (no class taken on any rank, the same on all), as CH4.1b
predicted and CH4.2b1 built; the model would take the form there. Errors 1.000× (φ and ∇φ); the
evaluation 53.0 s against 63.7 s (0.83×), of which the two fallback levels are 45.2 s. At 1 × 40
the same levels take the form and the evaluation is 14.3 s: one rank per node (or the signatures
per node, CH8.0) is what lets the coarse high-frequency levels use the form.

## 8. Many ranks per node

The time rule's constants are fitted at one rank per node (CH4.3b1). At one rank per core:

| machine | layout | case | evaluate off s | evaluate auto s | ÷ | build s off → auto (diagonal tables) |
| --- | --- | --- | ---: | ---: | ---: | --- |
| Kathleen | 40 × 1 | W1, kD = 500, 10⁻⁶ | 10.09 | 3.81 | 0.38 | 4.85 → 14.48 (10.05) |
| Kathleen | 40 × 1 | W1, kD = 500, 10⁻⁸ | 14.52 | 7.09 | 0.49 | 6.36 → 16.28 (10.92) |
| Kathleen | 40 × 1 | W3, kD = 500, 10⁻⁶ | 29.33 | 16.44 | 0.56 | 5.31 → 14.97 (9.89) |
| Kathleen | 40 × 1 | W3, kD = 500, 10⁻⁸ | 55.73 | 47.41 | 0.85 | 6.39 → 16.14 (9.79) |
| Kathleen | 40 × 1 | W1, kD = 1000, 10⁻⁶ | 63.74 | 53.04 | 0.83 | 42.34 → 45.66 (3.27) |
| locust | 72 × 1 | W1, kD = 500, 10⁻⁶ | 4.34 | 3.25 | 0.75 | 2.21 → 2.99 (1.05) |
| locust | 72 × 1 | W1, kD = 500, 10⁻⁸ | 6.29 | 5.16 | 0.82 | 2.64 → 3.43 (1.04) |
| locust | 72 × 1 | W3, kD = 500, 10⁻⁶ | 10.62 | 10.58 | 1.00 | 2.31 → 2.19 (0) |
| locust | 72 × 1 | W3, kD = 500, 10⁻⁸ | 18.77 | 18.64 | 0.99 | 2.78 → 2.64 (0) |

Errors 0.997–1.034× (§3). The choices (with the three costs, the first four rows per machine): 16
of 17 levels within 1.05× on Kathleen (17.97 s lost of 41.47 s), 19 of 20 on locust (0.20 s of
35.06 s). The two misses:

- **Kathleen, W3 at kD = 500, 10⁻⁸, level 4** (p_w = 44, k a 27.1): the waves 20.54 s, the form
  2.57 s, RCR 14.86 s; the default keeps the waves (8.0× the best). At 1 × 40 the same level's
  waves took 2.52 s: the waves pass 19.0 s per rank against 0.83 s (merge and shift 1.57 against
  1.69 s), at about the same chunks (6 against 4). The one-rank constant per box predicts the
  waves below the form there; at 40 ranks per node the waves on W3's dense levels cost about 8×
  as much, under the defaults too (so the defaults themselves lose 18 s there against the form).
  Not diagnosed here (each rank computes the waves of its distinct V sources, ghosts included; a
  rank's share of a dense level's sources is much more than its own boxes).
- **locust, W1 at kD = 500, 10⁻⁸, level 5** (p_w = 32): the waves 0.46 s, the form 0.26 s, 1.77×;
  the raised level of §5's misses.

Every rank builds every diagonal table serially (the interpolated build by ring, per rank): 9.8–
10.9 s at 40 × 1 on Kathleen, 1.0 s at 72 × 1 on locust. Build plus one evaluation at 40 × 1:
0.91–1.22× the defaults' (above 1 on W1 at kD = 500: 1.22× at 10⁻⁶, 1.12× at 10⁻⁸); from the
second evaluation on, the default wins on every case. The documented advice stays a few ranks per
node with threads (CH4.0's dense budget, CH4.1b); one copy of the tables per node is CH8.0.

## 9. Questions for sign-off

1. **The acceptance as read** (the plan's reading of the roadmap line): the errors within 1.1× on
   every case of W1 and W3 at 10⁻³–10⁻⁸ on Kathleen and locust (0.991–1.034×, §3); the evaluation
   and the build per level (§4, §7, §8); the order per level (§5); the whole-level fallback at
   40 × 1 on W1 at kD = 1000 (levels 2 and 3, §7); the Alsos: the dense-level misses listed (§5),
   the constants at many ranks per node measured (§8), W3 at kD = 200 repeated and separated (§6).
   *Recommended: accept.*
2. **`Auto` the default** (§2): `HelmholtzFmmBuilder`'s and the enum's default
   `HelmholtzDiagonal::Auto`, the operator's own default `Off`. *Recommended: adopt*: the
   evaluation 0.37–0.99× at matched accuracy on every case at one rank per node on both machines,
   build plus one evaluation at most 1.11×. *Alternative:* keep `Off` until question 4's task.
3. **D8 amended** (§5): under the default, a level that D8's rule gives the exponential M2L keeps
   it unless the time rule moves it (the form, where it beats the waves at the level's resolving
   degree, then the cheaper of the form and RCR). *Recommended:* amend D8's row and the proposal's
   §5.3 (the D8 bullet) and §6.4 with a line in Appendix A, in this PR.
4. **Many ranks per node and the dense levels** (§5, §8): a follow-up task, **CH4.3f The time
   rule's costs at many ranks per node**: *Done when:* the waves' cost per box at 40 × 1 on W3's
   dense levels explained (8× the one-rank constant at kD = 500, 10⁻⁸, level 4) and the costs
   fitted per ranks per node (or the rule's inputs given the rank's share), the dense-level misses
   at one rank (3.5 % and 2.2 % of the levels' time) refitted or closed, verified as CH4.3b2a on
   Kathleen and locust at 1 × N and N × 1. *Recommended:* add CH4.3f off the priority list,
   reconsidered at CH4.5 (which compares ours at 1 × N with the best split), with an Also on CH8.0
   (the diagonal tables per rank, about 10 s at 40 × 1). *Alternative:* on the list after CH4.3c3,
   before CH4.4.
5. **W3 at kD = 200 on Kathleen** (§6): closed as Kathleen's run-order drift (the same computation
   9–20 % slower the second time in a launch), not the form; no task. *Recommended: accept.*
