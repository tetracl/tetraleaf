# Helmholtz CH4.4b: the crosscheck mapped to levels and classes

Task plan and record: `docs/helmholtz/tasks/CH4.4b.md`. Proposal: §10 CH4.4 (the part "the
full-target crosscheck mapped to levels and classes, with no error concentrated at a class
boundary", which CH4.4a left), §9.1 (the box measure and its 100 ε bound, amended at CH4.0b), §3.6
(the rule per (level, class)). Machine: Kathleen, one node per job (2 × 20 Cascade Lake cores, Xeon
Gold 6248), `singlenode`, one rank × 40 threads (the box study runs on one rank), release, BLAS
threads 1, node-c11b-004 and node-c11b-013 excluded. Jobs 663781 (`low`), 663782 (`mid`), 663783
(`w3`), 663784 (`w1-1000`), 663788 (`pw1000`, `remedy`), 663789 (`w3-cells`); 663774 and 663780
the builds and checks.
Raw outputs: `~/projects/fmm-archive/kathleen/runs/ch4-4b/` (one Markdown file per case and ε, its
log beside it); the first runs without `--refine-faces` in `~/projects/fmm-archive/kathleen/runs/ch4-4b-superseded/` (§1).

## Summary

- **The mapping** (W1 to kD = 1000, W2 to 500, W3 resolved to 100 and at N = 10⁶ for kD = 200 and
  500, with the plane wave's trace and the constant density; 10⁻³, 10⁻⁶, 10⁻⁸; 47 launches on
  Kathleen under the library's defaults): the box errors are the levels' truncation (the levels'
  parts add up to within 0.4 % of Σ|e|² at 10⁻⁶ and 10⁻⁸, 4 % at 10⁻³),
  the finer levels' at 10⁻³ and the coarsest's at 10⁻⁸, spread over three to six levels at the worst
  boxes, which on W1 and W2 are the six tangent points' (CH4.0b's placement).
- **No error concentrated at a class boundary.** The diagonal form on RCR changes no worst box and
  no ℓ2 error (154 runs, 1.000, ℓ2 of φ within 1.001; at most 0.08 % of Σ|e|²); the classes beside the bound's threshold
  and the next decade on RCR the same (55 runs); Batched on RCR the same (11). Where the waves run,
  RCR at the same degree errs up to 1.19× more (30 runs; the default at most 1.023× RCR in the worst
  box, 1.034× in ℓ2), which is the waves' accuracy, not an error they add.
- **CH4.4a's 114.7 ε box** (W1, kD = 1000, the plane wave's trace, 10⁻⁸; 114.80 ε here, the same
  box) is the truncation of levels 2–5 at a tangent point (0.14–0.29 of its error each); its
  classes and the form add nothing. One level raised by two degrees leaves 89.4 ε; every level by
  one 48.4 ε (about 2–4 % more evaluation), by two 21.7 ε. Per pair the rule's degrees are 3.3–5.5 ε on
  every level, no outlier. Within D1's range (kD ≤ 500) every worst box is at most 68.2 ε.
  Recommended (question 2): §9.1's 100 ε bound for D1's range, beyond it reported; the
  recalibration is the measured alternative.
- **The tool.** `M2lOverride` (`tetraleaf-exec`, a measurement setting, empty by default; tested
  against the time rule's own RCR and against named degrees, bit for bit) and `helmholtz_boxes`'s
  studies under the defaults (`--study levels,m2l`, `--refine-faces`, `--raise-by`). The sampled
  reference reproduces CH4.4a's worst boxes on W1 to 0.8 % (§2).

## 1. Method

- **The cases.** CH4.4a's cases whose trees hold a level of k a ≥ 2.8 with the diagonal form or RCR
  under the defaults (its logs' "by level"): W1 at kD = 30 (10⁻⁶ only: Batched at 10⁻³, the waves at
  10⁻⁸), 100, 250, 500, 1000; the plane wave's trace at kD = 500 and 1000; the constant density on W1
  and W2 at kD = 30 (10⁻⁶) and 100; W2 at kD = 30 (10⁻⁶), 100, 250, 500; W3 resolved at kD = 30
  (10⁻⁶) and 100; W3 at N = 10⁶ for kD = 200 and 500. ε = 10⁻³, 10⁻⁶, 10⁻⁸; one complex charge
  vector, gradients on; 47 launches (one case and ε each; W3's 9 again with CH4.4a's worst boxes,
  below) and 4 remedy launches, every one exit 0.
- **The FMM.** The library's defaults (`Case::defaults`, CH4.4a's: the waves under Auto, the
  diagonal form's time rule, the refit degrees, c_leaf 4, the V, W and X rules). Every case's tree
  and each level's degree and M2L are CH4.4a's (its logs), at one rank × 40 threads against its
  4 × 10.
- **The reference.** `helmholtz_boxes`: the f64 direct sum (`tetraleaf_ref::helmholtz::p2p::
  direct_sum` × 1/(4π)) at every point to N = 2·10⁵; above at §9.1's 1,000 samples, two targets
  per box of the deepest level, and every point of the 8 worst boxes of φ and of ∇φ by that sample
  and of the boxes at the centres of the domain's faces, W1's six tangent points (`--refine 8
  --refine-faces`; up to 1.54·10⁵ targets at kD = 1000). The first runs had `--refine` alone: two
  targets per box missed CH4.4a's worst box at kD = 500 (29.1 against 40.4 ε), so every stage ran
  again with the faces (the first runs kept in `ch4-4b-superseded/`). The box measure (§9.1, CH3.6's): per box of the deepest
  level the largest error at one of its reference targets over the reference's RMS.
- **The levels** (CH4.0b's parts, under the defaults). Each level with V pairs raised to the rule's
  degree at ε/1000 (at least one above its degree as run) and run by RCR, every other level as the
  defaults choose, by a new measurement setting (`M2lOverride`, below); level l's part is the
  default's error minus that run's, its share Σ|part|² over Σ|error|² over the reference's targets;
  the rest what the parts leave. CH4.0b decomposed on all-Rotation, since named degrees change the
  waves' precision and take the form off; the override keeps the tolerance, the rules and the plan.
- **The M2L.** Each level whose M2L under the defaults is not RCR (the form, the waves, Batched) run
  by RCR at its degree as run (`OverrideM2l::Rcr`; the waves' resolving degree p_w where raised).
  "No error concentrated where the M2L changes" is read as: the default's worst box and ℓ2 at most
  1.1× that run's (`DIAGONAL_RATIO`, `PLANE_WAVE_RATIO`), one-sided (question 1).
- **The classes.** On each level of the form, the classes that pass the bound by decade of
  E(l, c)/(θ ε) below the threshold ((0.1, 1], "beside the threshold"; (0.01, 0.1]), each decade by
  RCR in one run (`OverrideM2l::RcrClasses`); on the remedy case also each class of the first decade
  alone. The classes that fail the bound run RCR in the default already.
- **The override** (`tetraleaf-exec`, a measurement setting, empty by default): `M2lOverride {
  level, degree, m2l }`, `HelmholtzFmmBuilder::m2l_overrides`. A degree replaces the rule's before
  every rule reads it (the builder's degrees, so the plan, the near rules and the waves' levels see
  it); `Rcr` takes the waves, the form and Batched's tables off the level; `RcrClasses` takes classes
  off the form after the bound and the time rule chose them. Tested against trusted paths
  (`tests/helmholtz/m2l_override.rs`; `check_m2l_override` in `tests/helmholtz_fmm.rs`: RCR on the
  form's level bit for bit the time rule's own RCR at free RCR costs, a named degree under Rotation
  bit for bit `Degrees::PerLevel`, 2 threads and the one-rank FMM bit for bit).
- **The remedy** (W1 with the plane wave's trace at kD = 500 and 1000, 10⁻⁸). FMM evidence: the level
  with the largest part at the worst ∇φ box raised by 1 and 2 under the defaults' M2L
  (`OverrideM2l::Setting`, the form at the new degree), and every level with V pairs raised at once
  by 1, 2 and 4. Pair evidence: spike S7's worst pair over the 316 offsets and three placements at
  each level's k a and degree (`spikes/truncation/results/c1_curves.csv`, read as CH4.3c2's pair
  evidence; `tools/helmholtz/ch4-4b-tables.py pairs`).

## 2. The reference against CH4.4a's

The base run of each launch is CH4.4a's FMM (the same tree, degrees and M2L per level; §9.1's
errors at the samples the same to every printed digit), measured against a sampled and refined
direct sum in place of FMM3D at 10⁻¹² at every point. The worst box over ε, φ / ∇φ:

| case | ε | CH4.4a φ / ∇φ | here φ / ∇φ |
| --- | --- | ---: | ---: |
| W1, kD = 30 | 1e-6 | 3.34 / 8.26 | 3.34 / 8.26 |
| W1, kD = 100 | 1e-3 | 5.43 / 7.20 | 5.43 / 7.20 |
| W1, kD = 100 | 1e-6 | 4.51 / 14.72 | 4.51 / 14.72 |
| W1, kD = 100 | 1e-8 | 5.04 / 8.67 | 5.04 / 8.67 |
| W1, kD = 250 | 1e-3 | 13.23 / 13.57 | 13.24 / 13.54 |
| W1, kD = 250 | 1e-6 | 8.85 / 16.98 | 8.84 / 16.90 |
| W1, kD = 250 | 1e-8 | 16.12 / 21.55 | 16.09 / 21.45 |
| W1, kD = 500 | 1e-3 | 21.28 / 22.99 | 21.43 / 23.01 |
| W1, kD = 500 | 1e-6 | 14.22 / 20.01 | 14.32 / 20.03 |
| W1, kD = 500 | 1e-8 | 27.38 / 40.42 | 27.59 / 40.47 |
| W1, kD = 1000 | 1e-3 | 19.96 / 20.26 | 20.03 / 20.25 |
| W1, kD = 1000 | 1e-6 | 25.94 / 28.00 | 26.05 / 27.99 |
| W1, kD = 1000 | 1e-8 | 41.84 / 43.31 | 42.00 / 43.30 |
| W1, constant, kD = 30 | 1e-6 | 0.76 / 5.45 | 0.76 / 5.45 |
| W1, constant, kD = 100 | 1e-3 | 2.92 / 8.24 | 2.92 / 8.24 |
| W1, constant, kD = 100 | 1e-6 | 5.55 / 15.02 | 5.55 / 15.02 |
| W1, constant, kD = 100 | 1e-8 | 22.48 / 34.85 | 22.48 / 34.85 |
| W1, planewave, kD = 500 | 1e-3 | 27.92 / 39.87 | 27.72 / 39.61 |
| W1, planewave, kD = 500 | 1e-6 | 26.02 / 42.19 | 26.14 / 42.42 |
| W1, planewave, kD = 500 | 1e-8 | 35.43 / 67.77 | 35.63 / 68.21 |
| W1, planewave, kD = 1000 | 1e-3 | 35.59 / 43.35 | 35.56 / 43.34 |
| W1, planewave, kD = 1000 | 1e-6 | 46.68 / 73.86 | 46.66 / 73.88 |
| W1, planewave, kD = 1000 | 1e-8 | 69.49 / 114.74 | 69.49 / 114.80 |
| W2, kD = 30 | 1e-6 | 2.34 / 7.99 | 2.34 / 7.99 |
| W2, kD = 100 | 1e-3 | 3.18 / 2.47 | 3.18 / 2.47 |
| W2, kD = 100 | 1e-6 | 10.59 / 11.31 | 10.59 / 11.31 |
| W2, kD = 100 | 1e-8 | 3.32 / 2.88 | 3.32 / 2.88 |
| W2, kD = 250 | 1e-3 | 10.11 / 5.36 | 15.91 / 10.84 |
| W2, kD = 250 | 1e-6 | 11.93 / 15.11 | 18.79 / 23.81 |
| W2, kD = 250 | 1e-8 | 10.84 / 8.47 | 17.08 / 17.18 |
| W2, kD = 500 | 1e-3 | 13.44 / 10.78 | 21.78 / 22.13 |
| W2, kD = 500 | 1e-6 | 11.57 / 15.72 | 18.76 / 32.34 |
| W2, kD = 500 | 1e-8 | 13.81 / 13.74 | 22.39 / 28.27 |
| W2, constant, kD = 30 | 1e-6 | 0.44 / 7.53 | 0.44 / 7.53 |
| W2, constant, kD = 100 | 1e-3 | 2.75 / 24.24 | 2.75 / 24.24 |
| W2, constant, kD = 100 | 1e-6 | 3.18 / 34.68 | 3.18 / 34.68 |
| W2, constant, kD = 100 | 1e-8 | 6.56 / 37.55 | 6.56 / 37.55 |
| W3 (N = 108849), kD = 30 | 1e-6 | 4.10 / 3.88 | 4.10 / 3.88 |
| W3 (N = 4031442), kD = 100 | 1e-3 | 11.05 / 6.34 | 11.02 / 3.85 |
| W3 (N = 4031442), kD = 100 | 1e-6 | 13.54 / 9.64 | 13.55 / 5.83 |
| W3 (N = 4031442), kD = 100 | 1e-8 | 9.72 / 7.97 | 9.72 / 4.83 |
| W3 (N = 1000000), kD = 200 | 1e-3 | 10.82 / 7.07 | 10.84 / 12.01 |
| W3 (N = 1000000), kD = 200 | 1e-6 | 10.21 / 10.58 | 10.23 / 17.97 |
| W3 (N = 1000000), kD = 200 | 1e-8 | 9.42 / 7.67 | 9.44 / 13.03 |
| W3 (N = 1000000), kD = 500 | 1e-3 | 16.43 / 15.86 | 16.45 / 19.03 |
| W3 (N = 1000000), kD = 500 | 1e-6 | 15.40 / 15.41 | 15.41 / 18.50 |
| W3 (N = 1000000), kD = 500 | 1e-8 | 24.59 / 24.76 | 24.61 / 29.71 |

- **W1** (and the constant density, the plane wave's trace): CH4.4a's figures to within 1 % (at
  most 0.8 % apart at kD ≥ 250, the same digits at kD ≤ 100; the 114.7 ε box 114.80 here against
  114.74, the same box, (127, 63, 64)); the worst boxes are the tangent points' (§1), which the
  reference holds whole. The differences are the reference's normaliser and FMM3D's own error at
  10⁻¹² against the direct sum here, and one rank against four (the same tree, degrees and M2L).
- **W2**: the same errors at the samples, and the worst ∇φ box CH4.4a's worst on 13 of 14
  launches (its second on the other, kD = 250 at 10⁻⁶; CH4.4a's logs' five worst), but larger figures at kD = 250
  and 500 (1.6–2.0×): the box measure divides by the reference's RMS over its own targets, and two
  targets per box weigh the sparse sphere more than the dense cap, whose |φ| and |∇φ| are larger. At
  kD ≤ 100 the reference holds every point and the figures are CH4.4a's.
- **W3**: without tangent points the sample's 8 worst boxes need not be the worst (the stage `w3`
  found lower figures at kD = 100: 7.5 against 11.1 ε in φ), so W3 ran again with CH4.4a's own five
  worst boxes of φ and of ∇φ per case and ε in the reference (`--refine-cells`, the stage
  `w3-cells`, job 663789; the rows above): φ's worst box is CH4.4a's to 0.3 % on every W3 case. ∇φ's
  differs by a factor fixed per case, the same at every ε (0.61 at kD = 100; 1.70 and 1.20 at
  N = 10⁶ for kD = 200 and 500): the normaliser, the reference's RMS of |∇φ| over the sample against
  over every point, which the cube's close pairs of points make heavy-tailed (N = 10⁶, kD = 200,
  10⁻⁸: CH4.4a's all-target ∇φ error 1.42·10⁻⁹ against 2.84·10⁻⁹ over this reference's targets).

So the gate's figures stay CH4.4a's (§9.1's box measure against every point); this task's measure
maps the error, and its shares and ratios (a run's error over the default's) do not depend on the
normaliser. On W1, where the bound binds, both measures agree.

## 3. The levels

Each level's share of Σ|e|² over the reference's targets is in §7 ("Levels"), with its M2L; the
level with the largest share of ∇φ's per launch in §7 ("Summary").

- **The parts add up.** The rest (what the levels' parts leave of the error) is at most 0.004 of
  Σ|e|² at 10⁻⁶ and 10⁻⁸ and at most 0.04 at 10⁻³ on every launch. On W3 the study raises the
  levels without V pairs too (`--all-levels`): at N = 10⁶, kD = 500, 10⁻⁸ the deepest level (p = 32,
  no V pairs: its translations and leaf operators) carries 0.06 (φ) and 0.12 (∇φ) of Σ|e|², which
  the first run's study, raising the levels with V pairs alone, left as a rest of 0.11 and 0.19.
  At 10⁻³ the rest (up to 0.04) is not explained by measurement; a likely part: a raised level's
  degree changes the plan's near pairs (the V, W and X crossovers read the degree, and the raised
  run's degree is the rule's at 10⁻⁶), so a raised run also converts pairs of other levels, which
  the decomposition does not separate. So the error is the levels' truncation, as CH4.0b found on
  all-Rotation trees, now under the defaults with the form, the waves and Batched on the levels.
- **Which levels.** At 10⁻³ the finer levels carry the most (W1: level 3 at kD = 100 and 250,
  level 6 at 500 and 1000, 19–63 % of ∇φ's Σ|e|²); at 10⁻⁸ the coarsest (W1: level 2 at
  kD = 100–1000, 30–55 %; level 3 with the plane wave's trace at kD = 500 and 1000; the constant density's
  level 2 at kD = 100, 68 %). At large kD it spreads: W1 at kD = 1000, 10⁻⁸, levels 2–6 carry 0.30,
  0.26, 0.11, 0.12 and 0.07 of ∇φ's. On W3 a middle level leads, the form's or the waves' (39–78 %),
  on W2 at 10⁻³ Batched's finest (levels 5–7, 24–38 %).
- **At the worst boxes.** Every worst box of W1 and W2 is one of the six tangent points' (CH4.0b's
  placement): its worst target lies near a corner of its box on every level at once (|v|/r 1.6–1.73
  on levels 2 to the deepest), so three to six levels' parts are each 0.11–0.69 of the error at the
  worst ∇φ target of W1 at kD ≥ 250 (§7, "Worst boxes"), never one level's alone. W3's worst
  boxes (CH4.4a's) lie inside the cube, and some are one level's (level 4 at N = 10⁶: 0.93–1.01 of
  the error at kD = 200, 10⁻³ and 10⁻⁸, and kD = 500, 10⁻⁶), the waves' or the form's, with the M2L's
  part there at most 0.023: that level's truncation, not its M2L.

## 4. The M2L and the classes

Runs with one level's M2L on RCR at its degree as run, over the default's (the worst box over all
boxes, the ℓ2 error over the reference's targets; a ratio above 1 is RCR erring more):

| the level's M2L under the defaults | runs | worst box φ | worst box ∇φ | ℓ2 φ | ℓ2 ∇φ | share of Σ|e|² φ, ∇φ |
| --- | ---: | --- | --- | --- | --- | --- |
| the diagonal form | 154 | 1.000–1.000 | 1.000–1.000 | 1.000–1.001 | 1.000–1.000 | ≤ 0.0008, ≤ 0.0001 |
| Batched | 11 | 1.000–1.000 | 1.000–1.000 | 1.000–1.000 | 1.000–1.000 | ≤ 0.0000, ≤ 0.0000 |
| the waves | 30 | 0.988–1.189 | 0.978–1.190 | 0.967–1.005 | 0.999–1.068 | ≤ 0.079, ≤ 0.033 |
| the form's classes by decade below the threshold, (0.1, 1] and (0.01, 0.1], each in one run | 55 | 1.000–1.000 | 1.000–1.000 | 1.000–1.001 | 1.000–1.000 | ≤ 0.0008, ≤ 0.0001 |

- **The form adds nothing measurable.** On every level the form runs (k a 5.4–217, 21 to 34
  classes, ε 10⁻³–10⁻⁸) RCR gives the same worst box to four digits and ℓ2 within 0.1 %: the form's
  rounding (E ≤ θ ε by the bound, CH4.3a's margin) is far below the level's truncation, which RCR
  and the form share. So no box's error comes from the form or from where it changes to RCR or the
  waves.
- **The classes beside the threshold** (E/(θ ε) in (0.1, 1]: classes 7 and 8 on W1's and W2's
  levels 4–6 at 10⁻⁸; classes 0 and 1, or 13–17, elsewhere) and the next decade: on RCR they change no
  worst box and no ℓ2 error (1.000). The bound's threshold is not a boundary of the error.
- **The waves.** RCR at the same degree errs up to 1.19× the waves' at the worst box (W2 at kD = 30,
  10⁻⁶, level 4) and as little as 0.978× (W2 at kD = 500, 10⁻⁸, level 7): the default at most 1.023×
  RCR in the worst box and 1.034× in ℓ2, within 1.1. The waves approximate the kernel at the FMM's
  tolerance, not RCR's truncation at p_w (CH6.1g), so they may err less than RCR there; that is no
  error they add.
- **Batched** is Dense's matrices, exact to rounding: RCR's output.

So, read with §9.1's box measure (question 1): **no error concentrates at a class boundary**: not
in the classes beside the bound's threshold, and not at the levels where the M2L changes (from the
form to RCR or the waves, from the waves to Batched); the box errors are the levels' truncation at
the tangent points, on every level at once.

## 5. The worst box of W1 at kD = 1000, the plane wave's trace, 10⁻⁸

CH4.4a's worst box, (127, 63, 64) on level 7, the sphere's tangent point on the cube's +x face, is
114.80 ε in ∇φ here (CH4.4a 114.74 ε), 69.49 ε in φ; the next worst 54.1 ε, the five next all at
the six tangent points.

- **Its levels.** At its worst target levels 2, 3, 4 and 5 carry 0.14, 0.29, 0.27 and 0.22 of the
  error (|v|/r 1.64–1.73 on levels 2–7: near a corner of its box on every level); over all targets
  levels 2–6 carry 0.22, 0.27, 0.17, 0.16 and 0.06 of ∇φ's Σ|e|². No level dominates.
- **Its classes.** Levels 2–6 run the form (34, 34, 34, 32 and 27 classes); every one of them on
  RCR, the decades below the threshold on RCR, and each of the two classes beside it (7 and 8 on
  level 6, E/(θ ε) in (0.1, 0.21]) alone on RCR: the worst box 114.80 ε in every run (1.000). The form
  and its classes add nothing; the box is the truncation of levels 2–5.
- **The same on the random charges and at kD = 500.** W1 with random charges at kD = 1000: 43.3 ε,
  levels 2–5 each 0.24–0.39 of the worst target's error; the plane wave's trace at kD = 500: 68.2 ε,
  levels 2–4 each 0.28–0.33. The plane wave's trace concentrates more (its far field is coherent at
  the tangent points, CH4.0b), and the measure grows with kD at the same ε (W1, 10⁻⁸: 8.7, 21.5,
  40.5 and 43.3 ε at kD = 100, 250, 500 and 1000 with random charges; 68.2 and 114.8 ε with the
  plane wave at 500 and 1000): the tangent point is a box corner on every level, and the tree has
  one level more each time kD doubles.

**The FMM evidence** (one launch each, Kathleen 1 × 40, the evaluations' seconds of the same launch
compared; the build changes too, by the tables of the new degrees):

| W1, the plane wave's trace, 10⁻⁸ | worst box φ, ∇φ ÷ ε | ℓ2 φ, ∇φ (the reference's targets) | evaluation s |
| --- | --- | --- | ---: |
| kD = 500, the defaults | 35.63, 68.21 | 6.18e-9, 1.14e-8 | 6.73 |
| kD = 500, level 3 at +1 | 29.88, 54.48 | 4.71e-9, 8.88e-9 | 6.93 |
| kD = 500, level 3 at +2 | 25.68, 49.07 | 4.45e-9, 8.65e-9 | 6.42 |
| kD = 500, levels 2–5 at +1 | 16.81, 25.02 | 2.30e-9, 4.14e-9 | 7.14 |
| kD = 500, levels 2–5 at +2 | 5.41, 10.36 | 9.64e-10, 1.72e-9 | 7.39 |
| kD = 500, levels 2–5 at +4 | 0.83, 1.73 | 1.78e-10, 3.11e-10 | 7.70 |
| kD = 1000, the defaults | 69.49, 114.80 | 5.43e-9, 9.17e-9 | 26.32 |
| kD = 1000, level 3 at +1 | 59.43, 97.10 | 4.76e-9, 8.12e-9 | 26.74 |
| kD = 1000, level 3 at +2 | 52.59, 89.42 | 4.57e-9, 7.99e-9 | 26.82 |
| kD = 1000, levels 2–6 at +1 | 29.41, 48.41 | 2.72e-9, 4.02e-9 | 27.35 |
| kD = 1000, levels 2–6 at +2 | 11.78, 21.72 | 1.32e-9, 2.03e-9 | 28.56 |
| kD = 1000, levels 2–6 at +4 | 2.90, 4.61 | 4.09e-10, 5.53e-10 | 29.84 |

(The level raised alone is the one with the largest part at the worst ∇φ box, level 3; its two
"defaults" rows come from different launches: 26.87 and 26.32 s, 6.93 and 6.73 s.) One level raised
takes off at most a third; every level with V pairs raised by one degree takes the worst box to
48.4 ε at kD = 1000 (2.4× less) for about 2–4 % more evaluation, by two to 21.7 ε for about 8 %
(one evaluation each, in one launch; the defaults' evaluation varies by 2 % between launches, 26.32
and 26.87 s; the build 14.2–14.8 s alike).

**The pair evidence** (spike S7's worst pair over the 316 offsets and three placements at the
levels' k a and degrees, at 10⁻⁸; levels 2–6 at kD = 1000; levels 2–5 at kD = 500 are the rows
of levels 3–6 at kD = 1000, the same (x, p)):

| x | p | worst pair ÷ ε at p | worst pair ÷ ε at p + 1 | worst pair ÷ ε at p + 2 | worst pair ÷ ε at p + 4 |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 216.51 | 248 | 3.97 | 2.31 | 1.35 | 0.367 |
| 108.25 | 134 | 4.18 | 1.74 | 0.744 | 0.155 |
| 54.13 | 75 | 3.25 | 1.28 | 0.563 | 0.0919 |
| 27.06 | 44 | 5.48 | 1.83 | 0.675 | 0.105 |
| 13.53 | 32 | 4.17 | 2.66 | 1.85 | 0.732 |

Per pair the rule's degrees are no outlier: 3.3–5.5 ε on every level, within the 19–42 ε per pair
that CH4.3c2's and CH4.3c3's FMM runs confirmed (the ℓ2 error here 0.09 ε at the samples), and one
degree more brings each to 1.3–2.7 ε. So no level's degree is short: the box exceeds 100 ε because
every level's truncation peaks at the same corner, not because one level is under-resolved.

**The remedy.** No single level's degree fixes it (one level at +2: 89.4 ε). Two remedies are
measured:
- the bound: §9.1's 100 ε holds on every case of D1's range (kD ≤ 500, signed off 2026-10-31: at
  most 68.2 ε, W1 with the plane wave's trace at 10⁻⁸), and the 114.8 ε lies at kD = 1000, beyond
  it. Stating the bound for D1's range, and kD = 1000 as reported against it, needs no code;
- every level with V pairs one degree up at 10⁻⁸ (a change of the degree rule's high-frequency
  branch at 10⁻⁸, x ≥ 13.5; the other ε not measured): 48.4 ε at kD = 1000 for about 2–4 % more
  evaluation, pair evidence 1.3–2.7 ε. It would be a recalibration (CH4.3c's scans and verification
  over D1, the constant density's included) for a gate bound, not for the tolerance: ε is §9.1's ℓ2
  tolerance, met with a margin here (0.09 ε at the samples), and the 100 ε bound is a check on the
  pointwise errors (CH4.0b, where a pointwise criterion cost 1.1–5.2× the evaluation). Scoping the
  bound to D1 leaves the one case beyond D1 above it, reported; that is question 2's choice.

Recommended: the first (question 2).

## 6. Questions for sign-off

1. **The acceptance as read** (plan, "Acceptance"): CH4.4a's cases with a level of the form or RCR
   at k a ≥ 2.8; the levels' parts under the defaults by `M2lOverride` (a level raised and on RCR);
   the M2L's parts by each non-RCR level on RCR, and the form's classes by decade below the bound's
   threshold; "no error concentrated at a class boundary" read one-sided, the default's worst box
   and ℓ2 at most 1.1× those of RCR on the level (amended in the plan from "either way" after the
   stage `low`: RCR errs up to 1.19× the waves', which is no error the waves add); the worst
   boxes' classes named by their effect, the classes' runs' part at each worst target (below 0.01
   everywhere), not by listing a box's V offsets' classes (they span nearly all 34); the sampled
   and refined reference, whose box figures equal CH4.4a's on W1 and differ on W2 and W3 by the
   normaliser (§2), CH4.4a's figures staying the gate's, with W3's worst boxes CH4.4a's own cells
   (`--refine-cells`, the stage `w3-cells`). Recommended: as read.
2. **The 114.7 ε box** (§5): its cause is the tangent-point placement on levels 2–5 at once (no
   level, class or M2L alone; the form, its classes and the waves add nothing). Recommended: §9.1's
   100 ε bound stated for D1's range (kD ≤ 500, where every worst box is at most 68.2 ε), and beyond
   it the worst box reported, not gated (proposal §9.1 amended; CH4.4 then closes with this
   finding). Alternatives: (b) a recalibration task, every level of the high-frequency branch one
   degree up where the box exceeds the bound (measured at 10⁻⁸: 114.8 → 48.4 ε at kD = 1000, about 2–4 % more evaluation; pair evidence
   1.3–2.7 ε), on CH4.3c's scans; (c) §9.1's bound raised (to 150 ε, say) for every kD.
3. **`M2lOverride` in `tetraleaf-exec`'s public API** as a measurement setting (documented, empty by
   default, no rule sets one), as `diagonal_costs` and `v_costs` are. Recommended: keep it public.
   Alternative: `#[doc(hidden)]`.

## 7. Tables

### Summary per launch

The worst box over ε (the refined reference), ℓ2 over the reference's targets, the level with the
largest share of ∇φ's Σ|e|² (its M2L), the rest of the levels' parts, the largest change of the worst
box and of ℓ2 (∇φ) over the M2L runs (RCR's over the default's) with its run, and the form's classes
beside the threshold by level. W3 beyond kD = 30 from the stage `w3-cells` (`-cells`).

| launch | worst box ÷ ε φ, ∇φ | ℓ2 φ, ∇φ | level (share of Σ|Δ∇φ|²) | rest φ, ∇φ | largest M2L change: worst box ∇φ (run) | ℓ2 ∇φ (run) | classes in (0.1, 1] |
| --- | --- | --- | --- | --- | --- | --- | --- |
| w1-kd100-constant-1e-3 | 2.92, 8.24 | 2.78e-4, 4.38e-4 | 3 (diagonal 34/34, 0.76) | 0.019, 0.025 | 1.000 (level 2 RCR (was diagonal 34/34)) | 1.000 (level 2 RCR (was diagonal 34/34)) | – |
| w1-kd100-constant-1e-6 | 5.55, 15.02 | 2.11e-7, 4.34e-7 | 4 (diagonal 34/34, 0.47) | 0.002, 0.002 | 1.000 (level 2 RCR (was diagonal 34/34)) | 1.000 (level 2 RCR (was diagonal 34/34)) | 4: [0, 1] |
| w1-kd100-constant-1e-8 | 22.48, 34.85 | 3.99e-9, 6.82e-9 | 2 (diagonal 32/34, 0.68) | 0.001, 0.001 | 1.002 (level 3 RCR (was waves)) | 1.000 (level 2 RCR (was diagonal 32/34)) | – |
| w1-kd100-random-1e-3 | 5.43, 7.20 | 1.68e-4, 2.42e-4 | 3 (diagonal 34/34, 0.63) | 0.015, 0.020 | 1.000 (level 2 RCR (was diagonal 34/34)) | 1.000 (level 2 RCR (was diagonal 34/34)) | – |
| w1-kd100-random-1e-6 | 4.51, 14.72 | 1.27e-7, 2.92e-7 | 4 (diagonal 34/34, 0.64) | 0.001, 0.001 | 1.000 (level 2 RCR (was diagonal 34/34)) | 1.000 (level 2 RCR (was diagonal 34/34)) | 4: [0, 1] |
| w1-kd100-random-1e-8 | 5.04, 8.67 | 7.64e-10, 1.35e-9 | 2 (diagonal 32/34, 0.55) | 0.001, 0.001 | 1.000 (level 2 RCR (was diagonal 32/34)) | 1.000 (level 2 RCR (was diagonal 32/34)) | – |
| w1-kd1000-planewave-1e-3 | 35.56, 43.34 | 5.85e-4, 9.72e-4 | 6 (diagonal 34/34, 0.30) | 0.035, 0.033 | 1.000 (level 2 RCR (was diagonal 34/34)) | 1.000 (level 2 RCR (was diagonal 34/34)) | – |
| w1-kd1000-planewave-1e-6 | 46.66, 73.88 | 4.31e-7, 7.79e-7 | 7 (diagonal 34/34, 0.28) | 0.003, 0.003 | 1.000 (level 2 RCR (was diagonal 34/34)) | 1.000 (level 2 RCR (was diagonal 34/34)) | – |
| w1-kd1000-planewave-1e-8 | 69.49, 114.80 | 5.43e-9, 9.17e-9 | 3 (diagonal 34/34, 0.27) | 0.001, 0.001 | 1.000 (level 2 RCR (was diagonal 34/34)) | 1.000 (level 2 RCR (was diagonal 34/34)) | 6: [7, 8] |
| w1-kd1000-random-1e-3 | 20.03, 20.25 | 6.08e-4, 8.31e-4 | 6 (diagonal 34/34, 0.33) | 0.033, 0.031 | 1.000 (level 2 RCR (was diagonal 34/34)) | 1.000 (level 2 RCR (was diagonal 34/34)) | – |
| w1-kd1000-random-1e-6 | 26.05, 27.99 | 5.28e-7, 8.64e-7 | 7 (diagonal 34/34, 0.50) | 0.003, 0.002 | 1.000 (level 2 RCR (was diagonal 34/34)) | 1.000 (level 2 RCR (was diagonal 34/34)) | – |
| w1-kd1000-random-1e-8 | 42.00, 43.30 | 6.66e-9, 8.48e-9 | 2 (diagonal 34/34, 0.30) | 0.001, 0.001 | 1.000 (level 2 RCR (was diagonal 34/34)) | 1.000 (level 2 RCR (was diagonal 34/34)) | 6: [7, 8] |
| w1-kd250-random-1e-3 | 13.24, 13.54 | 7.01e-4, 9.37e-4 | 3 (diagonal 34/34, 0.19) | 0.030, 0.038 | 1.000 (level 2 RCR (was diagonal 34/34)) | 1.000 (level 2 RCR (was diagonal 34/34)) | – |
| w1-kd250-random-1e-6 | 8.84, 16.90 | 4.74e-7, 9.15e-7 | 5 (diagonal 34/34, 0.49) | 0.002, 0.002 | 1.000 (level 2 RCR (was diagonal 34/34)) | 1.000 (level 2 RCR (was diagonal 34/34)) | – |
| w1-kd250-random-1e-8 | 16.09, 21.45 | 6.17e-9, 1.03e-8 | 2 (diagonal 34/34, 0.45) | 0.001, 0.001 | 1.000 (level 2 RCR (was diagonal 34/34)) | 1.000 (level 2 RCR (was diagonal 34/34)) | 4: [7, 8] |
| w1-kd30-constant-1e-6 | 0.76, 5.45 | 5.41e-8, 2.55e-7 | 2 (diagonal 34/34, 1.00) | 0.000, 0.000 | 1.000 (level 2 RCR (was diagonal 34/34)) | 1.000 (level 2 RCR (was diagonal 34/34)) | – |
| w1-kd30-random-1e-6 | 3.34, 8.26 | 1.39e-7, 3.26e-7 | 2 (diagonal 34/34, 1.00) | 0.000, 0.000 | 1.000 (level 2 RCR (was diagonal 34/34)) | 1.000 (level 2 RCR (was diagonal 34/34)) | – |
| w1-kd500-planewave-1e-3 | 27.72, 39.61 | 6.55e-4, 1.08e-3 | 6 (diagonal 34/34, 0.20) | 0.035, 0.040 | 1.000 (level 2 RCR (was diagonal 34/34)) | 1.000 (level 2 RCR (was diagonal 34/34)) | – |
| w1-kd500-planewave-1e-6 | 26.14, 42.42 | 4.78e-7, 9.43e-7 | 6 (diagonal 34/34, 0.26) | 0.003, 0.004 | 1.000 (level 2 RCR (was diagonal 34/34)) | 1.000 (level 2 RCR (was diagonal 34/34)) | – |
| w1-kd500-planewave-1e-8 | 35.63, 68.21 | 6.18e-9, 1.14e-8 | 3 (diagonal 34/34, 0.30) | 0.001, 0.001 | 1.000 (level 2 RCR (was diagonal 34/34)) | 1.000 (level 2 RCR (was diagonal 34/34)) | 5: [7, 8] |
| w1-kd500-random-1e-3 | 21.43, 23.01 | 7.30e-4, 9.44e-4 | 6 (diagonal 34/34, 0.22) | 0.036, 0.039 | 1.000 (level 2 RCR (was diagonal 34/34)) | 1.000 (level 2 RCR (was diagonal 34/34)) | – |
| w1-kd500-random-1e-6 | 14.32, 20.03 | 4.94e-7, 9.23e-7 | 6 (diagonal 34/34, 0.54) | 0.003, 0.002 | 1.000 (level 2 RCR (was diagonal 34/34)) | 1.000 (level 2 RCR (was diagonal 34/34)) | – |
| w1-kd500-random-1e-8 | 27.59, 40.47 | 6.20e-9, 9.28e-9 | 2 (diagonal 34/34, 0.33) | 0.001, 0.001 | 1.000 (level 2 RCR (was diagonal 34/34)) | 1.000 (level 2 RCR (was diagonal 34/34)) | 5: [7, 8] |
| w2-kd100-constant-1e-3 | 2.75, 24.24 | 1.24e-4, 6.39e-4 | 3 (diagonal 34/34, 0.46) | 0.011, 0.010 | 1.000 (level 2 RCR (was diagonal 34/34)) | 1.000 (level 2 RCR (was diagonal 34/34)) | – |
| w2-kd100-constant-1e-6 | 3.18, 34.68 | 1.05e-7, 7.66e-7 | 4 (diagonal 34/34, 0.57) | 0.001, 0.001 | 1.053 (level 5 RCR (was waves)) | 1.009 (level 5 RCR (was waves)) | 4: [0, 1] |
| w2-kd100-constant-1e-8 | 6.56, 37.55 | 8.05e-10, 3.76e-9 | 2 (diagonal 32/34, 0.51) | 0.001, 0.001 | 1.000 (level 2 RCR (was diagonal 32/34)) | 1.001 (level 3 RCR (was waves)) | – |
| w2-kd100-random-1e-3 | 3.18, 2.47 | 1.32e-4, 1.23e-4 | 5 (Batched, 0.38) | 0.011, 0.011 | 1.000 (level 2 RCR (was diagonal 34/34)) | 1.000 (level 2 RCR (was diagonal 34/34)) | – |
| w2-kd100-random-1e-6 | 10.59, 11.31 | 1.47e-7, 2.17e-7 | 4 (diagonal 34/34, 0.57) | 0.000, 0.000 | 0.995 (level 5 RCR (was waves)) | 1.015 (level 5 RCR (was waves)) | 4: [0, 1] |
| w2-kd100-random-1e-8 | 3.32, 2.88 | 5.45e-10, 4.70e-10 | 2 (diagonal 32/34, 0.46) | 0.001, 0.001 | 1.001 (level 3 RCR (was waves)) | 1.004 (level 4 RCR (was waves)) | – |
| w2-kd250-random-1e-3 | 15.91, 10.84 | 3.69e-4, 3.93e-4 | 6 (Batched, 0.27) | 0.025, 0.027 | 1.000 (level 2 RCR (was diagonal 34/34)) | 1.000 (level 2 RCR (was diagonal 34/34)) | – |
| w2-kd250-random-1e-6 | 18.79, 23.81 | 2.93e-7, 4.81e-7 | 5 (diagonal 34/34, 0.23) | 0.002, 0.001 | 1.100 (level 7 RCR (was waves)) | 1.027 (level 7 RCR (was waves)) | 6: [13, 14, 15, 16, 17] |
| w2-kd250-random-1e-8 | 17.08, 17.18 | 2.62e-9, 3.26e-9 | 3 (diagonal 32/34, 0.48) | 0.001, 0.001 | 0.990 (level 5 RCR (was waves)) | 1.002 (level 5 RCR (was waves)) | 4: [7, 8] |
| w2-kd30-constant-1e-6 | 0.44, 7.53 | 2.90e-8, 2.09e-7 | 2 (diagonal 34/34, 0.36) | 0.000, 0.000 | 1.190 (level 4 RCR (was waves)) | 1.068 (level 4 RCR (was waves)) | – |
| w2-kd30-random-1e-6 | 2.34, 7.99 | 1.00e-7, 2.08e-7 | 3 (RCR, 0.42) | 0.000, 0.000 | 1.147 (level 4 RCR (was waves)) | 1.064 (level 4 RCR (was waves)) | – |
| w2-kd500-random-1e-3 | 21.78, 22.13 | 3.89e-4, 4.27e-4 | 7 (Batched, 0.24) | 0.032, 0.034 | 1.000 (level 2 RCR (was diagonal 34/34)) | 1.000 (level 2 RCR (was diagonal 34/34)) | – |
| w2-kd500-random-1e-6 | 18.76, 32.34 | 2.37e-7, 4.31e-7 | 6 (diagonal 34/34, 0.33) | 0.003, 0.001 | 1.034 (level 8 RCR (was waves)) | 1.017 (level 8 RCR (was waves)) | 7: [13, 14, 15, 16, 17] |
| w2-kd500-random-1e-8 | 22.39, 28.27 | 2.47e-9, 2.83e-9 | 4 (diagonal 32/34, 0.32) | 0.001, 0.001 | 0.978 (level 7 RCR (was waves)) | 0.999 (level 6 RCR (was waves)) | 5: [7, 8] |
| w3-kd100-random-1e-3-cells | 11.02, 3.85 | 3.68e-4, 1.51e-4 | 3 (diagonal 34/34, 0.64) | 0.009, 0.014 | 1.000 (level 2 RCR (was diagonal 34/34)) | 1.000 (level 2 RCR (was diagonal 34/34)) | – |
| w3-kd100-random-1e-6-cells | 13.55, 5.83 | 2.29e-7, 1.35e-7 | 4 (waves, 0.32) | 0.001, 0.001 | 1.003 (level 4 RCR (was waves)) | 1.005 (level 5 RCR (was waves)) | – |
| w3-kd100-random-1e-8-cells | 9.72, 4.83 | 1.49e-9, 8.14e-10 | 3 (waves, 0.48) | 0.000, 0.001 | 1.021 (level 3 RCR (was waves)) | 1.001 (level 3 RCR (was waves)) | – |
| w3-kd30-random-1e-6 | 4.10, 3.88 | 1.28e-7, 1.16e-7 | 2 (diagonal 34/34, 0.59) | 0.000, 0.000 | 1.028 (level 3 RCR (was waves)) | 1.015 (level 3 RCR (was waves)) | – |
| w3:1000000-kd200-random-1e-3-cells | 10.84, 12.01 | 3.47e-4, 4.47e-4 | 4 (waves, 0.76) | 0.009, 0.011 | 1.000 (level 2 RCR (was diagonal 34/34)) | 1.000 (level 2 RCR (was diagonal 34/34)) | – |
| w3:1000000-kd200-random-1e-6-cells | 10.23, 17.97 | 1.95e-7, 3.00e-7 | 4 (waves, 0.49) | 0.001, 0.002 | 0.998 (level 4 RCR (was waves)) | 1.000 (level 2 RCR (was diagonal 34/34)) | 5: [0, 1] |
| w3:1000000-kd200-random-1e-8-cells | 9.44, 13.03 | 1.69e-9, 2.84e-9 | 4 (waves, 0.60) | 0.000, 0.000 | 0.999 (level 4 RCR (was waves)) | 1.000 (level 2 RCR (was diagonal 34/34)) | 2: [0, 1] |
| w3:1000000-kd500-random-1e-3-cells | 16.45, 19.03 | 3.96e-4, 5.38e-4 | 5 (diagonal 34/34, 0.49) | 0.017, 0.018 | 1.000 (level 2 RCR (was diagonal 34/34)) | 1.000 (level 2 RCR (was diagonal 34/34)) | – |
| w3:1000000-kd500-random-1e-6-cells | 15.41, 18.50 | 1.96e-7, 2.95e-7 | 4 (diagonal 34/34, 0.56) | 0.002, 0.003 | 1.000 (level 2 RCR (was diagonal 34/34)) | 1.000 (level 2 RCR (was diagonal 34/34)) | – |
| w3:1000000-kd500-random-1e-8-cells | 24.61, 29.71 | 2.70e-9, 4.24e-9 | 4 (waves, 0.47) | 0.001, 0.001 | 1.000 (level 2 RCR (was diagonal 34/34)) | 1.000 (level 2 RCR (was diagonal 34/34)) | – |

### Levels

Per launch and level: k a, the degree as run, the M2L, the share of Σ|e|² (φ, ∇φ) of the level's part,
and on the form's levels the classes by decade of E/(θ ε) ((0.1, 1] listed, (0.01, 0.1] counted)
and those that fail the bound.

| launch | level | k a | p | M2L | share φ | share ∇φ | classes (0.1, 1] | (0.01, 0.1] | failing |
| --- | ---: | ---: | ---: | --- | ---: | ---: | --- | ---: | ---: |
| w1-kd100-constant-1e-3 | 2 | 21.65 | 27 | diagonal 34/34 | 0.376 | 0.177 | [] | 0 | 0 |
| w1-kd100-constant-1e-3 | 3 | 10.83 | 15 | diagonal 34/34 | 0.608 | 0.755 | [] | 0 | 0 |
| w1-kd100-constant-1e-3 | 4 | 5.41 | 10 | diagonal 34/34 | 0.042 | 0.119 | [] | 0 | 0 |
| w1-kd100-constant-1e-6 | 2 | 21.65 | 34 | diagonal 34/34 | 0.672 | 0.423 | [] | 0 | 0 |
| w1-kd100-constant-1e-6 | 3 | 10.83 | 21 | diagonal 34/34 | 0.110 | 0.120 | [] | 0 | 0 |
| w1-kd100-constant-1e-6 | 4 | 5.41 | 15 | diagonal 34/34 | 0.141 | 0.474 | [0, 1] | 3 | 0 |
| w1-kd100-constant-1e-8 | 2 | 21.65 | 38 | diagonal 32/34 | 0.785 | 0.682 | [] | 3 | 2 |
| w1-kd100-constant-1e-8 | 3 | 10.83 | 32 | waves | 0.042 | 0.104 |  |  |  |
| w1-kd100-constant-1e-8 | 4 | 5.41 | 32 | – | – | – |  |  |  |
| w1-kd100-random-1e-3 | 2 | 21.65 | 27 | diagonal 34/34 | 0.169 | 0.102 | [] | 0 | 0 |
| w1-kd100-random-1e-3 | 3 | 10.83 | 15 | diagonal 34/34 | 0.598 | 0.627 | [] | 0 | 0 |
| w1-kd100-random-1e-3 | 4 | 5.41 | 10 | diagonal 34/34 | 0.073 | 0.152 | [] | 0 | 0 |
| w1-kd100-random-1e-6 | 2 | 21.65 | 34 | diagonal 34/34 | 0.321 | 0.119 | [] | 0 | 0 |
| w1-kd100-random-1e-6 | 3 | 10.83 | 21 | diagonal 34/34 | 0.225 | 0.161 | [] | 0 | 0 |
| w1-kd100-random-1e-6 | 4 | 5.41 | 15 | diagonal 34/34 | 0.369 | 0.642 | [0, 1] | 3 | 0 |
| w1-kd100-random-1e-8 | 2 | 21.65 | 38 | diagonal 32/34 | 0.706 | 0.553 | [] | 3 | 2 |
| w1-kd100-random-1e-8 | 3 | 10.83 | 32 | waves | 0.266 | 0.421 |  |  |  |
| w1-kd100-random-1e-8 | 4 | 5.41 | 32 | – | – | – |  |  |  |
| w1-kd1000-planewave-1e-3 | 2 | 216.51 | 225 | diagonal 34/34 | 0.044 | 0.019 | [] | 0 | 0 |
| w1-kd1000-planewave-1e-3 | 3 | 108.25 | 115 | diagonal 34/34 | 0.137 | 0.067 | [] | 0 | 0 |
| w1-kd1000-planewave-1e-3 | 4 | 54.13 | 60 | diagonal 34/34 | 0.080 | 0.050 | [] | 0 | 0 |
| w1-kd1000-planewave-1e-3 | 5 | 27.06 | 32 | diagonal 34/34 | 0.107 | 0.093 | [] | 0 | 0 |
| w1-kd1000-planewave-1e-3 | 6 | 13.53 | 18 | diagonal 34/34 | 0.219 | 0.297 | [] | 0 | 0 |
| w1-kd1000-planewave-1e-3 | 7 | 6.77 | 11 | diagonal 34/34 | 0.091 | 0.238 | [] | 0 | 0 |
| w1-kd1000-planewave-1e-6 | 2 | 216.51 | 240 | diagonal 34/34 | 0.296 | 0.145 | [] | 0 | 0 |
| w1-kd1000-planewave-1e-6 | 3 | 108.25 | 128 | diagonal 34/34 | 0.153 | 0.093 | [] | 0 | 0 |
| w1-kd1000-planewave-1e-6 | 4 | 54.13 | 70 | diagonal 34/34 | 0.136 | 0.103 | [] | 0 | 0 |
| w1-kd1000-planewave-1e-6 | 5 | 27.06 | 40 | diagonal 34/34 | 0.101 | 0.127 | [] | 0 | 0 |
| w1-kd1000-planewave-1e-6 | 6 | 13.53 | 25 | diagonal 34/34 | 0.019 | 0.048 | [] | 0 | 0 |
| w1-kd1000-planewave-1e-6 | 7 | 6.77 | 16 | diagonal 34/34 | 0.057 | 0.279 | [] | 2 | 0 |
| w1-kd1000-planewave-1e-8 | 2 | 216.51 | 248 | diagonal 34/34 | 0.328 | 0.215 | [] | 0 | 0 |
| w1-kd1000-planewave-1e-8 | 3 | 108.25 | 134 | diagonal 34/34 | 0.343 | 0.267 | [] | 0 | 0 |
| w1-kd1000-planewave-1e-8 | 4 | 54.13 | 75 | diagonal 34/34 | 0.145 | 0.171 | [] | 2 | 0 |
| w1-kd1000-planewave-1e-8 | 5 | 27.06 | 44 | diagonal 32/34 | 0.085 | 0.160 | [] | 0 | 2 |
| w1-kd1000-planewave-1e-8 | 6 | 13.53 | 32 | diagonal 27/34 | 0.019 | 0.064 | [7, 8] | 4 | 7 |
| w1-kd1000-planewave-1e-8 | 7 | 6.77 | 32 | – | – | – |  |  |  |
| w1-kd1000-random-1e-3 | 2 | 216.51 | 225 | diagonal 34/34 | 0.086 | 0.042 | [] | 0 | 0 |
| w1-kd1000-random-1e-3 | 3 | 108.25 | 115 | diagonal 34/34 | 0.158 | 0.081 | [] | 0 | 0 |
| w1-kd1000-random-1e-3 | 4 | 54.13 | 60 | diagonal 34/34 | 0.066 | 0.042 | [] | 0 | 0 |
| w1-kd1000-random-1e-3 | 5 | 27.06 | 32 | diagonal 34/34 | 0.114 | 0.093 | [] | 0 | 0 |
| w1-kd1000-random-1e-3 | 6 | 13.53 | 18 | diagonal 34/34 | 0.249 | 0.331 | [] | 0 | 0 |
| w1-kd1000-random-1e-3 | 7 | 6.77 | 11 | diagonal 34/34 | 0.108 | 0.256 | [] | 0 | 0 |
| w1-kd1000-random-1e-6 | 2 | 216.51 | 240 | diagonal 34/34 | 0.476 | 0.201 | [] | 0 | 0 |
| w1-kd1000-random-1e-6 | 3 | 108.25 | 128 | diagonal 34/34 | 0.114 | 0.067 | [] | 0 | 0 |
| w1-kd1000-random-1e-6 | 4 | 54.13 | 70 | diagonal 34/34 | 0.075 | 0.058 | [] | 0 | 0 |
| w1-kd1000-random-1e-6 | 5 | 27.06 | 40 | diagonal 34/34 | 0.058 | 0.058 | [] | 0 | 0 |
| w1-kd1000-random-1e-6 | 6 | 13.53 | 25 | diagonal 34/34 | 0.018 | 0.033 | [] | 0 | 0 |
| w1-kd1000-random-1e-6 | 7 | 6.77 | 16 | diagonal 34/34 | 0.142 | 0.503 | [] | 2 | 0 |
| w1-kd1000-random-1e-8 | 2 | 216.51 | 248 | diagonal 34/34 | 0.430 | 0.300 | [] | 0 | 0 |
| w1-kd1000-random-1e-8 | 3 | 108.25 | 134 | diagonal 34/34 | 0.271 | 0.258 | [] | 0 | 0 |
| w1-kd1000-random-1e-8 | 4 | 54.13 | 75 | diagonal 34/34 | 0.074 | 0.107 | [] | 2 | 0 |
| w1-kd1000-random-1e-8 | 5 | 27.06 | 44 | diagonal 32/34 | 0.056 | 0.124 | [] | 0 | 2 |
| w1-kd1000-random-1e-8 | 6 | 13.53 | 32 | diagonal 27/34 | 0.018 | 0.068 | [7, 8] | 4 | 7 |
| w1-kd1000-random-1e-8 | 7 | 6.77 | 32 | – | – | – |  |  |  |
| w1-kd250-random-1e-3 | 2 | 54.13 | 60 | diagonal 34/34 | 0.253 | 0.171 | [] | 0 | 0 |
| w1-kd250-random-1e-3 | 3 | 27.06 | 32 | diagonal 34/34 | 0.203 | 0.194 | [] | 0 | 0 |
| w1-kd250-random-1e-3 | 4 | 13.53 | 18 | diagonal 34/34 | 0.105 | 0.140 | [] | 0 | 0 |
| w1-kd250-random-1e-3 | 5 | 6.77 | 11 | diagonal 34/34 | 0.072 | 0.175 | [] | 0 | 0 |
| w1-kd250-random-1e-6 | 2 | 54.13 | 70 | diagonal 34/34 | 0.454 | 0.208 | [] | 0 | 0 |
| w1-kd250-random-1e-6 | 3 | 27.06 | 40 | diagonal 34/34 | 0.273 | 0.207 | [] | 0 | 0 |
| w1-kd250-random-1e-6 | 4 | 13.53 | 25 | diagonal 34/34 | 0.012 | 0.017 | [] | 0 | 0 |
| w1-kd250-random-1e-6 | 5 | 6.77 | 16 | diagonal 34/34 | 0.194 | 0.489 | [] | 2 | 0 |
| w1-kd250-random-1e-8 | 2 | 54.13 | 75 | diagonal 34/34 | 0.529 | 0.450 | [] | 2 | 0 |
| w1-kd250-random-1e-8 | 3 | 27.06 | 44 | diagonal 32/34 | 0.327 | 0.381 | [] | 0 | 2 |
| w1-kd250-random-1e-8 | 4 | 13.53 | 32 | diagonal 27/34 | 0.010 | 0.024 | [7, 8] | 4 | 7 |
| w1-kd250-random-1e-8 | 5 | 6.77 | 32 | – | – | – |  |  |  |
| w1-kd30-constant-1e-6 | 2 | 6.49 | 16 | diagonal 34/34 | 1.000 | 1.000 | [] | 2 | 0 |
| w1-kd30-random-1e-6 | 2 | 6.49 | 16 | diagonal 34/34 | 1.000 | 1.000 | [] | 2 | 0 |
| w1-kd500-planewave-1e-3 | 2 | 108.25 | 115 | diagonal 34/34 | 0.138 | 0.069 | [] | 0 | 0 |
| w1-kd500-planewave-1e-3 | 3 | 54.13 | 60 | diagonal 34/34 | 0.136 | 0.089 | [] | 0 | 0 |
| w1-kd500-planewave-1e-3 | 4 | 27.06 | 32 | diagonal 34/34 | 0.118 | 0.105 | [] | 0 | 0 |
| w1-kd500-planewave-1e-3 | 5 | 13.53 | 18 | diagonal 34/34 | 0.137 | 0.197 | [] | 0 | 0 |
| w1-kd500-planewave-1e-3 | 6 | 6.77 | 11 | diagonal 34/34 | 0.075 | 0.199 | [] | 0 | 0 |
| w1-kd500-planewave-1e-6 | 2 | 108.25 | 128 | diagonal 34/34 | 0.170 | 0.086 | [] | 0 | 0 |
| w1-kd500-planewave-1e-6 | 3 | 54.13 | 70 | diagonal 34/34 | 0.320 | 0.216 | [] | 0 | 0 |
| w1-kd500-planewave-1e-6 | 4 | 27.06 | 40 | diagonal 34/34 | 0.143 | 0.151 | [] | 0 | 0 |
| w1-kd500-planewave-1e-6 | 5 | 13.53 | 25 | diagonal 34/34 | 0.014 | 0.028 | [] | 0 | 0 |
| w1-kd500-planewave-1e-6 | 6 | 6.77 | 16 | diagonal 34/34 | 0.065 | 0.256 | [] | 2 | 0 |
| w1-kd500-planewave-1e-8 | 2 | 108.25 | 134 | diagonal 34/34 | 0.216 | 0.145 | [] | 0 | 0 |
| w1-kd500-planewave-1e-8 | 3 | 54.13 | 75 | diagonal 34/34 | 0.311 | 0.302 | [] | 2 | 0 |
| w1-kd500-planewave-1e-8 | 4 | 27.06 | 44 | diagonal 32/34 | 0.127 | 0.219 | [] | 0 | 2 |
| w1-kd500-planewave-1e-8 | 5 | 13.53 | 32 | diagonal 27/34 | 0.011 | 0.031 | [7, 8] | 4 | 7 |
| w1-kd500-planewave-1e-8 | 6 | 6.77 | 32 | – | – | – |  |  |  |
| w1-kd500-random-1e-3 | 2 | 108.25 | 115 | diagonal 34/34 | 0.240 | 0.136 | [] | 0 | 0 |
| w1-kd500-random-1e-3 | 3 | 54.13 | 60 | diagonal 34/34 | 0.118 | 0.083 | [] | 0 | 0 |
| w1-kd500-random-1e-3 | 4 | 27.06 | 32 | diagonal 34/34 | 0.094 | 0.092 | [] | 0 | 0 |
| w1-kd500-random-1e-3 | 5 | 13.53 | 18 | diagonal 34/34 | 0.144 | 0.210 | [] | 0 | 0 |
| w1-kd500-random-1e-3 | 6 | 6.77 | 11 | diagonal 34/34 | 0.083 | 0.221 | [] | 0 | 0 |
| w1-kd500-random-1e-6 | 2 | 108.25 | 128 | diagonal 34/34 | 0.286 | 0.106 | [] | 0 | 0 |
| w1-kd500-random-1e-6 | 3 | 54.13 | 70 | diagonal 34/34 | 0.269 | 0.138 | [] | 0 | 0 |
| w1-kd500-random-1e-6 | 4 | 27.06 | 40 | diagonal 34/34 | 0.140 | 0.110 | [] | 0 | 0 |
| w1-kd500-random-1e-6 | 5 | 13.53 | 25 | diagonal 34/34 | 0.017 | 0.023 | [] | 0 | 0 |
| w1-kd500-random-1e-6 | 6 | 6.77 | 16 | diagonal 34/34 | 0.199 | 0.538 | [] | 2 | 0 |
| w1-kd500-random-1e-8 | 2 | 108.25 | 134 | diagonal 34/34 | 0.451 | 0.325 | [] | 0 | 0 |
| w1-kd500-random-1e-8 | 3 | 54.13 | 75 | diagonal 34/34 | 0.251 | 0.286 | [] | 2 | 0 |
| w1-kd500-random-1e-8 | 4 | 27.06 | 44 | diagonal 32/34 | 0.118 | 0.175 | [] | 0 | 2 |
| w1-kd500-random-1e-8 | 5 | 13.53 | 32 | diagonal 27/34 | 0.016 | 0.044 | [7, 8] | 4 | 7 |
| w1-kd500-random-1e-8 | 6 | 6.77 | 32 | – | – | – |  |  |  |
| w2-kd100-constant-1e-3 | 2 | 21.65 | 27 | diagonal 34/34 | 0.231 | 0.072 | [] | 0 | 0 |
| w2-kd100-constant-1e-3 | 3 | 10.83 | 15 | diagonal 34/34 | 0.668 | 0.456 | [] | 0 | 0 |
| w2-kd100-constant-1e-3 | 4 | 5.41 | 10 | diagonal 34/34 | 0.089 | 0.143 | [] | 0 | 0 |
| w2-kd100-constant-1e-3 | 5 | 2.71 | 7 | Batched | 0.090 | 0.410 |  |  |  |
| w2-kd100-constant-1e-3 | 6 | 1.35 | 7 | Batched | 0.000 | 0.009 |  |  |  |
| w2-kd100-constant-1e-6 | 2 | 21.65 | 34 | diagonal 34/34 | 0.284 | 0.076 | [] | 0 | 0 |
| w2-kd100-constant-1e-6 | 3 | 10.83 | 21 | diagonal 34/34 | 0.310 | 0.248 | [] | 0 | 0 |
| w2-kd100-constant-1e-6 | 4 | 5.41 | 15 | diagonal 34/34 | 0.342 | 0.572 | [0, 1] | 3 | 0 |
| w2-kd100-constant-1e-6 | 5 | 2.71 | 15 | waves | 0.023 | 0.052 |  |  |  |
| w2-kd100-constant-1e-6 | 6 | 1.35 | 15 | waves | 0.004 | 0.008 |  |  |  |
| w2-kd100-constant-1e-8 | 2 | 21.65 | 38 | diagonal 32/34 | 0.704 | 0.513 | [] | 3 | 2 |
| w2-kd100-constant-1e-8 | 3 | 10.83 | 32 | waves | 0.111 | 0.278 |  |  |  |
| w2-kd100-constant-1e-8 | 4 | 5.41 | 32 | waves | 0.000 | 0.000 |  |  |  |
| w2-kd100-constant-1e-8 | 5 | 2.71 | 32 | waves | 0.000 | 0.000 |  |  |  |
| w2-kd100-constant-1e-8 | 6 | 1.35 | 32 | – | – | – |  |  |  |
| w2-kd100-random-1e-3 | 2 | 21.65 | 27 | diagonal 34/34 | 0.105 | 0.040 | [] | 0 | 0 |
| w2-kd100-random-1e-3 | 3 | 10.83 | 15 | diagonal 34/34 | 0.572 | 0.346 | [] | 0 | 0 |
| w2-kd100-random-1e-3 | 4 | 5.41 | 10 | diagonal 34/34 | 0.089 | 0.119 | [] | 0 | 0 |
| w2-kd100-random-1e-3 | 5 | 2.71 | 7 | Batched | 0.107 | 0.379 |  |  |  |
| w2-kd100-random-1e-3 | 6 | 1.35 | 7 | Batched | 0.003 | 0.046 |  |  |  |
| w2-kd100-random-1e-6 | 2 | 21.65 | 34 | diagonal 34/34 | 0.103 | 0.026 | [] | 0 | 0 |
| w2-kd100-random-1e-6 | 3 | 10.83 | 21 | diagonal 34/34 | 0.164 | 0.066 | [] | 0 | 0 |
| w2-kd100-random-1e-6 | 4 | 5.41 | 15 | diagonal 34/34 | 0.539 | 0.571 | [0, 1] | 3 | 0 |
| w2-kd100-random-1e-6 | 5 | 2.71 | 15 | waves | 0.050 | 0.221 |  |  |  |
| w2-kd100-random-1e-6 | 6 | 1.35 | 15 | waves | 0.005 | 0.095 |  |  |  |
| w2-kd100-random-1e-8 | 2 | 21.65 | 38 | diagonal 32/34 | 0.576 | 0.464 | [] | 3 | 2 |
| w2-kd100-random-1e-8 | 3 | 10.83 | 32 | waves | 0.202 | 0.287 |  |  |  |
| w2-kd100-random-1e-8 | 4 | 5.41 | 32 | waves | 0.000 | 0.000 |  |  |  |
| w2-kd100-random-1e-8 | 5 | 2.71 | 32 | waves | 0.000 | 0.000 |  |  |  |
| w2-kd100-random-1e-8 | 6 | 1.35 | 32 | – | – | – |  |  |  |
| w2-kd250-random-1e-3 | 2 | 54.13 | 60 | diagonal 34/34 | 0.134 | 0.057 | [] | 0 | 0 |
| w2-kd250-random-1e-3 | 3 | 27.06 | 32 | diagonal 34/34 | 0.201 | 0.108 | [] | 0 | 0 |
| w2-kd250-random-1e-3 | 4 | 13.53 | 18 | diagonal 34/34 | 0.246 | 0.215 | [] | 0 | 0 |
| w2-kd250-random-1e-3 | 5 | 6.77 | 11 | diagonal 34/34 | 0.106 | 0.158 | [] | 0 | 0 |
| w2-kd250-random-1e-3 | 6 | 3.38 | 7 | Batched | 0.098 | 0.272 |  |  |  |
| w2-kd250-random-1e-3 | 7 | 1.69 | 7 | Batched | 0.002 | 0.026 |  |  |  |
| w2-kd250-random-1e-6 | 2 | 54.13 | 70 | diagonal 34/34 | 0.216 | 0.066 | [] | 0 | 0 |
| w2-kd250-random-1e-6 | 3 | 27.06 | 40 | diagonal 34/34 | 0.344 | 0.154 | [] | 0 | 0 |
| w2-kd250-random-1e-6 | 4 | 13.53 | 25 | diagonal 34/34 | 0.023 | 0.017 | [] | 0 | 0 |
| w2-kd250-random-1e-6 | 5 | 6.77 | 16 | diagonal 34/34 | 0.158 | 0.233 | [] | 2 | 0 |
| w2-kd250-random-1e-6 | 6 | 3.38 | 15 | diagonal 21/34 | 0.034 | 0.178 | [13, 14, 15, 16, 17] | 4 | 13 |
| w2-kd250-random-1e-6 | 7 | 1.69 | 15 | waves | 0.009 | 0.222 |  |  |  |
| w2-kd250-random-1e-8 | 2 | 54.13 | 75 | diagonal 34/34 | 0.501 | 0.318 | [] | 2 | 0 |
| w2-kd250-random-1e-8 | 3 | 27.06 | 44 | diagonal 32/34 | 0.470 | 0.480 | [] | 0 | 2 |
| w2-kd250-random-1e-8 | 4 | 13.53 | 32 | diagonal 27/34 | 0.026 | 0.042 | [7, 8] | 4 | 7 |
| w2-kd250-random-1e-8 | 5 | 6.77 | 32 | waves | 0.000 | 0.000 |  |  |  |
| w2-kd250-random-1e-8 | 6 | 3.38 | 32 | waves | 0.000 | 0.000 |  |  |  |
| w2-kd250-random-1e-8 | 7 | 1.69 | 32 | – | – | – |  |  |  |
| w2-kd30-constant-1e-6 | 2 | 6.49 | 16 | diagonal 34/34 | 0.815 | 0.356 | [] | 2 | 0 |
| w2-kd30-constant-1e-6 | 3 | 3.25 | 15 | RCR | 0.145 | 0.323 |  |  |  |
| w2-kd30-constant-1e-6 | 4 | 1.62 | 15 | waves | 0.088 | 0.195 |  |  |  |
| w2-kd30-random-1e-6 | 2 | 6.49 | 16 | diagonal 34/34 | 0.736 | 0.234 | [] | 2 | 0 |
| w2-kd30-random-1e-6 | 3 | 3.25 | 15 | RCR | 0.249 | 0.418 |  |  |  |
| w2-kd30-random-1e-6 | 4 | 1.62 | 15 | waves | 0.076 | 0.418 |  |  |  |
| w2-kd500-random-1e-3 | 2 | 108.25 | 115 | diagonal 34/34 | 0.073 | 0.027 | [] | 0 | 0 |
| w2-kd500-random-1e-3 | 3 | 54.13 | 60 | diagonal 34/34 | 0.095 | 0.047 | [] | 0 | 0 |
| w2-kd500-random-1e-3 | 4 | 27.06 | 32 | diagonal 34/34 | 0.123 | 0.086 | [] | 0 | 0 |
| w2-kd500-random-1e-3 | 5 | 13.53 | 18 | diagonal 34/34 | 0.198 | 0.160 | [] | 0 | 0 |
| w2-kd500-random-1e-3 | 6 | 6.77 | 11 | diagonal 34/34 | 0.116 | 0.170 | [] | 0 | 0 |
| w2-kd500-random-1e-3 | 7 | 3.38 | 7 | Batched | 0.083 | 0.245 |  |  |  |
| w2-kd500-random-1e-3 | 8 | 1.69 | 7 | Batched | 0.001 | 0.020 |  |  |  |
| w2-kd500-random-1e-6 | 2 | 108.25 | 128 | diagonal 34/34 | 0.130 | 0.024 | [] | 0 | 0 |
| w2-kd500-random-1e-6 | 3 | 54.13 | 70 | diagonal 34/34 | 0.286 | 0.074 | [] | 0 | 0 |
| w2-kd500-random-1e-6 | 4 | 27.06 | 40 | diagonal 34/34 | 0.136 | 0.077 | [] | 0 | 0 |
| w2-kd500-random-1e-6 | 5 | 13.53 | 25 | diagonal 34/34 | 0.024 | 0.014 | [] | 0 | 0 |
| w2-kd500-random-1e-6 | 6 | 6.77 | 16 | diagonal 34/34 | 0.285 | 0.327 | [] | 2 | 0 |
| w2-kd500-random-1e-6 | 7 | 3.38 | 15 | diagonal 21/34 | 0.036 | 0.173 | [13, 14, 15, 16, 17] | 4 | 13 |
| w2-kd500-random-1e-6 | 8 | 1.69 | 15 | waves | 0.010 | 0.225 |  |  |  |
| w2-kd500-random-1e-8 | 2 | 108.25 | 134 | diagonal 34/34 | 0.529 | 0.229 | [] | 0 | 0 |
| w2-kd500-random-1e-8 | 3 | 54.13 | 75 | diagonal 34/34 | 0.347 | 0.305 | [] | 2 | 0 |
| w2-kd500-random-1e-8 | 4 | 27.06 | 44 | diagonal 32/34 | 0.211 | 0.320 | [] | 0 | 2 |
| w2-kd500-random-1e-8 | 5 | 13.53 | 32 | diagonal 27/34 | 0.029 | 0.051 | [7, 8] | 4 | 7 |
| w2-kd500-random-1e-8 | 6 | 6.77 | 32 | waves | 0.000 | 0.000 |  |  |  |
| w2-kd500-random-1e-8 | 7 | 3.38 | 32 | waves | 0.000 | 0.000 |  |  |  |
| w2-kd500-random-1e-8 | 8 | 1.69 | 32 | – | – | – |  |  |  |
| w3-kd100-random-1e-3-cells | 2 | 21.65 | 27 | diagonal 34/34 | 0.119 | 0.081 | [] | 0 | 0 |
| w3-kd100-random-1e-3-cells | 3 | 10.83 | 15 | diagonal 34/34 | 0.694 | 0.641 | [] | 0 | 0 |
| w3-kd100-random-1e-3-cells | 4 | 5.41 | 10 | Batched | 0.037 | 0.062 |  |  |  |
| w3-kd100-random-1e-3-cells | 5 | 2.71 | 7 | Batched | 0.014 | 0.061 |  |  |  |
| w3-kd100-random-1e-6-cells | 2 | 21.65 | 34 | diagonal 34/34 | 0.349 | 0.207 | [] | 0 | 0 |
| w3-kd100-random-1e-6-cells | 3 | 10.83 | 21 | diagonal 34/34 | 0.357 | 0.307 | [] | 0 | 0 |
| w3-kd100-random-1e-6-cells | 4 | 5.41 | 15 | waves | 0.178 | 0.323 |  |  |  |
| w3-kd100-random-1e-6-cells | 5 | 2.71 | 15 | waves | 0.012 | 0.084 |  |  |  |
| w3-kd100-random-1e-8-cells | 2 | 21.65 | 38 | diagonal 32/34 | 0.473 | 0.401 | [] | 3 | 2 |
| w3-kd100-random-1e-8-cells | 3 | 10.83 | 32 | waves | 0.412 | 0.483 |  |  |  |
| w3-kd100-random-1e-8-cells | 4 | 5.41 | 32 | waves | 0.000 | 0.000 |  |  |  |
| w3-kd100-random-1e-8-cells | 5 | 2.71 | 32 | waves | 0.000 | 0.000 |  |  |  |
| w3-kd30-random-1e-6 | 2 | 6.50 | 16 | diagonal 34/34 | 0.849 | 0.589 | [] | 2 | 0 |
| w3-kd30-random-1e-6 | 3 | 3.25 | 15 | waves | 0.150 | 0.410 |  |  |  |
| w3-kd30-random-1e-6 | 4 | 1.62 | 15 | – | – | – |  |  |  |
| w3:1000000-kd200-random-1e-3-cells | 2 | 43.30 | 49 | diagonal 34/34 | 0.089 | 0.056 | [] | 0 | 0 |
| w3:1000000-kd200-random-1e-3-cells | 3 | 21.65 | 27 | diagonal 34/34 | 0.110 | 0.080 | [] | 0 | 0 |
| w3:1000000-kd200-random-1e-3-cells | 4 | 10.83 | 15 | waves | 0.691 | 0.759 |  |  |  |
| w3:1000000-kd200-random-1e-3-cells | 5 | 5.41 | 10 | Batched | 0.009 | 0.020 |  |  |  |
| w3:1000000-kd200-random-1e-6-cells | 2 | 43.30 | 58 | diagonal 34/34 | 0.277 | 0.187 | [] | 0 | 0 |
| w3:1000000-kd200-random-1e-6-cells | 3 | 21.65 | 34 | diagonal 34/34 | 0.263 | 0.208 | [] | 0 | 0 |
| w3:1000000-kd200-random-1e-6-cells | 4 | 10.83 | 21 | waves | 0.376 | 0.493 |  |  |  |
| w3:1000000-kd200-random-1e-6-cells | 5 | 5.41 | 15 | diagonal 34/34 | 0.005 | 0.019 | [0, 1] | 3 | 0 |
| w3:1000000-kd200-random-1e-8-cells | 2 | 43.30 | 63 | diagonal 34/34 | 0.257 | 0.129 | [0, 1] | 0 | 0 |
| w3:1000000-kd200-random-1e-8-cells | 3 | 21.65 | 38 | diagonal 32/34 | 0.246 | 0.210 | [] | 3 | 2 |
| w3:1000000-kd200-random-1e-8-cells | 4 | 10.83 | 32 | waves | 0.453 | 0.599 |  |  |  |
| w3:1000000-kd200-random-1e-8-cells | 5 | 5.41 | 32 | – | 0.000 | 0.000 |  |  |  |
| w3:1000000-kd500-random-1e-3-cells | 2 | 108.25 | 115 | diagonal 34/34 | 0.026 | 0.014 | [] | 0 | 0 |
| w3:1000000-kd500-random-1e-3-cells | 3 | 54.13 | 60 | diagonal 34/34 | 0.091 | 0.057 | [] | 0 | 0 |
| w3:1000000-kd500-random-1e-3-cells | 4 | 27.06 | 32 | diagonal 34/34 | 0.341 | 0.283 | [] | 0 | 0 |
| w3:1000000-kd500-random-1e-3-cells | 5 | 13.53 | 18 | diagonal 34/34 | 0.388 | 0.491 | [] | 0 | 0 |
| w3:1000000-kd500-random-1e-6-cells | 2 | 108.25 | 128 | diagonal 34/34 | 0.052 | 0.025 | [] | 0 | 0 |
| w3:1000000-kd500-random-1e-6-cells | 3 | 54.13 | 70 | diagonal 34/34 | 0.250 | 0.186 | [] | 0 | 0 |
| w3:1000000-kd500-random-1e-6-cells | 4 | 27.06 | 40 | diagonal 34/34 | 0.487 | 0.559 | [] | 0 | 0 |
| w3:1000000-kd500-random-1e-6-cells | 5 | 13.53 | 25 | RCR | 0.062 | 0.118 |  |  |  |
| w3:1000000-kd500-random-1e-8-cells | 2 | 108.25 | 134 | diagonal 34/34 | 0.063 | 0.035 | [] | 0 | 0 |
| w3:1000000-kd500-random-1e-8-cells | 3 | 54.13 | 75 | diagonal 34/34 | 0.334 | 0.262 | [] | 2 | 0 |
| w3:1000000-kd500-random-1e-8-cells | 4 | 27.06 | 44 | waves | 0.374 | 0.475 |  |  |  |
| w3:1000000-kd500-random-1e-8-cells | 5 | 13.53 | 32 | – | 0.059 | 0.120 |  |  |  |

### Worst boxes

The worst boxes of ∇φ of W1 at kD ≥ 250 with the levels' parts (at least 0.1) and the M2L runs'
(at least 0.01; none reached it) at their worst targets.

**w1-kd1000-planewave-1e-3**

| # | cell | centre | error / ε | targets | boundary |  | v | / r, levels 2–7 | level parts | M2L parts |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | (127, 63, 64) | (0.996, 0.496, 0.504) | 43.34 | 154 | yes | 1.73 1.73 1.72 1.71 1.69 1.64 | 2: 0.10, 3: 0.20, 4: 0.17, 5: 0.15 |  |
| 2 | (63, 127, 64) | (0.496, 0.996, 0.504) | 37.83 | 155 | yes | 1.73 1.73 1.72 1.71 1.69 1.65 | 3: 0.27, 4: 0.15, 5: 0.13, 6: 0.11 |  |
| 3 | (127, 64, 63) | (0.996, 0.504, 0.496) | 33.21 | 154 | yes | 1.73 1.73 1.72 1.72 1.70 1.67 | 3: 0.21, 4: 0.25, 5: 0.24, 6: 0.17, 7: 0.12 |  |
| 4 | (86, 94, 116) | (0.676, 0.738, 0.910) | 19.13 | 1 | no | 0.99 0.94 0.73 0.99 1.66 1.59 | 6: 1.00 |  |
| 5 | (64, 127, 63) | (0.504, 0.996, 0.496) | 18.30 | 155 | yes | 1.73 1.73 1.73 1.72 1.71 1.69 | 3: 0.23, 4: 0.54, 5: 0.48, 6: 0.30, 7: 0.20 |  |
| 6 | (63, 64, 127) | (0.496, 0.504, 0.996) | 17.41 | 155 | yes | 1.73 1.73 1.72 1.71 1.69 1.64 | 2: 0.14, 3: 0.39, 4: 0.19, 6: 0.11 |  |
| 7 | (127, 64, 62) | (0.996, 0.504, 0.488) | 17.11 | 154 | yes | 1.66 1.60 1.49 1.40 1.70 1.67 | 6: 0.60, 7: 0.32 |  |
| 8 | (103, 19, 88) | (0.809, 0.152, 0.691) | 17.10 | 111 | no | 0.76 0.50 1.38 1.67 1.62 1.53 | 5: 0.52, 6: 0.25, 7: 0.11 |  |

**w1-kd1000-planewave-1e-6**

| # | cell | centre | error / ε | targets | boundary |  | v | / r, levels 2–7 | level parts | M2L parts |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | (127, 63, 64) | (0.996, 0.496, 0.504) | 73.88 | 154 | yes | 1.73 1.73 1.72 1.71 1.69 1.64 | 2: 0.16, 3: 0.23, 4: 0.25, 5: 0.21 |  |
| 2 | (127, 63, 63) | (0.996, 0.496, 0.496) | 40.22 | 156 | yes | 1.73 1.73 1.72 1.72 1.70 1.67 | 3: 0.14, 4: 0.32, 5: 0.47, 6: 0.15, 7: 0.18 |  |
| 3 | (63, 64, 127) | (0.496, 0.504, 0.996) | 32.28 | 155 | yes | 1.73 1.73 1.72 1.71 1.69 1.64 | 2: 0.50, 3: 0.18, 4: 0.23, 5: 0.19 |  |
| 4 | (63, 127, 63) | (0.496, 0.996, 0.496) | 23.08 | 154 | yes | 1.73 1.72 1.72 1.70 1.67 1.61 | 2: 0.14, 3: 0.33, 4: 0.61, 5: 0.38 |  |
| 5 | (64, 64, 127) | (0.504, 0.504, 0.996) | 23.00 | 154 | yes | 1.73 1.72 1.71 1.70 1.66 1.60 | 2: 0.30, 3: 0.26, 4: 0.44, 5: 0.24 |  |
| 6 | (63, 127, 64) | (0.496, 0.996, 0.504) | 22.43 | 155 | yes | 1.73 1.73 1.72 1.71 1.69 1.65 | 2: 0.81, 3: 0.29, 4: 0.26, 5: 0.31, 6: 0.10, 7: 0.13 |  |
| 7 | (63, 63, 127) | (0.496, 0.496, 0.996) | 19.98 | 154 | yes | 1.73 1.72 1.72 1.70 1.67 1.61 | 2: 0.21, 3: 0.60, 4: 0.36, 5: 0.46 |  |
| 8 | (64, 63, 127) | (0.504, 0.496, 0.996) | 17.83 | 156 | yes | 1.73 1.73 1.73 1.72 1.71 1.69 | 2: 0.62, 3: 0.38, 4: 0.53, 5: 0.49, 6: 0.11, 7: 0.18 |  |

**w1-kd1000-planewave-1e-8**

| # | cell | centre | error / ε | targets | boundary |  | v | / r, levels 2–7 | level parts | M2L parts |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | (127, 63, 64) | (0.996, 0.496, 0.504) | 114.80 | 154 | yes | 1.73 1.73 1.72 1.71 1.69 1.64 | 2: 0.14, 3: 0.29, 4: 0.27, 5: 0.22 |  |
| 2 | (127, 63, 63) | (0.996, 0.496, 0.496) | 54.08 | 156 | yes | 1.73 1.73 1.72 1.72 1.70 1.67 | 2: 0.11, 3: 0.16, 4: 0.44, 5: 0.61, 6: 0.18 |  |
| 3 | (63, 63, 127) | (0.496, 0.496, 0.996) | 47.38 | 154 | yes | 1.73 1.72 1.72 1.70 1.67 1.61 | 2: 0.15, 3: 0.73, 4: 0.26, 5: 0.22 |  |
| 4 | (63, 64, 127) | (0.496, 0.504, 0.996) | 36.42 | 155 | yes | 1.73 1.73 1.72 1.71 1.69 1.64 | 2: 0.45, 3: 0.49, 4: 0.33, 5: 0.20 |  |
| 5 | (63, 127, 63) | (0.496, 0.996, 0.496) | 35.05 | 154 | yes | 1.73 1.72 1.72 1.70 1.67 1.61 | 2: 0.13, 3: 0.47, 4: 0.80, 5: 0.42 |  |
| 6 | (63, 127, 64) | (0.496, 0.996, 0.504) | 31.89 | 155 | yes | 1.73 1.73 1.72 1.71 1.69 1.65 | 2: 1.15, 3: 0.35, 4: 0.19, 5: 0.24 |  |
| 7 | (64, 64, 127) | (0.504, 0.504, 0.996) | 31.67 | 154 | yes | 1.73 1.72 1.71 1.70 1.66 1.60 | 2: 0.51, 3: 0.49, 4: 0.38, 5: 0.36 |  |
| 8 | (127, 64, 63) | (0.996, 0.504, 0.496) | 19.66 | 154 | yes | 1.73 1.73 1.72 1.72 1.70 1.67 | 2: 0.35, 3: 1.20, 4: 1.72, 5: 1.58, 6: 0.41 |  |

**w1-kd1000-random-1e-3**

| # | cell | centre | error / ε | targets | boundary |  | v | / r, levels 2–7 | level parts | M2L parts |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | (63, 0, 63) | (0.496, 0.004, 0.496) | 20.25 | 153 | yes | 1.73 1.73 1.72 1.72 1.70 1.67 | 2: 0.17, 3: 0.20, 4: 0.15, 5: 0.26, 6: 0.23, 7: 0.13 |  |
| 2 | (63, 64, 0) | (0.496, 0.504, 0.004) | 19.05 | 154 | yes | 1.73 1.73 1.72 1.72 1.70 1.67 | 2: 0.13, 3: 0.26, 4: 0.17, 5: 0.29, 6: 0.18, 7: 0.16 |  |
| 3 | (64, 63, 127) | (0.504, 0.496, 0.996) | 17.84 | 156 | yes | 1.73 1.73 1.73 1.72 1.71 1.69 | 2: 0.18, 3: 0.28, 4: 0.28, 5: 0.11, 6: 0.29, 7: 0.22 |  |
| 4 | (0, 63, 64) | (0.004, 0.496, 0.504) | 16.77 | 156 | yes | 1.73 1.73 1.72 1.71 1.69 1.66 | 2: 0.20, 3: 0.40, 4: 0.30, 5: 0.30, 6: 0.19, 7: 0.18 |  |
| 5 | (63, 63, 127) | (0.496, 0.496, 0.996) | 14.92 | 154 | yes | 1.73 1.72 1.72 1.70 1.67 1.61 | 2: 0.15, 3: 0.34, 4: 0.13, 5: 0.33, 6: 0.18, 7: 0.20 |  |
| 6 | (64, 127, 63) | (0.504, 0.996, 0.496) | 13.58 | 155 | yes | 1.73 1.73 1.73 1.72 1.71 1.69 | 2: 0.43, 3: 0.19, 4: 0.37, 5: 0.36, 6: 0.35, 7: 0.27 |  |
| 7 | (63, 64, 127) | (0.496, 0.504, 0.996) | 13.23 | 155 | yes | 1.73 1.73 1.72 1.71 1.69 1.64 | 2: 0.19, 3: 0.56, 4: 0.27, 5: 0.48, 6: 0.35, 7: 0.12 |  |
| 8 | (0, 64, 64) | (0.004, 0.504, 0.504) | 12.65 | 154 | yes | 1.73 1.72 1.72 1.70 1.67 1.61 | 2: 0.14, 3: 0.32, 4: 0.12, 5: 0.32, 6: 0.32, 7: 0.10 |  |

**w1-kd1000-random-1e-6**

| # | cell | centre | error / ε | targets | boundary |  | v | / r, levels 2–7 | level parts | M2L parts |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | (64, 0, 64) | (0.504, 0.004, 0.504) | 27.99 | 155 | yes | 1.73 1.73 1.72 1.71 1.70 1.66 | 2: 0.69, 3: 0.15, 4: 0.29, 5: 0.13, 7: 0.20 |  |
| 2 | (64, 127, 63) | (0.504, 0.996, 0.496) | 26.03 | 155 | yes | 1.73 1.73 1.73 1.72 1.71 1.69 | 2: 0.37, 4: 0.50, 5: 0.25, 6: 0.13, 7: 0.25 |  |
| 3 | (64, 64, 0) | (0.504, 0.504, 0.004) | 25.88 | 154 | yes | 1.73 1.72 1.72 1.70 1.67 1.62 | 2: 0.38, 3: 0.40, 4: 0.21, 5: 0.27 |  |
| 4 | (64, 63, 127) | (0.504, 0.496, 0.996) | 24.02 | 156 | yes | 1.73 1.73 1.73 1.72 1.71 1.69 | 2: 0.58, 3: 0.25, 4: 0.46, 5: 0.33, 6: 0.14, 7: 0.35 |  |
| 5 | (0, 63, 63) | (0.004, 0.496, 0.496) | 21.89 | 154 | yes | 1.73 1.72 1.72 1.70 1.67 1.61 | 2: 0.41, 3: 0.39, 4: 0.38, 5: 0.17, 7: 0.12 |  |
| 6 | (127, 63, 64) | (0.996, 0.496, 0.504) | 20.45 | 154 | yes | 1.73 1.73 1.72 1.71 1.69 1.64 | 2: 0.34, 3: 0.10, 4: 0.44, 5: 0.32, 7: 0.13 |  |
| 7 | (0, 64, 63) | (0.004, 0.504, 0.496) | 19.78 | 155 | yes | 1.73 1.73 1.73 1.72 1.71 1.68 | 2: 0.82, 3: 0.29, 4: 0.40, 5: 0.56, 6: 0.10, 7: 0.46 |  |
| 8 | (64, 127, 64) | (0.504, 0.996, 0.504) | 19.00 | 155 | yes | 1.73 1.72 1.72 1.70 1.67 1.61 | 2: 0.41, 3: 0.42, 4: 0.30, 5: 0.25, 7: 0.14 |  |

**w1-kd1000-random-1e-8**

| # | cell | centre | error / ε | targets | boundary |  | v | / r, levels 2–7 | level parts | M2L parts |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | (63, 64, 127) | (0.496, 0.504, 0.996) | 43.30 | 155 | yes | 1.73 1.73 1.72 1.71 1.69 1.64 | 2: 0.24, 3: 0.36, 4: 0.37, 5: 0.39 |  |
| 2 | (63, 127, 64) | (0.496, 0.996, 0.504) | 38.52 | 155 | yes | 1.73 1.73 1.72 1.71 1.69 1.65 | 2: 0.20, 3: 0.45, 4: 0.28, 5: 0.30 |  |
| 3 | (64, 0, 63) | (0.504, 0.004, 0.496) | 35.14 | 155 | yes | 1.73 1.73 1.72 1.71 1.68 1.64 | 2: 0.26, 3: 0.38, 4: 0.32, 5: 0.32 |  |
| 4 | (63, 64, 0) | (0.496, 0.504, 0.004) | 34.63 | 154 | yes | 1.73 1.73 1.72 1.72 1.70 1.67 | 2: 0.29, 3: 0.37, 4: 0.40, 5: 0.39 |  |
| 5 | (63, 127, 63) | (0.496, 0.996, 0.496) | 34.20 | 154 | yes | 1.73 1.72 1.72 1.70 1.67 1.61 | 2: 0.49, 3: 0.28, 4: 0.34, 5: 0.21 |  |
| 6 | (64, 63, 127) | (0.504, 0.496, 0.996) | 32.96 | 156 | yes | 1.73 1.73 1.73 1.72 1.71 1.69 | 2: 0.75, 3: 0.64, 4: 0.26, 5: 0.48 |  |
| 7 | (0, 64, 63) | (0.004, 0.504, 0.496) | 32.74 | 155 | yes | 1.73 1.73 1.73 1.72 1.71 1.68 | 2: 0.28, 3: 0.58, 4: 0.43, 5: 0.59, 6: 0.12 |  |
| 8 | (64, 63, 0) | (0.504, 0.496, 0.004) | 32.67 | 155 | yes | 1.73 1.73 1.73 1.72 1.71 1.68 | 2: 0.73, 3: 0.47, 4: 0.46, 5: 0.36, 6: 0.17 |  |

**w1-kd250-random-1e-3**

| # | cell | centre | error / ε | targets | boundary |  | v | / r, levels 2–5 | level parts | M2L parts |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | (15, 15, 31) | (0.484, 0.484, 0.984) | 13.54 | 154 | yes | 1.72 1.70 1.67 1.61 | 2: 0.23, 3: 0.41, 4: 0.21, 5: 0.17 |  |
| 2 | (0, 15, 15) | (0.016, 0.484, 0.484) | 11.44 | 155 | yes | 1.72 1.71 1.69 1.64 | 2: 0.34, 3: 0.43, 4: 0.46, 5: 0.14 |  |
| 3 | (15, 0, 16) | (0.484, 0.016, 0.516) | 11.15 | 154 | yes | 1.72 1.71 1.69 1.65 | 2: 0.23, 3: 0.37, 4: 0.34, 5: 0.20 |  |
| 4 | (15, 16, 0) | (0.484, 0.516, 0.016) | 11.01 | 154 | yes | 1.72 1.71 1.68 1.64 | 2: 0.38, 3: 0.21, 4: 0.31, 5: 0.20 |  |
| 5 | (31, 16, 15) | (0.984, 0.516, 0.484) | 10.99 | 156 | yes | 1.72 1.71 1.70 1.66 | 2: 0.28, 3: 0.18, 4: 0.49, 5: 0.21 |  |
| 6 | (31, 16, 16) | (0.984, 0.516, 0.516) | 10.74 | 155 | yes | 1.72 1.71 1.68 1.64 | 2: 0.26, 3: 0.26, 4: 0.34, 5: 0.21 |  |
| 7 | (0, 16, 15) | (0.016, 0.516, 0.484) | 10.72 | 155 | yes | 1.72 1.71 1.68 1.63 | 2: 0.37, 3: 0.39, 4: 0.30, 5: 0.17 |  |
| 8 | (0, 16, 16) | (0.016, 0.516, 0.516) | 10.42 | 154 | yes | 1.73 1.72 1.72 1.70 | 2: 0.42, 3: 0.41, 4: 0.36, 5: 0.18 |  |

**w1-kd250-random-1e-6**

| # | cell | centre | error / ε | targets | boundary |  | v | / r, levels 2–5 | level parts | M2L parts |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | (31, 15, 16) | (0.984, 0.484, 0.516) | 16.90 | 156 | yes | 1.73 1.72 1.71 1.68 | 2: 0.32, 3: 0.50, 5: 0.32 |  |
| 2 | (0, 16, 16) | (0.016, 0.516, 0.516) | 14.77 | 154 | yes | 1.73 1.72 1.72 1.70 | 2: 0.31, 3: 0.56, 4: 0.16, 5: 0.90 |  |
| 3 | (15, 31, 15) | (0.484, 0.984, 0.484) | 13.63 | 154 | yes | 1.73 1.73 1.72 1.71 | 2: 0.44, 3: 0.48, 4: 0.22, 5: 0.74 |  |
| 4 | (15, 0, 16) | (0.484, 0.016, 0.516) | 13.15 | 154 | yes | 1.72 1.71 1.69 1.65 | 2: 0.33, 3: 0.63, 4: 0.17, 5: 0.54 |  |
| 5 | (31, 16, 15) | (0.984, 0.516, 0.484) | 12.58 | 156 | yes | 1.72 1.71 1.70 1.66 | 2: 0.59, 3: 0.31, 4: 0.15, 5: 0.31 |  |
| 6 | (0, 15, 15) | (0.016, 0.484, 0.484) | 11.97 | 155 | yes | 1.72 1.71 1.69 1.64 | 2: 0.42, 3: 0.78, 4: 0.14, 5: 0.33 |  |
| 7 | (15, 15, 0) | (0.484, 0.484, 0.016) | 10.79 | 156 | yes | 1.71 1.70 1.66 1.60 | 2: 0.39, 3: 0.39, 4: 0.14, 5: 0.33 |  |
| 8 | (16, 0, 15) | (0.516, 0.016, 0.484) | 9.96 | 156 | yes | 1.72 1.70 1.67 1.61 | 2: 0.35, 3: 0.49, 4: 0.13, 5: 0.40 |  |

**w1-kd250-random-1e-8**

| # | cell | centre | error / ε | targets | boundary |  | v | / r, levels 2–5 | level parts | M2L parts |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | (16, 0, 16) | (0.516, 0.016, 0.516) | 21.45 | 155 | yes | 1.73 1.72 1.71 1.68 | 2: 0.45, 3: 0.66, 4: 0.14 |  |
| 2 | (31, 16, 15) | (0.984, 0.516, 0.484) | 21.31 | 156 | yes | 1.72 1.71 1.70 1.66 | 2: 0.35, 3: 0.56, 4: 0.16 |  |
| 3 | (0, 15, 16) | (0.016, 0.484, 0.516) | 19.08 | 154 | yes | 1.71 1.70 1.66 1.60 | 2: 0.37, 3: 0.60 |  |
| 4 | (31, 15, 16) | (0.984, 0.484, 0.516) | 18.76 | 156 | yes | 1.73 1.72 1.71 1.68 | 2: 0.62, 3: 0.53, 4: 0.17 |  |
| 5 | (16, 16, 0) | (0.516, 0.516, 0.016) | 18.46 | 155 | yes | 1.73 1.72 1.71 1.69 | 2: 0.66, 3: 0.54, 4: 0.21 |  |
| 6 | (15, 0, 16) | (0.484, 0.016, 0.516) | 18.37 | 154 | yes | 1.72 1.70 1.67 1.61 | 2: 0.72, 3: 0.32 |  |
| 7 | (15, 16, 31) | (0.484, 0.516, 0.984) | 17.44 | 155 | yes | 1.72 1.71 1.69 1.64 | 2: 0.41, 3: 0.58, 4: 0.20 |  |
| 8 | (31, 16, 16) | (0.984, 0.516, 0.516) | 17.38 | 155 | yes | 1.72 1.71 1.68 1.64 | 2: 0.80, 3: 0.32, 4: 0.23 |  |

**w1-kd500-planewave-1e-3**

| # | cell | centre | error / ε | targets | boundary |  | v | / r, levels 2–6 | level parts | M2L parts |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | (63, 31, 32) | (0.992, 0.492, 0.508) | 39.61 | 155 | yes | 1.73 1.72 1.71 1.69 1.65 | 2: 0.12, 3: 0.19, 4: 0.22, 5: 0.14 |  |
| 2 | (63, 32, 31) | (0.992, 0.508, 0.492) | 26.93 | 156 | yes | 1.73 1.72 1.72 1.70 1.67 | 3: 0.23, 4: 0.33, 5: 0.24, 6: 0.14 |  |
| 3 | (31, 63, 32) | (0.492, 0.992, 0.508) | 24.44 | 155 | yes | 1.72 1.72 1.70 1.67 1.60 | 2: 0.30, 3: 0.23, 4: 0.18, 5: 0.10 |  |
| 4 | (63, 32, 32) | (0.992, 0.508, 0.508) | 16.01 | 154 | yes | 1.72 1.72 1.70 1.67 1.62 | 2: 0.15, 3: 0.37, 4: 0.59, 5: 0.42, 6: 0.17 |  |
| 5 | (31, 32, 63) | (0.492, 0.508, 0.992) | 16.00 | 155 | yes | 1.73 1.72 1.71 1.69 1.64 | 3: 0.32, 4: 0.31, 5: 0.18 |  |
| 6 | (32, 63, 32) | (0.508, 0.992, 0.508) | 15.81 | 155 | yes | 1.72 1.72 1.70 1.68 1.62 | 2: 0.34, 3: 0.56, 4: 0.60, 5: 0.32, 6: 0.18 |  |
| 7 | (31, 63, 31) | (0.492, 0.992, 0.492) | 15.39 | 155 | yes | 1.73 1.73 1.73 1.72 1.70 | 3: 0.36, 4: 0.84, 5: 0.72, 6: 0.53 |  |
| 8 | (31, 31, 63) | (0.492, 0.492, 0.992) | 15.27 | 154 | yes | 1.72 1.72 1.70 1.67 1.61 | 2: 0.22, 3: 0.52, 4: 0.22, 5: 0.14, 6: 0.19 |  |

**w1-kd500-planewave-1e-6**

| # | cell | centre | error / ε | targets | boundary |  | v | / r, levels 2–6 | level parts | M2L parts |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | (63, 31, 32) | (0.992, 0.492, 0.508) | 42.42 | 155 | yes | 1.73 1.72 1.71 1.69 1.66 | 2: 0.13, 3: 0.30, 4: 0.32, 6: 0.12 |  |
| 2 | (31, 63, 32) | (0.492, 0.992, 0.508) | 40.59 | 155 | yes | 1.72 1.72 1.70 1.67 1.60 | 2: 0.25, 3: 0.33, 4: 0.28 |  |
| 3 | (31, 63, 31) | (0.492, 0.992, 0.492) | 26.20 | 155 | yes | 1.73 1.73 1.73 1.72 1.70 | 2: 0.12, 3: 0.41, 4: 0.68, 5: 0.24, 6: 0.28 |  |
| 4 | (31, 31, 63) | (0.492, 0.492, 0.992) | 23.01 | 154 | yes | 1.72 1.72 1.70 1.67 1.61 | 2: 0.24, 3: 0.81, 4: 0.19, 5: 0.12, 6: 0.13 |  |
| 5 | (31, 32, 63) | (0.492, 0.508, 0.992) | 19.53 | 155 | yes | 1.73 1.72 1.71 1.69 1.64 | 2: 0.50, 3: 0.40, 4: 0.28, 6: 0.11 |  |
| 6 | (63, 32, 31) | (0.992, 0.508, 0.492) | 18.52 | 156 | yes | 1.73 1.72 1.72 1.70 1.67 | 2: 0.15, 3: 0.60, 4: 0.81, 5: 0.25, 6: 0.36 |  |
| 7 | (32, 31, 63) | (0.508, 0.492, 0.992) | 16.52 | 156 | yes | 1.73 1.73 1.72 1.71 1.70 | 2: 0.48, 3: 0.67, 4: 0.30, 6: 0.26 |  |
| 8 | (32, 32, 63) | (0.508, 0.508, 0.992) | 16.31 | 154 | yes | 1.72 1.71 1.70 1.66 1.60 | 2: 0.43, 3: 0.51, 4: 0.44, 5: 0.16 |  |

**w1-kd500-planewave-1e-8**

| # | cell | centre | error / ε | targets | boundary |  | v | / r, levels 2–6 | level parts | M2L parts |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | (31, 63, 32) | (0.492, 0.992, 0.508) | 68.21 | 155 | yes | 1.72 1.72 1.70 1.67 1.60 | 2: 0.30, 3: 0.33, 4: 0.28 |  |
| 2 | (31, 63, 31) | (0.492, 0.992, 0.492) | 56.37 | 155 | yes | 1.73 1.73 1.73 1.72 1.70 | 2: 0.10, 3: 0.29, 4: 0.58, 5: 0.18 |  |
| 3 | (63, 31, 32) | (0.992, 0.492, 0.508) | 54.96 | 155 | yes | 1.73 1.72 1.71 1.69 1.66 | 2: 0.20, 3: 0.38, 4: 0.37, 5: 0.10 |  |
| 4 | (31, 32, 63) | (0.492, 0.508, 0.992) | 30.47 | 155 | yes | 1.73 1.72 1.71 1.69 1.64 | 2: 0.41, 3: 0.50, 4: 0.47 |  |
| 5 | (31, 31, 63) | (0.492, 0.492, 0.992) | 30.16 | 154 | yes | 1.72 1.72 1.70 1.67 1.61 | 2: 0.17, 3: 0.72, 4: 0.30, 5: 0.13 |  |
| 6 | (63, 32, 31) | (0.992, 0.508, 0.492) | 25.21 | 156 | yes | 1.73 1.72 1.72 1.70 1.67 | 2: 0.22, 3: 0.65, 4: 0.91, 5: 0.26 |  |
| 7 | (32, 32, 63) | (0.508, 0.508, 0.992) | 22.40 | 154 | yes | 1.72 1.71 1.70 1.66 1.60 | 2: 0.50, 3: 0.75, 4: 0.50, 5: 0.16 |  |
| 8 | (32, 31, 63) | (0.508, 0.492, 0.992) | 20.17 | 156 | yes | 1.73 1.73 1.72 1.71 1.70 | 2: 0.70, 3: 0.92, 4: 0.55, 5: 0.17 |  |

**w1-kd500-random-1e-3**

| # | cell | centre | error / ε | targets | boundary |  | v | / r, levels 2–6 | level parts | M2L parts |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | (31, 31, 63) | (0.492, 0.492, 0.992) | 23.01 | 154 | yes | 1.72 1.72 1.70 1.67 1.61 | 2: 0.28, 3: 0.32, 4: 0.17, 5: 0.11 |  |
| 2 | (31, 63, 31) | (0.492, 0.992, 0.492) | 14.28 | 155 | yes | 1.73 1.73 1.73 1.72 1.70 | 2: 0.28, 3: 0.23, 4: 0.17, 5: 0.34, 6: 0.28 |  |
| 3 | (31, 32, 0) | (0.492, 0.508, 0.008) | 14.22 | 154 | yes | 1.72 1.71 1.69 1.66 1.58 | 2: 0.31, 3: 0.15, 4: 0.22, 5: 0.14 |  |
| 4 | (31, 31, 0) | (0.492, 0.492, 0.008) | 13.87 | 155 | yes | 1.73 1.72 1.71 1.69 1.65 | 2: 0.14, 3: 0.34, 4: 0.23, 5: 0.29, 6: 0.13 |  |
| 5 | (32, 63, 32) | (0.508, 0.992, 0.508) | 13.37 | 155 | yes | 1.73 1.72 1.70 1.68 1.63 | 2: 0.28, 3: 0.36, 4: 0.34, 5: 0.18, 6: 0.15 |  |
| 6 | (63, 31, 31) | (0.992, 0.492, 0.492) | 12.17 | 154 | yes | 1.72 1.72 1.70 1.67 1.61 | 2: 0.54, 3: 0.28, 4: 0.31, 5: 0.26, 6: 0.15 |  |
| 7 | (31, 0, 32) | (0.492, 0.008, 0.508) | 11.82 | 154 | yes | 1.73 1.73 1.72 1.71 1.69 | 2: 0.36, 3: 0.35, 4: 0.32, 5: 0.39, 6: 0.43 |  |
| 8 | (0, 31, 32) | (0.008, 0.492, 0.508) | 11.75 | 155 | yes | 1.73 1.73 1.73 1.72 1.71 | 2: 0.32, 3: 0.27, 4: 0.41, 5: 0.20, 6: 0.24 |  |

**w1-kd500-random-1e-6**

| # | cell | centre | error / ε | targets | boundary |  | v | / r, levels 2–6 | level parts | M2L parts |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | (31, 31, 0) | (0.492, 0.492, 0.008) | 20.03 | 155 | yes | 1.73 1.72 1.71 1.69 1.65 | 2: 0.18, 3: 0.31, 4: 0.31, 5: 0.11, 6: 0.40 |  |
| 2 | (31, 0, 32) | (0.492, 0.008, 0.508) | 19.91 | 154 | yes | 1.73 1.73 1.72 1.71 1.69 | 2: 0.29, 3: 0.27, 4: 0.59, 5: 0.11, 6: 0.23 |  |
| 3 | (0, 31, 32) | (0.008, 0.492, 0.508) | 19.89 | 155 | yes | 1.73 1.73 1.73 1.72 1.71 | 2: 0.29, 3: 0.47, 4: 0.67, 5: 0.12, 6: 0.30 |  |
| 4 | (32, 31, 63) | (0.508, 0.492, 0.992) | 19.85 | 156 | yes | 1.73 1.73 1.72 1.71 1.70 | 2: 0.36, 3: 0.38, 4: 0.45, 5: 0.17, 6: 0.32 |  |
| 5 | (63, 31, 32) | (0.992, 0.492, 0.508) | 18.79 | 155 | yes | 1.73 1.72 1.71 1.69 1.65 | 2: 0.28, 3: 0.38, 4: 0.45, 6: 0.27 |  |
| 6 | (63, 32, 32) | (0.992, 0.508, 0.508) | 17.14 | 154 | yes | 1.72 1.72 1.70 1.67 1.62 | 2: 0.38, 3: 0.29, 4: 0.38, 6: 0.17 |  |
| 7 | (63, 32, 31) | (0.992, 0.508, 0.492) | 15.05 | 156 | yes | 1.73 1.72 1.72 1.70 1.67 | 2: 0.26, 3: 0.11, 4: 0.45, 6: 0.48 |  |
| 8 | (32, 0, 31) | (0.508, 0.008, 0.492) | 14.52 | 154 | yes | 1.73 1.72 1.71 1.68 1.64 | 2: 0.44, 3: 0.50, 4: 0.53, 6: 0.32 |  |

**w1-kd500-random-1e-8**

| # | cell | centre | error / ε | targets | boundary |  | v | / r, levels 2–6 | level parts | M2L parts |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | (32, 31, 63) | (0.508, 0.492, 0.992) | 40.47 | 156 | yes | 1.73 1.73 1.72 1.71 1.70 | 2: 0.29, 3: 0.51, 4: 0.44 |  |
| 2 | (63, 32, 31) | (0.992, 0.508, 0.492) | 29.12 | 156 | yes | 1.73 1.72 1.72 1.70 1.67 | 2: 0.31, 3: 0.48, 4: 0.44, 5: 0.13 |  |
| 3 | (63, 31, 32) | (0.992, 0.492, 0.508) | 27.51 | 155 | yes | 1.73 1.72 1.71 1.69 1.66 | 2: 0.55, 3: 0.42, 4: 0.32, 5: 0.17 |  |
| 4 | (0, 31, 32) | (0.008, 0.492, 0.508) | 27.18 | 155 | yes | 1.73 1.73 1.73 1.72 1.71 | 2: 0.72, 3: 0.36, 4: 0.36, 5: 0.15 |  |
| 5 | (31, 63, 31) | (0.492, 0.992, 0.492) | 24.19 | 155 | yes | 1.73 1.73 1.73 1.72 1.70 | 2: 0.30, 3: 0.40, 4: 0.45, 5: 0.27 |  |
| 6 | (31, 31, 63) | (0.492, 0.492, 0.992) | 20.53 | 154 | yes | 1.72 1.72 1.70 1.67 1.61 | 2: 0.39, 3: 0.45, 4: 0.57, 5: 0.12 |  |
| 7 | (31, 0, 32) | (0.492, 0.008, 0.508) | 19.91 | 154 | yes | 1.72 1.71 1.69 1.64 1.56 | 2: 0.55, 3: 0.41, 4: 0.22 |  |
| 8 | (0, 31, 31) | (0.008, 0.492, 0.492) | 19.77 | 154 | yes | 1.72 1.72 1.70 1.67 1.60 | 2: 0.33, 3: 0.59, 4: 0.38 |  |
