# Helmholtz CH4.3d1: the transfer function's length, per pair

Helmholtz task CH4.3d1 ([docs/helmholtz/tasks/CH4.3d1.md](../../docs/helmholtz/tasks/CH4.3d1.md));
proposal [§3.6](../../docs/design/helmholtz-fmm-proposal.md) (the diagonal form, its rounding, "a
shorter T trading truncation for stability; CH4.3's to analyse"), §5.4 (the grid) and §10 CH4.3;
review 2026-11-12, M8 ([docs/reviews/2026-11-12-H6.md](../../docs/reviews/2026-11-12-H6.md) §1);
S6's report §6 ([spikes/diagonal/SPIKE_REPORT.md](../../spikes/diagonal/SPIKE_REPORT.md)).
Measured on 2026-11-22 on locust (GH200, 72 Grace cores, 70 threads; idle, load 0.04 before the
run) from branch `tbetcke/h-ch4-3d1`, release build, rustc 1.99.0; the tests and a smoke run on the
M6. Accuracy only: nothing here is timed.

## Summary

1. **2P − 2 stays; CH4.3d2 is not recommended.** A shorter transfer function is only worth a
   few per cent of a diagonal level's grid, and no length rule that can be computed before the
   run reaches even that safely.
2. **What the shortest admissible lengths would give** (§3; an oracle, read off the
   measurement): for real k the diagonal form stays within RCR's error down to 1–10 degrees
   below 2p at k a 22.6–32 and 6–28 below at k a ≥ 45. The grid that length allows has 0–16 %
   fewer samples than 2p's on W1's diagonal levels (3–14 % at k a ≥ 22.6). At k a 11.3–16 and
   10⁻⁶–10⁻⁷ some class needs 2p itself (no gain); at 10⁻⁸ there the nearest classes are on RCR
   and the others save up to 7 %. At 10⁻⁸
   a shorter T takes two classes from RCR on one level of W1 (k a = 27.1: 6 of an interior box's
   316 offsets) and none elsewhere. M8's expectation, "a few degrees fewer could move the nearest
   classes from RCR at 10⁻⁸", holds only there.
3. **The named lengths** (§4): ⌈3(P − 1)/2⌉ is within RCR's error on no class at k a ≥ 11.3 (S6
   saw it fail at high k a; at the FMM's degrees it fails everywhere the rule takes classes).
   The excess-bandwidth length 2ka + 1.8 (log₁₀ 1/ε)^{2/3} (2ka)^{1/3} (Song and Chew's constant)
   is admissible on every class for real k at k a ≥ 45. There it is at most 9 degrees below
   2p (0–8 % of the samples). At k a 11–32 it fails on 1–4 of the classes the bound takes, with
   up to 7× RCR's error per pair on the nearest class. The length needed is not of the form
   2ka + γ (log₁₀ 1/ε)^{2/3} (2ka)^{1/3} for one γ: the γ needed runs from 0.1 to 2.3 over k a
   and ε.
4. **A truncation model** (§5): Gegenbauer's tail at the pair's largest radius,
   τ_L = |k||t| Σ_{L<l≤2p} (2l + 1)|jₗ(2ka)||hₗ(k|t|)|, covers T's measured truncation with
   C_τ = 0.041 (the largest at 2×). It over-predicts by 67× in the median and up to 10⁸×, and
   more at larger k a. The rule it gives (the shortest L with 2 C_τ τ_L ≤ 0.1 ε) is safe on every
   one of the 4,356 classes the bound takes at k a ≥ 11. But it returns 2p on all but eight of
   134 rows: 0–3 % fewer samples. A model that knows the multipole's degree (the truncated
   expansions' coefficients, not the full Gegenbauer tail) might reach the oracle's 0–16 %. That
   is the most a further study could win.
5. **The bound holds at the shorter lengths** (§6): the f64 chain's rounding against the
   double-double chain stays within 2E (E at length L) on every row at L = 2p (largest 1.99 E,
   S6's own fitting point). At L < 2p six rows at k a 2.83, the rule's band edge, reach 2.0–2.7 E.
   At k a 4–8 the rounding is at most 1.1 E, at k a ≥ 11 at most 0.4 E.

## 1. The scan

S6's spike (`spikes/diagonal`) got a scan of T's length, `src/lengths.rs` and
`examples/lengths.rs`. Per unit (kind, k a, p), with p the FMM's degree (`DegreeRule::CALIBRATED`,
CH4.3c2's refit; the spike's copy `rules::refit` is pinned to the product's values by a test) at
ε = 10⁻³ … 10⁻⁸ below §5.5's pruning cut, every one of the 34 classes and S7's three placements
(uniform, corners, W1 patches: 128 sources, 64 targets), and the lengths L from
max(p, min(⌈3p/2⌉, ⌊2ka⌋ − ⌈(2ka)^{1/3}⌉)) to 2p in at most 25 steps (steps of 1–6 degrees, 6 at
p ≈ 250), with ⌈3p/2⌉ and the excess-bandwidth lengths of the unit's ε added exactly:

- **trunc**: S7's RCR at degree p against the direct sum, both in double-double;
- **dd**: the chain in double-double at length L against double-double RCR. This is T's
  truncation alone; it is computed at p ≤ 100 and is 10⁻²⁴–10⁻³¹ at L = 2p;
- **round**: the f64 chain (T rounded once from double-double, as the FMM builds it) against the
  double-double chain at the same L;
- **dev**, **total**: the f64 chain against double-double RCR and against the direct sum.

Every number is the worst over the placements of max(potential, gradient), S7's relative ℓ2 norm
over the targets. A length is **admissible** for a class at (k a, ε) when on every placement and
both norms the f64 chain's total is at most max(1.1 × RCR's, RCR's + 0.1 ε): as accurate as RCR
within CH4.3b's 1.1×, or within a tenth of ε more where RCR's own error is far below ε. The kinds are real k and k = κ(1 + 0.1i).
The k a values are S6's grid 2^{j/2} (j = 3..15, 2.83–181) plus W1's levels at kD = 500 and 1000
(108.3, 54.1, 27.1, 13.5 and 216.5). That makes 194 units and 387,804 rows, 26 minutes on locust
at 70 threads. The rows are archived in `~/projects/fmm-archive/locust/runs/ch4-3d1/` (62 MB, not
committed); `tools/helmholtz/ch4-3d1-tables.py` makes every table here from them.

Tests (`spikes/diagonal/tests/lengths.rs`, release, all pass):

- `the_scan_at_2p_is_s6s_chain`: at L = 2p the scan's f64 chain is S6's `f64:2P-2` variant bit for
  bit (its errors equal S6's `round` and `total`, RCR's truncation equal), on the 34 classes,
  real and damped, k a 4, p 10.
- `the_truncation_is_resolved`: the double-double chain at 2p within 10⁻²⁰ of double-double RCR
  (S6's identity through the scan; measured 10⁻²⁴–10⁻³¹), and T's truncation falling by at least
  two orders from the shortest length to 2p − 1 on the nearest and the farthest class. On the
  nearest it falls slowly, about Gegenbauer's (2a/|t|)^l = 0.87^l: from 1.1·10⁻³ to 7.6·10⁻⁸ at
  k a 8, p 18, and dropping only the last degree costs 7.6·10⁻⁸ to 6.9·10⁻⁶ over the placements.
- `the_tail_matches_double_double`: the model's sum in f64 against double-double within 10⁻¹⁰ in
  log₁₀, and finite where hₗ leaves f64's range.
- `the_lengths`, `the_rule_is_the_fmms` (the copy against `tetraleaf-exec`'s degrees at x = 1–216.5,
  10⁻³–10⁻⁸), `the_scan_is_deterministic` (1 and 3 threads, bit for bit).

## 2. What a shorter T costs

T's length L drops the (S|R) terms of degree n + n′ > L. Their weight decays like Gegenbauer's
series at the pair's largest radius 2a against |t|: geometrically, as (2a/|t|)^l = 0.87^l on the
nearest class (|t| = 2w, 2a = √3 w), once l passes 2ka, and faster on the far classes. At small k a
the decay is slow, so the last degrees below 2p carry errors of the size of RCR's own. At large k a,
where 2p is far above 2ka (the degree's excess over k a doubled), there is room below 2p.
Measured on the nearest class (|t|² = 4), worst placement: T's truncation (the double-double
chain against RCR at p ≤ 100; above, the f64 chain against RCR, whose rounding is far smaller
there) at a few lengths, beside RCR's own error against the direct sum (`ch4-3d1-tables.py`'s
data; the scan's steps are 1–3 degrees here):

| k a | ε | p | RCR's error | T's truncation | shortest admissible |
| --- | --- | --- | --- | --- | --- |
| 8.0 | 10⁻⁶ | 18 | 4.6·10⁻⁶ | 6.9·10⁻⁶ at 2p − 1, 1.8·10⁻⁵ at 2p − 2 | 2p (36) |
| 16.0 | 10⁻⁶ | 28 | 9.9·10⁻⁷ | 9.6·10⁻⁷ at 2p − 2 | 2p (56) |
| 27.1 | 10⁻⁶ | 41 | 5.7·10⁻⁷ | 4.0·10⁻⁸ at 2p − 2, 3.7·10⁻⁷ at 2p − 8 | 2p − 8 (74) |
| 45.3 | 10⁻⁶ | 61 | 8.0·10⁻⁷ | 2.6·10⁻¹⁰ at 2p − 1, 1.1·10⁻⁷ at 2p − 11 | 2p − 11 (111) |
| 54.1 | 10⁻⁸ | 76 | 3.7·10⁻⁹ | 4.1·10⁻¹² at 2p − 1, 2.1·10⁻⁹ at 2p − 16 (f64) | 2p − 16 (136) |
| 108.3 | 10⁻⁸ | 134 | 4.5·10⁻⁹ | (f64) 6.9·10⁻¹⁴ at 2p − 7, 9.0·10⁻¹² at 2p − 16 | 2p − 25 (243) |

## 3. The shortest admissible length, and what it would save

Per (kind, k a, ε): the shortest admissible length of the nearest class and the longest over the
classes (which sets a common grid), against 2p (from `ch4-3d1-tables.py lmin`; real k, W1's levels):

| k a | ε | p | 2p | nearest | longest | longest − 2p |
| --- | --- | --- | --- | --- | --- | --- |
| 13.5 | 10⁻³ | 20 | 40 | 34 | 34 | −6 |
| 13.5 | 10⁻⁶ | 25 | 50 | 50 | 50 | 0 |
| 13.5 | 10⁻⁸ | 27 | 54 | none | 51 | −3 |
| 27.1 | 10⁻³ | 34 | 68 | 60 | 61 | −7 |
| 27.1 | 10⁻⁶ | 41 | 82 | 74 | 74 | −8 |
| 27.1 | 10⁻⁸ | 44 | 88 | 86 | 86 | −2 |
| 54.1 | 10⁻³ | 61 | 122 | 114 | 114 | −8 |
| 54.1 | 10⁻⁶ | 71 | 142 | 127 | 127 | −15 |
| 54.1 | 10⁻⁸ | 76 | 152 | 136 | 136 | −16 |
| 108.3 | 10⁻³ | 115 | 230 | 218 | 224 | −6 |
| 108.3 | 10⁻⁶ | 128 | 256 | 234 | 237 | −19 |
| 108.3 | 10⁻⁸ | 134 | 268 | 243 | 246 | −22 |
| 216.5 | 10⁻³ | 224 | 448 | 431 | 436 | −12 |
| 216.5 | 10⁻⁶ | 240 | 480 | 450 | 455 | −25 |
| 216.5 | 10⁻⁸ | 248 | 496 | 462 | 468 | −28 |

Over the sweep, for real k the longest is 1–10 degrees below 2p at k a 22.6–32 and 6–28 below
at k a ≥ 45. At k a 11.3–16 it is 2p at 10⁻⁶ and 10⁻⁷ (one or two classes need 2p), 0–2 below at
10⁻⁵; at 10⁻⁸ the nearest classes are admissible at no length (on RCR by the bound) and the
others reach 2p − 3 ("none" in the table). At k a ≤ 8 with the plateau's degrees,
classes whose rounding already fails at 2p (on RCR by the bound) are admissible at no scanned
length or only at shorter ones (§4). For η = 0.1 the picture is the same, the longest lengths 5
degrees shorter to 14 longer than real k's at the same (k a, ε, p), and up to 19 classes per row
admissible at no length (at k a ≥ 11 all of them on RCR by the bound: §5's rule, safe on every
class the bound takes, finds a length for each of those).

**The grid.** A length L needs a grid exact for degree 2p + L: ⌈(2p + L + 1)/2⌉ Gauss–Legendre
rings and 2p + L + 1 azimuths padded to a multiple of 4, against 2p + 1 and 4p + 4 at L = 2p. The
level's conversions and pairs scale with the samples. With the longest admissible length over the
classes the bound takes at 2p (`rules`, column "measured"):

| k a | 10⁻³ | 10⁻⁶ | 10⁻⁸ |
| --- | --- | --- | --- |
| 13.5 | 0.839 | 1.000 | 0.929 |
| 27.1 | 0.888 | 0.906 | 0.870 |
| 54.1 | 0.936 | 0.892 | 0.898 |
| 108.3 | 0.970 | 0.924 | 0.916 |
| 216.5 | 0.973 | 0.945 | 0.945 |

So the oracle's grid is 0.84–1.00 of 2p's on W1's diagonal levels, 0.86–0.97 over the sweep's
real k at k a ≥ 22.6. The
diagonal levels are part of the M2L, so the evaluation would gain a fraction of that.

**Classes from RCR** (`gain`): classes the bound puts on RCR at 2p that some admissible length puts
within θ ε. There are two at k a = 27.1 at 10⁻⁸ (6 of an interior box's 316 offsets, about 1.9 %
of its V pairs; on W1's sphere the share per level differs and was not counted) and none on any
other W1 level at any ε. At k a = 13.5, 10⁻⁸ the two RCR classes need 2p and fail there.

## 4. The named lengths

`named`: the classes (of 34) on which each length is admissible.

- **⌈3(P − 1)/2⌉**: at k a ≥ 11.3, real and damped, on one class of one row (real k, k a 11.3,
  10⁻⁷) and on none elsewhere. At
  k a 2.8–8 with the plateau's degrees it is admissible on 0–27 classes, mostly the ones whose
  rounding fails at 2p: S6's "within a few times ε where 2P − 2 loses orders to rounding", below
  the rule's k a ≥ 2.8 band edge and on classes the bound sends to RCR.
- **2ka + 1.8 (log₁₀ 1/ε)^{2/3} (2ka)^{1/3}** (capped at 2p; `rules`, columns "1.8" and "1.8
  fails"): for real k it is admissible on every class the bound takes at k a ≥ 45. There it gives
  0.92–1.00 of 2p's samples (1.000 wherever it reaches 2p: at ε ≥ 10⁻⁴ to 10⁻⁵ from k a ≈ 54 up,
  since the FMM's degrees are already tight there). At k a 11.3–32 it fails on 1–4 of the
  classes the bound takes on 31 (kind, k a, ε) rows. The worst is the nearest class at
  k a 11.3–13.5, 10⁻⁵–10⁻⁷: up to 7.0× RCR's error per pair, 30–100 ε at 10⁻⁷. Overall the length
  is not admissible on 68 of the 4,356 classes the bound takes at k a ≥ 11.
- **One constant does not fit.** The γ that 2ka + γ (log₁₀ 1/ε)^{2/3} (2ka)^{1/3} needs to reach
  the longest admissible length is 0.1–1.4 at k a ≥ 90 and 0.5–2.3 at k a 11–64 (real k). It grows as
  ε tightens and k a falls, because the FMM's p already carries the ε-dependence.

## 5. A truncation model

The candidate of the plan: τ_L = Σ_{L<l≤2p} (2l + 1)|jₗ(2ka)||hₗ(k|t|)| / |h₀(k|t|)|, Gegenbauer's
tail at the radius 2a relative to the kernel at the centres' distance (for real k it is |k||t| Σ …,
E's form with |jₗ(2ka)| in place of u). It is computed in f64 (`transfer::log10_tail`) and
compared with T's measured truncation (dd) where the double-double chain resolves it (p ≤ 100, at
least 100× its own rounding: this drops k a ≤ 4 rows where T's terms reach 10³⁰). That leaves
93,564 rows (`model`):

| kind | k a | max | median | min |
| --- | --- | --- | --- | --- |
| real | 2.8 | −1.29 | −2.11 | −3.09 |
| real | 4.0 | −1.09 | −2.10 | −3.28 |
| real | 13.5 | −1.83 | −2.65 | −3.88 |
| real | 27.1 | −2.35 | −3.24 | −4.68 |
| real | 64.0 | −2.29 | −4.15 | −6.33 |
| η = 0.1 | 27.1 | −3.49 | −4.47 | −5.90 |
| η = 0.1 | 64.0 | −5.83 | −7.84 | −9.88 |

(log₁₀ of measured over τ_L; all rows in the table script's output.) C_τ = 0.041 puts the largest
ratio at 2×. The model is safe but loose: 67× in the median, and the looseness grows with k a and
with damping. That is where the room below 2p is. The rule it gives, the shortest scanned L
with 2 C_τ τ_L ≤ 0.1 ε, is admissible on all 4,356 classes the bound takes at k a ≥ 11. It returns
2p on 126 of 134 (kind, k a, ε) rows; the other eight (real k a ≥ 90, and η = 0.1 at k a 45.3
and 54.1, 10⁻⁷–10⁻⁸) keep 0.970–0.994 of the samples.

Why the model is loose: it bounds the tail of the full kernel at the radius 2a. The FMM's
multipole and local are truncated at degree p, and their coefficients decay beyond k a. A model
of the dropped terms Σ_{n+n′>L, n,n′≤p} |jₙ(ka)||jₙ′(ka)||h_{n+n′}(k|t|)| would be tighter. At most
it could reach the oracle's grid of §3 (3–14 % fewer samples at k a ≥ 22.6). Whether it would is
not measured here.

## 6. The rounding at shorter lengths

`round`: the f64 chain against the double-double chain at the same length, over E at that length
(proposal §3.6 with max_{j≤L}, C_T = 0.15, C₀ = 200, θ = 1/2; the damping factor for η = 0.1):

- L = 2p: 5,304 rows, at most 1.99 E (k a 2.83, p 7, the nearest class: S6's own point, where
  C_T was fitted to put the largest at 2), median 0.034 E;
- L < 2p: 94,962 rows, median 0.029 E; above 2 E six rows, all at k a 2.83 (p 9, L 13–17: at most
  2.73 E), the sweep's smallest k a, at the rule's band edge (k a ≥ 2.8), where S6 fitted C_T; at
  k a 4–8 at most 1.11 E (real k; 0.60 E for η = 0.1), at k a ≥ 11 at most 0.38 E.

So a shorter length lowers E with the rounding: the bound's form carries over, and it is not what
limits a shorter T.

## 7. Limits of this study

- **W3 and W1 at kD = 100 and 200.** Their diagonal levels' k a (86.6, 43.3, 21.65, 10.8, 5.4) are
  not on the sweep; the nearest sweep points (90.5, 45.3, 22.6, 11.3, 5.7) give the same picture
  (§3, `lmin` and `rules` over the sweep), so the conclusion covers them, but their per-level
  grids are not tabulated.
- **The table build** (the plan's gain (iii)): T's samples cost O(L) each, so a shorter length
  saves about (2p − L)/2p of a level's transfer-function work besides the smaller grid, a few per
  cent of a build part CH4.2c's interpolated build already made small; not measured.
- **Truncation above p = 100.** The double-double chain runs at p ≤ 100 only. The admissibility
  test reads the total against the direct sum and needs no separate truncation; the model of §5
  is fitted on the resolved rows (p ≤ 100) and applied above as an extrapolation, where its rule
  is still safe on every class (the total decides). Above p = 100 T's truncation is not reported
  separately.
- **The model's damping.** `log10_tail` divides by |h₀(k|t|)| instead of applying E's factor
  e^{Im(k)(|t| + 2a)}; for η = 0.1 it over-predicts by 10³–10⁶ (§5's table), so it is conservative.

## 8. Questions for sign-off

1. **The acceptance as read.** CH4.3d1 answers CH4.3d's pair questions: truncation and rounding
   per class against RCR, the excess-bandwidth length in the same sweep, the shortest lengths,
   a truncation model. "The FMM where one earns it" is CH4.3d2's. *Recommended*: accept.
2. **CH4.3d2.** *Recommended*: not run. Keep 2P − 2, and close CH4.3d2 unrun. The oracle's gain
   is 0–16 % of a diagonal level's samples, plus two classes on one level at 10⁻⁸. The one safe
   a-priori rule gains 0–3 %. Song and Chew's length is unsafe at k a 11–32. A tighter model is
   a study whose ceiling is the oracle's few per cent. Alternatives: (b) run CH4.3d2 with
   Song and Chew's length restricted to real k and k a ≥ 45, where it is safe (0–8 % of those
   levels' samples). That costs a new grid size per level in the tables, the conversions and the
   time rule's costs, then Kathleen runs. (c) First a study of the truncated-expansion model
   (§5), then decide.
3. **The proposal.** §3.6 says of ⌈3(P − 1)/2⌉ that it is "a shorter T trading truncation for
   stability; CH4.3's to analyse". *Recommended*: amend it to this report's finding:
   ⌈3(P − 1)/2⌉ is within RCR's error on no class at the FMM's degrees from k a 11.3; the shortest
   lengths within RCR's error save 0–16 % of a level's samples; no a-priori rule safely reaches
   them; 2P − 2 is kept (amended at CH4.3d1, with a line in Appendix A). Review M8 closes with it.

## 9. Reproduce

From `spikes/` (the spikes workspace; `spikes/README.md`):

```sh
cargo test -p tetraleaf-spike-diagonal --release --test lengths -- --nocapture
cargo run -p tetraleaf-spike-diagonal --release --example lengths -- --list
cargo run -p tetraleaf-spike-diagonal --release --example lengths -- --threads 70 --out <dir>/lengths.csv
```

Then, from the root, `python3 tools/helmholtz/ch4-3d1-tables.py <dir>/lengths.csv` (sections
`lmin`, `named`, `model`, `round`, `gain`, `rules`, `trunc`); their output is committed as
`spikes/diagonal/results/lengths-summary.md`. On locust the sweep ran from the tree
`tools/gh200/sync.sh` copied, writing outside it (`/data/ucahtbe/runs/ch4-3d1/`); its log is
archived with the rows.
