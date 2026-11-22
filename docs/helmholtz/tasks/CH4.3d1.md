# CH4.3d1 The transfer function's length, per pair

Started 2026-11-22 on `tbetcke/h-ch4-3d1`. Roadmap: CH4.3d1 "The transfer function's length, per
pair" (the first of CH4.3d1–CH4.3d2; CH4.3d split on 2026-11-22 on this branch by the size rule
of `/next-task` and the H6 gate's rule of pair and FMM evidence: a spike sweep on locust that may
recommend, then a change to the tables and the time rule with runs on Kathleen that confirms).
Proposal: §10 CH4.3 (the per-(level, class) rule, calibrated); §3.6 (the diagonal form, its
rounding, the bound E(l, c); "⌈3(P − 1)/2⌉ … a shorter T trading truncation for stability;
CH4.3's to analyse"); §5.4 (the grid (2P − 1)(4P − 3), exact for degree 4P − 4); §3.5 (the
calibrated degrees, CH4.3c2's branch); §2 (constraints 1, 2, 9). Review 2026-11-12, M8
(`docs/reviews/2026-11-12-H6.md` §1). S6's report §6 (`spikes/diagonal/SPIKE_REPORT.md`).

## Goal
Measure, per pair on S6's chain, what shortening the transfer function below 2P − 2 costs in
truncation and gains in rounding and grid size at the FMM's degrees, model the truncation, and
say whether a length rule is worth taking into the FMM (CH4.3d2).

## Acceptance (the roadmap's CH4.3d1, from CH4.3d's *Done when*)
CH4.3d's criterion: "(S6 §6) whether a shorter transfer function (⌈3(P − 1)/2⌉, whose rounding
stays far below 2P − 2's at k a 6–45, S6's report §6) earns a truncation analysis; (review
2026-11-12, M8) the excess-bandwidth length L ≈ 2ka + c (2ka)^{1/3} of the high-frequency FMM
literature, between S6's last two, in the same sweep: truncation and rounding per class against
RCR, and the FMM where one earns it; any adoption for sign-off."

CH4.3d1's part, as split: the pair sweep (truncation and rounding per class against RCR, the
total against the direct sum), the shortest length within RCR's error per (k a, ε, class), a
truncation model against the measurement, the gain of a length rule on W1's and W3's levels,
and whether CH4.3d2 runs, for sign-off. "The FMM where one earns it" is CH4.3d2's: by the H6
gate's rule a numerical rule changes only on pair and FMM evidence both, and this task may only
recommend.

M8 as the review wrote it: "S6 tested L = P − 1, ⌈3(P − 1)/2⌉ and 2P − 2. The excess-bandwidth
formula of the high-frequency FMM literature, L ≈ 2ka + c (2ka)^{1/3}, lies between the last
two. As the rounding grows exponentially in L − k|t|, a few degrees fewer could move the nearest
classes from RCR to the diagonal form at 10⁻⁸." Checked against the tree on 2026-11-22: the
length is 2p everywhere (`tetraleaf-ref`'s `Transfer`, `bound`'s `log10_hankel_term(z, 2 * p)`,
`AngleTable`; `tetraleaf-tables`' builds); the grid is (2p + 1) rings × 4p azimuths
(`diagonal::Grid`); no study since S6 has run another length.

What the items mean here (P = p + 1, so 2P − 2 = 2p):
- **Truncation.** In exact arithmetic the chain at length L < 2p drops the (S|R) terms of
  degree l = n + n′ > L; it is no longer the truncated M2L. Measured as the chain in
  double-double at length L against S7's double-double RCR at degree p, which resolves it to
  about 10⁻²⁰ (S6's `diagonal_is_rcr_in_double_double`).
- **Rounding.** The f64 chain (T rounded once from double-double, as the FMM builds it) at
  length L against the double-double chain at the same L, read against E at length L
  (proposal §3.6 with max_{j ≤ L}; C_T = 0.15, θ = 1/2): does the bound's form, fitted at 2p,
  still cover the rounding at shorter L within 2E?
- **Total.** The f64 chain at L against the direct sum in double-double, beside RCR's own
  truncation (RCR against the direct sum).
- **Within RCR's error.** A length L is admissible for class c at (k a, ε), p = the FMM's rule,
  when on every placement and both norms (potential and gradient, S7's relative ℓ2 over the
  targets) the f64 chain's total error at L is at most max(1.1 × RCR's, RCR's + 0.1 ε): the
  diagonal form at L as accurate as RCR within CH4.3b's 1.1×, or within a tenth of ε more where
  RCR's own error is far below ε (the far classes). The shortest admissible scanned length is
  L_min(k a, ε, c).
- **A truncation model.** The candidate, Gegenbauer's tail at the pair's largest radius 2a,
  τ_L = C_τ |k||t| Σ_{L < l ≤ 2p} (2l + 1)|j_l(2ka)||h_l(k|t|)| (times the damping factor of
  E for Im k > 0), the relative form of E's |k||t| (2j + 1)|h_j| with |j_l(2ka)| in place of u.
  Accepted as a model if a C_τ makes it cover the measured truncation within 2× (as S6's
  acceptance of E) over the sweep, and the length rule "the shortest L with τ_L ≤ θ_τ ε" is
  stated with it. If it does not, the report says how it fails, and an empirical rule is
  fitted to L_min instead (the excess-bandwidth form 2ka + γ (ln 1/ε)^{2/3} (2ka)^{1/3},
  γ fitted, capped at 2p), and the model is left as the open question it is.
- **The gain.** On the FMM's diagonal levels of W1 (kD = 100–1000) and W3 (kD = 200, 500) at
  10⁻³–10⁻⁸, from the measured L_min and the model or fit: (i) the grid a shorter length allows
  (exact for degree 2p + L: ⌈(2p + L)/2⌉ + 1 rings, 2p + L + 1 azimuths padded to a multiple of
  4, against 2p + 1 and 4p), set by the longest length of the level's diagonal classes, as a
  ratio of samples (the conversions and pairs scale with it); (ii) the classes the rule takes
  from RCR (E at L_c ≤ θ ε where E at 2p is not) and their share of the level's V pairs;
  (iii) the table build (O(L) per sample). The level's k a, p and pair counts from CH4.3b2c's
  report and runs.
- **Earns CH4.3d2.** The report recommends CH4.3d2 if the gain is material: the samples of
  some W1 or W3 diagonal level at least 5 % fewer, or a class of a level taken from RCR, at an
  ε of D1's range; otherwise it recommends dropping CH4.3d2 and keeping 2P − 2, which is the
  proposal's and S6's. Either way for sign-off.

## Read first
- `spikes/diagonal/SPIKE_REPORT.md` §1, §3, §6; `src/run.rs` (the unit, the variants, `Band`),
  `src/transfer.rs` (snapshots at several lengths in one pass; `log10_term`), `src/grid.rs`,
  `src/rules.rs`, `examples/sweep.rs`, `examples/summary.rs`.
- `spikes/truncation/src/{chain,radial}.rs` (S7's RCR, Bessel and Hankel functions in
  double-double).
- `tetraleaf-exec/src/helmholtz/degree.rs` (`DegreeRule::CALIBRATED`, CH4.3c2's branch).
- `tetraleaf-ref/src/helmholtz/diagonal/{bound,grid,transfer}.rs` (the product's length 2p).
- `tetraleaf-validate/results/helmholtz-ch4-3b2c.md` (the diagonal levels of W1 and W3, their
  classes on RCR and their pair counts).

## Changes
All in `spikes/diagonal` (the spikes workspace), plus the report and the records:
- `src/rules.rs`: CH4.3c2's `DegreeRule::CALIBRATED` (copied as CH3.3's was, with its
  constants), `refit`; CH3.3's `calibrated` stays (S6's sweep cites it).
- `src/lengths.rs` (new): the scan of lengths for one unit (kind, k a, p): the lengths
  (`lengths(ka, p, epsilons)`: from L_lo = max(p, min(⌈3p/2⌉, ⌊2ka⌋ − ⌈(2ka)^{1/3}⌉)) to 2p in
  steps of max(1, ⌈(2p − L_lo)/24⌉), plus ⌈3p/2⌉ and the excess-bandwidth lengths
  ⌈2ka + 1.8 (log₁₀ 1/ε)^{2/3} (2ka)^{1/3}⌉ (Song and Chew's constant) of the unit's ε, within
  [p, 2p]); per class and placement, S6's references (RCR and the direct sum in double-double),
  the f64 chain at each length (T's snapshots, one pass), and where p ≤ 100 the double-double
  chain at each length (`run::diagonal::<Dd>`); the errors (truncation, rounding, deviation of
  the f64 chain from RCR, total) and per length log₁₀ of E's term and of τ's sum. `run.rs`'s
  `run_unit` and its variants are not changed (S6's rows stay reproducible).
- `src/transfer.rs`: `log10_tail`, log₁₀ of |k||t| Σ_{L < l ≤ 2p} (2l + 1)|j_l(2ka)||h_l(k|t|)|
  in f64 with rescaling (as `log10_term`).
- `examples/lengths.rs` (new): the sweep, as `sweep.rs` (resumable CSV, shards, `--list`,
  `--kinds`, `--ka`), units = the distinct p of the refit rule at 10⁻³ … 10⁻⁸ per (kind, k a)
  below the pruning cut; k a = S6's grid 2^{j/2}, j = 3..15, plus 216.5 and the four levels of
  kD = 500 (108.3, 54.1, 27.1, 13.5); kinds real and η = 0.1.
- `examples/lengths_summary.rs` (new): the CSV reduced to the report's tables (L_min per
  (kind, k a, ε, class); truncation against τ_L; rounding against E_L; the gain per level of
  W1 and W3).
- `tests/lengths.rs` (new): the tests below. `lib.rs` and `SPIKE_REPORT.md` (a short section
  pointing at CH4.3d1's report), `spikes/README.md`'s row.
- Report `tetraleaf-validate/results/helmholtz-ch4-3d1.md`; the summary as
  `spikes/diagonal/results/lengths-summary.md`; raw rows archived in
  `~/projects/fmm-archive/locust/runs/ch4-3d1/`, not committed.
- No product crate, no public API of the root workspace, no convention, no dependency (the
  spike already depends on `truncation` alone; `spikes/Cargo.lock` unchanged).

## Tests that prove it
`cargo test -p tetraleaf-spike-diagonal --release` from `spikes/`:
1. **The scan is S6's chain at 2p**: at L = 2p the new path's f64 local equals `run_unit`'s
   `f64:2P-2` variant bit for bit (same samples, the same operations in the same order), and
   its errors the same, on a small unit (real k, k a 4, p 10, the 34 classes) and a damped one.
2. **Truncation resolved**: the double-double chain at L = 2p against double-double RCR within
   10⁻²⁰ relative (S6's identity, here through the new path); at L < 2p the truncation
   decreases (not strictly per step, but from L_lo to 2p by orders of magnitude) on the nearest
   and the farthest class.
3. **The tail term**: `log10_tail` against the double-double sum of S7's `radial` Bessel and
   Hankel functions within 10⁻¹⁰ in log₁₀ where both are in range (k a 3–64, |t|² = 4, 12, 27),
   and finite where hₗ leaves f64's range; 0 terms (L = 2p) give −∞.
4. **The lengths**: ascending, distinct, within [p, 2p], containing L_lo, ⌈3p/2⌉ where in range,
   the excess-bandwidth lengths in range, and 2p.
5. **The refit rule**: `refit` equals `DegreeRule::CALIBRATED`'s degrees at x = 1, 2.7, 13.5,
   21.65, 27.1, 54.1, 108.3, 216.5 and ε = 10⁻³ … 10⁻⁸, the expected values computed by
   `tetraleaf-exec`'s rule on this branch (a table pasted into the test, its source named).
6. **Determinism**: a unit's rows identical at 1 and 3 threads.
The sweep itself is the measurement, not a test; no timing is asserted.

## Checks
`/verify`: the branch changes `spikes/` and docs only, so the root workspace's checks run as the
skill picks them (fmt; no root crate changes), and from `spikes/`: `cargo fmt --all --check`,
`cargo clippy -p tetraleaf-spike-diagonal --all-targets -- -D warnings`, `cargo test -p
tetraleaf-spike-diagonal --release` (with `--include-ignored` once), `cargo check
--all-targets` (every spike still builds).

## Machine runs
- locust (accuracy only, no timing; a shared node, so `uptime` checked and the threads set
  below its idle cores): `tools/gh200/sync.sh`, then from `spikes/` `cargo run -p
  tetraleaf-spike-diagonal --release --example lengths -- --threads <n> --out <file>` in shards,
  detached, its log and CSV fetched before any re-sync. Estimate from S6's sweep log (k a ≤ 256,
  real k: 22–131 s per unit at p 279–318 on 72 threads with eight variants): about 2–4 hours of
  wall time for both kinds. It decides L_min, the model and the gain.
- The M6: the tests and a smoke sweep (k a ≤ 16), accuracy only (no timing, so no quiet run to
  ask for).

## Constraints and risks
- Constraint 1 (Laplace): not touched (a spike and docs).
- Constraint 9 (degrees and lengths from k, the domain, ε and the class): a length rule depends
  on k a, ε and the class only.
- Constraint 2 (determinism): the spike's units are deterministic per thread count (test 6).
- Risk: the double-double chain at every length costs O(p³) double-double work per length and
  placement; limited to p ≤ 100, where the rounding of the f64 chain can hide the truncation.
  Above p = 100 the truncation is read from the f64 chain's deviation from RCR where E at L is
  at least 10× below it, and marked unresolved otherwise.
- Risk: the scan's step at large p (up to 6 degrees at p ≈ 250) overstates L_min by up to a
  step; the excess-bandwidth lengths are in the scan exactly, and the report states the step.
- Risk: the gain is small (a priori: at k a ≥ 54 the nearest class's E is at its floor already
  at 2p, CH4.3c1, and the grid set by the longest length shrinks a few per cent; at k a ≤ 27
  Gegenbauer's tail converges like (2a/|t|)^l = 0.87^l on the nearest class, so a shorter T
  costs truncation fast). Then the outcome is "keep 2P − 2, drop CH4.3d2", which is an answer.
- Risk: the spike's degree rule copied from `tetraleaf-exec` drifts; test 5 pins it.

## Pending runs
None yet.

## Open questions
None yet.
