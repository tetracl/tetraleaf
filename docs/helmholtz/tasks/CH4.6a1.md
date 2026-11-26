# CH4.6a1 The trait and rustfft on the product's grid

Started 2026-11-26 on `tbetcke/h-ch4-6a1`. Roadmap: CH4.6a "The optional rustfft back end of the φ
step" (added at CH0.9, D7; item 13 of the priority list). Proposal: §10 CH4.6a; §6.7 ("Grid
sizes", "Plans per box", "Determinism and its preconditions", "Host back ends"), §5.4, §12 D7;
§2 constraints 1, 2, 6, 9. Spike: `spikes/fft-plans/SPIKE_REPORT.md` (S9), §3, §5–§7.

## Size: CH4.6a split

The §10 row is four pieces of work: a back end with its tests, rank check and CI job; times
per level on three machines; a second grid (the smooth N_φ, which changes the signatures, the
transfer functions and their tables, under a key of their own); and a rule choosing the grid per
level from measured times. That is more than one reviewable PR, more than two crates
(tetraleaf-tables, tetraleaf-exec, tetraleaf-ref for the second grid, CI) and an operator plus
measurements on several machines. Split on 2026-11-26, measuring before building the second grid (as CH4.1a
measured before CH4.2); the user: "Split a1/a2/a3 (Recommended)":

- **CH4.6a1 The trait and rustfft on the product's grid.** The φ step behind a trait internal to
  `tetraleaf-exec`, the product its default, rustfft behind an optional feature `rustfft`, on the
  product's grid (`azimuth_count`). *Done when:* against the product to round-off; per-box plans
  with a box's bits the same across batch positions, batch sizes and thread counts; at build
  rustfft's version, each level's N_φ and the dispatched CPU features all-gathered, and with one
  of two ranks made to differ (a test hook) both ranks run the product; the feature off by
  default, `cargo tree` of the default members without rustfft, one CI job with the feature.
  *Needs:* CH4.1b, CH0.9, D7.
- **CH4.6a2 The time per level against the product, both grids priced.** Each diagonal level's
  conversions and level call, the product against rustfft on the product's grid, on locust,
  Kathleen and the M6 (the M6 with the user's go-ahead), on W1 at kD = 250–1000 and W3 at
  kD = 500; the smooth grid priced from S9's per-box FFT times and the level's measured pair cost
  at (N_φ^smooth / N_φ) more samples. Ends with a question: build the smooth grid (CH4.6a3) or
  keep the product's grid for the FFT and close CH4.6a3 unrun. *Needs:* CH4.6a1. *Machines:*
  locust, Kathleen, the M6.
- **CH4.6a3 The smooth grid, chosen per level by its time.** If CH4.6a2's sign-off keeps it:
  `Grid` at a chosen N_φ (4P − 3 up to a 7-smooth multiple of 4), the conversions and transfer
  functions on it under their own table key, and a level's grid chosen by a rule fitted on
  CH4.6a2's times per machine class (constraint 9: a function of k, the geometry, ε and the
  calibration tables, never a time measured at run time). *Needs:* CH4.6a2 and its sign-off.

The rest of this plan is **CH4.6a1**.

## Goal (CH4.6a1)
A second, optional φ step for S → F and F → R on the diagonal form's levels: rustfft per box, on
the product's grid, selected by a builder option, falling back to the product where the ranks'
FFT builds differ; the product unchanged bit for bit and the default.

## Acceptance (proposal §10 CH4.6a, verbatim)
"against the product to round-off; per-box plans with bits identical across batch positions,
batch sizes and thread counts; at build rustfft's version, N_φ and the dispatched CPU features
all-gathered, and with one of two ranks made to differ (a test hook) both ranks run the product;
the feature off by default, one CI job with the feature; a level's grid (the product's or the
smooth one) chosen by its measured time; the time per level against the product on the M6,
locust and Kathleen"

CH4.6a1 takes all but the last two clauses (CH4.6a2, CH4.6a3).

## Read first
- `tetraleaf-tables/src/helmholtz/diagonal/mod.rs`: `Conversions` (`s2f_block`, `f2r_block`:
  the Legendre step per block of rings into per-order rows, then the φ step per residue μ mod 4
  over the base azimuths), `ConversionScratch`; `level.rs` (`Signatures`, `s2f_strides`,
  `f2r_from`).
- `tetraleaf-ref/src/helmholtz/diagonal/grid.rs`: φ_J = 2πJ/N_φ, `azimuth_count`.
- `tetraleaf-exec/src/helmholtz/operator_diagonal.rs` (the level's conversions, the threaded level
  call) and `fmm.rs` (the build's all-reduces, e.g. the agreed degrees).
- `spikes/fft-plans/src/lib.rs` (`RustPhi`: one plan per (N_φ, direction), executed ring by ring;
  the scatter of 2P − 1 modes; `simd_features`).

## Changes
- **Dependency** (signed off 2026-11-26, the user: "Exact =6.4.1 (Recommended)"): `rustfft` in `[workspace.dependencies]`, an optional dependency of
  `tetraleaf-exec` behind a feature `rustfft = ["dep:rustfft"]`, in no default set and enabled by
  no workspace crate. The pin `=6.4.1` (S9's version), so the version in the rank check is
  exact; the root `Cargo.lock` gains rustfft, `strength_reduce`, `transpose` and `primal-check`
  (num-complex is already there), all MIT OR Apache-2.0. `cargo deny` unchanged in config.
- **tetraleaf-tables** (`helmholtz::diagonal`, the Helmholtz diagonal form only; not on Laplace's
  paths): `Conversions` gains the two halves of a box's conversion around the φ step, public:
  - S → F: the Legendre step of one box into per-ring Fourier coefficients (ring i, order
    −p..p, complex: the four real planes [Re +μ, Im +μ, Re −μ, Im −μ] combined as
    c_{±μ} = (A_μ ∓ i B_μ)/2);
  - F → R: from per-ring Fourier coefficients (the analysis' 2P − 1 modes per ring, times the
    weights and 2π/N_φ as the product applies them) through the transposed Legendre step, added
    to the box's local.
  The existing `s2f`, `f2r`, `s2f_strides`, `f2r_from` stay as they are, bit for bit (the
  product's path does not go through the halves).
- **tetraleaf-exec** (Helmholtz only): a trait internal to the crate for a level's φ step, the
  product (the existing fused calls) and `RustFft` (feature `rustfft`): one plan per (level,
  N_φ, direction), built once per operator build, executed ring by ring per box with a
  per-thread scratch of two boxes (input zero outside the mode slots; analysis output), per S9.
  A builder option `HelmholtzPhiStep::{Product, RustFft}` (default `Product`; `RustFft` exists only
  with the feature). At build, per level with the diagonal form: a u64 fingerprint of
  (rustfft's version, N_φ, the CPU features rustfft dispatches on: AVX and FMA, SSE4.1 on
  x86_64, NEON on aarch64, scalar otherwise) all-reduced as min and max; if they differ, every
  rank runs that level's product. A hidden test hook perturbs one rank's fingerprint.
  Generic over `T: RealScalar` where rustfft's `FftNum` holds (f32 and f64).
- **CI**: one job, `run-tests-rustfft`, `cargo clippy -p tetraleaf-exec --all-targets --features
  rustfft -- -D warnings` and the Helmholtz tests with the feature, the rank test at 2 ranks
  through `run-mpi-tests.sh`; plus `cargo tree -e normal` of the default members shows no
  rustfft (in `run-licence-check`).
- Documentation: the crate docs' features section, the builder option, `tetraleaf-exec/CLAUDE.md`
  and `tetraleaf-tables/CLAUDE.md`; root `CLAUDE.md`'s CI section.
- Conventions: none (the samples and their order are §3.15's; only the arithmetic changes).
- Report scripts: none in CH4.6a1.

## Tests that prove it
- `tetraleaf-tables`: the halves composed with the product's φ step written as a plain sum equal
  `tetraleaf_ref::helmholtz::diagonal::{s2f, f2r}` within 256 u of the largest value (as
  `tests/helmholtz/diagonal.rs`); p ≤ 24 in debug, larger p ignored in release.
- `tetraleaf-exec`, feature on: S → F and F → R with `RustFft` against the product within
  256 u of the largest value, p = 6, 12, 24 (and 41, 72, 131 ignored, release), f64 and f32
  (f32 within its own u).
- Bits per box: a box's signature and local the same bits at every position in batches of 1, 2,
  7, 33 boxes, and on pools of 1, 2, 4 threads (and all threads, ignored).
- The FMM: the Helmholtz FMM with `RustFft` against the product's FMM on the accuracy tests'
  cases with the diagonal form on (the default `Auto` on a high-frequency tree), errors against
  the direct sum within 1.01× the product's, and bit for bit across 1 and 4 threads.
- Ranks (MPI, in `helmholtz_fmm`, feature on): at 2 ranks with `RustFft`, the output bit for bit
  the one-rank output; with the hook on one rank, both ranks report the product on every
  diagonal level and the output is the product's bit for bit.
- Feature off: the default build's outputs unchanged (the existing fixtures and bit-identity
  tests); `cargo tree -e normal` lists no rustfft.

## Checks
`/verify` (default, strict workspace, MPI list at 2 and 4 ranks, `cargo deny`), and with the
feature: `cargo clippy -p tetraleaf-exec --all-targets --features rustfft -- -D warnings`,
`RUST_MIN_STACK=8388608 cargo test -p tetraleaf-exec --features rustfft`, the ignored large-p
tests in release, `run-mpi-tests.sh` for `helmholtz_fmm` with the feature at 2 and 4 ranks.

## Machine runs
None in CH4.6a1 (correctness only; the feature's bits per box were measured on all three
machines by S9). Times are CH4.6a2's.

## Constraints and risks
- Constraint 1: no Laplace path is touched (tetraleaf-tables' `helmholtz::diagonal` and
  tetraleaf-exec's Helmholtz operator only); the product's fused path is unchanged bit for bit.
- Constraint 2: per-box plans, ring by ring; the fallback to the product is a collective
  decision at build on data-independent inputs (constraint 9).
- Risk: F → R's analysis as an FFT gives all N_φ modes; only 2P − 1 are kept (the grid is
  oversampled, so no aliasing into them beyond the product's own).

## Pending runs
None yet.

## Open questions
Asked before implementing, answered 2026-11-26: the dependency, "Exact =6.4.1 (Recommended)";
the split, "Split a1/a2/a3 (Recommended)".
