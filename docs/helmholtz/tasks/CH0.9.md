# CH0.9 S9: FFT plans, bits and host timing

Started 2026-11-26 on `tbetcke/h-ch0-9`. Roadmap: "CH0.9 S9: FFT plans, bits and host timing"
(H0; item 13 of the priority list, the FFT track's first task). Proposal: §10 CH0.9; §6.7
(Fourier transforms: grid sizes, plans per box, determinism, the host back ends), §12 D7, §8
R15.

## Goal
Find out, on the machines of record, whether per-box FFT plans of the φ step (FFTW through
RLST 0.12, and rustfft) give a box the same bits wherever and however it is transformed, and
what they cost per box against the product, so that D7 can be proposed per machine and §6.7's
plan layout, scratch and batch (2P − 1 or 2P rings) confirmed or revised.

## Acceptance (proposal §10, verbatim)
"a spike crate outside the workspace (a `[workspace]` table of its own, so nothing enters the
root lock file) on RLST 0.12's `rlst::fftw` with `fftw_system`, the code path of CH4.6, and
rustfft. Per-box plans of 2P − 1 rings and of 2P (one zero ring), interleaved scratch against
split plans, executed on boxes at different positions in a batch, in batches of different
sizes, from different threads and in fresh processes, with buffers in class 0: bits identical;
plan fingerprints compared across two ranks. For the record, whether plans over varying batch
or stride change bits. Time per box for the product, FFTW (`Rigor::Estimate`) and rustfft, for
P from 34 to 262 on smooth and power-of-two grids, on the M3 Max, locust and Kathleen, with
`rlst::fftw::version()` recorded per machine. Outputs: the plan layout and scratch storage of
Section 6.7 confirmed or revised, the batch (2P − 1 or 2P), D7 per machine"

The roadmap's *Done when* is the same, with "D7 proposed per machine" and the sign-off D7. The
M6 replaced the M3 Max as a machine of record on 2026-10-31 (D9's amendment); it has the same
Homebrew FFTW 3.3.11 without SIMD (`tools/m6/machine.md`).

## Read first
- Proposal §6.7 (all of it), §8 R15, §12 D7; `tetraleaf-validate/results/helmholtz-chr-4.md`
  (each machine's FFTW: M6/M3 Max scalar, locust NEON, Kathleen AVX-512);
  `tetraleaf-validate/results/helmholtz-ch4-1a.md` §1–§3 (the product per box; the φ step's
  share of the conversions' multiply-adds, about 80 %; an FFT bounded at 1–2 % of a run).
- `tetraleaf-tables/src/helmholtz/diagonal/mod.rs`: `Conversions` (`s2f`, `f2r`,
  `multiply_adds`), the product as shipped; `tetraleaf-ref/src/helmholtz/diagonal/grid.rs`:
  `azimuth_count` (N_φ = 4P − 3 rounded up to a multiple of 4), the signature order i N_φ + j.
- RLST 0.12.0's `src/fftw/` (`Layout::new`, `along_axis`, `split_along_axis`,
  `Planner::<f64>`, `Rigor::Estimate`, `Plan::execute`, `Plan::fingerprint`,
  `AlignmentClass`, `FFTW_ALIGNED`, `version()`), its `docs/fftw/` and the
  `fftw_determinism` test.

## Changes
- A new spike `spikes/fft-plans/` (package `tetraleaf-spike-fft-plans`), with a `[workspace]`
  table of its own, listed in `spikes/Cargo.toml`'s `exclude` so it is in neither the root's
  nor the spikes' workspace, and FFTW and rustfft enter neither `Cargo.lock` nor
  `spikes/Cargo.lock`. Its own lock file `spikes/fft-plans/Cargo.lock` is committed (the one
  crate-level lock file, the proposal's arrangement; root `CLAUDE.md`'s "Layout" gains a
  sentence). Target directory `target/spikes-fft-plans` (its own `.cargo/config.toml`, so it does
  not rebuild the spikes workspace's shared crates).
- Dependencies of the spike only: `rlst = { version = "0.12.0", features = ["fftw_system"] }`
  (the root's requirement plus the feature), `rustfft = "6.4"` (new; MIT OR Apache-2.0),
  `mpi = "0.8.2"` (the root's) for the fingerprint check across ranks, `rayon = "1"`,
  `tetraleaf-ref`, `tetraleaf-tables` and `tetraleaf-math` by path (N_φ, the product
  conversions). No product crate, manifest or lock file changes.
- `src/lib.rs`: the φ step's three back ends on one layout. A box's φ buffer is the
  interleaved complex [N_φ × R] array, ring after ring (R = 2P − 1, or 2P with one zero ring),
  φ contiguous, so its sample order is the signature's i N_φ + j (CONVENTIONS §3.15);
  boxes packed back to back in one 64-byte aligned allocation (16 N_φ R bytes per box, a
  multiple of 64, so every box in alignment class 0, checked).
  - FFTW: one plan per (N_φ, R, direction), R transforms of length N_φ, stride 1, distance
    N_φ, `Rigor::Estimate`, made once and executed box by box on the box's slice; the split
    variant (`Layout::split_along_axis` over two real planes) for the comparison.
  - rustfft: one `Arc<dyn Fft<f64>>` per (N_φ, direction) from `FftPlanner`, executed per
    ring with `process_with_scratch` and a per-thread scratch.
  - The φ step around the transform: S → F scatters each ring's 2P − 1 modes μ = −p..p into
    slots μ mod N_φ (zero elsewhere) and transforms; F → R transforms backward and gathers
    the 2P − 1 slots. A naive DFT (the φ step as a product, in f64 and double-double from
    `Grid`'s azimuths) is the trusted path for the accuracy check.
- `examples/record.rs`: the "for the record" variants (one plan over the batch, other ring
  strides, class 16, 2P rings, split, rustfft) against the per-box plan, and FFTW's plan
  descriptions. The fresh-process comparison is a test (`bits_in_fresh_processes`).
- `examples/times.rs`: time per box, one thread, for each p and back end; Markdown table.
- `examples/fingerprints.rs` (MPI): every rank plans every (p, grid, R, direction) and
  all-gathers the fingerprints; prints whether the ranks agree.
- `run.sh <machine>`: every run of one machine into `results/<machine>/` (`tests.txt`,
  `record.md`, `times-<variant>-<grid>.md`, `fingerprints.md`, `load.txt`), one times setting per
  launch.
- `SPIKE_REPORT.md`; `spikes/README.md` (the spike's row and its build line, the exception to
  "one spikes workspace").
- `tools/helmholtz/ch0-9-tables.py`: the report's time tables from one machine's times files,
  reading them with `common.tables` and `common.find_table`.
- Conventions: none.

## Tests that prove it
In `spikes/fft-plans/tests/` (run by hand, not in CI, like every spike):
- **Accuracy** against the naive DFT (as built: in f64, bound 10⁻¹³; the double-double
  reference and 10⁻¹⁴ planned here were not built), forward and backward, on random modes, every
  back end, p ∈ {33, 82, 261}, the product, smooth and power-of-two N_φ.
- **Bits per box, FFTW and rustfft** (`assert_eq!` on the `u64` bits): one box's input at
  every position of packed batches of B ∈ {1, 2, 7, 33} boxes; transformed by 1, 2, 4
  and all rayon threads (the box on different workers); R = 2P − 1 and 2P, interleaved and
  split, each against itself; every buffer's `AlignmentClass` 0.
- **Fresh processes**: the test runs its own binary twice as child processes (the ignored test
  `child_hashes`) and compares the printed hashes with each other and with its own.
- **For the record** (printed, not asserted): a plan over the whole batch (howmany = B · R)
  and over strides other than N_φ against the per-box plan: equal bits or not; 2P against
  2P − 1 rings on the live rings; split against interleaved.
- **Fingerprints**: `fingerprints` at 2 ranks agrees on every plan (asserted on every rank,
  all-reduced so it fails on every rank or none).

## Checks
- `cargo fmt --all -- --check` and `cargo clippy --all-targets -- -D warnings` in
  `spikes/fft-plans/` (with the spike's FFTW, on the M6 and on locust); `cargo test --release`
  there; `cargo check --all-targets` in `spikes/` (the spikes workspace still builds without
  it); in the root, `cargo tree -e normal` lists no fftw-sys and `cargo deny --locked check
  licenses` passes (nothing of the root changes; shown anyway).
- `/verify`: the branch changes no product crate; its section for documents and spikes.

## Machine runs
Each run builds in release, BLAS thread variables at 1, one thread for the times, and records
the machine, CPU, `rlst::fftw::version()` and rustfft's version.
- **locust** (NEON FFTW 3.3.11): `tools/gh200/sync.sh`, then
  `tools/gh200/remote.sh 'cd spikes/fft-plans && cargo test --release && cargo run --release --example bits && cargo run --release --example times && mpirun -n 2 target/.../fingerprints'`.
  Decides D7 on locust.
- **Kathleen** (AVX-512 FFTW 3.3.10): a job script `spikes/fft-plans/kathleen.sbatch` (the
  spike's, not `tools/kathleen/jobs/`): one node in `singlenode` for the tests, `bits` and
  `times` (one setting per launch for each back end's times, the 1.09–1.20× second-launch
  drift); `fingerprints` across 2 nodes, 1 rank each, in `test`. Decides D7 on Kathleen and
  the cross-node fingerprint check.
- **M6** (scalar Homebrew FFTW 3.3.11): the same commands, `HWLOC_COMPONENTS=-darwin` and the
  loopback flags for `mpirun -n 2`; **times only with the user's go-ahead, on a quiet
  machine**. Decides D7 on the M6.
- p ∈ {33, 41, 50, 66, 82, 100, 123, 135, 144, 200, 245, 261} (P from 34 to 262, the
  high-frequency degrees of CH4.0e–CH4.5's runs), N_φ smooth (`azimuth_count`, the product's)
  and the next power of two; per box: the product's S → F and F → R (`Conversions`, whole
  conversions, and its φ share by `multiply_adds`), FFTW interleaved and split at R = 2P − 1
  and 2P, rustfft; the scatter and gather counted in each FFT time.

## Constraints and risks
- Constraint 1 (Laplace): nothing of Laplace's is touched (a spike outside both workspaces).
- Constraint 2 and 9: the spike is what checks §6.7's route to them (per-box plans, fingerprints).
- The FFTW licence: GPL; the spike is outside both workspaces and never in CI or a default
  build, so the licence check of `deny.toml` is unaffected; `cargo tree` of the root shows it.
- The product's φ step is fused into `Conversions` and not timed apart; its time is taken as
  the conversion's time times the φ step's share of its multiply-adds (CH4.1a's method), and
  the report says so. The comparison D7 needs ("faster than the product") is then the FFT's φ
  step (with scatter and gather) against that share, and the whole conversion with the FFT in
  place of the φ step estimated as (1 − share) · product + FFT.
- rustfft's planner chooses algorithms by length and CPU features (§6.7): per box it is the
  same plan, so its bits per box should hold; across machines they differ, as FFTW's do.

## Pending runs
None. Done: locust (2026-11-26, `run.sh locust`), Kathleen jobs 672793 (one node, c11b-038) and
672796 (fingerprints across c11b-013 and c11b-031; 672794 before it could not show its nodes), the
M6 (tests, record, fingerprints; no times, the user). Raw outputs in `spikes/fft-plans/results/`.

## Open questions
1. The new dependencies (rustfft; rlst's `fftw_system`; mpi) in a spike with its own workspace
   and committed lock file `spikes/fft-plans/Cargo.lock`: the user, 2026-11-26, "Go ahead
   (Recommended)".
2. The M6 times: the user, 2026-11-26, "Skip the M6". Its tests, record and fingerprints ran.
3. For sign-off on the PR: D7, the batch and layout, the grid (`SPIKE_REPORT.md` §6). The user on
   #264, 2026-11-26: "Q1-3 as recommended". Recorded in the roadmap (D7, item 13, CH4.6a and
   CH4.6b) and the proposal (§6.7, §9, §10, §12, Appendix A).
