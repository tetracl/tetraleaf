# S9: FFT plans, bits and host timing (Helmholtz CH0.9)

Spike S9 of `docs/design/helmholtz-fmm-proposal.md` (§10 CH0.9, §6.7, §12 D7, §8 R15), run on
2026-11-26 on branch `tbetcke/h-ch0-9`. Plan: `docs/helmholtz/tasks/CH0.9.md`.

## Summary

- **Bits per box hold for both FFT libraries.** For FFTW (through RLST 0.12's `rlst::fftw`,
  per-box plans of `Rigor::Estimate`, interleaved or split, 2P − 1 or 2P rings) and for rustfft
  6.4, a box's output had the same bits at every position of a packed batch, in batches of 1, 2, 7
  and 33 boxes, serially and on 1, 2, 4 and all threads (72 on locust, 40 on Kathleen), and in two
  fresh processes. This held on the M6, locust and Kathleen.
- **For the record, nothing tried changed FFTW's bits on these builds.** These were: one plan over
  the whole batch (howmany = boxes · rings), a ring distance of N_φ + 1 or N_φ + 4, alignment class
  16, and 2P rings on the live rings. FFTW loops a one-ring plan over the rings (`dft-vrank>=1`) on
  NEON (locust) and AVX-512 (Kathleen) alike. Split storage (on locust and Kathleen; not on the
  scalar M6) and rustfft give other bits, within 8·10⁻¹⁶. The per-box plan of §6.7 stays: it is free
  and holds on any build.
- **The fingerprints agree across ranks.** All 288 plans per rank had equal fingerprints on 2 ranks:
  on the M6 (one node), on locust (one node), and on Kathleen across two nodes (c11b-013, c11b-031,
  job 672796).
- **2P rings are no faster than 2P − 1** (0.99–1.01× on both machines): the batch stays 2P − 1.
  Split plans take 0.92–1.43× (locust) and 0.96–2.30× (Kathleen) the interleaved ones' time, which confirms
  §6.7's interleaved scratch.
- **The FFT's φ step costs 0.04–0.25× the product's on both machines** (FFTW on the smooth grid; the
  product's φ step taken as its share, 0.80, of the measured conversions). With the faster FFT in
  place of the φ step, a whole conversion would take 0.28–0.40× the product's on locust and
  0.23–0.29× on Kathleen. CH4.1a bounds the effect on a run at about 1–2 % (2 % at kD = 1000),
  since the conversions are a small part of a run.
- **rustfft is about as fast as FFTW:** 0.88–2.16× FFTW's time on locust and 0.55–1.73× on
  Kathleen, where it is faster at p ≤ 82. It meets D7's condition (deterministic per box, faster than
  the product) on both machines, with no system library and no GPL.
- **An FFT needs the smooth grid.** The product's grid, `azimuth_count` (4P − 3 rounded up to a
  multiple of 4), is not the smooth grid that §6.7's "Grid sizes" describes: N_φ / 4 can have a large
  prime factor (580 = 4 · 5 · 29), and FFTW on it takes up to 8.7× (rustfft up to 3.9×) its time on
  the smooth N_φ (588). Powers of two cost 0.45–1.95×.
- **D7, proposed per machine:** the product stays the default and the reference everywhere. One
  optional FFT back end on its own smooth grid, rustfft rather than FFTW (question 1 of the PR).

## 1. Set-up

The spike is `spikes/fft-plans` (package `tetraleaf-spike-fft-plans`). It has a `[workspace]` of its
own, is excluded from the spikes workspace, and has its own committed `Cargo.lock`, so FFTW's GPL
bindings and rustfft are in no other lock file (signed off by the user, 2026-11-26: "Go ahead").
Its dependencies are rlst 0.12.0 with `fftw_system`, rustfft 6.4.1, mpi 0.8.2, rayon, and
tetraleaf-math, -ref and -tables by path.

**Layout.** A box's buffer holds the rings one after another, each of N_φ interleaved complex
samples, so sample j of ring i is at i N_φ + j (CONVENTIONS §3.15's signature order). Boxes are
packed back to back in one `AlignedDynArray` at `FFTW_ALIGNED` = 64 bytes. A box is
16 N_φ R bytes, so every box is in class 0, which the tests check.

**Synthesis and analysis.** Synthesis is S → F's φ step: the 2P − 1 modes of each ring go into
the slots μ mod N_φ, followed by FFTW's or rustfft's backward transform. Analysis is F → R's: the
forward transform, then a gather of the same slots.

**FFTW.** One plan per (N_φ, rings, direction): R transforms of length N_φ, stride 1, distance
N_φ, out of place and preserving. Each thread has a scratch of two boxes: the synthesis input,
zeroed once so only the mode slots are ever written, and the analysis output. A single scratch
box for both left analysis output in the slots that synthesis takes as zeros; the bit test found
that bug in the spike.

**rustfft.** One plan per (N_φ, direction) from a fresh `FftPlanner`, run in place over the live
rings. Synthesis zeroes and scatters into the box. Analysis copies the box first, since its input
is kept, as FFTW's is.

**The product.** This is `tetraleaf_tables::helmholtz::diagonal::Conversions`, the shipped S → F
and F → R at scale 1. The φ step is fused into it and is not timed apart. Its time is taken as the
measured conversion times its φ share of the multiply-adds (`Conversions::multiply_adds`, 0.80;
CH4.1a's method). The column "conversion with the faster FFT ÷ product" is
((1 − share) · product + FFT) / product, a model of the conversion with the FFT in place of the
φ step.

**Grids.** For each degree p = 33, 41, 50, 66, 82, 100, 123, 135, 144, 200, 245, 261 (P from 34 to
262), the times use three grids:
- product: `azimuth_count(p)`;
- smooth: the least multiple of 4 at or above 4P − 3 whose N_φ / 4 is 7-smooth;
- pow2: the least power of two at or above 4P − 3.

**Timing.** Each time is per box on one thread, including the scatter and the gather: the median
of 7 samples of at least 0.05 s each, over a batch of 2–16 boxes of distinct random modes. Each
setting (back end and grid) ran in its own launch (`run.sh`; review 2026-11-24, B12). The BLAS
thread variables were at 1, in release builds.

**Accuracy.** The trusted path is the DFT written as a sum (`naive_synthesis`,
`naive_analysis`).

| Machine | CPU, SIMD seen | FFTW (`rlst::fftw::version()`) | Runs |
| --- | --- | --- | --- |
| locust | Neoverse-V2, NEON, SVE | `fftw-3.3.11-neon` (spack) | `tools/gh200/remote.sh 'spikes/fft-plans/run.sh locust'`, 2026-11-26 00:33–00:40 UTC, load average ≤ 2.6 falling (own build), GPU idle |
| Kathleen | Xeon Gold 6248, AVX2, AVX-512 | `fftw-3.3.10-sse2-avx-avx2-avx2_128-avx512` (module) | job 672793 (node-c11b-038, `singlenode`, slow nodes c11b-004 and c11b-013 excluded): tests, record, times; job 672796 (`test`, node-c11b-013 and node-c11b-031, `mpirun -n 2 --map-by ppr:1:node`): fingerprints |
| M6 | Apple M6, NEON | `fftw-3.3.11` (Homebrew, no SIMD) | tests, record, fingerprints (2 ranks); **no times**, at the user's request (2026-11-26: "Skip the M6") |

Raw outputs are in `spikes/fft-plans/results/<machine>/`: `tests.txt`, `record.md`,
`fingerprints.md`, `times-<variant>-<grid>.md` and `load.txt`. The tables in §5 come from
`tools/helmholtz/ch0-9-tables.py spikes/fft-plans/results/<machine>`.

Kathleen's node ran the times at a 1-minute load average falling from 9.4 to 1.6. That load is the
decay of the job's own parallel build: Kathleen's nodes are exclusive.

A first fingerprint job (672794) printed `HOSTNAME`, which Slurm passes from the batch node to
every rank, so it could not show the two nodes. Its 288 plans agreed. The spike now asks the system
for the host name, and job 672796 is the record.

## 2. Accuracy

`tests/bits.rs::accuracy_against_the_sum` covers p = 33, 82 and 261, all three grids, 2P − 1 and 2P
rings, and FFTW interleaved, FFTW split and rustfft, against the sum. The largest relative ℓ2 error
is 9.2·10⁻¹⁶ for synthesis and 1.6·10⁻¹⁵ for analysis, on locust and Kathleen alike, against the
test's bound of 10⁻¹³.

## 3. Bits per box

`tests/bits.rs::bits_across_positions_batches_and_threads` checks, for p = 33 and 82, the three
grids, both ring counts and all three back ends (36 cases): every box of packed batches of 1, 2, 7
and 33 boxes, all with the same input, transformed serially and on 1, 2, 4 and all threads (rayon
`par_bridge`, a scratch per worker), against the single box serially. Synthesis and analysis are
compared bit for bit, and every buffer is checked to be in class 0. All 36 cases are the same bits
on the M6, locust (72 threads) and Kathleen (40).

`bits_in_fresh_processes` hashes the 36 cases' outputs in two child processes and in the test's own
process. All three are equal on every machine.

## 4. For the record: what changes a box's bits

`examples/record.rs` compares synthesis of 16 boxes with distinct modes, at p = 33, 82, 144 and the
three grids, against the per-box FFTW plan of 2P − 1 rings on packed class-0 buffers
(`results/<machine>/record.md`):

| Variant | M6 (scalar FFTW) | locust (NEON) | Kathleen (AVX-512) |
| --- | --- | --- | --- |
| one plan over the batch, 1, 2, 7, 16 boxes | same bits | same bits | same bits |
| ring distance N_φ + 1, N_φ + 4 | same bits | same bits | same bits |
| alignment class 16 (planned for it, buffers one complex off a line) | same bits | same bits | same bits |
| 2P rings, on the live rings | same bits | same bits | same bits |
| split storage | same bits | other bits, ≤ 2.8·10⁻¹⁶ | other bits, ≤ 2.5·10⁻¹⁶ |
| rustfft | other bits, ≤ 7.1·10⁻¹⁶ | other bits, ≤ 7.8·10⁻¹⁶ | other bits, ≤ 7.6·10⁻¹⁶ |
| rustfft one complex off a line, against rustfft aligned | same bits | same bits | same bits |

FFTW's plans at p = 82 (in `record.md`) loop a plan of one ring over the rings
(`dft-vrank>=1-x165/1`, then Cooley–Tukey, and Bluestein for the prime 83 of N_φ = 332). The batch
dimension therefore never enters a SIMD codelet on these builds. RLST's design saw a plan for 288
transforms differ from one for 289 at length 588, so this is a property of these builds and sizes,
not a guarantee. The per-box plan of §6.7 costs nothing and holds whatever FFTW chooses.

## 5. Time per box

### locust (one core, µs per box)

#### φ step per box, S → F + F → R, µs

| p | N_φ product, smooth, pow2 | product whole | φ share | product φ | FFTW smooth | FFTW smooth 2P | FFTW split smooth | FFTW pow2 | FFTW product grid | rustfft smooth | rustfft pow2 | rustfft product grid |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 33 | 136, 140, 256 | 340 | 0.804 | 273 | 58.3 | 59.0 | 79.8 | 95.9 | 94.6 | 88.8 | 118 | 68.8 |
| 41 | 168, 168, 256 | 625 | 0.804 | 502 | 113 | 114 | 110 | 120 | 113 | 99.5 | 147 | 99.8 |
| 50 | 204, 216, 256 | 943 | 0.800 | 755 | 188 | 190 | 173 | 150 | 208 | 195 | 184 | 162 |
| 66 | 268, 280, 512 | 1895 | 0.800 | 1516 | 233 | 235 | 328 | 426 | 802 | 424 | 536 | 1003 |
| 82 | 332, 336, 512 | 3187 | 0.800 | 2550 | 474 | 477 | 482 | 540 | 2652 | 609 | 684 | 2087 |
| 100 | 404, 420, 512 | 5341 | 0.800 | 4273 | 743 | 745 | 929 | 673 | 2260 | 1140 | 848 | 2485 |
| 123 | 496, 500, 512 | 9895 | 0.801 | 7926 | 1014 | 1018 | 1304 | 842 | 1768 | 2188 | 1068 | 1364 |
| 135 | 544, 560, 1024 | 12510 | 0.801 | 10021 | 1258 | 1261 | 1496 | 1990 | 1698 | 1954 | 2500 | 1461 |
| 144 | 580, 588, 1024 | 13380 | 0.800 | 10704 | 1555 | 1558 | 2225 | 2132 | 2356 | 2229 | 2652 | 2295 |
| 200 | 804, 840, 1024 | 35065 | 0.800 | 28052 | 3243 | 3247 | 4151 | 3052 | 7967 | 5028 | 3791 | 9507 |
| 245 | 984, 1000, 1024 | 57504 | 0.801 | 46060 | 5695 | 5669 | 6297 | 3976 | 14378 | 9072 | 4752 | 14240 |
| 261 | 1048, 1080, 2048 | 71714 | 0.801 | 57443 | 5727 | 5763 | 7489 | 8232 | 27911 | 8061 | 10016 | 16858 |

#### Ratios

| p | FFTW smooth ÷ product φ | rustfft smooth ÷ product φ | FFTW 2P ÷ 2P − 1 | FFTW split ÷ interleaved | FFTW pow2 ÷ smooth | FFTW product grid ÷ smooth | rustfft pow2 ÷ smooth | rustfft product grid ÷ smooth | conversion with the faster FFT ÷ product |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 33 | 0.21 | 0.32 | 1.01 | 1.37 | 1.64 | 1.62 | 1.32 | 0.77 | 0.37 |
| 41 | 0.23 | 0.20 | 1.01 | 0.97 | 1.06 | 1.00 | 1.48 | 1.00 | 0.36 |
| 50 | 0.25 | 0.26 | 1.01 | 0.92 | 0.80 | 1.10 | 0.95 | 0.83 | 0.40 |
| 66 | 0.15 | 0.28 | 1.01 | 1.41 | 1.83 | 3.45 | 1.26 | 2.37 | 0.32 |
| 82 | 0.19 | 0.24 | 1.00 | 1.02 | 1.14 | 5.59 | 1.12 | 3.43 | 0.35 |
| 100 | 0.17 | 0.27 | 1.00 | 1.25 | 0.91 | 3.04 | 0.74 | 2.18 | 0.34 |
| 123 | 0.13 | 0.28 | 1.00 | 1.29 | 0.83 | 1.74 | 0.49 | 0.62 | 0.30 |
| 135 | 0.13 | 0.20 | 1.00 | 1.19 | 1.58 | 1.35 | 1.28 | 0.75 | 0.30 |
| 144 | 0.15 | 0.21 | 1.00 | 1.43 | 1.37 | 1.52 | 1.19 | 1.03 | 0.32 |
| 200 | 0.12 | 0.18 | 1.00 | 1.28 | 0.94 | 2.46 | 0.75 | 1.89 | 0.29 |
| 245 | 0.12 | 0.20 | 1.00 | 1.11 | 0.70 | 2.52 | 0.52 | 1.57 | 0.30 |
| 261 | 0.10 | 0.14 | 1.01 | 1.31 | 1.44 | 4.87 | 1.24 | 2.09 | 0.28 |


### Kathleen (one core, node-c11b-038, job 672793, µs per box)

#### φ step per box, S → F + F → R, µs

| p | N_φ product, smooth, pow2 | product whole | φ share | product φ | FFTW smooth | FFTW smooth 2P | FFTW split smooth | FFTW pow2 | FFTW product grid | rustfft smooth | rustfft pow2 | rustfft product grid |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 33 | 136, 140, 256 | 814 | 0.804 | 654 | 92.6 | 93.8 | 110 | 63.6 | 112 | 74.4 | 80.6 | 104 |
| 41 | 168, 168, 256 | 1386 | 0.804 | 1114 | 150 | 149 | 148 | 82.7 | 146 | 85.7 | 103 | 82.2 |
| 50 | 204, 216, 256 | 2111 | 0.800 | 1689 | 240 | 243 | 231 | 107 | 255 | 132 | 136 | 233 |
| 66 | 268, 280, 512 | 4084 | 0.800 | 3267 | 391 | 394 | 481 | 384 | 2081 | 278 | 445 | 1071 |
| 82 | 332, 336, 512 | 6977 | 0.800 | 5581 | 600 | 603 | 636 | 522 | 2098 | 435 | 585 | 1527 |
| 100 | 404, 420, 512 | 11526 | 0.800 | 9220 | 589 | 581 | 1355 | 741 | 5148 | 750 | 819 | 2900 |
| 123 | 496, 500, 512 | 20094 | 0.801 | 16095 | 927 | 928 | 1967 | 1017 | 2557 | 1603 | 1152 | 3022 |
| 135 | 544, 560, 1024 | 24745 | 0.801 | 19821 | 1760 | 1764 | 2174 | 2330 | 2310 | 1645 | 2693 | 2378 |
| 144 | 580, 588, 1024 | 29333 | 0.800 | 23466 | 1382 | 1385 | 2768 | 2492 | 3415 | 1930 | 2910 | 4102 |
| 200 | 804, 840, 1024 | 74935 | 0.800 | 59948 | 2818 | 2828 | 5445 | 3638 | 19629 | 4040 | 4160 | 10965 |
| 245 | 984, 1000, 1024 | 136023 | 0.801 | 108954 | 4446 | 4431 | 7905 | 4592 | 32546 | 6345 | 5269 | 16318 |
| 261 | 1048, 1080, 2048 | 162984 | 0.801 | 130550 | 4836 | 4830 | 10025 | 9443 | 32936 | 6714 | 12937 | 23852 |

#### Ratios

| p | FFTW smooth ÷ product φ | rustfft smooth ÷ product φ | FFTW 2P ÷ 2P − 1 | FFTW split ÷ interleaved | FFTW pow2 ÷ smooth | FFTW product grid ÷ smooth | rustfft pow2 ÷ smooth | rustfft product grid ÷ smooth | conversion with the faster FFT ÷ product |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 33 | 0.14 | 0.11 | 1.01 | 1.19 | 0.69 | 1.21 | 1.08 | 1.40 | 0.29 |
| 41 | 0.13 | 0.08 | 0.99 | 0.99 | 0.55 | 0.97 | 1.20 | 0.96 | 0.26 |
| 50 | 0.14 | 0.08 | 1.01 | 0.96 | 0.45 | 1.06 | 1.03 | 1.77 | 0.26 |
| 66 | 0.12 | 0.09 | 1.01 | 1.23 | 0.98 | 5.32 | 1.60 | 3.85 | 0.27 |
| 82 | 0.11 | 0.08 | 1.00 | 1.06 | 0.87 | 3.49 | 1.35 | 3.51 | 0.26 |
| 100 | 0.06 | 0.08 | 0.99 | 2.30 | 1.26 | 8.73 | 1.09 | 3.87 | 0.25 |
| 123 | 0.06 | 0.10 | 1.00 | 2.12 | 1.10 | 2.76 | 0.72 | 1.89 | 0.25 |
| 135 | 0.09 | 0.08 | 1.00 | 1.23 | 1.32 | 1.31 | 1.64 | 1.45 | 0.27 |
| 144 | 0.06 | 0.08 | 1.00 | 2.00 | 1.80 | 2.47 | 1.51 | 2.13 | 0.25 |
| 200 | 0.05 | 0.07 | 1.00 | 1.93 | 1.29 | 6.97 | 1.03 | 2.71 | 0.24 |
| 245 | 0.04 | 0.06 | 1.00 | 1.78 | 1.03 | 7.32 | 0.83 | 2.57 | 0.23 |
| 261 | 0.04 | 0.05 | 1.00 | 2.07 | 1.95 | 6.81 | 1.93 | 3.55 | 0.23 |


Read per row: "FFTW smooth ÷ product φ" is the FFT against the φ step it replaces; "conversion with the faster FFT ÷ product" the whole S → F + F → R modelled with it. The product on Kathleen is 2.0–2.4× its time on locust per core; the FFTs differ less (0.66–1.68×).

## 6. The outputs §10 asks for

- **The plan layout and scratch of §6.7: confirmed.** One plan per (level, N_φ) of 2P − 1 rings,
  `Rigor::Estimate`, executed box by box on class-0 buffers packed back to back, with interleaved
  scratch. The per-thread scratch needs two boxes (synthesis input kept zero outside the mode
  slots, and analysis output), or a re-zeroing.
- **The batch: 2P − 1.** 2P rings ran 0.99–1.01× on both machines and changed no bits.
- **The grid.** An FFT back end needs N_φ / 4 7-smooth, which is §6.7's own rule. The product as
  built at CH4.1a uses `azimuth_count` (a multiple of 4 only), on which FFTW and rustfft lose up to
  8.7× and 3.9×. A level with an FFT back end therefore takes the smooth N_φ (at most 6 % more
  samples at these p: 140/136, 216/204, 588/580, 1080/1048), and with it its own transfer
  functions in the table cache. §6.7's "Grid sizes" should say that the product runs on the
  multiple of 4 (CH4.1a) and an FFT level on the smooth size (question 3).
- **`rlst::fftw::version()` per machine:** see §1. These are the same as CHR.4's records.
- **D7, per machine (question 1).** The proposal's recommendation and its condition for rustfft
  were met on both measured machines.

  | Machine | Product's φ step | FFTW (smooth) | rustfft (smooth) | Proposal |
  | --- | --- | --- | --- | --- |
  | locust | 1× | 0.10–0.25× | 0.14–0.32× | product default; rustfft the optional FFT |
  | Kathleen | 1× | 0.04–0.14× | 0.05–0.11× | product default; rustfft the optional FFT |
  | M6 | not timed (the user) | scalar FFTW (CHR.4) | not timed | product default; rustfft deterministic per box here too (§3); its time to be measured when the back end is built |

  The FFT saves about 1–2 % of a run (CH4.1a). The product therefore stays the default and the
  reference: it has no dependency, and its bits are those of every fixture and test today. For the
  optional back end, rustfft gives nearly all of FFTW's gain, without a system library, the GPL,
  pkg-config or the fftw-sys licence entry. FFTW's wrapper (CHR.3) and its back end (CH4.6) would
  then leave the priority list, and HR's gate would drop its third item.

  The alternatives:
  - as the proposal stands: FFTW behind `fftw` (CHR.3, CH4.6), rustfft not adopted;
  - both back ends behind their features;
  - neither: the FFT track leaves the list after S9, the product alone.

## 7. Caveats

- The product's φ step is not measured apart: 0.80 of the conversion is its share of the
  multiply-adds, not of the time. CH4.1a noted that the azimuths' fold, which an FFT also replaces,
  is not counted, so the real share is somewhat larger and the FFT's ratios somewhat smaller.
- Per box on one thread, with no other work on the node. A level's conversions run in parallel
  under memory bandwidth, where the FFT's larger buffers (N_φ R samples written and read) may cost
  more than here.
- rustfft has no plan description or fingerprint. A back end on it needs its own check across ranks,
  for example rustfft's version, N_φ and the CPU features it dispatches on (`simd_features`),
  all-gathered at build as §6.7 does for FFTW (question 1, its condition).
- The M6's times were not run (the user). Its bits (§3, §4) and fingerprints were.
