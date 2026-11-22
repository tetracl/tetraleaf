# Helmholtz CH4.3e: c_leaf re-decided with the time rule

Task: `docs/helmholtz/tasks/CH4.3e.md`. Roadmap: H4, CH4.3e (priority list item 10). Proposal §3.7
(leaves and tree depth, D11), §5.1, §5.7, §9.1; §2 (constraints 1, 2, 9). Builds on CH4.0d3
(c_leaf 2, `results/helmholtz-ch4-0d3.md` §4) and CH4.3b2c (the diagonal form's time rule the
builder's default, `results/helmholtz-ch4-3b2c.md`).

## 0. Summary

- **Leaves on high-frequency levels are the rule at every c_leaf.** A leaf of c_leaf λ has
  k a ≤ 5.4 c_leaf (a = √3 w/2), so even c_leaf 1's leaves (k a ≤ 5.4) reach the diagonal form's
  range, k a ≥ 2.8 (S6's, the time rule's floor). Counting every level at k a ≥ 2.8, whatever its
  M2L, 8–9 of Kathleen's 8–10 trees per ε and all 5 of locust's have occupied leaves there, at every
  c_leaf from 1 to 8 and every ε from 10⁻³ to 10⁻⁸ (§2; the exception the resolved W3 at kD = 100,
  whose leaves lie at k a 2.7). Their M2L is the form (the time rule's choice), RCR, the plane waves
  (D8's range reaches k a 27) or none (every V pair near). The census's first pass
  (`spikes/census`, D11's reason for c_leaf 1) classified by a bound without S6's constants; no
  c_leaf of this sweep keeps the leaves off such levels.
- **They are correct.** The leaf operators (P2M, L2P, P2L, M2P) do not depend on the M2L's
  regime: they run at the level's degree, which the degree rule covers to x = 216 (CH4.0c,
  CH4.3c2), and leaves of c_leaf 4 have x ≤ 21.8. Every one of the 456 runs (both machines,
  W1–W3, 10⁻³–10⁻⁸, c_leaf 1–8, both orders) erred at most 0.27 ε in φ and 0.44 ε in ∇φ (§3),
  and each order's output per c_leaf is the other's bit for bit. W2's trees put W and X lists on
  those leaves (M2P and P2L touching them: up to 57,048 and 70,728 interactions per tree, §2).
  The new test `check_form_leaves` holds the case on a small tree (§5).
- **c_leaf 4 the proposed default** (built here, for sign-off, question 2). Build plus one
  evaluation over c_leaf 2's, geometric means over the cases of both run orders (§4):

  | machine | cases | c_leaf 1 | 4 | 8 | without W3 at kD = 500: 1 | 4 | 8 |
  | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
  | Kathleen (1 × 40) | 37 | 1.417 | 0.874 | 0.849 | 1.247 | 0.995 | 0.987 |
  | locust (1 × 72) | 20 | 1.518 | 0.790 | 0.780 | 1.148 | 1.000 | 1.011 |

  The gain is the under-resolved W3 (N = 10⁶ at kD = 500), where c_leaf 2's leaves at level 6
  hold about 4 points: 0.29–0.32× on Kathleen and 0.26–0.36× on locust at c_leaf 4. Elsewhere
  c_leaf 4 changes only W1 and W2 at kD ≥ 500 on Kathleen (0.93–0.98×, a few coarser leaves where
  the sphere grazes a box) and W1 at kD = 500, 10⁻⁸ on locust (0.98×); on every other case the tree
  is c_leaf 2's and the ratio the run-order spread (0.94–1.11 on Kathleen, 0.99–1.03 on locust).
  c_leaf 8 changes only W3 at kD = 500 from c_leaf 4 (0.22–0.27× on Kathleen, 0.23× on locust at
  10⁻⁸), 1–3 % more on the means, at twice the phase range (Re(k) r ≤ 352) and leaf x (43.5);
  c_leaf 1 costs 1.4–1.5× (1.3–1.8× on W1 at kD ≥ 500, 3.7–4.8× on W3 at kD = 500).
- **The proposal's §3.7 and §5.1 are wrong as built** (question 3): "every leaf is on a
  low-frequency level, and so are the W- and X-list operators" holds at no c_leaf ≥ 1, and
  "under-resolved inputs get more, sparser leaves instead of high-frequency ones" is the price
  c_leaf 4 now declines (its gain is that case). The correction proposed: leaves may lie on levels
  of k a ≥ 2.8 whatever their M2L; the leaf operators stay correct there because they do not
  depend on the regime and the degree rule covers their arguments; c_leaf bounds the P2P's phase
  and the leaves' argument, not the regime.

## 1. Setup

- **Runs.** `helmholtz_accuracy` (release, BLAS threads 1, one rank × all cores), launched by
  `tools/helmholtz/ch4-3e-runs.sh` at fc3b38d: §9.1's complex charges and 1,000 samples,
  gradients, two evaluations (the rows the second's), the library's defaults named (`--m2l
  auto-pw`, `--diagonal auto --diagonal-build interpolated`, the V, W and X rules of the cost
  model, CH4.3b2d's `PerBox`, CH4.3c2's degrees, the leaf-size rule: 256 points on x86, 128 on
  NEON at ε ≥ 10⁻⁷ and 256 below), each case at `--c-leaf 1,2,4,8` (`fwd-*`) and `8,4,2,1`
  (`rev-*`) in one launch each. The ratios are the geometric mean of the two launches' ratios.
- **Cases.** W1 at kD = 100, 200, 500 and W2 at kD = 100, 200, 500 at 10⁻³, 10⁻⁶, 10⁻⁷, 10⁻⁸; W3 at
  N = 10⁶, kD = 200 and 500 at the same ε; on Kathleen also W1 at kD = 1000 at 10⁻³ and 10⁻⁶ and W3
  resolved at kD = 100 (N = 4.0·10⁶) at 10⁻³, 10⁻⁶, 10⁻⁸. Kathleen 37 cases, locust 20 (W1 and W3).
- **Machines.**
  - Kathleen, `singlenode`, one rank bound to the node's 40 cores (`--map-by slot:PE=40 --bind-to
    core`): job 663653 (node-c11b-012, 2026-11-22 14:19–14:53, W1 and W1 at kD = 1000), 663654
    (node-c11a-021, 14:20–14:34, W2), 663655 (node-c11a-021, 14:34–15:18, W3 and W3 resolved);
    every launch exit 0 (74).
  - locust (GH200, 72 cores), 14:21–14:40, W1 and W3, every launch exit 0 (40); idle before the
    start (1-minute load 0.03, no other process above 1.1 % CPU), the 1-minute load at most 44.7
    during the runs, below their own 72 threads (sampled every 30 s; the processes were not
    sampled during the runs).
- **Raw outputs**: `~/projects/fmm-archive/kathleen/runs/ch4-3e/` (`ch4-3e-w1`, `-w2`, `-w3`, the
  jobs' logs) and `~/projects/fmm-archive/locust/runs/ch4-3e/` (with `load.log`).
- **Tables**: `tools/helmholtz/ch4-3e-tables.py census|summary|cases|errors`.

## 2. The census: leaves on high-frequency levels

Per tree (c_leaf and ε), whether a level of k a = Re(k) √3 r_l ≥ 2.8 holds occupied leaves (the
table "Leaves"), whatever the M2L that ran there: the form (`diagonal l:c` in the column "M2L"),
RCR, the plane waves or Batched (the table "Time rule"), or none ("no V": W3 at kD = 500 at
c_leaf 1 and 2, whose leaves of 1–4 points make every V pair of the leaf level a near pair under
the V rule).

| machine | c_leaf | ε | trees | with leaves on a high-frequency level |
| --- | ---: | --- | ---: | ---: |
| kathleen | 1 | 1e-3 | 10 | 9 |
| kathleen | 1 | 1e-6 | 10 | 9 |
| kathleen | 1 | 1e-7 | 8 | 8 |
| kathleen | 1 | 1e-8 | 9 | 8 |
| kathleen | 2 | 1e-3 | 10 | 9 |
| kathleen | 2 | 1e-6 | 10 | 9 |
| kathleen | 2 | 1e-7 | 8 | 8 |
| kathleen | 2 | 1e-8 | 9 | 8 |
| kathleen | 4 | 1e-3 | 10 | 9 |
| kathleen | 4 | 1e-6 | 10 | 9 |
| kathleen | 4 | 1e-7 | 8 | 8 |
| kathleen | 4 | 1e-8 | 9 | 8 |
| kathleen | 8 | 1e-3 | 10 | 9 |
| kathleen | 8 | 1e-6 | 10 | 9 |
| kathleen | 8 | 1e-7 | 8 | 8 |
| kathleen | 8 | 1e-8 | 9 | 8 |
| locust | 1 | 1e-3 | 5 | 5 |
| locust | 1 | 1e-6 | 5 | 5 |
| locust | 1 | 1e-7 | 5 | 5 |
| locust | 1 | 1e-8 | 5 | 5 |
| locust | 2 | 1e-3 | 5 | 5 |
| locust | 2 | 1e-6 | 5 | 5 |
| locust | 2 | 1e-7 | 5 | 5 |
| locust | 2 | 1e-8 | 5 | 5 |
| locust | 4 | 1e-3 | 5 | 5 |
| locust | 4 | 1e-6 | 5 | 5 |
| locust | 4 | 1e-7 | 5 | 5 |
| locust | 4 | 1e-8 | 5 | 5 |
| locust | 8 | 1e-3 | 5 | 5 |
| locust | 8 | 1e-6 | 5 | 5 |
| locust | 8 | 1e-7 | 5 | 5 |
| locust | 8 | 1e-8 | 5 | 5 |

Selected trees on Kathleen (every row: `ch4-3e-tables.py census`): the leaf levels, the
high-frequency levels with leaves (the M2L, k a, the occupied leaves), the W interactions into
them (M2P at that level), the X interactions out of them (P2L at the next finer level), the phase
bound Re(k) × the farthest near pair's point distance (`Row::phase_bound`):

| machine | case | c_leaf | leaf levels (occupied) | levels of k a ≥ 2.8 with leaves: (M2L, k a): occupied | occupied there | W into them | X out of them | Re(k) r bound |
| --- | --- | ---: | --- | --- | ---: | ---: | ---: | ---: |
| kathleen | w1-1000-1e-6 | 1 | 8 | 8 (rcr, k a 3.4): 291849 | 291849 | 0 | 0 | 27.1 |
| kathleen | w1-1000-1e-6 | 2 | 7 | 7 (form 34, k a 6.8): 75295 | 75295 | 0 | 0 | 54.1 |
| kathleen | w1-1000-1e-6 | 4 | 6–7 | 6 (form 34, k a 13.5): 6252; 7 (form 34, k a 6.8): 64756 | 71008 | 0 | 0 | 108.3 |
| kathleen | w1-1000-1e-6 | 8 | 6–7 | 6 (form 34, k a 13.5): 6252; 7 (form 34, k a 6.8): 64756 | 71008 | 0 | 0 | 108.3 |
| kathleen | w2-100-1e-3 | 1 | 4–6 | 4 (form 34, k a 5.4): 1072 | 1072 | 9705 | 14248 | 43.3 |
| kathleen | w2-100-1e-3 | 2 | 3–6 | 3 (form 34, k a 10.8): 76; 4 (form 34, k a 5.4): 908 | 984 | 9705 | 14248 | 86.6 |
| kathleen | w2-100-1e-3 | 4 | 3–6 | 3 (form 34, k a 10.8): 76; 4 (form 34, k a 5.4): 908 | 984 | 9705 | 14248 | 86.6 |
| kathleen | w2-100-1e-3 | 8 | 3–6 | 3 (form 34, k a 10.8): 76; 4 (form 34, k a 5.4): 908 | 984 | 9705 | 14248 | 86.6 |
| kathleen | w2-500-1e-6 | 1 | 7–8 | 7 (rcr, k a 3.4): 70172 | 70172 | 0 | 0 | 27.1 |
| kathleen | w2-500-1e-6 | 2 | 6–8 | 6 (form 34, k a 6.8): 17556; 7 (form 21, k a 3.4): 2238 | 19794 | 34557 | 34557 | 54.1 |
| kathleen | w2-500-1e-6 | 4 | 5–8 | 5 (form 34, k a 13.5): 1433; 6 (form 34, k a 6.8): 15212; 7 (form 21, k a 3.4): 2238 | 18883 | 57048 | 57048 | 108.3 |
| kathleen | w2-500-1e-6 | 8 | 5–8 | 5 (form 34, k a 13.5): 1433; 6 (form 34, k a 6.8): 15212; 7 (form 21, k a 3.4): 2238 | 18883 | 57048 | 57048 | 108.3 |
| kathleen | w3-500-1e-6 | 1 | 7 | 7 (no V, k a 3.4): 795931 | 795931 | 0 | 0 | 27.1 |
| kathleen | w3-500-1e-6 | 2 | 6 | 6 (no V, k a 6.8): 256369 | 256369 | 0 | 0 | 54.1 |
| kathleen | w3-500-1e-6 | 4 | 5 | 5 (rcr, k a 13.5): 32768 | 32768 | 0 | 0 | 108.3 |
| kathleen | w3-500-1e-6 | 8 | 4–5 | 4 (form 34, k a 27.1): 3226; 5 (rcr, k a 13.5): 6960 | 10186 | 0 | 0 | 216.5 |
| kathleen | w3-500-1e-8 | 1 | 7 | 7 (no V, k a 3.4): 795931 | 795931 | 0 | 0 | 27.1 |
| kathleen | w3-500-1e-8 | 2 | 6 | 6 (no V, k a 6.8): 256369 | 256369 | 0 | 0 | 54.1 |
| kathleen | w3-500-1e-8 | 4 | 5 | 5 (no V, k a 13.5): 32768 | 32768 | 0 | 0 | 108.3 |
| kathleen | w3-500-1e-8 | 8 | 4–5 | 4 (waves, k a 27.1): 3226; 5 (no V, k a 13.5): 6960 | 10186 | 0 | 0 | 216.5 |

So: occupied leaves lie at k a 3.4–27 under every M2L the defaults run; on W1 the W and X pairs
there are near pairs (the W and X rule converts them) or absent, on W2's cap they are not, and M2P
and P2L run on the form's leaves (W2 at kD = 500, 10⁻⁶, c_leaf 4: 57,048 interactions each way).
The phase bound stays within the 44 c_leaf of converted V pairs (at most 43.3 c_leaf, W1 at
kD = 100 at c_leaf 1) and reaches 108.3 rad at c_leaf 4 and 216.5 at 8, inside the range the
Helmholtz P2P's sine and cosine are tested on (|θ| ≤ 2¹⁰, `tetraleaf-simd`'s `sincos`).

## 3. Errors

Relative ℓ2 errors at §9.1's samples, the worst over the runs, and each c_leaf's over c_leaf 2's
in the same (forward) launch; and whether every case's two launches gave each c_leaf the same
output (the checksum of the whole output):

| machine | ε | runs | φ ÷ ε max | ∇φ ÷ ε max | φ ÷ c_leaf 2's: c_leaf 1 | 4 | 8 | ∇φ ÷ c_leaf 2's: 1 | 4 | 8 | every order and run bit for bit the same output per c_leaf |
| --- | --- | ---: | ---: | ---: | --- | --- | --- | --- | --- | --- | --- |
| kathleen | 1e-3 | 80 | 0.12 | 0.16 | 1.000–1.048 | 0.999–1.000 | 0.982–1.000 | 1.000–1.346 | 0.969–1.000 | 0.939–1.000 | yes |
| kathleen | 1e-6 | 80 | 0.16 | 0.37 | 1.000–1.252 | 0.744–1.000 | 0.715–1.000 | 1.000–1.613 | 0.533–1.000 | 0.481–1.000 | yes |
| kathleen | 1e-7 | 64 | 0.13 | 0.22 | 0.938–1.000 | 0.991–1.010 | 0.647–1.010 | 0.842–1.000 | 0.976–1.000 | 0.647–1.000 | yes |
| kathleen | 1e-8 | 72 | 0.27 | 0.44 | 1.000–1.000 | 1.000–1.000 | 0.774–1.000 | 1.000–1.000 | 1.000–1.000 | 0.819–1.000 | yes |
| locust | 1e-3 | 40 | 0.12 | 0.16 | 1.000–1.036 | 1.000–1.000 | 1.000–1.000 | 1.000–1.237 | 0.985–1.000 | 0.985–1.000 | yes |
| locust | 1e-6 | 40 | 0.16 | 0.38 | 0.994–1.083 | 0.783–1.000 | 0.783–1.000 | 1.000–1.179 | 0.566–1.000 | 0.566–1.000 | yes |
| locust | 1e-7 | 40 | 0.12 | 0.18 | 0.987–1.000 | 1.000–1.052 | 1.000–1.052 | 0.877–1.000 | 1.000–1.071 | 1.000–1.071 | yes |
| locust | 1e-8 | 40 | 0.27 | 0.44 | 1.000–1.000 | 1.000–1.000 | 0.774–1.000 | 1.000–1.000 | 1.000–1.000 | 0.819–1.000 | yes |

Every run within ε by a factor of about 2.3 at least. Where a c_leaf changes the tree the errors
move both ways, 0.48–1.61× c_leaf 2's: c_leaf 1 up to 1.18–1.61× in ∇φ (W1 at kD = 500 and 1000,
W3 at kD = 200 and 500, at 10⁻³ and 10⁻⁶), c_leaf 4 and 8 down to 0.53 and 0.48× (W3 at kD = 500,
10⁻⁶: coarser leaves, more of the field exact in P2P).

## 4. Times

Build plus one evaluation and the evaluation alone over c_leaf 2's, geometric means over the cases
(the mean of both orders per case), and the cases whose output at c_leaf 4 or 8 differs from
c_leaf 2's (the tree changed):

| machine | ε | cases | build + evaluate ÷ c_leaf 2: c_leaf 1 | 4 | 8 | evaluate ÷ c_leaf 2: c_leaf 1 | 4 | 8 | cases where 4 or 8 changes the tree |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| kathleen | 1e-3 | 10 | 1.387 | 0.865 | 0.847 | 1.017 | 0.928 | 0.933 | 4 |
| kathleen | 1e-6 | 10 | 1.422 | 0.881 | 0.856 | 1.146 | 0.947 | 0.937 | 4 |
| kathleen | 1e-7 | 8 | 1.469 | 0.852 | 0.833 | 1.297 | 0.898 | 0.894 | 3 |
| kathleen | 1e-8 | 9 | 1.397 | 0.894 | 0.858 | 1.280 | 0.939 | 0.907 | 3 |
| kathleen | all | 37 | 1.417 | 0.874 | 0.849 | 1.171 | 0.929 | 0.919 | 14 |
| locust | 1e-3 | 5 | 1.443 | 0.760 | 0.772 | 1.065 | 0.981 | 0.982 | 1 |
| locust | 1e-6 | 5 | 1.471 | 0.779 | 0.787 | 1.186 | 0.979 | 0.978 | 1 |
| locust | 1e-7 | 5 | 1.455 | 0.809 | 0.813 | 1.293 | 0.982 | 0.980 | 1 |
| locust | 1e-8 | 5 | 1.722 | 0.814 | 0.749 | 1.380 | 0.992 | 0.926 | 2 |
| locust | all | 20 | 1.518 | 0.790 | 0.780 | 1.225 | 0.984 | 0.966 | 5 |

Per case (the build and evaluation of c_leaf 2 the geometric mean over both launches):

| machine | case | N | leaf size | leaf level at c_leaf 1, 2, 4, 8 | c_leaf 2: build s | evaluate s | build + evaluate ÷ c_leaf 2: 1 | 4 | 8 | evaluate ÷: 1 | 4 | 8 |
| --- | --- | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| kathleen | w1-100-1e-3 | 79578 | 256 | 4, 3, 2, 1 | 0.11 | 0.100 | 1.07 | 0.97 | 1.00 | 0.85 | 0.95 | 0.95 |
| kathleen | w1-100-1e-6 | 79578 | 256 | 4, 3, 2, 1 | 0.13 | 0.130 | 1.04 | 1.02 | 0.98 | 0.93 | 1.04 | 0.96 |
| kathleen | w1-100-1e-7 | 79578 | 256 | 4, 3, 2, 1 | 0.14 | 0.215 | 1.06 | 1.01 | 1.07 | 0.95 | 1.02 | 1.07 |
| kathleen | w1-100-1e-8 | 79578 | 256 | 4, 3, 2, 1 | 0.13 | 0.285 | 1.06 | 1.11 | 1.04 | 0.98 | 1.14 | 1.00 |
| kathleen | w1-200-1e-3 | 318310 | 256 | 5, 4, 3, 2 | 0.47 | 0.429 | 1.02 | 0.98 | 0.99 | 0.83 | 0.97 | 0.96 |
| kathleen | w1-200-1e-6 | 318310 | 256 | 5, 4, 3, 2 | 0.56 | 0.539 | 1.09 | 0.99 | 0.98 | 1.01 | 1.02 | 1.04 |
| kathleen | w1-200-1e-7 | 318310 | 256 | 5, 4, 3, 2 | 0.56 | 0.800 | 1.12 | 1.03 | 1.12 | 1.03 | 1.02 | 1.14 |
| kathleen | w1-200-1e-8 | 318310 | 256 | 5, 4, 3, 2 | 0.57 | 1.140 | 1.15 | 1.07 | 1.12 | 1.11 | 1.09 | 1.15 |
| kathleen | w1-500-1e-3 | 1989437 | 256 | 7, 6, 5, 4 | 4.76 | 3.005 | 1.66 | 0.93 | 0.93 | 0.94 | 0.97 | 1.02 |
| kathleen | w1-500-1e-6 | 1989437 | 256 | 7, 6, 5, 4 | 5.03 | 3.790 | 1.69 | 0.95 | 0.94 | 1.15 | 0.99 | 0.98 |
| kathleen | w1-500-1e-7 | 1989437 | 256 | 7, 6, 5, 4 | 5.05 | 5.005 | 1.77 | 0.98 | 0.94 | 1.46 | 1.01 | 0.98 |
| kathleen | w1-500-1e-8 | 1989437 | 256 | 7, 6, 5, 4 | 5.16 | 6.635 | 1.73 | 0.96 | 0.95 | 1.47 | 0.98 | 0.99 |
| kathleen | w1-1000-1e-3 | 7957748 | 256 | 8, 7, 6, 5 | 29.66 | 14.370 | 1.72 | 0.95 | 0.95 | 0.96 | 1.00 | 1.01 |
| kathleen | w1-1000-1e-6 | 7957748 | 256 | 8, 7, 6, 5 | 30.90 | 17.354 | 1.73 | 0.98 | 0.98 | 1.18 | 1.06 | 1.09 |
| kathleen | w2-100-1e-3 | 127554 | 256 | 4, 3, 2, 1 | 0.16 | 0.120 | 1.04 | 0.96 | 1.02 | 0.87 | 0.92 | 1.00 |
| kathleen | w2-100-1e-6 | 127554 | 256 | 4, 3, 2, 1 | 0.16 | 0.194 | 0.97 | 1.00 | 0.96 | 0.85 | 1.00 | 0.92 |
| kathleen | w2-100-1e-7 | 127554 | 256 | 4, 3, 2, 1 | 0.17 | 0.335 | 1.07 | 1.01 | 1.00 | 1.03 | 1.02 | 0.97 |
| kathleen | w2-100-1e-8 | 127554 | 256 | 4, 3, 2, 1 | 0.16 | 0.465 | 1.02 | 1.02 | 0.99 | 0.98 | 1.02 | 0.98 |
| kathleen | w2-200-1e-3 | 510214 | 256 | 5, 4, 3, 2 | 0.59 | 0.410 | 1.07 | 1.00 | 1.10 | 0.93 | 1.01 | 1.20 |
| kathleen | w2-200-1e-6 | 510214 | 256 | 5, 4, 3, 2 | 0.62 | 0.645 | 1.22 | 0.99 | 1.05 | 1.20 | 1.00 | 1.08 |
| kathleen | w2-200-1e-7 | 510214 | 256 | 5, 4, 3, 2 | 0.68 | 1.253 | 1.08 | 1.03 | 0.98 | 1.02 | 1.05 | 0.97 |
| kathleen | w2-200-1e-8 | 510214 | 256 | 5, 4, 3, 2 | 0.67 | 1.807 | 0.99 | 1.07 | 1.02 | 0.94 | 1.08 | 1.02 |
| kathleen | w2-500-1e-3 | 3188840 | 256 | 7, 6, 5, 4 | 4.54 | 2.905 | 1.47 | 0.97 | 0.93 | 0.93 | 1.04 | 0.98 |
| kathleen | w2-500-1e-6 | 3188840 | 256 | 7, 6, 5, 4 | 4.70 | 4.345 | 1.48 | 0.95 | 0.94 | 1.09 | 0.96 | 0.95 |
| kathleen | w2-500-1e-7 | 3188840 | 256 | 7, 6, 5, 4 | 4.83 | 6.345 | 1.50 | 0.96 | 0.96 | 1.27 | 0.98 | 0.98 |
| kathleen | w2-500-1e-8 | 3188840 | 256 | 7, 6, 5, 4 | 4.89 | 8.585 | 1.49 | 0.97 | 0.97 | 1.30 | 0.99 | 0.99 |
| kathleen | w3-200-1e-3 | 1000000 | 256 | 5, 4, 3, 2 | 0.90 | 1.918 | 1.43 | 1.01 | 0.93 | 1.05 | 1.02 | 0.89 |
| kathleen | w3-200-1e-6 | 1000000 | 256 | 5, 4, 3, 2 | 0.93 | 2.290 | 1.51 | 1.01 | 0.92 | 1.02 | 1.01 | 0.89 |
| kathleen | w3-200-1e-7 | 1000000 | 256 | 5, 4, 3, 2 | 0.94 | 2.710 | 1.37 | 0.94 | 0.93 | 0.96 | 0.93 | 0.90 |
| kathleen | w3-200-1e-8 | 1000000 | 256 | 5, 4, 3, 2 | 0.94 | 3.125 | 1.50 | 1.02 | 0.96 | 1.17 | 1.03 | 0.95 |
| kathleen | w3-500-1e-3 | 1000000 | 256 | 7, 6, 5, 4 | 17.14 | 7.250 | 3.74 | 0.29 | 0.22 | 2.55 | 0.57 | 0.55 |
| kathleen | w3-500-1e-6 | 1000000 | 256 | 7, 6, 5, 4 | 17.24 | 9.650 | 3.85 | 0.32 | 0.27 | 2.72 | 0.54 | 0.58 |
| kathleen | w3-500-1e-7 | 1000000 | 256 | 7, 6, 5, 4 | 17.66 | 14.045 | 4.38 | 0.29 | 0.24 | 4.38 | 0.41 | 0.41 |
| kathleen | w3-500-1e-8 | 1000000 | 256 | 7, 6, 5, 4 | 17.74 | 17.519 | 4.25 | 0.30 | 0.24 | 4.11 | 0.41 | 0.39 |
| kathleen | w3r-100-1e-3 | 4031442 | 256 | 4, 3, 2, 1 | 3.87 | 2.890 | 0.97 | 0.99 | 0.99 | 0.93 | 0.95 | 0.94 |
| kathleen | w3r-100-1e-6 | 4031442 | 256 | 4, 3, 2, 1 | 3.86 | 4.163 | 1.00 | 1.00 | 1.01 | 1.01 | 0.99 | 1.00 |
| kathleen | w3r-100-1e-8 | 4031442 | 256 | 4, 3, 2, 1 | 3.84 | 9.610 | 1.00 | 1.01 | 1.01 | 0.99 | 1.00 | 1.00 |
| locust | w1-100-1e-3 | 79578 | 128 | 4, 3, 2, 1 | 0.06 | 0.040 | 1.00 | 1.00 | 1.05 | 0.75 | 1.00 | 1.00 |
| locust | w1-100-1e-6 | 79578 | 128 | 4, 3, 2, 1 | 0.07 | 0.040 | 1.04 | 1.00 | 1.04 | 1.00 | 1.00 | 1.00 |
| locust | w1-100-1e-7 | 79578 | 128 | 4, 3, 2, 1 | 0.07 | 0.090 | 1.03 | 1.03 | 1.03 | 1.00 | 1.00 | 1.00 |
| locust | w1-100-1e-8 | 79578 | 256 | 4, 3, 2, 1 | 0.06 | 0.110 | 1.09 | 1.00 | 1.03 | 1.00 | 1.00 | 1.00 |
| locust | w1-200-1e-3 | 318310 | 128 | 5, 4, 3, 2 | 0.26 | 0.130 | 1.03 | 0.99 | 1.01 | 0.92 | 1.00 | 1.00 |
| locust | w1-200-1e-6 | 318310 | 128 | 5, 4, 3, 2 | 0.28 | 0.160 | 1.05 | 1.00 | 1.01 | 0.94 | 1.00 | 1.00 |
| locust | w1-200-1e-7 | 318310 | 128 | 5, 4, 3, 2 | 0.29 | 0.330 | 1.01 | 0.99 | 1.01 | 0.94 | 1.00 | 1.00 |
| locust | w1-200-1e-8 | 318310 | 256 | 5, 4, 3, 2 | 0.24 | 0.420 | 1.06 | 1.01 | 1.01 | 0.93 | 1.01 | 1.00 |
| locust | w1-500-1e-3 | 1989437 | 128 | 7, 6, 5, 4 | 2.74 | 0.775 | 1.34 | 0.99 | 0.99 | 0.81 | 1.00 | 1.00 |
| locust | w1-500-1e-6 | 1989437 | 128 | 7, 6, 5, 4 | 2.79 | 1.165 | 1.38 | 1.02 | 1.02 | 1.02 | 1.00 | 1.00 |
| locust | w1-500-1e-7 | 1989437 | 128 | 7, 6, 5, 4 | 2.86 | 1.720 | 1.32 | 1.00 | 1.01 | 1.03 | 1.01 | 1.01 |
| locust | w1-500-1e-8 | 1989437 | 256 | 7, 6, 5, 4 | 1.74 | 2.175 | 1.68 | 0.98 | 0.97 | 1.07 | 1.04 | 1.04 |
| locust | w3-200-1e-3 | 1000000 | 128 | 5, 4, 3, 2 | 1.25 | 0.580 | 1.01 | 1.00 | 1.01 | 1.00 | 1.00 | 1.00 |
| locust | w3-200-1e-6 | 1000000 | 128 | 5, 4, 3, 2 | 1.27 | 0.680 | 1.00 | 1.00 | 0.99 | 1.00 | 1.00 | 0.99 |
| locust | w3-200-1e-7 | 1000000 | 128 | 5, 4, 3, 2 | 1.37 | 1.385 | 1.00 | 1.00 | 1.00 | 1.00 | 1.01 | 1.00 |
| locust | w3-200-1e-8 | 1000000 | 256 | 5, 4, 3, 2 | 0.56 | 1.450 | 1.67 | 1.00 | 1.00 | 1.36 | 1.00 | 0.99 |
| locust | w3-500-1e-3 | 1000000 | 128 | 7, 6, 5, 4 | 9.70 | 1.655 | 4.52 | 0.26 | 0.26 | 2.44 | 0.91 | 0.91 |
| locust | w3-500-1e-6 | 1000000 | 128 | 7, 6, 5, 4 | 9.96 | 2.095 | 4.57 | 0.28 | 0.28 | 2.45 | 0.90 | 0.91 |
| locust | w3-500-1e-7 | 1000000 | 128 | 7, 6, 5, 4 | 9.80 | 3.025 | 4.79 | 0.34 | 0.34 | 3.73 | 0.90 | 0.89 |
| locust | w3-500-1e-8 | 1000000 | 256 | 7, 6, 5, 4 | 9.90 | 3.515 | 4.69 | 0.36 | 0.23 | 3.71 | 0.91 | 0.66 |

The leaf level in the table is the leaf rule's (`min_leaf_level`); the tree's leaves are split
further by points. Where c_leaf 4 builds c_leaf 2's tree (no checksum changes) its ratio is the
run-order spread: 0.94–1.11 on Kathleen (CH4.3b2c's report §6: 9–20 % between runs of one
launch), 0.99–1.03 on locust. The cases above 1.05 at c_leaf 4 are all of that kind (the tree
c_leaf 2's): W1 at kD = 100, 10⁻⁸ (1.11, 0.42 s in all), W1 and W2 at kD = 200, 10⁻⁸ (1.07); at
c_leaf 8 W1 at kD = 200, 10⁻⁷ and 10⁻⁸ (1.12), W2 at kD = 200, 10⁻³ (1.10) and W1 at kD = 100,
10⁻⁷ (1.07). By CH4.0d3's rule (another value where its geometric mean over a machine class's
cases is at least 5 % lower, errors within ε), c_leaf 4 and 8 both qualify on both classes; 1
does not. The whole margin is the under-resolved W3 (4 of 37 cases on Kathleen, 4 of 20 on
locust); without it c_leaf 4 is neutral (0.995 and 1.000). c_leaf 8 beats 4 by 2.9 % (Kathleen)
and 1.3 % (locust), under the rule's 5 %, on the same workload alone. Every run was one rank on
all cores; c_leaf 4 at many ranks per node (where the time rule's levels and the cut depend on the
ranks per node, CH4.3b2c) is not measured here (CH4.5's line, an Also).

## 5. The test

`check_form_leaves` (`tetraleaf-exec/tests/helmholtz_fmm.rs`): Gaussian clusters, N = 2,000,
kD = 40, ε = 10⁻⁴, complex charges with gradients, `DEFAULT_C_LEAF`, 64 points per leaf, the W and
X rules off; the time rule at free form costs and at the machine's. On the M6 the form takes
levels 2 and 3 (k a 8.7 and 4.3, all 34 classes); level 3 holds 84 occupied leaves, into which M2P
and out of which P2L compute 7,660 point interactions each (`HelmholtzFmm::far_interactions`); the
errors 3.875·10⁻⁶ (φ) and 8.560·10⁻⁶ (∇φ), 1.000× the form off's on the same tree; 2 threads and
on 2 and 4 ranks the one-rank FMM bit for bit. The test asserts occupied leaves, M2P and P2L work
on a level of the form. (Level 2 has 380 X entries naming leaves without points: the empty-box rule
keeps X entries of source-free leaves, whose P2L adds nothing; not followed up.)

## 6. Questions for sign-off

1. **The acceptance as read**: c_leaf re-decided on W1–W3 under the default (the time rule),
   "leaves on high-frequency levels" measured as built (levels of k a ≥ 2.8, with the M2L each ran),
   §3.7's "W and X lists off them" checked rather than enforced (the census, the errors, the
   test). Recommended.
2. **c_leaf 4 the builder's default** (`DEFAULT_C_LEAF`, built in this PR; D11 amended a second
   time, `C_LEAF` = 1 staying the leaf rule's and the census's constant): 0.874× and 0.790× c_leaf 2's
   build plus one evaluation, all of the margin from the under-resolved W3 (0.26–0.36×), the rest
   within the spread or better (0.93–0.98× where the tree changes), errors within ε; with an Also
   on CH4.5 (c_leaf 4 at its best split and at N × 1). Recommended: of the two that qualify, 4,
   because 8's further 2.9 % and 1.3 % come from that one workload, below the rule's 5 %, at
   twice the phase range (352 rad) and leaf argument (x ≤ 43.5). Alternatives: (a) 8 (0.849× and
   0.780×); (b) keep 2 (this PR's default reverted).
3. **Proposal §3.7, §5.1 and §5.7 amended** to the finding, with D11's row: in §3.7 the first
   consequence (leaves may lie on levels of k a ≥ 2.8 whatever their M2L, the form's included, with
   W and X lists; the leaf operators do not depend on the regime and the degree rule covers their
   arguments, x ≤ 5.4 c_leaf) and the third (under-resolved inputs get leaves up to c_leaf λ, 4 λ
   by default, not more and sparser ones), and the default 4; in §5.1 the leaves' arguments; in §5.7
   the phase range (22 c_leaf, 44 c_leaf with V pairs near, 176 rad at c_leaf 4, within the sine and
   cosine's tested 2¹⁰). Recommended.

## 7. Reproduce

```sh
cargo build --release -p tetraleaf-validate --example helmholtz_accuracy
tools/helmholtz/ch4-3e-runs.sh <dir> w1 w2 w3                 # on Kathleen inside run.sbatch
ONLY='w1-1000|w3r' tools/helmholtz/ch4-3e-runs.sh <dir> large
python3 tools/helmholtz/ch4-3e-tables.py census|summary|cases|errors kathleen=<dir>[,<dir>…] locust=<dir>
```

On locust the same script with `w1 w3` (`/data/ucahtbe/runs/ch4-3e/drive.sh` sampled the load).
