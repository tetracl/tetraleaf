# S9 for the record: what changes a box's bits

Machine: host `node-c11b-038.kathleen.ucl.ac.uk`, x86_64, features sse4.1, avx, fma, avx512f, FFTW `fftw-3.3.10-sse2-avx-avx2-avx2_128-avx512`, rustfft 6.4 (Cargo.lock)

Against the per-box FFTW plan of 2P − 1 rings on packed class-0 buffers, synthesis (backward) of boxes with distinct random modes.

| p | grid | N_φ | variant | boxes | boxes with other bits | max rel. diff. |
| --- | --- | --- | --- | --- | --- | --- |
| 33 | product | 136 | one plan over the batch | 1 | 0 | 0.0e0 |
| 33 | product | 136 | one plan over the batch | 2 | 0 | 0.0e0 |
| 33 | product | 136 | one plan over the batch | 7 | 0 | 0.0e0 |
| 33 | product | 136 | one plan over the batch | 16 | 0 | 0.0e0 |
| 33 | product | 136 | ring distance N_φ + 4 | 16 | 0 | 0.0e0 |
| 33 | product | 136 | ring distance N_φ + 1 | 16 | 0 | 0.0e0 |
| 33 | product | 136 | class 16 | 16 | 0 | 0.0e0 |
| 33 | product | 136 | 2P rings (live rings) | 16 | 0 | 0.0e0 |
| 33 | product | 136 | split | 16 | 16 | 1.4e-16 |
| 33 | product | 136 | rustfft | 16 | 16 | 2.6e-16 |
| 33 | product | 136 | rustfft one complex off a line, against rustfft aligned | 16 | 0 | 0.0e0 |
| 33 | smooth | 140 | one plan over the batch | 1 | 0 | 0.0e0 |
| 33 | smooth | 140 | one plan over the batch | 2 | 0 | 0.0e0 |
| 33 | smooth | 140 | one plan over the batch | 7 | 0 | 0.0e0 |
| 33 | smooth | 140 | one plan over the batch | 16 | 0 | 0.0e0 |
| 33 | smooth | 140 | ring distance N_φ + 4 | 16 | 0 | 0.0e0 |
| 33 | smooth | 140 | ring distance N_φ + 1 | 16 | 0 | 0.0e0 |
| 33 | smooth | 140 | class 16 | 16 | 0 | 0.0e0 |
| 33 | smooth | 140 | 2P rings (live rings) | 16 | 0 | 0.0e0 |
| 33 | smooth | 140 | split | 16 | 16 | 2.3e-16 |
| 33 | smooth | 140 | rustfft | 16 | 16 | 2.9e-16 |
| 33 | smooth | 140 | rustfft one complex off a line, against rustfft aligned | 16 | 0 | 0.0e0 |
| 33 | pow2 | 256 | one plan over the batch | 1 | 0 | 0.0e0 |
| 33 | pow2 | 256 | one plan over the batch | 2 | 0 | 0.0e0 |
| 33 | pow2 | 256 | one plan over the batch | 7 | 0 | 0.0e0 |
| 33 | pow2 | 256 | one plan over the batch | 16 | 0 | 0.0e0 |
| 33 | pow2 | 256 | ring distance N_φ + 4 | 16 | 0 | 0.0e0 |
| 33 | pow2 | 256 | ring distance N_φ + 1 | 16 | 0 | 0.0e0 |
| 33 | pow2 | 256 | class 16 | 16 | 0 | 0.0e0 |
| 33 | pow2 | 256 | 2P rings (live rings) | 16 | 0 | 0.0e0 |
| 33 | pow2 | 256 | split | 16 | 16 | 2.0e-16 |
| 33 | pow2 | 256 | rustfft | 16 | 16 | 2.7e-16 |
| 33 | pow2 | 256 | rustfft one complex off a line, against rustfft aligned | 16 | 0 | 0.0e0 |
| 82 | product | 332 | one plan over the batch | 1 | 0 | 0.0e0 |
| 82 | product | 332 | one plan over the batch | 2 | 0 | 0.0e0 |
| 82 | product | 332 | one plan over the batch | 7 | 0 | 0.0e0 |
| 82 | product | 332 | one plan over the batch | 16 | 0 | 0.0e0 |
| 82 | product | 332 | ring distance N_φ + 4 | 16 | 0 | 0.0e0 |
| 82 | product | 332 | ring distance N_φ + 1 | 16 | 0 | 0.0e0 |
| 82 | product | 332 | class 16 | 16 | 0 | 0.0e0 |
| 82 | product | 332 | 2P rings (live rings) | 16 | 0 | 0.0e0 |
| 82 | product | 332 | split | 16 | 16 | 1.1e-16 |
| 82 | product | 332 | rustfft | 16 | 16 | 7.6e-16 |
| 82 | product | 332 | rustfft one complex off a line, against rustfft aligned | 16 | 0 | 0.0e0 |
| 82 | smooth | 336 | one plan over the batch | 1 | 0 | 0.0e0 |
| 82 | smooth | 336 | one plan over the batch | 2 | 0 | 0.0e0 |
| 82 | smooth | 336 | one plan over the batch | 7 | 0 | 0.0e0 |
| 82 | smooth | 336 | one plan over the batch | 16 | 0 | 0.0e0 |
| 82 | smooth | 336 | ring distance N_φ + 4 | 16 | 0 | 0.0e0 |
| 82 | smooth | 336 | ring distance N_φ + 1 | 16 | 0 | 0.0e0 |
| 82 | smooth | 336 | class 16 | 16 | 0 | 0.0e0 |
| 82 | smooth | 336 | 2P rings (live rings) | 16 | 0 | 0.0e0 |
| 82 | smooth | 336 | split | 16 | 16 | 2.5e-16 |
| 82 | smooth | 336 | rustfft | 16 | 16 | 3.9e-16 |
| 82 | smooth | 336 | rustfft one complex off a line, against rustfft aligned | 16 | 0 | 0.0e0 |
| 82 | pow2 | 512 | one plan over the batch | 1 | 0 | 0.0e0 |
| 82 | pow2 | 512 | one plan over the batch | 2 | 0 | 0.0e0 |
| 82 | pow2 | 512 | one plan over the batch | 7 | 0 | 0.0e0 |
| 82 | pow2 | 512 | one plan over the batch | 16 | 0 | 0.0e0 |
| 82 | pow2 | 512 | ring distance N_φ + 4 | 16 | 0 | 0.0e0 |
| 82 | pow2 | 512 | ring distance N_φ + 1 | 16 | 0 | 0.0e0 |
| 82 | pow2 | 512 | class 16 | 16 | 0 | 0.0e0 |
| 82 | pow2 | 512 | 2P rings (live rings) | 16 | 0 | 0.0e0 |
| 82 | pow2 | 512 | split | 16 | 16 | 2.2e-16 |
| 82 | pow2 | 512 | rustfft | 16 | 16 | 2.8e-16 |
| 82 | pow2 | 512 | rustfft one complex off a line, against rustfft aligned | 16 | 0 | 0.0e0 |
| 144 | product | 580 | one plan over the batch | 1 | 0 | 0.0e0 |
| 144 | product | 580 | one plan over the batch | 2 | 0 | 0.0e0 |
| 144 | product | 580 | one plan over the batch | 7 | 0 | 0.0e0 |
| 144 | product | 580 | one plan over the batch | 16 | 0 | 0.0e0 |
| 144 | product | 580 | ring distance N_φ + 4 | 16 | 0 | 0.0e0 |
| 144 | product | 580 | ring distance N_φ + 1 | 16 | 0 | 0.0e0 |
| 144 | product | 580 | class 16 | 16 | 0 | 0.0e0 |
| 144 | product | 580 | 2P rings (live rings) | 16 | 0 | 0.0e0 |
| 144 | product | 580 | split | 16 | 16 | 2.0e-16 |
| 144 | product | 580 | rustfft | 16 | 16 | 3.4e-16 |
| 144 | product | 580 | rustfft one complex off a line, against rustfft aligned | 16 | 0 | 0.0e0 |
| 144 | smooth | 588 | one plan over the batch | 1 | 0 | 0.0e0 |
| 144 | smooth | 588 | one plan over the batch | 2 | 0 | 0.0e0 |
| 144 | smooth | 588 | one plan over the batch | 7 | 0 | 0.0e0 |
| 144 | smooth | 588 | one plan over the batch | 16 | 0 | 0.0e0 |
| 144 | smooth | 588 | ring distance N_φ + 4 | 16 | 0 | 0.0e0 |
| 144 | smooth | 588 | ring distance N_φ + 1 | 16 | 0 | 0.0e0 |
| 144 | smooth | 588 | class 16 | 16 | 0 | 0.0e0 |
| 144 | smooth | 588 | 2P rings (live rings) | 16 | 0 | 0.0e0 |
| 144 | smooth | 588 | split | 16 | 16 | 1.9e-16 |
| 144 | smooth | 588 | rustfft | 16 | 16 | 3.8e-16 |
| 144 | smooth | 588 | rustfft one complex off a line, against rustfft aligned | 16 | 0 | 0.0e0 |
| 144 | pow2 | 1024 | one plan over the batch | 1 | 0 | 0.0e0 |
| 144 | pow2 | 1024 | one plan over the batch | 2 | 0 | 0.0e0 |
| 144 | pow2 | 1024 | one plan over the batch | 7 | 0 | 0.0e0 |
| 144 | pow2 | 1024 | one plan over the batch | 16 | 0 | 0.0e0 |
| 144 | pow2 | 1024 | ring distance N_φ + 4 | 16 | 0 | 0.0e0 |
| 144 | pow2 | 1024 | ring distance N_φ + 1 | 16 | 0 | 0.0e0 |
| 144 | pow2 | 1024 | class 16 | 16 | 0 | 0.0e0 |
| 144 | pow2 | 1024 | 2P rings (live rings) | 16 | 0 | 0.0e0 |
| 144 | pow2 | 1024 | split | 16 | 16 | 2.3e-16 |
| 144 | pow2 | 1024 | rustfft | 16 | 16 | 3.2e-16 |
| 144 | pow2 | 1024 | rustfft one complex off a line, against rustfft aligned | 16 | 0 | 0.0e0 |

## FFTW's plans at p = 82 (synthesis)

- product N_φ 332 rings 165: fingerprint 447b5c5eef18204c

  ```text
  (dft-vrank>=1-x165/1
    (dft-ct-dit/4
      (dftw-direct-4/12 "t1buv_4_avx2")
      (dft-vrank>=1-x4/1
        (dft-bluestein-83/180
          (dft-ct-dit/5
            (dftw-direct-5/16 "t3fv_5_avx512")
            (dft-buffered-36-x5/5-2
              (dft-vrank>=1-x5/1
                (dft-ct-dit/9
                  (dftw-direct-9/64 "t1fv_9_avx512")
                  (dft-direct-4-x9 "n2fv_4_avx2_128")))
              (dft-r2hc-1
                (rdft-rank0-iter-ci/72-x5))
              (dft-nop)))))))
  ```

- product N_φ 332 rings 166: fingerprint 4803e798658719d5

  ```text
  (dft-vrank>=1-x166/1
    (dft-ct-dit/4
      (dftw-direct-4/12 "t1buv_4_avx2")
      (dft-vrank>=1-x4/1
        (dft-bluestein-83/180
          (dft-ct-dit/5
            (dftw-direct-5/16 "t3fv_5_avx512")
            (dft-buffered-36-x5/5-2
              (dft-vrank>=1-x5/1
                (dft-ct-dit/9
                  (dftw-direct-9/64 "t1fv_9_avx512")
                  (dft-direct-4-x9 "n2fv_4_avx2_128")))
              (dft-r2hc-1
                (rdft-rank0-iter-ci/72-x5))
              (dft-nop)))))))
  ```

- smooth N_φ 336 rings 165: fingerprint 0c48a1a2e8080268

  ```text
  (dft-vrank>=1-x165/1
    (dft-ct-dit/12
      (dftw-direct-12/88 "t1bv_12_avx512")
      (dft-vrank>=1-x12/1
        (dft-ct-dit/4
          (dftw-direct-4/12 "t1buv_4_avx2")
          (dft-direct-7-x4 "n1bv_7_avx512")))))
  ```

- smooth N_φ 336 rings 166: fingerprint 87f27eff966fb2e1

  ```text
  (dft-vrank>=1-x166/1
    (dft-ct-dit/12
      (dftw-direct-12/88 "t1bv_12_avx512")
      (dft-vrank>=1-x12/1
        (dft-ct-dit/4
          (dftw-direct-4/12 "t1buv_4_avx2")
          (dft-direct-7-x4 "n1bv_7_avx512")))))
  ```

- pow2 N_φ 512 rings 165: fingerprint 889d1664e0aba087

  ```text
  (dft-vrank>=1-x165/1
    (dft-ct-dit/16
      (dftw-direct-16/32 "t3bv_16_avx512")
      (dft-direct-32-x16 "n2bv_32_avx512")))
  ```

- pow2 N_φ 512 rings 166: fingerprint 5dca4f5cc3621850

  ```text
  (dft-vrank>=1-x166/1
    (dft-ct-dit/16
      (dftw-direct-16/32 "t3bv_16_avx512")
      (dft-direct-32-x16 "n2bv_32_avx512")))
  ```

