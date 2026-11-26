# S9 for the record: what changes a box's bits

Machine: host `locust.rc.ucl.ac.uk`, aarch64, features neon, sve, FFTW `fftw-3.3.11-neon`, rustfft 6.4 (Cargo.lock)

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
| 33 | product | 136 | split | 16 | 16 | 1.8e-16 |
| 33 | product | 136 | rustfft | 16 | 16 | 2.3e-16 |
| 33 | product | 136 | rustfft one complex off a line, against rustfft aligned | 16 | 0 | 0.0e0 |
| 33 | smooth | 140 | one plan over the batch | 1 | 0 | 0.0e0 |
| 33 | smooth | 140 | one plan over the batch | 2 | 0 | 0.0e0 |
| 33 | smooth | 140 | one plan over the batch | 7 | 0 | 0.0e0 |
| 33 | smooth | 140 | one plan over the batch | 16 | 0 | 0.0e0 |
| 33 | smooth | 140 | ring distance N_φ + 4 | 16 | 0 | 0.0e0 |
| 33 | smooth | 140 | ring distance N_φ + 1 | 16 | 0 | 0.0e0 |
| 33 | smooth | 140 | class 16 | 16 | 0 | 0.0e0 |
| 33 | smooth | 140 | 2P rings (live rings) | 16 | 0 | 0.0e0 |
| 33 | smooth | 140 | split | 16 | 16 | 2.2e-16 |
| 33 | smooth | 140 | rustfft | 16 | 16 | 2.5e-16 |
| 33 | smooth | 140 | rustfft one complex off a line, against rustfft aligned | 16 | 0 | 0.0e0 |
| 33 | pow2 | 256 | one plan over the batch | 1 | 0 | 0.0e0 |
| 33 | pow2 | 256 | one plan over the batch | 2 | 0 | 0.0e0 |
| 33 | pow2 | 256 | one plan over the batch | 7 | 0 | 0.0e0 |
| 33 | pow2 | 256 | one plan over the batch | 16 | 0 | 0.0e0 |
| 33 | pow2 | 256 | ring distance N_φ + 4 | 16 | 0 | 0.0e0 |
| 33 | pow2 | 256 | ring distance N_φ + 1 | 16 | 0 | 0.0e0 |
| 33 | pow2 | 256 | class 16 | 16 | 0 | 0.0e0 |
| 33 | pow2 | 256 | 2P rings (live rings) | 16 | 0 | 0.0e0 |
| 33 | pow2 | 256 | split | 16 | 16 | 1.9e-16 |
| 33 | pow2 | 256 | rustfft | 16 | 16 | 2.9e-16 |
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
| 82 | product | 332 | rustfft | 16 | 16 | 7.8e-16 |
| 82 | product | 332 | rustfft one complex off a line, against rustfft aligned | 16 | 0 | 0.0e0 |
| 82 | smooth | 336 | one plan over the batch | 1 | 0 | 0.0e0 |
| 82 | smooth | 336 | one plan over the batch | 2 | 0 | 0.0e0 |
| 82 | smooth | 336 | one plan over the batch | 7 | 0 | 0.0e0 |
| 82 | smooth | 336 | one plan over the batch | 16 | 0 | 0.0e0 |
| 82 | smooth | 336 | ring distance N_φ + 4 | 16 | 0 | 0.0e0 |
| 82 | smooth | 336 | ring distance N_φ + 1 | 16 | 0 | 0.0e0 |
| 82 | smooth | 336 | class 16 | 16 | 0 | 0.0e0 |
| 82 | smooth | 336 | 2P rings (live rings) | 16 | 0 | 0.0e0 |
| 82 | smooth | 336 | split | 16 | 16 | 2.4e-16 |
| 82 | smooth | 336 | rustfft | 16 | 16 | 3.0e-16 |
| 82 | smooth | 336 | rustfft one complex off a line, against rustfft aligned | 16 | 0 | 0.0e0 |
| 82 | pow2 | 512 | one plan over the batch | 1 | 0 | 0.0e0 |
| 82 | pow2 | 512 | one plan over the batch | 2 | 0 | 0.0e0 |
| 82 | pow2 | 512 | one plan over the batch | 7 | 0 | 0.0e0 |
| 82 | pow2 | 512 | one plan over the batch | 16 | 0 | 0.0e0 |
| 82 | pow2 | 512 | ring distance N_φ + 4 | 16 | 0 | 0.0e0 |
| 82 | pow2 | 512 | ring distance N_φ + 1 | 16 | 0 | 0.0e0 |
| 82 | pow2 | 512 | class 16 | 16 | 0 | 0.0e0 |
| 82 | pow2 | 512 | 2P rings (live rings) | 16 | 0 | 0.0e0 |
| 82 | pow2 | 512 | split | 16 | 16 | 2.1e-16 |
| 82 | pow2 | 512 | rustfft | 16 | 16 | 3.2e-16 |
| 82 | pow2 | 512 | rustfft one complex off a line, against rustfft aligned | 16 | 0 | 0.0e0 |
| 144 | product | 580 | one plan over the batch | 1 | 0 | 0.0e0 |
| 144 | product | 580 | one plan over the batch | 2 | 0 | 0.0e0 |
| 144 | product | 580 | one plan over the batch | 7 | 0 | 0.0e0 |
| 144 | product | 580 | one plan over the batch | 16 | 0 | 0.0e0 |
| 144 | product | 580 | ring distance N_φ + 4 | 16 | 0 | 0.0e0 |
| 144 | product | 580 | ring distance N_φ + 1 | 16 | 0 | 0.0e0 |
| 144 | product | 580 | class 16 | 16 | 0 | 0.0e0 |
| 144 | product | 580 | 2P rings (live rings) | 16 | 0 | 0.0e0 |
| 144 | product | 580 | split | 16 | 16 | 1.9e-16 |
| 144 | product | 580 | rustfft | 16 | 16 | 2.9e-16 |
| 144 | product | 580 | rustfft one complex off a line, against rustfft aligned | 16 | 0 | 0.0e0 |
| 144 | smooth | 588 | one plan over the batch | 1 | 0 | 0.0e0 |
| 144 | smooth | 588 | one plan over the batch | 2 | 0 | 0.0e0 |
| 144 | smooth | 588 | one plan over the batch | 7 | 0 | 0.0e0 |
| 144 | smooth | 588 | one plan over the batch | 16 | 0 | 0.0e0 |
| 144 | smooth | 588 | ring distance N_φ + 4 | 16 | 0 | 0.0e0 |
| 144 | smooth | 588 | ring distance N_φ + 1 | 16 | 0 | 0.0e0 |
| 144 | smooth | 588 | class 16 | 16 | 0 | 0.0e0 |
| 144 | smooth | 588 | 2P rings (live rings) | 16 | 0 | 0.0e0 |
| 144 | smooth | 588 | split | 16 | 16 | 2.8e-16 |
| 144 | smooth | 588 | rustfft | 16 | 16 | 4.2e-16 |
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
| 144 | pow2 | 1024 | rustfft | 16 | 16 | 3.5e-16 |
| 144 | pow2 | 1024 | rustfft one complex off a line, against rustfft aligned | 16 | 0 | 0.0e0 |

## FFTW's plans at p = 82 (synthesis)

- product N_φ 332 rings 165: fingerprint 546336ce8788ea4b

  ```text
  (dft-vrank>=1-x165/1
    (dft-ct-dit/4
      (dftw-direct-4/12 "t2bv_4_neon")
      (dft-vrank>=1-x4/1
        (dft-bluestein-83/180
          (dft-ct-dit/10
            (dftw-direct-10/6 "t3fv_10_neon")
            (dft-buffered-18-x10/10-4
              (dft-vrank>=1-x10/1
                (dft-ct-dit/3
                  (dftw-direct-3/4 "t1fv_3_neon")
                  (dft-direct-6-x3 "n2fv_6_neon")))
              (dft-r2hc-1
                (rdft-rank0-iter-ci/36-x10))
              (dft-nop)))))))
  ```

- product N_φ 332 rings 166: fingerprint 4672a4490cdeb43e

  ```text
  (dft-vrank>=1-x166/1
    (dft-ct-dit/4
      (dftw-direct-4/12 "t2bv_4_neon")
      (dft-vrank>=1-x4/1
        (dft-bluestein-83/180
          (dft-ct-dit/10
            (dftw-direct-10/6 "t3fv_10_neon")
            (dft-buffered-18-x10/10-4
              (dft-vrank>=1-x10/1
                (dft-ct-dit/3
                  (dftw-direct-3/4 "t1fv_3_neon")
                  (dft-direct-6-x3 "n2fv_6_neon")))
              (dft-r2hc-1
                (rdft-rank0-iter-ci/36-x10))
              (dft-nop)))))))
  ```

- smooth N_φ 336 rings 165: fingerprint 82c894394c4e0611

  ```text
  (dft-vrank>=1-x165/1
    (dft-ct-dit/16
      (dftw-direct-16/8 "t3bv_16_neon")
      (dft-vrank>=1-x16/1
        (dft-ct-dit/3
          (dftw-direct-3/4 "t1bv_3_neon")
          (dft-direct-7-x3 "n1bv_7_neon")))))
  ```

- smooth N_φ 336 rings 166: fingerprint 1166bb10098285c2

  ```text
  (dft-vrank>=1-x166/1
    (dft-ct-dit/16
      (dftw-direct-16/8 "t3bv_16_neon")
      (dft-vrank>=1-x16/1
        (dft-ct-dit/3
          (dftw-direct-3/4 "t1bv_3_neon")
          (dft-direct-7-x3 "n1bv_7_neon")))))
  ```

- pow2 N_φ 512 rings 165: fingerprint c5a66a57dd841368

  ```text
  (dft-vrank>=1-x165/1
    (dft-ct-dit/16
      (dftw-direct-16/8 "t3bv_16_neon")
      (dft-direct-32-x16 "n2bv_32_neon")))
  ```

- pow2 N_φ 512 rings 166: fingerprint b44c17a870d2b739

  ```text
  (dft-vrank>=1-x166/1
    (dft-ct-dit/16
      (dftw-direct-16/8 "t3bv_16_neon")
      (dft-direct-32-x16 "n2bv_32_neon")))
  ```

