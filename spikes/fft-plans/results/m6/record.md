# S9 for the record: what changes a box's bits

Machine: host `Skynet`, aarch64, features neon, FFTW `fftw-3.3.11`, rustfft 6.4 (Cargo.lock)

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
| 33 | product | 136 | split | 16 | 0 | 0.0e0 |
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
| 33 | smooth | 140 | split | 16 | 0 | 0.0e0 |
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
| 33 | pow2 | 256 | split | 16 | 0 | 0.0e0 |
| 33 | pow2 | 256 | rustfft | 16 | 16 | 2.8e-16 |
| 33 | pow2 | 256 | rustfft one complex off a line, against rustfft aligned | 16 | 0 | 0.0e0 |
| 82 | product | 332 | one plan over the batch | 1 | 0 | 0.0e0 |
| 82 | product | 332 | one plan over the batch | 2 | 0 | 0.0e0 |
| 82 | product | 332 | one plan over the batch | 7 | 0 | 0.0e0 |
| 82 | product | 332 | one plan over the batch | 16 | 0 | 0.0e0 |
| 82 | product | 332 | ring distance N_φ + 4 | 16 | 0 | 0.0e0 |
| 82 | product | 332 | ring distance N_φ + 1 | 16 | 0 | 0.0e0 |
| 82 | product | 332 | class 16 | 16 | 0 | 0.0e0 |
| 82 | product | 332 | 2P rings (live rings) | 16 | 0 | 0.0e0 |
| 82 | product | 332 | split | 16 | 0 | 0.0e0 |
| 82 | product | 332 | rustfft | 16 | 16 | 7.1e-16 |
| 82 | product | 332 | rustfft one complex off a line, against rustfft aligned | 16 | 0 | 0.0e0 |
| 82 | smooth | 336 | one plan over the batch | 1 | 0 | 0.0e0 |
| 82 | smooth | 336 | one plan over the batch | 2 | 0 | 0.0e0 |
| 82 | smooth | 336 | one plan over the batch | 7 | 0 | 0.0e0 |
| 82 | smooth | 336 | one plan over the batch | 16 | 0 | 0.0e0 |
| 82 | smooth | 336 | ring distance N_φ + 4 | 16 | 0 | 0.0e0 |
| 82 | smooth | 336 | ring distance N_φ + 1 | 16 | 0 | 0.0e0 |
| 82 | smooth | 336 | class 16 | 16 | 0 | 0.0e0 |
| 82 | smooth | 336 | 2P rings (live rings) | 16 | 0 | 0.0e0 |
| 82 | smooth | 336 | split | 16 | 0 | 0.0e0 |
| 82 | smooth | 336 | rustfft | 16 | 16 | 2.9e-16 |
| 82 | smooth | 336 | rustfft one complex off a line, against rustfft aligned | 16 | 0 | 0.0e0 |
| 82 | pow2 | 512 | one plan over the batch | 1 | 0 | 0.0e0 |
| 82 | pow2 | 512 | one plan over the batch | 2 | 0 | 0.0e0 |
| 82 | pow2 | 512 | one plan over the batch | 7 | 0 | 0.0e0 |
| 82 | pow2 | 512 | one plan over the batch | 16 | 0 | 0.0e0 |
| 82 | pow2 | 512 | ring distance N_φ + 4 | 16 | 0 | 0.0e0 |
| 82 | pow2 | 512 | ring distance N_φ + 1 | 16 | 0 | 0.0e0 |
| 82 | pow2 | 512 | class 16 | 16 | 0 | 0.0e0 |
| 82 | pow2 | 512 | 2P rings (live rings) | 16 | 0 | 0.0e0 |
| 82 | pow2 | 512 | split | 16 | 0 | 0.0e0 |
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
| 144 | product | 580 | split | 16 | 0 | 0.0e0 |
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
| 144 | smooth | 588 | split | 16 | 0 | 0.0e0 |
| 144 | smooth | 588 | rustfft | 16 | 16 | 4.1e-16 |
| 144 | smooth | 588 | rustfft one complex off a line, against rustfft aligned | 16 | 0 | 0.0e0 |
| 144 | pow2 | 1024 | one plan over the batch | 1 | 0 | 0.0e0 |
| 144 | pow2 | 1024 | one plan over the batch | 2 | 0 | 0.0e0 |
| 144 | pow2 | 1024 | one plan over the batch | 7 | 0 | 0.0e0 |
| 144 | pow2 | 1024 | one plan over the batch | 16 | 0 | 0.0e0 |
| 144 | pow2 | 1024 | ring distance N_φ + 4 | 16 | 0 | 0.0e0 |
| 144 | pow2 | 1024 | ring distance N_φ + 1 | 16 | 0 | 0.0e0 |
| 144 | pow2 | 1024 | class 16 | 16 | 0 | 0.0e0 |
| 144 | pow2 | 1024 | 2P rings (live rings) | 16 | 0 | 0.0e0 |
| 144 | pow2 | 1024 | split | 16 | 0 | 0.0e0 |
| 144 | pow2 | 1024 | rustfft | 16 | 16 | 3.5e-16 |
| 144 | pow2 | 1024 | rustfft one complex off a line, against rustfft aligned | 16 | 0 | 0.0e0 |

## FFTW's plans at p = 82 (synthesis)

- product N_φ 332 rings 165: fingerprint 47cb83f6ae3344c2

  ```text
  (dft-vrank>=1-x165/1
    (dft-ct-dit/4
      (dftw-direct-4/6 "t1_4")
      (dft-vrank>=1-x4/1
        (dft-generic-83))))
  ```

- product N_φ 332 rings 166: fingerprint 566a54fe24e82edb

  ```text
  (dft-vrank>=1-x166/1
    (dft-ct-dit/4
      (dftw-direct-4/6 "t1_4")
      (dft-vrank>=1-x4/1
        (dft-generic-83))))
  ```

- smooth N_φ 336 rings 165: fingerprint ce4b2ae4c2745e25

  ```text
  (dft-vrank>=1-x165/1
    (dft-ct-dit/16
      (dftw-direct-16/8 "t2_16")
      (dft-vrank>=1-x16/1
        (dft-ct-dit/3
          (dftw-direct-3/4 "t1_3")
          (dft-direct-7-x3 "n1_7")))))
  ```

- smooth N_φ 336 rings 166: fingerprint 3d4f1b7a523ebbf0

  ```text
  (dft-vrank>=1-x166/1
    (dft-ct-dit/16
      (dftw-direct-16/8 "t2_16")
      (dft-vrank>=1-x16/1
        (dft-ct-dit/3
          (dftw-direct-3/4 "t1_3")
          (dft-direct-7-x3 "n1_7")))))
  ```

- pow2 N_φ 512 rings 165: fingerprint 3f1188feca47aec6

  ```text
  (dft-vrank>=1-x165/1
    (dft-ct-dit/16
      (dftw-direct-16/8 "t2_16")
      (dft-direct-32-x16 "n1_32")))
  ```

- pow2 N_φ 512 rings 166: fingerprint 408e4de69542a21f

  ```text
  (dft-vrank>=1-x166/1
    (dft-ct-dit/16
      (dftw-direct-16/8 "t2_16")
      (dft-direct-32-x16 "n1_32")))
  ```

