//! S9's checks (Helmholtz CH0.9): the φ step's accuracy against the DFT as a sum, and a box's
//! bits across batch positions, batch sizes, threads and fresh processes, for FFTW
//! (interleaved and split, 2P − 1 and 2P rings) and rustfft. Run with `--release`.

use rayon::prelude::*;
use rlst::c64;
use tetraleaf_spike_fft_plans::{
    FftwPhi, FftwSplitPhi, Grid, Packed, RustPhi, Shape, hash, naive_analysis, naive_synthesis,
    random, relative_error, same_bits,
};

/// The bound of the accuracy check: the FFTs' round-off and the sum's, at N_φ ≤ 2048, are
/// both near 10⁻¹⁵.
const ACCURACY: f64 = 1e-13;

/// The back ends under test, by name.
const BACKENDS: [&str; 3] = ["fftw", "fftw-split", "rustfft"];

/// Synthesis then analysis of one box's modes by a back end: (samples of the live rings, the
/// analysis of the naive samples).
fn one_box(backend: &str, shape: Shape, modes: &[c64], samples_in: &[c64]) -> (Vec<c64>, Vec<c64>) {
    let live = shape.modes() * shape.azimuths;
    let mut out_modes = vec![c64::new(0.0, 0.0); shape.mode_len()];
    match backend {
        "fftw" => {
            let phi = FftwPhi::new(shape).unwrap();
            let mut scratch = phi.scratch();
            let mut boxes = Packed::<c64>::new(shape.box_len(), 1);
            phi.synthesize(modes, &mut scratch, boxes.get_mut(0));
            let samples = boxes.get(0)[..live].to_vec();
            boxes.get_mut(0)[..live].copy_from_slice(samples_in);
            phi.analyze(boxes.get(0), &mut scratch, &mut out_modes);
            (samples, out_modes)
        }
        "fftw-split" => {
            let phi = FftwSplitPhi::new(shape).unwrap();
            let mut scratch = phi.scratch();
            let mut boxes = Packed::<f64>::new(phi.box_reals(), 1);
            phi.synthesize(modes, &mut scratch, boxes.get_mut(0));
            let (re, im) = boxes.get(0).split_at(shape.box_len());
            let samples = (0..live).map(|i| c64::new(re[i], im[i])).collect();
            let b = boxes.get_mut(0);
            for (i, z) in samples_in.iter().enumerate() {
                b[i] = z.re;
                b[shape.box_len() + i] = z.im;
            }
            phi.analyze(boxes.get(0), &mut scratch, &mut out_modes);
            (samples, out_modes)
        }
        "rustfft" => {
            let phi = RustPhi::new(shape);
            let mut scratch = phi.scratch();
            let mut boxes = Packed::<c64>::new(shape.box_len(), 1);
            phi.synthesize(modes, &mut scratch, boxes.get_mut(0));
            let samples = boxes.get(0)[..live].to_vec();
            boxes.get_mut(0)[..live].copy_from_slice(samples_in);
            phi.analyze(boxes.get(0), &mut scratch, &mut out_modes);
            (samples, out_modes)
        }
        _ => unreachable!(),
    }
}

#[test]
fn accuracy_against_the_sum() {
    for p in [33, 82, 261] {
        for grid in Grid::ALL {
            for even in [false, true] {
                let shape = Shape::new(p, grid, even);
                let modes = random(shape.mode_len(), p as u64);
                let live = shape.modes() * shape.azimuths;
                let mut expected_samples = vec![c64::new(0.0, 0.0); live];
                naive_synthesis(&shape, &modes, &mut expected_samples);
                let samples_in = random(live, 1000 + p as u64);
                let mut expected_modes = vec![c64::new(0.0, 0.0); shape.mode_len()];
                naive_analysis(&shape, &samples_in, &mut expected_modes);
                for backend in BACKENDS {
                    let (samples, modes_out) = one_box(backend, shape, &modes, &samples_in);
                    let (es, ea) = (
                        relative_error(&samples, &expected_samples),
                        relative_error(&modes_out, &expected_modes),
                    );
                    println!(
                        "p {p} {} N_φ {} rings {} {backend}: synthesis {es:.2e}, analysis {ea:.2e}",
                        grid.name(),
                        shape.azimuths,
                        shape.rings
                    );
                    assert!(es < ACCURACY && ea < ACCURACY, "{backend} p {p} {grid:?}");
                }
            }
        }
    }
}

/// Every box of a batch of `boxes` holds the same input; the back end transforms each box
/// (`threads` workers, 0 for serial) and returns every box's output (synthesis samples, then
/// analysis modes), flattened per box.
fn batch_outputs(backend: &str, shape: Shape, boxes: usize, threads: usize) -> Vec<Vec<c64>> {
    let modes = random(shape.mode_len(), 11);
    let live = shape.modes() * shape.azimuths;
    let samples_in = random(live, 12);
    let run = |work: &(dyn Fn() -> Vec<Vec<c64>> + Sync)| {
        if threads == 0 {
            work()
        } else {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap()
                .install(work)
        }
    };
    match backend {
        "fftw" => {
            let phi = FftwPhi::new(shape).unwrap();
            run(&|| {
                let mut out = Packed::<c64>::new(shape.box_len(), boxes);
                let mut input = Packed::<c64>::new(shape.box_len(), boxes);
                for b in 0..boxes {
                    assert_eq!(out.box_class(b).offset(), 0);
                    input.get_mut(b)[..live].copy_from_slice(&samples_in);
                }
                let analysed: Vec<Vec<c64>> = out
                    .boxes_mut()
                    .zip(input.all().chunks_exact(shape.box_len()))
                    .par_bridge()
                    .map_init(
                        || phi.scratch(),
                        |scratch, (o, i)| {
                            phi.synthesize(&modes, scratch, o);
                            let mut m = vec![c64::new(0.0, 0.0); shape.mode_len()];
                            phi.analyze(i, scratch, &mut m);
                            let mut all = o.to_vec();
                            all.extend(m);
                            all
                        },
                    )
                    .collect();
                analysed
            })
        }
        "fftw-split" => {
            let phi = FftwSplitPhi::new(shape).unwrap();
            run(&|| {
                let reals = phi.box_reals();
                let mut out = Packed::<f64>::new(reals, boxes);
                let mut input = Packed::<f64>::new(reals, boxes);
                for b in 0..boxes {
                    assert_eq!(out.box_class(b).offset(), 0);
                    let x = input.get_mut(b);
                    for (k, z) in samples_in.iter().enumerate() {
                        x[k] = z.re;
                        x[shape.box_len() + k] = z.im;
                    }
                }
                out.boxes_mut()
                    .zip(input.all().chunks_exact(reals))
                    .par_bridge()
                    .map_init(
                        || phi.scratch(),
                        |scratch, (o, i)| {
                            phi.synthesize(&modes, scratch, o);
                            let mut m = vec![c64::new(0.0, 0.0); shape.mode_len()];
                            phi.analyze(i, scratch, &mut m);
                            let mut all: Vec<c64> = o.iter().map(|&x| c64::new(x, 0.0)).collect();
                            all.extend(m);
                            all
                        },
                    )
                    .collect()
            })
        }
        "rustfft" => {
            let phi = RustPhi::new(shape);
            run(&|| {
                let mut out = Packed::<c64>::new(shape.box_len(), boxes);
                let mut input = Packed::<c64>::new(shape.box_len(), boxes);
                for b in 0..boxes {
                    input.get_mut(b)[..live].copy_from_slice(&samples_in);
                }
                out.boxes_mut()
                    .zip(input.all().chunks_exact(shape.box_len()))
                    .par_bridge()
                    .map_init(
                        || phi.scratch(),
                        |scratch, (o, i)| {
                            phi.synthesize(&modes, scratch, o);
                            let mut m = vec![c64::new(0.0, 0.0); shape.mode_len()];
                            phi.analyze(i, scratch, &mut m);
                            let mut all = o.to_vec();
                            all.extend(m);
                            all
                        },
                    )
                    .collect()
            })
        }
        _ => unreachable!(),
    }
}

#[test]
fn bits_across_positions_batches_and_threads() {
    let max = std::thread::available_parallelism().map_or(4, |n| n.get());
    for p in [33, 82] {
        for grid in [Grid::Product, Grid::Smooth, Grid::PowerOfTwo] {
            for even in [false, true] {
                let shape = Shape::new(p, grid, even);
                for backend in BACKENDS {
                    let reference = &batch_outputs(backend, shape, 1, 0)[0];
                    for boxes in [1, 2, 7, 33] {
                        for threads in [0, 1, 2, 4, max] {
                            let outputs = batch_outputs(backend, shape, boxes, threads);
                            assert_eq!(outputs.len(), boxes);
                            for (b, out) in outputs.iter().enumerate() {
                                assert!(
                                    same_bits(out, reference),
                                    "{backend} p {p} {grid:?} even {even}: box {b} of {boxes} \
                                     on {threads} threads differs"
                                );
                            }
                        }
                    }
                    println!(
                        "p {p} {} rings {} {backend}: same bits at 1, 2, 7, 33 boxes, \
                         serial and 1, 2, 4, {max} threads",
                        grid.name(),
                        shape.rings
                    );
                }
            }
        }
    }
}

/// The hashes of one fixed set of transforms, one line per back end and shape.
fn hashes() -> Vec<String> {
    let mut lines = Vec::new();
    for p in [33, 82] {
        for grid in Grid::ALL {
            for even in [false, true] {
                let shape = Shape::new(p, grid, even);
                for backend in BACKENDS {
                    let out = &batch_outputs(backend, shape, 1, 0)[0];
                    lines.push(format!(
                        "HASH {backend}-p{p}-{}-{} {:016x}",
                        grid.name(),
                        shape.rings,
                        hash(out)
                    ));
                }
            }
        }
    }
    lines
}

/// Prints [`hashes`] for [`bits_in_fresh_processes`], which runs it in child processes.
#[test]
#[ignore = "the child of bits_in_fresh_processes"]
fn child_hashes() {
    for line in hashes() {
        println!("{line}");
    }
}

#[test]
fn bits_in_fresh_processes() {
    let child = || {
        let output = std::process::Command::new(std::env::current_exe().unwrap())
            .args([
                "--ignored",
                "--exact",
                "child_hashes",
                "--nocapture",
                "--test-threads=1",
            ])
            .output()
            .unwrap();
        assert!(output.status.success());
        let lines: Vec<String> = String::from_utf8(output.stdout)
            .unwrap()
            .lines()
            .filter_map(|l| l.find("HASH ").map(|at| l[at..].to_string()))
            .collect();
        lines
    };
    let (first, second) = (child(), child());
    let here = hashes();
    assert_eq!(first.len(), here.len());
    assert_eq!(first, second, "two fresh processes differ");
    assert_eq!(first, here, "a fresh process differs from this one");
    println!(
        "{} hashes equal in two fresh processes and this one",
        here.len()
    );
}
