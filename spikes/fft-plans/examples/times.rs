//! S9's time per box (Helmholtz CH0.9): one back end of the φ step on one grid per launch (on
//! Kathleen a second setting in one launch runs 1.09–1.20× slower; review 2026-11-24, B12),
//! one thread, for P = p + 1 from 34 to 262.
//!
//!     cargo run --release --example times -- <variant> <grid> [p,p,...]
//!
//! `<variant>`: `product` (the shipped `Conversions`, whole S → F and F → R, on the product
//! grid only, with the φ step's share of their multiply-adds), `fftw` (per-box plan of 2P − 1
//! rings, interleaved), `fftw-even` (2P rings), `fftw-split` (split storage, 2P − 1 rings),
//! `rustfft`. `<grid>`: `product`, `smooth`, `pow2`. Each FFT time includes the scatter or
//! gather of the modes. Prints one Markdown table; times are the median of 7 samples of at
//! least 0.05 s each, per box, in µs.

use rlst::c64;
use tetraleaf_spike_fft_plans::{
    DEGREES, FftwPhi, FftwSplitPhi, Grid, Packed, ProductBoxes, RustPhi, Shape, machine_line,
    random, time_per_call, timing_boxes,
};

const SAMPLES: usize = 7;
const MIN_SECONDS: f64 = 0.05;

/// (S → F µs, F → R µs) per box of an FFT back end.
fn fft_times(variant: &str, shape: Shape) -> (usize, f64, f64) {
    let boxes = timing_boxes(shape.box_len());
    let modes: Vec<Vec<c64>> = (0..boxes)
        .map(|b| random(shape.mode_len(), b as u64))
        .collect();
    let mut out_modes = vec![c64::new(0.0, 0.0); shape.mode_len()];
    let per_box = |t: f64| t / boxes as f64 * 1e6;
    let (s2f, f2r) = match variant {
        "fftw" | "fftw-even" => {
            let phi = FftwPhi::new(shape).unwrap();
            let mut scratch = phi.scratch();
            let mut samples = Packed::<c64>::new(shape.box_len(), boxes);
            let s2f = time_per_call(SAMPLES, MIN_SECONDS, || {
                for (b, m) in modes.iter().enumerate() {
                    phi.synthesize(m, &mut scratch, samples.get_mut(b));
                }
            });
            let f2r = time_per_call(SAMPLES, MIN_SECONDS, || {
                for b in 0..boxes {
                    phi.analyze(samples.get(b), &mut scratch, &mut out_modes);
                }
            });
            (s2f, f2r)
        }
        "fftw-split" => {
            let phi = FftwSplitPhi::new(shape).unwrap();
            let mut scratch = phi.scratch();
            let mut samples = Packed::<f64>::new(phi.box_reals(), boxes);
            let s2f = time_per_call(SAMPLES, MIN_SECONDS, || {
                for (b, m) in modes.iter().enumerate() {
                    phi.synthesize(m, &mut scratch, samples.get_mut(b));
                }
            });
            let f2r = time_per_call(SAMPLES, MIN_SECONDS, || {
                for b in 0..boxes {
                    phi.analyze(samples.get(b), &mut scratch, &mut out_modes);
                }
            });
            (s2f, f2r)
        }
        "rustfft" => {
            let phi = RustPhi::new(shape);
            let mut scratch = phi.scratch();
            let mut samples = Packed::<c64>::new(shape.box_len(), boxes);
            let s2f = time_per_call(SAMPLES, MIN_SECONDS, || {
                for (b, m) in modes.iter().enumerate() {
                    phi.synthesize(m, &mut scratch, samples.get_mut(b));
                }
            });
            let f2r = time_per_call(SAMPLES, MIN_SECONDS, || {
                for b in 0..boxes {
                    phi.analyze(samples.get(b), &mut scratch, &mut out_modes);
                }
            });
            (s2f, f2r)
        }
        _ => panic!("unknown variant {variant}"),
    };
    (boxes, per_box(s2f), per_box(f2r))
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let usage =
        "usage: times <product|fftw|fftw-even|fftw-split|rustfft> <product|smooth|pow2> [p,p,...]";
    let variant = args.first().expect(usage).as_str();
    let grid = Grid::parse(args.get(1).expect(usage)).expect(usage);
    let degrees: Vec<usize> = args.get(2).map_or(DEGREES.to_vec(), |s| {
        s.split(',').map(|p| p.parse().expect(usage)).collect()
    });
    println!("## times {variant} {}\n", grid.name());
    println!("Machine: {}\n", machine_line());
    if variant == "product" {
        assert_eq!(grid, Grid::Product, "the product runs on its own grid");
        println!("| p | N_φ | boxes | S → F µs | F → R µs | φ share | φ S → F µs | φ F → R µs |");
        println!("| --- | --- | --- | --- | --- | --- | --- | --- |");
        for &p in &degrees {
            let shape = Shape::new(p, grid, false);
            let boxes = timing_boxes(shape.modes() * shape.azimuths);
            let mut product = ProductBoxes::new(p, boxes);
            let s2f = time_per_call(SAMPLES, MIN_SECONDS, || product.s2f()) / boxes as f64 * 1e6;
            let f2r = time_per_call(SAMPLES, MIN_SECONDS, || product.f2r()) / boxes as f64 * 1e6;
            let share = product.phi_share();
            println!(
                "| {p} | {} | {boxes} | {s2f:.1} | {f2r:.1} | {share:.3} | {:.1} | {:.1} |",
                shape.azimuths,
                s2f * share,
                f2r * share
            );
        }
    } else {
        println!("| p | N_φ | rings | boxes | S → F µs | F → R µs |");
        println!("| --- | --- | --- | --- | --- | --- |");
        for &p in &degrees {
            let shape = Shape::new(p, grid, variant == "fftw-even");
            let (boxes, s2f, f2r) = fft_times(variant, shape);
            println!(
                "| {p} | {} | {} | {boxes} | {s2f:.1} | {f2r:.1} |",
                shape.azimuths, shape.rings
            );
        }
    }
}
