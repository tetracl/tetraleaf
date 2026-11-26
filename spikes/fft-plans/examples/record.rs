//! S9 "for the record" (Helmholtz CH0.9): which plan choices change a box's bits. Against the
//! per-box FFTW plan of 2P − 1 rings on packed class-0 buffers (the layout of proposal §6.7),
//! the synthesis of boxes with distinct modes by
//!
//! - one plan over the whole batch (howmany = boxes · rings), at several batch sizes;
//! - a per-box plan over rings `extra` elements further apart (another ring stride);
//! - a per-box plan for alignment class 16 on buffers one complex off a line;
//! - the per-box plan of 2P rings (one zero ring), on the live rings;
//! - the split plan; rustfft; rustfft on buffers one complex off a line.
//!
//! Prints Markdown: per variant the boxes whose bits differ and the largest relative
//! difference, then FFTW's plan descriptions at one degree.
//!
//!     cargo run --release --example record

use rlst::c64;
use rlst::fftw::{AlignmentClass, Planner};
use tetraleaf_spike_fft_plans::{
    FftwPhi, FftwSplitPhi, Grid, Packed, RustPhi, Shape, batch_layout, machine_line, random,
    relative_error, same_bits,
};

/// Each box's live samples from the per-box plan, packed, class 0.
fn per_box(shape: Shape, modes: &[Vec<c64>]) -> Vec<Vec<c64>> {
    let phi = FftwPhi::new(shape).unwrap();
    let mut scratch = phi.scratch();
    let mut out = Packed::<c64>::new(shape.box_len(), modes.len());
    let live = shape.modes() * shape.azimuths;
    (0..modes.len())
        .map(|b| {
            phi.synthesize(&modes[b], &mut scratch, out.get_mut(b));
            out.get(b)[..live].to_vec()
        })
        .collect()
}

/// One plan over every box of the batch.
fn whole_batch(shape: Shape, modes: &[Vec<c64>]) -> Vec<Vec<c64>> {
    let boxes = modes.len();
    let layout = batch_layout(&shape, boxes).unwrap();
    let plan = Planner::<f64>::new()
        .c2c(&layout, rlst::fftw::Direction::Backward)
        .unwrap();
    let mut input = Packed::<c64>::new(shape.box_len(), boxes);
    let mut out = Packed::<c64>::new(shape.box_len(), boxes);
    for (b, m) in modes.iter().enumerate() {
        shape.scatter(m, input.get_mut(b));
    }
    plan.execute(input.all(), out.all_mut()).unwrap();
    let live = shape.modes() * shape.azimuths;
    (0..boxes).map(|b| out.get(b)[..live].to_vec()).collect()
}

/// A per-box plan whose rings lie N_φ + `extra` apart.
fn strided(shape: Shape, modes: &[Vec<c64>], extra: usize) -> Vec<Vec<c64>> {
    let distance = shape.azimuths + extra;
    let phi = FftwPhi::with_planner(shape, &Planner::<f64>::new(), distance).unwrap();
    let len = distance * shape.rings;
    let (mut input, mut out) = (Packed::<c64>::new(len, 1), Packed::<c64>::new(len, 1));
    modes
        .iter()
        .map(|m| {
            let x = input.get_mut(0);
            x.fill(c64::new(0.0, 0.0));
            for (ring, row) in m.chunks_exact(shape.modes()).enumerate() {
                let mut line = vec![c64::new(0.0, 0.0); shape.azimuths];
                let single = Shape { rings: 1, ..shape };
                single.scatter(row, &mut line);
                x[ring * distance..ring * distance + shape.azimuths].copy_from_slice(&line);
            }
            phi.synthesis.execute(input.get(0), out.get_mut(0)).unwrap();
            let o = out.get(0);
            (0..shape.modes())
                .flat_map(|ring| o[ring * distance..ring * distance + shape.azimuths].to_vec())
                .collect()
        })
        .collect()
}

/// A per-box plan for class 16, on buffers one complex off a line.
fn class16(shape: Shape, modes: &[Vec<c64>]) -> Vec<Vec<c64>> {
    let class = AlignmentClass::from_offset(16);
    let planner = Planner::<f64>::new().alignment(class, class);
    let phi = FftwPhi::with_planner(shape, &planner, shape.azimuths).unwrap();
    let len = shape.box_len();
    let (mut input, mut out) = (
        Packed::<c64>::new(len + 4, 1),
        Packed::<c64>::new(len + 4, 1),
    );
    let live = shape.modes() * shape.azimuths;
    modes
        .iter()
        .map(|m| {
            let x = &mut input.get_mut(0)[1..=len];
            x.fill(c64::new(0.0, 0.0));
            shape.scatter(m, x);
            assert_eq!(AlignmentClass::of_ptr(x.as_ptr()).offset(), 16);
            phi.synthesis
                .execute(&input.get(0)[1..=len], &mut out.get_mut(0)[1..=len])
                .unwrap();
            out.get(0)[1..=live].to_vec()
        })
        .collect()
}

/// The split plan, read back as complex values.
fn split(shape: Shape, modes: &[Vec<c64>]) -> Vec<Vec<c64>> {
    let phi = FftwSplitPhi::new(shape).unwrap();
    let mut scratch = phi.scratch();
    let mut out = Packed::<f64>::new(phi.box_reals(), 1);
    let live = shape.modes() * shape.azimuths;
    modes
        .iter()
        .map(|m| {
            phi.synthesize(m, &mut scratch, out.get_mut(0));
            let (re, im) = out.get(0).split_at(shape.box_len());
            (0..live).map(|i| c64::new(re[i], im[i])).collect()
        })
        .collect()
}

/// rustfft, on buffers at `offset` complex values from a line.
fn rust(shape: Shape, modes: &[Vec<c64>], offset: usize) -> Vec<Vec<c64>> {
    let phi = RustPhi::new(shape);
    let mut scratch = phi.scratch();
    let len = shape.box_len();
    let mut out = Packed::<c64>::new(len + 4, 1);
    let live = shape.modes() * shape.azimuths;
    modes
        .iter()
        .map(|m| {
            let o = &mut out.get_mut(0)[offset..offset + len];
            phi.synthesize(m, &mut scratch, o);
            o[..live].to_vec()
        })
        .collect()
}

/// (boxes whose bits differ, the largest relative ℓ2 difference).
fn compare(a: &[Vec<c64>], b: &[Vec<c64>]) -> (usize, f64) {
    a.iter().zip(b).fold((0, 0.0_f64), |(n, e), (x, y)| {
        (
            n + usize::from(!same_bits(x, y)),
            e.max(relative_error(x, y)),
        )
    })
}

fn main() {
    println!("# S9 for the record: what changes a box's bits\n");
    println!("Machine: {}\n", machine_line());
    println!(
        "Against the per-box FFTW plan of 2P − 1 rings on packed class-0 buffers, synthesis \
         (backward) of boxes with distinct random modes.\n"
    );
    println!("| p | grid | N_φ | variant | boxes | boxes with other bits | max rel. diff. |");
    println!("| --- | --- | --- | --- | --- | --- | --- |");
    for p in [33, 82, 144] {
        for grid in Grid::ALL {
            let shape = Shape::new(p, grid, false);
            let modes: Vec<Vec<c64>> = (0..16).map(|b| random(shape.mode_len(), 100 + b)).collect();
            let reference = per_box(shape, &modes);
            let row = |variant: &str, boxes: usize, other: Vec<Vec<c64>>| {
                let (n, e) = compare(&reference[..boxes], &other);
                println!(
                    "| {p} | {} | {} | {variant} | {boxes} | {n} | {e:.1e} |",
                    grid.name(),
                    shape.azimuths
                );
            };
            for boxes in [1, 2, 7, 16] {
                row(
                    "one plan over the batch",
                    boxes,
                    whole_batch(shape, &modes[..boxes]),
                );
            }
            row("ring distance N_φ + 4", 16, strided(shape, &modes, 4));
            row("ring distance N_φ + 1", 16, strided(shape, &modes, 1));
            row("class 16", 16, class16(shape, &modes));
            let even = Shape::new(p, grid, true);
            row("2P rings (live rings)", 16, per_box(even, &modes));
            row("split", 16, split(shape, &modes));
            row("rustfft", 16, rust(shape, &modes, 0));
            let aligned = rust(shape, &modes, 0);
            let (n, e) = compare(&aligned, &rust(shape, &modes, 1));
            println!(
                "| {p} | {} | {} | rustfft one complex off a line, against rustfft aligned | 16 | {n} | {e:.1e} |",
                grid.name(),
                shape.azimuths
            );
        }
    }
    println!("\n## FFTW's plans at p = 82 (synthesis)\n");
    for grid in Grid::ALL {
        for even in [false, true] {
            let shape = Shape::new(82, grid, even);
            let phi = FftwPhi::new(shape).unwrap();
            println!(
                "- {} N_φ {} rings {}: fingerprint {:016x}\n\n  ```text\n  {}\n  ```\n",
                grid.name(),
                shape.azimuths,
                shape.rings,
                phi.synthesis.fingerprint(),
                phi.synthesis.description().replace('\n', "\n  ")
            );
        }
    }
}
