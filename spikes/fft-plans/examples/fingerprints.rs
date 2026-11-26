//! S9's fingerprints across ranks (Helmholtz CH0.9): every rank plans the φ step of every
//! degree, grid, ring count and storage (`Rigor::Estimate`, class 0, both directions), and the
//! ranks all-gather the plans' fingerprints (`Plan::fingerprint`, a hash of FFTW's
//! description of the plan). Every rank compares the same gathered values, so the check
//! fails on every rank or on none.
//!
//!     mpirun -n 2 target/spikes-fft-plans/release/examples/fingerprints

use mpi::traits::*;
use tetraleaf_spike_fft_plans::{DEGREES, FftwPhi, FftwSplitPhi, Grid, Shape, machine_line};

fn main() {
    let universe = mpi::initialize().unwrap();
    let world = universe.world();
    let (rank, size) = (world.rank(), world.size() as usize);
    let mut names = Vec::new();
    let mut prints: Vec<u64> = Vec::new();
    for p in DEGREES {
        for grid in Grid::ALL {
            for even in [false, true] {
                let shape = Shape::new(p, grid, even);
                let base = format!(
                    "p {p} {} N_φ {} rings {}",
                    grid.name(),
                    shape.azimuths,
                    shape.rings
                );
                for (kind, f) in [
                    ("interleaved", FftwPhi::new(shape).unwrap().fingerprints()),
                    ("split", FftwSplitPhi::new(shape).unwrap().fingerprints()),
                ] {
                    for (direction, value) in ["synthesis", "analysis"].iter().zip(f) {
                        names.push(format!("{base} {kind} {direction}"));
                        prints.push(value);
                    }
                }
            }
        }
    }
    let mut all = vec![0u64; prints.len() * size];
    world.all_gather_into(&prints[..], &mut all[..]);
    let machines = {
        let line = machine_line();
        let mut bytes = line.into_bytes();
        bytes.resize(256, b' ');
        let mut gathered = vec![0u8; 256 * size];
        world.all_gather_into(&bytes[..], &mut gathered[..]);
        gathered
            .as_chunks::<256>()
            .0
            .iter()
            .map(|c| String::from_utf8_lossy(c).trim_end().to_string())
            .collect::<Vec<_>>()
    };
    let differing: Vec<&String> = names
        .iter()
        .enumerate()
        .filter(|&(i, _)| (1..size).any(|r| all[r * prints.len() + i] != all[i]))
        .map(|(_, n)| n)
        .collect();
    if rank == 0 {
        println!("## fingerprints across {size} ranks\n");
        for (r, m) in machines.iter().enumerate() {
            println!("- rank {r}: {m}");
        }
        println!(
            "\n{} plans per rank (12 degrees × 3 grids × 2 ring counts × interleaved and split × 2 directions); {} with fingerprints that differ between ranks.",
            names.len(),
            differing.len()
        );
        let distinct: std::collections::BTreeSet<u64> = prints.iter().copied().collect();
        println!("{} distinct fingerprints on rank 0.", distinct.len());
        for n in &differing {
            println!("- differs: {n}");
        }
    }
    assert!(differing.is_empty(), "fingerprints differ across ranks");
}
