//! Where a Helmholtz leaf point's P2M and L2P time goes (Helmholtz CH4.9a1;
//! docs/helmholtz/tasks/CH4.9a1.md). Prints Markdown on stdout.
//!
//! `tetraleaf-exec`'s Helmholtz FMM runs P2M and L2P as one `tetraleaf_ref::helmholtz::leaf`
//! call per point, in the leaf's own frame: half-width 1, the points in [−1, 1)³, k̂ = k r
//! (`Densities::{p2m, l2p}`). Per point the reference forms the radial functions
//! (`tetraleaf_math::radial::bessel_j`, a backward recursion from a start order that depends
//! on the point, and with gradients their derivatives), the harmonics
//! (`tetraleaf_math::orthonormal::harmonics`, or `harmonics_derivatives` with gradients), the
//! damping factor e^{Im(k̂)ρ − β}, and then the sum over the (p + 1)² terms.
//!
//! For each degree of `--ps` and each k̂ of `--khat` (Re k̂; `--kim` Im k̂, default 0), on one
//! thread, over a batch of `--points` uniform random points (default 256) with complex charges:
//!
//! - **P2M**: the reference call over the batch (`leaf::p2m`), then its parts apart, each over
//!   the same points: the radial functions, the harmonics, the factor e (its `exp`), and the sum
//!   (the reference's own loop, `form`, restated: per degree the charge times e times the
//!   radial function, times each harmonic, added into the two planes). The parts are the same
//!   calls in the same order, so they reassemble the reference's multipole bit for bit; that is
//!   checked before any timing (the column "check").
//! - **L2P**: the reference call over the batch without and with gradients (`leaf::l2p`), and
//!   the radial functions without and with derivatives, the harmonics, and the harmonics with
//!   their derivatives (`harmonics_derivatives`); "rest" is the call with gradients less those
//!   two parts (the sum, the gradient's assembly, the call's overhead).
//!
//! Every figure is µs per point, the median over batches (`tetraleaf_validate::bench::
//! median_time_per_call`); "ns per term" divides by (p + 1)². One precision per launch
//! (`--precision f64|f32`, default f64). Run in release mode on one core of a quiet machine:
//!
//! ```text
//! cargo run --release -p tetraleaf-validate --example helmholtz_leaf -- --ps 17,32 --khat 4,16
//! ```
//!
//! Timings are wall times on this machine, reported and never asserted.

use std::hint::black_box;
use std::time::{Duration, Instant};

use tetraleaf_math::orthonormal::{DerivativeScratch, harmonics, harmonics_derivatives};
use tetraleaf_math::radial::bessel_j;
use tetraleaf_math::{Complex, RealScalar};
use tetraleaf_ref::Frame;
use tetraleaf_ref::helmholtz::{Level, Workspace, leaf};
use tetraleaf_validate::SplitMix64;
use tetraleaf_validate::bench::{cores, cpu_model, median_time_per_call, target};

/// The command line.
struct Arguments {
    ps: Vec<usize>,
    khat: Vec<f64>,
    kim: f64,
    points: usize,
    f32: bool,
}

/// Parses the options of the module documentation; exits with a message on anything else.
fn arguments() -> Arguments {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let usage = || -> ! {
        eprintln!(
            "usage: helmholtz_leaf [--ps p,…] [--khat x,…] [--kim y] [--points n] \
             [--precision f64|f32]; got {args:?}"
        );
        std::process::exit(2);
    };
    let mut parsed = Arguments {
        ps: vec![8, 12, 15, 17, 19, 26, 32, 41, 45],
        khat: vec![0.5, 2.0, 4.0, 8.0, 16.0],
        kim: 0.0,
        points: 256,
        f32: false,
    };
    let list = |v: &str| -> Option<Vec<f64>> { v.split(',').map(|x| x.parse().ok()).collect() };
    let mut i = 0;
    while i < args.len() {
        let value = args.get(i + 1).cloned().unwrap_or_else(|| usage());
        let ok = match args[i].as_str() {
            "--ps" => value
                .split(',')
                .map(|x| x.parse().ok())
                .collect::<Option<Vec<usize>>>()
                .map(|v| parsed.ps = v)
                .is_some(),
            "--khat" => list(&value)
                .filter(|v| v.iter().all(|x| *x > 0.0 && x.is_finite()))
                .map(|v| parsed.khat = v)
                .is_some(),
            "--kim" => value
                .parse()
                .ok()
                .filter(|x: &f64| *x >= 0.0 && x.is_finite())
                .map(|x| parsed.kim = x)
                .is_some(),
            "--points" => value
                .parse()
                .ok()
                .filter(|&n: &usize| n > 0)
                .map(|n| parsed.points = n)
                .is_some(),
            "--precision" => match value.as_str() {
                "f64" => Some(false),
                "f32" => Some(true),
                _ => None,
            }
            .map(|f| parsed.f32 = f)
            .is_some(),
            _ => false,
        };
        if !ok {
            usage();
        }
        i += 2;
    }
    parsed
}

/// One case's times, µs per point, and the check of the reassembled P2M.
struct Times {
    p2m: f64,
    radial: f64,
    harmonics: f64,
    exp: f64,
    sum: f64,
    check: String,
    l2p: f64,
    l2p_gradient: f64,
    radial_derivatives: f64,
    harmonics_derivatives: f64,
}

/// µs per point of `pass` over `n` points.
fn per_point(n: usize, mut pass: impl FnMut()) -> f64 {
    1e6 * median_time_per_call(|calls| {
        let start = Instant::now();
        for _ in 0..calls {
            pass();
        }
        start.elapsed()
    }) / n as f64
}

/// The radial functions, harmonics and factor of one point, as the reference's P2M forms
/// them.
fn point<T: RealScalar>(
    p: usize,
    level: &Level<T>,
    x: [T; 3],
    radial: &mut [Complex<T>],
    harm: &mut [T],
) -> T {
    let r = x[0].hypot(x[1]).hypot(x[2]);
    let k = level.k();
    bessel_j(k * r, level.scale(), radial, None);
    harmonics(p, x, harm);
    (k.im * r - level.damping()).exp()
}

/// The reference's sum (`leaf::form` with P2M's addition): per degree w = q e times the radial
/// function, times each harmonic of the degree, into the two planes.
fn sum<T: RealScalar>(
    p: usize,
    w: Complex<T>,
    radial: &[Complex<T>],
    harm: &[T],
    re: &mut [T],
    im: &mut [T],
) {
    for (n, &value) in radial[..=p].iter().enumerate() {
        let wn = w * value;
        let slots = n * n..(n + 1) * (n + 1);
        for (i, &y) in slots.clone().zip(&harm[slots]) {
            let t = wn * y;
            re[i] = re[i] + t.re;
            im[i] = im[i] + t.im;
        }
    }
}

fn run<T: RealScalar>(p: usize, khat: f64, kim: f64, n: usize) -> Times {
    let k = Complex::new(T::from_f64(khat), T::from_f64(kim));
    let level = Level::new(k, T::one());
    let frame = Frame::new([T::zero(); 3], T::one());
    let mut ws = Workspace::new(p);
    let mut rng = SplitMix64::new(0x4c45_4146 ^ p as u64);
    let mut uniform = || T::from_f64(rng.range(-1.0, 1.0));
    let points: Vec<[T; 3]> = (0..n).map(|_| [uniform(), uniform(), uniform()]).collect();
    let charges: Vec<Complex<T>> = (0..n).map(|_| Complex::new(uniform(), uniform())).collect();
    let len = (p + 1) * (p + 1);
    let mut radial = vec![Complex::new(T::zero(), T::zero()); p + 1];
    let mut derivatives = radial.clone();
    let mut harm = vec![T::zero(); len];
    let (mut d_theta, mut d_phi) = (harm.clone(), harm.clone());
    let mut scratch = DerivativeScratch::new(p);

    // The check: the parts, in the reference's order, reassemble its multipole bit for bit.
    let mut reference = vec![T::zero(); 2 * len];
    leaf::p2m(
        p,
        &level,
        &frame,
        &points,
        &charges,
        &mut ws,
        &mut reference,
    );
    let mut parts = vec![T::zero(); 2 * len];
    {
        let (re, im) = parts.split_at_mut(len);
        for (&x, &q) in points.iter().zip(&charges) {
            let e = point(p, &level, x, &mut radial, &mut harm);
            sum(p, q * e, &radial, &harm, re, im);
        }
    }
    let check = if parts == reference {
        "bit for bit".to_string()
    } else {
        let scale = reference
            .iter()
            .fold(0.0f64, |m, v| m.max(RealScalar::to_f64(*v).abs()));
        let diff = reference.iter().zip(&parts).fold(0.0f64, |m, (a, b)| {
            m.max((RealScalar::to_f64(*a) - RealScalar::to_f64(*b)).abs())
        });
        format!("{:.1e}", diff / scale)
    };

    let mut multipole = vec![T::zero(); 2 * len];
    let p2m = per_point(n, || {
        leaf::p2m(
            p,
            &level,
            &frame,
            &points,
            &charges,
            &mut ws,
            &mut multipole,
        );
        black_box(&multipole);
    });
    let radii: Vec<T> = points
        .iter()
        .map(|x| x[0].hypot(x[1]).hypot(x[2]))
        .collect();
    let radial_t = per_point(n, || {
        for &r in &radii {
            bessel_j(level.k() * r, level.scale(), &mut radial, None);
            black_box(&radial);
        }
    });
    let harmonics_t = per_point(n, || {
        for &x in &points {
            harmonics(p, x, &mut harm);
            black_box(&harm);
        }
    });
    let exp_t = per_point(n, || {
        for &r in &radii {
            black_box((level.k().im * r - level.damping()).exp());
        }
    });
    // The sum on one point's radial functions and harmonics, for every charge.
    let e = point(p, &level, points[0], &mut radial, &mut harm);
    let sum_t = per_point(n, || {
        let (re, im) = multipole.split_at_mut(len);
        for &q in &charges {
            sum(p, q * e, &radial, &harm, re, im);
        }
        black_box(&multipole);
    });

    let local: Vec<T> = (0..2 * len).map(|_| uniform()).collect();
    let mut potential = vec![Complex::new(T::zero(), T::zero()); n];
    let mut gradient = vec![[Complex::new(T::zero(), T::zero()); 3]; n];
    let l2p = per_point(n, || {
        leaf::l2p(
            p,
            &level,
            &frame,
            &local,
            &points,
            &mut ws,
            &mut potential,
            None,
        );
        black_box(&potential);
    });
    let l2p_gradient = per_point(n, || {
        leaf::l2p(
            p,
            &level,
            &frame,
            &local,
            &points,
            &mut ws,
            &mut potential,
            Some(&mut gradient),
        );
        black_box(&gradient);
    });
    let radial_derivatives = per_point(n, || {
        for &r in &radii {
            bessel_j(
                level.k() * r,
                level.scale(),
                &mut radial,
                Some(&mut derivatives),
            );
            black_box(&derivatives);
        }
    });
    let harmonics_derivatives_t = per_point(n, || {
        for &x in &points {
            harmonics_derivatives(p, x, &mut harm, &mut d_theta, &mut d_phi, &mut scratch);
            black_box(&d_phi);
        }
    });
    Times {
        p2m,
        radial: radial_t,
        harmonics: harmonics_t,
        exp: exp_t,
        sum: sum_t,
        check,
        l2p,
        l2p_gradient,
        radial_derivatives,
        harmonics_derivatives: harmonics_derivatives_t,
    }
}

fn main() {
    let arguments = arguments();
    let precision = if arguments.f32 { "f32" } else { "f64" };
    println!("# Helmholtz leaf operators per point, by part");
    println!();
    println!(
        "Machine: {}, {}, {}; one thread; {precision}; {} points per batch; Im k̂ = {}",
        cpu_model(),
        cores(),
        target(),
        arguments.points,
        arguments.kim
    );
    println!();
    println!(
        "| p | k̂ | P2M µs | radial µs | harmonics µs | exp µs | sum µs | parts ÷ P2M \
         | P2M ns per term | check | L2P µs | L2P ∇ µs | radial ′ µs | harmonics ′ µs \
         | rest ∇ µs | L2P ∇ ns per term |"
    );
    println!(
        "| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | ---: | ---: \
         | ---: | ---: | ---: | ---: |"
    );
    let start = Instant::now();
    for &p in &arguments.ps {
        for &khat in &arguments.khat {
            let t = if arguments.f32 {
                run::<f32>(p, khat, arguments.kim, arguments.points)
            } else {
                run::<f64>(p, khat, arguments.kim, arguments.points)
            };
            let terms = ((p + 1) * (p + 1)) as f64;
            let parts = t.radial + t.harmonics + t.exp + t.sum;
            let rest = t.l2p_gradient - t.radial_derivatives - t.harmonics_derivatives;
            println!(
                "| {p} | {khat} | {:.3} | {:.3} | {:.3} | {:.3} | {:.3} | {:.2} | {:.1} | {} \
                 | {:.3} | {:.3} | {:.3} | {:.3} | {:.3} | {:.1} |",
                t.p2m,
                t.radial,
                t.harmonics,
                t.exp,
                t.sum,
                parts / t.p2m,
                1e3 * t.p2m / terms,
                t.check,
                t.l2p,
                t.l2p_gradient,
                t.radial_derivatives,
                t.harmonics_derivatives,
                rest,
                1e3 * t.l2p_gradient / terms,
            );
        }
    }
    let elapsed: Duration = start.elapsed();
    println!();
    println!("Total {:.1} s.", elapsed.as_secs_f64());
}
