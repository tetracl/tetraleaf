//! Where a Helmholtz leaf point's P2M and L2P time goes (Helmholtz CH4.9a1;
//! docs/helmholtz/tasks/CH4.9a1.md; CH4.9a2a, CH4.9a2b and CH4.9b1 add the leaf routines' parts, and
//! CH4.9b2 the L2P kernel).
//! Prints Markdown on stdout.
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
//! - **The leaf routines** (Helmholtz CH4.9a2a; docs/helmholtz/tasks/CH4.9a2a.md): the same
//!   points `--lanes` B at a time (default 16, CH4.9a2a's fastest; 1, 4, 8 or 16) through
//!   `tetraleaf_math::orthonormal::LeafHarmonics`, once storing the full layout (`values`) and
//!   once column by column with every value added into an accumulator (`columns`, as a kernel
//!   would stream them), and `tetraleaf_math::radial::bessel_j_real` (real k only). Checked
//!   before timing: the radial values against `bessel_j`'s real parts (the values that differ
//!   are counted; each must be one the routine flushed), and the harmonics' largest difference
//!   from `harmonics` over (n + 1)² u √((2n + 1)/(4π)), the contract of its tests. "sum, leaf
//!   values" is the reference's sum loop on the leaf routines' own values (their harmonics, and
//!   for real k their radial values, both with under-range values flushed) of the batch's first
//!   16 points, cycled over the charges: the sum part above uses one point's reference values,
//!   whose cost in f32 depends on how many are subnormal (CH4.9a1, §2). "projected P2M" adds the
//!   new parts: the column-wise harmonics, the real radial functions (the reference's for
//!   complex k), the factor's `exp` and that sum; the call's overhead is left out (the
//!   reference's parts add up to 0.95–1.05 of its call in f64, CH4.9a1).
//! - **The P2M kernel** (Helmholtz CH4.9a2b; docs/helmholtz/tasks/CH4.9a2b.md): the batch as one
//!   leaf through `tetraleaf_tables::helmholtz::leaf::LeafP2m` (16 points at a time on the leaf
//!   routines, the sum vectorised over points), checked against the reference's multipole first
//!   ("kernel check": the largest difference over the largest value), its time over the
//!   reference's and the projection's, and "kernel sum": the kernel less the routines' parts
//!   (the columns, the radial functions, the factor's `exp`), its sum's share.
//! - **The derivatives** (Helmholtz CH4.9b1; docs/helmholtz/tasks/CH4.9b1.md): the harmonics with
//!   their derivatives column by column (`LeafHarmonics::derivative_columns`, every value summed),
//!   checked first against `harmonics_derivatives` ("derivatives check": the largest difference of
//!   ∂θ and (1/sin θ) ∂φ over (n + 1)³ u √((2n + 1)/(4π)), the contract of its tests), and for real
//!   k `bessel_j_real_derivatives` (checked against `bessel_j`'s derivatives: bit for bit, or in f32
//!   beside a flushed value within (n + 1) τ / s). "sum ∇, leaf values" is the reference's L2P sum
//!   with gradients (`evaluate`'s loop and the gradient's Cartesian assembly, restated) on the
//!   routines' values of the batch's first 16 points, cycled over the targets; "projected L2P ∇"
//!   adds the derivative columns, the radial functions with derivatives (the reference's for
//!   complex k), the factor's `exp` and that sum, over the reference's L2P with gradients.
//! - **The L2P kernel** (Helmholtz CH4.9b2; docs/helmholtz/tasks/CH4.9b2.md): the batch as one leaf
//!   through `tetraleaf_tables::helmholtz::leaf::LeafL2p` without and with gradients (16 targets at a
//!   time on the leaf routines, the sum vectorised over targets), checked against the reference's
//!   L2P with gradients first ("L2P kernel check": the largest difference over the largest value,
//!   of the potentials and of the gradients, the larger), its times over the reference's and the
//!   projection's, and "∇ kernel sum": the kernel with gradients less the routines' parts (the
//!   derivative columns, the radial functions with derivatives, the factor's `exp`), in ns per term
//!   and in GFLOP/s by the sum's operation count, (p + 1) (12.5 p + 8) per target
//!   (`l2p_sum_flops`), the ISA question's measure against the build's peak.
//!
//! Every figure is µs per point, the median over batches (`tetraleaf_validate::bench::
//! median_time_per_call`); "ns per term" divides by (p + 1)². One precision per launch
//! (`--precision f64|f32`, default f64). Run in release mode on one core of a quiet machine:
//!
//! ```text
//! cargo run --release -p tetraleaf-validate --example helmholtz_leaf -- --ps 17,32 --khat 4,16 --lanes 8
//! ```
//!
//! Timings are wall times on this machine, reported and never asserted.

use std::hint::black_box;
use std::time::{Duration, Instant};

use tetraleaf_math::orthonormal::{
    DerivativeScratch, LeafHarmonics, LeafScratch, harmonics, harmonics_derivatives,
};
use tetraleaf_math::radial::{bessel_j, bessel_j_real, bessel_j_real_derivatives};
use tetraleaf_math::{Complex, Layout, RealScalar};
use tetraleaf_ref::Frame;
use tetraleaf_ref::helmholtz::{Level, Workspace, leaf};
use tetraleaf_tables::helmholtz::leaf::{
    L2pScratch, LeafCharges, LeafL2p, LeafOutput, LeafP2m, P2mScratch,
};
use tetraleaf_validate::SplitMix64;
use tetraleaf_validate::bench::{cores, cpu_model, median_time_per_call, target};

/// The command line.
struct Arguments {
    ps: Vec<usize>,
    khat: Vec<f64>,
    kim: f64,
    points: usize,
    f32: bool,
    lanes: usize,
}

/// Parses the options of the module documentation; exits with a message on anything else.
fn arguments() -> Arguments {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let usage = || -> ! {
        eprintln!(
            "usage: helmholtz_leaf [--ps p,…] [--khat x,…] [--kim y] [--points n] \
             [--precision f64|f32] [--lanes 1|4|8|16]; got {args:?}"
        );
        std::process::exit(2);
    };
    let mut parsed = Arguments {
        ps: vec![8, 12, 15, 17, 19, 26, 32, 41, 45],
        khat: vec![0.5, 2.0, 4.0, 8.0, 16.0],
        kim: 0.0,
        points: 256,
        f32: false,
        lanes: 16,
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
            "--lanes" => value
                .parse()
                .ok()
                .filter(|b: &usize| [1, 4, 8, 16].contains(b))
                .map(|b| parsed.lanes = b)
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
    leaf: Leaf,
    /// The sum on the leaf routines' values.
    sum_leaf: f64,
    /// The batched P2M kernel over the batch as one leaf (Helmholtz CH4.9a2b).
    kernel: f64,
    /// The kernel's largest difference from the reference's multipole over its largest value.
    kernel_check: f64,
    /// The reference's L2P sum with gradients on the leaf routines' values (Helmholtz CH4.9b1).
    sum_gradient_leaf: f64,
    /// The batched L2P kernel over the batch as one leaf (Helmholtz CH4.9b2): potentials only,
    /// and with gradients.
    l2p_kernel: f64,
    l2p_kernel_gradient: f64,
    /// The kernel's largest difference from the reference's potentials and gradients, each over
    /// its largest value; the larger of the two.
    l2p_kernel_check: f64,
}

/// The leaf routines' times, µs per point, and their checks (module documentation).
struct Leaf {
    values: f64,
    columns: f64,
    /// `None` for complex k.
    radial: Option<f64>,
    /// The radial values that differ from `bessel_j`'s, each one flushed; `None` for complex k.
    radial_differ: Option<usize>,
    /// The harmonics' largest difference over (n + 1)² u √((2n + 1)/(4π)).
    harmonics_error: f64,
    /// The harmonics with their derivatives column by column (Helmholtz CH4.9b1).
    derivative_columns: f64,
    /// The derivatives' largest difference over (n + 1)³ u √((2n + 1)/(4π)).
    derivatives_error: f64,
    /// The real radial functions with their derivatives; `None` for complex k.
    radial_derivatives: Option<f64>,
    /// The radial derivatives that differ from `bessel_j`'s, each beside a flushed value;
    /// `None` for complex k.
    radial_derivatives_differ: Option<usize>,
}

/// One target of the L2P sum with gradients: the radial functions and their derivatives, the
/// harmonics with their derivatives (value, ∂θ, (1/sin θ) ∂φ), the angles of `Polar::of` and the
/// factor e.
struct Target<T> {
    r: T,
    angles: [T; 4],
    e: T,
    values: Vec<Complex<T>>,
    derivatives: Vec<Complex<T>>,
    harmonics: [Vec<T>; 3],
}

/// r and the angles [cos θ, sin θ, cos φ, sin φ] of `Polar::of` (by `hypot`; θ = 0 at the centre,
/// φ = 0 on the z-axis).
fn polar<T: RealScalar>(x: [T; 3]) -> (T, [T; 4]) {
    let rho = x[0].hypot(x[1]);
    let r = rho.hypot(x[2]);
    if r == T::zero() {
        return (r, [T::one(), T::zero(), T::one(), T::zero()]);
    }
    let (cp, sp) = if rho > T::zero() {
        (x[0] / rho, x[1] / rho)
    } else {
        (T::one(), T::zero())
    };
    (r, [x[2] / r, rho / r, cp, sp])
}

/// The reference's L2P sum with gradients at one target (`leaf::evaluate` with its gradient and
/// `Polar::cartesian`, restated as `sum` restates `form`): per degree the local's coefficients times
/// the harmonics and their derivatives, then the radial function, its derivative and the
/// tangential factor, the gradient in Cartesian components, all times e.
fn sum_gradient<T: RealScalar>(
    p: usize,
    level: &Level<T>,
    local: (&[T], &[T]),
    t: &Target<T>,
) -> (Complex<T>, [Complex<T>; 3]) {
    let zero = Complex::new(T::zero(), T::zero());
    let k = level.k();
    let at_centre = k * (T::one() / (T::from_f64(3.0) * level.scale()));
    let (mut phi, mut radial, mut polar, mut azimuthal) = (zero, zero, zero, zero);
    for n in 0..=p {
        let (mut value, mut theta, mut phi_sum) = (zero, zero, zero);
        for i in n * n..(n + 1) * (n + 1) {
            let c = Complex::new(local.0[i], local.1[i]);
            value = value + c * t.harmonics[0][i];
            theta = theta + c * t.harmonics[1][i];
            phi_sum = phi_sum + c * t.harmonics[2][i];
        }
        let tangential = if t.r > T::zero() {
            t.values[n] * (T::one() / t.r)
        } else if n == 1 {
            at_centre
        } else {
            zero
        };
        phi = phi + t.values[n] * value;
        radial = radial + t.derivatives[n] * value;
        polar = polar + tangential * theta;
        azimuthal = azimuthal + tangential * phi_sum;
    }
    let [ct, st, cp, sp] = t.angles;
    let kr = k * radial;
    let gradient = [
        kr * (st * cp) + polar * (ct * cp) + azimuthal * (-sp),
        kr * (st * sp) + polar * (ct * sp) + azimuthal * cp,
        kr * ct + polar * (-st),
    ];
    (phi * t.e, gradient.map(|g| g * t.e))
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

fn run<T: RealScalar>(p: usize, khat: f64, kim: f64, n: usize, lanes: usize) -> Times {
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

    // The same sum on the leaf routines' own values (`LeafHarmonics`, and for real k
    // `bessel_j_real`, both flushing under-range values) of 16 points, cycled over the charges:
    // the operands a kernel on the leaf routines would sum, and a mix of points rather than one.
    let sum_leaf = {
        let table = LeafHarmonics::<T>::new(p);
        let mut scratch = LeafScratch::<T, 1>::new(p);
        let mut lanes = vec![[T::zero(); 1]; len];
        let mut radial_lanes = vec![[T::zero(); 1]; p + 1];
        let data: Vec<(T, Vec<Complex<T>>, Vec<T>)> = points
            .iter()
            .take(16)
            .map(|&x| {
                let r = x[0].hypot(x[1]).hypot(x[2]);
                let k = level.k();
                let mut values = vec![Complex::new(T::zero(), T::zero()); p + 1];
                if k.im == T::zero() {
                    bessel_j_real(&[k.re * r], level.scale(), &mut radial_lanes);
                    for (v, l) in values.iter_mut().zip(&radial_lanes) {
                        *v = Complex::new(l[0], T::zero());
                    }
                } else {
                    bessel_j(k * r, level.scale(), &mut values, None);
                }
                table.values(p, &[x], &mut scratch, &mut lanes);
                let e = (k.im * r - level.damping()).exp();
                (e, values, lanes.iter().map(|v| v[0]).collect())
            })
            .collect();
        per_point(n, || {
            let (re, im) = multipole.split_at_mut(len);
            for (&q, (e, values, harm)) in charges.iter().zip(data.iter().cycle()) {
                sum(p, q * *e, values, harm, re, im);
            }
            black_box(&multipole);
        })
    };

    // The batched P2M kernel (Helmholtz CH4.9a2b): the batch as one leaf, checked against the
    // reference first.
    let (kernel, kernel_check) = {
        let table = LeafP2m::<T>::new(p);
        let mut scratch = P2mScratch::new(p, 1);
        let values: Vec<T> = charges.iter().flat_map(|q| [q.re, q.im]).collect();
        let leaf_charges = LeafCharges {
            values: &values,
            vectors: 1,
            complex: true,
        };
        let mut batched = vec![T::zero(); 2 * len];
        table.p2m(p, &level, &points, leaf_charges, &mut scratch, &mut batched);
        let scale = reference
            .iter()
            .fold(0.0f64, |m, v| m.max(RealScalar::to_f64(*v).abs()));
        let diff = reference.iter().zip(&batched).fold(0.0f64, |m, (a, b)| {
            m.max((RealScalar::to_f64(*a) - RealScalar::to_f64(*b)).abs())
        });
        let time = per_point(n, || {
            table.p2m(p, &level, &points, leaf_charges, &mut scratch, &mut batched);
            black_box(&batched);
        });
        (time, diff / scale)
    };

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
    // The L2P sum with gradients on the leaf routines' values (`LeafHarmonics::derivatives`, and
    // for real k `bessel_j_real_derivatives`) of 16 points, cycled over the targets.
    let sum_gradient_leaf = {
        let table = LeafHarmonics::<T>::new(p);
        let mut scratch = LeafScratch::<T, 1>::new(p);
        let mut out: [Vec<[T; 1]>; 3] = core::array::from_fn(|_| vec![[T::zero(); 1]; len]);
        let (mut radial_lanes, mut derivative_lanes) =
            (vec![[T::zero(); 1]; p + 1], vec![[T::zero(); 1]; p + 1]);
        let k = level.k();
        let data: Vec<Target<T>> = points
            .iter()
            .take(16)
            .map(|&x| {
                let (r, angles) = polar(x);
                let zero = Complex::new(T::zero(), T::zero());
                let (mut values, mut derivatives) = (vec![zero; p + 1], vec![zero; p + 1]);
                if k.im == T::zero() {
                    bessel_j_real_derivatives(
                        &[k.re * r],
                        level.scale(),
                        &mut radial_lanes,
                        &mut derivative_lanes,
                    );
                    for n in 0..=p {
                        values[n] = Complex::new(radial_lanes[n][0], T::zero());
                        derivatives[n] = Complex::new(derivative_lanes[n][0], T::zero());
                    }
                } else {
                    bessel_j(k * r, level.scale(), &mut values, Some(&mut derivatives));
                }
                {
                    let [v, t, f] = &mut out;
                    table.derivatives(p, &[x], &mut scratch, v, t, f);
                }
                Target {
                    r,
                    angles,
                    e: (k.im * r - level.damping()).exp(),
                    values,
                    derivatives,
                    harmonics: out.clone().map(|o| o.iter().map(|v| v[0]).collect()),
                }
            })
            .collect();
        let (re, im) = local.split_at(len);
        // The check: the restated sum on the reference's own values (`bessel_j` with derivatives,
        // `harmonics_derivatives`) is `leaf::l2p` with gradients, bit for bit, at the 16 points.
        {
            let zero = Complex::new(T::zero(), T::zero());
            let mut scratch = DerivativeScratch::new(p);
            for &x in points.iter().take(16) {
                let (r, angles) = polar(x);
                let (mut values, mut derivatives) = (vec![zero; p + 1], vec![zero; p + 1]);
                bessel_j(k * r, level.scale(), &mut values, Some(&mut derivatives));
                let mut harmonics: [Vec<T>; 3] = core::array::from_fn(|_| vec![T::zero(); len]);
                {
                    let [v, t, f] = &mut harmonics;
                    harmonics_derivatives(p, x, v, t, f, &mut scratch);
                }
                let target = Target {
                    r,
                    angles,
                    e: (k.im * r - level.damping()).exp(),
                    values,
                    derivatives,
                    harmonics,
                };
                let (phi, g) = sum_gradient(p, &level, (re, im), &target);
                let (mut want, mut want_g) = ([zero], [[zero; 3]]);
                leaf::l2p(
                    p,
                    &level,
                    &frame,
                    &local,
                    &[x],
                    &mut ws,
                    &mut want,
                    Some(&mut want_g),
                );
                let same = |a: Complex<T>, b: Complex<T>| a.re == b.re && a.im == b.im;
                assert!(
                    same(phi + zero, want[0]) && (0..3).all(|i| same(g[i] + zero, want_g[0][i])),
                    "the restated L2P sum with gradients differs from leaf::l2p"
                );
            }
        }
        per_point(n, || {
            for (out, t) in potential.iter_mut().zip(data.iter().cycle()) {
                let (phi, g) = sum_gradient(p, &level, (re, im), t);
                *out = *out + phi + g[0] + g[1] + g[2];
            }
            black_box(&potential);
        })
    };
    let leaf = match lanes {
        1 => leaf_parts::<T, 1>(p, &level, &points),
        4 => leaf_parts::<T, 4>(p, &level, &points),
        8 => leaf_parts::<T, 8>(p, &level, &points),
        _ => leaf_parts::<T, 16>(p, &level, &points),
    };
    // The batched L2P kernel (Helmholtz CH4.9b2): the batch as one leaf, without and with
    // gradients, checked against the reference first.
    let (l2p_kernel, l2p_kernel_gradient, l2p_kernel_check) = {
        let table = LeafL2p::<T>::new(p);
        let mut scratch = L2pScratch::new(p, 1);
        let zero = Complex::new(T::zero(), T::zero());
        let (mut phi, mut grad) = (vec![zero; n], vec![[zero; 3]; n]);
        leaf::l2p(
            p,
            &level,
            &frame,
            &local,
            &points,
            &mut ws,
            &mut phi,
            Some(&mut grad),
        );
        let mut out = vec![T::zero(); 8 * n];
        table.l2p(
            p,
            &level,
            &points,
            &local,
            &mut scratch,
            leaf_output(&mut out, true),
        );
        let largest = |values: &mut dyn Iterator<Item = f64>| values.fold(0.0f64, f64::max);
        let parts = |c: Complex<T>| [c.re.to_f64(), c.im.to_f64()];
        let phi_scale = largest(&mut phi.iter().flat_map(|&c| parts(c)).map(f64::abs));
        let phi_diff = largest(&mut (0..n).flat_map(|i| {
            let got = [out[i].to_f64(), out[n + i].to_f64()];
            let want = parts(phi[i]);
            [(got[0] - want[0]).abs(), (got[1] - want[1]).abs()]
        }));
        let grad_scale = largest(
            &mut grad
                .iter()
                .flat_map(|g| g.iter().flat_map(|&c| parts(c)))
                .map(f64::abs),
        );
        let grad_diff = largest(&mut (0..n).flat_map(|i| {
            let (out, grad) = (&out, &grad);
            (0..3).flat_map(move |a| {
                let got = [out[2 * n + 3 * i + a], out[2 * n + 3 * (n + i) + a]];
                let want = parts(grad[i][a]);
                [
                    (got[0].to_f64() - want[0]).abs(),
                    (got[1].to_f64() - want[1]).abs(),
                ]
            })
        }));
        let check = (phi_diff / phi_scale).max(grad_diff / grad_scale);
        let mut plain = vec![T::zero(); 2 * n];
        let time = per_point(n, || {
            table.l2p(
                p,
                &level,
                &points,
                &local,
                &mut scratch,
                leaf_output(&mut plain, false),
            );
            black_box(&plain);
        });
        let time_gradient = per_point(n, || {
            table.l2p(
                p,
                &level,
                &points,
                &local,
                &mut scratch,
                leaf_output(&mut out, true),
            );
            black_box(&out);
        });
        (time, time_gradient, check)
    };
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
        leaf,
        sum_leaf,
        kernel,
        kernel_check,
        sum_gradient_leaf,
        l2p_kernel,
        l2p_kernel_gradient,
        l2p_kernel_check,
    }
}

/// One vector's target output chunk for `LeafL2p::l2p`.
fn leaf_output<T>(values: &mut [T], gradients: bool) -> LeafOutput<'_, T> {
    LeafOutput {
        values,
        vectors: 1,
        gradients,
    }
}

/// The arithmetic of the L2P kernel's sum with gradients per target (Helmholtz CH4.9b2,
/// `LeafL2p::l2p`'s step 2): per order m ≥ 1 and degree n, m U once and on each of the two planes
/// A, B and the three sums (seven products and five sums), per degree n of order 0 two products
/// and two sums per plane: (p + 1) (12.5 p + 8) operations, about 12.5 per term. The degree loop
/// and the gradient's assembly (O(p) per target) are left out.
fn l2p_sum_flops(p: usize) -> f64 {
    let p = p as f64;
    (p + 1.0) * (12.5 * p + 8.0)
}

/// The points in batches of B lanes, the last batch padded with the first point.
fn batches<T: RealScalar, const B: usize>(points: &[[T; 3]]) -> Vec<[[T; 3]; B]> {
    points
        .chunks(B)
        .map(|chunk| core::array::from_fn(|j| *chunk.get(j).unwrap_or(&points[0])))
        .collect()
}

/// The leaf routines over the points, B lanes at a time (module documentation).
fn leaf_parts<T: RealScalar, const B: usize>(
    p: usize,
    level: &Level<T>,
    points: &[[T; 3]],
) -> Leaf {
    let n = points.len();
    let layout = Layout::new(p);
    let table = LeafHarmonics::<T>::new(p);
    let mut scratch = LeafScratch::<T, B>::new(p);
    let blocks = batches::<T, B>(points);
    let mut out = vec![[T::zero(); B]; layout.len()];

    // The harmonics' check against the reference, per point and term.
    let unit = T::epsilon().to_f64() / 2.0;
    let mut reference = vec![T::zero(); layout.len()];
    let mut harmonics_error = 0.0f64;
    for (b, block) in blocks.iter().enumerate() {
        table.values(p, block, &mut scratch, &mut out);
        for (j, &x) in block.iter().enumerate().take(n - (b * B).min(n)) {
            harmonics(p, x, &mut reference);
            for (degree, range) in layout.degrees() {
                let scale = ((degree + 1) * (degree + 1)) as f64
                    * unit
                    * ((2 * degree + 1) as f64 / (4.0 * std::f64::consts::PI)).sqrt();
                for i in range {
                    let d = (out[i][j].to_f64() - reference[i].to_f64()).abs();
                    harmonics_error = harmonics_error.max(d / scale);
                }
            }
        }
    }
    let values = per_point(n, || {
        for block in &blocks {
            table.values(p, block, &mut scratch, &mut out);
            black_box(&out);
        }
    });
    let mut sink = [T::zero(); B];
    let columns = per_point(n, || {
        for block in &blocks {
            table.columns(p, block, &mut scratch, |_, q, cos, sin| {
                for row in q {
                    for j in 0..B {
                        sink[j] = sink[j] + row[j] * cos[j] + row[j] * sin[j];
                    }
                }
            });
        }
        black_box(&sink);
    });

    let (radial, radial_differ) = if level.k().im == T::zero() {
        let k = level.k().re;
        let xs: Vec<[T; B]> = blocks
            .iter()
            .map(|block| {
                core::array::from_fn(|j| {
                    let [a, b, c] = block[j];
                    k * a.hypot(b).hypot(c)
                })
            })
            .collect();
        let mut lanes = vec![[T::zero(); B]; p + 1];
        let mut one = vec![Complex::new(T::zero(), T::zero()); p + 1];
        let mut differ = 0;
        for x in &xs {
            bessel_j_real(x, level.scale(), &mut lanes);
            for (j, &xj) in x.iter().enumerate() {
                bessel_j(Complex::new(xj, T::zero()), level.scale(), &mut one, None);
                for (order, (row, want)) in lanes.iter().zip(&one).enumerate() {
                    if row[j] != want.re {
                        // Only f32's flush may differ: a zero where the reference is below
                        // τ = 2⁻¹⁰⁰ at an order n ≥ x.
                        let flushed = row[j] == T::zero()
                            && want.re.abs().to_f64() < 2f64.powi(-100)
                            && order as f64 >= xj.to_f64();
                        assert!(flushed, "bessel_j_real differs from bessel_j");
                        differ += 1;
                    }
                }
            }
        }
        let time = per_point(n, || {
            for x in &xs {
                bessel_j_real(x, level.scale(), &mut lanes);
                black_box(&lanes);
            }
        });
        (Some(time), Some(differ))
    } else {
        (None, None)
    };
    // The derivatives (Helmholtz CH4.9b1): their check against `harmonics_derivatives` per point
    // and term, then the columns timed with every value summed.
    let mut derivatives_error = 0.0f64;
    {
        let mut out: [Vec<[T; B]>; 3] =
            core::array::from_fn(|_| vec![[T::zero(); B]; layout.len()]);
        let mut reference: [Vec<T>; 3] = core::array::from_fn(|_| vec![T::zero(); layout.len()]);
        let mut ref_scratch = DerivativeScratch::new(p);
        for (b, block) in blocks.iter().enumerate() {
            {
                let [v, t, f] = &mut out;
                table.derivatives(p, block, &mut scratch, v, t, f);
            }
            for (j, &x) in block.iter().enumerate().take(n - (b * B).min(n)) {
                {
                    let [v, t, f] = &mut reference;
                    harmonics_derivatives(p, x, v, t, f, &mut ref_scratch);
                }
                for (degree, range) in layout.degrees() {
                    let d = (degree + 1) as f64;
                    let scale = d
                        * d
                        * d
                        * unit
                        * ((2 * degree + 1) as f64 / (4.0 * std::f64::consts::PI)).sqrt();
                    for i in range {
                        for part in 1..3 {
                            let e = (out[part][i][j].to_f64() - reference[part][i].to_f64()).abs();
                            derivatives_error = derivatives_error.max(e / scale);
                        }
                    }
                }
            }
        }
    }
    // The contract of `tests/orthonormal/leaf.rs` (C = 10), checked before timing.
    assert!(
        derivatives_error <= 10.0,
        "the leaf derivatives outside their contract: {derivatives_error}"
    );
    let derivative_columns = per_point(n, || {
        for block in &blocks {
            table.derivative_columns(p, block, &mut scratch, |_, q, d, u, cos, sin| {
                for ((row, dr), ur) in q.iter().zip(d).zip(u) {
                    for j in 0..B {
                        sink[j] = sink[j] + row[j] * cos[j] + dr[j] * sin[j] + ur[j] * cos[j];
                    }
                }
            });
        }
        black_box(&sink);
    });
    let (radial_derivatives, radial_derivatives_differ) = if level.k().im == T::zero() {
        let k = level.k().re;
        let xs: Vec<[T; B]> = blocks
            .iter()
            .map(|block| {
                core::array::from_fn(|j| {
                    let [a, b, c] = block[j];
                    k * a.hypot(b).hypot(c)
                })
            })
            .collect();
        let (mut lanes, mut d_lanes) = (vec![[T::zero(); B]; p + 1], vec![[T::zero(); B]; p + 1]);
        let zero = Complex::new(T::zero(), T::zero());
        let (mut one, mut d_one) = (vec![zero; p + 1], vec![zero; p + 1]);
        let mut differ = 0;
        for x in &xs {
            bessel_j_real_derivatives(x, level.scale(), &mut lanes, &mut d_lanes);
            for (j, &xj) in x.iter().enumerate() {
                bessel_j(
                    Complex::new(xj, T::zero()),
                    level.scale(),
                    &mut one,
                    Some(&mut d_one),
                );
                for (order, (row, want)) in d_lanes.iter().zip(&d_one).enumerate() {
                    if row[j] != want.re {
                        // Only f32's flush may differ: beside a flushed neighbour, within
                        // (n + 1) τ / s of the reference.
                        let near = (row[j].to_f64() - want.re.to_f64()).abs()
                            <= (order + 1) as f64 * 2f64.powi(-100) / level.scale().to_f64();
                        assert!(near, "bessel_j_real_derivatives differs from bessel_j");
                        differ += 1;
                    }
                }
            }
        }
        let time = per_point(n, || {
            for x in &xs {
                bessel_j_real_derivatives(x, level.scale(), &mut lanes, &mut d_lanes);
                black_box(&d_lanes);
            }
        });
        (Some(time), Some(differ))
    } else {
        (None, None)
    };
    Leaf {
        values,
        columns,
        radial,
        radial_differ,
        harmonics_error,
        derivative_columns,
        derivatives_error,
        radial_derivatives,
        radial_derivatives_differ,
    }
}

fn main() {
    let arguments = arguments();
    let precision = if arguments.f32 { "f32" } else { "f64" };
    println!("# Helmholtz leaf operators per point, by part");
    println!();
    println!(
        "Machine: {}, {}, {}; one thread; {precision}; {} points per batch; Im k̂ = {}; \
         {} lanes",
        cpu_model(),
        cores(),
        target(),
        arguments.points,
        arguments.kim,
        arguments.lanes
    );
    println!();
    println!(
        "| p | k̂ | P2M µs | radial µs | harmonics µs | exp µs | sum µs | parts ÷ P2M \
         | P2M ns per term | check | L2P µs | L2P ∇ µs | radial ′ µs | harmonics ′ µs \
         | rest ∇ µs | L2P ∇ ns per term | leaf harmonics µs | leaf columns µs | real radial µs \
         | radial check | harmonics check | sum, leaf values µs | projected P2M µs \
         | projected ÷ P2M | P2M kernel µs | kernel ÷ P2M | kernel ÷ projected | kernel sum µs \
         | kernel check | leaf columns ′ µs | derivatives check | real radial ′ µs | radial ′ check \
         | sum ∇, leaf values µs | projected L2P ∇ µs | projected ÷ L2P ∇ | L2P kernel µs \
         | L2P kernel ÷ L2P | L2P ∇ kernel µs | ∇ kernel ÷ L2P ∇ | ∇ kernel ÷ projected \
         | ∇ kernel sum µs | ∇ kernel sum ns per term | ∇ kernel sum GFLOP/s | L2P kernel check |"
    );
    println!(
        "| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | ---: | ---: \
         | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: \
         | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | ---: | ---: | ---: \
         | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |"
    );
    let start = Instant::now();
    for &p in &arguments.ps {
        for &khat in &arguments.khat {
            let t = if arguments.f32 {
                run::<f32>(p, khat, arguments.kim, arguments.points, arguments.lanes)
            } else {
                run::<f64>(p, khat, arguments.kim, arguments.points, arguments.lanes)
            };
            let terms = ((p + 1) * (p + 1)) as f64;
            let parts = t.radial + t.harmonics + t.exp + t.sum;
            let rest = t.l2p_gradient - t.radial_derivatives - t.harmonics_derivatives;
            let leaf = &t.leaf;
            let radial = leaf.radial.unwrap_or(t.radial);
            let projected = leaf.columns + radial + t.exp + t.sum_leaf;
            let or_dash = |v: Option<String>| v.unwrap_or_else(|| "–".to_string());
            let radial_derivatives = leaf.radial_derivatives.unwrap_or(t.radial_derivatives);
            let projected_gradient =
                leaf.derivative_columns + radial_derivatives + t.exp + t.sum_gradient_leaf;
            // The L2P kernel with gradients less its routines' parts: its sum (with the
            // directions, the degree loop and the output).
            let kernel_sum =
                t.l2p_kernel_gradient - (leaf.derivative_columns + radial_derivatives + t.exp);
            println!(
                "| {p} | {khat} | {:.3} | {:.3} | {:.3} | {:.3} | {:.3} | {:.2} | {:.1} | {} \
                 | {:.3} | {:.3} | {:.3} | {:.3} | {:.3} | {:.1} | {:.3} | {:.3} | {} | {} \
                 | {:.2} | {:.3} | {:.3} | {:.2} | {:.3} | {:.2} | {:.2} | {:.3} | {:.1e} \
                 | {:.3} | {:.2} | {} | {} | {:.3} | {:.3} | {:.2} | {:.3} | {:.3} | {:.3} \
                 | {:.3} | {:.3} | {:.3} | {:.3} | {:.2} | {:.1e} |",
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
                leaf.values,
                leaf.columns,
                or_dash(leaf.radial.map(|r| format!("{r:.3}"))),
                or_dash(leaf.radial_differ.map(|d| if d == 0 {
                    "bit for bit".to_string()
                } else {
                    format!("{d} flushed")
                })),
                leaf.harmonics_error,
                t.sum_leaf,
                projected,
                projected / t.p2m,
                t.kernel,
                t.kernel / t.p2m,
                t.kernel / projected,
                t.kernel - (leaf.columns + radial + t.exp),
                t.kernel_check,
                leaf.derivative_columns,
                leaf.derivatives_error,
                or_dash(leaf.radial_derivatives.map(|r| format!("{r:.3}"))),
                or_dash(leaf.radial_derivatives_differ.map(|d| if d == 0 {
                    "bit for bit".to_string()
                } else {
                    format!("{d} beside flushed")
                })),
                t.sum_gradient_leaf,
                projected_gradient,
                projected_gradient / t.l2p_gradient,
                t.l2p_kernel,
                t.l2p_kernel / t.l2p,
                t.l2p_kernel_gradient,
                t.l2p_kernel_gradient / t.l2p_gradient,
                t.l2p_kernel_gradient / projected_gradient,
                kernel_sum,
                1e3 * kernel_sum / terms,
                l2p_sum_flops(p) / (1e3 * kernel_sum),
                t.l2p_kernel_check,
            );
        }
    }
    let elapsed: Duration = start.elapsed();
    println!();
    println!("Total {:.1} s.", elapsed.as_secs_f64());
}
