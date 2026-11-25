//! Where the Helmholtz FMM's errors concentrate (Helmholtz CH4.0b;
//! docs/helmholtz/tasks/CH4.0b.md). Prints Markdown on stdout.
//!
//! ```text
//! helmholtz_boxes [CASE] [--eps e,…] [--m2l rotation|default|defaults,…]
//!     [--study single|levels|m2l|levels,m2l] [--class-runs] [--refine n] [--raise-by d,…]
//!     [--raise-levels l,…] [--margin m] [--offset x,y,z] [--leaf-size n] [--raise-factor f]
//!     [--worst n] [--per-box s] [--full-max n]
//!     [--rule calibrated|ch3-3|ch4-3c2|ch4-3c3|γ,α,β,C,a,b,d[,cap[,margin]]] [--threads t]
//! ```
//!
//! Helmholtz CH4.4b (the crosscheck mapped to levels and classes): `--m2l defaults` runs
//! `Case::defaults`, every setting the library's (the diagonal form's time rule, the waves, the
//! near rules); there `--study levels` raises one level at a time by an override
//! (`M2lOverride`, its degree the rule's at ε × `--raise-factor` and its M2L RCR), the
//! tolerance and every other level's choice kept. `--study m2l` runs each level whose M2L is
//! not RCR on RCR at its degree as run, and on the form's levels the classes within a decade
//! below the bound's threshold (E/(θ ε) in (0.1, 1]) and the next decade on RCR, each group in
//! one run (`--class-runs`: each class of the first decade alone too); a table of the levels as
//! run and one of these runs' parts (their share of the error, their worst box and ℓ2 over the
//! first run's, their part at the worst ∇φ target). `--refine n` adds every point of the n worst
//! boxes of φ and of ∇φ (by the sampled reference) to the reference, so their maxima are
//! exact; `--refine-cells i,j,k;…` adds every point of the named boxes of the deepest level
//! (another run's worst boxes, CH4.4a's on W3); `--all-levels` makes `--study levels` raise
//! the levels without V pairs too (the deepest level's translations and leaf operators);
//! `--refine-faces` adds every point of the deepest level's boxes at the centres of the
//! domain's faces (W1's six tangent points with the points' bounding cube, where CH4.4a found
//! the worst boxes; a sample of two per box missed them). `--raise-by d,…` runs the level of `--raise-levels` (default: the level with the
//! largest part at the worst ∇φ box, from `--study levels`) at its degree as run plus d under
//! the setting's M2L, a remedy's FMM evidence (`--raise-levels v`: every level with V pairs;
//! `--raise-together`: one run per d with all of them raised at once). `--reference-cache dir` stores the sampled
//! reference (before `--refine`) in dir and reads it back in a launch of the same case and grid
//! (another ε), whose header must match; the refinement is then computed for that launch.
//!
//! `--rule ch3-3` (Helmholtz CH4.3c2) runs CH3.3's degree rule (`DegreeRule::CH3_3`) in place
//! of the library's default under the same tolerance, for the errors before the refit;
//! `--rule ch4-3c2` (Helmholtz CH4.3c3) CH4.3c2's rule (`DegreeRule::CH4_3C2`), the default
//! before CH4.3c3's refit; `--rule ch4-3c3` (Helmholtz CH4.8b) CH4.3c3's rule
//! (`DegreeRule::CH4_3C3`), the default before CH4.8b's margin;
//! `--rule γ,α,β,C,a,b,d[,cap[,margin]]` a `DegreeRule` of those constants (the refit's offset
//! varied; margin 0 unless named, so that CH4.8a's named rules repeat).
//!
//! One case on one rank, against the direct sum (`tetraleaf_validate::helmholtz_boxes`): at
//! every point when N ≤ `--full-max` (default 200,000), else at §9.1's samples plus
//! `--per-box` targets (default 2) per box of the FMM's deepest level. The FMM runs at each
//! ε of `--eps` (default 10⁻⁶) under each M2L of `--m2l` (`rotation`, the default here, or
//! `default`, the library's `Auto`), the first pair first. Per run: the degrees, §9.1's
//! relative ℓ2 errors at the 1,000 samples, the relative ℓ2 errors over the reference's
//! targets, the box measure of the deepest level (the largest error at one of a box's targets
//! over the reference's RMS) over the run's ε and over the ℓ2 error, the boxes above ε, the
//! build and the evaluation. For the first run: the `--worst` worst boxes of ∇φ and of φ
//! with their place (the box's cell and centre, whether it touches the domain's boundary, its
//! worst target's distance from its box's centre on every level from 2 in units of the
//! half-width, √3 at a corner), and the errors by the number of levels on which a target lies
//! within 0.1 of a face.
//!
//! `--study levels` also runs the first case once per level with M2L with that level's
//! degree raised to the rule's at ε × `--raise-factor` (default 10⁻³), and once with every
//! level raised; each level's part of the error is the first run's error minus the raised
//! run's (`helmholtz_boxes::level_shares`), over every target, at the worst boxes' worst
//! targets and by the target's distance from its box's centre on that level. Named degrees
//! change the exponential M2L's precision rule, so the study is read on `--m2l rotation`.
//!
//! CASE is `--workload w` (as `helmholtz_fmm3d`: `w1`, `w2`, `w3`, `w3:<n>`,
//! `w4-<distribution>:<n>`; default `w1`), `--kd x` (default 30), `--kim y` (default 0),
//! `--density random|planewave|constant` (default `random`, §9.1's) and `--charges
//! real|complex` (default complex). `--margin m` and `--offset x,y,z` set the FMM's domain to
//! [−m, 1 + m]³ shifted by the offset (default: the points' bounding cube); `--leaf-size n`
//! the points per leaf (`uniform` for N + 1: leaves on the leaf rule's level alone, no W or X
//! list).
//!
//! Run in release mode on one rank (`mpirun -n 1`), with one BLAS thread. Accuracy only; the
//! times are wall times, reported, never asserted.

use std::time::Instant;

use mpi::Threading;
use mpi::traits::*;
use tetraleaf_exec::helmholtz::{
    ChargeKind, DegreeRule, Degrees, HelmholtzM2l, M2lOverride, OverrideM2l, Wavenumber,
};
use tetraleaf_math::Complex;
use tetraleaf_validate::bench::{cores, cpu_model};
use tetraleaf_validate::fmm_accuracy::Distribution;
use tetraleaf_validate::helmholtz_accuracy::{
    Case, Density, Oracle, Problem, Row, Run, complex_k, run_output,
};
use tetraleaf_validate::helmholtz_boxes::{
    BoxError, Errors, Reference, box_sample, cell, class_ratios, face_closeness, level_parts_at,
    level_shares, offsets, threshold_decade, worst_boxes,
};
use tetraleaf_validate::helmholtz_fmm3d::Grid;
use tetraleaf_validate::helmholtz_workloads::{SEED, Workload};

/// The command line.
struct Arguments {
    workload: String,
    kd: f64,
    kim: f64,
    density: Density,
    charges: ChargeKind,
    eps: Vec<f64>,
    m2l: Vec<Setting>,
    levels: bool,
    /// `--study m2l` (Helmholtz CH4.4b): each level's M2L on RCR, the form's classes by decade.
    m2l_study: bool,
    /// `--class-runs`: each class of the threshold's decade alone on RCR.
    class_runs: bool,
    /// `--refine n`: every point of the n worst boxes of φ and of ∇φ in the reference.
    refine: usize,
    /// `--refine-faces`: every point of the deepest level's boxes at the centres of the
    /// domain's faces (W1's tangent points) in the reference.
    refine_faces: bool,
    /// `--refine-cells i,j,k;…`: every point of these boxes of the deepest level in the reference.
    refine_cells: Vec<[u64; 3]>,
    /// `--all-levels`: `--study levels` raises the levels without V pairs too.
    all_levels: bool,
    /// `--raise-by d,…`: the runs with a level's degree raised by d under the setting's M2L.
    raise_by: Vec<usize>,
    /// `--raise-levels l,…|v`: their levels (`v`, here an empty list: every level with V
    /// pairs; default: the level with the largest part at the worst ∇φ box, from `--study
    /// levels`).
    raise_levels: Option<Vec<usize>>,
    /// `--raise-together`: one run per d with every level of `--raise-levels` raised at once.
    raise_together: bool,
    /// `--reference-cache dir`: the sampled reference stored there, and read back by a launch
    /// of the same case and grid (another ε).
    reference_cache: Option<String>,
    margin: Option<f64>,
    offset: [f64; 3],
    leaf_size: Option<String>,
    raise_factor: f64,
    worst: usize,
    per_box: usize,
    full_max: usize,
    threads: usize,
    /// A degree rule in place of the default (`--rule ch3-3` or the constants, Helmholtz
    /// CH4.3c2).
    degree_rule: Option<DegreeRule>,
}

/// The FMM's settings of a run (`--m2l`).
#[derive(Clone, Copy, Debug, PartialEq)]
enum Setting {
    /// `Case::new` with Rotation (CH4.0b's reference).
    Rotation,
    /// `Case::new` with the library's M2L (`Auto`), the form off.
    Default,
    /// `Case::defaults`: every setting the library's (Helmholtz CH4.4b).
    Defaults,
}

impl Setting {
    fn name(self) -> &'static str {
        match self {
            Setting::Rotation => "Rotation",
            Setting::Default => "default",
            Setting::Defaults => "defaults",
        }
    }
}

/// The workload of a name at kD (W1 and W3 resolved depend on it).
fn workload(name: &str, kd: f64) -> Option<Workload> {
    let count = |n: &str| n.parse::<usize>().ok().filter(|&n| n > 0);
    match name.split_once(':') {
        None => match name {
            "w1" => Some(Workload::W1),
            "w2" => Some(Workload::W2),
            "w3" => Some(Workload::w3_resolved(kd)),
            _ => None,
        },
        Some(("w3", n)) => Some(Workload::W3 { n: count(n)? }),
        Some((w4, n)) => Some(Workload::W4 {
            distribution: Distribution::from_name(w4.strip_prefix("w4-")?)?,
            n: count(n)?,
        }),
    }
}

/// Parses the options of the module documentation; exits with a message on anything else.
fn arguments() -> Arguments {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let usage = || -> ! {
        eprintln!(
            "usage: helmholtz_boxes [--workload w] [--kd x] [--kim y] \
             [--density random|planewave|constant] [--charges real|complex] [--eps e] \
             [--m2l rotation|default|defaults] [--study single|levels|m2l|levels,m2l] \
             [--class-runs] [--refine n] [--refine-faces] [--refine-cells i,j,k;…] \
             [--all-levels] [--raise-by d,…] \
             [--raise-levels l,…|v] [--raise-together] \
             [--reference-cache dir] [--margin m] \
             [--offset x,y,z] [--leaf-size n|uniform] [--raise-factor f] [--worst n] \
             [--per-box s] [--full-max n] \
             [--rule calibrated|ch3-3|ch4-3c2|ch4-3c3|γ,α,β,C,a,b,d[,cap[,margin]]] [--threads t]; got {args:?}"
        );
        std::process::exit(2);
    };
    let mut parsed = Arguments {
        workload: "w1".into(),
        kd: 30.0,
        kim: 0.0,
        density: Density::Random,
        charges: ChargeKind::Complex,
        eps: vec![1e-6],
        m2l: vec![Setting::Rotation],
        levels: false,
        m2l_study: false,
        class_runs: false,
        refine: 0,
        refine_faces: false,
        refine_cells: Vec::new(),
        all_levels: false,
        raise_by: Vec::new(),
        raise_levels: None,
        raise_together: false,
        reference_cache: None,
        margin: None,
        offset: [0.0; 3],
        leaf_size: None,
        raise_factor: 1e-3,
        worst: 8,
        per_box: 2,
        full_max: 200_000,
        threads: 1,
        degree_rule: None,
    };
    let mut rest = &args[..];
    loop {
        if let [flag, tail @ ..] = rest
            && matches!(
                flag.as_str(),
                "--class-runs" | "--refine-faces" | "--raise-together" | "--all-levels"
            )
        {
            match flag.as_str() {
                "--class-runs" => parsed.class_runs = true,
                "--refine-faces" => parsed.refine_faces = true,
                "--all-levels" => parsed.all_levels = true,
                _ => parsed.raise_together = true,
            }
            rest = tail;
            continue;
        }
        let [flag, value, tail @ ..] = rest else {
            break;
        };
        let ok = match flag.as_str() {
            "--workload" => {
                parsed.workload = value.clone();
                workload(value, 1.0).is_some()
            }
            "--kd" => value.parse().map(|v| parsed.kd = v).is_ok(),
            "--kim" => value.parse().map(|v| parsed.kim = v).is_ok(),
            "--density" => Density::from_name(value)
                .map(|d| parsed.density = d)
                .is_some(),
            "--charges" => match value.as_str() {
                "real" => Some(ChargeKind::Real),
                "complex" => Some(ChargeKind::Complex),
                _ => None,
            }
            .map(|c| parsed.charges = c)
            .is_some(),
            "--eps" => value
                .split(',')
                .map(|v| v.parse::<f64>().ok().filter(|&e| e > 0.0 && e < 1.0))
                .collect::<Option<Vec<f64>>>()
                .map(|e| parsed.eps = e)
                .is_some(),
            "--m2l" => value
                .split(',')
                .map(|v| match v {
                    "rotation" => Some(Setting::Rotation),
                    "default" => Some(Setting::Default),
                    "defaults" => Some(Setting::Defaults),
                    _ => None,
                })
                .collect::<Option<Vec<Setting>>>()
                .map(|d| parsed.m2l = d)
                .is_some(),
            "--study" => value
                .split(',')
                .map(|v| match v {
                    "single" => Some((false, false)),
                    "levels" => Some((true, false)),
                    "m2l" => Some((false, true)),
                    _ => None,
                })
                .collect::<Option<Vec<(bool, bool)>>>()
                .map(|studies| {
                    parsed.levels = studies.iter().any(|s| s.0);
                    parsed.m2l_study = studies.iter().any(|s| s.1);
                })
                .is_some(),
            "--refine" => value.parse().map(|v| parsed.refine = v).is_ok(),
            "--refine-cells" => value
                .split(';')
                .map(|c| {
                    let v: Option<Vec<u64>> = c.split(',').map(|x| x.parse().ok()).collect();
                    v.filter(|v| v.len() == 3 && v.iter().all(|&x| x < 1 << 21))
                        .map(|v| [v[0], v[1], v[2]])
                })
                .collect::<Option<Vec<[u64; 3]>>>()
                .map(|cells| parsed.refine_cells = cells)
                .is_some(),
            "--reference-cache" => {
                parsed.reference_cache = Some(value.clone());
                true
            }
            "--raise-by" => value
                .split(',')
                .map(|v| v.parse::<usize>().ok().filter(|&d| d > 0))
                .collect::<Option<Vec<usize>>>()
                .map(|d| parsed.raise_by = d)
                .is_some(),
            "--raise-levels" if value == "v" => {
                parsed.raise_levels = Some(Vec::new());
                true
            }
            "--raise-levels" => value
                .split(',')
                .map(|v| v.parse::<usize>().ok().filter(|&l| l >= 2))
                .collect::<Option<Vec<usize>>>()
                .map(|l| parsed.raise_levels = Some(l))
                .is_some(),
            "--margin" => value
                .parse::<f64>()
                .ok()
                .filter(|&m| m >= 0.0)
                .map(|m| parsed.margin = Some(m))
                .is_some(),
            "--offset" => {
                let parts: Option<Vec<f64>> = value.split(',').map(|v| v.parse().ok()).collect();
                parts
                    .filter(|p| p.len() == 3)
                    .map(|p| parsed.offset = [p[0], p[1], p[2]])
                    .is_some()
            }
            "--leaf-size" => {
                parsed.leaf_size = Some(value.clone());
                value == "uniform" || value.parse::<usize>().is_ok_and(|n| n > 0)
            }
            "--raise-factor" => value
                .parse::<f64>()
                .ok()
                .filter(|&f| f > 0.0 && f < 1.0)
                .map(|f| parsed.raise_factor = f)
                .is_some(),
            "--worst" => value.parse().map(|v| parsed.worst = v).is_ok(),
            "--per-box" => value.parse().map(|v| parsed.per_box = v).is_ok(),
            "--rule" => match value.as_str() {
                "calibrated" => Some(None),
                "ch3-3" => Some(Some(DegreeRule::CH3_3)),
                "ch4-3c2" => Some(Some(DegreeRule::CH4_3C2)),
                "ch4-3c3" => Some(Some(DegreeRule::CH4_3C3)),
                constants => match constants
                    .split(',')
                    .map(|v| v.parse::<f64>().ok())
                    .collect::<Option<Vec<f64>>>()
                    .as_deref()
                {
                    Some(&[gamma, alpha, beta, c, a, b, d]) => Some(Some(DegreeRule {
                        gamma,
                        alpha,
                        beta,
                        c,
                        a,
                        b,
                        d,
                        cap: None,
                        margin: 0,
                    })),
                    Some(&[gamma, alpha, beta, c, a, b, d, cap]) => Some(Some(DegreeRule {
                        gamma,
                        alpha,
                        beta,
                        c,
                        a,
                        b,
                        d,
                        cap: Some(cap),
                        margin: 0,
                    })),
                    Some(&[gamma, alpha, beta, c, a, b, d, cap, margin])
                        if margin >= 0.0 && margin.fract() == 0.0 =>
                    {
                        Some(Some(DegreeRule {
                            gamma,
                            alpha,
                            beta,
                            c,
                            a,
                            b,
                            d,
                            cap: Some(cap),
                            margin: margin as usize,
                        }))
                    }
                    _ => None,
                },
            }
            .map(|r| parsed.degree_rule = r)
            .is_some(),
            "--full-max" => value.parse().map(|v| parsed.full_max = v).is_ok(),
            "--threads" => value.parse().map(|v| parsed.threads = v).is_ok(),
            _ => false,
        };
        if !ok {
            usage();
        }
        rest = tail;
    }
    if !rest.is_empty() || parsed.threads == 0 {
        usage();
    }
    if parsed.density == Density::PlaneWave && parsed.charges == ChargeKind::Real {
        eprintln!("helmholtz_boxes: a plane wave's trace needs complex charges");
        std::process::exit(2);
    }
    parsed
}

/// The degrees as a short string: runs of equal degrees as `p×count`.
fn degrees_string(degrees: &[usize]) -> String {
    let mut runs: Vec<(usize, usize)> = Vec::new();
    for &p in degrees {
        match runs.last_mut() {
            Some((q, count)) if *q == p => *count += 1,
            _ => runs.push((p, 1)),
        }
    }
    runs.iter()
        .map(|&(p, count)| {
            if count == 1 {
                p.to_string()
            } else {
                format!("{p}×{count}")
            }
        })
        .collect::<Vec<_>>()
        .join(", ")
}

/// One measured run.
struct Measured {
    label: String,
    eps: f64,
    run: Run,
    errors: Errors,
    boxes: Vec<BoxError>,
}

/// Runs `case` and measures it against `reference` on `grid` (the first run's deepest level,
/// or this run's when `None`).
#[allow(clippy::too_many_arguments)]
fn measure<C: Communicator>(
    label: String,
    eps: f64,
    case: &Case,
    problem: &Problem,
    oracle: &Oracle,
    reference: &Reference,
    grid: Option<&Grid>,
    threads: usize,
    world: &C,
) -> Measured {
    let start = Instant::now();
    let run = run_output(case, problem, oracle, threads, world).unwrap_or_else(|e| {
        eprintln!("helmholtz_boxes: {label}: {e}");
        std::process::exit(1);
    });
    let fmm = start.elapsed().as_secs_f64();
    let errors = Errors::of(&run.output, 0, reference);
    let own = Grid {
        domain: run.domain,
        level: run.row.levels - 1,
    };
    let boxes = errors.boxes(reference, &problem.points, grid.unwrap_or(&own));
    // The run as it finishes (the tables print after the last run): §9.1's errors and those
    // at the reference's targets, the build, the evaluation and its M2L by level, the rest
    // of `run_output` (the samples' errors) and the errors at the reference's targets.
    let row = &run.row;
    let m2l: Vec<String> = row
        .m2l_seconds
        .iter()
        .enumerate()
        .filter(|&(_, &s)| s > 0.0)
        .map(|(l, s)| format!("{l}: {s:.1}"))
        .collect();
    let all = errors.l2(reference);
    eprintln!(
        "helmholtz_boxes: {label}: p {}; §9.1 φ, ∇φ {:.2e}, {:.2e}, all {:.2e}, {:.2e}; \
         build {:.1} s, evaluate {:.1} s (M2L by level {}), the rest of the run {:.1} s, the \
         errors at {} targets {:.1} s",
        degrees_string(&row.degrees),
        row.errors[0][0],
        row.errors[0][1],
        all[0],
        all[1],
        row.build,
        row.evaluate,
        m2l.join(", "),
        fmm - row.build - row.evaluate,
        reference.targets.len(),
        start.elapsed().as_secs_f64() - fmm,
    );
    Measured {
        label,
        eps,
        run,
        errors,
        boxes,
    }
}

/// The row of one run in the runs table.
fn run_row(m: &Measured, reference: &Reference) {
    let (row, eps) = (&m.run.row, m.eps);
    let l2 = m.errors.l2(reference);
    let worst = m.boxes.iter().fold([0.0f64; 2], |w, b| {
        [w[0].max(b.error[0]), w[1].max(b.error[1])]
    });
    let above = |f: usize| m.boxes.iter().filter(|b| b.error[f] > eps).count();
    println!(
        "| {} | {} | {} | {:.2e}, {:.2e} | {:.2e}, {:.2e} | {:.2}, {:.2} | {:.0}, {:.0} | {}, {} of {} | {:.2}, {:.2} |",
        m.label,
        degrees_string(&row.degrees),
        row.leaves,
        row.errors[0][0],
        row.errors[0][1],
        l2[0],
        l2[1],
        worst[0] / eps,
        worst[1] / eps,
        worst[0] / l2[0].max(f64::MIN_POSITIVE),
        worst[1] / l2[1].max(f64::MIN_POSITIVE),
        above(0),
        above(1),
        m.boxes.len(),
        row.build,
        row.evaluate,
    );
}

/// Prints the worst boxes of field `field` of `base` with their place, and with the levels'
/// parts at their worst targets when `raised` is given, and the M2L runs' parts (Helmholtz
/// CH4.4b) when `m2l` is.
#[allow(clippy::too_many_arguments)]
fn worst_table(
    base: &Measured,
    raised: Option<&[Option<Errors>]>,
    m2l: Option<(&[Option<Errors>], &[&str])>,
    reference: &Reference,
    problem: &Problem,
    grid: &Grid,
    field: usize,
    count: usize,
    eps: f64,
) {
    let name = ["φ", "∇φ"][field];
    let deepest = grid.level;
    let cells = 1u64 << deepest;
    println!();
    println!(
        "Worst boxes of {name} (level {deepest}): the cell, its centre in the domain's unit \
         cube, the error over ε, the targets, whether a face lies on the domain's boundary, the \
         worst target's distance |v| from its box's centre on levels 2–{deepest} over the \
         half-width (√3 at a corner){}{}.",
        if raised.is_some() {
            ", and each level's part of its error (at least 0.1 of it)"
        } else {
            ""
        },
        if m2l.is_some() {
            ", and each M2L run's part (at least 0.01 of it)"
        } else {
            ""
        }
    );
    println!();
    println!(
        "| # | cell | centre | error / ε | targets | boundary | |v| / r, levels 2–{deepest} | level parts | M2L parts |"
    );
    println!("| ---: | --- | --- | ---: | ---: | --- | --- | --- | --- |");
    for (rank, b) in worst_boxes(&base.boxes, field, count).iter().enumerate() {
        let c = cell(b.key);
        let centre: Vec<String> = c
            .iter()
            .map(|&i| format!("{:.3}", (i as f64 + 0.5) / cells as f64))
            .collect();
        let boundary = c.iter().any(|&i| i == 0 || i == cells - 1);
        let j = b.worst[field];
        let x = problem.points[reference.targets[j]];
        let profile: Vec<String> = offsets(x, grid, deepest)
            .iter()
            .skip(2)
            .map(|o| format!("{:.2}", o.iter().map(|v| v * v).sum::<f64>().sqrt()))
            .collect();
        let parts = raised.map_or(String::new(), |raised| {
            level_parts_at(&base.errors, raised, j)
                .iter()
                .enumerate()
                .filter(|(_, p)| p[field] >= 0.1)
                .map(|(l, p)| format!("{l}: {:.2}", p[field]))
                .collect::<Vec<_>>()
                .join(", ")
        });
        let m2l_parts = m2l.map_or(String::new(), |(runs, labels)| {
            level_parts_at(&base.errors, runs, j)
                .iter()
                .zip(labels)
                .filter(|(p, _)| p[field] >= 0.01)
                .map(|(p, label)| format!("{label}: {:.3}", p[field]))
                .collect::<Vec<_>>()
                .join("; ")
        });
        println!(
            "| {} | ({}, {}, {}) | ({}) | {:.2} | {} | {} | {} | {} | {} |",
            rank + 1,
            c[0],
            c[1],
            c[2],
            centre.join(", "),
            b.error[field] / eps,
            b.targets,
            if boundary { "yes" } else { "no" },
            profile.join(" "),
            parts,
            m2l_parts,
        );
    }
}

/// The errors by the number of levels 2..=deepest on which a target lies within 0.1 of a
/// face of its box.
fn faces_table(base: &Measured, reference: &Reference, problem: &Problem, grid: &Grid, eps: f64) {
    let deepest = grid.level;
    let rms = reference.rms();
    let mut bins = vec![(0usize, [0.0f64; 2], [0.0f64; 2]); deepest.max(1)];
    for (j, &i) in reference.targets.iter().enumerate() {
        let closeness = face_closeness(&offsets(problem.points[i], grid, deepest));
        let faces = closeness.iter().skip(2).filter(|&&v| v >= 0.9).count();
        let e = base.errors.norms(j);
        let bin = &mut bins[faces.min(deepest.saturating_sub(1))];
        bin.0 += 1;
        for (f, e) in e.into_iter().enumerate() {
            bin.1[f] += e * e;
            bin.2[f] = bin.2[f].max(e);
        }
    }
    println!();
    println!(
        "Errors by the levels (of 2–{deepest}) on which the target lies within 0.1 of a face of \
         its box: the RMS error over the reference's RMS, and the largest, over ε."
    );
    println!();
    println!("| levels on a face | targets | RMS φ, ∇φ / ε | largest φ, ∇φ / ε |");
    println!("| ---: | ---: | --- | --- |");
    for (faces, (n, sums, max)) in bins.iter().enumerate() {
        if *n == 0 {
            continue;
        }
        let r = |f: usize| (sums[f] / *n as f64).sqrt() / rms[f] / eps;
        println!(
            "| {faces} | {n} | {:.3}, {:.3} | {:.2}, {:.2} |",
            r(0),
            r(1),
            max[0] / rms[0] / eps,
            max[1] / rms[1] / eps,
        );
    }
}

/// What ran on level `l` of `row`: the diagonal form with its classes, the exponential M2L,
/// the tables of Dense, Classes or Batched (Auto's tables are Batched's), RCR, or nothing (no
/// V pairs); as `helmholtz_fmm3d`'s log.
fn level_m2l(row: &Row, l: usize) -> String {
    let classes = row.diagonal_classes.get(l).copied().unwrap_or(0);
    if row.v_pairs.get(l).copied().unwrap_or(0) == 0 {
        "–".into()
    } else if classes > 0 {
        format!("diagonal {classes}/34")
    } else if row.plane_wave.get(l).is_some_and(Option::is_some) {
        "waves".into()
    } else if row.dense_levels.get(l).copied().unwrap_or(false) {
        match row.m2l {
            HelmholtzM2l::Auto => "Batched".into(),
            m2l => format!("{m2l:?}"),
        }
    } else {
        "RCR".into()
    }
}

/// A sampled reference stored as little-endian values after one header line that
/// names the case, the grid and the sample (`--reference-cache`): per target its position (u64),
/// φ (2 f64) and ∇φ (6 f64).
fn write_reference(path: &std::path::Path, header: &str, reference: &Reference) {
    let mut bytes = format!("{header}\n").into_bytes();
    for (j, &i) in reference.targets.iter().enumerate() {
        bytes.extend_from_slice(&(i as u64).to_le_bytes());
        let p = reference.potential[j];
        let g = reference.gradient[j];
        for v in [
            p.re, p.im, g[0].re, g[0].im, g[1].re, g[1].im, g[2].re, g[2].im,
        ] {
            bytes.extend_from_slice(&v.to_le_bytes());
        }
    }
    let partial = path.with_extension("partial");
    if let Err(e) = std::fs::write(&partial, &bytes).and_then(|()| std::fs::rename(&partial, path))
    {
        eprintln!(
            "helmholtz_boxes: the reference cache {}: {e}",
            path.display()
        );
    }
}

/// The reference of [`write_reference`] at `path` if its header is `header`, else `None`.
fn read_reference(path: &std::path::Path, header: &str) -> Option<Reference> {
    let bytes = std::fs::read(path).ok()?;
    let line = bytes.iter().position(|&b| b == b'\n')?;
    if &bytes[..line] != header.as_bytes() {
        return None;
    }
    let body = &bytes[line + 1..];
    let record = 8 * 9;
    if body.len() % record != 0 {
        return None;
    }
    let f = |at: &[u8]| f64::from_le_bytes(at.try_into().expect("8 bytes"));
    let mut reference = Reference {
        targets: Vec::new(),
        potential: Vec::new(),
        gradient: Vec::new(),
    };
    for r in body.chunks_exact(record) {
        let v: Vec<f64> = (1..9).map(|i| f(&r[8 * i..8 * i + 8])).collect();
        reference
            .targets
            .push(u64::from_le_bytes(r[..8].try_into().expect("8 bytes")) as usize);
        reference.potential.push(Complex::new(v[0], v[1]));
        reference.gradient.push([
            Complex::new(v[2], v[3]),
            Complex::new(v[4], v[5]),
            Complex::new(v[6], v[7]),
        ]);
    }
    Some(reference)
}

/// Level l's wavenumber in its units, k r_l, from the run's circumradius √3 r_l.
fn k_hat(row: &Row, k: Complex<f64>, l: usize) -> Complex<f64> {
    let r = row.circumradii[l] / 3f64.sqrt();
    Complex::new(k.re * r, k.im * r)
}

/// The diagonal classes of level l as run, by decade below the rule's threshold: the classes
/// of each decade (0: E/(θ ε) in (0.1, 1], 1: (0.01, 0.1], 2: the rest that pass), the classes
/// that fail, and the largest ratio among those that pass. `None` where the level did not run
/// the form.
fn class_decades(
    row: &Row,
    k: Complex<f64>,
    eps: f64,
    l: usize,
) -> Option<([Vec<usize>; 3], usize, f64)> {
    if row.diagonal_classes.get(l).copied().unwrap_or(0) == 0 {
        return None;
    }
    let ratios = class_ratios(
        row.degrees[l],
        k_hat(row, k, l),
        eps,
        row.diagonal_build_kind,
    );
    let mut decades: [Vec<usize>; 3] = Default::default();
    let mut failing = 0;
    let mut largest = 0.0f64;
    for (c, &ratio) in ratios.iter().enumerate() {
        match threshold_decade(ratio) {
            None => failing += 1,
            Some(d) => {
                decades[d.min(2)].push(c);
                largest = largest.max(ratio);
            }
        }
    }
    let passing: usize = decades.iter().map(Vec::len).sum();
    if passing != row.diagonal_classes[l] {
        eprintln!(
            "helmholtz_boxes: level {l}: {} classes ran, the bound passes {passing}",
            row.diagonal_classes[l]
        );
    }
    Some((decades, failing, largest))
}

/// The worst box of each field over the run's, and the ℓ2 errors over the reference's
/// targets, of `run` over `base`'s.
fn ratios_to(run: &Measured, base: &Measured, reference: &Reference) -> ([f64; 2], [f64; 2]) {
    let worst = |m: &Measured| {
        m.boxes.iter().fold([0.0f64; 2], |w, b| {
            [w[0].max(b.error[0]), w[1].max(b.error[1])]
        })
    };
    let (w, wb) = (worst(run), worst(base));
    let (l, lb) = (run.errors.l2(reference), base.errors.l2(reference));
    let r = |a: f64, b: f64| if b > 0.0 { a / b } else { 0.0 };
    (
        [r(w[0], wb[0]), r(w[1], wb[1])],
        [r(l[0], lb[0]), r(l[1], lb[1])],
    )
}

fn main() {
    let args = arguments();
    let (universe, _) =
        mpi::initialize_with_threading(Threading::Funneled).expect("MPI initialises");
    let world = universe.world();
    if world.size() != 1 {
        eprintln!("helmholtz_boxes: one rank only");
        std::process::exit(2);
    }
    let workload = workload(&args.workload, args.kd).expect("checked");
    let (eps, setting) = (args.eps[0], args.m2l[0]);
    let case_of = |setting: Setting, degrees: Degrees| {
        let mut case = match setting {
            Setting::Rotation => Case::new(workload, (args.kd, args.kim), degrees),
            Setting::Default => Case {
                m2l: HelmholtzM2l::Auto,
                ..Case::new(workload, (args.kd, args.kim), degrees)
            },
            Setting::Defaults => Case::defaults(workload, (args.kd, args.kim), degrees),
        };
        case.charges = args.charges;
        case.density = args.density;
        case.degree_rule = args.degree_rule;
        case
    };
    let mut case = case_of(setting, Degrees::Tolerance(eps));
    let setup = Instant::now();
    let problem = Problem::new(&case);
    let n = problem.len();
    if let Some(margin) = args.margin {
        let o = args.offset;
        case.domain = Some([
            -margin + o[0],
            -margin + o[1],
            -margin + o[2],
            1.0 + margin + o[0],
            1.0 + margin + o[1],
            1.0 + margin + o[2],
        ]);
    } else if args.offset != [0.0; 3] {
        eprintln!("helmholtz_boxes: --offset needs --margin");
        std::process::exit(2);
    }
    case.max_points_per_leaf = args.leaf_size.as_deref().map(|s| match s {
        "uniform" => n + 1,
        n => n.parse().expect("checked"),
    });
    let k = complex_k(&case);
    let oracle = Oracle::sharded(k, &problem, true, args.threads, &world);
    eprintln!(
        "helmholtz_boxes: the problem and the samples' oracle, {:.1} s",
        setup.elapsed().as_secs_f64()
    );

    println!(
        "# The box study: {} at kD = {}{}, ε = {:e}, {} density, {}",
        workload.name(),
        args.kd,
        if args.kim != 0.0 {
            format!(" + {}i", args.kim)
        } else {
            String::new()
        },
        eps,
        args.density.name(),
        setting.name(),
    );
    println!();
    println!(
        "Machine: {} ({} cores), 1 rank × {} threads. N = {n}. Domain: {}. Leaf size: {}.",
        cpu_model(),
        cores(),
        args.threads,
        case.domain
            .map_or("the points' bounding cube".into(), |d| format!("{d:?}")),
        args.leaf_size.as_deref().unwrap_or("the library's"),
    );

    // The first run fixes the grid of the box measure, and so the per-box sample.
    let first = measure(
        format!("ε = {eps:e}, {}", setting.name()),
        eps,
        &case,
        &problem,
        &oracle,
        &Reference::at(k, &problem, 0, problem.sample.clone(), args.threads),
        None,
        args.threads,
        &world,
    );
    let grid = Grid {
        domain: first.run.domain,
        level: first.run.row.levels - 1,
    };
    let start = Instant::now();
    // The sampled reference depends on the case, the grid and the sample alone, not on ε: a
    // launch of another ε reads it from the cache (Helmholtz CH4.4b).
    let header = format!(
        "tetraleaf helmholtz_boxes reference 1: {} kD {} {} density {} charges {:?} N {n} \
         domain {:?} level {} per box {} full max {} seed {SEED}",
        workload.name(),
        args.kd,
        args.kim,
        args.density.name(),
        args.charges,
        grid.domain,
        grid.level,
        args.per_box,
        args.full_max,
    );
    let cache = args.reference_cache.as_ref().map(|dir| {
        let name: String = header
            .chars()
            .map(|c| if c.is_ascii_alphanumeric() { c } else { '_' })
            .collect();
        let mut hash = 0xcbf2_9ce4_8422_2325u64;
        for b in name.bytes() {
            hash = (hash ^ u64::from(b)).wrapping_mul(0x100_0000_01b3);
        }
        std::path::Path::new(dir).join(format!("reference-{hash:016x}.bin"))
    });
    let cached = cache
        .as_ref()
        .and_then(|path| read_reference(path, &header));
    let from_cache = cached.is_some();
    let mut reference = cached.unwrap_or_else(|| {
        let reference = if n <= args.full_max {
            Reference::every(k, &problem, 0, args.threads)
        } else {
            let targets = box_sample(&problem.points, &problem.sample, &grid, args.per_box, SEED);
            Reference::at(k, &problem, 0, targets, args.threads)
        };
        if let Some(path) = &cache {
            write_reference(path, &header, &reference);
        }
        reference
    });
    let mut errors = Errors::of(&first.run.output, 0, &reference);
    let mut boxes = errors.boxes(&reference, &problem.points, &grid);
    let sampled = reference.targets.len();
    // Every point of the worst boxes of the sampled measure (Helmholtz CH4.4b), so that their
    // maxima are exact.
    let mut refined = 0;
    if (args.refine > 0 || args.refine_faces || !args.refine_cells.is_empty()) && sampled < n {
        let mut keys: std::collections::HashSet<u64> = [0, 1]
            .iter()
            .flat_map(|&f| worst_boxes(&boxes, f, args.refine))
            .map(|b| b.key)
            .collect();
        // The named boxes (`Grid::box_of`'s keys: i | j << 21 | k << 42).
        let named: std::collections::HashSet<u64> = args
            .refine_cells
            .iter()
            .map(|c| c[0] | (c[1] << 21) | (c[2] << 42))
            .collect();
        if !named.is_empty() {
            for &x in &problem.points {
                let key = grid.box_of(x);
                if named.contains(&key) {
                    keys.insert(key);
                }
            }
        }
        if args.refine_faces {
            // The boxes whose closure holds the centre of a face of the domain: the outermost
            // layer along the face's normal, the two middle cells along the other axes.
            let cells = 1u64 << grid.level;
            let middle = |i: u64| cells >= 2 && (i == cells / 2 - 1 || i == cells / 2);
            for &x in &problem.points {
                let key = grid.box_of(x);
                let c = cell(key);
                let at_face = (0..3).any(|a| {
                    (c[a] == 0 || c[a] == cells - 1)
                        && (0..3).filter(|&b| b != a).all(|b| middle(c[b]))
                });
                if at_face {
                    keys.insert(key);
                }
            }
        }
        let more: Vec<usize> = (0..n)
            .filter(|&i| keys.contains(&grid.box_of(problem.points[i])))
            .collect();
        reference = reference.extend(k, &problem, 0, &more, args.threads);
        refined = keys.len();
        errors = Errors::of(&first.run.output, 0, &reference);
        boxes = errors.boxes(&reference, &problem.points, &grid);
    }
    println!(
        "Reference: the direct sum at {} of {n} targets ({}{}), {:.1} s{}.",
        reference.targets.len(),
        if sampled == n {
            "every point".to_owned()
        } else {
            format!(
                "§9.1's samples and {} per box of level {}",
                args.per_box, grid.level
            )
        },
        if refined > 0 {
            format!(
                ", and every point of the {refined} worst boxes of φ and ∇φ by that sample{}, {} \
                 targets more",
                match (args.refine_faces, args.refine_cells.is_empty()) {
                    (true, true) => " and of the boxes at the centres of the domain's faces",
                    (true, false) => {
                        " and of the boxes at the centres of the domain's faces and those named"
                    }
                    (false, false) => " and of the boxes named",
                    (false, true) => "",
                },
                reference.targets.len() - sampled
            )
        } else {
            String::new()
        },
        start.elapsed().as_secs_f64(),
        if from_cache {
            " (the sample read from the cache)"
        } else {
            ""
        }
    );
    let base = Measured {
        boxes,
        errors,
        ..first
    };
    let row = base.run.row.clone();
    let p_run = row.degrees.clone();
    let measure_case = |label: String, case: &Case| {
        measure(
            label,
            eps,
            case,
            &problem,
            &oracle,
            &reference,
            Some(&grid),
            args.threads,
            &world,
        )
    };
    let with_overrides = |overrides: Vec<M2lOverride>| Case {
        m2l_overrides: overrides,
        ..case.clone()
    };

    let mut runs = Vec::new();
    for &other_eps in &args.eps {
        for &other in &args.m2l {
            if (other_eps, other) == (eps, setting) {
                continue;
            }
            let other_case = Case {
                domain: case.domain,
                max_points_per_leaf: case.max_points_per_leaf,
                ..case_of(other, Degrees::Tolerance(other_eps))
            };
            runs.push(measure(
                format!("ε = {other_eps:e}, {}", other.name()),
                other_eps,
                &other_case,
                &problem,
                &oracle,
                &reference,
                Some(&grid),
                args.threads,
                &world,
            ));
        }
    }
    let mut raised: Option<Vec<Option<Errors>>> = None;
    let mut raised_to = vec![0usize; p_run.len()];
    if args.levels {
        let wavenumber = Wavenumber::new(args.kd, args.kim).expect("a valid wavenumber");
        let rule = args.degree_rule.unwrap_or(DegreeRule::CALIBRATED);
        let tight = Degrees::Tolerance(eps * args.raise_factor)
            .resolve_with(&rule, &wavenumber, &base.run.domain, row.levels)
            .expect("a valid tolerance");
        let to = |l: usize| tight[l].max(p_run[l] + 1);
        // Under the library's defaults one level raised and on RCR by an override, the
        // tolerance and every other level's choice kept (Helmholtz CH4.4b); else CH4.0b's named
        // degrees (which change the waves' precision and take the form off).
        let raise_case = |levels: &[usize]| {
            if setting == Setting::Defaults {
                with_overrides(
                    levels
                        .iter()
                        .map(|&l| M2lOverride {
                            level: l,
                            degree: Some(to(l)),
                            m2l: OverrideM2l::Rcr,
                        })
                        .collect(),
                )
            } else {
                let mut named = p_run.clone();
                for &l in levels {
                    named[l] = to(l);
                }
                Case {
                    degrees: Degrees::PerLevel(named),
                    ..case.clone()
                }
            }
        };
        let mut parts = vec![None; p_run.len()];
        let levels: Vec<usize> = (2..p_run.len())
            .filter(|&l| args.all_levels || row.v_pairs[l] > 0)
            .collect();
        for &l in &levels {
            raised_to[l] = to(l);
            let m = measure_case(format!("level {l} at {}", to(l)), &raise_case(&[l]));
            parts[l] = Some(m.errors.clone());
            runs.push(m);
        }
        runs.push(measure_case(
            "every level raised".into(),
            &raise_case(&levels),
        ));
        raised = Some(parts);
    }

    // The M2L study (Helmholtz CH4.4b): each level whose M2L is not RCR on RCR at its degree
    // as run; on the form's levels the classes by decade below the threshold on RCR.
    let mut m2l_runs: Vec<Measured> = Vec::new();
    if args.m2l_study {
        for l in 2..p_run.len() {
            let ran = level_m2l(&row, l);
            if ran == "–" || ran == "RCR" {
                continue;
            }
            m2l_runs.push(measure_case(
                format!("level {l} RCR (was {ran})"),
                &with_overrides(vec![M2lOverride {
                    level: l,
                    degree: None,
                    m2l: OverrideM2l::Rcr,
                }]),
            ));
            let Some((decades, _, _)) = class_decades(&row, k, eps, l) else {
                continue;
            };
            for (d, classes) in decades.iter().enumerate().take(2) {
                if classes.is_empty() {
                    continue;
                }
                m2l_runs.push(measure_case(
                    format!("level {l} decade {d} {classes:?} RCR"),
                    &with_overrides(vec![M2lOverride {
                        level: l,
                        degree: None,
                        m2l: OverrideM2l::RcrClasses(classes.clone()),
                    }]),
                ));
            }
            if args.class_runs {
                for &c in &decades[0] {
                    m2l_runs.push(measure_case(
                        format!("level {l} class {c} RCR"),
                        &with_overrides(vec![M2lOverride {
                            level: l,
                            degree: None,
                            m2l: OverrideM2l::RcrClasses(vec![c]),
                        }]),
                    ));
                }
            }
        }
    }

    // The degree raised by d under the setting's M2L (Helmholtz CH4.4b: a remedy's FMM
    // evidence), on the levels named or the one with the largest part at the worst ∇φ box.
    if !args.raise_by.is_empty() {
        let levels = args.raise_levels.clone().unwrap_or_else(|| {
            let Some(raised) = &raised else {
                eprintln!("helmholtz_boxes: --raise-by needs --raise-levels or --study levels");
                std::process::exit(2);
            };
            let worst = worst_boxes(&base.boxes, 1, 1);
            let j = worst.first().map_or(0, |b| b.worst[1]);
            let parts = level_parts_at(&base.errors, raised, j);
            let level = (0..parts.len())
                .max_by(|&a, &b| parts[a][1].total_cmp(&parts[b][1]))
                .unwrap_or(2);
            vec![level]
        });
        let levels: Vec<usize> = if levels.is_empty() {
            (2..p_run.len()).filter(|&l| row.v_pairs[l] > 0).collect()
        } else {
            levels.into_iter().filter(|&l| l < p_run.len()).collect()
        };
        if args.raise_together {
            for &d in &args.raise_by {
                runs.push(measure_case(
                    format!("levels {levels:?} at +{d}, the setting's M2L"),
                    &with_overrides(
                        levels
                            .iter()
                            .map(|&l| M2lOverride {
                                level: l,
                                degree: Some(p_run[l] + d),
                                m2l: OverrideM2l::Setting,
                            })
                            .collect(),
                    ),
                ));
            }
        }
        for l in levels.into_iter().filter(|_| !args.raise_together) {
            for &d in &args.raise_by {
                let degree = p_run[l] + d;
                runs.push(measure_case(
                    format!("level {l} at {degree} (+{d}), the setting's M2L"),
                    &with_overrides(vec![M2lOverride {
                        level: l,
                        degree: Some(degree),
                        m2l: OverrideM2l::Setting,
                    }]),
                ));
            }
        }
    }

    println!();
    println!(
        "Runs: the degrees by level; the leaves; §9.1's relative ℓ2 errors at the 1,000 samples; \
         the relative ℓ2 errors over the reference's targets; the worst box of level {} over ε \
         and over the ℓ2 error; the boxes whose error exceeds ε; the build and the evaluation \
         in seconds.",
        grid.level
    );
    println!();
    println!(
        "| run | p by level | leaves | §9.1 φ, ∇φ | all φ, ∇φ | worst box / ε | worst / ℓ2 | boxes > ε φ, ∇φ | build, evaluate s |"
    );
    println!("| --- | --- | ---: | --- | --- | --- | --- | --- | --- |");
    run_row(&base, &reference);
    for m in runs.iter().chain(&m2l_runs) {
        run_row(m, &reference);
    }

    // The levels as run (Helmholtz CH4.4b).
    println!();
    println!(
        "Levels as run: k a (Re k √3 r_l), the degree, the V pairs, the M2L; on the form's levels \
         the classes that pass the bound by decade of E/(θ ε) below the threshold ((0.1, 1], \
         (0.01, 0.1], below), those that fail (on RCR), and the largest E/(θ ε) that passes."
    );
    println!();
    println!(
        "| level | k a | p | V pairs | M2L | classes (0.1, 1] | (0.01, 0.1] | below | failing | largest E/(θ ε) |"
    );
    println!("| ---: | ---: | ---: | ---: | --- | --- | ---: | ---: | ---: | ---: |");
    for (l, p) in p_run.iter().enumerate().skip(2) {
        let (decades, failing, largest) = class_decades(&row, k, eps, l).map_or(
            (String::new(), String::new(), String::new()),
            |(d, f, w)| {
                (
                    format!("{:?} | {} | {}", d[0], d[1].len(), d[2].len()),
                    f.to_string(),
                    format!("{w:.3}"),
                )
            },
        );
        println!(
            "| {l} | {:.2} | {} | {} | {} | {} | {failing} | {largest} |",
            k.re * row.circumradii[l],
            p,
            row.v_pairs[l],
            level_m2l(&row, l),
            if decades.is_empty() {
                " |  | ".to_owned()
            } else {
                decades
            },
        );
    }

    if let Some(raised) = &raised {
        let shares = level_shares(&base.errors, raised);
        println!();
        println!(
            "Each level's part of the error (the run minus the run with that level raised{}): its \
             share of Σ|e|² over the reference's targets, and the rest.",
            if setting == Setting::Defaults {
                " and on RCR, every other level as the defaults choose"
            } else {
                ""
            }
        );
        println!();
        println!("| level | p | M2L | raised to | V pairs | share φ | share ∇φ |");
        println!("| ---: | ---: | --- | ---: | ---: | ---: | ---: |");
        for (l, share) in shares.iter().enumerate().take(raised.len()) {
            if raised[l].is_none() {
                continue;
            }
            println!(
                "| {l} | {} | {} | {} | {} | {:.3} | {:.3} |",
                p_run[l],
                level_m2l(&row, l),
                raised_to[l],
                row.v_pairs[l],
                share[0],
                share[1]
            );
        }
        let rest = shares.last().expect("the rest");
        println!("| rest | | | | | {:.3} | {:.3} |", rest[0], rest[1]);
    }

    if !m2l_runs.is_empty() {
        let parts: Vec<Option<Errors>> = m2l_runs.iter().map(|m| Some(m.errors.clone())).collect();
        let shares = level_shares(&base.errors, &parts);
        let worst = worst_boxes(&base.boxes, 1, 1);
        let j = worst.first().map_or(0, |b| b.worst[1]);
        let at = level_parts_at(&base.errors, &parts, j);
        println!();
        println!(
            "The M2L's parts: each run's difference from the first run, its share of Σ|e|² over \
             the reference's targets; the run's worst box and ℓ2 error over the first run's; the \
             difference at the worst ∇φ box's worst target over the error there (φ, ∇φ)."
        );
        println!();
        println!(
            "| run | share φ | share ∇φ | worst box φ, ∇φ ÷ first | ℓ2 φ, ∇φ ÷ first | at the worst ∇φ target |"
        );
        println!("| --- | ---: | ---: | --- | --- | --- |");
        for (i, m) in m2l_runs.iter().enumerate() {
            let (w, l) = ratios_to(m, &base, &reference);
            println!(
                "| {} | {:.4} | {:.4} | {:.3}, {:.3} | {:.3}, {:.3} | {:.3}, {:.3} |",
                m.label, shares[i][0], shares[i][1], w[0], w[1], l[0], l[1], at[i][0], at[i][1]
            );
        }
    }

    let m2l_parts: Vec<Option<Errors>> = m2l_runs.iter().map(|m| Some(m.errors.clone())).collect();
    let m2l_labels: Vec<&str> = m2l_runs.iter().map(|m| m.label.as_str()).collect();
    for field in [1, 0] {
        worst_table(
            &base,
            raised.as_deref(),
            (!m2l_parts.is_empty()).then_some((&m2l_parts[..], &m2l_labels[..])),
            &reference,
            &problem,
            &grid,
            field,
            args.worst,
            eps,
        );
    }
    faces_table(&base, &reference, &problem, &grid, eps);
    if let Some(raised) = &raised {
        distance_table(&base, raised, &reference, &problem, &grid, eps);
    }
}

/// The edges of the bins of a target's distance from its box's centre over the half-width
/// (√3 at a corner).
const DISTANCE_BINS: [f64; 6] = [0.0, 0.5, 0.75, 1.0, 1.25, 1.5];

/// Each level's part of the error by the target's distance from the centre of its box on
/// that level, over the half-width: the RMS of the part over the reference's RMS, over ε.
fn distance_table(
    base: &Measured,
    raised: &[Option<Errors>],
    reference: &Reference,
    problem: &Problem,
    grid: &Grid,
    eps: f64,
) {
    let rms = reference.rms();
    let nbins = DISTANCE_BINS.len();
    println!();
    println!(
        "Each level's part of the ∇φ error (φ's in brackets) by the target's distance |v| from \
         the centre of its box on that level, over the half-width r: the RMS of the part over \
         the reference's RMS, over ε, and the targets per bin."
    );
    println!();
    let heads: Vec<String> = (0..nbins)
        .map(|b| match DISTANCE_BINS.get(b + 1) {
            Some(hi) => format!("{:.2}–{:.2}", DISTANCE_BINS[b], hi),
            None => format!("≥ {:.2}", DISTANCE_BINS[b]),
        })
        .collect();
    println!("| level | {} |", heads.join(" | "));
    println!("| ---: |{}", " ---: |".repeat(nbins));
    for (l, errors) in raised.iter().enumerate() {
        let Some(errors) = errors else { continue };
        let mut bins = vec![(0usize, [0.0f64; 2]); nbins];
        for (j, &i) in reference.targets.iter().enumerate() {
            let o = offsets(problem.points[i], grid, l);
            let distance = o[l].iter().map(|v| v * v).sum::<f64>().sqrt();
            let b = DISTANCE_BINS
                .iter()
                .rposition(|&e| distance >= e)
                .unwrap_or(0);
            let part = Errors {
                values: vec![std::array::from_fn(|c| {
                    base.errors.values[j][c] - errors.values[j][c]
                })],
            }
            .norms(0);
            bins[b].0 += 1;
            for (sum, e) in bins[b].1.iter_mut().zip(part) {
                *sum += e * e;
            }
        }
        let cells: Vec<String> = bins
            .iter()
            .map(|(n, sums)| {
                if *n == 0 {
                    "–".into()
                } else {
                    let r = |f: usize| (sums[f] / *n as f64).sqrt() / rms[f] / eps;
                    format!("{:.3} ({:.3}), {n}", r(1), r(0))
                }
            })
            .collect();
        println!("| {l} | {} |", cells.join(" | "));
    }
}
