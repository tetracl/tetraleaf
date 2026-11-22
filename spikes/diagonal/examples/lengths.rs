//! The scan of T's length (Helmholtz CH4.3d1, docs/helmholtz/tasks/CH4.3d1.md): every
//! distinct (kind, k a, p) that CH4.3c2's rule (the FMM's) gives at ε = 10⁻³ … 10⁻⁸ below the
//! pruning cut, every class and placement, at the lengths of `lengths::lengths`, as rows of a
//! CSV file.
//!
//! ```text
//! cargo run -p tetraleaf-spike-diagonal --release --example lengths -- \
//!     [--kinds real,complex0.1] [--ka 2.83,4,...] [--threads N] [--shard i/n] \
//!     [--out spikes/diagonal/results/lengths.csv] [--list] [--classes i,j,...] [--dd-max-p 100]
//! ```
//!
//! Units already in the output file are skipped, so a run can be resumed; `--shard i/n` takes
//! every n-th unit of the list ordered by cost (largest p first). `--list` prints the units
//! and their lengths and stops.

use std::collections::{BTreeMap, HashSet};
use std::fs::OpenOptions;
use std::io::{BufRead, BufReader, Write};

use tetraleaf_spike_diagonal::classes::classes;
use tetraleaf_spike_diagonal::lengths::{Config, lengths, run_lengths};
use tetraleaf_spike_diagonal::rules::refit;
use tetraleaf_spike_truncation::chain::Kind;
use tetraleaf_spike_truncation::geometry::Placement;

/// D1's ε, 10⁻³ … 10⁻⁸.
const EPSILONS: [f64; 6] = [1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8];

/// S6's half-octave grid 2^{j/2}, j = 3..15, the levels 2–5 of W1 at kD = 500 and level 2 at
/// kD = 1000 (κa = (√3/2) kD/2^l).
fn ka_grid() -> Vec<f64> {
    let mut v: Vec<f64> = (3..=15).map(|j| 2f64.powf(j as f64 / 2.0)).collect();
    let r = 3f64.sqrt() / 2.0;
    v.extend(
        [
            1000.0 / 4.0,
            500.0 / 4.0,
            500.0 / 8.0,
            500.0 / 16.0,
            500.0 / 32.0,
        ]
        .map(|x| x * r),
    );
    v.sort_by(f64::total_cmp);
    v
}

/// S6's pruning test (`examples/sweep.rs`).
fn below_cut(kind: Kind, ka: f64, eps: f64) -> bool {
    let im_ka = match kind {
        Kind::Real => 0.0,
        Kind::Damped(eta) => eta * ka,
        Kind::Imaginary => ka,
    };
    im_ka * 2.0 / 3f64.sqrt() <= (1.0 / eps).ln()
}

fn kind_of(name: &str) -> Kind {
    match name {
        "real" => Kind::Real,
        "complex0.1" => Kind::Damped(0.1),
        "complex1" => Kind::Damped(1.0),
        _ => panic!("unknown kind {name}"),
    }
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let get = |flag: &str| {
        args.iter()
            .position(|a| a == flag)
            .map(|i| args[i + 1].clone())
    };
    let kinds: Vec<Kind> = get("--kinds")
        .unwrap_or_else(|| "real,complex0.1".into())
        .split(',')
        .map(kind_of)
        .collect();
    let kas: Vec<f64> = get("--ka").map_or_else(ka_grid, |s| {
        s.split(',').map(|x| x.parse().unwrap()).collect()
    });
    let mut cfg = Config::default();
    if let Some(t) = get("--threads") {
        cfg.base.threads = t.parse().unwrap();
    }
    if let Some(t) = get("--dd-max-p") {
        cfg.dd_max_p = t.parse().unwrap();
    }
    let (shard, shards) = get("--shard").map_or((0, 1), |s| {
        let (a, b) = s.split_once('/').unwrap();
        (a.parse().unwrap(), b.parse().unwrap())
    });
    let out_path = get("--out")
        .unwrap_or_else(|| concat!(env!("CARGO_MANIFEST_DIR"), "/results/lengths.csv").to_string());
    // (kind, k a, p) → the ε that give p.
    let mut units: Vec<(Kind, f64, usize, Vec<f64>)> = Vec::new();
    for &kind in &kinds {
        for &ka in &kas {
            let mut by_p: BTreeMap<usize, Vec<f64>> = BTreeMap::new();
            for &eps in &EPSILONS {
                if below_cut(kind, ka, eps) {
                    by_p.entry(refit(ka, eps)).or_default().push(eps);
                }
            }
            units.extend(by_p.into_iter().map(|(p, e)| (kind, ka, p, e)));
        }
    }
    units.sort_by(|a, b| b.2.cmp(&a.2).then(a.1.total_cmp(&b.1)));
    let mine: Vec<_> = units
        .iter()
        .enumerate()
        .filter(|(i, _)| i % shards == shard)
        .map(|(_, u)| u.clone())
        .collect();
    if args.iter().any(|a| a == "--list") {
        let mut cost = 0.0;
        for (kind, ka, p, e) in &mine {
            let ls = lengths(*ka, *p, e);
            cost += (*p as f64).powi(3) * ls.len() as f64;
            println!("{} {ka:.3} {p} {e:?} {ls:?}", kind.name());
        }
        println!("{} units, Σ p³ × lengths = {cost:.3e}", mine.len());
        return;
    }
    let mut done: HashSet<(String, String, usize)> = HashSet::new();
    if let Ok(f) = std::fs::File::open(&out_path) {
        for line in BufReader::new(f).lines().skip(1) {
            let line = line.unwrap();
            let c: Vec<&str> = line.split(',').collect();
            done.insert((c[0].to_string(), c[1].to_string(), c[2].parse().unwrap()));
        }
    }
    let fresh = done.is_empty();
    let mut out = OpenOptions::new()
        .create(true)
        .append(true)
        .open(&out_path)
        .expect("output file");
    if fresh {
        writeln!(
            out,
            "kind,ka,p,eps,class,rep,d2,placement,l,log10_term,log10_tail,log10_max_t,\
             trunc_pot,trunc_grad,dd_pot,dd_grad,dev_pot,dev_grad,round_pot,round_grad,\
             total_pot,total_grad,seconds"
        )
        .unwrap();
    }
    let all = classes();
    let ids: Vec<usize> = get("--classes").map_or_else(
        || (0..all.len()).collect(),
        |s| s.split(',').map(|x| x.parse().unwrap()).collect(),
    );
    for (n, (kind, ka, p, eps)) in mine.iter().enumerate() {
        let key = (kind.name(), format!("{ka:.3}"), *p);
        if done.contains(&key) {
            continue;
        }
        let ls = lengths(*ka, *p, eps);
        let res = run_lengths(*kind, *ka, *p, &ls, &ids, &cfg);
        let eps_text = eps
            .iter()
            .map(|e| format!("{e:e}"))
            .collect::<Vec<_>>()
            .join(";");
        let mut text = String::new();
        for c in &res.classes {
            let rep = format!("{}:{}:{}", c.info.rep[0], c.info.rep[1], c.info.rep[2]);
            for (pl, &placement) in Placement::ALL.iter().enumerate() {
                for (li, &l) in res.lengths.iter().enumerate() {
                    let e = c.errors[pl][li];
                    text += &format!(
                        "{},{:.3},{},{},{},{},{},{},{},{:.4},{:.4},{:.4},{:.3e},{:.3e},{:.3e},{:.3e},{:.3e},{:.3e},{:.3e},{:.3e},{:.3e},{:.3e},{:.2}\n",
                        kind.name(),
                        ka,
                        p,
                        eps_text,
                        c.class,
                        rep,
                        c.info.d2,
                        placement.name(),
                        l,
                        c.log10_term[li],
                        c.log10_tail[li],
                        c.log10_max_t[li],
                        c.trunc[pl].0,
                        c.trunc[pl].1,
                        e.dd.0,
                        e.dd.1,
                        e.dev.0,
                        e.dev.1,
                        e.round.0,
                        e.round.1,
                        e.total.0,
                        e.total.1,
                        c.seconds
                    );
                }
            }
        }
        out.write_all(text.as_bytes()).unwrap();
        out.flush().unwrap();
        eprintln!(
            "[{}/{}] {} k a {ka:.3} p {p} ({} lengths): {:.1} s",
            n + 1,
            mine.len(),
            kind.name(),
            ls.len(),
            res.seconds
        );
    }
}
