//! The Helmholtz operator (Helmholtz CH3.1a): its sizes and errors, and its bodies pair by
//! pair (`PairOperator`) in chains against the direct sum
//! (`tetraleaf_ref::helmholtz::p2p::direct_sum`, 1/(4π) applied, through `Densities`' output
//! pass), relative L2 over the targets of the potential and of the gradient, for real and
//! complex k and real and complex charges; at named low degrees, where truncation
//! dominates, the same chains of `tetraleaf_ref::helmholtz` in absolute units; several
//! vectors bit for bit the one-vector operator; real charges against complex ones with
//! Im q = 0. The level calls through the drivers need a plan (MPI) and are CH3.1b's.
//! Yukawa (Helmholtz CH2.1b): the same chains against the direct sum at k = iκ, its complex
//! charges bit for bit the real path's two real vectors, Dense against Rotation, the M2L it
//! refuses, and the reference chains of `tetraleaf_ref`'s `_yukawa` forms.
//!
//! The geometry is the unit cube (dyadic, so every frame and leaf coordinate is exact),
//! with boxes on levels 2–4 of a five-level tree. The chains:
//!
//! - **far**: P2M on a leaf of level 4, M2M to level 2, M2L there, L2L to level 4, L2P;
//! - **v4**: P2M, M2L on level 4, L2P;
//! - **x**: P2L from a leaf of level 3 into a box of level 4 (its X list), L2P;
//! - **w**: P2M on a box of level 4, M2P at a leaf of level 3 (its W list);
//! - **near**: P2P from a neighbour and from the leaf itself.

use tetraleaf_exec::geometry::{self, Domain, radius};
use tetraleaf_exec::helmholtz::{
    ChargeKind, Charges, Densities, HelmholtzOperator, OperatorError, Output, Wavenumber,
    WavenumberKind, s7_degree,
};
use tetraleaf_octree::{MortonKey, PhysicalBox, morton};
use tetraleaf_plan::operator::{FmmSizes, PairOperator};
use tetraleaf_tables::geometry::m2l_offset_index;

use super::common::{C, Rng, exact, output_of, relative_l2};

const LEVELS: usize = 5;

fn domain() -> Domain {
    Domain::new(&PhysicalBox::new([0.0, 0.0, 0.0, 1.0, 1.0, 1.0])).unwrap()
}

fn key(index: [usize; 3], level: usize) -> MortonKey {
    morton::from_index_and_level(index, level)
}

/// n points inside `key`, away from its faces by a tenth of its half-width.
fn points_in(key: MortonKey, n: usize, rng: &mut Rng) -> Vec<[f64; 3]> {
    let domain = domain();
    let (c, r) = (
        geometry::centre(key, &domain),
        radius(morton::level(key) as usize, &domain),
    );
    (0..n)
        .map(|_| core::array::from_fn(|i| c[i] + 0.9 * r * rng.next()))
        .collect()
}

/// m charge vectors of `kind` at n points: (Re q, Im q), Im q zero for real charges.
fn charges(kind: ChargeKind, m: usize, n: usize, rng: &mut Rng) -> Vec<[Vec<f64>; 2]> {
    (0..m)
        .map(|_| {
            let re = (0..n).map(|_| rng.next()).collect();
            let im = match kind {
                ChargeKind::Real => vec![0.0; n],
                ChargeKind::Complex => (0..n).map(|_| rng.next()).collect(),
            };
            [re, im]
        })
        .collect()
}

fn as_charges(kind: ChargeKind, planes: &[[Vec<f64>; 2]]) -> Vec<Charges<'_, f64>> {
    planes
        .iter()
        .map(|[re, im]| match kind {
            ChargeKind::Real => Charges::Real(re),
            ChargeKind::Complex => Charges::Complex { re, im },
        })
        .collect()
}

/// A leaf's source chunk: leaf-scaled coordinates, then the strengths.
fn source_chunk(
    d: &Densities,
    leaf: MortonKey,
    points: &[[f64; 3]],
    q: &[Charges<'_, f64>],
) -> Vec<f64> {
    let n = points.len();
    let mut chunk = vec![0.0; n * d.source_point_size()];
    for (u, &x) in chunk.as_chunks_mut::<3>().0.iter_mut().zip(points) {
        *u = geometry::leaf_coordinates(x, leaf, &domain());
    }
    let indices: Vec<usize> = (0..n).collect();
    d.write_leaf_charges(q, &indices, &mut chunk[3 * n..]);
    chunk
}

fn target_input(leaf: MortonKey, points: &[[f64; 3]]) -> Vec<f64> {
    points
        .iter()
        .flat_map(|&x| geometry::leaf_coordinates::<f64>(x, leaf, &domain()))
        .collect()
}

/// A problem: k, the charge kind, m vectors, and the operator at S7's degrees at ε.
struct Problem {
    k: Wavenumber,
    kind: ChargeKind,
    densities: Densities,
    op: HelmholtzOperator<f64>,
}

impl Problem {
    fn new(k: Wavenumber, kind: ChargeKind, m: usize, epsilon: f64) -> Self {
        let densities = Densities::new(k.kind(), kind, m);
        // S7's per-pair degrees: the chains are single pairs, which S7's closed form bounds
        // (the default rule is calibrated on the whole FMM, Helmholtz CH3.3).
        let degrees: Vec<usize> = (0..LEVELS)
            .map(|l| match l {
                0 | 1 => 0,
                _ => s7_degree(&k, 3f64.sqrt() * radius(l, &domain()), epsilon),
            })
            .collect();
        // No stored M2L blocks: the few M2L calls form theirs (the same bits).
        let op =
            HelmholtzOperator::with_block_budget(k, densities, &domain(), &degrees, true, 64, 0)
                .unwrap();
        Self {
            k,
            kind,
            densities,
            op,
        }
    }

    fn multipole(&self, level: usize) -> Vec<f64> {
        vec![0.0; self.op.multipole_size(level)]
    }

    fn local(&self, level: usize) -> Vec<f64> {
        vec![0.0; self.op.local_size(level)]
    }

    fn output_chunk(&self, n: usize) -> Vec<f64> {
        vec![0.0; n * self.op.target_output_point_size()]
    }

    /// The physical φ and ∇φ of vector v of a leaf's target output chunk.
    fn physical(
        &self,
        leaf: MortonKey,
        chunk: &[f64],
        n: usize,
        v: usize,
    ) -> (Vec<C>, Vec<[C; 3]>) {
        let mut output = Output::new(&self.densities, n, true);
        let r = radius(morton::level(leaf) as usize, &domain());
        let indices: Vec<usize> = (0..n).collect();
        self.densities
            .write_leaf_output(chunk, r, true, &indices, &mut output);
        output_of(&output, v)
    }
}

/// The errors of a chain's output (vector v) against the direct sum of its sources.
fn errors(
    pb: &Problem,
    target: MortonKey,
    chunk: &[f64],
    targets: &[[f64; 3]],
    sources: &[[f64; 3]],
    q: &[[Vec<f64>; 2]],
    v: usize,
) -> (f64, f64) {
    let (phi, grad) = pb.physical(target, chunk, targets.len(), v);
    let qv: Vec<C> = (0..sources.len())
        .map(|j| C::new(q[v][0][j], q[v][1][j]))
        .collect();
    let (phi_e, grad_e) = exact(&pb.k, sources, &qv, targets);
    let flat = |g: &[[C; 3]]| g.iter().flatten().copied().collect::<Vec<C>>();
    (
        relative_l2(&phi, &phi_e),
        relative_l2(&flat(&grad), &flat(&grad_e)),
    )
}

/// How a chain draws its charges: `m` vectors of `kind`, of which it keeps all, or only
/// vector `pick`, and writes them as the problem's charge kind (real draws as complex
/// charges have Im q = 0).
///
/// With `split`, each drawn vector becomes two real ones, Re q then Im q (Yukawa's complex
/// charges as the real path's 2m vectors, CH2.1b).
#[derive(Clone, Copy)]
struct Draw {
    kind: ChargeKind,
    m: usize,
    pick: Option<usize>,
    split: bool,
}

impl Draw {
    fn all(kind: ChargeKind, m: usize) -> Self {
        Self {
            kind,
            m,
            pick: None,
            split: false,
        }
    }

    fn charges(self, n: usize, rng: &mut Rng) -> Vec<[Vec<f64>; 2]> {
        let q = charges(self.kind, self.m, n, rng);
        let q = match self.pick {
            None => q,
            Some(v) => vec![q[v].clone()],
        };
        if !self.split {
            return q;
        }
        q.into_iter()
            .flat_map(|[re, im]| [[re, vec![0.0; n]], [im, vec![0.0; n]]])
            .collect()
    }
}

/// A chain's target leaf, its target output chunk and its worst errors over the vectors.
struct Chain {
    name: &'static str,
    target: MortonKey,
    n: usize,
    output: Vec<f64>,
    errors: (f64, f64),
}

/// Runs every chain of the module documentation for `pb`.
fn chains(pb: &mut Problem, seed: u64, draw: Draw) -> Vec<Chain> {
    let m = pb.densities.vectors();
    let mut rng = Rng(seed);
    let mut out = Vec::new();
    let mut finish = |pb: &Problem,
                      name,
                      target,
                      output: Vec<f64>,
                      targets: &[[f64; 3]],
                      sources: &[[f64; 3]],
                      q: &[[Vec<f64>; 2]]| {
        let errors = (0..m).fold((0.0f64, 0.0f64), |w, v| {
            let e = errors(pb, target, &output, targets, sources, q, v);
            (w.0.max(e.0), w.1.max(e.1))
        });
        out.push(Chain {
            name,
            target,
            n: targets.len(),
            output,
            errors,
        });
    };
    let kind = pb.kind;

    // far: s4 = (1, 1, 1) and t4 = (12, 13, 11); their level-2 ancestors (0, 0, 0) and
    // (3, 3, 2) are a V-list pair, offset (3, 3, 2).
    {
        let (s4, t4) = (key([1, 1, 1], 4), key([12, 13, 11], 4));
        let (sources, targets) = (points_in(s4, 20, &mut rng), points_in(t4, 15, &mut rng));
        let q = draw.charges(sources.len(), &mut rng);
        let chunk = source_chunk(&pb.densities, s4, &sources, &as_charges(kind, &q));
        let up = [s4, morton::parent(s4).unwrap(), key([0, 0, 0], 2)];
        let down = [key([3, 3, 2], 2), morton::parent(t4).unwrap(), t4];
        let mut multipoles = [pb.multipole(4), pb.multipole(3), pb.multipole(2)];
        pb.op.p2m(up[0], &chunk, &mut multipoles[0]);
        for i in 0..2 {
            let (lo, hi) = multipoles.split_at_mut(i + 1);
            let o = morton::child_index(up[i]);
            pb.op.m2m(up[i], up[i + 1], o, &lo[i], &mut hi[0]);
        }
        let mut locals = [pb.local(2), pb.local(3), pb.local(4)];
        let d = m2l_offset_index([3, 3, 2]).unwrap();
        pb.op.m2l(up[2], down[0], d, &multipoles[2], &mut locals[0]);
        for i in 0..2 {
            let (lo, hi) = locals.split_at_mut(i + 1);
            let o = morton::child_index(down[i + 1]);
            pb.op.l2l(down[i], down[i + 1], o, &lo[i], &mut hi[0]);
        }
        let mut output = pb.output_chunk(targets.len());
        pb.op
            .l2p(t4, &locals[2], &target_input(t4, &targets), &mut output);
        finish(pb, "far", t4, output, &targets, &sources, &q);
    }

    // v4: s4 = (1, 1, 1), t4 = (3, 1, 1), offset (2, 0, 0) on level 4.
    {
        let (s4, t4) = (key([1, 1, 1], 4), key([3, 1, 1], 4));
        let (sources, targets) = (points_in(s4, 20, &mut rng), points_in(t4, 15, &mut rng));
        let q = draw.charges(sources.len(), &mut rng);
        let chunk = source_chunk(&pb.densities, s4, &sources, &as_charges(kind, &q));
        let mut multipole = pb.multipole(4);
        pb.op.p2m(s4, &chunk, &mut multipole);
        let mut local = pb.local(4);
        let d = m2l_offset_index([2, 0, 0]).unwrap();
        pb.op.m2l(s4, t4, d, &multipole, &mut local);
        let mut output = pb.output_chunk(targets.len());
        pb.op
            .l2p(t4, &local, &target_input(t4, &targets), &mut output);
        finish(pb, "v4", t4, output, &targets, &sources, &q);
    }

    // x: the leaf s3 = (0, 0, 0) on level 3 into t4 = (3, 0, 0) on level 4, whose parent
    // (1, 0, 0) is adjacent to s3 while t4 is not.
    {
        let (s3, t4) = (key([0, 0, 0], 3), key([3, 0, 0], 4));
        let (sources, targets) = (points_in(s3, 20, &mut rng), points_in(t4, 15, &mut rng));
        let q = draw.charges(sources.len(), &mut rng);
        let chunk = source_chunk(&pb.densities, s3, &sources, &as_charges(kind, &q));
        let mut local = pb.local(4);
        pb.op.p2l(s3, t4, &chunk, &mut local);
        let mut output = pb.output_chunk(targets.len());
        pb.op
            .l2p(t4, &local, &target_input(t4, &targets), &mut output);
        finish(pb, "x", t4, output, &targets, &sources, &q);
    }

    // w: the box s4 = (3, 0, 0) on level 4 at the leaf t3 = (0, 0, 0) on level 3.
    {
        let (s4, t3) = (key([3, 0, 0], 4), key([0, 0, 0], 3));
        let (sources, targets) = (points_in(s4, 20, &mut rng), points_in(t3, 15, &mut rng));
        let q = draw.charges(sources.len(), &mut rng);
        let chunk = source_chunk(&pb.densities, s4, &sources, &as_charges(kind, &q));
        let mut multipole = pb.multipole(4);
        pb.op.p2m(s4, &chunk, &mut multipole);
        let mut output = pb.output_chunk(targets.len());
        pb.op
            .m2p(s4, t3, &multipole, &target_input(t3, &targets), &mut output);
        finish(pb, "w", t3, output, &targets, &sources, &q);
    }

    // near: P2P from the neighbour (2, 1, 1) and from the leaf (1, 1, 1) itself.
    {
        let (s4, t4) = (key([2, 1, 1], 4), key([1, 1, 1], 4));
        let own = points_in(t4, 12, &mut rng);
        let other = points_in(s4, 12, &mut rng);
        let q_own = draw.charges(own.len(), &mut rng);
        let q_other = draw.charges(other.len(), &mut rng);
        let own_chunk = source_chunk(&pb.densities, t4, &own, &as_charges(kind, &q_own));
        let other_chunk = source_chunk(&pb.densities, s4, &other, &as_charges(kind, &q_other));
        let input = target_input(t4, &own);
        let mut output = pb.output_chunk(own.len());
        pb.op.p2p(s4, t4, &other_chunk, &input, &mut output);
        pb.op.p2p(t4, t4, &own_chunk, &input, &mut output);
        let sources: Vec<[f64; 3]> = other.iter().chain(&own).copied().collect();
        let q: Vec<[Vec<f64>; 2]> = q_other
            .iter()
            .zip(&q_own)
            .map(|(a, b)| core::array::from_fn(|p| [a[p].clone(), b[p].clone()].concat()))
            .collect();
        finish(pb, "near", t4, output, &own, &sources, &q);
    }
    out
}

/// The chains' errors within the tolerance ε at S7's degrees (the P2P to rounding, 1e-13),
/// for real k at
/// kD = 1, 20 and 100 and complex k at kD = 20 + 2i, real and complex charges, ε = 10⁻⁶
/// and 10⁻¹⁰.
#[test]
fn chains_against_the_direct_sum() {
    let ks = [(1.0, 0.0), (20.0, 0.0), (100.0, 0.0), (20.0, 2.0)];
    for (re, im) in ks {
        for kind in [ChargeKind::Real, ChargeKind::Complex] {
            for epsilon in [1e-6, 1e-10] {
                let k = Wavenumber::new(re, im).unwrap();
                let mut pb = Problem::new(k, kind, 1, epsilon);
                for chain in chains(&mut pb, 7, Draw::all(kind, 1)) {
                    let (name, (phi, grad)) = (chain.name, chain.errors);
                    println!(
                        "k = {re} + {im}i, {kind:?} charges, ε = {epsilon:e}, {name}: φ {phi:.1e}, ∇φ {grad:.1e}"
                    );
                    let bound = if name == "near" { 1e-13 } else { epsilon };
                    assert!(
                        phi <= bound && grad <= bound,
                        "k = {re} + {im}i, {kind:?}, ε = {epsilon:e}, {name}: φ {phi:e}, ∇φ {grad:e}"
                    );
                }
            }
        }
    }
}

/// Vector v of an operator of three vectors is the one-vector operator on vector v, bit
/// for bit, in every chain (the physical output, whose scaling is exact for equal input).
#[test]
fn each_vector_is_the_one_vector_operator() {
    for (re, im) in [(20.0, 0.0), (20.0, 2.0), (0.0, 20.0)] {
        for kind in [ChargeKind::Real, ChargeKind::Complex] {
            let k = Wavenumber::new(re, im).unwrap();
            let mut three = Problem::new(k, kind, 3, 1e-6);
            let many = chains(&mut three, 11, Draw::all(kind, 3));
            for v in 0..3 {
                let mut one = Problem::new(k, kind, 1, 1e-6);
                let pick = Draw {
                    kind,
                    m: 3,
                    pick: Some(v),
                    split: false,
                };
                for (a, b) in many.iter().zip(chains(&mut one, 11, pick)) {
                    let x = three.physical(a.target, &a.output, a.n, v);
                    let y = one.physical(b.target, &b.output, b.n, 0);
                    assert!(x == y, "k = {re} + {im}i, {kind:?}, {}, vector {v}", a.name);
                }
            }
        }
    }
}

/// Real charges against complex ones with Im q = 0: the same target output within 4 ε of
/// its largest value (printed: whether bit for bit).
#[test]
fn real_charges_are_complex_ones_with_zero_imaginary_parts() {
    for (re, im) in [(1.0, 0.0), (20.0, 0.0), (20.0, 2.0)] {
        let k = Wavenumber::new(re, im).unwrap();
        let mut real = Problem::new(k, ChargeKind::Real, 1, 1e-8);
        let mut complex = Problem::new(k, ChargeKind::Complex, 1, 1e-8);
        let draw = Draw::all(ChargeKind::Real, 1);
        let a = chains(&mut real, 3, draw);
        let b = chains(&mut complex, 3, draw);
        let mut bitwise = true;
        for (x, y) in a.iter().zip(&b) {
            let scale = y.output.iter().fold(0.0f64, |m, v| m.max(v.abs()));
            for (u, w) in x.output.iter().zip(&y.output) {
                bitwise &= u == w;
                assert!(
                    (u - w).abs() <= 4.0 * f64::EPSILON * scale,
                    "k = {re} + {im}i, {}: {u:e} against {w:e}",
                    x.name
                );
            }
        }
        println!("k = {re} + {im}i: real charges bit for bit complex ones: {bitwise}");
    }
}

#[test]
fn sizes_and_errors() {
    let domain = domain();
    let k = Wavenumber::new(10.0, 0.0).unwrap();
    for (charges, m) in [(ChargeKind::Real, 1), (ChargeKind::Complex, 3)] {
        let d = Densities::new(WavenumberKind::Real, charges, m);
        let op = HelmholtzOperator::<f64>::new(k, d, &domain, &[0, 0, 9, 7], false, 8).unwrap();
        for (level, p) in [(0, 0), (2, 9), (3, 7)] {
            assert_eq!(op.degree(level), p);
            assert_eq!(op.multipole_size(level), d.multipole_size(p));
            assert_eq!(op.local_size(level), d.local_size(p));
        }
        assert_eq!(op.source_point_size(), d.source_point_size());
        assert_eq!(op.target_input_point_size(), 3);
        assert_eq!(
            op.target_output_point_size(),
            d.target_output_point_size(false)
        );
        assert_eq!(op.nlevels(), 4);
    }
    // Yukawa (CH2.1b): (p + 1)² reals per expansion with real charges, two real ones per
    // complex vector.
    let yukawa = Wavenumber::new(0.0, 10.0).unwrap();
    for (charges, m) in [(ChargeKind::Real, 1), (ChargeKind::Complex, 3)] {
        let d = Densities::new(WavenumberKind::Yukawa, charges, m);
        let op = HelmholtzOperator::<f64>::new(yukawa, d, &domain, &[0, 0, 4], true, 8).unwrap();
        assert_eq!(op.multipole_size(2), d.multipole_size(4));
        assert_eq!(op.local_size(2), d.local_size(4));
        assert_eq!(
            op.target_output_point_size(),
            d.target_output_point_size(true)
        );
    }
    let d = Densities::new(WavenumberKind::Complex, ChargeKind::Real, 1);
    assert_eq!(
        HelmholtzOperator::<f64>::new(k, d, &domain, &[0, 0, 4], false, 8).err(),
        Some(OperatorError::KindMismatch {
            densities: WavenumberKind::Complex,
            k: WavenumberKind::Real
        })
    );
    let d = Densities::new(WavenumberKind::Real, ChargeKind::Real, 1);
    assert_eq!(
        HelmholtzOperator::<f64>::new(k, d, &domain, &[], false, 8).err(),
        Some(OperatorError::NoLevels)
    );
}

/// The operator at named low degrees, where truncation dominates the error, against the
/// same chains of `tetraleaf_ref::helmholtz` in absolute units (the trusted path: leaf
/// operators and `rcr` with the levels of half-widths r_l, the frames at the boxes'
/// centres): far, v4, x and w, potential and gradient, relative L2 within 2e-13 of each
/// other (1.2e-14 at most on the M6: the tables' exact-class angles and distances round
/// apart from `rcr`'s) while both are 10⁻⁵ to 10⁻¹ from the direct sum.
#[test]
fn low_degree_chains_are_the_reference_chains() {
    use tetraleaf_ref::Frame;
    use tetraleaf_ref::helmholtz::{Level, Workspace, leaf, rcr};

    let domain = domain();
    let frame = |key: MortonKey| {
        let level = morton::level(key) as usize;
        Frame::new(geometry::centre(key, &domain), radius(level, &domain))
    };
    for (re, im) in [(20.0, 0.0), (20.0, 2.0), (60.0, 0.0)] {
        for kind in [ChargeKind::Real, ChargeKind::Complex] {
            let k = Wavenumber::new(re, im).unwrap();
            let kc = C::new(re, im);
            let level = |l: usize| Level::new(kc, radius(l, &domain));
            let degrees = [0, 0, 8, 6, 5];
            let densities = Densities::new(k.kind(), kind, 1);
            let op = HelmholtzOperator::new(k, densities, &domain, &degrees, true, 64).unwrap();
            let mut pb = Problem {
                k,
                kind,
                densities,
                op,
            };
            let ours = chains(&mut pb, 5, Draw::all(kind, 1));
            let mut rng = Rng(5);
            let mut ws = Workspace::new(8);
            let zero = C::new(0.0, 0.0);
            // The reference chains, drawing the same points and charges as `chains`.
            for chain in ours.iter().filter(|c| c.name != "near") {
                let (s, t) = match chain.name {
                    "far" => (key([1, 1, 1], 4), key([12, 13, 11], 4)),
                    "v4" => (key([1, 1, 1], 4), key([3, 1, 1], 4)),
                    "x" => (key([0, 0, 0], 3), key([3, 0, 0], 4)),
                    _ => (key([3, 0, 0], 4), key([0, 0, 0], 3)),
                };
                let sources = points_in(s, 20, &mut rng);
                let targets = points_in(t, 15, &mut rng);
                let q = Draw::all(kind, 1).charges(sources.len(), &mut rng);
                let q: Vec<C> = (0..sources.len())
                    .map(|j| C::new(q[0][0][j], q[0][1][j]))
                    .collect();
                let (lt, ls) = (morton::level(t) as usize, morton::level(s) as usize);
                let (mut phi, mut grad) =
                    (vec![zero; targets.len()], vec![[zero; 3]; targets.len()]);
                let local_of = |ws: &mut Workspace<f64>| -> Vec<f64> {
                    match chain.name {
                        "far" => {
                            let up = [s, morton::parent(s).unwrap(), key([0, 0, 0], 2)];
                            let down = [key([3, 3, 2], 2), morton::parent(t).unwrap(), t];
                            let len = |l: usize| 2 * (degrees[l] + 1) * (degrees[l] + 1);
                            let mut m = vec![0.0; len(4)];
                            leaf::p2m(5, &level(4), &frame(s), &sources, &q, ws, &mut m);
                            for i in 0..2 {
                                let (lc, lp) = (4 - i, 3 - i);
                                let table = rcr::m2m_table(degrees[lc], degrees[lp], &level(lc));
                                let mut parent = vec![0.0; len(lp)];
                                rcr::m2m(
                                    degrees[lc],
                                    degrees[lp],
                                    &level(lc),
                                    &level(lp),
                                    &frame(up[i]),
                                    &frame(up[i + 1]),
                                    &table,
                                    ws,
                                    &m,
                                    &mut parent,
                                );
                                m = parent;
                            }
                            let table =
                                rcr::m2l_table(8, 8, &level(2), &frame(up[2]), &frame(down[0]));
                            let mut l = vec![0.0; len(2)];
                            rcr::m2l(
                                8,
                                8,
                                &level(2),
                                &frame(up[2]),
                                &frame(down[0]),
                                &table,
                                ws,
                                &m,
                                &mut l,
                            );
                            for i in 0..2 {
                                let (lp, lc) = (2 + i, 3 + i);
                                let table = rcr::l2l_table(degrees[lp], degrees[lc], &level(lc));
                                let mut child = vec![0.0; len(lc)];
                                rcr::l2l(
                                    degrees[lp],
                                    degrees[lc],
                                    &level(lp),
                                    &level(lc),
                                    &frame(down[i]),
                                    &frame(down[i + 1]),
                                    &table,
                                    ws,
                                    &l,
                                    &mut child,
                                );
                                l = child;
                            }
                            l
                        }
                        "v4" => {
                            let mut m = vec![0.0; 72];
                            leaf::p2m(5, &level(4), &frame(s), &sources, &q, ws, &mut m);
                            let table = rcr::m2l_table(5, 5, &level(4), &frame(s), &frame(t));
                            let mut l = vec![0.0; 72];
                            rcr::m2l(
                                5,
                                5,
                                &level(4),
                                &frame(s),
                                &frame(t),
                                &table,
                                ws,
                                &m,
                                &mut l,
                            );
                            l
                        }
                        _ => {
                            let mut l = vec![0.0; 72];
                            leaf::p2l(5, &level(4), &frame(t), &sources, &q, ws, &mut l);
                            l
                        }
                    }
                };
                if chain.name == "w" {
                    let mut m = vec![0.0; 72];
                    leaf::p2m(5, &level(ls), &frame(s), &sources, &q, &mut ws, &mut m);
                    leaf::m2p(
                        5,
                        &level(ls),
                        &frame(s),
                        &m,
                        &targets,
                        &mut ws,
                        &mut phi,
                        Some(&mut grad),
                    );
                } else {
                    let l = local_of(&mut ws);
                    leaf::l2p(
                        5,
                        &level(lt),
                        &frame(t),
                        &l,
                        &targets,
                        &mut ws,
                        &mut phi,
                        Some(&mut grad),
                    );
                }
                // φ̃ = 4πφ: the reference's fields carry the factor 4π.
                let scale = 1.0 / (4.0 * std::f64::consts::PI);
                let phi: Vec<C> = phi.iter().map(|&x| x * scale).collect();
                let grad: Vec<C> = grad.iter().flatten().map(|&x| x * scale).collect();
                let (our_phi, our_grad) = pb.physical(chain.target, &chain.output, chain.n, 0);
                let our_grad: Vec<C> = our_grad.iter().flatten().copied().collect();
                let (dp, dg) = (relative_l2(&our_phi, &phi), relative_l2(&our_grad, &grad));
                let (tp, tg) = chain.errors;
                println!(
                    "k = {re} + {im}i, {kind:?}, {}: against the reference φ {dp:.1e}, ∇φ {dg:.1e}; against the direct sum φ {tp:.1e}, ∇φ {tg:.1e}",
                    chain.name
                );
                assert!(dp <= 2e-13 && dg <= 2e-13, "{}: {dp:e}, {dg:e}", chain.name);
                assert!(
                    tp > 1e-9,
                    "{}: truncation not exercised ({tp:e})",
                    chain.name
                );
            }
        }
    }
}

/// P2P with `P2pChoice::Reference` and with `Isa::Scalar` (CH3.4c) against
/// `tetraleaf_ref::helmholtz::p2p::p2p` bit for bit: the neighbour's sources mapped into the
/// target leaf's units (ĉ(s|t) + r̂(s|t) u, exact on the unit cube), k̂ = k r_t, each
/// target's two output columns accumulated by the reference as complex values; then the
/// self pair.
#[test]
fn p2p_is_the_reference_bit_for_bit() {
    use tetraleaf_exec::operator::{Isa, P2pChoice};
    use tetraleaf_ref::helmholtz::p2p;

    let choices = [P2pChoice::Reference, P2pChoice::Isa(Isa::Scalar)];
    for ((re, im), choice) in [(20.0, 0.0), (20.0, 2.0)]
        .into_iter()
        .flat_map(|k| choices.map(|c| (k, c)))
    {
        let k = Wavenumber::new(re, im).unwrap();
        let mut pb = Problem::new(k, ChargeKind::Complex, 1, 1e-6);
        pb.op = pb.op.with_p2p(choice).unwrap();
        let (s, t) = (key([2, 1, 1], 4), key([1, 1, 1], 4));
        let mut rng = Rng(9);
        let (other, own) = (points_in(s, 10, &mut rng), points_in(t, 8, &mut rng));
        let q = Draw::all(ChargeKind::Complex, 1).charges(other.len() + own.len(), &mut rng);
        let q_other = vec![[q[0][0][..10].to_vec(), q[0][1][..10].to_vec()]];
        let q_own = vec![[q[0][0][10..].to_vec(), q[0][1][10..].to_vec()]];
        let other_chunk = source_chunk(&pb.densities, s, &other, &as_charges(pb.kind, &q_other));
        let own_chunk = source_chunk(&pb.densities, t, &own, &as_charges(pb.kind, &q_own));
        let input = target_input(t, &own);
        let mut output = pb.output_chunk(own.len());
        pb.op.p2p(s, t, &other_chunk, &input, &mut output);
        pb.op.p2p(t, t, &own_chunk, &input, &mut output);

        let r = radius(4, &domain());
        let k_hat = C::new(re * r, im * r);
        let frame = geometry::relative_frame::<f64>(s, t);
        let mapped: Vec<[f64; 3]> = other_chunk[..30]
            .as_chunks::<3>()
            .0
            .iter()
            .map(|u| core::array::from_fn(|i| frame.centre[i] + frame.radius * u[i]))
            .collect();
        let own_units: Vec<[f64; 3]> = own_chunk[..24].as_chunks::<3>().0.to_vec();
        let charges = |q: &[[Vec<f64>; 2]]| -> Vec<C> {
            (0..q[0][0].len())
                .map(|j| C::new(q[0][0][j], q[0][1][j]))
                .collect()
        };
        let zero = C::new(0.0, 0.0);
        let n = own.len();
        for i in 0..n {
            let (mut phi, mut grad) = ([zero], [[zero; 3]]);
            let x = &own_units[i..=i];
            p2p::p2p(
                k_hat,
                &mapped,
                &charges(&q_other),
                x,
                &mut phi,
                Some(&mut grad),
            );
            p2p::p2p(
                k_hat,
                &own_units,
                &charges(&q_own),
                x,
                &mut phi,
                Some(&mut grad),
            );
            let d = &pb.densities;
            assert_eq!(
                output[d.potential_index(n, 0, i)],
                phi[0].re,
                "{choice}, target {i}"
            );
            assert_eq!(
                output[d.potential_index(n, 1, i)],
                phi[0].im,
                "{choice}, target {i}"
            );
            for axis in 0..3 {
                assert_eq!(output[d.gradient_index(n, 0, i, axis)], grad[0][axis].re);
                assert_eq!(output[d.gradient_index(n, 1, i, axis)], grad[0][axis].im);
            }
        }
    }
}

/// Yukawa (CH2.1b): the chains' errors within ε at S7's degrees (the P2P to rounding, 1e-13),
/// at κD = 0.01, 1, 20 and 100 (κ a from about 2·10⁻⁴ to 22 on the chains' levels), real and
/// complex charges, ε = 10⁻⁶ and 10⁻¹⁰, against the direct sum at k = iκ. One exception: the
/// far chain at κD = 100 errs about 3·10⁻⁷ at every ε, round-off and not truncation. Its pair on
/// level 2 lies |t| ≈ 1.17 apart and its points ρ_s ≈ 0.05 and ρ_t ≈ 0.13 from their boxes'
/// centres, so the expansion's terms are about e^{−κ(|t| − ρ_s − ρ_t)} ≈ e^{−99} against a sum of
/// about e^{−119}: a cancellation of e^{20} ≈ 5·10⁸ times u, 10⁻⁸–10⁻⁷, the expansion's own
/// (any M2L of the pair shares it; `yukawa_low_degree_chains_are_the_reference_chains` shows the
/// chain is `rcr`'s to 3e-14 at κ = 100), the loss under strong damping of the roadmap's CH2.3
/// (the H1 gate's G3), which its pruning removes. It is held to `STRONG_DAMPING` (printed).
#[test]
fn yukawa_chains_against_the_direct_sum() {
    // The far chain at κD = 100 (above): 3.0e-7 measured on the M6; ten times that.
    const STRONG_DAMPING: f64 = 3e-6;
    for kappa in [0.01, 1.0, 20.0, 100.0] {
        for kind in [ChargeKind::Real, ChargeKind::Complex] {
            for epsilon in [1e-6, 1e-10] {
                let k = Wavenumber::new(0.0, kappa).unwrap();
                let mut pb = Problem::new(k, kind, 1, epsilon);
                for chain in chains(&mut pb, 7, Draw::all(kind, 1)) {
                    let (name, (phi, grad)) = (chain.name, chain.errors);
                    println!(
                        "Yukawa κ = {kappa}, {kind:?} charges, ε = {epsilon:e}, {name}: φ {phi:.1e}, ∇φ {grad:.1e}"
                    );
                    let bound = match name {
                        "near" => 1e-13,
                        "far" if kappa >= 100.0 => STRONG_DAMPING.max(epsilon),
                        _ => epsilon,
                    };
                    assert!(
                        phi <= bound && grad <= bound,
                        "Yukawa κ = {kappa}, {kind:?}, ε = {epsilon:e}, {name}: φ {phi:e}, ∇φ {grad:e}"
                    );
                }
            }
        }
    }
}

/// Yukawa's complex charges are its real path's two real vectors value for value
/// (`Densities`' chunks, CH2.1b): every chain's target output chunk with one complex vector
/// q is bit for bit that of the real path with the two vectors Re q and Im q, under both P2P
/// kernels.
#[test]
fn yukawa_complex_charges_are_two_real_vectors() {
    use tetraleaf_exec::operator::P2pChoice;

    for kappa in [1.0, 20.0] {
        for choice in [P2pChoice::Auto, P2pChoice::Reference] {
            let k = Wavenumber::new(0.0, kappa).unwrap();
            let mut complex = Problem::new(k, ChargeKind::Complex, 1, 1e-8);
            complex.op = complex.op.with_p2p(choice).unwrap();
            let mut real = Problem::new(k, ChargeKind::Real, 2, 1e-8);
            real.op = real.op.with_p2p(choice).unwrap();
            let a = chains(&mut complex, 13, Draw::all(ChargeKind::Complex, 1));
            let split = Draw {
                split: true,
                ..Draw::all(ChargeKind::Complex, 1)
            };
            let b = chains(&mut real, 13, split);
            for (x, y) in a.iter().zip(&b) {
                assert!(
                    x.output == y.output,
                    "κ = {kappa}, {choice}, {}: the chunks differ",
                    x.name
                );
            }
        }
    }
}

/// Yukawa's Dense M2L (CH2.1a's one-plane matrices) against Rotation in the far and v4
/// chains, at the named degrees [0, 0, 8, 6, 5] (so every level fits the dense budget): the
/// target outputs within 1e-14 relative L2 of each other; and the M2L the operator refuses on
/// Yukawa.
#[test]
fn yukawa_dense_against_rotation() {
    use tetraleaf_exec::helmholtz::{DEFAULT_DENSE_BUDGET, HelmholtzM2l};

    for kappa in [0.01, 1.0, 20.0] {
        for kind in [ChargeKind::Real, ChargeKind::Complex] {
            let k = Wavenumber::new(0.0, kappa).unwrap();
            let problem = || {
                let densities = Densities::new(k.kind(), kind, 2);
                let degrees = [0, 0, 8, 6, 5];
                let op =
                    HelmholtzOperator::new(k, densities, &domain(), &degrees, true, 64).unwrap();
                Problem {
                    k,
                    kind,
                    densities,
                    op,
                }
            };
            let (mut rotation, mut dense) = (problem(), problem());
            dense.op = dense.op.with_m2l(HelmholtzM2l::Dense, DEFAULT_DENSE_BUDGET);
            assert!(dense.op.dense_bytes() > 0, "κ = {kappa}: no dense level");
            let a = chains(&mut rotation, 17, Draw::all(kind, 2));
            let b = chains(&mut dense, 17, Draw::all(kind, 2));
            for (x, y) in a
                .iter()
                .zip(&b)
                .filter(|(x, _)| matches!(x.name, "far" | "v4"))
            {
                let diff: f64 = x
                    .output
                    .iter()
                    .zip(&y.output)
                    .map(|(u, w)| (u - w) * (u - w))
                    .sum();
                let norm: f64 = x.output.iter().map(|u| u * u).sum();
                let e = (diff / norm).sqrt();
                println!(
                    "Yukawa κ = {kappa}, {kind:?}, {}: Dense against Rotation {e:.1e}",
                    x.name
                );
                assert!(e <= 1e-14, "κ = {kappa}, {kind:?}, {}: {e:e}", x.name);
            }
        }
    }
    for m2l in [
        HelmholtzM2l::Batched,
        HelmholtzM2l::Classes,
        HelmholtzM2l::PlaneWave,
    ] {
        let refused = std::panic::catch_unwind(|| {
            let k = Wavenumber::new(0.0, 5.0).unwrap();
            Problem::new(k, ChargeKind::Real, 1, 1e-6)
                .op
                .with_m2l(m2l, DEFAULT_DENSE_BUDGET)
        });
        let message = refused.err().map(|e| {
            e.downcast_ref::<String>()
                .cloned()
                .or_else(|| e.downcast_ref::<&str>().map(|s| s.to_string()))
                .unwrap_or_default()
        });
        assert!(
            message
                .as_deref()
                .is_some_and(|m| m.contains("Yukawa") && m.contains("Rotation and Dense")),
            "{m2l:?}: {message:?}"
        );
    }
}

/// Yukawa's operator at named low degrees against the same chains of `tetraleaf_ref`'s
/// `_yukawa` forms in absolute units (as `low_degree_chains_are_the_reference_chains`; real
/// charges, one plane): far, v4, x and w within 2e-13 of each other, both 10⁻⁹ or more from
/// the direct sum.
#[test]
fn yukawa_low_degree_chains_are_the_reference_chains() {
    use tetraleaf_ref::Frame;
    use tetraleaf_ref::helmholtz::{Level, Workspace, leaf, rcr};

    let domain = domain();
    let frame = |key: MortonKey| {
        let level = morton::level(key) as usize;
        Frame::new(geometry::centre(key, &domain), radius(level, &domain))
    };
    let kind = ChargeKind::Real;
    for kappa in [0.5, 20.0, 60.0, 100.0] {
        let k = Wavenumber::new(0.0, kappa).unwrap();
        let kc = C::new(0.0, kappa);
        let level = |l: usize| Level::new(kc, radius(l, &domain));
        let degrees = [0, 0, 8, 6, 5];
        let len = |l: usize| (degrees[l] + 1) * (degrees[l] + 1);
        let densities = Densities::new(k.kind(), kind, 1);
        let op = HelmholtzOperator::new(k, densities, &domain, &degrees, true, 64).unwrap();
        let mut pb = Problem {
            k,
            kind,
            densities,
            op,
        };
        let ours = chains(&mut pb, 5, Draw::all(kind, 1));
        let mut rng = Rng(5);
        let mut ws = Workspace::new(8);
        for chain in ours.iter().filter(|c| c.name != "near") {
            let (s, t) = match chain.name {
                "far" => (key([1, 1, 1], 4), key([12, 13, 11], 4)),
                "v4" => (key([1, 1, 1], 4), key([3, 1, 1], 4)),
                "x" => (key([0, 0, 0], 3), key([3, 0, 0], 4)),
                _ => (key([3, 0, 0], 4), key([0, 0, 0], 3)),
            };
            let sources = points_in(s, 20, &mut rng);
            let targets = points_in(t, 15, &mut rng);
            let q = Draw::all(kind, 1).charges(sources.len(), &mut rng)[0][0].clone();
            let (lt, ls) = (morton::level(t) as usize, morton::level(s) as usize);
            let (mut phi, mut grad) = (vec![0.0; targets.len()], vec![[0.0; 3]; targets.len()]);
            let local_of = |ws: &mut Workspace<f64>| -> Vec<f64> {
                match chain.name {
                    "far" => {
                        let up = [s, morton::parent(s).unwrap(), key([0, 0, 0], 2)];
                        let down = [key([3, 3, 2], 2), morton::parent(t).unwrap(), t];
                        let mut m = vec![0.0; len(4)];
                        leaf::p2m_yukawa(5, &level(4), &frame(s), &sources, &q, ws, &mut m);
                        for i in 0..2 {
                            let (lc, lp) = (4 - i, 3 - i);
                            let table = rcr::m2m_table(degrees[lc], degrees[lp], &level(lc));
                            let mut parent = vec![0.0; len(lp)];
                            rcr::m2m_yukawa(
                                degrees[lc],
                                degrees[lp],
                                &level(lc),
                                &level(lp),
                                &frame(up[i]),
                                &frame(up[i + 1]),
                                &table,
                                ws,
                                &m,
                                &mut parent,
                            );
                            m = parent;
                        }
                        let table = rcr::m2l_table(8, 8, &level(2), &frame(up[2]), &frame(down[0]));
                        let mut l = vec![0.0; len(2)];
                        rcr::m2l_yukawa(
                            8,
                            8,
                            &level(2),
                            &frame(up[2]),
                            &frame(down[0]),
                            &table,
                            ws,
                            &m,
                            &mut l,
                        );
                        for i in 0..2 {
                            let (lp, lc) = (2 + i, 3 + i);
                            let table = rcr::l2l_table(degrees[lp], degrees[lc], &level(lc));
                            let mut child = vec![0.0; len(lc)];
                            rcr::l2l_yukawa(
                                degrees[lp],
                                degrees[lc],
                                &level(lp),
                                &level(lc),
                                &frame(down[i]),
                                &frame(down[i + 1]),
                                &table,
                                ws,
                                &l,
                                &mut child,
                            );
                            l = child;
                        }
                        l
                    }
                    "v4" => {
                        let mut m = vec![0.0; 36];
                        leaf::p2m_yukawa(5, &level(4), &frame(s), &sources, &q, ws, &mut m);
                        let table = rcr::m2l_table(5, 5, &level(4), &frame(s), &frame(t));
                        let mut l = vec![0.0; 36];
                        rcr::m2l_yukawa(
                            5,
                            5,
                            &level(4),
                            &frame(s),
                            &frame(t),
                            &table,
                            ws,
                            &m,
                            &mut l,
                        );
                        l
                    }
                    _ => {
                        let mut l = vec![0.0; 36];
                        leaf::p2l_yukawa(5, &level(4), &frame(t), &sources, &q, ws, &mut l);
                        l
                    }
                }
            };
            if chain.name == "w" {
                let mut m = vec![0.0; 36];
                leaf::p2m_yukawa(5, &level(ls), &frame(s), &sources, &q, &mut ws, &mut m);
                leaf::m2p_yukawa(
                    5,
                    &level(ls),
                    &frame(s),
                    &m,
                    &targets,
                    &mut ws,
                    &mut phi,
                    Some(&mut grad),
                );
            } else {
                let l = local_of(&mut ws);
                leaf::l2p_yukawa(
                    5,
                    &level(lt),
                    &frame(t),
                    &l,
                    &targets,
                    &mut ws,
                    &mut phi,
                    Some(&mut grad),
                );
            }
            // φ̃ = 4πφ: the reference's fields carry the factor 4π.
            let scale = 1.0 / (4.0 * std::f64::consts::PI);
            let phi: Vec<C> = phi.iter().map(|&x| C::new(x * scale, 0.0)).collect();
            let grad: Vec<C> = grad
                .iter()
                .flatten()
                .map(|&x| C::new(x * scale, 0.0))
                .collect();
            let (our_phi, our_grad) = pb.physical(chain.target, &chain.output, chain.n, 0);
            let our_grad: Vec<C> = our_grad.iter().flatten().copied().collect();
            let (dp, dg) = (relative_l2(&our_phi, &phi), relative_l2(&our_grad, &grad));
            let (tp, tg) = chain.errors;
            println!(
                "Yukawa κ = {kappa}, {}: against the reference φ {dp:.1e}, ∇φ {dg:.1e}; against the direct sum φ {tp:.1e}, ∇φ {tg:.1e}",
                chain.name
            );
            assert!(dp <= 2e-13 && dg <= 2e-13, "{}: {dp:e}, {dg:e}", chain.name);
            assert!(
                tp > 1e-9,
                "{}: truncation not exercised ({tp:e})",
                chain.name
            );
        }
    }
}

/// Yukawa's P2P with `P2pChoice::Reference` and `Isa::Scalar` against
/// `tetraleaf_ref::helmholtz::p2p::p2p_yukawa` bit for bit, per real column (two complex
/// vectors: four columns), as `p2p_is_the_reference_bit_for_bit`.
#[test]
fn yukawa_p2p_is_the_reference_bit_for_bit() {
    use tetraleaf_exec::operator::{Isa, P2pChoice};
    use tetraleaf_ref::helmholtz::p2p;

    for choice in [P2pChoice::Reference, P2pChoice::Isa(Isa::Scalar)] {
        let kappa = 20.0;
        let k = Wavenumber::new(0.0, kappa).unwrap();
        let mut pb = Problem::new(k, ChargeKind::Complex, 2, 1e-6);
        pb.op = pb.op.with_p2p(choice).unwrap();
        let (s, t) = (key([2, 1, 1], 4), key([1, 1, 1], 4));
        let mut rng = Rng(9);
        let (other, own) = (points_in(s, 10, &mut rng), points_in(t, 8, &mut rng));
        let q_other = Draw::all(ChargeKind::Complex, 2).charges(other.len(), &mut rng);
        let q_own = Draw::all(ChargeKind::Complex, 2).charges(own.len(), &mut rng);
        let other_chunk = source_chunk(&pb.densities, s, &other, &as_charges(pb.kind, &q_other));
        let own_chunk = source_chunk(&pb.densities, t, &own, &as_charges(pb.kind, &q_own));
        let input = target_input(t, &own);
        let mut output = pb.output_chunk(own.len());
        pb.op.p2p(s, t, &other_chunk, &input, &mut output);
        pb.op.p2p(t, t, &own_chunk, &input, &mut output);

        let kappa_hat = kappa * radius(4, &domain());
        let frame = geometry::relative_frame::<f64>(s, t);
        let mapped: Vec<[f64; 3]> = other_chunk[..30]
            .as_chunks::<3>()
            .0
            .iter()
            .map(|u| core::array::from_fn(|i| frame.centre[i] + frame.radius * u[i]))
            .collect();
        let own_units: Vec<[f64; 3]> = own_chunk[..24].as_chunks::<3>().0.to_vec();
        let n = own.len();
        let d = &pb.densities;
        for (column, (v, part)) in [(0, 0), (0, 1), (1, 0), (1, 1)].into_iter().enumerate() {
            let (mut phi, mut grad) = (vec![0.0; n], vec![[0.0; 3]; n]);
            p2p::p2p_yukawa(
                kappa_hat,
                &mapped,
                &q_other[v][part],
                &own_units,
                &mut phi,
                Some(&mut grad),
            );
            p2p::p2p_yukawa(
                kappa_hat,
                &own_units,
                &q_own[v][part],
                &own_units,
                &mut phi,
                Some(&mut grad),
            );
            assert_eq!(d.output_column(v, part), column);
            for i in 0..n {
                assert_eq!(
                    output[d.potential_index(n, column, i)],
                    phi[i],
                    "{choice}, column {column}, target {i}"
                );
                for axis in 0..3 {
                    assert_eq!(
                        output[d.gradient_index(n, column, i, axis)],
                        grad[i][axis],
                        "{choice}, column {column}, target {i}"
                    );
                }
            }
        }
    }
}
