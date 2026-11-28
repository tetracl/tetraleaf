//! The dense M2L tables of the Helmholtz family (Helmholtz CH3.2a;
//! docs/design/helmholtz-fmm-proposal.md §5.3, Dense): per level and V-list offset the
//! complex M2L as one real matrix, built in f64 from the level's RCR tables and rounded to T.
//!
//! # Layout
//!
//! The M2L of offset d on level l, of degree p = p_l and n = (p + 1)² coefficients, is a
//! complex n × n matrix A = A_r + i A_i from the multipole's coefficients to the local's
//! (CONVENTIONS §3.6 order in both). [`DenseM2l`] stores it as one 2n × n real matrix,
//! column-major, A_r in rows 0..n over A_i in rows n..2n: column j is the local of the unit
//! multipole e_j in two planes, real part then imaginary part, the layout of a box's
//! complex local (CONVENTIONS §3.15, "Storage"). The 316 matrices of a level lie one after
//! the other in offset order, 316 · 2n² values a level (21 MB at p = 7 in f64, 331 MB at
//! p = 15, 983 MB at p = 20).
//!
//! For Yukawa (k̂₀ = iκ̂₀, Helmholtz CH2.1a) the M2L is real, from a real multipole to a real
//! local (`HelmholtzTables`' "Yukawa"): one real n × n matrix A per offset, column j the
//! local of e_j in one plane, 316 n² values a level ([`DenseM2l::local_planes`] 1).
//!
//! # Build
//!
//! [`DenseM2l::build`] builds [`HelmholtzTables`] in f64 for the same k̂₀ and degrees, the
//! M2L tables of the chosen levels alone (since Helmholtz CH4.0h: no translations, no other
//! level's; the M2L of a level reads only its own tables), and takes column j of offset d on level l as their M2L of e_j (one real plane) into zeros:
//! the matrices hold the stored form of `HelmholtzTables` (its "The stored form, folded
//! in"), Rotation's coefficients applied to unit vectors, in f64, rounded to T entry by
//! entry. So Dense agrees with Rotation to the rounding of the two orders of operations,
//! and with `tetraleaf_ref::helmholtz::rcr` as Rotation does. The build is deterministic; it
//! costs n Rotation applications per offset, O(p⁵). [`DenseM2l::build_on`] computes its
//! parts (the f64 tables', then the matrices of a few offsets at a time) with a caller's
//! [`Parts`], the same bits (Helmholtz CH4.0h).
//!
//! # Levels and the budget
//!
//! Only levels ≥ 2 carry a V list. [`DenseM2l::select`] picks the levels whose tables fit a
//! byte budget, from the finest level up, as `HelmholtzTables::build` stores its y-blocks:
//! a level that does not fit what is left is skipped, and a coarser one may still fit. The
//! choice is a function of the degrees and the budget only (proposal §2, constraint 9).
//!
//! # Application
//!
//! [`DenseM2l::m2l`] adds the M2L of `vectors` multipoles to their locals, per vector the
//! one-vector operation, each product and sum rounded apart (as `MatrixSet::apply`):
//!
//! - one plane (the real multipole of real k and real charges, and Yukawa's):
//!   L\[r\] ← L\[r\] + A\[r, k\] M\[k\] for k = 0, 1, …, n − 1 in turn, over the 2n rows
//!   (n for Yukawa);
//! - two planes: in one pass over the table, for k = 0, 1, …, n − 1 in turn,
//!   L_r\[r\] ← (L_r\[r\] + A_r\[r, k\] M_r\[k\]) − A_i\[r, k\] M_i\[k\] and
//!   L_i\[r\] ← (L_i\[r\] + A_i\[r, k\] M_r\[k\]) + A_r\[r, k\] M_i\[k\], over the n rows.
//!
//! So the four real products of a complex M2L are combined pair by pair, with no per-target
//! epilogue (Helmholtz CH3.5, question 2), and each table is read once per pair: per pair
//! the application is bound by reading the table, 2n² values.

use tetraleaf_math::{Complex, Layout, RealScalar};

use super::{DEFAULT_BLOCK_BUDGET, HelmholtzTables, TableChoice};
use crate::geometry::M2L_OFFSET_COUNT;
use crate::parts::{Parts, Serial};

/// The dense M2L tables of the chosen levels of an octree (the
/// [module documentation](self)).
#[derive(Clone, Debug)]
pub struct DenseM2l<T: RealScalar> {
    /// Per level its degree and, if dense, its 316 matrices.
    levels: Vec<(usize, Option<Vec<T>>)>,
    /// The planes of a local: 2, or 1 for Yukawa.
    local_planes: usize,
}

impl<T: RealScalar> DenseM2l<T> {
    /// The bytes of one level's 316 matrices at degree p: 316 · 2 (p + 1)⁴ values of T
    /// ([`level_bytes_for`](Self::level_bytes_for) with two local planes).
    pub fn level_bytes(p: usize) -> usize {
        Self::level_bytes_for(p, 2)
    }

    /// The bytes of one level's 316 matrices at degree p into locals of `local_planes`
    /// planes (2, or 1 for Yukawa; Helmholtz CH2.1a): 316 · `local_planes` (p + 1)⁴ values
    /// of T.
    pub fn level_bytes_for(p: usize, local_planes: usize) -> usize {
        let n = Layout::new(p).len();
        M2L_OFFSET_COUNT * local_planes * n * n * size_of::<T>()
    }

    /// The levels whose tables fit in `budget` bytes, from the finest level up; never levels
    /// 0 and 1 (module documentation, "Levels and the budget"):
    /// [`select_for`](Self::select_for) with two local planes.
    pub fn select(degrees: &[usize], budget: usize) -> Vec<bool> {
        Self::select_for(degrees, budget, 2)
    }

    /// As [`select`](Self::select), for locals of `local_planes` planes (2, or 1 for
    /// Yukawa; Helmholtz CH2.1a), each level weighed by
    /// [`level_bytes_for`](Self::level_bytes_for).
    pub fn select_for(degrees: &[usize], budget: usize, local_planes: usize) -> Vec<bool> {
        let mut chosen = vec![false; degrees.len()];
        let mut left = budget;
        for (l, &p) in degrees.iter().enumerate().skip(2).rev() {
            let bytes = Self::level_bytes_for(p, local_planes);
            if bytes <= left {
                left -= bytes;
                chosen[l] = true;
            }
        }
        chosen
    }

    /// Builds the tables of the levels with `levels[l]` set, for the root's k̂₀ = k r₀ and
    /// the degrees p_l (module documentation, "Build"); with no level set, builds nothing.
    ///
    /// # Panics
    ///
    /// If `levels` and `degrees` differ in length, level 0 or 1 is chosen, or as
    /// [`HelmholtzTables::build`] panics.
    pub fn build(k_root: Complex<f64>, degrees: &[usize], levels: &[bool]) -> Self {
        Self::build_on(k_root, degrees, levels, &Serial)
    }

    /// As [`build`](Self::build), its parts computed by `parts` (module documentation,
    /// "Build"): the same tables bit for bit.
    ///
    /// # Panics
    ///
    /// As [`build`](Self::build).
    pub fn build_on(
        k_root: Complex<f64>,
        degrees: &[usize],
        levels: &[bool],
        parts: &impl Parts,
    ) -> Self {
        assert_eq!(
            levels.len(),
            degrees.len(),
            "one choice per level: {} choices, {} degrees",
            levels.len(),
            degrees.len()
        );
        assert!(
            levels.iter().take(2).all(|&on| !on),
            "no dense M2L on levels 0 and 1: they have no V list"
        );
        // As `HelmholtzTables`: Re k = 0 is Yukawa's real path, one local plane.
        let local_planes = if k_root.re == 0.0 { 1 } else { 2 };
        if !levels.contains(&true) {
            return Self {
                levels: degrees.iter().map(|&p| (p, None)).collect(),
                local_planes,
            };
        }
        let tables = m2l_tables(k_root, degrees, levels, parts);
        debug_assert_eq!(tables.local_planes(), local_planes);
        let all: Vec<usize> = (0..M2L_OFFSET_COUNT).collect();
        let levels = degrees
            .iter()
            .zip(levels)
            .enumerate()
            .map(|(l, (&p, &on))| (p, on.then(|| matrices_on(&tables, l, &all, parts))))
            .collect();
        Self {
            levels,
            local_planes,
        }
    }

    /// The planes of a local: 1 for Yukawa (a real n × n matrix per offset), else 2 (module
    /// documentation, "Layout").
    pub fn local_planes(&self) -> usize {
        self.local_planes
    }

    /// The number of levels.
    pub fn nlevels(&self) -> usize {
        self.levels.len()
    }

    /// The degree p_l of `level`.
    pub fn degree(&self, level: usize) -> usize {
        self.levels[level].0
    }

    /// Whether `level` has dense tables.
    pub fn is_dense(&self, level: usize) -> bool {
        self.levels[level].1.is_some()
    }

    /// The bytes of every level's tables.
    pub fn bytes(&self) -> usize {
        self.levels
            .iter()
            .filter_map(|(_, m)| m.as_ref())
            .map(|m| m.len() * size_of::<T>())
            .sum()
    }

    /// The 2n × n matrix of offset index `d` on `level`, column-major, or Yukawa's n × n
    /// (module documentation, "Layout").
    ///
    /// # Panics
    ///
    /// If `level` has no dense tables or d ≥ 316.
    pub fn matrix(&self, level: usize, d: usize) -> &[T] {
        assert!(d < M2L_OFFSET_COUNT, "offset index {d}");
        let (p, matrices) = &self.levels[level];
        let matrices = matrices
            .as_ref()
            .unwrap_or_else(|| panic!("no dense M2L on level {level}"));
        let size = self.local_planes * Layout::new(*p).len().pow(2);
        &matrices[d * size..(d + 1) * size]
    }

    /// M2L on `level` at the offset index `d`: adds the translations of `vectors`
    /// multipoles of `planes` planes (2, or 1 for the real multipoles of real k and real
    /// charges) to the target's `vectors` complex locals, or for Yukawa of `vectors` real
    /// multipoles to real locals (one plane each; module documentation, "Application").
    ///
    /// # Panics
    ///
    /// If `level` has no dense tables, d ≥ 316, `planes` is not 1 or 2 (1 for Yukawa), or a
    /// slice is not `vectors` expansions at the level's degree.
    pub fn m2l(
        &self,
        level: usize,
        d: usize,
        multipole: &[T],
        planes: usize,
        local: &mut [T],
        vectors: usize,
    ) {
        let a = self.matrix(level, d);
        assert!(matches!(planes, 1 | 2), "{planes} multipole planes");
        assert!(
            self.local_planes == 2 || planes == 1,
            "Yukawa's M2L takes one multipole plane, got {planes}"
        );
        let n = Layout::new(self.levels[level].0).len();
        let rows = self.local_planes * n;
        assert_eq!(multipole.len(), vectors * planes * n, "multipole length");
        assert_eq!(local.len(), vectors * rows, "local length");
        for (m, l) in multipole
            .chunks_exact(planes * n)
            .zip(local.chunks_exact_mut(rows))
        {
            apply(a, n, planes, m, l);
        }
    }
}

/// The f64 tables of the M2L of the chosen `levels` alone (no translations, no other
/// level's), the M2L blocks within [`DEFAULT_BLOCK_BUDGET`], built by `parts`: their M2L on
/// a chosen level is that of `HelmholtzTables::build`'s, bit for bit.
pub(super) fn m2l_tables(
    k_root: Complex<f64>,
    degrees: &[usize],
    levels: &[bool],
    parts: &impl Parts,
) -> HelmholtzTables<f64> {
    let choice = TableChoice {
        translations: false,
        m2l_levels: Some(levels),
        ..TableChoice::new(DEFAULT_BLOCK_BUDGET)
    };
    HelmholtzTables::build_by(k_root, degrees, &choice, parts)
}

/// The offsets of one part of [`matrices_on`].
const OFFSETS_PER_PART: usize = 8;

/// [`matrices`] of `offsets` on `level`, [`OFFSETS_PER_PART`] offsets per part of `parts`,
/// put together in the offsets' order: the same values.
pub(super) fn matrices_on<T: RealScalar>(
    tables: &HelmholtzTables<f64>,
    level: usize,
    offsets: &[usize],
    parts: &impl Parts,
) -> Vec<T> {
    let chunks: Vec<&[usize]> = offsets.chunks(OFFSETS_PER_PART).collect();
    parts
        .map(chunks.len(), &|i| {
            matrices::<T>(tables, level, chunks[i].iter().copied())
        })
        .concat()
}

/// The matrices of `offsets` on `level`, one after the other: column j of each Rotation's
/// M2L of e_j (one real plane) into zeros, in f64, rounded to T (module documentation,
/// "Build"); a column has the local's planes (two, or Yukawa's one).
pub(super) fn matrices<T: RealScalar>(
    tables: &HelmholtzTables<f64>,
    level: usize,
    offsets: impl ExactSizeIterator<Item = usize>,
) -> Vec<T> {
    let n = Layout::new(tables.degree(level)).len();
    let mut scratch = tables.scratch();
    let mut identity = vec![0.0; n * n];
    for j in 0..n {
        identity[j * n + j] = 1.0;
    }
    let size = tables.local_planes() * n * n;
    let mut columns = vec![0.0; size];
    let mut values = Vec::with_capacity(offsets.len() * size);
    for d in offsets {
        // Vector j of n is e_j; its local is column j, of the local's planes.
        columns.fill(0.0);
        tables.m2l(level, d, &identity, 1, &mut columns, n, &mut scratch);
        values.extend(columns.iter().map(|&v| T::from_f64(v)));
    }
    values
}

/// Adds the M2L of one multipole `m` of `planes` planes by the 2n × n matrix `a` to one
/// complex local `l` (module documentation, "Application"); with one plane, by a matrix of
/// `l.len()` rows (Yukawa's n × n into one real plane).
#[inline]
pub(super) fn apply<T: RealScalar>(a: &[T], n: usize, planes: usize, m: &[T], l: &mut [T]) {
    if planes == 1 {
        for (column, &x) in a.chunks_exact(l.len()).zip(m) {
            for (y, &e) in l.iter_mut().zip(column) {
                *y = *y + e * x;
            }
        }
        return;
    }
    let (m_r, m_i) = m.split_at(n);
    let (l_r, l_i) = l.split_at_mut(n);
    for ((column, &x_r), &x_i) in a.chunks_exact(2 * n).zip(m_r).zip(m_i) {
        let (a_r, a_i) = column.split_at(n);
        for ((y_r, y_i), (&e_r, &e_i)) in
            l_r.iter_mut().zip(l_i.iter_mut()).zip(a_r.iter().zip(a_i))
        {
            *y_r = *y_r + e_r * x_r - e_i * x_i;
            *y_i = *y_i + e_i * x_r + e_r * x_i;
        }
    }
}
