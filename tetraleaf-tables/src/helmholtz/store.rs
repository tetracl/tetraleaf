//! A store of the RCR tables' parts, keyed by what each part depends on, shared by the builds
//! of one process (Helmholtz CH2.2a; docs/design/helmholtz-fmm-proposal.md §6.3, "Table cache";
//! §8 R10).
//!
//! # Keys
//!
//! [`HelmholtzTables`](super::HelmholtzTables) is built in parts (its module documentation,
//! "Parts"), each a function of a few inputs alone. A [`TableStore`] keeps the parts under a
//! key of exactly those inputs:
//!
//! | Part | Key |
//! | --- | --- |
//! | M2M from a child level into its parent | k̂ of the child (the bits of its real and imaginary parts), the child's degree p_{l+1}, the parent's q_l |
//! | L2L from a parent level into a child | k̂ of the child, the parent's q_l, the child's p_{l+1} |
//! | M2L coaxial table of a level | k̂_l, q_l, the distance class \|d\|² (an integer) |
//! | y-blocks, the octants' and M2L's alike | the degree, the polar angle's integer key (sgn d_z, d_z²/\|d\|² reduced) |
//!
//! and the store's value type T. Whether a step is real (real k, Yukawa) follows from k̂ (Re or
//! Im zero), so the key holds it. A part taken from the store is the part a build computes from
//! the same inputs, so tables built with a store are those built without one, bit for bit,
//! whatever the store holds.
//!
//! Level l's tables depend on its level only through k̂_l = k̂₀ / 2^l, exact in floating point.
//! On one domain k̂₀(2k) = 2 k̂₀(k) exactly, so level l + 1 at 2k has level l's k̂ at k: across a
//! sweep of k in octaves every level's parts are reused but those of the new coarse levels
//! (proposal §6.3), where the degrees are a function of k̂_l and the tolerance (the rule of
//! `tetraleaf-exec`'s `Degrees`). The y-blocks depend on no k̂ at all: every level and every k of
//! one degree shares them.
//!
//! # Memory and eviction
//!
//! A part is held by an [`Arc`]: the tables of every build that took it share its memory with
//! the store. The store holds at most its cap in bytes (the parts' values, as
//! [`HelmholtzTables::bytes`](super::HelmholtzTables::bytes) counts them): after each build it
//! drops the parts used least recently (by build, ties in key order, so the same parts on every
//! run), until it fits. Dropping a part drops only the store's reference; tables that hold it
//! keep it. With a cap of 0 the store keeps nothing, and every build computes each of its distinct
//! parts once.
//!
//! # Threads and processes
//!
//! The store is behind a [`Mutex`], held only to look parts up and to put them in, never while
//! a part is computed, so builds on several threads may share one; two builds that compute the
//! same part keep the first one stored (the same bits). The store lives in memory, in one
//! process: each rank of an MPI run has its own.

use std::collections::BTreeMap;
use std::sync::{Arc, Mutex, MutexGuard};
use std::time::Duration;

use tetraleaf_math::RealScalar;

use super::Shift;

/// The key of a polar-angle class: sgn d_z and the reduced fraction d_z² / |d|² (the module
/// documentation of [`helmholtz`](super), "Shifts and their classes").
pub(super) type PolarKey = (i64, i64, i64);

/// What a part depends on (module documentation, "Keys").
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(super) enum PartKey {
    /// M2M from a child level: k̂ of the child as bits, degrees in and out.
    Up {
        k_hat: [u64; 2],
        p_in: usize,
        p_out: usize,
    },
    /// L2L into a child level: k̂ of the child as bits, degrees in and out.
    Down {
        k_hat: [u64; 2],
        p_in: usize,
        p_out: usize,
    },
    /// An M2L coaxial table: k̂ of the level as bits, its degree, |d|².
    M2l {
        k_hat: [u64; 2],
        degree: usize,
        distance: i64,
    },
    /// The y-blocks of one polar angle at one degree.
    Blocks { degree: usize, polar: PolarKey },
}

/// A stored part.
#[derive(Clone, Debug)]
pub(super) enum Stored<T> {
    Shift(Arc<Shift<T>>),
    Blocks(Arc<Vec<T>>),
}

impl<T> Stored<T> {
    /// The coaxial step of a shift's part.
    pub(super) fn shift(self) -> Arc<Shift<T>> {
        match self {
            Stored::Shift(shift) => shift,
            Stored::Blocks(_) => unreachable!("a shift's part"),
        }
    }

    /// The y-blocks of a blocks' part.
    pub(super) fn blocks(self) -> Arc<Vec<T>> {
        match self {
            Stored::Blocks(blocks) => blocks,
            Stored::Shift(_) => unreachable!("a blocks' part"),
        }
    }

    /// The values the part holds, as [`HelmholtzTables::bytes`](super::HelmholtzTables::bytes)
    /// counts them.
    pub(super) fn values(&self) -> usize {
        match self {
            Stored::Shift(shift) => shift.re.len() + shift.im.len(),
            Stored::Blocks(blocks) => blocks.len(),
        }
    }
}

/// A part in the store, with the build that used it last.
#[derive(Debug)]
struct Entry<T> {
    part: Stored<T>,
    bytes: usize,
    used: u64,
}

#[derive(Debug)]
struct Inner<T> {
    entries: BTreeMap<PartKey, Entry<T>>,
    bytes: usize,
    /// The builds begun so far, the clock of the least-recently-used order.
    builds: u64,
}

/// The parts of the RCR tables of the builds of one process (module documentation), within a
/// byte cap.
#[derive(Debug)]
pub struct TableStore<T: RealScalar> {
    cap: usize,
    inner: Mutex<Inner<T>>,
}

impl<T: RealScalar> TableStore<T> {
    /// An empty store that holds at most `cap` bytes of parts after each build (module
    /// documentation, "Memory and eviction").
    pub fn new(cap: usize) -> Self {
        Self {
            cap,
            inner: Mutex::new(Inner {
                entries: BTreeMap::new(),
                bytes: 0,
                builds: 0,
            }),
        }
    }

    /// An empty store without a cap: it keeps every part until [`clear`](Self::clear).
    pub fn unbounded() -> Self {
        Self::new(usize::MAX)
    }

    /// The cap in bytes.
    pub fn cap(&self) -> usize {
        self.cap
    }

    /// The bytes of the parts the store holds, each counted once.
    pub fn bytes(&self) -> usize {
        self.lock().bytes
    }

    /// The number of parts the store holds.
    pub fn entries(&self) -> usize {
        self.lock().entries.len()
    }

    /// Drops every part (tables that hold one keep it).
    pub fn clear(&self) {
        let mut inner = self.lock();
        inner.entries.clear();
        inner.bytes = 0;
    }

    fn lock(&self) -> MutexGuard<'_, Inner<T>> {
        // A build that panicked held the lock only between whole updates of the map.
        self.inner
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
    }

    /// Begins a build: its number in the least-recently-used order.
    pub(super) fn begin(&self) -> u64 {
        let mut inner = self.lock();
        inner.builds += 1;
        inner.builds
    }

    /// The parts of `keys` the store holds, each marked as used by build `build`.
    pub(super) fn get(&self, keys: &[PartKey], build: u64) -> Vec<Option<Stored<T>>> {
        let mut inner = self.lock();
        keys.iter()
            .map(|key| {
                inner.entries.get_mut(key).map(|entry| {
                    entry.used = build;
                    entry.part.clone()
                })
            })
            .collect()
    }

    /// Puts the parts that build `build` computed in the store (a key already there keeps its
    /// part, the same bits), then drops the least recently used until the store fits its cap.
    pub(super) fn insert(&self, parts: impl IntoIterator<Item = (PartKey, Stored<T>)>, build: u64) {
        let mut inner = self.lock();
        for (key, part) in parts {
            if let Some(entry) = inner.entries.get_mut(&key) {
                entry.used = entry.used.max(build);
                continue;
            }
            let bytes = part.values() * size_of::<T>();
            inner.bytes += bytes;
            inner.entries.insert(
                key,
                Entry {
                    part,
                    bytes,
                    used: build,
                },
            );
        }
        if inner.bytes > self.cap {
            let mut order: Vec<(u64, PartKey)> =
                inner.entries.iter().map(|(k, e)| (e.used, *k)).collect();
            order.sort_unstable();
            for (_, key) in order {
                if inner.bytes <= self.cap {
                    break;
                }
                let entry = inner.entries.remove(&key).expect("a stored key");
                inner.bytes -= entry.bytes;
            }
        }
    }
}

/// One level's share of a build of [`HelmholtzTables`](super::HelmholtzTables): the parts its
/// tables hold (the M2M, L2L and octant blocks of the pair of levels below it, its M2L tables
/// and stored M2L blocks), how many the build computed and how many it took from a store or
/// from an earlier part of the same build, the computed parts' time, and the level's bytes.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct LevelBuild {
    /// The parts computed.
    pub built: usize,
    /// The parts taken from the store, or from a part of the same key computed for another
    /// level of the build.
    pub reused: usize,
    /// The computed parts' time, each part's own clock summed (thread-seconds where the parts
    /// ran on several threads).
    pub seconds: Duration,
    /// The level's bytes, its shared parts counted ([`HelmholtzTables::level_bytes`](
    /// super::HelmholtzTables::level_bytes)).
    pub bytes: usize,
}

/// A build of [`HelmholtzTables`](super::HelmholtzTables) per level, and the store's size after
/// it ([`HelmholtzTables::build_report`](super::HelmholtzTables::build_report)).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct TableBuild {
    /// Per level, as [`LevelBuild`] says.
    pub levels: Vec<LevelBuild>,
    /// The store's bytes after the build, or `None` for a build without one.
    pub store_bytes: Option<usize>,
}

impl TableBuild {
    /// The parts computed over every level.
    pub fn built(&self) -> usize {
        self.levels.iter().map(|l| l.built).sum()
    }

    /// The parts reused over every level.
    pub fn reused(&self) -> usize {
        self.levels.iter().map(|l| l.reused).sum()
    }

    /// The computed parts' time over every level.
    pub fn seconds(&self) -> Duration {
        self.levels.iter().map(|l| l.seconds).sum()
    }
}
