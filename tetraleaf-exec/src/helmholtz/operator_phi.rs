//! The φ step of the diagonal form's conversions (Helmholtz CH4.6a1; proposal §6.7, D7).
//!
//! S → F and F → R are a Legendre step in θ and a φ step along each ring
//! (`tetraleaf_tables::helmholtz::diagonal`). The product (`Conversions::{s2f, f2r}`) does both
//! as products over a block of boxes; it is the default, always built, and the reference. With
//! the optional feature `rustfft` the φ step can be an FFT instead
//! ([`HelmholtzPhiStep::RustFft`]): `Conversions::s2f_modes` stops S → F after the Legendre step
//! with each ring's Fourier modes laid out as the input of an inverse DFT of length N_φ, written
//! into the box's own signature, which the level's plan transforms in place; F → R transforms an
//! accumulator in place with the forward plan and `Conversions::f2r_modes` takes it from there.
//! The two paths agree to rounding, not bit for bit.
//!
//! # The back ends
//!
//! An FFT back end is a [`PhiTransform`], the trait internal to this crate that proposal §6.7
//! asks for (RLST provides none): an unnormalised inverse and forward DFT of length N_φ run in
//! place over the 2p + 1 rings of one box with a thread's scratch. rustfft is the one back end
//! (D7, signed off at CH0.9); CH4.6b's FFTW would be a second.
//!
//! # Determinism
//!
//! One plan per (level, N_φ, direction), made once at the operator's build by a fresh planner,
//! executed box by box on the box's own samples by whichever worker holds the box: the plan never
//! sees the batch, so a box's bits do not depend on its position, the batch size or the thread
//! (spike S9 measured this on the M6, locust and Kathleen, `spikes/fft-plans/SPIKE_REPORT.md` §3,
//! and `tests/helmholtz/phi.rs` checks it). rustfft's planner chooses its algorithms by the
//! length and the CPU features it detects (AVX with FMA, or SSE4.1, on x86_64; NEON on aarch64),
//! so two ranks of one build on CPUs of different features may get different bits. The builder
//! therefore agrees each level's [`fingerprint`] (rustfft's version, the architecture, the
//! dispatched features, N_φ and the precision) over the ranks, and where they differ every rank
//! runs that level's product (`HelmholtzFmmBuilder::build`; constraints 2 and 9).

use tetraleaf_math::Complex;
use tetraleaf_simd::SimdScalar;

/// The φ step of the diagonal form's S → F and F → R (Helmholtz CH4.6a1; proposal §6.7, D7;
/// `HelmholtzOperator::with_phi_step`, `HelmholtzFmmBuilder::phi_step`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum HelmholtzPhiStep {
    /// The transform written as a product (`Conversions::{s2f, f2r}`): the default and the
    /// reference, always built.
    #[default]
    Product,
    /// rustfft's FFT, per box and ring (feature `rustfft`; off by default): the product to
    /// rounding, faster on the timed machines' high-frequency levels (S9). Where the ranks'
    /// rustfft builds differ (the build's agreement), a level runs the product.
    #[cfg(feature = "rustfft")]
    RustFft,
}

/// One level's φ step by FFT: an unnormalised inverse DFT (synthesis, e^{+2πi mJ/N}) and forward
/// DFT (analysis, e^{−2πi mJ/N}) of the level's length N_φ, each run in place over whole rings
/// (a multiple of N_φ values) with a thread's scratch of [`scratch_len`](Self::scratch_len).
pub(super) trait PhiTransform<T>: Send + Sync {
    /// The inverse DFT of every ring of `rings`, in place.
    fn synthesis(&self, rings: &mut [Complex<T>], scratch: &mut [Complex<T>]);
    /// The forward DFT of every ring of `rings`, in place.
    fn analysis(&self, rings: &mut [Complex<T>], scratch: &mut [Complex<T>]);
    /// The scratch either direction needs.
    fn scratch_len(&self) -> usize;
}

/// The level's φ step by FFT under `step` at length N_φ = `azimuths`: `None` for the product.
pub(super) fn transform<T: SimdScalar>(
    step: HelmholtzPhiStep,
    azimuths: usize,
) -> Option<Box<dyn PhiTransform<T>>> {
    match step {
        HelmholtzPhiStep::Product => {
            let _ = azimuths;
            None
        }
        #[cfg(feature = "rustfft")]
        HelmholtzPhiStep::RustFft => Some(rust::plans::<T>(azimuths)),
    }
}

/// The version of rustfft the feature builds, the workspace's exact pin (`Cargo.toml`; a unit
/// test holds the two together).
pub const RUSTFFT_VERSION: &str = "6.4.1";

/// The CPU features rustfft 6.4 dispatches on, detected at run time: on x86_64 `avx` and `fma`
/// (its AVX algorithms need both) and `sse4.1`, on aarch64 `neon`; none elsewhere.
pub fn rustfft_cpu_features() -> Vec<&'static str> {
    #[allow(unused_mut)]
    let mut found = Vec::new();
    #[cfg(target_arch = "x86_64")]
    {
        if std::arch::is_x86_feature_detected!("sse4.1") {
            found.push("sse4.1");
        }
        if std::arch::is_x86_feature_detected!("avx") {
            found.push("avx");
        }
        if std::arch::is_x86_feature_detected!("fma") {
            found.push("fma");
        }
    }
    #[cfg(target_arch = "aarch64")]
    {
        if std::arch::is_aarch64_feature_detected!("neon") {
            found.push("neon");
        }
    }
    found
}

/// A level's fingerprint under `step` (0 for the product): FNV-1a of the back end's name,
/// rustfft's version, the architecture, the CPU features it dispatches on, N_φ = `azimuths` and
/// the precision's `value_bytes`, xor `salt` (a test hook; 0 in every real build). Two ranks with equal
/// fingerprints make the same plans, so a box's bits agree across them.
pub(super) fn fingerprint(
    step: HelmholtzPhiStep,
    azimuths: usize,
    value_bytes: usize,
    salt: u64,
) -> u64 {
    match step {
        HelmholtzPhiStep::Product => {
            let _ = (azimuths, value_bytes, salt);
            0
        }
        #[cfg(feature = "rustfft")]
        HelmholtzPhiStep::RustFft => {
            let mut text = format!(
                "rustfft {RUSTFFT_VERSION} {} {} {azimuths} {}",
                std::env::consts::ARCH,
                rustfft_cpu_features().join(","),
                value_bytes
            );
            text.push(';');
            fnv1a(text.as_bytes()).max(1) ^ salt
        }
    }
}

/// FNV-1a, 64 bits.
#[cfg_attr(not(feature = "rustfft"), allow(dead_code))]
fn fnv1a(bytes: &[u8]) -> u64 {
    bytes.iter().fold(0xcbf2_9ce4_8422_2325, |h, &b| {
        (h ^ u64::from(b)).wrapping_mul(0x0000_0100_0000_01b3)
    })
}

#[cfg(feature = "rustfft")]
mod rust {
    //! rustfft's plans, one per direction of a level's N_φ.

    use super::PhiTransform;
    use rustfft::{Fft, FftNum, FftPlanner};
    use std::any::{Any, TypeId};
    use std::sync::Arc;
    use tetraleaf_math::Complex;
    use tetraleaf_simd::SimdScalar;

    /// The two plans of one length, from a fresh planner.
    struct Plans<F: FftNum> {
        synthesis: Arc<dyn Fft<F>>,
        analysis: Arc<dyn Fft<F>>,
    }

    impl<F: FftNum> Plans<F> {
        fn new(n: usize) -> Self {
            let mut planner = FftPlanner::<F>::new();
            Self {
                synthesis: planner.plan_fft_inverse(n),
                analysis: planner.plan_fft_forward(n),
            }
        }
    }

    impl<F: FftNum> PhiTransform<F> for Plans<F> {
        fn synthesis(&self, rings: &mut [Complex<F>], scratch: &mut [Complex<F>]) {
            let need = self.synthesis.get_inplace_scratch_len();
            self.synthesis
                .process_with_scratch(rings, &mut scratch[..need]);
        }

        fn analysis(&self, rings: &mut [Complex<F>], scratch: &mut [Complex<F>]) {
            let need = self.analysis.get_inplace_scratch_len();
            self.analysis
                .process_with_scratch(rings, &mut scratch[..need]);
        }

        fn scratch_len(&self) -> usize {
            self.synthesis
                .get_inplace_scratch_len()
                .max(self.analysis.get_inplace_scratch_len())
        }
    }

    /// The plans of length `n` for `T`, f32 or f64 (`SimdScalar` is sealed to the two): made for
    /// the concrete type and handed back through `Any`, so no bound beyond `SimdScalar` is needed.
    pub(super) fn plans<T: SimdScalar>(n: usize) -> Box<dyn PhiTransform<T>> {
        let any: Box<dyn Any> = if TypeId::of::<T>() == TypeId::of::<f64>() {
            Box::new(Box::new(Plans::<f64>::new(n)) as Box<dyn PhiTransform<f64>>)
        } else if TypeId::of::<T>() == TypeId::of::<f32>() {
            Box::new(Box::new(Plans::<f32>::new(n)) as Box<dyn PhiTransform<f32>>)
        } else {
            unreachable!("SimdScalar is f32 or f64")
        };
        *any.downcast::<Box<dyn PhiTransform<T>>>()
            .expect("the plans' type is T's")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The version in the fingerprint is the workspace's pin.
    #[test]
    fn the_version_is_the_pin() {
        let manifest = include_str!("../../../Cargo.toml");
        assert!(
            manifest.contains(&format!("rustfft = \"={RUSTFFT_VERSION}\"")),
            "Cargo.toml's rustfft pin is not {RUSTFFT_VERSION}"
        );
    }

    #[test]
    fn the_products_fingerprint_is_zero() {
        assert_eq!(fingerprint(HelmholtzPhiStep::Product, 580, 8, 0), 0);
        assert_eq!(HelmholtzPhiStep::default(), HelmholtzPhiStep::Product);
    }

    #[cfg(feature = "rustfft")]
    #[test]
    fn the_fingerprint_reads_the_length_the_precision_and_the_salt() {
        let f = |n, salt| fingerprint(HelmholtzPhiStep::RustFft, n, 8, salt);
        assert_ne!(f(580, 0), 0);
        assert_eq!(f(580, 0), f(580, 0));
        assert_ne!(f(580, 0), f(588, 0));
        assert_ne!(f(580, 0), f(580, 1));
        assert_ne!(f(580, 0), fingerprint(HelmholtzPhiStep::RustFft, 580, 4, 0));
    }
}
