pub mod alidisplay;
pub mod alphabet;
pub mod api;
pub mod bias;
pub mod domtblout;
pub mod dp;
pub mod forward;
pub mod hmmfile;
pub mod msa;
pub mod msv;
pub mod mpisupport;
pub mod ncbi;
pub mod omx;
pub(crate) mod optacc;
pub mod pli;
pub mod pfamtblout;
pub mod pipeline;
pub mod report;
pub mod seqio;
pub mod ssi;
pub mod tblout;
pub mod tracealign;
pub mod vitfilter;

mod easel;

/// Put the FPU in the same mode HMMER runs in.
///
/// C, impl_sse/impl_sse.h:355-374 (`impl_Init()`, called by `p7_Init()` at the top
/// of every HMMER program):
/// ```text
///   static inline void
///   impl_Init(void)
///   {
///   #ifdef HAVE_FLUSH_ZERO_MODE
///     /* In order to avoid the performance penalty dealing with sub-normal
///      * values in the floating point calculations, set the processor flag
///      * so sub-normals are "flushed" immediately to zero.
///      */
///     _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
///   #endif
///   #ifdef _PMMINTRIN_H_INCLUDED
///     /*
///      * FLUSH_ZERO doesn't necessarily work in non-SIMD calculations
///      * (yes on 64-bit, maybe not of 32-bit). This ensures that those
///      * scalar calculations will agree across architectures.
///      */
///     _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
///   #endif
///   }
/// ```
///
/// This is part of the arithmetic, not a tuning knob: with sub-normals flushed, an
/// `exp(lnP)` that underflows yields exactly `0.0`, which is why C prints `0` for
/// E-values below ~1e-308 where an IEEE-conforming run prints e.g. `1.3e-309`.
///
/// `MXCSR` is per-thread, so every worker thread must call this too. The CLI does it
/// for its rayon pool via `start_handler`; `HmmAnnotator::search` does it at the top
/// of each per-model task.
#[inline]
pub fn init() {
    #[cfg(target_arch = "x86_64")]
    {
        // Rust does not expose `_MM_SET_DENORMALS_ZERO_MODE` and has deprecated the
        // FTZ setter, so set the two MXCSR control bits directly:
        //   bit 15 (0x8000) = FTZ, flush-to-zero
        //   bit  6 (0x0040) = DAZ, denormals-are-zero
        #[allow(deprecated)]
        // SAFETY: `_mm_getcsr`/`_mm_setcsr` only read and write the current thread's
        // MXCSR; setting FTZ/DAZ is exactly what C's impl_Init() does.
        unsafe {
            use core::arch::x86_64::{_mm_getcsr, _mm_setcsr};
            _mm_setcsr(_mm_getcsr() | 0x8000 | 0x0040);
        }
    }
}
