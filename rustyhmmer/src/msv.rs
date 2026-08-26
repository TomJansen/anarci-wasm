













use crate::alphabet::{AMINO_FREQ, K, KP};
use crate::hmmfile::P7Hmm;





fn degen_set(code: usize) -> &'static [usize] {
    match code {
        21 => &[11, 2],                                             
        22 => &[7, 9],                                              
        23 => &[13, 3],                                             
        24 => &[8],                                                 
        25 => &[1],                                                 
        26 => &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], 
        _ => &[],
    }
}






#[inline]
fn unbiased_byteify(scale_b: f32, sc: f32) -> u8 {
    let v = -(scale_b * sc).round();
    if v >= 255.0 {
        255
    } else if v <= 0.0 {
        0
    } else {
        v as u8
    }
}






#[inline]
fn biased_byteify(scale_b: f32, bias_b: u8, sc: f32) -> u8 {
    let v = -(scale_b * sc).round();
    if v > (255 - bias_b as i32) as f32 {
        255
    } else {
        (v as i32 + bias_b as i32).clamp(0, 255) as u8
    }
}



#[derive(Debug, Clone)]
pub struct MsvProfile {
    pub m: usize,
    
    pub scale_b: f32,
    
    pub base_b: u8,
    
    pub bias_b: u8,
    
    
    pub tbm_b: u8,
    pub tec_b: u8,
    
    pub rbv: Vec<[u8; KP]>,
    
    
    
    
    
    pub rbv_t: Vec<u8>,
    
    
    
    
    
    pub sbv_t: Vec<i8>,
    /// Striped `om->sbv` for the banded SSV kernel: for residue code `x`, vector `q`
    /// is `sbv[(x*(nqb+P7O_EXTRA_SB) + q)*16 .. +16]`, and lane `z` of vector `q` is
    /// model position `k = q + 1 + z*nqb` (or C's k>M padding). The last
    /// `P7O_EXTRA_SB` vectors repeat the first ones so the kernel can index past
    /// `nqb` without a wrap test.
    pub sbv: Vec<i8>,
    /// `p7O_NQB(M)`.
    pub nqb: usize,
    
    pub mu: f64,
    pub lambda: f64,
}

impl MsvProfile {
    
    pub fn build(hmm: &P7Hmm) -> Self {
        let m = hmm.m;
        
        let scale_b = (3.0_f64 / std::f64::consts::LN_2) as f32; 
        let base_b: u8 = 190;

        
        
        
        let mut msc = vec![[f32::NEG_INFINITY; KP]; m + 1];
        let mut maxsc = 0.0f32;
        for k in 1..=m {
            let mut sc = [f32::NEG_INFINITY; KP];
            for x in 0..K {
                sc[x] = ((hmm.mat[k][x] as f64) / (AMINO_FREQ[x] as f32 as f64)).ln() as f32;
                if sc[x] > maxsc {
                    maxsc = sc[x];
                }
            }
            
            for code in K..KP {
                let set = degen_set(code);
                if set.is_empty() {
                    continue;
                }
                let mut left = 0.0f32;
                let mut right = 0.0f32;
                for &y in set {
                    left += (AMINO_FREQ[y] as f32) * sc[y];
                    right += AMINO_FREQ[y] as f32;
                }
                sc[code] = left / right;
            }
            msc[k] = sc;
        }

        
        let bias_b = unbiased_byteify(scale_b, -maxsc);
        let mut rbv = vec![[255u8; KP]; m + 1];
        for k in 1..=m {
            for x in 0..KP {
                rbv[k][x] = biased_byteify(scale_b, bias_b, msc[k][x]);
            }
        }
        
        let mut rbv_t = vec![255u8; KP * (m + 1)];
        for k in 1..=m {
            for x in 0..KP {
                rbv_t[x * (m + 1) + k] = rbv[k][x];
            }
        }
        
        
        
        
        
        
        let mut sbv_t = vec![0i8; KP * (m + 1)];
        let tmp: u8 = (127i16 + bias_b as i16) as u8; 
        for k in 1..=m {
            for x in 0..KP {
                sbv_t[x * (m + 1) + k] = (tmp.saturating_sub(rbv[k][x]) ^ 127) as i8;
            }
        }

        // C, impl_sse/p7_oprofile.c:754-758:
        //   for (x = 0; x < om->abc->Kp; x++)
        //     {
        //       for (q = 0;  q < nq;            q++) om->sbv[x][q] = _mm_xor_si128(_mm_subs_epu8(tmp, om->rbv[x][q]), tmp2);
        //       for (q = nq; q < nq + p7O_EXTRA_SB; q++) om->sbv[x][q] = om->sbv[x][q % nq];
        //     }
        // with om->rbv striped 16-way: lane z of vector q is model position
        // k = q+1 + z*nq, and 255 (i.e. -infinity) where k > M.
        let nqb = std::cmp::max(2, (m - 1) / 16 + 1);
        let vecs = nqb + P7O_EXTRA_SB;
        let mut sbv = vec![0i8; KP * vecs * 16];
        for x in 0..KP {
            for q in 0..nqb {
                for z in 0..16 {
                    let k = q + 1 + z * nqb;
                    let r = if k <= m { rbv[k][x] } else { 255 };
                    sbv[((x * vecs) + q) * 16 + z] = (tmp.saturating_sub(r) ^ 127) as i8;
                }
            }
            for q in nqb..vecs {
                let src = ((x * vecs) + (q % nqb)) * 16;
                let dst = ((x * vecs) + q) * 16;
                sbv.copy_within(src..src + 16, dst);
            }
        }

        let tbm_b = unbiased_byteify(scale_b, (2.0_f64 / (m as f64 * (m as f64 + 1.0))).ln() as f32);
        let tec_b = unbiased_byteify(scale_b, 0.5_f32.ln());

        MsvProfile {
            m,
            scale_b,
            base_b,
            bias_b,
            tbm_b,
            tec_b,
            rbv,
            rbv_t,
            sbv_t,
            sbv,
            nqb,
            mu: hmm.evparam.msv_mu,
            lambda: hmm.evparam.msv_lambda,
        }
    }
}


const LOG2: f64 = std::f64::consts::LN_2;

/// `MAX_BANDS` — impl_sse/ssvfilter.c:425. "Apparently, two registers are generally
/// used for something else, leaving 14 registers on 64 bit versions".
#[cfg_attr(not(target_arch = "x86_64"), allow(dead_code))] // used only by the x86_64 SSE SSV kernel
const MAX_BANDS: usize = 14;
/// `p7O_EXTRA_SB` — impl_sse/impl_sse.h:28.
const P7O_EXTRA_SB: usize = 17;

impl MsvProfile {
    
    
    
    
    
    
    
    
    
    pub fn msv_score(&self, dsq: &[u8], l: usize) -> f32 {
        
        
        let tjb_b = unbiased_byteify(self.scale_b, (3.0_f64 / (l as f64 + 3.0)).ln() as f32);
        
        let tjbm = tjb_b.wrapping_add(self.tbm_b);

        
        
        
        
        
        
        
        
        let m = self.m;
        let mut prev = vec![0u8; m + 1]; 
        let mut cur = vec![0u8; m + 1]; 
        let mut xj: u8 = 0;
        let mut xb: u8 = self.base_b.saturating_sub(tjbm);
        let bias = self.bias_b;

        for i in 1..=l {
            let res = dsq[i] as usize; 
            let base = res * (m + 1);
            let row = &self.rbv_t[base + 1..base + m + 1]; 
            let prev_km1 = &prev[0..m]; 
            let cur_k = &mut cur[1..m + 1]; 
            
            for ((c, &p), &r) in cur_k.iter_mut().zip(prev_km1.iter()).zip(row.iter()) {
                *c = p.max(xb).saturating_add(bias).saturating_sub(r);
            }
            
            let mut xe: u8 = 0;
            for &c in cur_k.iter() {
                if c > xe {
                    xe = c;
                }
            }
            
            if xe.saturating_add(bias) == 255 {
                return f32::INFINITY;
            }
            xe = xe.saturating_sub(self.tec_b); 
            if xj < xe {
                xj = xe;
            }
            xb = self.base_b.max(xj).saturating_sub(tjbm); 
            std::mem::swap(&mut prev, &mut cur); 
        }

        
        
        ((xj as f32 - tjb_b as f32) - self.base_b as f32) / self.scale_b - 3.0
    }

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    pub fn ssv_score(&self, dsq: &[u8], l: usize) -> Option<f32> {
        #[cfg(target_arch = "x86_64")]
        {
            if std::is_x86_feature_detected!("sse2") {
                // SAFETY: guarded by the sse2 feature check just above.
                return unsafe { ssv_score_striped(self, dsq, l) };
            }
        }
        self.ssv_score_scalar(dsq, l)
    }

    /// Scalar oracle for `ssv_score`, and the fallback off x86-64.
    pub fn ssv_score_scalar(&self, dsq: &[u8], l: usize) -> Option<f32> {
        
        
        let tjb_b = self.ssv_setup(l)?;

        let m = self.m;
        let mut prev = vec![-128i8; m + 1]; 
        let mut cur = vec![-128i8; m + 1]; 
        let mut raw_xe: u8 = 128; 

        for i in 1..=l {
            let res = dsq[i] as usize;
            let base = res * (m + 1);
            let srow = &self.sbv_t[base + 1..base + m + 1]; 
            let prev_km1 = &prev[0..m];
            let cur_k = &mut cur[1..m + 1];
            // Split the row update from the running-max reduction. Interleaving
            // the `if u > raw_xe` reduction with the store defeats auto-vectorization
            // (loop-carried scalar dependency); as two separate passes each loop
            // auto-vectorizes (psubsb / pmaxub) — the same shape msv_score already
            // uses. `raw_xe` is read only after the i-loop and max is associative,
            // so the result is bit-identical.
            for ((c, &p), &s) in cur_k.iter_mut().zip(prev_km1.iter()).zip(srow.iter()) {
                *c = p.saturating_sub(s);
            }
            for &c in cur_k.iter() {
                let u = c as u8;
                if u > raw_xe {
                    raw_xe = u;
                }
            }
            std::mem::swap(&mut prev, &mut cur);
        }

        
        
        

        
        
        self.ssv_finish(raw_xe, tjb_b)
    }

    /// The `tjb_b` and the early-out guard both kernels share.
    /// Returns `None` where the SSV shortcut is not applicable and the caller must
    /// fall back to the full MSV filter.
    #[inline]
    fn ssv_setup(&self, l: usize) -> Option<u8> {
        let tjb_b = unbiased_byteify(self.scale_b, (3.0_f64 / (l as f64 + 3.0)).ln() as f32);
        if tjb_b as i32 + self.tbm_b as i32 + self.tec_b as i32 + self.bias_b as i32 >= 127 {
            return None;
        }
        Some(tjb_b)
    }

    /// Everything after the DP: the overflow tests and the score conversion.
    #[inline]
    fn ssv_finish(&self, raw_xe: u8, tjb_b: u8) -> Option<f32> {
        if raw_xe as i32 >= 255 - self.bias_b as i32 {
            return None;
        }





        let xe: i32 =
            raw_xe as i32 + (self.base_b as i32 - tjb_b as i32 - self.tbm_b as i32) - 128;

        
        if xe >= 255 - self.bias_b as i32 {
            return None;
        }

        
        
        
        if xe < self.tec_b as i32 {
            return None;
        }
        let xj = xe - self.tec_b as i32;
        if xj > self.base_b as i32 {
            return None; 
        }

        
        
        
        let sc = ((xj as f32 - tjb_b as f32) - self.base_b as f32) / self.scale_b - 3.0;
        Some(sc)
    }

    
    
    
    pub fn pvalue(&self, msvsc_nats: f32, nullsc: f32) -> f64 {
        if msvsc_nats.is_infinite() {
            return 0.0; 
        }
        let seq_score = (msvsc_nats - nullsc) as f64 / LOG2;
        crate::easel::gumbel::esl_gumbel_surv(seq_score, self.mu, self.lambda)
    }
}



/// C, p7_bg.c (`p7_bg_SetLength` then `p7_bg_NullOne`):
/// ```text
///   p7_bg_SetLength(P7_BG *bg, int L)
///   {
///     bg->p1 = (float) L / (float) (L+1);
///     ...
///   }
///
///   p7_bg_NullOne(const P7_BG *bg, const ESL_DSQ *dsq, int L, float *ret_sc)
///   {
///     *ret_sc = (float) L * log(bg->p1) + log(1.-bg->p1);
///     return eslOK;
///   }
/// ```
///
/// `bg->p1` is a `float` field written by an f32 division, so both `log(bg->p1)` and
/// `1. - bg->p1` see the *f32-rounded* ratio; only the logs and the final sum are done
/// in double, and the sum is rounded once on store. `nullsc` is subtracted from every
/// reported score, so computing `p1` in f64 shifts all of them.
pub fn null_one(l: usize) -> f32 {
    let p1 = ((l as f32) / ((l + 1) as f32)) as f64;
    ((l as f32) as f64 * p1.ln() + (1.0 - p1).ln()) as f32
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seqio::read_fasta;

    fn globins() -> P7Hmm {
        let p = format!("{}/testdata/globins4.hmm", env!("CARGO_MANIFEST_DIR"));
        P7Hmm::read_all(&p).unwrap().pop().unwrap()
    }

    #[test]
    fn quantization_constants() {
        let mp = MsvProfile::build(&globins());
        assert_eq!(mp.base_b, 190);
        assert!((mp.scale_b - 4.328085).abs() < 1e-4, "scale_b={}", mp.scale_b);
        assert_eq!(mp.m, 149);
        
        
        
        for k in 1..=mp.m {
            for x in 0..K {
                let _ = mp.rbv[k][x]; 
            }
        }
    }

    #[test]
    fn tbm_matches_formula() {
        let mp = MsvProfile::build(&globins());
        let m = 149.0f64;
        let expect = unbiased_byteify(mp.scale_b, (2.0 / (m * (m + 1.0))).ln() as f32);
        assert_eq!(mp.tbm_b, expect);
        assert_eq!(mp.tec_b, unbiased_byteify(mp.scale_b, 0.5f32.ln()));
    }

    #[test]
    fn msv_passes_all_golden_hits() {
        
        
        let mp = MsvProfile::build(&globins());
        let seqs = read_fasta(&format!("{}/testdata/globins45.fa", env!("CARGO_MANIFEST_DIR")))
            .unwrap();
        let f1 = 0.02_f64;
        let mut n_pass = 0;
        let mut top_p = 1.0_f64;
        for s in &seqs {
            let l = s.len();
            let sc = mp.msv_score(&s.dsq, l);
            let p = mp.pvalue(sc, null_one(l));
            if p <= f1 {
                n_pass += 1;
            }
            top_p = top_p.min(p);
            
            if s.name == "MYG_ESCGI" {
                assert!(p <= f1, "MYG_ESCGI failed MSV F1 (P={p:.3e}, score={sc:.2})");
                assert!(p < 1e-10, "MYG_ESCGI MSV P unexpectedly weak: {p:.3e}");
            }
        }
        
        
        assert_eq!(n_pass, 45, "MSV F1 pass count {n_pass} != C's 45");
        assert!(top_p < 1e-10, "best MSV P too weak: {top_p:.3e}");
    }

    /// The SSE SSV kernel must be byte-identical to the scalar oracle across every
    /// M mod 16 boundary, a range of L, and many random profiles/sequences.
    #[test]
    fn sse_matches_scalar_ssv() {
        #[cfg(target_arch = "x86_64")]
        {
            if !std::is_x86_feature_detected!("sse2") {
                return;
            }
            let mut rng = crate::vitfilter::tests::Lcg(0x1234_5678_9abc_def0);
            let ms = [1usize, 2, 3, 15, 16, 17, 18, 31, 32, 33, 47, 48, 49, 64, 100, 149, 200, 301];
            let ls = [1usize, 2, 5, 17, 50, 200, 500];
            let mut n = 0u64;
            let mut n_some = 0u64;
            for &m in &ms {
                for _rep in 0..8 {
                    let hmm = crate::vitfilter::tests::random_hmm(&mut rng, m);
                    let mp = MsvProfile::build(&hmm);
                    for &l in &ls {
                        let dsq = crate::vitfilter::tests::random_dsq(&mut rng, l);
                        let a = mp.ssv_score_scalar(&dsq, l);
                        // Both SIMD kernels: the sequential-layout one and C's
                        // striped/banded one, each against the scalar oracle.
                        for (tag, b) in [
                            ("seq", unsafe { ssv_score_sse(&mp, &dsq, l) }),
                            ("striped", unsafe { ssv_score_striped(&mp, &dsq, l) }),
                        ] {
                            match (a, b) {
                                (None, None) => {}
                                (Some(x), Some(y)) => {
                                    assert_eq!(
                                        x.to_bits(),
                                        y.to_bits(),
                                        "SSV {tag} mismatch M={m} L={l}: scalar={x} simd={y}"
                                    );
                                    n_some += 1;
                                }
                                _ => panic!(
                                    "SSV {tag} applicability differs M={m} L={l}: {a:?} vs {b:?}"
                                ),
                            }
                            n += 1;
                        }
                    }
                }
            }
            assert!(n > 900, "expected many checks, got {n}");
            assert!(n_some > 100, "expected many scoring cases, got {n_some}");
        }
    }
}

/// Horizontal maximum of 16 unsigned bytes.
///
/// The same shuffle/`pmaxub` cascade `p7_MSVFilter()` uses to collapse `xEv`
/// (impl_sse/msvfilter.c:159-166).
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "sse2")]
// C, impl_sse/msvfilter.c:159-166:
//   tempv = _mm_shuffle_epi32(xEv, _MM_SHUFFLE(2, 3, 0, 1));
//   xEv = _mm_max_epu8(xEv, tempv);
//   tempv = _mm_shuffle_epi32(xEv, _MM_SHUFFLE(0, 1, 2, 3));
//   xEv = _mm_max_epu8(xEv, tempv);
//   tempv = _mm_shufflelo_epi16(xEv, _MM_SHUFFLE(2, 3, 0, 1));
//   xEv = _mm_max_epu8(xEv, tempv);
//   tempv = _mm_srli_si128(xEv, 1);
//   xEv = _mm_max_epu8(xEv, tempv);
unsafe fn hmax_epu8(a: core::arch::x86_64::__m128i) -> u8 {
    use core::arch::x86_64::*;
    // The constants are C's: _MM_SHUFFLE(2,3,0,1) = 0b10_11_00_01 and
    // _MM_SHUFFLE(0,1,2,3) = 0b00_01_10_11. After the two 32-bit shuffles every
    // word holds the byte-wise max of all four, so step three has to combine the
    // two 16-bit lanes *within* a word — it must swap lanes 0<->1, not 0<->2.
    let a = _mm_max_epu8(a, _mm_shuffle_epi32::<0b10_11_00_01>(a));
    let a = _mm_max_epu8(a, _mm_shuffle_epi32::<0b00_01_10_11>(a));
    let a = _mm_max_epu8(a, _mm_shufflelo_epi16::<0b10_11_00_01>(a));
    let a = _mm_max_epu8(a, _mm_srli_si128::<1>(a));
    _mm_cvtsi128_si32(a) as u8
}

/// SSE2 SSV kernel.
///
/// **Optimization divergence (DEVELOPMENT_PRINCIPLES §1-②, §3):**
///
/// (a) *What C does.* `p7_SSVFilter()` (impl_sse/ssvfilter.c) runs the same
///     recurrence on the **striped** `om->sbv` layout, with the whole DP row held in
///     up to 16 named `__m128i` registers and the Q loop fully macro-unrolled, so the
///     row never touches memory.
///
/// (b) *What this does instead, and why.* rustyhmmer keeps the **sequential**
///     `sbv_t[x*(M+1) + k]` layout, where the row must live in memory because
///     `prev[k-1]` is a one-byte-shifted read of the previous row. The previous
///     version wrote the row and then re-read it in a second pass, because fusing the
///     running-max into the store defeated LLVM's auto-vectorizer (loop-carried
///     scalar dependency). Written with explicit intrinsics the two phases fuse into
///     one pass — `psubsb` then `pmaxub` on the value still in a register — which
///     removes one full load of the DP row per residue, and the max accumulates in a
///     vector register instead of a scalar, collapsed once at the end.
///
/// (c) *Why it is still exact.* Every operation is integer: `_mm_subs_epi8` is the
///     same saturating subtract as `i8::saturating_sub`, and `_mm_max_epu8` the same
///     unsigned max as the scalar comparison. `raw_xe` is read only after the whole
///     sequence, and max is associative and commutative, so accumulating per lane and
///     collapsing at the end gives the identical byte. Verified by the
///     `sse_matches_scalar_ssv` differential test against `ssv_score_scalar` over
///     many random profiles/sequences and every M mod 16 boundary, and end to end by
///     `--tblout`/`--domtblout` md5 on 397,961 proteins x 117 profiles.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "sse2")]
// C, impl_sse/ssvfilter.c:431-433 — the single step this kernel repeats:
//   #define STEP_SINGLE(sv)                         \
//     sv   = _mm_subs_epi8(sv, *rsc); rsc++;        \
//     xEv  = _mm_max_epu8(xEv, sv);
// C keeps `sv00..sv15` (the whole DP row) live in registers across the sequence and
// macro-unrolls the Q loop; see the divergence note above for why this does not.
#[allow(dead_code)] // kept as the second oracle for the differential test
unsafe fn ssv_score_sse(f: &MsvProfile, dsq: &[u8], l: usize) -> Option<f32> { unsafe {
    use core::arch::x86_64::*;
    let tjb_b = f.ssv_setup(l)?;
    let m = f.m;

    let mut prev = vec![-128i8; m + 1];
    let mut cur = vec![-128i8; m + 1];
    // Four independent accumulators: a single `xev` would serialize the whole row
    // on one pmaxub dependency chain, which is what the auto-vectorized two-pass
    // form avoided by letting LLVM unroll the reduction with its own accumulators.
    let neutral = _mm_set1_epi8(-128); // 0x80 == 128u8, the scalar's initial raw_xe
    let mut x0 = neutral;
    let mut x1 = neutral;
    let mut x2 = neutral;
    let mut x3 = neutral;
    // The <16 cells that do not fill a vector are reduced scalar-side and folded in
    // once at the end; max is associative, so the split is exact.
    let mut tail_max: u8 = 128;
    let nvec = m / 16;

    for i in 1..=l {
        let base = (dsq[i] as usize) * (m + 1);
        let srow = f.sbv_t.as_ptr().add(base + 1);
        let pp = prev.as_ptr();
        let cp = cur.as_mut_ptr();

        let mut v = 0usize;
        while v + 4 <= nvec {
            let k = v * 16;
            let c0 = _mm_subs_epi8(
                _mm_loadu_si128(pp.add(k) as *const __m128i),
                _mm_loadu_si128(srow.add(k) as *const __m128i),
            );
            let c1 = _mm_subs_epi8(
                _mm_loadu_si128(pp.add(k + 16) as *const __m128i),
                _mm_loadu_si128(srow.add(k + 16) as *const __m128i),
            );
            let c2 = _mm_subs_epi8(
                _mm_loadu_si128(pp.add(k + 32) as *const __m128i),
                _mm_loadu_si128(srow.add(k + 32) as *const __m128i),
            );
            let c3 = _mm_subs_epi8(
                _mm_loadu_si128(pp.add(k + 48) as *const __m128i),
                _mm_loadu_si128(srow.add(k + 48) as *const __m128i),
            );
            _mm_storeu_si128(cp.add(k + 1) as *mut __m128i, c0);
            _mm_storeu_si128(cp.add(k + 17) as *mut __m128i, c1);
            _mm_storeu_si128(cp.add(k + 33) as *mut __m128i, c2);
            _mm_storeu_si128(cp.add(k + 49) as *mut __m128i, c3);
            x0 = _mm_max_epu8(x0, c0);
            x1 = _mm_max_epu8(x1, c1);
            x2 = _mm_max_epu8(x2, c2);
            x3 = _mm_max_epu8(x3, c3);
            v += 4;
        }
        while v < nvec {
            let k = v * 16;
            let c = _mm_subs_epi8(
                _mm_loadu_si128(pp.add(k) as *const __m128i),
                _mm_loadu_si128(srow.add(k) as *const __m128i),
            );
            _mm_storeu_si128(cp.add(k + 1) as *mut __m128i, c);
            x0 = _mm_max_epu8(x0, c);
            v += 1;
        }
        let mut k = nvec * 16;
        while k < m {
            let c = (*pp.add(k)).saturating_sub(*srow.add(k));
            *cp.add(k + 1) = c;
            let u = c as u8;
            if u > tail_max {
                tail_max = u;
            }
            k += 1;
        }
        std::mem::swap(&mut prev, &mut cur);
    }

    let xev = _mm_max_epu8(_mm_max_epu8(x0, x1), _mm_max_epu8(x2, x3));
    let raw_xe = hmax_epu8(xev).max(tail_max);
    f.ssv_finish(raw_xe, tjb_b)
}}

/// One band of `W` striped vectors, held in registers across the whole sequence.
///
/// Faithful transcription of `calc_band_<W>()`, which C generates from the `CALC`
/// macro (impl_sse/ssvfilter.c:668-713) with `STEP_SINGLE`/`STEP_BANDS_W` (431-505)
/// and `CONVERT_STEP`/`CONVERT_W` (513-590):
///
/// ```text
///   #define STEP_SINGLE(sv)                         \
///     sv   = _mm_subs_epi8(sv, *rsc); rsc++;        \
///     xEv  = _mm_max_epu8(xEv, sv);
///
///   #define CONVERT_STEP(step, length_check, label, sv, pos)  \
///     length_check(label)                                     \
///     rsc = om->sbv[dsq[i]] + pos;                            \
///     step()                                                  \
///     sv = _mm_slli_si128(sv, 1);                             \
///     sv = _mm_or_si128(sv, beginv);                          \
///     i++;
///
///   #define CALC(reset, step, convert, width)       \
///     int i; int i2; int Q = p7O_NQB(om->M); __m128i *rsc; int w = width;
///     dsq++;
///     reset()
///     for (i = 0; i < L && i < Q - q - w; i++)
///       { rsc = om->sbv[dsq[i]] + i + q;  step() }
///     i = Q - q - w;
///     convert(step, LENGTH_CHECK, done1)
///   done1:
///     for (i2 = Q - q; i2 < L - Q; i2 += Q)
///       {
///         for (i = 0; i < Q - w; i++) { rsc = om->sbv[dsq[i2 + i]] + i; step() }
///         i += i2;
///         convert(step, NO_CHECK, )
///       }
///     for (i = 0; i2 + i < L && i < Q - w; i++)
///       { rsc = om->sbv[dsq[i2 + i]] + i; step() }
///     i += i2;
///     convert(step, LENGTH_CHECK, done2)
///   done2:
///     return xEv;
/// ```
///
/// `CONVERT_W` expands to `CONVERT_STEP` for `sv[W-1]` at `pos = Q-W`, then
/// `sv[W-2]` at `Q-W+1`, ... down to `sv[0]` at `Q-1`.
///
/// The `_mm_or_si128(sv, beginv)` with `beginv = set1_epi8(-128)` looks like it would
/// corrupt all sixteen lanes, but every live value here has bit 7 set: `sv` starts at
/// -128 and the applicability guard (`tjb_b+tbm_b+tec_b+bias_b < 127`) keeps it in
/// [-128,-1]. So the OR only re-inserts -128 in lane 0 after the byte shift.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "sse2")]
#[inline]
unsafe fn calc_band<const W: usize>(
    f: &MsvProfile,
    dsq: &[u8],
    l: usize,
    q: usize,
    beginv: core::arch::x86_64::__m128i,
    mut xev: core::arch::x86_64::__m128i,
) -> core::arch::x86_64::__m128i { unsafe {
    use core::arch::x86_64::*;
    let qn = f.nqb;
    let vecs = qn + P7O_EXTRA_SB;
    // C does `dsq++` so its dsq[i] is residue i+1; rustyhmmer's dsq[0] is a sentinel.
    let row = |i: usize| -> *const __m128i {
        f.sbv.as_ptr().add((dsq[i + 1] as usize) * vecs * 16) as *const __m128i
    };

    let mut sv = [beginv; W];

    // STEP_BANDS_W: all W registers, in order sv[0]..sv[W-1].
    macro_rules! step {
        ($rsc:expr) => {{
            let rsc: *const __m128i = $rsc;
            for j in 0..W {
                sv[j] = _mm_subs_epi8(sv[j], _mm_loadu_si128(rsc.add(j)));
                xev = _mm_max_epu8(xev, sv[j]);
            }
        }};
    }

    // CONVERT_W. `check` selects LENGTH_CHECK vs NO_CHECK; on a length check firing,
    // C jumps to the corresponding `done` label, i.e. returns xEv as it stands.
    macro_rules! convert {
        ($i:expr, $check:expr) => {{
            for idx in 0..W {
                let j = W - 1 - idx;
                let pos = qn - W + idx;
                if $check && $i >= l {
                    return xev;
                }
                step!(row($i).add(pos));
                sv[j] = _mm_or_si128(_mm_slli_si128::<1>(sv[j]), beginv);
                $i += 1;
            }
        }};
    }

    let mut i = 0usize;
    while i < l && i < qn - q - W {
        step!(row(i).add(i + q));
        i += 1;
    }

    i = qn - q - W;
    convert!(i, true);

    let mut i2 = qn - q;
    while i2 + qn < l {
        let mut i = 0usize;
        while i < qn - W {
            step!(row(i2 + i).add(i));
            i += 1;
        }
        let mut ii = i + i2;
        convert!(ii, false);
        let _ = ii;
        i2 += qn;
    }

    let mut i = 0usize;
    while i2 + i < l && i < qn - W {
        step!(row(i2 + i).add(i));
        i += 1;
    }
    let mut ii = i + i2;
    convert!(ii, true);

    xev
}}

/// `p7_SSVFilter()`'s DP, in C's striped/banded form — impl_sse/ssvfilter.c:875-921.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "sse2")]
unsafe fn ssv_score_striped(f: &MsvProfile, dsq: &[u8], l: usize) -> Option<f32> { unsafe {
    let tjb_b = f.ssv_setup(l)?;
    f.ssv_finish(get_xe(f, dsq, l), tjb_b)
}}

/// `get_xE()` — impl_sse/ssvfilter.c:829-869.
///
/// ```text
///   beginv =  _mm_set1_epi8(-128);
///   xEv    =  beginv;
///   /* Use the highest number of bands but no more than MAX_BANDS */
///   bands = (Q + MAX_BANDS - 1) / MAX_BANDS;
///   for (i = 0; i < bands; i++) {
///     q = (Q * (i + 1)) / bands;
///     xEv = fs[q-last_q](dsq, L, om, last_q, beginv, xEv);
///     last_q = q;
///   }
///   return esl_sse_hmax_epu8(xEv);
/// ```
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "sse2")]
unsafe fn get_xe(f: &MsvProfile, dsq: &[u8], l: usize) -> u8 { unsafe {
    use core::arch::x86_64::*;
    let qn = f.nqb;
    let beginv = _mm_set1_epi8(-128);
    let mut xev = beginv;

    let bands = qn.div_ceil(MAX_BANDS);
    let mut last_q = 0usize;
    for i in 0..bands {
        let q = (qn * (i + 1)) / bands;
        let w = q - last_q;
        // C indexes a table of calc_band_1..calc_band_MAX_BANDS; const generics need
        // the width as a compile-time value, so dispatch on it.
        macro_rules! dispatch {
            ($($n:literal),*) => {
                match w {
                    $($n => xev = calc_band::<$n>(f, dsq, l, last_q, beginv, xev),)*
                    _ => unreachable!("band width {w} out of range 1..={MAX_BANDS}"),
                }
            };
        }
        dispatch!(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14);
        last_q = q;
    }
    hmax_epu8(xev)
}}
