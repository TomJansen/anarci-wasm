























use crate::alphabet::{AMINO_FREQ, K, KP};
use crate::hmmfile::{P7Hmm, TDD, TDM, TIM, TII, TMD, TMI, TMM};

const LOG2: f64 = std::f64::consts::LN_2;


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


const NEGINF: i16 = -32768;








#[inline]
fn wordify(scale_w: f32, sc: f32) -> i16 {
    if sc == f32::NEG_INFINITY {
        return NEGINF;
    }
    let v = (scale_w * sc).round();
    if v >= 32767.0 {
        32767
    } else if v <= -32768.0 {
        NEGINF
    } else {
        v as i16
    }
}






#[inline]
fn adds(a: i16, b: i16) -> i16 {
    a.saturating_add(b)
}


#[derive(Debug, Clone)]
pub struct VitFilter {
    pub m: usize,
    
    pub scale_w: f32,
    
    pub base_w: i16,

    #[allow(dead_code)]
    rwv: Vec<[i16; KP]>,
    
    
    
    
    rwv_t: Vec<i16>,
    
    
    
    tmm: Vec<i16>,
    tim: Vec<i16>,
    tdm: Vec<i16>,
    
    
    vbm: Vec<i16>,
    
    tmi: Vec<i16>,
    tii: Vec<i16>,
    
    tmd: Vec<i16>,
    tdd: Vec<i16>,
    
    xw_e_move: i16,
    xw_e_loop: i16,

    // ---- Striped SSE (Farrar) layout, built in `build` from the SAME quantized
    // i16 tables the scalar path uses, so the striped evaluation is byte-identical
    // (integer saturating max/add is associative). Faithful port of `vf_conversion`
    // (impl_sse/p7_oprofile.c) + `p7_ViterbiFilter` (impl_sse/vitfilter.c). ----
    /// Q = p7O_NQW(M) = ESL_MAX(2, (M-1)/8 + 1) — # of 8-wide i16 vectors.
    #[cfg_attr(not(target_arch = "x86_64"), allow(dead_code))]
    q: usize,
    /// Striped transitions `om->twv`: 7 interleaved vectors per q
    /// [BM,MM,IM,DM,MD,MI,II] at vector index `q*7+t`, then a DD block of Q vectors
    /// at `7*Q+q`. Flattened: vector `v` occupies `tw[v*8..v*8+8]`.
    #[cfg_attr(not(target_arch = "x86_64"), allow(dead_code))]
    tw: Vec<i16>,
    /// Striped match emissions `om->rwv`: residue x, vector q at `rw[(x*Q+q)*8..+8]`;
    /// slot z holds model position k=q+1+z*Q (or -inf if k>M).
    #[cfg_attr(not(target_arch = "x86_64"), allow(dead_code))]
    rw: Vec<i16>,
    /// Lazy-F DD short-circuit bound `om->ddbound_w`
    /// = max_{k=2..M-2}(TDD(k)+TDM(k+1)-TBM(k+2)). (vf_conversion:910-918)
    #[cfg_attr(not(target_arch = "x86_64"), allow(dead_code))]
    ddbound_w: i16,

    pub mu: f64,
    pub lambda: f64,
}

impl VitFilter {
    
    pub fn build(hmm: &P7Hmm) -> Self {
        let m = hmm.m;
        
        let scale_w = (500.0_f64 / LOG2) as f32; 
        let base_w: i16 = 12000;

        
        
        
        let mut rwv = vec![[NEGINF; KP]; m + 1];
        for k in 1..=m {
            let mut sc = [f32::NEG_INFINITY; KP];
            for x in 0..K {
                sc[x] = ((hmm.mat[k][x] as f64) / (AMINO_FREQ[x] as f32 as f64)).ln() as f32;
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
            for x in 0..KP {
                rwv[k][x] = wordify(scale_w, sc[x]); 
            }
        }
        
        let mut rwv_t = vec![NEGINF; KP * (m + 1)];
        for k in 1..=m {
            for x in 0..KP {
                rwv_t[x * (m + 1) + k] = rwv[k][x];
            }
        }

        
        
        
        
        
        let mut tmm = vec![NEGINF; m + 1];
        let mut tim = vec![NEGINF; m + 1];
        let mut tdm = vec![NEGINF; m + 1];
        let mut tmi = vec![NEGINF; m + 1];
        let mut tii = vec![NEGINF; m + 1];
        let mut tmd = vec![NEGINF; m + 1];
        let mut tdd = vec![NEGINF; m + 1];
        for j in 1..m {
            tmm[j] = wordify(scale_w, (hmm.t[j][TMM] as f64).ln() as f32).min(0);
            tim[j] = wordify(scale_w, (hmm.t[j][TIM] as f64).ln() as f32).min(0);
            tdm[j] = wordify(scale_w, (hmm.t[j][TDM] as f64).ln() as f32).min(0);
            tmi[j] = wordify(scale_w, (hmm.t[j][TMI] as f64).ln() as f32).min(0);
            tii[j] = wordify(scale_w, (hmm.t[j][TII] as f64).ln() as f32).min(-1);
            tmd[j] = wordify(scale_w, (hmm.t[j][TMD] as f64).ln() as f32).min(0);
            tdd[j] = wordify(scale_w, (hmm.t[j][TDD] as f64).ln() as f32).min(0);
        }

        
        
        let mut occ = vec![0.0f32; m + 1];
        if m >= 1 {
            occ[1] = hmm.t[0][TMI] + hmm.t[0][TMM];
        }
        for k in 2..=m {
            // C, p7_hmm.c:`p7_hmm_CalculateOccupancy` -- `(1.0-mocc[k-1])` uses a
            // double literal, so that term and the final addition are evaluated in
            // double and rounded to f32 once, on the store.
            occ[k] = ((occ[k - 1] * (hmm.t[k - 1][TMM] + hmm.t[k - 1][TMI])) as f64
                + (1.0 - occ[k - 1] as f64) * hmm.t[k - 1][TDM] as f64) as f32;
        }
        let mut z = 0.0f32;
        for k in 1..=m {
            z += occ[k] * (m - k + 1) as f32;
        }
        
        let mut vbm = vec![NEGINF; m]; 
        for k in 1..=m {
            let logv = ((occ[k] / z) as f64).ln() as f32;
            vbm[k - 1] = wordify(scale_w, logv).min(0);
        }

        
        
        let xw_e_move = wordify(scale_w, -(LOG2 as f32));
        let xw_e_loop = wordify(scale_w, -(LOG2 as f32));

        // ---- Striped SSE layout (vf_conversion), built from the SAME i16 tables
        // computed above. NB: rustyhmmer's `vbm` is 0-based — vbm[k-1] = B->M_k
        // (the scalar path reads vbm[j] for M at k=j+1). tmm/tim/tdm/tmi/tii/tmd/tdd
        // are 1-based [k] like C. ----
        // Q = ESL_MAX(2, (M-1)/8 + 1) (impl_sse.h:25). M>=1 always here.
        let q = std::cmp::max(2, (m - 1) / 8 + 1);

        // Striped match emissions rwv[x][q][z], k = q+1 + z*Q (vf_conversion:852-858),
        // k<=M else -inf.
        let mut rw = vec![NEGINF; KP * q * 8];
        for x in 0..KP {
            for qi in 0..q {
                for z in 0..8 {
                    let k = (qi + 1) + z * q; // model position
                    rw[(x * q + qi) * 8 + z] = if k <= m { rwv[k][x] } else { NEGINF };
                }
            }
        }

        // Striped transitions twv (vf_conversion:860-888): 7 interleaved vectors per q
        // [BM,MM,IM,DM,MD,MI,II] then a DD block of Q vectors. Into-M (BM/MM/IM/DM)
        // use k<=M reading vbm[k-1]/t[k-1]; the rest (MD/MI/II/DD) use k<M reading t[k].
        let mut tw = vec![NEGINF; (7 * q + q) * 8];
        for qi in 0..q {
            for z in 0..8 {
                let k = (qi + 1) + z * q; // model position
                let vbase = qi * 7;
                let (bm, mm, im, dm) = if k <= m {
                    (vbm[k - 1], tmm[k - 1], tim[k - 1], tdm[k - 1])
                } else {
                    (NEGINF, NEGINF, NEGINF, NEGINF)
                };
                tw[(vbase) * 8 + z] = bm;
                tw[(vbase + 1) * 8 + z] = mm;
                tw[(vbase + 2) * 8 + z] = im;
                tw[(vbase + 3) * 8 + z] = dm;
                let (md, mi, ii) = if k < m {
                    (tmd[k], tmi[k], tii[k])
                } else {
                    (NEGINF, NEGINF, NEGINF)
                };
                tw[(vbase + 4) * 8 + z] = md;
                tw[(vbase + 5) * 8 + z] = mi;
                tw[(vbase + 6) * 8 + z] = ii;
                tw[(7 * q + qi) * 8 + z] = if k < m { tdd[k] } else { NEGINF };
            }
        }

        // Lazy-F DD bound (vf_conversion:910-918): ddtmp = TDD(k)+TDM(k+1)-TBM(k+2),
        // k=2..M-2; TBM(k+2) = vbm[(k+2)-1] = vbm[k+1] in the 0-based array.
        let mut ddbound_w: i16 = NEGINF;
        let mut k = 2usize;
        while k < m.saturating_sub(1) {
            let ddtmp = tdd[k] as i32 + tdm[k + 1] as i32 - vbm[k + 1] as i32;
            ddbound_w = std::cmp::max(ddbound_w as i32, ddtmp) as i16;
            k += 1;
        }

        VitFilter {
            m,
            scale_w,
            base_w,
            rwv,
            rwv_t,
            tmm,
            tim,
            tdm,
            vbm,
            tmi,
            tii,
            tmd,
            tdd,
            xw_e_move,
            xw_e_loop,
            q,
            tw,
            rw,
            ddbound_w,
            mu: hmm.evparam.vit_mu,
            lambda: hmm.evparam.vit_lambda,
        }
    }

    
    
    
    
    
    /// Whole-sequence Viterbi filter score in bits. Runtime dispatcher: striped-SSE
    /// kernel [`vit_score_sse`] on x86_64 (byte-identical output, verified by the
    /// `sse_matches_scalar` differential test), else the scalar oracle.
    #[inline]
    pub fn vit_score(&self, dsq: &[u8], l: usize) -> f32 {
        #[cfg(target_arch = "x86_64")]
        {
            if std::is_x86_feature_detected!("sse2") {
                // SAFETY: guarded by the sse2 feature check just above.
                return unsafe { vit_score_sse(self, dsq, l) };
            }
        }
        self.vit_score_scalar(dsq, l)
    }

    /// Natural-order scalar oracle. Its full serial DD sweep is bit-identical to the
    /// striped SSE lazy-F convergence (integer saturating max/add is associative).
    pub fn vit_score_scalar(&self, dsq: &[u8], l: usize) -> f32 {
        let m = self.m;

        
        
        let pmove = 3.0_f32 / (l as f32 + 3.0);
        let xw_move = wordify(self.scale_w, pmove.ln()); 
        

        
        let mut mp = vec![NEGINF; m + 1]; 
        let mut ip = vec![NEGINF; m + 1]; 
        let mut dp = vec![NEGINF; m + 1]; 
        let mut mc = vec![NEGINF; m + 1]; 
        let mut ic = vec![NEGINF; m + 1]; 
        let mut dc = vec![NEGINF; m + 1]; 

        let mut xn: i16 = self.base_w; 
        let mut xb: i16 = adds(xn, xw_move); 
        let mut xj: i16 = NEGINF;
        let mut xc: i16 = NEGINF;

        for i in 1..=l {
            let x = dsq[i] as usize;
            mc[0] = NEGINF;
            ic[0] = NEGINF;
            dc[0] = NEGINF;

            
            
            
            
            
            
            
            

            
            
            {
                let base = x * (m + 1);
                let rwv_row = &self.rwv_t[base + 1..base + m + 1]; 
                let mpm = &mp[0..m]; 
                let ipm = &ip[0..m]; 
                let dpm = &dp[0..m]; 
                let vbm = &self.vbm[0..m];
                let tmm = &self.tmm[0..m];
                let tim = &self.tim[0..m];
                let tdm = &self.tdm[0..m];
                let mc_o = &mut mc[1..m + 1];
                for j in 0..m {
                    let sc = adds(xb, vbm[j])
                        .max(adds(mpm[j], tmm[j]))
                        .max(adds(ipm[j], tim[j]))
                        .max(adds(dpm[j], tdm[j]));
                    mc_o[j] = adds(sc, rwv_row[j]);
                }
            }
            
            {
                let mpk = &mp[1..m + 1];
                let ipk = &ip[1..m + 1];
                let tmi = &self.tmi[1..m + 1];
                let tii = &self.tii[1..m + 1];
                let ic_o = &mut ic[1..m + 1];
                for j in 0..m {
                    ic_o[j] = adds(mpk[j], tmi[j]).max(adds(ipk[j], tii[j]));
                }
            }
            
            let mut xe: i16 = NEGINF;
            for &v in &mc[1..m + 1] {
                if v > xe {
                    xe = v;
                }
            }
            
            
            for k in 1..=m {
                dc[k] = adds(mc[k - 1], self.tmd[k - 1]).max(adds(dc[k - 1], self.tdd[k - 1]));
            }

            
            
            if xe >= 32767 {
                return f32::INFINITY;
            }
            
            xn = adds(xn, 0);
            
            xc = xc.max(adds(xe, self.xw_e_move));
            
            xj = xj.max(adds(xe, self.xw_e_loop));
            
            xb = adds(xj, xw_move).max(adds(xn, xw_move));

            std::mem::swap(&mut mp, &mut mc);
            std::mem::swap(&mut ip, &mut ic);
            std::mem::swap(&mut dp, &mut dc);
        }

        
        if xc > NEGINF {
            let mut ret = xc as f32 + xw_move as f32 - self.base_w as f32; 
            ret /= self.scale_w;
            ret - 3.0 
        } else {
            f32::NEG_INFINITY
        }
    }

    
    
    
    pub fn pvalue(&self, vfsc: f32, filtersc: f32, _l: usize) -> f64 {
        if vfsc.is_infinite() {
            return if vfsc > 0.0 { 0.0 } else { 1.0 }; 
        }
        let seq_score = (vfsc - filtersc) as f64 / LOG2;
        crate::easel::gumbel::esl_gumbel_surv(seq_score, self.mu, self.lambda)
    }
}

/// Horizontal max of 8 i16 lanes — faithful `esl_sse_hmax_epi16` (esl_sse.h:75).
#[cfg(target_arch = "x86_64")]
#[inline]
#[target_feature(enable = "sse2")]
unsafe fn hmax_epi16(a: core::arch::x86_64::__m128i) -> i16 {
    use core::arch::x86_64::*;
    let a = _mm_max_epi16(a, _mm_shuffle_epi32::<0b01_00_11_10>(a));
    let a = _mm_max_epi16(a, _mm_shufflelo_epi16::<0b01_00_11_10>(a));
    let a = _mm_max_epi16(a, _mm_srli_epi32::<16>(a));
    _mm_cvtsi128_si32(a) as i16
}

/// Striped SIMD (Farrar) i16 Viterbi filter — faithful port of `p7_ViterbiFilter`
/// (impl_sse/vitfilter.c:82). Byte-identical to [`VitFilter::vit_score_scalar`]:
/// integer saturating max/add is associative, so the striped/interleaved order
/// yields the same xE/xC and thus the same final bit-score; the lazy-F DD loop
/// converges to the same D values the scalar full serial DD sweep computes.
///
/// SAFETY: caller must ensure the `sse2` target feature is available (guaranteed on
/// x86_64 baseline; checked in [`VitFilter::vit_score`]).
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "sse2")]
unsafe fn vit_score_sse(f: &VitFilter, dsq: &[u8], l: usize) -> f32 { unsafe {
    use core::arch::x86_64::*;
    let q = f.q;

    // p7_oprofile_ReconfigLength: xw[MOVE] = wordify(log(pmove)), pmove = 3/(L+3);
    // N/C/J LOOP costs stay 0.
    let pmove = 3.0_f32 / (l as f32 + 3.0);
    let xw_move = wordify(f.scale_w, pmove.ln());

    let neg_inf = _mm_set1_epi16(NEGINF);
    let neginfv = _mm_srli_si128::<14>(neg_inf);

    #[inline(always)]
    unsafe fn ld(base: *const i16, v: usize) -> core::arch::x86_64::__m128i { unsafe {
        core::arch::x86_64::_mm_loadu_si128(base.add(v * 8) as *const core::arch::x86_64::__m128i)
    }}
    let twp = f.tw.as_ptr();
    let rwp = f.rw.as_ptr();

    let mut mmx = vec![neg_inf; q];
    let mut imx = vec![neg_inf; q];
    let mut dmx = vec![neg_inf; q];

    let mut xn: i16 = f.base_w;
    let mut xb: i16 = adds(xn, xw_move);
    let mut xj: i16 = NEGINF;
    let mut xc: i16 = NEGINF;

    let dd_base = 7 * q;

    for i in 1..=l {
        let x = dsq[i] as usize;
        let rbase = x * q;
        let mut dcv = neg_inf;
        let mut xev = neg_inf;
        let mut dmaxv = neg_inf;
        let xbv = _mm_set1_epi16(xb);

        let mut mpv = _mm_or_si128(_mm_slli_si128::<2>(mmx[q - 1]), neginfv);
        let mut dpv = _mm_or_si128(_mm_slli_si128::<2>(dmx[q - 1]), neginfv);
        let mut ipv = _mm_or_si128(_mm_slli_si128::<2>(imx[q - 1]), neginfv);

        for qi in 0..q {
            let vb = qi * 7;
            // M(i,q) = max(B->M, M->M, I->M, D->M) + emission. (vitfilter.c:145-150)
            let mut sv = _mm_adds_epi16(xbv, ld(twp, vb)); // BM
            sv = _mm_max_epi16(sv, _mm_adds_epi16(mpv, ld(twp, vb + 1))); // MM
            sv = _mm_max_epi16(sv, _mm_adds_epi16(ipv, ld(twp, vb + 2))); // IM
            sv = _mm_max_epi16(sv, _mm_adds_epi16(dpv, ld(twp, vb + 3))); // DM
            sv = _mm_adds_epi16(sv, ld(rwp, rbase + qi)); // emission
            xev = _mm_max_epi16(xev, sv);

            mpv = mmx[qi];
            dpv = dmx[qi];
            ipv = imx[qi];
            mmx[qi] = sv;
            dmx[qi] = dcv;

            dcv = _mm_adds_epi16(sv, ld(twp, vb + 4)); // MD
            dmaxv = _mm_max_epi16(dcv, dmaxv);

            let sv_i = _mm_adds_epi16(mpv, ld(twp, vb + 5)); // MI
            imx[qi] = _mm_max_epi16(sv_i, _mm_adds_epi16(ipv, ld(twp, vb + 6))); // II
        }

        // Specials (vitfilter.c:175-181). Identical to the scalar oracle.
        let xe = hmax_epi16(xev);
        if xe >= 32767 {
            return f32::INFINITY; // eslERANGE — scalar returns +inf here too
        }
        xn = adds(xn, 0);
        xc = xc.max(adds(xe, f.xw_e_move));
        xj = xj.max(adds(xe, f.xw_e_loop));
        xb = adds(xj, xw_move).max(adds(xn, xw_move));

        // Lazy-F DD loop (vitfilter.c:197-231).
        let dmax = hmax_epi16(dmaxv);
        if (dmax as i32) + (f.ddbound_w as i32) > (xb as i32) {
            dcv = _mm_or_si128(_mm_slli_si128::<2>(dcv), neginfv);
            for qi in 0..q {
                let d = _mm_max_epi16(dcv, dmx[qi]);
                dmx[qi] = d;
                dcv = _mm_adds_epi16(d, ld(twp, dd_base + qi));
            }
            loop {
                dcv = _mm_or_si128(_mm_slli_si128::<2>(dcv), neginfv);
                let mut qi = 0;
                while qi < q {
                    if _mm_movemask_epi8(_mm_cmpgt_epi16(dcv, dmx[qi])) == 0 {
                        break;
                    }
                    let d = _mm_max_epi16(dcv, dmx[qi]);
                    dmx[qi] = d;
                    dcv = _mm_adds_epi16(d, ld(twp, dd_base + qi));
                    qi += 1;
                }
                if qi != q {
                    break;
                }
            }
        } else {
            dmx[0] = _mm_or_si128(_mm_slli_si128::<2>(dcv), neginfv);
        }
    }

    // C->T (vitfilter.c:239-247).
    if xc > NEGINF {
        let mut ret = xc as f32 + xw_move as f32 - f.base_w as f32;
        ret /= f.scale_w;
        ret - 3.0
    } else {
        f32::NEG_INFINITY
    }
}}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    fn globins() -> P7Hmm {
        let p = format!("{}/testdata/globins4.hmm", env!("CARGO_MANIFEST_DIR"));
        P7Hmm::read_all(&p).unwrap().pop().unwrap()
    }

    #[test]
    fn quantization_constants() {
        let vf = VitFilter::build(&globins());
        assert_eq!(vf.base_w, 12000);
        assert_eq!(vf.m, 149);
        
        assert!((vf.scale_w - 721.348).abs() < 0.01, "scale_w={}", vf.scale_w);
        
        assert_eq!(vf.xw_e_move, -500);
        assert_eq!(vf.xw_e_loop, -500);
    }

    /// Tiny deterministic LCG (Numerical Recipes) — no external RNG dependency.
    pub(crate) struct Lcg(pub(crate) u64);
    impl Lcg {
        fn next_u32(&mut self) -> u32 {
            self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (self.0 >> 33) as u32
        }
        /// Uniform f32 in [0.05, 1.05) — strictly positive so every log-score is finite.
        fn unit(&mut self) -> f32 {
            0.05 + (self.next_u32() as f32 / u32::MAX as f32)
        }
    }

    /// Build a random-but-valid protein p7 HMM of length M (only the fields the VF
    /// build reads: m, mat, t, evparam). t layout [MM,MI,MD,IM,II,DM,DD].
    pub(crate) fn random_hmm(rng: &mut Lcg, m: usize) -> P7Hmm {
        let mut mat = vec![[0.0f32; K]; m + 1];
        for k in 1..=m {
            let mut s = 0.0f32;
            for x in 0..K {
                mat[k][x] = rng.unit();
                s += mat[k][x];
            }
            for x in 0..K {
                mat[k][x] /= s;
            }
        }
        let mut t = vec![[0.0f32; 7]; m + 1];
        for k in 0..=m {
            let (mm, mi, md) = (rng.unit(), rng.unit(), rng.unit());
            let ms = mm + mi + md;
            let (im, ii) = (rng.unit(), rng.unit());
            let is = im + ii;
            let (dm, dd) = (rng.unit(), rng.unit());
            let ds = dm + dd;
            // [TMM,TMI,TMD,TIM,TII,TDM,TDD]
            t[k] = [mm / ms, mi / ms, md / ms, im / is, ii / is, dm / ds, dd / ds];
        }
        let mut h = P7Hmm {
            name: "rand".into(),
            acc: None,
            desc: None,
            m,
            mat,
            ins: vec![[0.0f32; K]; m + 1],
            t,
            compo: None,
            evparam: crate::hmmfile::EvParams::default(),
            cutoffs: crate::hmmfile::Cutoffs::default(),
            format: crate::hmmfile::FORMAT_3F,
            consensus: vec![b' '; m + 1],
            rf: None,
            mm: None,
            cs: None,
            map: None,
        };
        h.set_consensus();
        h
    }

    pub(crate) fn random_dsq(rng: &mut Lcg, l: usize) -> Vec<u8> {
        // Index 0 sentinel (unused by the 1..=l loop); codes 0..25 exercise
        // canonical + gap(20) + IUPAC-degenerate (21..25) emission rows.
        let mut dsq = vec![255u8; l + 2];
        for i in 1..=l {
            dsq[i] = (rng.next_u32() % 26) as u8;
        }
        dsq
    }

    /// The striped-SSE kernel must be byte-identical to the scalar oracle across
    /// every Q=ceil(M/8) boundary, a range of L, and many random profiles/seqs.
    #[test]
    fn sse_matches_scalar() {
        #[cfg(target_arch = "x86_64")]
        {
            if !std::is_x86_feature_detected!("sse2") {
                return;
            }
            let mut rng = Lcg(0x0f0e_0d0c_0b0a_0908);
            let ms = [1usize, 2, 3, 7, 8, 9, 15, 16, 17, 23, 24, 25, 31, 32, 40, 63, 64, 65, 100, 149, 200];
            let ls = [1usize, 2, 5, 17, 50, 200, 500];
            let mut n_checked = 0u64;
            for &m in &ms {
                for _rep in 0..8 {
                    let hmm = random_hmm(&mut rng, m);
                    let vf = VitFilter::build(&hmm);
                    for &l in &ls {
                        let dsq = random_dsq(&mut rng, l);
                        let sc_scalar = vf.vit_score_scalar(&dsq, l);
                        let sc_sse = unsafe { vit_score_sse(&vf, &dsq, l) };
                        assert_eq!(
                            sc_scalar.to_bits(),
                            sc_sse.to_bits(),
                            "VF mismatch M={m} L={l}: scalar={sc_scalar} sse={sc_sse}"
                        );
                        n_checked += 1;
                    }
                }
            }
            assert!(n_checked > 1000, "expected many checks, got {n_checked}");
        }
    }
}

