

















use crate::alphabet::{AMINO_FREQ, K, KP};
use crate::forward::ForwardFilter;
use crate::hmmfile::P7Hmm;
use crate::omx::{self, Omx, XFactors};
use std::sync::OnceLock;

const NEG_INF: f32 = f32::NEG_INFINITY;




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








const P7_LOGSUM_SCALE: f32 = 1000.0;
const P7_LOGSUM_TBL: usize = 16000;

fn flogsum_table() -> &'static [f32; P7_LOGSUM_TBL] {
    
    
    
    static TABLE: OnceLock<Box<[f32; P7_LOGSUM_TBL]>> = OnceLock::new();
    TABLE.get_or_init(|| {
        let mut t = Box::new([0.0f32; P7_LOGSUM_TBL]);
        for i in 0..P7_LOGSUM_TBL {
            t[i] = (1.0 + ((-(i as f64)) / P7_LOGSUM_SCALE as f64).exp()).ln() as f32;
        }
        t
    })
}






#[inline]
pub fn p7_flogsum(a: f32, b: f32) -> f32 {
    let max = if a > b { a } else { b };
    let min = if a > b { b } else { a };
    if min == NEG_INF || (max - min) >= 15.7f32 {
        max
    } else {
        max + flogsum_table()[((max - min) * P7_LOGSUM_SCALE) as usize]
    }
}






const P7P_MM: usize = 0;
const P7P_IM: usize = 1;
const P7P_DM: usize = 2;
const P7P_BM: usize = 3;
const P7P_MD: usize = 4;
const P7P_DD: usize = 5;
const P7P_MI: usize = 6;
const P7P_II: usize = 7;
const P7P_NTRANS: usize = 8;



const P7H_MM: usize = 0;
const P7H_MI: usize = 1;
const P7H_MD: usize = 2;
const P7H_IM: usize = 3;
const P7H_II: usize = 4;
const P7H_DM: usize = 5;
const P7H_DD: usize = 6;


const P7P_E: usize = 0;
const P7P_N: usize = 1;
const P7P_J: usize = 2;
const P7P_C: usize = 3;
const P7P_LOOP: usize = 0;
const P7P_MOVE: usize = 1;



#[derive(Clone)]
pub struct Gm {
    pub m: usize,
    pub tsc: Vec<f32>,      
    pub rsc: Vec<Vec<f32>>, 
    pub xsc: [[f32; 2]; 4],
    pub nj: f32,
    pub l: i32,
}

const P7P_NR: usize = 2;
const P7P_MSC: usize = 0;
const P7P_ISC: usize = 1;




pub fn build_local_profile(hmm: &P7Hmm, l: usize) -> Gm {
    let m = hmm.m;
    let mut tsc = vec![0.0f32; (m + 1) * P7P_NTRANS];
    let mut rsc = vec![vec![0.0f32; (m + 1) * P7P_NR]; KP];

    
    
    
    
    
    
    
    let mut occ = vec![0.0f32; m + 1];
    
    if m >= 1 {
        occ[1] = hmm.t[0][P7H_MI] + hmm.t[0][P7H_MM];
    }
    for k in 2..=m {
        // C, p7_hmm.c:`p7_hmm_CalculateOccupancy` -- `(1.0-mocc[k-1])` uses a double
        // literal, so that term and the final addition are evaluated in double and
        // rounded to f32 once, on the store.
        occ[k] = ((occ[k - 1] * (hmm.t[k - 1][P7H_MM] + hmm.t[k - 1][P7H_MI])) as f64
            + (1.0 - occ[k - 1] as f64) * hmm.t[k - 1][P7H_DM] as f64) as f32;
    }
    let mut z = 0.0f32;
    for k in 1..=m {
        z += occ[k] * (m - k + 1) as f32;
    }
    for k in 1..=m {
        
        tsc[(k - 1) * P7P_NTRANS + P7P_BM] = ((occ[k] / z) as f64).ln() as f32;
    }

    
    
    
    for k in 1..m {
        let tk = &hmm.t[k];
        let base = k * P7P_NTRANS;
        tsc[base + P7P_MM] = (tk[P7H_MM] as f64).ln() as f32;
        tsc[base + P7P_MI] = (tk[P7H_MI] as f64).ln() as f32;
        tsc[base + P7P_MD] = (tk[P7H_MD] as f64).ln() as f32;
        tsc[base + P7P_IM] = (tk[P7H_IM] as f64).ln() as f32;
        tsc[base + P7P_II] = (tk[P7H_II] as f64).ln() as f32;
        tsc[base + P7P_DM] = (tk[P7H_DM] as f64).ln() as f32;
        tsc[base + P7P_DD] = (tk[P7H_DD] as f64).ln() as f32;
    }

    
    
    
    
    
    
    for k in 1..=m {
        let mut sc = [0.0f32; KP];
        for x in 0..K {
            sc[x] = ((hmm.mat[k][x] as f64) / (AMINO_FREQ[x] as f32 as f64)).ln() as f32;
        }
        sc[K] = NEG_INF; 
        sc[KP - 2] = NEG_INF; 
        sc[KP - 1] = NEG_INF; 
        
        
        
        for x in (K + 1)..=(KP - 3) {
            let set = degen_set(x);
            let mut result = 0.0f32;
            let mut denom = 0.0f32;
            for &y in set {
                result += AMINO_FREQ[y] as f32 * sc[y];
                denom += AMINO_FREQ[y] as f32;
            }
            sc[x] = result / denom;
        }
        for x in 0..KP {
            rsc[x][k * P7P_NR + P7P_MSC] = sc[x];
        }
    }

    
    
    
    for x in 0..KP {
        for k in 1..m {
            rsc[x][k * P7P_NR + P7P_ISC] = 0.0;
        }
        rsc[x][m * P7P_NR + P7P_ISC] = NEG_INF;
    }
    for k in 1..=m {
        rsc[K][k * P7P_NR + P7P_ISC] = NEG_INF;
        rsc[KP - 2][k * P7P_NR + P7P_ISC] = NEG_INF;
        rsc[KP - 1][k * P7P_NR + P7P_ISC] = NEG_INF;
    }

    
    
    let mut xsc = [[0.0f32; 2]; 4];
    let log2 = std::f64::consts::LN_2 as f32;
    xsc[P7P_E][P7P_MOVE] = -log2;
    xsc[P7P_E][P7P_LOOP] = -log2;
    let nj = 1.0f32;

    let mut gm = Gm { m, tsc, rsc, xsc, nj, l: 0 };
    reconfig_length(&mut gm, l as i32);
    gm
}




pub fn reconfig_length(gm: &mut Gm, l: i32) {
    let pmove = (2.0f32 + gm.nj) / (l as f32 + 2.0f32 + gm.nj);
    let ploop = 1.0f32 - pmove;
    let ll = ploop.ln();
    let lm = pmove.ln();
    for &s in &[P7P_N, P7P_C, P7P_J] {
        gm.xsc[s][P7P_LOOP] = ll;
        gm.xsc[s][P7P_MOVE] = lm;
    }
    gm.l = l;
}





const P7G_M: usize = 0;
const P7G_I: usize = 1;
const P7G_D: usize = 2;
const P7G_NSCELLS: usize = 3;
const P7G_E: usize = 0;
const P7G_N: usize = 1;
const P7G_J: usize = 2;
const P7G_B: usize = 3;
const P7G_C: usize = 4;
const P7G_NXCELLS: usize = 5;

pub struct P7Gmx {
    pub m: usize,
    pub l: usize,
    pub dp: Vec<f32>,  
    pub xmx: Vec<f32>, 
}

impl P7Gmx {
    pub fn new(m: usize, l: usize) -> Self {
        P7Gmx {
            m,
            l,
            dp: vec![0.0; (l + 1) * (m + 1) * P7G_NSCELLS],
            xmx: vec![0.0; (l + 1) * P7G_NXCELLS],
        }
    }
    
    
    
    #[inline]
    fn mmx(&self, i: usize, k: usize) -> f32 {
        let idx = i * (self.m + 1) * P7G_NSCELLS + k * P7G_NSCELLS + P7G_M;
        debug_assert!(idx < self.dp.len());
        unsafe { *self.dp.get_unchecked(idx) }
    }
    #[inline]
    fn imx(&self, i: usize, k: usize) -> f32 {
        let idx = i * (self.m + 1) * P7G_NSCELLS + k * P7G_NSCELLS + P7G_I;
        debug_assert!(idx < self.dp.len());
        unsafe { *self.dp.get_unchecked(idx) }
    }
    #[inline]
    fn dmx(&self, i: usize, k: usize) -> f32 {
        let idx = i * (self.m + 1) * P7G_NSCELLS + k * P7G_NSCELLS + P7G_D;
        debug_assert!(idx < self.dp.len());
        unsafe { *self.dp.get_unchecked(idx) }
    }
    #[inline]
    fn xmx(&self, i: usize, s: usize) -> f32 {
        let idx = i * P7G_NXCELLS + s;
        debug_assert!(idx < self.xmx.len());
        unsafe { *self.xmx.get_unchecked(idx) }
    }
    #[inline]
    fn set_mmx(&mut self, i: usize, k: usize, v: f32) {
        let idx = i * (self.m + 1) * P7G_NSCELLS + k * P7G_NSCELLS + P7G_M;
        debug_assert!(idx < self.dp.len());
        unsafe { *self.dp.get_unchecked_mut(idx) = v; }
    }
    #[inline]
    fn set_imx(&mut self, i: usize, k: usize, v: f32) {
        let idx = i * (self.m + 1) * P7G_NSCELLS + k * P7G_NSCELLS + P7G_I;
        debug_assert!(idx < self.dp.len());
        unsafe { *self.dp.get_unchecked_mut(idx) = v; }
    }
    #[inline]
    fn set_dmx(&mut self, i: usize, k: usize, v: f32) {
        let idx = i * (self.m + 1) * P7G_NSCELLS + k * P7G_NSCELLS + P7G_D;
        debug_assert!(idx < self.dp.len());
        unsafe { *self.dp.get_unchecked_mut(idx) = v; }
    }
    #[inline]
    fn set_xmx(&mut self, i: usize, s: usize, v: f32) {
        let idx = i * P7G_NXCELLS + s;
        debug_assert!(idx < self.xmx.len());
        unsafe { *self.xmx.get_unchecked_mut(idx) = v; }
    }

    #[allow(dead_code)]
    fn clone_matrix(&self) -> (Vec<f32>, Vec<f32>) {
        (self.dp.clone(), self.xmx.clone())
    }
}

#[inline]
fn tsc(gm: &Gm, s: usize, k: usize) -> f32 {
    gm.tsc[k * P7P_NTRANS + s]
}








pub fn p7_gforward(dsq: &[u8], l: usize, gm: &Gm, gx: &mut P7Gmx) -> f32 {
    let mm = gm.m;
    
    
    let esc = 0.0f32;

    
    
    
    gx.set_xmx(0, P7G_N, 0.0);
    gx.set_xmx(0, P7G_B, gm.xsc[P7P_N][P7P_MOVE]);
    gx.set_xmx(0, P7G_E, NEG_INF);
    gx.set_xmx(0, P7G_C, NEG_INF);
    gx.set_xmx(0, P7G_J, NEG_INF);
    for k in 0..=mm {
        gx.set_mmx(0, k, NEG_INF);
        gx.set_imx(0, k, NEG_INF);
        gx.set_dmx(0, k, NEG_INF);
    }

    
    for i in 1..=l {
        let x = dsq[i] as usize;
        
        let rscx = &gm.rsc[x];

        gx.set_mmx(i, 0, NEG_INF);
        gx.set_imx(i, 0, NEG_INF);
        gx.set_dmx(i, 0, NEG_INF);
        gx.set_xmx(i, P7G_E, NEG_INF);

        for k in 1..mm {
            
            let sc = p7_flogsum(
                p7_flogsum(
                    gx.mmx(i - 1, k - 1) + tsc(gm, P7P_MM, k - 1),
                    gx.imx(i - 1, k - 1) + tsc(gm, P7P_IM, k - 1),
                ),
                p7_flogsum(
                    gx.xmx(i - 1, P7G_B) + tsc(gm, P7P_BM, k - 1),
                    gx.dmx(i - 1, k - 1) + tsc(gm, P7P_DM, k - 1),
                ),
            );
            gx.set_mmx(i, k, sc + rscx[k * P7P_NR + P7P_MSC]);

            
            let sc = p7_flogsum(
                gx.mmx(i - 1, k) + tsc(gm, P7P_MI, k),
                gx.imx(i - 1, k) + tsc(gm, P7P_II, k),
            );
            gx.set_imx(i, k, sc + rscx[k * P7P_NR + P7P_ISC]);

            
            let d = p7_flogsum(
                gx.mmx(i, k - 1) + tsc(gm, P7P_MD, k - 1),
                gx.dmx(i, k - 1) + tsc(gm, P7P_DD, k - 1),
            );
            gx.set_dmx(i, k, d);

            
            let e = p7_flogsum(
                p7_flogsum(gx.mmx(i, k) + esc, gx.dmx(i, k) + esc),
                gx.xmx(i, P7G_E),
            );
            gx.set_xmx(i, P7G_E, e);
        }

        
        let sc = p7_flogsum(
            p7_flogsum(
                gx.mmx(i - 1, mm - 1) + tsc(gm, P7P_MM, mm - 1),
                gx.imx(i - 1, mm - 1) + tsc(gm, P7P_IM, mm - 1),
            ),
            p7_flogsum(
                gx.xmx(i - 1, P7G_B) + tsc(gm, P7P_BM, mm - 1),
                gx.dmx(i - 1, mm - 1) + tsc(gm, P7P_DM, mm - 1),
            ),
        );
        gx.set_mmx(i, mm, sc + rscx[mm * P7P_NR + P7P_MSC]);
        gx.set_imx(i, mm, NEG_INF);

        
        let d = p7_flogsum(
            gx.mmx(i, mm - 1) + tsc(gm, P7P_MD, mm - 1),
            gx.dmx(i, mm - 1) + tsc(gm, P7P_DD, mm - 1),
        );
        gx.set_dmx(i, mm, d);

        
        let e = p7_flogsum(p7_flogsum(gx.mmx(i, mm), gx.dmx(i, mm)), gx.xmx(i, P7G_E));
        gx.set_xmx(i, P7G_E, e);

        
        let j = p7_flogsum(
            gx.xmx(i - 1, P7G_J) + gm.xsc[P7P_J][P7P_LOOP],
            gx.xmx(i, P7G_E) + gm.xsc[P7P_E][P7P_LOOP],
        );
        gx.set_xmx(i, P7G_J, j);

        
        let c = p7_flogsum(
            gx.xmx(i - 1, P7G_C) + gm.xsc[P7P_C][P7P_LOOP],
            gx.xmx(i, P7G_E) + gm.xsc[P7P_E][P7P_MOVE],
        );
        gx.set_xmx(i, P7G_C, c);

        
        let n = gx.xmx(i - 1, P7G_N) + gm.xsc[P7P_N][P7P_LOOP];
        gx.set_xmx(i, P7G_N, n);

        
        let b = p7_flogsum(
            gx.xmx(i, P7G_N) + gm.xsc[P7P_N][P7P_MOVE],
            gx.xmx(i, P7G_J) + gm.xsc[P7P_J][P7P_MOVE],
        );
        gx.set_xmx(i, P7G_B, b);
    }

    
    gx.m = mm;
    gx.l = l;
    gx.xmx(l, P7G_C) + gm.xsc[P7P_C][P7P_MOVE]
}







pub fn p7_gbackward(dsq: &[u8], l: usize, gm: &Gm, gx: &mut P7Gmx) -> f32 {
    let mm = gm.m;
    
    let esc = 0.0f32;

    
    gx.set_xmx(l, P7G_J, NEG_INF);
    gx.set_xmx(l, P7G_B, NEG_INF);
    gx.set_xmx(l, P7G_N, NEG_INF);
    gx.set_xmx(l, P7G_C, gm.xsc[P7P_C][P7P_MOVE]);
    gx.set_xmx(l, P7G_E, gx.xmx(l, P7G_C) + gm.xsc[P7P_E][P7P_MOVE]);

    gx.set_mmx(l, mm, gx.xmx(l, P7G_E));
    gx.set_dmx(l, mm, gx.xmx(l, P7G_E));
    gx.set_imx(l, mm, NEG_INF);
    for k in (1..mm).rev() {
        
        let mval = p7_flogsum(gx.xmx(l, P7G_E) + esc, gx.dmx(l, k + 1) + tsc(gm, P7P_MD, k));
        gx.set_mmx(l, k, mval);
        let dval = p7_flogsum(gx.xmx(l, P7G_E) + esc, gx.dmx(l, k + 1) + tsc(gm, P7P_DD, k));
        gx.set_dmx(l, k, dval);
        gx.set_imx(l, k, NEG_INF);
    }

    
    for i in (1..l).rev() {
        let xip1 = dsq[i + 1] as usize;
        let rscx = &gm.rsc[xip1];

        
        let mut b = gx.mmx(i + 1, 1) + tsc(gm, P7P_BM, 0) + rscx[1 * P7P_NR + P7P_MSC];
        for k in 2..=mm {
            b = p7_flogsum(b, gx.mmx(i + 1, k) + tsc(gm, P7P_BM, k - 1) + rscx[k * P7P_NR + P7P_MSC]);
        }
        gx.set_xmx(i, P7G_B, b);

        
        let j = p7_flogsum(
            gx.xmx(i + 1, P7G_J) + gm.xsc[P7P_J][P7P_LOOP],
            gx.xmx(i, P7G_B) + gm.xsc[P7P_J][P7P_MOVE],
        );
        gx.set_xmx(i, P7G_J, j);

        
        let c = gx.xmx(i + 1, P7G_C) + gm.xsc[P7P_C][P7P_LOOP];
        gx.set_xmx(i, P7G_C, c);

        
        let e = p7_flogsum(
            gx.xmx(i, P7G_J) + gm.xsc[P7P_E][P7P_LOOP],
            gx.xmx(i, P7G_C) + gm.xsc[P7P_E][P7P_MOVE],
        );
        gx.set_xmx(i, P7G_E, e);

        
        let n = p7_flogsum(
            gx.xmx(i + 1, P7G_N) + gm.xsc[P7P_N][P7P_LOOP],
            gx.xmx(i, P7G_B) + gm.xsc[P7P_N][P7P_MOVE],
        );
        gx.set_xmx(i, P7G_N, n);

        
        gx.set_mmx(i, mm, gx.xmx(i, P7G_E));
        gx.set_dmx(i, mm, gx.xmx(i, P7G_E));
        gx.set_imx(i, mm, NEG_INF);

        for k in (1..mm).rev() {
            
            
            let mval = p7_flogsum(
                p7_flogsum(
                    gx.mmx(i + 1, k + 1) + tsc(gm, P7P_MM, k) + rscx[(k + 1) * P7P_NR + P7P_MSC],
                    gx.imx(i + 1, k) + tsc(gm, P7P_MI, k) + rscx[k * P7P_NR + P7P_ISC],
                ),
                p7_flogsum(
                    gx.xmx(i, P7G_E) + esc,
                    gx.dmx(i, k + 1) + tsc(gm, P7P_MD, k),
                ),
            );
            gx.set_mmx(i, k, mval);

            
            let ival = p7_flogsum(
                gx.mmx(i + 1, k + 1) + tsc(gm, P7P_IM, k) + rscx[(k + 1) * P7P_NR + P7P_MSC],
                gx.imx(i + 1, k) + tsc(gm, P7P_II, k) + rscx[k * P7P_NR + P7P_ISC],
            );
            gx.set_imx(i, k, ival);

            
            let dval = p7_flogsum(
                gx.mmx(i + 1, k + 1) + tsc(gm, P7P_DM, k) + rscx[(k + 1) * P7P_NR + P7P_MSC],
                p7_flogsum(
                    gx.dmx(i, k + 1) + tsc(gm, P7P_DD, k),
                    gx.xmx(i, P7G_E) + esc,
                ),
            );
            gx.set_dmx(i, k, dval);
        }
    }

    
    let x1 = dsq[1] as usize;
    let rsc1 = &gm.rsc[x1];
    let mut b0 = gx.mmx(1, 1) + tsc(gm, P7P_BM, 0) + rsc1[1 * P7P_NR + P7P_MSC];
    for k in 2..=mm {
        b0 = p7_flogsum(b0, gx.mmx(1, k) + tsc(gm, P7P_BM, k - 1) + rsc1[k * P7P_NR + P7P_MSC]);
    }
    gx.set_xmx(0, P7G_B, b0);
    gx.set_xmx(0, P7G_J, NEG_INF);
    gx.set_xmx(0, P7G_C, NEG_INF);
    gx.set_xmx(0, P7G_E, NEG_INF);
    let n0 = p7_flogsum(
        gx.xmx(1, P7G_N) + gm.xsc[P7P_N][P7P_LOOP],
        gx.xmx(0, P7G_B) + gm.xsc[P7P_N][P7P_MOVE],
    );
    gx.set_xmx(0, P7G_N, n0);
    for k in 1..=mm {
        gx.set_mmx(0, k, NEG_INF);
        gx.set_imx(0, k, NEG_INF);
        gx.set_dmx(0, k, NEG_INF);
    }
    gx.xmx(0, P7G_N)
}









pub fn p7_reconfig_unihit(gm: &mut Gm, l: i32) {
    gm.xsc[P7P_E][P7P_MOVE] = 0.0;
    gm.xsc[P7P_E][P7P_LOOP] = NEG_INF;
    gm.nj = 0.0;
    reconfig_length(gm, l);
}



pub fn p7_reconfig_multihit(gm: &mut Gm, l: i32) {
    let log2 = std::f64::consts::LN_2 as f32;
    gm.xsc[P7P_E][P7P_MOVE] = -log2;
    gm.xsc[P7P_E][P7P_LOOP] = -log2;
    gm.nj = 1.0;
    reconfig_length(gm, l);
}





pub fn p7_gdomain_decoding(gm: &Gm, fwd: &P7Gmx, bck: &P7Gmx) -> (Vec<f32>, Vec<f32>, Vec<f32>) {
    let l = fwd.l;
    let overall_logp = fwd.xmx(l, P7G_C) + gm.xsc[P7P_C][P7P_MOVE];
    let mut btot = vec![0.0f32; l + 1];
    let mut etot = vec![0.0f32; l + 1];
    let mut mocc = vec![0.0f32; l + 1];
    for i in 1..=l {
        
        btot[i] = btot[i - 1]
            + ((fwd.xmx(i - 1, P7G_B) + bck.xmx(i - 1, P7G_B) - overall_logp) as f64).exp() as f32;
        etot[i] = etot[i - 1]
            + ((fwd.xmx(i, P7G_E) + bck.xmx(i, P7G_E) - overall_logp) as f64).exp() as f32;
        
        let mut njcp = (fwd.xmx(i - 1, P7G_N) + bck.xmx(i, P7G_N) + gm.xsc[P7P_N][P7P_LOOP] - overall_logp).exp();
        njcp += (fwd.xmx(i - 1, P7G_J) + bck.xmx(i, P7G_J) + gm.xsc[P7P_J][P7P_LOOP] - overall_logp).exp();
        njcp += (fwd.xmx(i - 1, P7G_C) + bck.xmx(i, P7G_C) + gm.xsc[P7P_C][P7P_LOOP] - overall_logp).exp();
        mocc[i] = 1.0 - njcp;
    }
    (btot, etot, mocc)
}



pub fn p7_gdecoding(gm: &Gm, fwd: &P7Gmx, bck: &P7Gmx, pp: &mut P7Gmx) {
    let l = fwd.l;
    let mm = gm.m;
    let overall_sc = fwd.xmx(l, P7G_C) + gm.xsc[P7P_C][P7P_MOVE];
    pp.m = mm;
    pp.l = l;

    
    pp.set_xmx(0, P7G_E, 0.0);
    pp.set_xmx(0, P7G_N, 0.0);
    pp.set_xmx(0, P7G_J, 0.0);
    pp.set_xmx(0, P7G_B, 0.0);
    pp.set_xmx(0, P7G_C, 0.0);
    for k in 0..=mm {
        pp.set_mmx(0, k, 0.0);
        pp.set_imx(0, k, 0.0);
        pp.set_dmx(0, k, 0.0);
    }

    for i in 1..=l {
        let mut denom = 0.0f32;
        pp.set_mmx(i, 0, 0.0);
        pp.set_imx(i, 0, 0.0);
        pp.set_dmx(i, 0, 0.0);
        for k in 1..mm {
            let mv = (fwd.mmx(i, k) + bck.mmx(i, k) - overall_sc).exp();
            pp.set_mmx(i, k, mv);
            denom += mv;
            let iv = (fwd.imx(i, k) + bck.imx(i, k) - overall_sc).exp();
            pp.set_imx(i, k, iv);
            denom += iv;
            pp.set_dmx(i, k, 0.0);
        }
        let mvm = (fwd.mmx(i, mm) + bck.mmx(i, mm) - overall_sc).exp();
        pp.set_mmx(i, mm, mvm);
        denom += mvm;
        pp.set_imx(i, mm, 0.0);
        pp.set_dmx(i, mm, 0.0);

        pp.set_xmx(i, P7G_E, 0.0);
        let n = (fwd.xmx(i - 1, P7G_N) + bck.xmx(i, P7G_N) + gm.xsc[P7P_N][P7P_LOOP] - overall_sc).exp();
        pp.set_xmx(i, P7G_N, n);
        let j = (fwd.xmx(i - 1, P7G_J) + bck.xmx(i, P7G_J) + gm.xsc[P7P_J][P7P_LOOP] - overall_sc).exp();
        pp.set_xmx(i, P7G_J, j);
        pp.set_xmx(i, P7G_B, 0.0);
        let c = (fwd.xmx(i - 1, P7G_C) + bck.xmx(i, P7G_C) + gm.xsc[P7P_C][P7P_LOOP] - overall_sc).exp();
        pp.set_xmx(i, P7G_C, c);
        denom += n + j + c;

        let inv = 1.0 / denom;
        for k in 1..mm {
            pp.set_mmx(i, k, pp.mmx(i, k) * inv);
            pp.set_imx(i, k, pp.imx(i, k) * inv);
        }
        pp.set_mmx(i, mm, pp.mmx(i, mm) * inv);
        pp.set_xmx(i, P7G_N, pp.xmx(i, P7G_N) * inv);
        pp.set_xmx(i, P7G_J, pp.xmx(i, P7G_J) * inv);
        pp.set_xmx(i, P7G_C, pp.xmx(i, P7G_C) * inv);
    }
}



pub fn p7_gnull2_by_expectation(gm: &Gm, pp: &mut P7Gmx) -> [f32; KP] {
    let mm = gm.m;
    let ld = pp.l;
    let ncell = (mm + 1) * P7G_NSCELLS;

    
    
    let (row0_start, row1_start) = (0usize, 1 * (mm + 1) * P7G_NSCELLS);
    for c in 0..ncell {
        pp.dp[row0_start + c] = pp.dp[row1_start + c];
    }
    for s in 0..P7G_NXCELLS {
        pp.xmx[s] = pp.xmx[P7G_NXCELLS + s];
    }
    
    for i in 2..=ld {
        let rowi = i * (mm + 1) * P7G_NSCELLS;
        for c in 0..ncell {
            pp.dp[row0_start + c] += pp.dp[rowi + c];
        }
        for s in 0..P7G_NXCELLS {
            pp.xmx[s] += pp.xmx[i * P7G_NXCELLS + s];
        }
    }
    
    
    
    
    
    
    
    
    let norm = 1.0f32 / ld as f32;
    for c in 0..ncell {
        pp.dp[row0_start + c] *= norm;
    }
    for s in 0..P7G_NXCELLS {
        pp.xmx[s] *= norm;
    }

    
    let xfactor = pp.xmx(0, P7G_N) + pp.xmx(0, P7G_C) + pp.xmx(0, P7G_J);

    let mut null2 = [0.0f32; KP];
    for x in 0..K {
        
        let mut sv = 0.0f32;
        for k in 1..mm {
            sv += pp.mmx(0, k) * gm.rsc[x][k * P7P_NR + P7P_MSC].exp();
            sv += pp.imx(0, k); 
        }
        sv += pp.mmx(0, mm) * gm.rsc[x][mm * P7P_NR + P7P_MSC].exp();
        null2[x] = sv + xfactor;
    }
    
    for x in (K + 1)..=(KP - 3) {
        let set = degen_set(x);
        let mut result = 0.0f32;
        let mut ndegen = 0.0f32;
        for &i in set {
            result += null2[i];
            ndegen += 1.0;
        }
        null2[x] = if ndegen > 0.0 { result / ndegen } else { 0.0 };
    }
    null2[K] = 1.0; 
    null2[KP - 2] = 1.0; 
    null2[KP - 1] = 1.0; 
    null2
}


#[derive(Clone, Debug)]
pub struct Domain {
    pub ienv: i64,
    pub jenv: i64,
    pub envsc: f32,
    pub domcorrection: f32,
    /// `dom->ad->hmmfrom` / `dom->ad->hmmto` — profile coordinates of the first and
    /// last M state of this domain's optimal-accuracy alignment. These are columns
    /// 16/17 of `--domtblout` (p7_tophits.c:1778-1779).
    pub hmm_from: i32,
    pub hmm_to: i32,
    /// `dom->iali` / `dom->jali` (= `ad->sqfrom`/`ad->sqto`, p7_domaindef.c:955-956),
    /// in original-dsq coordinates. Columns 18/19 of `--domtblout`.
    pub ali_from: i64,
    pub ali_to: i64,
    /// `dom->oasc` (p7_domaindef.c:960) — expected number of correctly aligned
    /// residues; the `acc` column is `oasc / (1 + |jenv - ienv|)`.
    pub oasc: f32,
    /// The optimal-accuracy traceback, in original-dsq coordinates. Feeds
    /// `p7_alidisplay_Create()`; the tabular outputs do not read it.
    pub trace: Vec<crate::optacc::TracePos>,
}


fn is_multidomain_region(btot: &[f32], etot: &[f32], i: usize, j: usize, rt3: f32) -> bool {
    let mut max = -1.0f32;
    for z in i..=j {
        let e = etot[z] - etot[i - 1];
        let b = btot[j] - btot[z - 1];
        let expected_n = if e < b { e } else { b };
        if expected_n > max {
            max = expected_n;
        }
    }
    max >= rt3
}






fn rescore_isolated_domain(
    ff: &ForwardFilter,
    dsq: &[u8],
    i: usize,
    j: usize,
    save_l: usize,
    do_null2: bool,
    precomputed_domcorrection: Option<f32>,
    n2sc: &mut [f32],
) -> Option<Domain> {
    
    
    
    
    
    
    
    // `None` is C's `eslERANGE` from `p7_Decoding`: the domain is abandoned as
    // "repetitive garbage" (p7_domaindef.c:846-852). The envelope stays counted.
    let er = crate::omx::envelope_rescore(ff, dsq, i, j, save_l)?;
    let (envsc, null2) = (er.envsc, er.null2);

    
    
    
    
    
    
    
    let domcorrection = if !do_null2 {
        0.0
    } else if let Some(dc) = precomputed_domcorrection {
        
        
        dc
    } else {
        
        
        
        let mut d = 0.0f32;
        for pos in i..=j {
            let v = null2[dsq[pos] as usize].ln();
            n2sc[pos] = v;
            d += v;
        }
        d
    };

    Some(Domain {
        ienv: i as i64,
        jenv: j as i64,
        envsc,
        domcorrection,
        hmm_from: er.hmm_from,
        hmm_to: er.hmm_to,
        // "hack the trace's sq coords to be correct w.r.t. original dsq"
        // (p7_domaindef.c:857-858): tr->i[z] += i-1, then dom->iali = ad->sqfrom.
        ali_from: er.ali_from + i as i64 - 1,
        ali_to: er.ali_to + i as i64 - 1,
        oasc: er.oasc,
        // "hack the trace's sq coords to be correct w.r.t. original dsq"
        // (p7_domaindef.c:856-858): tr->i[z] += i-1 for every emitting position.
        trace: {
            let mut t = er.trace;
            let off = i as i32 - 1;
            for p in t.iter_mut() {
                if p.2 > 0 {
                    p.2 += off;
                }
            }
            t
        },
    })
}
























struct EslRng {
    x: u32,
}

pub(crate) fn esl_mix3(mut a: u32, mut b: u32, mut c: u32) -> u32 {
    a = a.wrapping_sub(b); a = a.wrapping_sub(c); a ^= c >> 13;
    b = b.wrapping_sub(c); b = b.wrapping_sub(a); b ^= a << 8;
    c = c.wrapping_sub(a); c = c.wrapping_sub(b); c ^= b >> 13;
    a = a.wrapping_sub(b); a = a.wrapping_sub(c); a ^= c >> 12;
    b = b.wrapping_sub(c); b = b.wrapping_sub(a); b ^= a << 16;
    c = c.wrapping_sub(a); c = c.wrapping_sub(b); c ^= b >> 5;
    a = a.wrapping_sub(b); a = a.wrapping_sub(c); a ^= c >> 3;
    b = b.wrapping_sub(c); b = b.wrapping_sub(a); b ^= a << 10;
    c = c.wrapping_sub(a); c = c.wrapping_sub(b); c ^= b >> 15;
    c
}
impl EslRng {
    
    fn new(seed: u32) -> Self {
        let mut x = esl_mix3(seed, 87654321, 12345678);
        if x == 0 { x = 42; }
        EslRng { x }
    }
    
    fn knuth(&mut self) -> u32 {
        self.x = self.x.wrapping_mul(69069).wrapping_add(1);
        self.x
    }
    
    fn random(&mut self) -> f64 {
        (self.knuth() as f64) / 4294967296.0
    }
    
    fn fchoose(&mut self, p: &[f32]) -> usize {
        let mut norm = 0.0f64;
        let mut sum = 0.0f64;
        let roll = self.random();
        for &v in p { norm += v as f64; }
        for (i, &v) in p.iter().enumerate() {
            sum += v as f64;
            if roll < sum / norm { return i; }
        }
        p.len() - 1 
    }
}


pub(crate) fn esl_vec_fsum(v: &[f32]) -> f32 {
    let mut sum = 0.0f32;
    let mut c = 0.0f32;
    for &vi in v {
        let y = vi - c;
        let t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }
    sum
}

fn esl_vec_flogsum(v: &[f32]) -> f32 {
    let max = v.iter().cloned().fold(f32::NEG_INFINITY, f32::max);
    if max == f32::INFINITY { return f32::INFINITY; }
    let mut sum = 0.0f32;
    for &x in v {
        if x > max - 50.0 { sum += (x - max).exp(); }
    }
    sum.ln() + max
}


fn esl_vec_fnorm(v: &mut [f32]) {
    let sum = esl_vec_fsum(v);
    if sum != 0.0 {
        for x in v.iter_mut() {
            *x /= sum;
        }
    } else {
        let n = v.len() as f32;
        for x in v.iter_mut() {
            *x = 1.0 / n;
        }
    }
}

fn esl_vec_flognorm(v: &mut [f32]) {
    let denom = esl_vec_flogsum(v);
    for x in v.iter_mut() { *x += -denom; }
    for x in v.iter_mut() { *x = x.exp(); }
    let sum = esl_vec_fsum(v);
    if sum != 0.0 {
        for x in v.iter_mut() { *x /= sum; }
    } else {
        let n = v.len() as f32;
        for x in v.iter_mut() { *x = 1.0 / n; }
    }
}


pub const ST_T: u8 = 0;
pub const ST_C: u8 = 1;
pub const ST_E: u8 = 2;
pub const ST_M: u8 = 3;
pub const ST_D: u8 = 4;
pub const ST_I: u8 = 5;
pub const ST_N: u8 = 6;
pub const ST_B: u8 = 7;
pub const ST_J: u8 = 8;
pub const ST_S: u8 = 9;







#[allow(dead_code)]
fn p7_gstochastic_trace(rng: &mut EslRng, lr: usize, gm: &Gm, gx: &P7Gmx) -> Vec<(u8, i32, i32)> {
    let m = gm.m as i32;
    let mut tr: Vec<(u8, i32, i32)> = Vec::new();
    let mut k: i32 = 0;
    let mut i: i32 = lr as i32;
    tr.push((ST_T, k, i));
    tr.push((ST_C, k, i));
    let mut sprv = ST_C;
    while sprv != ST_S {
        let last = tr[tr.len() - 1].0;
        let iu = i as usize;
        let im1 = (i - 1) as usize;
        let ku = k as usize;
        let km1 = (k - 1) as usize;
        let mu = m as usize;
        let scur: u8;
        match last {
            ST_C => {
                
                let mut sc = [
                    gx.xmx(im1, P7G_C) + gm.xsc[P7P_C][P7P_LOOP],
                    gx.xmx(iu, P7G_E) + gm.xsc[P7P_E][P7P_MOVE],
                ];
                esl_vec_flognorm(&mut sc);
                scur = if rng.fchoose(&sc) == 0 { ST_C } else { ST_E };
            }
            ST_E => {
                
                
                
                
                let mut sc = vec![f32::NEG_INFINITY; 2 * mu + 1];
                for kk in 1..=mu { sc[kk] = gx.mmx(iu, kk); }
                for kk in 2..=mu { sc[kk + mu] = gx.dmx(iu, kk); }
                esl_vec_flognorm(&mut sc);
                let choice = rng.fchoose(&sc);
                if choice <= mu {
                    k = choice as i32;
                    scur = ST_M;
                } else {
                    k = (choice - mu) as i32;
                    scur = ST_D;
                }
            }
            ST_M => {
                
                let mut sc = [
                    gx.xmx(im1, P7G_B) + tsc(gm, P7P_BM, km1),
                    gx.mmx(im1, km1) + tsc(gm, P7P_MM, km1),
                    gx.imx(im1, km1) + tsc(gm, P7P_IM, km1),
                    gx.dmx(im1, km1) + tsc(gm, P7P_DM, km1),
                ];
                esl_vec_flognorm(&mut sc);
                scur = match rng.fchoose(&sc) {
                    0 => ST_B,
                    1 => ST_M,
                    2 => ST_I,
                    _ => ST_D,
                };
                k -= 1;
                i -= 1;
            }
            ST_D => {
                
                let mut sc = [
                    gx.mmx(iu, km1) + tsc(gm, P7P_MD, km1),
                    gx.dmx(iu, km1) + tsc(gm, P7P_DD, km1),
                ];
                esl_vec_flognorm(&mut sc);
                scur = if rng.fchoose(&sc) == 0 { ST_M } else { ST_D };
                k -= 1;
            }
            ST_I => {
                
                let mut sc = [
                    gx.mmx(im1, ku) + tsc(gm, P7P_MI, ku),
                    gx.imx(im1, ku) + tsc(gm, P7P_II, ku),
                ];
                esl_vec_flognorm(&mut sc);
                scur = if rng.fchoose(&sc) == 0 { ST_M } else { ST_I };
                i -= 1;
            }
            ST_N => {
                
                scur = if i == 0 { ST_S } else { ST_N };
            }
            ST_B => {
                
                let mut sc = [
                    gx.xmx(iu, P7G_N) + gm.xsc[P7P_N][P7P_MOVE],
                    gx.xmx(iu, P7G_J) + gm.xsc[P7P_J][P7P_MOVE],
                ];
                esl_vec_flognorm(&mut sc);
                scur = if rng.fchoose(&sc) == 0 { ST_N } else { ST_J };
            }
            ST_J => {
                
                let mut sc = [
                    gx.xmx(im1, P7G_J) + gm.xsc[P7P_J][P7P_LOOP],
                    gx.xmx(iu, P7G_E) + gm.xsc[P7P_E][P7P_LOOP],
                ];
                esl_vec_flognorm(&mut sc);
                scur = if rng.fchoose(&sc) == 0 { ST_J } else { ST_E };
            }
            _ => unreachable!("bogus state in stochastic traceback"),
        }
        tr.push((scur, k, i));
        
        if (scur == ST_N || scur == ST_J || scur == ST_C) && scur == sprv {
            i -= 1;
        }
        sprv = scur;
    }
    tr.reverse();
    tr
}


pub(crate) fn trace_index(tr: &[(u8, i32, i32)]) -> Vec<(i32, i32, i32, i32)> {
    let mut doms: Vec<(i32, i32, i32, i32)> = Vec::new(); 
    for &(st, k, i) in tr {
        match st {
            ST_B => doms.push((0, 0, 0, 0)),
            ST_M => {
                let d = doms.last_mut().unwrap();
                if d.0 == 0 { d.0 = i; }
                if d.2 == 0 { d.2 = k; }
                d.1 = i;
                d.3 = k;
            }
            _ => {}
        }
    }
    doms
}


#[derive(Clone, Copy)]
struct Seg { idx: i32, i: i32, j: i32, k: i32, m: i32 }


fn link_spsamples(h1: &Seg, h2: &Seg) -> bool {
    let min_overlap = 0.8f32;
    let max_diagdiff = 4i32;
    
    let nov = h1.j.min(h2.j) - h1.i.max(h2.i) + 1;
    let n = (h1.j - h1.i + 1).min(h2.j - h2.i + 1);
    if (nov as f32) / (n as f32) < min_overlap { return false; }
    
    let nov = h1.m.min(h2.m) - h1.k.max(h2.k);
    let n = (h1.m - h1.k + 1).min(h2.m - h2.k + 1);
    if (nov as f32) / (n as f32) < min_overlap { return false; }
    
    let (d1, d2) = (h1.i - h1.k, h2.i - h2.k);
    if (d1 - d2).abs() <= max_diagdiff { return true; }
    let (d1, d2) = (h1.j - h1.m, h2.j - h2.m);
    if (d1 - d2).abs() <= max_diagdiff { return true; }
    false
}


fn single_linkage(segs: &[Seg]) -> (Vec<usize>, usize) {
    let n = segs.len();
    let mut a: Vec<usize> = (0..n).map(|v| n - v - 1).collect(); 
    let mut na = n;
    let mut b = vec![0usize; n];
    let mut nb = 0usize;
    let mut c = vec![0usize; n];
    let mut nc = 0usize;
    while na > 0 {
        let v = a[na - 1]; na -= 1;
        b[nb] = v; nb += 1;
        while nb > 0 {
            let v = b[nb - 1]; nb -= 1;
            c[v] = nc;
            let mut i: i64 = na as i64 - 1;
            while i >= 0 {
                if link_spsamples(&segs[v], &segs[a[i as usize]]) {
                    let w = a[i as usize];
                    a[i as usize] = a[na - 1];
                    na -= 1;
                    b[nb] = w; nb += 1;
                }
                i -= 1;
            }
        }
        nc += 1;
    }
    (c, nc)
}

fn iargmax(v: &[i32]) -> usize {
    let mut best = 0usize;
    for i in 1..v.len() {
        if v[i] > v[best] { best = i; }
    }
    best
}





fn spensemble_cluster(segs: &[Seg], nsamples: usize) -> Vec<(i32, i32, f32)> {
    let min_posterior = 0.25f32;
    let min_endpointp = 0.02f32;
    if segs.is_empty() { return Vec::new(); }
    let (assignment, nc) = single_linkage(segs);
    let mut sigc: Vec<(i32, i32, i32, i32, f32)> = Vec::new(); 
    for c in 0..nc {
        
        let mut ninc = 0i32;
        let mut idx_of_last = -1i32;
        for h in 0..segs.len() {
            if assignment[h] == c {
                if segs[h].idx != idx_of_last { ninc += 1; }
                idx_of_last = segs[h].idx;
            }
        }
        if (ninc as f32) / (nsamples as f32) < min_posterior { continue; }
        
        let (mut imin, mut imax, mut jmin, mut jmax, mut kmin, mut kmax, mut mmin, mut mmax) =
            (0i32, 0i32, 0i32, 0i32, 0i32, 0i32, 0i32, 0i32);
        let mut started = false;
        for h in 0..segs.len() {
            if assignment[h] != c { continue; }
            let s = &segs[h];
            if !started {
                imin = s.i; imax = s.i; jmin = s.j; jmax = s.j;
                kmin = s.k; kmax = s.k; mmin = s.m; mmax = s.m;
                started = true;
            } else {
                imin = imin.min(s.i); imax = imax.max(s.i);
                jmin = jmin.min(s.j); jmax = jmax.max(s.j);
                kmin = kmin.min(s.k); kmax = kmax.max(s.k);
                mmin = mmin.min(s.m); mmax = mmax.max(s.m);
            }
        }
        let epc_threshold = ((ninc as f32) * min_endpointp).ceil() as i32;
        
        let mut epc = vec![0i32; (imax - imin + 1) as usize];
        for h in 0..segs.len() { if assignment[h] == c { epc[(segs[h].i - imin) as usize] += 1; } }
        let mut best_i = imin;
        while best_i <= imax { if epc[(best_i - imin) as usize] >= epc_threshold { break; } best_i += 1; }
        if best_i > imax { best_i = imin + iargmax(&epc) as i32; }
        
        let mut epc = vec![0i32; (kmax - kmin + 1) as usize];
        for h in 0..segs.len() { if assignment[h] == c { epc[(segs[h].k - kmin) as usize] += 1; } }
        let mut best_k = kmin;
        while best_k <= kmax { if epc[(best_k - kmin) as usize] >= epc_threshold { break; } best_k += 1; }
        if best_k > kmax { best_k = kmin + iargmax(&epc) as i32; }
        
        let mut epc = vec![0i32; (jmax - jmin + 1) as usize];
        for h in 0..segs.len() { if assignment[h] == c { epc[(segs[h].j - jmin) as usize] += 1; } }
        let mut best_j = jmax;
        while best_j >= jmin { if epc[(best_j - jmin) as usize] >= epc_threshold { break; } best_j -= 1; }
        if best_j < jmin { best_j = jmin + iargmax(&epc) as i32; }
        
        let mut epc = vec![0i32; (mmax - mmin + 1) as usize];
        for h in 0..segs.len() { if assignment[h] == c { epc[(segs[h].m - mmin) as usize] += 1; } }
        let mut best_m = mmax;
        while best_m >= mmin { if epc[(best_m - mmin) as usize] >= epc_threshold { break; } best_m -= 1; }
        if best_m < mmin { best_m = mmin + iargmax(&epc) as i32; }

        if best_i > best_j || best_k > best_m { continue; }
        sigc.push((best_i, best_j, best_k, best_m, ninc as f32 / nsamples as f32));
    }
    
    sigc.sort_by(|a, b| a.0.cmp(&b.0));
    
    let nsig = sigc.len();
    let mut dominated = vec![false; nsig];
    for d in 0..nsig {
        for d2 in (d + 1)..nsig {
            let nov = sigc[d].1.min(sigc[d2].1) - sigc[d].0.max(sigc[d2].0) + 1;
            if nov == 0 { break; }
            let n = (sigc[d].1 - sigc[d].0 + 1).min(sigc[d2].1 - sigc[d2].0 + 1);
            if (nov as f32) / (n as f32) >= 0.8 {
                if sigc[d].4 > sigc[d2].4 { dominated[d2] = true; } else { dominated[d] = true; }
            }
        }
    }
    sigc.iter().enumerate()
        .filter(|(d, _)| !dominated[*d])
        .map(|(_, s)| (s.0, s.1, s.4))
        .collect()
}











fn p7_stochastic_trace_odds(rng: &mut EslRng, ff: &ForwardFilter, xf: &XFactors, ox: &Omx) -> Vec<(u8, i32, i32)> {
    let ld = ox.ld as i32;
    let mut tr: Vec<(u8, i32, i32)> = Vec::new();
    let mut k: i32 = 0;
    let mut i: i32 = ld;
    tr.push((ST_T, k, i));
    tr.push((ST_C, k, i));
    let mut s0 = ST_C;
    while s0 != ST_S {
        let iu = i as usize;
        let im1 = (i - 1) as usize;
        let ku = k as usize;
        let s1: u8;
        match s0 {
            ST_M => {
                
                let mut path = [
                    ox.xb[im1] * ff.tbm[ku],
                    ox.mmx[im1][ku - 1] * ff.amm[ku],
                    ox.imx[im1][ku - 1] * ff.aim[ku],
                    ox.dmx[im1][ku - 1] * ff.adm[ku],
                ];
                esl_vec_fnorm(&mut path);
                s1 = [ST_B, ST_M, ST_I, ST_D][rng.fchoose(&path)];
                k -= 1;
                i -= 1;
            }
            ST_D => {
                
                let mut path = [ox.mmx[iu][ku - 1] * ff.tmd[ku - 1], ox.dmx[iu][ku - 1] * ff.tdd[ku - 1]];
                esl_vec_fnorm(&mut path);
                s1 = if rng.fchoose(&path) == 0 { ST_M } else { ST_D };
                k -= 1;
            }
            ST_I => {
                
                let mut path = [ox.mmx[im1][ku] * ff.tmi[ku], ox.imx[im1][ku] * ff.tii[ku]];
                esl_vec_fnorm(&mut path);
                s1 = if rng.fchoose(&path) == 0 { ST_M } else { ST_I };
                i -= 1;
            }
            ST_N => {
                
                s1 = if i == 0 { ST_S } else { ST_N };
            }
            ST_C => {
                
                let mut path = [ox.xc[im1] * xf.c_loop, ox.xe[iu] * xf.e_move * ox.scale[iu]];
                esl_vec_fnorm(&mut path);
                s1 = if rng.fchoose(&path) == 0 { ST_C } else { ST_E };
            }
            ST_J => {
                
                let mut path = [ox.xj[im1] * xf.j_loop, ox.xe[iu] * xf.e_loop * ox.scale[iu]];
                esl_vec_fnorm(&mut path);
                s1 = if rng.fchoose(&path) == 0 { ST_J } else { ST_E };
            }
            ST_E => {
                
                let (st, kk) = select_e_odds(rng, ox, iu);
                s1 = st;
                k = kk;
            }
            ST_B => {
                
                let mut path = [ox.xn[iu] * xf.n_move, ox.xj[iu] * xf.j_move];
                esl_vec_fnorm(&mut path);
                s1 = if rng.fchoose(&path) == 0 { ST_N } else { ST_J };
            }
            _ => unreachable!("bogus state in odds stochastic traceback"),
        }
        tr.push((s1, k, i));
        
        if (s1 == ST_N || s1 == ST_J || s1 == ST_C) && s1 == s0 {
            i -= 1;
        }
        s0 = s1;
    }
    tr.reverse();
    tr
}







fn select_e_odds(rng: &mut EslRng, ox: &Omx, i: usize) -> (u8, i32) {
    let m = ox.m;
    let q_n = (m + 3) / 4; 
    let norm = (1.0f64 / ox.xe[i] as f64) as f32; 
    let roll = rng.random(); 
    let mut sum = 0.0f64;
    
    
    for _pass in 0..4 {
        for q in 0..q_n {
            for r in 0..4 {
                let k = r * q_n + q + 1;
                let v = if k <= m { ox.mmx[i][k] } else { 0.0 };
                sum += (v * norm) as f64;
                if roll < sum {
                    return (ST_M, k as i32);
                }
            }
            for r in 0..4 {
                let k = r * q_n + q + 1;
                let v = if k <= m { ox.dmx[i][k] } else { 0.0 };
                sum += (v * norm) as f64;
                if roll < sum {
                    return (ST_D, k as i32);
                }
            }
        }
    }
    (ST_M, m as i32) 
}








fn compute_bytrace_null2(counter: &[f32], ld: f32, xsum: f32, rfv: &[[f32; KP]], m: usize) -> [f32; KP] {
    // C, impl_sse/null2.c (`p7_Null2_ByTrace`), the same striped sum
    // `p7_Null2_ByExpectation` uses:
    //   norm = 1.0 / (float) Ld;
    //   ...
    //   xfactor =  XMXo(0,p7X_N) + XMXo(0,p7X_C) + XMXo(0,p7X_J);
    //   for (x = 0; x < om->abc->K; x++)
    //     {
    //       sv = _mm_setzero_ps();
    //       rp = om->rfv[x];
    //       for (q = 0; q < Q; q++)
    //         {
    //           sv = _mm_add_ps(sv, _mm_mul_ps(wrk->dpf[0][q*3 + p7X_M], *rp)); rp++;
    //           sv = _mm_add_ps(sv,            wrk->dpf[0][q*3 + p7X_I]); /* insert emission odds implicitly 1.0 */
    //         }
    //       esl_sse_hsum_ps(sv, &(null2[x]));
    //       null2[x] += xfactor;
    //     }
    //
    // Four independent striped lanes, q ascending, collapsed as `(a0+a1) + (a2+a3)`.
    // The I cells stay zero because C bins *both* M and I counts into the M cell --
    // `q = p7X_NSCELLS * ((tr->k[z] - 1) % Q) + p7X_M`, with the state-type line left
    // commented out right above it -- so `counter` merging them matches, and the
    // emission odds end up applied to insert counts too.
    //
    // `xfactor` is each of the three normalized special counts added in N, C, J order.
    // Within a domain (`zstart`/`zend` are the trace's B and E) there are no N/C/J
    // emissions, so it is always 0 here; it is written this way to stay faithful.
    let norm = 1.0 / ld;
    let xfactor = xsum * norm;
    let mut null2 = [0.0f32; KP];
    let q_n = crate::optacc::nqf(m);
    for x in 0..K {
        let mut lane = [0.0f32; 4];
        for q in 0..q_n {
            for (z, l) in lane.iter_mut().enumerate() {
                let k = z * q_n + q + 1;
                if k <= m {
                    *l += (counter[k] * norm) * rfv[k][x];
                }
            }
        }
        null2[x] = ((lane[0] + lane[1]) + (lane[2] + lane[3])) + xfactor;
    }
    
    for x in (K + 1)..=(KP - 3) {
        let set = degen_set(x);
        let mut r = 0.0f32;
        let mut n = 0.0f32;
        for &y in set {
            r += null2[y];
            n += 1.0;
        }
        null2[x] = if n > 0.0 { r / n } else { 0.0 };
    }
    null2[K] = 1.0; 
    null2[KP - 2] = 1.0; 
    null2[KP - 1] = 1.0; 
    null2
}







fn trace_domain_null2s(tr: &[(u8, i32, i32)], rfv: &[[f32; KP]], m: usize) -> Vec<(i32, i32, [f32; KP])> {
    let mut out = Vec::new();
    let mut in_dom = false;
    let mut counter = vec![0.0f32; m + 1];
    let mut ld = 0.0f32;
    let mut xsum = 0.0f32; 
    let mut sqfrom = 0i32;
    let mut sqto = 0i32;
    for &(st, k, i) in tr {
        if st == ST_B {
            in_dom = true;
            for c in counter.iter_mut() {
                *c = 0.0;
            }
            ld = 0.0;
            xsum = 0.0;
            sqfrom = 0;
            sqto = 0;
        } else if in_dom {
            match st {
                ST_M => {
                    ld += 1.0;
                    if k > 0 {
                        counter[k as usize] += 1.0;
                    }
                    if sqfrom == 0 {
                        sqfrom = i;
                    }
                    sqto = i;
                }
                ST_I => {
                    if i > 0 {
                        ld += 1.0;
                        if k > 0 {
                            counter[k as usize] += 1.0;
                        }
                    }
                }
                ST_N | ST_C | ST_J => {
                    if i > 0 {
                        xsum += 1.0;
                    }
                }
                ST_E => {
                    let null2 = compute_bytrace_null2(&counter, ld, xsum, rfv, m);
                    out.push((sqfrom, sqto, null2));
                    in_dom = false;
                }
                _ => {}
            }
        }
    }
    out
}












fn region_trace_ensemble(
    ff: &ForwardFilter,
    dsq: &[u8],
    ri: usize,
    j: usize,
    save_l: i32,
    seed: u32,
) -> (Vec<(usize, usize)>, Vec<f32>) {
    let nsamples = 200usize;
    // "Normally, we reinitialize the RNG to the original seed every time we're about
    // to collect a stochastic trace ensemble. This eliminates run-to-run variability."
    // (p7_pipeline.c:120-123; the reinit itself is p7_domaindef.c:612-613.) A fresh
    // RNG per region is exactly that reseeding. `--seed 0` turns reseeding off in C
    // and carries RNG state across regions; that mode is non-reproducible by
    // construction, and `seed` then simply carries the one-time arbitrary seed.
    let lr = j - ri + 1;

    
    
    let xf = XFactors::multihit(save_l as usize);

    
    
    let mut region = vec![255u8];
    region.extend_from_slice(&dsq[ri..=j]);
    region.push(255u8);
    let ox = omx::forward(&ff, &xf, &region, lr);

    
    let mut n2sc = vec![0.0f32; lr + 1];

    
    
    let mut rng = EslRng::new(seed);
    let mut segs: Vec<Seg> = Vec::new();
    for t in 0..nsamples {
        let tr = p7_stochastic_trace_odds(&mut rng, &ff, &xf, &ox);

        
        let doms = trace_index(&tr);
        for (sqfrom, sqto, hmmfrom, hmmto) in doms {
            segs.push(Seg {
                idx: t as i32,
                i: sqfrom + ri as i32 - 1,
                j: sqto + ri as i32 - 1,
                k: hmmfrom,
                m: hmmto,
            });
        }

        
        
        
        
        let dnull2 = trace_domain_null2s(&tr, &ff.rfv, ff.m);
        let mut pos = 1i32;
        for (sqfrom, sqto, null2) in &dnull2 {
            while pos <= *sqfrom {
                n2sc[pos as usize] += 1.0;
                pos += 1;
            }
            while pos <= *sqto {
                let res = region[pos as usize] as usize; 
                n2sc[pos as usize] += null2[res];
                pos += 1;
            }
        }
        while pos <= lr as i32 {
            n2sc[pos as usize] += 1.0;
            pos += 1;
        }
    }

    
    for p in 1..=lr {
        n2sc[p] = (n2sc[p] / nsamples as f32).ln();
    }

    let clusters = spensemble_cluster(&segs, nsamples);
    let env = clusters
        .into_iter()
        .map(|(i2, j2, _p)| (i2 as usize, j2 as usize))
        .collect();
    (env, n2sc)
}







pub struct DomainDef {
    pub domains: Vec<Domain>,
    pub nexpected: f32, 
    pub nregions: i32,
    pub nclustered: i32,
    pub noverlaps: i32,
    pub nenvelopes: i32,
    
    
    
    
    pub n2sc: Vec<f32>,
}

pub fn p7_domaindef_local(
    ff: &ForwardFilter,
    gm: &mut Gm,
    dsq: &[u8],
    l: usize,
    do_null2: bool,
    rng_seed: u32,
) -> DomainDef {
    let rt1 = 0.25f32;
    let rt2 = 0.10f32;
    let rt3 = 0.20f32;
    let save_l = gm.l; 

    
    
    
    
    
    // C, p7_domaindef.c:400: `p7_DomainDecoding`'s eslERANGE propagates straight out of
    // `p7_domaindef_ByPosteriorHeuristics`, before any region is examined, so nothing
    // about this target is reported. An empty `DomainDef` gives the same result, since
    // `search_one` returns no hit when there are no domains.
    let Some((btot, etot, mocc)) = crate::omx::domain_decoding(ff, dsq, l) else {
        p7_reconfig_unihit(gm, save_l);
        p7_reconfig_multihit(gm, save_l);
        return DomainDef {
            domains: Vec::new(),
            nexpected: 0.0,
            nregions: 0,
            nclustered: 0,
            noverlaps: 0,
            nenvelopes: 0,
            n2sc: vec![0.0f32; l + 1],
        };
    };
    let nexpected = btot[l]; 

    
    p7_reconfig_unihit(gm, save_l);

    let mut domains = Vec::new();
    
    
    let mut n2sc_full = vec![0.0f32; l + 1];
    let (mut nregions, mut nclustered, mut noverlaps, mut nenvelopes) = (0i32, 0i32, 0i32, 0i32);
    
    let mut i: i64 = -1;
    let mut triggered = false;
    for j in 1..=l {
        if !triggered {
            if mocc[j] - (btot[j] - btot[j - 1]) < rt2 {
                i = j as i64;
            } else if i == -1 {
                i = j as i64;
            }
            if mocc[j] >= rt1 {
                triggered = true;
            }
        } else if mocc[j] - (etot[j] - etot[j - 1]) < rt2 {
            
            nregions += 1;
            let ri = i as usize;
            if is_multidomain_region(&btot, &etot, ri, j, rt3) {
                
                
                nclustered += 1; 

                
                
                
                let (clusters, region_n2sc) =
                    region_trace_ensemble(ff, dsq, ri, j, save_l, rng_seed);

                
                
                
                
                let lr = j - ri + 1;
                for p in 1..=lr {
                    n2sc_full[ri + p - 1] = region_n2sc[p];
                }

                
                let mut last_j2 = 0usize;
                for (i2, j2) in clusters {
                    if i2 <= last_j2 {
                        noverlaps += 1; 
                    }
                    nenvelopes += 1; 
                    
                    
                    let mut dc = 0.0f32;
                    for p in i2..=j2 {
                        dc += region_n2sc[p - ri + 1];
                    }
                    // C, p7_domaindef.c:471-472:
                    //   if (rescore_isolated_domain(...) == eslOK)
                    //        last_j2 = j2;
                    // so a domain lost to numeric overflow also leaves `last_j2` where
                    // it was, and the next cluster is not counted as overlapping it.
                    if let Some(dom) = rescore_isolated_domain(
                        ff, dsq, i2, j2, save_l as usize, do_null2, Some(dc), &mut n2sc_full,
                    ) {
                        domains.push(dom);
                        last_j2 = j2;
                    }
                }
            } else {
                
                
                // C, p7_domaindef.c:480-481: the single-domain path increments
                // `nenvelopes` and then ignores the rescore's return value, so an
                // overflowed domain leaves `env` > `dom` in the output counts.
                nenvelopes += 1;
                if let Some(dom) = rescore_isolated_domain(
                    ff, dsq, ri, j, save_l as usize, do_null2, None, &mut n2sc_full,
                ) {
                    domains.push(dom);
                }
            }
            i = -1;
            triggered = false;
        }
    }

    
    p7_reconfig_multihit(gm, save_l);
    DomainDef {
        domains,
        nexpected,
        nregions,
        nclustered,
        noverlaps,
        nenvelopes,
        n2sc: n2sc_full,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seqio::read_fasta;

    #[test]
    fn gforward_smoke() {
        let hmm = P7Hmm::read_all(&format!("{}/testdata/globins4.hmm", env!("CARGO_MANIFEST_DIR")))
            .unwrap()
            .pop()
            .unwrap();
        let seqs = read_fasta(&format!("{}/testdata/globins45.fa", env!("CARGO_MANIFEST_DIR")))
            .unwrap();
        let s = &seqs[0];
        let l = s.len();
        let gm = build_local_profile(&hmm, l);
        let mut gx = P7Gmx::new(gm.m, l);
        let sc = p7_gforward(&s.dsq, l, &gm, &mut gx);
        assert!(sc.is_finite(), "Forward score not finite: {sc}");
        assert!(sc > 0.0, "Forward score not positive-ish: {sc}");
    }

    
    
    
    #[test]
    fn fn3_7less_multidomain() {
        use crate::alphabet::AminoMap;
        use crate::pipeline::Model;
        use crate::seqio::Seq;

        
        let path = format!("{}/testdata/7LESS_DROME", env!("CARGO_MANIFEST_DIR"));
        let text = std::fs::read_to_string(&path).unwrap();
        let mut name = String::new();
        let mut residues: Vec<u8> = Vec::new();
        let mut in_sq = false;
        for line in text.lines() {
            if let Some(rest) = line.strip_prefix("ID   ") {
                name = rest.split_whitespace().next().unwrap_or("").to_string();
            } else if line.starts_with("SQ") {
                in_sq = true;
            } else if line.starts_with("//") {
                in_sq = false;
            } else if in_sq {
                residues.extend(line.bytes().filter(|b| b.is_ascii_alphabetic()));
            }
        }
        assert_eq!(name, "7LESS_DROME");
        assert_eq!(residues.len(), 2554, "expected 2554 AA");
        let seq = Seq {
            name,
            acc: String::new(),
            desc: String::new(),
            dsq: AminoMap::new().digitize(&residues),
            roff: 0,
            eoff: 0,
        };

        
        let hmm = P7Hmm::read_all(&format!("{}/testdata/fn3.hmm", env!("CARGO_MANIFEST_DIR")))
            .unwrap()
            .pop()
            .unwrap();
        let model = Model::new(hmm);
        let hit = model.search_one(&seq, &crate::pli::Pipeline::default()).expect("7LESS_DROME should be a hit");

        let z = 1.0_f64; 
        let e = hit.lnp.exp() * z;
        eprintln!(
            "fn3 × 7LESS_DROME: score={:.1} E={:.2e} ndom={} reg={} clu={} ov={} env={}",
            hit.score, e, hit.ndom, hit.nregions, hit.nclustered, hit.noverlaps, hit.nenvelopes
        );
        for (di, d) in hit.domains.iter().enumerate() {
            eprintln!("  dom {:>2}: env {}..{}  bits={:.1}", di + 1, d.ienv, d.jenv, d.bitscore);
        }

        
        
        assert!(hit.nregions >= 1, "expected ≥1 region");
        assert!(hit.ndom >= 1, "expected ≥1 domain");
    }
}
