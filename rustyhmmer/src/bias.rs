













use crate::alphabet::{AMINO_FREQ, K, KP};
use crate::hmmfile::P7Hmm;
use crate::msv::null_one;



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

pub struct BiasFilter {
    
    eo: Vec<[f32; 2]>,
    
    t1: [f32; 3],
    
    pi: [f32; 2],
}

impl BiasFilter {
    
    
    pub fn build(hmm: &P7Hmm) -> Option<Self> {
        let compo = hmm.compo?; 
        let l1 = hmm.m as f32 / 8.0; 

        
        
        let mut eo = vec![[1.0f32; 2]; KP];
        for x in 0..K {
            
            eo[x][0] = 1.0;
            eo[x][1] = compo[x] / AMINO_FREQ[x] as f32;
        }
        
        eo[K] = [1.0, 1.0];
        eo[KP - 2] = [1.0, 1.0];
        eo[KP - 1] = [1.0, 1.0];
        
        for code in (K + 1)..=(KP - 3) {
            let set = degen_set(code);
            
            let (mut n0, mut n1, mut den) = (0.0f32, 0.0f32, 0.0f32);
            for &y in set {
                n0 += AMINO_FREQ[y] as f32; 
                n1 += compo[y]; 
                den += AMINO_FREQ[y] as f32; 
            }
            eo[code][0] = if den > 0.0 { n0 / den } else { 0.0 };
            eo[code][1] = if den > 0.0 { n1 / den } else { 0.0 };
        }

        
        let t1 = [1.0 / (l1 + 1.0), l1 / (l1 + 1.0), 1.0];

        Some(BiasFilter { eo, t1, pi: [0.999, 0.001] })
    }

    
    
    pub fn filter_score(&self, dsq: &[u8], l: usize) -> f32 {
        if l == 0 {
            return 0.0;
        }
        
        let p1 = l as f32 / (l as f32 + 1.0);
        let t: [[f32; 3]; 2] = [[p1, 1.0 - p1, 1.0], self.t1];

        
        let x1 = dsq[1] as usize;
        let mut dp = [self.eo[x1][0] * self.pi[0], self.eo[x1][1] * self.pi[1]];
        let mut max = dp[0].max(dp[1]);
        dp[0] /= max;
        dp[1] /= max;
        let mut logsc = (max as f64).ln() as f32;

        
        for i in 2..=l {
            let x = dsq[i] as usize;
            let mut nd = [0.0f32; 2];
            for k in 0..2 {
                let mut s = 0.0f32;
                for m in 0..2 {
                    s += dp[m] * t[m][k];
                }
                nd[k] = s * self.eo[x][k];
            }
            max = nd[0].max(nd[1]);
            dp[0] = nd[0] / max;
            dp[1] = nd[1] / max;
            logsc += (max as f64).ln() as f32;
        }

        
        let end = dp[0] * t[0][2] + dp[1] * t[1][2];
        logsc += (end as f64).ln() as f32;

        
        
        logsc + null_one(l)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::msv::MsvProfile;
    use crate::seqio::read_fasta;

    #[test]
    fn bias_passes_all_golden_hits() {
        
        let hmm = P7Hmm::read_all(&format!("{}/testdata/globins4.hmm", env!("CARGO_MANIFEST_DIR")))
            .unwrap()
            .pop()
            .unwrap();
        let mp = MsvProfile::build(&hmm);
        let bf = BiasFilter::build(&hmm).expect("globins4 has COMPO");
        let seqs = read_fasta(&format!("{}/testdata/globins45.fa", env!("CARGO_MANIFEST_DIR")))
            .unwrap();
        let f1 = 0.02_f64;
        let mut n_pass = 0;
        for s in &seqs {
            let l = s.len();
            let usc = mp.msv_score(&s.dsq, l);
            
            if mp.pvalue(usc, crate::msv::null_one(l)) > f1 {
                continue;
            }
            let filtersc = bf.filter_score(&s.dsq, l);
            if mp.pvalue(usc, filtersc) <= f1 {
                n_pass += 1;
            }
        }
        assert_eq!(n_pass, 45, "bias filter pass count {n_pass} != C's 45");
    }
}
