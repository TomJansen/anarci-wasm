







use crate::bias::BiasFilter;
use crate::dp::{build_local_profile, p7_domaindef_local, p7_flogsum};
use crate::forward::ForwardFilter;
use crate::hmmfile::P7Hmm;
use crate::msv::{null_one, MsvProfile};
use crate::pli::{Pipeline, Stats};
use crate::seqio::Seq;
use crate::vitfilter::VitFilter;

const LOG2: f64 = std::f64::consts::LN_2;
const LOG2R: f64 = 1.0 / std::f64::consts::LN_2; 
const OMEGA: f64 = 1.0 / 256.0; 


#[inline]
fn exp_logsurv(x: f64, mu: f64, lambda: f64) -> f64 {
    if x < mu {
        0.0
    } else {
        -lambda * (x - mu)
    }
}


pub struct DomHit {
    pub ienv: i64,
    pub jenv: i64,
    pub bitscore: f32,
    pub dombias: f32,
    pub lnp: f64,
    /// Profile coordinates of this domain's optimal-accuracy alignment
    /// (`--domtblout` columns 16/17, `hmm from`/`hmm to`).
    pub hmm_from: i32,
    pub hmm_to: i32,
    /// Sequence coordinates of the alignment (`--domtblout` columns 18/19,
    /// `ali from`/`ali to`).
    pub ali_from: i64,
    pub ali_to: i64,
    /// `dom->oasc`; the `acc` column is `oasc / (1 + |jenv - ienv|)`.
    pub oasc: f32,
    /// Optimal-accuracy traceback in original-dsq coordinates, for
    /// `alidisplay::AliDisplay::create()`.
    pub trace: Vec<crate::optacc::TracePos>,
}


pub struct Hit {
    pub name: String,
    /// `hit->acc`, copied from the target's accession when it has one
    /// (p7_pipeline.c:843: `if (sq->acc[0] != '\0') esl_strdup(sq->acc, -1, &(hit->acc))`).
    pub acc: String,
    pub desc: String,
    /// Target sequence length (`sq->n`, reported as `ad->L`); the `tlen` column of
    /// `--domtblout`.
    pub tlen: usize,
    /// Index of the target in the searched sequence set. C keeps the `ESL_SQ` on the
    /// hit; we keep the index, because target names are not unique in practice (a
    /// concatenation of several proteomes repeats ORF names) and looking the
    /// sequence back up by name can return the wrong one.
    pub seq_idx: usize,
    pub score: f32,     
    pub pre_score: f32, 
    pub lnp: f64,
    pub nexpected: f32,
    pub nregions: i32,
    pub nclustered: i32,
    pub noverlaps: i32,
    pub nenvelopes: i32,
    pub ndom: i32,
    pub best: usize,
    pub domains: Vec<DomHit>,
}

impl Hit {
    
    pub fn bias(&self) -> f32 {
        self.pre_score - self.score
    }
}


pub struct Model {
    pub hmm: P7Hmm,
    pub msv: MsvProfile,
    pub vit: VitFilter,
    pub bias: Option<BiasFilter>,
    pub fwd: ForwardFilter,
}

impl Model {
    pub fn new(hmm: P7Hmm) -> Self {
        let msv = MsvProfile::build(&hmm);
        let vit = VitFilter::build(&hmm);
        let bias = BiasFilter::build(&hmm);
        let fwd = ForwardFilter::build(&hmm);
        Model {
            hmm,
            msv,
            vit,
            bias,
            fwd,
        }
    }

    
    
    pub fn search_one(&self, seq: &Seq, pli: &Pipeline) -> Option<Hit> {
        self.search_one_counted(seq, pli).0
    }

    /// As `search_one`, and also reports how far this target got through the
    /// acceleration pipeline, for `p7_pli_Statistics()`.
    pub fn search_one_counted(&self, seq: &Seq, pli: &Pipeline) -> (Option<Hit>, Stats) {
        let mut st = Stats { nseqs: 1, nres: seq.len() as u64, ..Default::default() };
        let hit = self.search_one_inner(seq, pli, &mut st);
        (hit, st)
    }

    fn search_one_inner(&self, seq: &Seq, pli: &Pipeline, st: &mut Stats) -> Option<Hit> {
        let dsq = &seq.dsq;
        let l = seq.len();
        if l == 0 {
            return None;
        }
        let nullsc = null_one(l) as f64;

        
        
        
        
        let usc = self
            .msv
            .ssv_score(dsq, l)
            .unwrap_or_else(|| self.msv.msv_score(dsq, l));
        if self.msv.pvalue(usc, nullsc as f32) > pli.f1 {
            return None;
        }
        st.n_past_msv += 1;
        
        let filtersc = if pli.do_biasfilter {
            let fs = self.bias.as_ref().unwrap().filter_score(dsq, l);
            if self.msv.pvalue(usc, fs) > pli.f1 {
                return None;
            }
            fs
        } else {
            nullsc as f32
        };
        st.n_past_bias += 1;
        
        
        
        let msv_p = self.msv.pvalue(usc, filtersc);
        if msv_p > pli.f2 {
            let vfsc = self.vit.vit_score(dsq, l);
            if self.vit.pvalue(vfsc, filtersc, l) > pli.f2 {
                return None;
            }
        }
        st.n_past_vit += 1;
        
        let fwdsc = self.fwd.score(dsq, l);
        if self.fwd.pvalue(fwdsc, filtersc) > pli.f3 {
            return None;
        }
        st.n_past_fwd += 1;

        
        
        
        
        
        
        let mut gm = build_local_profile(&self.hmm, l);
        let dd = p7_domaindef_local(&self.fwd, &mut gm, dsq, l, pli.do_null2, pli.rng_seed);
        let domains = &dd.domains;
        if domains.is_empty() {
            return None;
        }

        
        
        
        // C, p7_pipeline.c:700-708: `usc`, `vfsc`, `fwdsc`, `filtersc`, `nullsc`,
        // `seqbias`, `seq_score`, `sum_score`, `pre_score` and `pre2_score` are all
        // declared `float`. Only `eslCONST_LOG2` is a double, so each score is an f32
        // expression divided in f64 and rounded straight back to f32. Carrying these in
        // f64 accumulates differently and shifts the last printed digit.
        //
        // `log((float) sq->n / (float) (sq->n+3))` (p7_pipeline.c:820, :881) is likewise
        // an *f32* division widened only for `log`.
        let ln_ratio = (((l as f32) / ((l + 3) as f32)) as f64).ln();
        let nullsc = nullsc as f32;

        
        
        
        
        
        
        
        
        
        // C, p7_pipeline.c:779-785:
        //   seqbias = esl_vec_FSum(pli->ddef->n2sc, sq->n+1);
        //   seqbias = p7_FLogsum(0.0, log(bg->omega) + seqbias);
        //   ...
        //   pre_score =  (fwdsc - nullsc) / eslCONST_LOG2;
        //   seq_score =  (fwdsc - (nullsc + seqbias)) / eslCONST_LOG2;
        let seqbias_fwd = if pli.do_null2 {
            let n2sc_sum = crate::dp::esl_vec_fsum(&dd.n2sc);
            p7_flogsum(0.0, (OMEGA.ln() + n2sc_sum as f64) as f32)
        } else {
            0.0
        };

        let mut pre_score = ((fwdsc - nullsc) as f64 / LOG2) as f32;
        let mut seq_score = ((fwdsc - (nullsc + seqbias_fwd)) as f64 / LOG2) as f32;

        
        
        // C, p7_pipeline.c:792-822. `sum_score` and `seqbias` accumulate in f32.
        let mut sum_nats = 0.0f32;
        let mut ld = 0i64;
        let mut dombias_sum = 0.0f32;
        for d in domains {
            let significant = if pli.do_null2 {
                d.envsc - d.domcorrection > 0.0
            } else {
                d.envsc > 0.0
            };
            if significant {
                sum_nats += d.envsc;
                ld += d.jenv - d.ienv + 1;
                dombias_sum += d.domcorrection;
            }
        }
        let seqbias_sum = if pli.do_null2 {
            p7_flogsum(0.0, (OMEGA.ln() + dombias_sum as f64) as f32)
        } else {
            0.0
        };

        // `sum_score += (sq->n-Ld) * log(...)`: an f32 accumulator taking an f64 term,
        // so the sum rounds back to f32 here.
        sum_nats = (sum_nats as f64 + (l as i64 - ld) as f64 * ln_ratio) as f32;
        let pre2_score = ((sum_nats - nullsc) as f64 / LOG2) as f32;
        let sum_score = ((sum_nats - (nullsc + seqbias_sum)) as f64 / LOG2) as f32;
        if ld > 0 && sum_score > seq_score {
            seq_score = sum_score;
            pre_score = pre2_score;
        }

        let lnp = exp_logsurv(seq_score as f64, self.fwd.ftau, self.fwd.flambda);

        
        let mut dom_hits = Vec::with_capacity(domains.len());
        let mut best = 0usize;
        let mut best_bits = f32::NEG_INFINITY;
        for (di, d) in domains.iter().enumerate() {
            // C, p7_pipeline.c:880-884:
            //   Ld = hit->dcl[d].jenv - hit->dcl[d].ienv + 1;
            //   hit->dcl[d].bitscore = hit->dcl[d].envsc + (sq->n-Ld) * log((float) sq->n / (float) (sq->n+3)); /* NATS, for the moment... */
            //   hit->dcl[d].dombias  = (pli->do_null2 ? p7_FLogsum(0.0, log(bg->omega) + hit->dcl[d].domcorrection) : 0.0); /* NATS, and will stay so */
            //   hit->dcl[d].bitscore = (hit->dcl[d].bitscore - (nullsc + hit->dcl[d].dombias)) / eslCONST_LOG2; /* now BITS, as it should be */
            //
            // `bitscore` is a `float` field, so the nats value is rounded to f32 before
            // the bias is subtracted, and the subtraction itself is f32.
            let ld_d = d.jenv - d.ienv + 1;
            let bit_nats = (d.envsc as f64 + (l as i64 - ld_d) as f64 * ln_ratio) as f32;
            let dombias = if pli.do_null2 {
                p7_flogsum(0.0, (OMEGA.ln() + d.domcorrection as f64) as f32)
            } else {
                0.0
            };
            let bitscore = ((bit_nats - (nullsc + dombias)) as f64 / LOG2) as f32;
            let dlnp = exp_logsurv(
                bitscore as f64,
                self.fwd.ftau,
                self.fwd.flambda,
            );
            if bitscore > best_bits {
                best_bits = bitscore;
                best = di;
            }
            dom_hits.push(DomHit {
                ienv: d.ienv,
                jenv: d.jenv,
                bitscore,
                dombias,
                lnp: dlnp,
                hmm_from: d.hmm_from,
                hmm_to: d.hmm_to,
                ali_from: d.ali_from,
                ali_to: d.ali_to,
                oasc: d.oasc,
                trace: d.trace.clone(),
            });
        }

        Some(Hit {
            name: seq.name.clone(),
            acc: seq.acc.clone(),
            desc: seq.desc.clone(),
            tlen: l,
            seq_idx: usize::MAX,
            score: seq_score,
            pre_score,
            lnp,
            nexpected: dd.nexpected,
            nregions: dd.nregions,
            nclustered: dd.nclustered,
            noverlaps: dd.noverlaps,
            nenvelopes: dd.nenvelopes,
            ndom: domains.len() as i32,
            best,
            domains: dom_hits,
        })
    }
}


pub fn dombias_bits(dombias_nats: f32) -> f32 {
    (dombias_nats as f64 * LOG2R) as f32
}

/// The `acc` column of `--domtblout`: mean posterior probability of aligned residues.
///
/// `p7_tophits.c:1784` —
/// `(th->hit[h]->dcl[d].oasc / (1.0 + fabs((float) (dcl[d].jenv - dcl[d].ienv))))`.
/// The `(float)` cast is on the envelope length only; the division itself is done in
/// `double` because the literal `1.0` promotes it, so this is reproduced in `f64`.
pub fn dom_acc(d: &DomHit) -> f64 {
    d.oasc as f64 / (1.0 + ((d.jenv - d.ienv) as f32).abs() as f64)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seqio::read_fasta;

    #[test]
    fn globins_scores_match_golden() {
        let hmm = P7Hmm::read_all(&format!("{}/testdata/globins4.hmm", env!("CARGO_MANIFEST_DIR")))
            .unwrap()
            .pop()
            .unwrap();
        let model = Model::new(hmm);
        let seqs = read_fasta(&format!("{}/testdata/globins45.fa", env!("CARGO_MANIFEST_DIR")))
            .unwrap();
        let mut pli = Pipeline::default();
        pli.z = 45.0;
        pli.domz = pli.z;
        let z = 45.0_f64;

        
        let mut found = 0;
        for s in &seqs {
            if let Some(h) = model.search_one(s, &pli) {
                found += 1;
                if h.name == "MYG_ESCGI" {
                    let e = h.lnp.exp() * z;
                    eprintln!("MYG_ESCGI: score={:.1} bias={:.1} E={:.2e}", h.score, h.bias(), e);
                    assert!((h.score - 215.6).abs() < 0.15, "score {:.2} != 215.6", h.score);
                    assert!((h.bias() - 2.9).abs() < 0.15, "bias {:.2} != 2.9", h.bias());
                    
                    assert!((e / 8.7e-67).ln().abs() < 0.05, "E {:.3e} != 8.7e-67", e);
                }
            }
        }
        assert_eq!(found, 45, "expected 45 hits, got {found}");
    }
}
