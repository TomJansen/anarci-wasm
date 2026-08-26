use crate::alidisplay::AliDisplay;
use crate::pipeline::{dombias_bits, Hit, Model};
use crate::pli::{BitCutoffs, Pipeline};
use crate::seqio::Seq;
use rayon::prelude::*;

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Cutoff {
    GatheringGa,
    TrustedTc,
    NoiseNc,
    Evalue(f64),
}

/// One domain of a hit, carrying the same coordinates `--domtblout` reports.
///
/// `hmm_from`/`hmm_to` are the profile coordinates of the domain's optimal-accuracy
/// alignment (`--domtblout` columns 16/17), so HMM coverage is
/// `(hmm_to - hmm_from + 1) / model_len`.
#[derive(Clone, Debug)]
pub struct HmmDomain {
    pub score: f32,
    pub bias: f32,
    pub evalue: f64,
    pub c_evalue: f64,
    pub hmm_from: i32,
    pub hmm_to: i32,
    pub ali_from: i64,
    pub ali_to: i64,
    pub env_from: i64,
    pub env_to: i64,
    /// Mean posterior probability of aligned residues (`--domtblout` `acc` column).
    pub acc: f32,
    /// Consensus/model alignment columns for this domain, as printed in the
    /// human-readable alignment display.
    pub model_alignment: String,
    /// Target/query alignment columns for this domain, as printed in the
    /// human-readable alignment display.
    pub query_alignment: String,
}

impl HmmDomain {
    pub fn from_domain(
        domain: &crate::pipeline::DomHit,
        model: &Model,
        seq: &Seq,
        z: f64,
        domz: f64,
    ) -> Self {
        let (model_alignment, query_alignment) = domain_alignment(domain, model, seq);
        HmmDomain {
            score: domain.bitscore,
            bias: dombias_bits(domain.dombias),
            evalue: domain.lnp.exp() * z,
            c_evalue: domain.lnp.exp() * domz,
            hmm_from: domain.hmm_from,
            hmm_to: domain.hmm_to,
            ali_from: domain.ali_from,
            ali_to: domain.ali_to,
            env_from: domain.ienv,
            env_to: domain.jenv,
            acc: crate::pipeline::dom_acc(domain) as f32,
            model_alignment,
            query_alignment,
        }
    }
}

#[derive(Clone, Debug)]
pub struct HmmHit {
    pub query_name: String,
    pub query_acc: String,
    pub model_len: usize,
    pub target_name: String,
    pub seq_score: f32,
    pub seq_bias: f32,
    pub seq_evalue: f64,
    pub best_dom_score: f32,
    pub best_dom_bias: f32,
    pub best_dom_evalue: f64,
    pub env_from: i64,
    pub env_to: i64,
    /// Profile coordinates of the *best* domain's alignment, alongside the
    /// `env_from`/`env_to` above (which are also the best domain's).
    pub hmm_from: i32,
    pub hmm_to: i32,
    pub ali_from: i64,
    pub ali_to: i64,
    pub n_domains: i32,
    /// Every domain of this hit, in the order the domain-definition step found them.
    ///
    /// Unlike `--domtblout`, this list is **not** filtered by the per-domain
    /// reporting threshold, so it can contain domains C would not print (and it is
    /// therefore not indexed the same way as domtblout's `#` column, which counts
    /// only reported domains). Filter on `c_evalue` if you want domtblout's set:
    /// C reports a domain when `c_evalue <= 10`.
    pub domains: Vec<HmmDomain>,
}

pub struct HmmAnnotator {
    models: Vec<Model>,
    cutoff: Cutoff,
}

impl HmmAnnotator {
    pub fn from_hmm_file(path: &str) -> Result<Self, String> {
        let hmms = crate::hmmfile::P7Hmm::read_all(path).map_err(|e| e.to_string())?;
        if hmms.is_empty() {
            return Err(format!("no models in {path}"));
        }
        Ok(Self {
            models: hmms.into_iter().map(Model::new).collect(),
            cutoff: Cutoff::Evalue(10.0),
        })
    }

    pub fn from_models(models: Vec<Model>) -> Self {
        Self {
            models,
            cutoff: Cutoff::Evalue(10.0),
        }
    }

    pub fn with_cutoff(mut self, cutoff: Cutoff) -> Self {
        self.cutoff = cutoff;
        self
    }

    pub fn len(&self) -> usize {
        self.models.len()
    }
    pub fn is_empty(&self) -> bool {
        self.models.is_empty()
    }

    pub fn search(&self, seqs: &[Seq]) -> Vec<HmmHit> {
        let z = seqs.len() as f64;
        // Parallelize over PROFILES (coarse-grained), not over sequences within
        // each profile. A profile DB (~10^4-10^5 models) offers far more
        // parallelism than a single genome's sequence set, and most profiles are
        // rejected cheaply by the MSV filter — so a per-profile `par_iter` over
        // sequences paid rayon fork/join overhead on ~10^5 near-empty dispatches
        // (measured ~35% core usage on 16 threads). One coarse dispatch over the
        // profile DB — the scheme HMMER/pyhmmer use — keeps the cores saturated.
        // `par_iter().map(..).collect::<Vec<_>>()` is order-preserving
        // (IndexedParallelIterator), so the concatenated output is byte-identical
        // to the previous serial (model-then-sequence) order.
        let per_model: Vec<Vec<HmmHit>> = self
            .models
            .par_iter()
            .map(|model| {
                // MXCSR is per-thread and rayon's pool is not ours to configure;
                // setting it per task is a couple of instructions against a whole
                // model search.
                crate::init();
                // One P7_PIPELINE per query, as C does: `p7_pli_NewModelThresholds()`
                // installs this model's own GA/TC/NC pair (p7_pipeline.c:341-372).
                let mut pli = Pipeline::default();
                pli.z = z;
                let cuts = match self.cutoff {
                    Cutoff::Evalue(e) => {
                        pli.e = e;
                        None
                    }
                    Cutoff::GatheringGa => {
                        pli = pli.with_bit_cutoffs(BitCutoffs::Ga);
                        Some(model.hmm.cutoffs.ga)
                    }
                    Cutoff::TrustedTc => {
                        pli = pli.with_bit_cutoffs(BitCutoffs::Tc);
                        Some(model.hmm.cutoffs.tc)
                    }
                    Cutoff::NoiseNc => {
                        pli = pli.with_bit_cutoffs(BitCutoffs::Nc);
                        Some(model.hmm.cutoffs.nc)
                    }
                };
                let mut local = Vec::new();
                if let Some(c) = cuts {
                    // C fails the run here; the library API skips the model instead.
                    if pli.new_model_thresholds(c).is_err() {
                        return local;
                    }
                }

                let hits: Vec<Hit> = seqs
                    .iter()
                    .enumerate()
                    .filter_map(|(seq_idx, s)| {
                        model.search_one(s, &pli).map(|mut h| {
                            h.seq_idx = seq_idx;
                            h
                        })
                    })
                    .collect();

                let qname = &model.hmm.name;
                let qacc = model.hmm.acc.clone().unwrap_or_default();
                // `pli->domZ` is the number of reported targets for this query
                // (p7_tophits.c:963), and per-domain conditional E-values divide by it.
                let reported: Vec<&Hit> = hits
                    .iter()
                    .filter(|h| pli.target_reportable(h.score, h.lnp))
                    .collect();
                pli.domz = reported.len() as f64;
                let domz = pli.domz;
                for h in reported {
                    let bd = &h.domains[h.best];
                    local.push(HmmHit {
                        query_name: qname.clone(),
                        query_acc: qacc.clone(),
                        model_len: model.hmm.m,
                        target_name: h.name.clone(),
                        seq_score: h.score,
                        seq_bias: h.bias(),
                        seq_evalue: h.lnp.exp() * z,
                        best_dom_score: bd.bitscore,
                        best_dom_bias: dombias_bits(bd.dombias),
                        best_dom_evalue: bd.lnp.exp() * z,
                        env_from: bd.ienv,
                        env_to: bd.jenv,
                        hmm_from: bd.hmm_from,
                        hmm_to: bd.hmm_to,
                        ali_from: bd.ali_from,
                        ali_to: bd.ali_to,
                        n_domains: h.ndom,
                        domains: h
                            .domains
                            .iter()
                            .map(|d| HmmDomain::from_domain(d, model, &seqs[h.seq_idx], z, domz))
                            .collect(),
                    });
                }
                local
            })
            .collect();
        per_model.into_iter().flatten().collect()
    }
}

fn domain_alignment(
    domain: &crate::pipeline::DomHit,
    model: &Model,
    seq: &Seq,
) -> (String, String) {
    AliDisplay::create(&domain.trace, &model.hmm, &model.fwd, seq)
        .map(|ad| {
            (
                String::from_utf8_lossy(&ad.model).into_owned(),
                String::from_utf8_lossy(&ad.aseq).into_owned(),
            )
        })
        .unwrap_or_else(|| (String::new(), String::new()))
}
