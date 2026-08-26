//! Pipeline configuration: reporting/inclusion thresholds, search-space sizes,
//! acceleration-filter settings and RNG seeding.
//!
//! Faithful transcription of the corresponding parts of HMMER 3.4's `P7_PIPELINE`:
//!
//! | Rust                        | C (hmmer-3.4)                                  |
//! |-----------------------------|------------------------------------------------|
//! | `Pipeline::from_opts`       | `p7_pipeline_Create()`   p7_pipeline.c:108-235  |
//! | `Pipeline::new_model_thresholds` | `p7_pli_NewModelThresholds()` p7_pipeline.c:341-372 |
//! | `target_reportable`         | `p7_pli_TargetReportable()`   p7_pipeline.c:407-417 |
//! | `domain_reportable`         | `p7_pli_DomainReportable()`   p7_pipeline.c:427-436 |
//! | `target_includable`         | `p7_pli_TargetIncludable()`   p7_pipeline.c:446-457 |
//! | `domain_includable`         | `p7_pli_DomainIncludable()`   p7_pipeline.c:466-472 |
//!
//! Only the `long_targets == FALSE` (protein `hmmsearch`) branches are transcribed;
//! the nhmmer branches are marked where they are elided.

/// Which of the model's own bit-score cutoff pairs to threshold on
/// (`p7H_GA` / `p7H_TC` / `p7H_NC`, hmmer.h).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BitCutoffs {
    Ga,
    Tc,
    Nc,
}

/// `p7_ZSETBY_NTARGETS` / `p7_ZSETBY_OPTION` (hmmer.h).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ZSetBy {
    NTargets,
    Option,
}

#[derive(Clone, Debug)]
pub struct Pipeline {
    /* --- reporting thresholds (p7_pipeline.c:134-152) --- */
    pub by_e: bool,
    pub e: f64,
    pub t: f64,
    pub dom_by_e: bool,
    pub dom_e: f64,
    pub dom_t: f64,
    pub use_bit_cutoffs: Option<BitCutoffs>,

    /* --- inclusion thresholds (p7_pipeline.c:155-173) --- */
    pub inc_by_e: bool,
    pub inc_e: f64,
    pub inc_t: f64,
    pub incdom_by_e: bool,
    pub incdom_e: f64,
    pub incdom_t: f64,

    /* --- search space sizes (p7_pipeline.c:201-213) --- */
    pub z: f64,
    pub domz: f64,
    pub z_setby: ZSetBy,
    pub domz_setby: ZSetBy,

    /* --- acceleration pipeline (p7_pipeline.c:216-234) --- */
    pub do_max: bool,
    pub do_biasfilter: bool,
    pub do_null2: bool,
    pub f1: f64,
    pub f2: f64,
    pub f3: f64,

    /* --- RNG (p7_pipeline.c:124-130) --- */
    /// The seed as given; 0 means "one-time arbitrary seed".
    pub seed: u32,
    /// Effective seed actually handed to the stochastic-trace RNG. Equals `seed`
    /// unless `seed == 0`, in which case `choose_arbitrary_seed()` picked it.
    pub rng_seed: u32,
    pub do_reseeding: bool,

    /* --- display (p7_pipeline.c:204-205) --- */
    pub show_accessions: bool,
    pub show_alignments: bool,
}

impl Default for Pipeline {
    /// The defaults `p7_pipeline_Create()` installs when `go == NULL`, i.e. an
    /// `hmmsearch` run with no threshold options given.
    fn default() -> Self {
        Pipeline {
            by_e: true,
            e: 10.0,
            t: 0.0,
            dom_by_e: true,
            dom_e: 10.0,
            dom_t: 0.0,
            use_bit_cutoffs: None,

            inc_by_e: true,
            inc_e: 0.01,
            inc_t: 0.0,
            incdom_by_e: true,
            incdom_e: 0.01,
            incdom_t: 0.0,

            z: 0.0,
            domz: 0.0,
            z_setby: ZSetBy::NTargets,
            domz_setby: ZSetBy::NTargets,

            do_max: false,
            do_biasfilter: true,
            do_null2: true,
            f1: 0.02,
            f2: 1e-3,
            f3: 1e-5,

            seed: 42,
            rng_seed: 42,
            do_reseeding: true,

            show_accessions: false,
            show_alignments: true,
        }
    }
}

impl Pipeline {
    /// `--cut_ga` / `--cut_tc` / `--cut_nc` (p7_pipeline.c:176-199). Each zeroes both
    /// reporting and inclusion bit thresholds and switches all four off E-values; the
    /// actual values arrive per model in `new_model_thresholds`.
    // C, p7_pipeline.c:176-199 (one identical block per cutoff kind):
    //   if (go && esl_opt_GetBoolean(go, "--cut_ga"))
    //     {
    //       pli->T        = pli->domT        = 0.0;
    //       pli->by_E     = pli->dom_by_E    = FALSE;
    //       pli->incT     = pli->incdomT     = 0.0;
    //       pli->inc_by_E = pli->incdom_by_E = FALSE;
    //       pli->use_bit_cutoffs = p7H_GA;
    //     }
    pub fn with_bit_cutoffs(mut self, which: BitCutoffs) -> Self {
        self.t = 0.0;
        self.dom_t = 0.0;
        self.by_e = false;
        self.dom_by_e = false;
        self.inc_t = 0.0;
        self.incdom_t = 0.0;
        self.inc_by_e = false;
        self.incdom_by_e = false;
        self.use_bit_cutoffs = Some(which);
        self
    }

    /// `--max` (p7_pipeline.c:225-231): turn off the bias filter and open F1/F2/F3.
    /// Only the `long_targets == FALSE` value of F1 is used here.
    // C, p7_pipeline.c:225-231:
    //   if (go && esl_opt_GetBoolean(go, "--max"))
    //     {
    //       pli->do_max        = TRUE;
    //       pli->do_biasfilter = FALSE;
    //       pli->F2 = pli->F3 = 1.0;
    //       pli->F1 = (pli->long_targets ? 0.3 : 1.0);
    //     }
    pub fn with_max(mut self) -> Self {
        self.do_max = true;
        self.do_biasfilter = false;
        self.f2 = 1.0;
        self.f3 = 1.0;
        self.f1 = 1.0;
        self
    }

    /// `--seed <n>` (p7_pipeline.c:124-130). Seed 0 means "choose an arbitrary
    /// one-time seed and stop reseeding", so runs stop being reproducible.
    pub fn with_seed(mut self, seed: u32) -> Self {
        self.seed = seed;
        self.do_reseeding = seed != 0;
        self.rng_seed = if seed == 0 { choose_arbitrary_seed() } else { seed };
        self
    }

    /// `p7_pli_NewModelThresholds()` (p7_pipeline.c:341-372): for `--cut_ga/tc/nc`,
    /// install this model's own cutoff pair as both the reporting and the inclusion
    /// threshold. `cutoffs` is `(seq, domain)`; `None` means the model lacks them,
    /// which C reports as `eslEINVAL`.
    // C, p7_pipeline.c:345-368:
    //   if (pli->use_bit_cutoffs)
    //   {
    //     if (pli->use_bit_cutoffs == p7H_GA)
    //     {
    //       if (om->cutoff[p7_GA1] == p7_CUTOFF_UNSET)
    //         ESL_FAIL(eslEINVAL, pli->errbuf, "GA bit thresholds unavailable on model %s\n", om->name);
    //       pli->T    = pli->incT    = om->cutoff[p7_GA1];
    //       pli->domT = pli->incdomT = om->cutoff[p7_GA2];
    //     }
    //     else if  (pli->use_bit_cutoffs == p7H_TC) { ...p7_TC1/p7_TC2... }
    //     else if  (pli->use_bit_cutoffs == p7H_NC) { ...p7_NC1/p7_NC2... }
    //   }
    pub fn new_model_thresholds(&mut self, cutoffs: Option<(f64, f64)>) -> Result<(), ()> {
        if self.use_bit_cutoffs.is_none() {
            return Ok(());
        }
        match cutoffs {
            Some((c1, c2)) => {
                self.t = c1;
                self.inc_t = c1;
                self.dom_t = c2;
                self.incdom_t = c2;
                Ok(())
            }
            None => Err(()),
        }
    }

    /// `p7_pli_TargetReportable()` — p7_pipeline.c:407-417.
    // C, p7_pipeline.c:407-417:
    //   int
    //   p7_pli_TargetReportable(P7_PIPELINE *pli, float score, double lnP)
    //   {
    //     if      (  pli->by_E )
    //       {
    //         if ( !pli->long_targets  && exp(lnP) * pli->Z <= pli->E) return TRUE;
    //         if (  pli->long_targets  && exp(lnP) <= pli->E)          return TRUE;
    //       }
    //     else if (! pli->by_E   && score         >= pli->T) return TRUE;
    //     return FALSE;
    //   }
    pub fn target_reportable(&self, score: f32, lnp: f64) -> bool {
        if self.by_e {
            lnp.exp() * self.z <= self.e
        } else {
            (score as f64) >= self.t
        }
    }

    /// `p7_pli_DomainReportable()` — p7_pipeline.c:427-436.
    // C, p7_pipeline.c:427-436:
    //   if      (  pli->dom_by_E )
    //     {
    //       if ( !pli->long_targets  &&  exp(lnP) * pli->domZ <= pli->domE) return TRUE;
    //       if (  pli->long_targets  &&  exp(lnP) <= pli->domE) return TRUE;
    //     }
    //   else if (! pli->dom_by_E   && dom_score        >= pli->domT) return TRUE;
    //   return FALSE;
    pub fn domain_reportable(&self, dom_score: f32, lnp: f64) -> bool {
        if self.dom_by_e {
            lnp.exp() * self.domz <= self.dom_e
        } else {
            (dom_score as f64) >= self.dom_t
        }
    }

    /// `p7_pli_TargetIncludable()` — p7_pipeline.c:446-457.
    // C, p7_pipeline.c:446-457:
    //   if      (  pli->inc_by_E )
    //     {
    //       if ( !pli->long_targets && exp(lnP) * pli->Z <= pli->incE) return TRUE;
    //       if (  pli->long_targets && exp(lnP) <= pli->incE) return TRUE;
    //     }
    //   else if (! pli->inc_by_E   && score         >= pli->incT) return TRUE;
    //   return FALSE;
    pub fn target_includable(&self, score: f32, lnp: f64) -> bool {
        if self.inc_by_e {
            lnp.exp() * self.z <= self.inc_e
        } else {
            (score as f64) >= self.inc_t
        }
    }

    /// `p7_pli_DomainIncludable()` — p7_pipeline.c:466-472.
    // C, p7_pipeline.c:466-472:
    //   if      (  pli->incdom_by_E   && exp(lnP) * pli->domZ <= pli->incdomE) return TRUE;
    //   else if (! pli->incdom_by_E   && dom_score        >= pli->incdomT) return TRUE;
    //   else return FALSE;
    pub fn domain_includable(&self, dom_score: f32, lnp: f64) -> bool {
        if self.incdom_by_e {
            lnp.exp() * self.domz <= self.incdom_e
        } else {
            (dom_score as f64) >= self.incdom_t
        }
    }
}

/// The accounting `p7_pipeline_Create()` zeroes at p7_pipeline.c:236-249 and
/// `p7_pli_Statistics()` reports (p7_pipeline.c:1096-1160).
#[derive(Clone, Copy, Debug, Default)]
pub struct Stats {
    pub nmodels: u64,
    pub nnodes: u64,
    pub nseqs: u64,
    pub nres: u64,
    pub n_past_msv: u64,
    pub n_past_bias: u64,
    pub n_past_vit: u64,
    pub n_past_fwd: u64,
}

impl Stats {
    /// `p7_pipeline_Merge()` — p7_pipeline.c:544-573, minus the long-target fields.
    // C, p7_pipeline.c:559-573 (the non-long_targets counters):
    //   p1->nmodels     += p2->nmodels;
    //   p1->nseqs       += p2->nseqs;
    //   p1->nres        += p2->nres;
    //   p1->nnodes      += p2->nnodes;
    //   p1->n_past_msv  += p2->n_past_msv;
    //   p1->n_past_bias += p2->n_past_bias;
    //   p1->n_past_vit  += p2->n_past_vit;
    //   p1->n_past_fwd  += p2->n_past_fwd;
    pub fn merge(&mut self, o: &Stats) {
        self.nmodels += o.nmodels;
        self.nnodes += o.nnodes;
        self.nseqs += o.nseqs;
        self.nres += o.nres;
        self.n_past_msv += o.n_past_msv;
        self.n_past_bias += o.n_past_bias;
        self.n_past_vit += o.n_past_vit;
        self.n_past_fwd += o.n_past_fwd;
    }
}

/// `choose_arbitrary_seed()` — easel/esl_random.c. C mixes `time(NULL)`, `getpid()`
/// and `clock()`. std has no `getpid`/`clock`, so the process-identity and
/// CPU-time words come from the process's start-of-run instants instead; the
/// intent — decorrelating closely spaced invocations — and the `0 -> 42` guard are
/// preserved. Only reachable via `--seed 0`, which by definition asks for a
/// non-reproducible run.
fn choose_arbitrary_seed() -> u32 {
    use std::time::{SystemTime, UNIX_EPOCH};
    let now = SystemTime::now().duration_since(UNIX_EPOCH).unwrap_or_default();
    let a = now.as_secs() as u32;
    let b = now.subsec_nanos();
    let c = std::time::Instant::now().elapsed().subsec_nanos();
    let seed = crate::dp::esl_mix3(a, b, c);
    if seed == 0 {
        42
    } else {
        seed
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The `go == NULL` defaults of `p7_pipeline_Create()`.
    #[test]
    fn defaults_match_c() {
        let p = Pipeline::default();
        assert!(p.by_e && p.dom_by_e && p.inc_by_e && p.incdom_by_e);
        assert_eq!((p.e, p.dom_e, p.inc_e, p.incdom_e), (10.0, 10.0, 0.01, 0.01));
        assert_eq!((p.f1, p.f2, p.f3), (0.02, 1e-3, 1e-5));
        assert!(p.do_biasfilter && p.do_null2 && !p.do_max);
        assert_eq!(p.seed, 42);
        assert!(p.do_reseeding);
    }

    /// `--cut_ga` must switch reporting *and* inclusion onto bit scores, and
    /// `new_model_thresholds` must install GA1 for sequences and GA2 for domains.
    #[test]
    fn cut_ga_installs_both_thresholds() {
        let mut p = Pipeline::default().with_bit_cutoffs(BitCutoffs::Ga);
        assert!(!p.by_e && !p.dom_by_e && !p.inc_by_e && !p.incdom_by_e);
        p.new_model_thresholds(Some((116.25, 90.5))).unwrap();
        assert_eq!((p.t, p.inc_t), (116.25, 116.25));
        assert_eq!((p.dom_t, p.incdom_t), (90.5, 90.5));
        assert!(p.target_reportable(120.0, 0.0));
        assert!(!p.target_reportable(100.0, 0.0));
        assert!(p.domain_reportable(95.0, 0.0));
        assert!(!p.domain_reportable(80.0, 0.0));
        // A model without the cutoffs is an error, as in C.
        assert!(Pipeline::default()
            .with_bit_cutoffs(BitCutoffs::Ga)
            .new_model_thresholds(None)
            .is_err());
    }

    /// `--max` opens the filters and disables the bias filter.
    #[test]
    fn max_opens_filters() {
        let p = Pipeline::default().with_max();
        assert_eq!((p.f1, p.f2, p.f3), (1.0, 1.0, 1.0));
        assert!(p.do_max && !p.do_biasfilter);
    }

    /// Seed 0 stops reseeding and picks a nonzero arbitrary seed.
    #[test]
    fn seed_zero_is_arbitrary_and_stops_reseeding() {
        let p = Pipeline::default().with_seed(0);
        assert!(!p.do_reseeding);
        assert_ne!(p.rng_seed, 0);
        let q = Pipeline::default().with_seed(7);
        assert!(q.do_reseeding);
        assert_eq!(q.rng_seed, 7);
    }
}
