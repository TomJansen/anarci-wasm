use crate::{ViterbiHit, HMM_DATA};
use rustyhmmer::api::HmmDomain;
use rustyhmmer::hmmfile::P7Hmm;
use rustyhmmer::pipeline::Model;
use rustyhmmer::pli::Pipeline;
use rustyhmmer::seqio::Seq;
use rustyhmmer::tblout::fmt_g;
use std::sync::OnceLock;

static RUSTYHMMER_MODELS: OnceLock<Vec<Model>> = OnceLock::new();

#[derive(Clone, Debug)]
struct RustyAlignedHit {
    hmm_name: String,
    bit_score: f64,
    evalue: f64,
    evalue_text: String,
    seq_start: usize,
    seq_end: usize,
    hmm_start: usize,
    hmm_end: usize,
    model_alignment: String,
    query_alignment: String,
}

/// Use rustyhmmer for HMMER-compatible scoring and convert its structured
/// domain alignments into the state path ANARCI numbering consumes.
pub(crate) fn search_hits(
    sequence: &[u8],
    bit_score_threshold: f64,
) -> Result<Vec<ViterbiHit>, String> {
    let sequence =
        std::str::from_utf8(sequence).map_err(|err| format!("sequence is not UTF-8: {err}"))?;
    let models = get_models();
    let z = models.len() as f64;
    let seq = Seq::from_amino("query", "", sequence);
    search_sequence(models, &seq, z, bit_score_threshold)
        .map(|hits| hits.iter().map(rusty_hit_to_viterbi).collect())
}

pub(crate) fn num_models() -> usize {
    get_models().len()
}

fn get_models() -> &'static [Model] {
    RUSTYHMMER_MODELS.get_or_init(load_models)
}

fn load_models() -> Vec<Model> {
    let mut lines = HMM_DATA.lines().map(|line| Ok(line.to_string())).peekable();
    let mut hmms = Vec::new();
    while let Some(hmm) =
        P7Hmm::read_one(&mut lines).expect("embedded HMM database should parse successfully")
    {
        hmms.push(hmm);
    }
    hmms.into_iter().map(Model::new).collect()
}

fn search_sequence(
    models: &[Model],
    seq: &Seq,
    z: f64,
    bit_score_threshold: f64,
) -> Result<Vec<RustyAlignedHit>, String> {
    let seqs = [Seq::from_amino(&seq.name, &seq.desc, &seq_to_string(seq))];
    let mut hits = Vec::new();

    for model in models {
        rustyhmmer::init();
        let mut pli = Pipeline {
            z,
            ..Pipeline::default()
        };
        let Some(mut hit) = model.search_one(&seqs[0], &pli) else {
            continue;
        };
        if !pli.target_reportable(hit.score, hit.lnp) {
            continue;
        }
        hit.seq_idx = 0;

        pli.domz = 1.0;
        for parsed in hit.domains.iter().filter_map(|domain| {
            if !pli.domain_reportable(domain.bitscore, domain.lnp) {
                return None;
            }
            let (model_alignment, query_alignment) =
                rustyhmmer_domain_alignment(domain, model, &seqs[0]);
            Some(RustyAlignedHit {
                hmm_name: model.hmm.name.clone(),
                bit_score: domain.bitscore as f64,
                evalue: domain.lnp.exp() * z,
                evalue_text: fmt_g(domain.lnp.exp() * z, 2),
                seq_start: domain.ali_from as usize - 1,
                seq_end: domain.ali_to as usize,
                hmm_start: domain.hmm_from as usize,
                hmm_end: domain.hmm_to as usize,
                model_alignment,
                query_alignment,
            })
        }) {
            if parsed.bit_score >= bit_score_threshold {
                hits.push(parsed);
            }
        }
    }

    hits.sort_by(|a, b| b.bit_score.total_cmp(&a.bit_score));
    Ok(hits)
}

fn rustyhmmer_domain_alignment(
    domain: &rustyhmmer::pipeline::DomHit,
    model: &Model,
    seq: &Seq,
) -> (String, String) {
    let api_domain = HmmDomain::from_domain(domain, model, seq, 1.0, 1.0);
    (api_domain.model_alignment, api_domain.query_alignment)
}

fn rusty_hit_to_viterbi(hit: &RustyAlignedHit) -> ViterbiHit {
    ViterbiHit {
        hmm_name: hit.hmm_name.clone(),
        bit_score: hit.bit_score,
        evalue: hit.evalue,
        evalue_text: hit.evalue_text.clone(),
        seq_start: hit.seq_start,
        seq_end: hit.seq_end,
        hmm_start: hit.hmm_start,
        hmm_end: hit.hmm_end,
        path: alignment_to_path(hit),
    }
}

fn alignment_to_path(hit: &RustyAlignedHit) -> Vec<(usize, char)> {
    let hmm_states: Vec<usize> = (hit.hmm_start..=hit.hmm_end).collect();
    let mut h = 0usize;
    let mut path = Vec::new();

    for (model, query) in hit.model_alignment.bytes().zip(hit.query_alignment.bytes()) {
        let state_id = hmm_states
            .get(h)
            .copied()
            .unwrap_or_else(|| *hmm_states.last().unwrap_or(&hit.hmm_end));

        if query == b'-' {
            path.push((state_id, 'd'));
            h += 1;
        } else if model == b'.' {
            path.push((state_id, 'i'));
        } else {
            path.push((state_id, 'm'));
            h += 1;
        }
    }

    path
}

fn seq_to_string(seq: &Seq) -> String {
    seq.dsq[1..=seq.len()]
        .iter()
        .map(|&code| rustyhmmer::alphabet::SYM[code as usize] as char)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rustyhmmer_alignment_columns_convert_to_match_insert_delete_path() {
        let hit = RustyAlignedHit {
            hmm_name: "human_H".to_string(),
            bit_score: 100.0,
            evalue: 1e-20,
            evalue_text: "1e-20".to_string(),
            seq_start: 0,
            seq_end: 4,
            hmm_start: 1,
            hmm_end: 4,
            model_alignment: "AB.CD".to_string(),
            query_alignment: "A-CeD".to_string(),
        };

        assert_eq!(
            alignment_to_path(&hit),
            vec![(1, 'm'), (2, 'd'), (3, 'i'), (3, 'm'), (4, 'm')]
        );
    }
}
