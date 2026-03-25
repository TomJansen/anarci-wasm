use genomic_system_finder_hmm::alphabet;
use genomic_system_finder_hmm::hmm::{Plan7Hmm, TDD, TDM, TII, TIM, TMD, TMI, TMM};
use genomic_system_finder_hmm::parser;
use genomic_system_finder_hmm::pipeline::{self, SearchOptions, TargetSequence};
use std::io::Cursor;

const NEG_INF: f64 = f64::NEG_INFINITY;

const TB_START: u8 = 0;
const TB_FROM_M: u8 = 1;
const TB_FROM_I: u8 = 2;
const TB_FROM_D: u8 = 3;

/// A collection of profile HMMs loaded through `genomic_system_finder_hmm`.
pub struct HmmDatabase {
    pub profiles: Vec<Plan7Hmm>,
}

/// Result of an HMM hit adapted to the shape ANARCI's numbering code expects.
#[derive(Clone)]
pub struct ViterbiHit {
    pub hmm_name: String,
    pub bit_score: f64, // score in bits
    pub evalue: f64,
    pub evalue_text: String,
    pub seq_start: usize, // 0-based inclusive
    pub seq_end: usize,   // 0-based exclusive
    pub hmm_start: usize, // 1-based inclusive
    pub hmm_end: usize,   // 1-based inclusive
    /// State path: (hmm_state_1based, state_type) for each aligned position
    /// state_type: 'm' = match, 'i' = insert, 'd' = delete
    pub path: Vec<(usize, char)>,
}

impl HmmDatabase {
    /// Parse an embedded HMM database using `genomic_system_finder_hmm`.
    pub fn parse(text: &str) -> Self {
        let profiles = parser::parse_hmm(Cursor::new(text.as_bytes()))
            .expect("embedded HMM database should parse successfully");
        Self { profiles }
    }

    /// Run all profiles against a sequence, return hits above threshold sorted by score.
    pub fn search(&self, sequence: &[u8], bit_score_threshold: f64) -> Vec<ViterbiHit> {
        let targets = vec![TargetSequence {
            id: "query".to_string(),
            sequence: sequence.to_vec(),
        }];
        let options = SearchOptions {
            e_value: f64::INFINITY,
            i_evalue: f64::INFINITY,
            coverage_profile: 0.0,
            cut_ga: false,
            msv_evalue_threshold: None,
        };

        let mut hits = Vec::new();
        for hmm in &self.profiles {
            let result = pipeline::hmmsearch(hmm, &targets, &options);
            for hit in result.hits {
                if hit.score < bit_score_threshold {
                    continue;
                }
                if let Some(adapted) = adapt_hit(hmm, sequence, &hit) {
                    hits.push(adapted);
                }
            }
        }

        hits.sort_by(|a, b| b.bit_score.total_cmp(&a.bit_score));
        hits
    }
}

fn adapt_hit(hmm: &Plan7Hmm, sequence: &[u8], hit: &pipeline::DomainHit) -> Option<ViterbiHit> {
    let seq_start = hit.ali_from.checked_sub(1)? as usize;
    let seq_end = hit.ali_to as usize;
    let hmm_start = hit.hmm_from as usize;
    let hmm_end = hit.hmm_to as usize;

    if seq_start >= seq_end || seq_end > sequence.len() {
        return None;
    }
    if hmm_start == 0 || hmm_end < hmm_start || hmm_end > hmm.model_length {
        return None;
    }

    let path = traceback_path(hmm, &sequence[seq_start..seq_end], hmm_start, hmm_end)?;

    Some(ViterbiHit {
        hmm_name: hit.query_name.clone(),
        bit_score: hit.score,
        evalue: stable_evalue_from_bits(hit.score, hmm, 1),
        evalue_text: format_evalue_from_bits(hit.score, hmm, 1),
        seq_start,
        seq_end,
        hmm_start,
        hmm_end,
        path,
    })
}

fn traceback_path(
    hmm: &Plan7Hmm,
    sequence: &[u8],
    hmm_start: usize,
    hmm_end: usize,
) -> Option<Vec<(usize, char)>> {
    // `genomic_system_finder_hmm` currently returns domain boundaries and scores,
    // but ANARCI's numbering needs the per-state M/I/D path as well.
    let dsq = alphabet::encode_sequence(sequence);
    let l = dsq.len();
    let m = hmm_end.checked_sub(hmm_start)? + 1;

    if l == 0 || m == 0 {
        return None;
    }

    let stride = m + 1;
    let matrix_len = (l + 1) * stride;

    let mut dp_m = vec![NEG_INF; matrix_len];
    let mut dp_i = vec![NEG_INF; matrix_len];
    let mut dp_d = vec![NEG_INF; matrix_len];

    let mut tb_m = vec![TB_START; matrix_len];
    let mut tb_i = vec![TB_START; matrix_len];
    let mut tb_d = vec![TB_START; matrix_len];

    let idx = |j: usize, k: usize| j * stride + k;
    let abs_k = |k: usize| hmm_start + k - 1;

    dp_m[idx(1, 1)] = hmm.match_score(hmm_start, dsq[0]);

    for j in 1..=l {
        for k in 1..=m {
            if j == 1 && k == 1 {
                continue;
            }

            let cell = idx(j, k);
            let model_k = abs_k(k);

            if j > 1 && k > 1 {
                let prev = idx(j - 1, k - 1);
                let prev_model = abs_k(k - 1);

                let from_m = dp_m[prev] + hmm.transition(prev_model, TMM);
                let from_i = dp_i[prev] + hmm.transition(prev_model, TIM);
                let from_d = dp_d[prev] + hmm.transition(prev_model, TDM);

                let (best_prev, tb) = if from_m >= from_i && from_m >= from_d {
                    (from_m, TB_FROM_M)
                } else if from_i >= from_d {
                    (from_i, TB_FROM_I)
                } else {
                    (from_d, TB_FROM_D)
                };

                if best_prev.is_finite() {
                    dp_m[cell] = best_prev + hmm.match_score(model_k, dsq[j - 1]);
                    tb_m[cell] = tb;
                }
            }

            if j > 1 {
                let prev = idx(j - 1, k);
                let from_m = dp_m[prev] + hmm.transition(model_k, TMI);
                let from_i = dp_i[prev] + hmm.transition(model_k, TII);

                if from_m >= from_i {
                    if from_m.is_finite() {
                        dp_i[cell] = from_m + hmm.insert_score(model_k, dsq[j - 1]);
                        tb_i[cell] = TB_FROM_M;
                    }
                } else if from_i.is_finite() {
                    dp_i[cell] = from_i + hmm.insert_score(model_k, dsq[j - 1]);
                    tb_i[cell] = TB_FROM_I;
                }
            }

            if k > 1 {
                let prev = idx(j, k - 1);
                let prev_model = abs_k(k - 1);
                let from_m = dp_m[prev] + hmm.transition(prev_model, TMD);
                let from_d = dp_d[prev] + hmm.transition(prev_model, TDD);

                if from_m >= from_d {
                    if from_m.is_finite() {
                        dp_d[cell] = from_m;
                        tb_d[cell] = TB_FROM_M;
                    }
                } else if from_d.is_finite() {
                    dp_d[cell] = from_d;
                    tb_d[cell] = TB_FROM_D;
                }
            }
        }
    }

    let end_cell = idx(l, m);
    let mut state = if dp_m[end_cell] >= dp_d[end_cell] {
        'm'
    } else {
        'd'
    };

    if !dp_m[end_cell].is_finite() && !dp_d[end_cell].is_finite() {
        return None;
    }

    let mut j = l;
    let mut k = m;
    let mut path = Vec::with_capacity(l + m);

    loop {
        let cell = idx(j, k);
        let model_k = abs_k(k);

        match state {
            'm' => {
                path.push((model_k, 'm'));
                if j == 1 && k == 1 {
                    break;
                }
                let tb = tb_m[cell];
                if tb == TB_START || j == 0 || k == 0 {
                    return None;
                }
                j -= 1;
                k -= 1;
                state = match tb {
                    TB_FROM_M => 'm',
                    TB_FROM_I => 'i',
                    TB_FROM_D => 'd',
                    _ => return None,
                };
            }
            'i' => {
                path.push((model_k, 'i'));
                let tb = tb_i[cell];
                if tb == TB_START || j == 0 {
                    return None;
                }
                j -= 1;
                state = match tb {
                    TB_FROM_M => 'm',
                    TB_FROM_I => 'i',
                    _ => return None,
                };
            }
            'd' => {
                path.push((model_k, 'd'));
                let tb = tb_d[cell];
                if tb == TB_START || k == 0 {
                    return None;
                }
                k -= 1;
                state = match tb {
                    TB_FROM_M => 'm',
                    TB_FROM_D => 'd',
                    _ => return None,
                };
            }
            _ => return None,
        }
    }

    path.reverse();

    let emitted = path.iter().filter(|(_, state)| *state != 'd').count();
    if emitted != l {
        return None;
    }

    Some(path)
}

fn stable_evalue_from_bits(bit_score: f64, hmm: &Plan7Hmm, db_size: usize) -> f64 {
    let score_nats = bit_score * std::f64::consts::LN_2;
    let z = hmm.viterbi_gumbel.lambda * (score_nats - hmm.viterbi_gumbel.mu);
    if z > 700.0 {
        0.0
    } else {
        db_size as f64 * (1.0 - (-(-z).exp()).exp())
    }
}

fn format_evalue_from_bits(bit_score: f64, hmm: &Plan7Hmm, db_size: usize) -> String {
    let score_nats = bit_score * std::f64::consts::LN_2;
    let z = hmm.viterbi_gumbel.lambda * (score_nats - hmm.viterbi_gumbel.mu);

    if z <= 50.0 {
        return format!("{:.1e}", stable_evalue_from_bits(bit_score, hmm, db_size));
    }

    let ln_evalue = (db_size as f64).ln() - z;
    let log10_evalue = ln_evalue / std::f64::consts::LN_10;
    let exponent = log10_evalue.floor() as i32;
    let mantissa = 10f64.powf(log10_evalue - exponent as f64);

    if mantissa.is_finite() {
        format!("{:.1}e{:+03}", mantissa, exponent)
    } else {
        "0.0e+00".to_string()
    }
}
