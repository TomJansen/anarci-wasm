use crate::fasta;
use crate::schemes::resolve_scheme_name;
use crate::{
    number_domains_from_hits, parse_allowed_chains, parse_allowed_species, parse_hmm_name,
    NumberingResult, ViterbiHit,
};
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

#[derive(Clone)]
struct HmmerAlignedHit {
    domain_index: usize,
    hmm_name: String,
    bit_score: f64,
    evalue: f64,
    evalue_text: String,
    seq_start: usize,
    seq_end: usize,
    hmm_start: usize,
    hmm_end: usize,
    rf: String,
    query_alignment: String,
}

/// Native reference backend: run HMMER and number from its alignment columns.
///
/// This mirrors original ANARCI's protein-input path more closely than the WASM
/// backend can, because HMMER emits the RF/query alignment columns that ANARCI's
/// `_hmm_alignment_to_states` consumes directly. It is intentionally unavailable
/// for `wasm32`, where subprocesses and `hmmscan` are not available.
pub fn number_sequences_with_hmmer_cli(
    fasta_text: &str,
    scheme: &str,
    bit_score_threshold: f64,
    allowed_species_json: &str,
    restrict: &str,
) -> Result<Vec<NumberingResult>, String> {
    let scheme_name = resolve_scheme_name(scheme).ok_or_else(|| "unknown scheme".to_string())?;
    let allowed_species = parse_allowed_species(allowed_species_json);
    let allowed_chains = parse_allowed_chains(restrict, scheme_name);
    let mut results = Vec::new();
    for (name, sequence) in fasta::parse_fasta(fasta_text) {
        if !fasta::validate_protein_sequence(&sequence) {
            results.push(NumberingResult {
                id: name,
                sequence,
                input_type: "protein".to_string(),
                domains: Vec::new(),
            });
            continue;
        }

        let output = run_hmmscan_for_sequence(&name, &sequence)?;
        let mut hits = parse_hmmer_output(&output, bit_score_threshold)?;
        hits.retain(|hit| {
            let (_, chain_type) = parse_hmm_name(&hit.hmm_name);
            allowed_chains.contains(chain_type.as_str())
        });
        hits = filter_hmmer_hits_by_allowed_species(hits, allowed_species.as_deref());

        let hits: Vec<ViterbiHit> = select_non_overlapping_hmmer_hits(&hits)
            .into_iter()
            .map(|hit| hmmer_hit_to_viterbi(&hit))
            .collect();

        let domains = number_domains_from_hits(
            sequence.as_bytes(),
            &hits,
            scheme_name,
            allowed_species.as_deref(),
            None,
            None,
        );

        results.push(NumberingResult {
            id: name,
            sequence,
            input_type: "protein".to_string(),
            domains,
        });
    }

    Ok(results)
}

fn run_hmmscan_for_sequence(name: &str, sequence: &str) -> Result<String, String> {
    require_tool("hmmpress")?;
    require_tool("hmmscan")?;

    let temp_dir = unique_temp_dir();
    fs::create_dir_all(&temp_dir).map_err(|err| err.to_string())?;

    let hmm = temp_dir.join("ALL.hmm");
    let query = temp_dir.join("query.fa");
    let output = temp_dir.join("hmmscan.txt");

    fs::copy(
        Path::new(env!("CARGO_MANIFEST_DIR")).join("data/ALL.hmm"),
        &hmm,
    )
    .map_err(|err| err.to_string())?;
    fs::write(&query, format!(">{name}\n{sequence}\n")).map_err(|err| err.to_string())?;

    let press = Command::new("hmmpress")
        .arg(&hmm)
        .output()
        .map_err(|err| err.to_string())?;
    if !press.status.success() {
        let _ = fs::remove_dir_all(&temp_dir);
        return Err(String::from_utf8_lossy(&press.stderr).to_string());
    }

    let scan = Command::new("hmmscan")
        .arg("-o")
        .arg(&output)
        .arg(&hmm)
        .arg(&query)
        .output()
        .map_err(|err| err.to_string())?;
    if !scan.status.success() {
        let _ = fs::remove_dir_all(&temp_dir);
        return Err(String::from_utf8_lossy(&scan.stderr).to_string());
    }

    let text = fs::read_to_string(&output).map_err(|err| err.to_string());
    let _ = fs::remove_dir_all(&temp_dir);
    text
}

fn parse_hmmer_output(
    text: &str,
    bit_score_threshold: f64,
) -> Result<Vec<HmmerAlignedHit>, String> {
    let mut hits: Vec<HmmerAlignedHit> = Vec::new();
    let mut current_model: Option<String> = None;
    let mut current_domain: Option<usize> = None;
    let mut saw_model_alignment_line = false;

    for line in text.lines() {
        let trimmed = line.trim();

        if let Some(model) = trimmed.strip_prefix(">> ") {
            current_model = Some(model.split_whitespace().next().unwrap_or("").to_string());
            current_domain = None;
            saw_model_alignment_line = false;
            continue;
        }

        let Some(model) = current_model.as_ref() else {
            continue;
        };

        if let Some(domain_no) = parse_domain_alignment_header(trimmed) {
            current_domain = Some(domain_no);
            saw_model_alignment_line = false;
            continue;
        }

        let fields: Vec<&str> = trimmed.split_whitespace().collect();
        if fields.len() >= 12 && fields[1] == "!" {
            if let Some(hit) = parse_domain_table_row(model, &fields, bit_score_threshold) {
                hits.push(hit);
            }
            continue;
        }

        if trimmed.ends_with(" RF") {
            if let Some(hit) = current_hit_mut(&mut hits, model, current_domain) {
                let rf = trimmed.trim_end_matches(" RF").trim();
                hit.rf.push_str(rf);
            }
            saw_model_alignment_line = false;
            continue;
        }

        if is_alignment_sequence_line(&fields) {
            if fields[0] == model {
                saw_model_alignment_line = true;
            } else if saw_model_alignment_line {
                if let Some(hit) = current_hit_mut(&mut hits, model, current_domain) {
                    hit.query_alignment.push_str(fields[2]);
                }
                saw_model_alignment_line = false;
            }
        }
    }

    for hit in &hits {
        if hit.rf.len() != hit.query_alignment.len() {
            return Err(format!(
                "HMMER alignment length mismatch for {} domain {}: RF {}, query {}",
                hit.hmm_name,
                hit.domain_index,
                hit.rf.len(),
                hit.query_alignment.len()
            ));
        }
    }

    Ok(hits)
}

fn parse_domain_table_row(
    model: &str,
    fields: &[&str],
    bit_score_threshold: f64,
) -> Option<HmmerAlignedHit> {
    let bit_score = fields.get(2)?.parse().ok()?;
    if bit_score < bit_score_threshold {
        return None;
    }

    Some(HmmerAlignedHit {
        domain_index: fields.first()?.parse().ok()?,
        hmm_name: model.to_string(),
        bit_score,
        evalue: parse_evalue(fields.get(5)?),
        evalue_text: fields.get(5)?.to_string(),
        hmm_start: fields.get(6)?.parse().ok()?,
        hmm_end: fields.get(7)?.parse().ok()?,
        seq_start: fields.get(9)?.parse::<usize>().ok()?.saturating_sub(1),
        seq_end: fields.get(10)?.parse().ok()?,
        rf: String::new(),
        query_alignment: String::new(),
    })
}

fn parse_domain_alignment_header(trimmed: &str) -> Option<usize> {
    let rest = trimmed.strip_prefix("== domain ")?;
    rest.split_whitespace().next()?.parse().ok()
}

fn current_hit_mut<'a>(
    hits: &'a mut [HmmerAlignedHit],
    model: &str,
    current_domain: Option<usize>,
) -> Option<&'a mut HmmerAlignedHit> {
    let domain = current_domain?;
    hits.iter_mut()
        .find(|hit| hit.hmm_name == model && hit.domain_index == domain)
}

fn is_alignment_sequence_line(fields: &[&str]) -> bool {
    fields.len() >= 4
        && fields[1].parse::<usize>().is_ok()
        && fields
            .last()
            .and_then(|value| value.parse::<usize>().ok())
            .is_some()
}

fn hmmer_hit_to_viterbi(hit: &HmmerAlignedHit) -> ViterbiHit {
    ViterbiHit {
        hmm_name: hit.hmm_name.clone(),
        bit_score: hit.bit_score,
        evalue: hit.evalue,
        evalue_text: hit.evalue_text.clone(),
        seq_start: hit.seq_start,
        seq_end: hit.seq_end,
        hmm_start: hit.hmm_start,
        hmm_end: hit.hmm_end,
        path: hmmer_alignment_to_path(hit),
    }
}

fn hmmer_alignment_to_path(hit: &HmmerAlignedHit) -> Vec<(usize, char)> {
    let hmm_states: Vec<usize> = (hit.hmm_start..=hit.hmm_end).collect();
    let mut h = 0usize;
    let mut path = Vec::new();

    for (rf, query) in hit.rf.bytes().zip(hit.query_alignment.bytes()) {
        let state_id = hmm_states
            .get(h)
            .copied()
            .unwrap_or_else(|| *hmm_states.last().unwrap_or(&hit.hmm_end));

        if query == b'-' {
            path.push((state_id, 'd'));
            h += 1;
        } else if rf == b'x' || rf == b'X' {
            path.push((state_id, 'm'));
            h += 1;
        } else {
            path.push((state_id, 'i'));
        }
    }

    path
}

fn select_non_overlapping_hmmer_hits(hits: &[HmmerAlignedHit]) -> Vec<HmmerAlignedHit> {
    let mut domains = Vec::new();
    for hit in hits {
        let overlaps = domains.iter().any(|domain: &HmmerAlignedHit| {
            hit.seq_start < domain.seq_end && domain.seq_start < hit.seq_end
        });
        if !overlaps {
            domains.push(hit.clone());
        }
    }
    domains.sort_by_key(|hit| hit.seq_start);
    domains
}

fn filter_hmmer_hits_by_allowed_species(
    hits: Vec<HmmerAlignedHit>,
    allowed_species: Option<&[String]>,
) -> Vec<HmmerAlignedHit> {
    let Some(allowed_species) = allowed_species else {
        return hits;
    };

    if allowed_species.is_empty() {
        return hits;
    }

    let filtered: Vec<HmmerAlignedHit> = hits
        .iter()
        .filter(|hit| {
            let (species, _) = parse_hmm_name(&hit.hmm_name);
            allowed_species.iter().any(|allowed| allowed == &species)
        })
        .cloned()
        .collect();

    if filtered.is_empty() {
        hits
    } else {
        filtered
    }
}

fn parse_evalue(text: &str) -> f64 {
    text.parse().unwrap_or(0.0)
}

fn require_tool(name: &str) -> Result<(), String> {
    let status = Command::new("which")
        .arg(name)
        .status()
        .map_err(|err| err.to_string())?;
    if status.success() {
        Ok(())
    } else {
        Err(format!("{name} is required"))
    }
}

fn unique_temp_dir() -> PathBuf {
    let suffix = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    std::env::temp_dir().join(format!("anarci-hmmer-cli-{suffix}"))
}

#[cfg(test)]
mod tests {
    use super::{hmmer_alignment_to_path, HmmerAlignedHit};

    #[test]
    fn hmmer_alignment_columns_convert_to_match_insert_delete_path() {
        let hit = HmmerAlignedHit {
            domain_index: 1,
            hmm_name: "human_H".to_string(),
            bit_score: 1.0,
            evalue: 0.0,
            evalue_text: "0".to_string(),
            seq_start: 0,
            seq_end: 3,
            hmm_start: 10,
            hmm_end: 12,
            rf: "xx.x".to_string(),
            query_alignment: "ABc-".to_string(),
        };

        assert_eq!(
            hmmer_alignment_to_path(&hit),
            vec![(10, 'm'), (11, 'm'), (12, 'i'), (12, 'd')]
        );
    }
}
