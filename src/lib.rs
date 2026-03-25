mod alignment;
mod fasta;
mod germline;
mod hmm;
mod schemes;

use serde::Serialize;
use wasm_bindgen::prelude::*;

use alignment::{path_to_state_vector, StateVectorEntry};
use fasta::{InputSequenceType, TranslatedFrame};
use germline::{assign_germline, available_species as germline_species};
use hmm::{HmmDatabase, ViterbiHit};
use schemes::{chain_type_to_class, number_sequence, resolve_scheme_name};

use std::cmp::Ordering;
use std::collections::BTreeSet;
use std::sync::OnceLock;

/// Embedded HMM database (compiled into the WASM binary)
const HMM_DATA: &str = include_str!("../data/ALL.hmm");

/// Lazily parsed HMM database
static HMM_DB: OnceLock<HmmDatabase> = OnceLock::new();

fn get_db() -> &'static HmmDatabase {
    HMM_DB.get_or_init(|| HmmDatabase::parse(HMM_DATA))
}

#[derive(Serialize)]
pub struct NumberingResult {
    pub id: String,
    pub sequence: String,
    pub input_type: String,
    pub domains: Vec<DomainResult>,
}

#[derive(Serialize)]
pub struct DomainResult {
    pub domain_index: usize,
    pub species: String,
    pub chain_type: String,
    pub chain_class: String,
    pub identity_species: String,
    pub v_gene: String,
    pub v_identity: f64,
    pub j_gene: String,
    pub j_identity: f64,
    pub scheme: String,
    pub bit_score: f64,
    pub evalue: f64,
    pub evalue_text: String,
    pub seq_start: usize,
    pub seq_end: usize,
    pub translation_frame: String,
    pub nt_start: Option<usize>,
    pub nt_end: Option<usize>,
    pub numbering: Vec<NumberingEntry>,
}

#[derive(Serialize)]
pub struct NumberingEntry {
    pub position: i32,
    pub insertion: String,
    pub amino_acid: String,
}

#[wasm_bindgen]
pub fn number_sequences(fasta_text: &str, scheme: &str, bit_score_threshold: f64) -> String {
    number_sequences_with_options(
        fasta_text,
        scheme,
        bit_score_threshold,
        "[]",
        "",
        "protein",
    )
}

#[wasm_bindgen]
pub fn number_sequences_with_species(
    fasta_text: &str,
    scheme: &str,
    bit_score_threshold: f64,
    allowed_species_json: &str,
) -> String {
    number_sequences_with_options(
        fasta_text,
        scheme,
        bit_score_threshold,
        allowed_species_json,
        "",
        "protein",
    )
}

#[wasm_bindgen]
pub fn number_sequences_with_options(
    fasta_text: &str,
    scheme: &str,
    bit_score_threshold: f64,
    allowed_species_json: &str,
    restrict: &str,
    input_type: &str,
) -> String {
    let scheme_name = match resolve_scheme_name(scheme) {
        Some(s) => s,
        None => return serde_json::to_string(&Vec::<NumberingResult>::new()).unwrap(),
    };
    let input_type = parse_input_type(input_type);
    let allowed_species = parse_allowed_species(allowed_species_json);
    let allowed_chains = parse_allowed_chains(restrict, scheme_name);

    let db = get_db();
    let sequences = fasta::parse_fasta(fasta_text);
    let mut results: Vec<NumberingResult> = Vec::new();

    for (name, seq_str) in &sequences {
        let result = match input_type {
            InputSequenceType::Protein => process_protein_sequence(
                db,
                name,
                seq_str,
                scheme_name,
                bit_score_threshold,
                allowed_species.as_deref(),
                &allowed_chains,
            ),
            InputSequenceType::Dna => process_dna_sequence(
                db,
                name,
                seq_str,
                scheme_name,
                bit_score_threshold,
                allowed_species.as_deref(),
                &allowed_chains,
            ),
        };
        results.push(result);
    }

    serde_json::to_string(&results).unwrap_or_else(|_| "[]".to_string())
}

fn process_protein_sequence(
    db: &HmmDatabase,
    name: &str,
    sequence: &str,
    scheme_name: &str,
    bit_score_threshold: f64,
    allowed_species: Option<&[String]>,
    allowed_chains: &BTreeSet<String>,
) -> NumberingResult {
    if !fasta::validate_protein_sequence(sequence) {
        return NumberingResult {
            id: name.to_string(),
            sequence: sequence.to_string(),
            input_type: "protein".to_string(),
            domains: Vec::new(),
        };
    }

    let domains = number_domains(
        db,
        sequence.as_bytes(),
        scheme_name,
        bit_score_threshold,
        allowed_species,
        allowed_chains,
        None,
        None,
    );

    NumberingResult {
        id: name.to_string(),
        sequence: sequence.to_string(),
        input_type: "protein".to_string(),
        domains,
    }
}

fn process_dna_sequence(
    db: &HmmDatabase,
    name: &str,
    sequence: &str,
    scheme_name: &str,
    bit_score_threshold: f64,
    allowed_species: Option<&[String]>,
    allowed_chains: &BTreeSet<String>,
) -> NumberingResult {
    if !fasta::validate_dna_sequence(sequence) {
        return NumberingResult {
            id: name.to_string(),
            sequence: sequence.to_string(),
            input_type: "dna".to_string(),
            domains: Vec::new(),
        };
    }

    let translated = fasta::translate_six_frames(sequence);
    let mut best_translation: Option<&TranslatedFrame> = None;
    let mut best_domains: Vec<ViterbiHit> = Vec::new();
    let mut best_count = 0usize;
    let mut best_total = 0.0f64;
    let mut best_top = f64::NEG_INFINITY;

    for frame in &translated {
        if !fasta::validate_protein_sequence(&frame.sequence) {
            continue;
        }

        let hits = search_filtered_hits(
            db,
            frame.sequence.as_bytes(),
            bit_score_threshold,
            allowed_species,
            allowed_chains,
        );
        let domains = select_domains(&hits);
        if domains.is_empty() {
            continue;
        }

        let count = domains.len();
        let total = domains.iter().map(|hit| hit.bit_score).sum::<f64>();
        let top = domains
            .first()
            .map(|hit| hit.bit_score)
            .unwrap_or(f64::NEG_INFINITY);

        let ordering = count
            .cmp(&best_count)
            .then_with(|| total.total_cmp(&best_total))
            .then_with(|| top.total_cmp(&best_top));

        if ordering == Ordering::Greater {
            best_translation = Some(frame);
            best_domains = domains;
            best_count = count;
            best_total = total;
            best_top = top;
        }
    }

    let domains = best_translation
        .map(|frame| {
            number_domains_from_hits(
                db,
                frame.sequence.as_bytes(),
                &best_domains,
                scheme_name,
                allowed_species,
                Some(frame),
                Some(sequence.len()),
            )
        })
        .unwrap_or_default();

    NumberingResult {
        id: name.to_string(),
        sequence: sequence.to_string(),
        input_type: "dna".to_string(),
        domains,
    }
}

fn number_domains(
    db: &HmmDatabase,
    sequence: &[u8],
    scheme_name: &str,
    bit_score_threshold: f64,
    allowed_species: Option<&[String]>,
    allowed_chains: &BTreeSet<String>,
    translation: Option<&TranslatedFrame>,
    original_dna_len: Option<usize>,
) -> Vec<DomainResult> {
    let hits = search_filtered_hits(
        db,
        sequence,
        bit_score_threshold,
        allowed_species,
        allowed_chains,
    );
    let domains = select_domains(&hits);
    number_domains_from_hits(
        db,
        sequence,
        &domains,
        scheme_name,
        allowed_species,
        translation,
        original_dna_len,
    )
}

fn number_domains_from_hits(
    db: &HmmDatabase,
    sequence: &[u8],
    domains: &[ViterbiHit],
    scheme_name: &str,
    allowed_species: Option<&[String]>,
    translation: Option<&TranslatedFrame>,
    original_dna_len: Option<usize>,
) -> Vec<DomainResult> {
    let mut domain_results = Vec::new();

    for (di, hit) in domains.iter().enumerate() {
        let (species, chain_type) = parse_hmm_name(&hit.hmm_name);
        let state_vector = path_to_state_vector(hit, sequence.len(), di, domains.len());
        let state_vector = rescue_missing_j_region(db, sequence, &state_vector, domains.len());
        let germline = assign_germline(&state_vector, sequence, &chain_type, allowed_species);

        if let Some((numbered, start, end)) =
            number_sequence(&state_vector, sequence, scheme_name, &chain_type)
        {
            let entries: Vec<NumberingEntry> = numbered
                .iter()
                .map(|&((pos, ref ins), aa)| NumberingEntry {
                    position: pos,
                    insertion: ins.clone(),
                    amino_acid: aa.to_string(),
                })
                .collect();

            let (translation_frame, nt_start, nt_end) =
                if let (Some(frame), Some(dna_len), Some(start), Some(end)) =
                    (translation, original_dna_len, start, end)
                {
                    let (nt_start, nt_end) =
                        aa_hit_to_nt_range(frame, dna_len, start, end + 1);
                    (frame.label.clone(), Some(nt_start), Some(nt_end))
                } else {
                    (String::new(), None, None)
                };

            let seq_start = start.unwrap_or(hit.seq_start);
            let seq_end = end.map(|value| value + 1).unwrap_or(hit.seq_end);

            domain_results.push(DomainResult {
                domain_index: di,
                species: species.clone(),
                chain_type: chain_type.clone(),
                chain_class: chain_type_to_class(&chain_type).to_string(),
                identity_species: germline
                    .as_ref()
                    .map(|g| g.identity_species.clone())
                    .unwrap_or_default(),
                v_gene: germline
                    .as_ref()
                    .map(|g| g.v_gene.clone())
                    .unwrap_or_default(),
                v_identity: germline.as_ref().map(|g| g.v_identity).unwrap_or(0.0),
                j_gene: germline
                    .as_ref()
                    .map(|g| g.j_gene.clone())
                    .unwrap_or_default(),
                j_identity: germline.as_ref().map(|g| g.j_identity).unwrap_or(0.0),
                scheme: scheme_name.to_string(),
                bit_score: hit.bit_score,
                evalue: hit.evalue,
                evalue_text: hit.evalue_text.clone(),
                seq_start,
                seq_end,
                translation_frame,
                nt_start,
                nt_end,
                numbering: entries,
            });
        }
    }

    domain_results
}

fn search_filtered_hits(
    db: &HmmDatabase,
    sequence: &[u8],
    bit_score_threshold: f64,
    allowed_species: Option<&[String]>,
    allowed_chains: &BTreeSet<String>,
) -> Vec<ViterbiHit> {
    let hits: Vec<ViterbiHit> = db
        .search(sequence, bit_score_threshold)
        .into_iter()
        .filter(|hit| {
            let (_, chain_type) = parse_hmm_name(&hit.hmm_name);
            allowed_chains.contains(chain_type.as_str())
        })
        .collect();

    filter_hits_by_allowed_species(hits, allowed_species)
}

fn parse_input_type(input_type: &str) -> InputSequenceType {
    match input_type.trim().to_ascii_lowercase().as_str() {
        "dna" => InputSequenceType::Dna,
        _ => InputSequenceType::Protein,
    }
}

fn aa_hit_to_nt_range(
    translation: &TranslatedFrame,
    dna_len: usize,
    aa_start: usize,
    aa_end: usize,
) -> (usize, usize) {
    let translated_nt_start = translation.nt_offset + aa_start * 3;
    let translated_nt_end = translation.nt_offset + aa_end * 3;

    if translation.is_reverse {
        (dna_len - translated_nt_end, dna_len - translated_nt_start)
    } else {
        (translated_nt_start, translated_nt_end)
    }
}

fn parse_allowed_species(allowed_species_json: &str) -> Option<Vec<String>> {
    let parsed: Vec<String> = serde_json::from_str(allowed_species_json).unwrap_or_default();
    let normalized: Vec<String> = parsed
        .into_iter()
        .map(|species| species.trim().to_ascii_lowercase())
        .filter(|species| !species.is_empty())
        .collect();

    if normalized.is_empty() {
        None
    } else {
        Some(normalized)
    }
}

fn parse_allowed_chains(restrict: &str, scheme_name: &str) -> BTreeSet<String> {
    let mut allow: BTreeSet<String> = match restrict.trim().to_ascii_lowercase().as_str() {
        "" => ["H", "K", "L", "A", "B", "G", "D"]
            .into_iter()
            .map(String::from)
            .collect(),
        "ig" => ["H", "K", "L"].into_iter().map(String::from).collect(),
        "tr" => ["A", "B", "G", "D"].into_iter().map(String::from).collect(),
        "heavy" => ["H"].into_iter().map(String::from).collect(),
        "light" => ["K", "L"].into_iter().map(String::from).collect(),
        "h" => ["H"].into_iter().map(String::from).collect(),
        "k" => ["K"].into_iter().map(String::from).collect(),
        "l" => ["L"].into_iter().map(String::from).collect(),
        "a" => ["A"].into_iter().map(String::from).collect(),
        "b" => ["B"].into_iter().map(String::from).collect(),
        _ => ["H", "K", "L", "A", "B", "G", "D"]
            .into_iter()
            .map(String::from)
            .collect(),
    };

    if !matches!(scheme_name, "imgt" | "aho") {
        allow.retain(|chain| matches!(chain.as_str(), "H" | "K" | "L"));
    }

    allow
}

fn filter_hits_by_allowed_species(
    hits: Vec<ViterbiHit>,
    allowed_species: Option<&[String]>,
) -> Vec<ViterbiHit> {
    let Some(allowed_species) = allowed_species else {
        return hits;
    };

    if allowed_species.is_empty() {
        return hits;
    }

    let filtered: Vec<ViterbiHit> = hits
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

fn stitch_rescued_j_region(
    state_vector: &[StateVectorEntry],
    cys_alignment_index: usize,
    cys_seq_index: usize,
    j_region: &[StateVectorEntry],
) -> Vec<StateVectorEntry> {
    let Some(j_start_seq_index) = j_region.first().and_then(|entry| entry.1) else {
        return state_vector.to_vec();
    };

    let mut rescued = state_vector[..=cys_alignment_index].to_vec();
    let mut next_state = 105usize;
    for seq_index in (cys_seq_index + 1)..j_start_seq_index {
        let state = if next_state >= 116 {
            (116, 'i')
        } else {
            let state = (next_state, 'm');
            next_state += 1;
            state
        };
        rescued.push((state, Some(seq_index)));
    }
    rescued.extend_from_slice(j_region);
    rescued
}

fn rescue_missing_j_region(
    db: &HmmDatabase,
    sequence: &[u8],
    state_vector: &[StateVectorEntry],
    n_domains: usize,
) -> Vec<StateVectorEntry> {
    if n_domains != 1 || state_vector.is_empty() {
        return state_vector.to_vec();
    }

    let Some(&((last_state, _), Some(last_seq_index))) = state_vector
        .iter()
        .rev()
        .find(|entry| entry.1.is_some())
    else {
        return state_vector.to_vec();
    };

    if last_state >= 120 || last_seq_index + 30 >= sequence.len() {
        return state_vector.to_vec();
    }

    let Some(cys_seq_index) = state_vector
        .iter()
        .find_map(|entry| (entry.0 == (104, 'm')).then_some(entry.1).flatten())
    else {
        return state_vector.to_vec();
    };

    let Some(cys_alignment_index) = state_vector
        .iter()
        .position(|entry| *entry == ((104, 'm'), Some(cys_seq_index)))
    else {
        return state_vector.to_vec();
    };

    if cys_seq_index + 1 >= sequence.len() {
        return state_vector.to_vec();
    }

    let remaining = &sequence[cys_seq_index + 1..];
    let rescue_hits = select_domains(&db.search(remaining, 10.0));
    let Some(rescue_hit) = rescue_hits.first() else {
        return state_vector.to_vec();
    };

    let rescue_state_vector = path_to_state_vector(rescue_hit, remaining.len(), 0, rescue_hits.len());
    let Some(first_rescue_state) = rescue_state_vector.first().map(|entry| entry.0 .0) else {
        return state_vector.to_vec();
    };
    let Some(last_rescue_state) = rescue_state_vector.last().map(|entry| entry.0 .0) else {
        return state_vector.to_vec();
    };

    if first_rescue_state > 117 || last_rescue_state < 126 {
        return state_vector.to_vec();
    }

    let j_region: Vec<StateVectorEntry> = rescue_state_vector
        .iter()
        .filter_map(|&((state_id, state_type), seq_index)| {
            if state_id < 117 {
                return None;
            }
            let seq_index = seq_index?;
            Some(((state_id, state_type), Some(seq_index + cys_seq_index + 1)))
        })
        .collect();

    stitch_rescued_j_region(state_vector, cys_alignment_index, cys_seq_index, &j_region)
}

/// Select non-overlapping domains from hits (best score wins)
fn select_domains(hits: &[ViterbiHit]) -> Vec<ViterbiHit> {
    let mut domains: Vec<ViterbiHit> = Vec::new();
    for hit in hits {
        let overlaps = domains.iter().any(|d| {
            let (s1, e1) = (d.seq_start, d.seq_end);
            let (s2, e2) = (hit.seq_start, hit.seq_end);
            s2 < e1 && s1 < e2
        });
        if !overlaps {
            domains.push(hit.clone());
        }
    }
    domains.sort_by_key(|d| d.seq_start);
    domains
}

/// Parse HMM name like "human_H" into (species, chain_type)
fn parse_hmm_name(name: &str) -> (String, String) {
    if let Some((species, chain)) = name.split_once('_') {
        (species.to_string(), chain.to_string())
    } else {
        ("unknown".to_string(), "H".to_string())
    }
}

#[wasm_bindgen]
pub fn available_schemes() -> String {
    serde_json::to_string(&["imgt", "chothia", "kabat", "martin", "aho", "wolfguy"])
        .unwrap_or_else(|_| "[]".to_string())
}

#[wasm_bindgen]
pub fn available_species() -> String {
    serde_json::to_string(&germline_species()).unwrap_or_else(|_| "[]".to_string())
}

#[wasm_bindgen]
pub fn num_profiles() -> usize {
    get_db().profiles.len()
}

#[cfg(test)]
mod tests {
    use super::{
        filter_hits_by_allowed_species, number_sequences_with_options, stitch_rescued_j_region,
    };
    use crate::hmm::ViterbiHit;

    fn mock_hit(hmm_name: &str, seq_start: usize, seq_end: usize, bit_score: f64) -> ViterbiHit {
        ViterbiHit {
            hmm_name: hmm_name.to_string(),
            bit_score,
            evalue: 0.0,
            evalue_text: "0.0e+00".to_string(),
            seq_start,
            seq_end,
            hmm_start: 1,
            hmm_end: 128,
            path: Vec::new(),
        }
    }

    #[test]
    fn dna_vhh_sequence_is_detected_in_translated_frame() {
        let fasta = concat!(
            ">TPL1068_01_B06_M13rev-29_A01\n",
            "ATGAAAAAGACAGCTATCGCGATTGCAGTGGCACTGGCTGGTTTCGCTACCGTAGCGCAGGCCCAGGTGCAGCTGCAGGA\n",
            "GTCTGGAGGAGGCCTGGTGCAGCCCGGGGGGTCTCTGAGACTCTCCTGTGTAGTCTCTGGTGGACTCAGTACGGA\n",
            "TGATTACACCATGGGCTGGTTCCGCCAGGTTCCAGGGAAGGAGCGTGAGGGAATTGCATTAAGAA\n",
            "GTACGAGTGGTGCCGGTACCAGGTATGCGGACTCCGTGAAGGGCCGATTCACCATCTCCGGGGACAAGTCCAAGAA\n",
            "CACGGTGTATCTGCAAATGAACAACCTGAAACCAGAGGACACAGCCGTCTATTACTGTGCGGTAGGGCCGTGGTCGACAC\n",
            "ACCCGGGTCTAGGACAGTCTTTCTATGGCGACTGGGGCCAGGGGACCCAGGTCACCGTCTCCTCAGCG\n",
            "GCCGCAGACTACAAAGACCATGACGGTGATTATAAAGATCATGACATCGATTACAAGGATGACGATGACAAGGGTGCCGC\n",
            "ACACCACCATCACCATCACTAATAGAATTCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCGTTACCC\n",
            "AACTTAATCGCCTTGCAGCACATCCCCCTTTCGCCAGCTGGCGTAATAGCGAAGAGGCCCGCACCGATCGCCCTTCCCAA\n",
            "CAGTTGCGCAGCCTGAATGGCGAATGGCGCCTGATGCGGTATTTTCTCCTTACGCATCTGTGCGGTATTTCACACCGCAT\n",
            "ACGTCAAAGCAACCATAGTACGCGCCCTGTAGCGGCGCATTAAGCGCGGCGGGTGTGGTGGTTACGCGCAGCGTGACCGC\n",
            "TACACTTGCCAGCGCCTTAGCGCCCGCTCCTTTCGCTTTCTTCCCTTCCTATCTCGCCACGTTCGCCGGCTTTCCCCGTC\n",
            "AAGCTCTAAATCGCGGACTCCCTTTAGGGTTCCGATTAATGCTTTACGGC\n"
        );

        let json = number_sequences_with_options(fasta, "imgt", 80.0, "[]", "heavy", "dna");
        let results: serde_json::Value = serde_json::from_str(&json).unwrap();
        let domains = &results[0]["domains"];

        assert_eq!(domains.as_array().unwrap().len(), 1);
        assert_eq!(domains[0]["translation_frame"], "+1");
    }

    #[test]
    fn gapped_dna_vhh_sequence_is_detected_in_translated_frame() {
        let fasta = concat!(
            ">TPL1068_01_B06_M13rev-29_A01\n",
            "ATGAAAAAGACAGCTATCGCGATTGCAGTGGCACTGGCTGGTTTCGCTACCGTAGCGCAGGCCCAGGTGCAGCTGCAGGA\n",
            "GTCTGGAGGAGGCCTGGTGCAGCCCGGGGGGTCTCTGAGACTCTCCTGTGTAGTCTCTGGTGGACTCAGTACGGA\n",
            "TGATTACACCATGGGCTGGTTCCGCCAGGTTCCAGGGAAGGAGCGTGAGGGAATTGCATTAAGAA\n",
            "GTACGAGTGGTGCCGGTACCAGGTATGCGGACTCCGTGAAGGGCCGATTCACCATCTCCGGGGACAAGTCCAAGAA\n",
            "CACGGTGTATCTGCAAATGAACAACCTGAAACCAGAGGACACAGCCGTCTATTACTGTGCGGTAGGGCCGTGGTCGACAC\n",
            "ACCCGGGTCTAGGACAGTCTTTCTATGGCGACTGGGGCCAGGGGACCCAGGTCACCGTCTCCTCAGCG\n",
            "GCCGCAGACTACAAAGACCATGACGGTGATTATAAAGATCATGACATCGATTACAAGGATGACGATGACAAGGGTGCCGC\n",
            "ACACCACCATCACCATCACTAATAGAATTCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCGTTACCC\n",
            "AACTTAATCGCCTTGCAGCACATCCCCCTTTCGCCAGCTGGCGTAATAGCGAAGAGGCCCGCACCGATCGCCCTTCCCAA\n",
            "CAGTTGCGCAGCCTGAATGGCGAATGGCGCCTGATGCGGTATTTTCTCCTTACGCATCTGTGCGGTATTTCACACCGCAT\n",
            "ACGTCAAAGCAACCATAGTACGCGCCCTGTAGCGGCGCATTAAGCGCGGCGGGTGTGGTGGTTACGCGCAGCGTGACCGC\n",
            "TACACTTGCCAGCGCCTTAGCGCCCGCTCCTTTCGCTTTCTTCCCTTCCTATCTCGCCACGTTCGCCGGCTTTCCCCGTC\n",
            "AAGCTCTAAATCGCGGACTCCCTTTAGGGTTCCGATTAATGCTTTACGGC\n"
        );
        let gapped = fasta
            .lines()
            .map(|line| {
                if line.starts_with('>') {
                    line.to_string()
                } else {
                    let mut with_gaps = String::new();
                    for (index, ch) in line.chars().enumerate() {
                        if index > 0 && index % 17 == 0 {
                            with_gaps.push('-');
                        }
                        with_gaps.push(ch);
                    }
                    with_gaps
                }
            })
            .collect::<Vec<_>>()
            .join("\n");

        let json = number_sequences_with_options(&gapped, "imgt", 80.0, "[]", "heavy", "dna");
        let results: serde_json::Value = serde_json::from_str(&json).unwrap();
        let domains = &results[0]["domains"];

        assert_eq!(domains.as_array().unwrap().len(), 1);
        assert_eq!(domains[0]["translation_frame"], "+1");
    }

    #[test]
    fn allowed_species_prefers_requested_hits_when_available() {
        let hits = vec![
            mock_hit("human_H", 0, 120, 100.0),
            mock_hit("mouse_H", 0, 120, 95.0),
            mock_hit("human_K", 0, 110, 90.0),
        ];

        let filtered = filter_hits_by_allowed_species(hits, Some(&["mouse".to_string()]));

        assert_eq!(filtered.len(), 1);
        assert_eq!(filtered[0].hmm_name, "mouse_H");
    }

    #[test]
    fn allowed_species_falls_back_to_all_hits_when_no_match_exists() {
        let hits = vec![
            mock_hit("human_H", 0, 120, 100.0),
            mock_hit("mouse_H", 0, 120, 95.0),
        ];

        let filtered = filter_hits_by_allowed_species(hits.clone(), Some(&["rabbit".to_string()]));

        assert_eq!(filtered.len(), hits.len());
        assert_eq!(filtered[0].hmm_name, hits[0].hmm_name);
        assert_eq!(filtered[1].hmm_name, hits[1].hmm_name);
    }

    #[test]
    fn missing_j_region_rescue_stitches_v_and_j_segments() {
        let state_vector = vec![
            ((100, 'm'), Some(0)),
            ((101, 'm'), Some(1)),
            ((102, 'm'), Some(2)),
            ((103, 'm'), Some(3)),
            ((104, 'm'), Some(4)),
            ((110, 'm'), Some(5)),
        ];
        let j_region = vec![
            ((117, 'm'), Some(20)),
            ((118, 'm'), Some(21)),
            ((126, 'm'), Some(22)),
        ];

        let rescued = stitch_rescued_j_region(&state_vector, 4, 4, &j_region);

        assert_eq!(rescued[0], ((100, 'm'), Some(0)));
        assert_eq!(rescued[4], ((104, 'm'), Some(4)));
        assert_eq!(rescued[5], ((105, 'm'), Some(5)));
        assert_eq!(rescued[15], ((115, 'm'), Some(15)));
        assert_eq!(rescued[16], ((116, 'i'), Some(16)));
        assert_eq!(rescued[19], ((116, 'i'), Some(19)));
        assert_eq!(rescued[20], ((117, 'm'), Some(20)));
        assert_eq!(rescued[22], ((126, 'm'), Some(22)));
    }
}
