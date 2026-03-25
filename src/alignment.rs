/// State vector extraction and core numbering logic.
/// Ports the key functions from ANARCI's anarci.py and schemes.py.

use crate::hmm::ViterbiHit;

/// Insertion alphabet: A-Z, AA-ZZ, then space for "no insertion"
pub const ALPHABET: &[&str] = &[
    "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R",
    "S", "T", "U", "V", "W", "X", "Y", "Z", "AA", "BB", "CC", "DD", "EE", "FF", "GG", "HH",
    "II", "JJ", "KK", "LL", "MM", "NN", "OO", "PP", "QQ", "RR", "SS", "TT", "UU", "VV", "WW",
    "XX", "YY", "ZZ", " ",
];

/// A state in the alignment: (hmm_state_1based, state_type)
/// state_type: 'm' = match, 'i' = insert, 'd' = delete
pub type State = (usize, char);

/// A state vector entry: (state, sequence_index_or_none)
pub type StateVectorEntry = (State, Option<usize>);

/// A numbering annotation: ((position_number, insertion_code), amino_acid)
pub type NumberedResidue = ((i32, String), char);

/// Convert a Viterbi path into a state vector with sequence indices.
/// Equivalent to `_hmm_alignment_to_states` in ANARCI.
pub fn path_to_state_vector(hit: &ViterbiHit, seq_len: usize, domain_index: usize, n_domains: usize) -> Vec<StateVectorEntry> {
    let mut state_vector = Vec::new();
    let mut seq_idx = hit.seq_start;

    let hmm_start = hit.hmm_start;

    // Handle N-terminal extension for first domain
    // If the alignment doesn't start at state 1 and we're the first domain,
    // extend backwards to capture unmatched N-terminal residues
    if domain_index == 0 && hmm_start > 1 && hmm_start <= 5 {
        let n_extend = std::cmp::min(hmm_start - 1, hit.seq_start);
        let effective_hmm_start = hmm_start - n_extend;
        let effective_seq_start = hit.seq_start - n_extend;

        // Add N-terminal extension as match states
        for offset in 0..n_extend {
            state_vector.push((
                (effective_hmm_start + offset, 'm'),
                Some(effective_seq_start + offset),
            ));
        }
        seq_idx = hit.seq_start;
    }

    // Handle C-terminal extension for single domains
    let n_c_extend = if n_domains == 1 && hit.hmm_end < 128 && hit.hmm_end > 123 {
        let remaining_seq = seq_len - hit.seq_end;
        let remaining_hmm = 128 - hit.hmm_end;
        std::cmp::min(remaining_seq, remaining_hmm)
    } else {
        0
    };

    // Process the main path
    for &(state, state_type) in &hit.path {
        match state_type {
            'm' | 'i' => {
                state_vector.push(((state, state_type), Some(seq_idx)));
                seq_idx += 1;
            }
            'd' => {
                state_vector.push(((state, 'd'), None));
            }
            _ => {}
        }
    }

    // Add C-terminal extension
    for offset in 0..n_c_extend {
        let hmm_state = hit.hmm_end + 1 + offset;
        state_vector.push(((hmm_state, 'm'), Some(hit.seq_end + offset)));
    }

    state_vector
}

/// Smooth insertions in the state vector to correct HMM alignment artifacts.
/// Equivalent to `smooth_insertions` in ANARCI schemes.py.
pub fn smooth_insertions(state_vector: Vec<StateVectorEntry>) -> Vec<StateVectorEntry> {
    // Enforced patterns at CDR boundaries
    let enforced_patterns: Vec<Vec<(usize, char)>> = vec![
        vec![(25, 'm'), (26, 'm'), (27, 'm'), (28, 'i')],
        vec![(38, 'i'), (38, 'm'), (39, 'm'), (40, 'm')],
        vec![(54, 'm'), (55, 'm'), (56, 'm'), (57, 'i')],
        vec![(65, 'i'), (65, 'm'), (66, 'm'), (67, 'm')],
        vec![(103, 'm'), (104, 'm'), (105, 'm'), (106, 'i')],
        vec![(117, 'i'), (117, 'm'), (118, 'm'), (119, 'm')],
    ];

    let mut state_buffer: Vec<StateVectorEntry> = Vec::new();
    let mut sv: Vec<StateVectorEntry> = Vec::new();
    let mut reg: i32 = -1;

    for entry in &state_vector {
        let &((state_id, _state_type), _si) = entry;

        let in_region = if state_id < 23 {
            reg = -1;
            true
        } else if (25..28).contains(&state_id) {
            reg = 0;
            true
        } else if state_id > 37 && state_id <= 40 {
            reg = 1;
            true
        } else if (54..57).contains(&state_id) {
            reg = 2;
            true
        } else if state_id > 64 && state_id <= 67 {
            reg = 3;
            true
        } else if (103..106).contains(&state_id) {
            reg = 4;
            true
        } else if state_id > 116 && state_id <= 119 {
            reg = 5;
            true
        } else {
            false
        };

        if in_region {
            state_buffer.push(*entry);
        } else if !state_buffer.is_empty() {
            let nins = state_buffer.iter().filter(|s| s.0 .1 == 'i').count();

            if nins > 0 {
                if reg == -1 {
                    // FW1 handling
                    let mut nt_dels = state_buffer[0].0 .0 - 1;
                    for &((_, st), si_opt) in &state_buffer {
                        if st == 'd' || si_opt.is_none() {
                            nt_dels += 1;
                        } else {
                            break;
                        }
                    }

                    if nt_dels >= nins {
                        let match_states: Vec<(usize, char)> = state_buffer
                            .iter()
                            .filter(|s| s.0 .1 == 'm')
                            .map(|s| s.0)
                            .collect();
                        let first = match_states.first().map(|s| s.0).unwrap_or(1);

                        let filtered: Vec<StateVectorEntry> = state_buffer
                            .iter()
                            .filter(|s| s.0 .1 != 'd')
                            .copied()
                            .collect();

                        let add = filtered.len().saturating_sub(match_states.len());
                        let mut new_states: Vec<(usize, char)> = (first.saturating_sub(add)..first)
                            .map(|s| (s, 'm'))
                            .collect();
                        new_states.extend_from_slice(&match_states);

                        for i in 0..filtered.len() {
                            if i < new_states.len() {
                                sv.push((new_states[i], filtered[i].1));
                            }
                        }
                    } else {
                        sv.extend_from_slice(&state_buffer);
                    }
                } else {
                    // CDR boundary handling
                    let filtered: Vec<StateVectorEntry> = state_buffer
                        .iter()
                        .filter(|s| s.0 .1 != 'd')
                        .copied()
                        .collect();

                    let pattern = &enforced_patterns[reg as usize];
                    let buf_len = filtered.len();

                    let new_states: Vec<(usize, char)> = if reg % 2 == 1 {
                        // N-terminal FW
                        let prefix: Vec<(usize, char)> =
                            vec![pattern[0]; buf_len.saturating_sub(3)];
                        let skip = if buf_len >= 4 { 1 } else { 4 - buf_len.min(4) + 1 };
                        let suffix: Vec<(usize, char)> = pattern[skip.min(4)..].to_vec();
                        [prefix, suffix].concat()
                    } else {
                        // C-terminal FW
                        let prefix = pattern[..3].to_vec();
                        let extra: Vec<(usize, char)> =
                            vec![pattern[2]; buf_len.saturating_sub(3)];
                        [prefix, extra].concat()
                    };

                    for i in 0..filtered.len() {
                        if i < new_states.len() {
                            sv.push((new_states[i], filtered[i].1));
                        }
                    }
                }
            } else {
                sv.extend_from_slice(&state_buffer);
            }

            sv.push(*entry);
            state_buffer.clear();
        } else {
            sv.push(*entry);
        }
    }

    // Handle any remaining buffer
    if !state_buffer.is_empty() {
        sv.extend_from_slice(&state_buffer);
    }

    sv
}

/// Core region-based numbering function.
/// Equivalent to `_number_regions` in ANARCI schemes.py.
///
/// # Parameters
/// - `sequence`: the original sequence bytes
/// - `state_vector`: aligned state vector
/// - `state_string`: 128-char string, 'X' = scheme position, 'I' = insertion
/// - `region_string`: 128-char string, region labels
/// - `region_index_map`: maps region label char to region index
/// - `rels`: relative numbering offsets per region (mutated during processing)
/// - `n_regions`: number of regions
/// - `exclude_deletions`: region indices where deletions should be excluded
///
/// Returns (regions, start_index, end_index)
pub fn number_regions(
    sequence: &[u8],
    state_vector: &[StateVectorEntry],
    state_string: &[u8],  // 128 bytes, 'X' or 'I'
    region_string: &[u8],  // 128 bytes, region labels
    region_index_map: &std::collections::HashMap<u8, usize>,
    rels: &mut Vec<i32>,
    n_regions: usize,
    exclude_deletions: &[usize],
) -> (Vec<Vec<NumberedResidue>>, Option<usize>, Option<usize>) {
    let sv = smooth_insertions(state_vector.to_vec());

    let mut regions: Vec<Vec<NumberedResidue>> = (0..n_regions).map(|_| Vec::new()).collect();
    let mut insertion: i32 = -1;
    let mut previous_state_id: usize = 1;
    let mut _previous_state_type = 'd';
    let mut start_index: Option<usize> = None;
    let mut end_index: Option<usize> = None;
    let mut region: Option<usize> = None;

    for &((state_id, state_type), si) in &sv {
        if state_id == 0 || state_id > 128 {
            continue;
        }

        // Determine region
        if state_type != 'i' || region.is_none() {
            if let Some(&idx) = region_index_map.get(&region_string[state_id - 1]) {
                region = Some(idx);
            }
        }

        let reg = match region {
            Some(r) => r,
            None => continue,
        };

        match state_type {
            'm' => {
                if state_string[state_id - 1] == b'I' {
                    if _previous_state_type != 'd' {
                        insertion += 1;
                    }
                    rels[reg] -= 1;
                } else {
                    insertion = -1;
                }

                let ins_code = if insertion >= 0 && (insertion as usize) < ALPHABET.len() {
                    ALPHABET[insertion as usize].to_string()
                } else {
                    " ".to_string()
                };

                if let Some(seq_i) = si {
                    let aa = sequence[seq_i] as char;
                    let pos = state_id as i32 + rels[reg];
                    regions[reg].push(((pos, ins_code), aa));
                    if start_index.is_none() {
                        start_index = Some(seq_i);
                    }
                    end_index = Some(seq_i);
                }

                previous_state_id = state_id;
                _previous_state_type = 'm';
            }
            'i' => {
                insertion += 1;
                let ins_code = if insertion >= 0 && (insertion as usize) < ALPHABET.len() {
                    ALPHABET[insertion as usize].to_string()
                } else {
                    " ".to_string()
                };

                if let Some(seq_i) = si {
                    let aa = sequence[seq_i] as char;
                    let pos = previous_state_id as i32 + rels[reg];
                    regions[reg].push(((pos, ins_code), aa));
                    if start_index.is_none() {
                        start_index = Some(seq_i);
                    }
                    end_index = Some(seq_i);
                }

                _previous_state_type = 'i';
            }
            'd' => {
                _previous_state_type = 'd';
                if state_string[state_id - 1] == b'I' {
                    rels[reg] -= 1;
                    continue;
                }
                insertion = -1;
                previous_state_id = state_id;
            }
            _ => {}
        }

        // Reset insertion if it gets too high (for CDR regions)
        if insertion >= 25 && exclude_deletions.contains(&reg) {
            insertion = 0;
        }
    }

    (regions, start_index, end_index)
}

/// Fill in gaps for missing positions in the numbering.
/// Equivalent to `gap_missing` in ANARCI schemes.py.
pub fn gap_missing(numbering: &[Vec<NumberedResidue>]) -> Vec<NumberedResidue> {
    let flat: Vec<&NumberedResidue> = numbering.iter().flat_map(|r| r.iter()).collect();
    if flat.is_empty() {
        return Vec::new();
    }

    let mut result: Vec<NumberedResidue> = Vec::new();
    let mut last_pos: i32 = 0;

    for &&((pos, ref ins), aa) in &flat {
        if pos > last_pos + 1 {
            for gap_pos in (last_pos + 1)..pos {
                result.push(((gap_pos, " ".to_string()), '-'));
            }
        }
        result.push(((pos, ins.clone()), aa));
        last_pos = pos;
    }

    result
}

/// Get IMGT CDR annotations with symmetric gapping.
/// Equivalent to `get_imgt_cdr` in ANARCI schemes.py.
pub fn get_imgt_cdr(
    length: usize,
    max_length: usize,
    start: i32,
    end: i32,
) -> Vec<Option<(i32, String)>> {
    let size = std::cmp::max(length, max_length);
    let mut annotations: Vec<Option<(i32, String)>> = vec![None; size];

    if length == 0 {
        return annotations;
    }
    if length == 1 {
        annotations[0] = Some((start, " ".to_string()));
        return annotations;
    }

    let mut front: usize = 0;
    let mut back: isize = -1;

    let az: Vec<&str> = ALPHABET[..ALPHABET.len() - 1].to_vec();
    let za: Vec<&str> = az.iter().rev().copied().collect();

    for i in 0..std::cmp::min(length, max_length) {
        if i % 2 == 1 {
            let idx = (size as isize + back) as usize;
            annotations[idx] = Some((end + back as i32, " ".to_string()));
            back -= 1;
        } else {
            annotations[front] = Some((start + front as i32, " ".to_string()));
            front += 1;
        }
    }

    // Find centre point (None positions)
    let centre_points: Vec<usize> = annotations
        .iter()
        .enumerate()
        .filter(|(_, v)| v.is_none())
        .map(|(i, _)| i)
        .collect();

    if centre_points.is_empty() {
        return annotations;
    }

    let centre_left = annotations[*centre_points.first().unwrap() - 1]
        .as_ref()
        .unwrap()
        .0;
    let centre_right = annotations[*centre_points.last().unwrap() + 1]
        .as_ref()
        .unwrap()
        .0;

    let (frontfactor, backfactor) = if max_length % 2 == 0 {
        (max_length / 2, max_length / 2)
    } else {
        (max_length / 2 + 1, max_length / 2)
    };

    for i in 0..length.saturating_sub(max_length) {
        if i % 2 == 0 {
            let idx = (size as isize + back) as usize;
            let letter_offset = back + backfactor as isize;
            let letter_idx = if letter_offset >= 0 {
                letter_offset as usize
            } else {
                (za.len() as isize + letter_offset) as usize
            };
            if letter_idx < za.len() {
                annotations[idx] = Some((centre_right, za[letter_idx].to_string()));
            }
            back -= 1;
        } else {
            if front >= frontfactor && (front - frontfactor) < az.len() {
                annotations[front] = Some((centre_left, az[front - frontfactor].to_string()));
            }
            front += 1;
        }
    }

    annotations
}

#[cfg(test)]
mod tests {
    use super::get_imgt_cdr;

    #[test]
    fn imgt_cdr3_insertions_match_python_pattern() {
        let annotations = get_imgt_cdr(24, 13, 105, 118);
        let observed: Vec<(i32, String)> = annotations.into_iter().flatten().collect();
        let expected = vec![
            (105, " ".to_string()),
            (106, " ".to_string()),
            (107, " ".to_string()),
            (108, " ".to_string()),
            (109, " ".to_string()),
            (110, " ".to_string()),
            (111, " ".to_string()),
            (111, "A".to_string()),
            (111, "B".to_string()),
            (111, "C".to_string()),
            (111, "D".to_string()),
            (111, "E".to_string()),
            (112, "F".to_string()),
            (112, "E".to_string()),
            (112, "D".to_string()),
            (112, "C".to_string()),
            (112, "B".to_string()),
            (112, "A".to_string()),
            (112, " ".to_string()),
            (113, " ".to_string()),
            (114, " ".to_string()),
            (115, " ".to_string()),
            (116, " ".to_string()),
            (117, " ".to_string()),
        ];

        assert_eq!(observed, expected);
    }
}

/// Get CDR3 annotations for Chothia/Kabat schemes.
/// Equivalent to `get_cdr3_annotations` in ANARCI schemes.py.
pub fn get_cdr3_annotations(length: usize, scheme: &str, chain_type: &str) -> Vec<(i32, String)> {
    let az = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";

    if scheme == "imgt" {
        // Handled by get_imgt_cdr instead
        unreachable!("Use get_imgt_cdr for IMGT CDR3");
    }

    if (scheme == "chothia" || scheme == "kabat") && chain_type == "heavy" {
        let insertions = if length > 10 { length - 10 } else { 0 };
        let ordered: Vec<(i32, String)> = vec![
            (100, " ".into()), (99, " ".into()), (98, " ".into()), (97, " ".into()),
            (96, " ".into()), (95, " ".into()), (101, " ".into()), (102, " ".into()),
            (94, " ".into()), (93, " ".into()),
        ];
        let skip = if 10 > length { 10 - length } else { 0 };
        let mut anns: Vec<(i32, String)> = ordered[skip..].to_vec();
        for i in 0..insertions {
            anns.push((100, az[i..i + 1].to_string()));
        }
        anns.sort_by_key(|a| (a.0, a.1.clone()));
        anns
    } else if (scheme == "chothia" || scheme == "kabat") && chain_type == "light" {
        let insertions = if length > 9 { length - 9 } else { 0 };
        let ordered: Vec<(i32, String)> = vec![
            (95, " ".into()), (94, " ".into()), (93, " ".into()), (92, " ".into()),
            (91, " ".into()), (96, " ".into()), (97, " ".into()), (90, " ".into()),
            (89, " ".into()),
        ];
        let skip = if 9 > length { 9 - length } else { 0 };
        let mut anns: Vec<(i32, String)> = ordered[skip..].to_vec();
        for i in 0..insertions {
            anns.push((95, az[i..i + 1].to_string()));
        }
        anns.sort_by_key(|a| (a.0, a.1.clone()));
        anns
    } else {
        Vec::new()
    }
}
