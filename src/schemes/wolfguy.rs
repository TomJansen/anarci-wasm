/// Wolfguy numbering scheme implementation.
/// Port of number_wolfguy_heavy and number_wolfguy_light from Python ANARCI schemes.py.
///
/// Wolfguy uses large position numbers (H: 100-400, L: 500-800) and does not
/// use gap_missing — it returns the flat concatenation of regions instead.

use crate::alignment::*;
use std::collections::HashMap;

// ---------------------------------------------------------------------------
// BLOSUM62 scoring matrix
// ---------------------------------------------------------------------------

/// BLOSUM62 substitution score lookup.
/// Returns the BLOSUM62 score for a pair of uppercase amino acid characters.
/// Falls back to -1 for unknown pairs (consistent with 'X' row in BLOSUM62).
fn blosum62(a: u8, b: u8) -> i32 {
    // Store the full BLOSUM62 upper-triangle in a static HashMap built on first use.
    // The Python code stores (X,Y) pairs and falls back to (Y,X) on KeyError;
    // we pre-populate both directions.
    static SCORES: std::sync::LazyLock<HashMap<(u8, u8), i32>> = std::sync::LazyLock::new(|| {
        let entries: &[((u8, u8), i32)] = &[
            ((b'A', b'A'), 4),
            ((b'R', b'R'), 5),
            ((b'N', b'N'), 6),
            ((b'D', b'D'), 6),
            ((b'C', b'C'), 9),
            ((b'Q', b'Q'), 5),
            ((b'E', b'E'), 5),
            ((b'G', b'G'), 6),
            ((b'H', b'H'), 8),
            ((b'I', b'I'), 4),
            ((b'L', b'L'), 4),
            ((b'K', b'K'), 5),
            ((b'M', b'M'), 5),
            ((b'F', b'F'), 6),
            ((b'P', b'P'), 7),
            ((b'S', b'S'), 4),
            ((b'T', b'T'), 5),
            ((b'W', b'W'), 11),
            ((b'Y', b'Y'), 7),
            ((b'V', b'V'), 4),
            ((b'B', b'B'), 4),
            ((b'Z', b'Z'), 4),
            ((b'X', b'X'), -1),
            ((b'A', b'R'), -1),
            ((b'A', b'N'), -2),
            ((b'A', b'D'), -2),
            ((b'A', b'C'), 0),
            ((b'A', b'Q'), -1),
            ((b'A', b'E'), -1),
            ((b'A', b'G'), 0),
            ((b'A', b'H'), -2),
            ((b'A', b'I'), -1),
            ((b'A', b'L'), -1),
            ((b'A', b'K'), -1),
            ((b'A', b'M'), -1),
            ((b'A', b'F'), -2),
            ((b'A', b'P'), -1),
            ((b'A', b'S'), 1),
            ((b'A', b'T'), 0),
            ((b'A', b'W'), -3),
            ((b'A', b'Y'), -2),
            ((b'A', b'V'), 0),
            ((b'R', b'N'), 0),
            ((b'R', b'D'), -2),
            ((b'R', b'C'), -3),
            ((b'R', b'Q'), 1),
            ((b'R', b'E'), 0),
            ((b'R', b'G'), -2),
            ((b'R', b'H'), 0),
            ((b'R', b'I'), -3),
            ((b'R', b'L'), -2),
            ((b'R', b'K'), 2),
            ((b'R', b'M'), -1),
            ((b'R', b'F'), -3),
            ((b'R', b'P'), -2),
            ((b'R', b'S'), -1),
            ((b'R', b'T'), -1),
            ((b'R', b'W'), -3),
            ((b'R', b'Y'), -2),
            ((b'R', b'V'), -3),
            ((b'N', b'D'), 1),
            ((b'N', b'C'), -3),
            ((b'N', b'Q'), 0),
            ((b'N', b'E'), 0),
            ((b'N', b'G'), 0),
            ((b'N', b'H'), 1),
            ((b'N', b'I'), -3),
            ((b'N', b'L'), -3),
            ((b'N', b'K'), 0),
            ((b'N', b'M'), -2),
            ((b'N', b'F'), -3),
            ((b'N', b'P'), -2),
            ((b'N', b'S'), 1),
            ((b'N', b'T'), 0),
            ((b'N', b'W'), -4),
            ((b'N', b'Y'), -2),
            ((b'N', b'V'), -3),
            ((b'D', b'C'), -3),
            ((b'D', b'Q'), 0),
            ((b'D', b'E'), 2),
            ((b'D', b'G'), -1),
            ((b'D', b'H'), -1),
            ((b'D', b'I'), -3),
            ((b'D', b'L'), -4),
            ((b'D', b'K'), -1),
            ((b'D', b'M'), -3),
            ((b'D', b'F'), -3),
            ((b'D', b'P'), -1),
            ((b'D', b'S'), 0),
            ((b'D', b'T'), -1),
            ((b'D', b'W'), -4),
            ((b'D', b'Y'), -3),
            ((b'D', b'V'), -3),
            ((b'C', b'Q'), -3),
            ((b'C', b'E'), -4),
            ((b'C', b'G'), -3),
            ((b'C', b'H'), -3),
            ((b'C', b'I'), -1),
            ((b'C', b'L'), -1),
            ((b'C', b'K'), -3),
            ((b'C', b'M'), -1),
            ((b'C', b'F'), -2),
            ((b'C', b'P'), -3),
            ((b'C', b'S'), -1),
            ((b'C', b'T'), -1),
            ((b'C', b'W'), -2),
            ((b'C', b'Y'), -2),
            ((b'C', b'V'), -1),
            ((b'Q', b'E'), 2),
            ((b'Q', b'G'), -2),
            ((b'Q', b'H'), 0),
            ((b'Q', b'I'), -3),
            ((b'Q', b'L'), -2),
            ((b'Q', b'K'), 1),
            ((b'Q', b'M'), 0),
            ((b'Q', b'F'), -3),
            ((b'Q', b'P'), -1),
            ((b'Q', b'S'), 0),
            ((b'Q', b'T'), -1),
            ((b'Q', b'W'), -2),
            ((b'Q', b'Y'), -1),
            ((b'Q', b'V'), -2),
            ((b'E', b'G'), -2),
            ((b'E', b'H'), 0),
            ((b'E', b'I'), -3),
            ((b'E', b'L'), -3),
            ((b'E', b'K'), 1),
            ((b'E', b'M'), -2),
            ((b'E', b'F'), -3),
            ((b'E', b'P'), -1),
            ((b'E', b'S'), 0),
            ((b'E', b'T'), -1),
            ((b'E', b'W'), -3),
            ((b'E', b'Y'), -2),
            ((b'E', b'V'), -2),
            ((b'G', b'H'), -2),
            ((b'G', b'I'), -4),
            ((b'G', b'L'), -4),
            ((b'G', b'K'), -2),
            ((b'G', b'M'), -3),
            ((b'G', b'F'), -3),
            ((b'G', b'P'), -2),
            ((b'G', b'S'), 0),
            ((b'G', b'T'), -2),
            ((b'G', b'W'), -2),
            ((b'G', b'Y'), -3),
            ((b'G', b'V'), -3),
            ((b'H', b'I'), -3),
            ((b'H', b'L'), -3),
            ((b'H', b'K'), -1),
            ((b'H', b'M'), -2),
            ((b'H', b'F'), -1),
            ((b'H', b'P'), -2),
            ((b'H', b'S'), -1),
            ((b'H', b'T'), -2),
            ((b'H', b'W'), -2),
            ((b'H', b'Y'), 2),
            ((b'H', b'V'), -3),
            ((b'I', b'L'), 2),
            ((b'I', b'K'), -3),
            ((b'I', b'M'), 1),
            ((b'I', b'F'), 0),
            ((b'I', b'P'), -3),
            ((b'I', b'S'), -2),
            ((b'I', b'T'), -1),
            ((b'I', b'W'), -3),
            ((b'I', b'Y'), -1),
            ((b'I', b'V'), 3),
            ((b'L', b'K'), -2),
            ((b'L', b'M'), 2),
            ((b'L', b'F'), 0),
            ((b'L', b'P'), -3),
            ((b'L', b'S'), -2),
            ((b'L', b'T'), -1),
            ((b'L', b'W'), -2),
            ((b'L', b'Y'), -1),
            ((b'L', b'V'), 1),
            ((b'K', b'M'), -1),
            ((b'K', b'F'), -3),
            ((b'K', b'P'), -1),
            ((b'K', b'S'), 0),
            ((b'K', b'T'), -1),
            ((b'K', b'W'), -3),
            ((b'K', b'Y'), -2),
            ((b'K', b'V'), -2),
            ((b'M', b'F'), 0),
            ((b'M', b'P'), -2),
            ((b'M', b'S'), -1),
            ((b'M', b'T'), -1),
            ((b'M', b'W'), -1),
            ((b'M', b'Y'), -1),
            ((b'M', b'V'), 1),
            ((b'F', b'P'), -4),
            ((b'F', b'S'), -2),
            ((b'F', b'T'), -2),
            ((b'F', b'W'), 1),
            ((b'F', b'Y'), 3),
            ((b'F', b'V'), -1),
            ((b'P', b'S'), -1),
            ((b'P', b'T'), -1),
            ((b'P', b'W'), -4),
            ((b'P', b'Y'), -3),
            ((b'P', b'V'), -2),
            ((b'S', b'T'), 1),
            ((b'S', b'W'), -3),
            ((b'S', b'Y'), -2),
            ((b'S', b'V'), -2),
            ((b'T', b'W'), -2),
            ((b'T', b'Y'), -2),
            ((b'T', b'V'), 0),
            ((b'W', b'Y'), 2),
            ((b'W', b'V'), -3),
            ((b'Y', b'V'), -1),
            // B and Z ambiguity codes
            ((b'B', b'N'), 3),
            ((b'B', b'D'), 4),
            ((b'B', b'A'), -2),
            ((b'B', b'C'), -3),
            ((b'B', b'E'), 1),
            ((b'B', b'F'), -3),
            ((b'B', b'G'), -1),
            ((b'B', b'H'), 0),
            ((b'B', b'I'), -3),
            ((b'B', b'K'), 0),
            ((b'B', b'L'), -4),
            ((b'B', b'M'), -3),
            ((b'B', b'P'), -2),
            ((b'B', b'Q'), 0),
            ((b'B', b'R'), -1),
            ((b'B', b'S'), 0),
            ((b'B', b'T'), -1),
            ((b'B', b'V'), -3),
            ((b'B', b'W'), -4),
            ((b'B', b'Y'), -3),
            ((b'Z', b'A'), -1),
            ((b'Z', b'B'), 1),
            ((b'Z', b'C'), -3),
            ((b'Z', b'D'), 1),
            ((b'Z', b'E'), 4),
            ((b'Z', b'F'), -3),
            ((b'Z', b'G'), -2),
            ((b'Z', b'H'), 0),
            ((b'Z', b'I'), -3),
            ((b'Z', b'K'), 1),
            ((b'Z', b'L'), -3),
            ((b'Z', b'M'), -1),
            ((b'Z', b'N'), 0),
            ((b'Z', b'P'), -1),
            ((b'Z', b'Q'), 3),
            ((b'Z', b'R'), 0),
            ((b'Z', b'S'), 0),
            ((b'Z', b'T'), -1),
            ((b'Z', b'V'), -2),
            ((b'Z', b'W'), -3),
            ((b'Z', b'Y'), -2),
            // X wildcard
            ((b'X', b'A'), 0),
            ((b'X', b'B'), -1),
            ((b'X', b'C'), -2),
            ((b'X', b'D'), -1),
            ((b'X', b'E'), -1),
            ((b'X', b'F'), -1),
            ((b'X', b'G'), -1),
            ((b'X', b'H'), -1),
            ((b'X', b'I'), -1),
            ((b'X', b'K'), -1),
            ((b'X', b'L'), -1),
            ((b'X', b'M'), -1),
            ((b'X', b'N'), -1),
            ((b'X', b'P'), -2),
            ((b'X', b'Q'), -1),
            ((b'X', b'R'), -1),
            ((b'X', b'S'), 0),
            ((b'X', b'T'), 0),
            ((b'X', b'V'), -1),
            ((b'X', b'W'), -2),
            ((b'X', b'Y'), -1),
            ((b'X', b'Z'), -1),
        ];

        let mut map = HashMap::with_capacity(entries.len() * 2);
        for &((a, b), score) in entries {
            map.insert((a, b), score);
            if a != b {
                map.insert((b, a), score);
            }
        }
        map
    });

    let au = a.to_ascii_uppercase();
    let bu = b.to_ascii_uppercase();
    SCORES.get(&(au, bu)).copied().unwrap_or(-1)
}

// ---------------------------------------------------------------------------
// Wolfguy Heavy
// ---------------------------------------------------------------------------

/// Apply the Wolfguy numbering scheme for heavy chains.
///
/// Uses large position numbers (100-400 range). CDR regions are renumbered
/// with symmetric gapping patterns. Returns flat concatenation (no gap_missing).
pub fn number_wolfguy_heavy(
    state_vector: &[StateVectorEntry],
    sequence: &[u8],
) -> (Vec<NumberedResidue>, Option<usize>, Option<usize>) {
    let state_string = b"XXXXXXXXXIXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXIXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX";

    let region_string = b"11111111111111111111111111222222222222233333333333333344444444444444444444555555555555555555555555555555666666666666677777777777";

    let mut region_index_map: HashMap<u8, usize> = HashMap::new();
    region_index_map.insert(b'1', 0);
    region_index_map.insert(b'2', 1);
    region_index_map.insert(b'3', 2);
    region_index_map.insert(b'4', 3);
    region_index_map.insert(b'5', 4);
    region_index_map.insert(b'6', 5);
    region_index_map.insert(b'7', 6);

    let mut rels: Vec<i32> = vec![100, 124, 160, 196, 226, 244, 283];

    let n_regions = 7;
    let exclude_deletions: Vec<usize> = vec![1, 3, 5];

    let (regions, startindex, endindex) = number_regions(
        sequence,
        state_vector,
        state_string,
        region_string,
        &region_index_map,
        &mut rels,
        n_regions,
        &exclude_deletions,
    );

    let mut numbering: Vec<Vec<NumberedResidue>> = vec![
        regions[0].clone(), // FW1
        Vec::new(),         // CDRH1 - renumbered
        regions[2].clone(), // FW2
        Vec::new(),         // CDRH2 - renumbered
        regions[4].clone(), // FW3
        Vec::new(),         // CDRH3 - renumbered
        regions[6].clone(), // FW4
    ];

    // ========================
    // CDRH1 - region 1
    // ========================
    // Symmetric about 177. Delete right first.
    // ordered_deletions = [151] + interleaved pairs (152,199),(153,198),...,(175,176)
    let mut ordered_deletions: Vec<i32> = vec![151];
    for (p1, p2) in (152..176).zip((176..200).rev()) {
        ordered_deletions.push(p1);
        ordered_deletions.push(p2);
    }

    let length = regions[1].len();
    let mut annotations: Vec<i32> = ordered_deletions[..length.min(ordered_deletions.len())].to_vec();
    annotations.sort();
    numbering[1] = (0..length)
        .map(|i| ((annotations[i], " ".to_string()), regions[1][i].1))
        .collect();

    // ========================
    // CDRH2 - region 3
    // ========================
    // Symmetric about 271. Start with 251, interleave pairs (252,290),...,(270,272), then 271.
    // Prepend range(299,290,-1) = [299,298,...,291]
    let mut cdr2_inner: Vec<i32> = vec![251];
    for (p1, p2) in (252..271).zip((272..291).rev()) {
        cdr2_inner.push(p1);
        cdr2_inner.push(p2);
    }
    cdr2_inner.push(271);
    let prefix: Vec<i32> = (291..300).rev().collect(); // 299,298,...,291
    let mut ordered_deletions: Vec<i32> = Vec::new();
    ordered_deletions.extend_from_slice(&prefix);
    ordered_deletions.extend_from_slice(&cdr2_inner);

    let length = regions[3].len();
    let mut annotations: Vec<i32> = ordered_deletions[..length.min(ordered_deletions.len())].to_vec();
    annotations.sort();
    numbering[3] = (0..length)
        .map(|i| ((annotations[i], " ".to_string()), regions[3][i].1))
        .collect();

    // ========================
    // CDRH3 - region 5
    // ========================
    // Symmetric about 374. Pairs (356,391),(357,390),...,(373,375)
    // Prepend [354,394,355,393,392]
    // Prepend [331,332,399,398,351,352,397,353,396,395]
    let mut cdr3_inner: Vec<i32> = Vec::new();
    for (p1, p2) in (356..374).zip((374..392).rev()) {
        cdr3_inner.push(p1);
        cdr3_inner.push(p2);
    }
    let mid_prefix = vec![354, 394, 355, 393, 392];
    let outer_prefix = vec![331, 332, 399, 398, 351, 352, 397, 353, 396, 395];

    let mut ordered_deletions: Vec<i32> = Vec::new();
    ordered_deletions.extend_from_slice(&outer_prefix);
    ordered_deletions.extend_from_slice(&mid_prefix);
    ordered_deletions.extend_from_slice(&cdr3_inner);

    let length = regions[5].len();
    if length > ordered_deletions.len() {
        return (Vec::new(), startindex, endindex);
    }
    let mut annotations: Vec<i32> = ordered_deletions[..length].to_vec();
    annotations.sort();
    numbering[5] = (0..length)
        .map(|i| ((annotations[i], " ".to_string()), regions[5][i].1))
        .collect();

    // Flatten all regions (no gap_missing for wolfguy)
    let result: Vec<NumberedResidue> = numbering.into_iter().flatten().collect();
    (result, startindex, endindex)
}

// ---------------------------------------------------------------------------
// Wolfguy Light
// ---------------------------------------------------------------------------

/// Apply the Wolfguy numbering scheme for light chains.
///
/// Uses large position numbers (500-800 range). CDR regions are renumbered
/// with scheme-specific gapping and L1 canonical recognition.
/// Returns flat concatenation (no gap_missing).
pub fn number_wolfguy_light(
    state_vector: &[StateVectorEntry],
    sequence: &[u8],
) -> (Vec<NumberedResidue>, Option<usize>, Option<usize>) {
    let state_string = b"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXIXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX";

    let region_string = b"1111111AAABBBBBBBBBBBBB222222222222222223333333333333334444444444444455555555555666677777777777777777777888888888888899999999999";

    let mut region_index_map: HashMap<u8, usize> = HashMap::new();
    region_index_map.insert(b'1', 0);
    region_index_map.insert(b'A', 1);
    region_index_map.insert(b'B', 2);
    region_index_map.insert(b'2', 3);
    region_index_map.insert(b'3', 4);
    region_index_map.insert(b'4', 5);
    region_index_map.insert(b'5', 6);
    region_index_map.insert(b'6', 7);
    region_index_map.insert(b'7', 8);
    region_index_map.insert(b'8', 9);
    region_index_map.insert(b'9', 10);

    let mut rels: Vec<i32> = vec![500, 500, 500, 527, 560, 595, 631, 630, 630, 646, 683];

    let n_regions = 11;
    let exclude_deletions: Vec<usize> = vec![1, 3, 5, 7, 9];

    let (regions, startindex, endindex) = number_regions(
        sequence,
        state_vector,
        state_string,
        region_string,
        &region_index_map,
        &mut rels,
        n_regions,
        &exclude_deletions,
    );

    let mut numbering: Vec<Vec<NumberedResidue>> = vec![
        regions[0].clone(),  // 1
        Vec::new(),          // A - renumbered (FW gaps on 508)
        regions[2].clone(),  // B
        Vec::new(),          // 2 (CDRL1) - renumbered
        regions[4].clone(),  // 3
        Vec::new(),          // 4 (CDRL2) - renumbered
        regions[6].clone(),  // 5
        Vec::new(),          // 6 - renumbered (indel on 713/714)
        regions[8].clone(),  // 7
        Vec::new(),          // 8 (CDRL3) - renumbered
        regions[10].clone(), // 9
    ];

    // ========================
    // Region A (index 1) - FW gaps on 508
    // ========================
    let length = regions[1].len();
    {
        let mut base: Vec<(i32, String)> = vec![
            (510, " ".to_string()),
            (509, " ".to_string()),
            (508, " ".to_string()),
        ];
        base.truncate(length);
        // Add insertions on 508
        let extra = if length > 3 { length - 3 } else { 0 };
        for a in 0..extra {
            base.push((508, ALPHABET[a].to_string()));
        }
        base.sort_by(|a, b| a.0.cmp(&b.0).then_with(|| a.1.cmp(&b.1)));
        numbering[1] = (0..length)
            .map(|i| (base[i].clone(), regions[1][i].1))
            .collect();
    }

    // ========================
    // CDRL1 - Region 2 (index 3)
    // ========================
    let length = regions[3].len();
    {
        let annotations = get_wolfguy_l1(&regions[3], length);
        numbering[3] = (0..length)
            .map(|i| ((annotations[i], " ".to_string()), regions[3][i].1))
            .collect();
    }

    // ========================
    // CDRL2 - Region 4 (index 5)
    // ========================
    // Symmetric about 673, then right from 694. Maintain 651 as last deletion.
    let length = regions[5].len();
    {
        let mut ordered_deletions: Vec<i32> = Vec::new();
        for (p1, p2) in (652..673).zip((673..695).rev()) {
            ordered_deletions.push(p2);
            ordered_deletions.push(p1);
        }
        // Prepend: [651] + [699,698,697,696,695] + ordered_deletions + [673]
        let mut full: Vec<i32> = vec![651];
        full.extend((695..700).rev()); // 699,698,697,696,695
        full.extend_from_slice(&ordered_deletions);
        full.push(673);

        let mut annotations: Vec<i32> = full[..length.min(full.len())].to_vec();
        annotations.sort();
        numbering[5] = (0..length)
            .map(|i| ((annotations[i], " ".to_string()), regions[5][i].1))
            .collect();
    }

    // ========================
    // Region 6 (index 7) - indel on 713,714
    // ========================
    let length = regions[7].len();
    {
        let mut annotations: Vec<(i32, String)> = vec![
            (711, " ".to_string()),
            (712, " ".to_string()),
            (713, " ".to_string()),
            (714, " ".to_string()),
        ];
        annotations.truncate(length.min(4));
        let insertions = if length > 4 { length - 4 } else { 0 };
        for a in 0..insertions {
            annotations.push((714, ALPHABET[a].to_string()));
        }
        numbering[7] = (0..length)
            .map(|i| (annotations[i].clone(), regions[7][i].1))
            .collect();
    }

    // ========================
    // CDRL3 - Region 8 (index 9)
    // ========================
    // Symmetric about 775. Delete right first.
    let length = regions[9].len();
    {
        let mut ordered_deletions: Vec<i32> = Vec::new();
        for (p1, p2) in (751..775).zip((776..800).rev()) {
            ordered_deletions.push(p1);
            ordered_deletions.push(p2);
        }
        ordered_deletions.push(775);

        if length > ordered_deletions.len() {
            return (Vec::new(), startindex, endindex);
        }
        let mut annotations: Vec<i32> = ordered_deletions[..length].to_vec();
        annotations.sort();
        numbering[9] = (0..length)
            .map(|i| ((annotations[i], " ".to_string()), regions[9][i].1))
            .collect();
    }

    // Flatten all regions (no gap_missing for wolfguy)
    let result: Vec<NumberedResidue> = numbering.into_iter().flatten().collect();
    (result, startindex, endindex)
}

// ---------------------------------------------------------------------------
// Wolfguy L1 canonical recognition
// ---------------------------------------------------------------------------

/// Determine Wolfguy L1 annotation positions by matching against canonical
/// sequence patterns with BLOSUM62 scoring. For unknown lengths, uses
/// symmetric numbering about position 575.
fn get_wolfguy_l1(seq: &[NumberedResidue], length: usize) -> Vec<i32> {
    // Pre-defined canonical patterns for each CDR L1 length.
    // Each entry: (name, consensus_sequence, positions)
    struct Canonical {
        _name: &'static str,
        consensus: &'static [u8],
        positions: &'static [i32],
    }

    let l1_9: &[Canonical] = &[Canonical {
        _name: "9",
        consensus: b"XXXXXXXXX",
        positions: &[551, 552, 554, 556, 563, 572, 597, 598, 599],
    }];

    let l1_10: &[Canonical] = &[Canonical {
        _name: "10",
        consensus: b"XXXXXXXXXX",
        positions: &[551, 552, 553, 556, 561, 562, 571, 597, 598, 599],
    }];

    let l1_11: &[Canonical] = &[
        Canonical {
            _name: "11a",
            consensus: b"RASQDISSYLA",
            positions: &[551, 552, 553, 556, 561, 562, 571, 596, 597, 598, 599],
        },
        Canonical {
            _name: "11b",
            consensus: b"GGNNIGSKSVH",
            positions: &[551, 552, 554, 556, 561, 562, 571, 572, 597, 598, 599],
        },
        Canonical {
            _name: "11b.2",
            consensus: b"SGDQLPKKYAY",
            positions: &[551, 552, 554, 556, 561, 562, 571, 572, 597, 598, 599],
        },
    ];

    let l1_12: &[Canonical] = &[
        Canonical {
            _name: "12a",
            consensus: b"TLSSQHSTYTIE",
            positions: &[551, 552, 553, 554, 555, 556, 561, 563, 572, 597, 598, 599],
        },
        Canonical {
            _name: "12b",
            consensus: b"TASSSVSSSYLH",
            positions: &[551, 552, 553, 556, 561, 562, 571, 595, 596, 597, 598, 599],
        },
        Canonical {
            _name: "12c",
            consensus: b"RASQSVXNNYLA",
            positions: &[551, 552, 553, 556, 561, 562, 571, 581, 596, 597, 598, 599],
        },
        Canonical {
            _name: "12d",
            consensus: b"RSSHSIRRSRVH",
            positions: &[551, 552, 553, 556, 561, 562, 571, 581, 596, 597, 598, 599],
        },
    ];

    let l1_13: &[Canonical] = &[
        Canonical {
            _name: "13a",
            consensus: b"SGSSSNIGNNYVS",
            positions: &[551, 552, 554, 555, 556, 557, 561, 562, 571, 572, 597, 598, 599],
        },
        Canonical {
            _name: "13b",
            consensus: b"TRSSGSLANYYVQ",
            positions: &[551, 552, 553, 554, 556, 561, 562, 563, 571, 572, 597, 598, 599],
        },
    ];

    let l1_14: &[Canonical] = &[
        Canonical {
            _name: "14a",
            consensus: b"RSSTGAVTTSNYAN",
            positions: &[551, 552, 553, 554, 555, 561, 562, 563, 564, 571, 572, 597, 598, 599],
        },
        Canonical {
            _name: "14b",
            consensus: b"TGTSSDVGGYNYVS",
            positions: &[551, 552, 554, 555, 556, 557, 561, 562, 571, 572, 596, 597, 598, 599],
        },
    ];

    let l1_15: &[Canonical] = &[Canonical {
        _name: "15",
        consensus: b"XXXXXXXXXXXXXXX",
        positions: &[551, 552, 553, 556, 561, 562, 563, 581, 582, 594, 595, 596, 597, 598, 599],
    }];

    let l1_16: &[Canonical] = &[Canonical {
        _name: "16",
        consensus: b"XXXXXXXXXXXXXXXX",
        positions: &[
            551, 552, 553, 556, 561, 562, 563, 581, 582, 583, 594, 595, 596, 597, 598, 599,
        ],
    }];

    let l1_17: &[Canonical] = &[Canonical {
        _name: "17",
        consensus: b"XXXXXXXXXXXXXXXXX",
        positions: &[
            551, 552, 553, 556, 561, 562, 563, 581, 582, 583, 584, 594, 595, 596, 597, 598, 599,
        ],
    }];

    // Look up the canonical list for this length
    let canonicals: Option<&[Canonical]> = match length {
        9 => Some(l1_9),
        10 => Some(l1_10),
        11 => Some(l1_11),
        12 => Some(l1_12),
        13 => Some(l1_13),
        14 => Some(l1_14),
        15 => Some(l1_15),
        16 => Some(l1_16),
        17 => Some(l1_17),
        _ => None,
    };

    if let Some(canonicals) = canonicals {
        // Score each canonical and pick the best
        let mut best_positions: &[i32] = canonicals[0].positions;
        let mut best_score: i32 = i32::MIN;

        for canonical in canonicals {
            let mut sub_score: i32 = 0;
            for i in 0..length {
                let seq_aa = (seq[i].1 as u8).to_ascii_uppercase();
                let con_aa = canonical.consensus[i].to_ascii_uppercase();
                sub_score += blosum62(seq_aa, con_aa);
            }
            if sub_score > best_score {
                best_score = sub_score;
                best_positions = canonical.positions;
            }
        }

        best_positions.to_vec()
    } else {
        // Symmetric numbering about 575 for unknown lengths
        let mut ordered_deletions: Vec<i32> = Vec::new();
        for (p1, p2) in (551..575).zip((576..600).rev()) {
            ordered_deletions.push(p2);
            ordered_deletions.push(p1);
        }
        ordered_deletions.push(575);

        let mut result: Vec<i32> = ordered_deletions[..length.min(ordered_deletions.len())].to_vec();
        result.sort();
        result
    }
}
