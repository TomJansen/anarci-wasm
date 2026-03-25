/// Kabat numbering scheme for heavy and light chains.
/// Direct port of the Python numbering logic from ANARCI schemes.py.

use crate::alignment::*;
use std::collections::HashMap;

/// Apply the Kabat numbering scheme for heavy chains.
///
/// Regions:
///  0 (FW1)  - insertions on position 6
///  1        - simple mapping
///  2 (CDR1) - insertions on 35, delete forward from 35
///  3        - simple mapping
///  4 (CDR2) - insertions on 52
///  5        - simple mapping
///  6 (CDR3) - insertions on 100
///  7        - simple mapping
pub fn number_kabat_heavy(
    state_vector: &[StateVectorEntry],
    sequence: &[u8],
) -> (Vec<NumberedResidue>, Option<usize>, Option<usize>) {
    let state_string =
        b"XXXXXXXXXIXXXXXXXXXXXXXXXXXXXXIIIIXXXXXXXXXXXXXXXXXXXXXXXIXIIXXXXXXXXXXXIXXXXXXXXXXXXXXXXXXIIIXXXXXXXXXXXXXXXXXXIIIXXXXXXXXXXXXX";
    let region_string =
        b"11111111112222222222222333333333333333334444444444444455555555555666666666666666666666666666666666666666777777777777788888888888";

    let region_index_map: HashMap<u8, usize> = [
        (b'1', 0),
        (b'2', 1),
        (b'3', 2),
        (b'4', 3),
        (b'5', 4),
        (b'6', 5),
        (b'7', 6),
        (b'8', 7),
    ]
    .iter()
    .cloned()
    .collect();

    let mut rels: Vec<i32> = vec![0, -1, -1, -5, -5, -8, -12, -15];
    let n_regions = 8;
    let exclude_deletions: Vec<usize> = vec![2, 4, 6];

    let (regions, start_index, end_index) = number_regions(
        sequence,
        state_vector,
        state_string,
        region_string,
        &region_index_map,
        &mut rels,
        n_regions,
        &exclude_deletions,
    );

    // Build numbering: regions 1, 3, 5, 7 pass through; 0, 2, 4, 6 are renumbered
    let mut numbering: Vec<Vec<NumberedResidue>> = vec![
        Vec::new(),
        regions[1].clone(),
        Vec::new(),
        regions[3].clone(),
        Vec::new(),
        regions[5].clone(),
        Vec::new(),
        regions[7].clone(),
    ];

    // Region 0 (FW1): insertions on position 6
    let insertions_0 = regions[0].iter().filter(|r| r.0 .1 != " ").count();
    if insertions_0 > 0 {
        let start = regions[0][0].0 .0;
        let length = regions[0].len();
        let mut annotations: Vec<(i32, String)> = Vec::new();
        for pos in start..7 {
            annotations.push((pos, " ".into()));
        }
        for i in 0..insertions_0 {
            annotations.push((6, ALPHABET[i].to_string()));
        }
        annotations.push((7, " ".into()));
        annotations.push((8, " ".into()));
        annotations.push((9, " ".into()));
        numbering[0] = (0..length)
            .map(|i| (annotations[i].clone(), regions[0][i].1))
            .collect();
    } else {
        numbering[0] = regions[0].clone();
    }

    // Region 2 (CDR1): insertions on 35, delete forward from 35
    {
        let length = regions[2].len();
        let insertions = if length > 13 { length - 13 } else { 0 };
        // Base positions 23..=35 (13 positions), then insertions on 35
        let mut annotations: Vec<(i32, String)> = Vec::new();
        for pos in 23..36 {
            annotations.push((pos, " ".into()));
            if annotations.len() >= length {
                break;
            }
        }
        // Truncate to length if fewer than 13
        annotations.truncate(length);
        for i in 0..insertions {
            annotations.push((35, ALPHABET[i].to_string()));
        }
        numbering[2] = (0..length)
            .map(|i| (annotations[i].clone(), regions[2][i].1))
            .collect();
    }

    // Region 4 (CDR2): insertions on 52 (same as Chothia H)
    {
        let length = regions[4].len();
        let insertions = if length > 8 { length - 8 } else { 0 };
        // Delete order: 52, 51, 50, 53, 54, 55, 56, 57
        let front: Vec<(i32, String)> = vec![
            (50, " ".into()),
            (51, " ".into()),
            (52, " ".into()),
        ];
        let back: Vec<(i32, String)> = vec![
            (53, " ".into()),
            (54, " ".into()),
            (55, " ".into()),
            (56, " ".into()),
            (57, " ".into()),
        ];
        let front_take = if length > 5 { length - 5 } else { 0 }.min(3);
        let back_skip = if length >= 5 { 0 } else { 5 - length };
        let mut annotations: Vec<(i32, String)> = Vec::new();
        annotations.extend_from_slice(&front[..front_take]);
        for i in 0..insertions {
            annotations.push((52, ALPHABET[i].to_string()));
        }
        annotations.extend_from_slice(&back[back_skip..]);
        numbering[4] = (0..length)
            .map(|i| (annotations[i].clone(), regions[4][i].1))
            .collect();
    }

    // Region 6 (CDR3): insertions on 100
    {
        let length = regions[6].len();
        if length > 36 {
            return (Vec::new(), start_index, end_index);
        }
        let annotations = get_cdr3_annotations(length, "kabat", "heavy");
        numbering[6] = (0..length)
            .map(|i| (annotations[i].clone(), regions[6][i].1))
            .collect();
    }

    (gap_missing(&numbering), start_index, end_index)
}

/// Apply the Kabat numbering scheme for light chains.
///
/// Regions:
///  0        - simple mapping
///  1 (CDR1) - insertions on 27
///  2        - simple mapping
///  3 (CDR2) - insertions on 52
///  4        - simple mapping
///  5 (CDR3) - insertions on 95
///  6        - simple mapping
pub fn number_kabat_light(
    state_vector: &[StateVectorEntry],
    sequence: &[u8],
) -> (Vec<NumberedResidue>, Option<usize>, Option<usize>) {
    let state_string =
        b"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXIIIIIIXXXXXXXXXXXXXXXXXXXXXXIIIIIIIXXXXXXXXIXXXXXXXIIXXXXXXXXXXXXXXXXXXXXXXXXXXXIIIIXXXXXXXXXXXXXXX";
    let region_string =
        b"11111111111111111111111222222222222222223333333333333333444444444445555555555555555555555555555555555555666666666666677777777777";

    let region_index_map: HashMap<u8, usize> = [
        (b'1', 0),
        (b'2', 1),
        (b'3', 2),
        (b'4', 3),
        (b'5', 4),
        (b'6', 5),
        (b'7', 6),
    ]
    .iter()
    .cloned()
    .collect();

    let mut rels: Vec<i32> = vec![0, 0, -6, -6, -13, -16, -20];
    let n_regions = 7;
    let exclude_deletions: Vec<usize> = vec![1, 3, 5];

    let (regions, start_index, end_index) = number_regions(
        sequence,
        state_vector,
        state_string,
        region_string,
        &region_index_map,
        &mut rels,
        n_regions,
        &exclude_deletions,
    );

    // Build numbering: regions 0, 2, 4, 6 pass through; 1, 3, 5 are renumbered
    let mut numbering: Vec<Vec<NumberedResidue>> = vec![
        regions[0].clone(),
        Vec::new(),
        regions[2].clone(),
        Vec::new(),
        regions[4].clone(),
        Vec::new(),
        regions[6].clone(),
    ];

    // Region 1 (CDR1): insertions on 27
    {
        let length = regions[1].len();
        let insertions = if length > 11 { length - 11 } else { 0 };
        // 4 positions before insert point: 24, 25, 26, 27
        let front: Vec<(i32, String)> = vec![
            (24, " ".into()),
            (25, " ".into()),
            (26, " ".into()),
            (27, " ".into()),
        ];
        let back: Vec<(i32, String)> = vec![
            (28, " ".into()),
            (29, " ".into()),
            (30, " ".into()),
            (31, " ".into()),
            (32, " ".into()),
            (33, " ".into()),
            (34, " ".into()),
        ];
        let front_take = length.min(4);
        let back_skip = if length >= 11 { 0 } else { (11 - length).min(7) };
        let mut annotations: Vec<(i32, String)> = Vec::new();
        annotations.extend_from_slice(&front[..front_take]);
        for i in 0..insertions {
            annotations.push((27, ALPHABET[i].to_string()));
        }
        annotations.extend_from_slice(&back[back_skip..]);
        numbering[1] = (0..length)
            .map(|i| (annotations[i].clone(), regions[1][i].1))
            .collect();
    }

    // Region 3 (CDR2): insertions on 52 (same as Chothia L)
    {
        let length = regions[3].len();
        let insertions = if length > 4 { length - 4 } else { 0 };
        if insertions > 0 {
            let mut annotations: Vec<(i32, String)> = vec![
                (51, " ".into()),
                (52, " ".into()),
            ];
            for i in 0..insertions {
                annotations.push((52, ALPHABET[i].to_string()));
            }
            annotations.push((53, " ".into()));
            annotations.push((54, " ".into()));
            numbering[3] = (0..length)
                .map(|i| (annotations[i].clone(), regions[3][i].1))
                .collect();
        } else {
            // Let the alignment place gaps
            numbering[3] = regions[3].clone();
        }
    }

    // Region 5 (CDR3): insertions on 95
    {
        let length = regions[5].len();
        if length > 35 {
            return (Vec::new(), start_index, end_index);
        }
        let annotations = get_cdr3_annotations(length, "kabat", "light");
        numbering[5] = (0..length)
            .map(|i| (annotations[i].clone(), regions[5][i].1))
            .collect();
    }

    (gap_missing(&numbering), start_index, end_index)
}
