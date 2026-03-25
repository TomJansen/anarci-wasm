/// Martin (extended Chothia) numbering scheme for heavy and light chains.
/// Direct port of the Python numbering logic from ANARCI schemes.py.

use crate::alignment::*;
use std::collections::HashMap;

use super::chothia;

/// Apply the Martin numbering scheme for heavy chains.
///
/// Regions:
///  0 (FW1)  - insertions on position 8 (differs from Chothia which uses 6)
///  1        - simple mapping
///  2 (CDR1) - insertions on 31 (same as Chothia H CDR1)
///  3        - simple mapping
///  4 (CDR2) - insertions on 52 (same as Chothia H)
///  5 (FW3)  - insertions on 72 (Martin-specific)
///  6 (CDR3) - insertions on 100
///  7        - simple mapping
pub fn number_martin_heavy(
    state_vector: &[StateVectorEntry],
    sequence: &[u8],
) -> (Vec<NumberedResidue>, Option<usize>, Option<usize>) {
    let state_string =
        b"XXXXXXXXXIXXXXXXXXXXXXXXXXXXXXIIIIXXXXXXXXXXXXXXXXXXXXXXXIXIIXXXXXXXXXXXIXXXXXXXXIIIXXXXXXXXXXXXXXXXXXXXXXXXXXXXIIIXXXXXXXXXXXXX";
    let region_string =
        b"11111111112222222222222333333333333333444444444444444455555555555666666666666666666666666666666666666666777777777777788888888888";

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
    let exclude_deletions: Vec<usize> = vec![2, 4, 5, 6];

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

    // Region 0 (FW1): insertions on position 8 (Martin uses 8 instead of Chothia's 6)
    let insertions_0 = regions[0].iter().filter(|r| r.0 .1 != " ").count();
    if insertions_0 > 0 {
        let start = regions[0][0].0 .0;
        let length = regions[0].len();
        let mut annotations: Vec<(i32, String)> = Vec::new();
        for pos in start..9 {
            annotations.push((pos, " ".into()));
        }
        for i in 0..insertions_0 {
            annotations.push((8, ALPHABET[i].to_string()));
        }
        annotations.push((9, " ".into()));
        numbering[0] = (0..length)
            .map(|i| (annotations[i].clone(), regions[0][i].1))
            .collect();
    } else {
        numbering[0] = regions[0].clone();
    }

    // Region 2 (CDR1): insertions on 31 (same as Chothia H CDR1)
    {
        let length = regions[2].len();
        let insertions = if length > 11 { length - 11 } else { 0 };
        if insertions > 0 {
            let mut annotations: Vec<(i32, String)> = Vec::new();
            for pos in 23..32 {
                annotations.push((pos, " ".into()));
            }
            for i in 0..insertions {
                annotations.push((31, ALPHABET[i].to_string()));
            }
            annotations.push((32, " ".into()));
            annotations.push((33, " ".into()));
            numbering[2] = (0..length)
                .map(|i| (annotations[i].clone(), regions[2][i].1))
                .collect();
        } else {
            // When no insertions, take positions from 23 up to fill, then 32, 33
            let mut annotations: Vec<(i32, String)> = Vec::new();
            for pos in 23..32 {
                annotations.push((pos, " ".into()));
            }
            let front_take = if length >= 2 { length - 2 } else { 0 };
            annotations.truncate(front_take);
            let tail: Vec<(i32, String)> = vec![(32, " ".into()), (33, " ".into())];
            let tail_take = length.min(2);
            annotations.extend_from_slice(&tail[..tail_take]);
            numbering[2] = (0..length)
                .map(|i| (annotations[i].clone(), regions[2][i].1))
                .collect();
        }
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

    // Region 5 (FW3): insertions on 72 (Martin-specific)
    // 35 base positions (58..=92), insertions placed on 72
    {
        let length = regions[5].len();
        let insertions = if length > 35 { length - 35 } else { 0 };
        if insertions > 0 {
            let mut annotations: Vec<(i32, String)> = Vec::new();
            for pos in 58..73 {
                annotations.push((pos, " ".into()));
            }
            for i in 0..insertions {
                annotations.push((72, ALPHABET[i].to_string()));
            }
            for pos in 73..93 {
                annotations.push((pos, " ".into()));
            }
            numbering[5] = (0..length)
                .map(|i| (annotations[i].clone(), regions[5][i].1))
                .collect();
        } else {
            // Deletions - let alignment place them.
            // Python source: _numbering[4] = _regions[4] (passes CDR2 through from alignment)
            numbering[4] = regions[4].clone();
        }
    }

    // Region 6 (CDR3): insertions on 100
    {
        let length = regions[6].len();
        if length > 36 {
            return (Vec::new(), start_index, end_index);
        }
        let annotations = get_cdr3_annotations(length, "chothia", "heavy");
        numbering[6] = (0..length)
            .map(|i| (annotations[i].clone(), regions[6][i].1))
            .collect();
    }

    (gap_missing(&numbering), start_index, end_index)
}

/// Apply the Martin numbering scheme for light chains.
///
/// Martin light is identical to Chothia light. The Martin and Chothia specifications
/// for light chains are very similar. Martin is more explicit in the location of
/// indels but unlike the heavy chain these are additions instead of changes to the
/// Chothia scheme. Thus Chothia light is implemented as Martin light.
pub fn number_martin_light(
    state_vector: &[StateVectorEntry],
    sequence: &[u8],
) -> (Vec<NumberedResidue>, Option<usize>, Option<usize>) {
    chothia::number_chothia_light(state_vector, sequence)
}
