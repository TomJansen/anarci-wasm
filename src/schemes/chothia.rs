/// Chothia numbering scheme implementation for heavy and light chains.
/// Port of number_chothia_heavy and number_chothia_light from Python ANARCI schemes.py.

use crate::alignment::*;
use std::collections::HashMap;

/// Apply the Chothia numbering scheme for heavy chains.
///
/// 8 regions are defined:
///  0 - FW1: insertions placed at Chothia position 6
///  1 - Simple mapping
///  2 - CDRH1: 23-33, insertions on 31
///  3 - Simple mapping
///  4 - CDRH2: 50-57, insertions on 52
///  5 - Simple mapping
///  6 - CDRH3: insertions on 100, uses get_cdr3_annotations
///  7 - Simple mapping (FW4)
pub fn number_chothia_heavy(
    state_vector: &[StateVectorEntry],
    sequence: &[u8],
) -> (Vec<NumberedResidue>, Option<usize>, Option<usize>) {
    let state_string = b"XXXXXXXXXIXXXXXXXXXXXXXXXXXXXXIIIIXXXXXXXXXXXXXXXXXXXXXXXIXIIXXXXXXXXXXXIXXXXXXXXXXXXXXXXXXIIIXXXXXXXXXXXXXXXXXXIIIXXXXXXXXXXXXX";

    let region_string = b"11111111112222222222222333333333333333444444444444444455555555555666666666666666666666666666666666666666777777777777788888888888";

    let mut region_index_map: HashMap<u8, usize> = HashMap::new();
    region_index_map.insert(b'1', 0);
    region_index_map.insert(b'2', 1);
    region_index_map.insert(b'3', 2);
    region_index_map.insert(b'4', 3);
    region_index_map.insert(b'5', 4);
    region_index_map.insert(b'6', 5);
    region_index_map.insert(b'7', 6);
    region_index_map.insert(b'8', 7);

    let mut rels: Vec<i32> = vec![0, -1, -1, -5, -5, -8, -12, -15];

    let n_regions = 8;
    let exclude_deletions: Vec<usize> = vec![0, 2, 4, 6];

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

    // Build renumbered output: regions 0, 2, 4, 6 are renumbered; 1, 3, 5, 7 pass through
    let mut numbering: Vec<Vec<NumberedResidue>> = vec![
        Vec::new(),         // region 0: FW1 (renumbered)
        regions[1].clone(), // region 1: pass through
        Vec::new(),         // region 2: CDR1 (renumbered)
        regions[3].clone(), // region 3: pass through
        Vec::new(),         // region 4: CDR2 (renumbered)
        regions[5].clone(), // region 5: pass through
        Vec::new(),         // region 6: CDR3 (renumbered)
        regions[7].clone(), // region 7: FW4 pass through
    ];

    // Region 0: FW1 - insertions placed at Chothia position 6
    let insertions_r0 = regions[0].iter().filter(|r| r.0 .1 != " ").count();
    if insertions_r0 > 0 {
        let start = regions[0][0].0 .0; // Starting Chothia number
        let length = regions[0].len();
        let mut annotations: Vec<(i32, String)> = Vec::new();
        for pos in start..7 {
            annotations.push((pos, " ".to_string()));
        }
        for i in 0..insertions_r0 {
            annotations.push((6, ALPHABET[i].to_string()));
        }
        annotations.push((7, " ".to_string()));
        annotations.push((8, " ".to_string()));
        annotations.push((9, " ".to_string()));
        numbering[0] = (0..length)
            .map(|i| (annotations[i].clone(), regions[0][i].1))
            .collect();
    } else {
        numbering[0] = regions[0].clone();
    }

    // Region 2: CDR1 - range 23-33 (11 base positions), insertions on 31
    let length_r2 = regions[2].len();
    let insertions_r2 = if length_r2 > 11 { length_r2 - 11 } else { 0 };
    let annotations_r2: Vec<(i32, String)> = if insertions_r2 > 0 {
        let mut anns: Vec<(i32, String)> = (23..32).map(|p| (p, " ".to_string())).collect();
        for i in 0..insertions_r2 {
            anns.push((31, ALPHABET[i].to_string()));
        }
        anns.push((32, " ".to_string()));
        anns.push((33, " ".to_string()));
        anns
    } else {
        // No insertions: take positions 23..32 truncated, then 32,33
        let base_count = if length_r2 >= 2 { length_r2 - 2 } else { 0 };
        let mut anns: Vec<(i32, String)> = (23..32)
            .take(base_count)
            .map(|p| (p, " ".to_string()))
            .collect();
        let tail: Vec<(i32, String)> = vec![(32, " ".to_string()), (33, " ".to_string())];
        let tail_take = length_r2.min(tail.len() + anns.len()) - anns.len();
        anns.extend_from_slice(&tail[..tail_take]);
        anns
    };
    numbering[2] = (0..length_r2)
        .map(|i| (annotations_r2[i].clone(), regions[2][i].1))
        .collect();

    // Region 4: CDR2 - range 50-57 (8 base positions), insertions on 52
    // Delete order: 52, 51, 50, 53, 54, 55, 56, 57
    let length_r4 = regions[4].len();
    let insertions_r4 = if length_r4 > 8 { length_r4 - 8 } else { 0 };
    let annotations_r4: Vec<(i32, String)> = {
        // Front positions: (50, 51, 52) up to max(0, length-5) items
        let front_count = if length_r4 > 5 { length_r4 - 5 } else { 0 };
        let front_positions = vec![
            (50, " ".to_string()),
            (51, " ".to_string()),
            (52, " ".to_string()),
        ];
        let mut anns: Vec<(i32, String)> = front_positions[..front_count.min(3)].to_vec();

        // Insertion annotations on 52
        for i in 0..insertions_r4 {
            anns.push((52, ALPHABET[i].to_string()));
        }

        // Tail positions: (53, 54, 55, 56, 57) skipping abs(min(0, length-5)) items
        let tail_positions = vec![
            (53, " ".to_string()),
            (54, " ".to_string()),
            (55, " ".to_string()),
            (56, " ".to_string()),
            (57, " ".to_string()),
        ];
        let skip = if length_r4 < 5 { 5 - length_r4 } else { 0 };
        anns.extend_from_slice(&tail_positions[skip..]);
        anns
    };
    numbering[4] = (0..length_r4)
        .map(|i| (annotations_r4[i].clone(), regions[4][i].1))
        .collect();

    // Region 6: CDR3 - insertions on 100, max length 36
    let length_r6 = regions[6].len();
    if length_r6 > 36 {
        return (Vec::new(), startindex, endindex);
    }
    let annotations_r6 = get_cdr3_annotations(length_r6, "chothia", "heavy");
    numbering[6] = (0..length_r6)
        .map(|i| (annotations_r6[i].clone(), regions[6][i].1))
        .collect();

    (gap_missing(&numbering), startindex, endindex)
}

/// Apply the Chothia numbering scheme for light chains.
///
/// 7 regions are defined:
///  0 - FW1: simple mapping
///  1 - CDRL1: 24-34, insertions on 30
///  2 - Simple mapping
///  3 - CDRL2: 51-54, insertions on 52
///  4 - FW3: insertions on 68
///  5 - CDRL3: insertions on 95, uses get_cdr3_annotations
///  6 - FW4: simple mapping
pub fn number_chothia_light(
    state_vector: &[StateVectorEntry],
    sequence: &[u8],
) -> (Vec<NumberedResidue>, Option<usize>, Option<usize>) {
    let state_string = b"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXIIIIIIXXXXXXXXXXXXXXXXXXXXXXIIIIIIIXXXXXXXXIXXXXXXXIIXXXXXXXXXXXXXXXXXXXXXXXXXXXIIIIXXXXXXXXXXXXXXX";

    let region_string = b"11111111111111111111111222222222222222223333333333333333444444444445555555555555555555555555555555555555666666666666677777777777";

    let mut region_index_map: HashMap<u8, usize> = HashMap::new();
    region_index_map.insert(b'1', 0);
    region_index_map.insert(b'2', 1);
    region_index_map.insert(b'3', 2);
    region_index_map.insert(b'4', 3);
    region_index_map.insert(b'5', 4);
    region_index_map.insert(b'6', 5);
    region_index_map.insert(b'7', 6);

    let mut rels: Vec<i32> = vec![0, 0, -6, -6, -13, -16, -20];

    let n_regions = 7;
    let exclude_deletions: Vec<usize> = vec![1, 3, 4, 5];

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
        regions[0].clone(), // FW1: pass through
        Vec::new(),         // CDR1 (renumbered)
        regions[2].clone(), // FW2: pass through
        Vec::new(),         // CDR2 (renumbered)
        Vec::new(),         // FW3 (renumbered for insertions on 68)
        Vec::new(),         // CDR3 (renumbered)
        regions[6].clone(), // FW4: pass through
    ];

    // Region 1: CDR1 - range 24-34 (11 base positions), insertions on 30
    let length_r1 = regions[1].len();
    let insertions_r1 = if length_r1 > 11 { length_r1 - 11 } else { 0 };
    let annotations_r1: Vec<(i32, String)> = {
        // Front positions: 24-30 (up to length items)
        let front_count = length_r1.min(7);
        let mut anns: Vec<(i32, String)> = (24..31)
            .take(front_count)
            .map(|p| (p, " ".to_string()))
            .collect();

        // Insertion annotations on 30
        for i in 0..insertions_r1 {
            anns.push((30, ALPHABET[i].to_string()));
        }

        // Tail positions: 31, 32, 33, 34
        let tail = vec![
            (31, " ".to_string()),
            (32, " ".to_string()),
            (33, " ".to_string()),
            (34, " ".to_string()),
        ];
        // Skip abs(min(0, length-11)) items from tail
        let skip = if length_r1 < 11 { 11 - length_r1 } else { 0 };
        // Only take from tail[skip..] to ensure we don't exceed annotations needed
        if skip < tail.len() {
            anns.extend_from_slice(&tail[skip..]);
        }
        anns
    };
    numbering[1] = (0..length_r1)
        .map(|i| (annotations_r1[i].clone(), regions[1][i].1))
        .collect();

    // Region 3: CDR2 - range 51-54 (4 base positions), insertions on 52
    let length_r3 = regions[3].len();
    let insertions_r3 = if length_r3 > 4 { length_r3 - 4 } else { 0 };
    if insertions_r3 > 0 {
        let mut annotations_r3: Vec<(i32, String)> = vec![
            (51, " ".to_string()),
            (52, " ".to_string()),
        ];
        for i in 0..insertions_r3 {
            annotations_r3.push((52, ALPHABET[i].to_string()));
        }
        annotations_r3.push((53, " ".to_string()));
        annotations_r3.push((54, " ".to_string()));
        numbering[3] = (0..length_r3)
            .map(|i| (annotations_r3[i].clone(), regions[3][i].1))
            .collect();
    } else {
        // Let alignment handle gapping when no insertions
        numbering[3] = regions[3].clone();
    }

    // Region 4: FW3 - range 55-88 (34 base positions), insertions on 68, first deletion on 68
    let length_r4 = regions[4].len();
    let insertions_r4 = if length_r4 > 34 { length_r4 - 34 } else { 0 };
    if insertions_r4 > 0 {
        // Insertions on 68
        let mut annotations_r4: Vec<(i32, String)> =
            (55..69).map(|p| (p, " ".to_string())).collect();
        for i in 0..insertions_r4 {
            annotations_r4.push((68, ALPHABET[i].to_string()));
        }
        let tail: Vec<(i32, String)> = (69..89).map(|p| (p, " ".to_string())).collect();
        annotations_r4.extend(tail);
        numbering[4] = (0..length_r4)
            .map(|i| (annotations_r4[i].clone(), regions[4][i].1))
            .collect();
    } else if length_r4 == 33 {
        // First deletion on 68
        let mut annotations_r4: Vec<(i32, String)> =
            (55..68).map(|p| (p, " ".to_string())).collect();
        let tail: Vec<(i32, String)> = (69..89).map(|p| (p, " ".to_string())).collect();
        annotations_r4.extend(tail);
        numbering[4] = (0..length_r4)
            .map(|i| (annotations_r4[i].clone(), regions[4][i].1))
            .collect();
    } else {
        // More deletions - allow alignment to place them
        numbering[4] = regions[4].clone();
    }

    // Region 5: CDR3 - insertions on 95, max length 35
    let length_r5 = regions[5].len();
    if length_r5 > 35 {
        return (Vec::new(), startindex, endindex);
    }
    let annotations_r5 = get_cdr3_annotations(length_r5, "chothia", "light");
    numbering[5] = (0..length_r5)
        .map(|i| (annotations_r5[i].clone(), regions[5][i].1))
        .collect();

    (gap_missing(&numbering), startindex, endindex)
}
