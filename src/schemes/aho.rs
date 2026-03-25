/// Aho numbering scheme implementation.
/// Port of number_aho from Python ANARCI schemes.py.

use crate::alignment::*;
use std::collections::HashMap;

/// Apply the Aho numbering scheme for any chain type.
///
/// The Aho scheme uses 128 HMM states (all 'X'), 11 regions (A-K),
/// and renumbers CDR1, CDR2, FW3, and CDR3 with scheme-specific
/// gap placement rules that depend on chain_type.
pub fn number_aho(
    state_vector: &[StateVectorEntry],
    sequence: &[u8],
    chain_type: &str,
) -> (Vec<NumberedResidue>, Option<usize>, Option<usize>) {
    // State string - all X's, every position maps directly
    let state_string = b"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX";

    // Region string - 11 regions (A=empty now in B, B=FW1, C, D=CDR1, E, F=CDR2, G=empty now in H, H=FW3, I, J=CDR3, K=FW4)
    let region_string = b"BBBBBBBBBBCCCCCCCCCCCCCCDDDDDDDDDDDDDDDDEEEEEEEEEEEEEEEFFFFFFFFFFFFFFFFFFFFHHHHHHHHHHHHHHHHIIIIIIIIIIIIIJJJJJJJJJJJJJKKKKKKKKKKK";

    let mut region_index_map: HashMap<u8, usize> = HashMap::new();
    region_index_map.insert(b'A', 0);
    region_index_map.insert(b'B', 1);
    region_index_map.insert(b'C', 2);
    region_index_map.insert(b'D', 3);
    region_index_map.insert(b'E', 4);
    region_index_map.insert(b'F', 5);
    region_index_map.insert(b'G', 6);
    region_index_map.insert(b'H', 7);
    region_index_map.insert(b'I', 8);
    region_index_map.insert(b'J', 9);
    region_index_map.insert(b'K', 10);

    let mut rels: Vec<i32> = vec![0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 21];

    let n_regions = 11;
    let exclude_deletions: Vec<usize> = vec![1, 3, 4, 5, 7, 9];

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

    // Build the numbering assembly: regions 3,5,7,9 are renumbered (replaced with [])
    let mut numbering: Vec<Vec<NumberedResidue>> = vec![
        regions[0].clone(),  // region A (empty)
        regions[1].clone(),  // region B (FW1) - will be renumbered below
        regions[2].clone(),  // region C
        Vec::new(),          // region D (CDR1) - renumbered
        regions[4].clone(),  // region E
        Vec::new(),          // region F (CDR2) - renumbered
        regions[6].clone(),  // region G (empty)
        Vec::new(),          // region H (FW3) - renumbered
        regions[8].clone(),  // region I
        regions[9].clone(),  // region J (CDR3) - will be replaced below
        regions[10].clone(), // region K
    ];

    // ========================
    // FW1 - Region B (index 1)
    // ========================
    // Place indels on position 8
    let length = regions[1].len();
    if length > 0 {
        let start = regions[1][0].0 .0;
        let stretch_len = 10 - (start - 1);
        if length as i32 > stretch_len {
            // Insertions are present. Place on 8
            let num_insertions = length as i32 - stretch_len;
            let mut annotations: Vec<(i32, String)> = Vec::new();
            for p in start..9 {
                annotations.push((p, " ".to_string()));
            }
            for a in 0..num_insertions as usize {
                annotations.push((8, ALPHABET[a].to_string()));
            }
            annotations.push((9, " ".to_string()));
            annotations.push((10, " ".to_string()));
            numbering[1] = annotations
                .iter()
                .enumerate()
                .take(length)
                .map(|(i, ann)| (ann.clone(), regions[1][i].1))
                .collect();
        } else {
            // Deletions: ordered_deletions starting with 8, then rest
            let mut ordered_deletions: Vec<(i32, String)> = vec![(8, " ".to_string())];
            for p in start..11 {
                if p != 8 {
                    ordered_deletions.push((p, " ".to_string()));
                }
            }
            let skip = (stretch_len - length as i32).max(0) as usize;
            let mut annotations: Vec<(i32, String)> = ordered_deletions[skip..].to_vec();
            annotations.sort_by(|a, b| a.0.cmp(&b.0).then_with(|| a.1.cmp(&b.1)));
            numbering[1] = annotations
                .iter()
                .enumerate()
                .take(length)
                .map(|(i, ann)| (ann.clone(), regions[1][i].1))
                .collect();
        }
    }

    // ========================
    // CDR1 - Region D (index 3)
    // ========================
    // Gap ordering depends on chain type
    let cdr1_deletions: Vec<i32> = match chain_type {
        "L" => vec![28, 36, 35, 37, 34, 38, 27, 29, 33, 39, 32, 40, 26, 30, 25, 31, 41, 42],
        "K" => vec![28, 27, 36, 35, 37, 34, 38, 33, 39, 32, 40, 29, 26, 30, 25, 31, 41, 42],
        "H" => vec![28, 36, 35, 37, 34, 38, 27, 33, 39, 32, 40, 29, 26, 30, 25, 31, 41, 42],
        "A" => vec![28, 36, 35, 37, 34, 38, 33, 39, 27, 32, 40, 29, 26, 30, 25, 31, 41, 42],
        "B" => vec![28, 36, 35, 37, 34, 38, 33, 39, 27, 32, 40, 29, 26, 30, 25, 31, 41, 42],
        "D" => vec![28, 36, 35, 37, 34, 38, 27, 33, 39, 32, 40, 29, 26, 30, 25, 31, 41, 42],
        "G" => vec![28, 36, 35, 37, 34, 38, 27, 33, 39, 32, 40, 29, 26, 30, 25, 31, 41, 42],
        _ => vec![28, 36, 35, 37, 34, 38, 27, 33, 39, 32, 40, 29, 26, 30, 25, 31, 41, 42],
    };

    let length = regions[3].len();
    let skip = if 18 > length { 18 - length } else { 0 };
    let mut annotations: Vec<(i32, String)> = {
        let mut selected: Vec<i32> = cdr1_deletions[skip..].to_vec();
        selected.sort();
        selected.iter().map(|&i| (i, " ".to_string())).collect()
    };

    let insertions = if length > 18 { length - 18 } else { 0 };
    if insertions > 26 {
        return (Vec::new(), startindex, endindex);
    } else if insertions > 0 {
        // Insert on position 36
        let insertat = annotations.iter().position(|a| a.0 == 36 && a.1 == " ").unwrap() + 1;
        let ins_annotations: Vec<(i32, String)> = (0..insertions)
            .map(|a| (36, ALPHABET[a].to_string()))
            .collect();
        annotations.splice(insertat..insertat, ins_annotations);
    }

    numbering[3] = (0..length)
        .map(|i| (annotations[i].clone(), regions[3][i].1))
        .collect();

    // ========================
    // CDR2 - Region F (index 5)
    // ========================
    let cdr2_deletions: Vec<i32> = if chain_type == "A" {
        vec![74, 73, 63, 62, 64, 61, 65, 60, 66, 59, 67, 58, 68, 69, 70, 71, 72, 75, 76, 77]
    } else {
        vec![63, 62, 64, 61, 65, 60, 66, 59, 67, 58, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77]
    };

    let length = regions[5].len();
    let skip = if 20 > length { 20 - length } else { 0 };
    let mut annotations: Vec<(i32, String)> = {
        let mut selected: Vec<i32> = cdr2_deletions[skip..].to_vec();
        selected.sort();
        selected.iter().map(|&i| (i, " ".to_string())).collect()
    };

    let insertions = if length > 20 { length - 20 } else { 0 };
    if insertions > 26 {
        return (Vec::new(), startindex, endindex);
    } else if insertions > 0 {
        let insertat = annotations.iter().position(|a| a.0 == 63 && a.1 == " ").unwrap() + 1;
        let ins_annotations: Vec<(i32, String)> = (0..insertions)
            .map(|a| (63, ALPHABET[a].to_string()))
            .collect();
        annotations.splice(insertat..insertat, ins_annotations);
    }

    numbering[5] = (0..length)
        .map(|i| (annotations[i].clone(), regions[5][i].1))
        .collect();

    // ========================
    // FW3 - Region H (index 7)
    // ========================
    let fw3_deletions: Vec<i32> = vec![86, 85, 87, 84, 88, 83, 89, 82, 90, 81, 91, 80, 92, 79, 93, 78];
    let length = regions[7].len();
    let skip = if 16 > length { 16 - length } else { 0 };
    let mut annotations: Vec<(i32, String)> = {
        let mut selected: Vec<i32> = fw3_deletions[skip..].to_vec();
        selected.sort();
        selected.iter().map(|&i| (i, " ".to_string())).collect()
    };

    let insertions = if length > 16 { length - 16 } else { 0 };
    if insertions > 26 {
        return (Vec::new(), startindex, endindex);
    } else if insertions > 0 {
        let insertat = annotations.iter().position(|a| a.0 == 85 && a.1 == " ").unwrap() + 1;
        let ins_annotations: Vec<(i32, String)> = (0..insertions)
            .map(|a| (85, ALPHABET[a].to_string()))
            .collect();
        annotations.splice(insertat..insertat, ins_annotations);
    }

    numbering[7] = (0..length)
        .map(|i| (annotations[i].clone(), regions[7][i].1))
        .collect();

    // ========================
    // CDR3 - Region J (index 9)
    // ========================
    let cdr3_deletions: Vec<i32> = vec![
        123, 124, 122, 125, 121, 126, 120, 127, 119, 128, 118, 129, 117, 130,
        116, 131, 115, 132, 114, 133, 113, 134, 112, 135, 111, 136, 110, 137,
        109, 138, 108, 107,
    ];
    let length = regions[9].len();
    let skip = if 32 > length { 32 - length } else { 0 };
    let mut annotations: Vec<(i32, String)> = {
        let mut selected: Vec<i32> = cdr3_deletions[skip..].to_vec();
        selected.sort();
        selected.iter().map(|&i| (i, " ".to_string())).collect()
    };

    let insertions = if length > 32 { length - 32 } else { 0 };
    if insertions > 26 {
        return (Vec::new(), startindex, endindex);
    } else if insertions > 0 {
        let insertat = annotations.iter().position(|a| a.0 == 123 && a.1 == " ").unwrap() + 1;
        let ins_annotations: Vec<(i32, String)> = (0..insertions)
            .map(|a| (123, ALPHABET[a].to_string()))
            .collect();
        annotations.splice(insertat..insertat, ins_annotations);
    }

    numbering[9] = (0..length)
        .map(|i| (annotations[i].clone(), regions[9][i].1))
        .collect();

    // ========================
    // Final assembly with gap_missing
    // ========================
    let mut result = gap_missing(&numbering);

    // AHo includes one extra position (149) for light chains if position 148 is occupied
    if !result.is_empty() {
        let last = &result[result.len() - 1];
        if last.0 .0 == 148
            && last.0 .1 == " "
            && last.1 != '-'
            && endindex.is_some()
            && endindex.unwrap() + 1 < sequence.len()
        {
            let next_aa = sequence[endindex.unwrap() + 1] as char;
            result.push(((149, " ".to_string()), next_aa));
            return (result, startindex, Some(endindex.unwrap() + 1));
        }
    }

    (result, startindex, endindex)
}
