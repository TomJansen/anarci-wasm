/// IMGT numbering scheme implementation.
/// Port of number_imgt from Python ANARCI schemes.py.

use crate::alignment::*;
use std::collections::HashMap;

/// Apply the IMGT numbering scheme for heavy or light chains.
///
/// All 128 HMM states map directly to IMGT positions (state_string is all 'X').
/// 7 regions are defined, with CDR1 (region 1), CDR2 (region 3), and CDR3 (region 5)
/// renumbered using symmetric IMGT gapping rules.
pub fn number_imgt(
    state_vector: &[StateVectorEntry],
    sequence: &[u8],
) -> (Vec<NumberedResidue>, Option<usize>, Option<usize>) {
    // State string - all X's, every IMGT position maps directly
    let state_string = b"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX";

    // Region string - 7 regions
    let region_string = b"11111111111111111111111111222222222222333333333333333334444444444555555555555555555555555555555555555555666666666666677777777777";

    let mut region_index_map: HashMap<u8, usize> = HashMap::new();
    region_index_map.insert(b'1', 0);
    region_index_map.insert(b'2', 1);
    region_index_map.insert(b'3', 2);
    region_index_map.insert(b'4', 3);
    region_index_map.insert(b'5', 4);
    region_index_map.insert(b'6', 5);
    region_index_map.insert(b'7', 6);

    let mut rels: Vec<i32> = vec![0, 0, 0, 0, 0, 0, 0];

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

    // Build the renumbered output
    let mut numbering: Vec<Vec<NumberedResidue>> = vec![
        regions[0].clone(), // FW1
        Vec::new(),         // CDR1 (renumbered below)
        regions[2].clone(), // FW2
        Vec::new(),         // CDR2 (renumbered below)
        regions[4].clone(), // FW3
        Vec::new(),         // CDR3 (renumbered below)
        regions[6].clone(), // FW4
    ];

    // CDR1: range 27 (inc) to 39 (exc), max length 12
    let cdr1seq: Vec<u8> = regions[1]
        .iter()
        .filter(|x| x.1 != '-')
        .map(|x| x.1 as u8)
        .collect();
    let cdr1length = cdr1seq.len();
    let mut si: usize = 0;
    let mut prev_state: i32 = 26;
    for ann in get_imgt_cdr(cdr1length, 12, 27, 39) {
        match ann {
            None => {
                numbering[1].push(((prev_state + 1, " ".to_string()), '-'));
                prev_state += 1;
            }
            Some((pos, ins)) => {
                numbering[1].push(((pos, ins), cdr1seq[si] as char));
                prev_state = pos;
                si += 1;
            }
        }
    }

    // CDR2: range 56 (inc) to 66 (exc), max length 10
    let cdr2seq: Vec<u8> = regions[3]
        .iter()
        .filter(|x| x.1 != '-')
        .map(|x| x.1 as u8)
        .collect();
    let cdr2length = cdr2seq.len();
    si = 0;
    prev_state = 55;
    for ann in get_imgt_cdr(cdr2length, 10, 56, 66) {
        match ann {
            None => {
                numbering[3].push(((prev_state + 1, " ".to_string()), '-'));
                prev_state += 1;
            }
            Some((pos, ins)) => {
                numbering[3].push(((pos, ins), cdr2seq[si] as char));
                prev_state = pos;
                si += 1;
            }
        }
    }

    // CDR3: range 105 (inc) to 118 (exc), max length 13
    // Maximum practical length is 117 (13 positions + 52 insertions)
    let cdr3seq: Vec<u8> = regions[5]
        .iter()
        .filter(|x| x.1 != '-')
        .map(|x| x.1 as u8)
        .collect();
    let cdr3length = cdr3seq.len();
    if cdr3length > 117 {
        // Too many insertions; do not apply numbering
        return (Vec::new(), startindex, endindex);
    }
    si = 0;
    let mut previous_state_id: i32 = 104;
    for ann in get_imgt_cdr(cdr3length, 13, 105, 118) {
        match ann {
            None => {
                numbering[5].push(((previous_state_id + 1, " ".to_string()), '-'));
                previous_state_id += 1;
            }
            Some((pos, ins)) => {
                numbering[5].push(((pos, ins), cdr3seq[si] as char));
                previous_state_id = pos;
                si += 1;
            }
        }
    }

    // Fill in gaps for missing positions and return
    (gap_missing(&numbering), startindex, endindex)
}
