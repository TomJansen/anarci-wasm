/// Numbering scheme implementations.
/// Each scheme converts a state vector (IMGT HMM alignment) to a numbered sequence.

pub mod imgt;
pub mod chothia;
pub mod kabat;
pub mod martin;
pub mod aho;
pub mod wolfguy;

use crate::alignment::{NumberedResidue, StateVectorEntry};

/// Apply a numbering scheme to a state vector and sequence.
pub fn number_sequence(
    state_vector: &[StateVectorEntry],
    sequence: &[u8],
    scheme: &str,
    chain_type: &str,
) -> Option<(Vec<NumberedResidue>, Option<usize>, Option<usize>)> {
    match scheme {
        "imgt" => Some(imgt::number_imgt(state_vector, sequence)),
        "chothia" => match chain_type {
            "H" => Some(chothia::number_chothia_heavy(state_vector, sequence)),
            "K" | "L" => Some(chothia::number_chothia_light(state_vector, sequence)),
            _ => None,
        },
        "kabat" => match chain_type {
            "H" => Some(kabat::number_kabat_heavy(state_vector, sequence)),
            "K" | "L" => Some(kabat::number_kabat_light(state_vector, sequence)),
            _ => None,
        },
        "martin" => match chain_type {
            "H" => Some(martin::number_martin_heavy(state_vector, sequence)),
            "K" | "L" => Some(martin::number_martin_light(state_vector, sequence)),
            _ => None,
        },
        "aho" => Some(aho::number_aho(state_vector, sequence, chain_type)),
        "wolfguy" => match chain_type {
            "H" => Some(wolfguy::number_wolfguy_heavy(state_vector, sequence)),
            "K" | "L" => Some(wolfguy::number_wolfguy_light(state_vector, sequence)),
            _ => None,
        },
        _ => None,
    }
}

/// Resolve scheme name aliases
pub fn resolve_scheme_name(scheme: &str) -> Option<&'static str> {
    match scheme.to_lowercase().as_str() {
        "imgt" | "i" => Some("imgt"),
        "chothia" | "c" => Some("chothia"),
        "kabat" | "k" => Some("kabat"),
        "martin" | "m" => Some("martin"),
        "aho" | "a" => Some("aho"),
        "wolfguy" | "w" => Some("wolfguy"),
        _ => None,
    }
}

/// Map chain type character to chain class
pub fn chain_type_to_class(chain_type: &str) -> &'static str {
    match chain_type {
        "H" => "H",
        "K" | "L" => "L",
        "A" => "A",
        "B" => "B",
        "G" => "G",
        "D" => "D",
        _ => "H",
    }
}
