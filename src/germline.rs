use crate::alignment::StateVectorEntry;
use serde::de::{Deserializer, MapAccess, Visitor};
use serde::Deserialize;
use std::collections::BTreeSet;
use std::fmt;
use std::marker::PhantomData;
use std::sync::OnceLock;

/// A JSON object deserialized into a Vec of (key, value) pairs, preserving the
/// original insertion order. Used so that `get_hmm_length` can pick the *first*
/// J gene exactly as Python ANARCI does (`serde_json`'s default map types sort
/// or are unordered without the `preserve_order` feature).
struct Ordered<V>(Vec<(String, V)>);

impl<'de, V: Deserialize<'de>> Deserialize<'de> for Ordered<V> {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        struct OrderedVisitor<V>(PhantomData<V>);

        impl<'de, V: Deserialize<'de>> Visitor<'de> for OrderedVisitor<V> {
            type Value = Ordered<V>;

            fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
                f.write_str("a JSON object")
            }

            fn visit_map<M: MapAccess<'de>>(self, mut map: M) -> Result<Self::Value, M::Error> {
                let mut entries = Vec::with_capacity(map.size_hint().unwrap_or(0));
                while let Some((key, value)) = map.next_entry::<String, V>()? {
                    entries.push((key, value));
                }
                Ok(Ordered(entries))
            }
        }

        deserializer.deserialize_map(OrderedVisitor(PhantomData))
    }
}

// Order-preserving view of the germline JSON, used only for `get_hmm_length`.
type OrderedGeneMap = Ordered<String>;
type OrderedSpeciesMap = Ordered<OrderedGeneMap>;
type OrderedChainMap = Ordered<OrderedSpeciesMap>;

const GERMLINE_DATA: &str = include_str!("../data/germlines.json");

static GERMLINES: OnceLock<GermlineDb> = OnceLock::new();

#[derive(Deserialize)]
struct GermlineDb {
    /// V germlines preserving JSON order. Germline assignment iterates species
    /// and genes in this order so that ties resolve exactly as Python ANARCI's
    /// `max(seq_ids, key=...)` does (first maximum in insertion order).
    #[serde(rename = "V")]
    v: OrderedChainMap,
    /// J germlines preserving JSON order, for assignment and `get_hmm_length`.
    #[serde(rename = "J")]
    j: OrderedChainMap,
}

impl Ordered<OrderedSpeciesMap> {
    fn chain(&self, chain_type: &str) -> Option<&OrderedSpeciesMap> {
        self.0.iter().find(|(c, _)| c == chain_type).map(|(_, v)| v)
    }
}

impl Ordered<OrderedGeneMap> {
    fn species(&self, species: &str) -> Option<&OrderedGeneMap> {
        self.0.iter().find(|(s, _)| s == species).map(|(_, v)| v)
    }
}

#[derive(Clone)]
pub struct GermlineAssignment {
    pub identity_species: String,
    pub v_gene: String,
    pub v_identity: f64,
    pub j_gene: String,
    pub j_identity: f64,
}

#[derive(Clone)]
struct BestGene {
    species: String,
    gene: String,
    identity: f64,
}

fn germlines() -> &'static GermlineDb {
    GERMLINES.get_or_init(|| {
        serde_json::from_str(GERMLINE_DATA)
            .expect("embedded germline database should parse successfully")
    })
}

/// The canonical species iteration order used when no `allowed_species` filter
/// is given. Python ANARCI derives this from the JSON key order of `V.H`
/// (`all_species`) and applies it to every chain type.
fn canonical_species_order() -> Vec<&'static str> {
    germlines()
        .v
        .chain("H")
        .map(|species_map| species_map.0.iter().map(|(s, _)| s.as_str()).collect())
        .unwrap_or_default()
}

pub fn assign_germline(
    state_vector: &[StateVectorEntry],
    sequence: &[u8],
    chain_type: &str,
    allowed_species: Option<&[String]>,
) -> Option<GermlineAssignment> {
    let db = germlines();
    let state_sequence = build_state_sequence(state_vector, sequence);

    let v_chain = db.v.chain(chain_type)?;

    // Determine the species iteration order, matching Python ANARCI:
    // - With an explicit filter, iterate the caller's list in order. If any
    //   requested species is absent for this chain, Python returns no assignment.
    // - Otherwise iterate the canonical `all_species` order (JSON order of V.H).
    let species_order: Vec<&str> = match allowed_species {
        Some(allowed) => {
            if allowed.iter().any(|sp| v_chain.species(sp).is_none()) {
                return None;
            }
            allowed.iter().map(|s| s.as_str()).collect()
        }
        None => canonical_species_order(),
    };

    // V gene: highest identity wins; ties resolve to the first encountered in
    // (species_order, then JSON gene order) — exactly Python's first `max`.
    let v_hit = best_gene(&state_sequence, v_chain, &species_order)?;

    // J gene uses the species assigned to the V gene.
    let j_hit = db
        .j
        .chain(chain_type)
        .and_then(|species_map| species_map.species(&v_hit.species))
        .and_then(|genes| best_gene_within_species(&state_sequence, &v_hit.species, genes));

    Some(GermlineAssignment {
        identity_species: v_hit.species.clone(),
        v_gene: v_hit.gene,
        v_identity: v_hit.identity,
        j_gene: j_hit.as_ref().map(|hit| hit.gene.clone()).unwrap_or_default(),
        j_identity: j_hit.map(|hit| hit.identity).unwrap_or(0.0),
    })
}

/// Length of the HMM for a given species/chain type.
/// Mirrors Python ANARCI's `get_hmm_length`: the number of non-trailing-gap
/// positions in the first J germline for the chain/species, defaulting to 128.
pub fn get_hmm_length(species: &str, chain_type: &str) -> usize {
    let db = germlines();
    db.j
        .chain(chain_type)
        .and_then(|species_map| species_map.species(species))
        .and_then(|genes| genes.0.first())
        .map(|(_, seq)| seq.trim_end_matches('-').chars().count())
        .unwrap_or(128)
}

pub fn available_species() -> Vec<String> {
    let db = germlines();
    let mut species = BTreeSet::new();
    for (_, species_map) in &db.v.0 {
        for (species_name, _) in &species_map.0 {
            species.insert(species_name.clone());
        }
    }
    species.into_iter().collect()
}

fn build_state_sequence(state_vector: &[StateVectorEntry], sequence: &[u8]) -> [u8; 128] {
    let mut state_sequence = [b'-'; 128];
    for &((state_id, state_type), seq_index) in state_vector {
        if state_type != 'm' {
            continue;
        }
        if !(1..=128).contains(&state_id) {
            continue;
        }
        if let Some(seq_index) = seq_index {
            if let Some(&aa) = sequence.get(seq_index) {
                state_sequence[state_id - 1] = aa.to_ascii_uppercase();
            }
        }
    }
    state_sequence
}

/// Best V gene across the given species, visited in `species_order`. Within each
/// species, genes are visited in JSON order. The first maximum-identity hit wins,
/// matching Python's `max(seq_ids, key=...)` over an insertion-ordered dict.
fn best_gene(
    state_sequence: &[u8; 128],
    species_map: &OrderedSpeciesMap,
    species_order: &[&str],
) -> Option<BestGene> {
    let mut best: Option<BestGene> = None;
    for &species in species_order {
        let Some(genes) = species_map.species(species) else {
            continue;
        };
        if let Some(hit) = best_gene_within_species(state_sequence, species, genes) {
            match &best {
                Some(current) if current.identity >= hit.identity => {}
                _ => best = Some(hit),
            }
        }
    }
    best
}

fn best_gene_within_species(
    state_sequence: &[u8; 128],
    species: &str,
    genes: &OrderedGeneMap,
) -> Option<BestGene> {
    let mut best: Option<BestGene> = None;
    for (gene, germline_sequence) in &genes.0 {
        let identity = get_identity(state_sequence, germline_sequence.as_bytes());
        match &best {
            Some(current) if current.identity >= identity => {}
            _ => {
                best = Some(BestGene {
                    species: species.to_string(),
                    gene: gene.clone(),
                    identity,
                })
            }
        }
    }
    best
}

fn get_identity(state_sequence: &[u8; 128], germline_sequence: &[u8]) -> f64 {
    let mut n = 0usize;
    let mut m = 0usize;

    for i in 0..128 {
        let Some(&germline_aa) = germline_sequence.get(i) else {
            break;
        };
        if germline_aa == b'-' {
            continue;
        }
        if state_sequence[i] == germline_aa {
            m += 1;
        }
        n += 1;
    }

    if n == 0 {
        0.0
    } else {
        m as f64 / n as f64
    }
}

#[cfg(test)]
mod tests {
    use super::{assign_germline, canonical_species_order, get_hmm_length};

    #[test]
    fn canonical_species_order_follows_v_heavy_json_order() {
        // Must match Python's `all_species = list(all_germlines['V']['H'].keys())`,
        // which is JSON insertion order, not alphabetical.
        assert_eq!(
            canonical_species_order(),
            vec!["human", "mouse", "rat", "rabbit", "rhesus", "pig", "alpaca", "cow"]
        );
    }

    #[test]
    fn germline_ties_resolve_to_first_in_canonical_then_json_order() {
        // An empty state vector yields an all-gap state sequence, so every germline
        // scores identity 0.0 — a total tie. Python's first-`max` then picks the
        // first species in canonical order (human) and its first V gene in JSON
        // order (IGHV1-18*01).
        let assignment = assign_germline(&[], &[], "H", None).unwrap();
        assert_eq!(assignment.identity_species, "human");
        assert_eq!(assignment.v_gene, "IGHV1-18*01");
        assert_eq!(assignment.v_identity, 0.0);
    }

    #[test]
    fn germline_ties_respect_allowed_species_list_order() {
        // With an explicit filter the species are visited in the caller's order,
        // so on a total tie the first listed species wins regardless of canonical
        // order. mouse before human => mouse is assigned.
        let mouse_first =
            assign_germline(&[], &[], "H", Some(&["mouse".into(), "human".into()])).unwrap();
        assert_eq!(mouse_first.identity_species, "mouse");

        let human_first =
            assign_germline(&[], &[], "H", Some(&["human".into(), "mouse".into()])).unwrap();
        assert_eq!(human_first.identity_species, "human");
    }

    #[test]
    fn germline_assignment_aborts_if_requested_species_missing_for_chain() {
        // Python returns no germline assignment when any requested species is not
        // present for the chain type. Alpha chains only have human and mouse.
        assert!(assign_germline(&[], &[], "A", Some(&["alpaca".into()])).is_none());
        assert!(assign_germline(&[], &[], "A", Some(&["human".into()])).is_some());
    }

    #[test]
    fn hmm_length_uses_first_j_gene_in_json_order() {
        // For these species/chains the J germlines disagree on length (127 vs 128).
        // Python ANARCI takes the *first* gene in JSON order; matching that here
        // requires preserving insertion order rather than sorting the gene keys.
        assert_eq!(get_hmm_length("human", "L"), 128);
        assert_eq!(get_hmm_length("human", "G"), 128);
        assert_eq!(get_hmm_length("alpaca", "H"), 128);
    }

    #[test]
    fn hmm_length_defaults_to_128_for_unknown_chain_or_species() {
        assert_eq!(get_hmm_length("nonexistent", "H"), 128);
        assert_eq!(get_hmm_length("human", "Z"), 128);
    }
}
