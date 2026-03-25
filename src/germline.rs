use crate::alignment::StateVectorEntry;
use serde::Deserialize;
use std::collections::{BTreeMap, BTreeSet};
use std::sync::OnceLock;

type GeneMap = BTreeMap<String, String>;
type SpeciesMap = BTreeMap<String, GeneMap>;
type ChainMap = BTreeMap<String, SpeciesMap>;

const GERMLINE_DATA: &str = include_str!("../data/germlines.json");

static GERMLINES: OnceLock<GermlineDb> = OnceLock::new();

#[derive(Deserialize)]
struct GermlineDb {
    #[serde(rename = "V")]
    v: ChainMap,
    #[serde(rename = "J")]
    j: ChainMap,
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

pub fn assign_germline(
    state_vector: &[StateVectorEntry],
    sequence: &[u8],
    chain_type: &str,
    allowed_species: Option<&[String]>,
) -> Option<GermlineAssignment> {
    let db = germlines();
    let state_sequence = build_state_sequence(state_vector, sequence);

    let v_hit = best_gene(&state_sequence, db.v.get(chain_type)?, allowed_species)?;
    let j_hit = db
        .j
        .get(chain_type)
        .and_then(|species_map| species_map.get(&v_hit.species))
        .and_then(|genes| best_gene_within_species(&state_sequence, &v_hit.species, genes));

    Some(GermlineAssignment {
        identity_species: v_hit.species.clone(),
        v_gene: v_hit.gene,
        v_identity: v_hit.identity,
        j_gene: j_hit.as_ref().map(|hit| hit.gene.clone()).unwrap_or_default(),
        j_identity: j_hit.map(|hit| hit.identity).unwrap_or(0.0),
    })
}

pub fn available_species() -> Vec<String> {
    let db = germlines();
    let mut species = BTreeSet::new();
    for chain_map in db.v.values() {
        for species_name in chain_map.keys() {
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

fn best_gene(
    state_sequence: &[u8; 128],
    species_map: &SpeciesMap,
    allowed_species: Option<&[String]>,
) -> Option<BestGene> {
    let mut best: Option<BestGene> = None;
    let allowed_species: Option<BTreeSet<&str>> = allowed_species
        .map(|species| species.iter().map(|s| s.as_str()).collect());

    for (species, genes) in species_map {
        if let Some(allowed_species) = &allowed_species {
            if !allowed_species.contains(species.as_str()) {
                continue;
            }
        }
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
    genes: &GeneMap,
) -> Option<BestGene> {
    let mut best: Option<BestGene> = None;
    for (gene, germline_sequence) in genes {
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
