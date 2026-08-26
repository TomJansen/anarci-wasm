use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

use anarci_wasm::{number_sequences_with_hmmer_cli, number_sequences_with_options};
use rustyhmmer::api::HmmAnnotator;
use rustyhmmer::seqio::Seq;

const HEAVY_CHAIN_EXAMPLE: &str = concat!(
    ">test\n",
    "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS\n",
);

const BROAD_GERMLINE_CASES: &[(&str, &str, &str)] = &[
    (
        "human_igh",
        "h",
        "QVQLVQSGA-EVKKPGASVKVSCKASGYTF----TSYGISWVRQAPGQGLEWMGWISAY--NGNTNYAQKLQ-GRVTMTTDTSTSTAYMELRSLRSDDTAVYYCAR----------------------",
    ),
    (
        "human_igk",
        "k",
        "DIQMTQSPSSVSASVGDRVTITCRASQGI------SSWLAWYQQKPGKAPKLLIYAA-------SSLQSGVP-SRFSGSG--SGTDFTLTISSLQPEDFATYYCQQAN--------------------",
    ),
    (
        "human_igl",
        "l",
        "QSVLTQPPS-VSEAPRQRVTISCSGSSSNI----GNNAVNWYQQLPGKAPKLLIYYD-------DLLPSGVS-DRFSGSK--SGTSASLAISGLQSEDEADYYCAAWD--------------------",
    ),
    (
        "human_tcra",
        "a",
        "GQSLEQ-PSEVTAVEGAIVQINCTYQTSG------FYGLSWYQQHDGGAPTFLSYNAL----DGLEET-----GRFSSFLSRSDSYGYLLLQELQMKDSASYFCAVR---------------------",
    ),
    (
        "human_tcrb",
        "b",
        "DAEITQSPRHKITETGRQVTLACHQTWNH-------NNMFWYRQDLGHGLRLIHYSYG----VQDTNKGEVS-DGYSVSRS-NTEDLPLTLESAASSQTSVYFCASSE--------------------",
    ),
    (
        "human_tcrg",
        "g",
        "SSNLEGRTKSVIRQTGSSAEITCDLAEGS------NGYIHWYLHQEGKAPQRLQYYDSY--NSKVVLESGVSPGKYYTYAS-TRNNLRLILRNLIENDSGVYYCATWD--------------------",
    ),
    (
        "human_tcrd",
        "d",
        "AQKITQTQPGMFVQEKEAVTLDCTYDTSDP-----SYGLFWYKQPSSGEMIFLIYQGSY--DQQNATE-----GRYSLNFQKARKSANLVISASQLGDSAMYFCAMRE--------------------",
    ),
];

#[derive(Debug)]
struct HmmerDomain {
    species: String,
    chain_type: String,
    score: f64,
    i_evalue: String,
    ali_from: usize,
    ali_to: usize,
}

#[test]
#[ignore = "requires hmmscan/hmmpress"]
fn default_rust_output_matches_hmmer_for_heavy_chain_example() {
    let hmmer = run_hmmer(HEAVY_CHAIN_EXAMPLE);
    let rust = run_rust(HEAVY_CHAIN_EXAMPLE);

    assert_eq!(rust["species"].as_str().unwrap(), hmmer.species);
    assert_eq!(rust["chain_type"].as_str().unwrap(), hmmer.chain_type);
    assert_eq!(
        rust["seq_start"].as_u64().unwrap() as usize,
        hmmer.ali_from - 1
    );
    assert_eq!(rust["seq_end"].as_u64().unwrap() as usize, hmmer.ali_to);

    let rust_score = rust["bit_score"].as_f64().unwrap();
    assert!(
        (rust_score - hmmer.score).abs() <= 0.1,
        "Rust score {rust_score:.3} differs from HMMER score {:.3}",
        hmmer.score
    );
    assert_eq!(rust["evalue_text"].as_str().unwrap(), hmmer.i_evalue);
}

#[test]
#[ignore = "requires hmmscan/hmmpress"]
fn hmmer_cli_backend_matches_hmmer_for_heavy_chain_example() {
    let hmmer = run_hmmer(HEAVY_CHAIN_EXAMPLE);
    let results =
        number_sequences_with_hmmer_cli(HEAVY_CHAIN_EXAMPLE, "imgt", 80.0, "[]", "heavy").unwrap();
    let domain = &results[0].domains[0];

    assert_eq!(domain.species, hmmer.species);
    assert_eq!(domain.chain_type, hmmer.chain_type);
    assert_eq!(domain.seq_start, hmmer.ali_from - 1);
    assert_eq!(domain.seq_end, hmmer.ali_to);
    assert!((domain.bit_score - hmmer.score).abs() <= 0.1);
    assert_eq!(domain.evalue_text, hmmer.i_evalue);
}

#[test]
#[ignore = "requires rustyhmmer and hmmscan/hmmpress"]
fn rustyhmmer_scores_match_hmmer_for_heavy_chain_example() {
    let hmmer = run_hmmer(HEAVY_CHAIN_EXAMPLE);
    let sequence = HEAVY_CHAIN_EXAMPLE
        .lines()
        .filter(|line| !line.starts_with('>'))
        .collect::<String>();
    let annotator =
        HmmAnnotator::from_hmm_file("data/ALL.hmm").expect("rustyhmmer should load ANARCI HMMs");
    let hits = annotator.search(&[Seq::from_amino("test", "", &sequence)]);
    let best = hits
        .iter()
        .max_by(|left, right| left.best_dom_score.total_cmp(&right.best_dom_score))
        .expect("rustyhmmer should find a heavy-chain hit");

    let (species, chain_type) = best
        .query_name
        .split_once('_')
        .map(|(species, chain)| (species.to_string(), chain.to_string()))
        .unwrap();

    assert_eq!(species, hmmer.species);
    assert_eq!(chain_type, hmmer.chain_type);
    assert_eq!(best.ali_from as usize, hmmer.ali_from);
    assert_eq!(best.ali_to as usize, hmmer.ali_to);
    assert!(
        ((best.best_dom_score as f64) - hmmer.score).abs() <= 0.1,
        "rustyhmmer score {:.3} differs from HMMER score {:.3}",
        best.best_dom_score,
        hmmer.score
    );
    assert_eq!(
        format_evalue(best.best_dom_evalue * annotator.len() as f64),
        hmmer.i_evalue
    );
}

#[test]
#[ignore = "requires rustyhmmer and hmmscan/hmmpress"]
fn default_backend_matches_hmmer_cli_numbering_for_heavy_chain_example() {
    let hmmer =
        number_sequences_with_hmmer_cli(HEAVY_CHAIN_EXAMPLE, "imgt", 80.0, "[]", "heavy").unwrap();
    let rust = run_default_results(HEAVY_CHAIN_EXAMPLE, "imgt", 80.0, "[]", "heavy", "protein");

    let hmmer_domain = &hmmer[0].domains[0];
    let rust_domain = &rust[0].domains[0];

    assert_eq!(rust_domain.species, hmmer_domain.species);
    assert_eq!(rust_domain.chain_type, hmmer_domain.chain_type);
    assert_eq!(rust_domain.seq_start, hmmer_domain.seq_start);
    assert_eq!(rust_domain.seq_end, hmmer_domain.seq_end);
    assert!((rust_domain.bit_score - hmmer_domain.bit_score).abs() <= 0.1);
    assert_eq!(rust_domain.evalue_text, hmmer_domain.evalue_text);
    assert_eq!(
        numbering_tuples(&rust_domain.numbering),
        numbering_tuples(&hmmer_domain.numbering)
    );
}

#[test]
#[ignore = "requires rustyhmmer and hmmscan/hmmpress"]
fn default_backend_matches_hmmer_cli_across_receptor_chain_classes() {
    for (name, restrict, aligned_sequence) in BROAD_GERMLINE_CASES {
        let fasta = format!(">{name}\n{}\n", aligned_sequence.replace('-', ""));
        let hmmer = number_sequences_with_hmmer_cli(&fasta, "imgt", 20.0, "[]", restrict).unwrap();
        let rust = run_default_results(&fasta, "imgt", 20.0, "[]", restrict, "protein");

        assert_numbering_results_match(name, &rust, &hmmer);
    }
}

#[test]
#[ignore = "requires rustyhmmer and hmmscan/hmmpress"]
fn default_backend_matches_hmmer_cli_species_filter_fallback() {
    let hmmer = number_sequences_with_hmmer_cli(
        HEAVY_CHAIN_EXAMPLE,
        "imgt",
        80.0,
        r#"["not_a_species"]"#,
        "heavy",
    )
    .unwrap();
    let rust = run_default_results(
        HEAVY_CHAIN_EXAMPLE,
        "imgt",
        80.0,
        r#"["not_a_species"]"#,
        "heavy",
        "protein",
    );

    assert_numbering_results_match("species_filter_fallback", &rust, &hmmer);
}

fn run_rust(fasta: &str) -> serde_json::Value {
    let json = number_sequences_with_options(fasta, "imgt", 80.0, "[]", "heavy", "protein");
    let results: serde_json::Value = serde_json::from_str(&json).unwrap();
    results[0]["domains"][0].clone()
}

fn run_default_results(
    fasta: &str,
    scheme: &str,
    bit_score_threshold: f64,
    allowed_species_json: &str,
    restrict: &str,
    input_type: &str,
) -> Vec<anarci_wasm::NumberingResult> {
    let json = number_sequences_with_options(
        fasta,
        scheme,
        bit_score_threshold,
        allowed_species_json,
        restrict,
        input_type,
    );
    serde_json::from_str(&json).unwrap()
}

fn assert_numbering_results_match(
    case_name: &str,
    rusty: &[anarci_wasm::NumberingResult],
    hmmer: &[anarci_wasm::NumberingResult],
) {
    assert_eq!(rusty.len(), hmmer.len(), "{case_name}: sequence count");
    assert_eq!(rusty[0].id, hmmer[0].id, "{case_name}: id");
    assert_eq!(rusty[0].sequence, hmmer[0].sequence, "{case_name}: sequence");
    assert!(
        !hmmer[0].domains.is_empty(),
        "{case_name}: HMMER backend found no domains"
    );
    assert_eq!(
        rusty[0].domains.len(),
        hmmer[0].domains.len(),
        "{case_name}: domain count"
    );

    for (domain_index, (rusty_domain, hmmer_domain)) in rusty[0]
        .domains
        .iter()
        .zip(&hmmer[0].domains)
        .enumerate()
    {
        assert_eq!(
            rusty_domain.species, hmmer_domain.species,
            "{case_name} domain {domain_index}: species"
        );
        assert_eq!(
            rusty_domain.chain_type, hmmer_domain.chain_type,
            "{case_name} domain {domain_index}: chain type"
        );
        assert_eq!(
            rusty_domain.seq_start, hmmer_domain.seq_start,
            "{case_name} domain {domain_index}: seq_start"
        );
        assert_eq!(
            rusty_domain.seq_end, hmmer_domain.seq_end,
            "{case_name} domain {domain_index}: seq_end"
        );
        assert!(
            (rusty_domain.bit_score - hmmer_domain.bit_score).abs() <= 0.1,
            "{case_name} domain {domain_index}: rusty score {:.3} differs from HMMER score {:.3}",
            rusty_domain.bit_score,
            hmmer_domain.bit_score
        );
        assert_eq!(
            rusty_domain.evalue_text, hmmer_domain.evalue_text,
            "{case_name} domain {domain_index}: evalue"
        );
        assert_eq!(
            numbering_tuples(&rusty_domain.numbering),
            numbering_tuples(&hmmer_domain.numbering),
            "{case_name} domain {domain_index}: numbering"
        );
    }
}

fn format_evalue(value: f64) -> String {
    if value == 0.0 {
        "0".to_string()
    } else if !(0.001..10000.0).contains(&value) {
        format!("{value:.1e}")
    } else {
        format!("{value:.3}")
    }
}

fn numbering_tuples(numbering: &[anarci_wasm::NumberingEntry]) -> Vec<(i32, String, String)> {
    numbering
        .iter()
        .map(|entry| {
            (
                entry.position,
                entry.insertion.clone(),
                entry.amino_acid.clone(),
            )
        })
        .collect()
}

fn run_hmmer(fasta: &str) -> HmmerDomain {
    require_tool("hmmpress");
    require_tool("hmmscan");

    let temp_dir = unique_temp_dir();
    fs::create_dir_all(&temp_dir).unwrap();

    let hmm = temp_dir.join("ALL.hmm");
    let query = temp_dir.join("query.fa");
    let output = temp_dir.join("hmmscan.txt");

    fs::copy(
        Path::new(env!("CARGO_MANIFEST_DIR")).join("data/ALL.hmm"),
        &hmm,
    )
    .unwrap();
    fs::write(&query, fasta).unwrap();

    let press = Command::new("hmmpress")
        .arg(&hmm)
        .output()
        .expect("failed to run hmmpress");
    assert!(
        press.status.success(),
        "hmmpress failed: {}",
        String::from_utf8_lossy(&press.stderr)
    );

    let scan = Command::new("hmmscan")
        .arg("-o")
        .arg(&output)
        .arg(&hmm)
        .arg(&query)
        .output()
        .expect("failed to run hmmscan");
    assert!(
        scan.status.success(),
        "hmmscan failed: {}",
        String::from_utf8_lossy(&scan.stderr)
    );

    let text = fs::read_to_string(&output).unwrap();
    let _ = fs::remove_dir_all(&temp_dir);
    parse_best_domain(&text)
}

fn parse_best_domain(text: &str) -> HmmerDomain {
    let mut current_model: Option<String> = None;

    for line in text.lines() {
        if let Some(model) = line.trim_start().strip_prefix(">> ") {
            current_model = Some(model.split_whitespace().next().unwrap().to_string());
            continue;
        }

        let Some(model) = current_model.as_ref() else {
            continue;
        };
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 12 || fields[0] != "1" || fields[1] != "!" {
            continue;
        }

        let (species, chain_type) = model
            .split_once('_')
            .map(|(species, chain)| (species.to_string(), chain.to_string()))
            .unwrap();

        return HmmerDomain {
            species,
            chain_type,
            score: fields[2].parse().unwrap(),
            i_evalue: fields[5].to_string(),
            ali_from: fields[9].parse().unwrap(),
            ali_to: fields[10].parse().unwrap(),
        };
    }

    panic!("no HMMER domain row found");
}

fn require_tool(name: &str) {
    let status = Command::new("which")
        .arg(name)
        .status()
        .unwrap_or_else(|_| panic!("failed to check for {name}"));
    assert!(status.success(), "{name} is required for this ignored test");
}

fn unique_temp_dir() -> PathBuf {
    let suffix = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    std::env::temp_dir().join(format!("anarci-hmmer-parity-{suffix}"))
}
