// Validate the HmmAnnotator API against the byte-parity binary.
//   cargo run --release --example annotate -- <db.hmm> <proteins.fa>
use rustyhmmer::api::{Cutoff, HmmAnnotator};
use rustyhmmer::seqio;

fn main() {
    let mut a = std::env::args().skip(1);
    let hmm = a.next().expect("usage: annotate <db.hmm> <proteins.fa>");
    let fa = a.next().expect("need proteins.fa");
    let cutoff = match a.next().as_deref() {
        Some("ga") => Cutoff::GatheringGa,
        Some("tc") => Cutoff::TrustedTc,
        _ => Cutoff::Evalue(10.0),
    };

    let seqs = seqio::read_fasta(&fa).expect("read fasta");
    let ann = HmmAnnotator::from_hmm_file(&hmm).expect("load hmm").with_cutoff(cutoff);
    let hits = ann.search(&seqs);
    // One line per domain: the same coordinates `--domtblout` reports, plus the
    // HMM coverage they make computable.
    for h in &hits {
        for (i, d) in h.domains.iter().enumerate() {
            let coverage = (d.hmm_to - d.hmm_from + 1) as f64 / h.model_len as f64;
            println!(
                "{}\t{}\t{}\t{}/{}\tscore={:.1}\tE={:.2e}\thmm={}..{}\tali={}..{}\tenv={}..{}\tacc={:.2}\tcov={:.3}",
                h.query_name, h.query_acc, h.target_name,
                i + 1, h.n_domains,
                d.score, d.evalue,
                d.hmm_from, d.hmm_to, d.ali_from, d.ali_to, d.env_from, d.env_to,
                d.acc, coverage,
            );
        }
    }
    eprintln!("{} hits", hits.len());
}
