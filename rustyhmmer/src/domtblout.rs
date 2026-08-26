//! `--domtblout`: parseable table of per-domain hits.
//!
//! Faithful transcription of `p7_tophits_TabularDomains()`
//! (hmmer-3.4 src/p7_tophits.c:1714-1789), including its three header lines and the
//! exact `printf` field widths, so the output is byte-identical to C's.
//!
//! Domain reporting/inclusion flags come from `p7_tophits_Threshold()`
//! (p7_tophits.c:928-1000) for E-value thresholds, or from the pipeline's
//! `use_bit_cutoffs` block (p7_pipeline.c:913-931) for `--cut_ga/tc/nc`.

use crate::pipeline::{dom_acc, dombias_bits, Hit};
use crate::pli::Pipeline;
use crate::tblout::fmt_g2;

/// Column widths, as computed at the top of `p7_tophits_TabularDomains()`
/// (p7_tophits.c:1717-1720).
#[derive(Clone, Copy)]
pub struct DomWidths {
    pub tnamew: usize,
    pub taccw: usize,
    pub qnamew: usize,
    pub qaccw: usize,
}

impl DomWidths {
    /// `th` holds exactly the reportable targets — the pipeline only creates a hit
    /// for a target that passed `p7_pli_TargetReportable()` (p7_pipeline.c:838) — so
    /// the max-name/accession scans run over the reported set.
    pub fn compute(qname: &str, qacc: &str, reported: &[&Hit]) -> Self {
        DomWidths {
            tnamew: 20.max(reported.iter().map(|h| h.name.len()).max().unwrap_or(0)),
            // C, p7_tophits.c:1720: `ESL_MAX(10, p7_tophits_GetMaxAccessionLength(th))`.
            // FASTA carries no accession, but the EMBL/UniProt `AC` and GenBank/DDBJ
            // `VERSION` lines do.
            taccw: 10.max(reported.iter().map(|h| h.acc.len()).max().unwrap_or(0)),
            qnamew: 20.max(qname.len()),
            qaccw: if qacc.is_empty() { 10 } else { 10.max(qacc.len()) },
        }
    }
}

/// The three `#` header lines.
///
/// C, p7_tophits.c:1726-1735:
///   if (fprintf(ofp, "#%*s %22s %40s %11s %11s %11s\n", tnamew+qnamew-1+15+taccw+qaccw, "",
///        "--- full sequence ---", "-------------- this domain -------------",
///        "hmm coord", "ali coord", "env coord") < 0) ...
///   if (fprintf(ofp, "#%-*s %-*s %5s %-*s %-*s %5s %9s %6s %5s %3s %3s %9s %9s %6s %5s %5s %5s %5s %5s %5s %5s %4s %s\n",
///        tnamew-1, " target name", taccw, "accession", "tlen", qnamew, "query name",
///        qaccw, "accession", "qlen", "E-value", "score", "bias", "#", "of",
///        "c-Evalue", "i-Evalue", "score", "bias", "from", "to", "from", "to",
///        "from", "to", "acc", "description of target") < 0) ...
///   if (fprintf(ofp, "#%*s %*s %5s %*s %*s %5s %9s %6s %5s %3s %3s %9s %9s %6s %5s %5s %5s %5s %5s %5s %5s %4s %s\n",
///        tnamew-1, "-------------------", taccw, "----------", "-----", qnamew,
///        "--------------------", qaccw, "----------", "-----", "---------", "------",
///        "-----", "---", "---", "---------", "---------", "------", "-----", "-----",
///        "-----", "-----", "-----", "-----", "-----", "----", "---------------------") < 0) ...
pub fn header(w: &DomWidths) -> String {
    let mut s = String::new();

    // "#%*s %22s %40s %11s %11s %11s\n" with width tnamew+qnamew-1+15+taccw+qaccw
    let grp = w.tnamew + w.qnamew - 1 + 15 + w.taccw + w.qaccw;
    s.push('#');
    s.push_str(&format!(
        "{:grp$} {:>22} {:>40} {:>11} {:>11} {:>11}\n",
        "",
        "--- full sequence ---",
        "-------------- this domain -------------",
        "hmm coord",
        "ali coord",
        "env coord",
        grp = grp,
    ));

    // "#%-*s %-*s %5s %-*s %-*s %5s %9s %6s %5s %3s %3s %9s %9s %6s %5s %5s %5s %5s %5s %5s %5s %4s %s\n"
    s.push('#');
    s.push_str(&format!(
        "{:<tn$} {:<ta$} {:>5} {:<qn$} {:<qa$} {:>5} {:>9} {:>6} {:>5} {:>3} {:>3} {:>9} {:>9} \
         {:>6} {:>5} {:>5} {:>5} {:>5} {:>5} {:>5} {:>5} {:>4} {}\n",
        " target name",
        "accession",
        "tlen",
        "query name",
        "accession",
        "qlen",
        "E-value",
        "score",
        "bias",
        "#",
        "of",
        "c-Evalue",
        "i-Evalue",
        "score",
        "bias",
        "from",
        "to",
        "from",
        "to",
        "from",
        "to",
        "acc",
        "description of target",
        tn = w.tnamew - 1,
        ta = w.taccw,
        qn = w.qnamew,
        qa = w.qaccw,
    ));

    // "#%*s %*s %5s %*s %*s %5s %9s ... %s\n" (all right-justified)
    s.push('#');
    s.push_str(&format!(
        "{:>tn$} {:>ta$} {:>5} {:>qn$} {:>qa$} {:>5} {:>9} {:>6} {:>5} {:>3} {:>3} {:>9} {:>9} \
         {:>6} {:>5} {:>5} {:>5} {:>5} {:>5} {:>5} {:>5} {:>4} {}\n",
        "-------------------",
        "----------",
        "-----",
        "--------------------",
        "----------",
        "-----",
        "---------",
        "------",
        "-----",
        "---",
        "---",
        "---------",
        "---------",
        "------",
        "-----",
        "-----",
        "-----",
        "-----",
        "-----",
        "-----",
        "-----",
        "----",
        "---------------------",
        tn = w.tnamew - 1,
        ta = w.taccw,
        qn = w.qnamew,
        qa = w.qaccw,
    ));

    s
}

/// All reported-domain rows for one target hit (p7_tophits.c:1738-1787).
///
/// `tlen` is the target sequence length and `qlen` the model length: in hmmsearch
/// (`p7_SEARCH_SEQS`) `qlen = ad->M` and `tlen = ad->L` (p7_tophits.c:1753).
#[allow(clippy::too_many_arguments)]
pub fn format_rows(
    hit: &Hit,
    tlen: usize,
    qname: &str,
    qacc: &str,
    qlen: usize,
    pli: &Pipeline,
    w: &DomWidths,
) -> String {
    let (z, domz) = (pli.z, pli.domz);
    let nreported = hit
        .domains
        .iter()
        .filter(|d| pli.domain_reportable(d.bitscore, d.lnp))
        .count();

    let desc = if hit.desc.is_empty() { "-" } else { &hit.desc };
    let qacc_s = if qacc.is_empty() { "-" } else { qacc };
    let e_full = hit.lnp.exp() * z;

    let mut out = String::new();
    let mut nd = 0usize;
    for d in &hit.domains {
        if !pli.domain_reportable(d.bitscore, d.lnp) {
            continue;
        }
        nd += 1;
        // C, p7_tophits.c:1757-1786:
        //   if (fprintf(ofp, "%-*s %-*s %5d %-*s %-*s %5d %9.2g %6.1f %5.1f %3d %3d %9.2g %9.2g %6.1f %5.1f %5d %5d %5" PRId64 " %5" PRId64 " %5" PRId64 " %5" PRId64 " %4.2f %s\n",
        //     tnamew, th->hit[h]->name,
        //     taccw,  th->hit[h]->acc ? th->hit[h]->acc : "-",
        //     tlen,
        //     qnamew, qname,
        //     qaccw,  ( (qacc != NULL && qacc[0] != '\0') ? qacc : "-"),
        //     qlen,
        //     exp(th->hit[h]->lnP) * pli->Z,
        //     th->hit[h]->score,
        //     th->hit[h]->pre_score - th->hit[h]->score, /* bias correction */
        //     nd,
        //     th->hit[h]->nreported,
        //     exp(th->hit[h]->dcl[d].lnP) * pli->domZ,
        //     exp(th->hit[h]->dcl[d].lnP) * pli->Z,
        //     th->hit[h]->dcl[d].bitscore,
        //     th->hit[h]->dcl[d].dombias * eslCONST_LOG2R, /* NATS to BITS at last moment */
        //     th->hit[h]->dcl[d].ad->hmmfrom,
        //     th->hit[h]->dcl[d].ad->hmmto,
        //     th->hit[h]->dcl[d].ad->sqfrom,
        //     th->hit[h]->dcl[d].ad->sqto,
        //     th->hit[h]->dcl[d].ienv,
        //     th->hit[h]->dcl[d].jenv,
        //     (th->hit[h]->dcl[d].oasc / (1.0 + fabs((float) (th->hit[h]->dcl[d].jenv - th->hit[h]->dcl[d].ienv)))),
        //     (th->hit[h]->desc ?  th->hit[h]->desc : "-")) < 0) ...
        out.push_str(&format!(
            "{tname:<tnamew$} {tacc:<taccw$} {tlen:>5} {qname:<qnamew$} {qacc:<qaccw$} \
             {qlen:>5} {efull:>9} {sc:6.1} {bi:5.1} {nd:>3} {of:>3} {ce:>9} {ie:>9} \
             {dsc:6.1} {dbi:5.1} {hf:>5} {ht:>5} {af:>5} {at:>5} {ef:>5} {et:>5} \
             {acc:4.2} {desc}\n",
            tname = hit.name,
            tacc = if hit.acc.is_empty() { "-" } else { &hit.acc },
            tlen = tlen,
            qname = qname,
            qacc = qacc_s,
            qlen = qlen,
            tnamew = w.tnamew,
            taccw = w.taccw,
            qnamew = w.qnamew,
            qaccw = w.qaccw,
            efull = fmt_g2(e_full),
            sc = hit.score,
            bi = hit.bias(),
            nd = nd,
            of = nreported,
            // c-Evalue is conditional on the reported targets (domZ); i-Evalue is
            // independent, over the whole database (Z). p7_tophits.c:1773-1774.
            ce = fmt_g2(d.lnp.exp() * domz),
            ie = fmt_g2(d.lnp.exp() * z),
            dsc = d.bitscore,
            // NATS to BITS at last moment (p7_tophits.c:1777)
            dbi = dombias_bits(d.dombias),
            hf = d.hmm_from,
            ht = d.hmm_to,
            af = d.ali_from,
            at = d.ali_to,
            ef = d.ienv,
            et = d.jenv,
            acc = dom_acc(d),
            desc = desc,
        ));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hmmfile::P7Hmm;
    use crate::pipeline::Model;
    use crate::seqio::read_fasta;
    use crate::tblout::sort_hits;

    /// Every `--domtblout` data row must be byte-identical to HMMER 3.4's.
    #[test]
    fn globins_domtblout_rows_byte_identical() {
        let base = env!("CARGO_MANIFEST_DIR");
        let hmm = P7Hmm::read_all(&format!("{base}/testdata/globins4.hmm"))
            .unwrap()
            .pop()
            .unwrap();
        let qname = hmm.name.clone();
        let qlen = hmm.m;
        let model = Model::new(hmm);
        let seqs = read_fasta(&format!("{base}/testdata/globins45.fa")).unwrap();
        let mut pli = Pipeline::default();
        pli.z = seqs.len() as f64;

        let mut hits: Vec<Hit> = seqs.iter().filter_map(|s| model.search_one(s, &pli)).collect();
        sort_hits(&mut hits, &pli);
        let reported: Vec<&Hit> = hits.iter().collect();
        // domZ = number of reported targets (p7_tophits.c:963)
        pli.domz = reported.len() as f64;
        let w = DomWidths::compute(&qname, "", &reported);
        let mine: Vec<String> = reported
            .iter()
            .flat_map(|h| {
                format_rows(h, h.tlen, &qname, "", qlen, &pli, &w)
                    .lines()
                    .map(|l| l.trim_end().to_string())
                    .collect::<Vec<_>>()
            })
            .collect();

        // Reference rows from HMMER 3.4:
        //   hmmsearch --domtblout globins4.domtblout testdata/globins4.hmm \
        //             testdata/globins45.fa
        // (`#` comment/trailer lines stripped, trailing blanks trimmed).
        let golden =
            std::fs::read_to_string(format!("{base}/testdata/globins4.domtblout")).unwrap();
        let gold_rows: Vec<String> = golden
            .lines()
            .filter(|l| !l.starts_with('#'))
            .map(|l| l.trim_end().to_string())
            .collect();

        assert_eq!(mine.len(), gold_rows.len(), "row count");
        let mut nbad = 0;
        for (m, g) in mine.iter().zip(gold_rows.iter()) {
            if m != g {
                nbad += 1;
                if nbad <= 5 {
                    eprintln!("MINE : {m:?}");
                    eprintln!("GOLD : {g:?}");
                }
            }
        }
        assert_eq!(nbad, 0, "{nbad}/{} rows differ from golden", gold_rows.len());
    }

    /// The three `#` header lines must match C's byte-for-byte.
    #[test]
    fn globins_domtblout_header_byte_identical() {
        let base = env!("CARGO_MANIFEST_DIR");
        let golden =
            std::fs::read_to_string(format!("{base}/testdata/globins4.domtblout.hdr")).unwrap();
        let w = DomWidths {
            tnamew: 20,
            taccw: 10,
            qnamew: 20,
            qaccw: 10,
        };
        let mine: Vec<String> = header(&w).lines().map(|l| l.trim_end().to_string()).collect();
        let gold: Vec<String> = golden.lines().map(|l| l.trim_end().to_string()).collect();
        assert_eq!(mine, gold);
    }
}
