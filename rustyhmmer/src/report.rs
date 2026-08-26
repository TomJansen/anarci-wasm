//! The default human-readable `hmmsearch` report — what C writes to stdout or `-o`.
//!
//! Faithful transcription of HMMER 3.4:
//!
//! | Rust                  | C (hmmer-3.4)                                        |
//! |-----------------------|------------------------------------------------------|
//! | `query_header`        | `hmmsearch.c:487-489`                                |
//! | `targets`             | `p7_tophits_Targets()`     p7_tophits.c:1046-1180    |
//! | `domains`             | `p7_tophits_Domains()`     p7_tophits.c:1195-1420    |
//! | `pli_statistics`      | `p7_pli_Statistics()`      p7_pipeline.c:1096-1160   |
//!
//! Only the `long_targets == FALSE`, `mode == p7_SEARCH_SEQS` branches are
//! transcribed; the nhmmer/hmmscan branches are noted where elided.

use crate::alidisplay::AliDisplay;
use crate::forward::ForwardFilter;
use crate::hmmfile::P7Hmm;
use crate::pipeline::{dom_acc, dombias_bits, Hit};
use crate::pli::{Pipeline, Stats, ZSetBy};
use crate::seqio::Seq;
use crate::tblout::{domain_counts, fmt_g};

/// `p7_tophits_GetMaxNameLength()` / `GetMaxShownLength()` (p7_tophits.c:836-880).
/// Sequence targets from FASTA carry no accession, so the `--acc` variant falls
/// back to names.
/// The name `--acc` shows: the accession when the target has one, else the name.
///
/// C, p7_tophits.c:1186-1192:
/// ```text
///   if (pli->show_accessions)
///     {   /* the --acc option: report accessions rather than names if possible */
///       if (th->hit[h]->acc != NULL && th->hit[h]->acc[0] != '\0') showname = th->hit[h]->acc;
///       else                                                       showname = th->hit[h]->name;
///     }
///   else
///     showname = th->hit[h]->name;
/// ```
fn shown_name(h: &Hit, show_accessions: bool) -> &str {
    if show_accessions && !h.acc.is_empty() {
        &h.acc
    } else {
        &h.name
    }
}

/// `p7_tophits_GetMaxShownLength()` (p7_tophits.c:627-644): the same fallback, maximized.
/// C picks between this and `GetMaxNameLength` on `pli->show_accessions`
/// (p7_tophits.c:1133-1134).
fn max_shown_len(reported: &[&Hit], show_accessions: bool) -> usize {
    reported
        .iter()
        .map(|h| shown_name(h, show_accessions).len())
        .max()
        .unwrap_or(0)
}

/// `Query:` / `Accession:` / `Description:` — hmmsearch.c:487-489.
pub fn query_header(out: &mut String, hmm: &P7Hmm) {
    out.push_str(&format!("Query:       {}  [M={}]\n", hmm.name, hmm.m));
    if let Some(a) = &hmm.acc {
        out.push_str(&format!("Accession:   {a}\n"));
    }
    if let Some(d) = &hmm.desc {
        out.push_str(&format!("Description: {d}\n"));
    }
}

/// `p7_tophits_Targets()` — p7_tophits.c:1046-1180.
///
/// `textw` is C's `textw`: > 0 truncates descriptions, <= 0 means unlimited.
pub fn targets(out: &mut String, reported: &[&Hit], th: &[&Hit], pli: &Pipeline, textw: i32) {
    // namew is measured over C's `th` (p7_tophits.c:1055-1057); see tblout::Widths
    // for why the reported set stands in for it.
    let namew = 8.max(max_shown_len(th, pli.show_accessions));
    // "61 chars excluding desc is from the format: 2 + 22+2 +22+2 +8+2 +<name>+1"
    let descw = if textw > 0 {
        32.max(textw - namew as i32 - 61) as usize
    } else {
        0
    };

    // C, p7_tophits.c:1090-1099:
    //   fprintf(ofp, "Scores for complete sequence%s (score includes all domains):\n", pli->mode == p7_SEARCH_SEQS ? "s" : "");
    //   fprintf(ofp, "  %22s  %22s  %8s\n", " --- full sequence ---", " --- best 1 domain ---", "-#dom-");
    //   fprintf(ofp, "  %9s %6s %5s  %9s %6s %5s  %5s %2s  %-*s %s\n",
    //     "E-value", " score", " bias", "E-value", " score", " bias", "  exp",  "N", namew,
    //     (pli->mode == p7_SEARCH_SEQS ? "Sequence":"Model"), "Description");
    //   fprintf(ofp, "  %9s %6s %5s  %9s %6s %5s  %5s %2s  %-*s %s\n",
    //     "-------", "------", "-----", "-------", "------", "-----", " ----", "--", namew,
    //     "--------", "-----------");
    out.push_str("Scores for complete sequences (score includes all domains):\n");
    out.push_str(&format!(
        "  {:>22}  {:>22}  {:>8}\n",
        " --- full sequence ---", " --- best 1 domain ---", "-#dom-"
    ));
    out.push_str(&format!(
        "  {:>9} {:>6} {:>5}  {:>9} {:>6} {:>5}  {:>5} {:>2}  {:<namew$} {}\n",
        "E-value", " score", " bias", "E-value", " score", " bias", "  exp", "N", "Sequence",
        "Description"
    ));
    out.push_str(&format!(
        "  {:>9} {:>6} {:>5}  {:>9} {:>6} {:>5}  {:>5} {:>2}  {:<namew$} {}\n",
        "-------", "------", "-----", "-------", "------", "-----", " ----", "--", "--------",
        "-----------"
    ));

    let mut have_printed_incthresh = false;
    for h in reported {
        let bd = &h.domains[h.best];
        if !pli.target_includable(h.score, h.lnp) && !have_printed_incthresh {
            out.push_str("  ------ inclusion threshold ------\n");
            have_printed_incthresh = true;
        }
        // C, p7_tophits.c:1145-1156. `newness` is ' ' outside jackhmmer's iterative
        // mode (p7_IS_NEW / p7_IS_DROPPED are never set by hmmsearch):
        //   fprintf(ofp, "%c %9.2g %6.1f %5.1f  %9.2g %6.1f %5.1f  %5.1f %2d  %-*s ",
        //     newness,
        //     exp(th->hit[h]->lnP) * pli->Z,
        //     th->hit[h]->score,
        //     th->hit[h]->pre_score - th->hit[h]->score, /* bias correction */
        //     exp(th->hit[h]->dcl[d].lnP) * pli->Z,
        //     th->hit[h]->dcl[d].bitscore,
        //     eslCONST_LOG2R * th->hit[h]->dcl[d].dombias, /* convert NATS to BITS at last moment */
        //     th->hit[h]->nexpected,
        //     th->hit[h]->nreported,
        //     namew, showname);
        let (nrep, _ninc) = domain_counts(h, pli);
        out.push_str(&format!(
            "{} {:>9} {:6.1} {:5.1}  {:>9} {:6.1} {:5.1}  {:5.1} {:2}  {:<namew$} ",
            ' ',
            fmt_g(h.lnp.exp() * pli.z, 2),
            h.score,
            h.bias(),
            fmt_g(bd.lnp.exp() * pli.z, 2),
            bd.bitscore,
            dombias_bits(bd.dombias),
            h.nexpected,
            nrep,
            shown_name(h, pli.show_accessions),
        ));
        if textw > 0 {
            let d: String = h.desc.chars().take(descw).collect();
            out.push_str(&format!(" {d}\n"));
        } else {
            out.push_str(&format!(" {}\n", h.desc));
        }
    }

    if reported.is_empty() {
        out.push_str("\n   [No hits detected that satisfy reporting thresholds]\n");
    }
}

/// `p7_tophits_Domains()` — p7_tophits.c:1195-1420.
///
/// `seqs` supplies the target sequence for the alignment display, indexed by
/// `Hit::seq_idx`; C keeps the `P7_ALIDISPLAY` on the hit instead, built back when
/// the domain was scored.
#[allow(clippy::too_many_arguments)]
pub fn domains(
    out: &mut String,
    reported: &[&Hit],
    pli: &Pipeline,
    textw: i32,
    hmm: &P7Hmm,
    ff: &ForwardFilter,
    seqs: &[Seq],
) {
    out.push_str(&format!(
        "Domain annotation for each sequence{}:\n",
        if pli.show_alignments {
            " (and alignments)"
        } else {
            ""
        }
    ));

    for h in reported {
        // C, p7_tophits.c:1292-1300: `--acc` swaps the name here too, and `namew` becomes
        // that string's own length rather than a column width.
        let shown = shown_name(h, pli.show_accessions);
        let namew = shown.len();
        if textw > 0 {
            let descw = 32.max(textw - namew as i32 - 5) as usize;
            let d: String = h.desc.chars().take(descw).collect();
            out.push_str(&format!(">> {}  {}\n", shown, d));
        } else {
            out.push_str(&format!(">> {}  {}\n", shown, h.desc));
        }

        let (nrep, _) = domain_counts(h, pli);
        if nrep == 0 {
            out.push_str(
                "   [No individual domains that satisfy reporting thresholds (although complete target did)]\n\n",
            );
            continue;
        }

        out.push_str(&format!(
            " {:>3}   {:>6} {:>5} {:>9} {:>9} {:>7} {:>7} {:>2} {:>7} {:>7} {:>2} {:>7} {:>7} {:>2} {:>4}\n",
            "#", "score", "bias", "c-Evalue", "i-Evalue", "hmmfrom", "hmm to", "  ", "alifrom",
            "ali to", "  ", "envfrom", "env to", "  ", "acc"
        ));
        out.push_str(&format!(
            " {:>3}   {:>6} {:>5} {:>9} {:>9} {:>7} {:>7} {:>2} {:>7} {:>7} {:>2} {:>7} {:>7} {:>2} {:>4}\n",
            "---", "------", "-----", "---------", "---------", "-------", "-------", "  ",
            "-------", "-------", "  ", "-------", "-------", "  ", "----"
        ));

        let seq_included = pli.target_includable(h.score, h.lnp);
        let l = h.tlen as i64;
        let m = hmm.m as i32;

        let mut nd = 0;
        for d in &h.domains {
            if !pli.domain_reportable(d.bitscore, d.lnp) {
                continue;
            }
            nd += 1;
            let included = seq_included && pli.domain_includable(d.bitscore, d.lnp);
            // C, p7_tophits.c:1314-1345 (four fprintf calls building one row):
            //   fprintf(ofp, " %3d %c %6.1f %5.1f %9.2g %9.2g %7d %7d %c%c",
            //     nd, th->hit[h]->dcl[d].is_included ? '!' : '?',
            //     th->hit[h]->dcl[d].bitscore,
            //     th->hit[h]->dcl[d].dombias * eslCONST_LOG2R,
            //     exp(th->hit[h]->dcl[d].lnP) * pli->domZ,
            //     exp(th->hit[h]->dcl[d].lnP) * pli->Z,
            //     th->hit[h]->dcl[d].ad->hmmfrom, th->hit[h]->dcl[d].ad->hmmto,
            //     (th->hit[h]->dcl[d].ad->hmmfrom == 1) ? '[' : '.',
            //     (th->hit[h]->dcl[d].ad->hmmto   == th->hit[h]->dcl[d].ad->M ) ? ']' : '.');
            //   fprintf(ofp, " %7"PRId64" %7"PRId64" %c%c", ad->sqfrom, ad->sqto,
            //     (ad->sqfrom == 1) ? '[' : '.', (ad->sqto == ad->L) ? ']' : '.');
            //   fprintf(ofp, " %7"PRId64" %7"PRId64" %c%c", dcl[d].ienv, dcl[d].jenv,
            //     (dcl[d].ienv == 1) ? '[' : '.', (dcl[d].jenv == ad->L) ? ']' : '.');
            //   fprintf(ofp, " %4.2f\n", (dcl[d].oasc / (1.0 + fabs((float)(dcl[d].jenv - dcl[d].ienv)))));
            out.push_str(&format!(
                " {:>3} {} {:6.1} {:5.1} {:>9} {:>9} {:>7} {:>7} {}{}",
                nd,
                if included { '!' } else { '?' },
                d.bitscore,
                dombias_bits(d.dombias),
                fmt_g(d.lnp.exp() * pli.domz, 2),
                fmt_g(d.lnp.exp() * pli.z, 2),
                d.hmm_from,
                d.hmm_to,
                if d.hmm_from == 1 { '[' } else { '.' },
                if d.hmm_to == m { ']' } else { '.' },
            ));
            out.push_str(&format!(
                " {:>7} {:>7} {}{}",
                d.ali_from,
                d.ali_to,
                if d.ali_from == 1 { '[' } else { '.' },
                if d.ali_to == l { ']' } else { '.' },
            ));
            out.push_str(&format!(
                " {:>7} {:>7} {}{}",
                d.ienv,
                d.jenv,
                if d.ienv == 1 { '[' } else { '.' },
                if d.jenv == l { ']' } else { '.' },
            ));
            out.push_str(&format!(" {:4.2}\n", dom_acc(d)));
        }

        if pli.show_alignments {
            out.push_str("\n  Alignments for each domain:\n");
            let mut nd = 0;
            for d in &h.domains {
                if !pli.domain_reportable(d.bitscore, d.lnp) {
                    continue;
                }
                nd += 1;
                out.push_str(&format!("  == domain {nd}"));
                out.push_str(&format!("  score: {:.1} bits", d.bitscore));
                out.push_str(&format!(
                    ";  conditional E-value: {}\n",
                    fmt_g(d.lnp.exp() * pli.domz, 2)
                ));
                if let Some(sq) = seqs.get(h.seq_idx) {
                    if let Some(ad) = AliDisplay::create(&d.trace, hmm, ff, sq) {
                        // p7_alidisplay_Print(ofp, ad, 40, textw, pli)
                        ad.print(out, 40, textw, pli.show_accessions);
                    }
                }
                out.push('\n');
            }
        } else {
            out.push('\n');
        }
    }

    if reported.is_empty() {
        out.push_str("\n   [No targets detected that satisfy reporting thresholds]\n");
    }
}

/// `p7_pli_Statistics()` — p7_pipeline.c:1096-1160, `p7_SEARCH_SEQS` branch.
///
/// The trailing `# CPU time:` / `# Mc/sec:` lines that C emits from its stopwatch
/// are inherently run-specific; `elapsed` supplies them.
pub fn pli_statistics(out: &mut String, pli: &Pipeline, st: &Stats, elapsed: Option<f64>) {
    out.push_str("Internal pipeline statistics summary:\n");
    out.push_str("-------------------------------------\n");
    out.push_str(&format!(
        "Query model(s):              {:>15}  ({} nodes)\n",
        st.nmodels, st.nnodes
    ));
    out.push_str(&format!(
        "Target sequences:            {:>15}  ({} residues searched)\n",
        st.nseqs, st.nres
    ));
    let ntargets = st.nseqs as f64;

    let row = |name: &str, n: u64, f: f64| {
        format!(
            "{name}{:>15}  ({}); expected {:.1} ({})\n",
            n,
            fmt_g(n as f64 / ntargets, 6),
            f * ntargets,
            fmt_g(f, 6)
        )
    };
    // C, p7_pipeline.c:1138-1160:
    //   fprintf(ofp, "Passed MSV filter:           %15" PRId64 "  (%.6g); expected %.1f (%.6g)\n",
    //       pli->n_past_msv, (double) pli->n_past_msv / ntargets, pli->F1 * ntargets, pli->F1);
    //   ... same shape for bias/Vit/Fwd with F1/F2/F3 ...
    //   fprintf(ofp, "Initial search space (Z):    %15.0f  %s\n", pli->Z,
    //       pli->Z_setby    == p7_ZSETBY_OPTION ? "[as set by --Z on cmdline]"    : "[actual number of targets]");
    //   fprintf(ofp, "Domain search space  (domZ): %15.0f  %s\n", pli->domZ,
    //       pli->domZ_setby == p7_ZSETBY_OPTION ? "[as set by --domZ on cmdline]" : "[number of targets reported over threshold]");
    out.push_str(&row("Passed MSV filter:           ", st.n_past_msv, pli.f1));
    out.push_str(&row("Passed bias filter:          ", st.n_past_bias, pli.f1));
    out.push_str(&row("Passed Vit filter:           ", st.n_past_vit, pli.f2));
    out.push_str(&row("Passed Fwd filter:           ", st.n_past_fwd, pli.f3));

    out.push_str(&format!(
        "Initial search space (Z):    {:>15.0}  {}\n",
        pli.z,
        if pli.z_setby == ZSetBy::Option {
            "[as set by --Z on cmdline]"
        } else {
            "[actual number of targets]"
        }
    ));
    out.push_str(&format!(
        "Domain search space  (domZ): {:>15.0}  {}\n",
        pli.domz,
        if pli.domz_setby == ZSetBy::Option {
            "[as set by --domZ on cmdline]"
        } else {
            "[number of targets reported over threshold]"
        }
    ));

    if let Some(secs) = elapsed {
        // esl_stopwatch_Display(ofp, w, "# CPU time: ") — easel/esl_stopwatch.c.
        out.push_str(&format!(
            "# CPU time: {:.2}u {:.2}s 00:00:{:05.2} Elapsed: 00:00:{:05.2}\n",
            secs, 0.0, secs, secs
        ));
        let mcs = st.nres as f64 * st.nnodes as f64 / (secs * 1.0e6);
        out.push_str(&format!("# Mc/sec: {mcs:.2}\n"));
    }
    out.push_str("//\n");
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pipeline::Model;
    use crate::tblout::sort_hits;

    /// The whole human-readable report — banner, score table, per-domain table and
    /// the alignment display with its CS/consensus/match/sequence/PP lines — must be
    /// byte-identical to `hmmsearch`'s. The reference was produced by
    ///   hmmsearch -E 1e-10 testdata/fn3.hmm testdata/7LESS_DROME.fa
    /// with only the run-specific lines removed (`# CPU time:`, `# Mc/sec:`, and the
    /// two input-path lines).
    #[test]
    fn fn3_7less_report_byte_identical() {
        let base = env!("CARGO_MANIFEST_DIR");
        let hmm = P7Hmm::read_all(&format!("{base}/testdata/fn3.hmm"))
            .unwrap()
            .pop()
            .unwrap();
        let seqs = crate::seqio::read_fasta(&format!("{base}/testdata/7LESS_DROME.fa")).unwrap();

        let mut pli = Pipeline::default();
        pli.e = 1e-10;
        pli.z = seqs.len() as f64;
        let model = Model::new(hmm);

        let mut hits: Vec<Hit> = Vec::new();
        let mut idx = 0usize;
        let mut st = Stats {
            nmodels: 1,
            nnodes: model.hmm.m as u64,
            ..Default::default()
        };
        for s in &seqs {
            let (h, one) = model.search_one_counted(s, &pli);
            st.merge(&one);
            if let Some(mut h) = h {
                h.seq_idx = idx;
                hits.push(h);
            }
            idx += 1;
        }
        sort_hits(&mut hits, &pli);
        let reported: Vec<&Hit> = hits
            .iter()
            .filter(|h| pli.target_reportable(h.score, h.lnp))
            .collect();
        pli.domz = reported.len() as f64;

        let mut got = String::new();
        got.push_str("# hmmsearch :: search profile(s) against a sequence database\n");
        got.push_str("# HMMER 3.4 (Aug 2023); http://hmmer.org/\n");
        got.push_str("# Copyright (C) 2023 Howard Hughes Medical Institute.\n");
        got.push_str("# Freely distributed under the BSD open source license.\n");
        got.push_str("# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
        got.push_str("# sequence reporting threshold:    E-value <= 1e-10\n");
        got.push_str("# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");
        query_header(&mut got, &model.hmm);
        targets(&mut got, &reported, &reported, &pli, 120);
        got.push_str("\n\n");
        domains(&mut got, &reported, &pli, 120, &model.hmm, &model.fwd, &seqs);
        got.push_str("\n\n");
        pli_statistics(&mut got, &pli, &st, None);
        got.push_str("[ok]\n");

        let want = std::fs::read_to_string(format!("{base}/testdata/fn3_7LESS.out")).unwrap();
        let g: Vec<&str> = got.lines().map(|l| l.trim_end()).collect();
        let w: Vec<&str> = want.lines().map(|l| l.trim_end()).collect();
        let mut nbad = 0;
        for (i, (a, b)) in g.iter().zip(w.iter()).enumerate() {
            if a != b {
                nbad += 1;
                if nbad <= 6 {
                    eprintln!("line {}:\n  MINE: {a:?}\n  GOLD: {b:?}", i + 1);
                }
            }
        }
        assert_eq!(g.len(), w.len(), "line count");
        assert_eq!(nbad, 0, "{nbad}/{} report lines differ", w.len());
    }
}
