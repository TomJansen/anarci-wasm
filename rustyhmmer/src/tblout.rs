


use crate::pipeline::{dombias_bits, Hit};
use crate::pli::Pipeline;





pub fn fmt_g2(x: f64) -> String {
    fmt_g(x, 2)
}

/// C's `%.*g`. Chooses `%e` or `%f` by exponent, then strips trailing zeros.
pub fn fmt_g(x: f64, prec: i32) -> String {
    let prec = prec.max(1);
    // glibc's `%g` spells the non-finite values out, with the sign bit honoured. This is
    // reachable: an empty target database makes p7_pli_Statistics() divide 0 by 0, and
    // hmmsearch prints "Passed MSV filter: 0  (-nan)".
    if !x.is_finite() {
        let body = if x.is_nan() { "nan" } else { "inf" };
        return if x.is_sign_negative() {
            format!("-{body}")
        } else {
            body.to_string()
        };
    }
    if x == 0.0 {
        return "0".to_string();
    }
    let neg = x < 0.0;
    let ax = x.abs();

    let sci = format!("{:.*e}", (prec - 1) as usize, ax);
    let (mant, exp_str) = sci.split_once('e').unwrap();
    let exp: i32 = exp_str.parse().unwrap(); 

    let body = if exp >= -4 && exp < prec {
        
        
        let dec = ((prec - 1) - exp).max(0) as usize;
        let f = format!("{:.*}", dec, ax);
        if f.contains('.') {
            f.trim_end_matches('0').trim_end_matches('.').to_string()
        } else {
            f
        }
    } else {
        
        let mant = if mant.contains('.') {
            mant.trim_end_matches('0').trim_end_matches('.')
        } else {
            mant
        };
        let (sign, digits) = match exp_str.strip_prefix('-') {
            Some(d) => ("-", d),
            None => ("+", exp_str.trim_start_matches('+')),
        };
        format!("{}e{}{:0>2}", mant, sign, digits)
    };
    if neg {
        format!("-{body}")
    } else {
        body
    }
}





#[derive(Clone, Copy)]
pub struct Widths {
    pub tnamew: usize,
    pub taccw: usize,
    pub qnamew: usize,
    pub qaccw: usize,
}

impl Widths {
    
    
    
    pub fn compute(qname: &str, qacc: &str, reported: &[&Hit]) -> Self {
        let tnamew = 20.max(reported.iter().map(|h| h.name.len()).max().unwrap_or(0));
        // C, p7_tophits.c:1615: `taccw = ESL_MAX(10, p7_tophits_GetMaxAccessionLength(th))`,
        // and `GetMaxAccessionLength` returns 0 when no hit has an accession -- which is
        // every FASTA run, since FASTA carries none. EMBL/UniProt and GenBank/DDBJ do.
        let taccw = 10.max(reported.iter().map(|h| h.acc.len()).max().unwrap_or(0));
        let qnamew = 20.max(qname.len());
        let qaccw = if qacc.is_empty() { 10 } else { 10.max(qacc.len()) };
        Widths { tnamew, taccw, qnamew, qaccw }
    }
}




pub fn header(w: &Widths) -> String {
    let grp = w.tnamew + w.qnamew + w.taccw + w.qaccw + 2;
    let mut s = String::new();
    
    s.push('#');
    s.push_str(&format!(
        "{:grp$} {:>22} {:>22} {:>33}\n",
        "", "--- full sequence ----", "--- best 1 domain ----", "--- domain number estimation ----",
        grp = grp,
    ));
    
    s.push('#');
    s.push_str(&format!(
        "{:<tn$} {:<ta$} {:<qn$} {:<qa$} {:>9} {:>6} {:>5} {:>9} {:>6} {:>5} {:>5} {:>3} {:>3} {:>3} {:>3} {:>3} {:>3} {:>3} {}\n",
        " target name", "accession", "query name", "accession",
        "  E-value", " score", " bias", "  E-value", " score", " bias",
        "exp", "reg", "clu", " ov", "env", "dom", "rep", "inc", "description of target",
        tn = w.tnamew - 1, ta = w.taccw, qn = w.qnamew, qa = w.qaccw,
    ));
    
    // C, p7_tophits.c:1637-1638:
    //   if (fprintf(ofp, "#%*s %*s %*s %*s %9s %6s %5s %9s %6s %5s %5s %3s %3s %3s %3s %3s %3s %3s %s\n",
    //     tnamew-1, "-------------------", taccw, "----------", qnamew, "--------------------", qaccw, "----------", "---------", "------", "-----", "---------", "------", "-----", "---", "---", "---", "---", "---", "---", "---", "---", "---------------------") < 0)
    //
    // Count the dash strings, not the conversions: there are six for the two
    // score/bias triples and then *eight* `"---"` for the eight count columns. So the
    // `exp` column's `%5s` is fed `"---"`, printing `  ---`, not five dashes. (The
    // header's own `exp` label is `"exp"` in `%5s` for the same reason.)
    s.push('#');
    s.push_str(&format!(
        "{:>tn$} {:>ta$} {:>qn$} {:>qa$} {:>9} {:>6} {:>5} {:>9} {:>6} {:>5} {:>5} {:>3} {:>3} {:>3} {:>3} {:>3} {:>3} {:>3} {}\n",
        "-------------------", "----------", "--------------------", "----------",
        "---------", "------", "-----", "---------", "------", "-----", "---",
        "---", "---", "---", "---", "---", "---", "---", "---------------------",
        tn = w.tnamew - 1, ta = w.taccw, qn = w.qnamew, qa = w.qaccw,
    ));
    s
}






pub fn format_row(hit: &Hit, qname: &str, qacc: &str, pli: &Pipeline, w: &Widths) -> String {
    let z = pli.z;
    
    
    let e_full = hit.lnp.exp() * z;
    let bd = &hit.domains[hit.best];
    let e_dom = bd.lnp.exp() * z;

    
    
    
    
    // p7_tophits_Threshold() (p7_tophits.c:968-999): a domain can only be reported
    // if its sequence is, and only be included if its sequence is included too.
    let (nrep, ninc) = domain_counts(hit, pli);

    let desc = if hit.desc.is_empty() { "-" } else { &hit.desc };
    let qacc = if qacc.is_empty() { "-" } else { qacc };

    format!(
        "{tname:<tnamew$} {tacc:<taccw$} {qname:<qnamew$} {qacc:<qaccw$} {efull:>9} {sc:6.1} {bi:5.1} \
         {edom:>9} {dsc:6.1} {dbi:5.1} {exp:5.1} {reg:3} {clu:3} {ov:3} {env:3} {dom:3} {rep:3} {inc:3} {desc}",
        tname = hit.name,
        // `th->hit[h]->acc ? th->hit[h]->acc : "-"`
        tacc = if hit.acc.is_empty() { "-" } else { &hit.acc },
        qname = qname,
        qacc = qacc,
        tnamew = w.tnamew,
        taccw = w.taccw,
        qnamew = w.qnamew,
        qaccw = w.qaccw,
        efull = fmt_g2(e_full),
        sc = hit.score,
        bi = hit.bias(),
        edom = fmt_g2(e_dom),
        dsc = bd.bitscore,
        dbi = dombias_bits(bd.dombias),
        exp = hit.nexpected,
        reg = hit.nregions,
        clu = hit.nclustered,
        ov = hit.noverlaps,
        env = hit.nenvelopes,
        dom = hit.ndom,
        rep = nrep,
        inc = ninc,
        desc = desc,
    )
}



/// `nreported` / `nincluded` for one hit, per `p7_tophits_Threshold()`
/// (p7_tophits.c:968-999). The hit is assumed reportable (only reportable targets
/// reach the output stage; p7_pipeline.c:838).
pub fn domain_counts(hit: &Hit, pli: &Pipeline) -> (i32, i32) {
    let seq_included = pli.target_includable(hit.score, hit.lnp);
    let mut nrep = 0;
    let mut ninc = 0;
    for d in &hit.domains {
        if pli.domain_reportable(d.bitscore, d.lnp) {
            nrep += 1;
        }
        if seq_included && pli.domain_includable(d.bitscore, d.lnp) {
            ninc += 1;
        }
    }
    (nrep, ninc)
}

/// C's `th`: the hit list the pipeline accumulated, which is what the output column
/// widths (`namew`, `tnamew`) are measured over — not the finally-reported subset.
///
/// `p7_Pipeline()` creates a hit for every target that passes
/// `p7_pli_TargetReportable()` (p7_pipeline.c:838), but `p7_pli_NewSeq()` has just
/// set `pli->Z = pli->nseqs` (p7_pipeline.c:576-582), so the k-th target is tested
/// against a *running lower bound* Z = k rather than the final Z. A target early in
/// the file can therefore enter `th` on a lower-bound E-value and then fail the real
/// threshold in `p7_tophits_Threshold()`.
///
/// Because the bound depends only on the target's position in the file, this is
/// deterministic and independent of the order rayon happens to evaluate targets in.
/// It reproduces C's `--cpu 1` behaviour exactly; C's worker threads count their own
/// sequences, so their bound is looser, but that was verified not to change the
/// output here (identical for `--cpu` 1, 2 and 8 on bac120 x E. coli).
pub fn hitlist<'a>(hits: &'a [Hit], pli: &Pipeline) -> Vec<&'a Hit> {
    let mut during = pli.clone();
    hits.iter()
        .filter(|h| {
            if during.z_setby == crate::pli::ZSetBy::NTargets {
                // seq_idx is 0-based; pli->nseqs is 1-based at this point.
                during.z = (h.seq_idx as f64) + 1.0;
            }
            during.target_reportable(h.score, h.lnp)
        })
        .collect()
}

/// Sort the target hit list the way C's `p7_tophits_SortBySortkey()` does.
///
/// C, p7_tophits.c:310-331 (`hit_sorter_by_sortkey`, the comparator `qsort` is handed
/// at p7_tophits.c:399):
/// ```text
///   static int
///   hit_sorter_by_sortkey(const void *vh1, const void *vh2)
///   {
///     P7_HIT *h1 = *((P7_HIT **) vh1);  /* don't ask. don't change. Don't Panic. */
///     P7_HIT *h2 = *((P7_HIT **) vh2);
///     int     c;
///
///     if      (h1->sortkey < h2->sortkey) return  1;
///     else if (h1->sortkey > h2->sortkey) return -1;
///     else {
///       if ( (c = strcmp(h1->name, h2->name)) != 0) return c;
///
///       /* if on different strand, the positive strand goes first, else use position */
///       int dir1 = (h1->dcl[0].iali < h1->dcl[0].jali ? 1 : -1);
///       int dir2 = (h2->dcl[0].iali < h2->dcl[0].jali ? 1 : -1);
///       if (dir1 != dir2) return dir2;
///       else {
///         if     (h1->dcl[0].iali > h2->dcl[0].iali) return  1;
///         else if(h1->dcl[0].iali < h2->dcl[0].iali) return -1;
///         else                                       return  0;
///       }
///     }
///   }
/// ```
///
/// The sortkey itself is `hit->sortkey = pli->inc_by_E ? -lnP : seq_score`
/// (p7_pipeline.c:862), so descending `-lnP` is ascending `lnP`.
///
/// The tie-break is not cosmetic. Real target sets produce exact score ties — a
/// concatenated set of seven proteomes has several per profile — and without the
/// `strcmp` step a stable sort leaves them in target-file order, which is a different
/// order than C's and shows up as transposed rows in `--tblout`/`--domtblout`.
pub fn sort_hits(hits: &mut [Hit], pli: &Pipeline) {
    // `hit->sortkey = pli->inc_by_E ? -lnP : seq_score` -- bigger is better, so the
    // comparator sorts *descending* on it.
    let key = |h: &Hit| -> f64 {
        if pli.inc_by_e {
            -h.lnp
        } else {
            h.score as f64
        }
    };
    // `dir = (iali < jali ? 1 : -1)`. hmmsearch's envelopes always run left to right,
    // so both dirs are 1 and the branch never fires; nhmmer's reverse-strand hits are
    // what it exists for.
    let dir = |h: &Hit| -> i32 {
        let d = h.domains.first();
        match d {
            Some(d) if d.ali_from < d.ali_to => 1,
            Some(_) => -1,
            None => 1,
        }
    };
    let iali = |h: &Hit| -> i64 { h.domains.first().map(|d| d.ali_from).unwrap_or(0) };

    hits.sort_by(|a, b| {
        key(b)
            .partial_cmp(&key(a))
            .unwrap()
            // `strcmp` compares as `unsigned char`; Rust's byte-slice ordering matches.
            .then_with(|| a.name.as_bytes().cmp(b.name.as_bytes()))
            // `if (dir1 != dir2) return dir2;` -- positive strand first.
            .then_with(|| dir(b).cmp(&dir(a)))
            .then_with(|| iali(a).cmp(&iali(b)))
    });
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hmmfile::P7Hmm;
    use crate::pipeline::Model;
    use crate::seqio::read_fasta;

    /// glibc's `%g` for the non-finite values, which an empty target database reaches:
    /// C prints "Passed MSV filter: 0  (-nan)" there, because `0/0` on x86 sets the
    /// NaN's sign bit.
    #[test]
    fn fmt_g_spells_out_non_finite() {
        let zero_over_zero = 0.0f64 / std::hint::black_box(0.0f64);
        assert!(zero_over_zero.is_sign_negative());
        assert_eq!(fmt_g(zero_over_zero, 6), "-nan");
        assert_eq!(fmt_g(f64::NAN, 6), "nan");
        assert_eq!(fmt_g(f64::INFINITY, 6), "inf");
        assert_eq!(fmt_g(f64::NEG_INFINITY, 6), "-inf");
    }

    /// The three `#` header lines, against C's own output. The dash rule in particular
    /// is easy to get wrong: C feeds `"---"` to the `exp` column's `%5s`, so it prints
    /// `  ---`, and only the data rows were covered before.
    #[test]
    fn globins_tblout_header_byte_identical() {
        let base = env!("CARGO_MANIFEST_DIR");
        let golden =
            std::fs::read_to_string(format!("{base}/testdata/globins4.tblout.hdr")).unwrap();
        let w = Widths { tnamew: 20, taccw: 10, qnamew: 20, qaccw: 10 };
        let mine: Vec<String> = header(&w).lines().map(|l| l.trim_end().to_string()).collect();
        let gold: Vec<String> = golden.lines().map(|l| l.trim_end().to_string()).collect();
        assert_eq!(mine, gold);
    }

    #[test]
    fn globins_tblout_rows_byte_identical() {
        let base = env!("CARGO_MANIFEST_DIR");
        let hmm = P7Hmm::read_all(&format!("{base}/testdata/globins4.hmm"))
            .unwrap()
            .pop()
            .unwrap();
        let qname = hmm.name.clone();
        let model = Model::new(hmm);
        let seqs = read_fasta(&format!("{base}/testdata/globins45.fa")).unwrap();
        let mut pli = Pipeline::default();
        pli.z = seqs.len() as f64;
        pli.domz = pli.z;

        let mut hits: Vec<Hit> = seqs.iter().filter_map(|s| model.search_one(s, &pli)).collect();
        sort_hits(&mut hits, &pli);
        let reported: Vec<&Hit> = hits.iter().collect();
        let w = Widths::compute(&qname, "", &reported);
        let mine: Vec<String> = hits
            .iter()
            .map(|h| format_row(h, &qname, "", &pli, &w))
            .collect();

        


        // Reference rows from HMMER 3.4:
        //   hmmsearch --tblout globins4.tblout testdata/globins4.hmm testdata/globins45.fa
        // (`#` comment/trailer lines stripped, trailing blanks trimmed).
        let golden = std::fs::read_to_string(format!("{base}/testdata/globins4.tblout")).unwrap();
        let gold_rows: Vec<String> = golden
            .lines()
            .filter(|l| !l.starts_with('#'))
            .map(|l| l.trim_end().to_string())
            .collect();

        assert_eq!(mine.len(), gold_rows.len(), "row count");
        let mut nbad = 0;
        for (m, g) in mine.iter().zip(gold_rows.iter()) {
            if m.trim_end() != g.trim_end() {
                nbad += 1;
                if nbad <= 5 {
                    eprintln!("MINE : {:?}", m.trim_end());
                    eprintln!("GOLD : {:?}", g);
                }
            }
        }
        assert_eq!(nbad, 0, "{nbad}/{} rows differ from golden", gold_rows.len());
    }
}
