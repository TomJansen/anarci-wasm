//! `--pfamtblout`: the Pfam-style tabular format.
//!
//! C, p7_tophits.c:`p7_tophits_TabularXfam`. `hmmsearch` calls it once per query
//! (hmmsearch.c:549), inside the query loop, and then writes the usual tabular tail once
//! at the end (hmmsearch.c:599) -- so the two blocks below repeat for every query.
//!
//! Only the `!pli->long_targets` branch is implemented; the other is nhmmer's.

use crate::pipeline::Hit;
use crate::pli::Pipeline;
use crate::tblout::fmt_g2;

/// `tnamew`, the only width the sequence-search branch uses.
///
/// C, p7_tophits.c: `int tnamew = ESL_MAX(20, p7_tophits_GetMaxNameLength(th));`
/// measured over `th`, the accumulated hit list -- the same basis as `--tblout`'s, not
/// the finally-reported subset. (`taccw` and `qnamew` are computed there too but only
/// the long-target branch prints them.)
pub fn tnamew(th: &[&Hit]) -> usize {
    20.max(th.iter().map(|h| h.name.len()).max().unwrap_or(0))
}

/// One query's two blocks.
///
/// ```text
///   if (fprintf(ofp, "# Sequence scores\n# ---------------\n#\n") < 0) ...
///   if (fprintf(ofp, "# %-*s %6s %9s %3s %5s %5s    %s\n",
///   tnamew-1, "name",  " bits", "  E-value", "n",  "exp", " bias", "description") < 0) ...
///   if (fprintf(ofp, "# %*s %6s %9s %3s %5s %5s    %s\n",
///   tnamew-1, "-------------------",  "------", "---------","---", "-----",  "-----", "---------------------") < 0) ...
///
///   for (h = 0; h < th->N; h++)
///   {
///     if (th->hit[h]->flags & p7_IS_REPORTED)
///     {
///       if (fprintf(ofp, "%-*s  %6.1f %9.2g %3d %5.1f %5.1f    %s\n",
///       tnamew, th->hit[h]->name,
///       th->hit[h]->score,
///       exp(th->hit[h]->lnP) * pli->Z,
///       th->hit[h]->ndom,
///       th->hit[h]->nexpected,
///       th->hit[h]->pre_score - th->hit[h]->score, /* bias correction */
///       (th->hit[h]->desc == NULL ? "-" : th->hit[h]->desc)) < 0) ...
///     }
///   }
///   if (fprintf(ofp, "\n") < 0) ...
/// ```
///
/// then, after re-sorting the reported domains as pseudo-hits:
///
/// ```text
///   if (fprintf(ofp, "# Domain scores\n# -------------\n#\n") < 0) ...
///   if (fprintf(ofp, "# %-*s %6s %9s %5s %5s %6s %6s %6s %6s %6s %6s     %s\n",
///   tnamew-1, " name",  "bits", "E-value", "hit", "bias",      "env-st",  "env-en",  "ali-st",  "ali-en",  "hmm-st",  "hmm-en",   "description") < 0) ...
///   if (fprintf(ofp, "# %*s %6s %9s %5s %5s %6s %6s %6s %6s %6s %6s      %s\n",
///   tnamew-1, "-------------------",  "------", "---------", "-----", "-----", "------", "------", "------", "------", "------", "------", "---------------------") < 0) ...
///
///   for (h = 0; h < domHitlist->N; h++)
///   {
///     domhit = domHitlist->hit[h];
///     if (fprintf(ofp, "%-*s  %6.1f %9.2g %5d %5.1f %6" PRId64 " %6" PRId64 " %6" PRId64 " %6" PRId64 " %6d %6d     %s\n",
///           tnamew, domHitlist->hit[h]->name,
///           domhit->dcl[0].bitscore,
///           exp(domhit->dcl[0].lnP) * pli->Z, //i-Evalue
///           domhit->ndom,
///           domhit->dcl[0].dombias * eslCONST_LOG2R, // NATS to BITS at last moment
///           domhit->dcl[0].ienv,
///           domhit->dcl[0].jenv,
///           domhit->dcl[0].ad->sqfrom,
///           domhit->dcl[0].ad->sqto,
///           domhit->dcl[0].ad->hmmfrom,
///           domhit->dcl[0].ad->hmmto,
///           (domhit->desc ?  domhit->desc : "-")) < 0) ...
///   }
/// ```
///
/// Note the two dash rules are *not* the same: the sequence one has four spaces before
/// the description column, the domain one has six.
///
/// There is no blank line after the domain rows, so the next query's `# Sequence scores`
/// follows immediately.
pub fn query_blocks(out: &mut String, reported: &[&Hit], pli: &Pipeline, w: usize) {
    let tn = w - 1;

    out.push_str("# Sequence scores\n# ---------------\n#\n");
    out.push_str(&format!(
        "# {:<tn$} {:>6} {:>9} {:>3} {:>5} {:>5}    {}\n",
        "name", " bits", "  E-value", "n", "exp", " bias", "description",
    ));
    out.push_str(&format!(
        "# {:>tn$} {:>6} {:>9} {:>3} {:>5} {:>5}    {}\n",
        "-------------------", "------", "---------", "---", "-----", "-----",
        "---------------------",
    ));
    for h in reported {
        out.push_str(&format!(
            "{name:<w$}  {sc:6.1} {ev:>9} {n:3} {exp:5.1} {bias:5.1}    {desc}\n",
            name = h.name,
            sc = h.score,
            ev = fmt_g2(h.lnp.exp() * pli.z),
            n = h.ndom,
            exp = h.nexpected,
            bias = h.bias(),
            desc = if h.desc.is_empty() { "-" } else { &h.desc },
        ));
    }
    out.push('\n');

    // "Need to sort the domains. One way to do this is to re-use the hit sorting
    // machinery, so we create one 'hit' for each domain, then hand it off to the sorter"
    // -- and `domhit->ndom` is overloaded to carry the domain's 1-based ordinal within
    // its original hit, which is the `hit` column.
    struct Row<'a> {
        name: &'a str,
        desc: &'a str,
        ordinal: i32,
        dom: &'a crate::pipeline::DomHit,
        sortkey: f64,
    }
    let mut rows: Vec<Row> = Vec::new();
    for h in reported {
        let mut ordinal = 0;
        for d in &h.domains {
            if pli.domain_reportable(d.bitscore, d.lnp) {
                ordinal += 1;
                rows.push(Row {
                    name: &h.name,
                    desc: if h.desc.is_empty() { "-" } else { &h.desc },
                    ordinal,
                    dom: d,
                    // `pli->inc_by_E ? -1.0 * dcl[d].lnP : dcl[d].bitscore`
                    sortkey: if pli.inc_by_e {
                        -d.lnp
                    } else {
                        d.bitscore as f64
                    },
                });
            }
        }
    }
    // `p7_tophits_SortBySortkey(domHitlist)` -- the same comparator the target list
    // uses, so the same tie-break applies: descending sortkey, then `strcmp` of the
    // name, then strand, then `iali`. See `tblout::sort_hits`.
    rows.sort_by(|a, b| {
        let dir = |r: &Row| -> i32 {
            if r.dom.ali_from < r.dom.ali_to {
                1
            } else {
                -1
            }
        };
        b.sortkey
            .partial_cmp(&a.sortkey)
            .unwrap()
            .then_with(|| a.name.as_bytes().cmp(b.name.as_bytes()))
            .then_with(|| dir(b).cmp(&dir(a)))
            .then_with(|| a.dom.ali_from.cmp(&b.dom.ali_from))
    });

    out.push_str("# Domain scores\n# -------------\n#\n");
    out.push_str(&format!(
        "# {:<tn$} {:>6} {:>9} {:>5} {:>5} {:>6} {:>6} {:>6} {:>6} {:>6} {:>6}     {}\n",
        " name", "bits", "E-value", "hit", "bias", "env-st", "env-en", "ali-st", "ali-en",
        "hmm-st", "hmm-en", "description",
    ));
    out.push_str(&format!(
        "# {:>tn$} {:>6} {:>9} {:>5} {:>5} {:>6} {:>6} {:>6} {:>6} {:>6} {:>6}      {}\n",
        "-------------------", "------", "---------", "-----", "-----", "------", "------",
        "------", "------", "------", "------", "---------------------",
    ));
    for r in &rows {
        out.push_str(&format!(
            "{name:<w$}  {sc:6.1} {ev:>9} {hit:5} {bias:5.1} {es:6} {ee:6} {as_:6} {ae:6} {hs:6} {he:6}     {desc}\n",
            name = r.name,
            sc = r.dom.bitscore,
            ev = fmt_g2(r.dom.lnp.exp() * pli.z),
            hit = r.ordinal,
            bias = crate::pipeline::dombias_bits(r.dom.dombias),
            es = r.dom.ienv,
            ee = r.dom.jenv,
            as_ = r.dom.ali_from,
            ae = r.dom.ali_to,
            hs = r.dom.hmm_from,
            he = r.dom.hmm_to,
            desc = r.desc,
        ));
    }
}
