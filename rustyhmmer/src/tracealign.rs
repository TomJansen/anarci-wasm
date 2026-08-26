//! Turning traces back into a multiple alignment: `tracealign.c`'s text path.
//!
//! Reached from `hmmsearch -A` via `p7_tophits_Alignment()`, which back-converts every
//! *included* domain's alidisplay into a faux sequence + trace and hands the array to
//! `p7_tracealign_Seqs()`. `hmmsearch` passes `p7_ALL_CONSENSUS_COLS` and not
//! `p7_DIGITIZE` or `p7_TRIM`, so only that combination is implemented; the `p7T_X`
//! fragment states `make_text_msa` handles cannot occur in a profile trace.

use crate::alidisplay::AliDisplay;
use crate::msa::Msa;
use crate::optacc::TracePos;

use crate::dp::{ST_C, ST_D, ST_E, ST_I, ST_M, ST_N};

/// `esl_abc_DigitizeSymbol`: one character to its code, with anything unrecognized
/// becoming X, matching [`crate::alphabet::AminoMap::digitize`]'s per-character rule.
fn digitize_one(map: &crate::alphabet::AminoMap, c: u8) -> u8 {
    let code = map.code(c);
    if code == crate::alphabet::ILLEGAL {
        crate::alphabet::AMINO_X
    } else {
        code
    }
}

/// One back-converted domain: the faux subsequence's name/desc/acc plus its trace.
pub struct FauxSeq {
    pub name: String,
    pub desc: Option<String>,
    pub acc: Option<String>,
    /// `(state, k, i, pp)` as [`crate::optacc::TracePos`], with `i` indexing the faux
    /// subsequence 1..subL rather than the original target.
    pub tr: Vec<TracePos>,
    /// `sq->dsq[1..=subL]`, digital.
    pub dsq: Vec<u8>,
}

/// C, p7_alidisplay.c:`p7_alidisplay_Backconvert`: rebuild a subsequence and a trace from
/// the alignment display.
///
/// ```text
///   k = ad->hmmfrom - 1;   // -1 so the first M causes k to advance to <hmmfrom>.
///   i = 1;                 //    ... which assumes <ad> always starts with M; currently true, all alis are local alis.
///   for (a = 0; a < ad->N; a++)
///     {
///       /* here, anything that could appear in the original input
///        * sequence, as opposed to an alignment gap, needs to be
///        * reconstructed.  So, do not test for IsResidue(), because that
///        * will fail to reconstruct * chars. use ! IsGap() instead.
///        * [xref iss#135]
///        */
///       if   (! esl_abc_CIsGap(abc, ad->model[a])) { k++; s = (! esl_abc_CIsGap(abc, ad->aseq[a]) ? p7T_M : p7T_D); }
///       else s = p7T_I;
///       ...
///       switch (s) {
///       case p7T_M: sq->dsq[i] = esl_abc_DigitizeSymbol(abc, ad->aseq[a]); i++; break;
///       case p7T_I: sq->dsq[i] = esl_abc_DigitizeSymbol(abc, ad->aseq[a]); i++; break;
///       case p7T_D:                                                             break;
///       }
///     }
/// ```
///
/// The trace is bracketed by S-N-B ... E-C-T, and the name is
/// `<sqname>/<sqfrom>-<sqto>` with the description `[subseq from] <desc-or-name>`.
///
/// The three sanity checks C makes (`tr->N == ad->N + 6`, `k == ad->hmmto`,
/// `i == subL+1`) are kept as debug assertions: they check this transcription, not the
/// input.
pub fn backconvert(ad: &AliDisplay) -> FauxSeq {
    let map = crate::alphabet::AminoMap::new();
    // Easel's `esl_abc_CIsGap` is `strchr(a->gapsym, c)`, and the amino alphabet's
    // `gapsym` is "-_.~" -- so `*` is *not* a gap here. That is the whole point of
    // iss#135, which C's comment flags: testing `IsResidue` instead would drop `*`.
    let is_gap = |c: u8| matches!(c, b'-' | b'_' | b'.' | b'~');

    let mut tr: Vec<TracePos> = Vec::with_capacity(ad.n + 6);
    let mut dsq: Vec<u8> = vec![255]; // dsq[0] = eslDSQ_SENTINEL
    tr.push((crate::dp::ST_S, 0, 0, 0.0));
    tr.push((ST_N, 0, 0, 0.0));
    tr.push((crate::dp::ST_B, 0, 0, 0.0));

    let mut k = ad.hmmfrom - 1;
    let mut i: i32 = 1;
    for a in 0..ad.n {
        let s = if !is_gap(ad.model[a]) {
            k += 1;
            if !is_gap(ad.aseq[a]) {
                ST_M
            } else {
                ST_D
            }
        } else {
            ST_I
        };
        let pp = ad
            .ppline
            .as_ref()
            .map(|p| crate::alidisplay::decode_postprob(p[a]))
            .unwrap_or(0.0);
        tr.push((s, k, i, pp));
        if s == ST_M || s == ST_I {
            dsq.push(digitize_one(&map, ad.aseq[a]));
            i += 1;
        }
    }
    tr.push((ST_E, 0, 0, 0.0));
    tr.push((ST_C, 0, 0, 0.0));
    tr.push((crate::dp::ST_T, 0, 0, 0.0));

    debug_assert_eq!(tr.len(), ad.n + 6, "backconverted trace has unexpected size");
    debug_assert_eq!(k, ad.hmmto, "backconverted trace didn't end at hmmto");

    FauxSeq {
        name: format!("{}/{}-{}", ad.sqname, ad.sqfrom, ad.sqto),
        desc: Some(format!(
            "[subseq from] {}",
            if ad.sqdesc.is_empty() { &ad.sqname } else { &ad.sqdesc }
        )),
        acc: if ad.sqacc.is_empty() { None } else { Some(ad.sqacc.clone()) },
        tr,
        dsq,
    }
}

/// `matuse` / `matmap` / `inscount` / `alen`.
///
/// C, tracealign.c:`map_new_msa`:
/// ```text
///   if (optflags & p7_ALL_CONSENSUS_COLS) esl_vec_ISet(matuse+1, M, TRUE);
///   else                                  esl_vec_ISet(matuse+1, M, FALSE);
///   for (idx = 0; idx < nseq; idx++)
///     {
///       esl_vec_ISet(insnum, M+1, 0);
///       for (z = 1; z < tr[idx]->N; z++)
///         {
///           switch (tr[idx]->st[z]) {
///           case p7T_I:                                insnum[tr[idx]->k[z]]++; break;
///           case p7T_N: if (tr[idx]->st[z-1] == p7T_N) insnum[0]++;             break;
///           case p7T_C: if (tr[idx]->st[z-1] == p7T_C) insnum[M]++;             break;
///           case p7T_M: matuse[tr[idx]->k[z]] = TRUE;                           break;
///           ...
///         }
///       for (k = 0; k <= M; k++)
///         inscount[k] = ESL_MAX(inscount[k], insnum[k]);
///     }
///   ...
///   alen      = inscount[0];
///   for (k = 1; k <= M; k++) {
///     if (matuse[k]) { matmap[k] = alen+1; alen += 1+inscount[k]; }
///     else           { matmap[k] = alen;   alen +=   inscount[k]; }
///   }
/// ```
///
/// The N/C cases count *emitting* N and C states, which is what the `st[z-1] == st[z]`
/// test detects (the first N and the first C do not emit). Back-converted traces have
/// exactly one N and one C, neither emitting, so `inscount[0]` and `inscount[M]` stay 0
/// here -- but the general form is transcribed because `hmmalign` will feed real traces
/// through this.
fn map_new_msa(traces: &[Vec<TracePos>], m: usize) -> (Vec<usize>, Vec<bool>, Vec<usize>, usize) {
    let mut inscount = vec![0usize; m + 1];
    // p7_ALL_CONSENSUS_COLS: every consensus column appears.
    let mut matuse = vec![true; m + 1];
    matuse[0] = false;
    let mut matmap = vec![0usize; m + 1];

    let mut insnum = vec![0usize; m + 1];
    for tr in traces {
        insnum.iter_mut().for_each(|v| *v = 0);
        for z in 1..tr.len() {
            let (st, k, _, _) = tr[z];
            match st {
                ST_I => insnum[k as usize] += 1,
                ST_N if tr[z - 1].0 == ST_N => insnum[0] += 1,
                ST_C if tr[z - 1].0 == ST_C => insnum[m] += 1,
                ST_M => matuse[k as usize] = true,
                _ => {}
            }
        }
        for k in 0..=m {
            inscount[k] = inscount[k].max(insnum[k]);
        }
    }

    let mut alen = inscount[0];
    for k in 1..=m {
        if matuse[k] {
            matmap[k] = alen + 1;
            alen += 1 + inscount[k];
        } else {
            matmap[k] = alen;
            alen += inscount[k];
        }
    }
    (inscount, matuse, matmap, alen)
}

/// C, tracealign.c:`p7_tracealign_Seqs` for the text, all-consensus-columns case:
/// `map_new_msa` -> `make_text_msa` -> `annotate_rf` -> `annotate_mm` ->
/// `annotate_posterior_probability` -> `rejustify_insertions_text`.
pub fn seqs(seqs: &[FauxSeq], m: usize, mm: Option<&[u8]>) -> Msa {
    let traces: Vec<Vec<TracePos>> = seqs.iter().map(|s| s.tr.clone()).collect();
    let (inscount, matuse, matmap, alen) = map_new_msa(&traces, m);

    let sym = |d: u8| crate::alphabet::SYM[d as usize];

    // make_text_msa: every row starts as all '.', with '-' written into each consensus
    // column, and is then overwritten from the trace.
    let mut msa = Msa {
        alen,
        ..Default::default()
    };
    for s in seqs {
        let mut row = vec![b'.'; alen];
        for k in 1..=m {
            if matuse[k] {
                row[matmap[k] - 1] = b'-';
            }
        }
        let mut apos = 0usize;
        for &(st, k, i, _) in &s.tr {
            match st {
                ST_M => {
                    row[matmap[k as usize] - 1] = sym(s.dsq[i as usize]).to_ascii_uppercase();
                    // "i.e. one past the match column. remember, text mode is 0..alen-1"
                    apos = matmap[k as usize];
                }
                ST_D => {
                    // "bug #h77: if all column is deletes, do nothing; do NOT overwrite
                    // a column"
                    if matuse[k as usize] {
                        row[matmap[k as usize] - 1] = b'-';
                    }
                    apos = matmap[k as usize];
                }
                ST_I => {
                    row[apos] = sym(s.dsq[i as usize]).to_ascii_lowercase();
                    apos += 1;
                }
                ST_N | ST_C => {
                    if i > 0 {
                        row[apos] = sym(s.dsq[i as usize]).to_ascii_lowercase();
                        apos += 1;
                    }
                }
                ST_E => apos = matmap[m], // "set position for C-terminal tail"
                _ => {}
            }
        }
        msa.sqname.push(s.name.clone());
        msa.sqacc.push(s.acc.clone());
        msa.sqdesc.push(s.desc.clone());
        msa.aseq.push(row);
    }

    // annotate_rf: '.' everywhere, 'x' on the consensus columns.
    let mut rf = vec![b'.'; alen];
    for k in 1..=m {
        if matuse[k] {
            rf[matmap[k] - 1] = b'x';
        }
    }
    msa.rf = Some(rf);

    // annotate_mm: only when the HMM carries an MM line.
    //
    // C, tracealign.c:`annotate_mm`:
    //   for (k = 0; k < hmm->M; k++)
    //     if (matuse[k])
    //       msa->mm[matmap[k]-1] = hmm->mm[k];
    // Note the loop bounds and indexing are C's, off-by-one relative to every other
    // `matuse`/`matmap` walk here (`0 <= k < M` rather than `1 <= k <= M`), so node k's
    // mask lands in node k's column only if you read `hmm->mm` the same shifted way.
    // Transcribed as written rather than "corrected".
    if let Some(mmline) = mm {
        let mut v = vec![b'.'; alen];
        for k in 0..m {
            if matuse[k] {
                v[matmap[k] - 1] = mmline[k];
            }
        }
        msa.mm = Some(v);
    }

    // annotate_posterior_probability.
    let mut totp = vec![0.0f64; alen];
    let mut nppc = vec![0usize; alen];
    for s in seqs {
        let mut pp = vec![b'.'; alen];
        let mut apos = 0usize;
        for &(st, k, _, ppv) in &s.tr {
            match st {
                // Note C's fall-through: p7T_M sets the character *and* falls into
                // p7T_D's `apos = matmap[k]`.
                ST_M => {
                    pp[matmap[k as usize] - 1] = crate::alidisplay::encode_postprob(ppv);
                    totp[matmap[k as usize] - 1] += ppv as f64;
                    nppc[matmap[k as usize] - 1] += 1;
                    apos = matmap[k as usize];
                }
                ST_D => apos = matmap[k as usize],
                ST_I => {
                    pp[apos] = crate::alidisplay::encode_postprob(ppv);
                    apos += 1;
                }
                ST_N | ST_C => {}
                ST_E => apos = matmap[m],
                _ => {}
            }
        }
        msa.pp.push(Some(pp));
    }
    let mut ppc = vec![b'.'; alen];
    for apos in 0..alen {
        if nppc[apos] > 0 {
            ppc[apos] = crate::alidisplay::encode_postprob((totp[apos] / nppc[apos] as f64) as f32);
        }
    }
    msa.pp_cons = Some(ppc);

    rejustify_insertions_text(&mut msa, &inscount, &matmap, &matuse, m);
    msa
}

/// C, tracealign.c:`rejustify_insertions_text`. Insert runs longer than one column are
/// split: half the residues left-justified, half right, except at the N-terminus which is
/// right-justified entirely.
///
/// ```text
///   for (k = 0; k < M; k++)
///     if (inserts[k] > 1)
///       {
///         for (nins = 0, apos = matmap[k]; apos < matmap[k+1]-matuse[k+1]; apos++)
///           if (esl_abc_CIsResidue(abc, msa->aseq[idx][apos])) nins++;
///
///         if (k == 0) nins = 0;    /* N-terminus is right justified */
///         else        nins /= 2;   /* split in half; nins now = # of residues left left-justified  */
///
///         opos = npos = -1+matmap[k+1]-matuse[k+1];
///         while (opos >= matmap[k]+nins) {
///           if (esl_abc_CIsGap(abc, msa->aseq[idx][opos])) opos--;
///           else {
///             msa->aseq[idx][npos] = msa->aseq[idx][opos];
///             if (msa->pp != NULL && msa->pp[idx] != NULL) msa->pp[idx][npos] = msa->pp[idx][opos];
///             npos--;
///             opos--;
///           }
///         }
///         while (npos >= matmap[k]+nins) {
///           msa->aseq[idx][npos] = '.';
///           if (msa->pp != NULL && msa->pp[idx] != NULL) msa->pp[idx][npos] = '.';
///           npos--;
///         }
///       }
/// ```
///
/// `matuse` is used as an integer here (`matmap[k+1]-matuse[k+1]`), which is why it stays
/// a truth value rather than becoming a bitset.
fn rejustify_insertions_text(
    msa: &mut Msa,
    inserts: &[usize],
    matmap: &[usize],
    matuse: &[bool],
    m: usize,
) {
    let is_gap = |c: u8| c == b'-' || c == b'.' || c == b'_' || c == b'~';
    for idx in 0..msa.nseq() {
        for k in 0..m {
            if inserts[k] <= 1 {
                continue;
            }
            let hi = matmap[k + 1] - usize::from(matuse[k + 1]);
            let mut nins = 0usize;
            for apos in matmap[k]..hi {
                if !is_gap(msa.aseq[idx][apos]) {
                    nins += 1;
                }
            }
            if k == 0 {
                nins = 0; // N-terminus is right justified
            } else {
                nins /= 2;
            }

            let lo = matmap[k] + nins;
            // C walks with ints that can go below `lo`; the loop guard is `>= lo`, so
            // signed arithmetic is needed to express "one before lo".
            let mut opos = hi as isize - 1;
            let mut npos = hi as isize - 1;
            while opos >= lo as isize {
                if is_gap(msa.aseq[idx][opos as usize]) {
                    opos -= 1;
                } else {
                    msa.aseq[idx][npos as usize] = msa.aseq[idx][opos as usize];
                    if let Some(pp) = &mut msa.pp[idx] {
                        pp[npos as usize] = pp[opos as usize];
                    }
                    npos -= 1;
                    opos -= 1;
                }
            }
            while npos >= lo as isize {
                msa.aseq[idx][npos as usize] = b'.';
                if let Some(pp) = &mut msa.pp[idx] {
                    pp[npos as usize] = b'.';
                }
                npos -= 1;
            }
        }
    }
}
