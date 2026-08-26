//! A text-mode multiple alignment, and the Stockholm/Pfam writer.
//!
//! Only what `hmmsearch -A` needs: `esl_msafile_Write()` is called with
//! `eslMSAFILE_STOCKHOLM` or `eslMSAFILE_PFAM` (hmmsearch.c:566-567), both of which are
//! `stockholm_write()` differing only in characters per line, and
//! `p7_tophits_Alignment()` builds the alignment with `p7_ALL_CONSENSUS_COLS` and
//! *without* `p7_DIGITIZE` -- so this is the text path, and the digital one
//! (`make_digital_msa`, `rejustify_insertions_digital`) is not needed.

/// The subset of `ESL_MSA` that the text path fills in.
#[derive(Default)]
pub struct Msa {
    pub name: Option<String>,
    pub acc: Option<String>,
    pub desc: Option<String>,
    pub au: Option<String>,
    /// `sqname[idx]`, and the aligned text `aseq[idx]` (`alen` chars).
    pub sqname: Vec<String>,
    pub sqacc: Vec<Option<String>>,
    pub sqdesc: Vec<Option<String>>,
    pub aseq: Vec<Vec<u8>>,
    /// `#=GR <seq> PP`
    pub pp: Vec<Option<Vec<u8>>>,
    /// `#=GC PP_cons` / `#=GC RF` / `#=GC MM`
    pub pp_cons: Option<Vec<u8>>,
    pub rf: Option<Vec<u8>>,
    pub mm: Option<Vec<u8>>,
    pub alen: usize,
}

impl Msa {
    pub fn nseq(&self) -> usize {
        self.sqname.len()
    }
}

/// Characters per line: `eslMSAFILE_STOCKHOLM` wraps at 200,
/// `eslMSAFILE_PFAM` writes one block (esl_msafile.c:358-359):
/// ```text
///   case eslMSAFILE_PFAM:       return stockholm_write(fp, msa, msa->alen);
///   case eslMSAFILE_STOCKHOLM:  return stockholm_write(fp, msa, 200);
/// ```
/// `hmmsearch` picks between them on `textw > 0`, i.e. `--notextw` gives Pfam.
pub const STOCKHOLM_CPL: usize = 200;

/// C, easel/esl_msafile_stockholm.c:`stockholm_write`.
///
/// The left-margin arithmetic is the fiddly part, and it is quoted in full because
/// getting it wrong shifts every alignment row:
///
/// ```text
///   /* The left margin of an alignment block can be composed of:
///    *
///    * <seqname>                      max length: uniqwidth + maxname + 1
///    * #=GC <gc_tag>                  max length: 4 + 1 + maxgc + 1
///    * #=GR <seqname> <gr_tag>        max length: 4 + 1 + uniqwidth + maxname + 1 + maxgr + 1
///    *
///    * <margin> is the max of these. It is the total length of the
///    * left margin that we need to leave, inclusive of the last space.
///    *
///    * Then when we output, we do:
///    * name:  <leftmargin-uniqwidth-1>
///    * gc:    #=GC <leftmargin-6>
///    * gr:    #=GR <uniqwidth><maxname> <leftmargin-maxname-uniqwidth-7>
///    */
///   maxname = esl_str_GetMaxWidth(msa->sqname, msa->nseq);
///   maxgf   = esl_str_GetMaxWidth(msa->gf_tag, msa->ngf);
///   if (maxgf < 2) maxgf = 2;
///   maxgc   = esl_str_GetMaxWidth(msa->gc_tag, msa->ngc);
///   if (msa->rf      && maxgc < 2) maxgc = 2;
///   if (msa->mm      && maxgc < 2) maxgc = 2;
///   ...
///   if (msa->pp_cons && maxgc < 7) maxgc = 7;
///   maxgr   = esl_str_GetMaxWidth(msa->gr_tag, msa->ngr);
///   ...
///   if (msa->pp && maxgr < 2) maxgr = 2;
///
///   margin = uniqwidth + maxname + 1;
///   if (maxgc > 0 && maxgc+6 > margin)                   margin = maxgc+6;
///   if (maxgr > 0 && uniqwidth+maxname+maxgr+7 > margin) margin = uniqwidth+maxname+maxgr+7;
/// ```
///
/// `uniqwidth` is 0 unless the names collide, which `esl_msa_CheckUniqueNames` decides;
/// `-A`'s names are `<seqname>/<from>-<to>`, so two domains of the same target differ,
/// but nothing guarantees uniqueness in general -- see [`write`]'s handling.
///
/// `esl_str_GetMaxWidth` over an empty list returns 0, so with no `#=GF` tags of our own
/// `maxgf` is the floor of 2, which is what pads `ID`/`AC`/`DE`/`AU`.
pub fn write(out: &mut String, msa: &Msa, cpl: usize) {
    let nseq = msa.nseq();

    // C, esl_msa_CheckUniqueNames: if any two names match, every name gets a
    // `<seq#>|` prefix and a warning line is emitted.
    let mut sorted: Vec<&String> = msa.sqname.iter().collect();
    sorted.sort();
    let make_uniquenames = sorted.windows(2).any(|w| w[0] == w[1]);
    let uniqwidth = if make_uniquenames {
        let mut w = 0;
        let mut n = nseq;
        while n != 0 {
            w += 1;
            n /= 10;
        }
        w + 1 // includes the '|'
    } else {
        0
    };

    let maxname = msa.sqname.iter().map(|s| s.len()).max().unwrap_or(0);
    let maxgf = 2;
    let mut maxgc = 0usize;
    if msa.rf.is_some() {
        maxgc = maxgc.max(2);
    }
    if msa.mm.is_some() {
        maxgc = maxgc.max(2);
    }
    if msa.pp_cons.is_some() {
        maxgc = maxgc.max(7);
    }
    let mut maxgr = 0usize;
    if msa.pp.iter().any(|p| p.is_some()) {
        maxgr = maxgr.max(2);
    }

    let mut margin = uniqwidth + maxname + 1;
    if maxgc > 0 && maxgc + 6 > margin {
        margin = maxgc + 6;
    }
    if maxgr > 0 && uniqwidth + maxname + maxgr + 7 > margin {
        margin = uniqwidth + maxname + maxgr + 7;
    }

    out.push_str("# STOCKHOLM 1.0\n");
    if make_uniquenames {
        out.push_str(
            "# WARNING: seq names have been made unique by adding a prefix of \"<seq#>|\"\n",
        );
    }

    // GF section. `#=GF %-*s %s\n` with maxgf.
    let mut any_gf = false;
    for (tag, val) in [
        ("ID", &msa.name),
        ("AC", &msa.acc),
        ("DE", &msa.desc),
        ("AU", &msa.au),
    ] {
        if let Some(v) = val {
            out.push_str(&format!("#=GF {tag:<maxgf$} {v}\n"));
            any_gf = true;
        }
    }
    let _ = any_gf; // C emits the separating blank line unconditionally.
    out.push('\n');

    // GS section. Only DE and AC are ever set here; each present block is followed by a
    // blank line, and C emits that blank line whenever the *array* exists -- which for
    // `-A` it always does, since `p7_alidisplay_Backconvert` sets a description on every
    // faux sequence.
    // C keeps the uniquifying prefix and the name in *separate* conversions, so the
    // padding applies to the name alone:
    //   `#=GS %0*d|%-*s DE %s\n` with (uniqwidth-1, i, maxname, sqname[i])
    // Padding the already-prefixed string to `maxname` would be a different width.
    let prefix = |i: usize| -> String {
        if make_uniquenames {
            format!("{:0>w$}|", i, w = uniqwidth - 1)
        } else {
            String::new()
        }
    };
    if msa.sqacc.iter().any(|a| a.is_some()) {
        for i in 0..nseq {
            if let Some(a) = &msa.sqacc[i] {
                out.push_str(&format!(
                    "#=GS {}{:<maxname$} AC {}\n",
                    prefix(i),
                    msa.sqname[i],
                    a
                ));
            }
        }
        out.push('\n');
    }
    if msa.sqdesc.iter().any(|d| d.is_some()) {
        for i in 0..nseq {
            if let Some(d) = &msa.sqdesc[i] {
                out.push_str(&format!(
                    "#=GS {}{:<maxname$} DE {}\n",
                    prefix(i),
                    msa.sqname[i],
                    d
                ));
            }
        }
        out.push('\n');
    }

    // Alignment blocks.
    // `%-*s %s\n` with `margin-1`, or `%0*d|%-*s %s\n` with `margin-uniqwidth-1`.
    let seqw = margin - uniqwidth - 1;
    let grw = margin - maxname - uniqwidth - 7;
    let gcw = margin - 6;
    let mut currpos = 0usize;
    while currpos < msa.alen {
        let acpl = (msa.alen - currpos).min(cpl);
        if currpos > 0 {
            out.push('\n');
        }
        let slice = |v: &[u8]| -> String {
            String::from_utf8_lossy(&v[currpos..currpos + acpl]).into_owned()
        };
        for i in 0..nseq {
            out.push_str(&format!(
                "{}{:<seqw$} {}\n",
                prefix(i),
                msa.sqname[i],
                slice(&msa.aseq[i])
            ));
            if let Some(pp) = &msa.pp[i] {
                out.push_str(&format!(
                    "#=GR {}{:<maxname$} {:<grw$} {}\n",
                    prefix(i),
                    msa.sqname[i],
                    "PP",
                    slice(pp),
                ));
            }
        }
        if let Some(v) = &msa.pp_cons {
            out.push_str(&format!("#=GC {:<gcw$} {}\n", "PP_cons", slice(v)));
        }
        if let Some(v) = &msa.rf {
            out.push_str(&format!("#=GC {:<gcw$} {}\n", "RF", slice(v)));
        }
        if let Some(v) = &msa.mm {
            out.push_str(&format!("#=GC {:<gcw$} {}\n", "MM", slice(v)));
        }
        currpos += cpl;
    }
    out.push_str("//\n");
}
