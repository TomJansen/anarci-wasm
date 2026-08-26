//! Human-readable alignment display, i.e. what `hmmsearch` prints under
//! "Alignments for each domain".
//!
//! Faithful transcription of HMMER 3.4 `src/p7_alidisplay.c`:
//!
//! | Rust                        | C (hmmer-3.4)                                      |
//! |-----------------------------|----------------------------------------------------|
//! | `AliDisplay::create`        | `p7_alidisplay_Create()`   p7_alidisplay.c:60-260   |
//! | `AliDisplay::print`         | `p7_nontranslated_alidisplay_Print()` :1162-1229    |
//! | `encode_postprob`           | `p7_alidisplay_EncodePostProb()`   :1046-1050      |
//! | `integer_textwidth`         | `integer_textwidth()`              :1140-1145      |
//!
//! Only the protein (`long_targets == FALSE`, non-translated) path is transcribed;
//! the nhmmer/DNA reverse-complement branches are marked where elided.

use crate::alphabet::{AminoMap, SYM};
use crate::forward::ForwardFilter;
use crate::hmmfile::P7Hmm;
use crate::optacc::TracePos;
use crate::seqio::Seq;

use crate::dp::{ST_B, ST_D, ST_E, ST_I, ST_M};

/// `p7_alidisplay_EncodePostProb()` — p7_alidisplay.c:1046-1050.
#[inline]
// C, p7_alidisplay.c:1046-1050:
//   char
//   p7_alidisplay_EncodePostProb(float p)
//   {
//     return (p + 0.05 >= 1.0) ? '*' :  (char) ((p + 0.05) * 10.0) + '0';
//   }
pub fn encode_postprob(p: f32) -> u8 {
    // `p` is a `float` but `0.05` and `10.0` are **double** literals, so both the
    // addition and the multiply happen in double and only the truncation is integral.
    //
    // This is load-bearing at the boundary. For `p = 0.949999988079071f32` the exact
    // sum is 0.999999988079071: in double that is `< 1.0`, so the result is
    // `(char)9.99999988 = '9'`. Adding `0.05f32` instead rounds the sum *up* to exactly
    // `1.0f32` -- the nearest f32 to 0.99999998808 is 1.0, not 0.99999994 -- and the
    // result becomes `'*'`. `#=GC PP_cons` in `-A` output hits this.
    let q = p as f64 + 0.05;
    if q >= 1.0 {
        b'*'
    } else {
        // (char) ((p + 0.05) * 10.0) + '0' — C truncates toward zero.
        ((q * 10.0) as u8) + b'0'
    }
}

/// `p7_alidisplay_DecodePostProb()` — p7_alidisplay.c:1102-1109.
///
/// C:
/// ```text
///   float
///   p7_alidisplay_DecodePostProb(char pc)
///   {
///     if      (pc == '0') return 0.01;
///     else if (pc == '*') return 1.0;
///     else if (pc == '.') return 0.0;
///     else                return ((float) (pc - '0') / 10.);
///   }
/// ```
///
/// Deliberately crude: `pc` has already been discretized. `'0'` decodes to 0.01 rather
/// than 0 "just to avoid any possible absorbing-zero artifacts", and the round trip
/// `EncodePostProb(DecodePostProb(pc)) == pc` is required to hold.
#[inline]
pub fn decode_postprob(pc: u8) -> f32 {
    match pc {
        b'0' => 0.01,
        b'*' => 1.0,
        b'.' => 0.0,
        _ => (pc - b'0') as f32 / 10.0,
    }
}

/// `integer_textwidth()` — p7_alidisplay.c:1140-1145.
// C, p7_alidisplay.c:1140-1145:
//   static int
//   integer_textwidth(long n)
//   {
//     int w = (n < 0)? 1 : 0;
//     while (n != 0) { n /= 10; w++; }
//     return w;
//   }
fn integer_textwidth(mut n: i64) -> usize {
    let mut w = if n < 0 { 1 } else { 0 };
    while n != 0 {
        n /= 10;
        w += 1;
    }
    w
}

/// `P7_ALIDISPLAY`, restricted to the fields the protein path fills.
#[derive(Clone, Debug)]
pub struct AliDisplay {
    pub rfline: Option<Vec<u8>>,
    pub mmline: Option<Vec<u8>>,
    pub csline: Option<Vec<u8>>,
    pub model: Vec<u8>,
    pub mline: Vec<u8>,
    pub aseq: Vec<u8>,
    pub ppline: Option<Vec<u8>>,
    /// `ad->N` — number of alignment columns.
    pub n: usize,

    pub hmmname: String,
    pub hmmacc: String,
    pub hmmdesc: String,
    pub hmmfrom: i32,
    pub hmmto: i32,
    /// `ad->M` — the model's length.
    pub m: i32,

    pub sqname: String,
    pub sqacc: String,
    pub sqdesc: String,
    pub sqfrom: i64,
    pub sqto: i64,
    /// `ad->L` — the target sequence's length.
    pub l: i64,
}

impl AliDisplay {
    /// `p7_alidisplay_Create(tr, which=0, om, sq, ntsq=NULL)` — p7_alidisplay.c:60-260.
    ///
    /// `tr` is the optimal-accuracy traceback with sequence coordinates already
    /// shifted to the original `dsq` (as C does at p7_domaindef.c:857-858). Returns
    /// `None` where C returns `NULL`: a trace with no M state in the domain.
    pub fn create(tr: &[TracePos], hmm: &P7Hmm, ff: &ForwardFilter, sq: &Seq) -> Option<Self> {
        // C, p7_alidisplay.c:113-119:
        //   } else {  /* without an index, we can still do it fine:    */
        //     for (z1 = 0; which >= 0 && z1 < tr->N; z1++) if (tr->st[z1] == p7T_B) which--; /* find the right B state */
        //     if (z1 == tr->N) return NULL;
        //     for (; z1 < tr->N; z1++) if (tr->st[z1] == p7T_M) break;   /* find next M state */
        //     if (z1 == tr->N) return NULL;
        //     for (z2 = z1; z2 < tr->N; z2++) if (tr->st[z2] == p7T_E) break; /* find the next E state */
        //     for (; z2 >= 0;    z2--) if (tr->st[z2] == p7T_M) break;   /* find prev M state */
        //     if (z2 == -1) return NULL;
        //   }
        let mut z1 = 0usize;
        while z1 < tr.len() && tr[z1].0 != ST_B {
            z1 += 1;
        }
        while z1 < tr.len() && tr[z1].0 != ST_M {
            z1 += 1;
        }
        if z1 == tr.len() {
            return None;
        }
        let mut z2 = z1;
        while z2 < tr.len() && tr[z2].0 != ST_E {
            z2 += 1;
        }
        loop {
            if tr[z2].0 == ST_M {
                break;
            }
            if z2 == 0 {
                return None;
            }
            z2 -= 1;
        }

        let n = z2 - z1 + 1;
        let map = AminoMap::new();

        // The annotation lines exist only if the model carries them (om->rf[0] != 0
        // etc., p7_alidisplay.c:153-156). C 3.4 hardwires mmline to NULL.
        let mut rfline = hmm.rf.as_ref().map(|_| vec![0u8; n]);
        let csline = hmm.cs.as_ref().map(|_| vec![0u8; n]);
        let mmline: Option<Vec<u8>> = None;
        let mut csline = csline;
        let mut ppline = vec![0u8; n];
        let mut model = vec![0u8; n];
        let mut mline = vec![0u8; n];
        let mut aseq = vec![0u8; n];

        for (o, z) in (z1..=z2).enumerate() {
            let (s, k, i, pp) = tr[z];
            let k = k as usize;

            if let Some(rf) = rfline.as_mut() {
                rf[o] = if s == ST_I { b'.' } else { hmm.rf.as_ref().unwrap()[k] };
            }
            if let Some(cs) = csline.as_mut() {
                cs[o] = if s == ST_I { b'.' } else { hmm.cs.as_ref().unwrap()[k] };
            }
            ppline[o] = if s == ST_D { b'.' } else { encode_postprob(pp) };

            // C, p7_alidisplay.c:211-241:
            //   for (z = z1; z <= z2; z++)
            //     {
            //       k = tr->k[z];  i = tr->i[z];  x = sq->dsq[i];  s = tr->st[z];
            //       switch (s) {
            //       case p7T_M:
            //         ad->model[z-z1] = om->consensus[k];
            //         if      (x == esl_abc_DigitizeSymbol(om->abc, om->consensus[k])) ad->mline[z-z1] = ad->model[z-z1];
            //         else if (p7_oprofile_FGetEmission(om, k, x) > 1.0)               ad->mline[z-z1] = '+'; /* >1 not >0; om has odds ratios, not scores */
            //         else                                                             ad->mline[z-z1] = ' ';
            //         ad->aseq  [z-z1] = toupper(Alphabet[x]);
            //         break;
            //       case p7T_I:
            //         ad->model [z-z1] = '.';
            //         ad->mline [z-z1] = ' ';
            //         ad->aseq  [z-z1] = tolower(Alphabet[x]);
            //         break;
            //       case p7T_D:
            //         ad->model [z-z1] = om->consensus[k];
            //         ad->mline [z-z1] = ' ';
            //         ad->aseq  [z-z1] = '-';
            //         break;
            //       }
            //     }
            match s {
                ST_M => {
                    let x = sq.dsq[i as usize] as usize;
                    model[o] = hmm.consensus[k];
                    mline[o] = if x == map.code(hmm.consensus[k]) as usize {
                        model[o]
                    } else if ff.rfv[k][x] > 1.0 {
                        // "> 1 not > 0; om has odds ratios, not scores"
                        b'+'
                    } else {
                        b' '
                    };
                    aseq[o] = SYM[x].to_ascii_uppercase();
                }
                ST_I => {
                    let x = sq.dsq[i as usize] as usize;
                    model[o] = b'.';
                    mline[o] = b' ';
                    aseq[o] = SYM[x].to_ascii_lowercase();
                }
                ST_D => {
                    model[o] = hmm.consensus[k];
                    mline[o] = b' ';
                    aseq[o] = b'-';
                }
                _ => return None, // "invalid state in trace: not M,D,I"
            }
        }

        Some(AliDisplay {
            rfline,
            mmline,
            csline,
            model,
            mline,
            aseq,
            ppline: Some(ppline),
            n,
            hmmname: hmm.name.clone(),
            hmmacc: hmm.acc.clone().unwrap_or_default(),
            hmmdesc: hmm.desc.clone().unwrap_or_default(),
            hmmfrom: tr[z1].1,
            hmmto: tr[z2].1,
            m: hmm.m as i32,
            sqname: sq.name.clone(),
            // C, p7_alidisplay.c:174: `strcpy(ad->sqacc, sq->acc)`. Empty for FASTA
            // targets, set for EMBL/UniProt and GenBank/DDBJ ones -- which is what
            // `--acc` shows instead of the name, and what `-A` writes as `#=GS ... AC`.
            sqacc: sq.acc.clone(),
            sqdesc: sq.desc.clone(),
            sqfrom: tr[z1].2 as i64,
            sqto: tr[z2].2 as i64,
            l: sq.len() as i64,
        })
    }

    /// `p7_nontranslated_alidisplay_Print()` — p7_alidisplay.c:1162-1229.
    ///
    /// `linewidth` is C's `textw` (0 or negative means "unlimited": lay the whole
    /// alignment out on one line).
    pub fn print(&self, out: &mut String, min_aliwidth: usize, linewidth: i32, show_accessions: bool) {
        // "implement the --acc option for preferring accessions over names"
        let show_hmmname = if show_accessions && !self.hmmacc.is_empty() {
            &self.hmmacc
        } else {
            &self.hmmname
        };
        let show_seqname = if show_accessions && !self.sqacc.is_empty() {
            &self.sqacc
        } else {
            &self.sqname
        };

        let namewidth = show_hmmname.len().max(show_seqname.len());
        let coordwidth = integer_textwidth(self.hmmfrom as i64)
            .max(integer_textwidth(self.hmmto as i64))
            .max(integer_textwidth(self.sqfrom).max(integer_textwidth(self.sqto)));

        let mut aliwidth = if linewidth > 0 {
            (linewidth as isize) - namewidth as isize - 2 * coordwidth as isize - 5
        } else {
            self.n as isize
        };
        if aliwidth < self.n as isize && aliwidth < min_aliwidth as isize {
            aliwidth = min_aliwidth as isize;
        }
        let aliwidth = aliwidth.max(1) as usize;

        // C, p7_alidisplay.c:1192-1225:
        //   i1 = ad->sqfrom;
        //   k1 = ad->hmmfrom;
        //   for (pos = 0; pos < ad->N; pos += aliwidth)
        //     {
        //       if (pos > 0) fprintf(fp, "\n"); /* blank line betweeen blocks */
        //       ni = nk = 0;
        //       for (z = pos; z < pos + aliwidth && z < ad->N; z++) {
        //         if (ad->model[z] != '.') nk++; /* k advances except on insert states */
        //         if (ad->aseq[z]  != '-') ni++; /* i advances except on delete states */
        //       }
        //       k2 = k1+nk-1;
        //       if (ad->sqfrom < ad->sqto) i2 = i1+ni-1;
        //       else                       i2 = i1-ni+1; // revcomp hit for DNA
        //
        //       if (ad->csline != NULL) { strncpy(buf, ad->csline+pos, aliwidth); fprintf(fp, "  %*s %s CS\n", namewidth+coordwidth+1, "", buf); }
        //       if (ad->rfline != NULL) { strncpy(buf, ad->rfline+pos, aliwidth); fprintf(fp, "  %*s %s RF\n", namewidth+coordwidth+1, "", buf); }
        //       if (ad->mmline != NULL) { strncpy(buf, ad->mmline+pos, aliwidth); fprintf(fp, "  %*s %s MM\n", namewidth+coordwidth+1, "", buf); }
        //
        //       strncpy(buf, ad->model+pos, aliwidth); fprintf(fp, "  %*s %*d %s %-*d\n", namewidth, show_hmmname, coordwidth, k1, buf, coordwidth, k2);
        //       strncpy(buf, ad->mline+pos, aliwidth); fprintf(fp, "  %*s %s\n", namewidth+coordwidth+1, " ", buf);
        //
        //       if (ni > 0) { strncpy(buf, ad->aseq+pos, aliwidth); fprintf(fp, "  %*s %*ld %s %-*ld\n", namewidth, show_seqname, coordwidth, i1,  buf, coordwidth, i2);  }
        //       else        { strncpy(buf, ad->aseq+pos, aliwidth); fprintf(fp, "  %*s %*s %s %*s\n",    namewidth, show_seqname, coordwidth, "-", buf, coordwidth, "-"); }
        //
        //       if (ad->ppline != NULL)  { strncpy(buf, ad->ppline+pos, aliwidth);  fprintf(fp, "  %*s %s PP\n", namewidth+coordwidth+1, "", buf);  }
        //
        //       k1 += nk;
        //       if   (ad->sqfrom < ad->sqto)  i1 += ni;
        //       else                          i1 -= ni;  // revcomp hit for DNA
        //     }
        let mut i1 = self.sqfrom;
        let mut k1 = self.hmmfrom;
        let mut pos = 0usize;
        // C writes `strncpy(buf, line+pos, aliwidth)` into a buffer whose terminator
        // sits at buf[aliwidth]: a short final block keeps whatever the previous
        // block left after its own data. Slicing to the remaining length reproduces
        // that only because every line is the same length and the final block is the
        // tail of each — so each line's own tail is what would remain.
        while pos < self.n {
            if pos > 0 {
                out.push('\n');
            }
            let end = (pos + aliwidth).min(self.n);

            let mut ni = 0i64;
            let mut nk = 0i32;
            for z in pos..end {
                if self.model[z] != b'.' {
                    nk += 1;
                }
                if self.aseq[z] != b'-' {
                    ni += 1;
                }
            }
            let k2 = k1 + nk - 1;
            // The `sqfrom > sqto` (revcomp) branch is nhmmer-only.
            let i2 = i1 + ni - 1;

            let lead = namewidth + coordwidth + 1;
            let sl = |v: &[u8]| String::from_utf8_lossy(&v[pos..end]).into_owned();

            if let Some(cs) = &self.csline {
                out.push_str(&format!("  {:>lead$} {} CS\n", "", sl(cs)));
            }
            if let Some(rf) = &self.rfline {
                out.push_str(&format!("  {:>lead$} {} RF\n", "", sl(rf)));
            }
            if let Some(mm) = &self.mmline {
                out.push_str(&format!("  {:>lead$} {} MM\n", "", sl(mm)));
            }

            out.push_str(&format!(
                "  {:>namewidth$} {:>coordwidth$} {} {:<coordwidth$}\n",
                show_hmmname,
                k1,
                sl(&self.model),
                k2
            ));
            out.push_str(&format!("  {:>lead$} {}\n", " ", sl(&self.mline)));

            if ni > 0 {
                out.push_str(&format!(
                    "  {:>namewidth$} {:>coordwidth$} {} {:<coordwidth$}\n",
                    show_seqname,
                    i1,
                    sl(&self.aseq),
                    i2
                ));
            } else {
                out.push_str(&format!(
                    "  {:>namewidth$} {:>coordwidth$} {} {:>coordwidth$}\n",
                    show_seqname,
                    "-",
                    sl(&self.aseq),
                    "-"
                ));
            }

            if let Some(pp) = &self.ppline {
                out.push_str(&format!("  {:>lead$} {} PP\n", "", sl(pp)));
            }

            k1 += nk;
            i1 += ni;
            pos += aliwidth;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `p7_alidisplay_EncodePostProb()`, against values printed by C itself.
    ///
    /// The `0.95` case is the one that matters and the one an earlier version of this
    /// test got wrong: `0.05` and `10.0` are *double* literals in C, so `0.95f`
    /// (= 0.949999988079071) sums to `0.999999988 < 1.0` and encodes as `'9'`. Doing the
    /// addition in f32 rounds the sum up to exactly `1.0f` and yields `'*'`. The
    /// boundary sits between `0.949998975` and `0.950000107`, both checked here.
    #[test]
    fn postprob_encoding_matches_c() {
        assert_eq!(encode_postprob(0.0), b'0');
        assert_eq!(encode_postprob(0.04), b'0');
        assert_eq!(encode_postprob(0.05), b'1');
        assert_eq!(encode_postprob(0.5), b'5');
        assert_eq!(encode_postprob(0.94), b'9');
        assert_eq!(encode_postprob(0.95), b'9');
        assert_eq!(encode_postprob(1.0), b'*');
        assert_eq!(encode_postprob(0.949998975), b'9');
        assert_eq!(encode_postprob(0.950000107), b'*');
        // The two averages that `#=GC PP_cons` disagreed on before the fix.
        assert_eq!(encode_postprob(0.949999988079071), b'9');
        assert_eq!(encode_postprob(0.6499999761581421), b'6');
    }

    /// `p7_alidisplay_DecodePostProb()` must round-trip through the encoder, which is the
    /// property C documents.
    #[test]
    fn postprob_round_trips() {
        for pc in b"0123456789*" {
            assert_eq!(
                encode_postprob(decode_postprob(*pc)),
                *pc,
                "round trip failed for {:?}",
                *pc as char
            );
        }
        assert_eq!(decode_postprob(b'.'), 0.0);
    }

    #[test]
    fn integer_widths() {
        assert_eq!(integer_textwidth(0), 0);
        assert_eq!(integer_textwidth(1), 1);
        assert_eq!(integer_textwidth(9), 1);
        assert_eq!(integer_textwidth(10), 2);
        assert_eq!(integer_textwidth(2554), 4);
        assert_eq!(integer_textwidth(-5), 2);
    }
}
