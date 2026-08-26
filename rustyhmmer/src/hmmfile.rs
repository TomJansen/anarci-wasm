












use crate::alphabet::K;


pub const TMM: usize = 0;
pub const TMI: usize = 1;
pub const TMD: usize = 2;
pub const TIM: usize = 3;
pub const TII: usize = 4;
pub const TDM: usize = 5;
pub const TDD: usize = 6;



#[derive(Debug, Clone, Copy, Default)]
pub struct EvParams {
    pub msv_mu: f64,
    pub msv_lambda: f64,
    pub vit_mu: f64,
    pub vit_lambda: f64,
    pub fwd_tau: f64,
    pub fwd_lambda: f64,
}



#[derive(Debug, Clone, Copy, Default)]
pub struct Cutoffs {
    pub ga: Option<(f64, f64)>,
    pub tc: Option<(f64, f64)>,
    pub nc: Option<(f64, f64)>,
}





#[derive(Debug, Clone)]
pub struct P7Hmm {
    pub name: String,
    pub acc: Option<String>,
    pub desc: Option<String>,
    pub m: usize,
    
    pub mat: Vec<[f32; K]>,
    
    pub ins: Vec<[f32; K]>,
    
    pub t: Vec<[f32; 7]>,
    
    pub compo: Option<[f32; K]>,
    pub evparam: EvParams,
    pub cutoffs: Cutoffs,
    /// ASCII format generation of the record this was read from
    /// (`p7_HMMFILE_3a`..`3f` -> 0..5). Decides which optional fields the match
    /// lines carry (p7_hmmfile.c:1520-1537).
    pub format: u8,
    /// `hmm->consensus[0..=M]`, index 0 is `' '`. Always present: for pre-3e files
    /// C synthesizes it with `p7_hmm_SetConsensus()` (p7_hmmfile.c:1557).
    pub consensus: Vec<u8>,
    /// `hmm->rf` / `hmm->mm` / `hmm->cs`, each `[0..=M]` with index 0 unused.
    /// `None` when the corresponding header flag is `no`.
    pub rf: Option<Vec<u8>>,
    pub mm: Option<Vec<u8>>,
    pub cs: Option<Vec<u8>>,
    /// `hmm->map[1..=M]`, the alignment column each match state came from.
    pub map: Option<Vec<i32>>,
}

/// `p7_HMMFILE_3a` .. `p7_HMMFILE_3f` (p7_hmmfile.h), as a comparable generation
/// number: the suffix letter's offset from 'a'.
pub const FORMAT_3E: u8 = 4;
pub const FORMAT_3F: u8 = 5;


#[derive(Debug)]
pub struct ParseError(pub String);

impl std::fmt::Display for ParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "hmm parse error: {}", self.0)
    }
}
impl std::error::Error for ParseError {}


#[inline]
fn score_to_prob(tok: &str) -> f32 {
    if tok == "*" {
        0.0
    } else {
        (-tok.parse::<f32>().unwrap_or(f32::INFINITY)).exp()
    }
}

impl P7Hmm {
    
    
    
    pub fn read_one<I>(lines: &mut std::iter::Peekable<I>) -> Result<Option<Self>, ParseError>
    where
        I: Iterator<Item = std::io::Result<String>>,
    {
        
        let magic = loop {
            match lines.next() {
                None => return Ok(None), 
                Some(Ok(l)) if l.trim().is_empty() => continue,
                Some(Ok(l)) => break l,
                Some(Err(e)) => return Err(ParseError(format!("io: {e}"))),
            }
        };
        if !magic.starts_with("HMMER3/") {
            return Err(ParseError(format!("expected HMMER3/f header, got {magic:?}")));
        }
        // p7_hmmfile.c:190-196 maps the suffix letter to p7_HMMFILE_3a..3f. Unlike C,
        // which latches the first record's format for the whole file (and so aborts on
        // a mixed-format file), each record is read in its own format.
        let r#gen = magic.as_bytes().get(7).copied().unwrap_or(b'f');
        if !(b'a'..=b'f').contains(&r#gen) {
            return Err(ParseError(format!("unknown HMMER3 format in {magic:?}")));
        }
        let format = r#gen - b'a';

        let mut name = String::new();
        let mut acc = None;
        let mut desc = None;
        let mut m: usize = 0;
        let mut alph = String::new();
        let mut ev = EvParams::default();
        let mut cut = Cutoffs::default();
        // p7H_RF / p7H_MMASK / p7H_CONS / p7H_CS / p7H_MAP header flags.
        let (mut has_rf, mut has_mm, mut has_cons, mut has_cs, mut has_map) =
            (false, false, false, false, false);

        
        loop {
            let line = next_line(lines)?
                .ok_or_else(|| ParseError("EOF in header".into()))?;
            let t = line.trim();
            if t.starts_with("HMM ") || t == "HMM" {
                break; 
            }
            let p: Vec<&str> = t.split_whitespace().collect();
            if p.is_empty() {
                continue;
            }
            match p[0] {
                "NAME" => name = p.get(1).unwrap_or(&"").to_string(),
                "ACC" => acc = p.get(1).map(|s| s.to_string()),
                "DESC" => desc = Some(t[p[0].len()..].trim().to_string()),
                "LENG" => m = p.get(1).and_then(|s| s.parse().ok()).unwrap_or(0),
                "ALPH" => alph = p.get(1).unwrap_or(&"").to_string(),
                
                // C, hmmer.h:188 and :269 -- `P7_HMM.evparam` and `P7_PROFILE.evparam`
                // are `float[p7_NEVPARAM]`, and `p7_hmmfile.c` parses the `STATS LOCAL`
                // line straight into them. So the model's mu/lambda/tau are the *f32*
                // roundings of the decimals in the file, promoted back to double when
                // `esl_gumbel_surv`/`esl_exp_surv` (which take doubles) are called.
                //
                // Keeping the full f64 value of the decimal string instead shifts every
                // E-value in the last printed digit: `lnP = -lambda*(score-tau)` is
                // linear in both, so an f32-vs-f64 difference in lambda of ~1e-8
                // relative moves a lnP of ~-250 by ~2e-6, which is enough to cross a
                // `%.2g` boundary.
                "STATS" if p.len() >= 5 && p[1] == "LOCAL" => {
                    let a = p[3].parse::<f32>().unwrap_or(0.0) as f64;
                    let b = p[4].parse::<f32>().unwrap_or(0.0) as f64;
                    match p[2] {
                        "MSV" => { ev.msv_mu = a; ev.msv_lambda = b; }
                        "VITERBI" => { ev.vit_mu = a; ev.vit_lambda = b; }
                        "FORWARD" => { ev.fwd_tau = a; ev.fwd_lambda = b; }
                        _ => {}
                    }
                }
                
                "RF" => has_rf = p.get(1).is_some_and(|v| v.eq_ignore_ascii_case("yes")),
                "MM" => has_mm = p.get(1).is_some_and(|v| v.eq_ignore_ascii_case("yes")),
                "CONS" => has_cons = p.get(1).is_some_and(|v| v.eq_ignore_ascii_case("yes")),
                "CS" => has_cs = p.get(1).is_some_and(|v| v.eq_ignore_ascii_case("yes")),
                "MAP" => has_map = p.get(1).is_some_and(|v| v.eq_ignore_ascii_case("yes")),
                "GA" if p.len() >= 3 => cut.ga = Some((pf(p[1]), pf(p[2]))),
                "TC" if p.len() >= 3 => cut.tc = Some((pf(p[1]), pf(p[2]))),
                "NC" if p.len() >= 3 => cut.nc = Some((pf(p[1]), pf(p[2]))),
                _ => {}
            }
        }

        if !alph.eq_ignore_ascii_case("amino") {
            return Err(ParseError(format!(
                "rustyhmmer is amino-only; model {name:?} has ALPH {alph:?}"
            )));
        }
        if m == 0 {
            return Err(ParseError(format!("model {name:?} has LENG 0")));
        }

        let mut hmm = P7Hmm {
            name,
            acc,
            desc,
            m,
            mat: vec![[0.0; K]; m + 1],
            ins: vec![[0.0; K]; m + 1],
            t: vec![[0.0; 7]; m + 1],
            compo: None,
            evparam: ev,
            cutoffs: cut,
            format,
            consensus: vec![b' '; m + 1],
            rf: has_rf.then(|| vec![b' '; m + 1]),
            mm: has_mm.then(|| vec![b' '; m + 1]),
            cs: has_cs.then(|| vec![b' '; m + 1]),
            map: has_map.then(|| vec![0i32; m + 1]),
        };

        
        let _ = next_line(lines)?;

        
        
        
        
        
        
        let first = next_line(lines)?
            .ok_or_else(|| ParseError("EOF before COMPO/inserts".into()))?;
        let fp: Vec<&str> = first.split_whitespace().collect();
        let node0_ins_line = if fp.first() == Some(&"COMPO") {
            let mut c = [0.0f32; K];
            for i in 0..K {
                c[i] = score_to_prob(fp.get(i + 1).copied().unwrap_or("*"));
            }
            hmm.compo = Some(c);
            next_line(lines)?.ok_or_else(|| ParseError("EOF at node-0 inserts".into()))?
        } else {
            first
        };
        parse_emission(&node0_ins_line, &mut hmm.ins[0])?; 
        let node0_t = next_line(lines)?
            .ok_or_else(|| ParseError("EOF at node-0 transitions".into()))?;
        parse_transition(&node0_t, &mut hmm.t[0])?;

        
        for k in 1..=m {
            let mline = next_line(lines)?
                .ok_or_else(|| ParseError(format!("EOF at node {k} match")))?;
            if mline.trim_start().starts_with("//") {
                return Err(ParseError(format!("model ended early at node {k}")));
            }
            
            let mp: Vec<&str> = mline.split_whitespace().collect();
            for i in 0..K {
                hmm.mat[k][i] = score_to_prob(mp.get(i + 1).copied().unwrap_or("*"));
            }
            // C, p7_hmmfile.c:1519-1537:
            //   if ((status = esl_fileparser_GetTokenOnLine(hfp->efp, &tok1, NULL)) != eslOK) ...
            //   if (hmm->flags & p7H_MAP) hmm->map[k] = atoi(tok1);
            //
            //   if (hfp->format >= p7_HMMFILE_3e) {
            //     if ((status = esl_fileparser_GetTokenOnLine(hfp->efp, &tok1, NULL)) != eslOK) ...
            //     if (hmm->flags & p7H_CONS) hmm->consensus[k] = *tok1;
            //   }
            //   if ((status = esl_fileparser_GetTokenOnLine(hfp->efp, &tok1, NULL)) != eslOK) ...
            //   if (hmm->flags & p7H_RF) hmm->rf[k]   = *tok1;
            //
            //   if (hfp->format >= p7_HMMFILE_3f) {
            //     if ((status = esl_fileparser_GetTokenOnLine(hfp->efp, &tok1, NULL)) != eslOK) ...
            //     if (hmm->flags & p7H_MMASK) hmm->mm[k] = *tok1;
            //   }
            //   if ((status = esl_fileparser_GetTokenOnLine(hfp->efp, &tok1, NULL)) != eslOK) ...
            //   if (hmm->flags & p7H_CS) hmm->cs[k]   = *tok1;
            //
            // Trailing annotation fields, in C's order:
            //   MAP, [CONS if fmt >= 3e], RF, [MM if fmt >= 3f], CS.
            // Every field is always present on the line (as "-" when unset); only
            // storing it is conditional on the header flag.
            {
                let mut f = K + 1;
                let take = |f: &mut usize| -> u8 {
                    let t = mp.get(*f).copied().unwrap_or("-");
                    *f += 1;
                    t.as_bytes().first().copied().unwrap_or(b'-')
                };
                if let Some(v) = hmm.map.as_mut() {
                    v[k] = mp.get(f).and_then(|t| t.parse().ok()).unwrap_or(0);
                }
                f += 1;
                if format >= FORMAT_3E {
                    let c = take(&mut f);
                    if has_cons {
                        hmm.consensus[k] = c;
                    }
                }
                let c = take(&mut f);
                if let Some(v) = hmm.rf.as_mut() {
                    v[k] = c;
                }
                if format >= FORMAT_3F {
                    let c = take(&mut f);
                    if let Some(v) = hmm.mm.as_mut() {
                        v[k] = c;
                    }
                }
                let c = take(&mut f);
                if let Some(v) = hmm.cs.as_mut() {
                    v[k] = c;
                }
            }
            let iline = next_line(lines)?
                .ok_or_else(|| ParseError(format!("EOF at node {k} insert")))?;
            parse_emission(&iline, &mut hmm.ins[k])?;
            let tline = next_line(lines)?
                .ok_or_else(|| ParseError(format!("EOF at node {k} transition")))?;
            parse_transition(&tline, &mut hmm.t[k])?;
        }

        // "legacy issues" (p7_hmmfile.c:1556-1557): formats before 3e carry no CONS
        // field, so the consensus is synthesized from the emission probabilities.
        if format < FORMAT_3E || !has_cons {
            hmm.set_consensus();
        }

        if let Some(Ok(l)) = lines.peek() {
            if l.trim() == "//" {
                let _ = lines.next();
            }
        }

        Ok(Some(hmm))
    }

    /// `p7_hmm_SetConsensus(hmm, NULL)` — p7_hmm.c:702-731. The most likely residue
    /// at each match state, upper-cased when its probability reaches the threshold
    /// (0.5 for amino alphabets).
    // C, p7_hmm.c:713-726:
    //   if      (hmm->abc->type == eslAMINO) mthresh = 0.5;
    //   else if (hmm->abc->type == eslDNA)   mthresh = 0.9;
    //   else if (hmm->abc->type == eslRNA)   mthresh = 0.9;
    //   else                                 mthresh = 0.5;
    //
    //   hmm->consensus[0] = ' ';
    //   for (k = 1; k <= hmm->M; k++)
    //     {
    //       x = (sq ?  sq->dsq[k] : esl_vec_FArgMax(hmm->mat[k], hmm->abc->K));
    //       hmm->consensus[k] = ((hmm->mat[k][x] >= mthresh) ? toupper(hmm->abc->sym[x]) : tolower(hmm->abc->sym[x]));
    //     }
    //   hmm->consensus[hmm->M+1] = '\0';
    pub(crate) fn set_consensus(&mut self) {
        const MTHRESH: f32 = 0.5;
        self.consensus[0] = b' ';
        for k in 1..=self.m {
            // esl_vec_FArgMax: first maximum wins ties (esl_vectorops.c:403).
            let mut x = 0usize;
            for i in 1..K {
                if self.mat[k][i] > self.mat[k][x] {
                    x = i;
                }
            }
            let sym = crate::alphabet::SYM[x];
            self.consensus[k] = if self.mat[k][x] >= MTHRESH {
                sym.to_ascii_uppercase()
            } else {
                sym.to_ascii_lowercase()
            };
        }
    }

    /// Read every model in `path`, in either the ASCII or the `hmmpress` binary format.
    ///
    /// C, p7_hmmfile.c:`open_engine` sniffs the first four bytes for one of the six
    /// binary magic numbers and installs `read_bin30hmm` as the parser; otherwise it
    /// falls through to the ASCII `read_asc30hmm`. Pfam and TIGRFAM ship `hmmpress`ed
    /// files (`.h3m`), so this is the common case in practice, not an exotic one.
    /// `path` may be `-`, meaning standard input: C opens the query file through
    /// `p7_hmmfile_Open()`, which routes `"-"` to stdin
    /// (p7_hmmfile.c:`open_engine`: `else if (strcmp(filename, "-") == 0) { hfp->f = stdin; hfp->do_stdin = TRUE; ... }`).
    /// Both the ASCII and the binary format are sniffed the same way either way.
    pub fn read_all(path: &str) -> Result<Vec<Self>, ParseError> {
        use std::io::{BufRead, Read};
        let src: Box<dyn Read> = if path == "-" {
            Box::new(std::io::stdin())
        } else {
            Box::new(std::fs::File::open(path).map_err(|e| ParseError(format!("open {path}: {e}")))?)
        };
        // Sniff the magic through the BufReader rather than slurping: an ASCII
        // `Pfam-A.hmm` is well over a gigabyte, and it streams line by line.
        let mut r = std::io::BufReader::with_capacity(1 << 16, src);
        let is_binary = {
            let head = r.fill_buf().map_err(|e| ParseError(format!("read {path}: {e}")))?;
            head.len() >= 4 && magic_format(&head[..4]).is_some()
        };
        if is_binary {
            // Pressed files are an order of magnitude smaller than the ASCII form and
            // the records are positional, so this one is read whole.
            let mut bytes = Vec::new();
            r.read_to_end(&mut bytes)
                .map_err(|e| ParseError(format!("read {path}: {e}")))?;
            return Self::read_all_binary(&bytes);
        }
        let mut lines = r.lines().peekable();
        let mut out = Vec::new();
        while let Some(hmm) = Self::read_one(&mut lines)? {
            out.push(hmm);
        }
        Ok(out)
    }

    /// Every model in an `hmmpress` binary stream.
    ///
    /// C latches the format at open (p7_hmmfile.c:381-386) and then re-checks the magic
    /// at the head of every subsequent record against *that* format
    /// (p7_hmmfile.c:1611-1616), so a file mixing generations is rejected — the same
    /// latching that makes C abort partway through a mixed ASCII file.
    fn read_all_binary(bytes: &[u8]) -> Result<Vec<Self>, ParseError> {
        let format = magic_format(&bytes[..4]).expect("caller checked the magic");
        let mut r = Bin { b: bytes, at: 0 };
        let mut out = Vec::new();
        while r.at < bytes.len() {
            let magic = r.u32("magic")?;
            match magic_format(&magic.to_le_bytes()) {
                Some(f) if f == format => {}
                _ => return Err(ParseError("bad magic number at start of HMM".into())),
            }
            out.push(Self::read_one_binary(&mut r, format)?);
        }
        Ok(out)
    }

    /// C, p7_hmmfile.c:`read_bin30hmm`. The field order below is that function's,
    /// exactly; every optional field is gated on the same `hmm->flags` bit or format
    /// generation C gates it on, because the fields are positional and a wrong guess
    /// desynchronizes the rest of the record rather than failing cleanly.
    fn read_one_binary(r: &mut Bin, format: u8) -> Result<Self, ParseError> {
        // hmmer.h:109-127.
        const P7H_DESC: i32 = 1 << 1;
        const P7H_RF: i32 = 1 << 2;
        const P7H_CS: i32 = 1 << 3;
        const P7H_MAP: i32 = 1 << 8;
        const P7H_ACC: i32 = 1 << 9;
        const P7H_GA: i32 = 1 << 10;
        const P7H_TC: i32 = 1 << 11;
        const P7H_NC: i32 = 1 << 12;
        const P7H_CA: i32 = 1 << 13;
        const P7H_COMPO: i32 = 1 << 14;
        const P7H_CONS: i32 = 1 << 16;
        const P7H_MMASK: i32 = 1 << 17;

        let flags = r.i32("flags")?;
        let m_i = r.i32("model size M")?;
        if m_i <= 0 {
            return Err(ParseError(format!("nonsensical model size M={m_i}")));
        }
        let m = m_i as usize;
        // esl_alphabet.h:16-18 -- eslRNA 1, eslDNA 2, eslAMINO 3.
        let alphabet_type = r.i32("alphabet_type")?;
        if alphabet_type != 3 {
            return Err(ParseError(format!(
                "alphabet type {alphabet_type} unsupported (only amino acid is implemented)"
            )));
        }

        // Core model probabilities: mat[1..=M], ins[0..=M], t[0..=M].
        let mut mat = vec![[0.0f32; K]; m + 1];
        for k in 1..=m {
            r.f32s(&mut mat[k], "mat")?;
        }
        let mut ins = vec![[0.0f32; K]; m + 1];
        for k in 0..=m {
            r.f32s(&mut ins[k], "ins")?;
        }
        let mut t = vec![[0.0f32; 7]; m + 1];
        for k in 0..=m {
            r.f32s(&mut t[k], "t")?;
        }

        let name = r.string("name")?.unwrap_or_default();
        let acc = if flags & P7H_ACC != 0 { r.string("acc")? } else { None };
        let desc = if flags & P7H_DESC != 0 { r.string("desc")? } else { None };
        // Each annotation line is M+2 chars: index 0 unused, 1..M, and a trailing NUL.
        let rf = if flags & P7H_RF != 0 { Some(r.line(m, "rf")?) } else { None };
        let mm = if flags & P7H_MMASK != 0 { Some(r.line(m, "mm")?) } else { None };
        let cons = if flags & P7H_CONS != 0 { Some(r.line(m, "consensus")?) } else { None };
        let cs = if flags & P7H_CS != 0 { Some(r.line(m, "cs")?) } else { None };
        // `ca` (surface accessibility) is read and dropped: nothing downstream uses it,
        // but the bytes have to be consumed or every later field shifts.
        if flags & P7H_CA != 0 {
            r.line(m, "ca")?;
        }
        // Likewise comlog / nseq / eff_nseq / max_length / ctime / checksum. `hmmpress`
        // and `hmmscan` will need these preserved when they land; `hmmsearch` prints
        // none of them.
        r.string("comlog")?;
        r.i32("nseq")?;
        r.f32("eff_nseq")?;
        if format >= FORMAT_3C {
            r.i32("max_length")?;
        }
        r.string("ctime")?;
        let map = if flags & P7H_MAP != 0 {
            // hmm->map is int[M+1], index 0 unused.
            let mut v = vec![0i32; m + 1];
            for e in v.iter_mut() {
                *e = r.i32("map")?;
            }
            Some(v)
        } else {
            None
        };
        r.u32("checksum")?;

        // E-value parameters, in p7_evparams_e order (hmmer.h:72).
        let mut ev = [0.0f32; 6];
        if format >= FORMAT_3B {
            r.f32s(&mut ev, "statistical params")?;
        } else {
            // C, p7_hmmfile.c:1676-1684: "3/a files stored 3 floats: LAMBDA, MU, TAU.
            // Read 3 #'s and carefully copy/rearrange them into new 6 format".
            let mut old = [0.0f32; 3];
            r.f32s(&mut old, "statistical params")?;
            ev[5] = old[0]; // p7_FLAMBDA
            ev[4] = old[2]; // p7_FTAU
            ev[3] = old[0]; // p7_VLAMBDA
            ev[2] = old[1]; // p7_VMU
            ev[1] = ev[3]; // p7_MLAMBDA
            ev[0] = ev[2]; // p7_MMU
        }
        let evparam = EvParams {
            msv_mu: ev[0] as f64,
            msv_lambda: ev[1] as f64,
            vit_mu: ev[2] as f64,
            vit_lambda: ev[3] as f64,
            fwd_tau: ev[4] as f64,
            fwd_lambda: ev[5] as f64,
        };

        // Pfam cutoffs, in p7_cutoffs_e order (hmmer.h:73). Always six floats on disk;
        // the flags say which pairs are meaningful.
        let mut cut = [0.0f32; 6];
        r.f32s(&mut cut, "Pfam score cutoffs")?;
        let pair = |lo: usize, on: bool| -> Option<(f64, f64)> {
            if on {
                Some((cut[lo] as f64, cut[lo + 1] as f64))
            } else {
                None
            }
        };
        let cutoffs = Cutoffs {
            ga: pair(0, flags & P7H_GA != 0),
            tc: pair(2, flags & P7H_TC != 0),
            nc: pair(4, flags & P7H_NC != 0),
        };

        let compo = if flags & P7H_COMPO != 0 {
            let mut c = [0.0f32; K];
            r.f32s(&mut c, "model composition")?;
            Some(c)
        } else {
            None
        };

        let mut hmm = P7Hmm {
            name,
            acc,
            desc,
            m,
            mat,
            ins,
            t,
            compo,
            evparam,
            cutoffs,
            format,
            consensus: cons.unwrap_or_else(|| vec![b' '; m + 1]),
            rf,
            mm,
            cs,
            map,
        };
        // C, p7_hmmfile.c:1694: pre-3e files carry no consensus line, so C synthesizes
        // one. (`p7H_CONS` did not exist before 3e, which is why the flag alone is a
        // sufficient test above.)
        if format < FORMAT_3E {
            hmm.set_consensus();
        }
        Ok(hmm)
    }
}

/// `p7_HMMFILE_3b` / `3c` as this crate's generation numbers (suffix letter minus 'a').
const FORMAT_3B: u8 = 1;
const FORMAT_3C: u8 = 2;

/// The binary magic number at the head of every `hmmpress` record, mapped to this
/// crate's generation number.
///
/// C, p7_hmmfile.c:47-52:
/// ```text
///   static uint32_t  v3a_magic = 0xe8ededb6; /* 3/a binary: "hmm6" + 0x80808080 */
///   static uint32_t  v3b_magic = 0xe8ededb7; /* 3/b binary: "hmm7" + 0x80808080 */
///   static uint32_t  v3c_magic = 0xe8ededb8; /* 3/c binary: "hmm8" + 0x80808080 */
///   static uint32_t  v3d_magic = 0xe8ededb9; /* 3/d binary: "hmm9" + 0x80808080 */
///   static uint32_t  v3e_magic = 0xe8ededb0; /* 3/e binary: "hmm0" + 0x80808080 */
///   static uint32_t  v3f_magic = 0xe8ededba; /* 3/f binary: "hmma" + 0x80808080 */
/// ```
///
/// Note 3/e breaks the sequence (`...b0`, not `...ba`), so this is a lookup, not
/// arithmetic.
fn magic_format(b: &[u8]) -> Option<u8> {
    match u32::from_le_bytes([b[0], b[1], b[2], b[3]]) {
        0xe8ed_edb6 => Some(0), // 3/a
        0xe8ed_edb7 => Some(1), // 3/b
        0xe8ed_edb8 => Some(2), // 3/c
        0xe8ed_edb9 => Some(3), // 3/d
        0xe8ed_edb0 => Some(4), // 3/e
        0xe8ed_edba => Some(5), // 3/f
        _ => None,
    }
}

/// A cursor over an `hmmpress` binary stream.
///
/// **Endianness:** C writes these files with `fwrite` of raw `int`/`float`/`uint32_t`,
/// so they are host-byte-order and not portable — an `.h3m` written on a big-endian
/// machine is unreadable by C on a little-endian one, and vice versa. Reading
/// little-endian therefore matches C on every platform HMMER is actually used on, and
/// on a big-endian host neither implementation would accept the other's file anyway.
struct Bin<'a> {
    b: &'a [u8],
    at: usize,
}

impl Bin<'_> {
    fn take(&mut self, n: usize, what: &str) -> Result<&[u8], ParseError> {
        let end = self.at.checked_add(n).filter(|e| *e <= self.b.len());
        match end {
            Some(e) => {
                let s = &self.b[self.at..e];
                self.at = e;
                Ok(s)
            }
            None => Err(ParseError(format!("failed to read {what}: truncated file"))),
        }
    }
    fn i32(&mut self, what: &str) -> Result<i32, ParseError> {
        let b = self.take(4, what)?;
        Ok(i32::from_le_bytes([b[0], b[1], b[2], b[3]]))
    }
    fn u32(&mut self, what: &str) -> Result<u32, ParseError> {
        Ok(self.i32(what)? as u32)
    }
    fn f32(&mut self, what: &str) -> Result<f32, ParseError> {
        Ok(f32::from_bits(self.u32(what)?))
    }
    fn f32s(&mut self, out: &mut [f32], what: &str) -> Result<(), ParseError> {
        let b = self.take(out.len() * 4, what)?;
        for (o, c) in out.iter_mut().zip(b.chunks_exact(4)) {
            *o = f32::from_bits(u32::from_le_bytes([c[0], c[1], c[2], c[3]]));
        }
        Ok(())
    }
    /// C's `read_bin_string`: a length, then that many bytes *including* the trailing
    /// NUL (`write_bin_string` stores `strlen(s)+1`). A length of 0 means the field was
    /// `NULL`.
    fn string(&mut self, what: &str) -> Result<Option<String>, ParseError> {
        let len = self.i32(what)?;
        if len <= 0 {
            return Ok(None);
        }
        let b = self.take(len as usize, what)?;
        let s = b.strip_suffix(&[0]).unwrap_or(b);
        Ok(Some(String::from_utf8_lossy(s).into_owned()))
    }
    /// A per-node annotation line: `M+2` bytes on disk (index 0, then 1..M, then a NUL),
    /// kept as `[0..=M]` the way the ASCII reader keeps it.
    fn line(&mut self, m: usize, what: &str) -> Result<Vec<u8>, ParseError> {
        let b = self.take(m + 2, what)?;
        Ok(b[..=m].to_vec())
    }
}

#[inline]
/// C's `atof()`, which `read_asc30hmm()` uses for the GA/TC/NC cutoffs
/// (p7_hmmfile.c:1431/1436 and 1837-1838).
///
/// `atof` converts the longest leading numeric prefix and ignores whatever follows,
/// and returns 0.0 if no conversion can be performed. That matters here: Pfam and
/// TIGRFAM write the cutoff line terminated with a semicolon —
///
/// ```text
/// GA    116.25 116.25;
/// ```
///
/// so the *second* token is `"116.25;"`. Rust's `str::parse::<f64>()` rejects it
/// outright, which previously silently produced a per-domain cutoff of 0.0 (every
/// domain then passes `--cut_ga`/`--cut_tc`/`--cut_nc`), while the per-sequence
/// cutoff — an unsuffixed token — parsed fine. All 117 bac120 models and all of Pfam
/// are written this way.
///
/// Only the decimal forms `atof` can see in an HMM file are handled; hex floats and
/// `inf`/`nan` do not occur in `GA`/`TC`/`NC` lines.
fn pf(tok: &str) -> f64 {
    let s = tok.trim_start();
    let b = s.as_bytes();
    let mut end = 0usize;
    if end < b.len() && (b[end] == b'+' || b[end] == b'-') {
        end += 1;
    }
    while end < b.len() && b[end].is_ascii_digit() {
        end += 1;
    }
    if end < b.len() && b[end] == b'.' {
        end += 1;
        while end < b.len() && b[end].is_ascii_digit() {
            end += 1;
        }
    }
    // An exponent counts only if at least one digit follows it; otherwise `atof`
    // stops before the 'e' (e.g. "1.5e" converts as 1.5).
    if end < b.len() && (b[end] == b'e' || b[end] == b'E') {
        let mut e = end + 1;
        if e < b.len() && (b[e] == b'+' || b[e] == b'-') {
            e += 1;
        }
        if e < b.len() && b[e].is_ascii_digit() {
            while e < b.len() && b[e].is_ascii_digit() {
                e += 1;
            }
            end = e;
        }
    }
    s[..end].parse().unwrap_or(0.0)
}

fn next_line<I>(lines: &mut std::iter::Peekable<I>) -> Result<Option<String>, ParseError>
where
    I: Iterator<Item = std::io::Result<String>>,
{
    match lines.next() {
        None => Ok(None),
        Some(Ok(l)) => Ok(Some(l)),
        Some(Err(e)) => Err(ParseError(format!("io: {e}"))),
    }
}


fn parse_emission(line: &str, out: &mut [f32; K]) -> Result<(), ParseError> {
    let p: Vec<&str> = line.split_whitespace().collect();
    if p.len() < K {
        return Err(ParseError(format!("emission line has {} fields, need {K}", p.len())));
    }
    for i in 0..K {
        out[i] = score_to_prob(p[i]);
    }
    Ok(())
}


fn parse_transition(line: &str, out: &mut [f32; 7]) -> Result<(), ParseError> {
    let p: Vec<&str> = line.split_whitespace().collect();
    if p.len() < 7 {
        return Err(ParseError(format!("transition line has {} fields, need 7", p.len())));
    }
    for i in 0..7 {
        out[i] = score_to_prob(p[i]);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn td(name: &str) -> String {
        format!("{}/testdata/{}", env!("CARGO_MANIFEST_DIR"), name)
    }

    /// The `hmmpress` binary form of a model must parse to the same thing as its ASCII
    /// form. `globins4.h3m` was produced by C's own `hmmpress` from `globins4.hmm`.
    ///
    /// This is the whole point of the binary reader: Pfam and TIGRFAM ship pressed
    /// files, and before this `hmmsearch x.hmm.h3m` failed with "stream did not contain
    /// valid UTF-8".
    #[test]
    fn binary_and_ascii_agree() {
        let asc = P7Hmm::read_all(&td("globins4.hmm")).unwrap();
        let bin = P7Hmm::read_all(&td("globins4.h3m")).unwrap();
        assert_eq!(asc.len(), bin.len());
        for (a, b) in asc.iter().zip(bin.iter()) {
            assert_eq!(a.name, b.name);
            assert_eq!(a.m, b.m);
            assert_eq!(a.acc, b.acc);
            assert_eq!(a.desc, b.desc);
            assert_eq!(a.format, b.format);
            assert_eq!(a.consensus, b.consensus);
            assert_eq!(a.rf, b.rf);
            assert_eq!(a.mm, b.mm);
            assert_eq!(a.cs, b.cs);
            assert_eq!(a.map, b.map);
            // The ASCII form stores negative natural-log probabilities to 5 decimals,
            // so it is the *lossy* one; the binary form holds the f32 the model was
            // built with. Only agreement to the ASCII file's own precision is testable.
            for k in 1..=a.m {
                for x in 0..K {
                    let (p, q) = (a.mat[k][x], b.mat[k][x]);
                    assert!(
                        (p - q).abs() <= 1e-5 * q.max(1e-6),
                        "mat[{k}][{x}]: ascii {p:e} vs binary {q:e}"
                    );
                }
                for i in 0..7 {
                    let (p, q) = (a.t[k][i], b.t[k][i]);
                    assert!(
                        (p - q).abs() <= 1e-5 * q.max(1e-6),
                        "t[{k}][{i}]: ascii {p:e} vs binary {q:e}"
                    );
                }
            }
            // E-value parameters and cutoffs are exact in both: five decimals in ASCII,
            // and the f32 they round to is the same one the binary form stores.
            assert_eq!(a.evparam.fwd_tau, b.evparam.fwd_tau);
            assert_eq!(a.evparam.fwd_lambda, b.evparam.fwd_lambda);
            assert_eq!(a.evparam.msv_mu, b.evparam.msv_mu);
            assert_eq!(a.cutoffs.ga, b.cutoffs.ga);
        }
    }

    /// A file that is neither an `HMMER3/` header nor a known binary magic must be
    /// rejected as such, not misparsed.
    #[test]
    fn rejects_unknown_magic() {
        let dir = std::env::temp_dir().join("rustyhmmer-bad-magic.hmm");
        std::fs::write(&dir, b"\x01\x02\x03\x04not an hmm").unwrap();
        let e = P7Hmm::read_all(dir.to_str().unwrap()).unwrap_err();
        assert!(format!("{e}").contains("expected HMMER3/f header"), "{e}");
        let _ = std::fs::remove_file(&dir);
    }

    #[test]
    fn parse_globins4() {
        let hmms = P7Hmm::read_all(&td("globins4.hmm")).unwrap();
        assert_eq!(hmms.len(), 1);
        let h = &hmms[0];
        assert_eq!(h.name, "globins4");
        assert_eq!(h.m, 149);
        
        assert!((h.evparam.msv_mu - -9.9014).abs() < 1e-4);
        assert!((h.evparam.msv_lambda - 0.70957).abs() < 1e-5);
        assert!((h.evparam.fwd_tau - -4.1637).abs() < 1e-4);
        assert!(h.cutoffs.ga.is_none()); 
        
        assert!(h.compo.is_some());
        for k in 1..=h.m {
            let s: f32 = h.mat[k].iter().sum();
            assert!((s - 1.0).abs() < 1e-2, "node {k} emission sum {s}");
        }
    }

    #[test]
    fn parse_fn3_cutoffs() {
        let hmms = P7Hmm::read_all(&td("fn3.hmm")).unwrap();
        let h = &hmms[0];
        
        let (ga_seq, ga_dom) = h.cutoffs.ga.expect("fn3 has GA");
        assert!((ga_seq - 8.00).abs() < 1e-6);
        assert!((ga_dom - 7.20).abs() < 1e-6);
    }

    /// `pf` must behave like C's `atof`. The semicolon-terminated form is what Pfam
    /// and TIGRFAM actually write (`GA    116.25 116.25;`), and getting it wrong
    /// silently zeroes the per-domain cutoff.
    #[test]
    fn cutoff_tokens_parse_like_atof() {
        assert_eq!(pf("116.25"), 116.25);
        assert_eq!(pf("116.25;"), 116.25);
        assert_eq!(pf("  7.20;"), 7.20);
        assert_eq!(pf("-3.5;"), -3.5);
        assert_eq!(pf("25"), 25.0);
        assert_eq!(pf("1.5e2;"), 150.0);
        assert_eq!(pf("1.5e;"), 1.5); // exponent with no digits: atof stops at 'e'
        assert_eq!(pf(";"), 0.0); // no conversion possible
        assert_eq!(pf(""), 0.0);
    }

    /// `set_consensus()` must reproduce exactly the `CONS` column that
    /// `hmmbuild` wrote, because C generated that column with the very same
    /// `p7_hmm_SetConsensus()` (p7_builder.c:976). This is what a pre-3e file
    /// relies on, where the column is absent and we synthesize it
    /// (p7_hmmfile.c:1557).
    #[test]
    fn set_consensus_reproduces_the_file_cons_line() {
        for f in ["globins4.hmm", "fn3.hmm"] {
            let hmms = P7Hmm::read_all(&td(f)).unwrap();
            assert!(!hmms.is_empty());
            for h in &hmms {
                let parsed = h.consensus.clone();
                let mut recomputed = h.clone();
                recomputed.consensus = vec![b' '; h.m + 1];
                recomputed.set_consensus();
                assert_eq!(
                    String::from_utf8_lossy(&parsed[1..]),
                    String::from_utf8_lossy(&recomputed.consensus[1..]),
                    "{f}: consensus for model {}", h.name
                );
            }
        }
    }

    /// The trailing annotation fields must land in the right slots: fn3 has
    /// `MAP yes / CONS yes / RF no / MM no / CS yes`, so the match line carries
    /// MAP, CONS, RF, MM, CS and only MAP/CONS/CS are stored.
    #[test]
    fn match_line_annotation_fields() {
        let h = &P7Hmm::read_all(&td("fn3.hmm")).unwrap()[0];
        assert_eq!(h.format, FORMAT_3F);
        assert!(h.rf.is_none(), "fn3 has RF no");
        assert!(h.mm.is_none(), "fn3 has MM no");
        assert!(h.cs.is_some(), "fn3 has CS yes");
        assert!(h.map.is_some(), "fn3 has MAP yes");
        // fn3.hmm match-line tails, as "MAP CONS RF MM CS":
        //   k=1 "1 p - - -"   k=4 "5 P - - -"   k=5 "7 e - - C"   k=6 "8 n - - E"
        let (map, cs) = (h.map.as_ref().unwrap(), h.cs.as_ref().unwrap());
        assert_eq!((map[1], h.consensus[1], cs[1]), (1, b'p', b'-'));
        assert_eq!((map[4], h.consensus[4], cs[4]), (5, b'P', b'-'));
        assert_eq!((map[5], h.consensus[5], cs[5]), (7, b'e', b'C'));
        assert_eq!((map[6], h.consensus[6], cs[6]), (8, b'n', b'E'));
    }
}
