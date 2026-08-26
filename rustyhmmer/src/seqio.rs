

use crate::alphabet::AminoMap;




#[derive(Debug)]
pub struct Seq {
    pub name: String,
    /// `sq->acc`. FASTA carries no accession, so it stays empty there; the EMBL/UniProt
    /// `AC` line and the GenBank/DDBJ `VERSION` line set it, and C then propagates it to
    /// `hit->acc` (p7_pipeline.c:843) -- which widens `--tblout`'s `accession` column and
    /// produces `#=GS <name> AC <acc>` lines in `-A`.
    pub acc: String,
    pub desc: String,
    pub dsq: Vec<u8>,
    /// `sq->roff`, the record's byte offset in the sequence file -- what an SSI index
    /// stores and what `esl_sqfile_PositionByKey` seeks to. Left at 0 by the readers
    /// that do not track it; only the FASTA family does, which is also the only family
    /// `esl_sqfile_Position` can seek in.
    pub roff: u64,
    /// `sq->eoff`, the record's last byte. Only the formats that set `roff` set this;
    /// `--mpi`'s block formation is the one thing that reads it.
    pub eoff: u64,
}

impl Seq {
    
    pub fn len(&self) -> usize {
        self.dsq.len().saturating_sub(2)
    }
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn from_amino(name: impl Into<String>, desc: impl Into<String>, aa: &str) -> Self {
        Seq {
            name: name.into(),
            acc: String::new(),
            desc: desc.into(),
            dsq: AminoMap::new().digitize(aa.as_bytes()),
            roff: 0,
    eoff: 0,
        }
    }
}

pub fn digitize_amino<'a, I>(records: I) -> Vec<Seq>
where
    I: IntoIterator<Item = (&'a str, &'a str)>,
{
    let map = AminoMap::new();
    records
        .into_iter()
        .map(|(name, aa)| Seq {
            name: name.to_string(),
            acc: String::new(),
            desc: String::new(),
            dsq: map.digitize(aa.as_bytes()),
            roff: 0,
    eoff: 0,
        })
        .collect()
}


/// Why a sequence file could not be read, in the terms hmmsearch reports it: C
/// distinguishes a failure to *open* the file (one message, from `esl_sqfile_Open`'s
/// status) from a failure to *parse* a record (a different message, carrying
/// `esl_sqfile_GetErrorBuf()`), so the reader has to say which happened.
#[derive(Debug)]
pub enum SeqError {
    /// `esl_sqfile_Open` returned `eslENOTFOUND`.
    NotFound,
    /// `esl_sqfile_Open` returned `eslEFORMAT`. Open-time format errors do fill
    /// Easel's `errbuf`, but hmmsearch.c never prints it for this case.
    OpenFormat,
    /// `esl_sqfile_Open` returned some other non-`eslOK` status. In practice this is
    /// only `eslEOF` (3), which `fileheader_hmmpgmd` can return.
    OpenStatus(i32),
    /// A read-time `eslEFORMAT`, carrying `errbuf` verbatim. Held as bytes because the
    /// non-ASCII-character message interpolates a raw byte with `%c`.
    Parse(Vec<u8>),
    /// A genuine I/O failure. Easel would report these as system errors.
    Io(std::io::Error),
}

impl std::fmt::Display for SeqError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SeqError::NotFound => write!(f, "no such file"),
            SeqError::OpenFormat => write!(f, "empty or misformatted"),
            SeqError::OpenStatus(s) => write!(f, "unexpected open status {s}"),
            SeqError::Parse(m) => write!(f, "{}", String::from_utf8_lossy(m)),
            SeqError::Io(e) => write!(f, "{e}"),
        }
    }
}

fn io_to_seqerror(e: std::io::Error) -> SeqError {
    if e.kind() == std::io::ErrorKind::NotFound {
        SeqError::NotFound
    } else {
        SeqError::Io(e)
    }
}

/// `easel.h:355-359`:
/// ```text
///   #define eslDSQ_SENTINEL 255     /* sentinel bytes 0,L+1 in a dsq */
///   #define eslDSQ_ILLEGAL  254     /* input symbol is unmapped and unexpected */
///   #define eslDSQ_IGNORED  253     /* input symbol is unmapped and ignored */
///   #define eslDSQ_EOL      252     /* input symbol marks end of a line */
///   #define eslDSQ_EOD      251     /* input symbol marks end of a seq record */
/// ```
const DSQ_SENTINEL: u8 = 255;
const DSQ_ILLEGAL: u8 = 254;
const DSQ_IGNORED: u8 = 253;
const DSQ_EOL: u8 = 252;
const DSQ_EOD: u8 = 251;

/// `esl_sqio.h:122`: `#define eslREADBUFSIZE 4096`.
const READBUFSIZE: usize = 4096;

/// The digital amino alphabet's input map, `abc->inmap`.
///
/// C, easel/esl_alphabet.c:`create_amino`:
/// ```text
///   a = esl_alphabet_CreateCustom("ACDEFGHIKLMNPQRSTVWY-BJZOUX*~", 20, 29);
///   esl_alphabet_SetEquiv(a, '_', '-');       /* allow _ as a gap too */
///   esl_alphabet_SetEquiv(a, '.', '-');       /* allow . as a gap too */
///   esl_alphabet_SetCaseInsensitive(a);       /* allow lower case input */
/// ```
/// and `esl_alphabet_CreateCustom`:
/// ```text
///   for (c = 0; c < 128; c++)   a->inmap[c]               = eslDSQ_ILLEGAL;
///   for (x = 0; x < a->Kp; x++) a->inmap[(int) a->sym[x]] = x;
/// ```
/// So `.` and `_` map to the *gap* code and stay legal in a sequence body, while the
/// sequence-file map below overrides `-` to illegal. `SetDegeneracy` touches only the
/// degeneracy tables, not `inmap`.
fn abc_amino_inmap() -> [u8; 256] {
    let mut m = [DSQ_ILLEGAL; 256];
    for (code, &c) in crate::alphabet::SYM.iter().enumerate() {
        m[c as usize] = code as u8;
    }
    m[b'_' as usize] = crate::alphabet::GAP;
    m[b'.' as usize] = crate::alphabet::GAP;
    // esl_alphabet_SetCaseInsensitive: "for every letter that is mapped in either lower
    // or upper case, map the other case to the same internal code". Only upper case is
    // in `sym`, so this is the one direction.
    for lc in b'a'..=b'z' {
        let uc = lc.to_ascii_uppercase();
        if m[uc as usize] != DSQ_ILLEGAL {
            m[lc as usize] = m[uc as usize];
        }
    }
    m
}

/// The sequence-file input map, `sqfp->inmap`, for the FASTA family in digital mode.
///
/// C, easel/esl_sqio_ascii.c:`inmap_fasta` (and `inmap_daemon`, which differs only in
/// the end-of-data character):
/// ```text
///   if (abc_inmap != NULL) {
///     for (x = 0; x < 128; x++) sqfp->inmap[x] = abc_inmap[x];
///     sqfp->inmap['-']  = eslDSQ_ILLEGAL;   // don't let the abc_inmap map the gap char; this is an ungapped file format
///   } else { ... }
///   sqfp->inmap['*']  = '*';         /* accept * as a nonresidue/stop codon character */
///   sqfp->inmap[' ']  = eslDSQ_IGNORED;
///   sqfp->inmap['\t'] = eslDSQ_IGNORED;
///   sqfp->inmap['\r'] = eslDSQ_IGNORED;/* DOS eol compatibility */
///   sqfp->inmap['\n'] = eslDSQ_EOL;
///   sqfp->inmap['>']  = eslDSQ_EOD;
/// ```
/// `hmmsearch` always reaches this through `esl_sqfile_OpenDigital`, so `abc_inmap` is
/// the amino map above.
///
/// One deliberate simplification: C's `sqfp->inmap['*'] = '*'` stores the literal
/// character 42, not the alphabet's code. `seebuf()` only asks whether the value is
/// `<= 127` -- 42 and `p7_NONRESIDUE` (27) both are -- and `addbuf()` digitizes through
/// `sq->abc->inmap` rather than this table, so leaving `*` at the alphabet's code gives
/// identical results from one table instead of two.
fn sq_amino_inmap(eod: u8) -> [u8; 256] {
    let mut m = abc_amino_inmap();
    m[b'-' as usize] = DSQ_ILLEGAL;
    m[b' ' as usize] = DSQ_IGNORED;
    m[b'\t' as usize] = DSQ_IGNORED;
    m[b'\r' as usize] = DSQ_IGNORED;
    m[b'\n' as usize] = DSQ_EOL;
    m[eod as usize] = DSQ_EOD;
    m
}

/// C's `isspace()` in the default locale.
#[inline]
fn is_space(c: u8) -> bool {
    matches!(c, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
}

fn parse_err(msg: impl Into<Vec<u8>>) -> SeqError {
    SeqError::Parse(msg.into())
}

/// The `ESL_SQASCII_DATA` block cursor: `buf`, `nc`, `bpos`, `linenumber`.
///
/// C reads a *block* of `eslREADBUFSIZE` bytes at a time for the FASTA family
/// (`config_fasta`/`config_daemon` both set `is_linebased = FALSE`), and one of the
/// parsers -- `end_daemon` -- behaves differently depending on whether the record it is
/// terminating happens to end at a block boundary. So the block size is part of the
/// observable behaviour and is reproduced here rather than abstracted away.
struct SqBlocks {
    src: Box<dyn std::io::Read>,
    buf: Vec<u8>,
    nc: usize,
    bpos: usize,
    /// `ascii->boff`, the file offset of the block now in `buf`. `header_fasta` records
    /// `sq->roff = ascii->boff + ascii->bpos` at the `>`, which is the offset an SSI
    /// index stores for that record.
    boff: u64,
    /// `ascii->linenumber`, which `esl_sqascii_Open` initializes to 1.
    linenumber: i64,
}

impl SqBlocks {
    fn new(src: Box<dyn std::io::Read>) -> Self {
        SqBlocks {
            // Two bytes of slack past eslREADBUFSIZE: `end_daemon` reads two bytes
            // without a bounds check, which in C can run one past `nc` when a record's
            // first `/` is the block's last byte. C reads whatever stale byte is there;
            // reading a zero here takes the same branch (the "did not find //
            // terminator" error) without a panic.
            src,
            buf: vec![0u8; READBUFSIZE + 2],
            nc: 0,
            bpos: 0,
            boff: 0,
            linenumber: 1,
        }
    }

    /// `loadbuf()` for `is_linebased == FALSE`, which is `loadmem()`'s single
    /// `fread(mem, 1, eslREADBUFSIZE, fp)` handed out as the whole buffer:
    /// ```text
    ///   ascii->buf    = ascii->mem  + ascii->mpos;
    ///   ascii->bpos   = 0;
    ///   ascii->nc     = ascii->mn - ascii->mpos;
    ///   ...
    ///   return (ascii->nc == 0 ? eslEOF : eslOK);
    /// ```
    /// Returns `false` for `eslEOF`.
    fn loadbuf(&mut self) -> Result<bool, SeqError> {
        // fread() returns a short count only at end of file, so fill the block.
        let mut n = 0;
        while n < READBUFSIZE {
            match std::io::Read::read(&mut self.src, &mut self.buf[n..READBUFSIZE]) {
                Ok(0) => break,
                Ok(k) => n += k,
                Err(ref e) if e.kind() == std::io::ErrorKind::Interrupted => {}
                Err(e) => return Err(SeqError::Io(e)),
            }
        }
        self.boff += self.nc as u64;
        self.nc = n;
        self.bpos = 0;
        Ok(n != 0)
    }

    /// `nextchar()`:
    /// ```text
    ///   ascii->bpos++;
    ///   if (ascii->nc == ascii->bpos && (status = loadbuf(sqfp)) != eslOK) return status;
    ///   *ret_c = ascii->buf[ascii->bpos];
    /// ```
    /// Note that `bpos` ends up *at* the returned character, not past it. `None` is
    /// `eslEOF`, and C leaves `c` at its previous value in that case, so callers here
    /// must do the same.
    fn nextchar(&mut self) -> Result<Option<u8>, SeqError> {
        self.bpos += 1;
        if self.bpos == self.nc && !self.loadbuf()? {
            return Ok(None);
        }
        Ok(Some(self.buf[self.bpos]))
    }
}

/// Which of the three FASTA-family formats we are parsing. They share `header_fasta`
/// and differ only in the end-of-data character, whether EOF may end a record, and
/// whether the file carries a header line.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum FastaKind {
    /// `config_fasta`: EOD `>`, `eof_is_ok = TRUE`, `parse_end = end_fasta`.
    Fasta,
    /// `config_daemon`: EOD `/`, `eof_is_ok = FALSE`, `parse_end = end_daemon`.
    Daemon,
    /// `eslSQFILE_HMMPGMD` uses `config_fasta` plus `fileheader_hmmpgmd`.
    Hmmpgmd,
}

/// `header_fasta()`, which parses one `>name description` line.
///
/// ```text
///   c =  ascii->buf[ascii->bpos];
///   while (status == eslOK && isspace(c)) status = nextchar(sqfp, &c); /* skip space (including \n) */
///   if (status == eslEOF) return eslEOF;
///   if (status == eslOK && c == '>') {    /* accept the > */
///     sq->roff = ascii->boff + ascii->bpos;
///     status = nextchar(sqfp, &c);
///   } else if (c != '>') ESL_FAIL(eslEFORMAT, ascii->errbuf, "Line %" PRId64 ": unexpected char %c; expected FASTA to start with >", ascii->linenumber, c);
///   while (status == eslOK && (c == '\t' || c == ' ')) status = nextchar(sqfp, &c);
///   pos = 0;
///   while (status == eslOK && ! isspace(c)) { sq->name[pos++] = c; ... status = nextchar(sqfp, &c); }
///   if (pos == 0) ESL_FAIL(eslEFORMAT, ascii->errbuf, "Line %" PRId64 ": no FASTA name found", ascii->linenumber);
///   while (status == eslOK &&  (c == '\t' || c == ' ')) status = nextchar(sqfp, &c);
///   pos = 0;
///   while (status == eslOK && c != '\n' && c != '\r' && c != 1) { sq->desc[pos++] = c; ... }
///   while (status == eslOK && c != '\n' && c != '\r') status = nextchar(sqfp, &c);
///   while (status == eslOK && (c == '\n' || c == '\r')) status = nextchar(sqfp, &c);
///   ascii->linenumber++;
/// ```
/// The description is delimited by end of line *or* ctrl-A, "Patched to deal with NCBI
/// NR desclines [SRE:H1/82]", and is stored verbatim -- trailing spaces included, which
/// is visible in `--tblout`'s description column.
///
/// `None` is `eslEOF`: no more records.
fn header_fasta(b: &mut SqBlocks) -> Result<Option<(String, String, u64)>, SeqError> {
    // "make sure there are characters in the buffer"
    if b.nc == b.bpos && !b.loadbuf()? {
        return Ok(None);
    }
    let mut c = b.buf[b.bpos];
    let mut eof = false;
    while is_space(c) {
        match b.nextchar()? {
            Some(x) => c = x,
            None => {
                eof = true;
                break;
            }
        }
    }
    if eof {
        return Ok(None);
    }
    let roff = b.boff + b.bpos as u64; // sq->roff = ascii->boff + ascii->bpos
    if c == b'>' {
        match b.nextchar()? {
            Some(x) => c = x,
            None => eof = true,
        }
    } else {
        return Err(parse_err(format!(
            "Line {}: unexpected char {}; expected FASTA to start with >",
            b.linenumber, c as char
        )));
    }
    macro_rules! advance {
        () => {
            match b.nextchar()? {
                Some(x) => c = x,
                None => {
                    eof = true;
                }
            }
        };
    }
    while !eof && (c == b'\t' || c == b' ') {
        advance!();
    }
    let mut name = Vec::new();
    while !eof && !is_space(c) {
        name.push(c);
        advance!();
    }
    if name.is_empty() {
        return Err(parse_err(format!(
            "Line {}: no FASTA name found",
            b.linenumber
        )));
    }
    while !eof && (c == b'\t' || c == b' ') {
        advance!();
    }
    let mut desc = Vec::new();
    while !eof && c != b'\n' && c != b'\r' && c != 1 {
        desc.push(c);
        advance!();
    }
    while !eof && c != b'\n' && c != b'\r' {
        advance!();
    }
    while !eof && (c == b'\n' || c == b'\r') {
        advance!();
    }
    b.linenumber += 1;
    Ok(Some((
        String::from_utf8_lossy(&name).into_owned(),
        String::from_utf8_lossy(&desc).into_owned(),
        roff,
    )))
}

/// `fileheader_hmmpgmd()`: skip the `#<nseq> <nres>` line that an hmmpgmd database
/// starts with; the remainder of the file is FASTA.
///
/// ```text
///   c =  ascii->buf[ascii->bpos];
///   while (status == eslOK && isspace(c)) status = nextchar(sqfp, &c); /* skip space (including \n, \r) */
///   if (status == eslEOF) return eslEOF;
///   if (c != '#') ESL_FAIL(eslEFORMAT, ascii->errbuf, "hmmpgmd file expected to start with #");
///   /* skip first line; remainder of file is FASTA format */
///   while (status == eslOK && (c != '\n' && c != '\r')) status = nextchar(sqfp, &c);
///   if (status == eslEOF) return eslEOF;
/// ```
/// Both `eslEOF` returns propagate out of `esl_sqascii_Open` unchanged, which is why
/// hmmsearch reports "Unexpected error 3 opening sequence file" for a database that is
/// only a header line with no terminating newline. Note that this does not increment
/// `linenumber`, so line numbers in later error messages are one less than the file's.
///
/// `Ok(false)` is `eslEOF`.
fn fileheader_hmmpgmd(b: &mut SqBlocks) -> Result<bool, SeqError> {
    let mut c = b.buf[b.bpos];
    while is_space(c) {
        match b.nextchar()? {
            Some(x) => c = x,
            None => return Ok(false),
        }
    }
    if c != b'#' {
        return Err(SeqError::OpenFormat);
    }
    while c != b'\n' && c != b'\r' {
        match b.nextchar()? {
            Some(x) => c = x,
            None => return Ok(false),
        }
    }
    Ok(true)
}

/// `seebuf()`: scan the rest of the block, counting residues and validating them
/// against the inmap, and stop early at the end-of-data character.
///
/// ```text
///   for (bpos = ascii->bpos; nres < maxn && bpos < ascii->nc; bpos++)
///   {
///       sym = ascii->buf[bpos];
///       if (!isascii(sym)) ESL_FAIL(eslEFORMAT, ascii->errbuf, "Line %" PRId64 ": non-ASCII character %c in sequence", ascii->linenumber, sym);
///       x   = sqfp->inmap[sym];
///       if      (x <= 127) nres++;
///       else if (x == eslDSQ_EOL)     { ... if (ascii->linenumber != -1) ascii->linenumber++; }
///       else if (x == eslDSQ_ILLEGAL) ESL_FAIL(eslEFORMAT, ascii->errbuf, "Line %" PRId64 ": illegal character %c", ascii->linenumber, sym);
///       else if (x == eslDSQ_EOD)     { status = eslEOD; break; }
///       else if (x != eslDSQ_IGNORED) ESL_FAIL(eslEFORMAT, ascii->errbuf, "inmap corruption?");
///   }
/// ```
/// Returns `(nres, endpos, saw_eod)`. C's `maxn` limit exists for `ReadWindow`, which
/// hmmsearch does not use, so this is always the `-1` (whole block) case. C's
/// `curbpl/currpl/prvbpl/prvrpl/bpl/rpl` bookkeeping is dropped with it: those only
/// feed SSI index *creation*, which hmmsearch never does.
fn seebuf(b: &mut SqBlocks, inmap: &[u8; 256]) -> Result<(usize, usize, bool), SeqError> {
    let mut nres = 0usize;
    let mut bpos = b.bpos;
    let mut eod = false;
    while bpos < b.nc {
        let sym = b.buf[bpos];
        if sym >= 128 {
            let mut msg = format!("Line {}: non-ASCII character ", b.linenumber).into_bytes();
            msg.push(sym);
            msg.extend_from_slice(b" in sequence");
            return Err(SeqError::Parse(msg));
        }
        let x = inmap[sym as usize];
        if x <= 127 {
            nres += 1;
        } else if x == DSQ_EOL {
            b.linenumber += 1;
        } else if x == DSQ_ILLEGAL {
            return Err(parse_err(format!(
                "Line {}: illegal character {}",
                b.linenumber, sym as char
            )));
        } else if x == DSQ_EOD {
            eod = true;
            break;
        }
        bpos += 1;
    }
    Ok((nres, bpos, eod))
}

/// `addbuf()`: re-walk the residues `seebuf` just counted and digitize them.
///
/// ```text
///   while (nres) {
///     x  = sq->abc->inmap[(int) ascii->buf[ascii->bpos++]];
///     if (x <= 127) { nres--; sq->dsq[++sq->n] = x; }
///   } /* we skipped IGNORED, EOL. EOD, ILLEGAL don't occur; seebuf() already checked  */
/// ```
/// The map here is the *alphabet's*, not the sequence file's -- the two differ on `*`
/// (see [`sq_amino_inmap`]) and on the characters the file map marks ignorable, which
/// the alphabet marks illegal and this loop therefore skips.
fn addbuf(b: &mut SqBlocks, dsq: &mut Vec<u8>, abc: &[u8; 256], mut nres: usize) {
    while nres > 0 {
        let x = abc[b.buf[b.bpos] as usize];
        b.bpos += 1;
        if x <= 127 {
            nres -= 1;
            dsq.push(x);
        }
    }
}

/// `end_fasta()`, which only asserts that the reader stopped where it should have:
/// ```text
///   if (ascii->bpos < ascii->nc) {
///     if (ascii->buf[ascii->bpos] != '>') ESL_FAIL(eslEFORMAT, ascii->errbuf, "Whoops, FASTA reader is corrupted");
///     sq->eoff = ascii->boff + ascii->bpos - 1; /* this puts eoff at the last \n */
///   } /* else, EOF, and we don't have to do anything. */
/// ```
fn end_fasta(b: &SqBlocks) -> Result<(), SeqError> {
    if b.bpos < b.nc && b.buf[b.bpos] != b'>' {
        return Err(parse_err(&b"Whoops, FASTA reader is corrupted"[..]));
    }
    Ok(())
}

/// `end_daemon()`, which consumes the `//` that terminates a daemon record.
///
/// ```text
///   if (ascii->nc < 3) ESL_FAIL(eslEFORMAT, ascii->errbuf, "Whoops, daemon input stream is corrupted");
///   c =  ascii->buf[ascii->bpos++];
///   if (c != '/') ESL_FAIL(eslEFORMAT, ascii->errbuf, "Line %" PRId64 ": did not find // terminator at end of seq record", ascii->linenumber);
///   c =  ascii->buf[ascii->bpos++];
///   if (c != '/') ESL_FAIL(eslEFORMAT, ascii->errbuf, "Line %" PRId64 ": did not find // terminator at end of seq record", ascii->linenumber);
///   /* skip to end of line */
///   while (c != '\n' && c != '\r' && ascii->bpos < ascii->nc) c =  ascii->buf[ascii->bpos++];
///   /* skip past end of line */
///   while ((c == '\n' || c == '\r') && ascii->bpos < ascii->nc) c =  ascii->buf[ascii->bpos++];
/// ```
/// The second loop overshoots by one character whenever the block does not end at the
/// `//` line: it reads the newline, sees it is a newline, and reads *again*, consuming
/// the `>` that opens the next record. So C only accepts a multi-record daemon stream
/// if there is a filler character between the newline and the next `>` -- which is
/// consistent with the format's purpose, a pipe carrying one record at a time, where
/// the block ends at the `//` and `bpos < nc` stops the loop. Transcribed as-is:
/// rejecting the same files C rejects is the point.
fn end_daemon(b: &mut SqBlocks) -> Result<(), SeqError> {
    let terminator = |b: &SqBlocks| {
        parse_err(format!(
            "Line {}: did not find // terminator at end of seq record",
            b.linenumber
        ))
    };
    if b.nc < 3 {
        return Err(parse_err(
            &b"Whoops, daemon input stream is corrupted"[..],
        ));
    }
    let mut c = b.buf[b.bpos];
    b.bpos += 1;
    if c != b'/' {
        return Err(terminator(b));
    }
    c = b.buf[b.bpos];
    b.bpos += 1;
    if c != b'/' {
        return Err(terminator(b));
    }
    while c != b'\n' && c != b'\r' && b.bpos < b.nc {
        c = b.buf[b.bpos];
        b.bpos += 1;
    }
    while (c == b'\n' || c == b'\r') && b.bpos < b.nc {
        c = b.buf[b.bpos];
        b.bpos += 1;
    }
    Ok(())
}

/// `sqascii_Read`'s main case, for the three FASTA-family formats.
///
/// ```text
///   if (ascii->nc == 0) return eslEOF;
///   if ((status = ascii->parse_header(sqfp, sq)) != eslOK) return status;
///   do {
///     if ((status = seebuf(sqfp, -1, &n, &epos)) == eslEFORMAT) return status;
///     if (esl_sq_GrowTo(sq, sq->n + n) != eslOK) return eslEMEM;
///     addbuf(sqfp, sq, n);
///     ascii->L   += n;
///     sq->eoff   = ascii->boff + epos - 1;
///     if (status == eslEOD)     break;
///   } while ((status = loadbuf(sqfp)) == eslOK);
///   if      (status == eslEOF)
///     {
///       if (! ascii->eof_is_ok) ESL_FAIL(eslEFORMAT, ascii->errbuf, "Unexpected EOF; file truncated?");
///       if ((status = ascii->parse_end(sqfp, sq)) != eslOK) return status;
///     }
///   else if (status == eslEOD)
///     {
///       ascii->bpos = epos;
///       if ((status = ascii->parse_end(sqfp, sq)) != eslOK) return status;
///     }
/// ```
/// hmmsearch reads the whole database through this loop, interleaved with searching;
/// here it is read up front, so a parse error surfaces before any output rather than
/// after a partial report. In practice C's threaded reader loads a block of sequences
/// before searching any of them, so the error precedes output there too.
fn read_fasta_family(b: &mut SqBlocks, kind: FastaKind) -> Result<Vec<Seq>, SeqError> {
    let eod = if kind == FastaKind::Daemon { b'/' } else { b'>' };
    let eof_is_ok = kind != FastaKind::Daemon;
    let sq_inmap = sq_amino_inmap(eod);
    let abc = abc_amino_inmap();
    let parse_end = |b: &mut SqBlocks| -> Result<(), SeqError> {
        if kind == FastaKind::Daemon {
            end_daemon(b)
        } else {
            end_fasta(b)
        }
    };

    let mut out = Vec::new();
    loop {
        if b.nc == 0 {
            break;
        }
        let (name, desc, roff) = match header_fasta(b)? {
            Some(h) => h,
            None => break,
        };
        let mut dsq = vec![DSQ_SENTINEL];
        // `sq->eoff = ascii->boff + epos - 1`, reassigned on every buffer.
        let mut eoff;
        loop {
            let (n, epos, saw_eod) = seebuf(b, &sq_inmap)?;
            eoff = b.boff + epos as u64;
            // C's `esl_sq_GrowTo(sq, sq->n + n)`, which is why seebuf() returns a count
            // before addbuf() writes anything.
            dsq.reserve(n);
            addbuf(b, &mut dsq, &abc, n);
            if saw_eod {
                b.bpos = epos;
                parse_end(b)?;
                break;
            }
            if !b.loadbuf()? {
                if !eof_is_ok {
                    return Err(parse_err(&b"Unexpected EOF; file truncated?"[..]));
                }
                parse_end(b)?;
                break;
            }
        }
        dsq.push(DSQ_SENTINEL);
        out.push(Seq {
            name,
            acc: String::new(),
            desc,
            dsq,
            roff,
            // `sq->eoff = ascii->boff + epos - 1` from the read loop, which end_fasta
            // then confirms; at end of file the last block's end is the file's end.
            eoff: eoff.saturating_sub(1),
        });
    }
    Ok(out)
}

/// `esl_sqascii_Open`'s work for the FASTA family: open the stream and preload the
/// first block, plus the file header for hmmpgmd.
///
/// ```text
///   if (strcmp(filename, "-") == 0) { ascii->fp = stdin; ascii->do_stdin = TRUE; }
///   else if ((ascii->fp = fopen(filename, "r")) == NULL) { status = eslENOTFOUND; goto ERROR; }
///   ...
///   /* Preload the first line or chunk of file. */
///   status = loadbuf(sqfp);
///   if      (status == eslEOF) { status = eslEFORMAT; goto ERROR; }
///   else if (status != eslOK)  { goto ERROR; }
///   switch (format) {
///   case eslSQFILE_HMMPGMD:   status = fileheader_hmmpgmd(sqfp); break;
///   default:                  status = eslOK;                    break;
///   }
///   if (status != eslOK) goto ERROR;
/// ```
/// `.gz` files, which C handles by `popen`ing `gzip -dc`, are not supported here.
fn open_fasta_family(path: &str, kind: FastaKind) -> Result<SqBlocks, SeqError> {
    let src: Box<dyn std::io::Read> = if path == "-" {
        Box::new(std::io::stdin())
    } else {
        Box::new(std::fs::File::open(path).map_err(|_| SeqError::NotFound)?)
    };
    open_fasta_stream(src, kind)
}

/// `esl_sqfile_Open`'s NCBI branch, from easel/esl_sqio.c:`sqfile_open`:
/// ```text
///   if (strcmp(filename, "-") == 0) { /* stdin special case */
///     if ((status = esl_strdup(filename, -1, &path)) != eslOK) goto ERROR;
///     if ((status = esl_sqascii_Open(path, sqfp->format, sqfp)) != eslOK) goto ERROR;
///   } else {
///     status = eslENOTFOUND;
///     if (format == eslSQFILE_NCBI && status == eslENOTFOUND)
///       status = esl_sqncbi_Open(sqfp->filename, sqfp->format, sqfp);
///     if (status == eslENOTFOUND)
///       status = esl_sqascii_Open(sqfp->filename, sqfp->format, sqfp);
///     /* if it's not there, then check in directory list provided by <env>. */
///     ...
/// ```
/// A pressed database cannot come from stdin: `esl_sqascii_Open` opens with
/// `if (format == eslSQFILE_NCBI) return eslENOTFOUND;`, so the stdin path fails, and so
/// does the ASCII fallback after `esl_sqncbi_Open` declines.
///
/// The `BLASTDB` directory search that follows is not implemented, because it does not
/// work in C either. `sqfile_open` builds its copy of the variable with
/// ```text
///   ESL_ALLOC(list, sizeof(char) * (strlen(s1) + 1));
///   strcpy(list + 2, s1);
/// ```
/// -- two bytes into a buffer sized for none, leaving `list[0]` and `list[1]`
/// uninitialized in front of the first directory and writing two bytes past the end. Every
/// candidate path it then builds is prefixed with that garbage. Verified: with `BLASTDB`
/// pointing at a directory holding a valid database, C HMMER 3.4 still reports "Failed to
/// open sequence file %s for reading", for one directory and for several.
fn open_ncbi(path: &str) -> Result<crate::ncbi::NcbiDb, SeqError> {
    if path == "-" {
        return Err(SeqError::NotFound);
    }
    crate::ncbi::open(path)
}

/// [`open_fasta_family`] once the stream is open.
fn open_fasta_stream(src: Box<dyn std::io::Read>, kind: FastaKind) -> Result<SqBlocks, SeqError> {
    let mut b = SqBlocks::new(src);
    if !b.loadbuf()? {
        return Err(SeqError::OpenFormat);
    }
    if kind == FastaKind::Hmmpgmd && !fileheader_hmmpgmd(&mut b)? {
        return Err(SeqError::OpenStatus(3)); // eslEOF, straight out of Open
    }
    Ok(b)
}

/// An opened sequence database, between `esl_sqfile_Open` and the first
/// `esl_sqio_Read`.
///
/// The split matters for more than tidiness: C opens the target database *before* the
/// query HMM file but does not parse a single record until the search loop, so a
/// database whose first record is malformed still gets past the open, and it is the HMM
/// file's error -- or the report header -- that comes out first.
pub struct SeqSource {
    inner: SourceInner,
}

enum SourceInner {
    FastaFamily { blocks: SqBlocks, kind: FastaKind },
    Ncbi(crate::ncbi::NcbiDb),
    /// The flat-file and alignment readers here are not split into open and read
    /// phases; they re-open the path and parse it in one pass. Only the `fopen` half of
    /// C's open is reproduced for them, which is what distinguishes `eslENOTFOUND`.
    Deferred { path: String, fmt: Format },
}

/// `esl_sqfile_Open`. Everything this can detect, C detects at open time too.
pub fn open_format(path: &str, fmt: Format) -> Result<SeqSource, SeqError> {
    let inner = match fmt {
        Format::Fasta => SourceInner::FastaFamily {
            blocks: open_fasta_family(path, FastaKind::Fasta)?,
            kind: FastaKind::Fasta,
        },
        Format::Daemon => SourceInner::FastaFamily {
            blocks: open_fasta_family(path, FastaKind::Daemon)?,
            kind: FastaKind::Daemon,
        },
        Format::Hmmpgmd => SourceInner::FastaFamily {
            blocks: open_fasta_family(path, FastaKind::Hmmpgmd)?,
            kind: FastaKind::Hmmpgmd,
        },
        // esl_sqio_ascii.c:`esl_sqascii_Open` dispatches the non-alignment formats with
        // `default: status = eslEFORMAT; goto ERROR;`, and `eslSQFILE_FMINDEX` has no
        // case -- so `esl_sqfile_Open` fails for every file, and hmmsearch.c turns that
        // into "Sequence file %s is empty or misformatted".
        Format::Fmindex => return Err(SeqError::OpenFormat),
        // esl_sqio.c:`sqfile_open` tries `esl_sqncbi_Open` first for this format, and
        // falls back to `esl_sqascii_Open` -- which returns eslENOTFOUND outright for
        // eslSQFILE_NCBI -- and then to the BLASTDB directory list.
        Format::Ncbi => SourceInner::Ncbi(open_ncbi(path)?),
        _ => {
            if path != "-" {
                std::fs::File::open(path).map_err(|_| SeqError::NotFound)?;
            }
            SourceInner::Deferred {
                path: path.to_string(),
                fmt,
            }
        }
    };
    Ok(SeqSource { inner })
}

impl SeqSource {
    /// Read the database to the end, as hmmsearch's search loop does.
    pub fn read(self) -> Result<Vec<Seq>, SeqError> {
        match self.inner {
            SourceInner::FastaFamily { mut blocks, kind } => {
                read_fasta_family(&mut blocks, kind)
            }
            SourceInner::Ncbi(db) => db.read(),
            SourceInner::Deferred { path, fmt } => read_deferred(&path, fmt),
        }
    }
}

/// Read a FASTA file, or standard input when `path` is `-`.
///
/// C, hmmsearch.c: the target database goes through `esl_sqfile_Open()`, and Easel's
/// `esl_sqfile_Open` treats the filename `"-"` as stdin
/// (`esl_sqio_ascii.c`: `if (strcmp(filename, "-") == 0) { ascii->fp = stdin; ... }`).
/// The header lines print the name as given, i.e. the literal `-`.
pub fn read_fasta(path: &str) -> Result<Vec<Seq>, SeqError> {
    read_fasta_family(&mut open_fasta_family(path, FastaKind::Fasta)?, FastaKind::Fasta)
}

/// A sequence-file format code, as `--tformat` names them.
///
/// C, easel/esl_sqio.c:`esl_sqio_EncodeFormat`, which falls through to
/// easel/esl_msafile.c:`esl_msafile_EncodeFormat` for the aligned formats:
/// ```text
///   if (strcasecmp(fmtstring, "fasta")     == 0) return eslSQFILE_FASTA;
///   if (strcasecmp(fmtstring, "embl")      == 0) return eslSQFILE_EMBL;
///   if (strcasecmp(fmtstring, "genbank")   == 0) return eslSQFILE_GENBANK;
///   if (strcasecmp(fmtstring, "ddbj")      == 0) return eslSQFILE_DDBJ;
///   if (strcasecmp(fmtstring, "uniprot")   == 0) return eslSQFILE_UNIPROT;
///   if (strcasecmp(fmtstring, "ncbi")      == 0) return eslSQFILE_NCBI;
///   if (strcasecmp(fmtstring, "daemon")    == 0) return eslSQFILE_DAEMON;
///   if (strcasecmp(fmtstring, "hmmpgmd")   == 0) return eslSQFILE_HMMPGMD;
///   if (strcasecmp(fmtstring, "fmindex")   == 0) return eslSQFILE_FMINDEX;
///   return esl_msafile_EncodeFormat(fmtstring);
/// ```
/// The comparison is case-insensitive.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    Fasta,
    /// EMBL and UniProtKB share one parser (`config_embl`); they differ only in which
    /// header lines are present, and neither is parsed beyond ID/AC/DE.
    Embl,
    UniProt,
    /// GenBank and DDBJ likewise share `config_genbank`.
    GenBank,
    Ddbj,
    /// Stockholm and its one-block variant Pfam share `esl_msafile_stockholm_Read`;
    /// reading an alignment as a sequence database dealigns each row.
    Stockholm,
    Pfam,
    /// Aligned FASTA: FASTA whose rows carry gaps.
    Afa,
    /// PSI-BLAST's blocked output.
    PsiBlast,
    /// Clustal, and `clustallike` which is Clustal without the `CLUSTAL` word in the
    /// header line (MUSCLE, MAFFT).
    Clustal,
    ClustalLike,
    /// A2M: FASTA-like, with lowercase/`.` marking insert columns.
    A2m,
    /// SELEX: blocks of `<name> <text>` rows plus `#=RF`/`#=CS`/`#=SS`/`#=SA` markup.
    Selex,
    /// PHYLIP interleaved, and `phylips` sequential.
    Phylip,
    PhylipS,
    /// A pressed NCBI BLAST database, version 4: `<base>.pin`/`.phr`/`.psq`, or the
    /// nucleotide `.nin`/`.nhr`/`.nsq`, or a `.pal`/`.nal` alias listing volumes.
    Ncbi,
    /// "Farrar format, hmmpgmd queries: fasta + // terminator" (esl_sqio.h:113).
    Daemon,
    /// "Farrar hmmpgmd database format: fasta + # header" (esl_sqio.h:114).
    Hmmpgmd,
    /// "Pressed FM-index format used in HMMER" (esl_sqio.h:115). `nhmmer` reads these
    /// through its own `--fmindex` path, not through `esl_sqfile_Open`, which has no
    /// case for `eslSQFILE_FMINDEX` at all -- so asserting this format always fails.
    Fmindex,
}

impl Format {
    /// `esl_sqio_IsAlignment(fmt)`, which is `fmt >= 100` -- the aligned formats are the
    /// ones Easel numbers in the `esl_msafile` range and reads through the MSA module.
    ///
    /// It matters beyond dispatch: `esl_sqascii_Open` only takes its stdin branch for a
    /// non-alignment format, so `ascii->do_stdin` -- and with it
    /// `esl_sqfile_IsRewindable` -- is false for an alignment read from standard input.
    pub fn is_alignment(self) -> bool {
        matches!(
            self,
            Format::Stockholm
                | Format::Pfam
                | Format::Afa
                | Format::PsiBlast
                | Format::Clustal
                | Format::ClustalLike
                | Format::A2m
                | Format::Selex
                | Format::Phylip
                | Format::PhylipS
        )
    }
}

/// `esl_sqio_EncodeFormat`, restricted to what is implemented. Unknown *and*
/// not-yet-implemented names are distinguished by the caller so the error message can
/// say which it is.
pub fn encode_format(s: &str) -> Option<Format> {
    let l = s.to_ascii_lowercase();
    match l.as_str() {
        "fasta" => Some(Format::Fasta),
        "embl" => Some(Format::Embl),
        "uniprot" => Some(Format::UniProt),
        "genbank" => Some(Format::GenBank),
        "ddbj" => Some(Format::Ddbj),
        "stockholm" => Some(Format::Stockholm),
        "pfam" => Some(Format::Pfam),
        "afa" => Some(Format::Afa),
        "psiblast" => Some(Format::PsiBlast),
        "clustal" => Some(Format::Clustal),
        "clustallike" => Some(Format::ClustalLike),
        "a2m" => Some(Format::A2m),
        "selex" => Some(Format::Selex),
        "phylip" => Some(Format::Phylip),
        "phylips" => Some(Format::PhylipS),
        "ncbi" => Some(Format::Ncbi),
        "daemon" => Some(Format::Daemon),
        "hmmpgmd" => Some(Format::Hmmpgmd),
        "fmindex" => Some(Format::Fmindex),
        _ => None,
    }
}

/// Open and read a sequence file in one step, or standard input when `path` is `-`.
pub fn read_format(path: &str, fmt: Format) -> Result<Vec<Seq>, SeqError> {
    open_format(path, fmt)?.read()
}

/// The formats whose reader does its own opening.
fn read_deferred(path: &str, fmt: Format) -> Result<Vec<Seq>, SeqError> {
    let ascii = |r: std::io::Result<Vec<Seq>>| r.map_err(io_to_seqerror);
    match fmt {
        Format::Fasta | Format::Daemon | Format::Hmmpgmd | Format::Fmindex | Format::Ncbi => {
            unreachable!("handled by open_format")
        }
        Format::Embl | Format::UniProt => ascii(read_embl(path)),
        Format::GenBank | Format::Ddbj => ascii(read_genbank(path)),
        Format::Stockholm | Format::Pfam => ascii(read_stockholm(path)),
        Format::Afa => ascii(read_afa(path)),
        Format::PsiBlast => ascii(read_blocked(path, Blocked::PsiBlast)),
        Format::Clustal => ascii(read_blocked(path, Blocked::Clustal)),
        Format::ClustalLike => ascii(read_blocked(path, Blocked::ClustalLike)),
        // A2M is FASTA-like; C's reader is more elaborate only because it reconstructs an
        // RF line from the case/dot convention, which a sequence database does not need.
        // The *residues* are the letters either way, upper and lower, so dealigning gives
        // the same sequence.
        Format::A2m => ascii(read_afa(path)),
        Format::Selex => ascii(read_blocked(path, Blocked::Selex)),
        Format::Phylip => ascii(read_phylip(path, false)),
        Format::PhylipS => ascii(read_phylip(path, true)),
    }
}

/// PHYLIP, interleaved (`phylip`) or sequential (`phylips`).
///
/// C, easel/esl_msafile_phylip.c. The header is `<nseq> <alen>`, and the name occupies a
/// fixed-width field -- `namewidth = (afp->fmtd.namewidth ? afp->fmtd.namewidth : 10)`,
/// "default: strict PHYLIP, namewidth 10" -- which is why long names come back truncated.
///
/// `phylip_rectify_input_name` trims the field and turns interior spaces into `_`:
/// ```text
///   for (endpos = n-1; endpos > 0;    endpos--) if (p[endpos] != ' ') break;
///   for (pos    = 0;   pos <= endpos; pos++)    if (p[pos]    != ' ') break;
///   for (          ;   pos <= endpos; pos++)
///     {
///       if (! isgraph(p[pos]) && p[pos] != ' ') return eslEINVAL;
///       namebuf[npos++] = (p[pos] == ' ' ? '_' : p[pos]);
///     }
/// ```
///
/// Interleaved: the name field appears only in the first block (`if (nblocks == 0)`), and
/// each block holds `nseq` lines. Sequential: each sequence's first line carries the name
/// and its lines continue until `alen_stated` residues have accumulated.
pub fn read_phylip(path: &str, sequential: bool) -> std::io::Result<Vec<Seq>> {
    use std::io::BufRead;
    const NAMEWIDTH: usize = 10;
    let map = AminoMap::new();
    let rectify = |field: &str| -> String {
        field.trim_matches(' ').replace(' ', "_")
    };

    let mut lines = open_reader(path)?.lines();
    // Header: `<nseq> <alen>`, after any leading blank lines.
    let mut nseq = 0usize;
    let mut alen = 0usize;
    for line in lines.by_ref() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let mut it = line.split_whitespace();
        nseq = it.next().and_then(|t| t.parse().ok()).unwrap_or(0);
        alen = it.next().and_then(|t| t.parse().ok()).unwrap_or(0);
        break;
    }
    if nseq == 0 {
        return Ok(Vec::new());
    }

    let mut names: Vec<String> = Vec::new();
    let mut rows: Vec<Vec<u8>> = vec![Vec::new(); nseq];
    let mut idx = 0usize;
    let mut first_block = true;

    for line in lines {
        let line = line?;
        if line.trim().is_empty() {
            // "tolerate blank lines only at the end of a block"
            if !sequential && idx >= nseq {
                idx = 0;
                first_block = false;
            }
            continue;
        }
        let take_name = if sequential {
            rows.get(idx).map(|r| r.is_empty()).unwrap_or(false)
        } else {
            first_block
        };
        let text = if take_name {
            if line.len() < NAMEWIDTH {
                continue; // "PHYLIP line too short to find sequence name"
            }
            names.push(rectify(&line[..NAMEWIDTH]));
            &line[NAMEWIDTH..]
        } else {
            &line[..]
        };
        if idx < nseq {
            rows[idx].extend(text.split_whitespace().flat_map(|w| w.bytes()));
        }
        if sequential {
            // A sequence ends once it has `alen` aligned columns.
            if rows[idx].len() >= alen {
                idx += 1;
                if idx >= nseq {
                    break;
                }
            }
        } else {
            idx += 1;
            if idx >= nseq {
                idx = 0;
                first_block = false;
            }
        }
    }

    Ok(names
        .into_iter()
        .zip(rows)
        .map(|(name, row)| Seq {
            name,
            acc: String::new(),
            desc: String::new(),
            dsq: map.digitize(&dealign(&row)),
            roff: 0,
    eoff: 0,
        })
        .collect())
}

/// Aligned FASTA.
///
/// C, easel/esl_msafile_afa.c:`esl_msafile_afa_Read`: the FASTA grammar exactly, with the
/// rows kept aligned. Leading whitespace and blank lines inside a record are tolerated and
/// skipped, and every row must end up the same length.
///
/// It differs from `--tformat fasta` only in that gaps are expected and then removed, so
/// this shares [`read_fasta`]'s shape and adds [`dealign`].
pub fn read_afa(path: &str) -> std::io::Result<Vec<Seq>> {
    use std::io::BufRead;
    let map = AminoMap::new();
    let mut out = Vec::new();
    let (mut name, mut desc) = (String::new(), String::new());
    let mut row: Vec<u8> = Vec::new();
    let mut have = false;
    for line in open_reader(path)?.lines() {
        let line = line?;
        let t = line.trim_start();
        if let Some(hdr) = t.strip_prefix('>') {
            if have {
                out.push(Seq {
                    name: std::mem::take(&mut name),
                    acc: String::new(),
                    desc: std::mem::take(&mut desc),
                    dsq: map.digitize(&dealign(&row)),
                    roff: 0,
    eoff: 0,
                });
                row.clear();
            }
            have = true;
            let hdr = hdr.trim();
            match hdr.split_once(char::is_whitespace) {
                Some((n, d)) => {
                    name = n.to_string();
                    desc = d.trim().to_string();
                }
                None => {
                    name = hdr.to_string();
                    desc = String::new();
                }
            }
        } else {
            row.extend(t.split_whitespace().flat_map(|w| w.bytes()));
        }
    }
    if have {
        out.push(Seq {
            name,
            acc: String::new(),
            desc,
            dsq: map.digitize(&dealign(&row)),
            roff: 0,
    eoff: 0,
        });
    }
    Ok(out)
}

/// Which of the interleaved `<name> <aligned text>` formats is being read.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Blocked {
    PsiBlast,
    Clustal,
    ClustalLike,
    Selex,
}

/// PSI-BLAST and Clustal: blocks of `<name> <aligned text>` lines.
///
/// C, easel/esl_msafile_psiblast.c and esl_msafile_clustal.c. Both tokenize each line the
/// same way -- first non-space run is the name, next non-space run starts the aligned text
/// -- and both concatenate blocks in first-block sequence order:
///
/// ```text
///   for (pos = 0;     pos < n; pos++) if (! isspace(p[pos])) break;
///   name_start = pos;
///   for (pos = pos+1; pos < n; pos++) if (  isspace(p[pos])) break;
///   name_len   = pos - name_start;
///   for (pos = pos+1; pos < n; pos++) if (! isspace(p[pos])) break;
///   seq_start  = pos;
/// ```
///
/// They differ at the end of the row and the end of the block:
///
/// * PSI-BLAST takes the text to the last non-space character of the line
///   (`for (pos = afp->n-1; pos > 0; pos--) if (! isspace(afp->line[pos])) break;`), and
///   blocks are separated by blank lines.
/// * Clustal takes only the *first* whitespace-delimited run ("expect one block; ignore
///   trailing stuff, inc. optional coords"), and each block is closed by a consensus line,
///   detected as a line made only of `" .:*"`:
///   `while (esl_memspn(afp->line, afp->n, " .:*") < afp->n);`
///
/// Clustal also requires a header line whose first token starts with `CLUSTAL` and which
/// contains the word "alignment"; `clustallike` drops the `CLUSTAL` requirement so MUSCLE
/// and MAFFT headers pass.
///
/// Note the consensus-line test means a Clustal *sequence* row consisting only of gaps and
/// `*`/`.`/`:` characters would be mistaken for a consensus line -- that is C's behaviour,
/// not an artifact of this transcription.
fn read_blocked(path: &str, kind: Blocked) -> std::io::Result<Vec<Seq>> {
    use std::io::BufRead;
    let map = AminoMap::new();
    let mut order: Vec<String> = Vec::new();
    let mut rows: std::collections::HashMap<String, Vec<u8>> = std::collections::HashMap::new();
    let clustal = matches!(kind, Blocked::Clustal | Blocked::ClustalLike);
    // Clustal's header line is consumed before any alignment data. SELEX and PSI-BLAST
    // have none.
    let mut need_header = clustal;

    for line in open_reader(path)?.lines() {
        let line = line?;
        if need_header {
            if line.trim().is_empty() {
                continue;
            }
            need_header = false;
            continue;
        }
        let t = line.trim_start();
        if t.is_empty() {
            continue; // block separator
        }
        // Clustal's consensus line: only " .:*" characters.
        if clustal && line.bytes().all(|c| matches!(c, b' ' | b'.' | b':' | b'*')) {
            continue;
        }
        // SELEX markup and comments. C classifies `#=RF`/`#=CS`/`#=SS`/`#=SA` as markup
        // (esl_msafile_selex.c:552-556) and treats any other `#` line as a comment.
        if kind == Blocked::Selex && t.starts_with('#') {
            continue;
        }
        let mut it = t.splitn(2, char::is_whitespace);
        let Some(name) = it.next() else { continue };
        let Some(rest) = it.next() else { continue };
        let text: &str = if clustal {
            // "expect one block; ignore trailing stuff, inc. optional coords"
            rest.trim_start().split_whitespace().next().unwrap_or("")
        } else {
            rest.trim()
        };
        if text.is_empty() {
            continue;
        }
        let e = rows.entry(name.to_string()).or_insert_with(|| {
            order.push(name.to_string());
            Vec::new()
        });
        if clustal {
            e.extend(text.bytes());
        } else {
            e.extend(text.split_whitespace().flat_map(|w| w.bytes()));
        }
    }

    let mut out = Vec::new();
    for name in order {
        if let Some(row) = rows.remove(&name) {
            out.push(Seq {
                name,
                acc: String::new(),
                desc: String::new(),
                dsq: map.digitize(&dealign(&row)),
                roff: 0,
    eoff: 0,
            });
        }
    }
    Ok(out)
}

/// Strip the alignment gaps from one aligned row.
///
/// C reads an alignment as a sequence database by parsing the MSA and then calling
/// `esl_sq_FetchFromMSA` per row (easel/esl_sqio_ascii.c:770), whose digital branch ends
/// with
/// ```text
///   esl_abc_XDealign(sq->abc, sq->dsq,  sq->dsq, &(sq->n));
/// ```
/// and `esl_abc_XDealign` drops exactly the gap and missing symbols:
/// ```text
///   if (! esl_abc_XIsGap(abc, ref_ax[apos]) && ! esl_abc_XIsMissing(abc, ref_ax[apos]) )
///     x[n++] = x[apos];
/// ```
/// For the amino alphabet `"ACDEFGHIKLMNPQRSTVWY-BJZOUX*~"` (K=20, Kp=29) that is code 20
/// (`-`, the gap) and code 28 (`~`, missing). `create_amino` additionally maps `_` and `.`
/// onto the gap (`esl_alphabet_SetEquiv(a, '_', '-')`, `SetEquiv(a, '.', '-')`), so all
/// four go. `*` is code 27, a nonresidue, and is **kept**.
///
/// The dealigning happens before digitizing rather than after, which is equivalent and
/// necessary here: this crate's `AminoMap` has no entry for `.` or `_`, so digitizing
/// first would turn them into X residues instead of dropping them.
fn dealign(row: &[u8]) -> Vec<u8> {
    row.iter()
        .copied()
        .filter(|c| !matches!(c, b'-' | b'_' | b'.' | b'~'))
        .collect()
}

/// Stockholm, and its one-block variant Pfam.
///
/// C, easel/esl_msafile_stockholm.c:`esl_msafile_stockholm_Read`. Only what a sequence
/// database needs is kept: the per-row name, `#=GS ... AC`, `#=GS ... DE`, and the
/// aligned text. `#=GF`, `#=GC`, `#=GR`, `#=GS ... WT` and free comments are recognized
/// and skipped.
///
/// ```text
///   /* Skip leading blank lines in file. EOF here is a normal EOF return. */
///   do {
///     if ( ( status = esl_msafile_GetLine(afp, &p, &n)) != eslOK) goto ERROR;
///   } while (esl_memspn(afp->line, afp->n, " \t") == afp->n ||                  /* skip blank lines             */
///            (esl_memstrpfx(afp->line, afp->n, "#")                             /* and skip comment lines       */
///             && ! esl_memstrpfx(afp->line, afp->n, "# STOCKHOLM")));           /* but stop on Stockholm header */
///
///   if (! esl_memstrpfx(afp->line, afp->n, "# STOCKHOLM 1."))  ESL_XFAIL(eslEFORMAT, afp->errmsg, "missing Stockholm header");
/// ```
///
/// The block structure: a blank line or `//` closes a block, `//` ends the record, and
/// rows of later blocks append to the row of the same name from the first block (C
/// enforces that the order and count match; here a name lookup is enough, since a
/// malformed file would only mis-concatenate rather than corrupt the search).
///
/// A file may hold several records; `esl_sqio_Read` reads the next MSA when the current
/// one is exhausted (easel/esl_sqio_ascii.c:754-757), so all records' rows are returned
/// in order.
pub fn read_stockholm(path: &str) -> std::io::Result<Vec<Seq>> {
    use std::io::BufRead;
    let map = AminoMap::new();
    let mut out = Vec::new();
    // Per record: insertion-ordered rows.
    let mut order: Vec<String> = Vec::new();
    let mut rows: std::collections::HashMap<String, (String, String, Vec<u8>)> =
        std::collections::HashMap::new();
    let mut in_record = false;

    let flush = |order: &mut Vec<String>,
                     rows: &mut std::collections::HashMap<String, (String, String, Vec<u8>)>,
                     out: &mut Vec<Seq>| {
        for name in order.drain(..) {
            if let Some((acc, desc, row)) = rows.remove(&name) {
                out.push(Seq { name, acc, desc, dsq: map.digitize(&dealign(&row)),
    roff: 0,
    eoff: 0,
});
            }
        }
        rows.clear();
    };

    for line in open_reader(path)?.lines() {
        let line = line?;
        let t = line.trim_start_matches([' ', '\t']);
        if !in_record {
            // Blank lines and comments are skipped until the magic header.
            if t.starts_with("# STOCKHOLM 1.") {
                in_record = true;
            }
            continue;
        }
        if t.is_empty() {
            continue; // end of block; rows continue in the next one
        }
        if t.starts_with("//") {
            flush(&mut order, &mut rows, &mut out);
            in_record = false;
            continue;
        }
        if t.starts_with('#') {
            if let Some(rest) = t.strip_prefix("#=GS") {
                // `#=GS <seqname> <tag> <annotation>`. C tokenizes the first three
                // fields with `esl_memtok` and then, for DE, hands the *remainder of the
                // line* to `esl_msa_SetSeqDescription(msa, seqidx, p, n)` verbatim -- so
                // runs of spaces inside a description are preserved. (UniProt
                // descriptions really do contain them: `... sevenless;          EC=...`.)
                // Splitting on whitespace and rejoining collapses them and changes the
                // description in every output that prints it.
                let after_name = rest.trim_start();
                let (name, after_name) = match after_name.split_once(char::is_whitespace) {
                    Some(v) => v,
                    None => continue,
                };
                let after_tag = after_name.trim_start();
                let (tag, val) = match after_tag.split_once(char::is_whitespace) {
                    Some((t, v)) => (t, v.trim_start()),
                    None => (after_tag, ""),
                };
                let e = rows.entry(name.to_string()).or_insert_with(|| {
                    order.push(name.to_string());
                    (String::new(), String::new(), Vec::new())
                });
                match tag {
                    // "should have only one field, the accession"
                    "AC" => e.0 = val.split_whitespace().next().unwrap_or("").to_string(),
                    "DE" => e.1 = val.trim_end().to_string(),
                    _ => {}
                }
            }
            // #=GF / #=GC / #=GR / free comments carry nothing a sequence db needs.
            continue;
        }
        // A sequence line: `<name> <aligned text>`.
        let mut it = t.splitn(2, char::is_whitespace);
        let Some(name) = it.next() else { continue };
        let Some(text) = it.next() else { continue };
        let e = rows.entry(name.to_string()).or_insert_with(|| {
            order.push(name.to_string());
            (String::new(), String::new(), Vec::new())
        });
        e.2.extend(text.split_whitespace().flat_map(|w| w.bytes()));
    }
    // C fails with "missing // terminator after MSA"; be lenient and return what parsed,
    // since a truncated database is more useful than none.
    flush(&mut order, &mut rows, &mut out);
    Ok(out)
}

/// The residue-character rules the EMBL and GenBank parsers share.
///
/// C, easel/esl_sqio_ascii.c:`inmap_embl` / `inmap_genbank` (identical apart from a
/// comment):
/// ```text
///   for (x = '0'; x <= '9'; x++)
///     sqfp->inmap[x] = eslDSQ_IGNORED;    /* EMBL DNA sequence format puts coordinates after each line */
///   sqfp->inmap['*']  = '*';              /* accept * as a nonresidue/stop codon character */
///   sqfp->inmap[' ']  = eslDSQ_IGNORED;
///   sqfp->inmap['\t'] = eslDSQ_IGNORED;
///   sqfp->inmap['\n'] = eslDSQ_IGNORED;
///   sqfp->inmap['\r'] = eslDSQ_IGNORED;  /* DOS eol compatibility */
///   sqfp->inmap['/']  = eslDSQ_EOD;
/// ```
/// and, in digital mode, `sqfp->inmap['-'] = eslDSQ_ILLEGAL` -- "don't let the abc_inmap
/// map the gap char; this is an ungapped file format". So digits and whitespace are
/// dropped, `*` is kept, and `-` is *not* accepted.
fn push_flatfile_residues(buf: &mut Vec<u8>, line: &str) {
    for c in line.bytes() {
        match c {
            b'0'..=b'9' | b' ' | b'\t' | b'\r' | b'\n' => {}
            _ => buf.push(c),
        }
    }
}

/// Lines, plus the byte offset each one starts at.
///
/// The line-based readers need offsets for one reason: `sq->roff` is what an SSI index
/// stores, and `--restrictdb_stkey` resolves a key to that offset. `BufRead::lines()`
/// cannot supply it, because it has already discarded whether the terminator was `\n` or
/// `\r\n`.
struct OffsetLines<R> {
    r: R,
    off: u64,
    buf: Vec<u8>,
}

impl<R: std::io::BufRead> OffsetLines<R> {
    fn new(r: R) -> Self {
        OffsetLines { r, off: 0, buf: Vec::new() }
    }

    /// The next line with its terminator stripped, and the offset it began at.
    fn next_line(&mut self) -> std::io::Result<Option<(String, u64)>> {
        self.buf.clear();
        let start = self.off;
        let n = self.r.read_until(b'\n', &mut self.buf)?;
        if n == 0 {
            return Ok(None);
        }
        self.off += n as u64;
        while matches!(self.buf.last(), Some(b'\n') | Some(b'\r')) {
            self.buf.pop();
        }
        Ok(Some((
            String::from_utf8_lossy(&self.buf).into_owned(),
            start,
        )))
    }
}

fn open_reader(path: &str) -> std::io::Result<std::io::BufReader<Box<dyn std::io::Read>>> {
    let src: Box<dyn std::io::Read> = if path == "-" {
        Box::new(std::io::stdin())
    } else {
        Box::new(std::fs::File::open(path)?)
    };
    Ok(std::io::BufReader::new(src))
}

/// EMBL / UniProtKB flat file.
///
/// C, easel/esl_sqio_ascii.c:`header_embl`. Records begin with an `ID` line and end with
/// `//`; only three header lines are parsed:
///
/// ```text
///   /* ID line is defined as:
///    *     ID   ENTRY_NAME DATA_CLASS; MOLECULE_TYPE; SEQUENCE_LENGTH.
///    * We're only after the ENTRY_NAME.
///    * Examples:
///    *  ID   SNRPA_DROME    STANDARD;      PRT;   216 AA.
///    *  ID   SNRPA_DROME             Reviewed;         216 AA.
///    *  ID   X06347; SV 1; linear; mRNA; STD; HUM; 1209 BP.
///    */
///   s = ascii->buf+5;
///   esl_strtok(&s, " ;", &tok)
/// ```
/// -- so the name is the first token delimited by space *or* semicolon.
///
/// ```text
///   /*   AC   P43332; Q9W4D7;
///    *  Note that Easel only stores primary accessions.
///    *  Because there can be more than one accession line, we check to
///    *  see if the accession is already set before storing a line.
///    */
///   if (strncmp(ascii->buf, "AC   ", 5) == 0 && sq->acc[0] == '\0')
/// ```
///
/// ```text
///   /* "...In cases where more than one DE line is required, the text is
///    * only divided between words and only the last DE line is
///    * terminated by a period."
///    * We'll make no attempt to parse the structured UniProt description header,
///    * for the moment.
///    */
///   if (strncmp(ascii->buf, "DE   ", 5) == 0) { ... esl_sq_AppendDesc(sq, s) }
/// ```
/// `esl_sq_AppendDesc` joins successive DE lines with a single space.
///
/// Sequence data starts after the `SQ   ` line: "We don't parse this line; we just look
/// for it as the last line before the sequence starts."
pub fn read_embl(path: &str) -> std::io::Result<Vec<Seq>> {
    let map = AminoMap::new();
    let mut out = Vec::new();
    let mut name = String::new();
    let mut acc = String::new();
    let mut desc = String::new();
    let mut buf: Vec<u8> = Vec::new();
    // Three states: outside a record, in its header, in its sequence data.
    let (mut in_rec, mut in_seq) = (false, false);
    // `header_embl`: `sq->roff = ascii->boff;  /* record the offset of the ID line */`
    let mut roff = 0u64;
    let mut lines = OffsetLines::new(open_reader(path)?);
    while let Some((line, off)) = lines.next_line()? {
        if !in_rec {
            if let Some(rest) = line.strip_prefix("ID   ") {
                roff = off;
                // first token delimited by space or ';'
                name = rest
                    .split([' ', ';'])
                    .find(|t| !t.is_empty())
                    .unwrap_or("")
                    .to_string();
                acc.clear();
                desc.clear();
                buf.clear();
                in_rec = true;
                in_seq = false;
            }
            continue;
        }
        if !in_seq {
            // `strncmp(buf,"AC   ",5)==0 && sq->acc[0]=='\0'` -- only the primary
            // accession, i.e. the first `;`-delimited token of the *first* AC line.
            if acc.is_empty() {
                if let Some(rest) = line.strip_prefix("AC   ") {
                    acc = rest.split(';').next().unwrap_or("").trim().to_string();
                }
            }
            if let Some(rest) = line.strip_prefix("DE   ") {
                // esl_sq_AppendDesc: successive lines joined by a single space.
                if !desc.is_empty() {
                    desc.push(' ');
                }
                desc.push_str(rest.trim_end());
            } else if line.starts_with("SQ   ") {
                in_seq = true;
            }
            continue;
        }
        // `sqfp->inmap['/'] = eslDSQ_EOD` -- the `//` terminator ends the data.
        if line.starts_with("//") {
            out.push(Seq {
                name: std::mem::take(&mut name),
                acc: std::mem::take(&mut acc),
                desc: std::mem::take(&mut desc),
                dsq: map.digitize(&buf),
                roff,
                eoff: lines.off.saturating_sub(1),
            });
            buf.clear();
            in_rec = false;
            in_seq = false;
            continue;
        }
        push_flatfile_residues(&mut buf, &line);
    }
    Ok(out)
}

/// GenBank / DDBJ flat file.
///
/// C, easel/esl_sqio_ascii.c:`header_genbank`. The record starts at `LOCUS   ` (note
/// three spaces in the test but the name is read from *column 12*), the accession comes
/// from `VERSION`, the description from `DEFINITION`, and the sequence follows `ORIGIN`:
///
/// ```text
///   while (strncmp(ascii->buf, "LOCUS   ", 8) != 0) { ... }
///   s = ascii->buf+12;
///   esl_strtok(&s, " ", &tok)
///   ...
///   if (strncmp(ascii->buf, "VERSION   ", 10) == 0) { s = ascii->buf+12; esl_strtok(&s, " \t\n", &tok); ... }
///   if (strncmp(ascii->buf, "DEFINITION ", 11) == 0) { s = ascii->buf+12; esl_strchop(...); esl_sq_AppendDesc(sq, s); }
///   } while (strncmp(ascii->buf, "ORIGIN", 6) != 0);
/// ```
///
/// The fixed column-12 offsets are C's, and are why a `LOCUS` line whose name starts
/// earlier or later than column 13 parses differently than one might expect. Transcribed
/// as written.
pub fn read_genbank(path: &str) -> std::io::Result<Vec<Seq>> {
    let map = AminoMap::new();
    let mut out = Vec::new();
    let mut name = String::new();
    let mut acc = String::new();
    let mut desc = String::new();
    let mut buf: Vec<u8> = Vec::new();
    let (mut in_rec, mut in_seq) = (false, false);
    // `ascii->buf+12` on a line shorter than 12 chars would read past the end in C; here
    // it yields an empty field.
    let from12 = |l: &str| -> String { l.chars().skip(12).collect() };
    // `header_genbank`: `sq->roff = ascii->boff;  /* record the disk offset to the LOCUS line */`
    let mut roff = 0u64;
    let mut lines = OffsetLines::new(open_reader(path)?);
    while let Some((line, off)) = lines.next_line()? {
        if !in_rec {
            if line.starts_with("LOCUS   ") {
                roff = off;
                let rest = from12(&line);
                name = rest
                    .split(' ')
                    .find(|t| !t.is_empty())
                    .unwrap_or("")
                    .to_string();
                acc.clear();
                desc.clear();
                buf.clear();
                in_rec = true;
                in_seq = false;
            }
            continue;
        }
        if !in_seq {
            // "Optional VERSION line is parsed as \"accession\"." -- first
            // whitespace-delimited token from column 12.
            if line.starts_with("VERSION   ") {
                let rest = from12(&line);
                acc = rest.split_whitespace().next().unwrap_or("").to_string();
            }
            if line.starts_with("DEFINITION ") {
                let rest = from12(&line);
                if !desc.is_empty() {
                    desc.push(' ');
                }
                desc.push_str(rest.trim_end());
            } else if line.starts_with("ORIGIN") {
                in_seq = true;
            }
            continue;
        }
        if line.starts_with("//") {
            out.push(Seq {
                name: std::mem::take(&mut name),
                acc: std::mem::take(&mut acc),
                desc: std::mem::take(&mut desc),
                dsq: map.digitize(&buf),
                roff,
                eoff: lines.off.saturating_sub(1),
            });
            buf.clear();
            in_rec = false;
            in_seq = false;
            continue;
        }
        push_flatfile_residues(&mut buf, &line);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn read_bytes(data: &[u8], kind: FastaKind) -> Result<Vec<Seq>, SeqError> {
        let mut b = open_fasta_stream(Box::new(std::io::Cursor::new(data.to_vec())), kind)?;
        read_fasta_family(&mut b, kind)
    }

    fn residues(s: &Seq) -> Vec<u8> {
        s.dsq[1..s.dsq.len() - 1].to_vec()
    }

    /// Expectations here come from running C HMMER 3.4 on the same inputs, not from
    /// reasoning about the transcription: `hmmsearch testdata/globins4.hmm <file>` reports
    /// "8 residues searched" for each of `*`, `.`, `_` and `~`, and
    /// "Parse failed ... Line 2: illegal character" for `1` and `-`.
    #[test]
    fn fasta_inmap_matches_c() {
        use crate::alphabet::{GAP, MISSING, NONRESIDUE};
        let one = |body: &str| -> Result<Vec<u8>, SeqError> {
            read_bytes(format!(">s1\n{body}\n").as_bytes(), FastaKind::Fasta)
                .map(|v| residues(&v[0]))
        };
        // `*` is a residue (p7_NONRESIDUE), and `.`/`_`/`~` are legal too -- the file
        // inmap only overrides `-`.
        assert_eq!(one("AC*").unwrap(), vec![0, 1, NONRESIDUE]);
        assert_eq!(one("AC.").unwrap(), vec![0, 1, GAP]);
        assert_eq!(one("AC_").unwrap(), vec![0, 1, GAP]);
        assert_eq!(one("AC~").unwrap(), vec![0, 1, MISSING]);
        // Lower case maps to the same codes; spaces and tabs are ignored.
        assert_eq!(one("ac \tD").unwrap(), vec![0, 1, 2]);
        // A digit and the gap character are format errors.
        for bad in ["AC1", "AC-"] {
            match one(bad) {
                Err(SeqError::Parse(m)) => assert_eq!(
                    String::from_utf8_lossy(&m),
                    format!("Line 2: illegal character {}", &bad[2..])
                ),
                other => panic!("{bad}: expected a parse error, got {other:?}"),
            }
        }
    }

    /// C stores the description verbatim from the first non-blank character to the end of
    /// the line, so trailing spaces survive into `--tblout`'s description column.
    #[test]
    fn fasta_header_keeps_trailing_description_space() {
        let seqs = read_bytes(b">s1  desc  with   trail   \nAC\n", FastaKind::Fasta).unwrap();
        assert_eq!(seqs[0].name, "s1");
        assert_eq!(seqs[0].desc, "desc  with   trail   ");
    }

    /// `header_fasta` delimits the description at ctrl-A as well as at end of line
    /// ("Patched to deal with NCBI NR desclines [SRE:H1/82]").
    #[test]
    fn fasta_header_stops_description_at_ctrl_a() {
        let seqs = read_bytes(b">s1 first\x01second\nAC\n", FastaKind::Fasta).unwrap();
        assert_eq!(seqs[0].desc, "first");
    }

    /// `end_daemon`'s second skip loop reads one character past the newline whenever the
    /// block does not end there, so C accepts a lone record but rejects two back-to-back
    /// -- and accepts them again with a filler character in between. Verified against C:
    /// `hmmsearch --tformat daemon` reports 1 sequence, then
    /// "Line 3: unexpected char s; expected FASTA to start with >", then 2 sequences.
    #[test]
    fn daemon_terminator_overshoots_by_one() {
        let one = read_bytes(b">s1 d1\nACDEFGH\n//\n", FastaKind::Daemon).unwrap();
        assert_eq!(one.len(), 1);
        assert_eq!(residues(&one[0]).len(), 7);

        match read_bytes(b">s1 d1\nAC\n//\n>s2 d2\nMK\n//\n", FastaKind::Daemon) {
            Err(SeqError::Parse(m)) => assert_eq!(
                String::from_utf8_lossy(&m),
                "Line 3: unexpected char s; expected FASTA to start with >"
            ),
            other => panic!("expected a parse error, got {other:?}"),
        }

        let two = read_bytes(b">s1 d1\nAC\n//\nX>s2 d2\nMK\n//\n", FastaKind::Daemon).unwrap();
        assert_eq!(two.len(), 2);
        assert_eq!(two[1].name, "s2");
    }

    /// A daemon record with no `//` is a truncated file, because `config_daemon` sets
    /// `eof_is_ok = FALSE`.
    #[test]
    fn daemon_requires_the_terminator() {
        match read_bytes(b">s1 d1\nACDEFGH\n", FastaKind::Daemon) {
            Err(SeqError::Parse(m)) => {
                assert_eq!(String::from_utf8_lossy(&m), "Unexpected EOF; file truncated?")
            }
            other => panic!("expected a parse error, got {other:?}"),
        }
    }

    /// hmmpgmd is FASTA behind one `#` line. C tolerates leading whitespace before the
    /// `#`, fails the *open* if the first non-space character is something else, and
    /// returns `eslEOF` out of the open if the header line has no terminating newline.
    #[test]
    fn hmmpgmd_header_matches_c() {
        let seqs = read_bytes(b"\n\n  #2 4\n>s1 d1\nAC\n>s2 d2\nDE\n", FastaKind::Hmmpgmd).unwrap();
        assert_eq!(seqs.len(), 2);
        assert_eq!(seqs[1].name, "s2");
        // A header line and nothing else is an empty database, not an error.
        assert_eq!(read_bytes(b"#2 4\n", FastaKind::Hmmpgmd).unwrap().len(), 0);
        assert!(matches!(
            read_bytes(b">s1\nAC\n", FastaKind::Hmmpgmd),
            Err(SeqError::OpenFormat)
        ));
        assert!(matches!(
            read_bytes(b"#2 4", FastaKind::Hmmpgmd),
            Err(SeqError::OpenStatus(3))
        ));
        assert!(matches!(
            read_bytes(b"  \n \n", FastaKind::Hmmpgmd),
            Err(SeqError::OpenStatus(3))
        ));
    }

    /// The reader's block size is observable, so a record that straddles one has to come
    /// out whole: `eslREADBUFSIZE` is 4096, and 60-residue lines put the boundary inside
    /// a sequence.
    #[test]
    fn records_straddle_block_boundaries() {
        let mut file = String::new();
        for i in 0..100 {
            file.push_str(&format!(">s{i} desc {i}\n"));
            for _ in 0..5 {
                file.push_str("ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY\n");
            }
        }
        assert!(file.len() > 4 * READBUFSIZE);
        let seqs = read_bytes(file.as_bytes(), FastaKind::Fasta).unwrap();
        assert_eq!(seqs.len(), 100);
        for (i, s) in seqs.iter().enumerate() {
            assert_eq!(s.name, format!("s{i}"));
            assert_eq!(s.desc, format!("desc {i}"));
            assert_eq!(residues(s).len(), 300);
        }
    }

    /// An empty file is an open-time format error, which hmmsearch reports as
    /// "Sequence file %s is empty or misformatted".
    #[test]
    fn empty_file_fails_the_open() {
        assert!(matches!(
            read_bytes(b"", FastaKind::Fasta),
            Err(SeqError::OpenFormat)
        ));
    }
}
