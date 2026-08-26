//! NCBI BLAST database input, `--tformat ncbi`.
//!
//! C, easel/esl_sqio_ncbi.c. A BLAST database is three files sharing a base name: the
//! index `.pin`, the headers `.phr` and the packed sequences `.psq` (`.nin`/`.nhr`/`.nsq`
//! for nucleotide). The index gives the record boundaries in the other two; the headers
//! are BER-encoded ASN.1 `Blast-def-line-set` structures; the sequences are raw
//! NCBIstdaa codes with a zero byte after each.
//!
//! Only version 4 is supported, exactly as C is: `sqncbi_ParseIndexFile` fails on
//! anything else, so a database built by a modern `makeblastdb` without
//! `-blastdb_version 4` is rejected by C too.
//!
//! Nucleotide databases are read too, because C reads them: `sqncbi_SetDigital` ignores
//! the alphabet it is handed (`return eslOK;` and nothing else), so hmmsearch searches a
//! `.nin` database with its amino model and finds whatever the DNA codes happen to mean
//! as amino codes. That is deterministic, and matching it is the point.

use crate::alphabet;
use crate::seqio::{Seq, SeqError};

/// `#define NCBI_VERSION_4 4`
const VERSION_4: u32 = 4;
/// `#define NCBI_DNA_DB 0` / `#define NCBI_AMINO_DB 1`
const DNA_DB: u32 = 0;
const AMINO_DB: u32 = 1;
/// `#define MAX_DB_VOLUMES 100` (esl_sqio.h).
const MAX_DB_VOLUMES: usize = 100;

/// A parsed `.pin` index plus the two data files it describes.
struct Volume {
    /// `.phr` contents.
    hdr: Vec<u8>,
    /// `.psq` contents.
    seq: Vec<u8>,
    /// The header index table, `num_seq + 1` offsets into `.phr`.
    hdr_indexes: Vec<u32>,
    /// The sequence index table, `num_seq + 1` offsets into `.psq`.
    seq_indexes: Vec<u32>,
    /// The ambiguity index table, present only in a nucleotide index:
    /// ```text
    ///   if (ncbi->alphatype == eslDNA) {
    ///     ncbi->amb_off = ncbi->seq_off + sizeof(uint32_t) * (ncbi->num_seq + 1);
    /// ```
    amb_indexes: Vec<u32>,
    num_seq: usize,
}

/// `sqncbi_ParseIndexFile`, which reads the fixed header of a `.pin`/`.nin` and locates
/// the two index tables:
/// ```text
///   if (fread(&info[0], sizeof(uint32_t), 3, ncbi->fppin) != 3) status = eslFAIL;
///   if (htobe32(info[0]) != NCBI_VERSION_4)                     status = eslEFORMAT;
///   if (htobe32(info[1]) != dbtype)                             status = eslEUNIMPLEMENTED;
///   ...
///   /* read the database title */
///   len = htobe32(info[2]);
///   ...
///   /* read the database time stamp */
///   ...
///   /* read in database stats */
///   if (fread(&info[0], sizeof(uint32_t), 4, ncbi->fppin) != 4) { status = eslFAIL; goto ERROR; }
///   ncbi->num_seq   = htobe32(info[0]);
///   memcpy(&ncbi->total_res, info+1, sizeof(uint64_t));
///   ncbi->max_seq   = htobe32(info[3]);
///
///   /* save the offsets to the index tables */
///   ncbi->hdr_off = ftell(ncbi->fppin);
///   ncbi->seq_off = ncbi->hdr_off + sizeof(uint32_t) * (ncbi->num_seq + 1);
/// ```
/// The counts are big-endian but `total_res` is `memcpy`'d without a swap -- and that is
/// right, because `makeblastdb` writes that one field in host order. It is only used for
/// the volume accounting that hmmsearch never prints, so it is not kept here.
///
/// The two integer statuses are distinct and both surface: `eslEFORMAT` (wrong version)
/// becomes hmmsearch's "is empty or misformatted", `eslEUNIMPLEMENTED` (24) becomes
/// "Unexpected error 24 opening sequence file".
fn parse_index_file(pin: &[u8], dbtype: u32) -> Result<(usize, usize, usize), SeqError> {
    // Returns (num_seq, hdr_off, seq_off).
    let be32 = |o: usize| -> Option<u32> {
        pin.get(o..o + 4)
            .map(|b| u32::from_be_bytes([b[0], b[1], b[2], b[3]]))
    };
    // fread of 3 uint32_t; a short read is eslFAIL, which hmmsearch reports as
    // "Unexpected error 1 opening sequence file".
    let (version, ftype, title_len) = match (be32(0), be32(4), be32(8)) {
        (Some(a), Some(b), Some(c)) => (a, b, c),
        _ => return Err(SeqError::OpenStatus(1)), // eslFAIL
    };
    if version != VERSION_4 {
        return Err(SeqError::OpenFormat);
    }
    if ftype != dbtype {
        return Err(SeqError::OpenStatus(24)); // eslEUNIMPLEMENTED
    }
    let mut o = 12usize;
    // The title and the timestamp are read but only used by esl-sfetch's -n output.
    if pin.len() < o + title_len as usize {
        return Err(SeqError::OpenStatus(1));
    }
    o += title_len as usize;
    let stamp_len = match be32(o) {
        Some(v) => v as usize,
        None => return Err(SeqError::OpenStatus(1)),
    };
    o += 4;
    if pin.len() < o + stamp_len {
        return Err(SeqError::OpenStatus(1));
    }
    o += stamp_len;
    // num_seq, total_res (8 bytes, host order), max_seq
    if pin.len() < o + 16 {
        return Err(SeqError::OpenStatus(1));
    }
    let num_seq = be32(o).unwrap() as usize;
    o += 16;
    let hdr_off = o;
    let seq_off = hdr_off + 4 * (num_seq + 1);
    Ok((num_seq, hdr_off, seq_off))
}

/// `sqncbi_DbOpen`: open `<base>.Xin`, `<base>.Xhr` and `<base>.Xsq`, where `X` is `p`
/// for protein and `n` for nucleotide, then parse the index.
///
/// ```text
///   strcpy(name+len, ".Xin");
///   name[len+1] = (dbtype == NCBI_DNA_DB) ? 'n' : 'p';
///   if ((ncbi->fppin = fopen(name, "rb")) == NULL) { status = eslENOTFOUND; goto ERROR; }
///   ... .Xhr ... .Xsq ...
///   if ((status = sqncbi_ParseIndexFile(ncbi, dbtype)) != eslOK) goto ERROR;
/// ```
/// C keeps the files open and seeks per record; this reads each whole, which is the same
/// for a search that visits every record once.
fn db_open(base: &str, dbtype: u32) -> Result<Volume, SeqError> {
    let x = if dbtype == DNA_DB { 'n' } else { 'p' };
    let slurp = |ext: &str| -> Result<Vec<u8>, SeqError> {
        std::fs::read(format!("{base}.{x}{ext}")).map_err(|_| SeqError::NotFound)
    };
    let pin = slurp("in")?;
    let hdr = slurp("hr")?;
    let seq = slurp("sq")?;
    let (num_seq, hdr_off, seq_off) = parse_index_file(&pin, dbtype)?;

    let table = |off: usize| -> Result<Vec<u32>, SeqError> {
        let end = off + 4 * (num_seq + 1);
        let b = pin.get(off..end).ok_or(SeqError::OpenFormat)?;
        Ok(b.chunks_exact(4)
            .map(|c| u32::from_be_bytes([c[0], c[1], c[2], c[3]]))
            .collect())
    };
    Ok(Volume {
        hdr,
        seq,
        hdr_indexes: table(hdr_off)?,
        seq_indexes: table(seq_off)?,
        amb_indexes: if dbtype == DNA_DB {
            table(seq_off + 4 * (num_seq + 1))?
        } else {
            Vec::new()
        },
        num_seq,
    })
}

/// `sqncbi_AliasOpen`: a `.pal`/`.nal` alias file lists the volumes of a split database
/// on a `DBLIST` line.
///
/// ```text
///   /* find the DBLIST directive */
///   ...
///   if (strncmp(buffer, "DBLIST", 6) == 0 && isspace(buffer[6])) { done = 1; ptr = buffer + 6; }
///   ...
///   /* if the name in the DBLIST directive was ablsolute, do not
///    * use the working directory as a prefix. */
///   if (dbname[len] == eslDIRSLASH) { dbptr = dbname + len; } else { dbptr = dbname; }
///   status = sqncbi_DbOpen(ncbi, dbptr, dbtype);
/// ```
/// Volume names are resolved against the alias file's own directory unless they are
/// absolute. C reads the `DBLIST` line through an 80-byte `fgets` buffer and continues
/// across refills, so a name may straddle the boundary; reading the whole line here is
/// equivalent because the tokens are the same either way.
fn alias_open(base: &str, dbtype: u32) -> Result<Vec<Volume>, SeqError> {
    let x = if dbtype == DNA_DB { 'n' } else { 'p' };
    let text =
        std::fs::read_to_string(format!("{base}.{x}al")).map_err(|_| SeqError::NotFound)?;

    // "remove the filename keeping the path"
    let dir = match base.rfind('/') {
        Some(i) => &base[..=i],
        None => "",
    };
    let mut list: Option<&str> = None;
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix("DBLIST") {
            if rest.starts_with([' ', '\t']) {
                list = Some(rest);
                break;
            }
        }
    }
    let list = list.ok_or(SeqError::OpenFormat)?;

    let mut vols = Vec::new();
    for tok in list.split_whitespace() {
        if vols.len() >= MAX_DB_VOLUMES {
            break;
        }
        let name = if tok.starts_with('/') {
            tok.to_string()
        } else {
            format!("{dir}{tok}")
        };
        vols.push(db_open(&name, dbtype)?);
    }
    if vols.is_empty() {
        return Err(SeqError::OpenFormat);
    }
    Ok(vols)
}

/// `sqncbi_Open`'s search order:
/// ```text
///   /* first try to open a single protein database */
///   status = sqncbi_DbOpen(ncbi, filename, NCBI_AMINO_DB);
///   /* if the database was not found, look for protein volume */
///   if (status == eslENOTFOUND) status = sqncbi_AliasOpen(ncbi, filename, NCBI_AMINO_DB);
///   /* if nothing so far, try a dna database */
///   if (status == eslENOTFOUND) status = sqncbi_DbOpen(ncbi, filename, NCBI_DNA_DB);
///   /* still nothing, look for dna volume */
///   if (status == eslENOTFOUND) status = sqncbi_AliasOpen(ncbi, filename, NCBI_DNA_DB);
/// ```
/// Only `eslENOTFOUND` falls through to the next attempt, so a `.pin` that exists but is
/// the wrong version stops the search with that error rather than being retried as
/// nucleotide.
pub struct NcbiDb {
    volumes: Vec<Volume>,
    /// `ncbi->alphatype`. hmmsearch only ever wants amino, but a nucleotide database has
    /// to be *read* far enough to report what C reports.
    is_dna: bool,
}

pub fn open(base: &str) -> Result<NcbiDb, SeqError> {
    let attempt = |r: Result<Vec<Volume>, SeqError>| -> Option<Result<Vec<Volume>, SeqError>> {
        match r {
            Err(SeqError::NotFound) => None,
            other => Some(other),
        }
    };
    if let Some(r) = attempt(db_open(base, AMINO_DB).map(|v| vec![v])) {
        return r.map(|volumes| NcbiDb { volumes, is_dna: false });
    }
    if let Some(r) = attempt(alias_open(base, AMINO_DB)) {
        return r.map(|volumes| NcbiDb { volumes, is_dna: false });
    }
    if let Some(r) = attempt(db_open(base, DNA_DB).map(|v| vec![v])) {
        return r.map(|volumes| NcbiDb { volumes, is_dna: true });
    }
    if let Some(r) = attempt(alias_open(base, DNA_DB)) {
        return r.map(|volumes| NcbiDb { volumes, is_dna: true });
    }
    Err(SeqError::NotFound)
}

/// `inmap_ncbi_amino`, which maps NCBIstdaa codes to Easel's amino codes:
/// ```text
///   const char *ncbisym = "-ABCDEFGHIKLMNPQRSTVWXYZU*OJ";
///   for (x =  0;  x < 128;  x++) sqfp->inmap[x] = eslDSQ_ILLEGAL;
///   /* for each letter in the ncbi alphabet, find that letter in the
///    * hmmer alphabet and map the translation. */
///   for (x = 0; x < strlen(ncbisym); ++x) {
///     for (y = 0; y < strlen(abc->sym); ++y) {
///       if (ncbisym[x] == abc->sym[y]) { sqfp->inmap[x] = y; break; }
///     }
///     ...
/// ```
/// Note that the index into the map is `x`, the *position* in `ncbisym` -- not the
/// character. That is deliberate: a `.psq` stores NCBIstdaa numeric codes, and `x` is
/// the code. Codes 28 and above stay `eslDSQ_ILLEGAL`, and `read_amino` copies that
/// straight into the dsq without checking.
const NCBI_STDAA: &[u8; 28] = b"-ABCDEFGHIKLMNPQRSTVWXYZU*OJ";

fn amino_inmap() -> [u8; 256] {
    let mut m = [254u8; 256]; // eslDSQ_ILLEGAL
    for (code, &c) in NCBI_STDAA.iter().enumerate() {
        if let Some(y) = alphabet::SYM.iter().position(|&s| s == c) {
            m[code] = y as u8;
        }
    }
    m
}

/// `inmap_ncbi_dna`, the same construction against NCBI's 4-bit nucleotide alphabet and
/// Easel's DNA alphabet:
/// ```text
///   const char *ncbisym = "-ACMGRSVTWYHKDBN";
///   if ((abc = esl_alphabet_Create(eslDNA)) == NULL) return eslEMEM;
/// ```
/// with `esl_alphabet.c:create_dna` supplying
/// `esl_alphabet_CreateCustom("ACGT-RYMKSWHBVDN*~", 4, 18)`.
///
/// hmmsearch never revisits this: `sqncbi_SetDigital` does not check the alphabet, so the
/// codes this produces are *DNA* codes that the amino pipeline then reads as amino codes.
const NCBI_NA4: &[u8; 16] = b"-ACMGRSVTWYHKDBN";
const ESL_DNA_SYM: &[u8; 18] = b"ACGT-RYMKSWHBVDN*~";

fn dna_inmap() -> [u8; 256] {
    let mut m = [254u8; 256];
    for (code, &c) in NCBI_NA4.iter().enumerate() {
        if let Some(y) = ESL_DNA_SYM.iter().position(|&s| s == c) {
            m[code] = y as u8;
        }
    }
    m
}

/// The string form of a saved header field: an offset and length into the header buffer.
///
/// C keeps raw pointers into `ncbi->hdr_buf` and zero-terminates them in place at the end
/// of `parse_def_line` (`ncbi->name_ptr[ncbi->name_size] = '\0'`). Those writes only ever
/// land on the first byte of an ASN.1 terminator the parser has already consumed, so
/// slicing by offset and length is equivalent and leaves the buffer intact.
#[derive(Clone, Copy)]
struct Field {
    off: usize,
    len: usize,
}

/// The parser state of `ESL_SQNCBI_DATA` that `reset_header_values` clears:
/// ```text
///   ncbi->name_ptr    = NULL;  ncbi->name_size   = 0;
///   ncbi->acc_ptr     = NULL;  ncbi->acc_size    = 0;
///   ncbi->int_id      = -1;
///   ncbi->str_id_ptr  = NULL;  ncbi->str_id_size = 0;
/// ```
/// It is reset once per *header*, not per def line, which is how "only the information
/// from the first usable defition will be used" is implemented.
struct HdrParse<'a> {
    buf: &'a [u8],
    pos: usize,
    name: Option<Field>,
    acc: Option<Field>,
    str_id: Option<Field>,
    int_id: i32,
}

impl<'a> HdrParse<'a> {
    fn new(buf: &'a [u8]) -> Self {
        HdrParse {
            buf,
            pos: 0,
            name: None,
            acc: None,
            str_id: None,
            int_id: -1,
        }
    }

    fn text(&self, f: Field) -> String {
        String::from_utf8_lossy(&self.buf[f.off..f.off + f.len]).into_owned()
    }

    /// `parse_expect`: the next `str.len()` bytes must match, and are consumed.
    fn expect(&mut self, s: &[u8]) -> Result<(), SeqError> {
        if self.pos + s.len() > self.buf.len() || &self.buf[self.pos..self.pos + s.len()] != s {
            return Err(SeqError::OpenFormat);
        }
        self.pos += s.len();
        Ok(())
    }

    /// `parse_accept`: consume the bytes only if they all match.
    fn accept(&mut self, s: &[u8]) -> bool {
        if self.pos + s.len() > self.buf.len() || &self.buf[self.pos..self.pos + s.len()] != s {
            return false;
        }
        self.pos += s.len();
        true
    }

    /// `parse_peek`.
    fn peek(&self) -> Result<u8, SeqError> {
        self.buf.get(self.pos).copied().ok_or(SeqError::OpenFormat)
    }

    /// `parse_consume` of a single byte.
    fn consume(&mut self) -> Result<u8, SeqError> {
        let c = self.peek()?;
        self.pos += 1;
        Ok(c)
    }

    /// `parse_advance`.
    fn advance(&mut self, len: usize) -> Result<(), SeqError> {
        if self.pos + len > self.buf.len() {
            return Err(SeqError::OpenFormat);
        }
        self.pos += len;
        Ok(())
    }

    /// `parse_string`, a BER `VisibleString` (tag `0x1a`) with a definite length that is
    /// either short-form (< 128) or long-form (low 7 bits give the byte count):
    /// ```text
    ///   if (parse_expect(ncbi, "\x1a", 1) != eslOK)  return eslEFORMAT;
    ///   if (parse_consume(ncbi, &c, 1) != eslOK)     return eslEFORMAT;
    ///   if (c < 128) { n = c; }
    ///   else {
    ///     c = c & 0x7f;
    ///     if (c > sizeof(n)) return eslEFORMAT;
    ///     n = 0;
    ///     while (c > 0) { ...; n = (n << 8) + (unsigned int) x; --c; }
    ///   }
    ///   ptr = ncbi->hdr_ptr;
    ///   if (parse_advance(ncbi, n) != eslOK) return eslEFORMAT;
    /// ```
    fn string(&mut self) -> Result<Field, SeqError> {
        self.expect(b"\x1a")?;
        let c = self.consume()?;
        let n: usize = if c < 128 {
            c as usize
        } else {
            let mut c = c & 0x7f;
            // "if (c > sizeof(n))": C's int is 4 bytes here.
            if c as usize > std::mem::size_of::<i32>() {
                return Err(SeqError::OpenFormat);
            }
            let mut n: u32 = 0;
            while c > 0 {
                n = (n << 8) + self.consume()? as u32;
                c -= 1;
            }
            n as usize
        };
        let off = self.pos;
        self.advance(n)?;
        Ok(Field { off, len: n })
    }

    /// `parse_integer`, a BER `INTEGER` (tag `0x02`) read big-endian into a C `int`, so
    /// "if the integer is more bytes than the native int format, the most significant
    /// bytes will be lost":
    /// ```text
    ///   if (parse_expect(ncbi, "\x02", 1) != eslOK) return eslEFORMAT;
    ///   if (parse_peek(ncbi, &c) != eslOK)          return eslEFORMAT;
    ///   ptr = ncbi->hdr_ptr + 1;
    ///   if (parse_advance(ncbi, c + 1) != eslOK)    return eslEFORMAT;
    ///   n = 0;
    ///   while (c > 0) { n = (n << 8) + (unsigned int) *ptr++; --c; }
    /// ```
    /// The shift is on a signed `int` and wraps in practice; `wrapping_shl`/`wrapping_add`
    /// reproduce that without a debug-build panic.
    fn integer(&mut self) -> Result<i32, SeqError> {
        self.expect(b"\x02")?;
        let c = self.peek()?;
        let start = self.pos + 1;
        self.advance(c as usize + 1)?;
        let mut n: i32 = 0;
        for i in 0..c as usize {
            n = n.wrapping_shl(8).wrapping_add(self.buf[start + i] as i32);
        }
        Ok(n)
    }

    /// `ignore_sequence_of_integer`.
    fn ignore_integers(&mut self) -> Result<(), SeqError> {
        self.expect(b"\x30\x80")?;
        while self.peek()? == 0x02 {
            self.integer()?;
        }
        self.expect(b"\x00\x00")
    }

    /// `parse_textseq_id`:
    /// ```text
    /// Textseq-id ::= SEQUENCE {
    ///     name      VisibleString OPTIONAL ,
    ///     accession VisibleString OPTIONAL ,
    ///     release   VisibleString OPTIONAL ,
    ///     version   INTEGER       OPTIONAL
    /// }
    /// ```
    /// with C's precedence for which def line's identifiers win:
    /// ```text
    ///   if (acc != NULL && name != NULL) {
    ///     if (ncbi->name_ptr == NULL || ncbi->acc_ptr == NULL) { ...save both... }
    ///   } else if (ncbi->name_ptr == NULL && ncbi->acc_ptr == NULL) {
    ///     /* if neither the accession or name have been set, and the
    ///      * header supplied one, save it off. */
    ///     if (acc  != NULL) { ...save acc...  }
    ///     if (name != NULL) { ...save name... }
    ///   }
    /// ```
    fn textseq_id(&mut self) -> Result<(), SeqError> {
        self.expect(b"\x30\x80")?;
        let mut name = None;
        let mut acc = None;
        if self.accept(b"\xa0\x80") {
            name = Some(self.string()?);
            self.expect(b"\x00\x00")?;
        }
        if self.accept(b"\xa1\x80") {
            acc = Some(self.string()?);
            self.expect(b"\x00\x00")?;
        }
        if self.accept(b"\xa2\x80") {
            self.string()?;
            self.expect(b"\x00\x00")?;
        }
        if self.accept(b"\xa3\x80") {
            self.integer()?;
            self.expect(b"\x00\x00")?;
        }
        self.expect(b"\x00\x00")?;

        if acc.is_some() && name.is_some() {
            if self.name.is_none() || self.acc.is_none() {
                self.name = name;
                self.acc = acc;
            }
        } else if self.name.is_none() && self.acc.is_none() {
            if acc.is_some() {
                self.acc = acc;
            }
            if name.is_some() {
                self.name = name;
            }
        }
        Ok(())
    }

    /// `parse_object_id`:
    /// ```text
    /// Object-id ::= CHOICE { id INTEGER , str VisibleString }
    /// ```
    fn object_id(&mut self) -> Result<(), SeqError> {
        let mut id_str = None;
        let mut id = -1i32;
        if self.accept(b"\xa0\x80") {
            id = self.integer()?;
        } else if self.accept(b"\xa1\x80") {
            id_str = Some(self.string()?);
        } else {
            return Err(SeqError::OpenFormat);
        }
        self.expect(b"\x00\x00")?;
        if self.int_id == -1 && self.str_id.is_none() {
            if id_str.is_some() {
                self.str_id = id_str;
            } else if id != -1 {
                self.int_id = id;
            }
        }
        Ok(())
    }

    /// `parse_dbtag`:
    /// ```text
    /// Dbtag ::= SEQUENCE { db VisibleString , tag Object-id }
    /// ```
    /// "it looks like the dbtag is used when formatdb is run without parsing sequence ids
    /// (ie -o F). if that is the case, the id is equal to the sequence number in the
    /// database. so for dbtag headers, nothing will be saved. to do this lets create a
    /// bogus id value and restore it after dbtag is parsed."
    fn dbtag(&mut self) -> Result<(), SeqError> {
        self.expect(b"\x30\x80")?;
        self.expect(b"\xa0\x80")?;
        self.string()?;
        self.expect(b"\x00\x00")?;

        let temp_id = self.int_id;
        self.int_id = 1;
        self.expect(b"\xa1\x80")?;
        self.object_id()?;
        self.expect(b"\x00\x00")?;
        self.int_id = temp_id;

        self.expect(b"\x00\x00")
    }

    /// `parse_giimport_id`:
    /// ```text
    /// Giimport-id ::= SEQUENCE {
    ///     id INTEGER, db VisibleString OPTIONAL, release VisibleString OPTIONAL }
    /// ```
    fn giimport_id(&mut self) -> Result<(), SeqError> {
        self.expect(b"\x30\x80")?;
        self.expect(b"\xa0\x80")?;
        let id = self.integer()?;
        if self.accept(b"\xa1\x80") {
            self.string()?;
            self.expect(b"\x00\x00")?;
        }
        if self.accept(b"\xa2\x80") {
            self.string()?;
            self.expect(b"\x00\x00")?;
        }
        self.expect(b"\x00\x00")?;
        if self.int_id == -1 && self.str_id.is_none() {
            self.int_id = id;
        }
        Ok(())
    }

    /// `parse_id_pat`:
    /// ```text
    /// Id-pat ::= SEQUENCE {
    ///     country  VisibleString ,
    ///     id       CHOICE { number VisibleString , app-number VisibleString } ,
    ///     doc-type VisibleString OPTIONAL
    /// }
    /// ```
    fn id_pat(&mut self) -> Result<(), SeqError> {
        self.expect(b"\x30\x80")?;
        self.expect(b"\xa0\x80")?;
        self.string()?;
        self.expect(b"\xa1\x80")?;
        self.expect(b"\x30\x80")?;
        if self.accept(b"\xa0\x80") || self.accept(b"\xa1\x80") {
            self.string()?;
        } else {
            return Err(SeqError::OpenFormat);
        }
        self.expect(b"\x00\x00")?;
        if self.accept(b"\xa3\x80") {
            self.string()?;
        }
        self.expect(b"\x00\x00")
    }

    /// `parse_patent_seq_id`:
    /// ```text
    /// Patent-seq-id ::= SEQUENCE { seqid INTEGER , cit Id-pat }
    /// ```
    fn patent_seq_id(&mut self) -> Result<(), SeqError> {
        self.expect(b"\x30\x80")?;
        self.expect(b"\xa0\x80")?;
        let id = self.integer()?;
        self.expect(b"\xa1\x80")?;
        self.id_pat()?;
        self.expect(b"\x00\x00")?;
        if self.int_id == -1 && self.str_id.is_none() {
            self.int_id = id;
        }
        Ok(())
    }

    /// `parse_date_std`, every field of which is discarded.
    fn date_std(&mut self) -> Result<(), SeqError> {
        self.expect(b"\x30\x80")?;
        self.expect(b"\xa0\x80")?;
        self.integer()?;
        self.expect(b"\x00\x00")?;
        for tag in [
            &b"\xa1\x80"[..],
            &b"\xa2\x80"[..],
            &b"\xa3\x80"[..],
            &b"\xa4\x80"[..],
            &b"\xa5\x80"[..],
            &b"\xa6\x80"[..],
        ] {
            if self.accept(tag) {
                // The season (0xa3) is a string; the rest are integers.
                if tag == b"\xa3\x80" {
                    self.string()?;
                } else {
                    self.integer()?;
                }
                self.expect(b"\x00\x00")?;
            }
        }
        self.expect(b"\x00\x00")
    }

    /// `parse_pdb_seq_id`:
    /// ```text
    /// PDB-seq-id ::= SEQUENCE {
    ///     mol   PDB-mol-id , chain INTEGER , rel Date OPTIONAL }
    /// ```
    /// The molecule name is saved as the string id unconditionally -- unlike the other
    /// branches, C does not guard this one on the id being unset before overwriting the
    /// local `id`, only on the save itself.
    fn pdb_seq_id(&mut self) -> Result<(), SeqError> {
        self.expect(b"\x30\x80")?;
        self.expect(b"\xa0\x80")?;
        let id = self.string()?;
        self.expect(b"\x00\x00")?;
        if self.accept(b"\xa1\x80") {
            self.integer()?;
            self.expect(b"\x00\x00")?;
        }
        if self.accept(b"\xa2\x80") {
            if self.accept(b"\xa0\x80") {
                self.string()?;
            } else if self.accept(b"\xa1\x80") {
                self.date_std()?;
            } else {
                return Err(SeqError::OpenFormat);
            }
            self.expect(b"\x00\x00")?;
            self.expect(b"\x00\x00")?;
        }
        if self.int_id == -1 && self.str_id.is_none() {
            self.str_id = Some(id);
        }
        self.expect(b"\x00\x00")
    }

    /// `parse_seq_id`, the `Seq-id` CHOICE, dispatched on the context tag:
    /// ```text
    ///   while (c != 0x00) {
    ///     if (parse_expect(ncbi, "\x80", 1) != eslOK) return eslEFORMAT;
    ///     switch (c) {
    ///     case 0xa0: /* LOCAL     */ status = parse_object_id(ncbi);      break;
    ///     case 0xa1: /* GIBBSQ    */ id_ptr = (ncbi->int_id != -1) ? NULL : &ncbi->int_id;
    ///                                status = parse_integer(ncbi, id_ptr); break;
    ///     case 0xa2: /* GIBBMT    */ status = parse_integer(ncbi, NULL);  break;
    ///     case 0xa3: /* GIIM      */ status = parse_giimport_id(ncbi);    break;
    ///     case 0xa4: /* GENBANK   */
    ///     case 0xa5: /* EMBL      */
    ///     case 0xa6: /* PIR       */
    ///     case 0xa7: /* SWISSPROT */ status = parse_textseq_id(ncbi);     break;
    ///     case 0xa8: /* PATENT    */ status = parse_patent_seq_id(ncbi);  break;
    ///     case 0xa9: /* OTHER     */ status = parse_textseq_id(ncbi);     break;
    ///     case 0xaa: /* GENERAL   */ status = parse_dbtag(ncbi);          break;
    ///     case 0xab: /* GI        */ status = parse_integer(ncbi, NULL);  break;
    ///     case 0xac: /* DDBJ      */
    ///     case 0xad: /* PRF       */ status = parse_textseq_id(ncbi);     break;
    ///     case 0xae: /* PDB       */ status = parse_pdb_seq_id(ncbi);     break;
    ///     case 0xaf ... 0xb3:        status = parse_textseq_id(ncbi);     break;
    ///     default:                   status = eslEFORMAT;
    ///     }
    ///     ...
    /// ```
    fn seq_id(&mut self) -> Result<(), SeqError> {
        self.expect(b"\x30\x80")?;
        let mut c = self.consume()?;
        while c != 0x00 {
            self.expect(b"\x80")?;
            match c {
                0xa0 => self.object_id()?,
                0xa1 => {
                    // GIBBSQ: only the first integer id is kept.
                    let id = self.integer()?;
                    if self.int_id == -1 {
                        self.int_id = id;
                    }
                }
                0xa2 => {
                    self.integer()?;
                }
                0xa3 => self.giimport_id()?,
                0xa4..=0xa7 => self.textseq_id()?,
                0xa8 => self.patent_seq_id()?,
                0xa9 => self.textseq_id()?,
                0xaa => self.dbtag()?,
                0xab => {
                    self.integer()?;
                }
                0xac | 0xad => self.textseq_id()?,
                0xae => self.pdb_seq_id()?,
                0xaf..=0xb3 => self.textseq_id()?,
                _ => return Err(SeqError::OpenFormat),
            }
            self.expect(b"\x00\x00")?;
            c = self.consume()?;
        }
        self.expect(b"\x00")
    }

    /// `parse_def_line`:
    /// ```text
    /// Blast-def-line ::= SEQUENCE {
    ///     title       VisibleString       OPTIONAL,  -- simple title
    ///     seqid       SEQUENCE OF Seq-id,            -- Regular NCBI Seq-Id
    ///     taxid       INTEGER             OPTIONAL,  -- taxonomy id
    ///     memberships SEQUENCE OF INTEGER OPTIONAL,  -- bit arrays
    ///     links       SEQUENCE OF INTEGER OPTIONAL,  -- bit arrays
    ///     other-info  SEQUENCE OF INTEGER OPTIONAL   -- for future use
    /// }
    /// ```
    /// then the name/accession/description assignment, whose three-way fallback is what
    /// decides how a `>name description` FASTA header comes back out:
    /// ```text
    ///   if (ncbi->name_ptr != NULL || ncbi->acc_ptr != NULL) {
    ///     if (ncbi->name_ptr != NULL) {
    ///       esl_sq_SetName(sq, ncbi->name_ptr);
    ///       if (ncbi->acc_ptr != NULL) esl_sq_SetAccession(sq, ncbi->acc_ptr);
    ///     } else { esl_sq_SetName(sq, ncbi->acc_ptr); }
    ///     if (title != NULL) esl_sq_SetDesc(sq, title);
    ///   } else if (ncbi->str_id_ptr != NULL || ncbi->int_id != -1) {
    ///     if (ncbi->str_id_ptr != NULL) esl_sq_SetName(sq, ncbi->str_id_ptr);
    ///     else { char id[32]; snprintf(id, 32, "%d", ncbi->int_id); esl_sq_SetName(sq, id); }
    ///     if (title != NULL) esl_sq_SetDesc(sq, title);
    ///   } else if (title != NULL) {
    ///     /* take the first word of the title and use that for the name. the
    ///      * remaining portion of the title will be used for the description. */
    ///     ...
    ///   }
    /// ```
    /// `sq->desc[0] = 0` at the top of every def line, and the assignment runs for every
    /// def line in the header, so with several def lines the *last* one's title decides
    /// the description while the identifiers stay those of the first usable one.
    fn def_line(&mut self, name: &mut String, acc: &mut String, desc: &mut String) -> Result<(), SeqError> {
        self.expect(b"\x30\x80")?;

        desc.clear();
        let mut title = None;
        if self.accept(b"\xa0\x80") {
            title = Some(self.string()?);
            self.expect(b"\x00\x00")?;
        }

        self.expect(b"\xa1\x80")?;
        self.seq_id()?;
        self.expect(b"\x00\x00")?;

        if self.accept(b"\xa2\x80") {
            self.integer()?; // taxid, which hmmsearch never prints
            self.expect(b"\x00\x00")?;
        }
        for tag in [&b"\xa3\x80"[..], &b"\xa4\x80"[..], &b"\xa5\x80"[..]] {
            if self.accept(tag) {
                self.ignore_integers()?;
                self.expect(b"\x00\x00")?;
            }
        }
        self.expect(b"\x00\x00")?;

        if self.name.is_some() || self.acc.is_some() {
            if let Some(n) = self.name {
                *name = self.text(n);
                if let Some(a) = self.acc {
                    *acc = self.text(a);
                }
            } else {
                *name = self.text(self.acc.unwrap());
            }
            if let Some(t) = title {
                *desc = self.text(t);
            }
        } else if self.str_id.is_some() || self.int_id != -1 {
            *name = match self.str_id {
                Some(s) => self.text(s),
                None => format!("{}", self.int_id),
            };
            if let Some(t) = title {
                *desc = self.text(t);
            }
        } else if let Some(t) = title {
            let text = &self.buf[t.off..t.off + t.len];
            let mut i = text.len();
            for (j, &b) in text.iter().enumerate() {
                if is_space(b) {
                    i = j;
                    break;
                }
            }
            *name = String::from_utf8_lossy(&text[..i]).into_owned();
            i += 1;
            while i < text.len() && is_space(text[i]) {
                i += 1;
            }
            if i < text.len() {
                *desc = String::from_utf8_lossy(&text[i..]).into_owned();
            }
        }
        Ok(())
    }
}

fn is_space(c: u8) -> bool {
    matches!(c, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
}

/// `parse_header`, which walks the `Blast-def-line-set`:
/// ```text
///   reset_header_values(ncbi);
///   ...
///   /* verify we are at the beginning of a structure */
///   if (parse_expect(ncbi, "\x30\x80", 2) != eslOK)  return eslEFORMAT;
///   if (parse_peek(ncbi, &c) != eslOK)               return eslEFORMAT;
///   /* parse the different seq id structures */
///   while (c != 0x00) {
///     if ((status = parse_def_line(ncbi, sq)) != eslOK) return status;
///     if (parse_peek(ncbi, &c) != eslOK)             return eslEFORMAT;
///   }
///   /* verify we are at the end of the structure */
///   if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)  return eslEFORMAT;
/// ```
fn parse_header(buf: &[u8]) -> Result<(String, String, String), SeqError> {
    let mut p = HdrParse::new(buf);
    let (mut name, mut acc, mut desc) = (String::new(), String::new(), String::new());
    p.expect(b"\x30\x80")?;
    while p.peek()? != 0x00 {
        p.def_line(&mut name, &mut acc, &mut desc)?;
    }
    p.expect(b"\x00\x00")?;
    Ok((name, acc, desc))
}

/// `read_amino`, which is a straight table lookup over the record's bytes:
/// ```text
///   size = sq->eoff - sq->doff;
///   ESL_DSQ *ptr = sq->dsq + 1;
///   if (fread(ptr, sizeof(char), size, ncbi->fppsq) != size) return eslEFORMAT;
///   for (inx = 0; inx < size - 1; ++inx) { *ptr = sqfp->inmap[(int) *ptr]; ++ptr; }
///   *ptr = eslDSQ_SENTINEL;
///   ...
///   sq->n = size - 1;
/// ```
/// The record's trailing zero byte is not translated; it is overwritten by the dsq's end
/// sentinel, which is why the length is `size - 1`.
fn read_amino(seq: &[u8], doff: usize, eoff: usize, inmap: &[u8; 256]) -> Vec<u8> {
    let size = eoff - doff;
    let mut dsq = Vec::with_capacity(size + 1);
    dsq.push(255u8); // eslDSQ_SENTINEL
    for &b in &seq[doff..doff + size - 1] {
        dsq.push(inmap[b as usize]);
    }
    dsq.push(255u8);
    dsq
}

/// `read_dna`: unpack the 2-bit-per-base sequence, then overlay the ambiguity table.
///
/// ```text
///   size = sq->eoff - sq->doff;
///   ...
///   ssize     = ncbi->seq_apos - sq->doff - 1;
///   remainder = *(ncbi->hdr_buf + ssize) & 0x03;
///   length    = ssize * 4 + remainder;
///   ...
///   for (inx = 0; inx < ssize; ++inx) {
///     c = ncbi->hdr_buf[inx];
///     n = 1 << ((c >> 6) & 0x03);  *ptr = sqfp->inmap[n]; ++ptr;
///     n = 1 << ((c >> 4) & 0x03);  *ptr = sqfp->inmap[n]; ++ptr;
///     n = 1 << ((c >> 2) & 0x03);  *ptr = sqfp->inmap[n]; ++ptr;
///     n = 1 << ((c >> 0) & 0x03);  *ptr = sqfp->inmap[n]; ++ptr;
///   }
///   /* handle the remainder */
///   c = ncbi->hdr_buf[inx];
///   for (inx = 0; inx < remainder; ++inx) {
///     n = 1 << ((c >> (6 - inx * 2)) & 0x03);  *ptr = sqfp->inmap[n]; ++ptr;
///   }
///   *ptr = eslDSQ_SENTINEL;
/// ```
/// The last packed byte's low two bits hold the count of bases in it, so an exact
/// multiple of four ends with a byte of zero and `remainder == 0`.
///
/// Then the ambiguity table, whose entry width is announced by the top bit of its first
/// byte -- "we need to look that the first by the the ambiguity table to see if the
/// entries are 32 or 64 bit entries":
/// ```text
///   amb32 = 0;
///   if (ncbi->seq_apos - sq->doff < size) {
///     amb32 = ((ncbi->hdr_buf[ncbi->seq_apos - sq->doff] & 0x80) == 0);
///   }
///   /* skip past the count and start processing the abmiguity table */
///   ssize = ncbi->seq_apos - sq->doff + 4;
///   while (ssize < size) {
///     n = ((ncbi->hdr_buf[ssize] >> 4) & 0x0f);
///     c = sqfp->inmap[n];
///     if (amb32) {
///       cnt = (ncbi->hdr_buf[ssize] & 0x0f); cnt += 1;
///       off = ncbi->hdr_buf[ssize+1];
///       off = (off << 8) | ncbi->hdr_buf[ssize+2];
///       off = (off << 8) | ncbi->hdr_buf[ssize+3];
///       for (inx = 0; inx < cnt; ++inx) ptr[off+inx] = c;
///       ssize += 4;
///     } else {
///       cnt = (ncbi->hdr_buf[ssize] & 0x0f);
///       cnt = (cnt << 8) | ncbi->hdr_buf[ssize+1];  cnt += 1;
///       off = ncbi->hdr_buf[ssize+2];
///       ... six bytes ...
///       for (inx = 0; inx < cnt; ++inx) ptr[off+inx] = c;
///       ssize += 8;
///     }
///   }
/// ```
fn read_dna(
    seq: &[u8],
    doff: usize,
    eoff: usize,
    seq_apos: usize,
    inmap: &[u8; 256],
) -> Result<Vec<u8>, SeqError> {
    let bad = || SeqError::Parse(b"Error reading sequence index".to_vec());
    let size = eoff - doff;
    let rec = &seq[doff..eoff];
    if seq_apos <= doff || seq_apos - doff > size {
        return Err(bad());
    }
    let ssize = seq_apos - doff - 1;
    let remainder = (rec[ssize] & 0x03) as usize;
    let length = ssize * 4 + remainder;

    let mut dsq = Vec::with_capacity(length + 2);
    dsq.push(255u8);
    for &c in &rec[..ssize] {
        for shift in [6, 4, 2, 0] {
            dsq.push(inmap[1usize << ((c >> shift) & 0x03)]);
        }
    }
    let c = rec[ssize];
    for inx in 0..remainder {
        dsq.push(inmap[1usize << ((c >> (6 - inx * 2)) & 0x03)]);
    }
    dsq.push(255u8);

    // The ambiguity table overwrites already-unpacked residues in place, so it works on
    // the dsq's residue span (C's `ptr`, which is `sq->dsq + 1`).
    let mut amb32 = false;
    if seq_apos - doff < size {
        amb32 = (rec[seq_apos - doff] & 0x80) == 0;
    }
    let mut pos = seq_apos - doff + 4;
    while pos < size {
        let n = ((rec[pos] >> 4) & 0x0f) as usize;
        let c = inmap[n];
        let (cnt, off, step) = if amb32 {
            if pos + 4 > size {
                return Err(bad());
            }
            let cnt = (rec[pos] & 0x0f) as usize + 1;
            let off = ((rec[pos + 1] as usize) << 16)
                | ((rec[pos + 2] as usize) << 8)
                | rec[pos + 3] as usize;
            (cnt, off, 4)
        } else {
            if pos + 8 > size {
                return Err(bad());
            }
            let cnt = ((((rec[pos] & 0x0f) as usize) << 8) | rec[pos + 1] as usize) + 1;
            let mut off = 0usize;
            for i in 2..8 {
                off = (off << 8) | rec[pos + i] as usize;
            }
            (cnt, off, 8)
        };
        for inx in 0..cnt {
            // C writes through `ptr[off+inx]` with no bounds check; a corrupt table would
            // scribble past the sequence. Stopping at the end of the dsq is the one place
            // this reader is deliberately safer than its original.
            let at = 1 + off + inx;
            if at >= dsq.len() - 1 {
                break;
            }
            dsq[at] = c;
        }
        pos += step;
    }
    Ok(dsq)
}

impl NcbiDb {
    /// Read the whole database, as hmmsearch's search loop does through
    /// `sqncbi_ReadBlock` -- which, for anything that is not a long nucleotide target,
    /// "just read[s] in a sequence at a time" via `sqncbi_Read`.
    ///
    /// `sqncbi_Read` positions from the index tables, reads the sequence, then parses the
    /// header:
    /// ```text
    ///   if (ncbi->index >= ncbi->num_seq) return eslEOF;
    ///   if ((status = pos_sequence(ncbi, ncbi->index)) != eslOK) return status;
    ///   sq->idx  = ncbi->index;  sq->roff = ncbi->roff;  sq->doff = ncbi->doff;
    ///   sq->hoff = ncbi->hoff;   sq->eoff = ncbi->eoff;
    ///   if (ncbi->alphatype == eslAMINO) status = read_amino(sqfp, sq);
    ///   else                             status = read_dna(sqfp, sq);
    ///   if (status != eslOK) return status;
    ///   if ((status = parse_header(ncbi, sq)) != eslOK) return status;
    /// ```
    /// and `read_amino` takes the record's length straight from the index table, with the
    /// `.psq` record's trailing zero byte becoming the dsq's end sentinel:
    /// ```text
    ///   size = sq->eoff - sq->doff;
    ///   ESL_DSQ *ptr = sq->dsq + 1;
    ///   if (fread(ptr, sizeof(char), size, ncbi->fppsq) != size) return eslEFORMAT;
    ///   for (inx = 0; inx < size - 1; ++inx) { *ptr = sqfp->inmap[(int) *ptr]; ++ptr; }
    ///   *ptr = eslDSQ_SENTINEL;
    ///   ...
    ///   sq->n = size - 1;
    /// ```
    pub fn read(&self) -> Result<Vec<Seq>, SeqError> {
        let inmap = if self.is_dna { dna_inmap() } else { amino_inmap() };
        let mut out = Vec::new();
        for vol in &self.volumes {
            for inx in 0..vol.num_seq {
                let doff = vol.seq_indexes[inx] as usize;
                let eoff = vol.seq_indexes[inx + 1] as usize;
                let roff = vol.hdr_indexes[inx] as usize;
                let hoff = vol.hdr_indexes[inx + 1] as usize;
                if eoff <= doff || eoff > vol.seq.len() || hoff < roff || hoff > vol.hdr.len() {
                    return Err(SeqError::Parse(b"Error reading sequence index".to_vec()));
                }
                let dsq = if self.is_dna {
                    read_dna(&vol.seq, doff, eoff, vol.amb_indexes[inx] as usize, &inmap)?
                } else {
                    read_amino(&vol.seq, doff, eoff, &inmap)
                };
                let (name, acc, desc) = parse_header(&vol.hdr[roff..hoff])?;
                out.push(Seq { name, acc, desc, dsq, roff: 0, eoff: 0 });
            }
        }
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The NCBIstdaa-to-Easel translation, spot-checked against the two alphabet strings
    /// C intersects. The index is the *code* stored in a `.psq`, not a character.
    #[test]
    fn ncbi_amino_inmap_translates_stdaa_codes() {
        let m = amino_inmap();
        assert_eq!(m[0], alphabet::GAP); // '-'
        assert_eq!(m[1], 0); // 'A'
        assert_eq!(m[2], alphabet::AMINO_B); // 'B'
        assert_eq!(m[3], 1); // 'C'
        assert_eq!(m[21], alphabet::AMINO_X); // 'X'
        assert_eq!(m[22], 19); // 'Y'
        assert_eq!(m[23], alphabet::AMINO_Z); // 'Z'
        assert_eq!(m[24], alphabet::AMINO_U); // 'U'
        assert_eq!(m[25], alphabet::NONRESIDUE); // '*'
        assert_eq!(m[26], alphabet::AMINO_O); // 'O'
        assert_eq!(m[27], alphabet::AMINO_J); // 'J'
        // Codes past the alphabet stay eslDSQ_ILLEGAL, and read_amino copies them
        // through unchecked.
        assert_eq!(m[28], 254);
        assert_eq!(m[255], 254);
    }

    /// A minimal `Blast-def-line-set` with one title and one `local str` id, encoded by
    /// hand from the ASN.1 in esl_sqio_ncbi.c's comments.
    fn defline(title: Option<&[u8]>, local: &[u8]) -> Vec<u8> {
        let mut v = Vec::new();
        let string = |v: &mut Vec<u8>, s: &[u8]| {
            v.push(0x1a);
            v.push(s.len() as u8);
            v.extend_from_slice(s);
        };
        v.extend_from_slice(b"\x30\x80"); // Blast-def-line-set
        v.extend_from_slice(b"\x30\x80"); // Blast-def-line
        if let Some(t) = title {
            v.extend_from_slice(b"\xa0\x80");
            string(&mut v, t);
            v.extend_from_slice(b"\x00\x00");
        }
        v.extend_from_slice(b"\xa1\x80"); // seqid
        v.extend_from_slice(b"\x30\x80"); // SEQUENCE OF Seq-id
        v.extend_from_slice(b"\xa0\x80"); // local
        v.extend_from_slice(b"\xa1\x80"); // Object-id str
        string(&mut v, local);
        v.extend_from_slice(b"\x00\x00"); // end Object-id
        v.extend_from_slice(b"\x00\x00"); // end local
        // The SEQUENCE OF Seq-id terminator is two bytes: parse_seq_id consumes one into
        // its loop variable and then expects the other.
        v.extend_from_slice(b"\x00\x00");
        v.extend_from_slice(b"\x00\x00"); // end seqid
        v.extend_from_slice(b"\x00\x00"); // end Blast-def-line
        v.extend_from_slice(b"\x00\x00"); // end set
        v
    }

    #[test]
    fn header_uses_the_local_id_as_the_name() {
        let (name, acc, desc) = parse_header(&defline(Some(b"a description"), b"SEQ1")).unwrap();
        assert_eq!(name, "SEQ1");
        assert_eq!(acc, "");
        assert_eq!(desc, "a description");
    }

    /// With no identifier at all, C splits the title into a name and a description at the
    /// first run of whitespace.
    #[test]
    fn header_falls_back_to_splitting_the_title() {
        let mut v = Vec::new();
        v.extend_from_slice(b"\x30\x80\x30\x80\xa0\x80");
        v.push(0x1a);
        let t = b"NAME1   rest of  it";
        v.push(t.len() as u8);
        v.extend_from_slice(t);
        v.extend_from_slice(b"\x00\x00");
        // An empty but present seqid sequence, then the four two-byte terminators for
        // the SEQUENCE OF, the seqid, the def line and the set.
        v.extend_from_slice(b"\xa1\x80\x30\x80\x00\x00\x00\x00\x00\x00\x00\x00");
        let (name, acc, desc) = parse_header(&v).unwrap();
        assert_eq!(name, "NAME1");
        assert_eq!(acc, "");
        assert_eq!(desc, "rest of  it");
    }

    /// A long-form string length: the tag byte's high bit means "the low 7 bits are a
    /// byte count for the real length".
    #[test]
    fn long_form_string_lengths_parse() {
        let long = "x".repeat(300);
        let mut v = Vec::new();
        v.extend_from_slice(b"\x30\x80\x30\x80\xa0\x80");
        v.push(0x1a);
        v.push(0x82); // two length bytes
        v.extend_from_slice(&[0x01, 0x2c]); // 300
        v.extend_from_slice(long.as_bytes());
        v.extend_from_slice(b"\x00\x00");
        v.extend_from_slice(b"\xa1\x80\x30\x80\x00\x00\x00\x00\x00\x00\x00\x00");
        let (name, _, desc) = parse_header(&v).unwrap();
        assert_eq!(name, long);
        assert_eq!(desc, "");
    }

    /// A `.pin` whose version is not 4 is `eslEFORMAT`, and a protein read of a
    /// nucleotide index is `eslEUNIMPLEMENTED` (24) -- two different messages out of
    /// hmmsearch.
    #[test]
    fn index_version_and_type_are_checked() {
        let mut pin = Vec::new();
        pin.extend_from_slice(&5u32.to_be_bytes()); // version 5
        pin.extend_from_slice(&1u32.to_be_bytes());
        pin.extend_from_slice(&0u32.to_be_bytes());
        assert!(matches!(
            parse_index_file(&pin, AMINO_DB),
            Err(SeqError::OpenFormat)
        ));

        let mut pin = Vec::new();
        pin.extend_from_slice(&4u32.to_be_bytes());
        pin.extend_from_slice(&0u32.to_be_bytes()); // a nucleotide index
        pin.extend_from_slice(&0u32.to_be_bytes());
        assert!(matches!(
            parse_index_file(&pin, AMINO_DB),
            Err(SeqError::OpenStatus(24))
        ));
    }
}
