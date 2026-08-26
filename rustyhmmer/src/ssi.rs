//! SSI sequence-index input, for `--restrictdb_stkey` and `--ssifile`.
//!
//! C, easel/esl_ssi.c. Only the lookup half is transcribed: hmmsearch reads an index that
//! `esl-sfetch --index` wrote, and never creates one.
//!
//! Every field is stored in network byte order. `esl_ssi_Open` accepts a byteswapped
//! magic number as well as the native one but then reads the rest big-endian regardless,
//! so a byteswapped index would parse to garbage; that is reproduced by accepting both
//! magics and reading big-endian.

/// `esl_ssi.c:19-20`:
/// ```text
///   static uint32_t v30magic = 0xd3d3c9b3; /* SSI 3.0: "ssi3" + 0x80808080 */
///   static uint32_t v30swap  = 0xb3c9d3d3; /* byteswapped */
/// ```
const V30MAGIC: u32 = 0xd3d3_c9b3;
const V30SWAP: u32 = 0xb3c9_d3d3;

/// Why an index could not be opened, mapped to the Easel status codes hmmsearch
/// distinguishes.
#[derive(Debug, PartialEq, Eq)]
pub enum SsiError {
    /// `eslENOTFOUND`: no such file.
    NotFound,
    /// `eslEFORMAT`: not an SSI file, or truncated.
    Format,
    /// `eslERANGE`: an offset size this build cannot represent.
    Range,
}

/// An open SSI index, held in memory. C seeks within the open file for each lookup; the
/// indexes hmmsearch reads are small enough that reading them whole is simpler and gives
/// the same answers.
#[derive(Debug)]
pub struct Ssi {
    buf: Vec<u8>,
    offsz: usize,
    nprimary: u64,
    nsecondary: u64,
    plen: usize,
    slen: usize,
    precsize: usize,
    srecsize: usize,
    poffset: u64,
    soffset: u64,
}

fn be32(b: &[u8], o: usize) -> Option<u32> {
    b.get(o..o + 4)
        .map(|c| u32::from_be_bytes([c[0], c[1], c[2], c[3]]))
}
fn be16(b: &[u8], o: usize) -> Option<u16> {
    b.get(o..o + 2).map(|c| u16::from_be_bytes([c[0], c[1]]))
}
fn be64(b: &[u8], o: usize) -> Option<u64> {
    b.get(o..o + 8).map(|c| {
        u64::from_be_bytes([c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7]])
    })
}

/// `esl_fread_offset`, which reads either 4 or 8 bytes depending on the index's
/// `offsz`.
fn offset(b: &[u8], o: usize, offsz: usize) -> Option<u64> {
    match offsz {
        8 => be64(b, o),
        4 => be32(b, o).map(u64::from),
        _ => None,
    }
}

impl Ssi {
    /// `esl_ssi_Open`:
    /// ```text
    ///   status = eslENOTFOUND;
    ///   if ((ssi->fp = fopen(filename, "rb")) == NULL) goto ERROR;
    ///   status = eslEFORMAT;
    ///   if (esl_fread_u32(ssi->fp, &magic)        != eslOK) goto ERROR;
    ///   if (magic != v30magic && magic != v30swap)          goto ERROR;
    ///   if (esl_fread_u32(ssi->fp, &(ssi->flags)) != eslOK) goto ERROR;
    ///   if (esl_fread_u32(ssi->fp, &(ssi->offsz)) != eslOK) goto ERROR;
    ///   status = eslERANGE;
    ///   if (ssi->offsz != 4 && ssi->offsz != 8) goto ERROR;
    ///   if (ssi->offsz > sizeof(off_t))         goto ERROR;
    ///   status = eslEFORMAT;
    ///   if (esl_fread_u16(ssi->fp, &(ssi->nfiles))     != eslOK) goto ERROR;
    ///   if (esl_fread_u64(ssi->fp, &(ssi->nprimary))   != eslOK) goto ERROR;
    ///   if (esl_fread_u64(ssi->fp, &(ssi->nsecondary)) != eslOK) goto ERROR;
    ///   if (esl_fread_u32(ssi->fp, &(ssi->flen))       != eslOK) goto ERROR;
    ///   if (esl_fread_u32(ssi->fp, &(ssi->plen))       != eslOK) goto ERROR;
    ///   if (esl_fread_u32(ssi->fp, &(ssi->slen))       != eslOK) goto ERROR;
    ///   if (esl_fread_u32(ssi->fp, &(ssi->frecsize))   != eslOK) goto ERROR;
    ///   if (esl_fread_u32(ssi->fp, &(ssi->precsize))   != eslOK) goto ERROR;
    ///   if (esl_fread_u32(ssi->fp, &(ssi->srecsize))   != eslOK) goto ERROR;
    ///   if (esl_fread_offset(ssi->fp, ssi->offsz, &(ssi->foffset)) != eslOK) goto ERROR;
    ///   if (esl_fread_offset(ssi->fp, ssi->offsz, &(ssi->poffset)) != eslOK) goto ERROR;
    ///   if (esl_fread_offset(ssi->fp, ssi->offsz, &(ssi->soffset)) != eslOK) goto ERROR;
    ///   ...
    ///   if (ssi->nfiles == 0) goto ERROR;
    /// ```
    /// The per-file records are read too, but only `esl_ssi_FileInfo` uses them and
    /// hmmsearch never calls it.
    pub fn open(path: &str) -> Result<Ssi, SsiError> {
        let buf = std::fs::read(path).map_err(|_| SsiError::NotFound)?;
        let magic = be32(&buf, 0).ok_or(SsiError::Format)?;
        if magic != V30MAGIC && magic != V30SWAP {
            return Err(SsiError::Format);
        }
        let _flags = be32(&buf, 4).ok_or(SsiError::Format)?;
        let offsz = be32(&buf, 8).ok_or(SsiError::Format)? as usize;
        if offsz != 4 && offsz != 8 {
            return Err(SsiError::Range);
        }
        let nfiles = be16(&buf, 12).ok_or(SsiError::Format)?;
        let nprimary = be64(&buf, 14).ok_or(SsiError::Format)?;
        let nsecondary = be64(&buf, 22).ok_or(SsiError::Format)?;
        let _flen = be32(&buf, 30).ok_or(SsiError::Format)? as usize;
        let plen = be32(&buf, 34).ok_or(SsiError::Format)? as usize;
        let slen = be32(&buf, 38).ok_or(SsiError::Format)? as usize;
        let _frecsize = be32(&buf, 42).ok_or(SsiError::Format)? as usize;
        let precsize = be32(&buf, 46).ok_or(SsiError::Format)? as usize;
        let srecsize = be32(&buf, 50).ok_or(SsiError::Format)? as usize;
        let mut o = 54;
        let _foffset = offset(&buf, o, offsz).ok_or(SsiError::Format)?;
        o += offsz;
        let poffset = offset(&buf, o, offsz).ok_or(SsiError::Format)?;
        o += offsz;
        let soffset = offset(&buf, o, offsz).ok_or(SsiError::Format)?;
        if nfiles == 0 {
            return Err(SsiError::Format);
        }
        Ok(Ssi {
            buf,
            offsz,
            nprimary,
            nsecondary,
            plen,
            slen,
            precsize,
            srecsize,
            poffset,
            soffset,
        })
    }

    /// `binary_search`, which compares fixed-width NUL-padded key fields with `strcmp`:
    /// ```text
    ///   if (maxidx == 0) return eslENOTFOUND; /* special case: empty index */
    ///   left  = 0;
    ///   right = maxidx-1;
    ///   while (1) {
    ///     mid   = (left+right) / 2;
    ///     if (fseeko(ssi->fp, base + recsize*mid, SEEK_SET) != 0)    goto ERROR;
    ///     if (fread(name, sizeof(char), klen, ssi->fp)      != klen) goto ERROR;
    ///     cmp = strcmp(name, key);
    ///     if      (cmp == 0) break;
    ///     else if (left >= right) goto ERROR;
    ///     else if (cmp < 0)       left  = mid+1;
    ///     else if (cmp > 0) {
    ///       if (mid == 0) goto ERROR;
    ///       else right = mid-1;
    ///     }
    ///   }
    /// ```
    /// Returns the offset just past the matched key field, which is where C leaves the
    /// stream for the caller to read the record's remaining fields.
    fn binary_search(&self, key: &str, klen: usize, base: u64, recsize: usize, maxidx: u64) -> Option<usize> {
        if maxidx == 0 {
            return None;
        }
        let mut left: u64 = 0;
        let mut right: u64 = maxidx - 1;
        loop {
            let mid = (left + right) / 2;
            let at = (base + recsize as u64 * mid) as usize;
            let field = self.buf.get(at..at + klen)?;
            // strcmp() stops at the first NUL of the padded field.
            let name = match field.iter().position(|&b| b == 0) {
                Some(i) => &field[..i],
                None => field,
            };
            match name.cmp(key.as_bytes()) {
                std::cmp::Ordering::Equal => return Some(at + klen),
                _ if left >= right => return None,
                std::cmp::Ordering::Less => left = mid + 1,
                std::cmp::Ordering::Greater => {
                    if mid == 0 {
                        return None;
                    }
                    right = mid - 1;
                }
            }
        }
    }

    /// `esl_ssi_FindName`: look `key` up as a primary key, and failing that as a
    /// secondary key, which maps to a primary key that is then looked up again.
    ///
    /// ```text
    ///   status = binary_search(ssi, key, ssi->plen, ssi->poffset, ssi->precsize, ssi->nprimary);
    ///   if (status == eslOK)
    ///     { /* We found it as a primary key; get our data & return. */
    ///       if (esl_fread_u16(ssi->fp, ret_fh)                  != eslOK) goto ERROR;
    ///       if (esl_fread_offset(ssi->fp, ssi->offsz, ret_roff) != eslOK) goto ERROR;
    ///       ...
    ///     }
    ///   else if (status == eslENOTFOUND)
    ///     { /* Not in the primary keys? OK, try the secondary keys. */
    ///       if (ssi->nsecondary > 0) {
    ///         if ((status = binary_search(ssi, key, ssi->slen, ssi->soffset, ssi->srecsize, ssi->nsecondary)) != eslOK) goto ERROR;
    ///         /* We have the secondary key; flip to its primary key, then look that up. */
    ///         if (fread(pkey, sizeof(char), ssi->plen, ssi->fp) != ssi->plen) goto ERROR;
    ///         if ((status = esl_ssi_FindName(ssi, pkey, ret_fh, ret_roff, &doff, &L)) != eslOK) goto ERROR;
    ///       } else goto ERROR;
    ///     } else goto ERROR;
    /// ```
    /// Returns the record offset in the indexed file.
    pub fn find_name(&self, key: &str) -> Option<u64> {
        if let Some(at) = self.binary_search(key, self.plen, self.poffset, self.precsize, self.nprimary) {
            // fh (u16), then roff.
            return offset(&self.buf, at + 2, self.offsz);
        }
        if self.nsecondary == 0 {
            return None;
        }
        let at = self.binary_search(key, self.slen, self.soffset, self.srecsize, self.nsecondary)?;
        let field = self.buf.get(at..at + self.plen)?;
        let pkey = match field.iter().position(|&b| b == 0) {
            Some(i) => &field[..i],
            None => field,
        };
        let pkey = std::str::from_utf8(pkey).ok()?;
        self.find_name(pkey)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build an index by hand from the layout `esl_ssi_Open` reads, with 32-bit offsets,
    /// one file, `nprimary` primary keys and `nsecondary` secondary keys.
    fn make_ssi(primary: &[(&str, u32)], secondary: &[(&str, &str)]) -> Vec<u8> {
        let flen = 16usize;
        let plen = 16usize;
        let slen = 16usize;
        let offsz = 4usize;
        let frecsize = flen + 16;
        let precsize = plen + 2 + offsz + offsz + 8;
        let srecsize = slen + plen;
        let header = 54 + 3 * offsz;
        let foffset = header;
        let poffset = foffset + frecsize;
        let soffset = poffset + precsize * primary.len();

        let mut v = Vec::new();
        v.extend_from_slice(&V30MAGIC.to_be_bytes());
        v.extend_from_slice(&0u32.to_be_bytes()); // flags
        v.extend_from_slice(&(offsz as u32).to_be_bytes());
        v.extend_from_slice(&1u16.to_be_bytes()); // nfiles
        v.extend_from_slice(&(primary.len() as u64).to_be_bytes());
        v.extend_from_slice(&(secondary.len() as u64).to_be_bytes());
        v.extend_from_slice(&(flen as u32).to_be_bytes());
        v.extend_from_slice(&(plen as u32).to_be_bytes());
        v.extend_from_slice(&(slen as u32).to_be_bytes());
        v.extend_from_slice(&(frecsize as u32).to_be_bytes());
        v.extend_from_slice(&(precsize as u32).to_be_bytes());
        v.extend_from_slice(&(srecsize as u32).to_be_bytes());
        v.extend_from_slice(&(foffset as u32).to_be_bytes());
        v.extend_from_slice(&(poffset as u32).to_be_bytes());
        v.extend_from_slice(&(soffset as u32).to_be_bytes());
        assert_eq!(v.len(), header);

        let padded = |s: &str, n: usize| {
            let mut b = s.as_bytes().to_vec();
            b.resize(n, 0);
            b
        };
        // one file record
        v.extend_from_slice(&padded("seq.fa", flen));
        for _ in 0..4 {
            v.extend_from_slice(&0u32.to_be_bytes());
        }
        // primary key records, which the index requires to be sorted
        for (k, roff) in primary {
            v.extend_from_slice(&padded(k, plen));
            v.extend_from_slice(&0u16.to_be_bytes()); // fh
            v.extend_from_slice(&roff.to_be_bytes()); // roff
            v.extend_from_slice(&0u32.to_be_bytes()); // doff
            v.extend_from_slice(&0i64.to_be_bytes()); // L
        }
        for (s, p) in secondary {
            v.extend_from_slice(&padded(s, slen));
            v.extend_from_slice(&padded(p, plen));
        }
        v
    }

    fn write_tmp(name: &str, bytes: &[u8]) -> String {
        let p = std::env::temp_dir().join(format!("rustyhmmer-ssi-test-{name}"));
        std::fs::write(&p, bytes).unwrap();
        p.to_string_lossy().into_owned()
    }

    #[test]
    fn finds_primary_keys_by_binary_search() {
        // Sorted by strcmp, as esl-sfetch writes them.
        let keys: Vec<(String, u32)> = (0..64)
            .map(|i| (format!("k{i:04}"), (i as u32 + 1) * 100))
            .collect();
        let refs: Vec<(&str, u32)> = keys.iter().map(|(k, o)| (k.as_str(), *o)).collect();
        let path = write_tmp("primary", &make_ssi(&refs, &[]));
        let ssi = Ssi::open(&path).unwrap();
        for (k, o) in &refs {
            assert_eq!(ssi.find_name(k), Some(u64::from(*o)), "key {k}");
        }
        assert_eq!(ssi.find_name("k9999"), None);
        assert_eq!(ssi.find_name("a"), None);
        let _ = std::fs::remove_file(&path);
    }

    /// `esl_ssi_FindName` falls back to the secondary keys and then looks the primary key
    /// they name up again -- which is how an accession resolves to a record.
    #[test]
    fn secondary_keys_resolve_through_their_primary() {
        let path = write_tmp(
            "secondary",
            &make_ssi(&[("aaa", 10), ("bbb", 20)], &[("ACC_A", "aaa"), ("ACC_B", "bbb")]),
        );
        let ssi = Ssi::open(&path).unwrap();
        assert_eq!(ssi.find_name("bbb"), Some(20));
        assert_eq!(ssi.find_name("ACC_A"), Some(10));
        assert_eq!(ssi.find_name("ACC_B"), Some(20));
        assert_eq!(ssi.find_name("ACC_C"), None);
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn rejects_a_file_that_is_not_an_index() {
        let path = write_tmp("bogus", b"not an ssi file at all");
        assert!(matches!(Ssi::open(&path), Err(SsiError::Format)));
        let _ = std::fs::remove_file(&path);
        assert!(matches!(Ssi::open("/nonexistent/x.ssi"), Err(SsiError::NotFound)));
    }

    /// `if (ssi->offsz != 4 && ssi->offsz != 8) goto ERROR;` with `status = eslERANGE`.
    #[test]
    fn rejects_an_unsupported_offset_size() {
        let mut bytes = make_ssi(&[("a", 1)], &[]);
        bytes[8..12].copy_from_slice(&2u32.to_be_bytes());
        let path = write_tmp("offsz", &bytes);
        assert!(matches!(Ssi::open(&path), Err(SsiError::Range)));
        let _ = std::fs::remove_file(&path);
    }
}
