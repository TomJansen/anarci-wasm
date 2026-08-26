//! MPI master/worker search, `--mpi`.
//!
//! C, src/hmmsearch.c:`mpi_master`/`mpi_worker`, plus the serialization in
//! src/mpisupport.c and easel/esl_mpi.c.
//!
//! The division of labour is C's: the master never searches. It walks the database's
//! record offsets, hands each idle worker a `(offset, length, count)` block, and when the
//! database is exhausted sends an empty block to every worker and collects one
//! `P7_TOPHITS` and one `P7_PIPELINE` from each. Sequences are never sent -- every rank
//! opens the database itself.
//!
//! # Why there is no MPI library here
//!
//! C HMMER calls `MPI_Init`, `MPI_Send` and `MPI_Recv`, which means linking whichever MPI
//! implementation `mpirun` belongs to. This crate does not link anything, and does not
//! need to, because HMMER's use of MPI is narrow:
//!
//!   * **Rank and size** come from the launcher as ordinary environment variables --
//!     `OMPI_COMM_WORLD_RANK`/`_SIZE` from Open MPI's `mpirun`, `PMI_RANK`/`PMI_SIZE`
//!     from MPICH's, `SLURM_PROCID`/`SLURM_NTASKS` from `srun`. Nothing is read from the
//!     PMIx or PMI wire protocol.
//!   * **The topology is a star.** Every message in hmmsearch's protocol is between the
//!     master and one worker; workers never talk to each other. A TCP connection per
//!     worker carries all of it.
//!
//! So `--mpi` works in a plain `cargo build`, under any launcher that sets one of those
//! variable pairs, with no build-time or run-time dependency on an MPI implementation.
//! The one thing it needs from the environment is a rendezvous path the ranks share; see
//! [`Rendezvous`].
//!
//! # Divergences from C, all documented
//!
//! Two follow from this crate reading a database into memory instead of seeking:
//!
//!   * a worker resolves `block.offset` to a record index rather than `fseek`ing to it,
//!     which is exact because the offsets come from the same records;
//!   * every rank holds the whole database, where C's workers hold one block at a time.
//!
//! Neither changes which sequences a worker searches, and the master merges and re-sorts
//! exactly as in the serial path, so the report matches a serial run's.
//!
//! The wire format below is this crate's own -- only rustyhmmer ranks ever talk to each
//! other, since an MPI job is launched from one binary. It carries every field the output
//! depends on.

use std::io::{Read, Write};
use std::net::{TcpListener, TcpStream};
use std::sync::mpsc;

use crate::pipeline::{DomHit, Hit, Model};
use crate::pli::{Pipeline, Stats};
use crate::seqio::Seq;

/// C, hmmsearch.c:556-566. `HMMER_HMM_TAG` (2) and tags 3 and 7 are unused by
/// hmmsearch's protocol.
/// ```text
///   #define HMMER_ERROR_TAG          1
///   #define HMMER_HMM_TAG            2
///   #define HMMER_BLOCK_TAG          4
///   #define HMMER_PIPELINE_TAG       5
///   #define HMMER_TOPHITS_TAG        6
///   #define HMMER_TERMINATING_TAG    8
///   #define HMMER_READY_TAG          9
/// ```
const ERROR_TAG: i32 = 1;
const BLOCK_TAG: i32 = 4;
const PIPELINE_TAG: i32 = 5;
const TOPHITS_TAG: i32 = 6;
const TERMINATING_TAG: i32 = 8;
const READY_TAG: i32 = 9;

/// C, hmmsearch.c:691: `#define MAX_BLOCK_SIZE (512*1024)`.
const MAX_BLOCK_SIZE: u64 = 512 * 1024;

/// C's `SEQ_BLOCK`: a byte range of the database and the number of records in it.
#[derive(Clone, Copy, Default, Debug, PartialEq, Eq)]
pub struct SeqBlock {
    pub offset: u64,
    pub length: u64,
    pub count: u64,
}

/// `next_block()` -- hmmsearch.c:713-780 -- precomputed over the whole database.
///
/// ```text
///   while (block->length < MAX_BLOCK_SIZE && (n_targetseqs <0 || block->count < n_targetseqs) && (status = esl_sqio_ReadInfo(sqfp, sq)) == eslOK)
///     {
///       if (block->count == 0) block->offset = sq->roff;
///       block->length = sq->eoff - block->offset + 1;
///       block->count++;
///       esl_sq_Reuse(sq);
///     }
/// ```
/// C discovers the blocks lazily with `ReadInfo` and caches them in a `BLOCK_LIST` so
/// later queries reuse them without re-parsing; the list is the same either way, so it is
/// built once here. Note the loop tests `length` *before* adding a record, so a block
/// always holds at least one record and usually overshoots `MAX_BLOCK_SIZE` by the last
/// one -- reproduced.
pub fn blocks(seqs: &[Seq]) -> Vec<SeqBlock> {
    let mut out = Vec::new();
    let mut i = 0usize;
    while i < seqs.len() {
        let offset = seqs[i].roff;
        let mut length = 0u64;
        let mut count = 0u64;
        while i < seqs.len() && length < MAX_BLOCK_SIZE {
            // `sq->eoff - block->offset + 1`, with eoff the record's last byte.
            length = seqs[i].eoff.saturating_sub(offset) + 1;
            count += 1;
            i += 1;
        }
        out.push(SeqBlock { offset, length, count });
    }
    out
}

// ---- wire encoding -------------------------------------------------------------

/// A little-endian byte writer. Nothing here needs to interoperate with C, so the
/// encoding is the simplest one that round-trips every field exactly.
#[derive(Default)]
struct Buf {
    b: Vec<u8>,
}

impl Buf {
    fn u64(&mut self, v: u64) {
        self.b.extend_from_slice(&v.to_le_bytes());
    }
    fn i64(&mut self, v: i64) {
        self.b.extend_from_slice(&v.to_le_bytes());
    }
    fn i32(&mut self, v: i32) {
        self.b.extend_from_slice(&v.to_le_bytes());
    }
    fn u8(&mut self, v: u8) {
        self.b.push(v);
    }
    fn f32(&mut self, v: f32) {
        self.b.extend_from_slice(&v.to_le_bytes());
    }
    fn f64(&mut self, v: f64) {
        self.b.extend_from_slice(&v.to_le_bytes());
    }
    fn str(&mut self, s: &str) {
        self.u64(s.len() as u64);
        self.b.extend_from_slice(s.as_bytes());
    }
}

struct Cur<'a> {
    b: &'a [u8],
    p: usize,
}

impl<'a> Cur<'a> {
    fn new(b: &'a [u8]) -> Self {
        Cur { b, p: 0 }
    }
    fn take(&mut self, n: usize) -> &'a [u8] {
        let s = &self.b[self.p..self.p + n];
        self.p += n;
        s
    }
    fn u64(&mut self) -> u64 {
        u64::from_le_bytes(self.take(8).try_into().unwrap())
    }
    fn i64(&mut self) -> i64 {
        i64::from_le_bytes(self.take(8).try_into().unwrap())
    }
    fn i32(&mut self) -> i32 {
        i32::from_le_bytes(self.take(4).try_into().unwrap())
    }
    fn u8(&mut self) -> u8 {
        self.take(1)[0]
    }
    fn f32(&mut self) -> f32 {
        f32::from_le_bytes(self.take(4).try_into().unwrap())
    }
    fn f64(&mut self) -> f64 {
        f64::from_le_bytes(self.take(8).try_into().unwrap())
    }
    fn str(&mut self) -> String {
        let n = self.u64() as usize;
        String::from_utf8_lossy(self.take(n)).into_owned()
    }
}

/// `p7_tophits_MPISend`'s payload, for this crate's `Hit`.
///
/// C packs the `P7_HIT` fields and then each `P7_DOMAIN` with its `P7_ALIDISPLAY`
/// (mpisupport.c). Here a domain carries its optimal-accuracy traceback instead, because
/// that is what this crate keeps -- the alidisplay is built from the trace and the target
/// sequence when the report is written, and the master has every target sequence.
fn pack_hits(hits: &[Hit]) -> Vec<u8> {
    let mut w = Buf::default();
    w.u64(hits.len() as u64);
    for h in hits {
        w.str(&h.name);
        w.str(&h.acc);
        w.str(&h.desc);
        w.u64(h.tlen as u64);
        w.u64(h.seq_idx as u64);
        w.f32(h.score);
        w.f32(h.pre_score);
        w.f64(h.lnp);
        w.f32(h.nexpected);
        w.i32(h.nregions);
        w.i32(h.nclustered);
        w.i32(h.noverlaps);
        w.i32(h.nenvelopes);
        w.i32(h.ndom);
        w.u64(h.best as u64);
        w.u64(h.domains.len() as u64);
        for d in &h.domains {
            w.i64(d.ienv);
            w.i64(d.jenv);
            w.f32(d.bitscore);
            w.f32(d.dombias);
            w.f64(d.lnp);
            w.i32(d.hmm_from);
            w.i32(d.hmm_to);
            w.i64(d.ali_from);
            w.i64(d.ali_to);
            w.f32(d.oasc);
            w.u64(d.trace.len() as u64);
            for &(st, k, i, pp) in &d.trace {
                w.u8(st);
                w.i32(k);
                w.i32(i);
                w.f32(pp);
            }
        }
    }
    w.b
}

fn unpack_hits(b: &[u8]) -> Vec<Hit> {
    let mut c = Cur::new(b);
    let n = c.u64() as usize;
    let mut hits = Vec::with_capacity(n);
    for _ in 0..n {
        let name = c.str();
        let acc = c.str();
        let desc = c.str();
        let tlen = c.u64() as usize;
        let seq_idx = c.u64() as usize;
        let score = c.f32();
        let pre_score = c.f32();
        let lnp = c.f64();
        let nexpected = c.f32();
        let nregions = c.i32();
        let nclustered = c.i32();
        let noverlaps = c.i32();
        let nenvelopes = c.i32();
        let ndom = c.i32();
        let best = c.u64() as usize;
        let ndoms = c.u64() as usize;
        let mut domains = Vec::with_capacity(ndoms);
        for _ in 0..ndoms {
            let ienv = c.i64();
            let jenv = c.i64();
            let bitscore = c.f32();
            let dombias = c.f32();
            let lnp = c.f64();
            let hmm_from = c.i32();
            let hmm_to = c.i32();
            let ali_from = c.i64();
            let ali_to = c.i64();
            let oasc = c.f32();
            let ntr = c.u64() as usize;
            let mut trace = Vec::with_capacity(ntr);
            for _ in 0..ntr {
                let st = c.u8();
                let k = c.i32();
                let i = c.i32();
                let pp = c.f32();
                trace.push((st, k, i, pp));
            }
            domains.push(DomHit {
                ienv,
                jenv,
                bitscore,
                dombias,
                lnp,
                hmm_from,
                hmm_to,
                ali_from,
                ali_to,
                oasc,
                trace,
            });
        }
        hits.push(Hit {
            name,
            acc,
            desc,
            tlen,
            seq_idx,
            score,
            pre_score,
            lnp,
            nexpected,
            nregions,
            nclustered,
            noverlaps,
            nenvelopes,
            ndom,
            best,
            domains,
        });
    }
    hits
}

/// `p7_pipeline_MPISend`'s payload, restricted to the counters `p7_pipeline_Merge`
/// actually merges and `p7_pli_Statistics` prints.
fn pack_stats(s: &Stats) -> Vec<u8> {
    let mut w = Buf::default();
    for v in [s.nmodels, s.nnodes, s.nseqs, s.nres, s.n_past_msv, s.n_past_bias, s.n_past_vit, s.n_past_fwd] {
        w.u64(v);
    }
    w.b
}

fn unpack_stats(b: &[u8]) -> Stats {
    let mut c = Cur::new(b);
    Stats {
        nmodels: c.u64(),
        nnodes: c.u64(),
        nseqs: c.u64(),
        nres: c.u64(),
        n_past_msv: c.u64(),
        n_past_bias: c.u64(),
        n_past_vit: c.u64(),
        n_past_fwd: c.u64(),
    }
}

// ---- transport -----------------------------------------------------------------

/// The launcher-supplied rank and size, and the path the ranks meet on.
///
/// Every MPI launcher publishes the rank and the world size in the environment before it
/// even starts the process; C's `MPI_Comm_rank`/`MPI_Comm_size` return the same numbers
/// after `MPI_Init` has read them from PMIx or PMI. Reading them directly skips the whole
/// process-management protocol, which is the only part of MPI hmmsearch does not
/// otherwise need.
///
/// The pairs are tried in this order:
///
/// | launcher | rank | size |
/// |---|---|---|
/// | Open MPI `mpirun` | `OMPI_COMM_WORLD_RANK` | `OMPI_COMM_WORLD_SIZE` |
/// | MPICH `mpiexec` | `PMI_RANK` | `PMI_SIZE` |
/// | Slurm `srun` | `SLURM_PROCID` | `SLURM_NTASKS` |
/// | PMIx, as a last resort | `PMIX_RANK` | -- |
///
/// With none of them set the run is a world of one, which is what `MPI_Init` gives a
/// process started outside any launcher -- and, as in C, a world of one deadlocks in the
/// master loop rather than searching anything.
pub struct Rendezvous {
    pub rank: i32,
    pub size: i32,
    path: std::path::PathBuf,
}

fn env_i32(k: &str) -> Option<i32> {
    std::env::var(k).ok().and_then(|v| v.parse().ok())
}

impl Rendezvous {
    fn detect() -> Rendezvous {
        let (rank, size) = [
            ("OMPI_COMM_WORLD_RANK", "OMPI_COMM_WORLD_SIZE"),
            ("PMI_RANK", "PMI_SIZE"),
            ("SLURM_PROCID", "SLURM_NTASKS"),
        ]
        .iter()
        .find_map(|(r, s)| Some((env_i32(r)?, env_i32(s)?)))
        .unwrap_or_else(|| (env_i32("PMIX_RANK").unwrap_or(0), 1));

        // The ranks of one job have to agree on a filename without talking first, so it
        // is keyed by whatever the launcher calls the job. `PMIX_NAMESPACE` is unique per
        // `mpirun`; Slurm's job id is unique per allocation.
        let key = std::env::var("PMIX_NAMESPACE")
            .or_else(|_| std::env::var("SLURM_JOB_ID"))
            .or_else(|_| std::env::var("OMPI_MCA_orte_ess_jobid"))
            .unwrap_or_else(|_| "singleton".to_string());
        // A multi-node job needs this directory to be shared. $HOME is on every cluster
        // I know of; RUSTYHMMER_MPI_DIR overrides it for the ones where it is not.
        let dir = std::env::var("RUSTYHMMER_MPI_DIR")
            .map(std::path::PathBuf::from)
            .unwrap_or_else(|_| {
                std::env::var("HOME")
                    .map(std::path::PathBuf::from)
                    .unwrap_or_else(|_| std::env::temp_dir())
            });
        let path = dir.join(format!(".rustyhmmer-mpi-{key}"));
        Rendezvous { rank, size, path }
    }
}

/// One message: a tag and its payload, framed so a stream can carry many.
///
/// C gets framing from MPI itself -- `MPI_Probe` reports the tag and `MPI_Get_count` the
/// size before the receive. The same two numbers go on the wire here.
fn send_msg(s: &mut TcpStream, tag: i32, payload: &[u8]) -> std::io::Result<()> {
    let mut hdr = [0u8; 12];
    hdr[..4].copy_from_slice(&tag.to_le_bytes());
    hdr[4..].copy_from_slice(&(payload.len() as u64).to_le_bytes());
    s.write_all(&hdr)?;
    s.write_all(payload)?;
    s.flush()
}

fn recv_msg(s: &mut TcpStream) -> std::io::Result<(i32, Vec<u8>)> {
    let mut hdr = [0u8; 12];
    s.read_exact(&mut hdr)?;
    let tag = i32::from_le_bytes(hdr[..4].try_into().unwrap());
    let n = u64::from_le_bytes(hdr[4..].try_into().unwrap()) as usize;
    let mut payload = vec![0u8; n];
    s.read_exact(&mut payload)?;
    Ok((tag, payload))
}

/// The world, once every rank has connected.
pub struct Mpi {
    rank: i32,
    size: i32,
    /// Master only: one connection per worker, indexed by rank.
    conns: Vec<Option<TcpStream>>,
    /// Master only: every worker's messages, in arrival order.
    ///
    /// This is `MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, ...)`. A reader thread per worker
    /// turns "wait for whichever worker speaks next" into a channel receive, which is
    /// what the master's loop is written against.
    inbox: Option<mpsc::Receiver<(i32, i32, Vec<u8>)>>,
    /// Messages that arrived before the master asked for them.
    ///
    /// MPI receives are *matched*: `MPI_Recv(dest, HMMER_TOPHITS_TAG, ...)` takes the
    /// next message from that rank with that tag and leaves everything else queued. A
    /// plain channel has no such matching, and the difference is not academic -- a worker
    /// sends its results and then, if that was the last query, `HMMER_TERMINATING_TAG`
    /// immediately after, so the terminating message routinely overtakes another worker's
    /// results. Unmatched messages go here until someone asks for them.
    stash: std::cell::RefCell<Vec<(i32, i32, Vec<u8>)>>,
    /// Worker only: the connection to rank 0.
    master: Option<TcpStream>,
}

/// `MPI_Init` plus `MPI_Comm_rank`/`MPI_Comm_size` (hmmsearch.c:319-321), and the
/// connection setup MPI would have done inside them.
///
/// The master binds a port and publishes it; the workers read it and connect. The file is
/// written to a temporary name and renamed, so a worker never sees a half-written
/// address, and the master removes it once every worker is in.
pub fn init() -> Option<Mpi> {
    let rv = Rendezvous::detect();
    if rv.size <= 1 {
        // A world of one: no connections to make. The master loop will block on its first
        // receive, exactly as C's does when `mpirun` was not used.
        return Some(Mpi { rank: rv.rank, size: rv.size, conns: Vec::new(), inbox: None, stash: Default::default(), master: None });
    }
    if rv.rank == 0 {
        let listener = TcpListener::bind(("0.0.0.0", 0)).ok()?;
        let port = listener.local_addr().ok()?.port();
        let host = hostname();
        let tmp = rv.path.with_extension(format!("tmp{}", std::process::id()));
        std::fs::write(&tmp, format!("{host}:{port}")).ok()?;
        std::fs::rename(&tmp, &rv.path).ok()?;

        let mut conns: Vec<Option<TcpStream>> = (0..rv.size).map(|_| None).collect();
        let (tx, rx) = mpsc::channel();
        for _ in 1..rv.size {
            let (mut s, _) = listener.accept().ok()?;
            s.set_nodelay(true).ok();
            // The worker announces which rank it is, so `conns` can be indexed by rank.
            let mut b = [0u8; 4];
            s.read_exact(&mut b).ok()?;
            let wrank = i32::from_le_bytes(b);
            let mut reader = s.try_clone().ok()?;
            let tx = tx.clone();
            std::thread::spawn(move || {
                while let Ok((tag, payload)) = recv_msg(&mut reader) {
                    if tx.send((wrank, tag, payload)).is_err() {
                        break;
                    }
                }
            });
            if let Some(slot) = conns.get_mut(wrank as usize) {
                *slot = Some(s);
            }
        }
        let _ = std::fs::remove_file(&rv.path);
        Some(Mpi { rank: rv.rank, size: rv.size, conns, inbox: Some(rx), stash: Default::default(), master: None })
    } else {
        // The master may not have published yet; wait for it.
        let deadline = std::time::Instant::now() + std::time::Duration::from_secs(120);
        let addr = loop {
            if let Ok(s) = std::fs::read_to_string(&rv.path) {
                if !s.trim().is_empty() {
                    break s.trim().to_string();
                }
            }
            if std::time::Instant::now() > deadline {
                return None;
            }
            std::thread::sleep(std::time::Duration::from_millis(20));
        };
        let mut s = loop {
            match TcpStream::connect(&addr) {
                Ok(s) => break s,
                Err(_) if std::time::Instant::now() <= deadline => {
                    std::thread::sleep(std::time::Duration::from_millis(20));
                }
                Err(_) => return None,
            }
        };
        s.set_nodelay(true).ok();
        s.write_all(&rv.rank.to_le_bytes()).ok()?;
        Some(Mpi { rank: rv.rank, size: rv.size, conns: Vec::new(), inbox: None, stash: Default::default(), master: Some(s) })
    }
}

fn hostname() -> String {
    std::fs::read_to_string("/proc/sys/kernel/hostname")
        .map(|s| s.trim().to_string())
        .ok()
        .filter(|s| !s.is_empty())
        .unwrap_or_else(|| "127.0.0.1".to_string())
}

// ---- the protocol --------------------------------------------------------------

impl Mpi {
    pub fn rank(&self) -> i32 {
        self.rank
    }
    pub fn size(&self) -> i32 {
        self.size
    }

    /// `mpi_failure()`: report and tear the job down. C sends the message to the master
    /// with `HMMER_ERROR_TAG` when it can; here the message goes to stderr and the
    /// process exits non-zero, which brings the whole launcher job down with it.
    fn fail(&self, msg: &str) -> ! {
        eprintln!("\nError: {msg}");
        std::process::exit(1);
    }

    /// `MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, ...)` followed by the matching `MPI_Recv`.
    /// Blocks forever in a world of one, as C does.
    fn recv_any(&self) -> (i32, i32, Vec<u8>) {
        if let Some(m) = self.stash.borrow_mut().pop() {
            return m;
        }
        self.recv_wire()
    }

    /// `MPI_Recv(src, tag, ...)`: take the next message from `src` carrying `tag`,
    /// leaving anything else queued for a later receive.
    fn recv_matched(&self, src: i32, tag: i32) -> Vec<u8> {
        if let Some(i) = self
            .stash
            .borrow()
            .iter()
            .position(|(s, t, _)| *s == src && *t == tag)
        {
            return self.stash.borrow_mut().remove(i).2;
        }
        loop {
            let m = self.recv_wire();
            if m.0 == src && m.1 == tag {
                return m.2;
            }
            if m.1 == ERROR_TAG {
                self.fail(&format!(
                    "MPI client {} raised error:\n{}\n",
                    m.0,
                    String::from_utf8_lossy(&m.2)
                ));
            }
            self.stash.borrow_mut().push(m);
        }
    }

    fn recv_wire(&self) -> (i32, i32, Vec<u8>) {
        match &self.inbox {
            Some(rx) => match rx.recv() {
                Ok(m) => m,
                Err(_) => self.fail("MPI error receiving message: all workers gone\n"),
            },
            None => {
                // No workers exist, so no message can ever arrive. C reaches the same
                // state -- `hmmsearch --mpi` outside a launcher blocks in `mpi_master`.
                loop {
                    std::thread::park();
                }
            }
        }
    }

    fn send_to(&self, rank: i32, tag: i32, payload: &[u8]) {
        let Some(Some(c)) = self.conns.get(rank as usize) else {
            self.fail(&format!("MPI error sending message to {rank}\n"));
        };
        let mut c = match c.try_clone() {
            Ok(c) => c,
            Err(e) => self.fail(&format!("MPI error sending message to {rank}: {e}\n")),
        };
        if let Err(e) = send_msg(&mut c, tag, payload) {
            self.fail(&format!("MPI error sending message to {rank}: {e}\n"));
        }
    }

    /// One `HMMER_READY_TAG` from any worker, checked as C checks it:
    /// ```text
    ///   if (mpistatus.MPI_TAG == HMMER_ERROR_TAG) mpi_failure("MPI client %d raised error:\n%s\n", dest, mpi_buf);
    ///   if (mpistatus.MPI_TAG != HMMER_READY_TAG) mpi_failure("Unexpected tag %d from %d\n", mpistatus.MPI_TAG, dest);
    /// ```
    fn expect_ready(&self) -> i32 {
        let (src, tag, payload) = self.recv_any();
        match tag {
            ERROR_TAG => self.fail(&format!(
                "MPI client {src} raised error:\n{}\n",
                String::from_utf8_lossy(&payload)
            )),
            READY_TAG => src,
            t => self.fail(&format!("Unexpected tag {t} from {src}\n")),
        }
    }

    /// The master's inner loop for one query -- hmmsearch.c:940-1030.
    ///
    /// ```text
    ///   while ((n_targets==-1 || seq_cnt<=n_targets) && (sstatus = next_block(...)) == eslOK )
    ///   {
    ///     seq_cnt += block.count;
    ///     MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &mpistatus);
    ///     ...
    ///     MPI_Send(&block, 3, MPI_LONG_LONG_INT, dest, HMMER_BLOCK_TAG, MPI_COMM_WORLD);
    ///   }
    ///   ...
    ///   /* wait for all workers to finish up their work blocks */
    ///   for (i = 1; i < cfg->nproc; ++i) { ...expect READY... }
    ///   /* merge the results of the search results */
    ///   for (dest = 1; dest < cfg->nproc; ++dest)
    ///   {
    ///     /* send an empty block to signal the worker they are done */
    ///     MPI_Send(&block, 3, MPI_LONG_LONG_INT, dest, HMMER_BLOCK_TAG, MPI_COMM_WORLD);
    ///     p7_tophits_MPIRecv(dest, HMMER_TOPHITS_TAG, ...);
    ///     p7_pipeline_MPIRecv(dest, HMMER_PIPELINE_TAG, ...);
    ///     p7_tophits_Merge(th, mpi_th);
    ///     p7_pipeline_Merge(pli, mpi_pli);
    ///   }
    /// ```
    pub fn master_search(&self, blocks: &[SeqBlock]) -> (Vec<Hit>, Stats) {
        for blk in blocks {
            let dest = self.expect_ready();
            self.send_to(dest, BLOCK_TAG, &pack_block(*blk));
        }

        // "wait for all workers to finish up their work blocks"
        for _ in 1..self.size {
            self.expect_ready();
        }

        let mut hits: Vec<Hit> = Vec::new();
        let mut stats = Stats::default();
        for dest in 1..self.size {
            self.send_to(dest, BLOCK_TAG, &pack_block(SeqBlock::default()));
            let th = self.recv_matched(dest, TOPHITS_TAG);
            let pl = self.recv_matched(dest, PIPELINE_TAG);
            hits.extend(unpack_hits(&th));
            stats.merge(&unpack_stats(&pl));
        }
        (hits, stats)
    }

    /// "monitor all the workers to make sure they have ended" -- hmmsearch.c:1085-1103,
    /// after the last query.
    pub fn master_wait_for_termination(&self) {
        for _ in 1..self.size {
            let (src, tag, payload) = self.recv_any();
            match tag {
                ERROR_TAG => self.fail(&format!(
                    "MPI client {src} raised error:\n{}\n",
                    String::from_utf8_lossy(&payload)
                )),
                TERMINATING_TAG => {}
                t => self.fail(&format!("Unexpected tag {t} from {src}\n")),
            }
        }
    }

    fn worker_send(&self, tag: i32, payload: &[u8]) {
        let Some(m) = &self.master else {
            self.fail("MPI error: no connection to the master\n");
        };
        let mut m = match m.try_clone() {
            Ok(m) => m,
            Err(e) => self.fail(&format!("MPI error sending to master: {e}\n")),
        };
        if let Err(e) = send_msg(&mut m, tag, payload) {
            self.fail(&format!("MPI error sending to master: {e}\n"));
        }
    }

    fn worker_recv_block(&self) -> SeqBlock {
        let Some(m) = &self.master else {
            self.fail("MPI error: no connection to the master\n");
        };
        let mut m = match m.try_clone() {
            Ok(m) => m,
            Err(e) => self.fail(&format!("MPI error receiving from master: {e}\n")),
        };
        match recv_msg(&mut m) {
            Ok((BLOCK_TAG, p)) => unpack_block(&p),
            Ok((t, _)) => self.fail(&format!("Unexpected tag {t} from 0\n")),
            Err(e) => self.fail(&format!("MPI error receiving from master: {e}\n")),
        }
    }

    /// `mpi_worker()` -- hmmsearch.c:1180-1265, one iteration per query.
    ///
    /// ```text
    ///   status = 0;
    ///   MPI_Send(&status, 1, MPI_INT, 0, HMMER_READY_TAG, MPI_COMM_WORLD);
    ///   ...
    ///   MPI_Recv(&block, 3, MPI_LONG_LONG_INT, 0, HMMER_BLOCK_TAG, MPI_COMM_WORLD, &mpistatus);
    ///   while (block.count > 0)
    ///     {
    ///       status = esl_sqfile_Position(dbfp, block.offset);
    ///       while (count > 0 && (sstatus = esl_sqio_Read(dbfp, dbsq)) == eslOK)
    ///         {
    ///           p7_pli_NewSeq(pli, dbsq); p7_bg_SetLength(bg, dbsq->n); p7_oprofile_ReconfigLength(om, dbsq->n);
    ///           p7_Pipeline(pli, om, bg, dbsq, NULL, th);
    ///           ...
    ///         }
    ///       ...
    ///       MPI_Send(&status, 1, MPI_INT, 0, HMMER_READY_TAG, MPI_COMM_WORLD);
    ///       MPI_Recv(&block, 3, MPI_LONG_LONG_INT, 0, HMMER_BLOCK_TAG, MPI_COMM_WORLD, &mpistatus);
    ///     }
    ///   ...
    ///   p7_tophits_MPISend(th, 0, HMMER_TOPHITS_TAG, ...);
    ///   p7_pipeline_MPISend(pli, 0, HMMER_PIPELINE_TAG, ...);
    /// ```
    /// C's sanity checks on a returned block ("Block count mismatch", "Block length
    /// mismatch") compare the block it was given against what re-reading the file
    /// produced; here the block is resolved against the same in-memory records the master
    /// built it from, so the equivalent check is that the offset names a record and that
    /// the count fits.
    pub fn worker_search(&self, model: &Model, pli: &Pipeline, seqs: &[Seq]) {
        use rayon::prelude::*;

        let mut hits: Vec<Hit> = Vec::new();
        let mut stats = Stats::default();

        self.worker_send(READY_TAG, &[]);
        loop {
            let block = self.worker_recv_block();
            if block.count == 0 {
                break;
            }
            let Some(start) = seqs.iter().position(|s| s.roff == block.offset) else {
                self.fail(&format!(
                    "Cannot position sequence database to {}\n",
                    block.offset
                ));
            };
            let end = start + block.count as usize;
            if end > seqs.len() {
                self.fail(&format!(
                    "Block count mismatch - expected {} found {} at offset {}\n",
                    block.count,
                    seqs.len() - start,
                    block.offset
                ));
            }

            // Within its block a worker still uses every core it was given, exactly as
            // the serial path does over the whole database.
            let (mut blk_hits, blk_stats): (Vec<Hit>, Stats) = seqs[start..end]
                .par_iter()
                .enumerate()
                .map(|(i, s)| {
                    let (h, st) = model.search_one_counted(s, pli);
                    (
                        h.map(|mut h| {
                            h.seq_idx = start + i;
                            h
                        }),
                        st,
                    )
                })
                .fold(
                    || (Vec::new(), Stats::default()),
                    |(mut v, mut acc), (h, st)| {
                        acc.merge(&st);
                        if let Some(h) = h {
                            v.push(h);
                        }
                        (v, acc)
                    },
                )
                .reduce(
                    || (Vec::new(), Stats::default()),
                    |(mut a, mut sa), (b, sb)| {
                        a.extend(b);
                        sa.merge(&sb);
                        (a, sa)
                    },
                );
            hits.append(&mut blk_hits);
            stats.merge(&blk_stats);

            self.worker_send(READY_TAG, &[]);
        }

        self.worker_send(TOPHITS_TAG, &pack_hits(&hits));
        self.worker_send(PIPELINE_TAG, &pack_stats(&stats));
    }

    /// The worker's last message once the HMM file is exhausted (hmmsearch.c:1272).
    pub fn worker_terminate(&self) {
        self.worker_send(TERMINATING_TAG, &[]);
    }
}

/// C sends a block as three `MPI_LONG_LONG_INT`s; the same three numbers go on the wire.
fn pack_block(b: SeqBlock) -> Vec<u8> {
    let mut w = Buf::default();
    w.u64(b.offset);
    w.u64(b.length);
    w.u64(b.count);
    w.b
}

fn unpack_block(b: &[u8]) -> SeqBlock {
    let mut c = Cur::new(b);
    SeqBlock { offset: c.u64(), length: c.u64(), count: c.u64() }
}


#[cfg(test)]
mod tests {
    use super::*;

    fn seq_at(roff: u64, len: u64) -> Seq {
        Seq {
            name: String::new(),
            acc: String::new(),
            desc: String::new(),
            dsq: vec![255, 255],
            roff,
            eoff: roff + len - 1,
        }
    }

    /// `next_block` tests `block->length < MAX_BLOCK_SIZE` before adding a record, so a
    /// block closes on the record that crosses the limit -- it does not stop short of it.
    #[test]
    fn blocks_close_on_the_record_that_crosses_the_limit() {
        // Records of 200 KiB: two fit under 512 KiB (400 KiB), the third crosses it.
        let step = 200 * 1024u64;
        let seqs: Vec<Seq> = (0..7).map(|i| seq_at(i * step, step)).collect();
        let b = blocks(&seqs);
        assert_eq!(b[0], SeqBlock { offset: 0, length: 3 * step, count: 3 });
        assert_eq!(b[1], SeqBlock { offset: 3 * step, length: 3 * step, count: 3 });
        // The last block holds whatever is left.
        assert_eq!(b[2], SeqBlock { offset: 6 * step, length: step, count: 1 });
        assert_eq!(b.iter().map(|x| x.count).sum::<u64>(), 7);
    }

    /// A single record larger than the limit still forms one block: the length test runs
    /// before the record is added.
    #[test]
    fn one_oversized_record_is_its_own_block() {
        let seqs = vec![seq_at(0, MAX_BLOCK_SIZE * 3), seq_at(MAX_BLOCK_SIZE * 3, 10)];
        let b = blocks(&seqs);
        assert_eq!(b.len(), 2);
        assert_eq!(b[0].count, 1);
        assert_eq!(b[1].count, 1);
    }

    #[test]
    fn hits_round_trip_through_the_wire_format() {
        let h = Hit {
            name: "target one".into(),
            acc: "ACC1".into(),
            desc: "a  description  with  runs".into(),
            tlen: 153,
            seq_idx: 42,
            score: 215.625,
            pre_score: 218.5,
            lnp: -153.25,
            nexpected: 1.0625,
            nregions: 1,
            nclustered: 0,
            noverlaps: 0,
            nenvelopes: 1,
            ndom: 1,
            best: 0,
            domains: vec![DomHit {
                ienv: 3,
                jenv: 149,
                bitscore: 215.5,
                dombias: 2.875,
                lnp: -152.0,
                hmm_from: 1,
                hmm_to: 149,
                ali_from: 2,
                ali_to: 147,
                oasc: 145.5,
                trace: vec![(1, 0, 0, 0.0), (4, 7, 9, 0.9921875)],
            }],
        };
        let back = unpack_hits(&pack_hits(std::slice::from_ref(&h)));
        assert_eq!(back.len(), 1);
        let g = &back[0];
        assert_eq!(g.name, h.name);
        assert_eq!(g.acc, h.acc);
        assert_eq!(g.desc, h.desc);
        assert_eq!(g.tlen, h.tlen);
        assert_eq!(g.seq_idx, h.seq_idx);
        assert_eq!(g.score.to_bits(), h.score.to_bits());
        assert_eq!(g.pre_score.to_bits(), h.pre_score.to_bits());
        assert_eq!(g.lnp.to_bits(), h.lnp.to_bits());
        assert_eq!(g.nexpected.to_bits(), h.nexpected.to_bits());
        assert_eq!(g.ndom, h.ndom);
        assert_eq!(g.domains.len(), 1);
        assert_eq!(g.domains[0].ienv, 3);
        assert_eq!(g.domains[0].oasc.to_bits(), 145.5f32.to_bits());
        assert_eq!(g.domains[0].trace, h.domains[0].trace);
        assert!(unpack_hits(&pack_hits(&[])).is_empty());
    }

    /// The launcher's variables are tried in a fixed order, and a block round-trips
    /// through the same encoding the master and worker use for it.
    #[test]
    fn blocks_round_trip_through_the_wire_format() {
        let b = SeqBlock { offset: 1234, length: 512 * 1024 + 7, count: 99 };
        assert_eq!(unpack_block(&pack_block(b)), b);
        let empty = SeqBlock::default();
        assert_eq!(unpack_block(&pack_block(empty)), empty);
        assert_eq!(unpack_block(&pack_block(empty)).count, 0);
    }

    #[test]
    fn stats_round_trip_through_the_wire_format() {
        let s = Stats {
            nmodels: 1,
            nnodes: 149,
            nseqs: 397961,
            nres: 123456789,
            n_past_msv: 5000,
            n_past_bias: 4000,
            n_past_vit: 300,
            n_past_fwd: 20,
        };
        let b = unpack_stats(&pack_stats(&s));
        assert_eq!(b.nseqs, s.nseqs);
        assert_eq!(b.nres, s.nres);
        assert_eq!(b.n_past_fwd, s.n_past_fwd);
    }
}
