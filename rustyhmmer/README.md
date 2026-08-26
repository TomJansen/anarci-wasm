# rustyhmmer

Pure-Rust reimplementation of [HMMER3](http://hmmer.org/)'s `hmmsearch`,
producing **byte-identical output** to HMMER 3.4 for protein profile-HMM
search — the human-readable report including the alignment display, and every
output file: `--tblout`, `--domtblout`, `--pfamtblout` and `-A`.

No dependency but [rayon](https://crates.io/crates/rayon): no C library, no
`*-sys` crate, no build script. That holds for `--mpi` too (see below).

HMMER is pronounced "hammer"; this is a rusty one.

## Status

- Amino-acid `hmmsearch` pipeline: MSV/SSV (F1) → bias (F3b) → Viterbi (F2)
  → Forward (F3) → domain definition → optimal-accuracy alignment →
  per-sequence / per-domain scoring → reporting.
- Every command-line option of C `hmmsearch` is implemented.
- Output is validated byte-for-byte against C HMMER 3.4 on 117 bac120 profiles
  × 397,961 proteins, comparing the report, `--tblout`, `--domtblout`, `-A` and
  `--pfamtblout` together.
- SIMD through auto-vectorization plus an SSE2 path, so it runs on both x86-64
  and AArch64.

## Build

```sh
cargo build --release
```

## Usage

```sh
rustyhmmer-hmmsearch [options] <query.hmm> <target.faa>
```

The query may be an ASCII HMMER3 `.hmm` file or an `hmmpress`ed binary `.h3m`;
either argument may be `-` to read from standard input. The target database is
FASTA by default, and `--tformat` accepts every format C accepts, including
pressed NCBI BLAST databases:

```sh
rustyhmmer-hmmsearch --tformat ncbi      query.hmm blastdb_basename
rustyhmmer-hmmsearch --tformat stockholm query.hmm alignment.sto
```

### Running across machines with `--mpi`

`--mpi` works under `mpirun` or `srun` exactly as C HMMER's does, but this
binary links no MPI library — check with `ldd`. It needs only two things from
the launcher, and every launcher provides them: the rank and world size, which
arrive as ordinary environment variables (`OMPI_COMM_WORLD_RANK`/`_SIZE`,
`PMI_RANK`/`PMI_SIZE`, or `SLURM_PROCID`/`SLURM_NTASKS`). HMMER's protocol is a
star — the master hands blocks to workers and merges what they return, and
workers never talk to each other — so one TCP connection per worker carries it.

```sh
mpirun -np 8 rustyhmmer-hmmsearch --mpi query.hmm target.faa
```

The ranks find each other through a file the master writes and then removes,
by default under `$HOME`. A job spread over several nodes needs that directory
to be shared; set `RUSTYHMMER_MPI_DIR` if `$HOME` is not.

### Per-domain profile coordinates from the library API

`--domtblout`'s `hmm from`/`hmm to` columns are also available in-process, per
domain, so consumers can compute HMM coverage without shelling out:

```rust
use rustyhmmer::api::{Cutoff, HmmAnnotator};

let ann = HmmAnnotator::from_hmm_file("bac120.hmm")?
    .with_cutoff(Cutoff::Evalue(1e-10));
for hit in ann.search(&seqs) {
    for d in &hit.domains {
        let coverage =
            (d.hmm_to - d.hmm_from + 1) as f64 / hit.model_len as f64;
        // d.ali_from / d.ali_to / d.env_from / d.env_to / d.acc are also here
    }
}
```

## Scope and limitations

- Protein (amino-acid) models only.
- `hmmsearch` only (no `hmmscan` / `nhmmer` / `phmmer` / `jackhmmer`).
- Local, multihit mode (HMMER's default).
- Byte-parity is verified for both E-value (`-E`) and model-specific
  (`--cut_ga` / `--cut_tc` / `--cut_nc`) reporting thresholds.
- `.gz` databases are not read. C pipes them through `gzip -dc`; this does not.
- `--seed 0` asks for an arbitrary seed, so its results are not reproducible —
  in C either.

## Changelog

### 0.1.4

Every `hmmsearch` command-line option is now implemented, and parity is
verified for every output file at once rather than the tabular ones alone.

- **`--mpi`, with no MPI library linked.** The master/worker protocol of
  `mpi_master`/`mpi_worker` is transcribed: blocks of database offsets go out to
  idle workers, hits and pipeline counters come back and are merged. Rank and
  world size are read from the launcher's environment variables and the star
  topology is carried over TCP, so `cargo build` gains no dependency and the
  binary links no MPI. Verified four ways at once — C serial, C under `mpirun`,
  rustyhmmer serial, rustyhmmer under `mpirun` — on 117 profiles × 397,961
  proteins at `-np 8`. `--stall` too.
- **`--tformat` accepts all 19 of C's format names.** Added pressed NCBI BLAST
  databases (the version-4 index, header and sequence files, the ASN.1
  `Blast-def-line-set` parser, alias files listing volumes, and nucleotide
  databases, which C reads with an amino model and so does this), Farrar's
  `daemon` and `hmmpgmd`, and `fmindex`, which C rejects for every file.
- **`--restrictdb_stkey` / `--restrictdb_n` / `--ssifile`**, including the SSI
  index reader — primary keys, and accessions through the secondary-key table.
- **The FASTA reader is now Easel's, block for block.** The previous
  line-oriented reader silently diverged in three ways that only wide data
  exposes: trailing spaces in a description were trimmed where C keeps them,
  illegal characters were mapped to `X` where C fails the parse, and
  descriptions were not delimited at ctrl-A as NCBI headers require.
- **Command-line parsing is `esl_getopts`.** Unambiguous prefix abbreviation
  (`--tbl`), `--opt=value`, `-E1`, clustered flags, and every diagnostic
  verbatim — including the 24-character truncation and the lower-case "option"
  in one of the range messages. `-h` and the usage failures are C's text.
- **Failure paths match too**, which turns out to be most of the surface: which
  output files exist and what is in them when a run dies mid-search, the four
  places Easel aborts with a fatal exception (same message, same signal), and
  the three ways a multi-query run fails that a single-query run cannot reach.
- Fixed `fmt_g` on non-finite values — an empty database makes C print
  `(-nan)`, and this used to panic.
- Fixed `--cut_ga` / `--cut_tc` / `--cut_nc` against a model that carries no
  such cutoff. C names the kind (`GA bit thresholds unavailable on model %s`)
  and fails after the query header is written, leaving the report ending at the
  `Query:` line and the tabular files present and empty; this printed a message
  of its own shape and wrote nothing.
- **Rust 2024 edition**, MSRV `1.85`. The declared MSRV had been `1.74` while
  rayon already required `1.80`; cargo's MSRV-aware resolver silently picks
  older dependencies to satisfy a false declaration rather than failing, so
  this was not visible in a build.
- **Match HMMER's FPU mode (flush-to-zero, denormals-are-zero).** C calls
  `impl_Init()` at startup (`impl_sse/impl_sse.h:355-374`), which sets FTZ and DAZ
  in `MXCSR`; its own comment notes this is what makes even the scalar arithmetic
  agree across architectures. rustyhmmer ran IEEE-conforming, so an `exp(lnP)`
  that underflows produced a sub-normal where C produced exactly `0.0`. This was
  most of the remaining divergence: against C on 397,961 proteins x 117 profiles,
  differing `--tblout` rows fell from 31 to 18 of 1,453 (`--domtblout` 38 to 25 of
  1,989; with `--cut_ga`, 24 to 11 and 27 to 14), and on the E. coli set from 3 to
  1 of 227 rows. New `rustyhmmer::init()`; the CLI calls it and sets it on every
  rayon worker, since `MXCSR` is per-thread.
- **The Forward filter's E-state sum now uses C's striped lane order.**
  `omx::forward` already did; `forward::score` was a plain sequential sum, and
  floating-point addition is not associative, so the order is part of the
  transcription. (Measured: no output change on the test data — but it was a real
  structural divergence, and that score feeds the reported bit score.)
- **Striped, register-resident SSV filter kernel.** The SSV filter is ~50% of the
  run time (every model x sequence pair goes through it). It now runs HMMER's
  striped/banded form (`p7_SSVFilter`, `impl_sse/ssvfilter.c`): the DP row is held
  in SSE registers, up to 14 vectors at a time, and marched diagonally, so it
  never touches memory. **-43.7% instructions retired** end to end; output is
  **byte-identical** (differential test against the scalar oracle, plus md5 on
  397,961 proteins x 117 profiles).
- **`--cpu` now defaults to 2, as C does** (`p7_NCPU`, `p7_config.h:40`), capped at
  the core count (`hmmsearch.c:423`). rustyhmmer previously used every core, which
  on a many-core machine multiplied the per-worker DP matrices: peak resident
  memory on a 128-core host dropped from 2.5 GB to 219 MB (C: 119 MB). Pass
  `--cpu <n>` for more, exactly as with C.
- **The default human-readable report.** `hmmsearch`'s stdout output is now
  produced in full: banner, per-query score table, per-domain table and the
  alignment display with its CS/RF/consensus/match/sequence/PP lines, plus the
  internal pipeline statistics. New `-o`, `--textw`/`--notextw`. Previously
  rustyhmmer printed `--tblout` to stdout, which alone broke drop-in use.
- **HMM annotation lines are now parsed** (`RF`/`MM`/`CONS`/`CS`/`MAP`), with the
  per-format field layout C uses (`p7_hmmfile.c:1519-1537`) and the legacy
  `p7_hmm_SetConsensus()` synthesis for pre-3e files that carry no `CONS` column.
  Unlike C, the format is read per record rather than latched from the first one,
  so a file mixing `HMMER3/b` and `HMMER3/f` records reads through.
- **The optimal-accuracy traceback now carries per-position posterior
  probabilities** (`get_postprob()`, `p7_trace_Reverse()`'s N/C/J pull-back).
  These feed the alignment `PP` line. Verified by C's own invariant:
  `p7_trace_GetExpectedAccuracy(tr) == oasc`.
- **Full reporting/inclusion threshold system.** `p7_pipeline.c`'s `P7_PIPELINE`
  threshold, search-space and acceleration configuration is now transcribed
  (`src/pli.rs`), and the CLI exposes it: `-T --domE --domT --incE --incT
  --incdomE --incdomT -Z --domZ --max --F1 --F2 --F3 --nobias --nonull2 --seed
  --acc`. Reporting and inclusion counts come from `p7_tophits_Threshold()`
  instead of the previous ad-hoc rules.
- **Fix: `-T` was parsed and then ignored** — it is advertised in `--help` but
  never reached the reporting decision. Unrecognized options are now a usage
  error instead of a silent no-op.
- **Per-domain profile coordinates, and `--domtblout`.** Ported HMMER's
  optimal-accuracy alignment stage — `p7_OptimalAccuracy()` and `p7_OATrace()`
  from `impl_sse/optacc.c`, in the striped (Farrar) layout the C uses — so each
  domain now carries the profile coordinates of its alignment. New
  `--domtblout <f>` writes the full per-domain table
  (`p7_tophits_TabularDomains()`); `dp::Domain`, `pipeline::DomHit` and the new
  `api::HmmDomain` expose `hmm_from`/`hmm_to`, `ali_from`/`ali_to` and `oasc`,
  and `api::HmmHit` gained a `domains` list.
- Verified against C HMMER 3.4 on 397,961 proteins × 117 profiles (1,989 domain
  rows): every `hmm from`/`hmm to`/`ali from`/`ali to`/`env from`/`env to`/`acc`
  field is byte-identical, and `--tblout` is unchanged (same md5 as 0.1.3).
- **Fix: `--cut_ga`/`--cut_tc`/`--cut_nc` per-domain cutoffs.** Pfam and TIGRFAM
  terminate the cutoff line with a semicolon (`GA    116.25 116.25;`), so the
  second token was `"116.25;"`, which `str::parse::<f64>()` rejects — the
  per-domain cutoff silently became `0.0` and every domain passed. C uses
  `atof()`, which converts the leading numeric prefix; the parser now does the
  same. Per-*sequence* cutoffs were never affected. This changes the `rep`/`inc`
  columns of `--tblout` for `--cut_ga/tc/nc` runs on such files (E-value runs are
  unaffected, bit-for-bit).
- **Fix:** the into-M transitions for model node 0 (`amm`/`aim`/`adm` at k=1)
  held the raw HMM node-0 probabilities instead of C's `-eslINFINITY`
  (`p7_profile.c:84`, "node 0 nonexistent, has no transitions"). Forward and
  Backward multiply those by a zero cell so the values never mattered there, but
  the optimal-accuracy DP *selects* on `transition > 0`. `--tblout` output is
  unchanged.

### 0.1.3
- **Striped-SSE (Farrar) Viterbi filter.** The Viterbi filter score now runs a
  striped-SIMD i16 kernel (`p7_ViterbiFilter`, `impl_sse/vitfilter.c`) on
  SSE2-capable x86-64, built from the same quantized tables as the scalar path
  and dispatched at runtime. Output is **byte-identical** to the scalar oracle
  (verified by the `sse_matches_scalar` differential test); pure throughput.

### 0.1.2
- **Profile-parallel `search()`.** The in-process library `search()` now
  parallelizes over *profiles* (coarse-grained) instead of over sequences
  within each profile. On a genome-sized query set most profiles are rejected
  cheaply by the MSV filter, so the old per-profile sequence `par_iter` paid
  rayon fork/join overhead on ~10^5 near-empty dispatches (~35% core use on 16
  threads). The coarse dispatch is an order-preserving `IndexedParallelIterator`,
  so the concatenated output is **byte-identical** to before — pure throughput.

### 0.1.1
- Synced hmmsearch / Easel / DP code from the dev tree; full `--help`, graceful
  thread-pool error handling.

### 0.1.0
- Initial release: pure-Rust HMMER3 `hmmsearch` (`--tblout`), byte-parity with
  C HMMER 3.4; in-process `HmmAnnotator` library API.

## License

BSD-3-Clause. rustyhmmer is a derivative work of HMMER and its Easel library
and is distributed under the same terms. See [LICENSE](LICENSE) for the full
HMMER and Easel copyright notices.

## Author

Sunju Kim <n.e.coli.1822@gmail.com>
