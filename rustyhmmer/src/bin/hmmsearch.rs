










use rayon::prelude::*;
use std::process::exit;

#[derive(Default)]
struct Args {
    hmm: String,
    seq: String,
    tblout: Option<String>,
    domtblout: Option<String>,
    pfamtblout: Option<String>,
    afile: Option<String>,
    tformat: Option<String>,
    cut_ga: bool,
    cut_tc: bool,
    cut_nc: bool,
    // reporting thresholds
    evalue: Option<f64>,
    bitscore: Option<f64>,
    dom_evalue: Option<f64>,
    dom_bitscore: Option<f64>,
    // inclusion thresholds
    inc_evalue: Option<f64>,
    inc_bitscore: Option<f64>,
    incdom_evalue: Option<f64>,
    incdom_bitscore: Option<f64>,
    // search space
    z: Option<f64>,
    domz: Option<f64>,
    // acceleration
    max: bool,
    f1: Option<f64>,
    f2: Option<f64>,
    f3: Option<f64>,
    nobias: bool,
    nonull2: bool,
    // misc
    seed: Option<u32>,
    acc: bool,
    noali: bool,
    out: Option<String>,
    textw: Option<i32>,
    notextw: bool,
    cpu: Option<usize>,
    /// `--stall`: "arrest after start: for debugging MPI under gdb".
    stall: bool,
    /// `--mpi`: "run as an MPI parallel program".
    mpi: bool,
    /// `--restrictdb_stkey`: "Search starts at the sequence with name <s>".
    restrictdb_stkey: Option<String>,
    /// `--restrictdb_n`: "Search <j> target sequences (starting at --restrictdb_stkey)".
    restrictdb_n: Option<i64>,
    /// `--ssifile`: "restrictdb_x values require ssi file. Override default to <s>".
    ssifile: Option<String>,
    /// The literal string given for each value-bearing option, as typed.
    ///
    /// C's header lines are guarded by `esl_opt_IsUsed()`, and `esl_opt_IsDefault()`
    /// (easel/esl_getopts.c) decides by comparing *strings*:
    /// ```text
    ///   if (g->setby[opti] == eslARG_SETBY_DEFAULT)               return TRUE;
    ///   if (g->val[opti] == NULL && g->opt[opti].defval == NULL)  return TRUE;
    ///   if (esl_strcmp(g->opt[opti].defval, g->val[opti]) == 0)   return TRUE; /* option may have been set but restored to original default value */
    /// ```
    /// so `--cpu 2` against the default `"2"` is *not* "used" and prints no header
    /// line, while `--cpu 2.0` would. Reproducing that needs the token, not the
    /// parsed number.
    raw: std::collections::HashMap<&'static str, String>,
    /// Every option seen on the command line, flags mapped to `None`. Needed for
    /// [`Args::spoof_cmdline`], which unlike the header lines includes options given at
    /// their default value.
    set: std::collections::HashMap<&'static str, Option<String>>,
}

impl Args {
    /// C, easel/esl_getopts.c:`esl_opt_SpoofCmdline` -- `argv[0]`, then every option
    /// whose `setby` is not `eslARG_SETBY_DEFAULT` *in option-table order*, then the
    /// positional arguments; each token followed by a space, including the last.
    ///
    /// ```text
    ///   ntot = strlen(g->argv[0]) + 1;               // +1 for the space
    ///   snprintf(cmdline, ntot+1, "%s ", g->argv[0]);
    ///   for (i = 0; i < g->nopts; i++)
    ///     if (g->setby[i] != eslARG_SETBY_DEFAULT)
    ///       {
    ///         if (g->opt[i].type == eslARG_NONE) snprintf(cmdline + ntot, n+1, "%s ",    g->opt[i].name);
    ///         else                               snprintf(cmdline + ntot, n+1, "%s %s ", g->opt[i].name, g->val[i]);
    ///         ...
    ///       }
    ///   for (j = g->optind; j < g->argc; j++) { snprintf(cmdline + ntot, n+1, "%s ", g->argv[j]); ... }
    /// ```
    ///
    /// Note this is `setby != DEFAULT`, not `IsUsed`: an option given at its default
    /// value (`--cpu 2`) appears here even though it prints no header line.
    fn spoof_cmdline(&self) -> String {
        let mut out = std::env::args().next().unwrap_or_default();
        out.push(' ');
        for sp in OPTIONS {
            let name = sp.name;
            if let Some(val) = self.set.get(name) {
                out.push_str(name);
                if let Some(v) = val {
                    out.push(' ');
                    out.push_str(v);
                }
                out.push(' ');
            }
        }
        out.push_str(&self.hmm);
        out.push(' ');
        out.push_str(&self.seq);
        out.push(' ');
        out
    }

    /// `esl_opt_IsUsed()`: the option was given, and not at its default string.
    fn is_used(&self, flag: &str) -> bool {
        match self.raw.get(flag) {
            None => false,
            Some(given) => spec(flag)
                .and_then(|s| s.default)
                .is_none_or(|d| given != d),
        }
    }
}

/// C, src/errors.c:`p7_Fail`:
/// ```text
///   fprintf(stderr, "\nError: ");
///   vfprintf(stderr, format, argp);
///   fprintf(stderr, "\n");
///   fflush(stderr);
///   exit(1);
/// ```
/// `msg` is the already-formatted text, including whatever trailing newline the C format
/// string carries -- so the usual case ends in a blank line.
fn p7_fail(msg: &str) -> ! {
    eprint!("\nError: {msg}\n");
    exit(1);
}

/// C, easel/easel.c:`esl_fatal`, which writes the message and one newline, with no
/// `Error:` prefix. Takes bytes because one of the messages it carries -- Easel's
/// non-ASCII-character complaint -- interpolates a raw byte with `%c`.
fn esl_fatal(msg: &[u8]) -> ! {
    use std::io::Write;
    let mut err = std::io::stderr();
    let _ = err.write_all(msg);
    let _ = err.write_all(b"\n");
    let _ = err.flush();
    exit(1);
}

/// C's `hmmsearch -h`, verbatim.
///
/// hmmsearch.c builds it from its option table with `p7_banner`, `esl_usage` and seven
/// `esl_opt_DisplayHelp` calls (groups 1, 2, 4, 5, 6, 7, 12); group 99 -- the three
/// `--restrictdb_*`/`--ssifile` options -- is deliberately not shown. The result is fixed
/// text apart from the program name, which `p7_banner` and `esl_usage` take as the
/// *basename* of argv[0]; this crate reports itself as `hmmsearch` throughout so that its
/// output is identical whatever the binary is called.
const HELP: &str = "# hmmsearch :: search profile(s) against a sequence database
# HMMER 3.4 (Aug 2023); http://hmmer.org/
# Copyright (C) 2023 Howard Hughes Medical Institute.
# Freely distributed under the BSD open source license.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Usage: hmmsearch [options] <hmmfile> <seqdb>

Basic options:
  -h : show brief help on version and usage

Options directing output:
  -o <f>           : direct output to file <f>, not stdout
  -A <f>           : save multiple alignment of all hits to file <f>
  --tblout <f>     : save parseable table of per-sequence hits to file <f>
  --domtblout <f>  : save parseable table of per-domain hits to file <f>
  --pfamtblout <f> : save table of hits and domains to file, in Pfam format <f>
  --acc            : prefer accessions over names in output
  --noali          : don't output alignments, so output is smaller
  --notextw        : unlimit ASCII text output line width
  --textw <n>      : set max width of ASCII text output lines  [120]  (n>=120)

Options controlling reporting thresholds:
  -E <x>     : report sequences <= this E-value threshold in output  [10.0]  (x>0)
  -T <x>     : report sequences >= this score threshold in output
  --domE <x> : report domains <= this E-value threshold in output  [10.0]  (x>0)
  --domT <x> : report domains >= this score cutoff in output

Options controlling inclusion (significance) thresholds:
  --incE <x>    : consider sequences <= this E-value threshold as significant
  --incT <x>    : consider sequences >= this score threshold as significant
  --incdomE <x> : consider domains <= this E-value threshold as significant
  --incdomT <x> : consider domains >= this score threshold as significant

Options controlling model-specific thresholding:
  --cut_ga : use profile's GA gathering cutoffs to set all thresholding
  --cut_nc : use profile's NC noise cutoffs to set all thresholding
  --cut_tc : use profile's TC trusted cutoffs to set all thresholding

Options controlling acceleration heuristics:
  --max    : Turn all heuristic filters off (less speed, more power)
  --F1 <x> : Stage 1 (MSV) threshold: promote hits w/ P <= F1  [0.02]
  --F2 <x> : Stage 2 (Vit) threshold: promote hits w/ P <= F2  [1e-3]
  --F3 <x> : Stage 3 (Fwd) threshold: promote hits w/ P <= F3  [1e-5]
  --nobias : turn off composition bias filter

Other expert options:
  --nonull2     : turn off biased composition score corrections
  -Z <x>        : set # of comparisons done, for E-value calculation
  --domZ <x>    : set # of significant seqs, for domain E-value calculation
  --seed <n>    : set RNG seed to <n> (if 0: one-time arbitrary seed)  [42]
  --tformat <s> : assert target <seqfile> is in format <s>: no autodetection
  --cpu <n>     : number of parallel CPU workers to use for multithreads  [2]
  --stall       : arrest after start: for debugging MPI under gdb
  --mpi         : run as an MPI parallel program
";

/// C's `process_commandline` FAILURE label: the message, then usage, then the "basic"
/// option group, then a pointer at `-h`. All on stdout; exit 1.
fn usage_failure(msg: &str) -> ! {
    let argv0 = std::env::args().next().unwrap_or_default();
    println!("{msg}");
    println!("Usage: hmmsearch [options] <hmmfile> <seqdb>");
    println!("\nwhere most common options are:");
    println!("  -h : show brief help on version and usage");
    println!("\nTo see more help on available options, do {argv0} -h\n");
    exit(1);
}


/// The type of an option's argument, from `ESL_OPTIONS`'s `type` field.
#[derive(Clone, Copy, PartialEq, Eq)]
enum OptType {
    /// `eslARG_NONE`: a flag.
    None,
    /// `eslARG_INT`.
    Int,
    /// `eslARG_REAL`.
    Real,
    /// `eslARG_STRING` and `eslARG_OUTFILE`, which differ only in what `-h` prints.
    Str,
}

/// One row of hmmsearch.c's `static ESL_OPTIONS options[]`.
struct OptSpec {
    name: &'static str,
    ty: OptType,
    /// `defval`, as the literal string C stores, which is what `esl_opt_IsDefault`
    /// compares against.
    default: Option<&'static str>,
    envvar: Option<&'static str>,
    /// `range`, in Easel's `x>0` / `n>=120` notation.
    range: Option<&'static str>,
    /// `required_opts`.
    reqs: Option<&'static str>,
    /// `incompat_opts`.
    incompat: Option<&'static str>,
}

/// hmmsearch.c's option table, verbatim and in order. `esl_opt_SpoofCmdline` walks it,
/// `esl_opt_VerifyConfig` walks it, and `esl_opt_IsDefault` compares against the
/// defaults here, so one table drives all three.
///
/// ```text
///   #define REPOPTS     "-E,-T,--cut_ga,--cut_nc,--cut_tc"
///   #define DOMREPOPTS  "--domE,--domT,--cut_ga,--cut_nc,--cut_tc"
///   #define INCOPTS     "--incE,--incT,--cut_ga,--cut_nc,--cut_tc"
///   #define INCDOMOPTS  "--incdomE,--incdomT,--cut_ga,--cut_nc,--cut_tc"
///   #define THRESHOPTS  "-E,-T,--domE,--domT,--incE,--incT,--incdomE,--incdomT,--cut_ga,--cut_nc,--cut_tc"
///   #define CPUOPTS     "--mpi"
///   #define MPIOPTS     "--cpu"
/// ```
/// The `--restrictdb_*`/`--ssifile` rows are in docgroup 99, which `-h` never displays.
const OPTIONS: &[OptSpec] = &[
    o("-h", OptType::None, None, None, None),
    o("-o", OptType::Str, None, None, None),
    o("-A", OptType::Str, None, None, None),
    o("--tblout", OptType::Str, None, None, None),
    o("--domtblout", OptType::Str, None, None, None),
    o("--pfamtblout", OptType::Str, None, None, None),
    o("--acc", OptType::None, None, None, None),
    o("--noali", OptType::None, None, None, None),
    o("--notextw", OptType::None, None, None, Some("--textw")),
    OptSpec { name: "--textw", ty: OptType::Int, default: Some("120"), envvar: None, range: Some("n>=120"), reqs: None, incompat: Some("--notextw") },
    OptSpec { name: "-E", ty: OptType::Real, default: Some("10.0"), envvar: None, range: Some("x>0"), reqs: None, incompat: Some(REPOPTS) },
    OptSpec { name: "-T", ty: OptType::Real, default: None, envvar: None, range: None, reqs: None, incompat: Some(REPOPTS) },
    OptSpec { name: "--domE", ty: OptType::Real, default: Some("10.0"), envvar: None, range: Some("x>0"), reqs: None, incompat: Some(DOMREPOPTS) },
    OptSpec { name: "--domT", ty: OptType::Real, default: None, envvar: None, range: None, reqs: None, incompat: Some(DOMREPOPTS) },
    OptSpec { name: "--incE", ty: OptType::Real, default: Some("0.01"), envvar: None, range: Some("x>0"), reqs: None, incompat: Some(INCOPTS) },
    OptSpec { name: "--incT", ty: OptType::Real, default: None, envvar: None, range: None, reqs: None, incompat: Some(INCOPTS) },
    OptSpec { name: "--incdomE", ty: OptType::Real, default: Some("0.01"), envvar: None, range: Some("x>0"), reqs: None, incompat: Some(INCDOMOPTS) },
    OptSpec { name: "--incdomT", ty: OptType::Real, default: None, envvar: None, range: None, reqs: None, incompat: Some(INCDOMOPTS) },
    o("--cut_ga", OptType::None, None, None, Some(THRESHOPTS)),
    o("--cut_nc", OptType::None, None, None, Some(THRESHOPTS)),
    o("--cut_tc", OptType::None, None, None, Some(THRESHOPTS)),
    o("--max", OptType::None, None, None, Some("--F1,--F2,--F3")),
    OptSpec { name: "--F1", ty: OptType::Real, default: Some("0.02"), envvar: None, range: None, reqs: None, incompat: Some("--max") },
    OptSpec { name: "--F2", ty: OptType::Real, default: Some("1e-3"), envvar: None, range: None, reqs: None, incompat: Some("--max") },
    OptSpec { name: "--F3", ty: OptType::Real, default: Some("1e-5"), envvar: None, range: None, reqs: None, incompat: Some("--max") },
    o("--nobias", OptType::None, None, None, Some("--max")),
    o("--nonull2", OptType::None, None, None, None),
    OptSpec { name: "-Z", ty: OptType::Real, default: None, envvar: None, range: Some("x>0"), reqs: None, incompat: None },
    OptSpec { name: "--domZ", ty: OptType::Real, default: None, envvar: None, range: Some("x>0"), reqs: None, incompat: None },
    OptSpec { name: "--seed", ty: OptType::Int, default: Some("42"), envvar: None, range: Some("n>=0"), reqs: None, incompat: None },
    o("--tformat", OptType::Str, None, None, None),
    // p7_config.h:40 -- `#define p7_NCPU "2"`.
    OptSpec { name: "--cpu", ty: OptType::Int, default: Some("2"), envvar: Some("HMMER_NCPU"), range: Some("n>=0"), reqs: None, incompat: Some("--mpi") },
    o("--stall", OptType::None, None, Some("--mpi"), None),
    o("--mpi", OptType::None, None, None, Some("--cpu")),
    OptSpec { name: "--restrictdb_stkey", ty: OptType::Str, default: Some("0"), envvar: None, range: None, reqs: None, incompat: None },
    OptSpec { name: "--restrictdb_n", ty: OptType::Int, default: Some("-1"), envvar: None, range: None, reqs: None, incompat: None },
    o("--ssifile", OptType::Str, None, None, None),
];

const REPOPTS: &str = "-E,-T,--cut_ga,--cut_nc,--cut_tc";
const DOMREPOPTS: &str = "--domE,--domT,--cut_ga,--cut_nc,--cut_tc";
const INCOPTS: &str = "--incE,--incT,--cut_ga,--cut_nc,--cut_tc";
const INCDOMOPTS: &str = "--incdomE,--incdomT,--cut_ga,--cut_nc,--cut_tc";
const THRESHOPTS: &str =
    "-E,-T,--domE,--domT,--incE,--incT,--incdomE,--incdomT,--cut_ga,--cut_nc,--cut_tc";

/// A table row with no default, no env var and no range, which most of them are.
const fn o(
    name: &'static str,
    ty: OptType,
    range: Option<&'static str>,
    reqs: Option<&'static str>,
    incompat: Option<&'static str>,
) -> OptSpec {
    OptSpec { name, ty, default: None, envvar: None, range, reqs, incompat }
}

fn spec(name: &str) -> Option<&'static OptSpec> {
    OPTIONS.iter().find(|s| s.name == name)
}

/// C's `%.24s`, which every `esl_getopts` message uses.
fn t24(s: &str) -> String {
    s.chars().take(24).collect()
}

/// `esl_opt_ProcessCmdline`'s failures all end here: the message is prefixed and the
/// usage block follows.
fn cmdline_failure(msg: String) -> ! {
    usage_failure(&format!("Failed to parse command line: {msg}"));
}

/// `get_optidx_abbrev`: an option may be abbreviated to any unambiguous prefix, and an
/// exact match always wins.
///
/// ```text
///   ESL_FAIL(eslESYNTAX, g->errbuf, "Abbreviated option \"%.24s\" is ambiguous.", g->argv[g->optind]);
///   ESL_FAIL(eslESYNTAX, g->errbuf, "No such option \"%.24s\".", g->argv[g->optind]);
/// ```
fn match_long(tok: &str, name: &str) -> &'static OptSpec {
    if let Some(s) = spec(name) {
        return s;
    }
    let hits: Vec<&OptSpec> = OPTIONS.iter().filter(|s| s.name.starts_with(name)).collect();
    match hits.len() {
        1 => hits[0],
        0 => cmdline_failure(format!("No such option \"{}\".", t24(tok))),
        _ => cmdline_failure(format!("Abbreviated option \"{}\" is ambiguous.", t24(tok))),
    }
}

/// `verify_type_and_range`. The lower-case "option" in the integer-range message and the
/// upper-case "Option" everywhere else are C's, and are reproduced as they are.
///
/// ```text
///   if (! esl_str_IsInteger(val))
///     ESL_FAIL(..., "Option %.24s takes integer arg; got %.24s %s", ...);
///   if (verify_integer_range(val, g->opt[i].range) != eslOK)
///     ESL_FAIL(..., "option %.24s takes integer arg in range %.24s; got %.24s %s", ...);
///   if (! esl_str_IsReal(val))
///     ESL_FAIL(..., "Option %.24s takes real-valued arg; got %.24s %s", ...);
///   if (verify_real_range(val, g->opt[i].range) != eslOK)
///     ESL_FAIL(..., "Option %.24s takes real-valued arg in range %.24s; got %.24s %s", ...);
/// ```
/// `where` is "on cmdline" or "in env".
fn verify_type_and_range(sp: &OptSpec, val: &str, where_: &str) -> Result<(), String> {
    let in_range = |x: f64, range: &str| -> bool {
        // Easel's range grammar, restricted to the forms hmmsearch's table uses.
        let (op, num) = match range.find(|c: char| c.is_ascii_digit() || c == '-' || c == '.') {
            Some(i) => (&range[1..i], &range[i..]),
            None => return true,
        };
        let n: f64 = match num.parse() {
            Ok(v) => v,
            Err(_) => return true,
        };
        match op {
            ">" => x > n,
            ">=" => x >= n,
            "<" => x < n,
            "<=" => x <= n,
            _ => true,
        }
    };
    match sp.ty {
        OptType::None | OptType::Str => Ok(()),
        OptType::Int => {
            let Ok(v) = val.parse::<i64>() else {
                return Err(format!(
                    "Option {} takes integer arg; got {} {where_}",
                    t24(sp.name),
                    t24(val)
                ));
            };
            match sp.range {
                Some(r) if !in_range(v as f64, r) => Err(format!(
                    "option {} takes integer arg in range {}; got {} {where_}",
                    t24(sp.name),
                    t24(r),
                    t24(val)
                )),
                _ => Ok(()),
            }
        }
        OptType::Real => {
            let Ok(v) = val.parse::<f64>() else {
                return Err(format!(
                    "Option {} takes real-valued arg; got {} {where_}",
                    t24(sp.name),
                    t24(val)
                ));
            };
            match sp.range {
                Some(r) if !in_range(v, r) => Err(format!(
                    "Option {} takes real-valued arg in range {}; got {} {where_}",
                    t24(sp.name),
                    t24(r),
                    t24(val)
                )),
                _ => Ok(()),
            }
        }
    }
}

/// `esl_opt_VerifyConfig`, run after the command line is parsed:
/// ```text
///   /* For all options that are set ... verify that all their required_opts are set. */
///   ESL_FAIL(eslESYNTAX, g->errbuf, "Option %.24s requires (or has no effect without) option(s) %.24s", ...);
///   /* ... verify that no incompatible options are set to non-default values */
///   ESL_FAIL(eslESYNTAX, g->errbuf, "Option %.24s is incompatible with option(s) %.24s", ...);
/// ```
/// Both loops walk the table in order, and the message quotes the whole opt-list string,
/// not the option that actually clashed. The incompatibility check skips the option
/// itself, which is why a list may name it.
fn verify_config(set: &std::collections::HashMap<&'static str, Option<String>>) {
    for sp in OPTIONS {
        if !set.contains_key(sp.name) {
            continue;
        }
        if let Some(reqs) = sp.reqs {
            for r in reqs.split(',') {
                if !set.contains_key(r) {
                    cmdline_failure(format!(
                        "Option {} requires (or has no effect without) option(s) {}",
                        t24(sp.name),
                        t24(reqs)
                    ));
                }
            }
        }
    }
    for sp in OPTIONS {
        if !set.contains_key(sp.name) {
            continue;
        }
        if let Some(list) = sp.incompat {
            for other in list.split(',') {
                if other != sp.name && set.contains_key(other) {
                    cmdline_failure(format!(
                        "Option {} is incompatible with option(s) {}",
                        t24(sp.name),
                        t24(list)
                    ));
                }
            }
        }
    }
}

/// `esl_opt_ProcessCmdline` and `esl_opt_ProcessEnvironment`, for hmmsearch's table.
///
/// Option processing stops at the first argument that is not an option, so anything
/// after the two file names is counted as an extra argument rather than parsed -- which
/// is why `hmmsearch a.hmm b.fa -o out` is "Incorrect number of command line arguments".
fn parse_args() -> Args {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let mut set: std::collections::HashMap<&'static str, Option<String>> =
        std::collections::HashMap::new();
    let mut raw: std::collections::HashMap<&'static str, String> =
        std::collections::HashMap::new();

    let mut i = 0usize;
    let take_value = |sp: &'static OptSpec,
                          attached: Option<String>,
                          argv: &Vec<String>,
                          i: &mut usize|
     -> Option<String> {
        if sp.ty == OptType::None {
            if attached.is_some() {
                cmdline_failure(format!(
                    "Option {} does not take an argument.",
                    t24(sp.name)
                ));
            }
            return None;
        }
        let v = match attached {
            Some(v) => v,
            None => {
                if *i >= argv.len() {
                    cmdline_failure(format!("Option {} requires an argument.", t24(sp.name)));
                }
                let v = argv[*i].clone();
                *i += 1;
                // "Watch out for options masquerading as missing arguments."
                if sp.ty == OptType::Str && v.starts_with('-') {
                    cmdline_failure(format!(
                        "Arg looks like option? Use {}={} if you really mean it.",
                        t24(sp.name),
                        t24(&v)
                    ));
                }
                v
            }
        };
        if let Err(e) = verify_type_and_range(sp, &v, "on cmdline") {
            cmdline_failure(e);
        }
        Some(v)
    };

    while i < argv.len() {
        let tok = argv[i].clone();
        if tok == "-" || !tok.starts_with('-') {
            break; // start of the positional arguments
        }
        if let Some(body) = tok.strip_prefix("--") {
            // `--foo=arg` is split without touching argv.
            let (name, attached) = match body.split_once('=') {
                Some((n, v)) => (n.to_string(), Some(v.to_string())),
                None => (body.to_string(), None),
            };
            let sp = match_long(&tok, &format!("--{name}"));
            i += 1;
            let v = take_value(sp, attached, &argv, &mut i);
            if let Some(v) = &v {
                raw.insert(sp.name, v.clone());
            }
            set.insert(sp.name, v);
        } else {
            // `process_stdopt`: one-char options may be concatenated, and only the last
            // one in the string may take an argument -- either the rest of this element
            // or the next one.
            let chars: Vec<char> = tok[1..].chars().collect();
            i += 1;
            let mut c = 0usize;
            while c < chars.len() {
                let name = format!("-{}", chars[c]);
                let Some(sp) = spec(&name) else {
                    cmdline_failure(format!("No such option \"{}\".", t24(&name)));
                };
                c += 1;
                let attached = if sp.ty != OptType::None && c < chars.len() {
                    let rest: String = chars[c..].iter().collect();
                    c = chars.len();
                    Some(rest)
                } else {
                    None
                };
                let v = take_value(sp, attached, &argv, &mut i);
                if let Some(v) = &v {
                    raw.insert(sp.name, v.clone());
                }
                set.insert(sp.name, v);
            }
        }
    }

    // `-h` is handled before the argument count is checked.
    if set.contains_key("-h") {
        print!("{HELP}");
        exit(0);
    }

    // esl_opt_ProcessEnvironment runs *before* the command line in C, but only options
    // the command line did not set keep the environment's value, so checking it here is
    // the same -- except that its type errors carry a different prefix.
    for sp in OPTIONS {
        let Some(env) = sp.envvar else { continue };
        if set.contains_key(sp.name) {
            continue;
        }
        if let Ok(v) = std::env::var(env) {
            if let Err(e) = verify_type_and_range(sp, &v, "in env") {
                usage_failure(&format!("Failed to process environment: {e}"));
            }
            raw.insert(sp.name, v);
        }
    }

    verify_config(&set);

    let pos = &argv[i..];
    // C, hmmsearch.c:195 and its FAILURE label:
    //   if (esl_opt_ArgNumber(go) != 2) { puts("Incorrect number of command line arguments."); goto FAILURE; }
    //   ...
    //  FAILURE:  /* all errors handled here are user errors, so be polite.  */
    //   esl_usage(stdout, argv[0], usage);
    //   puts("\nwhere most common options are:");
    //   esl_opt_DisplayHelp(stdout, go, 1, 2, 80);
    //   printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
    //   exit(1);
    if pos.len() != 2 {
        usage_failure("Incorrect number of command line arguments.");
    }
    if pos[0] == "-" && pos[1] == "-" {
        usage_failure("Either <hmmfile> or <seqdb> may be '-' (to read from stdin), but not both.");
    }

    let num = |n: &str| -> Option<f64> { raw.get(n).and_then(|v| v.parse().ok()) };
    let int = |n: &str| -> Option<i64> { raw.get(n).and_then(|v| v.parse().ok()) };
    let string = |n: &str| -> Option<String> { raw.get(n).cloned() };
    let on = |n: &str| -> bool { set.contains_key(n) };

    Args {
        hmm: pos[0].clone(),
        seq: pos[1].clone(),
        tblout: string("--tblout"),
        domtblout: string("--domtblout"),
        pfamtblout: string("--pfamtblout"),
        afile: string("-A"),
        tformat: string("--tformat"),
        cut_ga: on("--cut_ga"),
        cut_tc: on("--cut_tc"),
        cut_nc: on("--cut_nc"),
        evalue: num("-E"),
        bitscore: num("-T"),
        dom_evalue: num("--domE"),
        dom_bitscore: num("--domT"),
        inc_evalue: num("--incE"),
        inc_bitscore: num("--incT"),
        incdom_evalue: num("--incdomE"),
        incdom_bitscore: num("--incdomT"),
        z: num("-Z"),
        domz: num("--domZ"),
        max: on("--max"),
        f1: num("--F1"),
        f2: num("--F2"),
        f3: num("--F3"),
        nobias: on("--nobias"),
        nonull2: on("--nonull2"),
        seed: int("--seed").map(|v| v as u32),
        acc: on("--acc"),
        noali: on("--noali"),
        out: string("-o"),
        textw: int("--textw").map(|v| v as i32),
        notextw: on("--notextw"),
        cpu: int("--cpu").map(|v| v as usize),
        stall: on("--stall"),
        mpi: on("--mpi"),
        restrictdb_stkey: string("--restrictdb_stkey"),
        restrictdb_n: int("--restrictdb_n"),
        ssifile: string("--ssifile"),
        raw,
        set,
    }
}


fn main() {
    // p7_Init() -> impl_Init(): flush-to-zero / denormals-are-zero, as C does.
    rustyhmmer::init();
    let args = parse_args();

    
    
    
    
    
    // hmmsearch.c:423 -- ncpus = ESL_MIN(--cpu, esl_threads_GetCPUCount()), and the
    // default for --cpu is p7_NCPU = 2 (p7_config.h:40), *not* every core. Matching
    // that matters for drop-in use: on a many-core machine the difference is a large
    // amount of resident memory, since each worker holds its own DP matrices.
    {
        let ncores = std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(1);
        // C's option table gives `--cpu` the environment variable `HMMER_NCPU`
        // (hmmsearch.c: `{ "--cpu", eslARG_INT, p7_NCPU, "HMMER_NCPU", "n>=0", ... }`),
        // so an unset `--cpu` falls back to that before falling back to p7_NCPU.
        let env_ncpu = std::env::var("HMMER_NCPU").ok().and_then(|v| v.parse().ok());
        let ncpus = args.cpu.or(env_ncpu).unwrap_or(2).min(ncores).max(1);
        rayon::ThreadPoolBuilder::new()
            .num_threads(ncpus)
            // MXCSR is per-thread; every worker needs C's FPU mode too.
            .start_handler(|_| rustyhmmer::init())
            .build_global()
            .unwrap_or_else(|e| {
                eprintln!("error: cannot start {ncpus} worker threads: {e}");
                exit(1);
            });
    }

    // C, hmmsearch.c:310-314:
    //   /* pause the execution of the programs execution until the user has a
    //    * chance to attach with a debugger and send a signal to resume execution
    //    * i.e. (gdb) signal SIGCONT
    //    */
    //   if (esl_opt_GetBoolean(go, "--stall")) pause();
    // The option table makes `--stall` require `--mpi`, so this is only reachable
    // together with it.
    if args.stall {
        loop {
            std::thread::park();
        }
    }

    // C, hmmsearch.c:316-330, guarded by `#ifdef HMMER_MPI`:
    //   if (esl_opt_GetBoolean(go, "--mpi"))
    //     {
    //       cfg.do_mpi     = TRUE;
    //       MPI_Init(&argc, &argv);
    //       MPI_Comm_rank(MPI_COMM_WORLD, &(cfg.my_rank));
    //       MPI_Comm_size(MPI_COMM_WORLD, &(cfg.nproc));
    //       if (cfg.my_rank > 0)  status = mpi_worker(go, &cfg);
    //       else 		    status = mpi_master(go, &cfg);
    //       MPI_Finalize();
    //     }
    //
    // Implemented in `rustyhmmer::mpisupport`, which takes the rank and world size from
    // the launcher's environment variables and carries HMMER's master/worker protocol
    // over TCP -- so no MPI library is linked. See that module for why that is enough.
    // As in C, a world of one deadlocks in the master loop rather than searching.
    let mpi = if args.mpi {
        match rustyhmmer::mpisupport::init() {
            Some(m) => Some(m),
            None => p7_fail("MPI_Init failed\n"),
        }
    } else {
        None
    };

    // C, hmmsearch.c:296-303, before anything is opened:
    //   if (esl_opt_IsUsed(go, "--restrictdb_stkey") )
    //     if ((cfg.firstseq_key = esl_opt_GetString(go, "--restrictdb_stkey")) == NULL)  p7_Fail("Failure capturing --restrictdb_stkey\n");
    //   if (esl_opt_IsUsed(go, "--restrictdb_n") )
    //     cfg.n_targetseq = esl_opt_GetInteger(go, "--restrictdb_n");
    //   if ( cfg.n_targetseq != -1 && cfg.n_targetseq < 1 )
    //     p7_Fail("--restrictdb_n must be >= 1\n");
    let firstseq_key: Option<&str> = if args.is_used("--restrictdb_stkey") {
        args.restrictdb_stkey.as_deref()
    } else {
        None
    };
    let n_targetseq: i64 = if args.is_used("--restrictdb_n") {
        args.restrictdb_n.unwrap_or(-1)
    } else {
        -1
    };
    if n_targetseq != -1 && n_targetseq < 1 {
        p7_fail("--restrictdb_n must be >= 1\n");
    }

    // C, hmmsearch.c:386-389:
    //   if (esl_opt_IsOn(go, "--tformat")) {
    //     dbfmt = esl_sqio_EncodeFormat(esl_opt_GetString(go, "--tformat"));
    //     if (dbfmt == eslSQFILE_UNKNOWN) p7_Fail("%s is not a recognized sequence database file format\n", ...);
    //   }
    // Without --tformat the format is autodetected (`eslSQFILE_UNKNOWN`); we only
    // autodetect FASTA, which is what every previous release assumed.
    let fmt = match &args.tformat {
        None => rustyhmmer::seqio::Format::Fasta,
        Some(name) => match rustyhmmer::seqio::encode_format(name) {
            Some(f) => f,
            None => p7_fail(&format!(
                "{name} is not a recognized sequence database file format\n"
            )),
        },
    };
    // C opens the target database *before* the query HMM file (hmmsearch.c:392 vs :408),
    // so when both are unusable it is the database's message that gets reported.
    //
    //   status = esl_sqfile_Open(cfg->dbfile, dbfmt, p7_SEQDBENV, &dbfp);
    //   if      (status == eslENOTFOUND) p7_Fail("Failed to open sequence file %s for reading\n",          cfg->dbfile);
    //   else if (status == eslEFORMAT)   p7_Fail("Sequence file %s is empty or misformatted\n",            cfg->dbfile);
    //   else if (status == eslEINVAL)    p7_Fail("Can't autodetect format of a stdin or .gz seqfile");
    //   else if (status != eslOK)        p7_Fail("Unexpected error %d opening sequence file %s\n", status, cfg->dbfile);
    //
    // The records themselves are not touched until the search loop below.
    let dbfp = {
        use rustyhmmer::seqio::SeqError;
        match rustyhmmer::seqio::open_format(&args.seq, fmt) {
            Ok(s) => s,
            Err(SeqError::NotFound) => p7_fail(&format!(
                "Failed to open sequence file {} for reading\n",
                args.seq
            )),
            Err(SeqError::OpenFormat) => p7_fail(&format!(
                "Sequence file {} is empty or misformatted\n",
                args.seq
            )),
            Err(SeqError::OpenStatus(s)) => p7_fail(&format!(
                "Unexpected error {s} opening sequence file {}\n",
                args.seq
            )),
            Err(e) => {
                eprintln!("error: cannot read {}: {e}", args.seq);
                exit(1);
            }
        }
    };

    // C, hmmsearch.c:399-404, immediately after the database opens:
    //   if (esl_opt_IsUsed(go, "--restrictdb_stkey") || esl_opt_IsUsed(go, "--restrictdb_n")) {
    //     if (esl_opt_IsUsed(go, "--ssifile"))
    //       esl_sqfile_OpenSSI(dbfp, esl_opt_GetString(go, "--ssifile"));
    //     else
    //       esl_sqfile_OpenSSI(dbfp, NULL);
    //   }
    //
    // The return value is discarded, so a missing index is not reported here -- it only
    // surfaces later, as an exception out of `esl_sqfile_PositionByKey`. Three of
    // `sqascii_OpenSSI`'s guards are `ESL_EXCEPTION`s, though, and those abort:
    //
    //   if (ascii->do_gzip)     ESL_EXCEPTION(eslEINVAL, "can't open an SSI index for a .gz compressed seq file");
    //   if (ascii->do_stdin)    ESL_EXCEPTION(eslEINVAL, "can't open an SSI index for standard input");
    //   if (ascii->afp != NULL) ESL_EXCEPTION(eslEINVAL, "can't open an SSI index for sequential input from an MSA");
    //
    // and the default SSI file name is the sequence file's with `.ssi` appended.
    let ssi = if firstseq_key.is_some() || n_targetseq != -1 {
        if fmt == rustyhmmer::seqio::Format::Ncbi {
            // C has no SSI entry points for a pressed BLAST database: `esl_sqncbi_Open`
            // never assigns `sqfp->open_ssi`, which `sqfile_open` left NULL, so
            // `esl_sqfile_OpenSSI` calls through a null pointer. Verified: C HMMER 3.4
            // segfaults here for either restrictdb option. Refusing is the one place
            // this program deliberately does something else.
            p7_fail(&format!(
                "Can't open an SSI index for the NCBI database {}\n",
                args.seq
            ));
        }
        if args.seq == "-" {
            esl_fatal_bytes_no_newline(
                b"Fatal exception (source file esl_sqio_ascii.c, line 1744):\n\
                  can't open an SSI index for standard input\n",
            );
        }
        if fmt.is_alignment() {
            esl_fatal_bytes_no_newline(
                b"Fatal exception (source file esl_sqio_ascii.c, line 1745):\n\
                  can't open an SSI index for sequential input from an MSA\n",
            );
        }
        let path = match &args.ssifile {
            Some(p) => p.clone(),
            None => format!("{}.ssi", args.seq),
        };
        rustyhmmer::ssi::Ssi::open(&path).ok()
    } else {
        None
    };

    // C, hmmsearch.c:408-410:
    //   status = p7_hmmfile_Open(cfg->hmmfile, NULL, &hfp, errbuf);
    //   if      (status == eslENOTFOUND) p7_Fail("File existence/permissions problem in trying to open HMM file %s.\n%s\n", cfg->hmmfile, errbuf);
    //   else if (status == eslEFORMAT)   p7_Fail("File format problem in trying to open HMM file %s.\n%s\n",                cfg->hmmfile, errbuf);
    // with p7_hmmfile.c:339 supplying the errbuf for the missing-file case.
    if args.hmm != "-" && !std::path::Path::new(&args.hmm).exists() {
        p7_fail(&format!(
            "File existence/permissions problem in trying to open HMM file {}.\nHMM file {} not found (nor an .h3m binary of it)\n",
            args.hmm, args.hmm
        ));
    }
    let hmms = match rustyhmmer::hmmfile::P7Hmm::read_all(&args.hmm) {
        Ok(h) if !h.is_empty() => h,
        Ok(_) => {
            eprintln!("error: no models in {}", args.hmm);
            exit(1);
        }
        Err(e) => {
            eprintln!("error: {e}");
            exit(1);
        }
    };
    // C, hmmsearch.c:414-419, "Open the results output files". These are `fopen(..., "w")`
    // calls, so every requested file exists -- and is truncated -- from this point on,
    // whether or not the run ever writes to it. A run that dies mid-search leaves the
    // ones it had not reached empty, and that is observable.
    //
    //   if (esl_opt_IsOn(go, "-o"))          { if ((ofp      = fopen(...)) == NULL) p7_Fail("Failed to open output file %s for writing\n", ...); }
    //   if (esl_opt_IsOn(go, "-A"))          { if ((afp      = fopen(...)) == NULL) p7_Fail("Failed to open alignment file %s for writing\n", ...); }
    //   if (esl_opt_IsOn(go, "--tblout"))    { if ((tblfp    = fopen(...)) == NULL) esl_fatal("Failed to open tabular per-seq output file %s for writing\n", ...); }
    //   if (esl_opt_IsOn(go, "--domtblout")) { if ((domtblfp = fopen(...)) == NULL) esl_fatal("Failed to open tabular per-dom output file %s for writing\n", ...); }
    //   if (esl_opt_IsOn(go, "--pfamtblout")){ if ((pfamtblfp= fopen(...)) == NULL) esl_fatal("Failed to open pfam-style tabular output file %s for writing\n", ...); }
    //
    // Note the two message styles: `-o` and `-A` go through p7_Fail, the tabular three
    // through esl_fatal.
    for (path, msg, is_p7_fail) in [
        (&args.out, "Failed to open output file", true),
        (&args.afile, "Failed to open alignment file", true),
        (&args.tblout, "Failed to open tabular per-seq output file", false),
        (&args.domtblout, "Failed to open tabular per-dom output file", false),
        (
            &args.pfamtblout,
            "Failed to open pfam-style tabular output file",
            false,
        ),
    ] {
        if let Some(p) = path {
            if std::fs::File::create(p).is_err() {
                let text = format!("{msg} {p} for writing\n");
                if is_p7_fail {
                    p7_fail(&text);
                } else {
                    esl_fatal(text.as_bytes());
                }
            }
        }
    }

    use rustyhmmer::domtblout::{
        format_rows as format_dom_rows, header as domtblout_header, DomWidths,
    };
    use rustyhmmer::pli::{BitCutoffs, Pipeline, Stats, ZSetBy};
    use rustyhmmer::report::{domains, pli_statistics, query_header, targets};
    use rustyhmmer::tblout::{format_row, header, hitlist, sort_hits, Widths};

    // p7_pipeline_Create() (p7_pipeline.c:108-235). The option-to-field mapping and
    // the order of the overrides (--cut_* after -T/--domT, --max after --F1/2/3) are
    // C's.
    let mut base = Pipeline::default();
    if let Some(v) = args.evalue { base.e = v; }
    if let Some(v) = args.dom_evalue { base.dom_e = v; }
    if let Some(v) = args.bitscore { base.t = v; base.by_e = false; }
    if let Some(v) = args.dom_bitscore { base.dom_t = v; base.dom_by_e = false; }
    if let Some(v) = args.inc_evalue { base.inc_e = v; }
    if let Some(v) = args.incdom_evalue { base.incdom_e = v; }
    if let Some(v) = args.inc_bitscore { base.inc_t = v; base.inc_by_e = false; }
    if let Some(v) = args.incdom_bitscore { base.incdom_t = v; base.incdom_by_e = false; }
    if args.cut_ga { base = base.with_bit_cutoffs(BitCutoffs::Ga); }
    if args.cut_nc { base = base.with_bit_cutoffs(BitCutoffs::Nc); }
    if args.cut_tc { base = base.with_bit_cutoffs(BitCutoffs::Tc); }
    let mut rep = String::new();
    output_header(&mut rep, &args);

    // hmmsearch.c reads the database inside the per-query search loop, and by the time
    // the first read happens it has written output_header() *and* the first query's
    // `Query:`/`Accession:`/`Description:` lines (hmmsearch.c:485-487). So that is
    // exactly what is on disk when a malformed record aborts the run.
    let seqs = match dbfp.read() {
        Ok(s) => s,
        Err(e) => {
            let mut prefix = rep.clone();
            if let Some(first) = hmms.first() {
                rustyhmmer::report::query_header(&mut prefix, first);
            }
            fail_on_seq_parse(&args, &prefix, e)
        }
    };

    // C, hmmsearch.c:481-485, once per query and just before the `Query:` line:
    //   if ( cfg->firstseq_key != NULL ) {
    //     sstatus = esl_sqfile_PositionByKey(dbfp, cfg->firstseq_key);
    //     if (sstatus != eslOK)
    //       p7_Fail("Failure setting restrictdb_stkey to %d\n", cfg->firstseq_key);
    //   }
    // and then `serial_loop`/`thread_loop` stop after `cfg->n_targetseq` sequences.
    //
    // The checks are the same for every query, so evaluating them once here is exact:
    // whichever way they go, they go that way on the first query.
    let seqs = {
        let mut seqs = seqs;
        let mut start = 0usize;
        if let Some(key) = firstseq_key {
            //   if (ascii->ssi == NULL) ESL_EXCEPTION(eslEINVAL,"Need an open SSI index to call esl_sqfile_PositionByKey()");
            // -- an exception, so this aborts with nothing flushed.
            let Some(ssi) = &ssi else {
                esl_fatal_bytes_no_newline(
                    b"Fatal exception (source file esl_sqio_ascii.c, line 1802):\n\
                      Need an open SSI index to call esl_sqfile_PositionByKey()\n",
                );
            };
            // p7_Fail's format string here is "%d" applied to a `char *`, so C prints
            // the low bits of a pointer -- a different number on every run, under ASLR.
            // There is no byte to match, so this prints the key that C meant to print.
            let Some(roff) = ssi.find_name(key) else {
                emit_report(&args, &rep);
                p7_fail(&format!("Failure setting restrictdb_stkey to {key}\n"));
            };
            match seqs.iter().position(|s| s.roff == roff) {
                Some(i) => start = i,
                None => {
                    // The index disagrees with the file. C would `fseek` to the offset
                    // and parse whatever bytes are there, with a message that depends on
                    // them; refusing is clearer than guessing which one.
                    emit_report(&args, &rep);
                    p7_fail(&format!(
                        "SSI index offset {roff} for {key} is not a record boundary in {}\n",
                        args.seq
                    ));
                }
            }
        }
        let end = if n_targetseq == -1 {
            seqs.len()
        } else {
            seqs.len().min(start + n_targetseq as usize)
        };
        seqs.truncate(end);
        seqs.drain(..start);
        seqs
    };

    // A worker never writes any output: it loops over the same queries the master does,
    // searching whatever blocks it is handed, and exits. `mpi_worker` (hmmsearch.c:1134)
    // is a separate function in C only because C's master duplicates the reporting code;
    // here the master simply continues down this function.
    if let Some(m) = &mpi {
        if m.rank() > 0 {
            let mut wbase = base.clone();
            wbase.z = seqs.len() as f64;
            if let Some(v) = args.z { wbase.z = v; wbase.z_setby = ZSetBy::Option; }
            for hmm in hmms {
                let model = rustyhmmer::pipeline::Model::new(hmm);
                let mut pli = wbase.clone();
                let cut = match pli.use_bit_cutoffs {
                    Some(BitCutoffs::Ga) => Some(model.hmm.cutoffs.ga),
                    Some(BitCutoffs::Tc) => Some(model.hmm.cutoffs.tc),
                    Some(BitCutoffs::Nc) => Some(model.hmm.cutoffs.nc),
                    None => None,
                };
                if let Some(c) = cut {
                    if pli.new_model_thresholds(c).is_err() {
                        // Same failure as the master's, below; a worker writes no
                        // report, so it only reports and dies.
                        let kind = match pli.use_bit_cutoffs {
                            Some(BitCutoffs::Ga) => "GA",
                            Some(BitCutoffs::Tc) => "TC",
                            Some(BitCutoffs::Nc) => "NC",
                            None => unreachable!("a cutoff pair without a cutoff kind"),
                        };
                        p7_fail(&format!(
                            "{kind} bit thresholds unavailable on model {}\n",
                            model.hmm.name
                        ));
                    }
                }
                m.worker_search(&model, &pli, &seqs);
            }
            m.worker_terminate();
            return;
        }
    }

    // The block list C builds lazily in `next_block` and caches in its `BLOCK_LIST`
    // for reuse across queries.
    let mpi_blocks = match &mpi {
        Some(_) => rustyhmmer::mpisupport::blocks(&seqs),
        None => Vec::new(),
    };

    // "Configure search space sizes for E value calculations" (p7_pipeline.c:201-213).
    // Unset means p7_ZSETBY_NTARGETS: Z is the sequence count, domZ the reported count.
    base.z = seqs.len() as f64;
    if let Some(v) = args.z { base.z = v; base.z_setby = ZSetBy::Option; }
    if let Some(v) = args.domz { base.domz = v; base.domz_setby = ZSetBy::Option; }
    // ESL_MIN(1.0, x) as C does (p7_pipeline.c:219-221).
    if let Some(v) = args.f1 { base.f1 = v.min(1.0); }
    if let Some(v) = args.f2 { base.f2 = v.min(1.0); }
    if let Some(v) = args.f3 { base.f3 = v.min(1.0); }
    if args.max { base = base.with_max(); }
    if args.nonull2 { base.do_null2 = false; }
    if args.nobias { base.do_biasfilter = false; }
    if let Some(v) = args.seed { base = base.with_seed(v); }
    base.show_accessions = args.acc;
    base.show_alignments = !args.noali;

    
    
    
    let mut out = String::new();
    let mut header_written = false;
    let mut domout = String::new();
    let mut pfamout = String::new();
    let mut aout = String::new();
    let mut dom_header_written = false;

    // hmmsearch.c: textw defaults to 120; --notextw makes it unlimited (0).
    let textw: i32 = if args.notextw { 0 } else { args.textw.unwrap_or(120) };

    for (nquery, hmm) in hmms.into_iter().enumerate() {
        let qt0 = std::time::Instant::now();
        let model = rustyhmmer::pipeline::Model::new(hmm);
        let qname = model.hmm.name.clone();
        let qacc = model.hmm.acc.clone().unwrap_or_default();

        // C rereads the whole database for every query after the first, and two things
        // can go wrong there (hmmsearch.c:470-476):
        //
        //   if (nquery > 1)
        //   {
        //     if (! esl_sqfile_IsRewindable(dbfp))
        //       esl_fatal("Target sequence file %s isn't rewindable; can't search it with multiple queries", cfg->dbfile);
        //     if (! esl_opt_IsUsed(go, "--restrictdb_stkey") )
        //       esl_sqfile_Position(dbfp, 0);
        //   }
        //
        // Reading the database once cannot reproduce either by accident, so both are
        // reproduced explicitly. Everything written so far is flushed first, as `exit()`
        // does for C's still-open output files -- so the tabular files keep the first
        // query's rows and lose their trailers.
        if nquery > 0 {
            // An alignment read from standard input reaches neither branch cleanly: C's
            // `IsRewindable` says yes (an MSA file never sets `do_stdin`), so it calls
            // `esl_sqfile_Position`, which closes the msafile and reopens `-` -- a stream
            // that is already at end of file. `esl_buffer.c:479` throws "invalid stream"
            // and Easel's default exception handler calls `abort()`, losing every stdio
            // buffer, so all the output files stay empty and the shell reports 134.
            // Reproduced with the same message and the same signal.
            if args.seq == "-" && fmt.is_alignment() {
                esl_fatal_bytes_no_newline(
                    b"Fatal exception (source file esl_buffer.c, line 479):\ninvalid stream\n",
                );
            }
            // `sqascii_IsRewindable` is `! (ascii->do_gzip || ascii->do_stdin)`, and
            // `do_stdin` is only ever set for a non-alignment format.
            if args.seq == "-" {
                flush_partial(&args, &rep, &out, &domout, &pfamout, &aout);
                esl_fatal(
                    format!(
                        "Target sequence file {} isn't rewindable; can't search it with multiple queries",
                        args.seq
                    )
                    .as_bytes(),
                );
            }
            // `esl_sqfile_Position(dbfp, 0)` rewinds to byte zero and reloads the buffer,
            // but nothing re-runs `fileheader_hmmpgmd` -- that only happens inside
            // `esl_sqascii_Open`. So the second query's first read finds the `#` header
            // line where a FASTA record should be, and every hmmpgmd database with more
            // than one query dies here. The message is fixed: `header_fasta` reports the
            // first non-space character, which the open already required to be `#`, and
            // `Position` resets `linenumber` to 1.
            if fmt == rustyhmmer::seqio::Format::Hmmpgmd {
                rustyhmmer::report::query_header(&mut rep, &model.hmm);
                flush_partial(&args, &rep, &out, &domout, &pfamout, &aout);
                let mut msg =
                    format!("Parse failed (sequence file {}):\n", args.seq).into_bytes();
                msg.extend_from_slice(
                    b"Line 1: unexpected char #; expected FASTA to start with >\n",
                );
                esl_fatal(&msg);
            }
        }

        // p7_pli_NewModelThresholds() (p7_pipeline.c:341-372): --cut_* thresholds are
        // per model, so each query gets its own pipeline.
        let mut pli = base.clone();
        let cut = match pli.use_bit_cutoffs {
            Some(BitCutoffs::Ga) => Some(model.hmm.cutoffs.ga),
            Some(BitCutoffs::Tc) => Some(model.hmm.cutoffs.tc),
            Some(BitCutoffs::Nc) => Some(model.hmm.cutoffs.nc),
            None => None,
        };
        if let Some(c) = cut {
            if pli.new_model_thresholds(c).is_err() {
                // C, p7_pipeline.c:543/550/557 -- one block per cutoff kind, and the
                // message names which:
                //   if (om->cutoff[p7_GA1] == p7_CUTOFF_UNSET)
                //     ESL_FAIL(eslEINVAL, pli->errbuf, "GA bit thresholds unavailable on model %s\n", om->name);
                // which hmmsearch.c:504 hands straight to p7_Fail:
                //   status = p7_pli_NewModel(info[i].pli, info[i].om, info[i].bg);
                //   if (status == eslEINVAL) p7_Fail(info->pli->errbuf);
                //
                // This sits *after* the query header is written, so the report on disk
                // ends with the `Query:` line; exit() then flushes it, and the tabular
                // files stay as `fopen` left them -- present and empty.
                let kind = match pli.use_bit_cutoffs {
                    Some(BitCutoffs::Ga) => "GA",
                    Some(BitCutoffs::Tc) => "TC",
                    Some(BitCutoffs::Nc) => "NC",
                    None => unreachable!("a cutoff pair without a cutoff kind"),
                };
                let mut prefix = rep.clone();
                query_header(&mut prefix, &model.hmm);
                emit_report(&args, &prefix);
                p7_fail(&format!(
                    "{kind} bit thresholds unavailable on model {qname}\n"
                ));
            }
        }

        // With --mpi the master does not search: it hands blocks to the workers and
        // merges what comes back (hmmsearch.c:940-1030). Without it, every target is
        // searched here.
        let mpi_result = mpi.as_ref().map(|m| m.master_search(&mpi_blocks));

        let (mut hits, st): (Vec<rustyhmmer::pipeline::Hit>, Stats) = match mpi_result {
            Some(r) => r,
            None => seqs
            .par_iter()
            .enumerate()
            .map(|(idx, s)| {
                let (h, st) = model.search_one_counted(s, &pli);
                (
                    h.map(|mut h| {
                        h.seq_idx = idx;
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
            ),
        };
        // Counters are per query: hmmsearch.c creates a fresh P7_PIPELINE inside the
        // loop, so p7_pli_Statistics() reports this query only (hmmsearch.c:541-553).
        let mut qstats = st;
        qstats.nmodels = 1;
        qstats.nnodes = model.hmm.m as u64;
        sort_hits(&mut hits, &pli);

        // C's `th` (see tblout::hitlist) is what the column widths are measured over;
        // `reported` is what p7_tophits_Threshold() then actually prints.
        let th = hitlist(&hits, &pli);
        let reported: Vec<&rustyhmmer::pipeline::Hit> = hits
            .iter()
            .filter(|h| pli.target_reportable(h.score, h.lnp))
            .collect();

        // "Now we can determine domZ, the effective search space in which additional
        // domains are found" (p7_tophits.c:962-963).
        if pli.domz_setby == ZSetBy::NTargets {
            pli.domz = reported.len() as f64;
        }

        // Human-readable report for this query (hmmsearch.c:487-553).
        query_header(&mut rep, &model.hmm);
        targets(&mut rep, &reported, &th, &pli, textw);
        rep.push_str("\n\n");
        domains(
            &mut rep,
            &reported,
            &pli,
            textw,
            &model.hmm,
            &model.fwd,
            &seqs,
        );
        rep.push_str("\n\n");

        pli_statistics(&mut rep, &pli, &qstats, Some(qt0.elapsed().as_secs_f64()));

        // C, hmmsearch.c:555-574: the -A block comes *after* `p7_pli_Statistics()` and
        // after the `//` that closes the query, so its report line is the last thing in
        // the per-query block -- past the `//`, not before the statistics.
        if args.afile.is_some() {
            match write_alignment(&mut aout, &reported, &pli, &model, &seqs, textw) {
                Some(nseq) => rep.push_str(&format!(
                    "# Alignment of {nseq} hits satisfying inclusion thresholds saved to: {}\n",
                    args.afile.as_deref().unwrap_or("")
                )),
                None => rep
                    .push_str("# No hits satisfy inclusion thresholds; no alignment saved\n"),
            }
        }

        let w = Widths::compute(&qname, &qacc, &th);
        if !header_written {
            out.push_str(&header(&w));
            header_written = true;
        }
        for h in &reported {
            out.push_str(&format_row(h, &qname, &qacc, &pli, &w));
            out.push('\n');
        }

        if args.domtblout.is_some() {
            let dw = DomWidths::compute(&qname, &qacc, &th);
            if !dom_header_written {
                domout.push_str(&domtblout_header(&dw));
                dom_header_written = true;
            }
            let qlen = model.hmm.m;
            for h in &reported {
                domout.push_str(&format_dom_rows(
                    h, h.tlen, &qname, &qacc, qlen, &pli, &dw,
                ));
            }
        }

        // Unlike --tblout, C emits the whole Pfam-style block per query
        // (hmmsearch.c:549 is inside the query loop and `p7_tophits_TabularXfam` has no
        // show-header guard), so these repeat.
        if args.pfamtblout.is_some() {
            let pw = rustyhmmer::pfamtblout::tnamew(&th);
            rustyhmmer::pfamtblout::query_blocks(&mut pfamout, &reported, &pli, pw);
        }
    }
    // "monitor all the workers to make sure they have ended" (hmmsearch.c:1085).
    if let Some(m) = &mpi {
        m.master_wait_for_termination();
    }

    if !header_written {
        out.push_str(&header(&Widths { tnamew: 20, taccw: 10, qnamew: 20, qaccw: 10 }));
    }
    write_tblout_trailer(&mut out, &args);

    if let Some(path) = &args.tblout {
        if let Err(e) = std::fs::write(path, &out) {
            eprintln!("error: cannot write {path}: {e}");
            exit(1);
        }
    }

    if let Some(path) = &args.afile {
        if let Err(e) = std::fs::write(path, &aout) {
            eprintln!("error: cannot write {path}: {e}");
            exit(1);
        }
    }

    if let Some(path) = &args.pfamtblout {
        // hmmsearch.c:599 writes the standard tabular tail here too.
        write_tblout_trailer(&mut pfamout, &args);
        if let Err(e) = std::fs::write(path, &pfamout) {
            eprintln!("error: cannot write {path}: {e}");
            exit(1);
        }
    }

    {
        rep.push_str("[ok]\n");
        emit_report(&args, &rep);
    }

    if let Some(path) = &args.domtblout {
        if !dom_header_written {
            domout.push_str(&domtblout_header(&DomWidths {
                tnamew: 20,
                taccw: 10,
                qnamew: 20,
                qaccw: 10,
            }));
        }
        write_tblout_trailer(&mut domout, &args);
        if let Err(e) = std::fs::write(path, &domout) {
            eprintln!("error: cannot write {path}: {e}");
            exit(1);
        }
    }
}






/// C, hmmsearch.c:517-527, reached when `esl_sqio_ReadBlock` fails mid-search:
/// ```text
///   case eslEFORMAT:
///     esl_fatal("Parse failed (sequence file %s):\n%s\n",
///         dbfp->filename, esl_sqfile_GetErrorBuf(dbfp));
///     break;
///   case eslEOF:
///     /* do nothing */
///     break;
///   default:
///     esl_fatal("Unexpected error %d reading sequence file %s", sstatus, dbfp->filename);
/// ```
/// Whatever has already been written to `ofp` is flushed by `exit()`, so the report
/// buffer goes out first.
fn fail_on_seq_parse(args: &Args, rep: &str, e: rustyhmmer::seqio::SeqError) -> ! {
    use rustyhmmer::seqio::SeqError;
    emit_report(args, rep);
    match e {
        SeqError::Parse(errbuf) => {
            let mut msg = format!("Parse failed (sequence file {}):\n", args.seq).into_bytes();
            msg.extend_from_slice(&errbuf);
            msg.push(b'\n');
            esl_fatal(&msg);
        }
        other => {
            eprintln!("error: cannot read {}: {other}", args.seq);
            exit(1);
        }
    }
}

/// Easel's default exception handler, easel.c:`esl_exception`, which ends in `abort()`.
/// Unlike `esl_fatal` it adds no newline of its own, and unlike `exit()` it flushes
/// nothing -- so any output file the run had opened is left as it was on disk.
fn esl_fatal_bytes_no_newline(msg: &[u8]) -> ! {
    use std::io::Write;
    let mut err = std::io::stderr();
    let _ = err.write_all(msg);
    let _ = err.flush();
    std::process::abort();
}

/// Write every output file with exactly what the run produced before it aborted.
///
/// C holds all of these open and writes as it goes, so `exit()` flushes partial content:
/// the tabular files keep the completed queries' rows but never get their `#` trailer,
/// and the report never gets its `# [ok]`.
fn flush_partial(args: &Args, rep: &str, tbl: &str, dom: &str, pfam: &str, ali: &str) {
    emit_report(args, rep);
    for (path, body) in [
        (&args.tblout, tbl),
        (&args.domtblout, dom),
        (&args.pfamtblout, pfam),
        (&args.afile, ali),
    ] {
        if let Some(p) = path {
            let _ = std::fs::write(p, body);
        }
    }
}

/// Send the human-readable report to `-o`'s file, or to stdout when it is unset.
///
/// C writes to `ofp` as it goes, and `exit()` flushes it; so does the abort path in
/// [`fail_on_seq_parse`], where the header is already in the buffer when the first read
/// fails.
fn emit_report(args: &Args, rep: &str) {
    match &args.out {
        Some(path) => {
            if let Err(e) = std::fs::write(path, rep) {
                eprintln!("error: cannot write {path}: {e}");
                exit(1);
            }
        }
        None => print!("{rep}"),
    }
}

/// `output_header()` — hmmsearch.c:220-262. Only the options actually given are
/// echoed, in C's order.
fn output_header(out: &mut String, args: &Args) {
    out.push_str("# hmmsearch :: search profile(s) against a sequence database\n");
    out.push_str("# HMMER 3.4 (Aug 2023); http://hmmer.org/\n");
    out.push_str("# Copyright (C) 2023 Howard Hughes Medical Institute.\n");
    out.push_str("# Freely distributed under the BSD open source license.\n");
    out.push_str("# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
    out.push_str(&format!("# query HMM file:                  {}\n", args.hmm));
    out.push_str(&format!("# target sequence database:        {}\n", args.seq));
    let g = |x: f64| rustyhmmer::tblout::fmt_g(x, 6);
    if let Some(v) = &args.out { out.push_str(&format!("# output directed to file:         {v}\n")); }
    if let Some(v) = &args.afile { out.push_str(&format!("# MSA of all hits saved to file:   {v}\n")); }
    if let Some(v) = &args.tblout { out.push_str(&format!("# per-seq hits tabular output:     {v}\n")); }
    if let Some(v) = &args.domtblout { out.push_str(&format!("# per-dom hits tabular output:     {v}\n")); }
    if let Some(v) = &args.pfamtblout { out.push_str(&format!("# pfam-style tabular hit output:   {v}\n")); }
    if args.acc { out.push_str("# prefer accessions over names:    yes\n"); }
    if args.noali { out.push_str("# show alignments in output:       no\n"); }
    if args.notextw { out.push_str("# max ASCII text line length:      unlimited\n"); }
    if let Some(v) = args.textw { if args.is_used("--textw") { out.push_str(&format!("# max ASCII text line length:      {v}\n")); } }
    if let Some(v) = args.evalue { if args.is_used("-E") { out.push_str(&format!("# sequence reporting threshold:    E-value <= {}\n", g(v))); } }
    if let Some(v) = args.bitscore { out.push_str(&format!("# sequence reporting threshold:    score >= {}\n", g(v))); }
    if let Some(v) = args.dom_evalue { if args.is_used("--domE") { out.push_str(&format!("# domain reporting threshold:      E-value <= {}\n", g(v))); } }
    if let Some(v) = args.dom_bitscore { out.push_str(&format!("# domain reporting threshold:      score >= {}\n", g(v))); }
    if let Some(v) = args.inc_evalue { if args.is_used("--incE") { out.push_str(&format!("# sequence inclusion threshold:    E-value <= {}\n", g(v))); } }
    if let Some(v) = args.inc_bitscore { out.push_str(&format!("# sequence inclusion threshold:    score >= {}\n", g(v))); }
    if let Some(v) = args.incdom_evalue { if args.is_used("--incdomE") { out.push_str(&format!("# domain inclusion threshold:      E-value <= {}\n", g(v))); } }
    if let Some(v) = args.incdom_bitscore { out.push_str(&format!("# domain inclusion threshold:      score >= {}\n", g(v))); }
    if args.cut_ga { out.push_str("# model-specific thresholding:     GA cutoffs\n"); }
    if args.cut_nc { out.push_str("# model-specific thresholding:     NC cutoffs\n"); }
    if args.cut_tc { out.push_str("# model-specific thresholding:     TC cutoffs\n"); }
    if args.max { out.push_str("# Max sensitivity mode:            on [all heuristic filters off]\n"); }
    if let Some(v) = args.f1 { if args.is_used("--F1") { out.push_str(&format!("# MSV filter P threshold:       <= {}\n", g(v))); } }
    if let Some(v) = args.f2 { if args.is_used("--F2") { out.push_str(&format!("# Vit filter P threshold:       <= {}\n", g(v))); } }
    if let Some(v) = args.f3 { if args.is_used("--F3") { out.push_str(&format!("# Fwd filter P threshold:       <= {}\n", g(v))); } }
    if args.nobias { out.push_str("# biased composition HMM filter:   off\n"); }
    if let Some(v) = &args.restrictdb_stkey { if args.is_used("--restrictdb_stkey") { out.push_str(&format!("# Restrict db to start at seq key: {v}\n")); } }
    if let Some(v) = args.restrictdb_n { if args.is_used("--restrictdb_n") { out.push_str(&format!("# Restrict db to # target seqs:    {v}\n")); } }
    if let Some(v) = &args.ssifile { out.push_str(&format!("# Override ssi file to:            {v}\n")); }
    if args.nonull2 { out.push_str("# null2 bias corrections:          off\n"); }
    if let Some(v) = args.z { out.push_str(&format!("# sequence search space set to:    {v:.0}\n")); }
    if let Some(v) = args.domz { out.push_str(&format!("# domain search space set to:      {v:.0}\n")); }
    if let Some(v) = args.seed {
        if args.is_used("--seed") {
            // C, hmmsearch.c:256-259:
            //   if (esl_opt_IsUsed(go, "--seed"))  {
            //     if (esl_opt_GetInteger(go, "--seed") == 0 && fprintf(ofp, "# random number seed:              one-time arbitrary\n") < 0) ESL_EXCEPTION_SYS(...);
            //     else if (                               fprintf(ofp, "# random number seed set to:       %d\n", esl_opt_GetInteger(go, "--seed")) < 0) ESL_EXCEPTION_SYS(...);
            //   }
            // The `&&` binds the write into the *condition*, so when seed is 0 and the
            // write succeeds the condition is false and the `else if` runs too: `--seed 0`
            // prints both lines. Verified against C HMMER 3.4.
            if v == 0 { out.push_str("# random number seed:              one-time arbitrary\n"); }
            out.push_str(&format!("# random number seed set to:       {v}\n"));
        }
    }
    if let Some(v) = &args.tformat { out.push_str(&format!("# targ <seqfile> format asserted:  {v}\n")); }
    if let Some(v) = args.cpu { if args.is_used("--cpu") { out.push_str(&format!("# number of worker threads:        {v}\n")); } }
    // hmmsearch.c:267 -- the last of the header lines, after --cpu's.
    if args.mpi { out.push_str("# MPI:                             on\n"); }
    // C ends the header with "...\n\n" (hmmsearch.c:269), so the blank line
    // belongs to the header, not to each query block.
    out.push_str("# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");
}

/// C, p7_tophits.c:`p7_tophits_TabularTail`:
/// ```text
///   if (fprintf(ofp, "#\n") < 0) ...
///   if (fprintf(ofp, "# Program:         %s\n",      (progname == NULL) ? "[none]" : progname) < 0) ...
///   if (fprintf(ofp, "# Version:         %s (%s)\n", HMMER_VERSION, HMMER_DATE) < 0) ...
///   if (fprintf(ofp, "# Pipeline mode:   %s\n",      modestamp) < 0) ...
///   if (fprintf(ofp, "# Query file:      %s\n",      (qfile    == NULL) ? "[none]" : qfile) < 0) ...
///   if (fprintf(ofp, "# Target file:     %s\n",      (tfile    == NULL) ? "[none]" : tfile) < 0) ...
///   if (fprintf(ofp, "# Option settings: %s\n",      spoof_cmd) < 0) ...
///   if (fprintf(ofp, "# Current dir:     %s\n",      (cwd      == NULL) ? "[unknown]" : cwd) < 0) ...
///   if (fprintf(ofp, "# Date:            %s",        timestamp) < 0) /* timestamp ends in \n */ ...
///   if (fprintf(ofp, "# [ok]\n") < 0) ...
/// ```
///
/// The last three lines were missing entirely. Their *contents* are in
/// DEVELOPMENT_PRINCIPLES §7's documented-variable set (command line, path, date), but
/// their absence shifted the line count of every tabular file.
fn write_tblout_trailer(out: &mut String, args: &Args) {
    out.push_str("#\n");
    out.push_str("# Program:         hmmsearch\n");
    out.push_str("# Version:         3.4 (Aug 2023)\n");
    out.push_str("# Pipeline mode:   SEARCH\n");
    out.push_str(&format!("# Query file:      {}\n", args.hmm));
    out.push_str(&format!("# Target file:     {}\n", args.seq));
    out.push_str(&format!("# Option settings: {}\n", args.spoof_cmdline()));
    let cwd = std::env::current_dir()
        .map(|p| p.display().to_string())
        .unwrap_or_else(|_| "[unknown]".to_string());
    out.push_str(&format!("# Current dir:     {cwd}\n"));
    out.push_str(&format!("# Date:            {}", asctime_now()));
    out.push_str("# [ok]\n");
}

/// `ctime_r()`'s text, which is `asctime()`'s: `"Thu Jul 30 11:24:57 2026\n"`.
///
/// glibc's `asctime` format is `"%.3s %.3s%3d %.2d:%.2d:%.2d %d\n"` -- note the day is
/// `%3d`, so it supplies the separating space itself and single-digit days get two
/// spaces (`"Jul  5"`).
///
/// Local time, like `ctime_r()`: `std` has no calendar or time zone support, so the UTC
/// offset comes from [`local_utc_offset`].
fn asctime_now() -> String {
    const WDAY: [&str; 7] = ["Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat"];
    const MON: [&str; 12] = [
        "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
    ];
    let utc = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_secs() as i64)
        .unwrap_or(0);
    let secs = utc + local_utc_offset(utc);
    let days = secs.div_euclid(86_400);
    let tod = secs.rem_euclid(86_400);
    // 1970-01-01 was a Thursday.
    let wday = (days + 4).rem_euclid(7) as usize;
    // civil_from_days (Howard Hinnant), shifting the era to start in March so that the
    // leap day lands at the end of the year.
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z.rem_euclid(146_097);
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = if m <= 2 { y + 1 } else { y };
    format!(
        "{} {}{:>3} {:02}:{:02}:{:02} {}\n",
        WDAY[wday],
        MON[(m - 1) as usize],
        d,
        tod / 3600,
        (tod % 3600) / 60,
        tod % 60,
        y
    )
}

/// Seconds to add to a UTC epoch second to get local time, read from the platform time
/// zone database.
///
/// `ctime_r()` renders local time, and `std` exposes neither a calendar nor the time
/// zone, so the offset is taken from the TZif file (RFC 8536) that `TZ` or
/// `/etc/localtime` names. Returns 0 if none can be read, which degrades to UTC rather
/// than failing.
///
/// Only the `utoff` of the type in effect at `utc` is needed; DST designations, leap
/// seconds and the POSIX footer string are skipped.
fn local_utc_offset(utc: i64) -> i64 {
    fn parse_tzif(buf: &[u8], utc: i64) -> Option<i64> {
        let be32 = |b: &[u8]| -> u32 { u32::from_be_bytes([b[0], b[1], b[2], b[3]]) };

        // RFC 8536 §3.1: "TZif", version, 15 reserved, then six counts.
        let hdr = |at: usize| -> Option<[u32; 6]> {
            let h = buf.get(at..at + 44)?;
            if &h[0..4] != b"TZif" {
                return None;
            }
            let mut c = [0u32; 6];
            for (i, v) in c.iter_mut().enumerate() {
                *v = be32(&h[20 + i * 4..]);
            }
            Some(c)
        };
        // §3.2 data block size, for `time_size` bytes per transition time.
        let block = |c: &[u32; 6], time_size: usize| -> usize {
            c[3] as usize * (time_size + 1)   // transition times + type indices
                + c[4] as usize * 6           // local time type records
                + c[5] as usize               // time zone designations
                + c[2] as usize * (time_size + 4) // leap-second records
                + c[1] as usize               // standard/wall indicators
                + c[0] as usize // UT/local indicators
        };

        let c1 = hdr(0)?;
        let version = buf[4];
        // Version 2+ repeats everything with 64-bit transition times; that block wins.
        let (at, c, tsz) = if version >= b'2' {
            let at = 44 + block(&c1, 4);
            (at + 44, hdr(at)?, 8)
        } else {
            (44, c1, 4)
        };

        let (timecnt, typecnt) = (c[3] as usize, c[4] as usize);
        if typecnt == 0 {
            return None;
        }
        let times = buf.get(at..at + timecnt * tsz)?;
        let idxs = buf.get(at + timecnt * tsz..at + timecnt * (tsz + 1))?;
        let types_at = at + timecnt * (tsz + 1);

        let at_time = |i: usize| -> i64 {
            let b = &times[i * tsz..];
            if tsz == 8 {
                i64::from_be_bytes([b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7]])
            } else {
                be32(b) as i32 as i64
            }
        };
        // Last transition at or before `utc`; before the first, type 0 applies.
        let mut ty = 0usize;
        for i in 0..timecnt {
            if at_time(i) > utc {
                break;
            }
            ty = idxs[i] as usize;
        }
        if ty >= typecnt {
            return None;
        }
        let rec = buf.get(types_at + ty * 6..types_at + ty * 6 + 6)?;
        Some(be32(rec) as i32 as i64)
    }

    let path = match std::env::var("TZ") {
        // A leading ':' is a common way to force a file path; a bare name is relative to
        // the zoneinfo directory. Anything else (a POSIX rule string) we do not parse.
        Ok(tz) if !tz.is_empty() => {
            let tz = tz.strip_prefix(':').unwrap_or(&tz);
            if tz.starts_with('/') {
                tz.to_string()
            } else {
                format!("/usr/share/zoneinfo/{tz}")
            }
        }
        _ => "/etc/localtime".to_string(),
    };
    std::fs::read(&path)
        .ok()
        .and_then(|buf| parse_tzif(&buf, utc))
        .unwrap_or(0)
}

/// Build and append one query's `-A` alignment; returns the number of sequences in it, or
/// `None` when nothing satisfied the inclusion thresholds.
///
/// C, p7_tophits.c:`p7_tophits_Alignment` then hmmsearch.c:559-567:
/// ```text
///   for (h = 0; h < th->N; h++)
///     if (th->hit[h]->flags & p7_IS_INCLUDED)
///       for (d = 0; d < th->hit[h]->ndom; d++)
///         if (th->hit[h]->dcl[d].is_included)
///           {
///             if (M == 0) M = th->hit[h]->dcl[d].ad->M;
///             ndom++;
///           }
///   ...
///   if (inc_n+ndom == 0) { status = eslFAIL; goto ERROR; }
///   ...
///           p7_alidisplay_Backconvert(th->hit[h]->dcl[d].ad, abc, &(sqarr[y]), &(trarr[y]))
///   ...
///   p7_tracealign_Seqs(sqarr, trarr, inc_n+ndom, M, optflags, NULL, &msa)
/// ```
/// ```text
///   esl_msa_SetName     (msa, hmm->name, -1);
///   esl_msa_SetAccession(msa, hmm->acc,  -1);
///   esl_msa_SetDesc     (msa, hmm->desc, -1);
///   esl_msa_FormatAuthor(msa, "hmmsearch (HMMER %s)", HMMER_VERSION);
///
///   if (textw > 0) esl_msafile_Write(afp, msa, eslMSAFILE_STOCKHOLM);
///   else           esl_msafile_Write(afp, msa, eslMSAFILE_PFAM);
/// ```
///
/// `M` is taken from the first included domain's alidisplay, not from the HMM -- C has no
/// `hmm` here ("we don't have hmm, but every ali has a copy"). They agree, but the
/// alidisplay is the documented source.
fn write_alignment(
    out: &mut String,
    reported: &[&rustyhmmer::pipeline::Hit],
    pli: &rustyhmmer::pli::Pipeline,
    model: &rustyhmmer::pipeline::Model,
    seqs: &[rustyhmmer::seqio::Seq],
    textw: i32,
) -> Option<usize> {
    use rustyhmmer::alidisplay::AliDisplay;

    let mut faux = Vec::new();
    let mut m = 0usize;
    for h in reported {
        if !pli.target_includable(h.score, h.lnp) {
            continue;
        }
        for d in &h.domains {
            if !pli.domain_includable(d.bitscore, d.lnp) {
                continue;
            }
            let Some(sq) = seqs.get(h.seq_idx) else { continue };
            let Some(ad) = AliDisplay::create(&d.trace, &model.hmm, &model.fwd, sq) else {
                continue;
            };
            if m == 0 {
                m = ad.m as usize;
            }
            faux.push(rustyhmmer::tracealign::backconvert(&ad));
        }
    }
    if faux.is_empty() {
        return None;
    }

    let mut msa = rustyhmmer::tracealign::seqs(&faux, m, model.hmm.mm.as_deref());
    msa.name = Some(model.hmm.name.clone());
    msa.acc = model.hmm.acc.clone();
    msa.desc = model.hmm.desc.clone();
    msa.au = Some("hmmsearch (HMMER 3.4)".to_string());

    let cpl = if textw > 0 {
        rustyhmmer::msa::STOCKHOLM_CPL
    } else {
        msa.alen
    };
    let n = msa.nseq();
    rustyhmmer::msa::write(out, &msa, cpl);
    Some(n)
}
