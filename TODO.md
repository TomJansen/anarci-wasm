# TODO: RustyHMMER hmmscan Workaround

The current Rust backend is intentionally RustyHMMER-only, but RustyHMMER does
not expose a direct in-memory `hmmscan` API that returns the same structured
domain alignment data ANARCI needs.

To match original ANARCI behavior, [src/rustyhmmer_engine.rs](src/rustyhmmer_engine.rs)
currently emulates `hmmscan` over the embedded `data/ALL.hmm` database:

- Parse the embedded HMMER3 profiles with `rustyhmmer::hmmfile::P7Hmm`.
- Convert every HMM into a RustyHMMER `Model`.
- For each query sequence, run `Model::search_one` against each model.
- Set the effective database size `Z` to the number of ANARCI HMMs so E-values
  line up with `hmmscan` over `ALL.hmm`.
- Ask RustyHMMER to render the domain report text.
- Parse the HMMER-style domain table and alignment blocks back into structured
  hit data.
- Convert the displayed model/query alignment columns into the match/insert/delete
  path used by ANARCI numbering.

This works and is covered by parity tests against the original `hmmscan` CLI,
but it is still a workaround. The fragile part is the
text parsing needed to recover alignment/path information because RustyHMMER's
public API does not currently expose the full structured domain trace we need.

## Desired Upstream Shape

If RustyHMMER grows a proper scan API, replace the workaround with something
closer to:

- Load HMMs from an embedded reader or byte slice.
- Run an in-memory `hmmscan` equivalent over all models.
- Pass explicit `Z`/`domZ` settings matching HMMER semantics.
- Receive structured hits, domains, coordinates, scores, E-values, and trace or
  alignment columns directly.
- Build ANARCI's match/insert/delete path without parsing rendered text.

## Keep Before Replacing

Any replacement must continue passing:

- `cargo test`
- `cargo test --test hmmer_parity -- --ignored --nocapture`
- `wasm-pack build --target web --release`

The important parity target is not just "finds the same top chain"; it must
match original ANARCI/HMMER on chain class, species, sequence coordinates, bit
score tolerance, E-value text, and final numbering positions.
