# ANARCI-WASM

Antibody Numbering and Receptor ClassIfication — compiled to WebAssembly.

A browser-based port of [ANARCI](https://github.com/oxpig/ANARCI) that runs entirely client-side. No server, no Python, no HMMER installation required.
The browser UI uses Web Workers to keep numbering responsive and to parallelize multi-sequence jobs across available CPU cores.

## Supported numbering schemes

- IMGT
- Chothia (IG only)
- Kabat (IG only)
- Martin / Enhanced Chothia (IG only)
- AHo
- Wolfguy (IG only)

## Prerequisites

- [Rust](https://rustup.rs/) (1.70+)
- [wasm-pack](https://rustwasm.github.io/wasm-pack/installer/)

```bash
# Install wasm-pack
curl https://rustwasm.github.io/wasm-pack/installer/init.sh -sSf | sh
```

## Build

```bash
cd anarci-wasm

# Build the WASM package
wasm-pack build --target web --release

# The output goes to pkg/
```

## Run locally

After building, serve the `www/` directory with any static file server:

```bash
# Option 1: Python
cd www && ln -sf ../pkg . && python3 -m http.server 8080

# Option 2: npx
cd www && ln -sf ../pkg . && npx serve .
```

Then open http://localhost:8080 in your browser.

## Usage

1. Open the web page
2. Paste FASTA sequences or drag-and-drop a `.fasta` file
3. Select a numbering scheme (IMGT, Chothia, Kabat, Martin, AHo, Wolfguy)
4. Click "Number sequences"
5. Results appear inline with CDR highlighting; download as CSV

## HMM data

The HMM profiles (`data/ALL.hmm`) are embedded into the WASM binary at compile time. These are the same HMMER3 profiles used by the original ANARCI, built from IMGT germline sequences.

To rebuild HMMs from scratch, use the `build_pipeline/` in the main ANARCI repository.

## API

The WASM module exports:

```javascript
// Number sequences from FASTA text, returns JSON
number_sequences(fasta_text: string, scheme: string, bit_score_threshold: number): string

// List available schemes
available_schemes(): string

// Number of loaded HMM profiles
num_profiles(): number
```

## License

GPL v3.0
