// SVG chromatogram renderer for .ab1 traces.
//
// Given a parsed trace (from ab1.js) and the numbering result for that sequence,
// renders an expandable panel showing:
//   - the 4-channel chromatogram (A/C/G/T) as SVG paths
//   - the base-call row aligned to peak locations
//   - the best-scoring frame's translation (one AA centred under each codon)
// Stop codons read through as `X` (from the Rust side) are highlighted in both
// the base and AA rows so the user can jump to the peak and judge the base call.

const CHANNEL_COLORS = { A: '#2e7d32', C: '#1565c0', G: '#000000', T: '#c62828' };
const TRACE_HEIGHT = 120;
const BASE_ROW_Y = TRACE_HEIGHT + 16;
const AA_ROW_Y = TRACE_HEIGHT + 46;
const PX_PER_SAMPLE = 2.2; // horizontal scale of the chromatogram
const SVG_NS = 'http://www.w3.org/2000/svg';

// Build the collapsible trace panel. `res` is the numbering result for this
// sequence; `escHtml` is passed in to avoid duplicating the app's escaper.
export function buildTracePanel(trace, res, escHtml) {
  const wrap = document.createElement('div');
  wrap.className = 'trace-panel';

  const btn = document.createElement('button');
  btn.type = 'button';
  btn.className = 'trace-toggle';
  btn.textContent = 'Show trace';
  btn.setAttribute('aria-expanded', 'false');

  const viewport = document.createElement('div');
  viewport.className = 'trace-viewport hidden';

  let scroller = null; // the rendered .trace-scroll element, once expanded

  // Render lazily on first expand — SVG for long reads is non-trivial.
  const ensureRendered = () => {
    if (scroller) return;
    scroller = renderTraceSvg(trace, res, escHtml);
    viewport.appendChild(scroller);
  };

  const setOpen = (open, scrollX) => {
    viewport.classList.toggle('hidden', !open);
    btn.textContent = open ? 'Hide trace' : 'Show trace';
    btn.setAttribute('aria-expanded', open ? 'true' : 'false');
    if (open) {
      ensureRendered();
      const targetX = scrollX != null ? scrollX : scroller._scrollToX;
      if (targetX != null) {
        requestAnimationFrame(() => { scroller.scrollLeft = Math.max(0, targetX); });
      }
    }
  };

  btn.addEventListener('click', () => setOpen(viewport.classList.contains('hidden')));

  wrap.appendChild(btn);
  wrap.appendChild(viewport);

  // Imperative API used by the numbering grid: open the trace (if needed),
  // centre the k-th annotated codon in the viewport, and briefly flash it.
  wrap.scrollToCodon = k => {
    ensureRendered();
    // The grid addresses codons by domain-residue index; the trace's codon
    // array has the front-flank codons prepended, so shift by that offset.
    const codonIdx = (scroller._domainCodonOffset || 0) + k;
    const cx = scroller._codonCenters && scroller._codonCenters[codonIdx];
    if (cx == null) return;
    setOpen(true);
    // Centre after layout is final, so clientWidth is correct (the trace may
    // have just been rendered/unhidden this frame).
    requestAnimationFrame(() => {
      scroller.scrollLeft = Math.max(0, cx - scroller.clientWidth / 2);
      flashCodon(scroller, codonIdx);
    });
  };

  return wrap;
}

// Briefly overlay a light-gray band on codon `k` to mark the click target.
function flashCodon(scroller, k) {
  const svg = scroller.querySelector('svg.trace-svg');
  const span = svg && svg._codonSpans && svg._codonSpans[k];
  if (!span) return;
  const band = document.createElementNS(SVG_NS, 'rect');
  band.setAttribute('x', span.left.toFixed(1));
  band.setAttribute('y', '0');
  band.setAttribute('width', (span.right - span.left).toFixed(1));
  band.setAttribute('height', svg.getAttribute('height'));
  band.setAttribute('class', 'trace-flash-band');
  // Insert behind the chromatogram/glyphs so it tints rather than covers them.
  svg.insertBefore(band, svg.firstChild);
  setTimeout(() => band.remove(), 1500);
}

function renderTraceSvg(trace, res, escHtml) {
  const { channels, peakLocations, traceLength, bases } = trace;
  const width = Math.max(traceLength * PX_PER_SAMPLE, 200);
  const height = AA_ROW_Y + 18;

  const svg = document.createElementNS(SVG_NS, 'svg');
  svg.setAttribute('class', 'trace-svg');
  svg.setAttribute('width', String(Math.round(width)));
  svg.setAttribute('height', String(height));
  svg.setAttribute('viewBox', `0 0 ${Math.round(width)} ${height}`);

  // Scale traces to fit TRACE_HEIGHT.
  let maxVal = 1;
  for (const base of Object.keys(channels)) {
    for (const v of channels[base]) if (v > maxVal) maxVal = v;
  }
  const yScale = (TRACE_HEIGHT - 4) / maxVal;
  const xOf = sampleIdx => sampleIdx * PX_PER_SAMPLE;
  const yOf = v => TRACE_HEIGHT - v * yScale;

  // One path per channel.
  for (const base of ['A', 'C', 'G', 'T']) {
    const data = channels[base];
    if (!data) continue;
    const path = document.createElementNS(SVG_NS, 'path');
    let d = '';
    // Subsample for very long traces to keep the path light.
    const step = data.length > 4000 ? 2 : 1;
    for (let i = 0; i < data.length; i += step) {
      d += `${i === 0 ? 'M' : 'L'}${xOf(i).toFixed(1)} ${yOf(data[i]).toFixed(1)} `;
    }
    path.setAttribute('d', d);
    path.setAttribute('fill', 'none');
    path.setAttribute('stroke', CHANNEL_COLORS[base]);
    path.setAttribute('stroke-width', '1');
    path.setAttribute('opacity', '0.85');
    svg.appendChild(path);
  }

  const domain = res.domains[0];
  const frame = domain ? domain.translation_frame : '';
  const { codons, firstCodon, firstDomainIdx } = mapTranslation(res, bases.length);

  const xOfBase = i =>
    peakLocations[i] != null ? xOf(peakLocations[i]) : (i / bases.length) * width;

  // Map each base-call index to the codon (k) it belongs to, so base glyphs can
  // be tagged and flashed alongside their amino acid.
  const codonOfBase = new Map();
  codons.forEach((codon, k) => {
    for (const idx of codon.baseIdxs) codonOfBase.set(idx, k);
  });

  // Centre x of each codon's AA glyph, indexed by k — used for click-to-scroll.
  const codonCenters = codons.map(codon => {
    const xs = codon.baseIdxs.map(xOfBase);
    return xs.reduce((a, b) => a + b, 0) / xs.length;
  });

  // Pixel span (left/right x) of each codon, used by the transient click flash.
  const codonSpans = codons.map(codon => {
    const xs = codon.baseIdxs.map(xOfBase);
    return { left: Math.min(...xs) - 6, right: Math.max(...xs) + 6 };
  });
  svg._codonSpans = codonSpans;

  // Codon grouping: a short bracket under each triplet that funnels into the AA
  // glyph, tying the three bases to their amino acid.
  codons.forEach((codon, k) => {
    const { left, right } = codonSpans[k];
    const cx = (left + right) / 2;
    const bracketY = BASE_ROW_Y + 5;
    const bracket = document.createElementNS(SVG_NS, 'path');
    bracket.setAttribute(
      'd',
      `M${(left + 2).toFixed(1)} ${bracketY} L${(right - 2).toFixed(1)} ${bracketY} ` +
        `M${cx.toFixed(1)} ${bracketY} L${cx.toFixed(1)} ${bracketY + 5}`
    );
    bracket.setAttribute(
      'class',
      `trace-codon-bracket${codon.isFlank ? ' trace-flank' : ''}`
    );
    svg.appendChild(bracket);
  });

  // Base-call row, aligned to peak locations.
  for (let i = 0; i < bases.length; i++) {
    const x = xOfBase(i);
    const base = bases[i];

    const t = document.createElementNS(SVG_NS, 'text');
    t.setAttribute('x', x.toFixed(1));
    t.setAttribute('y', String(BASE_ROW_Y));
    t.setAttribute('text-anchor', 'middle');
    t.setAttribute('class', 'trace-base');
    t.setAttribute('fill', CHANNEL_COLORS[base] || '#555');
    if (codonOfBase.has(i)) t.setAttribute('data-codon', String(codonOfBase.get(i)));
    t.textContent = base;
    svg.appendChild(t);
  }

  // AA row: one glyph per codon, centred between the codon's three bases.
  codons.forEach((codon, k) => {
    const cx = codonCenters[k];

    if (codon.isStop) {
      // Boxed glyph so X / stop read-throughs stand out from normal residues.
      const box = document.createElementNS(SVG_NS, 'rect');
      box.setAttribute('x', (cx - 7).toFixed(1));
      box.setAttribute('y', String(AA_ROW_Y - 11));
      box.setAttribute('width', '14');
      box.setAttribute('height', '15');
      box.setAttribute('rx', '2');
      box.setAttribute('class', 'trace-aa-stopbox');
      box.setAttribute('data-codon', String(k));
      svg.appendChild(box);
    }

    const a = document.createElementNS(SVG_NS, 'text');
    a.setAttribute('x', cx.toFixed(1));
    a.setAttribute('y', String(AA_ROW_Y));
    a.setAttribute('text-anchor', 'middle');
    a.setAttribute(
      'class',
      `trace-aa${codon.isStop ? ' trace-stop' : ''}${codon.isFlank ? ' trace-flank' : ''}`
    );
    a.setAttribute('data-codon', String(k));
    a.textContent = codon.residue;
    svg.appendChild(a);
  });

  const container = document.createElement('div');
  container.className = 'trace-scroll';
  const caption = document.createElement('div');
  caption.className = 'trace-caption';
  caption.innerHTML = domain
    ? `Frame <strong>${escHtml(frame || '?')}</strong> &middot; ${escHtml(domain.chain_type)} &middot; score ${domain.bit_score.toFixed(1)}`
    : 'No antibody domain detected — showing best-effort base calls';
  container.appendChild(caption);
  container.appendChild(svg);
  container._codonCenters = codonCenters;
  // Grid residue index k maps to codon (firstDomainIdx + k) here, since the
  // flank codons are prepended to the combined array.
  container._domainCodonOffset = firstDomainIdx;

  // Jump to the first annotated residue on open, rather than the read start
  // (which is often messy primer / leader sequence).
  if (firstCodon) {
    const firstX = firstCodon.baseIdxs.map(xOfBase).reduce((a, b) => a + b, 0) /
      firstCodon.baseIdxs.length;
    container._scrollToX = Math.max(0, firstX - 60);
  }
  return container;
}

// Map each base-call index to the amino acid of the best-scoring frame, and
// determine which base indices correspond to stop-codon read-throughs.
// nt coordinates come from the domain (nt_start/nt_end) and the frame label
// gives strand + offset. We align AA glyphs to the *middle* base of each codon.
function mapTranslation(res, numBases) {
  // One annotation per codon: the residue, whether it is a stop, and the three
  // base-call indices the codon spans (so the caller can centre the AA glyph
  // between them and locate the flash band).
  const codons = [];
  const domain = res.domains[0];
  if (!domain || domain.nt_start == null) {
    return { codons, firstCodon: null };
  }

  const frame = domain.translation_frame || '+1';
  const isReverse = frame.startsWith('-');
  const numbered = domain.numbering.filter(e => e.amino_acid !== '-');
  const stopSet = new Set(domain.stop_positions || []);
  const flankBefore = domain.flank_before || '';
  const flankAfter = domain.flank_after || '';

  // Codons run from nt_start, every 3 bases, across the domain span. For the
  // forward strand the base index is straightforward; for the reverse strand
  // the codon's bases are mirrored onto the displayed (forward) trace.
  const baseIndexOf = nt => (isReverse ? numBases - 1 - nt : nt);

  // Build one codon annotation for a residue whose codon starts at `codonStartNt`
  // (nt coordinate in the frame). Returns null if it falls off the read.
  const makeCodon = (residue, codonStartNt, isStop, isFlank) => {
    const baseIdxs = [0, 1, 2]
      .map(b => baseIndexOf(codonStartNt + b))
      .filter(idx => idx >= 0 && idx < numBases);
    if (baseIdxs.length === 0) return null;
    return { residue, isStop, isFlank, baseIdxs };
  };

  // Front flank: the N-terminal residues (from the start codon up to the domain)
  // sit immediately before nt_start, one codon per residue, rendered lighter.
  for (let i = 0; i < flankBefore.length; i++) {
    const codonStartNt = domain.nt_start - (flankBefore.length - i) * 3;
    const c = makeCodon(flankBefore[i], codonStartNt, false, true);
    if (c) codons.push(c);
  }

  const firstDomainIdx = codons.length; // first non-flank codon, for initial scroll
  for (let k = 0; k < numbered.length; k++) {
    const residue = numbered[k].amino_acid;
    const codonStartNt = domain.nt_start + k * 3;
    const isStop = stopSet.has(domain.seq_start + k) || residue === 'X';
    const c = makeCodon(residue, codonStartNt, isStop, false);
    if (c) codons.push(c);
  }

  // Back flank: the C-terminal residues after the domain, up to the first stop.
  const ntAfter = domain.nt_start + numbered.length * 3;
  for (let i = 0; i < flankAfter.length; i++) {
    const c = makeCodon(flankAfter[i], ntAfter + i * 3, false, true);
    if (c) codons.push(c);
  }

  return {
    codons,
    firstCodon: codons[firstDomainIdx] || codons[0] || null,
    // Index of the first domain codon in `codons` (flanks are prepended). The
    // numbering grid addresses codons by domain-residue index, so callers add
    // this offset to map a grid index onto the combined codon array.
    firstDomainIdx,
  };
}
