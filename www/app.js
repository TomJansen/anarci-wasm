import init, { number_sequences_with_options, available_schemes, available_species, compute_pim } from './pkg/anarci_wasm.js';
import { parseAb1 } from './ab1.js';
import { buildTracePanel } from './trace.js';

// Trace data extracted from dropped .ab1 files, keyed by FASTA sequence id.
window._ab1Traces = window._ab1Traces || new Map();

let wasmReady = false;
window._showUniqueOnly = true;
// Summary panels are collapsed by default.
window._duplicateMapCollapsed = true;
window._stopCodonMapCollapsed = true;
window._noDomainMapCollapsed = true;
const MAX_PARALLEL_WORKERS = 4;
let workerPoolState = null;
let activeRunToken = 0;

async function initWasm() {
  try {
    await init();
    populateSpeciesFilter();
    document.getElementById('status').textContent = '';
    document.getElementById('status').className = '';
    document.getElementById('runBtn').disabled = false;
    wasmReady = true;
  } catch(e) {
    document.getElementById('status').textContent = 'Failed to load WASM: ' + e.message;
    document.getElementById('status').className = 'error';
  }
}

function populateSpeciesFilter() {
  const speciesFilter = document.getElementById('speciesFilter');
  const species = JSON.parse(available_species());
  speciesFilter.innerHTML = '<option value="">All species</option>';
  for (const entry of species) {
    const option = document.createElement('option');
    option.value = entry;
    option.textContent = entry;
    speciesFilter.appendChild(option);
  }
}

// File handling
const fileDrop = document.getElementById('fileDrop');
const fileInput = document.getElementById('fileInput');
const fastaInput = document.getElementById('fastaInput');

fileDrop.addEventListener('click', () => fileInput.click());
fileDrop.addEventListener('dragover', e => { e.preventDefault(); fileDrop.classList.add('dragover'); });
fileDrop.addEventListener('dragleave', () => fileDrop.classList.remove('dragover'));
fileDrop.addEventListener('drop', e => {
  e.preventDefault();
  fileDrop.classList.remove('dragover');
  if (e.dataTransfer.files.length) readFiles(e.dataTransfer.files);
});
fileInput.addEventListener('change', e => { if (e.target.files.length) readFiles(e.target.files); });

async function readFiles(fileList) {
  const files = Array.from(fileList);
  const ab1Files = files.filter(f => /\.ab1$/i.test(f.name));
  const textFiles = files.filter(f => !/\.ab1$/i.test(f.name));

  const textRecords = await Promise.all(textFiles.map(f => f.text()));
  const ab1Records = await Promise.all(ab1Files.map(readAb1File));

  const fastaChunks = [
    ...textRecords.map(c => c.trim()).filter(Boolean),
    ...ab1Records.filter(Boolean).map(r => r.fasta),
  ];
  fastaInput.value = fastaChunks.join('\n\n');

  // If any .ab1 was loaded, switch to DNA mode so base calls are translated.
  if (ab1Records.some(Boolean)) {
    document.getElementById('inputType').value = 'dna';
  }
}

// Parse one .ab1 file: returns a FASTA record for its base calls and stashes
// the chromatogram trace keyed by the sequence id (the file's base name).
async function readAb1File(file) {
  const buffer = await file.arrayBuffer();
  const trace = parseAb1(buffer);
  if (!trace) {
    setStatus(`Could not parse ${escHtml(file.name)} as an .ab1 file.`, 'error');
    return null;
  }
  const id = file.name.replace(/\.ab1$/i, '');
  window._ab1Traces.set(id, trace);
  return { fasta: `>${id}\n${trace.bases}` };
}

function splitFastaRecords(fastaText) {
  const normalized = fastaText.replace(/\r\n?/g, '\n').trim();
  if (!normalized) return [];

  const records = [];
  let current = [];

  for (const line of normalized.split('\n')) {
    if (line.startsWith('>')) {
      if (current.length > 0) {
        records.push(current.join('\n'));
      }
      current = [line];
    } else if (current.length > 0) {
      current.push(line);
    }
  }

  if (current.length > 0) {
    records.push(current.join('\n'));
  }

  return records;
}

function getWorkerCount(taskCount) {
  const cores = navigator.hardwareConcurrency || 2;
  return Math.max(1, Math.min(MAX_PARALLEL_WORKERS, cores, taskCount));
}

function shutdownWorkerPool() {
  if (!workerPoolState) return;
  for (const worker of workerPoolState.workers) {
    worker.terminate();
  }
  workerPoolState = null;
}

function getWorkerPool(size) {
  if (workerPoolState && workerPoolState.size === size) {
    return workerPoolState;
  }

  shutdownWorkerPool();

  workerPoolState = {
    size,
    workers: Array.from({ length: size }, () => new Worker('./worker.js', { type: 'module' })),
  };
  return workerPoolState;
}

window.addEventListener('beforeunload', shutdownWorkerPool);

function numberSequencesOnMainThread({
  fasta,
  scheme,
  threshold,
  allowedSpeciesJson,
  restrict,
  inputType,
}) {
  const json = number_sequences_with_options(
    fasta,
    scheme,
    threshold,
    allowedSpeciesJson,
    restrict,
    inputType,
  );
  return JSON.parse(json);
}

async function numberSequencesInWorkers({
  fasta,
  scheme,
  threshold,
  allowedSpeciesJson,
  restrict,
  inputType,
  statusPrefix,
  runToken,
}) {
  const records = splitFastaRecords(fasta);
  if (records.length === 0 || typeof Worker !== 'function') {
    return {
      results: numberSequencesOnMainThread({
        fasta,
        scheme,
        threshold,
        allowedSpeciesJson,
        restrict,
        inputType,
      }),
      workerCount: 0,
    };
  }

  const workerCount = getWorkerCount(records.length);
  const pool = getWorkerPool(workerCount);
  const results = new Array(records.length).fill(null);
  let nextTaskIndex = 0;
  let completed = 0;

  const updateProgress = () => {
    if (runToken !== activeRunToken) return;
    setStatus(
      `<span class="loading"></span>${statusPrefix}Numbering ${records.length} ${pluralize(records.length, 'sequence')} using ${workerCount} ${pluralize(workerCount, 'worker')}... (${completed}/${records.length})`,
      ''
    );
  };

  updateProgress();

  return new Promise((resolve, reject) => {
    let settled = false;

    const cleanup = () => {
      for (const worker of pool.workers) {
        worker.onmessage = null;
        worker.onerror = null;
      }
    };

    const fail = error => {
      if (settled) return;
      settled = true;
      cleanup();
      shutdownWorkerPool();
      reject(error instanceof Error ? error : new Error(String(error)));
    };

    const maybeFinish = () => {
      if (completed !== records.length || settled) {
        return false;
      }
      settled = true;
      cleanup();
      resolve({
        results: results.filter(result => result !== null),
        workerCount,
      });
      return true;
    };

    const dispatch = worker => {
      if (settled) return;

      const taskId = nextTaskIndex;
      if (taskId >= records.length) {
        maybeFinish();
        return;
      }

      nextTaskIndex += 1;
      worker.postMessage({
        taskId,
        fasta: records[taskId],
        scheme,
        threshold,
        allowedSpeciesJson,
        restrict,
        inputType,
      });
    };

    for (const worker of pool.workers) {
      worker.onmessage = event => {
        if (settled) return;

        const message = event.data || {};
        if (!message.ok) {
          fail(message.error || 'Worker failed to number sequences.');
          return;
        }

        let parsed;
        try {
          parsed = JSON.parse(message.json);
        } catch (error) {
          fail(error);
          return;
        }

        const taskId = Number(message.taskId);
        results[taskId] = Array.isArray(parsed) ? (parsed[0] ?? null) : null;
        completed += 1;
        updateProgress();

        if (maybeFinish()) {
          return;
        }
        dispatch(worker);
      };

      worker.onerror = event => {
        fail(event.message || 'Worker crashed.');
      };

      dispatch(worker);
    }
  });
}

// Run numbering
document.getElementById('runBtn').addEventListener('click', runNumbering);
document.getElementById('resetBtn').addEventListener('click', resetPage);
document.getElementById('outputFormat').addEventListener('change', updateDownloadButton);
document.getElementById('uniqueOnly').addEventListener('change', toggleUniqueOnly);
// Omitting stops affects the download only — no re-render needed, but refresh
// the download button label/state to reflect the changed count.
document.getElementById('omitStopsDownload').addEventListener('change', updateDownloadButton);
document.getElementById('omitNoDomainDownload').addEventListener('change', updateDownloadButton);
window.addEventListener('resize', () => {
  if (window._lastResults) rerenderResults();
});

async function runNumbering() {
  if (!wasmReady) return;
  const fasta = fastaInput.value.trim();
  if (!fasta) { setStatus('Please enter FASTA sequences.', 'error'); return; }

  const scheme = document.getElementById('scheme').value;
  const threshold = parseFloat(document.getElementById('threshold').value) || 120;
  const inputTypeSelect = document.getElementById('inputType');
  let inputType = inputTypeSelect.value;
  if (inputType === 'protein' && fastaLooksLikeDnaOnly(fasta)) {
    inputType = 'dna';
    inputTypeSelect.value = 'dna';
  }
  const selectedSpecies = document.getElementById('speciesFilter').value;
  const allowedSpecies = selectedSpecies ? [selectedSpecies] : [];
  const allowedSpeciesJson = JSON.stringify(allowedSpecies);
  const restrict = document.getElementById('restrictFilter').value;

  const statusPrefix = inputTypeSelect.value === 'dna' && fastaLooksLikeDnaOnly(fasta)
    ? 'Detected DNA FASTA. '
    : '';
  setStatus(`<span class="loading"></span>${statusPrefix}Numbering sequences...`, '');
  document.getElementById('runBtn').disabled = true;
  document.getElementById('resetBtn').disabled = true;

  const runToken = ++activeRunToken;

  // Yield once so the loading state paints before any fallback main-thread work.
  await new Promise(resolve => setTimeout(resolve, 0));

  try {
    const t0 = performance.now();
    let results;
    let workerCount = 0;

    try {
      const parallel = await numberSequencesInWorkers({
        fasta,
        scheme,
        threshold,
        allowedSpeciesJson,
        restrict,
        inputType,
        statusPrefix,
        runToken,
      });
      results = parallel.results;
      workerCount = parallel.workerCount;
    } catch (error) {
      console.warn('Falling back to main-thread numbering:', error);
      setStatus(`<span class="loading"></span>${statusPrefix}Worker setup failed, retrying on the main thread...`, '');
      await new Promise(resolve => setTimeout(resolve, 0));
      results = numberSequencesOnMainThread({
        fasta,
        scheme,
        threshold,
        allowedSpeciesJson,
        restrict,
        inputType,
      });
    }

    if (runToken !== activeRunToken) {
      return;
    }

    const elapsed = ((performance.now() - t0) / 1000).toFixed(2);
    const nDomains = results.reduce((s, r) => s + r.domains.length, 0);
    const parallelNote = workerCount > 0
      ? ` using ${workerCount} ${pluralize(workerCount, 'worker')}`
      : '';

    setStatus(
      `Numbered ${nDomains} ${pluralize(nDomains, 'domain')} in ${results.length} ${pluralize(results.length, 'sequence')} in ${elapsed}s${parallelNote}.`,
      'success'
    );
    window._lastResults = results;
    window._lastScheme = scheme;
    rerenderResults();
  } catch(e) {
    if (runToken === activeRunToken) {
      setStatus('Error: ' + e.message, 'error');
    }
  } finally {
    if (runToken === activeRunToken) {
      document.getElementById('runBtn').disabled = false;
      document.getElementById('resetBtn').disabled = false;
    }
  }
}

function resetPage() {
  fastaInput.value = '';
  fileInput.value = '';
  document.getElementById('scheme').value = 'imgt';
  document.getElementById('threshold').value = '120';
  document.getElementById('inputType').value = 'protein';
  document.getElementById('speciesFilter').value = '';
  document.getElementById('restrictFilter').value = '';
  document.getElementById('outputFormat').value = 'csv';
  document.getElementById('correctResults').innerHTML = '';
  document.getElementById('summaryBar').innerHTML = '';
  document.getElementById('summaryBar').classList.add('hidden');
  document.getElementById('duplicateMap').innerHTML = '';
  document.getElementById('duplicateMap').classList.add('hidden');
  document.getElementById('stopCodonMap').innerHTML = '';
  document.getElementById('stopCodonMap').classList.add('hidden');
  document.getElementById('noDomainMap').innerHTML = '';
  document.getElementById('noDomainMap').classList.add('hidden');
  document.getElementById('downloadBtn').classList.add('hidden');
  document.getElementById('omitStopsDownload').checked = true;
  document.getElementById('omitNoDomainDownload').checked = true;
  window._lastResults = null;
  window._lastScheme = null;
  window._showUniqueOnly = true;
  window._duplicateMapCollapsed = true;
  window._stopCodonMapCollapsed = true;
  window._noDomainMapCollapsed = true;
  updateUniqueToggleButton();

  if (wasmReady) {
    setStatus('', '');
  } else {
    setStatus('Loading WASM module...', '');
  }
}

function setStatus(html, cls) {
  const el = document.getElementById('status');
  el.innerHTML = html;
  el.className = cls;
}

function pluralize(count, singular, plural = `${singular}s`) {
  return count === 1 ? singular : plural;
}

function fastaLooksLikeDnaOnly(fastaText) {
  const seq = fastaText
    .split('\n')
    .filter(line => line && !line.startsWith('>'))
    .join('')
    .replace(/\s+/g, '')
    .replace(/[-.]/g, '')
    .toUpperCase();

  if (seq.length < 30) return false;
  return /^[ACGTUNRYKMSWBDHV]+$/.test(seq);
}

// IMGT CDR definitions for highlighting
const CDR_RANGES = {
  imgt: { H: [[27,38],[56,65],[105,117]], L: [[27,38],[56,65],[105,117]] },
  chothia: { H: [[26,32],[52,56],[95,102]], L: [[24,34],[50,56],[89,97]] },
  kabat: { H: [[31,35],[50,65],[95,102]], L: [[24,34],[50,56],[89,97]] },
  martin: { H: [[26,32],[52,56],[95,102]], L: [[24,34],[50,56],[89,97]] },
  aho: { H: [[25,42],[58,77],[107,138]], L: [[25,42],[58,77],[107,138]] },
  wolfguy: { H: [[151,199],[251,299],[331,399]], L: [[551,599],[651,699],[751,799]] },
};

function cdrClassForPosition(pos, scheme, chainClass) {
  const ranges = (CDR_RANGES[scheme] || {})[chainClass] || [];
  for (let i = 0; i < ranges.length; i += 1) {
    const [start, end] = ranges[i];
    if (pos >= start && pos <= end) {
      return `cdr${i + 1}`;
    }
  }
  return '';
}

function resultLabel(res, index) {
  return `${res.id} [${index + 1}]`;
}

function domainResidueSequence(domain) {
  return (domain.numbering || [])
    .map(entry => entry.amino_acid)
    .filter(aa => aa && aa !== '-')
    .join('');
}

function duplicateKeyForResult(res) {
  if (res.input_type === 'dna' && Array.isArray(res.domains) && res.domains.length > 0) {
    return res.domains
      .map(domain => `${domain.chain_type}:${domainResidueSequence(domain)}`)
      .join('|')
      .toUpperCase();
  }

  return (res.sequence || '').toUpperCase();
}

// A result contains an internal stop codon if any of its domains read one
// through (as X). Only DNA/.ab1 inputs ever populate stop_positions.
function resultHasInternalStop(res) {
  return (res.domains || []).some(d => (d.stop_positions || []).length > 0);
}

function summarizeResults(results) {
  const seen = new Map();
  const uniqueResults = [];
  const duplicateMappings = [];

  results.forEach((res, index) => {
    const key = duplicateKeyForResult(res);
    if (!seen.has(key)) {
      seen.set(key, { canonicalIndex: index, canonical: res });
      uniqueResults.push(res);
      return;
    }

    const match = seen.get(key);
    duplicateMappings.push({
      duplicateIndex: index,
      duplicate: res,
      canonicalIndex: match.canonicalIndex,
      canonical: match.canonical,
    });
  });

  const stopResults = results.filter(resultHasInternalStop);
  const noDomainResults = results.filter(res => (res.domains || []).length === 0);

  // "Correct" = has a domain and no internal stop codon. These are the cards
  // shown in the Correct sequences section.
  const correctResults = results.filter(
    res => (res.domains || []).length > 0 && !resultHasInternalStop(res)
  );
  // "Usable" = correct AND unique — the deduplicated end set you take forward.
  const uniqueSet = new Set(uniqueResults);
  const usableResults = correctResults.filter(res => uniqueSet.has(res));

  return {
    uniqueResults,
    duplicateMappings,
    stopResults,
    noDomainResults,
    correctResults,
    usableResults,
  };
}

function groupDuplicateMappings(duplicateMappings) {
  const groups = new Map();

  duplicateMappings.forEach(mapping => {
    const key = mapping.canonicalIndex;
    if (!groups.has(key)) {
      groups.set(key, {
        canonicalIndex: mapping.canonicalIndex,
        canonical: mapping.canonical,
        duplicates: [],
      });
    }
    groups.get(key).duplicates.push(mapping);
  });

  return Array.from(groups.values()).sort((a, b) => a.canonicalIndex - b.canonicalIndex);
}

function getRenderedResults() {
  const results = window._lastResults || [];
  const summary = summarizeResults(results);
  return {
    results: window._showUniqueOnly ? summary.uniqueResults : results,
    summary,
  };
}

// Results to write to the downloaded file: the currently-shown set, minus
// sequences with internal stop codons when that option is checked. This is a
// download-only filter and does not affect what is rendered on the page.
function getResultsForDownload() {
  const { results } = getRenderedResults();
  const omitStops = document.getElementById('omitStopsDownload').checked;
  const omitNoDomain = document.getElementById('omitNoDomainDownload').checked;
  return results.filter(res =>
    (!omitStops || !resultHasInternalStop(res)) &&
    (!omitNoDomain || (res.domains || []).length > 0)
  );
}

function updateUniqueToggleButton() {
  // Keep the checkbox in sync with state (e.g. after a reset).
  const box = document.getElementById('uniqueOnly');
  if (box) box.checked = window._showUniqueOnly;
}

function updateDuplicateMapToggle(container) {
  const body = container.querySelector('.duplicate-map-body');
  const toggle = container.querySelector('.duplicate-map-toggle');
  if (!body || !toggle) return;

  body.classList.toggle('hidden', window._duplicateMapCollapsed);
  toggle.setAttribute('aria-expanded', String(!window._duplicateMapCollapsed));
  toggle.textContent = window._duplicateMapCollapsed ? 'Expand' : 'Collapse';
}

function toggleUniqueOnly() {
  window._showUniqueOnly = document.getElementById('uniqueOnly').checked;
  if (!window._lastResults) return;
  rerenderResults();
}

// Top-of-results summary. Leads with the count that actually matters — the
// usable set (correct + unique) — then a logical breakdown of where the rest
// of the input went.
function renderSummaryBar(summary) {
  const container = document.getElementById('summaryBar');
  const total = (window._lastResults || []).length;
  if (total === 0) {
    container.innerHTML = '';
    container.classList.add('hidden');
    return;
  }

  const usable = summary.usableResults.length;
  const correct = summary.correctResults.length;
  const unique = summary.uniqueResults.length;
  const duplicate = summary.duplicateMappings.length;
  const stop = summary.stopResults.length;
  const noDomain = summary.noDomainResults.length;

  const breakdown = [
    ['Total input', total],
    ['Correct', correct],
    ['Unique', unique],
    ['Duplicates', duplicate],
    ['With stop codon', stop],
    ['No domain', noDomain],
  ]
    .map(([label, value]) => `<span class="dup-stat"><strong>${label}:</strong> ${value}</span>`)
    .join('');

  container.innerHTML = `
    <div class="summary-headline">
      Usable sequences <span class="summary-sub">(correct &amp; unique)</span>:
      <span class="summary-number">${usable}</span>
      <span class="summary-sub">of ${total}</span>
    </div>
    <div class="dup-stats">${breakdown}</div>
  `;
  container.classList.remove('hidden');
}

function renderDuplicateMap(summary) {
  const container = document.getElementById('duplicateMap');
  const { uniqueResults, duplicateMappings } = summary;

  if (!window._lastResults || window._lastResults.length === 0) {
    container.innerHTML = '';
    container.classList.add('hidden');
    return;
  }

  const duplicateCount = duplicateMappings.length;
  const groups = groupDuplicateMappings(duplicateMappings);

  const rows = groups
    .map((group, groupIndex) => {
      const detailId = `duplicate-detail-${groupIndex}`;
      const duplicateCount = group.duplicates.length;
      const amount = duplicateCount + 1;
      const detailItems = group.duplicates
        .map(
          mapping =>
            `<li>${escHtml(resultLabel(mapping.duplicate, mapping.duplicateIndex))}</li>`
        )
        .join('');

      return `
        <tr>
          <td><strong>${escHtml(resultLabel(group.canonical, group.canonicalIndex))}</strong></td>
          <td>
            <button
              type="button"
              class="duplicate-expand-btn"
              aria-expanded="false"
              aria-controls="${detailId}"
              data-detail-id="${detailId}"
            >
              Show ${duplicateCount} ${pluralize(duplicateCount, 'sequence')}
            </button>
          </td>
          <td>${amount}</td>
        </tr>
        <tr id="${detailId}" class="duplicate-detail-row hidden">
          <td colspan="3" class="duplicate-detail-cell">
            <ul class="duplicate-detail-list">${detailItems}</ul>
          </td>
        </tr>
      `;
    })
    .join('');

  const collapsed = window._duplicateMapCollapsed;
  container.innerHTML = `
    <div class="duplicate-map-header" role="button" tabindex="0" aria-controls="duplicateMapBody">
      <h2>Duplicate Sequence Map</h2>
      <button type="button" class="duplicate-map-toggle" aria-expanded="${!collapsed}">${collapsed ? 'Expand' : 'Collapse'}</button>
    </div>
    <div id="duplicateMapBody" class="duplicate-map-body${collapsed ? ' hidden' : ''}">
      ${duplicateMappings.length === 0
        ? '<p class="dup-note">No duplicate sequences.</p>'
        : `<table class="duplicate-map-table">
        <thead>
          <tr>
            <th>Unique</th>
            <th>Non unique</th>
            <th>Amount</th>
          </tr>
        </thead>
        <tbody>${rows}</tbody>
      </table>`}
    </div>
  `;
  container.classList.remove('hidden');
  updateDuplicateMapToggle(container);

  const header = container.querySelector('.duplicate-map-header');
  const toggleMap = () => {
    window._duplicateMapCollapsed = !window._duplicateMapCollapsed;
    updateDuplicateMapToggle(container);
  };
  header.addEventListener('click', event => {
    if (event.target.closest('.duplicate-expand-btn')) return;
    if (event.target.closest('.duplicate-map-toggle')) return;
    toggleMap();
  });
  header.addEventListener('keydown', event => {
    if (event.key === 'Enter' || event.key === ' ') {
      event.preventDefault();
      toggleMap();
    }
  });

  container.querySelector('.duplicate-map-toggle').addEventListener('click', event => {
    event.stopPropagation();
    toggleMap();
  });

  container.querySelectorAll('.duplicate-expand-btn').forEach(btn => {
    btn.addEventListener('click', () => {
      const detailId = btn.dataset.detailId;
      const row = container.querySelector(`#${detailId}`);
      if (!row) return;
      const isHidden = row.classList.toggle('hidden');
      btn.setAttribute('aria-expanded', String(!isHidden));
      btn.textContent = isHidden
        ? btn.textContent.replace(/^Hide/, 'Show')
        : btn.textContent.replace(/^Show/, 'Hide');
    });
  });
}

function rerenderResults() {
  if (!window._lastResults) {
    document.getElementById('correctResults').innerHTML = '';
    document.getElementById('summaryBar').innerHTML = '';
    document.getElementById('summaryBar').classList.add('hidden');
    document.getElementById('duplicateMap').innerHTML = '';
    document.getElementById('duplicateMap').classList.add('hidden');
    document.getElementById('stopCodonMap').innerHTML = '';
    document.getElementById('stopCodonMap').classList.add('hidden');
    document.getElementById('noDomainMap').innerHTML = '';
    document.getElementById('noDomainMap').classList.add('hidden');
    updateUniqueToggleButton();
    updateDownloadButton();
    return;
  }

  const scheme = window._lastScheme || document.getElementById('scheme').value;
  const summary = summarizeResults(window._lastResults);

  renderSummaryBar(summary);
  renderDuplicateMap(summary);
  renderStopCodonMap(summary, scheme);
  renderNoDomainMap(summary, scheme);
  renderResults(applyUniqueFilter(summary.correctResults, summary), scheme);
  updateUniqueToggleButton();
  updateDownloadButton();
}

// Generic collapsible panel that renders full seq-result cards. Used for the
// stop-codon and no-domain sections. `collapseKey` is the window state flag.
function renderCardPanel({ containerId, title, note, results, scheme, collapseKey }) {
  const container = document.getElementById(containerId);
  if (results.length === 0) {
    container.innerHTML = '';
    container.classList.add('hidden');
    return;
  }

  const collapsed = window[collapseKey];
  container.innerHTML = `
    <div class="duplicate-map-header" role="button" tabindex="0">
      <h2>${escHtml(title)} (${results.length})</h2>
      <button type="button" class="duplicate-map-toggle" aria-expanded="${!collapsed}">${collapsed ? 'Expand' : 'Collapse'}</button>
    </div>
    <div class="duplicate-map-body${collapsed ? ' hidden' : ''}">
      ${note ? `<p class="dup-note">${escHtml(note)}</p>` : ''}
    </div>
  `;

  const body = container.querySelector('.duplicate-map-body');
  for (const res of results) {
    body.appendChild(buildSeqResultCard(res, scheme));
  }
  container.classList.remove('hidden');

  const header = container.querySelector('.duplicate-map-header');
  const toggle = container.querySelector('.duplicate-map-toggle');
  const toggleMap = () => {
    window[collapseKey] = !window[collapseKey];
    body.classList.toggle('hidden', window[collapseKey]);
    toggle.setAttribute('aria-expanded', String(!window[collapseKey]));
    toggle.textContent = window[collapseKey] ? 'Expand' : 'Collapse';
  };
  header.addEventListener('click', event => {
    if (event.target.closest('.duplicate-map-toggle')) return;
    toggleMap();
  });
  header.addEventListener('keydown', event => {
    if (event.key === 'Enter' || event.key === ' ') {
      event.preventDefault();
      toggleMap();
    }
  });
  toggle.addEventListener('click', event => {
    event.stopPropagation();
    toggleMap();
  });
}

// Restrict a result list to the deduplicated set when "show unique only" is on.
function applyUniqueFilter(results, summary) {
  if (!window._showUniqueOnly) return results;
  const uniqueSet = new Set(summary.uniqueResults);
  return results.filter(res => uniqueSet.has(res));
}

function renderStopCodonMap(summary, scheme) {
  renderCardPanel({
    containerId: 'stopCodonMap',
    title: 'Sequences with internal stop codons',
    note: 'These are excluded from the download when "Omit sequences with internal stop codons" is checked. Each card shows the stop codon (e.g. TGA) at its numbered position.',
    results: applyUniqueFilter(summary.stopResults || [], summary),
    scheme,
    collapseKey: '_stopCodonMapCollapsed',
  });
}

function renderNoDomainMap(summary, scheme) {
  renderCardPanel({
    containerId: 'noDomainMap',
    title: 'No domain found',
    note: 'No antibody/receptor domain was detected in these sequences.',
    results: applyUniqueFilter(summary.noDomainResults || [], summary),
    scheme,
    collapseKey: '_noDomainMapCollapsed',
  });
}

// Render the "Correct sequences" section: a labelled heading plus a card for
// each correct (has-domain, stop-free) result, honouring the unique-only toggle.
function renderResults(results, scheme) {
  const container = document.getElementById('correctResults');
  container.innerHTML = '';

  const heading = document.createElement('h2');
  heading.className = 'results-section-heading';
  heading.textContent = `Correct sequences (${results.length})`;
  container.appendChild(heading);

  if (results.length === 0) {
    const empty = document.createElement('p');
    empty.className = 'dup-note';
    empty.textContent = 'No correct sequences.';
    container.appendChild(empty);
    return;
  }

  for (const res of results) {
    container.appendChild(buildSeqResultCard(res, scheme));
  }
}

// Build one collapsible seq-result card (header + domains + optional trace).
function buildSeqResultCard(res, scheme) {
    const div = document.createElement('div');
    div.className = 'seq-result';

    const nDom = res.domains.length;
    const header = document.createElement('div');
    header.className = 'seq-header';
    header.innerHTML = `<span>${escHtml(res.id)}</span><span>${nDom ? `${nDom} ${pluralize(nDom, 'domain')}` : 'No domains found'}</span>`;

    const body = document.createElement('div');

    const domTrace = window._ab1Traces.get(res.id);

    for (const dom of res.domains) {
      // Domain info bar — ANARCI scores, then (separated) the .ab1 read quality.
      const info = document.createElement('div');
      info.className = 'domain-info';
      const badge = `chain-${dom.chain_class}`;
      // Phred quality summary over the domain's base calls, when this came from
      // an .ab1 with a quality track. Shown apart from the ANARCI scores.
      const qstats = dom === res.domains[0] ? domainQualityStats(dom, domTrace) : null;
      info.innerHTML = `
        <span class="chain-badge ${badge}">${dom.chain_type}</span>
        <span>Species: ${dom.species}</span>
        <span>Score: ${dom.bit_score.toFixed(1)}</span>
        <span>E-value: ${dom.evalue_text || dom.evalue.toExponential(1)}</span>
        <span>${dom.translation_frame ? 'AA region' : 'Region'}: ${dom.seq_start+1}-${dom.seq_end}</span>
        ${dom.translation_frame ? `<span>Frame: ${dom.translation_frame}</span>` : ''}
        ${dom.nt_start != null && dom.nt_end != null ? `<span>NT region: ${dom.nt_start + 1}-${dom.nt_end}</span>` : ''}
        ${qstats != null ? `<span class="domain-quality" title="Phred quality over this domain, averaged in error-probability space (the standard 'mean quality'). Q20% = fraction of bases at ≤1% error.">Quality: <strong>Q${qstats.meanQ}</strong> &middot; ${qstats.q20Pct}% &ge;Q20</span>` : ''}
      `;
      body.appendChild(info);

      // Numbering table (chunked into rows of 40)
      const tableDiv = document.createElement('div');
      tableDiv.className = 'numbering-table';

      const residues = dom.numbering.filter(e => e.amino_acid !== '-');
      const grid = document.createElement('div');
      grid.className = 'numbering-grid';

      // The k-th non-gap residue maps to the k-th codon in the trace, but only
      // for the best-scoring (first) domain, which is what the trace renders.
      const traceLinked = window._ab1Traces.has(res.id) && dom === res.domains[0];

      residues.forEach((e, k) => {
        const cell = document.createElement('div');
        const cdr = cdrClassForPosition(e.position, scheme, dom.chain_class);
        const cdrClass = cdr ? ` ${cdr}` : '';
        const gap = e.amino_acid === '-' ? ' gap' : '';
        // Stop read-throughs (X) get the boxed-red treatment, matching the trace.
        const isStop = e.amino_acid === 'X';
        const stop = isStop ? ' residue-stop' : '';
        const linkClass = traceLinked ? ' trace-linked' : '';
        cell.className = `residue-cell${cdrClass}${gap}${stop}${linkClass}`;
        if (traceLinked) cell.dataset.codon = String(k);

        // For a stop, show the actual codon (TGA/TAA/TAG) in place of the
        // position number, so the user sees which stop it is. Only the
        // best-scoring domain maps cleanly onto the trace's codons.
        const codon = isStop && dom === res.domains[0]
          ? codonNucleotides(dom, domTrace, k)
          : null;
        const label = codon
          ? `<span class="residue-label residue-codon" title="Stop codon at position ${escHtml(`${e.position}${e.insertion.trim() || ''}`)}">${escHtml(codon)}</span>`
          : `<span class="residue-label">${escHtml(`${e.position}${e.insertion.trim() || ''}`)}</span>`;
        cell.innerHTML = `
          <span class="residue-aa">${escHtml(e.amino_acid)}</span>
          ${label}
        `;
        grid.appendChild(cell);
      });

      tableDiv.appendChild(grid);
      body.appendChild(tableDiv);
    }

    // Chromatogram trace panel, if this sequence came from an .ab1 file.
    const trace = window._ab1Traces.get(res.id);
    if (trace) {
      const panel = buildTracePanel(trace, res, escHtml);
      body.appendChild(panel);
      // Clicking a residue cell of the best-scoring domain scrolls the trace to
      // that codon. Delegated so it survives the grid being (re)built.
      body.addEventListener('click', event => {
        const cell = event.target.closest('.residue-cell.trace-linked');
        if (!cell || !body.contains(cell)) return;
        const k = Number(cell.dataset.codon);
        if (Number.isInteger(k)) panel.scrollToCodon(k);
      });
    }

    // Toggle collapse
    let collapsed = false;
    header.addEventListener('click', () => {
      collapsed = !collapsed;
      body.style.display = collapsed ? 'none' : '';
    });

    div.appendChild(header);
    div.appendChild(body);
    return div;
}

function escHtml(s) {
  return s.replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;');
}

// Phred quality summary over the base calls spanned by a domain, or null if
// there is no trace / quality track / nt mapping. Quality is indexed like the
// forward base-call string; a reverse-strand frame maps its nt span onto
// mirrored indices (same convention as the trace renderer).
//
// Phred is a log scale (Q = -10 log10 P_error), so a plain arithmetic mean of
// Q-values overweights the good bases and understates the real error rate. We
// instead average in *probability* space — mean error probability, converted
// back to Q — which is the standard way QC tools report an "average quality".
// We also return the Q20 fraction (bases at <=1% error), the other headline
// QC metric (FastQC/Illumina style).
function domainQualityStats(dom, trace) {
  if (!trace || !trace.quality || trace.quality.length === 0) return null;
  if (dom.nt_start == null || dom.nt_end == null) return null;
  const n = trace.quality.length;
  const isReverse = (dom.translation_frame || '').startsWith('-');
  let errSum = 0;
  let count = 0;
  let q20 = 0;
  for (let nt = dom.nt_start; nt < dom.nt_end; nt++) {
    const idx = isReverse ? n - 1 - nt : nt;
    if (idx >= 0 && idx < n) {
      const q = trace.quality[idx];
      errSum += Math.pow(10, -q / 10);
      if (q >= 20) q20++;
      count++;
    }
  }
  if (!count) return null;
  const meanErr = errSum / count;
  // Mean error of 0 (all extremely high Q) → clamp to a sane ceiling.
  const meanQ = meanErr > 0 ? Math.round(-10 * Math.log10(meanErr)) : 60;
  return { meanQ, q20Pct: Math.round((q20 / count) * 100) };
}

const COMPLEMENT = { A: 'T', T: 'A', G: 'C', C: 'G', U: 'A', N: 'N' };

// The coding-strand nucleotide triplet for the k-th codon of a domain, so a
// stop reads as TGA/TAA/TAG. Returns null without a trace / nt mapping. For a
// reverse-strand frame the codon's bases are mirrored on the displayed trace,
// so we read them right-to-left and complement to recover the coding triplet.
function codonNucleotides(dom, trace, k) {
  if (!trace || !trace.bases || dom.nt_start == null) return null;
  const bases = trace.bases;
  const n = bases.length;
  const isReverse = (dom.translation_frame || '').startsWith('-');
  const codonStartNt = dom.nt_start + k * 3;
  let out = '';
  for (let b = 0; b < 3; b++) {
    const idx = isReverse ? n - 1 - (codonStartNt + b) : codonStartNt + b;
    if (idx < 0 || idx >= n) return null;
    const base = bases[idx];
    out += isReverse ? COMPLEMENT[base] || 'N' : base;
  }
  return out;
}

function csvEscape(value) {
  const s = String(value ?? '');
  if (/[",\n]/.test(s)) {
    return `"${s.replace(/"/g, '""')}"`;
  }
  return s;
}

function positionKey(position, insertion) {
  return `${position}|${insertion}`;
}

function positionLabel(position, insertion) {
  return `${position}${(insertion || '').trim()}`;
}

function collectPythonStylePositions(results) {
  const posRanks = new Map();
  const allPositions = new Map();

  for (const res of results) {
    for (const dom of res.domains) {
      let lastPos = null;
      let rank = 0;

      for (const entry of dom.numbering) {
        const pos = entry.position;
        const ins = entry.insertion || ' ';
        if (pos !== lastPos) {
          lastPos = pos;
          rank = 0;
        } else {
          rank += 1;
        }

        const key = positionKey(pos, ins);
        posRanks.set(key, Math.max(rank, posRanks.get(key) ?? rank));
        allPositions.set(key, { position: pos, insertion: ins });
      }
    }
  }

  return Array.from(allPositions.values()).sort((a, b) => {
    if (a.position !== b.position) return a.position - b.position;
    const rankA = posRanks.get(positionKey(a.position, a.insertion)) ?? 0;
    const rankB = posRanks.get(positionKey(b.position, b.insertion)) ?? 0;
    return rankA - rankB;
  });
}

function buildPythonStyleCsv(results) {
  const positions = collectPythonStylePositions(results);
  const fields = [
    'Id',
    'domain_no',
    'hmm_species',
    'chain_type',
    'e-value',
    'score',
    'seqstart_index',
    'seqend_index',
    'identity_species',
    'v_gene',
    'v_identity',
    'j_gene',
    'j_identity',
    ...positions.map(p => positionLabel(p.position, p.insertion)),
  ];

  const lines = [fields.map(csvEscape).join(',')];

  for (const res of results) {
    for (const dom of res.domains) {
      const numbered = new Map(
        dom.numbering.map(e => [positionKey(e.position, e.insertion || ' '), e.amino_acid])
      );

      const row = [
        res.id,
        dom.domain_index,
        dom.species,
        dom.chain_type,
        dom.evalue_text || dom.evalue,
        dom.bit_score.toFixed(1),
        dom.seq_start,
        dom.seq_end > dom.seq_start ? dom.seq_end - 1 : dom.seq_end,
        dom.identity_species ?? '',
        dom.v_gene ?? '',
        dom.v_identity ? dom.v_identity.toFixed(2) : '',
        dom.j_gene ?? '',
        dom.j_identity ? dom.j_identity.toFixed(2) : '',
        ...positions.map(p => numbered.get(positionKey(p.position, p.insertion)) ?? '-'),
      ];

      lines.push(row.map(csvEscape).join(','));
    }
  }

  return lines.join('\n') + '\n';
}

function buildAnarciText(results) {
  const lines = [];

  for (const res of results) {
    lines.push(`# ${res.id}`);

    if (res.domains.length > 0) {
      lines.push('# ANARCI numbered');
    }

    for (const [index, dom] of res.domains.entries()) {
      const seqEnd = dom.seq_end > dom.seq_start ? dom.seq_end - 1 : dom.seq_end;
      lines.push(`# Domain ${index + 1} of ${res.domains.length}`);
      lines.push('# Most significant HMM hit');
      lines.push('#|species|chain_type|e-value|score|seqstart_index|seqend_index|');
      lines.push(
        `#|${dom.species}|${dom.chain_type}|${dom.evalue_text || dom.evalue}|${dom.bit_score.toFixed(1)}|${dom.seq_start}|${seqEnd}|`
      );

      if (dom.identity_species || dom.v_gene || dom.j_gene || dom.v_identity || dom.j_identity) {
        lines.push('# Most sequence-identical germlines');
        lines.push('#|species|v_gene|v_identity|j_gene|j_identity|');
        lines.push(
          `#|${dom.identity_species || ''}|${dom.v_gene || 'unknown'}|${(dom.v_identity || 0).toFixed(2)}|${dom.j_gene || 'unknown'}|${(dom.j_identity || 0).toFixed(2)}|`
        );
      }

      lines.push(`# Scheme = ${dom.scheme}`);
      if (!dom.numbering.length) {
        lines.push(`# Warning: ${dom.scheme} scheme could not be applied to this sequence.`);
      }

      for (const entry of dom.numbering) {
        const insertion = entry.insertion || ' ';
        lines.push(`${dom.chain_class} ${String(entry.position).padEnd(5)} ${insertion} ${entry.amino_acid}`);
      }
    }

    lines.push('//');
  }

  return lines.join('\n') + '\n';
}

// Wrap a sequence to fixed-width lines for FASTA output.
function fastaWrap(seq, width = 60) {
  const lines = [];
  for (let i = 0; i < seq.length; i += width) {
    lines.push(seq.slice(i, i + width));
  }
  return lines.join('\n');
}

// The FASTA header label for one domain. Suffixes the chain type, and the
// domain index when a sequence has more than one domain, to keep ids unique.
function fastaDomainLabel(res, dom, index, domainCount) {
  const suffix = domainCount > 1 ? `_d${index + 1}` : '';
  const chain = dom.chain_type ? `_${dom.chain_type}` : '';
  return `${res.id}${chain}${suffix}`;
}

// One FASTA record per domain, using the numbered amino-acid sequence (the
// non-gap residues, including any X stop read-throughs). Works for both protein
// and DNA inputs since both populate domain.numbering.
function buildFastaAminoAcids(results) {
  const records = [];
  for (const res of results) {
    const domains = res.domains || [];
    for (const [index, dom] of domains.entries()) {
      const seq = domainResidueSequence(dom);
      if (!seq) continue;
      records.push(`>${fastaDomainLabel(res, dom, index, domains.length)}\n${fastaWrap(seq)}`);
    }
  }
  return records.join('\n') + (records.length ? '\n' : '');
}

// The coding-strand nucleotide span for a domain, recovered from the input DNA.
// For a reverse-strand frame the span is read right-to-left and complemented,
// matching the convention used by the trace / codon helpers.
function domainCodingDna(res, dom) {
  if (dom.nt_start == null || dom.nt_end == null) return '';
  const dna = (res.sequence || '').toUpperCase();
  const n = dna.length;
  const isReverse = (dom.translation_frame || '').startsWith('-');
  let out = '';
  for (let nt = dom.nt_start; nt < dom.nt_end; nt++) {
    const idx = isReverse ? n - 1 - nt : nt;
    if (idx < 0 || idx >= n) continue;
    const base = dna[idx];
    out += isReverse ? COMPLEMENT[base] || 'N' : base;
  }
  return out;
}

// One FASTA record per domain, using the coding DNA span. Only DNA inputs carry
// nucleotide coordinates; protein-input domains have no DNA and are skipped.
function buildFastaDna(results) {
  const records = [];
  for (const res of results) {
    if (res.input_type !== 'dna') continue;
    const domains = res.domains || [];
    for (const [index, dom] of domains.entries()) {
      const seq = domainCodingDna(res, dom);
      if (!seq) continue;
      records.push(`>${fastaDomainLabel(res, dom, index, domains.length)}\n${fastaWrap(seq)}`);
    }
  }
  return records.join('\n') + (records.length ? '\n' : '');
}

function updateDownloadButton() {
  const btn = document.getElementById('downloadBtn');
  if (!window._lastResults) {
    btn.classList.add('hidden');
    return;
  }

  const outputFormat = document.getElementById('outputFormat').value;
  btn.classList.remove('hidden');
  const base =
    outputFormat === 'anarci'
      ? 'Download ANARCI text'
      : outputFormat === 'fasta-aa'
      ? 'Download FASTA (AA)'
      : outputFormat === 'fasta-dna'
      ? 'Download FASTA (DNA)'
      : outputFormat.startsWith('pim-')
      ? 'Download PIM'
      : 'Download CSV';

  // Show how much the download will contain, so the omit filters' effect is
  // visible before clicking. FASTA exports emit one record per domain (DNA only
  // for the nucleotide variant); other formats emit one row per sequence.
  const results = getResultsForDownload();
  let count;
  if (outputFormat === 'fasta-aa') {
    count = results.reduce((n, res) => n + (res.domains || []).length, 0);
  } else if (outputFormat === 'fasta-dna') {
    count = results.reduce(
      (n, res) => n + (res.input_type === 'dna' ? (res.domains || []).length : 0),
      0
    );
  } else {
    count = results.length;
  }
  btn.textContent = `${base} (${count})`;
}

function hierarchicalClusterOrder(n, distFn) {
  if (n <= 1) return Array.from({length: n}, (_, i) => i);

  const dist = Array.from({length: n}, (_, i) =>
    Array.from({length: n}, (_, j) => i === j ? 0 : distFn(i, j))
  );
  const sizes = new Array(n).fill(1);
  // Each entry is either a leaf index (number) or a nested [left, right] pair
  const subtrees = Array.from({length: n}, (_, i) => i);
  const active = new Set(Array.from({length: n}, (_, i) => i));

  while (active.size > 1) {
    const activeArr = [...active];
    let minDist = Infinity, minI = -1, minJ = -1;
    for (let ai = 0; ai < activeArr.length; ai++) {
      for (let aj = ai + 1; aj < activeArr.length; aj++) {
        const i = activeArr[ai], j = activeArr[aj];
        if (dist[i][j] < minDist) { minDist = dist[i][j]; minI = i; minJ = j; }
      }
    }
    const sI = sizes[minI], sJ = sizes[minJ];
    subtrees[minI] = [subtrees[minI], subtrees[minJ]];
    for (const k of active) {
      if (k === minI || k === minJ) continue;
      dist[minI][k] = dist[k][minI] =
        (dist[minI][k] * sI + dist[minJ][k] * sJ) / (sI + sJ);
    }
    sizes[minI] = sI + sJ;
    active.delete(minJ);
  }

  const order = [];
  function traverse(node) {
    if (typeof node === 'number') order.push(node);
    else { traverse(node[0]); traverse(node[1]); }
  }
  traverse(subtrees[[...active][0]]);
  return order;
}

function extractDomainPositionMap(domain, posFilter) {
  const map = new Map();
  for (const entry of domain.numbering) {
    if (entry.amino_acid === '-') continue;
    if (posFilter !== null && !posFilter.has(entry.position)) continue;
    map.set(positionKey(entry.position, entry.insertion || ' '), entry.amino_acid);
  }
  return map;
}

function computePairwiseIdentity(mapA, mapB) {
  const allKeys = new Set([...mapA.keys(), ...mapB.keys()]);
  if (allKeys.size === 0) return null;
  let matches = 0;
  for (const key of allKeys) {
    const a = mapA.get(key);
    const b = mapB.get(key);
    if (a && b && a === b) matches++;
  }
  return (matches / allKeys.size) * 100;
}

function buildPimCsv(results, scope, scheme) {
  const byChain = new Map();

  for (const res of results) {
    const chainSeqCount = new Map();
    for (const dom of res.domains) {
      const cc = dom.chain_class;
      const idx = chainSeqCount.get(cc) || 0;
      chainSeqCount.set(cc, idx + 1);
      const label = idx > 0 ? `${res.id}_${cc}${idx + 1}` : res.id;

      // K (kappa) and L (lambda) share the same CDR position ranges
      const cdrCC = (cc === 'K' || cc === 'L') ? 'L' : cc;
      const ranges = (CDR_RANGES[scheme] || {})[cdrCC] || [];

      let posFilter = null;
      if (scope === 'cdr1' || scope === 'cdr2' || scope === 'cdr3') {
        const cdrIdx = scope === 'cdr1' ? 0 : scope === 'cdr2' ? 1 : 2;
        posFilter = new Set();
        if (cdrIdx < ranges.length) {
          const [start, end] = ranges[cdrIdx];
          for (let p = start; p <= end; p++) posFilter.add(p);
        }
      } else if (scope === 'cdrs') {
        posFilter = new Set();
        for (const [start, end] of ranges) {
          for (let p = start; p <= end; p++) posFilter.add(p);
        }
      }

      const posMap = extractDomainPositionMap(dom, posFilter);
      if (!byChain.has(cc)) byChain.set(cc, []);
      byChain.get(cc).push({ label, posMap });
    }
  }

  if (byChain.size === 0) return '';

  const scopeLabel = scope === 'full' ? 'Full sequence'
    : scope === 'cdr1' ? 'CDR1'
    : scope === 'cdr2' ? 'CDR2'
    : scope === 'cdr3' ? 'CDR3'
    : 'All CDRs';

  const lines = [];
  for (const [chainClass, entries] of byChain) {
    if (lines.length > 0) lines.push('');
    const n = entries.length;

    // Precompute full identity matrix
    const identMatrix = Array.from({length: n}, (_, i) =>
      Array.from({length: n}, (_, j) => {
        if (i === j) return 100;
        return computePairwiseIdentity(entries[i].posMap, entries[j].posMap) ?? 0;
      })
    );

    // Reorder by hierarchical clustering (UPGMA, distance = 100 - identity)
    const order = hierarchicalClusterOrder(n, (i, j) => 100 - identMatrix[i][j]);
    const reordered = order.map(i => entries[i]);

    lines.push(`# Chain: ${chainClass} | Scope: ${scopeLabel} | Scheme: ${scheme} | N=${n}`);
    lines.push(['', ...reordered.map(e => csvEscape(e.label))].join(','));
    for (let i = 0; i < n; i++) {
      const row = [csvEscape(reordered[i].label)];
      for (let j = 0; j < n; j++) {
        row.push(identMatrix[order[i]][order[j]].toFixed(2));
      }
      lines.push(row.join(','));
    }
  }

  return lines.join('\n') + '\n';
}

document.getElementById('downloadBtn').addEventListener('click', () => {
  if (!window._lastResults) return;
  const scheme = window._lastScheme || document.getElementById('scheme').value;
  const results = getResultsForDownload();
  const outputFormat = document.getElementById('outputFormat').value;

  let text, mimeType, filename;
  if (outputFormat === 'anarci') {
    text = buildAnarciText(results);
    mimeType = 'text/plain';
    filename = `anarci_${scheme}.anarci`;
  } else if (outputFormat === 'fasta-aa') {
    text = buildFastaAminoAcids(results);
    mimeType = 'text/plain';
    filename = `anarci_${scheme}_aa.fasta`;
  } else if (outputFormat === 'fasta-dna') {
    text = buildFastaDna(results);
    mimeType = 'text/plain';
    filename = `anarci_${scheme}_dna.fasta`;
  } else if (outputFormat.startsWith('pim-')) {
    const scope = outputFormat.slice(4);
    text = compute_pim(JSON.stringify(results), scope, scheme);
    mimeType = 'text/csv';
    filename = `anarci_pim_${scope}_${scheme}.csv`;
  } else {
    text = buildPythonStyleCsv(results);
    mimeType = 'text/csv';
    filename = `anarci_${scheme}.csv`;
  }

  const blob = new Blob([text], { type: mimeType });
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = filename;
  a.click();
});

initWasm();
