import init, { number_sequences_with_options } from './pkg/anarci_wasm.js';

let initPromise = null;

async function ensureWasmReady() {
  if (!initPromise) {
    initPromise = init();
  }
  await initPromise;
}

self.onmessage = async event => {
  const {
    taskId,
    fasta,
    scheme,
    threshold,
    allowedSpeciesJson,
    restrict,
    inputType,
  } = event.data || {};

  try {
    await ensureWasmReady();
    const json = number_sequences_with_options(
      fasta,
      scheme,
      threshold,
      allowedSpeciesJson,
      restrict,
      inputType,
    );

    self.postMessage({ taskId, ok: true, json });
  } catch (error) {
    self.postMessage({
      taskId,
      ok: false,
      error: error instanceof Error ? error.message : String(error),
    });
  }
};
