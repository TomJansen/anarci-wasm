// Minimal ABIF (.ab1) parser — extracts base calls, peak locations, and the
// four trace channels. No external dependencies.
//
// ABIF is a tagged binary format (Applied Biosystems). Layout:
//   - 4-byte magic "ABIF"
//   - version (int16)
//   - a root directory entry (28 bytes) at offset 6, whose dataoffset points to
//     an array of 28-byte directory entries.
// Each directory entry names a tag (e.g. "PBAS", "PLOC", "DATA") + tag number,
// an element type, count, and either an inline value or an offset to the data.
//
// We pull:
//   PBAS.2 — edited base calls (ASCII)               -> sequence string
//   PLOC.2 — peak location index per base (int16[])  -> x-position of each base
//   DATA.9..12 — the 4 processed trace channels (int16[])
//   FWO_.1 — base order for the 4 channels, e.g. "GATC"
//
// Returns null if the bytes are not a valid ABIF file.

const ABIF_TYPES = {
  1: 'byte', // unsigned byte
  2: 'char', // signed char (we read ASCII)
  3: 'word', // unsigned short
  4: 'short', // signed short
  5: 'long', // signed int32
  7: 'float',
  8: 'double',
  10: 'date',
  11: 'time',
  18: 'pString',
  19: 'cString',
};

const ENTRY_SIZE = 28;

export function parseAb1(arrayBuffer) {
  const dv = new DataView(arrayBuffer);
  if (dv.byteLength < 6) return null;

  // Magic "ABIF"
  if (
    dv.getUint8(0) !== 0x41 ||
    dv.getUint8(1) !== 0x42 ||
    dv.getUint8(2) !== 0x49 ||
    dv.getUint8(3) !== 0x46
  ) {
    return null;
  }

  // The root entry sits at offset 6 and describes the directory. Within a
  // 28-byte entry: name@0(4) num@4(4) etype@8(2) esize@10(2) numelements@12(4)
  // datasize@16(4) dataoffset@20(4).
  const rootNumEntries = dv.getInt32(6 + 12, false);
  const rootDataOffset = dv.getInt32(6 + 20, false);

  const entries = new Map();
  for (let i = 0; i < rootNumEntries; i++) {
    const entry = readDirEntry(dv, rootDataOffset + i * ENTRY_SIZE);
    if (entry) entries.set(`${entry.name}.${entry.number}`, entry);
  }

  const baseOrder = readString(dv, entries.get('FWO_.1')) || 'GATC';
  const bases = readString(dv, entries.get('PBAS.2') || entries.get('PBAS.1'));
  const peakLocations = readShortArray(dv, entries.get('PLOC.2') || entries.get('PLOC.1'));
  // PCON.2 = per-base Phred quality (one byte per base). May be absent.
  const quality = readByteArray(dv, entries.get('PCON.2') || entries.get('PCON.1'));

  // Channels DATA.9..12 are the processed traces, one per base in baseOrder.
  const rawChannels = [9, 10, 11, 12].map(n =>
    readShortArray(dv, entries.get(`DATA.${n}`))
  );

  if (!bases || rawChannels.some(c => !c)) return null;

  // Map each channel to its base via FWO_ order so callers get {A,C,G,T}.
  const channels = {};
  for (let i = 0; i < 4 && i < baseOrder.length; i++) {
    channels[baseOrder[i]] = rawChannels[i];
  }

  return {
    bases: bases.toUpperCase(),
    peakLocations: peakLocations || [],
    channels, // { A: int16[], C: ..., G: ..., T: ... }
    quality: quality || [], // per-base Phred quality, empty if absent
    baseOrder,
    traceLength: rawChannels[0].length,
  };
}

function readDirEntry(dv, offset) {
  if (offset + ENTRY_SIZE > dv.byteLength) return null;
  const name = String.fromCharCode(
    dv.getUint8(offset),
    dv.getUint8(offset + 1),
    dv.getUint8(offset + 2),
    dv.getUint8(offset + 3)
  );
  const number = dv.getInt32(offset + 4, false);
  const elementType = dv.getInt16(offset + 8, false);
  const elementSize = dv.getInt16(offset + 10, false);
  const numElements = dv.getInt32(offset + 12, false);
  const dataSize = dv.getInt32(offset + 16, false);
  const dataOffset = dv.getInt32(offset + 20, false);
  return {
    name,
    number,
    type: ABIF_TYPES[elementType] || elementType,
    elementSize,
    numElements,
    dataSize,
    // Values <= 4 bytes are stored inline in the dataOffset field itself.
    dataOffset: dataSize <= 4 ? offset + 20 : dataOffset,
  };
}

function readString(dv, entry) {
  if (!entry) return null;
  let start = entry.dataOffset;
  let len = entry.numElements;
  // pString: first byte is the length; cString: NUL-terminated.
  if (entry.type === 'pString') {
    len = dv.getUint8(start);
    start += 1;
  }
  let s = '';
  for (let i = 0; i < len; i++) {
    const c = dv.getUint8(start + i);
    if (entry.type === 'cString' && c === 0) break;
    s += String.fromCharCode(c);
  }
  return s;
}

function readShortArray(dv, entry) {
  if (!entry) return null;
  const out = new Array(entry.numElements);
  for (let i = 0; i < entry.numElements; i++) {
    out[i] = dv.getInt16(entry.dataOffset + i * 2, false);
  }
  return out;
}

function readByteArray(dv, entry) {
  if (!entry) return null;
  const out = new Array(entry.numElements);
  for (let i = 0; i < entry.numElements; i++) {
    out[i] = dv.getUint8(entry.dataOffset + i);
  }
  return out;
}
