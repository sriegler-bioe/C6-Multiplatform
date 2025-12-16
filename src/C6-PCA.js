import { revcomp, gccontent } from './C6-Seq.js';

/**
 * @file C6-PCA.js
 * @author Spencer Riegler with Claude Sonnet 4.5 and Christopher Anderson
 * @copyright 2025 University of California, Berkeley
 * @license See the LICENSE file included in the repository
 * @version 1.0.0
 * @module C6-PCA
 * @description
 * Deterministic Polymerase Chain Assembly (PCA) simulator.
 *
 * Models PCA as overlap-driven anneal/extend cycles under idealized conditions:
 * - Exact-match overlaps only
 * - Automatic orientation (forward + reverse-complement considered)
 * - Deterministic enumeration of feasible products (no yields/kinetics)
 *
 * Includes evidence-grounded overlap validation:
 * - Hard cutoff: overlap < 15 bp (lowest clear observed value in the literature, Sequiera et al. 2016) 
 * - Warning: 15–19 bp (> 20 is ideal)
 *
 * Includes approximate Tm screening (Wallace rule) for overlap regions:
 * - Used only for warnings and guidance, not to gate assembly.
 */

/**
 * Normalize and validate a DNA sequence.
 * Accepts mixed case; removes whitespace; converts to uppercase.
 * Only A/C/G/T are allowed.
 *
 * @param {string} seq
 * @returns {string}
 */
function normalizeDNA(seq) {
  if (typeof seq !== 'string') throw new Error('Sequence must be a string.');
  const s = seq.replace(/\s+/g, '').toUpperCase();
  if (s.length === 0) throw new Error('Sequence is empty after normalization.');
  if (!/^[ACGT]+$/.test(s)) {
    throw new Error(`Invalid DNA sequence (allowed A,C,G,T only): ${seq}`);
  }
  return s;
}

/**
 * Wallace-rule approximate Tm (°C): 2*(A+T) + 4*(G+C)
 * Reasonable as a lightweight screening metric for short oligos/overlaps.
 *
 * @param {string} seq DNA (ACGT), already normalized
 * @returns {number}
 */
function tmWallace(seq) {
  const s = normalizeDNA(seq);
  const a = (s.match(/A/g) || []).length;
  const t = (s.match(/T/g) || []).length;
  const g = (s.match(/G/g) || []).length;
  const c = (s.match(/C/g) || []).length;
  return 2 * (a + t) + 4 * (g + c);
}

/**
 * Generate oriented candidates for each oligo (forward + reverse-complement).
 *
 * @param {{name:string, sequence:string}[]} oligos
 * @returns {{name:string, sequence:string, orientation:'F'|'R'}[]}
 */
function orientOligos(oligos) {
  if (!Array.isArray(oligos) || oligos.length === 0) {
    throw new Error('runPCA requires a non-empty oligos array.');
  }
  const out = [];
  for (const o of oligos) {
    if (!o || typeof o.name !== 'string') throw new Error('Each oligo must have a string "name".');
    const seqF = normalizeDNA(o.sequence);
    out.push({ name: o.name, sequence: seqF, orientation: 'F' });

    const seqR = normalizeDNA(revcomp(seqF));
    // Avoid duplicating palindromic sequences (F == R)
    if (seqR !== seqF) out.push({ name: o.name, sequence: seqR, orientation: 'R' });
  }
  return out;
}

/**
 * Find the best (longest) suffix→prefix exact overlap between A and B.
 * Returns 0 if none with length >= minOverlap.
 *
 * @param {string} a
 * @param {string} b
 * @param {number} minOverlap
 * @returns {number} overlap length
 */
function bestSuffixPrefixOverlap(a, b, minOverlap) {
  const A = normalizeDNA(a);
  const B = normalizeDNA(b);
  const maxK = Math.min(A.length, B.length);
  for (let k = maxK; k >= minOverlap; k--) {
    if (A.slice(A.length - k) === B.slice(0, k)) return k;
  }
  return 0;
}

/**
 * Deterministically extend A by B using overlap length k:
 * A + B[k:].
 *
 * @param {string} a
 * @param {string} b
 * @param {number} k
 * @returns {string}
 */
function extendByOverlap(a, b, k) {
  const A = normalizeDNA(a);
  const B = normalizeDNA(b);
  if (k <= 0) throw new Error('extendByOverlap requires k > 0.');
  return A + B.slice(k);
}

/**
 * Build overlap edges among a pool of oriented sequences.
 *
 * @param {{name:string, sequence:string, orientation:'F'|'R'}[]} pool
 * @param {number} minOverlap
 * @returns {Array<{from:number,to:number,k:number,overlapSeq:string,tm:number,gc:number}>}
 */
function buildOverlapEdges(pool, minOverlap) {
  const edges = [];
  for (let i = 0; i < pool.length; i++) {
    for (let j = 0; j < pool.length; j++) {
      if (i === j) continue;
      const a = pool[i].sequence;
      const b = pool[j].sequence;
      const k = bestSuffixPrefixOverlap(a, b, minOverlap);
      if (k >= minOverlap) {
        const overlapSeq = a.slice(a.length - k);
        const tm = tmWallace(overlapSeq);
        const gc = gccontent(overlapSeq);
        edges.push({ from: i, to: j, k, overlapSeq, tm, gc });
      }
    }
  }
  return edges;
}

/**
 * Run a deterministic PCA simulation.
 *
 * @param {{name:string, sequence:string}[]} oligos
 * @param {Object} [options]
 * @param {number} [options.cycles] default: oligos.length
 * @param {number} [options.minOverlap] default: 15 (hard cutoff)
 * @param {number} [options.warnOverlap] default: 20
 * @param {number} [options.longProductBp] default: 1000
 * @returns {{
 *   products: Array<{sequence:string,length:number,supports:string[],circularizable:boolean}>,
 *   warnings: string[],
 *   tmSummary: {min:number|null,max:number|null,spread:number|null},
 *   protocolAdvice: {recommendedPCA:'one-step'|'two-step', rationale:string}
 * }}
 */
function runPCA(oligos, options = {}) {
  const minOverlap = (options.minOverlap ?? 15);
  const warnOverlap = (options.warnOverlap ?? 20);
  const longProductBp = (options.longProductBp ?? 1000);

  if (!Number.isInteger(minOverlap) || minOverlap < 1) throw new Error('minOverlap must be a positive integer.');
  if (minOverlap !== 15) {
    throw new Error('This PCA simulator v1.0 fixes minOverlap at 15 bp).');
  }

  const oriented = orientOligos(oligos);
  const cycles = (options.cycles ?? oligos.length);
  if (!Number.isInteger(cycles) || cycles < 1) throw new Error('cycles must be a positive integer.');

  const warnings = [];

  // Precompute base overlap graph among oriented oligos (for early fatal detection + Tm summary).
  const baseEdges = buildOverlapEdges(oriented, minOverlap);
  if (baseEdges.length === 0) {
    throw new Error(`No overlaps ≥ ${minOverlap} bp found among oligos even after considering reverse complements.`);
  }

  // Overlap warnings and Tm summary from baseEdges (screening).
  let tmMin = null, tmMax = null;
  for (const e of baseEdges) {
    if (tmMin === null || e.tm < tmMin) tmMin = e.tm;
    if (tmMax === null || e.tm > tmMax) tmMax = e.tm;

    if (e.k < warnOverlap) {
      warnings.push(
        `Short overlap (${e.k} bp) between ${oriented[e.from].name}(${oriented[e.from].orientation}) → ` +
        `${oriented[e.to].name}(${oriented[e.to].orientation}) (below recommended ${warnOverlap} bp).`
      );
    } else if (e.k > 40) {
      warnings.push(
        `Long overlap (${e.k} bp) between ${oriented[e.from].name}(${oriented[e.from].orientation}) → ` +
        `${oriented[e.to].name}(${oriented[e.to].orientation}); may increase mispriming/secondary-structure risk.`
      );
    }
  }
  const tmSpread = (tmMin === null || tmMax === null) ? null : (tmMax - tmMin);

  // Product representation:
  // sequence => { sequence, supports:Set<string> }
  let pool = new Map();
  for (const o of oriented) {
    pool.set(o.sequence, { sequence: o.sequence, supports: new Set([o.name]) });
  }

  // Deterministic cycle expansion
  for (let c = 0; c < cycles; c++) {
    const entries = Array.from(pool.values());
    const newProducts = new Map();

    for (let i = 0; i < entries.length; i++) {
      for (let j = 0; j < entries.length; j++) {
        if (i === j) continue;
        const A = entries[i].sequence;
        const B = entries[j].sequence;

        const k = bestSuffixPrefixOverlap(A, B, minOverlap);
        if (k >= minOverlap) {
          const seq = extendByOverlap(A, B, k);
          const supports = new Set(entries[i].supports);
          for (const s of entries[j].supports) supports.add(s);

          // Deduplicate: keep version with more supports if sequence already exists
          const existing = newProducts.get(seq) || pool.get(seq);
          if (!existing) {
            newProducts.set(seq, { sequence: seq, supports });
          } else {
            const existingCount = existing.supports ? existing.supports.size : 0;
            if (supports.size > existingCount) {
              newProducts.set(seq, { sequence: seq, supports });
            }
          }
        }
      }
    }

    // If no expansion, stop early
    if (newProducts.size === 0) break;

    // Merge new products into pool
    for (const [seq, obj] of newProducts.entries()) {
      const prev = pool.get(seq);
      if (!prev || obj.supports.size > prev.supports.size) {
        pool.set(seq, obj);
      }
    }
  }

  // Final classification
  const products = Array.from(pool.values()).map(p => {
    const seq = p.sequence;
    const circularizable = (bestSuffixPrefixOverlap(seq, seq, minOverlap) >= minOverlap);
    return {
      sequence: seq,
      length: seq.length,
      supports: Array.from(p.supports).sort(),
      circularizable
    };
  });

  // Sort: prefer longer, then more supports, then lexicographic
  products.sort((x, y) => {
    if (y.length !== x.length) return y.length - x.length;
    if (y.supports.length !== x.supports.length) return y.supports.length - x.supports.length;
    return x.sequence.localeCompare(y.sequence);
  });

  // Long product warning
  if (products.length > 0 && products[0].length > longProductBp) {
    warnings.push(
      `Largest assembled product is ${products[0].length} bp (> ${longProductBp} bp). ` +
      `Large multi-oligo assemblies can be limited by competitive annealing and oligo synthesis errors; ` +
      `consider two-step PCA, subpools, or post-assembly cleanup strategies.`
    );
  }

  // Protocol advice (heuristic)
  const best = products[0];
  let recommendedPCA = 'one-step';
  let rationale = 'Heuristic: shorter/simple assemblies typically work in one-step PCA.';
  if (best) {
    const nOligos = new Set(best.supports).size;
    if (best.length > 1000 || nOligos > 4) {
      recommendedPCA = 'two-step';
      rationale =
        `Heuristic: assembled length ${best.length} bp and/or ${nOligos} oligos suggests increased complexity; ` +
        `two-step PCA (assembly then outer-primer PCR) is often more robust for longer constructs.`;
    } else {
      recommendedPCA = 'one-step';
      rationale =
        `Heuristic: assembled length ${best.length} bp using ${nOligos} oligos is relatively simple; ` +
        `one-step PCA is typically sufficient.`;
    }
  }

  return {
    products,
    warnings: dedupeStrings(warnings),
    tmSummary: { min: tmMin, max: tmMax, spread: tmSpread },
    protocolAdvice: { recommendedPCA, rationale }
  };
}

/**
 * Dedupe strings while preserving order.
 * @param {string[]} arr
 * @returns {string[]}
 */
function dedupeStrings(arr) {
  const seen = new Set();
  const out = [];
  for (const s of arr) {
    if (!seen.has(s)) { seen.add(s); out.push(s); }
  }
  return out;
}

export {
  runPCA,
  normalizeDNA,
  tmWallace,
  bestSuffixPrefixOverlap,
  extendByOverlap,
  buildOverlapEdges
};
