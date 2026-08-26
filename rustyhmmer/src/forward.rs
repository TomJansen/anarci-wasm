












use crate::alphabet::{AMINO_FREQ, K, KP};
use crate::hmmfile::{P7Hmm, TDD, TDM, TIM, TII, TMD, TMI, TMM};

fn degen_set(code: usize) -> &'static [usize] {
    match code {
        21 => &[11, 2],
        22 => &[7, 9],
        23 => &[13, 3],
        24 => &[8],
        25 => &[1],
        26 => &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19],
        _ => &[],
    }
}


pub struct ForwardFilter {
    pub m: usize,
    /// `p7O_NQF(M)`: the number of striped f32 vectors the DD relaxation walks.
    pub(crate) q_n: usize,
    /// tDD in C's striped vector-major order, `4*q_n` long. See `dd_relax`.
    pub(crate) tdd_v: Vec<f32>,
    
    
    pub(crate) rfv: Vec<[f32; KP]>,
    
    
    
    
    
    
    pub(crate) rfv_t: Vec<f32>,
    
    pub(crate) tbm: Vec<f32>,
    
    pub(crate) amm: Vec<f32>,
    pub(crate) aim: Vec<f32>,
    pub(crate) adm: Vec<f32>,
    
    pub(crate) tmi: Vec<f32>,
    pub(crate) tii: Vec<f32>,
    pub(crate) tmd: Vec<f32>,
    pub(crate) tdd: Vec<f32>,
    pub ftau: f64,
    pub flambda: f64,
}

/// The D-state row of one Forward pass, computed the way C's striped kernel computes it.
///
/// C, impl_sse/fwdback.c:344-395. After the main `q` loop, `DMO(dpc,q)` holds only the
/// delayed M->D store and `dcv` holds the one that fell off the end; the DD paths are
/// then relaxed in four passes:
///
/// ```text
///   dcv        = esl_sse_rightshiftz_float(dcv);
///   DMO(dpc,0) = zerov;
///   tp         = om->tfv + 7*Q;  /* set tp to start of the DD's */
///   for (q = 0; q < Q; q++)
///     {
///       DMO(dpc,q) = _mm_add_ps(dcv, DMO(dpc,q));
///       dcv        = _mm_mul_ps(DMO(dpc,q), *tp); tp++; /* extend DMO(q), so we include M->D and D->D paths */
///     }
///
///   if (om->M < 100)
///     {                  /* Fully serialized version */
///       for (j = 1; j < 4; j++)
///         {
///           dcv = esl_sse_rightshiftz_float(dcv);
///           tp  = om->tfv + 7*Q;
///           for (q = 0; q < Q; q++)
///             { /* note, extend dcv, not DMO(q); only adding DD paths now */
///               DMO(dpc,q) = _mm_add_ps(dcv, DMO(dpc,q));
///               dcv        = _mm_mul_ps(dcv, *tp);   tp++;
///             }
///         }
///     }
///   else
///     {                  /* Slightly parallelized version, but which incurs some overhead */
///       for (j = 1; j < 4; j++)
///         {
///           register __m128 cv;  /* keeps track of whether any DD's change DMO(q) */
///
///           dcv = esl_sse_rightshiftz_float(dcv);
///           tp  = om->tfv + 7*Q;
///           cv  = zerov;
///           for (q = 0; q < Q; q++)
///             { /* using cmpgt below tests if DD changed any DMO(q) *without* conditional branch */
///               sv         = _mm_add_ps(dcv, DMO(dpc,q));
///               cv         = _mm_or_ps(cv, _mm_cmpgt_ps(sv, DMO(dpc,q)));
///               DMO(dpc,q) = sv;                                     /* store new DMO(q) */
///               dcv        = _mm_mul_ps(dcv, *tp);   tp++;            /* note, extend dcv, not DMO(q) */
///             }
///           if (! _mm_movemask_ps(cv)) break; /* DD's didn't change any DMO(q)? Then done, break out. */
///         }
///     }
/// ```
///
/// # Why this replaces the serial recurrence
///
/// The obvious form, `dc[k] = mc[k-1]*tmd[k-1] + dc[k-1]*tdd[k-1]`, is the *same sum*
/// as the four passes above — they differ only in how the terms are associated, and
/// four passes over four lanes is exactly enough to carry a D path from lane 0 through
/// to lane 3, so both are mathematically exact. But f32 addition is not associative,
/// so they are not the same float. The Forward score depends on this at the ±1 ULP
/// level, which is enough to move the last printed digit of a score or E-value.
///
/// # The lane mapping
///
/// Striped vector `q` lane `z` holds model position `k = q + 1 + z*q_n`
/// (`p7_oprofile.c:948`, `k + z*nq`). `esl_sse_rightshiftz_float` (lane `z` takes lane
/// `z-1`, zero shifted in) is therefore the `k -> k+1` carry between stripes, which is
/// why four passes over four lanes is exactly enough.
///
/// The seed is the delayed M->D store: the slot for `k` holds `mc[k-1] * tmd[k-1]` —
/// the same M->D term as the serial recurrence. C clears all four lanes of vector 0
/// (`DMO(dpc,0) = zerov`) and re-supplies lanes 1..3 of it from the shifted `dcv`,
/// which carries exactly `mc[z*q_n]*tmd[z*q_n]`; lane 0 of vector 0 is `k = 1`, where
/// `dc[1] = 0`.
///
/// The `4*q_n` slots cover `k = 1..=4*q_n`, which is `>= m`; the tail beyond `m` is the
/// model's zero padding, and `tdd_v` being zero there is what stops the DD tail from
/// leaking off the end of the model into the E-state sum.
///
/// # Layout
///
/// `d` is scratch in C's *vector-major* order — `d[q*4 + z]` is vector `q` lane `z`, so
/// each vector's four lanes are contiguous and one `__m128` load. `tdd_v` is tDD in the
/// same order. `mc`/`dc` stay in `k` order (that is what the rest of the DP and the
/// posterior decoding index by), so the seed gathers and the result scatters across the
/// stride. `d` is fully overwritten on entry.
///
/// **Optimization divergence (DEVELOPMENT_PRINCIPLES §3):** (a) C holds the row in
/// striped `__m128` registers/memory throughout, so no gather or scatter is needed at
/// all. (b) A first cut of this function instead kept `d` in `k` order and walked the
/// four lanes as `d[z*q_n + q]`, which is the same arithmetic but turns each of the
/// `16*q_n` inner steps into a strided scalar load/store: it measured **+10.9%
/// instructions** over the whole search (63.17G -> 70.06G, `perf stat -e instructions`,
/// 30k proteins x 117 profiles, `--cpu 1`, pinned). Vector-major with the same SSE2
/// operations C uses brings the passes down to `4*q_n` vector ops. (c) The gather/scatter
/// against `k` order is the price of not having ported the whole row to striped layout;
/// the arithmetic inside the passes is unchanged, and parity is re-verified below.
pub(crate) fn dd_relax(
    m: usize,
    q_n: usize,
    mc: &[f32],
    tmd: &[f32],
    tdd_v: &[f32],
    d: &mut [f32],
    dc: &mut [f32],
) {
    // Seed: the delayed M->D stores, with vector 0's four lanes cleared.
    //
    // `dcv` after C's main q loop is the store that fell off the end,
    // `mc[(z+1)*q_n] * tmd[(z+1)*q_n]`; `esl_sse_rightshiftz_float` then makes lane z
    // carry `mc[z*q_n] * tmd[z*q_n]` and clears lane 0.
    let mut dcv = [0.0f32; 4];
    for z in 1..4 {
        let s = z * q_n;
        dcv[z] = if s <= m { mc[s] * tmd[s] } else { 0.0 };
    }
    d[0..4].fill(0.0);
    for q in 1..q_n {
        for z in 0..4 {
            let s = z * q_n + q;
            d[q * 4 + z] = if s <= m { mc[s] * tmd[s] } else { 0.0 };
        }
    }

    #[cfg(target_arch = "x86_64")]
    // SAFETY: SSE2 is part of the x86_64 baseline, so the intrinsics below are always
    // available. Every load/store is a full 16-byte vector inside `d`/`tdd_v`, both of
    // which are `4*q_n` long, and `q < q_n`.
    unsafe {
        dd_relax_sse2(m, q_n, tdd_v, d, &mut dcv)
    }
    #[cfg(not(target_arch = "x86_64"))]
    dd_relax_scalar(m, q_n, tdd_v, d, &mut dcv);

    dc[0] = 0.0;
    for q in 0..q_n {
        for z in 0..4 {
            let k = z * q_n + q + 1;
            if k <= m {
                dc[k] = d[q * 4 + z];
            }
        }
    }
}

/// The four relaxation passes, vector-major, using the same SSE2 operations C uses.
#[cfg(target_arch = "x86_64")]
#[inline]
unsafe fn dd_relax_sse2(m: usize, q_n: usize, tdd_v: &[f32], d: &mut [f32], dcv: &mut [f32; 4]) { unsafe {
    use core::arch::x86_64::*;

    let dp = d.as_mut_ptr();
    let tp0 = tdd_v.as_ptr();
    let mut dv = _mm_loadu_ps(dcv.as_ptr());

    // Pass 1: extend DMO(q), so this pass carries both the M->D and the D->D paths.
    for q in 0..q_n {
        let dq = dp.add(q * 4);
        let sv = _mm_add_ps(dv, _mm_loadu_ps(dq));
        _mm_storeu_ps(dq, sv);
        dv = _mm_mul_ps(sv, _mm_loadu_ps(tp0.add(q * 4)));
    }

    // Passes 2..4: extend `dcv`, not DMO(q) — only DD paths are being added now.
    //
    // The two branches differ only in C's early exit, which is data-dependent and so
    // cannot be assumed to be a no-op: once a pass adds nothing anywhere, the *next*
    // pass still right-shifts, and a `dcv` that was negligible against lane z's D
    // values is not necessarily negligible against lane z+1's. Both branches are
    // transcribed so the choice of branch reproduces C's for every M.
    let zerov = _mm_setzero_ps();
    let track_changes = m >= 100;
    for _ in 1..4 {
        // esl_sse_rightshiftz_float: [4 8 12 x] -> [x 4 8 12], zeros shifted on.
        dv = _mm_castsi128_ps(_mm_slli_si128(_mm_castps_si128(dv), 4));
        let mut cv = zerov;
        for q in 0..q_n {
            let dq = dp.add(q * 4);
            let dmo = _mm_loadu_ps(dq);
            let sv = _mm_add_ps(dv, dmo);
            cv = _mm_or_ps(cv, _mm_cmpgt_ps(sv, dmo));
            _mm_storeu_ps(dq, sv);
            dv = _mm_mul_ps(dv, _mm_loadu_ps(tp0.add(q * 4)));
        }
        if track_changes && _mm_movemask_ps(cv) == 0 {
            break;
        }
    }
    _mm_storeu_ps(dcv.as_mut_ptr(), dv);
}}

/// Scalar equivalent of [`dd_relax_sse2`] for non-x86_64 targets. Kept in lockstep with
/// it lane by lane so both produce the same floats.
#[cfg(not(target_arch = "x86_64"))]
fn dd_relax_scalar(m: usize, q_n: usize, tdd_v: &[f32], d: &mut [f32], dcv: &mut [f32; 4]) {
    let mut dv = *dcv;
    for q in 0..q_n {
        for z in 0..4 {
            let i = q * 4 + z;
            d[i] += dv[z];
            dv[z] = d[i] * tdd_v[i];
        }
    }
    let track_changes = m >= 100;
    for _ in 1..4 {
        dv = [0.0, dv[0], dv[1], dv[2]];
        let mut changed = false;
        for q in 0..q_n {
            for z in 0..4 {
                let i = q * 4 + z;
                let sv = d[i] + dv[z];
                changed |= sv > d[i];
                d[i] = sv;
                dv[z] *= tdd_v[i];
            }
        }
        if track_changes && !changed {
            break;
        }
    }
    *dcv = dv;
}

/// The D-state row of one Backward pass, the way C's striped kernel computes it.
///
/// This is [`dd_relax`] mirrored: the Backward recurrence pulls `D(k)` from `D(k+1)`, so
/// the passes walk `q` downwards and the inter-stripe carry is a *left* shift.
///
/// C, impl_sse/fwdback.c:616-641 (the per-residue row; phase 3 does one DD step
/// together with the `{MD}->E` paths, phase 4 finishes the other three):
///
/// ```text
///   /* phase 3: {MD}->E paths and one step of the D->D paths */
///   tp  = om->tfv + 8*Q - 1;      /* <*tp> now the [4 8 12 x] TDD quad */
///   dpv = _mm_add_ps(DMO(dpc,0), xEv);
///   dpv = _mm_move_ss(dpv, zerov);
///   dpv = _mm_shuffle_ps(dpv, dpv, _MM_SHUFFLE(0,3,2,1));
///   for (q = Q-1; q >= 0; q--)
///     {
///       dcv        = _mm_mul_ps(dpv, *tp); tp--;
///       DMO(dpc,q) = _mm_add_ps(DMO(dpc,q), _mm_add_ps(dcv, xEv));
///       dpv        = DMO(dpc,q);
///       MMO(dpc,q) = _mm_add_ps(MMO(dpc,q), xEv);
///     }
///
///   /* phase 4: finish extending the DD paths */
///   /* fully serialized for now */
///   for (j = 1; j < 4; j++)   /* three passes: we've already done 1 segment, we need 4 total */
///     {
///       dcv = _mm_move_ss(dcv, zerov);
///       dcv = _mm_shuffle_ps(dcv, dcv, _MM_SHUFFLE(0,3,2,1));
///       tp  = om->tfv + 8*Q - 1;      /* <*tp> now the [4 8 12 x] TDD quad */
///       for (q = Q-1; q >= 0; q--)
///         {
///           dcv        = _mm_mul_ps(dcv, *tp); tp--;
///           DMO(dpc,q) = _mm_add_ps(DMO(dpc,q), dcv);
///         }
///     }
/// ```
///
/// Note there is no early exit here — Backward is "fully serialized for now" for every
/// M, unlike Forward's `M >= 100` branch.
///
/// The row-L initialisation (impl_sse/fwdback.c:504-522) is the same shape with
/// `xEv` already sitting in every `DMO` lane and so no `xEv` addend in the pass:
///
/// ```text
///   for (q = 0; q < Q; q++) MMO(dpc,q) = DMO(dpc,q) = xEv;
///   ...
///   dpv = _mm_move_ss(DMO(dpc,Q-1), zerov);
///   dpv = _mm_shuffle_ps(dpv, dpv, _MM_SHUFFLE(0,3,2,1));
///   for (q = Q-1; q >= 0; q--)
///     {
///       dcv        = _mm_mul_ps(dpv, *tp);      tp--;
///       DMO(dpc,q) = _mm_add_ps(DMO(dpc,q), dcv);
///       dpv        = DMO(dpc,q);
///     }
/// ```
///
/// so both callers share this function, distinguished only by `addend`: `xE` for a
/// per-residue row, `0.0` for row L. (`x + (y + 0.0) == x + y` exactly for the
/// non-negative values in a probability-space Backward matrix.)
///
/// On entry `d` holds the phase-1 D seed in the vector-major layout [`dd_relax`]
/// documents; on return `d` is relaxed and `dc[1..=m]` holds it in `k` order.
pub(crate) fn bck_dd_relax(
    q_n: usize,
    m: usize,
    addend: f32,
    tdd_v: &[f32],
    d: &mut [f32],
    dc: &mut [f32],
) {
    #[cfg(target_arch = "x86_64")]
    // SAFETY: SSE2 is baseline on x86_64. Every access is a full 16-byte vector at
    // `q*4` inside `d`/`tdd_v`, both `4*q_n` long, with `q < q_n`.
    unsafe {
        use core::arch::x86_64::*;
        let dp = d.as_mut_ptr();
        let tp0 = tdd_v.as_ptr();
        let xev = _mm_set1_ps(addend);

        // phase 3: dpv = leftshift(DMO(dpc,0) + xEv)
        let mut dpv = _mm_add_ps(_mm_loadu_ps(dp), xev);
        dpv = _mm_castsi128_ps(_mm_srli_si128(_mm_castps_si128(dpv), 4));
        let mut dcv = _mm_setzero_ps();
        for q in (0..q_n).rev() {
            let dq = dp.add(q * 4);
            dcv = _mm_mul_ps(dpv, _mm_loadu_ps(tp0.add(q * 4)));
            dpv = _mm_add_ps(_mm_loadu_ps(dq), _mm_add_ps(dcv, xev));
            _mm_storeu_ps(dq, dpv);
        }
        // phase 4: three more passes, extending dcv only
        for _ in 1..4 {
            dcv = _mm_castsi128_ps(_mm_srli_si128(_mm_castps_si128(dcv), 4));
            for q in (0..q_n).rev() {
                let dq = dp.add(q * 4);
                dcv = _mm_mul_ps(dcv, _mm_loadu_ps(tp0.add(q * 4)));
                _mm_storeu_ps(dq, _mm_add_ps(_mm_loadu_ps(dq), dcv));
            }
        }
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        let mut dpv = [
            d[1] + addend,
            d[2] + addend,
            d[3] + addend,
            0.0,
        ];
        let mut dcv = [0.0f32; 4];
        for q in (0..q_n).rev() {
            for z in 0..4 {
                let i = q * 4 + z;
                dcv[z] = dpv[z] * tdd_v[i];
                d[i] += dcv[z] + addend;
                dpv[z] = d[i];
            }
        }
        for _ in 1..4 {
            dcv = [dcv[1], dcv[2], dcv[3], 0.0];
            for q in (0..q_n).rev() {
                for z in 0..4 {
                    let i = q * 4 + z;
                    dcv[z] *= tdd_v[i];
                    d[i] += dcv[z];
                }
            }
        }
    }

    for q in 0..q_n {
        for z in 0..4 {
            let k = z * q_n + q + 1;
            if k <= m {
                dc[k] = d[q * 4 + z];
            }
        }
    }
}

impl ForwardFilter {

    pub fn build(hmm: &P7Hmm) -> Self {
        let m = hmm.m;

        
        
        let mut rfv = vec![[0.0f32; KP]; m + 1];
        for k in 1..=m {
            let mut sc = [f32::NEG_INFINITY; KP];
            for x in 0..K {
                // C, modelconfig.c:144 (`p7_ProfileConfig`):
                //   sc[x] = log((double)hmm->mat[k][x] / bg->f[x]);
                // `bg->f` is `float *` (p7_bg.c:`ESL_ALLOC(bg->f, sizeof(float)*abc->K)`),
                // so the divisor is the *f32-rounded* Swiss-Prot frequency widened back to
                // double -- not the f64 literal. `sc` is `float`, so the double `log` is
                // rounded once on store.
                sc[x] = ((hmm.mat[k][x] as f64) / (AMINO_FREQ[x] as f32 as f64)).ln() as f32;
            }
            for code in K..KP {
                let set = degen_set(code);
                if set.is_empty() {
                    continue;
                }
                let mut left = 0.0f32;
                let mut right = 0.0f32;
                for &y in set {
                    left += (AMINO_FREQ[y] as f32) * sc[y];
                    right += AMINO_FREQ[y] as f32;
                }
                sc[code] = left / right;
            }
            for x in 0..KP {
                // C, impl_sse/p7_oprofile.c:947-953 (`fb_conversion`):
                //   for (k = 1, q = 0; q < nq; q++, k++)
                //     {
                //       for (z = 0; z < 4; z++) tmp.x[z] = (k+ z*nq <= M) ? p7P_MSC(gm, k+z*nq, x) : -eslINFINITY;
                //       om->rfv[x][q] = esl_sse_expf(tmp.v);
                //     }
                // `esl_sse_expf`, not `expf` -- see easel::sse.
                rfv[k][x] = crate::easel::sse::expf(sc[x]);
            }
        }
        
        let mut rfv_t = vec![0.0f32; KP * (m + 1)];
        for k in 1..=m {
            for x in 0..KP {
                rfv_t[x * (m + 1) + k] = rfv[k][x];
            }
        }

        
        let mut occ = vec![0.0f32; m + 1];
        if m >= 1 {
            occ[1] = hmm.t[0][TMI] + hmm.t[0][TMM];
        }
        for k in 2..=m {
            // C, p7_hmm.c:`p7_hmm_CalculateOccupancy`:
            //   mocc[k] = mocc[k-1] * (hmm->t[k-1][p7H_MM] + hmm->t[k-1][p7H_MI]) +
            //     (1.0-mocc[k-1]) * hmm->t[k-1][p7H_DM];
            //
            // `1.0` is a *double* literal, so C's usual arithmetic conversions compute
            // `1.0 - mocc[k-1]` in double, promote `t[k-1][p7H_DM]` to double for the
            // multiply, promote the (float) first term for the addition, and round to
            // f32 only on the store. Doing it all in f32 drifts by an ULP or two, and
            // because `occ` feeds `tbm` through `log(occ[k]/Z)` the drift lands in the
            // B->Mk entry probability the DP multiplies by on every row.
            occ[k] = ((occ[k - 1] * (hmm.t[k - 1][TMM] + hmm.t[k - 1][TMI])) as f64
                + (1.0 - occ[k - 1] as f64) * hmm.t[k - 1][TDM] as f64) as f32;
        }
        
        
        let mut z = 0.0f32;
        for k in 1..=m {
            z += occ[k] * (m - k + 1) as f32;
        }
        // C, modelconfig.c:96-97:
        //   for (k = 1; k <= hmm->M; k++)
        //     p7P_TSC(gm, k-1, p7P_BM) = log(occ[k] / Z);
        // then fb_conversion() exponentiates it with esl_sse_expf(). `occ[k] / Z` is an
        // f32 division (both are `float`); `log` is the double libm call, rounded once
        // when stored into the `float` gm->tsc.
        let mut tbm = vec![0.0f32; m + 1];
        for k in 1..=m {
            tbm[k] = crate::easel::sse::expf(((occ[k] / z) as f64).ln() as f32);
        }

        
        // Into-M transitions, held rotated by -1 like C's `gm->tsc` (so index k is the
        // transition *into* M_k, read from HMM node k-1).
        //
        // Node 0 is special: `p7_profile_Create()` (p7_profile.c:84) does
        //   esl_vec_FSet(gm->tsc, p7P_NTRANS, -eslINFINITY);  /* node 0 nonexistent,
        //                                                        has no transitions */
        // and `p7_ProfileConfig()`'s transition loop starts at k=1
        // (modelconfig.c:126), so only p7P_TSC(gm,0,p7P_BM) is ever written for node 0
        // — local-mode entry comes exclusively through `tbm`. So MM/IM/DM into M_1 are
        // -eslINFINITY in C, i.e. 0.0 after fb_conversion()'s expf().
        //
        // Forward/Backward/stochastic-traceback multiply these by cell M(i-1,0) /
        // I(i-1,0) / D(i-1,0), which are identically 0.0, so the raw HMM node-0
        // probabilities were harmless there. The optimal-accuracy DP is different: it
        // *selects* on `transition > 0` rather than multiplying, so a nonzero entry
        // here makes M(1,1) look reachable from M(0,0) and the OA traceback walks off
        // the matrix. Store what C stores.
        //
        // The stored value is *not* the HMM's transition probability. C round-trips it
        // through log space and back:
        //
        //   modelconfig.c:126-134 (`p7_ProfileConfig`):
        //     for (k = 1; k < gm->M; k++) {
        //       tp = gm->tsc + k * p7P_NTRANS;
        //       tp[p7P_MM] = log(hmm->t[k][p7H_MM]);
        //       ...
        //       tp[p7P_DD] = log(hmm->t[k][p7H_DD]);
        //     }
        //
        //   impl_sse/p7_oprofile.c:955-981 (`fb_conversion`):
        //     for (z = 0; z < 4; z++) tmp.x[z] = (kb+z*nq < M) ? p7P_TSC(gm, kb+z*nq, tg) : -eslINFINITY;
        //     om->tfv[j++] = esl_sse_expf(tmp.v);
        //
        // So the DP multiplies by `esl_sse_expf((float)log((double)p))`, which is a
        // couple of ULP away from `p`. Every one of the ~8M transitions is affected and
        // the DP multiplies L of them together per path, so skipping the round trip
        // moves the Forward score far more than any summation-order detail does.
        let rt = |p: f32| crate::easel::sse::expf((p as f64).ln() as f32);
        let mut amm = vec![0.0f32; m + 1];
        let mut aim = vec![0.0f32; m + 1];
        let mut adm = vec![0.0f32; m + 1];
        for k in 2..=m {
            amm[k] = rt(hmm.t[k - 1][TMM]);
            aim[k] = rt(hmm.t[k - 1][TIM]);
            adm[k] = rt(hmm.t[k - 1][TDM]);
        }

        let mut tmi = vec![0.0f32; m + 1];
        let mut tii = vec![0.0f32; m + 1];
        let mut tmd = vec![0.0f32; m + 1];
        let mut tdd = vec![0.0f32; m + 1];
        for k in 1..m {
            tmi[k] = rt(hmm.t[k][TMI]);
            tii[k] = rt(hmm.t[k][TII]);
            tmd[k] = rt(hmm.t[k][TMD]);
            tdd[k] = rt(hmm.t[k][TDD]);
        }

        // tDD in C's striped order, padded out to the full 4*Q slots the DD relaxation
        // walks. Slot `s` carries the transition D_{s+1} -> D_{s+2}, which exists only
        // for `s+1 < m`; C zeroes the rest (`(k+z*nq < M) ? ... : -eslINFINITY`, then
        // `esl_sse_expf(-inf) == 0`). Those zero slots are what stops the DD tail from
        // leaking off the end of the model into the E-state sum.
        let q_n = crate::optacc::nqf(m);
        let mut tdd_v = vec![0.0f32; 4 * q_n];
        for q in 0..q_n {
            for z in 0..4 {
                let s = z * q_n + q;
                if s + 1 < m {
                    tdd_v[q * 4 + z] = tdd[s + 1];
                }
            }
        }

        ForwardFilter {
            m,
            q_n,
            tdd_v,
            rfv,
            rfv_t,
            tbm,
            amm,
            aim,
            adm,
            tmi,
            tii,
            tmd,
            tdd,
            ftau: hmm.evparam.fwd_tau,
            flambda: hmm.evparam.fwd_lambda,
        }
    }

    
    
    
    pub fn score(&self, dsq: &[u8], l: usize) -> f32 {
        let m = self.m;

        
        let nj = 1.0f32;
        let pmove = (2.0 + nj) / (l as f32 + 2.0 + nj);
        let ploop = 1.0 - pmove;
        let (xf_e_move, xf_e_loop) = (0.5f32, 0.5f32); 
        let (xf_n_loop, xf_n_move) = (ploop, pmove);
        let (xf_c_loop, xf_c_move) = (ploop, pmove);
        let (xf_j_loop, xf_j_move) = (ploop, pmove);

        
        let mut mp = vec![0.0f32; m + 1];
        let mut ip = vec![0.0f32; m + 1];
        let mut dp = vec![0.0f32; m + 1];
        let mut mc = vec![0.0f32; m + 1];
        let mut ic = vec![0.0f32; m + 1];
        let mut dc = vec![0.0f32; m + 1];
        // Scratch for the striped DD relaxation; `dd_relax` overwrites it every row.
        let mut dd = vec![0.0f32; 4 * self.q_n];

        let mut xn = 1.0f32;
        let mut xj = 0.0f32;
        let mut xb = xf_n_move;
        let mut xc = 0.0f32;
        // C, impl_sse/impl_sse.h:203 -- `totscale` is a *float* field of P7_OMX:
        //   float     totscale;    /* log of the product of all scale factors (0.0 if unscaled)   */
        // so each `ox->totscale += log(xE)` (fwdback.c:432) rounds the running sum
        // back down to f32. Accumulating in f64 and rounding once at the end is a
        // different number.
        let mut totscale = 0.0f32;

        for i in 1..=l {
            let x = dsq[i] as usize;
            mc[0] = 0.0;
            ic[0] = 0.0;
            dc[0] = 0.0;

            
            
            
            
            {
                let base = x * (m + 1);
                let rfv_row = &self.rfv_t[base + 1..base + m + 1]; 
                let mpm = &mp[0..m]; 
                let ipm = &ip[0..m]; 
                let dpm = &dp[0..m]; 
                let tbm = &self.tbm[1..m + 1];
                let amm = &self.amm[1..m + 1];
                let aim = &self.aim[1..m + 1];
                let adm = &self.adm[1..m + 1];
                let mc_o = &mut mc[1..m + 1];
                for j in 0..m {
                    let sv = xb * tbm[j] + mpm[j] * amm[j] + ipm[j] * aim[j] + dpm[j] * adm[j];
                    mc_o[j] = sv * rfv_row[j];
                }
            }
            {
                let mpk = &mp[1..m + 1];
                let ipk = &ip[1..m + 1];
                let tmi = &self.tmi[1..m + 1];
                let tii = &self.tii[1..m + 1];
                let ic_o = &mut ic[1..m + 1];
                for j in 0..m {
                    ic_o[j] = mpk[j] * tmi[j] + ipk[j] * tii[j];
                }
            }

            dd_relax(m, self.q_n, &mc, &self.tmd, &self.tdd_v, &mut dd, &mut dc);

            
            // C, impl_sse/fwdback.c:290-291 (M's, inside the q loop) and 343:
            //   xEv  = _mm_add_ps(xEv, sv);          /* per striped vector q */
            //   ...
            //   /* Add D's to xEv */
            //   for (q = 0; q < Q; q++) xEv = _mm_add_ps(DMO(dpc,q), xEv);
            //   /* horizontal sum of xEv's elements */
            //   xEv = _mm_add_ps(xEv, _mm_shuffle_ps(xEv, xEv, _MM_SHUFFLE(0, 3, 2, 1)));
            //   xEv = _mm_add_ps(xEv, _mm_shuffle_ps(xEv, xEv, _MM_SHUFFLE(1, 0, 3, 2)));
            //   _mm_store_ss(&xE, xEv);
            //
            // Floating-point addition is not associative, so the *order* is part of
            // the transcription: xE accumulates in four striped lanes (lane r takes
            // vector q at model position k = q+1+r*Q, in q order, M's first and then
            // D's) and the horizontal sum collapses as (a0+a1) + (a2+a3). A plain
            // `for k in 1..=m { xe += mc[k] }` is a different order and gives a
            // different float. `omx::forward` already did this correctly.
            let q_n = crate::optacc::nqf(m);
            let mut lane = [0.0f32; 4];
            for (r, l) in lane.iter_mut().enumerate() {
                for q in 0..q_n {
                    let k = r * q_n + q + 1;
                    if k <= m {
                        *l += mc[k];
                    }
                }
            }
            for (r, l) in lane.iter_mut().enumerate() {
                for q in 0..q_n {
                    let k = r * q_n + q + 1;
                    if k <= m {
                        *l += dc[k];
                    }
                }
            }
            let xe = (lane[0] + lane[1]) + (lane[2] + lane[3]);

            
            xn *= xf_n_loop;
            xc = (xc * xf_c_loop) + (xe * xf_e_move);
            xj = (xj * xf_j_loop) + (xe * xf_e_loop);
            xb = (xj * xf_j_move) + (xn * xf_n_move);

            
            // C, impl_sse/fwdback.c:418-433:
            //   if (xE > 1.0e4)      /* that's a little less than e^10, ~10% of our dynamic range */
            //     {
            //       xN  = xN / xE;
            //       xC  = xC / xE;
            //       xJ  = xJ / xE;
            //       xB  = xB / xE;
            //       xEv = _mm_set1_ps(1.0 / xE);
            //       for (q = 0; q < Q; q++)
            //         {
            //           MMO(dpc,q) = _mm_mul_ps(MMO(dpc,q), xEv);
            //           DMO(dpc,q) = _mm_mul_ps(DMO(dpc,q), xEv);
            //           IMO(dpc,q) = _mm_mul_ps(IMO(dpc,q), xEv);
            //         }
            //       ox->xmx[i*p7X_NXCELLS+p7X_SCALE] = xE;
            //       ox->totscale += log(xE);
            //       xE = 1.0;
            //     }
            //
            // The asymmetry is deliberate transcription, not an oversight of ours:
            // the four *specials* are **divided** by xE, while the DP cells are
            // **multiplied** by the reciprocal. `x / e` and `x * (1/e)` round
            // differently, so using one form for both diverges from C at every
            // rescaling event (and long sequences trigger dozens of them).
            if xe > 1.0e4 {
                xn /= xe;
                xc /= xe;
                xj /= xe;
                xb /= xe;
                let inv = 1.0 / xe;
                for k in 0..=m {
                    mc[k] *= inv;
                    dc[k] *= inv;
                    ic[k] *= inv;
                }
                // `log()` is the C double-precision libm call, but the accumulator is
                // f32: `(float)((double)totscale + log((double)xE))`.
                totscale = (totscale as f64 + (xe as f64).ln()) as f32;
            }

            std::mem::swap(&mut mp, &mut mc);
            std::mem::swap(&mut ip, &mut ic);
            std::mem::swap(&mut dp, &mut dc);
        }

        
        // C, impl_sse/fwdback.c:461:
        //   if (opt_sc != NULL) *opt_sc = ox->totscale + log(xC * om->xf[p7O_C][p7O_MOVE]);
        //
        // `xC * om->xf[...]` is an f32 product (both operands are float, and on
        // x86-64 SSE FLT_EVAL_METHOD is 0 so it is not evaluated wider); only then is
        // it widened for `log`. Widening both factors first and multiplying in f64 is
        // a different value.
        (totscale as f64 + ((xc * xf_c_move) as f64).ln()) as f32
    }

    
    
    pub fn pvalue(&self, fwdsc: f32, filtersc: f32) -> f64 {
        let seq_score = (fwdsc - filtersc) as f64 / std::f64::consts::LN_2;
        crate::easel::exponential::esl_exp_surv(seq_score, self.ftau, self.flambda)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bias::BiasFilter;
    use crate::msv::{null_one, MsvProfile};
    use crate::seqio::read_fasta;

    #[test]
    fn forward_filter_passes_all_golden_hits() {
        
        let hmm = P7Hmm::read_all(&format!("{}/testdata/globins4.hmm", env!("CARGO_MANIFEST_DIR")))
            .unwrap()
            .pop()
            .unwrap();
        let mp = MsvProfile::build(&hmm);
        let bf = BiasFilter::build(&hmm).unwrap();
        let ff = ForwardFilter::build(&hmm);
        let seqs = read_fasta(&format!("{}/testdata/globins45.fa", env!("CARGO_MANIFEST_DIR")))
            .unwrap();
        let (f1, f3) = (0.02_f64, 1e-5_f64);
        let mut n_pass = 0;
        for s in &seqs {
            let l = s.len();
            let usc = mp.msv_score(&s.dsq, l);
            if mp.pvalue(usc, null_one(l)) > f1 {
                continue;
            }
            let filtersc = bf.filter_score(&s.dsq, l);
            if mp.pvalue(usc, filtersc) > f1 {
                continue;
            }
            
            let fwdsc = ff.score(&s.dsq, l);
            if ff.pvalue(fwdsc, filtersc) <= f3 {
                n_pass += 1;
            }
        }
        assert_eq!(n_pass, 45, "Fwd filter pass count {n_pass} != C's 45");
    }
}
