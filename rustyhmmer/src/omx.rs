






















use crate::alphabet::{K, KP};
use crate::forward::ForwardFilter;



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




pub(crate) struct XFactors {
    
    pub(crate) n_loop: f32,
    pub(crate) n_move: f32,
    pub(crate) c_loop: f32,
    pub(crate) c_move: f32,
    pub(crate) j_loop: f32,
    pub(crate) j_move: f32,
    
    pub(crate) e_move: f32,
    pub(crate) e_loop: f32,
}

impl XFactors {
    
    
    
    
    
    fn unihit(save_l: usize) -> Self {
        let nj = 0.0f32;
        let pmove = (2.0 + nj) / (save_l as f32 + 2.0 + nj);
        let ploop = 1.0 - pmove;
        XFactors {
            n_loop: ploop,
            n_move: pmove,
            c_loop: ploop,
            c_move: pmove,
            j_loop: ploop,
            j_move: pmove,
            e_move: 1.0,
            e_loop: 0.0,
        }
    }

    
    
    
    
    pub(crate) fn multihit(save_l: usize) -> Self {
        let nj = 1.0f32;
        let pmove = (2.0 + nj) / (save_l as f32 + 2.0 + nj);
        let ploop = 1.0 - pmove;
        XFactors {
            n_loop: ploop,
            n_move: pmove,
            c_loop: ploop,
            c_move: pmove,
            j_loop: ploop,
            j_move: pmove,
            e_move: 0.5,
            e_loop: 0.5,
        }
    }
}



/// One `(ld+1) x (m+1)` DP plane, stored contiguously.
///
/// **Optimization divergence (DEVELOPMENT_PRINCIPLES §1-②, §3):** (a) C keeps the
/// whole matrix in one allocation inside `P7_OMX` (`ox->dpf`, impl_sse/p7_omx.c),
/// striped, and the pipeline reuses it. (b) rustyhmmer had `Vec<Vec<f32>>`, i.e. one
/// allocation *per row* — `3*(Ld+1)` of them per `Omx`, so a 400-residue envelope
/// cost ~1200 allocations and every `mmx[i][k]` paid a double indirection with rows
/// scattered across the heap. This stores each plane as one buffer with stride
/// `m+1`; `Index`/`IndexMut` still hand out `&[f32]` rows, so every `mmx[i][k]` call
/// site is unchanged. (c) Purely a storage change — the same values are written and
/// read in the same order, so the arithmetic is bit-identical by construction;
/// re-verified by --tblout/--domtblout md5 on 397,961 proteins x 117 profiles.
#[derive(Clone)]
pub(crate) struct Rows {
    data: Vec<f32>,
    stride: usize,
}

impl Rows {
    #[inline]
    fn new(rows: usize, stride: usize) -> Self {
        Rows { data: vec![0.0f32; rows * stride], stride }
    }
    /// Two distinct rows, mutably. Replaces the `split_at_mut` dance the
    /// row-of-Vecs layout needed.
    #[inline]
    fn pair_mut(&mut self, a: usize, b: usize) -> (&mut [f32], &mut [f32]) {
        debug_assert!(a < b);
        let (lo, hi) = self.data.split_at_mut(b * self.stride);
        (
            &mut lo[a * self.stride..(a + 1) * self.stride],
            &mut hi[..self.stride],
        )
    }
}

impl std::ops::Index<usize> for Rows {
    type Output = [f32];
    #[inline]
    fn index(&self, i: usize) -> &[f32] {
        &self.data[i * self.stride..(i + 1) * self.stride]
    }
}

impl std::ops::IndexMut<usize> for Rows {
    #[inline]
    fn index_mut(&mut self, i: usize) -> &mut [f32] {
        &mut self.data[i * self.stride..(i + 1) * self.stride]
    }
}

pub(crate) struct Omx {
    pub(crate) m: usize,
    pub(crate) ld: usize,
    
    pub(crate) mmx: Rows,
    pub(crate) imx: Rows,
    pub(crate) dmx: Rows,
    
    pub(crate) xe: Vec<f32>,
    pub(crate) xn: Vec<f32>,
    pub(crate) xj: Vec<f32>,
    pub(crate) xb: Vec<f32>,
    pub(crate) xc: Vec<f32>,
    
    pub(crate) scale: Vec<f32>,
}

impl Omx {
    fn new(m: usize, ld: usize) -> Self {
        Omx {
            m,
            ld,
            mmx: Rows::new(ld + 1, m + 1),
            imx: Rows::new(ld + 1, m + 1),
            dmx: Rows::new(ld + 1, m + 1),
            xe: vec![0.0f32; ld + 1],
            xn: vec![0.0f32; ld + 1],
            xj: vec![0.0f32; ld + 1],
            xb: vec![0.0f32; ld + 1],
            xc: vec![0.0f32; ld + 1],
            scale: vec![1.0f32; ld + 1],
        }
    }
}





pub(crate) fn forward(ff: &ForwardFilter, xf: &XFactors, sub: &[u8], ld: usize) -> Omx {
    let m = ff.m;
    let mut ox = Omx::new(m, ld);

    
    
    ox.xe[0] = 0.0;
    ox.xn[0] = 1.0;
    ox.xj[0] = 0.0;
    ox.xb[0] = xf.n_move; 
    ox.xc[0] = 0.0;
    ox.scale[0] = 1.0;

    
    let mut xn = 1.0f32;
    let mut xj = 0.0f32;
    let mut xb = xf.n_move;
    let mut xc = 0.0f32;

    
    let mut mp = vec![0.0f32; m + 1];
    let mut ip = vec![0.0f32; m + 1];
    let mut dp = vec![0.0f32; m + 1];
    // Scratch for the striped DD relaxation; overwritten every row.
    let mut dd = vec![0.0f32; 4 * ff.q_n];

    for i in 1..=ld {
        let x = sub[i] as usize;
        let mc = &mut ox.mmx[i];
        let ic = &mut ox.imx[i];
        let dc = &mut ox.dmx[i];
        mc[0] = 0.0;
        ic[0] = 0.0;
        dc[0] = 0.0;

        
        
        
        
        
        {
            let base = x * (m + 1);
            let rfv_row = &ff.rfv_t[base + 1..base + m + 1]; 
            let mpm = &mp[0..m]; 
            let ipm = &ip[0..m];
            let dpm = &dp[0..m];
            let tbm = &ff.tbm[1..m + 1];
            let amm = &ff.amm[1..m + 1];
            let aim = &ff.aim[1..m + 1];
            let adm = &ff.adm[1..m + 1];
            let mc_o = &mut mc[1..m + 1];
            for j in 0..m {
                let sv = xb * tbm[j] + mpm[j] * amm[j] + ipm[j] * aim[j] + dpm[j] * adm[j];
                mc_o[j] = sv * rfv_row[j];
            }
        }
        {
            let mpk = &mp[1..m + 1];
            let ipk = &ip[1..m + 1];
            let tmi = &ff.tmi[1..m + 1];
            let tii = &ff.tii[1..m + 1];
            let ic_o = &mut ic[1..m + 1];
            for j in 0..m {
                ic_o[j] = mpk[j] * tmi[j] + ipk[j] * tii[j];
            }
        }

        // Same engine as the F3 filter -- C runs both through `forward_engine()` with
        // only the `do_full` flag differing (impl_sse/fwdback.c:99, :141) -- so the DD
        // paths are relaxed the striped way here too. See `forward::dd_relax`.
        crate::forward::dd_relax(m, ff.q_n, mc, &ff.tmd, &ff.tdd_v, &mut dd, dc);

        
        
        
        
        
        
        
        // `p7O_NQF(M) = ESL_MAX(2, ((((M)-1) / 4) + 1))` (impl_sse/impl_sse.h). The
        // `max(2, ...)` floor matters for M <= 4, where `(m+3)/4` gives 1 vector.
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
        let mut xe = (lane[0] + lane[1]) + (lane[2] + lane[3]);

        
        
        xn *= xf.n_loop;
        xc = (xc * xf.c_loop) + (xe * xf.e_move);
        xj = (xj * xf.j_loop) + (xe * xf.e_loop);
        xb = (xj * xf.j_move) + (xn * xf.n_move);

        
        // C, impl_sse/fwdback.c:418-433. The four *specials* are **divided** by xE
        // while the DP cells are **multiplied** by `1.0/xE`:
        //     xN  = xN / xE;   xC = xC / xE;   xJ = xJ / xE;   xB = xB / xE;
        //     xEv = _mm_set1_ps(1.0 / xE);
        //     for (q = 0; q < Q; q++) { MMO(dpc,q) = _mm_mul_ps(MMO(dpc,q), xEv); ... }
        // `x / e` and `x * (1/e)` round differently, so the asymmetry is part of the
        // transcription.
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
            ox.scale[i] = xe;
            xe = 1.0;
        } else {
            ox.scale[i] = 1.0;
        }

        
        ox.xe[i] = xe;
        ox.xn[i] = xn;
        ox.xj[i] = xj;
        ox.xb[i] = xb;
        ox.xc[i] = xc;

        
        mp.copy_from_slice(&ox.mmx[i]);
        ip.copy_from_slice(&ox.imx[i]);
        dp.copy_from_slice(&ox.dmx[i]);
    }

    ox
}






/// `xB` for one Backward row: the posterior-weighted sum of the B->Mk entries.
///
/// C, impl_sse/fwdback.c:591-597 (inside the `q` loop, `q` descending) and 604-607:
/// ```text
///   mpv        = _mm_mul_ps(MMO(dpp,q), *rp);  rp--;  /* obtain mpv for next q. i+1's MMO(q) is freed  */
///   ...
///   xBv = _mm_add_ps(xBv, _mm_mul_ps(mpv, *tp)); tp--;
///   ...
///   /* this incantation is a horiz sum of xBv elements: (_mm_hadd_ps() would require SSE3) */
///   xBv = _mm_add_ps(xBv, _mm_shuffle_ps(xBv, xBv, _MM_SHUFFLE(0, 3, 2, 1)));
///   xBv = _mm_add_ps(xBv, _mm_shuffle_ps(xBv, xBv, _MM_SHUFFLE(1, 0, 3, 2)));
///   _mm_store_ss(&xB, xBv);
/// ```
///
/// Both the association and the order are part of the value: the emission is folded in
/// first (`(M * e) * tBM`, since `mpv` is already `MMO * rp` when it meets `*tp`), the
/// four striped lanes accumulate independently with `q` running *down*, and the
/// horizontal sum collapses as `(a0+a1) + (a2+a3)`.
///
/// Model positions past `m` are skipped: their emission odds are `esl_sse_expf(-inf)`,
/// i.e. exactly `0.0`, so C's padding lanes contribute an exact zero.
/// `descending` picks the `q` order, which is *not* the same in both places C does this
/// sum and therefore is not the same float:
///
/// * per-residue rows walk `q` down (`for (q = Q-1; q >= 0; q--)`, fwdback.c:583)
/// * the `i=0` termination walks `q` up (`for (q = 0; q < Q; q++)`, fwdback.c:700)
fn bck_xb(ff: &ForwardFilter, x: usize, mm_next: &[f32], descending: bool) -> f32 {
    let (m, q_n) = (ff.m, ff.q_n);
    let mut lane = [0.0f32; 4];
    let acc = |q: usize, lane: &mut [f32; 4]| {
        for (z, l) in lane.iter_mut().enumerate() {
            let k = z * q_n + q + 1;
            if k <= m {
                *l += (mm_next[k] * ff.rfv[k][x]) * ff.tbm[k];
            }
        }
    };
    if descending {
        for q in (0..q_n).rev() {
            acc(q, &mut lane);
        }
    } else {
        for q in 0..q_n {
            acc(q, &mut lane);
        }
    }
    (lane[0] + lane[1]) + (lane[2] + lane[3])
}

/// One Backward row: C's five phases, in C's order.
///
/// C, impl_sse/fwdback.c:563-651. Phase 1 builds `I` and the partial `M`/`D` from row
/// `i+1`; phase 2 does the specials (the caller has already done that, since it needs
/// `xE`); phase 3 adds the `{MD}->E` paths and one DD step; phase 4 finishes the DD
/// paths; phase 5 adds the M->D paths into `M`.
///
/// The order matters twice over. `D` has to be *fully* relaxed before `M` reads
/// `dc[k+1]` in phase 5, and the additions have to be grouped C's way: phase 1 leaves
/// `DMO = mpv*tDM`, then phase 3 stores `DMO + (dcv + xE)` — not `(DMO + dcv) + xE`.
#[allow(clippy::too_many_arguments)]
fn bck_row(
    ff: &ForwardFilter,
    x: usize,
    xe: f32,
    mm_next: &[f32],
    im_next: &[f32],
    me: &mut [f32],
    dd: &mut [f32],
    mc: &mut [f32],
    dc: &mut [f32],
    ic: &mut [f32],
) {
    let (m, q_n) = (ff.m, ff.q_n);

    // `mpv` = M(i+1,k+1) * e(k+1), the one value phases 1 and 5 both need. C keeps it
    // in a register as the q loop walks down; we cache it per k because the D seed is
    // gathered in striped order while M and I are written in k order.
    //
    // C, impl_sse/fwdback.c:580-582 and :596:
    //   mpv = _mm_mul_ps(MMO(dpp,0), om->rfv[dsq[i+1]][0]); /* precalc M(i+1,k+1) * e(M_k+1, x_{i+1}) */
    //   ... mpv = _mm_mul_ps(MMO(dpp,q), *rp);  rp--;
    // At k = m the emission is a padding lane, i.e. exactly 0.0.
    me[m] = 0.0;
    for k in 1..m {
        me[k] = mm_next[k + 1] * ff.rfv[k + 1][x];
    }

    // Phase 1, D only: DMO(dpc,q) = _mm_mul_ps(mpv, tdmv)   (fwdback.c:589)
    for q in 0..q_n {
        for z in 0..4 {
            let k = z * q_n + q + 1;
            dd[q * 4 + z] = if k < m { me[k] * ff.adm[k + 1] } else { 0.0 };
        }
    }
    // Phases 3 and 4.
    crate::forward::bck_dd_relax(q_n, m, xe, &ff.tdd_v, dd, dc);

    // Phase 1 for I and M, plus phase 3's `MMO += xEv` and phase 5's M->D.
    //
    // C, impl_sse/fwdback.c:586-588, :634, :644-649:
    //   IMO(dpc,q) = _mm_add_ps(_mm_mul_ps(ipv, *tp), _mm_mul_ps(mpv, timv));   tp--;
    //   mcv        = _mm_add_ps(_mm_mul_ps(ipv, *tp), _mm_mul_ps(mpv, tmmv));   tp-= 2;
    //   ...
    //   MMO(dpc,q) = _mm_add_ps(MMO(dpc,q), xEv);
    //   ...
    //   MMO(dpc,q) = _mm_add_ps(MMO(dpc,q), _mm_mul_ps(dcv, *tp)); tp -= 7;
    //   dcv        = DMO(dpc,q);
    //
    // giving `((ipv*tMI + mpv*tMM) + xE) + D(k+1)*tMD`. `dcv` there is the fully
    // relaxed `DMO(dpc,q+1)`, which is what `dc[k+1]` is now.
    //
    // At k = m every transition C would use is a padding lane (`tMI`, `tII` and `tMD`
    // are guarded by `kb+z*nq < M`) and `me[m]` is 0, so the row ends at `xE`.
    mc[m] = xe;
    ic[m] = 0.0;
    for k in (1..m).rev() {
        mc[k] = (me[k] * ff.amm[k + 1] + im_next[k] * ff.tmi[k] + xe) + dc[k + 1] * ff.tmd[k];
        ic[k] = me[k] * ff.aim[k + 1] + im_next[k] * ff.tii[k];
    }
}

/// The row-L Backward initialisation: every M and D lane starts at `xE`, I at 0, then
/// the DD paths are relaxed and the M->D paths folded in.
///
/// C, impl_sse/fwdback.c:504-522 (see [`crate::forward::bck_dd_relax`] for the passes):
/// ```text
///   for (q = 0; q < Q; q++) MMO(dpc,q) = DMO(dpc,q) = xEv;
///   for (q = 0; q < Q; q++) IMO(dpc,q) = zerov;
/// ```
fn bck_row_l(
    ff: &ForwardFilter,
    xe: f32,
    dd: &mut [f32],
    mc: &mut [f32],
    dc: &mut [f32],
    ic: &mut [f32],
) {
    let (m, q_n) = (ff.m, ff.q_n);
    dd[..4 * q_n].fill(xe);
    // No `xEv` addend: it is already in every DMO lane.
    crate::forward::bck_dd_relax(q_n, m, 0.0, &ff.tdd_v, dd, dc);
    mc[m] = xe;
    ic[m] = 0.0;
    for k in (1..m).rev() {
        mc[k] = xe + dc[k + 1] * ff.tmd[k];
        ic[k] = 0.0;
    }
}

/// Returns the matrix and whether it ended up using its own scale factors.
///
/// C, impl_sse/fwdback.c:648-662 -- one `backward_engine` serves both `p7_Backward` and
/// `p7_BackwardParser`, so *both* can switch off `fwd`'s scale factors mid-sequence:
/// ```text
///   /* In rare cases [J3/119] scale factors from <fwd> are
///    * insufficient and backwards will overflow. In this case, we
///    * switch on the fly to using our own scale factors, different
///    * from those in <fwd>. This will complicate subsequent
///    * posterior decoding routines.
///    */
///   if (xB > 1.0e16) bck->has_own_scales = TRUE;
///
///   if      (bck->has_own_scales)  bck->xmx[i*p7X_NXCELLS+p7X_SCALE] = (xB > 1.0e4) ? xB : 1.0;
///   else                           bck->xmx[i*p7X_NXCELLS+p7X_SCALE] = fwd->xmx[i*p7X_NXCELLS+p7X_SCALE];
/// ```
/// This used to be implemented only in the domain-definition Backward, which meant the
/// envelope-rescoring one silently kept using `fwd`'s factors past the point where they
/// overflow -- and `decoding` then had no reason to advance `scaleproduct`.
fn backward(
    ff: &ForwardFilter,
    xf: &XFactors,
    sub: &[u8],
    ld: usize,
    fwd: &Omx,
) -> (Omx, bool) {
    let m = ff.m;
    let mut bck = Omx::new(m, ld);
    let mut has_own_scales = false;

    
    let mut xj = 0.0f32;
    let mut xb = 0.0f32;
    let mut xn = 0.0f32;
    let mut xc = xf.c_move; 
    let mut xe = xc * xf.e_move; 
    // Scratch for the striped DD relaxation, reused by every row below.
    let mut dd = vec![0.0f32; 4 * ff.q_n];
    {
        let (mmx, dmx, imx) = (&mut bck.mmx, &mut bck.dmx, &mut bck.imx);
        bck_row_l(ff, xe, &mut dd, &mut mmx[ld], &mut dmx[ld], &mut imx[ld]);
    }
    
    let scl = fwd.scale[ld];
    if scl > 1.0 {
        // C, impl_sse/fwdback.c:535-548 / 664-676: specials divided by the scale
        // factor, DP cells multiplied by its reciprocal.
        xe /= scl;
        xn /= scl;
        xj /= scl;
        xb /= scl;
        xc /= scl;
        let inv = 1.0 / scl;
        for k in 0..=m {
            bck.mmx[ld][k] *= inv;
            bck.dmx[ld][k] *= inv;
            bck.imx[ld][k] *= inv;
        }
    }
    bck.scale[ld] = scl;
    bck.xe[ld] = xe;
    bck.xn[ld] = xn;
    bck.xj[ld] = xj;
    bck.xb[ld] = xb;
    bck.xc[ld] = xc;

    
    // Per-row scratch, hoisted out of the residue loop.
    //
    // **Optimization divergence (DEVELOPMENT_PRINCIPLES §1-②, §3):** (a) C writes the
    // Backward row straight into its preallocated, reused `P7_OMX`
    // (impl_sse/fwdback.c's backward_engine) and never allocates per row. (b) These
    // three vectors were declared *inside* the `for i` loop, so a 400-residue
    // envelope cost 1200 allocations and 1200 memsets per Backward call. (c) Purely
    // where the buffer lives: indices 1..=m are fully rewritten every row (`mc[m]`,
    // `dc[m]`, `ic[m]` then the `k in (1..m).rev()` loop), and index 0 is never
    // written by anyone, so it keeps the 0.0 it is allocated with — exactly as
    // before. The values and their order are untouched.
    let mut mc = vec![0.0f32; m + 1];
    let mut dc = vec![0.0f32; m + 1];
    let mut ic = vec![0.0f32; m + 1];
    let mut me = vec![0.0f32; m + 1];
    for i in (1..ld).rev() {
        let x = sub[i + 1] as usize; 

        xb = bck_xb(ff, x, &bck.mmx[i + 1], true);

        
        
        xc = xc * xf.c_loop;
        xj = (xb * xf.j_move) + (xj * xf.j_loop);
        xn = (xb * xf.n_move) + (xn * xf.n_loop);
        xe = (xc * xf.e_move) + (xj * xf.e_loop);

        
        
        {
            let (mm_next, im_next) = (&bck.mmx[i + 1], &bck.imx[i + 1]);
            bck_row(
                ff, x, xe, mm_next, im_next, &mut me, &mut dd, &mut mc, &mut dc, &mut ic,
            );
            bck.mmx[i].copy_from_slice(&mc);
            bck.dmx[i].copy_from_slice(&dc);
            bck.imx[i].copy_from_slice(&ic);
        }

        
        
        if xb > 1.0e16 {
            has_own_scales = true;
        }
        let scl = if has_own_scales {
            if xb > 1.0e4 {
                xb
            } else {
                1.0
            }
        } else {
            fwd.scale[i]
        };
        if scl > 1.0 {
            // C, impl_sse/fwdback.c:664-676: specials divided, cells multiplied.
            xe /= scl;
            xn /= scl;
            xj /= scl;
            xb /= scl;
            xc /= scl;
            let inv = 1.0 / scl;
            for k in 0..=m {
                bck.mmx[i][k] *= inv;
                bck.dmx[i][k] *= inv;
                bck.imx[i][k] *= inv;
            }
        }
        bck.scale[i] = scl;
        bck.xe[i] = xe;
        bck.xn[i] = xn;
        bck.xj[i] = xj;
        bck.xb[i] = xb;
        bck.xc[i] = xc;
    }

    
    // C, impl_sse/fwdback.c:695-708 -- row 0's B sum is striped exactly like the
    // per-residue one:
    //   for (q = 0; q < Q; q++)
    //     {
    //       mpv = _mm_mul_ps(MMO(dpp,q), *rp);  rp++;
    //       mpv = _mm_mul_ps(mpv,        *tp);  tp += 7;
    //       xBv = _mm_add_ps(xBv,        mpv);
    //     }
    //   /* horizontal sum of xBv */
    //   ...
    // A sequential `sum_k M*e*tBM` is a different float, and `bck.xn[0]` becomes
    // `scaleproduct = 1.0/bck.xn[0]` in the decoding, which scales every posterior.
    let x1 = sub[1] as usize;
    let b0 = bck_xb(ff, x1, &bck.mmx[1], false);
    
    xn = (b0 * xf.n_move) + (xn * xf.n_loop);
    bck.xb[0] = b0;
    bck.xc[0] = 0.0;
    bck.xj[0] = 0.0;
    bck.xn[0] = xn;
    bck.xe[0] = 0.0;
    bck.scale[0] = 1.0;

    (bck, has_own_scales)
}









fn backward_full(
    ff: &ForwardFilter,
    xf: &XFactors,
    sub: &[u8],
    ld: usize,
    fwd: &Omx,
) -> (Omx, bool) {
    let m = ff.m;
    let mut bck = Omx::new(m, ld);
    let mut has_own_scales = false; 

    
    let mut xj = 0.0f32;
    let mut xb = 0.0f32;
    let mut xn = 0.0f32;
    let mut xc = xf.c_move; 
    let mut xe = xc * xf.e_move; 
    // Scratch for the striped DD relaxation, reused by every row below.
    let mut dd = vec![0.0f32; 4 * ff.q_n];
    {
        let (mmx, dmx, imx) = (&mut bck.mmx, &mut bck.dmx, &mut bck.imx);
        bck_row_l(ff, xe, &mut dd, &mut mmx[ld], &mut dmx[ld], &mut imx[ld]);
    }
    
    
    let scl = fwd.scale[ld];
    if scl > 1.0 {
        // C, impl_sse/fwdback.c:535-548 / 664-676: specials divided by the scale
        // factor, DP cells multiplied by its reciprocal.
        xe /= scl;
        xn /= scl;
        xj /= scl;
        xb /= scl;
        xc /= scl;
        let inv = 1.0 / scl;
        for k in 0..=m {
            bck.mmx[ld][k] *= inv;
            bck.dmx[ld][k] *= inv;
            bck.imx[ld][k] *= inv;
        }
    }
    bck.scale[ld] = scl;
    bck.xe[ld] = xe;
    bck.xn[ld] = xn;
    bck.xj[ld] = xj;
    bck.xb[ld] = xb;
    bck.xc[ld] = xc;

    
    // Per-row scratch, hoisted out of the residue loop.
    //
    // **Optimization divergence (DEVELOPMENT_PRINCIPLES §1-②, §3):** (a) C writes the
    // Backward row straight into its preallocated, reused `P7_OMX`
    // (impl_sse/fwdback.c's backward_engine) and never allocates per row. (b) These
    // three vectors were declared *inside* the `for i` loop, so a 400-residue
    // envelope cost 1200 allocations and 1200 memsets per Backward call. (c) Purely
    // where the buffer lives: indices 1..=m are fully rewritten every row (`mc[m]`,
    // `dc[m]`, `ic[m]` then the `k in (1..m).rev()` loop), and index 0 is never
    // written by anyone, so it keeps the 0.0 it is allocated with — exactly as
    // before. The values and their order are untouched.
    let mut mc = vec![0.0f32; m + 1];
    let mut dc = vec![0.0f32; m + 1];
    let mut ic = vec![0.0f32; m + 1];
    let mut me = vec![0.0f32; m + 1];
    for i in (1..ld).rev() {
        let x = sub[i + 1] as usize; 

        xb = bck_xb(ff, x, &bck.mmx[i + 1], true);

        
        xc = xc * xf.c_loop;
        xj = (xb * xf.j_move) + (xj * xf.j_loop);
        xn = (xb * xf.n_move) + (xn * xf.n_loop);
        xe = (xc * xf.e_move) + (xj * xf.e_loop);

        
        {
            let (mm_next, im_next) = (&bck.mmx[i + 1], &bck.imx[i + 1]);
            bck_row(
                ff, x, xe, mm_next, im_next, &mut me, &mut dd, &mut mc, &mut dc, &mut ic,
            );
            bck.mmx[i].copy_from_slice(&mc);
            bck.dmx[i].copy_from_slice(&dc);
            bck.imx[i].copy_from_slice(&ic);
        }

        
        
        
        if xb > 1.0e16 {
            has_own_scales = true;
        }
        let scl = if has_own_scales {
            if xb > 1.0e4 {
                xb
            } else {
                1.0
            }
        } else {
            fwd.scale[i]
        };
        
        if scl > 1.0 {
            // C, impl_sse/fwdback.c:664-676: specials divided, cells multiplied.
            xe /= scl;
            xn /= scl;
            xj /= scl;
            xb /= scl;
            xc /= scl;
            let inv = 1.0 / scl;
            for k in 0..=m {
                bck.mmx[i][k] *= inv;
                bck.dmx[i][k] *= inv;
                bck.imx[i][k] *= inv;
            }
        }
        bck.scale[i] = scl;
        bck.xe[i] = xe;
        bck.xn[i] = xn;
        bck.xj[i] = xj;
        bck.xb[i] = xb;
        bck.xc[i] = xc;
    }

    
    // C, impl_sse/fwdback.c:695-708 -- row 0's B sum is striped exactly like the
    // per-residue one:
    //   for (q = 0; q < Q; q++)
    //     {
    //       mpv = _mm_mul_ps(MMO(dpp,q), *rp);  rp++;
    //       mpv = _mm_mul_ps(mpv,        *tp);  tp += 7;
    //       xBv = _mm_add_ps(xBv,        mpv);
    //     }
    //   /* horizontal sum of xBv */
    //   ...
    // A sequential `sum_k M*e*tBM` is a different float, and `bck.xn[0]` becomes
    // `scaleproduct = 1.0/bck.xn[0]` in the decoding, which scales every posterior.
    let x1 = sub[1] as usize;
    let b0 = bck_xb(ff, x1, &bck.mmx[1], false);
    xn = (b0 * xf.n_move) + (xn * xf.n_loop);
    bck.xb[0] = b0;
    bck.xc[0] = 0.0;
    bck.xj[0] = 0.0;
    bck.xn[0] = xn; 
    bck.xe[0] = 0.0;
    bck.scale[0] = 1.0;

    (bck, has_own_scales)
}







/// `None` on the numeric overflow C reports as `eslERANGE`.
///
/// C, impl_sse/decoding.c:`p7_DomainDecoding`, last line:
/// ```text
///   if (isinf(scaleproduct)) return eslERANGE;
///   else                     return eslOK;
/// ```
/// `p7_domaindef_ByPosteriorHeuristics` returns it immediately (p7_domaindef.c:400) and
/// `p7_Pipeline` turns it into `ESL_FAIL(status, ..., "domain definition workflow
/// failure")` -- which `hmmsearch.c` does not check (it calls `p7_Pipeline` and discards
/// the status, hmmsearch.c:1222/1304). So the observable behaviour is that the target is
/// silently skipped and contributes no hit.
pub(crate) fn domain_decoding(
    ff: &ForwardFilter,
    dsq: &[u8],
    l: usize,
) -> Option<(Vec<f32>, Vec<f32>, Vec<f32>)> {
    
    
    
    let xf = XFactors::multihit(l);

    
    
    
    let fwd = forward(ff, &xf, dsq, l);
    let (bck, has_own_scales) = backward_full(ff, &xf, dsq, l, &fwd);

    
    let mut btot = vec![0.0f32; l + 1];
    let mut etot = vec![0.0f32; l + 1];
    let mut mocc = vec![0.0f32; l + 1];
    
    let mut scaleproduct = 1.0f32 / bck.xn[0];
    
    for i in 1..=l {
        
        btot[i] = btot[i - 1]
            + (fwd.xb[i - 1] * bck.xb[i - 1] * fwd.scale[i - 1] * scaleproduct);
        
        if has_own_scales {
            scaleproduct *= fwd.scale[i - 1] / bck.scale[i - 1];
        }
        
        etot[i] = etot[i - 1] + (fwd.xe[i] * bck.xe[i] * fwd.scale[i] * scaleproduct);
        
        let mut njcp = fwd.xn[i - 1] * bck.xn[i] * xf.n_loop * scaleproduct;
        njcp += fwd.xj[i - 1] * bck.xj[i] * xf.j_loop * scaleproduct;
        njcp += fwd.xc[i - 1] * bck.xc[i] * xf.c_loop * scaleproduct;
        mocc[i] = 1.0 - njcp;
    }
    if scaleproduct.is_infinite() {
        return None;
    }
    Some((btot, etot, mocc))
}





/// Returns `None` on the numeric overflow C reports as `eslERANGE`.
///
/// C, impl_sse/decoding.c:`p7_Decoding`, last line:
/// ```text
///   if (isinf(scaleproduct)) return eslERANGE;
///   else                     return eslOK;
/// ```
/// `rescore_isolated_domain` (p7_domaindef.c:846-852) turns that into `eslFAIL` and
/// abandons the domain: "rare: numeric overflow; domain is assumed to be repetitive
/// garbage [J3/119-121]". The envelope has already been counted by then, so such a
/// target reports `env` greater than `dom`.
fn decoding(xf: &XFactors, fwd: &Omx, bck: &Omx, bck_has_own: bool) -> Option<Omx> {
    let m = fwd.m;
    let ld = fwd.ld;
    let mut pp = Omx::new(m, ld);

    
    let mut scaleproduct = 1.0 / bck.xn[0];

    for i in 1..=ld {
        
        let totrv = scaleproduct * fwd.scale[i];
        
        
        {
            let fm = &fwd.mmx[i][1..m + 1];
            let bm = &bck.mmx[i][1..m + 1];
            let pm = &mut pp.mmx[i][1..m + 1];
            for j in 0..m {
                pm[j] = fm[j] * bm[j] * totrv;
            }
        }
        {
            let fi = &fwd.imx[i][1..m + 1];
            let bi = &bck.imx[i][1..m + 1];
            let pi = &mut pp.imx[i][1..m + 1];
            for j in 0..m {
                pi[j] = fi[j] * bi[j] * totrv;
            }
        }
        
        pp.xe[i] = 0.0;
        pp.xn[i] = fwd.xn[i - 1] * bck.xn[i] * xf.n_loop * scaleproduct;
        pp.xj[i] = fwd.xj[i - 1] * bck.xj[i] * xf.j_loop * scaleproduct;
        pp.xc[i] = fwd.xc[i - 1] * bck.xc[i] * xf.c_loop * scaleproduct;
        pp.xb[i] = 0.0;

        
        if bck_has_own {
            scaleproduct *= fwd.scale[i] / bck.scale[i];
        }
    }

    if scaleproduct.is_infinite() {
        return None;
    }
    Some(pp)
}










/// Everything `rescore_isolated_domain()` (p7_domaindef.c:814-963) extracts from one
/// envelope's Forward/Backward/Decoding/OptimalAccuracy pass.
pub struct EnvRescore {
    /// `envsc` — Forward score of the envelope, in nats.
    pub envsc: f32,
    /// null2 odds ratios from `p7_Null2_ByExpectation()`.
    pub null2: [f32; KP],
    /// `ad->hmmfrom` / `ad->hmmto` — model coordinates of the first/last M state of
    /// the optimal-accuracy alignment (p7_alidisplay.c:178-179).
    pub hmm_from: i32,
    pub hmm_to: i32,
    /// `ad->sqfrom` / `ad->sqto` (p7_alidisplay.c:182-183), **relative to the
    /// envelope** (1..Ld); the caller adds `i-1` to reach original dsq coordinates,
    /// as C does at p7_domaindef.c:857-858.
    pub ali_from: i64,
    pub ali_to: i64,
    /// `oasc` — expected number of correctly aligned residues.
    pub oasc: f32,
    /// The optimal-accuracy traceback, `(state, k, i, pp)` in S..T order with
    /// sequence coordinates relative to the envelope. `p7_alidisplay_Create()`
    /// reads this; `--tblout`/`--domtblout` do not.
    pub trace: Vec<crate::optacc::TracePos>,
}

/// `None` when the posterior decoding overflows -- see [`decoding`]. The caller must
/// drop the domain, as C's `rescore_isolated_domain` does.
pub fn envelope_rescore(
    ff: &ForwardFilter,
    dsq: &[u8],
    i: usize,
    j: usize,
    save_l: usize,
) -> Option<EnvRescore> {



    let xf = XFactors::unihit(save_l);
    let m = ff.m;



    let ld = j - i + 1;
    let mut sub = vec![255u8];
    sub.extend_from_slice(&dsq[i..=j]);
    sub.push(255u8);


    let fwd = forward(ff, &xf, &sub, ld);


    // C, impl_sse/fwdback.c:432 and :461:
    //   ox->totscale += log(xE);                                     /* per rescaling event */
    //   ...
    //   if (opt_sc != NULL) *opt_sc = ox->totscale + log(xC * om->xf[p7O_C][p7O_MOVE]);
    //
    // `ox->totscale` is `float` (impl_sse/impl_sse.h:203), so the running sum is rounded
    // to f32 on every store, and `xC * om->xf[...]` is an f32 product widened only for
    // `log`. This is the score the domain is *reported* with, so both details are
    // load-bearing. Rows that were not rescaled have `scale[r] == 1.0` and contribute
    // an exact `log(1.0) == 0.0`, which is why summing over every row here matches C
    // adding only at the rescaling events.
    let mut totscale = 0.0f32;
    for r in 1..=ld {
        totscale = (totscale as f64 + (fwd.scale[r] as f64).ln()) as f32;
    }
    let envsc = (totscale as f64 + ((fwd.xc[ld] * xf.c_move) as f64).ln()) as f32;

    let (bck, bck_has_own) = backward(ff, &xf, &sub, ld, &fwd);

    let mut pp = decoding(&xf, &fwd, &bck, bck_has_own)?;
    // Nothing below reads the Forward/Backward matrices again — only `pp`. Release
    // them before the optimal-accuracy matrix is allocated so the three don't have to
    // be live at once. (C reuses two preallocated P7_OMXs, `ox1`/`ox2`, for the same
    // reason: p7_domaindef.c:851-857 overwrites ox1 with the OA scores.)
    drop(fwd);
    drop(bck);

    // ---- Optimal accuracy alignment (p7_domaindef.c:856-858):
    //        p7_OptimalAccuracy(om, ox2, ox1, &oasc);
    //        p7_OATrace        (om, ox2, ox1, ddef->tr);
    //      This must run BEFORE the p7_Null2_ByExpectation() block below, which
    //      reuses row 0 of `pp` as its accumulator — C keeps them apart by writing
    //      the OA scores into a separate matrix (ox1) and leaving `pp` (ox2) intact.
    let oatrans = crate::optacc::OaTrans::build(ff);
    let (oax, oasc) = crate::optacc::optimal_accuracy(ff, &xf, &oatrans, &pp, ld);
    let tr = crate::optacc::oa_trace(&oatrans, &xf, &pp, &oax);
    // p7_alidisplay_Create() takes hmmfrom/hmmto/sqfrom/sqto from the first and last
    // M state of the (single, because the model is unilocal here) domain in the
    // trace — exactly what trace_index() computes (p7_alidisplay.c:107-119, 178-183).
    let idx: Vec<(u8, i32, i32)> = tr.iter().map(|&(s, k, i, _)| (s, k, i)).collect();
    let (ali_from, ali_to, hmm_from, hmm_to) = crate::dp::trace_index(&idx)
        .into_iter()
        .find(|d| d.2 != 0)
        .unwrap_or((0, 0, 0, 0));


    
    
    {
        let (r0, r1) = pp.mmx.pair_mut(0, 1);
        r0.copy_from_slice(r1);
    }
    {
        let (r0, r1) = pp.imx.pair_mut(0, 1);
        r0.copy_from_slice(r1);
    }
    pp.xn[0] = pp.xn[1];
    pp.xc[0] = pp.xc[1];
    pp.xj[0] = pp.xj[1];
    
    for r in 2..=ld {
        for k in 0..=m {
            pp.mmx[0][k] += pp.mmx[r][k];
            pp.imx[0][k] += pp.imx[r][k];
        }
        pp.xn[0] += pp.xn[r];
        pp.xc[0] += pp.xc[r];
        pp.xj[0] += pp.xj[r];
    }
    
    let norm = 1.0f32 / ld as f32;
    for k in 0..=m {
        pp.mmx[0][k] *= norm;
        pp.imx[0][k] *= norm;
    }
    pp.xn[0] *= norm;
    pp.xc[0] *= norm;
    pp.xj[0] *= norm;

    
    let xfactor = pp.xn[0] + pp.xc[0] + pp.xj[0];

    
    
    // C, impl_sse/null2.c (`p7_Null2_ByExpectation`):
    //   xfactor = XMXo(0, p7X_N) + XMXo(0, p7X_C) + XMXo(0, p7X_J);
    //   for (x = 0; x < om->abc->K; x++)
    //     {
    //       sv = _mm_setzero_ps();
    //       rp = om->rfv[x];
    //       for (q = 0; q < Q; q++)
    //         {
    //           sv = _mm_add_ps(sv, _mm_mul_ps(pp->dpf[0][q*3 + p7X_M], *rp)); rp++;
    //           sv = _mm_add_ps(sv,            pp->dpf[0][q*3 + p7X_I]);              /* insert odds implicitly 1.0 */
    //         }
    //       esl_sse_hsum_ps(sv, &(null2[x]));
    //       null2[x] += xfactor;
    //     }
    //
    // `esl_sse_hsum_ps` (easel/esl_sse.h:121-126) collapses as `(a0+a1) + (a2+a3)`, so
    // the four striped lanes accumulate independently -- with M and I alternating within
    // each lane, q ascending -- before being combined. A single accumulator over
    // ascending k is a different sum.
    //
    // `pp.imx[0][m]` is included, as C's lane for k = m is: it is exactly 0.0, because
    // I_M's transitions are padding lanes in the Forward matrix it came from.
    let mut null2 = [0.0f32; KP];
    let q_n = ff.q_n;
    for x in 0..K {
        let mut lane = [0.0f32; 4];
        for q in 0..q_n {
            for (z, l) in lane.iter_mut().enumerate() {
                let k = z * q_n + q + 1;
                if k <= m {
                    *l += pp.mmx[0][k] * ff.rfv[k][x];
                    *l += pp.imx[0][k];
                }
            }
        }
        null2[x] = ((lane[0] + lane[1]) + (lane[2] + lane[3])) + xfactor;
    }
    
    for x in (K + 1)..=(KP - 3) {
        let set = degen_set(x);
        let mut result = 0.0f32;
        let mut ndegen = 0.0f32;
        for &y in set {
            result += null2[y];
            ndegen += 1.0;
        }
        null2[x] = if ndegen > 0.0 { result / ndegen } else { 0.0 };
    }
    null2[K] = 1.0;
    null2[KP - 2] = 1.0;
    null2[KP - 1] = 1.0;
    Some(EnvRescore {
        envsc,
        null2,
        hmm_from,
        hmm_to,
        ali_from: ali_from as i64,
        ali_to: ali_to as i64,
        oasc,
        trace: tr,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hmmfile::P7Hmm;
    use crate::seqio::read_fasta;

    #[test]
    fn envelope_null2_smoke() {
        let hmm = P7Hmm::read_all(&format!("{}/testdata/globins4.hmm", env!("CARGO_MANIFEST_DIR")))
            .unwrap()
            .pop()
            .unwrap();
        let seqs = read_fasta(&format!("{}/testdata/globins45.fa", env!("CARGO_MANIFEST_DIR")))
            .unwrap();
        let s = seqs
            .iter()
            .find(|s| s.name == "MYG_ESCGI")
            .expect("MYG_ESCGI not found");
        let l = s.len();
        assert_eq!(l, 153, "MYG_ESCGI length");

        
        let (i, j, save_l) = (1usize, 147usize, 153usize);
        let ff = ForwardFilter::build(&hmm);
        let null2 = envelope_rescore(&ff, &s.dsq, i, j, save_l).unwrap().null2;

        
        for x in 0..K {
            assert!(null2[x].is_finite(), "null2[{x}] not finite: {}", null2[x]);
            assert!(null2[x] > 0.0, "null2[{x}] not > 0: {}", null2[x]);
        }

        
        let mut domcorrection = 0.0f32;
        for pos in i..=j {
            domcorrection += null2[s.dsq[pos] as usize].ln();
        }
        println!("domcorrection (odds-space null2) = {domcorrection}");
        assert!(
            domcorrection > 6.5 && domcorrection < 8.5,
            "domcorrection {domcorrection} out of sane range"
        );
    }

    /// The OA traceback's per-position posterior probabilities must sum back to the
    /// OA score: `p7_trace_GetExpectedAccuracy(tr) == oasc`, the invariant C's own
    /// optacc unit test checks (impl_sse/optacc.c:670). This is what validates the
    /// `get_postprob()` transcription, since `oasc` comes from the independent DP
    /// fill (`XMXo(L, p7X_C)`).
    #[test]
    fn oa_trace_pp_sums_to_oasc() {
        use crate::alphabet::AminoMap;
        use crate::seqio::read_fasta;
        let base = env!("CARGO_MANIFEST_DIR");
        let hmm = P7Hmm::read_all(&format!("{base}/testdata/globins4.hmm"))
            .unwrap()
            .pop()
            .unwrap();
        let ff = ForwardFilter::build(&hmm);
        let seqs = read_fasta(&format!("{base}/testdata/globins45.fa")).unwrap();
        let _ = AminoMap::new();

        let mut n = 0usize;
        for s in seqs.iter().take(12) {
            let l = s.len();
            // A handful of envelopes per sequence, including the full-length one.
            for &(i, j) in &[(1usize, l), (1, l / 2 + 1), (l / 4 + 1, l)] {
                if j <= i {
                    continue;
                }
                let er = envelope_rescore(&ff, &s.dsq, i, j, l).expect("no decoding overflow expected here");
                let acc = crate::optacc::expected_accuracy(&er.trace);
                let tol = 1e-3 * er.oasc.abs().max(1.0);
                assert!(
                    (acc - er.oasc).abs() <= tol,
                    "{}: sum(tr.pp)={acc} != oasc={} (env {i}..{j})",
                    s.name,
                    er.oasc
                );
                // The trace must start at S and end at T.
                assert_eq!(er.trace.first().unwrap().0, crate::dp::ST_S);
                assert_eq!(er.trace.last().unwrap().0, crate::dp::ST_T);
                n += 1;
            }
        }
        assert!(n >= 24, "expected to check many envelopes, checked {n}");
    }
}
