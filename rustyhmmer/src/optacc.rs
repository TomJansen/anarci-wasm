//! Optimal accuracy alignment (sensu Kall05) — DP fill and traceback.
//!
//! Faithful transcription of HMMER 3.4 `src/impl_sse/optacc.c`, which is the
//! implementation `hmmsearch` actually runs (`rescore_isolated_domain()` calls
//! `p7_OptimalAccuracy()`/`p7_OATrace()`, not the `generic_optacc.c` pair):
//!
//! | Rust (this file)      | C (hmmer-3.4)                                    |
//! |-----------------------|--------------------------------------------------|
//! | `OaTrans::build`      | `fb_conversion()`   impl_sse/p7_oprofile.c:929-981 |
//! | `optimal_accuracy`    | `p7_OptimalAccuracy()` impl_sse/optacc.c:59-166  |
//! | `oa_trace`            | `p7_OATrace()`      impl_sse/optacc.c:216-262    |
//! | `select_m`            | `select_m()`        impl_sse/optacc.c:283-317    |
//! | `select_d`            | `select_d()`        impl_sse/optacc.c:320-344    |
//! | `select_i`            | `select_i()`        impl_sse/optacc.c:347-360    |
//! | `select_n`            | `select_n()`        impl_sse/optacc.c:363-367    |
//! | `select_c`            | `select_c()`        impl_sse/optacc.c:370-377    |
//! | `select_j`            | `select_j()`        impl_sse/optacc.c:380-388    |
//! | `select_e`            | `select_e()`        impl_sse/optacc.c:392-412    |
//! | `select_b`            | `select_b()`        impl_sse/optacc.c:416-424    |
//!
//! The striped (Farrar) *layout* is reproduced literally rather than collapsed to
//! a sequential k=1..M recurrence. That is deliberate (DEVELOPMENT_PRINCIPLES §2):
//! the SSE code's observable behaviour depends on the layout in ways a "obviously
//! equivalent" sequential rewrite would not reproduce —
//!   (a) lanes for k>M are padding whose transitions are `expf(-inf) == 0.0`, so
//!       they contribute `0.0` (not -inf) to `sv`, and therefore participate in
//!       both the `xEv` reduction and `select_e()`'s argmax over k;
//!   (b) the DP fill shifts `-eslINFINITY` on (`esl_sse_rightshift_ps`) while the
//!       traceback's `select_m`/`select_d` shift `0.0` on
//!       (`esl_sse_rightshiftz_float`) — a real asymmetry in the C;
//!   (c) the lazy-F D->D resolution runs one "collect from DMO" pass followed by
//!       three "propagate dcv" passes, which is not the plain sequential
//!       `D(i,k) = max(D(i,k), D(i,k-1))` recurrence.
//! Each of these is transcribed as-is; no equivalence argument is relied upon.
//!
//! A "vector" below is 4 consecutive `f32` lanes, matching `__m128`.

use crate::dp::{ST_B, ST_C, ST_D, ST_E, ST_I, ST_J, ST_M, ST_N, ST_S, ST_T};
use crate::forward::ForwardFilter;
use crate::omx::{Omx, XFactors};

const NEG: f32 = f32::NEG_INFINITY;

/// `p7X_M` / `p7X_D` / `p7X_I` — `enum p7x_scells_e`, impl_sse/impl_sse.h:167.
const PX_M: usize = 0;
const PX_D: usize = 1;
const PX_I: usize = 2;

/// `p7O_NQF(M)` — impl_sse/impl_sse.h:26. Number of 4-wide f32 vectors per row.
#[inline]
pub(crate) fn nqf(m: usize) -> usize {
    debug_assert!(m >= 1, "p7O_NQF is only defined for M >= 1");
    std::cmp::max(2, (m - 1) / 4 + 1)
}

/* ---- scalar transcriptions of the SSE intrinsics used by optacc.c ---- */

/// `_mm_max_ps(a, b)` — per-lane `a > b ? a : b`.
#[inline]
fn v_max(a: [f32; 4], b: [f32; 4]) -> [f32; 4] {
    let mut r = [0.0f32; 4];
    for z in 0..4 {
        r[z] = if a[z] > b[z] { a[z] } else { b[z] };
    }
    r
}

/// `_mm_and_ps(_mm_cmpgt_ps(t, zerov), v)` — keep `v` where `t > 0`, else `+0.0`
/// (the IEEE754 all-zero bit pattern the mask leaves behind).
#[inline]
fn v_and_gt0(t: &[f32], v: [f32; 4]) -> [f32; 4] {
    let mut r = [0.0f32; 4];
    for z in 0..4 {
        r[z] = if t[z] > 0.0 { v[z] } else { 0.0 };
    }
    r
}

/// `_mm_add_ps(a, b)`
#[inline]
fn v_add(a: [f32; 4], b: [f32; 4]) -> [f32; 4] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2], a[3] + b[3]]
}

/// `esl_sse_rightshift_ps(a, b)` — easel/esl_sse.h:212. Returns `{b[0],a[0],a[1],a[2]}`;
/// call sites pass `infv`, so `-eslINFINITY` is shifted on.
#[inline]
fn rightshift(a: [f32; 4], b: f32) -> [f32; 4] {
    [b, a[0], a[1], a[2]]
}

/// `esl_sse_rightshiftz_float(a)` — easel/esl_sse.h:181. Shifts a **zero** on, not -inf.
#[inline]
fn rightshiftz(a: [f32; 4]) -> [f32; 4] {
    [0.0, a[0], a[1], a[2]]
}

/// `esl_sse_hmax_ps(a, &max)` — easel/esl_sse.h:93. Two shuffle/max rounds, then lane 0.
#[inline]
fn hmax(a: [f32; 4]) -> f32 {
    // a = _mm_max_ps(a, _mm_shuffle_ps(a, a, _MM_SHUFFLE(0, 3, 2, 1)));
    let s = [a[1], a[2], a[3], a[0]];
    let a = v_max(a, s);
    // a = _mm_max_ps(a, _mm_shuffle_ps(a, a, _MM_SHUFFLE(1, 0, 3, 2)));
    let s = [a[2], a[3], a[0], a[1]];
    let a = v_max(a, s);
    a[0]
}

/// `esl_vec_FArgMax()` — easel/esl_vectorops.c:403. First maximum wins ties.
#[inline]
fn argmax(path: &[f32]) -> usize {
    let mut best = 0usize;
    for i in 1..path.len() {
        if path[i] > path[best] {
            best = i;
        }
    }
    best
}

/* ------------------------- striped transitions ------------------------- */

/// Striped Forward/Backward transition probabilities, i.e. C's `om->tfv`.
///
/// Faithful port of `fb_conversion()` (impl_sse/p7_oprofile.c:953-980), except that
/// rustyhmmer's `ForwardFilter` already stores odds ratios (= `expf` of `gm`'s tsc),
/// so the `esl_sse_expf()` step is already done and an out-of-range lane is `0.0`
/// (= `expf(-eslINFINITY)`) rather than `-eslINFINITY`.
pub(crate) struct OaTrans {
    q: usize,
    /// 7 interleaved vectors per q — [BM,MM,IM,DM,MD,MI,II] at vector `7*q+t` —
    /// followed by the DD block of Q vectors at `7*Q+q`. Vector `v` is `tfv[v*4..v*4+4]`.
    tfv: Vec<f32>,
}

/// `enum p7o_tsc_e` — impl_sse/impl_sse.h:73.
const O_BM: usize = 0;
const O_MM: usize = 1;
const O_IM: usize = 2;
const O_DM: usize = 3;
const O_MD: usize = 4;
const O_MI: usize = 5;
const O_II: usize = 6;

impl OaTrans {
    pub(crate) fn build(ff: &ForwardFilter) -> Self {
        let m = ff.m;
        let q = nqf(m);
        let mut tfv = vec![0.0f32; 8 * q * 4];

        // "Transition scores, all but the DD's" (p7_oprofile.c:953-973).
        // BM/MM/IM/DM read gm at kb = k-1 (gm stores those rotated by -1), so in
        // rustyhmmer's 1-based arrays they are index kk = k + z*nq, valid for kk <= M.
        // MD/MI/II read gm "straight up" at kb = k, in range while kb + z*nq < M.
        for qi in 0..q {
            let k = qi + 1;
            for t in 0..7 {
                let v = qi * 7 + t;
                for z in 0..4 {
                    let kk = k + z * q;
                    tfv[v * 4 + z] = match t {
                        O_BM => {
                            if kk <= m {
                                ff.tbm[kk]
                            } else {
                                0.0
                            }
                        }
                        O_MM => {
                            if kk <= m {
                                ff.amm[kk]
                            } else {
                                0.0
                            }
                        }
                        O_IM => {
                            if kk <= m {
                                ff.aim[kk]
                            } else {
                                0.0
                            }
                        }
                        O_DM => {
                            if kk <= m {
                                ff.adm[kk]
                            } else {
                                0.0
                            }
                        }
                        O_MD => {
                            if kk < m {
                                ff.tmd[kk]
                            } else {
                                0.0
                            }
                        }
                        O_MI => {
                            if kk < m {
                                ff.tmi[kk]
                            } else {
                                0.0
                            }
                        }
                        _ => {
                            debug_assert_eq!(t, O_II);
                            if kk < m {
                                ff.tii[kk]
                            } else {
                                0.0
                            }
                        }
                    };
                }
            }
        }

        // "And finally the DD's, which are at the end of the optimized tfv vector"
        // (p7_oprofile.c:976-980).
        for qi in 0..q {
            let k = qi + 1;
            let v = 7 * q + qi;
            for z in 0..4 {
                let kk = k + z * q;
                tfv[v * 4 + z] = if kk < m { ff.tdd[kk] } else { 0.0 };
            }
        }

        OaTrans { q, tfv }
    }

    /// Vector `v` of `om->tfv`, as 4 lanes.
    #[inline]
    fn v(&self, v: usize) -> &[f32] {
        &self.tfv[v * 4..v * 4 + 4]
    }

    #[inline]
    fn vv(&self, v: usize) -> [f32; 4] {
        let s = self.v(v);
        [s[0], s[1], s[2], s[3]]
    }
}

/* ---------------------------- the OA matrix ---------------------------- */

/// Striped OA score matrix — C's `P7_OMX` used as `ox1` by `rescore_isolated_domain()`.
pub(crate) struct OaMx {
    q: usize,
    l: usize,
    /// Row `i` holds `3*Q` vectors; cell `(q, s)` is vector `q*3 + s`, matching the
    /// `MMO`/`DMO`/`IMO` macros (impl_sse/impl_sse.h:222-224).
    dp: Vec<f32>,
    /// `XMXo(i, s)` for s = E,N,J,B,C (impl_sse/impl_sse.h:171).
    xe: Vec<f32>,
    xn: Vec<f32>,
    xj: Vec<f32>,
    xb: Vec<f32>,
    xc: Vec<f32>,
}

impl OaMx {
    fn new(q: usize, l: usize) -> Self {
        OaMx {
            q,
            l,
            dp: vec![0.0f32; (l + 1) * 3 * q * 4],
            xe: vec![0.0f32; l + 1],
            xn: vec![0.0f32; l + 1],
            xj: vec![0.0f32; l + 1],
            xb: vec![0.0f32; l + 1],
            xc: vec![0.0f32; l + 1],
        }
    }

    #[inline]
    fn off(&self, i: usize, qi: usize, s: usize) -> usize {
        ((i * self.q + qi) * 3 + s) * 4
    }

    /// Read cell `(i, q, s)` as 4 lanes.
    #[inline]
    fn get(&self, i: usize, qi: usize, s: usize) -> [f32; 4] {
        let o = self.off(i, qi, s);
        [self.dp[o], self.dp[o + 1], self.dp[o + 2], self.dp[o + 3]]
    }

    #[inline]
    fn set(&mut self, i: usize, qi: usize, s: usize, v: [f32; 4]) {
        let o = self.off(i, qi, s);
        self.dp[o..o + 4].copy_from_slice(&v);
    }
}

/* --------------------------- 1. OA DP fill ----------------------------- */

/// `p7_OptimalAccuracy()` — impl_sse/optacc.c:59-166.
///
/// `pp` is the posterior decoding matrix (C's `ox2` after `p7_Decoding`), `ld` its
/// length `pp->L`. Returns the filled OA matrix and `*ret_e` (the "oasc": expected
/// number of correctly decoded positions).
pub(crate) fn optimal_accuracy(
    ff: &ForwardFilter,
    xf: &XFactors,
    t: &OaTrans,
    pp: &Omx,
    ld: usize,
) -> (OaMx, f32) {
    let m = ff.m;
    let q = t.q;
    let mut ox = OaMx::new(q, ld);
    let infv = [NEG; 4];

    // Posterior probabilities for the 4 lanes of vector `qi` on row `i`. Lanes with
    // k > M are padding: C's pp holds fwd*bck*totrv there, and both fwd and bck
    // padding lanes are 0.0, so the product is 0.0.
    let ppm = |i: usize, qi: usize| -> [f32; 4] {
        let mut v = [0.0f32; 4];
        for (z, slot) in v.iter_mut().enumerate() {
            let k = qi + 1 + z * q;
            if k <= m {
                *slot = pp.mmx[i][k];
            }
        }
        v
    };
    let ppi = |i: usize, qi: usize| -> [f32; 4] {
        let mut v = [0.0f32; 4];
        for (z, slot) in v.iter_mut().enumerate() {
            let k = qi + 1 + z * q;
            if k <= m {
                *slot = pp.imx[i][k];
            }
        }
        v
    };

    // for (q = 0; q < Q; q++) MMO(dpc,q) = IMO(dpc,q) = DMO(dpc,q) = infv;
    for qi in 0..q {
        ox.set(0, qi, PX_M, infv);
        ox.set(0, qi, PX_I, infv);
        ox.set(0, qi, PX_D, infv);
    }
    ox.xe[0] = NEG;
    ox.xn[0] = 0.0;
    ox.xj[0] = NEG;
    ox.xb[0] = 0.0;
    ox.xc[0] = NEG;

    // C, impl_sse/optacc.c:88-125:
    //   for (i = 1; i <= pp->L; i++)
    //     {
    //       dpp = dpc;             /* previous DP row in OA matrix */
    //       dpc = ox->dpf[i];      /* current DP row in OA matrix  */
    //       ppp = pp->dpf[i];      /* current row in the posterior probabilities  */
    //       tp  = om->tfv;         /* transition probabilities */
    //       dcv = infv;
    //       xEv = infv;
    //       xBv = _mm_set1_ps(XMXo(i-1, p7X_B));
    //
    //       mpv = esl_sse_rightshift_ps(MMO(dpp,Q-1), infv);
    //       dpv = esl_sse_rightshift_ps(DMO(dpp,Q-1), infv);
    //       ipv = esl_sse_rightshift_ps(IMO(dpp,Q-1), infv);
    //       for (q = 0; q < Q; q++)
    //         {
    //           sv  =                _mm_and_ps(_mm_cmpgt_ps(*tp, zerov), xBv);  tp++;
    //           sv  = _mm_max_ps(sv, _mm_and_ps(_mm_cmpgt_ps(*tp, zerov), mpv)); tp++;
    //           sv  = _mm_max_ps(sv, _mm_and_ps(_mm_cmpgt_ps(*tp, zerov), ipv)); tp++;
    //           sv  = _mm_max_ps(sv, _mm_and_ps(_mm_cmpgt_ps(*tp, zerov), dpv)); tp++;
    //           sv  = _mm_add_ps(sv, *ppp);                                      ppp += 2;
    //           xEv = _mm_max_ps(xEv, sv);
    //
    //           mpv = MMO(dpp,q);
    //           dpv = DMO(dpp,q);
    //           ipv = IMO(dpp,q);
    //
    //           MMO(dpc,q) = sv;
    //           DMO(dpc,q) = dcv;
    //
    //           dcv = _mm_and_ps(_mm_cmpgt_ps(*tp, zerov), sv); tp++;
    //
    //           sv         =                _mm_and_ps(_mm_cmpgt_ps(*tp, zerov), mpv);   tp++;
    //           sv         = _mm_max_ps(sv, _mm_and_ps(_mm_cmpgt_ps(*tp, zerov), ipv));  tp++;
    //           IMO(dpc,q) = _mm_add_ps(sv, *ppp);                                       ppp++;
    //         }
    for i in 1..=ld {
        let mut dcv = infv;
        let mut xev = infv;
        let xbv = [ox.xb[i - 1]; 4];

        // Right shifts by 4 bytes. 4,8,12,x becomes x,4,8,12.
        let mut mpv = rightshift(ox.get(i - 1, q - 1, PX_M), NEG);
        let mut dpv = rightshift(ox.get(i - 1, q - 1, PX_D), NEG);
        let mut ipv = rightshift(ox.get(i - 1, q - 1, PX_I), NEG);

        for qi in 0..q {
            let tp = qi * 7;
            let mut sv = v_and_gt0(t.v(tp + O_BM), xbv);
            sv = v_max(sv, v_and_gt0(t.v(tp + O_MM), mpv));
            sv = v_max(sv, v_and_gt0(t.v(tp + O_IM), ipv));
            sv = v_max(sv, v_and_gt0(t.v(tp + O_DM), dpv));
            sv = v_add(sv, ppm(i, qi));
            xev = v_max(xev, sv);

            mpv = ox.get(i - 1, qi, PX_M);
            dpv = ox.get(i - 1, qi, PX_D);
            ipv = ox.get(i - 1, qi, PX_I);

            ox.set(i, qi, PX_M, sv);
            ox.set(i, qi, PX_D, dcv);

            dcv = v_and_gt0(t.v(tp + O_MD), sv);

            let mut sv2 = v_and_gt0(t.v(tp + O_MI), mpv);
            sv2 = v_max(sv2, v_and_gt0(t.v(tp + O_II), ipv));
            ox.set(i, qi, PX_I, v_add(sv2, ppi(i, qi)));
        }

        // C, impl_sse/optacc.c:126-146:
        //   /* dcv has carried through from end of q loop above; store it
        //    * in first pass, we add M->D and D->D path into DMX
        //    */
        //   dcv = esl_sse_rightshift_ps(dcv, infv);
        //   tp  = om->tfv + 7*Q;        /* set tp to start of the DD's */
        //   for (q = 0; q < Q; q++)
        //     {
        //       DMO(dpc, q) = _mm_max_ps(dcv, DMO(dpc, q));
        //       dcv         = _mm_and_ps(_mm_cmpgt_ps(*tp, zerov), DMO(dpc,q));   tp++;
        //     }
        //
        //   /* fully serialized D->D; can optimize later */
        //   for (j = 1; j < 4; j++)
        //     {
        //       dcv = esl_sse_rightshift_ps(dcv, infv);
        //       tp  = om->tfv + 7*Q;
        //       for (q = 0; q < Q; q++)
        //         {
        //           DMO(dpc, q) = _mm_max_ps(dcv, DMO(dpc, q));
        //           dcv         = _mm_and_ps(_mm_cmpgt_ps(*tp, zerov), dcv);   tp++;
        //         }
        //     }
        // NB the second and later passes propagate `dcv`, not DMO -- that is C's,
        // not a typo here.
        dcv = rightshift(dcv, NEG);
        for qi in 0..q {
            let d = v_max(dcv, ox.get(i, qi, PX_D));
            ox.set(i, qi, PX_D, d);
            dcv = v_and_gt0(t.v(7 * q + qi), d);
        }

        // Fully serialized D->D.
        for _j in 1..4 {
            dcv = rightshift(dcv, NEG);
            for qi in 0..q {
                let d = v_max(dcv, ox.get(i, qi, PX_D));
                ox.set(i, qi, PX_D, d);
                dcv = v_and_gt0(t.v(7 * q + qi), dcv);
            }
        }

        // D->E paths
        for qi in 0..q {
            xev = v_max(xev, ox.get(i, qi, PX_D));
        }

        // C, impl_sse/optacc.c:148-166:
        //   /* D->E paths */
        //   for (q = 0; q < Q; q++) xEv = _mm_max_ps(xEv, DMO(dpc,q));
        //
        //   /* Specials */
        //   esl_sse_hmax_ps(xEv, &(XMXo(i,p7X_E)));
        //
        //   t1 = ( (om->xf[p7O_J][p7O_LOOP] == 0.0) ? 0.0 : ox->xmx[(i-1)*p7X_NXCELLS+p7X_J] + pp->xmx[i*p7X_NXCELLS+p7X_J]);
        //   t2 = ( (om->xf[p7O_E][p7O_LOOP] == 0.0) ? 0.0 : ox->xmx[   i *p7X_NXCELLS+p7X_E]);
        //   ox->xmx[i*p7X_NXCELLS+p7X_J] = ESL_MAX(t1, t2);
        //
        //   t1 = ( (om->xf[p7O_C][p7O_LOOP] == 0.0) ? 0.0 : ox->xmx[(i-1)*p7X_NXCELLS+p7X_C] + pp->xmx[i*p7X_NXCELLS+p7X_C]);
        //   t2 = ( (om->xf[p7O_E][p7O_MOVE] == 0.0) ? 0.0 : ox->xmx[   i *p7X_NXCELLS+p7X_E]);
        //   ox->xmx[i*p7X_NXCELLS+p7X_C] = ESL_MAX(t1, t2);
        //
        //   ox->xmx[i*p7X_NXCELLS+p7X_N] = ((om->xf[p7O_N][p7O_LOOP] == 0.0) ? 0.0 : ox->xmx[(i-1)*p7X_NXCELLS+p7X_N] + pp->xmx[i*p7X_NXCELLS+p7X_N]);
        //
        //   t1 = ( (om->xf[p7O_N][p7O_MOVE] == 0.0) ? 0.0 : ox->xmx[i*p7X_NXCELLS+p7X_N]);
        //   t2 = ( (om->xf[p7O_J][p7O_MOVE] == 0.0) ? 0.0 : ox->xmx[i*p7X_NXCELLS+p7X_J]);
        //   ox->xmx[i*p7X_NXCELLS+p7X_B] = ESL_MAX(t1, t2);
        ox.xe[i] = hmax(xev);

        let t1 = if xf.j_loop == 0.0 {
            0.0
        } else {
            ox.xj[i - 1] + pp.xj[i]
        };
        let t2 = if xf.e_loop == 0.0 { 0.0 } else { ox.xe[i] };
        ox.xj[i] = if t1 > t2 { t1 } else { t2 };

        let t1 = if xf.c_loop == 0.0 {
            0.0
        } else {
            ox.xc[i - 1] + pp.xc[i]
        };
        let t2 = if xf.e_move == 0.0 { 0.0 } else { ox.xe[i] };
        ox.xc[i] = if t1 > t2 { t1 } else { t2 };

        ox.xn[i] = if xf.n_loop == 0.0 {
            0.0
        } else {
            ox.xn[i - 1] + pp.xn[i]
        };

        let t1 = if xf.n_move == 0.0 { 0.0 } else { ox.xn[i] };
        let t2 = if xf.j_move == 0.0 { 0.0 } else { ox.xj[i] };
        ox.xb[i] = if t1 > t2 { t1 } else { t2 };
    }

    let ret_e = ox.xc[ld];
    (ox, ret_e)
}

/* --------------------------- 2. OA traceback --------------------------- */

/// M(i,k) is reached from B(i-1), M(i-1,k-1), D(i-1,k-1), or I(i-1,k-1).
/// `select_m()` — impl_sse/optacc.c:283-317.
// C, impl_sse/optacc.c:284-317:
//   static inline int
//   select_m(const P7_OPROFILE *om, const P7_OMX *ox, int i, int k)
//   {
//     int     Q     = p7O_NQF(ox->M);
//     int     q     = (k-1) % Q;
//     int     r     = (k-1) / Q;
//     __m128 *tp    = om->tfv + 7*q;
//     __m128  xBv   = _mm_set1_ps(ox->xmx[(i-1)*p7X_NXCELLS+p7X_B]);
//     __m128  mpv, dpv, ipv;
//     union { __m128 v; float p[4]; } u, tv;
//     float   path[4];
//     int     state[4] = { p7T_M, p7T_I, p7T_D, p7T_B };
//
//     if (q > 0) {
//       mpv = ox->dpf[i-1][(q-1)*3 + p7X_M];
//       dpv = ox->dpf[i-1][(q-1)*3 + p7X_D];
//       ipv = ox->dpf[i-1][(q-1)*3 + p7X_I];
//     } else {
//       mpv = esl_sse_rightshiftz_float(ox->dpf[i-1][(Q-1)*3 + p7X_M]);
//       dpv = esl_sse_rightshiftz_float(ox->dpf[i-1][(Q-1)*3 + p7X_D]);
//       ipv = esl_sse_rightshiftz_float(ox->dpf[i-1][(Q-1)*3 + p7X_I]);
//     }
//
//     /* paths are numbered so that most desirable choice in case of tie is first. */
//     u.v = xBv;  tv.v = *tp;  path[3] = ((tv.p[r] == 0.0) ?  -eslINFINITY : u.p[r]);  tp++;
//     u.v = mpv;  tv.v = *tp;  path[0] = ((tv.p[r] == 0.0) ?  -eslINFINITY : u.p[r]);  tp++;
//     u.v = ipv;  tv.v = *tp;  path[1] = ((tv.p[r] == 0.0) ?  -eslINFINITY : u.p[r]);  tp++;
//     u.v = dpv;  tv.v = *tp;  path[2] = ((tv.p[r] == 0.0) ?  -eslINFINITY : u.p[r]);
//     return state[esl_vec_FArgMax(path, 4)];
//   }
// NB `esl_sse_rightshiftz_float` shifts a ZERO on, where the DP fill above shifts
// -eslINFINITY on. That asymmetry is C's.
fn select_m(t: &OaTrans, ox: &OaMx, i: usize, k: usize) -> u8 {
    let q = ox.q;
    let qi = (k - 1) % q;
    let r = (k - 1) / q;
    let tp = 7 * qi;
    let xbv = [ox.xb[i - 1]; 4];

    let (mpv, dpv, ipv) = if qi > 0 {
        (
            ox.get(i - 1, qi - 1, PX_M),
            ox.get(i - 1, qi - 1, PX_D),
            ox.get(i - 1, qi - 1, PX_I),
        )
    } else {
        // NB: `esl_sse_rightshiftz_float` shifts a ZERO on here, unlike the DP fill.
        (
            rightshiftz(ox.get(i - 1, q - 1, PX_M)),
            rightshiftz(ox.get(i - 1, q - 1, PX_D)),
            rightshiftz(ox.get(i - 1, q - 1, PX_I)),
        )
    };

    // Paths are numbered so that most desirable choice in case of tie is first.
    let mut path = [0.0f32; 4];
    path[3] = if t.v(tp + O_BM)[r] == 0.0 { NEG } else { xbv[r] };
    path[0] = if t.v(tp + O_MM)[r] == 0.0 { NEG } else { mpv[r] };
    path[1] = if t.v(tp + O_IM)[r] == 0.0 { NEG } else { ipv[r] };
    path[2] = if t.v(tp + O_DM)[r] == 0.0 { NEG } else { dpv[r] };
    const STATE: [u8; 4] = [ST_M, ST_I, ST_D, ST_B];
    STATE[argmax(&path)]
}

/// D(i,k) is reached from M(i,k-1) or D(i,k-1). `select_d()` — impl_sse/optacc.c:320-344.
// C, impl_sse/optacc.c:321-344:
//   if (q > 0) {
//     mpv.v  = ox->dpf[i][(q-1)*3 + p7X_M];
//     dpv.v  = ox->dpf[i][(q-1)*3 + p7X_D];
//     tmdv.v = om->tfv[7*(q-1) + p7O_MD];
//     tddv.v = om->tfv[7*Q + (q-1)];
//   } else {
//     mpv.v  = esl_sse_rightshiftz_float(ox->dpf[i][(Q-1)*3 + p7X_M]);
//     dpv.v  = esl_sse_rightshiftz_float(ox->dpf[i][(Q-1)*3 + p7X_D]);
//     tmdv.v = esl_sse_rightshiftz_float(om->tfv[7*(Q-1) + p7O_MD]);
//     tddv.v = esl_sse_rightshiftz_float(om->tfv[8*Q-1]);
//   }
//   path[0] = ((tmdv.p[r] == 0.0) ? -eslINFINITY : mpv.p[r]);
//   path[1] = ((tddv.p[r] == 0.0) ? -eslINFINITY : dpv.p[r]);
//   return  ((path[0] >= path[1]) ? p7T_M : p7T_D);
fn select_d(t: &OaTrans, ox: &OaMx, i: usize, k: usize) -> u8 {
    let q = ox.q;
    let qi = (k - 1) % q;
    let r = (k - 1) / q;

    let (mpv, dpv, tmdv, tddv) = if qi > 0 {
        (
            ox.get(i, qi - 1, PX_M),
            ox.get(i, qi - 1, PX_D),
            t.vv(7 * (qi - 1) + O_MD),
            t.vv(7 * q + (qi - 1)),
        )
    } else {
        (
            rightshiftz(ox.get(i, q - 1, PX_M)),
            rightshiftz(ox.get(i, q - 1, PX_D)),
            rightshiftz(t.vv(7 * (q - 1) + O_MD)),
            rightshiftz(t.vv(8 * q - 1)),
        )
    };

    let p0 = if tmdv[r] == 0.0 { NEG } else { mpv[r] };
    let p1 = if tddv[r] == 0.0 { NEG } else { dpv[r] };
    if p0 >= p1 {
        ST_M
    } else {
        ST_D
    }
}

/// I(i,k) is reached from M(i-1,k) or I(i-1,k). `select_i()` — impl_sse/optacc.c:347-360.
// C, impl_sse/optacc.c:348-360:
//   __m128 *tp   = om->tfv + 7*q + p7O_MI;
//   mpv.v = ox->dpf[i-1][q*3 + p7X_M]; tv.v = *tp;  path[0] = ((tv.p[r] == 0.0) ? -eslINFINITY : mpv.p[r]);  tp++;
//   ipv.v = ox->dpf[i-1][q*3 + p7X_I]; tv.v = *tp;  path[1] = ((tv.p[r] == 0.0) ? -eslINFINITY : ipv.p[r]);
//   return  ((path[0] >= path[1]) ? p7T_M : p7T_I);
fn select_i(t: &OaTrans, ox: &OaMx, i: usize, k: usize) -> u8 {
    let q = ox.q;
    let qi = (k - 1) % q;
    let r = (k - 1) / q;
    let tp = 7 * qi + O_MI;

    let mpv = ox.get(i - 1, qi, PX_M);
    let p0 = if t.v(tp)[r] == 0.0 { NEG } else { mpv[r] };
    let ipv = ox.get(i - 1, qi, PX_I);
    let p1 = if t.v(tp + 1)[r] == 0.0 { NEG } else { ipv[r] };
    if p0 >= p1 {
        ST_M
    } else {
        ST_I
    }
}

/// N(i) must come from N(i-1) for i>0; else it comes from S.
/// `select_n()` — impl_sse/optacc.c:363-367.
#[inline]
fn select_n(i: usize) -> u8 {
    if i == 0 {
        ST_S
    } else {
        ST_N
    }
}

/// C(i) is reached from E(i) or C(i-1). `select_c()` — impl_sse/optacc.c:370-377.
fn select_c(xf: &XFactors, pp: &Omx, ox: &OaMx, i: usize) -> u8 {
    let p0 = if xf.c_loop == 0.0 {
        NEG
    } else {
        ox.xc[i - 1] + pp.xc[i]
    };
    let p1 = if xf.e_move == 0.0 { NEG } else { ox.xe[i] };
    if p0 > p1 {
        ST_C
    } else {
        ST_E
    }
}

/// J(i) is reached from E(i) or J(i-1). `select_j()` — impl_sse/optacc.c:380-388.
fn select_j(xf: &XFactors, pp: &Omx, ox: &OaMx, i: usize) -> u8 {
    let p0 = if xf.j_loop == 0.0 {
        NEG
    } else {
        ox.xj[i - 1] + pp.xj[i]
    };
    let p1 = if xf.e_loop == 0.0 { NEG } else { ox.xe[i] };
    if p0 > p1 {
        ST_J
    } else {
        ST_E
    }
}

/// E(i) is reached from any M(i,k=1..M) or D(i,k=2..M); assumes all M_k->E, D_k->E are 1.0.
/// `select_e()` — impl_sse/optacc.c:392-412.
// C, impl_sse/optacc.c:392-412:
//   /* precedence rules in case of ties here are a little tricky: M beats D: note the >= max!  */
//   for (q = 0; q < Q; q++)
//     {
//       u.v   = *dp; dp++;  for (r = 0; r < 4; r++) if (u.p[r] >= max) { max = u.p[r]; smax = p7T_M; kmax = r*Q + q + 1; }
//       u.v   = *dp; dp+=2; for (r = 0; r < 4; r++) if (u.p[r] > max)  { max = u.p[r]; smax = p7T_D; kmax = r*Q + q + 1; }
//     }
//   *ret_k = kmax;
//   return smax;
fn select_e(ox: &OaMx, i: usize) -> (u8, usize) {
    let q = ox.q;
    let mut max = NEG;
    // C leaves smax/kmax uninitialized until the first comparison succeeds; with
    // max starting at -eslINFINITY the very first `>= max` always fires, so these
    // initial values are never observed.
    let mut smax = ST_M;
    let mut kmax = 0usize;

    // Precedence rules in case of ties here are a little tricky: M beats D — note the >= max!
    for qi in 0..q {
        let u = ox.get(i, qi, PX_M);
        for (r, &val) in u.iter().enumerate() {
            if val >= max {
                max = val;
                smax = ST_M;
                kmax = r * q + qi + 1;
            }
        }
        let u = ox.get(i, qi, PX_D);
        for (r, &val) in u.iter().enumerate() {
            if val > max {
                max = val;
                smax = ST_D;
                kmax = r * q + qi + 1;
            }
        }
    }
    (smax, kmax)
}

/// B(i) is reached from N(i) or J(i). `select_b()` — impl_sse/optacc.c:416-424.
fn select_b(xf: &XFactors, ox: &OaMx, i: usize) -> u8 {
    let p0 = if xf.n_move == 0.0 { NEG } else { ox.xn[i] };
    let p1 = if xf.j_move == 0.0 { NEG } else { ox.xj[i] };
    if p0 > p1 {
        ST_N
    } else {
        ST_J
    }
}

/// Posterior probability annotation for one traceback position.
/// `get_postprob()` — impl_sse/optacc.c:265-281.
///
/// C's `switch` deliberately falls through for N/C/J: each case returns only when
/// `sprv == scur`, otherwise control reaches `default` and yields 0.0.
// C, impl_sse/optacc.c:265-281:
//   static inline float
//   get_postprob(const P7_OMX *pp, int scur, int sprv, int k, int i)
//   {
//     int     Q     = p7O_NQF(pp->M);
//     int     q     = (k-1) % Q;
//     int     r     = (k-1) / Q;
//     union { __m128 v; float p[4]; } u;
//
//     switch (scur) {
//     case p7T_M: u.v = MMO(pp->dpf[i], q); return u.p[r];
//     case p7T_I: u.v = IMO(pp->dpf[i], q); return u.p[r];
//     case p7T_N: if (sprv == scur) return pp->xmx[i*p7X_NXCELLS+p7X_N];
//     case p7T_C: if (sprv == scur) return pp->xmx[i*p7X_NXCELLS+p7X_C];
//     case p7T_J: if (sprv == scur) return pp->xmx[i*p7X_NXCELLS+p7X_J];
//     default:    return 0.0;
//     }
//   }
fn get_postprob(pp: &Omx, scur: u8, sprv: u8, k: usize, i: usize) -> f32 {
    let m = pp.m;
    match scur {
        ST_M => {
            if k >= 1 && k <= m {
                pp.mmx[i][k]
            } else {
                0.0
            }
        }
        ST_I => {
            if k >= 1 && k <= m {
                pp.imx[i][k]
            } else {
                0.0
            }
        }
        ST_N if sprv == scur => pp.xn[i],
        ST_C if sprv == scur => pp.xc[i],
        ST_J if sprv == scur => pp.xj[i],
        _ => 0.0,
    }
}

/// One traceback position: `(state, k, i, pp)`, matching C's parallel
/// `tr->st[] / tr->k[] / tr->i[] / tr->pp[]` arrays.
pub(crate) type TracePos = (u8, i32, i32, f32);

/// `p7_OATrace()` — impl_sse/optacc.c:216-262.
///
/// Returns the traceback in S..T order. Sequence coordinates are relative to the
/// envelope, i.e. 1..`ox.l`. The `pp` field carries `get_postprob()`, which
/// `p7_alidisplay_Create()` turns into the alignment's `PP` line and which
/// `p7_trace_GetExpectedAccuracy()` sums back to the OA score.
pub(crate) fn oa_trace(t: &OaTrans, xf: &XFactors, pp: &Omx, ox: &OaMx) -> Vec<TracePos> {
    let mut i = ox.l;
    let mut k = 0usize;
    let mut tr: Vec<TracePos> = Vec::new();

    append(&mut tr, ST_T, k as i32, i as i32, 0.0);
    append(&mut tr, ST_C, k as i32, i as i32, 0.0);

    // C, impl_sse/optacc.c:235-260:
    //   if ((status = p7_trace_AppendWithPP(tr, p7T_T, k, i, 0.0)) != eslOK) return status;
    //   if ((status = p7_trace_AppendWithPP(tr, p7T_C, k, i, 0.0)) != eslOK) return status;
    //   s0 = tr->st[tr->N-1];
    //   while (s0 != p7T_S)
    //     {
    //       switch (s0) {
    //       case p7T_M: s1 = select_m(om,     ox, i, k);  k--; i--; break;
    //       case p7T_D: s1 = select_d(om,     ox, i, k);  k--;      break;
    //       case p7T_I: s1 = select_i(om,     ox, i, k);       i--; break;
    //       case p7T_N: s1 = select_n(i);                           break;
    //       case p7T_C: s1 = select_c(om, pp, ox, i);               break;
    //       case p7T_J: s1 = select_j(om, pp, ox, i);               break;
    //       case p7T_E: s1 = select_e(om,     ox, i, &k);           break;
    //       case p7T_B: s1 = select_b(om,     ox, i);               break;
    //       }
    //       postprob = get_postprob(pp, s1, s0, k, i);
    //       if ((status = p7_trace_AppendWithPP(tr, s1, k, i, postprob)) != eslOK) return status;
    //       if ( (s1 == p7T_N || s1 == p7T_J || s1 == p7T_C) && s1 == s0) i--;
    //       s0 = s1;
    //     }
    //   tr->M = om->M;  tr->L = ox->L;
    //   return p7_trace_Reverse(tr);
    let mut s0 = ST_C;
    while s0 != ST_S {
        let s1 = match s0 {
            ST_M => {
                let s = select_m(t, ox, i, k);
                k -= 1;
                i -= 1;
                s
            }
            ST_D => {
                let s = select_d(t, ox, i, k);
                k -= 1;
                s
            }
            ST_I => {
                let s = select_i(t, ox, i, k);
                i -= 1;
                s
            }
            ST_N => select_n(i),
            ST_C => select_c(xf, pp, ox, i),
            ST_J => select_j(xf, pp, ox, i),
            ST_E => {
                let (s, kk) = select_e(ox, i);
                k = kk;
                s
            }
            ST_B => select_b(xf, ox, i),
            _ => unreachable!("bogus state in OA traceback"),
        };

        let postprob = get_postprob(pp, s1, s0, k, i);
        append(&mut tr, s1, k as i32, i as i32, postprob);

        if (s1 == ST_N || s1 == ST_J || s1 == ST_C) && s1 == s0 {
            i -= 1;
        }
        s0 = s1;
    }

    trace_reverse(&mut tr);
    tr
}

/// `p7_trace_AppendWithPP()` — p7_trace.c:1016-1056.
///
/// Which of `k`, `i` and `pp` a position actually stores depends on the state:
/// emit-on-transition states (N/C/J) keep the residue only when the *previously
/// appended* position is the same state, D keeps `k` but no residue, and the
/// states outside the main model keep nothing.
// C, p7_trace.c:1020-1052:
//   switch (st) {
//     /* Emit-on-transition states: */
//   case p7T_N:
//   case p7T_C:
//   case p7T_J:
//     if (tr->st[tr->N-1] == st)
//       { tr->i[tr->N]  = i;  tr->pp[tr->N] = pp; }
//     else
//       { tr->i[tr->N]  = 0;  tr->pp[tr->N] = 0.0; }
//     tr->k[tr->N] = 0;
//     break;
//     /* Nonemitting states, outside main model: */
//   case p7T_X: case p7T_S: case p7T_B: case p7T_E:
//   case p7T_T: tr->i[tr->N] = 0; tr->pp[tr->N] = 0.0; tr->k[tr->N] = 0; break;
//     /* Nonemitting, but in main model (k valid) */
//   case p7T_D: tr->i[tr->N] = 0; tr->pp[tr->N] = 0.0; tr->k[tr->N] = k; break;
//     /* Emitting states, with valid k position in model: */
//   case p7T_M:
//   case p7T_I: tr->i[tr->N] = i; tr->pp[tr->N] = pp;  tr->k[tr->N] = k; break;
//   }
fn append(tr: &mut Vec<TracePos>, st: u8, k: i32, i: i32, pp: f32) {
    let pos = match st {
        ST_N | ST_C | ST_J => {
            if tr.last().map(|p| p.0) == Some(st) {
                (st, 0, i, pp)
            } else {
                (st, 0, 0, 0.0)
            }
        }
        ST_S | ST_B | ST_E | ST_T => (st, 0, 0, 0.0),
        ST_D => (st, k, 0, 0.0),
        ST_M | ST_I => (st, k, i, pp),
        _ => unreachable!("no such state; can't append"),
    };
    tr.push(pos);
}

/// `p7_trace_Reverse()` — p7_trace.c:1109-1151.
///
/// "For emit-on-transition states N,C,J, traces always obey the C-,Cx,Cx,Cx
/// convention even when they were constructed backwards; so we make them
/// Cx,Cx,Cx,C- by pulling residues backwards by one, just before reversing them."
/// The residue *and* its posterior probability move together.
// C, p7_trace.c:1121-1148:
//   for (z = 0; z < tr->N; z++)
//     {
//       if ( (tr->st[z] == p7T_N && tr->st[z+1] == p7T_N) ||
//            (tr->st[z] == p7T_C && tr->st[z+1] == p7T_C) ||
//            (tr->st[z] == p7T_J && tr->st[z+1] == p7T_J))
//         {
//           if (tr->i[z] == 0 && tr->i[z+1] > 0)
//             {
//               tr->i[z]   = tr->i[z+1];
//               tr->i[z+1] = 0;
//               if (tr->pp != NULL) { tr->pp[z] = tr->pp[z+1]; tr->pp[z+1] = 0.0; }
//             }
//         }
//     }
//   /* Reverse the trace in place. */
//   for (z = 0; z < tr->N/2; z++) { ...swap st/k/i/pp... }
fn trace_reverse(tr: &mut [TracePos]) {
    for z in 0..tr.len().saturating_sub(1) {
        let (s, s_next) = (tr[z].0, tr[z + 1].0);
        let emit_on_transition =
            (s == ST_N && s_next == ST_N) || (s == ST_C && s_next == ST_C) || (s == ST_J && s_next == ST_J);
        if emit_on_transition && tr[z].2 == 0 && tr[z + 1].2 > 0 {
            tr[z].2 = tr[z + 1].2;
            tr[z + 1].2 = 0;
            tr[z].3 = tr[z + 1].3;
            tr[z + 1].3 = 0.0;
        }
    }
    tr.reverse();
}

/// `p7_trace_GetExpectedAccuracy()` — p7_trace.c:1085-1093: the sum of the
/// per-position posterior probabilities, which must come back to the OA score.
///
/// Currently exercised only by the `oa_trace_pp_sums_to_oasc` invariant test; C also
/// uses it in `hmmalign`/`tracealign`, which is not ported yet.
#[allow(dead_code)]
pub(crate) fn expected_accuracy(tr: &[TracePos]) -> f32 {
    let mut acc = 0.0f32;
    for p in tr {
        acc += p.3;
    }
    acc
}
