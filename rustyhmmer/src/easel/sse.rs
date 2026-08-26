//! Easel's vectorized transcendentals.
//!
//! HMMER does **not** use libm for the probability-space profile it hands to
//! Forward/Backward. `fb_conversion()` (`impl_sse/p7_oprofile.c:929`) exponentiates
//! every striped emission and transition score with `esl_sse_expf()`, a Cephes/Pommier
//! polynomial approximation:
//!
//! ```text
//!   om->rfv[x][q] = esl_sse_expf(tmp.v);          /* match emissions   */
//!   om->tfv[j++]  = esl_sse_expf(tmp.v);          /* transitions       */
//! ```
//!
//! It is *not* `expf()`: it differs from a correctly-rounded exponential by a couple
//! of ULP, in a value-dependent direction. Since it is applied to every one of the
//! `Kp*M` emissions and `8*M` transitions that the DP then multiplies together, using
//! libm here shifts the Forward score by far more than the DP's own summation order
//! does. It is part of the model, not a numerical detail.
//!
//! Only the specials go through libm (`om->xf[...] = expf(gm->xsc[...])`,
//! `p7_oprofile.c:986-994`), so those stay on `f32::exp`.

/// `r = exp(x)`, bit-identical to Easel's `esl_sse_expf()`.
///
/// C, easel/esl_sse.c:182-236. The routine is lane-independent, so a scalar
/// transcription is exact — every operation below is the single-precision operation
/// the corresponding `_mm_*_ps` intrinsic performs on one lane.
///
/// ```text
///   __m128
///   esl_sse_expf(__m128 x)
///   {
///     static float cephes_p[6] = { 1.9875691500E-4f, 1.3981999507E-3f, 8.3334519073E-3f,
///                                  4.1665795894E-2f, 1.6666665459E-1f, 5.0000001201E-1f };
///     static float cephes_c[2] = { 0.693359375f,    -2.12194440e-4f };
///     static float maxlogf     =  88.3762626647949f;
///     static float minlogf     = -88.3762626647949f;
///     ...
///     maxmask = _mm_cmpgt_ps(x, _mm_set1_ps(maxlogf));
///     minmask = _mm_cmple_ps(x, _mm_set1_ps(minlogf));
///
///     /* range reduction: exp(x) = 2^k e^f = exp(f + k log 2); k = floorf(0.5 + x / log2): */
///     fx = _mm_mul_ps(x,  _mm_set1_ps(eslCONST_LOG2R));
///     fx = _mm_add_ps(fx, _mm_set1_ps(0.5f));
///
///     /* floorf() with SSE:  */
///     k    = _mm_cvttps_epi32(fx);
///     tmp  = _mm_cvtepi32_ps(k);
///     mask = _mm_cmpgt_ps(tmp, fx);
///     mask = _mm_and_ps(mask, _mm_set1_ps(1.0f));
///     fx   = _mm_sub_ps(tmp, mask);
///     k    = _mm_cvttps_epi32(fx);
///
///     /* polynomial approx for e^f for f in range [-0.5, 0.5] */
///     tmp = _mm_mul_ps(fx, _mm_set1_ps(cephes_c[0]));
///     z   = _mm_mul_ps(fx, _mm_set1_ps(cephes_c[1]));
///     x   = _mm_sub_ps(x, tmp);
///     x   = _mm_sub_ps(x, z);
///     z   = _mm_mul_ps(x, x);
///
///     y =               _mm_set1_ps(cephes_p[0]);    y = _mm_mul_ps(y, x);
///     y = _mm_add_ps(y, _mm_set1_ps(cephes_p[1]));   y = _mm_mul_ps(y, x);
///     y = _mm_add_ps(y, _mm_set1_ps(cephes_p[2]));   y = _mm_mul_ps(y, x);
///     y = _mm_add_ps(y, _mm_set1_ps(cephes_p[3]));   y = _mm_mul_ps(y, x);
///     y = _mm_add_ps(y, _mm_set1_ps(cephes_p[4]));   y = _mm_mul_ps(y, x);
///     y = _mm_add_ps(y, _mm_set1_ps(cephes_p[5]));   y = _mm_mul_ps(y, z);
///     y = _mm_add_ps(y, x);
///     y = _mm_add_ps(y, _mm_set1_ps(1.0f));
///
///     /* build 2^k by hand, by creating a IEEE754 float */
///     k  = _mm_add_epi32(k, _mm_set1_epi32(127));
///     k  = _mm_slli_epi32(k, 23);
///     fx = _mm_castsi128_ps(k);
///
///     /* put 2^k e^f together (fx = 2^k,  y = e^f) and we're done */
///     y = _mm_mul_ps(y, fx);
///
///     /* special/range cleanup */
///     y = esl_sse_select_ps(y, _mm_set1_ps(eslINFINITY), maxmask);
///     y = esl_sse_select_ps(y, _mm_set1_ps(0.0f),        minmask);
///     return y;
///   }
/// ```
///
/// Two transcription notes:
///
/// * `_mm_cvttps_epi32` truncates toward zero and returns `0x80000000` for NaN or for
///   magnitudes that do not fit in an `i32`. Rust's `as i32` float cast saturates
///   instead of wrapping, so the SSE behaviour is spelled out. In practice the masks
///   below overwrite every input that could reach that range, but the intermediate
///   `2^k` construction is done from the raw bits either way, so it must match.
/// * The masks are applied *after* the arithmetic, exactly as C does, because a NaN
///   input must fall through to the `maxmask`/`minmask` selects (both compare false
///   for NaN, so `esl_sse_expf(NaN)` returns whatever the polynomial produced).
#[inline]
pub fn expf(x_in: f32) -> f32 {
    const CEPHES_P: [f32; 6] = [
        1.987_569_15E-4,
        1.398_199_95E-3,
        8.333_451_9E-3,
        4.166_579_6E-2,
        1.666_666_5E-1,
        5.000_000_1E-1,
    ];
    const CEPHES_C: [f32; 2] = [0.693_359_375, -2.121_944_4E-4];
    const MAXLOGF: f32 = 88.376_262_664_794_9;
    const MINLOGF: f32 = -88.376_262_664_794_9;
    // eslCONST_LOG2R = 1/log(2), as the f32 Easel spells it (esl_sse.c uses the
    // `eslCONST_LOG2R` macro from easel.h: 1.44269504088896341).
    const LOG2R: f32 = 1.442_695_040_888_963_4;

    let maxmask = x_in > MAXLOGF;
    let minmask = x_in <= MINLOGF;

    // Range reduction.
    let mut fx = x_in * LOG2R;
    fx += 0.5;

    // floorf() the SSE way.
    let mut k = cvtt(fx);
    let tmp0 = k as f32;
    let mask = if tmp0 > fx { 1.0f32 } else { 0.0f32 };
    fx = tmp0 - mask;
    k = cvtt(fx);

    // Polynomial for e^f, f in [-0.5, 0.5].
    let mut x = x_in - fx * CEPHES_C[0];
    x -= fx * CEPHES_C[1];
    let z = x * x;

    let mut y = CEPHES_P[0];
    y *= x;
    y += CEPHES_P[1];
    y *= x;
    y += CEPHES_P[2];
    y *= x;
    y += CEPHES_P[3];
    y *= x;
    y += CEPHES_P[4];
    y *= x;
    y += CEPHES_P[5];
    y *= z;
    y += x;
    y += 1.0;

    // 2^k, assembled bitwise.
    let two_k = f32::from_bits(((k.wrapping_add(127)) << 23) as u32);
    y *= two_k;

    if maxmask {
        f32::INFINITY
    } else if minmask {
        0.0
    } else {
        y
    }
}

/// `_mm_cvttps_epi32` on one lane: truncate toward zero, or the "integer indefinite"
/// value `0x80000000` when the result is not representable (NaN included).
#[inline]
fn cvtt(v: f32) -> i32 {
    if v.is_nan() || v >= 2147483648.0 || v < -2147483648.0 {
        i32::MIN
    } else {
        v as i32
    }
}

#[cfg(test)]
mod tests {
    use super::expf;

    /// Reference values captured from Easel's own `esl_sse_expf()` (built from
    /// hmmer-3.4/easel and printed as raw f32 bits), covering the range the profile
    /// actually exercises: emission log-odds and transition logs.
    #[test]
    fn matches_easel_bit_for_bit() {
        // (input bits, expected output bits)
        const CASES: &[(u32, u32)] = &[
            (0x00000000, 0x3f800000), // 0              -> 1
            (0xbf800000, 0x3ebc5ab2), // -1             -> 0.36787945
            (0x3f800000, 0x402df854), // 1              -> 2.71828175
            (0xc0490fdb, 0x3d310113), // -3.141593      -> 0.0432139151
            (0xbf000000, 0x3f1b4598), // -0.5           -> 0.606530666
            (0x3f000000, 0x3fd3094c), // 0.5            -> 1.64872122
            (0xc1200000, 0x383e6bce), // -10            -> 4.5399931e-05
            (0xc2c80000, 0x00000000), // -100           -> 0 (minmask)
            (0x41200000, 0x46ac14ee), // 10             -> 22026.4648
            (0xc0000000, 0x3e0a9555), // -2             -> 0.135335281
            (0xbfc00000, 0x3e647c3c), // -1.5           -> 0.223130167
            (0xc07fffff, 0x3c960ab0), // -3.99999976    -> 0.0183156431
            (0x3dcccccd, 0x3f8d763e), // 0.1            -> 1.10517097
            (0xbdcccccd, 0x3f67a36d), // -0.1           -> 0.90483743
            (0xc2b00000, 0x00000000), // -88            -> 0 (2^-127 builds IEEE 0)
            (0x42b00000, 0x7ef882b7), // 88             -> 1.65163627e+38
            (0xc2b0f000, 0x00000000), // -88.46875      -> 0 (minmask)
            (0x42b0f000, 0x7f800000), // 88.46875       -> inf (maxmask)
        ];
        for &(xi, yi) in CASES {
            let x = f32::from_bits(xi);
            let got = expf(x);
            assert_eq!(
                got.to_bits(),
                yi,
                "expf({x:e}) = {got:e} ({:08x}), want {:08x}",
                got.to_bits(),
                yi
            );
        }
    }

    /// Out-of-range clamping, per the `maxmask`/`minmask` selects.
    #[test]
    fn clamps_out_of_range() {
        assert_eq!(expf(89.0), f32::INFINITY);
        assert_eq!(expf(-89.0), 0.0);
        // Verified against Easel: 0x7f800000 -> 0x7f800000, 0xff800000 -> 0x00000000.
        assert_eq!(expf(f32::NEG_INFINITY), 0.0);
        assert_eq!(expf(f32::INFINITY), f32::INFINITY);
    }
}
