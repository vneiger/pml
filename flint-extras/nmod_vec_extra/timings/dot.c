#include <flint/flint.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_mat.h>
#include <flint/nmod.h>

#include "dot.h"

ulong
_nmod_vec_flintdot(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod, int nlimbs)
{
    ulong res;
    slong i;
    FLINT_NMOD_VEC_DOT(res, i, len, vec1[i], vec2[i], mod, nlimbs);
    return res;
}

ulong
_nmod_vec_flintdot_rev(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod, int nlimbs)
{
    ulong res;
    slong i;
    FLINT_NMOD_VEC_DOT(res, i, len, vec1[len-1-i], vec2[len-1-i], mod, nlimbs);
    return res;
}

ulong
_nmod_vec_newdot(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod, int nlimbs, ulong red_pow)
{
    ulong res;
    slong i;
    _NMOD_VEC_DOT(res, i, len, vec1[i], vec2[i], mod, nlimbs, red_pow);
    return res;
}

ulong
_nmod_vec_newdot2(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod, dot_params_t params)
{
    ulong res = UWORD(0);   /* covers _DOT0 */
    slong i;

    if (params.method == _DOT1)
#ifdef HAVE_AVX2   // TODO
    {
        vec4n dp = vec4n_zero();

        slong i = 0;
        for ( ; i+31 < len; i += 32)
        {
            dp = vec4n_add(dp, vec4n_mul(vec4n_load_unaligned(vec1+i+ 0), vec4n_load_unaligned(vec2+i+ 0)));
            dp = vec4n_add(dp, vec4n_mul(vec4n_load_unaligned(vec1+i+ 4), vec4n_load_unaligned(vec2+i+ 4)));
            dp = vec4n_add(dp, vec4n_mul(vec4n_load_unaligned(vec1+i+ 8), vec4n_load_unaligned(vec2+i+ 8)));
            dp = vec4n_add(dp, vec4n_mul(vec4n_load_unaligned(vec1+i+12), vec4n_load_unaligned(vec2+i+12)));
            dp = vec4n_add(dp, vec4n_mul(vec4n_load_unaligned(vec1+i+16), vec4n_load_unaligned(vec2+i+16)));
            dp = vec4n_add(dp, vec4n_mul(vec4n_load_unaligned(vec1+i+20), vec4n_load_unaligned(vec2+i+20)));
            dp = vec4n_add(dp, vec4n_mul(vec4n_load_unaligned(vec1+i+24), vec4n_load_unaligned(vec2+i+24)));
            dp = vec4n_add(dp, vec4n_mul(vec4n_load_unaligned(vec1+i+28), vec4n_load_unaligned(vec2+i+28)));
        }

        for ( ; i + 3 < len; i += 4)
            dp = vec4n_add(dp, vec4n_mul(vec4n_load_unaligned(vec1+i), vec4n_load_unaligned(vec2+i)));

        res = dp[0] + dp[1] + dp[2] + dp[3];

        for (; i < len; i++)
            res += vec1[i] * vec2[i];

        NMOD_RED(res, res, mod);
    }
#else // ifdef HAVE_AVX2
        _NMOD_VEC_DOT1(res, i, len, vec1[i], vec2[i], mod)
#endif // ifdef HAVE_AVX2
    else if (params.method == _DOT2_32_SPLIT)
#ifdef HAVE_AVX2
    {
        const vec4n low_bits = vec4n_set_n(DOT_SPLIT_MASK);
        vec4n dp_lo = vec4n_zero();
        vec4n dp_hi = vec4n_zero();

        slong i = 0;
        for ( ; i+31 < len; i += 32)
        {
            dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(vec1+i+ 0), vec4n_load_unaligned(vec2+i+ 0)));
            dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(vec1+i+ 4), vec4n_load_unaligned(vec2+i+ 4)));
            dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(vec1+i+ 8), vec4n_load_unaligned(vec2+i+ 8)));
            dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(vec1+i+12), vec4n_load_unaligned(vec2+i+12)));
            dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(vec1+i+16), vec4n_load_unaligned(vec2+i+16)));
            dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(vec1+i+20), vec4n_load_unaligned(vec2+i+20)));
            dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(vec1+i+24), vec4n_load_unaligned(vec2+i+24)));
            dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(vec1+i+28), vec4n_load_unaligned(vec2+i+28)));

            dp_hi = vec4n_add(dp_hi, vec4n_bit_shift_right(dp_lo, DOT_SPLIT_BITS));
            dp_lo = vec4n_bit_and(dp_lo, low_bits);
        }

        for ( ; i + 3 < len; i += 4)
            dp_lo = vec4n_add(dp_lo, vec4n_mul(vec4n_load_unaligned(vec1+i), vec4n_load_unaligned(vec2+i)));

        dp_hi = vec4n_add(dp_hi, vec4n_bit_shift_right(dp_lo, DOT_SPLIT_BITS));
        dp_lo = vec4n_bit_and(dp_lo, low_bits);

        ulong hsum_lo = dp_lo[0] + dp_lo[1] + dp_lo[2] + dp_lo[3];
        const ulong hsum_hi = dp_hi[0] + dp_hi[1] + dp_hi[2] + dp_hi[3] + (hsum_lo >> DOT_SPLIT_BITS);
        hsum_lo &= DOT_SPLIT_MASK;

        for (; i < len; i++)
            hsum_lo += vec1[i] * vec2[i];

        NMOD_RED(res, params.pow2_precomp * hsum_hi + hsum_lo, mod);
    }
#else // ifdef HAVE_AVX2
        _NMOD_VEC_DOT2_32_SPLIT(res, i, len, vec1[i], vec2[i], mod,
                params.pow2_precomp)
#endif // ifdef HAVE_AVX2
    else if (params.method == _DOT2_32)
        _NMOD_VEC_DOT2_32(res, i, len, vec1[i], vec2[i], mod)
    else if (params.method == _DOT2)
        _NMOD_VEC_DOT2(res, i, len, vec1[i], vec2[i], mod)
    else if (params.method == _DOT3_UNROLL)
        _NMOD_VEC_DOT3_UNROLL(res, i, len, vec1[i], vec2[i], mod)
    else if (params.method == _DOT3)
        _NMOD_VEC_DOT3(res, i, len, vec1[i], vec2[i], mod)

    return res;
}

ulong
_nmod_vec_newdot_rev(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod, int nlimbs, ulong red_pow)
{
    // TODO
    ulong res;
    slong i;
    _NMOD_VEC_DOT(res, i, len, vec1[len-1-i], vec2[len-1-i], mod, nlimbs, red_pow);
    return res;
}

dot_params_t _nmod_vec_dot_params(ulong len, nmod_t mod)
{
    ulong t2, t1, t0, u1, u0;

    umul_ppmm(t1, t0, mod.n - 1, mod.n - 1);
    umul_ppmm(t2, t1, t1, len);
    umul_ppmm(u1, u0, t0, len);
    add_ssaaaa(t2, t1, t2, t1, UWORD(0), u1);

    dot_params_t params = {_DOT0, UWORD(0)};

    // TODO rethink where to test what; and adjust tests for unrolling

    if (t2 != 0) // three limbs
    {
        if (mod.n <= 6521908912666391107L)
            params.method = _DOT3_UNROLL;
        else
            params.method = _DOT3;
    }

    else if (t1 != 0) // two limbs
    {
        if (mod.n <= UWORD(1515531528) && len <= WORD(134744072)) // TODO update using long
        {
            params.method = _DOT2_32_SPLIT;
            NMOD_RED(params.pow2_precomp, (1L << DOT_SPLIT_BITS), mod);
        }
        else if (mod.n <= (UWORD(1) << (FLINT_BITS / 2)))
            params.method = _DOT2_32;
        else
            params.method = _DOT2;
    }

    // single limb
    else if (u0 != 0)
        params.method = _DOT1;

    // remaining case: u0 == 0
    // <=> mod.n == 1 or len == 0
    // => dot product is zero, _DOT0 already set

    return params;
}

// matmul, C does not alias A or B
void nmod_mat_mul_flint(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
{
    const int nlimbs = _nmod_vec_dot_bound_limbs(A->c, A->mod);
    for (slong i = 0; i < A->r; i++)
        for (slong j = 0; j < B->r; j++)
            C->rows[i][j] = _nmod_vec_flintdot(A->rows[i], B->rows[j], A->c, A->mod, nlimbs);
}

void nmod_mat_mul_newdot1(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
{
    const int nlimbs = _nmod_vec_dot_bound_limbs(A->c, A->mod);
    const ulong red_pow = (1L << DOT_SPLIT_BITS) % A->mod.n;
    for (slong i = 0; i < A->r; i++)
        for (slong j = 0; j < B->r; j++)
            C->rows[i][j] = _nmod_vec_newdot(A->rows[i], B->rows[j], A->c, A->mod, nlimbs, red_pow);
}

#define LOOP(meth)                                                                                         \
do                                                                                                         \
{                                                                                                          \
    if (params.method == (meth))                                                                           \
    {                                                                                                      \
        for (slong i = 0; i < A->r; i++)                                                                   \
            for (slong j = 0; j < B->r; j++)                                                               \
                C->rows[i][j] = _nmod_vec_newdot2(A->rows[i], B->rows[j], A->c, A->mod, params);           \
        return;                                                                                            \
    }                                                                                                      \
} while (0);
void nmod_mat_mul_newdot(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
{
    const dot_params_t params = _nmod_vec_dot_params(A->c, A->mod);
    LOOP(_DOT2);
    LOOP(_DOT1);
    LOOP(_DOT2_32_SPLIT);
    LOOP(_DOT3_UNROLL);
    LOOP(_DOT3);
    LOOP(_DOT2_32);
    LOOP(_DOT0); // could be simplified
}

/*--------------------------------------------------------------*/
/* this is a tiny bit slower ! (not negligible for small sizes) */
/*--------------------------------------------------------------*/

void _nmod_mat_mul_newdot(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B, const dot_params_t params)
{
    for (slong i = 0; i < A->r; i++)
        for (slong j = 0; j < B->r; j++)
            C->rows[i][j] = _nmod_vec_newdot2(A->rows[i], B->rows[j], A->c, A->mod, params);
}

inline void nmod_mat_mul_newdot2(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
{
    const dot_params_t params = _nmod_vec_dot_params(A->c, A->mod);
    _nmod_mat_mul_newdot(C, A, B, params);
}
