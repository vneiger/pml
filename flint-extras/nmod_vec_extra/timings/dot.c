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
    ulong res;
    slong i;
    _NMOD_VEC_DOT_NEW(res, i, len, vec1[i], vec2[i], mod, params);
    return res;
}

ulong
_nmod_vec_newdot_rev(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod, int nlimbs, ulong red_pow)
{
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
