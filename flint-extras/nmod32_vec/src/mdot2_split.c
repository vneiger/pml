#include <flint/nmod_vec.h>

#include "nmod32_vec.h"

void _nmod32_vec_mdot2_split(n32_ptr mv, n32_srcptr mat, n32_srcptr vec,
                             slong nrows, slong len, slong stride, nmod_t mod)
{
    ulong pow2_precomp;
    NMOD_RED(pow2_precomp, (UWORD(1) << DOT_SPLIT_BITS), mod);

    slong i = 0;
    for ( ; i+1 < nrows; i+=2)
        _nmod32_vec_dot2_split(mv+i, mv+i+1,
                               vec, mat + i*stride, mat + (i+1)*stride,
                               len, mod, pow2_precomp);

    if (i == nrows - 1)
        mv[i] = _nmod32_vec_dot_split(vec, mat + i*stride, len, mod, pow2_precomp);
}

void _nmod32_vec_mdot2_split_avx2(n32_ptr mv, n32_srcptr mat, n32_srcptr vec,
                                  slong nrows, slong len, slong stride, nmod_t mod)
{
    ulong pow2_precomp;
    NMOD_RED(pow2_precomp, (UWORD(1) << DOT_SPLIT_BITS), mod);

    slong i = 0;
    for ( ; i+1 < nrows; i+=2)
        _nmod32_vec_dot2_split_avx2(mv+i, mv+i+1,
                                    vec, mat + i*stride, mat + (i+1)*stride,
                                    len, mod, pow2_precomp);

    if (i == nrows - 1)
        mv[i] = _nmod32_vec_dot_split_avx2(vec, mat + i*stride, len, mod, pow2_precomp);
}

void _nmod32_vec_mdot2_split_avx512(n32_ptr mv, n32_srcptr mat, n32_srcptr vec,
                                    slong nrows, slong len, slong stride, nmod_t mod)
{
    ulong pow2_precomp;
    NMOD_RED(pow2_precomp, (UWORD(1) << DOT_SPLIT_BITS), mod);

    slong i = 0;
    for ( ; i+1 < nrows; i+=2)
        _nmod32_vec_dot2_split(mv+i, mv+i+1,
                               vec, mat + i*stride, mat + (i+1)*stride,
                               len, mod, pow2_precomp);

    if (i == nrows - 1)
        mv[i] = _nmod32_vec_dot_split(vec, mat + i*stride, len, mod, pow2_precomp);
}

