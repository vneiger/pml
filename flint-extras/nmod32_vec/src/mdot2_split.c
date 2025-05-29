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
        _nmod32_vec_dot2_split_avx512(mv+i, mv+i+1,
                                      vec, mat + i*stride, mat + (i+1)*stride,
                                      len, mod, pow2_precomp);

    if (i == nrows - 1)
        mv[i] = _nmod32_vec_dot_split_avx512(vec, mat + i*stride, len, mod, pow2_precomp);
}

void _nmod32_vec_mdot3_split_avx2(n32_ptr mv, n32_srcptr mat, n32_srcptr vec,
                                  slong nrows, slong len, slong stride, nmod_t mod)
{
    ulong pow2_precomp;
    NMOD_RED(pow2_precomp, (UWORD(1) << DOT_SPLIT_BITS), mod);

    slong i = 0;

    for ( ; i+2 < nrows; i+=3)
        _nmod32_vec_dot3_split_avx2(mv+i, mv+i+1, mv+i+2,
                                    vec, mat + i*stride, mat + (i+1)*stride, mat + (i+2)*stride,
                                    len, mod, pow2_precomp);

    for ( ; i < nrows; i++)
        mv[i] = _nmod32_vec_dot_split_avx512(vec, mat + i*stride, len, mod, pow2_precomp);
}

void _nmod32_vec_mdot4_split_avx512(n32_ptr mv, n32_srcptr mat, n32_srcptr vec,
                                    slong nrows, slong len, slong stride, nmod_t mod)
{
    ulong pow2_precomp;
    NMOD_RED(pow2_precomp, (UWORD(1) << DOT_SPLIT_BITS), mod);

    slong i = 0;

    for ( ; i+3 < nrows; i+=4)
        _nmod32_vec_dot4_split_avx512(mv+i, mv+i+1, mv+i+2, mv+i+3, vec,
                                      mat + i*stride,
                                      mat + (i+1)*stride,
                                      mat + (i+2)*stride,
                                      mat + (i+3)*stride,
                                      len, mod, pow2_precomp);

    for ( ; i < nrows; i++)
        mv[i] = _nmod32_vec_dot_split_avx512(vec, mat + i*stride, len, mod, pow2_precomp);
}

