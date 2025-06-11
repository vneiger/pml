/*
    Copyright (C) 2025 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

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

#if PML_HAVE_AVX2

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

void _nmod32_vec_mdot3_split_avx2(n32_ptr mv, n32_srcptr mat, n32_srcptr vec,
                                  slong nrows, slong len, slong stride, nmod_t mod)
{
    ulong pow2_precomp;
    NMOD_RED(pow2_precomp, (UWORD(1) << DOT_SPLIT_BITS), mod);

    slong i = 0;

    for ( ; i+2 < nrows; i+=3)
        _nmod32_vec_dot3_split_avx2(mv+i, vec,
                                    mat + i*stride,
                                    mat + (i+1)*stride,
                                    mat + (i+2)*stride,
                                    len, mod, pow2_precomp);

    if (nrows - i == 2)
    {
        _nmod32_vec_dot2_split_avx2(mv+i, mv+i+1,
                                    vec, mat + i*stride, mat + (i+1)*stride,
                                    len, mod, pow2_precomp);
    }
    else if (nrows - i == 1)
        mv[i] = _nmod32_vec_dot_split_avx2(vec, mat + i*stride, len, mod, pow2_precomp);
}

#endif  /* PML_HAVE_AVX2 */

#if PML_HAVE_AVX512

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

void _nmod32_vec_mdot3_split_avx512(n32_ptr mv, n32_srcptr mat, n32_srcptr vec,
                                    slong nrows, slong len, slong stride, nmod_t mod)
{
    ulong pow2_precomp;
    NMOD_RED(pow2_precomp, (UWORD(1) << DOT_SPLIT_BITS), mod);

    slong i = 0;

    for ( ; i+2 < nrows; i+=3)
        _nmod32_vec_dot3_split_avx512(mv+i,
                                      vec,
                                      mat + i*stride,
                                      mat + (i+1)*stride,
                                      mat + (i+2)*stride,
                                      len, mod, pow2_precomp);

    if (nrows - i == 2)
    {
        _nmod32_vec_dot2_split_avx512(mv+i, mv+i+1,
                                      vec, mat + i*stride, mat + (i+1)*stride,
                                      len, mod, pow2_precomp);
    }

    else if (nrows - i == 1)
        mv[i] = _nmod32_vec_dot_split_avx512(vec, mat + i*stride, len, mod, pow2_precomp);
}

#endif  /* PML_HAVE_AVX512 */
