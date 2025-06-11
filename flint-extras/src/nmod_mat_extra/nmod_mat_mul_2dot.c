/*
    Copyright (C) 2025 Vincent Neiger, Éric Schost

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/nmod_mat.h>
#include <flint/nmod_vec.h>

#include "nmod_vec_extra.h"

#if FLINT_BITS == 64
/** matrix multiplication using AVX2 instructions for moduli less than 2^30 */
void nmod_mat_mul_2dot(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
{
    // C aliases A or B: use a temporary
    if (C == A || C == B)
    {
        nmod_mat_t T;
        nmod_mat_init(T, A->r, B->c, A->mod.n);
        nmod_mat_mul_2dot(T, A, B);
        nmod_mat_swap_entrywise(C, T);
        nmod_mat_clear(T);
        return;
    }

    ulong power_two;
    NMOD_RED(power_two, UWORD(1)<<DOT_SPLIT_BITS, A->mod);

    // transpose of B
    nmod_mat_t BT;
    nmod_mat_init(BT, B->c, B->r, B->mod.n);
    nmod_mat_transpose(BT, B);

    // let's go
    ulong res[2];
    slong i = 0;
    for (; i+1 < A->r; i += 2)
        for (slong j = 0; j < BT->r; j++)
        {
#if __FLINT_VERSION < 3 || (__FLINT_VERSION == 3 && __FLINT_VERSION_MINOR < 3)
            _nmod_vec_2dot2_split(res, A->rows[i], A->rows[i+1], BT->rows[j],
                                  A->c, A->mod, power_two);
            C->rows[i][j] = res[0];
            C->rows[i+1][j] = res[1];
#else
            _nmod_vec_2dot2_split(res,
                                  nmod_mat_entry_ptr(A, i, 0),
                                  nmod_mat_entry_ptr(A, i+1, 0),
                                  nmod_mat_entry_ptr(BT, j, 0),
                                  A->c, A->mod, power_two);
            nmod_mat_entry(C, i, j) = res[0];
            nmod_mat_entry(C, i+1, j) = res[1];
#endif
        }

    for (; i < A->r; i++)
        for (slong j = 0; j < BT->r; j++)
#if __FLINT_VERSION < 3 || (__FLINT_VERSION == 3 && __FLINT_VERSION_MINOR < 3)
            C->rows[i][j] = _nmod_vec_dot2_split(A->rows[i], BT->rows[j], A->c, A->mod, power_two);
#else
            nmod_mat_entry(C, i, j) = _nmod_vec_dot2_split(nmod_mat_entry_ptr(A, i, 0),
                                                           nmod_mat_entry_ptr(BT, j, 0),
                                                           A->c, A->mod, power_two);
#endif

    nmod_mat_clear(BT);
}

void nmod_mat_mul_nmod_vec_2dot(nn_ptr v, const nmod_mat_t A, nn_srcptr u, ulong len)
{
    ulong power_two;
    NMOD_RED(power_two, UWORD(1)<<DOT_SPLIT_BITS, A->mod);

    // let's go
    slong i = 0;
    for (; i+1 < A->r; i += 2)
    {
#if __FLINT_VERSION < 3 || (__FLINT_VERSION == 3 && __FLINT_VERSION_MINOR < 3)
        _nmod_vec_2dot2_split(v+i, A->rows[i], A->rows[i+1], u, len, A->mod, power_two);
#else
        _nmod_vec_2dot2_split(v+i,
                              nmod_mat_entry_ptr(A, i, 0),
                              nmod_mat_entry_ptr(A, i+1, 0),
                              u,
                              len, A->mod, power_two);
#endif
    }

    for (; i < A->r; i++)
#if __FLINT_VERSION < 3 || (__FLINT_VERSION == 3 && __FLINT_VERSION_MINOR < 3)
        v[i] = _nmod_vec_dot2_split(A->rows[i], u, len, A->mod, power_two);
#else
        v[i] = _nmod_vec_dot2_split(nmod_mat_entry_ptr(A, i, 0), u,
                                    len, A->mod, power_two);
#endif
}

#endif  /* FLINT_BITS == 64 */
