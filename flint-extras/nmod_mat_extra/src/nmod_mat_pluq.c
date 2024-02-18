#include <flint/nmod_vec.h> // for _nmod_vec_set

#include "nmod_mat_extra.h"





void _nmod_mat_l_u_from_compactlu_zeroin(nmod_mat_t L,
                                        nmod_mat_t U,
                                        const nmod_mat_t LU,
                                        slong rk)
{
    // LU->r == 0 -> L is 0 x 0 and U is 0 x LU->c  -> nothing to do
    if (LU->r == 0)
        return;

    // L: fill entries row by row, from lower triangular part of LU
    nmod_mat_entry(L, 0, 0) = UWORD(1);  // note L->r == LU->r > 0
    for (slong i = 1; i < rk; i++)
    {
        _nmod_vec_set(L->rows[i], LU->rows[i], i);
        nmod_mat_entry(L, i, i) = UWORD(1);
    }
    for (slong i = rk; i < LU->r; i++)
    {
        _nmod_vec_set(L->rows[i], LU->rows[i], rk);
        nmod_mat_entry(L, i, i) = UWORD(1);
    }

    // U: if U aliases LU, just zero out entries of L
    if (U == LU)
    {
        for (slong i = 1; i < rk; i++)
            _nmod_vec_zero(U->rows[i], i);
        for (slong i = rk; i < LU->r; i++)
            _nmod_vec_zero(U->rows[i], rk);
    }
    // U: else, U == 0, copy entries from upper triangular part of LU
    else
        for (slong i = 0; i < rk; i++)
            _nmod_vec_set(U->rows[i]+i, LU->rows[i]+i, U->c - i);
}

void nmod_mat_l_u_from_compactlu(nmod_mat_t L,
                                 nmod_mat_t U,
                                 const nmod_mat_t LU,
                                 slong rk)
{
    // L: set to zero all entries above diagonal or beyond rk-th column
    for (slong i = 1; i < rk; i++)
        _nmod_vec_zero(L->rows[i]+i+1, L->c - i-1);
    for (slong i = rk; i < L->r; i++)
        _nmod_vec_zero(L->rows[i]+rk, L->c - rk);

    // U: set to zero all entries below diagonal or beyond rk-th row
    if (U != LU)
    {
        for (slong i = 1; i < rk; i++)
            _nmod_vec_zero(U->rows[i], i);
        for (slong i = rk; i < U->r; i++)
            _nmod_vec_zero(U->rows[i], U->c);
    }

    // call "zeroin" variant
    _nmod_mat_l_u_from_compactlu_zeroin(L, U, LU, rk);
}

slong nmod_mat_pluq(slong * P, nmod_mat_t A, slong * Q)
{
    slong nullity = 0;
    slong rank = 0;

    while (rank + nullity < A->r)
    {
        // seek pivot: on row rank, first nonzero entry with column index >= rank
        slong pivot = rank;
        while (pivot < A->c && nmod_mat_entry(A, rank, pivot) == 0)
            pivot++;

        // if no pivot was found: increase rank defect + perform row rotation
        // (current row becomes last row, others rows in between are moved up)
        if (pivot == A->c)
        {
            // TODO row rotation
            nullity++;
        }
        else
        {
            // TODO column rotation
            
            rank++;
        }
    }
    return rank;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
