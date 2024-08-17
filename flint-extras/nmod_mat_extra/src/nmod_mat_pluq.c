#include <flint/flint.h>
#include <flint/nmod.h> // for nmod_div
#include <flint/nmod_mat.h>
#include <flint/ulong_extras.h> // for n_invmod
#include <flint/nmod_vec.h> // for _nmod_vec_set

#include "nmod_vec_extra.h"
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

slong nmod_mat_pluq(nmod_mat_t A, slong * P, slong * Q)
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
            _nmod_mat_rotate_rows_upward(A, P, rank, A->r - 1);
            nullity++;
        }
        else
        {
            // pivot found, move its column to rank-th column
            _nmod_mat_rotate_columns_rightward(A, Q, rank, pivot);
            // perform elimination
            ulong inv_pivot = n_invmod(nmod_mat_entry(A, rank, rank), A->mod.n);
            for (slong i = rank+1; i < A->r - nullity; i++)
            {
                nmod_mat_entry(A, i, rank) = nmod_mul(nmod_mat_entry(A, i, rank), inv_pivot, A->mod);
                _nmod_vec_scalar_addmul_nmod(A->rows[i]+rank+1,
                                             A->rows[rank]+rank+1,
                                             A->c - rank - 1,
                                             nmod_neg(nmod_mat_entry(A, i, rank), A->mod),
                                             A->mod);
            }
            rank++;
        }
    }
    return rank;
}

slong nmod_mat_pluq_crout(nmod_mat_t A, slong * P, slong * Q)
{
    // nbits of modulus, for delayed reduction
    const ulong nbits = FLINT_BIT_COUNT(A->mod.n);

    // nullity and rank that will be continuously updated
    slong nullity = 0;
    slong rank = 0;
    // will store parts of columns of A, maximum length final rank
    ulong * A_i_piv = flint_malloc(FLINT_MIN(A->r, A->c) * sizeof(ulong));
    // will store results of parts of A multiplied by A_i_piv
    ulong * A_vec = flint_malloc(A->r * sizeof(ulong));

    // we will use some windows of A, for 1st VERSION and OTHERVERSION
    //nmod_mat_t Awin;
    //// other things, only for OTHERVERSION
    //nmod_mat_t Arowwin;
    //nmod_mat_t matvec;
    // tmp vec, only for YETOTHERVERSION
    nn_ptr vecmat = _nmod_vec_init(A->c);
    nn_ptr * Arows_rk = flint_malloc(A->r * sizeof(nn_ptr)); // TODO could be avoided by providing starting column index to dot_product_multi

    while (rank + nullity < A->r)
    {
        // update row rank
        // A[rank, rank:] <-- A[rank, rank:] - A[rank, :rank] * A[:rank, rank:]
        //// 1stVERSION, TESTED:
        //for (slong i = 0; i < rank; i++)
        //    _nmod_vec_scalar_addmul_nmod(A->rows[rank] + rank,
        //                                 A->rows[i] + rank,
        //                                 A->c - rank,
        //                                 nmod_neg(nmod_mat_entry(A, rank, i), A->mod),
        //                                 A->mod);
        //// OTHERVERSION, TESTED: with mul, which actually seems faster than nmod_mat_nmod_vec_mul...
        //nmod_mat_window_init(Arowwin, A, rank, 0, rank+1, rank);
        //nmod_mat_window_init(Awin, A, 0, rank, rank, A->c);
        //nmod_mat_init(matvec, 1, A->c - rank, A->mod.n);
        //nmod_mat_mul(matvec, Arowwin, Awin);
        //_nmod_vec_sub(A->rows[rank]+rank, A->rows[rank]+rank, matvec->rows[0], A->c - rank, A->mod);
        //nmod_mat_clear(matvec);
        //nmod_mat_window_clear(Awin);
        //nmod_mat_window_clear(Arowwin);
        // YETOTHERVERSION, TESTED: with new nmod_mat_nmod_vec_mul (via dot_product_multi)
        for (slong i = 0; i < rank; i++)
            Arows_rk[i] = A->rows[i]+rank;
        nmod_vec_dot_product_multi(vecmat, A->rows[rank], (nn_srcptr *) Arows_rk, rank, A->c - rank, nbits, nbits, A->mod);
        _nmod_vec_sub(A->rows[rank]+rank, A->rows[rank]+rank, vecmat, A->c - rank, A->mod);

        // seek pivot: on row rank, first nonzero entry with column index >= rank
        slong pivot = rank;
        while (pivot < A->c && nmod_mat_entry(A, rank, pivot) == 0)
            pivot++;

        // if no pivot was found: increase rank defect + perform row rotation
        // (current row becomes last row, others rows in between are moved up)
        if (pivot == A->c)
        {
            // TODO rank-sensitive: delay this, or handle zero rows (beyond A->r - nullity) in a smarter way
            _nmod_mat_rotate_rows_upward(A, P, rank, A->r - 1);
            nullity++;
        }
        else
        {
            // precompute inverse of pivot
            ulong inv_pivot = n_invmod(nmod_mat_entry(A, rank, pivot), A->mod.n);
            // compute matrix-vector product A[rank+1:A->r-nullity, :rank] * A[:rank, pivot]
            // 1stVERSION, TESTED:
            //for (slong i = 0; i < rank; i++)
            //    A_i_piv[i] = nmod_mat_entry(A, i, pivot);
            //nmod_mat_window_init(Awin, A, rank+1, 0, A->r - nullity, rank);
            //nmod_mat_mul_nmod_vec(A_vec, Awin, A_i_piv, rank);
            //nmod_mat_window_clear(Awin);
            // 2ndVERSION, TESTED: with new nmod_vec_dot
            const dot_params_t params = _nmod_vec_dot_params(rank, A->mod);
            for (slong i = 0; i < rank; i++)
                A_i_piv[i] = nmod_mat_entry(A, i, pivot);
            for (slong i = 0; i < A->r - nullity - rank - 1; i++)
                A_vec[i] = _nmod_vec_dot(A->rows[i+rank+1], A_i_piv, rank, A->mod, params);
            // update:  from A[rank+1:A->r-nullity, pivot], subtract A_vec and divide by pivot
            for (slong i = rank+1; i < A->r - nullity; i++)
            {
                nmod_mat_entry(A, i, pivot) = _nmod_sub(nmod_mat_entry(A, i, pivot), A_vec[i-rank-1], A->mod);
                nmod_mat_entry(A, i, pivot) = nmod_mul(nmod_mat_entry(A, i, pivot), inv_pivot, A->mod);
            }
            // move pivot column to rank-th column by rotation
            // TODO rank-sensitive: only do this on rows up to A->r - nullity
            _nmod_mat_rotate_columns_rightward(A, Q, rank, pivot);
            rank++;
        }
    }
    flint_free(A_vec);
    flint_free(A_i_piv);

    // for YETOTHERVERSION
    _nmod_vec_clear(vecmat);
    flint_free(Arows_rk);
    return rank;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
