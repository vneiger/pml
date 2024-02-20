#include "nmod_vec_extra.h"
#include "nmod_mat_extra.h"
#include <flint/flint.h>
#include <flint/nmod_mat.h>
#include <flint/perm.h>
#include <stdlib.h> // qsort

static int _slong_comparator(const void * a, const void * b)
{
    if ( *(slong *)a  <  *(slong *)b )
        return -1;
    if ( *(slong *)a  == *(slong *)b )
        return 0;
    else // if ( *(slong *)a  >  *(slong *)b )
        return 1;
}

// uniform random
void nmod_mat_rand(nmod_mat_t mat, flint_rand_t state)
{
    _nmod_vec_rand(mat->entries, state, mat->r * mat->c, mat->mod);
}

// random dense with specified rank
void nmod_mat_randrank_dense(nmod_mat_t mat,
                             flint_rand_t state,
                             slong rank)
{
    nmod_mat_randrank(mat, state, rank);
    nmod_mat_randops(mat, state, (mat->r+mat->c)*(mat->r+mat->c));
    // heuristic: number (nrows + ncols)**2 taken through a few experiments to
    // make matrix look rather dense
}

// random lower row echelon form
void nmod_mat_rand_lref(nmod_mat_t mat,
                        flint_rand_t state,
                        slong rank,
                        int unit)
{
    // check given rank is acceptable
    if (rank < 0 || rank > mat->r || rank > mat->c)
    {
        flint_printf("Exception (nmod_mat_rand_lref). Impossible rank.\n");
        flint_abort();
    }

    // there will be rank nonzero rows, we need to fix their respective pivot
    // indices: take a random permutation of length mat->c and keep the first
    // rank entries as pivots
    slong * pivots = (slong *) flint_malloc(rank * sizeof(slong));
    slong * tmp = _perm_init(mat->c);
    _perm_randtest(tmp, mat->c, state);
    for (slong i = 0; i < rank; i++)
        pivots[i] = tmp[i];
    _perm_clear(tmp);
    qsort(pivots, rank, sizeof(slong), _slong_comparator);

    slong i, j;

    // fill nonzero rows
    for (i = 0; i < rank; i++)
    {
        for (j = 0; j < pivots[i]; j++)
            nmod_mat_entry(mat, i, j) = n_randint(state, mat->mod.n);

        nmod_mat_entry(mat, i, pivots[i]) = n_randint(state, mat->mod.n);
        if (unit || nmod_mat_entry(mat, i, pivots[i]) == UWORD(0))
            nmod_mat_entry(mat, i, pivots[i]) = UWORD(1);

        for (j = pivots[i]+1; j < mat->c; j++)
            nmod_mat_entry(mat, i, j) = UWORD(0);
    }

    // make sure remaining rows are zero
    for (i = rank; i < mat->r; i++)
        for (j = 0; j < mat->c; j++)
            nmod_mat_entry(mat, i, j) = UWORD(0);

    flint_free(pivots);
}

// random upper row echelon form
void nmod_mat_rand_uref(nmod_mat_t mat,
                        flint_rand_t state,
                        slong rank,
                        int unit)
{
    // check given rank is acceptable
    if (rank < 0 || rank > mat->r || rank > mat->c)
    {
        flint_printf("Exception (nmod_mat_rand_uref). Impossible rank.\n");
        flint_abort();
    }

    // there will be rank nonzero rows, we need to fix their respective pivot
    // indices: take a random permutation of length mat->c and keep the first
    // rank entries as pivots
    slong * pivots = (slong *) flint_malloc(rank * sizeof(slong));
    slong * tmp = _perm_init(mat->c);
    _perm_randtest(tmp, mat->c, state);
    for (slong i = 0; i < rank; i++)
        pivots[i] = tmp[i];
    _perm_clear(tmp);
    qsort(pivots, rank, sizeof(slong), _slong_comparator);

    slong i, j;

    // fill nonzero rows
    for (i = 0; i < rank; i++)
    {
        for (j = 0; j < pivots[i]; j++)
            nmod_mat_entry(mat, i, j) = UWORD(0);

        nmod_mat_entry(mat, i, pivots[i]) = n_randint(state, mat->mod.n);
        if (unit || nmod_mat_entry(mat, i, pivots[i]) == UWORD(0))
            nmod_mat_entry(mat, i, pivots[i]) = UWORD(1);


        for (j = pivots[i]+1; j < mat->c; j++)
            nmod_mat_entry(mat, i, j) = n_randint(state, mat->mod.n);
    }

    // make sure remaining rows are zero
    for (i = rank; i < mat->r; i++)
        for (j = 0; j < mat->c; j++)
            nmod_mat_entry(mat, i, j) = UWORD(0);

    flint_free(pivots);
}

// random lower reduced row echelon form
void nmod_mat_rand_lrref(nmod_mat_t mat,
                         flint_rand_t state,
                         slong rank)
{
    // check given rank is acceptable
    if (rank < 0 || rank > mat->r || rank > mat->c)
    {
        flint_printf("Exception (nmod_mat_rand_lrref). Impossible rank.\n");
        flint_abort();
    }

    // there will be rank nonzero rows, we need to fix their respective pivot
    // indices: take a random permutation of length mat->c and keep the first
    // rank entries as pivots
    slong * pivots = (slong *) flint_malloc(rank * sizeof(slong));
    slong * tmp = _perm_init(mat->c);
    _perm_randtest(tmp, mat->c, state);
    for (slong i = 0; i < rank; i++)
        pivots[i] = tmp[i];
    _perm_clear(tmp);
    qsort(pivots, rank, sizeof(slong), _slong_comparator);

    slong i, j;

    // fill nonzero rows
    for (i = rank-1; i >= 0; i--)
    {
        for (j = 0; j < pivots[i]; j++)
            nmod_mat_entry(mat, i, j) = n_randint(state, mat->mod.n);

        nmod_mat_entry(mat, i, pivots[i]) = UWORD(1);

        for (j = pivots[i]+1; j < mat->c; j++)
            nmod_mat_entry(mat, i, j) = UWORD(0);

        for (j = i+1; j < rank; j++)
            nmod_mat_entry(mat, j, pivots[i]) = UWORD(0);
    }

    // make sure remaining rows are zero
    for (i = rank; i < mat->r; i++)
        for (j = 0; j < mat->c; j++)
            nmod_mat_entry(mat, i, j) = UWORD(0);

    flint_free(pivots);
}

// random upper reduced row echelon form
void nmod_mat_rand_urref(nmod_mat_t mat,
                         flint_rand_t state,
                         slong rank)
{
    // check given rank is acceptable
    if (rank < 0 || rank > mat->r || rank > mat->c)
    {
        flint_printf("Exception (nmod_mat_rand_urref). Impossible rank.\n");
        flint_abort();
    }

    // there will be rank nonzero rows, we need to fix their respective pivot
    // indices: take a random permutation of length mat->c and keep the first
    // rank entries as pivots
    slong * pivots = (slong *) flint_malloc(rank * sizeof(slong));
    slong * tmp = _perm_init(mat->c);
    _perm_randtest(tmp, mat->c, state);
    for (slong i = 0; i < rank; i++)
        pivots[i] = tmp[i];
    _perm_clear(tmp);
    qsort(pivots, rank, sizeof(slong), _slong_comparator);

    slong i, j;

    // fill nonzero rows
    for (i = 0; i < rank; i++)
    {
        for (j = 0; j < pivots[i]; j++)
            nmod_mat_entry(mat, i, j) = UWORD(0);

        nmod_mat_entry(mat, i, pivots[i]) = UWORD(1);

        for (j = pivots[i]+1; j < mat->c; j++)
            nmod_mat_entry(mat, i, j) = n_randint(state, mat->mod.n);

        for (j = 0; j < i; j++)
            nmod_mat_entry(mat, j, pivots[i]) = UWORD(0);
    }

    // make sure remaining rows are zero
    for (i = rank; i < mat->r; i++)
        for (j = 0; j < mat->c; j++)
            nmod_mat_entry(mat, i, j) = UWORD(0);

    flint_free(pivots);
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
