#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>
#include <flint/perm.h>
#include <flint/ulong_extras.h>
#include "nmod_poly_mat_utils.h"

/*******************
*  tested moduli  *
*******************/

// very small 2, 3, 5; then i*10 bits for i = 1..6
const slong _test_collection_nb_primes = 9;
slong _test_collection_primes[]  = {
    2, 3, 5,    // very small ones
    521, 524309, 536870923,  // 10, 20, 30 bits
    549755813911, 562949953421381, 576460752303423619  // 40, 50, 60 bits
 };

/***********************************
*  tested dimensions and degrees  *
***********************************/

const slong _test_collection_nb_dim = 7;
slong _test_collection_rdims[] = {1, 2, 3, 5, 10, 15, 23};
slong _test_collection_cdims[] = {1, 2, 3, 5, 11, 17, 21};

const slong _test_collection_nb_smalldegs = 13;
slong _test_collection_smalldegs[] = {0, 1, 2, 3, 4, 5, 10, 15, 25, 50, 75, 100, 125};

const slong _test_collection_nb_degs = 16;
slong _test_collection_degs[] = {0, 1, 2, 3, 4, 5, 10, 15, 25, 50, 75, 100, 125, 200, 400, 1000};

/*******************
*  tested shifts  *
*******************/

// uniform [0,...,0]
inline void _test_collection_shift_uniform(slong * shift, slong cdim)
{
	for (long i = 0; i < cdim; ++i)
        shift[i] = 0;
}

// increasing [0,1,2,..,rdim-1]
inline void _test_collection_shift_increasing(slong * shift, slong cdim)
{
	for (long i = 0; i < cdim; ++i)
        shift[i] = i;
}

// decreasing [rdim,..,3,2,1]
inline void _test_collection_shift_decreasing(slong * shift, slong cdim)
{
	for (long i = 0; i < cdim; ++i)
        shift[i] = cdim-i;
}

// random shuffle of [0,1,...,rdim-1]
inline void _test_collection_shift_shuffle(slong * shift, slong cdim, flint_rand_s * state)
{
    _perm_randtest(shift, cdim, state);
}

// Hermite shift [0, cdim*deg, 2*cdim*deg, ..., (cdim-1)*cdim*deg]
inline void _test_collection_shift_hermite(slong * shift, slong cdim, slong deg)
{
    for (long i = 0; i < cdim; ++i)
        shift[i] = i*cdim*deg;
}

// reverse Hermite shift [cdim*cdim*deg, ..., 2*cdim*deg, cdim*deg]
inline void _test_collection_shift_rhermite(slong * shift, slong cdim, slong deg)
{
    for (long i = 0; i < cdim; ++i)
        shift[i] = (cdim-i)*cdim*deg;
}

// plateau shift   [0 ... 0  cdim*deg ... cdim*deg]
inline void _test_collection_shift_plateau(slong * shift, slong cdim, slong deg)
{
    for (long i = 0; i < cdim/2; ++i)
        shift[i] = 0;
    for (long i = cdim/2; i < cdim; ++i)
        shift[i] = cdim*deg;
}


// reverse plateau shift   [0 ... 0  cdim*deg ... cdim*deg]
inline void _test_collection_shift_rplateau(slong * shift, slong cdim, slong deg)
{
    for (long i = 0; i < cdim/2; ++i)
        shift[i] = cdim*deg;
    for (long i = cdim/2; i < cdim; ++i)
        shift[i] = 0;
}


/*********************
*  tested matrices  *
*********************/

// TODO : non-uniform degree profiles

// zero matrices
inline void _test_collection_mat_zero(nmod_poly_mat_t mat)
{
    nmod_poly_mat_zero(mat);
}

// uniformly random matrices, uniform degree
inline void _test_collection_mat_uniform(nmod_poly_mat_t mat, slong deg, flint_rand_s * state)
{
    nmod_poly_mat_rand(mat, state, deg+1);
}

// randtest matrices
inline void _test_collection_mat_test(nmod_poly_mat_t mat, slong deg, flint_rand_s * state)
{
    nmod_poly_mat_randtest(mat, state, deg+1);
}

// sparse matrices
inline void _test_collection_mat_sparse(nmod_poly_mat_t mat, slong deg, flint_rand_s * state)
{
    nmod_poly_mat_randtest_sparse(mat, state, deg+1, 0.05);
}

// rank-deficient matrices
void _test_collection_mat_rkdef(nmod_poly_mat_t mat, slong deg, flint_rand_s * state)
{
    // pick random row index and column index
    const long ii = n_randint(state, mat->r);
    const long jj = n_randint(state, mat->c);

    // first fill with random entries
    nmod_poly_mat_randtest(mat, state, deg+1);

    // build random polynomial linear combination of all rows except ii-th
    nmod_poly_mat_t comb_rows;
    nmod_poly_mat_init(comb_rows, 1, mat->r, mat->modulus);
    nmod_poly_zero(nmod_poly_mat_entry(comb_rows, 0, ii));
    nmod_poly_mat_mul(comb_rows, comb_rows, mat);

    // replace ii-th row by the linear combination
    for (long j = 0; j < mat->c; ++j)
        nmod_poly_swap(nmod_poly_mat_entry(mat, ii, j), nmod_poly_mat_entry(comb_rows, 0, j));

    // build random polynomial linear combination of all columns except jj-th
    nmod_poly_mat_t comb_cols;
    nmod_poly_mat_init(comb_cols, mat->c, 1, mat->modulus);
    nmod_poly_zero(nmod_poly_mat_entry(comb_cols, jj, 0));
    nmod_poly_mat_mul(comb_cols, mat, comb_cols);

    // replace jj-th column by the linear combination
    for (long i = 0; i < mat->r; ++i)
        nmod_poly_swap(nmod_poly_mat_entry(mat, i, jj), nmod_poly_mat_entry(comb_rows, i, 0));
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
