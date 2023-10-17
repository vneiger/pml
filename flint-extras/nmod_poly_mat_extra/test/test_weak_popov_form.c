#include <flint/flint.h>
#include <flint/nmod_poly.h>
#include <stdlib.h>

#include "nmod_poly_mat_forms.h"
#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_io.h"
#include "sagemath_extra.h"

#include "testing_collection.h"
#include <flint/nmod_poly_mat.h>
#include <flint/ulong_extras.h>
#include <flint/profiler.h>

//#define NOTRANS

// verify weak Popov form
int verify_weak_popov_form(const nmod_poly_mat_t wpf, const slong * shift, const nmod_poly_mat_t tsf, slong rk, const nmod_poly_mat_t mat, flint_rand_t state)
{
    if (! tsf)
    {
        printf("~~~ upper, rowwise ~~~ transformation not provided, skipping verification\n");
        return 1;
    }

    // testing correctness, from fastest to slowest:
    // 0. dimensions
    if (mat->r != wpf->r || mat->c != wpf->c || mat->r != tsf->r || mat->r != tsf->c)
    {
            printf("~~~ upper, rowwise ~~~ INCORRECT: dimension mismatch\n");
            return 0;
    }

    // empty matrices are fine (checked by is_weak_popov, but done here otherwise
    // tsf does not pass test on tsf_deg >= 0)
    if (mat->r == 0 || mat->c == 0)
        return 1;

    // 1. unimodular transformation
    // upper bound on degree of determinant
    slong tsf_deg = nmod_poly_mat_degree(tsf);
    if (tsf_deg < 0)
    {
        printf("~~~ upper, rowwise ~~~ INCORRECT: transformation is zero\n");
        return 0;
    }
    ulong degdet = tsf->r * (ulong)tsf_deg;
    if (mat->modulus > 100*degdet)
    {
        if (!nmod_poly_mat_is_unimodular_randomized(tsf, state))
        {
            printf("~~~ upper, rowwise ~~~ INCORRECT: transformation is not unimodular (randomized test)\n");
            return 0;
        }
    }
    else
    {
        if (!nmod_poly_mat_is_unimodular(tsf))
        {
            printf("~~~ upper, rowwise ~~~ INCORRECT: transformation is not unimodular\n");
            return 0;
        }
    }

    // 2. tsf * mat == wpf
    nmod_poly_mat_t tmp;
    nmod_poly_mat_init(tmp, mat->r, mat->c, mat->modulus);
    nmod_poly_mat_mul(tmp, tsf, mat);
    if (! nmod_poly_mat_equal(tmp, wpf))
    {
        printf("~~~ upper, rowwise ~~~ INCORRECT: tsf * mat != wpf\n");
        nmod_poly_mat_clear(tmp);
        return 0;
    }

    // 3. wpf is in weak Popov form
    // first prune zero rows since this is a requirement of is_weak_popov
    // (note that nonzero rows are not necessarily the ones at rk ... rdim-1)
    nmod_poly_mat_clear(tmp);
    nmod_poly_mat_init(tmp, rk, wpf->c, wpf->modulus);
    slong current_row = 0;
    for (slong i = 0; i < wpf->r; i++)
    {
        slong jj = 0;
        while (jj < mat->c && nmod_poly_is_zero(nmod_poly_mat_entry(wpf, i, jj)))
            jj++;
        if (jj < mat->c) // row i is nonzero, jj is first nonzero
        {
            for (slong j = jj; j < wpf->c; j++)
                nmod_poly_set(nmod_poly_mat_entry(tmp, current_row, j), nmod_poly_mat_entry(wpf, i, j));
            current_row++;
        }
    }
    if (! nmod_poly_mat_is_weak_popov_rowwise(tmp, shift))
    {
        printf("~~~ upper, rowwise ~~~ INCORRECT: wpf not in shifted weak Popov form\n");
        nmod_poly_mat_clear(tmp);
        return 0;
    }

    nmod_poly_mat_clear(tmp);
    return 1;
}

// test one given input for weak Popov form
int core_test_weak_popov_form(const nmod_poly_mat_t mat, const slong * shift, int time, flint_rand_t state)
{
    // init copy of mat
    nmod_poly_mat_t wpf;
    nmod_poly_mat_init(wpf, mat->r, mat->c, mat->modulus);

    // init unimodular transformation tsf such that wpf = tsf * mat
    nmod_poly_mat_t tsf;
    nmod_poly_mat_init(tsf, mat->r, mat->r, mat->modulus);

    { // Mulders and Storjohann's algorithm
        nmod_poly_mat_set(wpf, mat);
        nmod_poly_mat_one(tsf);
        timeit_t timer;
        timeit_start(timer);
#ifndef NOTRANS
        slong rk = nmod_poly_mat_weak_popov_mulders_storjohann_upper_rowwise(wpf, shift, tsf);
#else
        slong rk = nmod_poly_mat_weak_popov_mulders_storjohann_upper_rowwise(wpf, NULL);
#endif /* ifndef NOTRANS */
        timeit_stop(timer);
        if (time)
            flint_printf("-- time (Mulders-Storjohann): %wd ms\n", timer->wall);
        timeit_start(timer);
#ifndef NOTRANS
        if (! verify_weak_popov_form(wpf, shift, tsf, rk, mat, state))
#else
        if (! verify_weak_popov_form(wpf, shift, NULL, rk, mat, state))
#endif /* ifndef NOTRANS */
        {
            printf("weak Popov form -- Mulders-Storjohann -- failure.\n");
            return 0;
        }
        timeit_stop(timer);
        if (time)
            flint_printf("-- time (verif): %wd ms\n", timer->wall);
    }

    nmod_poly_mat_clear(wpf);
    nmod_poly_mat_clear(tsf);

    return 1;
}

//** Test against the whole testing collection */
int collection_test_weak_popov_form(slong iter, flint_rand_t state)
{
    // input matrix
    nmod_poly_mat_t mat;
    // TODO add shifts

    long total_nb_tests =
            iter // number of iterations
            * 5 // number of mats (currently 5)
            * _test_collection_nb_primes
            * _test_collection_nb_minidims
            * _test_collection_nb_minidims
            * _test_collection_nb_minidegs;

    printf("Launching testing collection (%ld cases)\n", total_nb_tests);

    for (slong i_primes = 0; i_primes < _test_collection_nb_primes; i_primes++)
        for (slong i_rdims = 0; i_rdims < _test_collection_nb_minidims; i_rdims++)
            for (slong i_cdims = 0; i_cdims < _test_collection_nb_minidims; i_cdims++)
                for (slong i_degs = 0; i_degs < _test_collection_nb_minidegs; i_degs++)
                    for (slong it = 0; it < iter; it++)
                    {
                        const long prime = _test_collection_primes[i_primes];
                        const long rdim = _test_collection_minirdims[i_rdims];
                        const long cdim = _test_collection_minicdims[i_cdims];
                        const long len = _test_collection_minidegs[i_degs];
                        printf("prime %ld, rdim %ld, cdim %ld, length %ld.\n", prime, rdim, cdim, len);

                        nmod_poly_mat_init(mat, rdim, cdim, prime);

                        _test_collection_mat_zero(mat);
                        if (! core_test_weak_popov_form(mat, NULL, 0, state))
                        { printf("failed %s -- %s,\n...exiting\n", "uniform", "zero"); return 0; }

                        _test_collection_mat_uniform(mat, len-1, state);
                        if (! core_test_weak_popov_form(mat, NULL, 0, state))
                        { printf("failed %s -- %s,\n...exiting\n", "uniform", "uniform"); return 0; }

                        _test_collection_mat_test(mat, len-1, state);
                        if (! core_test_weak_popov_form(mat, NULL, 0, state))
                        { printf("failed %s -- %s,\n...exiting\n", "uniform", "test"); return 0; }

                        _test_collection_mat_sparse(mat, len-1, state);
                        if (! core_test_weak_popov_form(mat, NULL, 0, state))
                        { printf("failed %s -- %s,\n...exiting\n", "uniform", "sparse"); return 0; }

                        _test_collection_mat_rkdef(mat, len-1, state);
                        if (! core_test_weak_popov_form(mat, NULL, 0, state))
                        { printf("failed %s -- %s,\n...exiting\n", "uniform", "rkdef"); return 0; }

                        nmod_poly_mat_clear(mat);
                    }

    printf("--> Successful\n");
    return 1;
}

int main(int argc, char ** argv)
{
    printf("Usage: %s OR %s [nbits] [rdim] [cdim] [order]\n--\n", argv[0], argv[0]);

    // disable line buffering
    setbuf(stdout, NULL);

    if (argc != 1 && argc != 5)
        return 1;

    flint_rand_t state;
    flint_randinit(state);
    srand(time(NULL));
    flint_randseed(state, rand(), rand());

    int res = 0;

    if (argc == 1)
    {
        printf("launching test collection...\n");
        res = collection_test_weak_popov_form(10, state);
    }
    else if (argc == 5)
    {
        slong nbits = atoi(argv[1]);
        slong rdim = atoi(argv[2]);
        slong cdim = atoi(argv[3]);
        slong order = atoi(argv[4]);

        slong prime = n_randprime(state, nbits, 1);
        printf("Launching test with\n\tprime = %ld,\n\trdim = %ld,\n\tcdim = %ld,\
               \n\tlen = %ld...\n",prime,rdim,cdim,order);

        nmod_poly_mat_t mat;
        nmod_poly_mat_init(mat, rdim, cdim, prime);
        nmod_poly_mat_rand(mat, state, order);

        res = core_test_weak_popov_form(mat, NULL, 1, state);

        nmod_poly_mat_t ref;
        nmod_poly_mat_init(ref, rdim, cdim, prime);
        nmod_poly_t den;
        nmod_poly_init(den, prime);
        timeit_t timer;
        timeit_start(timer);
        nmod_poly_mat_rref(ref, den, mat);
        timeit_stop(timer);
        nmod_poly_mat_clear(ref);
        nmod_poly_clear(den);
        flint_printf("-- time FFLU: %wd ms (for comparison)\n", timer->wall);

        nmod_poly_mat_clear(mat);
    }

    flint_randclear(state);

    if (res == 0)
    {
        printf("FAILURE\n");
        return 1;
    }
    else
    {
        printf("SUCCESS\n");
        return 0;
    }
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
