#include <stdlib.h>
#include <flint/nmod_types.h>
#include <flint/profiler.h>

#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_extra.h"  // det_iter currently in _extra.h

#include "testing_collection.h"

// verify determinant
int verify_determinant(const nmod_poly_t det, const nmod_poly_mat_t mat, flint_rand_t state)
{
    // checking dimensions
    if (mat->r != mat->c)
    {
            printf("~~~ verify determinant ~~~ INCORRECT: dimension mismatch\n");
            return 0;
    }

    // compute determinant via Flint's native method
    nmod_poly_t det_correct;
    nmod_poly_init(det_correct, mat->modulus);
    nmod_poly_mat_det(det_correct, mat);

    if (! nmod_poly_equal(det, det_correct))
    {
        printf("~~~ verify determinant ~~~ INCORRECT: determinant is wrong\n");
        nmod_poly_print_pretty(det, "X"); printf("\n");
        nmod_poly_print_pretty(det_correct, "X"); printf("\n");
        nmod_poly_clear(det_correct);
        return 0;
    }
    else
    {
        nmod_poly_clear(det_correct);
        return 1;
    }
}

// test one given input
int core_test_determinant(const nmod_poly_mat_t mat, int time, flint_rand_t state)
{
    nmod_poly_t det;
    nmod_poly_init(det, mat->modulus);
    nmod_poly_t det_correct;
    nmod_poly_init(det_correct, mat->modulus);
    // init copy of mat
    nmod_poly_mat_t copy_mat;
    nmod_poly_mat_init(copy_mat, mat->r, mat->c, mat->modulus);

    // verification of determinant
    int verif_det;

    { // Mulders and Storjohann's algorithm, row by row variant
        nmod_poly_mat_set(copy_mat, mat);
        timeit_t timer;
        timeit_start(timer);
        nmod_poly_mat_det_iter(det, copy_mat);
        timeit_stop(timer);
        if (time)
            flint_printf("-- time (iter - rowbyrow): %wd ms\n", timer->wall);
        timeit_start(timer);
        verif_det = verify_determinant(det, mat, state);
        timeit_stop(timer);
        if (! verif_det)
            printf("determinant -- iter-rowbyrow -- determinant failure.\n");
        else
            nmod_poly_set(det_correct, det);
        if (time)
            flint_printf("-- time (verif): %wd ms\n", timer->wall);
    }

    { // Mulders and Storjohann's algorithm, row by row variant
        nmod_poly_mat_set(copy_mat, mat);
        timeit_t timer;
        timeit_start(timer);
        nmod_poly_mat_det_iter(det, copy_mat);
        timeit_stop(timer);
        if (time)
            flint_printf("-- time (iter - rowbyrow): %wd ms\n", timer->wall);
        timeit_start(timer);
        if (verif_det)
            verif_det = nmod_poly_equal(det_correct, det);
        else
            verif_det = verify_determinant(det, mat, state);
        if (! verif_det)
            printf("determinant -- iter-rowbyrow-bis -- determinant failure.\n");
        timeit_stop(timer);
        if (time)
            flint_printf("-- time (verif): %wd ms\n", timer->wall);
    }


    nmod_poly_mat_clear(copy_mat);
    nmod_poly_clear(det);

    return verif_det;
}

//** Test against the whole testing collection */
int collection_test_determinant(slong iter, flint_rand_t state)
{
    // input matrix
    nmod_poly_mat_t mat;

    long total_nb_tests =
            iter // number of iterations
            * 5 // number of mats (currently 5)
            * _test_collection_nb_primes
            * _test_collection_nb_smalldims
            * _test_collection_nb_minidegs;

    printf("Launching testing collection (%ld cases)\n", total_nb_tests);

    for (slong i_primes = 0; i_primes < _test_collection_nb_primes; i_primes++)
        for (slong i_rdims = 0; i_rdims < _test_collection_nb_minidims; i_rdims++)
            for (slong i_degs = 0; i_degs < _test_collection_nb_minidegs; i_degs++)
                for (slong it = 0; it < iter; it++)
                {
                    const long prime = _test_collection_primes[i_primes];
                    const long rdim = _test_collection_minirdims[i_rdims];
                    const long cdim = rdim;
                    const long len = _test_collection_minidegs[i_degs];
                    printf("prime %ld, rdim %ld, cdim %ld, length %ld.\n", prime, rdim, cdim, len);

                    nmod_poly_mat_init(mat, rdim, cdim, prime);

                    _test_collection_mat_zero(mat);
                    if (! core_test_determinant(mat, 0, state))
                    { printf("failed %s -- %s,\n...exiting\n", "uniform", "zero"); return 0; }

                    _test_collection_mat_uniform(mat, len-1, state);
                    if (! core_test_determinant(mat, 0, state))
                    { printf("failed %s -- %s,\n...exiting\n", "uniform", "uniform"); return 0; }

                    _test_collection_mat_test(mat, len-1, state);
                    if (! core_test_determinant(mat, 0, state))
                    { printf("failed %s -- %s,\n...exiting\n", "uniform", "test"); return 0; }

                    _test_collection_mat_sparse(mat, len-1, state);
                    if (! core_test_determinant(mat, 0, state))
                    { printf("failed %s -- %s,\n...exiting\n", "uniform", "sparse"); return 0; }

                    _test_collection_mat_rkdef(mat, len-1, state);
                    if (! core_test_determinant(mat, 0, state))
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
    flint_rand_init(state);
    srand(time(NULL));
    flint_randseed(state, rand(), rand());

    int res = 0;

    if (argc == 1)
    {
        printf("launching test collection...\n");
        res = collection_test_determinant(100, state);
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

        res = core_test_determinant(mat, 1, state);



        nmod_poly_mat_clear(mat);
    }

    flint_rand_clear(state);

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
