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

// verify Hermite form

int verify_hermite_form(const nmod_poly_mat_t hnf, const nmod_poly_mat_t tsf, slong rk, const nmod_poly_mat_t mat)
{
    if (! tsf)
    {
        printf("~~~ upper, rowwise ~~~ transformation not provided, skipping verification\n");
        return 1;
    }

    // testing correctness, from fastest to slowest:
    // 0. dimensions
    if (mat->r != hnf->r || mat->c != hnf->c || mat->r != tsf->r || mat->r != tsf->c)
    {
            printf("~~~ upper, rowwise ~~~ INCORRECT: dimension mismatch\n");
            return 0;
    }

    // 1. unimodular transformation
    if (mat->r != mat->c)
    {
        if (! nmod_poly_mat_is_unimodular(tsf))
        {
            printf("~~~ upper, rowwise ~~~ INCORRECT: transformation is not unimodular\n");
            return 0;
        }
    }
    else
    {
        // NOTE the above is slow due to likely large degrees in transformation
        // --> circumvent it in the square case
        nmod_poly_t det1;
        nmod_poly_init(det1, mat->modulus);
        nmod_poly_mat_det(det1, mat);
        nmod_poly_make_monic(det1, det1);

        nmod_poly_t det2;
        nmod_poly_init(det2, mat->modulus);
        nmod_poly_one(det2);
        for (slong i = 0; i < hnf->r; i++)
            nmod_poly_mul(det2, det2, nmod_poly_mat_entry(hnf, i, i));
        nmod_poly_make_monic(det2, det2);
        if (! nmod_poly_equal(det1, det2))
        {
            printf("~~~ upper, rowwise ~~~ INCORRECT: transformation is not unimodular\n");
            nmod_poly_clear(det1);
            nmod_poly_clear(det2);
            return 0;
        }
        nmod_poly_clear(det1);
        nmod_poly_clear(det2);
    }

    // 2. tsf * mat == hnf
    nmod_poly_mat_t tmp;
    nmod_poly_mat_init(tmp, mat->r, mat->c, mat->modulus);
    nmod_poly_mat_mul(tmp, tsf, mat);
    if (! nmod_poly_mat_equal(tmp, hnf))
    {
        printf("~~~ upper, rowwise ~~~ INCORRECT: tsf * mat != hnf\n");
        nmod_poly_mat_clear(tmp);
        return 0;
    }

    // 3. hnf is in Hermite form
    // first prune zero rows since this is a requirement of is_hermite
    // (note that nonzero rows are necessarily the ones at rk ... rdim-1)
    nmod_poly_mat_clear(tmp);
    nmod_poly_mat_init(tmp, rk, hnf->c, hnf->modulus);
    for (slong i = 0; i < rk; i++)
        for (slong j = 0; j < hnf->c; j++)
            nmod_poly_set(nmod_poly_mat_entry(tmp, i, j), nmod_poly_mat_entry(hnf, i, j));
    if (! nmod_poly_mat_is_uhermite_rowwise(tmp))
    {
        printf("~~~ upper, rowwise ~~~ INCORRECT: hnf not in Hermite form\n");
        nmod_poly_mat_clear(tmp);
        return 0;
    }

    nmod_poly_mat_clear(tmp);
    return 1;
}

// test one given input for hermite form
int core_test_hermite_form(nmod_poly_mat_t mat, int time)
{
    // init copy of mat
    nmod_poly_mat_t hnf;
    nmod_poly_mat_init(hnf, mat->r, mat->c, mat->modulus);

    // init unimodular transformation tsf such that hnf = tsf * mat
    nmod_poly_mat_t tsf;
    nmod_poly_mat_init(tsf, mat->r, mat->r, mat->modulus);

    { // Rosser's algorithm
        nmod_poly_mat_set(hnf, mat);
        nmod_poly_mat_one(tsf);
        timeit_t timer;
        timeit_start(timer);
        slong rk = nmod_poly_mat_upper_hermite_form_rowwise_rosser(hnf, tsf);
        //slong rk = nmod_poly_mat_upper_hermite_form_rowwise_rosser(hnf, NULL);
        timeit_stop(timer);
        if (time)
            flint_printf("-- time (Rosser): %wd ms\n", timer->wall);
        timeit_start(timer);
        if (! verify_hermite_form(hnf, tsf, rk, mat))
        //if (! verify_hermite_form(hnf, NULL, rk, mat))
            printf("Hermite form -- Rosser -- failure.\n");
        timeit_stop(timer);
        if (time)
            flint_printf("-- time (verif): %wd ms\n", timer->wall);
    }

    { // Bradley's algorithm
        nmod_poly_mat_set(hnf, mat);
        nmod_poly_mat_one(tsf);
        timeit_t timer;
        timeit_start(timer);
        slong rk = nmod_poly_mat_upper_hermite_form_rowwise_rosser(hnf, tsf);
        //slong rk = nmod_poly_mat_upper_hermite_form_rowwise_bradley(hnf, NULL);
        timeit_stop(timer);
        if (time)
            flint_printf("-- time (Bradley): %wd ms\n", timer->wall);
        timeit_start(timer);
        if (! verify_hermite_form(hnf, tsf, rk, mat))
        //if (! verify_hermite_form(hnf, NULL, rk, mat))
            printf("Hermite form -- Bradley -- failure.\n");
        timeit_stop(timer);
        if (time)
            flint_printf("-- time (verif): %wd ms\n", timer->wall);
    }

    nmod_poly_mat_clear(hnf);
    nmod_poly_mat_clear(tsf);

    return 1;
}

//** Test against the whole testing collection */
int collection_test_hermite_form(slong iter)
{
    flint_rand_t state;
    flint_randinit(state);
    srand(time(NULL));
    flint_randseed(state, rand(), rand());

    // input matrix
    nmod_poly_mat_t mat;

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
                        if (! core_test_hermite_form(mat, 0))
                        { printf("failed %s -- %s,\n...exiting\n", "uniform", "zero"); return 0; }

                        _test_collection_mat_uniform(mat, len-1, state);
                        if (! core_test_hermite_form(mat, 0))
                        { printf("failed %s -- %s,\n...exiting\n", "uniform", "uniform"); return 0; }

                        _test_collection_mat_test(mat, len-1, state);
                        if (! core_test_hermite_form(mat, 0))
                        { printf("failed %s -- %s,\n...exiting\n", "uniform", "test"); return 0; }

                        _test_collection_mat_sparse(mat, len-1, state);
                        if (! core_test_hermite_form(mat, 0))
                        { printf("failed %s -- %s,\n...exiting\n", "uniform", "sparse"); return 0; }

                        _test_collection_mat_rkdef(mat, len-1, state);
                        if (! core_test_hermite_form(mat, 0))
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

    if (argc == 1)
    {
        printf("launching test collection...\n");
        collection_test_hermite_form(10);
    }
    else if (argc == 5)
    {
        slong nbits = atoi(argv[1]);
        slong rdim = atoi(argv[2]);
        slong cdim = atoi(argv[3]);
        slong order = atoi(argv[4]);

        flint_rand_t state;
        flint_randinit(state);
        srand(time(NULL));
        flint_randseed(state, rand(), rand());
            
        slong prime = n_randprime(state, nbits, 1);
        printf("Launching test with\n\tprime = %ld,\n\trdim = %ld,\n\tcdim = %ld,\
               \n\tlen = %ld...\n",prime,rdim,cdim,order);

        nmod_poly_mat_t mat;
        nmod_poly_mat_init(mat, rdim, cdim, prime);
        nmod_poly_mat_rand(mat, state, order);

        core_test_hermite_form(mat, 1);

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

    return EXIT_SUCCESS;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
