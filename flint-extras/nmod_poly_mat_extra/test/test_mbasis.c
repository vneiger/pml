#include <time.h>
#include <stdlib.h>
#include <flint/nmod_types.h>

#include "nmod_poly_mat_approximant.h"
#include "nmod_poly_mat_io.h"
#include "sagemath_extra.h"

#include "testing_collection.h"

int shift_equal(slong * shift1, slong * shift2, slong length)
{
    for (slong i = 0; i < length; i++)
        if (shift1[i] != shift2[i])
            return 0;
    return 1;
}

// test one given input for mbasis
int core_test_mbasis(nmod_poly_mat_t mat, slong order, slong * shift)
{
    const slong rdim = mat->r;
    nmod_poly_mat_t appbas;

    slong cshift[rdim];
    for (slong i = 0; i < rdim; i++)
        cshift[i] = shift[i];

    nmod_poly_mat_init(appbas, rdim, rdim, mat->modulus);
    nmod_poly_mat_mbasis(appbas, cshift, mat, order);

    // testing correctness of nmod_poly_mat_mbasis
    if (!nmod_poly_mat_is_approximant_basis(appbas, mat, order, shift, ROW_LOWER))
    {
        printf("nmod_poly_mat_mbasis output is not a minimal approximant basis\n");
        printf("Input matrix:\n");
        nmod_poly_mat_print_pretty(mat, "X");
        printf("Output matrix:\n");
        nmod_poly_mat_print(appbas,"X");
        printf("Residual matrix:\n");
        nmod_poly_mat_t res;
        nmod_poly_mat_init(res, appbas->r, mat->c, mat->modulus);
        nmod_poly_mat_mul(res, appbas, mat);
        nmod_poly_mat_print_pretty(res,"X");
        nmod_poly_mat_clear(res);
        printf("Degree matrix:\n");
        nmod_poly_mat_degree_matrix_print_pretty(appbas);
        printf("Input shift:\t");
        slongvec_print_sagemath(shift, rdim);
        printf("Output shift:\t");
        slongvec_print_sagemath(cshift, rdim);
        printf("\n");
        return 0;
    }

    //mbasis(appbas, oshift, mat, order, shift);

    //// testing correctness of mbasis
    //if (!nmod_poly_mat_is_approximant_basis(appbas, mat, order, shift, ROW_LOWER))
    //{
    //    printf("mbasis output is not a minimal approximant basis\n");
    //    printf("Input matrix:\n");
    //    nmod_poly_mat_print(mat, "X");
    //    printf("Output matrix:\n");
    //    nmod_poly_mat_print(appbas,"X");
    //    printf("Residual matrix:\n");
    //    nmod_poly_mat_t res;
    //    nmod_poly_mat_init(res, appbas->r, mat->c, mat->modulus);
    //    nmod_poly_mat_mul(res, appbas, mat);
    //    nmod_poly_mat_print(res,"X");
    //    printf("Degree matrix:\n");
    //    nmod_poly_mat_degree_matrix_print_pretty(appbas);
    //    printf("Input shift:\t");
    //    slongvec_print_sagemath(shift, rdim);
    //    printf("Output shift:\t");
    //    slongvec_print_sagemath(oshift, rdim);
    //    printf("\n");
    //    return 0;
    //}

    //// testing other functions return the same
    //mbasisII(appbas2, oshift2, mat, order, shift);
    //if (!nmod_poly_mat_equal(appbas, appbas2))
    //{
    //    printf("mbasis and mbasisII don't return the same basis\n");
    //    nmod_poly_mat_print(appbas,"X");
    //    nmod_poly_mat_print(appbas2,"X");
    //    return 0;
    //}
    //if (!shift_equal(oshift, oshift2, rdim))
    //{
    //    printf("mbasis and mbasisII don't return the same shifts\n");
    //    return 0;
    //}

    //mbasisIII(appbas2, oshift2, mat, order, shift);
    //if (!nmod_poly_mat_equal(appbas, appbas2))
    //{
    //    printf("mbasis and mbasisIII don't return the same basis\n");
    //    nmod_poly_mat_print(appbas,"X");
    //    nmod_poly_mat_print(appbas2,"X");
    //    return 0;
    //}
    //if (!shift_equal(oshift, oshift2, rdim))
    //{
    //    printf("mbasis and mbasisIII don't return the same shifts\n");
    //    return 0;
    //}

    //mbasisIV(appbas2, oshift2, mat, order, shift);
    //if (!nmod_poly_mat_equal(appbas, appbas2))
    //{
    //    printf("mbasis and mbasisIV don't return the same basis\n");
    //    nmod_poly_mat_print(appbas,"X");
    //    nmod_poly_mat_print(appbas2,"X");
    //    return 0;
    //}
    //if (!shift_equal(oshift, oshift2, rdim))
    //{
    //    printf("mbasis and mbasisIV don't return the same shifts\n");
    //    return 0;
    //}


    //mbasisV(appbas2, oshift2, mat, order, shift);
    //if (!nmod_poly_mat_equal(appbas, appbas2))
    //{
    //    printf("mbasis and mbasisV don't return the same basis\n");
    //    nmod_poly_mat_print(appbas,"X");
    //    nmod_poly_mat_print(appbas2,"X");
    //    return 0;
    //}
    //if (!shift_equal(oshift, oshift2, rdim))
    //{
    //    printf("mbasis and mbasisV don't return the same shifts\n");
    //    return 0;
    //}

    nmod_poly_mat_clear(appbas);

    return 1;
}

/** Test with specified parameters, uniform shift */
int one_test_mbasis(slong prime, slong rdim, slong cdim, slong order, slong len, slong iter)
{
    flint_rand_t state;
    flint_randinit(state);
    srand(time(NULL));
    flint_randseed(state, rand(), rand());

    // random matrix
    nmod_poly_mat_t mat;
    nmod_poly_mat_init(mat, rdim, cdim, prime);

    // shift random uniform entries
    slong shift[rdim];
    _test_collection_shift_uniform(shift, rdim);

    for (slong i = 0; i < iter; i++)
    {
        //nmod_poly_mat_randtest(mat, state, len);
        nmod_poly_mat_rand(mat, state, len);

        //_perm_randtest(shift, rdim, state);
        //for (slong i = 0; i < rdim; i++)
        //    shift[i] = rand() % len - len/2;
        int res = core_test_mbasis(mat, order, shift);
        if (res == 0)
        {
            printf("failed at iteration %ld... exiting\n", i);
            return 0;
        }
    }

    printf("All %ld iterations went fine! exiting\n", iter);
    nmod_poly_mat_clear(mat);
    flint_randclear(state);

    return 1;
}

/** Test against the whole testing collection */
int collection_test_mbasis(slong iter)
{
    flint_rand_t state;
    flint_randinit(state);
    srand(time(NULL));
    flint_randseed(state, rand(), rand());

    // input matrix for approximation
    nmod_poly_mat_t mat;

    // input shift
    slong * shift;

    long total_nb_tests =
            iter // number of iterations
            * 40 // number of mats (currently 5) x number of shifts (currently 8)
            * _test_collection_nb_primes
            * _test_collection_nb_dims
            * _test_collection_nb_dims
            * _test_collection_nb_smalldegs;
    printf("Launching testing collection (%ld cases)\n", total_nb_tests);

    for (slong i_primes = 0; i_primes < _test_collection_nb_primes; i_primes++)
        for (slong i_rdims = 0; i_rdims < _test_collection_nb_dims; i_rdims++)
            for (slong i_cdims = 0; i_cdims < _test_collection_nb_dims; i_cdims++)
                for (slong i_degs = 0; i_degs < _test_collection_nb_smalldegs; i_degs++)
                    for (slong it = 0; it < iter; it++)
                    {
                        const long prime = _test_collection_primes[i_primes];
                        const long rdim = _test_collection_rdims[i_rdims];
                        const long cdim = _test_collection_cdims[i_cdims];
                        const long order = _test_collection_smalldegs[i_degs];
                        //printf("prime %ld, rdim %ld, cdim %ld, order %ld.\n", prime, rdim, cdim, order);

                        shift = (slong *) flint_malloc(rdim * sizeof(slong));
                        nmod_poly_mat_init(mat, rdim, cdim, prime);

                        _test_collection_shift_uniform(shift, rdim);
                        _test_collection_mat_zero(mat);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "uniform", "zero"); return 0; }

                        _test_collection_shift_uniform(shift, rdim);
                        _test_collection_mat_uniform(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "uniform", "uniform"); return 0; }

                        _test_collection_shift_uniform(shift, rdim);
                        _test_collection_mat_test(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "uniform", "test"); return 0; }

                        _test_collection_shift_uniform(shift, rdim);
                        _test_collection_mat_sparse(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "uniform", "sparse"); return 0; }

                        _test_collection_shift_uniform(shift, rdim);
                        _test_collection_mat_rkdef(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "uniform", "rkdef"); return 0; }


                        _test_collection_shift_increasing(shift, rdim);
                        _test_collection_mat_zero(mat);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "increasing", "zero"); return 0; }

                        _test_collection_shift_increasing(shift, rdim);
                        _test_collection_mat_uniform(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "increasing", "uniform"); return 0; }

                        _test_collection_shift_increasing(shift, rdim);
                        _test_collection_mat_test(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "increasing", "test"); return 0; }

                        _test_collection_shift_increasing(shift, rdim);
                        _test_collection_mat_sparse(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "increasing", "sparse"); return 0; }

                        _test_collection_shift_increasing(shift, rdim);
                        _test_collection_mat_rkdef(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "increasing", "rkdef"); return 0; }


                        _test_collection_shift_decreasing(shift, rdim);
                        _test_collection_mat_zero(mat);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "decreasing", "zero"); return 0; }

                        _test_collection_shift_decreasing(shift, rdim);
                        _test_collection_mat_uniform(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "decreasing", "uniform"); return 0; }

                        _test_collection_shift_decreasing(shift, rdim);
                        _test_collection_mat_test(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "decreasing", "test"); return 0; }

                        _test_collection_shift_decreasing(shift, rdim);
                        _test_collection_mat_sparse(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "decreasing", "sparse"); return 0; }

                        _test_collection_shift_decreasing(shift, rdim);
                        _test_collection_mat_rkdef(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "decreasing", "rkdef"); return 0; }


                        _test_collection_shift_shuffle(shift, rdim, state);
                        _test_collection_mat_zero(mat);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "shuffle", "zero"); return 0; }

                        _test_collection_shift_shuffle(shift, rdim, state);
                        _test_collection_mat_uniform(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "shuffle", "uniform"); return 0; }

                        _test_collection_shift_shuffle(shift, rdim, state);
                        _test_collection_mat_test(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "shuffle", "test"); return 0; }

                        _test_collection_shift_shuffle(shift, rdim, state);
                        _test_collection_mat_sparse(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "shuffle", "sparse"); return 0; }

                        _test_collection_shift_shuffle(shift, rdim, state);
                        _test_collection_mat_rkdef(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "shuffle", "rkdef"); return 0; }


                        _test_collection_shift_hermite(shift, rdim, order);
                        _test_collection_mat_zero(mat);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "hermite", "zero"); return 0; }

                        _test_collection_shift_hermite(shift, rdim, order);
                        _test_collection_mat_uniform(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "hermite", "uniform"); return 0; }

                        _test_collection_shift_hermite(shift, rdim, order);
                        _test_collection_mat_test(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "hermite", "test"); return 0; }

                        _test_collection_shift_hermite(shift, rdim, order);
                        _test_collection_mat_sparse(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "hermite", "sparse"); return 0; }

                        _test_collection_shift_hermite(shift, rdim, order);
                        _test_collection_mat_rkdef(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "hermite", "rkdef"); return 0; }


                        _test_collection_shift_rhermite(shift, rdim, order);
                        _test_collection_mat_zero(mat);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "rhermite", "zero"); return 0; }

                        _test_collection_shift_rhermite(shift, rdim, order);
                        _test_collection_mat_uniform(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "rhermite", "uniform"); return 0; }

                        _test_collection_shift_rhermite(shift, rdim, order);
                        _test_collection_mat_test(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "rhermite", "test"); return 0; }

                        _test_collection_shift_rhermite(shift, rdim, order);
                        _test_collection_mat_sparse(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "rhermite", "sparse"); return 0; }

                        _test_collection_shift_rhermite(shift, rdim, order);
                        _test_collection_mat_rkdef(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "rhermite", "rkdef"); return 0; }


                        _test_collection_shift_plateau(shift, rdim, order);
                        _test_collection_mat_zero(mat);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "plateau", "zero"); return 0; }

                        _test_collection_shift_plateau(shift, rdim, order);
                        _test_collection_mat_uniform(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "plateau", "uniform"); return 0; }

                        _test_collection_shift_plateau(shift, rdim, order);
                        _test_collection_mat_test(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "plateau", "test"); return 0; }

                        _test_collection_shift_plateau(shift, rdim, order);
                        _test_collection_mat_sparse(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "plateau", "sparse"); return 0; }

                        _test_collection_shift_plateau(shift, rdim, order);
                        _test_collection_mat_rkdef(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "plateau", "rkdef"); return 0; }


                        _test_collection_shift_rplateau(shift, rdim, order);
                        _test_collection_mat_zero(mat);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "rplateau", "zero"); return 0; }

                        _test_collection_shift_rplateau(shift, rdim, order);
                        _test_collection_mat_uniform(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "rplateau", "uniform"); return 0; }

                        _test_collection_shift_rplateau(shift, rdim, order);
                        _test_collection_mat_test(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "rplateau", "test"); return 0; }

                        _test_collection_shift_rplateau(shift, rdim, order);
                        _test_collection_mat_sparse(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "rplateau", "sparse"); return 0; }

                        _test_collection_shift_rplateau(shift, rdim, order);
                        _test_collection_mat_rkdef(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "rplateau", "rkdef"); return 0; }

                        nmod_poly_mat_clear(mat);
                        flint_free(shift);
                    }

    printf("--> Successful\n");
    return 1;
}

int main(int argc, char ** argv)
{
    printf("Usage: %s OR %s [nbits] [rdim] [cdim] [order]\n--\n", argv[0], argv[0]);

    if (argc != 1 && argc != 5)
        return 1;

    if (argc == 1)
    {
        printf("launching test collection...\n");
        collection_test_mbasis(10);
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
               \n\torder = %ld,\n\tlen = %ld...\n",prime,rdim,cdim,order,order);

        one_test_mbasis(prime, rdim, cdim, order, order, 10000);
    }

    return EXIT_SUCCESS;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
