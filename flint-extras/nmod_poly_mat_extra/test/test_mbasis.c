#include "nmod_poly_mat_approximant.h"
#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_io.h"
#include "sagemath_extra.h"

#include "testing_collection.h"
#include <flint/nmod_poly_mat.h>
#include <flint/ulong_extras.h>

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
    nmod_poly_mat_init(appbas, rdim, rdim, mat->modulus);

    slong cshift[rdim];
    for (slong i = 0; i < rdim; i++)
        cshift[i] = shift[i];

    nmod_poly_mat_init(appbas, rdim, rdim, mat->modulus);
    nmod_poly_mat_mbasis(appbas, cshift, mat, order);

    printf("Degree matrix:\n");
    nmod_poly_mat_degree_matrix_print_pretty(appbas);

    // testing correctness of nmod_poly_mat_mbasis
    if (!nmod_poly_mat_is_approximant_basis(appbas, mat, order, shift, ROW_WISE))
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
    //if (!nmod_poly_mat_is_approximant_basis(appbas, mat, order, shift, ROW_WISE))
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

    // shift random small entries
    slong shift[rdim];
    _test_collection_shift_uniform(shift, rdim);

    for (slong i = 0; i < iter; i++)
    {
        //nmod_poly_mat_randtest(mat, state, len);
        nmod_poly_mat_rand(mat, state, len);

        //_perm_randtest(shift, rdim, state);
        //for (slong i = 0; i < rdim; i++)
        //    shift[i] = rand() % len - len/2;

        if (core_test_mbasis(mat, order, shift) == 0)
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

    // shift random small entries
    slong shift;

    // TODO finish
}

int main(int argc, char ** argv)
{

    if (argc != 1 && argc != 5)
    {
        printf("Usage: %s OR %s [nbits] [rdim] [cdim] [order]\n", argv[0], argv[0]);
        return 1;
    }


    if (argc == 1)
    {
        printf("launching test collection...\n");
        collection_test_mbasis(1000);
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
