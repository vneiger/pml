#include "nmod_poly_mat_approximant.h"
#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_io.h"
#include "sagemath_extra.h"

#include "testing_collection.h"
#include <flint/nmod_poly_mat.h>

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
    nmod_poly_mat_t appbas, appbas2;
    nmod_poly_mat_init(appbas, rdim, rdim, mat->modulus);
    nmod_poly_mat_init(appbas2, rdim, rdim, mat->modulus);

    slong oshift[rdim], oshift2[rdim];

    mbasis(appbas, oshift, mat, order, shift);

    // testing correctness of mbasis
    if (!nmod_poly_mat_is_approximant_basis(appbas, mat, order, shift, ROW_WISE))
    {
        printf("mbasis output is not a minimal approximant basis\n");
        nmod_poly_mat_print(mat, "X");
        nmod_poly_mat_print(appbas,"X");
        slongvec_print_sagemath(shift, rdim);
        slongvec_print_sagemath(oshift, rdim);
        printf("\n");
        return 0;
    }

    // testing other functions return the same
    mbasisII(appbas2, oshift2, mat, order, shift);
    if (!nmod_poly_mat_equal(appbas, appbas2))
    {
        printf("mbasis and mbasisII don't return the same basis\n");
        nmod_poly_mat_print(appbas,"X");
        nmod_poly_mat_print(appbas2,"X");
        return 0;
    }
    if (!shift_equal(oshift, oshift2, rdim))
    {
        printf("mbasis and mbasisII don't return the same shifts\n");
        return 0;
    }

    mbasisIII(appbas2, oshift2, mat, order, shift);
    if (!nmod_poly_mat_equal(appbas, appbas2))
    {
        printf("mbasis and mbasisIII don't return the same basis\n");
        nmod_poly_mat_print(appbas,"X");
        nmod_poly_mat_print(appbas2,"X");
        return 0;
    }
    if (!shift_equal(oshift, oshift2, rdim))
    {
        printf("mbasis and mbasisIII don't return the same shifts\n");
        return 0;
    }

    mbasisIV(appbas2, oshift2, mat, order, shift);
    if (!nmod_poly_mat_equal(appbas, appbas2))
    {
        printf("mbasis and mbasisIV don't return the same basis\n");
        nmod_poly_mat_print(appbas,"X");
        nmod_poly_mat_print(appbas2,"X");
        return 0;
    }
    if (!shift_equal(oshift, oshift2, rdim))
    {
        printf("mbasis and mbasisIV don't return the same shifts\n");
        return 0;
    }


    mbasisV(appbas2, oshift2, mat, order, shift);
    if (!nmod_poly_mat_equal(appbas, appbas2))
    {
        printf("mbasis and mbasisV don't return the same basis\n");
        nmod_poly_mat_print(appbas,"X");
        nmod_poly_mat_print(appbas2,"X");
        return 0;
    }
    if (!shift_equal(oshift, oshift2, rdim))
    {
        printf("mbasis and mbasisV don't return the same shifts\n");
        return 0;
    }

    nmod_poly_mat_clear(appbas);
    nmod_poly_mat_clear(appbas2);

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

    // shift shift random small entries
    slong shift[rdim];

    for (slong i = 0; i < iter; i++)
    {
        //nmod_poly_mat_randtest(mat, state, len);
        nmod_poly_mat_rand(mat, state, len);

        _perm_randtest(shift, rdim, state);
        for (slong i = 0; i < rdim; i++)
            shift[i] = rand() % len - len/2;

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

int main(void)
{
    slong prime = 1125899906842679;
    slong rdim = 3, cdim = 1, order = 4, len = 4;
    slong iter = 100;
    one_test_mbasis(prime, rdim, cdim, order, len, iter);

    return EXIT_SUCCESS;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
