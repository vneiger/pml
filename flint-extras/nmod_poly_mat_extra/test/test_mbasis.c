#include "nmod_poly_mat_approximant.h"
#include "nmod_poly_mat_io.h"

// TODO make random choice for given prime (or for random prime of given size)
#define PRIME_30_BITS 536870923
#define PRIME_60_BITS 576460752303423619

// TODO make more robust and reusable framework test for approximant basis algorithms
// and use it for mbasis1

// TODO general: test for memory leaks

int shift_equal(slong * shift1, slong * shift2, slong length)
{
    for (slong i = 0; i < length; i++)
        if (shift1[i] != shift2[i])
            return 0;
    return 1;
}

/** Verify if all versions of mbasis give the same results */
int test_mbasis(void)
{
    nmod_poly_mat_t mat, appbas1, appbas2, appbas3, appbas4, appbas5;
    slong rdim = 3, cdim = 1, prime = 9001, order = 4, len = 4;
    slong shifts[rdim], shift1[rdim], shift2[rdim];
    slong shift3[rdim], shift4[rdim], shift5[rdim];

    flint_rand_t state;

    /** init */
    nmod_poly_mat_init(mat, rdim, cdim, prime);

    nmod_poly_mat_init(appbas1, rdim, rdim, prime);
    nmod_poly_mat_init(appbas2, rdim, rdim, prime);
    nmod_poly_mat_init(appbas3, rdim, rdim, prime);
    nmod_poly_mat_init(appbas4, rdim, rdim, prime);
    nmod_poly_mat_init(appbas5, rdim, rdim, prime);

    flint_randinit(state);
    srand(time(NULL));
    flint_randseed(state, rand(), rand());

    nmod_poly_mat_randtest(mat, state, len);

    for (slong i = 0; i < 100000; i++)
    {
        nmod_poly_mat_randtest(mat, state, len);

        //_perm_randtest(shifts, rdim, state);
        //for (slong i = 0; i < rdim; i++)
        //    shifts[i] = rand() % len - 10;

        mbasis(appbas1, shift1, mat, order, shifts);

        mbasisII(appbas2, shift2, mat, order, shifts);

        mbasisIII(appbas3, shift3, mat, order, shifts);

        mbasisIV(appbas4, shift4, mat, order, shifts);

        mbasisV(appbas5, shift5, mat, order, shifts);

        // testing correctness of mbasis
        if (!nmod_poly_mat_is_approximant_basis(appbas1, mat, order, shifts, ROW_WISE))
        {
            printf("mbasis output is not a minimal approximant basis\n");
            nmod_poly_mat_print(mat, "X");
            nmod_poly_mat_print(appbas1,"X");
            printf("\n");
            return 0;
        }

        // testing other functions return the same
        if (!nmod_poly_mat_equal(appbas1, appbas2))
        {
            printf("mbasis and mbasisII don't return the same result\n");
            nmod_poly_mat_print(appbas1,"X");
            nmod_poly_mat_print(appbas2,"X");
            return 0;
        }

        if (!nmod_poly_mat_equal(appbas3, appbas1))
        {
            printf("mbasis and mbasisIII don't return the same result\n");
            nmod_poly_mat_print(appbas1,"X");
            nmod_poly_mat_print(appbas3,"X");
            return 0;
        }

        if (!nmod_poly_mat_equal(appbas4, appbas1))
        {
            printf("mbasis and mbasisIV don't return the same result\n");
            nmod_poly_mat_print(appbas1,"X");
            nmod_poly_mat_print(appbas4,"X");
            return 0;
        }

        if (!nmod_poly_mat_equal(appbas5, appbas1))
        {
            printf("mbasis and mbasisV don't return the same result\n");
            nmod_poly_mat_print(appbas1,"X");
            nmod_poly_mat_print(appbas5,"X");
            return 0;
        }

        if (!shift_equal(shift1, shift2, rdim))
        {
            printf("mbasis and mbasisII don't return the same shifts result\n");
            return 0;
        }

        if (!shift_equal(shift1, shift3, rdim))
        {
            printf("mbasis and mbasisIII don't return the same shifts result\n");
            return 0;
        }

        if (!shift_equal(shift1, shift4, rdim))
        {
            printf("mbasis and mbasisIV don't return the same shifts result\n");
            return 0;
        }

        if (!shift_equal(shift1, shift5, rdim))
        {
            printf("mbasis and mbasisV don't return the same shifts result\n");
            return 0;
        }
    }
    /** clear **/
    nmod_poly_mat_clear(mat);
    nmod_poly_mat_clear(appbas1);
    nmod_poly_mat_clear(appbas2);
    nmod_poly_mat_clear(appbas3);
    nmod_poly_mat_clear(appbas4);
    nmod_poly_mat_clear(appbas5);

    flint_randclear(state);
    return 0;
}

int main(void)
{
    test_mbasis();
    return EXIT_SUCCESS;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
