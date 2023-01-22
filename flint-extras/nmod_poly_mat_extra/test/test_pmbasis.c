#include "nmod_poly_mat_approximant.h"
#include "nmod_poly_mat_extra.h"
#include "sagemath_extra.h"

#include <flint/profiler.h>

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

int test_pmbasis(void)
{
    nmod_poly_mat_t mat, appbas1, appbas2;
    slong rdim = 20, cdim = 10, prime = 101 , order = 400, len = 400;
    slong shifts[rdim], shift1[rdim], shift2[rdim];
    flint_rand_t state;

    /** init **/
    nmod_poly_mat_init(mat, rdim, cdim, prime);

    nmod_poly_mat_init(appbas1, rdim, rdim, prime);
    nmod_poly_mat_init(appbas2, rdim, rdim, prime);

    flint_randinit(state);
    srand(time(NULL));
    flint_randseed(state, rand(), rand());

    for (slong i = 0; i < rdim; i++)
        shifts[i] = rand() % 10 - 5;

    nmod_poly_mat_randtest(mat, state, len);

    mbasisIV(appbas1, shift1, mat, order, shifts);

    pmbasis(appbas2, shift2, mat, order, shifts);

    if (!shift_equal(shift1, shift2, rdim))
        printf("pmbasis and mbasis don't return the same shift\n");

    if (!nmod_poly_mat_equal(appbas1, appbas2))
        printf("pmbasis and mbasis don't return the same approximant basis\n");

    nmod_poly_mat_clear(mat);
    nmod_poly_mat_clear(appbas1);
    nmod_poly_mat_clear(appbas2);

    flint_randclear(state);
    return 0;
}

int main(void)
{
    test_pmbasis();
    return EXIT_SUCCESS;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
