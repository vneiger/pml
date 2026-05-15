/*
    Copyright (C) 2025 Vincent Neiger, Éric Schost

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/flint.h>
#include <flint/ulong_extras.h>
#include <flint/nmod_poly_mat.h>
#include <flint/test_helpers.h>

#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_multiply.h"

int test_mat_mulmid(ulong prime, nmod_poly_mat_t A, nmod_poly_mat_t B, slong nlo, slong nhi)
{
    nmod_poly_mat_t C, D;
    nmod_poly_mat_init(C, A->r, B->c, prime);
    nmod_poly_mat_init(D, A->r, B->c, prime);

    /* most naive way */
    /* flint_printf("preparing test\n"); */
    if (nlo >= nhi)
        nmod_poly_mat_zero(C);
    else
    {
        nmod_poly_mat_mul(C, A, B);
        nmod_poly_mat_shift_right(C, C, nlo);
        nmod_poly_mat_truncate(C, nhi - nlo);
    }

    /* mulmid */
    /* flint_printf("compute\n"); */
    nmod_poly_mat_mulmid(D, A, B, nlo, nhi);

    /* flint_printf("test equal\n"); */
    int res = nmod_poly_mat_equal(C, D);
    
    nmod_poly_mat_clear(C);
    nmod_poly_mat_clear(D);
   
    return res;
}

TEST_FUNCTION_START(nmod_poly_mat_mulmid, state)
{
    int i, result;
    srand(time(NULL));
    flint_rand_set_seed(state, rand(), 12);

    for (i = 0; i < 300 * flint_test_multiplier(); i++)
    {
        /* TODO geometric variant supposed to fail if no suitable geometric
         * progression.. (i.e., if field too small) -> make test more robust */
        ulong bits = 16 + n_randint(state, 49);
        ulong prime = n_randprime(state, bits, 1);

        ulong m = 1 + n_randint(state, 20);
        ulong n = 1 + n_randint(state, 20);
        ulong p = 1 + n_randint(state, 20);

        /* constraints: */
        slong lenA, lenB, nlo, nhi;
        /* if (i < 100 * flint_test_multiplier()) */
        if (0)
        {
            /* tests on random parameters */
            nlo = n_randint(state, 50);
            nhi = n_randint(state, 100);
            lenA = n_randint(state, 50);
            lenB = n_randint(state, 50);
        }
        else if (i < 10000 * flint_test_multiplier())
        {
            /* tests on "classical" parameters for transposed multiplication: */
            /* nlo < nhi, lenA <= nlo+1, lenB <= nhi */
            nlo = n_randint(state, 50);
            nhi = nlo+1 + n_randint(state, 75);
            lenA = n_randint(state, nlo+2);
            lenB = n_randint(state, nhi+1);
        }
        else
        {
            /* same as above, permuted: lenA <= nhi, lenB <= nlo+1 */
            nlo = n_randint(state, 50);
            nhi = n_randint(state, 100);
            lenA = n_randint(state, nhi+1);
            lenB = n_randint(state, nlo+2);
        }

        /* flint_printf("NEW -- prime = %wu, m = %wu, n = %wu, p = %wu\n" */
        /*              "bits = %wu, lenA = %wd, lenB = %wd, nlo = %wd, nhi = %wd\n", */
        /*              prime, m, n, p, bits, lenA, lenB, nlo, nhi); */

        nmod_poly_mat_t A, B;
        nmod_poly_mat_init(A, m, n, prime);
        if (lenA > 0)
            nmod_poly_mat_randtest(A, state, lenA);  /* TODO use randtest */
        nmod_poly_mat_init(B, n, p, prime);
        if (lenB > 0)
            nmod_poly_mat_randtest(B, state, lenB);  /* TODO use randtest */


        /* flint_printf("let's go\n"); */
        result = test_mat_mulmid(prime, A, B, nlo, nhi);

        if (!result)
            TEST_FUNCTION_FAIL(
                    "prime = %wu, m = %wu, n = %wu, p = %wu\n"
                    "bits = %wu, lenA = %wd, lenB = %wd, nlo = %wd, nhi = %wd\n",
                    prime, m, n, p, bits, lenA, lenB, nlo, nhi);

        nmod_poly_mat_clear(A);
        nmod_poly_mat_clear(B);
    }

    TEST_FUNCTION_END(state);
}
