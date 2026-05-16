/*
    Copyright (C) 2026 Vincent Neiger, Éric Schost

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/nmod_poly_mat.h>
#include <flint/test_helpers.h>
#include <flint/ulong_extras.h>

#include "nmod_poly_mat_multiply.h"

/* TODO more primes, more variations of degrees, etc. */

/*--------------------------------------------------------------*/
/* middle product using different implementations               */
/*--------------------------------------------------------------*/
int test_mat_mul_vandermonde(ulong bits, ulong m, ulong n, ulong p, slong len, flint_rand_t state)
{
    ulong prime = n_randprime(state, bits, 1);
    ulong prod_len = 2*len - 1;

    /* for vandermonde1, we need 0, ..., prod_len-1 to be distinct points in Z/prime Z */
    if (prime <= prod_len)
        return 0;

    nmod_poly_mat_t A, B, C1, C2;

    nmod_poly_mat_init(A, m, n, prime);
    nmod_poly_mat_init(B, n, p, prime);
    nmod_poly_mat_init(C1, m, p, prime);
    nmod_poly_mat_init(C2, m, p, prime);

    nmod_poly_mat_randtest(A, state, len);
    nmod_poly_mat_randtest(B, state, len);
    nmod_poly_mat_randtest(C1, state, len);
    nmod_poly_mat_randtest(C2, state, len);

    nmod_poly_mat_mul(C1, A, B);

    nmod_poly_mat_mul_vandermonde1(C2, A, B);
    int res1 = nmod_poly_mat_equal(C1, C2);

    int res2 = 1;
    /* for vandermonde2, we need 1**2, 2**2, ..., prod_len**2 to be distinct points in Z/prime Z */
    if (prime > prod_len*prod_len)
    {
        nmod_poly_mat_mul_vandermonde2(C2, A, B);
        res2 = nmod_poly_mat_equal(C1, C2);
    }

    nmod_poly_mat_clear(C1);
    nmod_poly_mat_clear(C2);
    nmod_poly_mat_clear(B);
    nmod_poly_mat_clear(A);

    if (!res1)
        return 1;
    if (!res2)
        return -1;
    return 0;
}

TEST_FUNCTION_START(nmod_poly_mat_mul_vandermonde, state)
{
    int i, result;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        ulong bits = 2 + n_randint(state, 63);
        ulong m = 1 + n_randint(state, 25);
        ulong n = 1 + n_randint(state, 25);
        ulong p = 1 + n_randint(state, 25);
        ulong len = 0 + n_randint(state, 100);

        result = test_mat_mul_vandermonde(bits, m, n, p, len, state);

        if (result)
            TEST_FUNCTION_FAIL(
                    "error code = %d, m = %wu, n = %wu, p = %wu\n"
                    "n_bits = %wu\n",
                    result, m, n, p, bits);
    }

    TEST_FUNCTION_END(state);
}

