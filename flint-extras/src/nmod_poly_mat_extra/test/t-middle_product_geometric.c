/*
    Copyright (C) 2025 Vincent Neiger, Éric Schost

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/nmod_poly_mat.h>
#include <flint/test_helpers.h>

#include "nmod_poly_mat_utils.h" // for rand
#include "nmod_poly_mat_multiply.h"

/*--------------------------------------------------------------*/
/* middle product using different implementations               */
/*--------------------------------------------------------------*/
int test_mat_middle_product_geometric(ulong bits, ulong m, ulong n, ulong p, ulong deg, flint_rand_t state)
{
    ulong prime = n_randprime(state, bits, 1);

    nmod_poly_mat_t A, B, C1, C2;

    nmod_poly_mat_init(A, m, n, prime);
    nmod_poly_mat_init(B, n, p, prime);
    nmod_poly_mat_init(C1, m, p, prime);
    nmod_poly_mat_init(C2, m, p, prime);

    nmod_poly_mat_rand(A, state, deg);
    nmod_poly_mat_rand(B, state, 2 * deg - 1);
    nmod_poly_mat_rand(C1, state, deg);
    nmod_poly_mat_rand(C2, state, deg);

    nmod_poly_mat_middle_product_naive(C1, A, B, deg-1, deg-1);
    nmod_poly_mat_middle_product_geometric(C2, A, B, deg-1, deg-1);

    int res = nmod_poly_mat_equal(C1, C2);
    
    nmod_poly_mat_clear(C1);
    nmod_poly_mat_clear(C2);
    nmod_poly_mat_clear(B);
    nmod_poly_mat_clear(A);
   
    return res;
}

TEST_FUNCTION_START(nmod_poly_mat_middle_product_geometric, state)
{
    int i, result;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        /* FIXME field must be "large enough" : add proper test/fallback */
        ulong bits = 16 + n_randint(state, 49);
        ulong m = 1 + n_randint(state, 50);
        ulong n = 1 + n_randint(state, 50);
        ulong p = 1 + n_randint(state, 50);
        ulong deg = 0 + n_randint(state, 50);

        /* flint_printf("m = %wu, n = %wu, p = %wu\n" */
        /*              "deg = %wu, n_bits = %wu\n", */
        /*              m, n, p, deg, bits); */

        result = test_mat_middle_product_geometric(bits, m, n, p, deg, state);

        if (!result)
            TEST_FUNCTION_FAIL(
                    "m = %wu, n = %wu, p = %wu\n"
                    "deg = %wu, n_bits = %wu\n",
                    m, n, p, deg, bits);
    }

    TEST_FUNCTION_END(state);
}
