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
int test_mat_mulmid_geometric(ulong bits, ulong m, ulong n, ulong p, ulong len, flint_rand_t state)
{
    ulong prime = n_randprime(state, bits, 1);

    nmod_poly_mat_t A, B, C1, C2;

    nmod_poly_mat_init(A, m, n, prime);
    nmod_poly_mat_init(B, n, p, prime);
    nmod_poly_mat_init(C1, m, p, prime);
    nmod_poly_mat_init(C2, m, p, prime);

    nmod_poly_mat_rand(A, state, len);
    nmod_poly_mat_rand(B, state, 2 * len - 1);
    nmod_poly_mat_rand(C1, state, len);
    nmod_poly_mat_rand(C2, state, len);

    nmod_poly_mat_mulmid_naive(C1, A, B, len-1, 2*len-1);
    nmod_poly_mat_mulmid_geometric(C2, A, B, len-1, 2*len-1);

    int res = nmod_poly_mat_equal(C1, C2);
    
    nmod_poly_mat_clear(C1);
    nmod_poly_mat_clear(C2);
    nmod_poly_mat_clear(B);
    nmod_poly_mat_clear(A);
   
    return res;
}

TEST_FUNCTION_START(nmod_poly_mat_mulmid_geometric, state)
{
    int i, result;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        /* FIXME field must be "large enough" : add proper test/fallback */
        ulong bits = 16 + n_randint(state, 49);
        ulong m = 1 + n_randint(state, 50);
        ulong n = 1 + n_randint(state, 50);
        ulong p = 1 + n_randint(state, 50);
        ulong len = 1 + n_randint(state, 50);

        /* flint_printf("m = %wu, n = %wu, p = %wu\n" */
        /*              "len = %wu, n_bits = %wu\n", */
        /*              m, n, p, len, bits); */

        result = test_mat_mulmid_geometric(bits, m, n, p, len, state);

        if (!result)
            TEST_FUNCTION_FAIL(
                    "m = %wu, n = %wu, p = %wu\n"
                    "len = %wu, n_bits = %wu\n",
                    m, n, p, len, bits);
    }

    TEST_FUNCTION_END(state);
}
