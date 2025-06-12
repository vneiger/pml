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
/* multiplies matrices using different implementations          */
/*--------------------------------------------------------------*/
int one_test_nmod_poly_mat_mul_waksman(ulong bits, ulong m, ulong n, ulong p, ulong deg, flint_rand_t state)
{
    ulong prime = n_randprime(state, bits, 1);

    nmod_poly_mat_t A, B, C1, C2;

    flint_rand_init(state);

    nmod_poly_mat_init(A, m, n, prime);
    nmod_poly_mat_init(B, n, p, prime);
    nmod_poly_mat_init(C1, m, p, prime);
    nmod_poly_mat_init(C2, m, p, prime);

    nmod_poly_mat_randtest(A, state, deg);
    nmod_poly_mat_randtest(B, state, deg);
    nmod_poly_mat_randtest(C1, state, deg);
    nmod_poly_mat_randtest(C2, state, deg);

    nmod_poly_mat_mul(C1, A, B);
    nmod_poly_mat_mul_waksman(C2, A, B);

    int result = nmod_poly_mat_equal(C1, C2);
    
    nmod_poly_mat_clear(C1);
    nmod_poly_mat_clear(C2);
    nmod_poly_mat_clear(B);
    nmod_poly_mat_clear(A);

    return result;
}

TEST_FUNCTION_START(nmod_poly_mat_mul_waksman, state)
{
    int i, result;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        ulong bits = 2 + n_randint(state, 63);
        ulong m = 1 + n_randint(state, 50);
        ulong n = 1 + n_randint(state, 50);
        ulong p = 1 + n_randint(state, 50);
        ulong deg = 0 + n_randint(state, 50);

        result = one_test_nmod_poly_mat_mul_waksman(bits, m, n, p, deg, state);

        if (!result)
            TEST_FUNCTION_FAIL(
                    "m = %wu, n = %wu, p = %wu\n"
                    "n_bits = %wu\n",
                    m, n, p);
    }

    TEST_FUNCTION_END(state);
}
