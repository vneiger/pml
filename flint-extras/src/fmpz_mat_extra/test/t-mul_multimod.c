/*
    Copyright (C) 2025 Vincent Neiger, Éric Schost

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/fmpz_mat.h>
#include <flint/test_helpers.h>

#include "fmpz_mat_extra.h"

#if PML_HAVE_AVX2

int test_fmpz_mat_mul(ulong m, ulong n, ulong p, ulong n_bits, flint_rand_t state)
{
    fmpz_mat_t A, B, C1, C2;

    fmpz_mat_init(A, m, n);
    fmpz_mat_init(B, n, p);
    fmpz_mat_init(C1, m, p);
    fmpz_mat_init(C2, m, p);

    fmpz_mat_randbits(A, state, n_bits);
    fmpz_mat_randbits(B, state, n_bits);
    fmpz_mat_randbits(C1, state, n_bits);
    fmpz_mat_randbits(C2, state, n_bits);

    fmpz_mat_mul(C1, A, B);
    fmpz_mat_mul_multimod(C2, A, B);

    int ret = fmpz_mat_equal(C1, C2);

    fmpz_mat_clear(C1);
    fmpz_mat_clear(C2);
    fmpz_mat_clear(B);
    fmpz_mat_clear(A);

    return ret;
}

TEST_FUNCTION_START(fmpz_mat_mul_multimod, state)
{
    int i;
    int FLINT_SET_BUT_UNUSED(result);

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        ulong factor = 1 + n_randint(state, 50);
        ulong n_bits = factor * 200;
        ulong m = n_randint(state, 200 / factor);
        ulong n = n_randint(state, 200 / factor);
        ulong p = n_randint(state, 200 / factor);

        result = test_fmpz_mat_mul(m, n, p, n_bits, state);;

        if (!result)
            TEST_FUNCTION_FAIL(
                    "m = %wu, n = %wu, p = %wu\n"
                    "n_bits = %wu\n",
                    m, n, p, n_bits);
    }

    TEST_FUNCTION_END(state);
}

#else  /* PML_HAVE_AVX2 */

/* just to make sure to have at least one test in main.c */
TEST_FUNCTION_START(fmpz_mat_mul_multimod, state)
{
    int i, result;

    for (i = 0; i < flint_test_multiplier(); i++)
    {
        result = 1;
    }

    TEST_FUNCTION_END(state);
}

#endif  /* PML_HAVE_AVX2 */
