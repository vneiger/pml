/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/nmod_mat.h>
#include <flint/nmod_poly_mat.h>
#include <flint/test_helpers.h>

#include "nmod_mat_poly.h"
#include "nmod_poly_mat_arith.h"
#include "nmod_poly_mat_utils.h"

/* pmat <- the constant polynomial matrix cmat */
static void
_mul_nmod_mat_constant(nmod_poly_mat_t pmat, const nmod_mat_t cmat)
{
    nmod_mat_poly_t matp;
    nmod_mat_poly_init_set_from_nmod_mat(matp, cmat);
    nmod_poly_mat_set_from_mat_poly(pmat, matp);
    nmod_mat_poly_clear(matp);
}

/* Compares the products by a constant matrix, on the left and on the right,
 * with nmod_poly_mat_mul; also with the output aliasing the polynomial
 * input. */
TEST_FUNCTION_START(nmod_poly_mat_mul_nmod_mat, state)
{
    int i, result;

    for (i = 0; i < 300 * flint_test_multiplier(); i++)
    {
        const ulong prime = n_randprime(state, 2 + n_randint(state, 63), 1);
        const slong m = n_randint(state, 10);
        const slong n = n_randint(state, 10);
        const slong k = n_randint(state, 10);
        const slong len = n_randint(state, 20);

        nmod_mat_t A, B;
        nmod_poly_mat_t P, Ap, Bp, C, ref;

        nmod_poly_mat_init(P, m, n, prime);
        nmod_poly_mat_randtest(P, state, len);

        /* left: A * P, with A constant k x m */
        nmod_mat_init(A, k, m, prime);
        nmod_mat_randtest(A, state);
        nmod_poly_mat_init(Ap, k, m, prime);
        _mul_nmod_mat_constant(Ap, A);
        nmod_poly_mat_init(ref, k, n, prime);
        nmod_poly_mat_mul(ref, Ap, P);

        nmod_poly_mat_init(C, k, n, prime);
        nmod_poly_mat_randtest(C, state, len + 2);
        nmod_poly_mat_mul_left_nmod_mat(C, A, P);
        result = nmod_poly_mat_equal(C, ref);

        if (result && k == m)
        {
            nmod_poly_mat_set(C, P);
            nmod_poly_mat_mul_left_nmod_mat(C, A, C);
            result = nmod_poly_mat_equal(C, ref);
        }

        if (!result)
            TEST_FUNCTION_FAIL("left, prime = %wu, k = %wd, m = %wd, n = %wd, len = %wd\n",
                               prime, k, m, n, len);

        nmod_poly_mat_clear(C);
        nmod_poly_mat_clear(ref);

        /* right: P * B, with B constant n x k */
        nmod_mat_init(B, n, k, prime);
        nmod_mat_randtest(B, state);
        nmod_poly_mat_init(Bp, n, k, prime);
        _mul_nmod_mat_constant(Bp, B);
        nmod_poly_mat_init(ref, m, k, prime);
        nmod_poly_mat_mul(ref, P, Bp);

        nmod_poly_mat_init(C, m, k, prime);
        nmod_poly_mat_randtest(C, state, len + 2);
        nmod_poly_mat_mul_right_nmod_mat(C, P, B);
        result = nmod_poly_mat_equal(C, ref);

        if (result && k == n)
        {
            nmod_poly_mat_set(C, P);
            nmod_poly_mat_mul_right_nmod_mat(C, C, B);
            result = nmod_poly_mat_equal(C, ref);
        }

        if (!result)
            TEST_FUNCTION_FAIL("right, prime = %wu, m = %wd, n = %wd, k = %wd, len = %wd\n",
                               prime, m, n, k, len);

        nmod_poly_mat_clear(C);
        nmod_poly_mat_clear(ref);
        nmod_poly_mat_clear(Bp);
        nmod_poly_mat_clear(Ap);
        nmod_mat_clear(B);
        nmod_mat_clear(A);
        nmod_poly_mat_clear(P);
    }

    TEST_FUNCTION_END(state);
}
