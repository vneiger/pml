/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>
#include <flint/test_helpers.h>
#include <flint/ulong_extras.h>

#include "nmod_poly_mat_linsolve.h"
#include "linsolve_test_utils.h"

/* the solver under test: nmod_poly_mat_linsolve_via_series, or
 * nmod_poly_mat_linsolve_via_high_order_lifting if hol is nonzero */
static int _lsseriestest_solve(nmod_poly_mat_t u, nmod_poly_t den,
                               const nmod_poly_mat_t A, const nmod_poly_mat_t b,
                               linsolve_den_t den_method, int hol)
{
    return hol ? nmod_poly_mat_linsolve_via_high_order_lifting(u, den, A, b, den_method)
               : nmod_poly_mat_linsolve_via_series(u, den, A, b, den_method);
}

/* fixed examples: [[x, 0], [0, 1]] (A(0) singular), whose solution for b =
 * [1, 1] is [1, x] / x; and over GF(2), [[x^2 + x]] whose determinant vanishes
 * at all elements of the field (via_kernel is used), with b = [1] */
static int _lsseriestest_examples(linsolve_den_t den_method, int hol)
{
    int ok;
    nmod_poly_mat_t A, b, u, v;
    nmod_poly_t den;

    nmod_poly_mat_init(A, 2, 2, 101);
    nmod_poly_mat_init(b, 1, 2, 101);
    nmod_poly_mat_init(u, 1, 2, 101);
    nmod_poly_mat_init(v, 1, 2, 101);
    nmod_poly_init(den, 101);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(A, 0, 0), 1, 1);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(A, 1, 1), 0, 1);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(b, 0, 0), 0, 1);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(b, 0, 1), 0, 1);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(v, 0, 0), 0, 1);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(v, 0, 1), 1, 1);
    ok = _lsseriestest_solve(u, den, A, b, den_method, hol)
         && nmod_poly_mat_equal(u, v) && nmod_poly_degree(den) == 1
         && den->coeffs[0] == 0 && den->coeffs[1] == 1;
    nmod_poly_mat_clear(A);
    nmod_poly_mat_clear(b);
    nmod_poly_mat_clear(u);
    nmod_poly_mat_clear(v);
    nmod_poly_clear(den);

    nmod_poly_mat_init(A, 1, 1, 2);
    nmod_poly_mat_init(b, 1, 1, 2);
    nmod_poly_mat_init(u, 1, 1, 2);
    nmod_poly_init(den, 2);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(A, 0, 0), 1, 1);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(A, 0, 0), 2, 1);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(b, 0, 0), 0, 1);
    ok = ok && _lsseriestest_solve(u, den, A, b, den_method, hol)
         && nmod_poly_is_one(nmod_poly_mat_entry(u, 0, 0))
         && nmod_poly_equal(den, nmod_poly_mat_entry(A, 0, 0));
    nmod_poly_mat_clear(A);
    nmod_poly_mat_clear(b);
    nmod_poly_mat_clear(u);
    nmod_poly_clear(den);

    return ok;
}

TEST_FUNCTION_START(nmod_poly_mat_linsolve_via_series, state)
{
    int i, result;

    for (int hol = 0; hol < 2; hol++)
        if (!_lsseriestest_examples(NMOD_POLY_MAT_LINSOLVE_RECONSTRUCTION, hol)
            || !_lsseriestest_examples(NMOD_POLY_MAT_LINSOLVE_DETERMINANT, hol))
            TEST_FUNCTION_FAIL("fixed examples, hol = %d\n", hol);

    for (i = 0; i < 300 * flint_test_multiplier(); i++)
    {
        /* moduli from 2 bits (2 or 3) to 64 bits */
        const ulong prime = n_randprime(state, 2 + n_randint(state, 63), 1);
        /* kind of A: 0 random, 1 random with A(0) singular, 2 singular, 3 a
         * few rows of large degree; kind of b: 0 random, 1 zero, 2 v A for
         * some random v (then den = 1 if A is nonsingular) */
        const slong kindA = n_randint(state, 4);
        const slong kindb = _linsolvetest_rand_kindb(state);
        slong n, len;
        _linsolvetest_rand_dims(&n, &len, state);

        nmod_poly_mat_t A, b, u, u0, v, w;
        nmod_poly_t den, den0, f;
        nmod_poly_mat_init(A, n, n, prime);
        nmod_poly_mat_init(b, 1, n, prime);
        nmod_poly_mat_init(u, 1, n, prime);
        nmod_poly_mat_init(u0, 1, n, prime);
        nmod_poly_mat_init(v, 1, n, prime);
        nmod_poly_mat_init(w, 1, n, prime);
        nmod_poly_init(den, prime);
        nmod_poly_init(den0, prime);
        nmod_poly_init(f, prime);

        _linsolvetest_rand_A(A, state, len, kindA == 2);

        if (kindA == 1 && n > 0)
        {
            /* column j multiplied by x, and sometimes by x + c */
            const slong j = n_randint(state, n);
            for (slong k = 0; k < n; k++)
                nmod_poly_shift_left(nmod_poly_mat_entry(A, k, j), nmod_poly_mat_entry(A, k, j), 1);
            if (n_randint(state, 2))
            {
                nmod_poly_set_coeff_ui(f, 1, 1);
                nmod_poly_set_coeff_ui(f, 0, n_randint(state, prime));
                for (slong k = 0; k < n; k++)
                    nmod_poly_mul(nmod_poly_mat_entry(A, n - 1 - k, k), nmod_poly_mat_entry(A, n - 1 - k, k), f);
            }
        }

        if (kindA == 3 && n > 0)
            _linsolvetest_long_rows(A, state, len);

        const slong lenb = _linsolvetest_rand_b(b, v, A, state, kindb, len);

        const int nonsingular = (nmod_poly_mat_rank(A) == n);

        /* reference: via_kernel */
        const int expected = nmod_poly_mat_linsolve_via_kernel(w, den0, A, b);

        /* dm: den_method; hol: series by Dixon's lifting or by high-order
         * lifting */
        for (int k = 0; k < 4; k++)
        {
            const int dm = k & 1;
            const int hol = k >> 1;
            const linsolve_den_t den_method = (dm == 0) ? NMOD_POLY_MAT_LINSOLVE_RECONSTRUCTION
                                                        : NMOD_POLY_MAT_LINSOLVE_DETERMINANT;

            /* return value, correctness, and u and den not modified in case
             * of failure; if A is nonsingular, the solution is unique */
            nmod_poly_mat_randtest(u, state, 3);
            nmod_poly_randtest(den, state, 3);
            nmod_poly_mat_set(u0, u);
            nmod_poly_set(f, den);
            result = (_lsseriestest_solve(u, den, A, b, den_method, hol) == expected);
            if (result && expected)
                result = _linsolvetest_check(u, den, A, b);
            if (result && expected && nonsingular)
                result = nmod_poly_mat_equal(u, w) && nmod_poly_equal(den, den0);
            if (result && expected && nonsingular && kindb == 2)
                result = nmod_poly_mat_equal(u, v) && nmod_poly_is_one(den);
            if (result && !expected)
                result = nmod_poly_mat_equal(u, u0) && nmod_poly_equal(den, f);

            /* aliasing u and b */
            if (result)
            {
                nmod_poly_mat_set(u0, b);
                result = (_lsseriestest_solve(u0, f, A, u0, den_method, hol) == expected)
                         && nmod_poly_mat_equal(u0, expected ? u : b)
                         && nmod_poly_equal(f, den);
            }

            if (!result)
                TEST_FUNCTION_FAIL("prime = %wu, kindA = %wd, kindb = %wd, n = %wd, len = %wd, lenb = %wd, den_method = %d, hol = %d, nonsingular = %d\n",
                                   prime, kindA, kindb, n, len, lenb, dm, hol, nonsingular);
        }

        nmod_poly_mat_clear(A);
        nmod_poly_mat_clear(b);
        nmod_poly_mat_clear(u);
        nmod_poly_mat_clear(u0);
        nmod_poly_mat_clear(v);
        nmod_poly_mat_clear(w);
        nmod_poly_clear(den);
        nmod_poly_clear(den0);
        nmod_poly_clear(f);
    }

    TEST_FUNCTION_END(state);
}
