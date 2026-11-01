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
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>
#include <flint/test_helpers.h>

#include "nmod_poly_mat_forms.h"
#include "nmod_poly_mat_utils.h"

/* checks that rmat is the reverse of pmat with bounds rhi[i] + chi[j] (NULL
 * meaning zero), truncated at order (no truncation if order < 0), directly
 * from the definition */
static int _revtest_check(const nmod_poly_mat_t rmat, const nmod_poly_mat_t pmat,
                          const slong * rhi, const slong * chi, slong order)
{
    for (slong i = 0; i < pmat->r; i++)
        for (slong j = 0; j < pmat->c; j++)
        {
            const slong hi = (rhi ? rhi[i] : 0) + (chi ? chi[j] : 0);
            const nmod_poly_struct * p = nmod_poly_mat_entry(pmat, i, j);
            const nmod_poly_struct * r = nmod_poly_mat_entry(rmat, i, j);
            /* bound on the length of r */
            slong len = FLINT_MAX(0, hi + 1);
            if (order >= 0 && order < len)
                len = order;
            if (r->length > len)
                return 0;
            for (slong k = 0; k < len; k++)
                if (nmod_poly_get_coeff_ui(r, k) != nmod_poly_get_coeff_ui(p, hi - k))
                    return 0;
        }
    return 1;
}

TEST_FUNCTION_START(nmod_poly_mat_reverse, state)
{
    int i, result;

    for (i = 0; i < 500 * flint_test_multiplier(); i++)
    {
        const ulong prime = n_randprime(state, 2 + n_randint(state, 62), 1);
        const slong m = n_randint(state, 6);
        const slong n = n_randint(state, 6);
        const slong len = n_randint(state, 10);

        nmod_poly_mat_t A, R, B;
        nmod_poly_mat_init(A, m, n, prime);
        nmod_poly_mat_init(R, m, n, prime);
        nmod_poly_mat_init(B, m, n, prime);
        nmod_poly_mat_randtest(A, state, len);

        slong * rhi = flint_malloc(FLINT_MAX(1, m) * sizeof(slong));
        slong * chi = flint_malloc(FLINT_MAX(1, n) * sizeof(slong));
        for (slong k = 0; k < m; k++)
            rhi[k] = (slong) n_randint(state, 2 * len + 4) - 3;
        for (slong k = 0; k < n; k++)
            chi[k] = (slong) n_randint(state, 2 * len + 4) - 3;
        const slong hi = (slong) n_randint(state, 2 * len + 4) - 3;
        const slong order = (slong) n_randint(state, len + 4) - 1;

        /* uniform bound */
        nmod_poly_mat_reverse(R, A, hi);
        slong * hirow = flint_malloc(FLINT_MAX(1, m) * sizeof(slong));
        for (slong k = 0; k < m; k++)
            hirow[k] = hi;
        result = _revtest_check(R, A, hirow, NULL, -1);

        /* with respect to the degree; and reversing twice with the same
         * bound hi >= deg(A) gives A back */
        if (result)
        {
            nmod_poly_mat_reverse_degree(R, A);
            const slong deg = nmod_poly_mat_degree(A);
            for (slong k = 0; k < m; k++)
                hirow[k] = deg;
            result = _revtest_check(R, A, hirow, NULL, -1);
        }
        if (result && hi >= nmod_poly_mat_degree(A))
        {
            nmod_poly_mat_reverse(R, A, hi);
            nmod_poly_mat_reverse(B, R, hi);
            result = nmod_poly_mat_equal(A, B);
        }

        /* row-wise and column-wise bounds, with aliasing */
        if (result)
        {
            nmod_poly_mat_set(B, A);
            nmod_poly_mat_reverse_rowwise(B, B, rhi);
            result = _revtest_check(B, A, rhi, NULL, -1);
        }
        if (result)
        {
            nmod_poly_mat_reverse_columnwise(R, A, chi);
            result = _revtest_check(R, A, NULL, chi, -1);
        }

        /* general version, truncated, with and without aliasing */
        if (result)
        {
            nmod_poly_mat_reverse_trunc(R, A, rhi, chi, order);
            result = _revtest_check(R, A, rhi, chi, order);
            nmod_poly_mat_set(B, A);
            nmod_poly_mat_reverse_trunc(B, B, rhi, chi, order);
            result = result && nmod_poly_mat_equal(B, R);
        }

        /* shifted reversal: constant coefficient is the shifted leading
         * matrix (on nonzero rows) */
        if (result && m > 0 && n > 0)
        {
            slong * d = flint_malloc(m * sizeof(slong));
            slong * ms = flint_malloc(n * sizeof(slong));
            nmod_poly_mat_row_degree(d, A, chi);
            for (slong k = 0; k < n; k++)
                ms[k] = - chi[k];
            nmod_poly_mat_reverse_trunc(R, A, d, ms, 1);
            nmod_mat_t L0, L;
            nmod_mat_init(L0, m, n, prime);
            nmod_mat_init(L, m, n, prime);
            nmod_poly_mat_get_coeff_mat(L0, R, 0);
            nmod_poly_mat_leading_matrix(L, A, chi, ROW_LOWER);
            for (slong ii = 0; result && ii < m; ii++)
            {
                int zero_row = 1;
                for (slong jj = 0; jj < n; jj++)
                    if (!nmod_poly_is_zero(nmod_poly_mat_entry(A, ii, jj)))
                        zero_row = 0;
                for (slong jj = 0; !zero_row && jj < n; jj++)
                    if (nmod_mat_entry(L0, ii, jj) != nmod_mat_entry(L, ii, jj))
                        result = 0;
            }
            /* truncation at order 1 leaves constant polynomials */
            result = result && (nmod_poly_mat_degree(R) <= 0);
            nmod_mat_clear(L0);
            nmod_mat_clear(L);
            flint_free(ms);
            flint_free(d);
        }

        if (!result)
            TEST_FUNCTION_FAIL("prime = %wu, m = %wd, n = %wd, len = %wd, hi = %wd, order = %wd\n",
                               prime, m, n, len, hi, order);

        flint_free(hirow);
        flint_free(rhi);
        flint_free(chi);
        nmod_poly_mat_clear(A);
        nmod_poly_mat_clear(R);
        nmod_poly_mat_clear(B);
    }

    TEST_FUNCTION_END(state);
}
