/*
    Copyright (C) 2026 Vincent Neiger

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

#include "nmod_poly_mat_multiply.h"

TEST_FUNCTION_START(nmod_poly_mat_mulmid_precomp_right, state)
{
    int i, result = 1;

    for (i = 0; i < 300 * flint_test_multiplier(); i++)
    {
        const ulong prime = n_randprime(state, 2 + n_randint(state, 63), 1);

        /* dimensions: mostly small, some row vectors, sometimes zero */
        const slong m = (n_randint(state, 4) == 0) ? 1 : (slong) n_randint(state, 12);
        const slong k = n_randint(state, 12);
        const slong n = n_randint(state, 12);

        /* lengths and window: random, or the shapes of Dixon's lifting
         * (low product by a long operand, and middle product by a short one
         * of a long left operand) */
        slong len1, len2, nlo, nhi;
        switch (n_randint(state, 3))
        {
            case 0:
                len1 = n_randint(state, 40);
                len2 = n_randint(state, 40);
                nlo = n_randint(state, 40);
                nhi = n_randint(state, 80);
                break;
            case 1:
                len2 = 1 + n_randint(state, 60);
                len1 = 1 + n_randint(state, 20);
                nlo = 0;
                nhi = len2;
                break;
            default:
                len1 = 1 + n_randint(state, 60);
                len2 = 1 + n_randint(state, 20);
                nlo = len1 - 1;
                nhi = len1 + len2 - 1;
        }
        const int method = (int) n_randint(state, 3) - 1;

        nmod_poly_mat_t pmat1, pmat2, res, expected;
        nmod_poly_mat_init(pmat2, k, n, prime);
        nmod_poly_mat_randtest(pmat2, state, len2);

        nmod_poly_mat_mulmid_precomp_right_t pre;
        _nmod_poly_mat_mulmid_precomp_right_init(pre, pmat2, len1, nlo, nhi, method);

        /* several products with the same precomputation: left operands
         * possibly shorter than len1, or longer (fallback), windows
         * possibly shorter, or longer than at initialization (reduced),
         * and aliasing res and pmat1 */
        for (slong j = 0; j < 4 && result; j++)
        {
            const slong l1 = (n_randint(state, 6) == 0) ? len1 + 1 + (slong) n_randint(state, 5)
                                                        : (slong) n_randint(state, len1 + 1);
            slong h = (n_randint(state, 4) == 0) ? nhi + (slong) n_randint(state, 5)
                                                 : nlo + (slong) n_randint(state, FLINT_MAX(nhi - nlo, 0) + 2);
            const slong hh = FLINT_MIN(h, pre->nhi);

            nmod_poly_mat_init(pmat1, m, k, prime);
            nmod_poly_mat_init(res, m, n, prime);
            nmod_poly_mat_init(expected, m, n, prime);
            nmod_poly_mat_randtest(pmat1, state, l1);
            nmod_poly_mat_randtest(res, state, 3);

            /* reference: full product, then shift and truncation (not
             * nmod_poly_mat_mulmid, which the products without stored
             * transforms call) */
            _nmod_poly_mat_mulmid_naive(expected, pmat1, nmod_poly_mat_max_length(pmat1),
                                        pmat2, nmod_poly_mat_max_length(pmat2), nlo, hh);

            if (n_randint(state, 4) == 0 && k == n)
            {
                nmod_poly_mat_set(res, pmat1);
                nmod_poly_mat_mulmid_precomp_right(res, res, pre, h);
            }
            else
                nmod_poly_mat_mulmid_precomp_right(res, pmat1, pre, h);
            result = nmod_poly_mat_equal(res, expected);

            if (!result)
                TEST_FUNCTION_FAIL("prime = %wu, m = %wd, k = %wd, n = %wd, len1 = %wd, len2 = %wd, "
                                   "nlo = %wd, nhi = %wd, method = %d, fft = %d, l1 = %wd, h = %wd\n",
                                   prime, m, k, n, len1, len2, nlo, nhi, method,
                                   nmod_poly_mat_mulmid_precomp_right_is_fft(pre), l1, h);

            nmod_poly_mat_clear(pmat1);
            nmod_poly_mat_clear(res);
            nmod_poly_mat_clear(expected);
        }

        nmod_poly_mat_mulmid_precomp_right_clear(pre);
        nmod_poly_mat_clear(pmat2);
    }

    TEST_FUNCTION_END(state);
}
