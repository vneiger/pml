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

#include "nmod_poly_mat_utils.h"

/* checks that all entries are normalized (a nonzero leading coefficient) */
static int _shifttrunctest_is_normalised(const nmod_poly_mat_t c)
{
    for (slong i = 0; i < c->r; i++)
        for (slong j = 0; j < c->c; j++)
        {
            const nmod_poly_struct * p = nmod_poly_mat_entry(c, i, j);
            if (p->length > 0 && p->coeffs[p->length - 1] == 0)
                return 0;
        }
    return 1;
}

TEST_FUNCTION_START(nmod_poly_mat_shift_right_trunc, state)
{
    int i, result;

    for (i = 0; i < 2000 * flint_test_multiplier(); i++)
    {
        const ulong prime = n_randprime(state, 2 + n_randint(state, 63), 1);
        const slong r = n_randint(state, 4);
        const slong s = n_randint(state, 4);
        const slong len_in = n_randint(state, 30);
        /* shift: possibly beyond the length; truncation length: possibly
         * nonpositive, or beyond the remaining length */
        const slong k = n_randint(state, 35);
        const slong len = (slong) n_randint(state, 35) - 3;

        nmod_poly_mat_t a, b, c;
        nmod_poly_mat_init(a, r, s, prime);
        nmod_poly_mat_init(b, r, s, prime);
        nmod_poly_mat_init(c, r, s, prime);
        nmod_poly_mat_randtest(a, state, len_in);
        /* sometimes zero coefficients just below the truncation order, to
         * check normalization */
        if (n_randint(state, 4) == 0 && len > 0)
            for (slong ii = 0; ii < r; ii++)
                for (slong jj = 0; jj < s; jj++)
                    for (slong d = 0; d < 3 && d < len; d++)
                        nmod_poly_set_coeff_ui(nmod_poly_mat_entry(a, ii, jj), k + len - 1 - d, 0);

        /* reference: shift, then truncate */
        nmod_poly_mat_shift_right(b, a, k);
        nmod_poly_mat_truncate(b, FLINT_MAX(len, 0));

        nmod_poly_mat_randtest(c, state, 5);
        nmod_poly_mat_shift_right_trunc(c, a, k, len);
        result = nmod_poly_mat_equal(b, c) && _shifttrunctest_is_normalised(c);

        /* aliasing */
        if (result)
        {
            nmod_poly_mat_set(c, a);
            nmod_poly_mat_shift_right_trunc(c, c, k, len);
            result = nmod_poly_mat_equal(b, c) && _shifttrunctest_is_normalised(c);
        }

        if (!result)
            TEST_FUNCTION_FAIL("prime = %wu, r = %wd, s = %wd, len_in = %wd, k = %wd, len = %wd\n",
                               prime, r, s, len_in, k, len);

        nmod_poly_mat_clear(a);
        nmod_poly_mat_clear(b);
        nmod_poly_mat_clear(c);
    }

    TEST_FUNCTION_END(state);
}
