/*
    Copyright (C) 2025 Vincent Neiger, Éric Schost

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

int test_mat_mulmid(ulong prime, nmod_poly_mat_t pmat1, nmod_poly_mat_t pmat2, slong nlo, slong nhi)
{
    nmod_poly_mat_t res_true, res;
    nmod_poly_mat_init(res_true, pmat1->r, pmat2->c, prime);
    nmod_poly_mat_init(res, pmat1->r, pmat2->c, prime);

    /* most naive way */
    if (nlo >= nhi)
        nmod_poly_mat_zero(res_true);
    else
    {
        nmod_poly_mat_mul(res_true, pmat1, pmat2);
        nmod_poly_mat_shift_right(res_true, res_true, nlo);
        nmod_poly_mat_truncate(res_true, nhi - nlo);
    }

    /* mulmid */
    nmod_poly_mat_mulmid(res, pmat1, pmat2, nlo, nhi);

    int check = nmod_poly_mat_equal(res_true, res);

    nmod_poly_mat_clear(res_true);
    nmod_poly_mat_clear(res);

    return check;
}

TEST_FUNCTION_START(nmod_poly_mat_mulmid, state)
{
    int i, result;

    for (i = 0; i < 300 * flint_test_multiplier(); i++)
    {
        ulong bits = 2 + n_randint(state, 63);
        ulong prime = n_randprime(state, bits, 1);

        ulong rdim = 1 + n_randint(state, 20);
        ulong idim = 1 + n_randint(state, 20);
        ulong cdim = 1 + n_randint(state, 20);

        /* constraints: */
        slong len1, len2, nlo, nhi;
        if (i < 100 * flint_test_multiplier())
        {
            /* tests on random parameters */
            nlo = n_randint(state, 50);
            nhi = n_randint(state, 100);
            len1 = n_randint(state, 50);
            len2 = n_randint(state, 50);
        }
        else if (i < 200 * flint_test_multiplier())
        {
            /* tests on "classical" parameters for transposed multiplication: */
            /* nlo < nhi, len1 <= nlo+1, len2 <= nhi */
            nlo = n_randint(state, 50);
            nhi = nlo+1 + n_randint(state, 50);
            len1 = n_randint(state, nlo+2);
            len2 = n_randint(state, nhi+1);
        }
        else
        {
            /* same as above, permuted: len1 <= nhi, len2 <= nlo+1 */
            nlo = n_randint(state, 50);
            nhi = nlo+1 + n_randint(state, 50);
            len1 = n_randint(state, nhi+1);
            len2 = n_randint(state, nlo+2);
        }

        nmod_poly_mat_t pmat1, pmat2;
        nmod_poly_mat_init(pmat1, rdim, idim, prime);
        nmod_poly_mat_randtest(pmat1, state, len1);
        nmod_poly_mat_init(pmat2, idim, cdim, prime);
        nmod_poly_mat_randtest(pmat2, state, len2);

        result = test_mat_mulmid(prime, pmat1, pmat2, nlo, nhi);

        if (!result)
            TEST_FUNCTION_FAIL(
                    "prime = %wu, rdim = %wu, idim = %wu, cdim = %wu\n"
                    "bits = %wu, len1 = %wd, len2 = %wd, nlo = %wd, nhi = %wd\n",
                    prime, rdim, idim, cdim, bits, len1, len2, nlo, nhi);

        nmod_poly_mat_clear(pmat1);
        nmod_poly_mat_clear(pmat2);
    }

    TEST_FUNCTION_END(state);
}
