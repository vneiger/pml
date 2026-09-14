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

#include "nmod_extra.h"   /* for nmod_find_root */
#include "nmod_poly_mat_multiply.h"
#include "nmod_poly_mat_extra/impl.h"

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

    /* the geometric variant directly, with a memory bound small enough to
       force the grouping of rows and columns, several threads, and
       entries of unequal lengths with zero entries, rows and columns
       sprinkled in -- what the approximant basis algorithms feed it */
    {
        const slong nthreads_save = flint_get_num_threads();

        for (i = 0; i < 200 * flint_test_multiplier(); i++)
        {
            const ulong bits = 16 + n_randint(state, 49);
            const ulong prime = n_randprime(state, bits, 1);
            const slong rdim = 1 + n_randint(state, 20);
            const slong idim = 1 + n_randint(state, 20);
            const slong cdim = 1 + n_randint(state, 20);
            const slong nlo = n_randint(state, 40);
            const slong nhi = nlo + 1 + n_randint(state, 40);
            const int swap = n_randint(state, 2);
            /* one operand fits below nlo+1, the other below nhi */
            const slong lenX = 1 + n_randint(state, nlo + 1);
            const slong lenY = 1 + n_randint(state, nhi);
            const slong len1 = swap ? lenY : lenX;
            const slong len2 = swap ? lenX : lenY;
            const ulong membytes = (n_randint(state, 3) == 0)
                                   ? 0 : 8 * (1 + n_randint(state, 4000));
            slong r, c;

            flint_set_num_threads(1 + n_randint(state, 4));

            nmod_poly_mat_t pmat1, pmat2, res, res_true;
            nmod_poly_mat_init(pmat1, rdim, idim, prime);
            nmod_poly_mat_init(pmat2, idim, cdim, prime);
            nmod_poly_mat_init(res, rdim, cdim, prime);
            nmod_poly_mat_init(res_true, rdim, cdim, prime);

            for (r = 0; r < rdim; r++)
                for (c = 0; c < idim; c++)
                    nmod_poly_randtest(nmod_poly_mat_entry(pmat1, r, c), state,
                                       n_randint(state, len1 + 1));
            for (r = 0; r < idim; r++)
                for (c = 0; c < cdim; c++)
                    nmod_poly_randtest(nmod_poly_mat_entry(pmat2, r, c), state,
                                       n_randint(state, len2 + 1));
            if (n_randint(state, 3) == 0)   /* a zero row of pmat1 */
                for (c = 0; c < idim; c++)
                    nmod_poly_zero(nmod_poly_mat_entry(pmat1, n_randint(state, rdim), c));
            if (n_randint(state, 3) == 0)   /* a zero column of pmat2 */
                for (r = 0; r < idim; r++)
                    nmod_poly_zero(nmod_poly_mat_entry(pmat2, r, n_randint(state, cdim)));
            nmod_poly_mat_randtest(res, state, 5);   /* stale data in the output */

            nmod_poly_mat_mul(res_true, pmat1, pmat2);
            nmod_poly_mat_shift_right(res_true, res_true, nlo);
            nmod_poly_mat_truncate(res_true, nhi - nlo);

            {
                nmod_t mod;
                nmod_geometric_progression_t G;
                ulong w;

                nmod_init(&mod, prime);
                w = nmod_find_root(2 * nhi, mod);
                _nmod_geometric_progression_init_function(G, w, nhi, mod, UWORD(3));
                _nmod_poly_mat_mulmid_geometric_precomp_bounded(res, pmat1, len1,
                                                    pmat2, len2, nlo, nhi, G, membytes);
                nmod_geometric_progression_clear(G);
            }

            result = nmod_poly_mat_equal(res_true, res);

            if (!result)
                TEST_FUNCTION_FAIL(
                        "geometric: prime = %wu, rdim = %wd, idim = %wd, cdim = %wd\n"
                        "len1 = %wd, len2 = %wd, nlo = %wd, nhi = %wd, "
                        "membytes = %wu, threads = %wd\n",
                        prime, rdim, idim, cdim, len1, len2, nlo, nhi, membytes,
                        flint_get_num_threads());

            nmod_poly_mat_clear(pmat1);
            nmod_poly_mat_clear(pmat2);
            nmod_poly_mat_clear(res);
            nmod_poly_mat_clear(res_true);
        }

        flint_set_num_threads(nthreads_save);
    }

    TEST_FUNCTION_END(state);
}
