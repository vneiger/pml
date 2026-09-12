/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

/**
 * Tests the nmod_poly_mat -> nmod_mat_poly conversion.
 *
 * For each random instance, the blocked implementation is run with every
 * option accepted by _nmod_mat_poly_set_trunc_from_poly_mat, plus the default
 * option in nmod_mat_poly_set_trunc_from_poly_mat, and each result is compared
 * with the straightforward conversion (*_naive below).
 *
 * Unavailable kernels fall back to an available one, so the loop below is
 * portable; on a given build several columns may test the same code.
 *
 * The instances deliberately cover:
 *   - dimensions and orders on both sides of the block sizes 4 and 8, and
 *     the degenerate 0 x c, r x 0, order 0 cases;
 *   - input entries of unequal lengths, including zero entries (so that the
 *     zero-padding above each entry's length is exercised, both inside a
 *     block and in the residual bands);
 *   - a destination that already holds data, of larger and of smaller length
 *     than the requested order, so that stale coefficients must be
 *     overwritten;
 *   - truncation orders below, at, and above the length of the input;
 *   - a window of a polynomial matrix as input, whose row stride differs
 *     from its number of columns;
 *   - a destination whose own row stride is larger than its number of
 *     columns, which takes the row-by-row path.
 */

#include <flint/test_helpers.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>

#include "nmod_mat_poly.h"
#include "nmod_mat_poly_extra/impl.h"

/* reference point for checks */
static void _set_trunc_from_poly_mat_naive(nmod_mat_poly_t matp,
                                           const nmod_poly_mat_t pmat,
                                           slong order)
{
    const slong len = nmod_poly_mat_max_length(pmat);
    if (order > len)
        order = len;

    nmod_mat_poly_fit_length(matp, order);
    _nmod_mat_poly_set_length(matp, order);

    for (slong k = 0; k < order; k++)
        for (slong i = 0; i < matp->r; i++)
            for (slong j = 0; j < matp->c; j++)
                nmod_mat_poly_entry(matp, k, i, j) = nmod_poly_get_coeff_ui(nmod_poly_mat_entry(pmat, i, j), k);

    if (order < len)
        _nmod_mat_poly_normalise(matp);
}

/* full equality, including length */
static int _matp_equal(const nmod_mat_poly_t a, const nmod_mat_poly_t b)
{
    nmod_mat_t ca, cb;
    if (a->length != b->length || a->r != b->r || a->c != b->c)
        return 0;
    for (slong k = 0; k < a->length; k++)
    {
        nmod_mat_poly_coeff_attach(ca, a, k);
        nmod_mat_poly_coeff_attach(cb, b, k);
        if (! nmod_mat_equal(ca, cb))
            return 0;
    }
    return 1;
}

/* also check against the definition, entry by entry */
static int _matp_matches_pmat(const nmod_mat_poly_t matp,
                              const nmod_poly_mat_t pmat,
                              slong order)
{
    for (slong k = 0; k < matp->length; k++)
        for (slong i = 0; i < matp->r; i++)
            for (slong j = 0; j < matp->c; j++)
                if (nmod_mat_poly_entry(matp, k, i, j)
                        != nmod_poly_get_coeff_ui(nmod_poly_mat_entry(pmat, i, j), k))
                    return 0;
    /* nothing of degree >= length, up to the truncation order */
    for (slong k = matp->length; k < order; k++)
        for (slong i = 0; i < matp->r; i++)
            for (slong j = 0; j < matp->c; j++)
                if (nmod_poly_get_coeff_ui(nmod_poly_mat_entry(pmat, i, j), k) != UWORD(0))
                    return 0;
    return 1;
}

static int core_test_set_from_poly_mat(const nmod_poly_mat_t pmat,
                                       slong order,
                                       flint_rand_t state)
{
    const slong r = nmod_poly_mat_nrows(pmat);
    const slong c = nmod_poly_mat_ncols(pmat);
    const ulong prime = nmod_poly_mat_modulus(pmat);
    int ok = 1;

    nmod_mat_poly_t ref, res;
    nmod_mat_poly_init(ref, r, c, prime);
    nmod_mat_poly_init(res, r, c, prime);

    _set_trunc_from_poly_mat_naive(ref, pmat, order);

    if (! _matp_matches_pmat(ref, pmat, order))
        ok = 0;

    for (int kern = -1; ok && kern <= NMOD_MAT_POLY_CONV_VEC8; kern++)
        for (int dmaj = -1; ok && dmaj <= 1; dmaj++)
        {
            /* leave some stale data behind, of a length unrelated to order */
            nmod_mat_poly_rand(res, state, n_randint(state, 2 * order + 3));

            _nmod_mat_poly_set_trunc_from_poly_mat(res, pmat, order, kern, dmaj);

            if (! _matp_equal(ref, res))
            {
                flint_printf("failure with kern = %d, dmaj = %d\n", kern, dmaj);
                ok = 0;
            }
        }

    /* the dispatcher, and the non-truncated entry point */
    if (ok)
    {
        nmod_mat_poly_rand(res, state, n_randint(state, 2 * order + 3));
        nmod_mat_poly_set_trunc_from_poly_mat(res, pmat, order);
        ok = _matp_equal(ref, res);
    }

    if (ok && order >= nmod_poly_mat_max_length(pmat))
    {
        nmod_mat_poly_rand(res, state, n_randint(state, 2 * order + 3));
        nmod_mat_poly_set_from_poly_mat(res, pmat);
        ok = _matp_equal(ref, res);
    }

    /* same, with a destination whose row stride is not the number of columns:
       this takes the row-by-row path of the implementation */
    if (ok && r > 0 && c > 0)
    {
        nmod_mat_poly_t wide;
        nmod_mat_poly_init(wide, r, c, prime);
        wide->stride = c + 1 + n_randint(state, 8);
        nmod_mat_poly_set_trunc_from_poly_mat(wide, pmat, order);
        for (slong k = 0; ok && k < wide->length; k++)
            for (slong i = 0; ok && i < r; i++)
                for (slong j = 0; ok && j < c; j++)
                    if (nmod_mat_poly_entry(wide, k, i, j) != nmod_mat_poly_entry(ref, k, i, j))
                        ok = 0;
        if (ok && wide->length != ref->length)
            ok = 0;
        nmod_mat_poly_clear(wide);
    }

    nmod_mat_poly_clear(ref);
    nmod_mat_poly_clear(res);
    return ok;
}

TEST_FUNCTION_START(nmod_mat_poly_set_from_poly_mat, state)
{
    int result;

    for (slong iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        const slong nbits = 2 + n_randint(state, 63);
        const ulong prime = n_randprime(state, nbits, 1);

        /* dimensions and lengths straddling the block sizes 4 and 8 */
        const slong rdim = n_randint(state, 12);
        const slong cdim = n_randint(state, 12);
        const slong len = n_randint(state, 26);

        nmod_poly_mat_t pmat;
        nmod_poly_mat_init(pmat, rdim, cdim, prime);

        /* entries of unequal lengths, some of them zero */
        for (slong i = 0; i < rdim; i++)
            for (slong j = 0; j < cdim; j++)
            {
                if (n_randint(state, 6) == 0)
                    nmod_poly_zero(nmod_poly_mat_entry(pmat, i, j));
                else
                    nmod_poly_randtest(nmod_poly_mat_entry(pmat, i, j), state,
                                       1 + n_randint(state, len + 1));
            }

        /* order below, at, and above the length of the input */
        const slong maxlen = nmod_poly_mat_max_length(pmat);
        const slong order = n_randint(state, maxlen + 3);

        result = core_test_set_from_poly_mat(pmat, order, state);

        if (result && rdim > 1 && cdim > 1)
        {
            /* a window: its row stride differs from its column count */
            nmod_poly_mat_t win;
            nmod_poly_mat_window_init(win, pmat, 0, 0, rdim - 1, cdim - 1);
            result = core_test_set_from_poly_mat(win, order, state);
            nmod_poly_mat_window_clear(win);
        }

        if (!result)
            TEST_FUNCTION_FAIL("prime = %wu, rdim = %wd, cdim = %wd, len = %wd, order = %wd\n",
                               prime, rdim, cdim, len, order);

        nmod_poly_mat_clear(pmat);
    }

    /* larger instances, to reach the blocked path with both schedules */
    {
        const ulong prime = UWORD(1099511627791);
        const slong shapes[4][2] = {{9, 9}, {5, 40}, {40, 5}, {16, 16}};
        const slong orders[3] = {7, 33, 130};

        for (int s = 0; s < 4; s++)
            for (int o = 0; o < 3; o++)
            {
                nmod_poly_mat_t pmat;
                nmod_poly_mat_init(pmat, shapes[s][0], shapes[s][1], prime);
                nmod_poly_mat_rand(pmat, state, orders[o]);
                /* make a few entries shorter, including in the last block */
                for (slong i = 0; i < shapes[s][0]; i++)
                    for (slong j = 0; j < shapes[s][1]; j++)
                        if ((i + j) % 5 == 0)
                            nmod_poly_truncate(nmod_poly_mat_entry(pmat, i, j),
                                               n_randint(state, orders[o] + 1));

                result = core_test_set_from_poly_mat(pmat, orders[o], state)
                      && core_test_set_from_poly_mat(pmat, orders[o] / 2, state);

                if (!result)
                    TEST_FUNCTION_FAIL("large instance %wd x %wd, order %wd\n",
                                       shapes[s][0], shapes[s][1], orders[o]);

                nmod_poly_mat_clear(pmat);
            }
    }

    TEST_FUNCTION_END(state);
}
