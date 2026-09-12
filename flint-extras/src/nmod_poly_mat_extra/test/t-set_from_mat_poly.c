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
 * Tests the nmod_mat_poly -> nmod_poly_mat conversion.
 *
 * For each random instance, the blocked implementation is run with every
 * option accepted by _nmod_poly_mat_set_trunc_from_mat_poly, plus the default
 * options in nmod_poly_mat_set_trunc_from_mat_poly and
 * nmod_poly_mat_set_from_mat_poly, and each result is compared with the
 * straightforward conversion (*_naive below).
 *
 * Unavailable kernels fall back to an available one, so the loop below is
 * portable; on a given build several of its iterations test the same code.
 *
 * The instances deliberately cover:
 *   - dimensions and orders on both sides of the block sizes 4 and 8, and the
 *     degenerate 0 x c, r x 0, order 0 and length 0 cases;
 *   - inputs whose trailing coefficients are zero, entirely or in part, so
 *     that the output entries come out with unequal lengths and the
 *     normalisation is exercised;
 *   - a destination that already holds data, of larger and of smaller length
 *     than the requested order, so that stale coefficients must be overwritten
 *     and stale lengths corrected;
 *   - truncation orders below, at, and above the length of the input;
 *   - an input whose row stride is larger than its number of columns, which
 *     takes the row-by-row path of the implementation;
 *   - a window of a polynomial matrix as destination, whose row stride
 *     differs from its number of columns.
 */

#include <flint/test_helpers.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>

#include "nmod_mat_poly.h"
#include "nmod_poly_mat_utils.h"
#include "nmod_mat_poly_extra/impl.h"
#include "nmod_poly_mat_extra/impl.h"

/* reference point for checks */
static void _set_trunc_from_mat_poly_naive(nmod_poly_mat_t pmat,
                                           const nmod_mat_poly_t matp,
                                           slong order)
{
    if (order > matp->length)
        order = matp->length;

    for (slong i = 0; i < pmat->r; i++)
        for (slong j = 0; j < pmat->c; j++)
            nmod_poly_fit_length(nmod_poly_mat_entry(pmat, i, j), order);

    for (slong k = 0; k < order; k++)
        for (slong i = 0; i < pmat->r; i++)
            for (slong j = 0; j < pmat->c; j++)
                nmod_poly_mat_entry(pmat, i, j)->coeffs[k] = nmod_mat_poly_entry(matp, k, i, j);

    for (slong i = 0; i < pmat->r; i++)
        for (slong j = 0; j < pmat->c; j++)
        {
            _nmod_poly_set_length(nmod_poly_mat_entry(pmat, i, j), order);
            _nmod_poly_normalise(nmod_poly_mat_entry(pmat, i, j));
        }
}

/* full equality, including the length of every entry */
static int _pmat_equal(const nmod_poly_mat_t a, const nmod_poly_mat_t b)
{
    if (a->r != b->r || a->c != b->c)
        return 0;
    for (slong i = 0; i < a->r; i++)
        for (slong j = 0; j < a->c; j++)
        {
            const nmod_poly_struct * pa = nmod_poly_mat_entry(a, i, j);
            const nmod_poly_struct * pb = nmod_poly_mat_entry(b, i, j);
            if (pa->length != pb->length || ! nmod_poly_equal(pa, pb))
                return 0;
        }
    return 1;
}

/* also check against the definition, entry by entry */
static int _pmat_matches_matp(const nmod_poly_mat_t pmat,
                              const nmod_mat_poly_t matp,
                              slong order)
{
    if (order > matp->length)
        order = matp->length;

    for (slong i = 0; i < pmat->r; i++)
        for (slong j = 0; j < pmat->c; j++)
        {
            const nmod_poly_struct * p = nmod_poly_mat_entry(pmat, i, j);
            if (p->length > order)
                return 0;
            for (slong k = 0; k < order; k++)
                if (nmod_poly_get_coeff_ui(p, k) != nmod_mat_poly_entry(matp, k, i, j))
                    return 0;
        }
    return 1;
}

/* fill a destination with stale data of a length unrelated to `order` */
static void _stale(nmod_poly_mat_t pmat, slong order, flint_rand_t state)
{
    for (slong i = 0; i < pmat->r; i++)
        for (slong j = 0; j < pmat->c; j++)
            nmod_poly_randtest(nmod_poly_mat_entry(pmat, i, j), state,
                               n_randint(state, 2 * order + 4));
}

static int core_test_set_from_mat_poly(const nmod_mat_poly_t matp,
                                       slong order,
                                       flint_rand_t state)
{
    const slong r = matp->r;
    const slong c = matp->c;
    const ulong prime = matp->mod.n;
    int ok = 1;

    nmod_poly_mat_t ref, res;
    nmod_poly_mat_init(ref, r, c, prime);
    nmod_poly_mat_init(res, r, c, prime);

    _stale(ref, order, state);
    _set_trunc_from_mat_poly_naive(ref, matp, order);

    if (! _pmat_matches_matp(ref, matp, order))
        ok = 0;

    for (int kern = -1; ok && kern <= NMOD_MAT_POLY_CONV_VEC8; kern++)
        for (int dmaj = -1; ok && dmaj <= 1; dmaj++)
        {
            _stale(res, order, state);

            _nmod_poly_mat_set_trunc_from_mat_poly(res, matp, order, kern, dmaj);

            if (! _pmat_equal(ref, res))
            {
                flint_printf("failure with kern = %d, dmaj = %d\n", kern, dmaj);
                ok = 0;
            }
        }

    /* the dispatcher, and the non-truncated entry point */
    if (ok)
    {
        _stale(res, order, state);
        nmod_poly_mat_set_trunc_from_mat_poly(res, matp, order);
        ok = _pmat_equal(ref, res);
    }

    if (ok && order >= matp->length)
    {
        _stale(res, order, state);
        nmod_poly_mat_set_from_mat_poly(res, matp);
        ok = _pmat_equal(ref, res);
    }

    /* same, into a window, whose row stride is not the number of columns */
    if (ok && r > 1 && c > 1)
    {
        nmod_poly_mat_t big, win;
        nmod_poly_mat_init(big, r + 1, c + 2, prime);
        _stale(big, order, state);
        nmod_poly_mat_window_init(win, big, 1, 2, 1 + r, 2 + c);
        nmod_poly_mat_set_trunc_from_mat_poly(win, matp, order);
        ok = _pmat_equal(ref, win);
        nmod_poly_mat_window_clear(win);
        nmod_poly_mat_clear(big);
    }

    nmod_poly_mat_clear(ref);
    nmod_poly_mat_clear(res);
    return ok;
}

/* random matrix polynomial; `zeros` asks for some coefficients to be zero,
   including trailing ones, so that the output entries get unequal lengths */
static void _rand_matp(nmod_mat_poly_t matp, flint_rand_t state, slong len, int zeros)
{
    nmod_mat_poly_rand(matp, state, len);
    if (zeros)
        for (slong k = 0; k < matp->length; k++)
            for (slong i = 0; i < matp->r; i++)
                for (slong j = 0; j < matp->c; j++)
                    if (n_randint(state, 3) == 0 || (k + 1 == matp->length && (i + j) % 2))
                        nmod_mat_poly_entry(matp, k, i, j) = UWORD(0);
}

TEST_FUNCTION_START(nmod_poly_mat_set_from_mat_poly, state)
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

        nmod_mat_poly_t matp;
        nmod_mat_poly_init(matp, rdim, cdim, prime);
        _rand_matp(matp, state, len, n_randint(state, 2));

        /* order below, at, and above the length of the input */
        const slong order = n_randint(state, matp->length + 3);

        result = core_test_set_from_mat_poly(matp, order, state);

        /* same input with a row stride larger than its number of columns:
           this takes the row-by-row path of the implementation */
        if (result && rdim > 0 && cdim > 0)
        {
            nmod_mat_poly_t wide;
            nmod_mat_poly_init(wide, rdim, cdim, prime);
            wide->stride = cdim + 1 + n_randint(state, 8);
            nmod_mat_poly_fit_length(wide, matp->length);
            _nmod_mat_poly_set_length(wide, matp->length);
            for (slong k = 0; k < matp->length; k++)
                for (slong i = 0; i < rdim; i++)
                    for (slong j = 0; j < cdim; j++)
                        nmod_mat_poly_entry(wide, k, i, j) = nmod_mat_poly_entry(matp, k, i, j);
            result = core_test_set_from_mat_poly(wide, order, state);
            nmod_mat_poly_clear(wide);
        }

        if (!result)
            TEST_FUNCTION_FAIL("prime = %wu, rdim = %wd, cdim = %wd, len = %wd, order = %wd\n",
                               prime, rdim, cdim, len, order);

        nmod_mat_poly_clear(matp);
    }

    /* larger instances, to reach the blocked path with both schedules */
    {
        const ulong prime = UWORD(1099511627791);
        const slong shapes[4][2] = {{9, 9}, {5, 40}, {40, 5}, {16, 16}};
        const slong orders[3] = {7, 33, 130};

        for (int s = 0; s < 4; s++)
            for (int o = 0; o < 3; o++)
            {
                nmod_mat_poly_t matp;
                nmod_mat_poly_init(matp, shapes[s][0], shapes[s][1], prime);
                _rand_matp(matp, state, orders[o], 1);

                result = core_test_set_from_mat_poly(matp, orders[o], state)
                      && core_test_set_from_mat_poly(matp, orders[o] / 2, state);

                if (!result)
                    TEST_FUNCTION_FAIL("large instance %wd x %wd, order %wd\n",
                                       shapes[s][0], shapes[s][1], orders[o]);

                nmod_mat_poly_clear(matp);
            }
    }

    TEST_FUNCTION_END(state);
}
