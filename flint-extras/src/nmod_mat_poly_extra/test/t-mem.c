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
 * Checks the storage invariants of nmod_mat_poly:
 *   - the row stride is the number of columns, and the entry arrays of all
 *     coefficients are aligned on NMOD_MAT_POLY_ALIGN bytes;
 *   - a coefficient that comes into existence, whether the length grows from
 *     scratch or after a shrink, is zero;
 *   - growing the length, which reallocates the array of coefficient
 *     pointers, leaves the already existing coefficients untouched;
 *   - the nmod_mat_poly_entry macro and the nmod_mat_t view returned by
 *     nmod_mat_poly_coeff_attach address the same words;
 *   - truncate, realloc down and up, zero, one, and clear behave.
 */

#include <stdint.h>

#include <flint/test_helpers.h>
#include <flint/nmod_mat.h>

#include "nmod_mat_poly.h"

static int _check_shape(const nmod_mat_poly_t matp)
{
    if (matp->stride != matp->c)
        return 0;
    for (slong k = 0; k < matp->length; k++)
    {
        if (matp->r == 0 || matp->c == 0)
        {
            if (matp->coeffs[k] != NULL)
                return 0;
        }
        else if (matp->coeffs[k] == NULL
                 || ((uintptr_t) matp->coeffs[k]) % NMOD_MAT_POLY_ALIGN != 0)
            return 0;
    }
    return 1;
}

static int _check_zero_from(const nmod_mat_poly_t matp, slong from)
{
    for (slong k = from; k < matp->length; k++)
        for (slong i = 0; i < matp->r; i++)
            for (slong j = 0; j < matp->c; j++)
                if (nmod_mat_poly_entry(matp, k, i, j) != UWORD(0))
                    return 0;
    return 1;
}

/* the entry macro and the nmod_mat_t view must agree */
static int _check_view(const nmod_mat_poly_t matp)
{
    nmod_mat_t cmat;
    for (slong k = 0; k < matp->length; k++)
    {
        nmod_mat_poly_coeff_attach(cmat, matp, k);
        if (cmat->r != matp->r || cmat->c != matp->c
                || cmat->stride != matp->stride || cmat->mod.n != matp->mod.n)
            return 0;
        if (nmod_mat_poly_coeff_ptr(matp, k) != cmat->entries)
            return 0;
        for (slong i = 0; i < matp->r; i++)
            for (slong j = 0; j < matp->c; j++)
                if (&nmod_mat_entry(cmat, i, j) != nmod_mat_poly_entry_ptr(matp, k, i, j))
                    return 0;
    }
    return 1;
}

TEST_FUNCTION_START(nmod_mat_poly_mem, state)
{
    for (slong iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        const ulong prime = n_randprime(state, 2 + n_randint(state, 63), 1);
        const slong r = n_randint(state, 10);
        const slong c = n_randint(state, 10);
        const slong len = n_randint(state, 20);

        nmod_mat_poly_t matp;
        nmod_mat_poly_init(matp, r, c, prime);

        if (matp->stride != c || matp->alloc != 0 || matp->length != 0)
            TEST_FUNCTION_FAIL("init: r = %wd, c = %wd\n", r, c);

        /* growing from scratch: everything must be zero and aligned */
        nmod_mat_poly_fit_length(matp, len);
        _nmod_mat_poly_set_length(matp, len);
        if (! _check_shape(matp) || ! _check_zero_from(matp, 0) || ! _check_view(matp))
            TEST_FUNCTION_FAIL("grow from scratch: r = %wd, c = %wd, len = %wd\n", r, c, len);

        /* fill with data, remember it, then grow again: the reallocation of
           the pointer array must not disturb the existing coefficients */
        nmod_mat_poly_rand(matp, state, len);
        const slong len1 = matp->length;
        nmod_mat_poly_t copy;
        nmod_mat_poly_init(copy, r, c, prime);
        nmod_mat_poly_set(copy, matp);

        const slong len2 = len1 + 1 + n_randint(state, 30);
        nmod_mat_poly_fit_length(matp, len2);
        _nmod_mat_poly_set_length(matp, len2);
        if (! _check_shape(matp) || ! _check_view(matp))
            TEST_FUNCTION_FAIL("grow again: shape, len1 = %wd, len2 = %wd\n", len1, len2);
        if (! _check_zero_from(matp, len1))
            TEST_FUNCTION_FAIL("grow again: new coefficients not zero\n");
        for (slong k = 0; k < len1; k++)
            for (slong i = 0; i < r; i++)
                for (slong j = 0; j < c; j++)
                    if (nmod_mat_poly_entry(matp, k, i, j) != nmod_mat_poly_entry(copy, k, i, j))
                        TEST_FUNCTION_FAIL("grow again: old coefficients disturbed\n");

        /* shrink then regrow: the coefficients that come back must be zero */
        const slong len3 = n_randint(state, len2 + 1);
        _nmod_mat_poly_set_length(matp, len3);
        _nmod_mat_poly_set_length(matp, len2);
        if (! _check_shape(matp) || ! _check_zero_from(matp, len3))
            TEST_FUNCTION_FAIL("shrink and regrow: len3 = %wd, len2 = %wd\n", len3, len2);

        /* truncate and normalise */
        nmod_mat_poly_rand(matp, state, len2);
        const slong len4 = n_randint(state, len2 + 1);
        nmod_mat_poly_truncate(matp, len4);
        if (matp->length > len4 || ! _check_shape(matp))
            TEST_FUNCTION_FAIL("truncate: len4 = %wd, length = %wd\n", len4, matp->length);
        if (matp->length > 0)
        {
            nmod_mat_t lead;
            nmod_mat_poly_lead_attach(lead, matp);
            if (nmod_mat_is_zero(lead) || nmod_mat_poly_lead_ptr(matp) != lead->entries)
                TEST_FUNCTION_FAIL("truncate: not normalised\n");
        }

        /* realloc down to zero and back up */
        nmod_mat_poly_realloc(matp, 0);
        if (matp->alloc != 0 || matp->length != 0 || matp->coeffs != NULL)
            TEST_FUNCTION_FAIL("realloc to zero\n");
        nmod_mat_poly_fit_length(matp, 5);
        _nmod_mat_poly_set_length(matp, 5);
        if (! _check_shape(matp) || ! _check_zero_from(matp, 0))
            TEST_FUNCTION_FAIL("regrow after realloc to zero\n");

        /* zero and one */
        nmod_mat_poly_zero(matp);
        if (! nmod_mat_poly_is_zero(matp) || matp->length != 0)
            TEST_FUNCTION_FAIL("zero\n");
        if (r == c && r > 0)
        {
            nmod_mat_poly_one(matp);
            if (! nmod_mat_poly_is_one(matp) || ! _check_shape(matp))
                TEST_FUNCTION_FAIL("one\n");
        }

        nmod_mat_poly_clear(copy);
        nmod_mat_poly_clear(matp);
    }

    TEST_FUNCTION_END(state);
}
