/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/nmod_poly_mat.h>
#include <flint/test_helpers.h>
#include <flint/ulong_extras.h>

#include "nmod_poly_mat_utils.h"

/* reference: c = a + x^k b (sub == 0) or a - x^k b (sub == 1), via a left
 * shift followed by an addition or subtraction */
static void _addshifttest_ref(nmod_poly_mat_t c, const nmod_poly_mat_t a,
                              const nmod_poly_mat_t b, slong k, int sub)
{
    nmod_poly_mat_t t;
    nmod_poly_mat_init(t, b->r, b->c, b->modulus);
    nmod_poly_mat_shift_left(t, b, k);
    if (sub)
        nmod_poly_mat_sub(c, a, t);
    else
        nmod_poly_mat_add(c, a, t);
    nmod_poly_mat_clear(t);
}

static void _addshifttest_call(nmod_poly_mat_t c, const nmod_poly_mat_t a,
                               const nmod_poly_mat_t b, slong k, int sub)
{
    if (sub)
        nmod_poly_mat_sub_shift_left(c, a, b, k);
    else
        nmod_poly_mat_add_shift_left(c, a, b, k);
}

/* checks that all entries are normalized (a nonzero leading coefficient) */
static int _addshifttest_is_normalised(const nmod_poly_mat_t c)
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

TEST_FUNCTION_START(nmod_poly_mat_add_shift_left, state)
{
    int i, result;

    for (i = 0; i < 2000 * flint_test_multiplier(); i++)
    {
        const ulong prime = n_randprime(state, 2 + n_randint(state, 63), 1);
        const slong r = n_randint(state, 4);
        const slong s = n_randint(state, 4);
        const slong k = n_randint(state, 20);
        /* lengths of a: below, around, or above k (+ length of b) */
        const slong lena = n_randint(state, 2*k + 20);
        const slong lenb = n_randint(state, 20);
        const int sub = n_randint(state, 2);

        nmod_poly_mat_t a, b, c, d;
        nmod_poly_mat_init(a, r, s, prime);
        nmod_poly_mat_init(b, r, s, prime);
        nmod_poly_mat_init(c, r, s, prime);
        nmod_poly_mat_init(d, r, s, prime);

        nmod_poly_mat_randtest(a, state, lena);
        nmod_poly_mat_randtest(b, state, lenb);
        if (n_randint(state, 4) == 0)
        {
            /* b = +-(a div x^k), so that x^k b cancels the high part of a;
             * sometimes the constant coefficient of one entry of b is
             * perturbed, so that the cancellation is only partial */
            nmod_poly_mat_shift_right(b, a, k);
            if (! sub)
                nmod_poly_mat_neg(b, b);
            if (n_randint(state, 2) && r > 0 && s > 0)
            {
                nmod_poly_struct * p = nmod_poly_mat_entry(b, 0, 0);
                if (p->length > 1)
                    nmod_poly_set_coeff_ui(p, 0, n_randint(state, prime));
            }
        }
        nmod_poly_mat_randtest(c, state, 5);

        _addshifttest_ref(d, a, b, k, sub);

        /* no aliasing */
        _addshifttest_call(c, a, b, k, sub);
        result = nmod_poly_mat_equal(c, d) && _addshifttest_is_normalised(c);

        /* c aliases a */
        if (result)
        {
            nmod_poly_mat_set(c, a);
            _addshifttest_call(c, c, b, k, sub);
            result = nmod_poly_mat_equal(c, d) && _addshifttest_is_normalised(c);
        }

        /* c aliases b */
        if (result)
        {
            nmod_poly_mat_set(c, b);
            _addshifttest_call(c, a, c, k, sub);
            result = nmod_poly_mat_equal(c, d) && _addshifttest_is_normalised(c);
        }

        /* c aliases a and b */
        if (result)
        {
            _addshifttest_ref(d, a, a, k, sub);
            nmod_poly_mat_set(c, a);
            _addshifttest_call(c, c, c, k, sub);
            result = nmod_poly_mat_equal(c, d) && _addshifttest_is_normalised(c);
        }

        /* a aliases b, c does not */
        if (result)
        {
            _addshifttest_call(c, a, a, k, sub);
            result = nmod_poly_mat_equal(c, d) && _addshifttest_is_normalised(c);
        }

        if (!result)
            TEST_FUNCTION_FAIL("prime = %wu, r = %wd, s = %wd, k = %wd, lena = %wd, lenb = %wd, sub = %d\n",
                               prime, r, s, k, lena, lenb, sub);

        nmod_poly_mat_clear(a);
        nmod_poly_mat_clear(b);
        nmod_poly_mat_clear(c);
        nmod_poly_mat_clear(d);
    }

    TEST_FUNCTION_END(state);
}
