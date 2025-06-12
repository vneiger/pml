/*
    Copyright (C) 2025 Vincent Neiger, Éric Schost

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/nmod_vec.h>
#include <flint/test_helpers.h>

#include "nmod_vec_extra.h"


void _nmod_vec_randtest_not_zero(nn_ptr vec, flint_rand_t state, slong len, nmod_t mod)
{
    slong i, sparseness;

    if (n_randint(state, 2))
    {
        for (i = 0; i < len; i++)
            vec[i] = 1 + (n_randtest(state) % (mod.n - 1));
    }
    else
    {
        sparseness = 1 + n_randint(state, FLINT_MAX(2, len));

        for (i = 0; i < len; i++)
        {
            if (n_randint(state, sparseness))
                vec[i] = 1;
            else
                vec[i] = 1 + (n_randtest(state) % (mod.n - 1));
        }
    }
}

void _nmod_vec_inv_naive(nn_ptr res, nn_srcptr vec, ulong len, nmod_t mod)
{
    for (ulong k = 0; k < len; k++)
        res[k] = nmod_inv(vec[k], mod);
}


int check_nmod_vec_inv(slong len, flint_bitcnt_t bit_len, flint_rand_t state)
{
    ulong p;
    nn_ptr vec, res1, res2, res3;
    nmod_t mod;

    ulong res = 0;

    p = n_randprime(state, bit_len, 1);
    nmod_init(&mod, p);

    vec = _nmod_vec_init(len);
    _nmod_vec_randtest_not_zero(vec, state, len, mod);

    res1 = _nmod_vec_init(len);
    _nmod_vec_inv(res1, vec, len, mod);

    res2 = _nmod_vec_init(len);
    _nmod_vec_inv_naive(res2, vec, len, mod);

    res3 = _nmod_vec_init(len);
    _nmod_vec_inv2(res3, vec, len, mod);

    if (!_nmod_vec_equal(res1, res2, len))
        res = 1;

    if (!_nmod_vec_equal(res3, res2, len))
        res = 2;

    _nmod_vec_clear(res1);
    _nmod_vec_clear(res2);
    _nmod_vec_clear(res3);
    _nmod_vec_clear(vec);

    return res;
}

TEST_FUNCTION_START(nmod_vec_inv, state)
{
    int i, result;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        ulong bits = 2 + n_randint(state, 62);  /* FIXME is 64 allowed? */
        ulong len = 1 + n_randint(state, 1000);

        result = check_nmod_vec_inv(len, bits, state);

        if (result)
            TEST_FUNCTION_FAIL("bits = %wu, len = %wu, ret code = %wu\n",
                    bits, len, result);
    }

    TEST_FUNCTION_END(state);
}
