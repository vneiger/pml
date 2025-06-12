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


/*--------------------------------------------------------------*/
/* computes an integer dot-product, length len, bitsize bit_len */
/* checks modulo a prime                                        */
/*--------------------------------------------------------------*/
int check_nmod_vec_integer_dot_product(slong len, flint_bitcnt_t bit_len, flint_rand_t state)
{
    ulong p, res1, res2;
    nn_ptr v1, v2, res;
    nmod_t mod;

    p = n_randprime(state, bit_len, 0);
    nmod_init(&mod, p);

    v1 = _nmod_vec_init(len);
    v2 = _nmod_vec_init(len);

    _nmod_vec_randtest(v1, state, len, mod);
    _nmod_vec_randtest(v2, state, len, mod);

    res = _nmod_vec_init(3);
    nmod_vec_integer_dot_product(res, v1, v2, len, bit_len, bit_len);
    NMOD_RED3(res1, res[2], res[1], res[0], mod);
    dot_params_t params = {_DOT3, UWORD(0)};
    res2 = _nmod_vec_dot(v1, v2, len, mod, params);

    int result = (res1 == res2);

    _nmod_vec_clear(res);
    _nmod_vec_clear(v1);
    _nmod_vec_clear(v2);

    return result;
}

TEST_FUNCTION_START(nmod_vec_integer_dot_product, state)
{
    int i, result;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        ulong bits = 2 + n_randint(state, 63);
        ulong len = 1 + n_randint(state, 2000);

        result = check_nmod_vec_integer_dot_product(len, bits,state);

        if (!result)
            TEST_FUNCTION_FAIL("bits = %wu, len = %wu\n", bits, len);
    }

    TEST_FUNCTION_END(state);
}
