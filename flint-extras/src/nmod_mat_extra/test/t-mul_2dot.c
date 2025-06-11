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
#include <flint/nmod.h>
#include <flint/nmod_vec.h>
#include <flint/test_helpers.h>

#include "nmod_mat_extra.h"
 

/*--------------------------------------------------------------*/
/* computes a square matrix product in size len modulo n        */
/*--------------------------------------------------------------*/
int check_nmod_mat_mul_2dot(ulong len, ulong n)
{
    flint_rand_t state;
    nmod_mat_t a, b, c1, c2;
    nmod_t mod;

    flint_rand_init(state);
    nmod_init(&mod, n);

    nmod_mat_init(a, len, len, mod.n);
    nmod_mat_init(b, len, len, mod.n);
    nmod_mat_init(c1, len, len, mod.n);
    nmod_mat_init(c2, len, len, mod.n);

    nmod_mat_rand(a, state);
    nmod_mat_rand(b, state);
    nmod_mat_rand(c1, state);
    nmod_mat_rand(c2, state);
    
    nmod_mat_mul(c1, a, b);
    nmod_mat_mul_2dot(c2, a, b);

    int res = nmod_mat_equal(c1, c2);
    
    nmod_mat_clear(a);
    nmod_mat_clear(b);
    nmod_mat_clear(c1);
    nmod_mat_clear(c2);
    flint_rand_clear(state);

    return res;
}

TEST_FUNCTION_START(nmod_mat_mul_2dot, state)
{
    flint_rand_set_seed(state, time(NULL), time(NULL)+12984125L);

    int i;

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        slong len = n_randint(state, 100) + 1;

        ulong n = n_randint(state, UWORD(1)<<30);

        int res = check_nmod_mat_mul_2dot(len, n);

        if (!res)
            TEST_FUNCTION_FAIL("mul_2dot, len = %wu, n = %wd\n", len, n);
    }

    TEST_FUNCTION_END(state);
}
