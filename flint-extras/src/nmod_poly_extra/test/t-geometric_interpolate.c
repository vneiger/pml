/*
    Copyright (C) 2025 Vincent Neiger, Éric Schost

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/test_helpers.h>

#include "nmod_extra.h"
#include "nmod_poly_extra.h"

TEST_FUNCTION_START(nmod_poly_geometric_interpolate, state)
{
    int i, result;

    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        slong ii;
        ulong nn;
        nmod_t mod;

        /* no find-root test yet: must be sufficiently large */
        ulong bits = 16 + n_randint(state, 48);
        nn = n_randprime(state, bits, 1);
        nmod_init(&mod, nn);

        for (ii = 100; ii < 5000; ii+=500)
        {
            ulong r;
            nmod_geometric_progression_t G;
            nmod_poly_t P, P2;
            nn_ptr val;

            r = nmod_find_root(2*ii, mod);
            nmod_geometric_progression_init_set(G, r, ii, mod);

            nmod_poly_init(P, nn);
            nmod_poly_init(P2, nn);
            nmod_poly_randtest(P, state, ii);

            val = _nmod_vec_init(ii);
            nmod_geometric_progression_evaluate(val, P, G);
            nmod_geometric_progression_interpolate(P2, val, G);

            result = (nmod_poly_equal(P, P2));

            _nmod_vec_clear(val);
            nmod_poly_clear(P2);
            nmod_poly_clear(P);
            nmod_geometric_progression_clear(G);

            if (!result)
                TEST_FUNCTION_FAIL("mod = %wu, ii = %wu, r = %wu\n", nn, ii, r);
        }
    }

    TEST_FUNCTION_END(state);
}
