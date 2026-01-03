/*
    Copyright (C) 2025 Gilles Villard

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
#include <flint/nmod_types.h>
#include <flint/profiler.h>
#include <flint/test_helpers.h>

#include "nmod_poly_mat_kernel.h"

TEST_FUNCTION_START(nmod_poly_mat_kernel_via_approx, state)
{
    int i, result;

    for (i = 0; i < 16 * flint_test_multiplier(); i++)
    {
        ulong nbits = 2 + n_randint(state, 63);
        slong rdim = n_randint(state, 30);
        slong cdim = n_randint(state, 30);
        ulong len = 1 + n_randint(state, 100);
        flint_printf("TEST iter %ld -- %ld, %ld, %ld\n", i, rdim, cdim, len);

        slong * shift = FLINT_ARRAY_ALLOC(rdim, slong);
        for (slong i = 0; i < rdim; i++)
            shift[i] = 0;

        slong * pivind = FLINT_ARRAY_ALLOC(rdim, slong);
        slong * rdeg = FLINT_ARRAY_ALLOC(rdim, slong);
        for (slong i = 0; i < rdim; i++)
            rdeg[i] = shift[i];

        ulong prime = n_randprime(state, nbits, 1);

        nmod_poly_mat_t pmat;
        nmod_poly_mat_init(pmat, rdim, cdim, prime);
        nmod_poly_mat_randtest(pmat, state, len);

        nmod_poly_mat_t ker;
        nmod_poly_mat_init(ker, rdim, rdim, prime);
        /* TODO test weak Popov + pivind */
        /* TODO currently does not test generation */
        slong nz = nmod_poly_mat_kernel_via_approx(ker, pivind, rdeg, pmat);
        flint_printf("kernel computed, now testing...\n");
        result = nmod_poly_mat_is_kernel(ker, nz, shift, pmat, ROW_LOWER);

        nmod_poly_mat_clear(pmat);
        flint_free(shift);
        flint_free(rdeg);
        flint_free(pivind);

        if (!result)
        {
            TEST_FUNCTION_FAIL("rdim = %wu, cdim = %wu, length = %wu, p = %wu\n", \
                rdim, cdim, len, prime);
        }
    }

    TEST_FUNCTION_END(state);
}
