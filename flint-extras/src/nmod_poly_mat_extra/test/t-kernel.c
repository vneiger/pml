/*
    Copyright (C) 2025 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/flint.h>
#include <flint/nmod_vec.h>
#include <flint/ulong_extras.h>
#include <flint/nmod_poly_mat.h>
#include <flint/nmod_types.h>
#include <flint/profiler.h>
#include <flint/test_helpers.h>

#include "nmod_poly_mat_forms.h"
#include "nmod_poly_mat_kernel.h"

TEST_FUNCTION_START(nmod_poly_mat_kernel, state)
{
    int i, result;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        ulong nbits = 2 + n_randint(state, 63);
        slong rdim = n_randint(state, 12);
        slong cdim = n_randint(state, 12);
        ulong len = 1 + n_randint(state, 150);

        slong * shift = FLINT_ARRAY_ALLOC(rdim, slong);
        for (slong i = 0; i < rdim; i++)
            shift[i] = n_randint(state, len);

        slong * pivind = FLINT_ARRAY_ALLOC(rdim, slong);
        slong * pivind_check = FLINT_ARRAY_ALLOC(rdim, slong);
        slong * rdeg = FLINT_ARRAY_ALLOC(rdim, slong);
        slong * rdeg_check = FLINT_ARRAY_ALLOC(rdim, slong);

        ulong prime = n_randprime(state, nbits, 1);

        nmod_poly_mat_t pmat;
        nmod_poly_mat_init(pmat, rdim, cdim, prime);
        nmod_poly_mat_randtest(pmat, state, len);

        nmod_poly_mat_t ker;
        nmod_poly_mat_t kernz;
        nmod_poly_mat_init(ker, rdim, rdim, prime);

        {
            for (slong i = 0; i < rdim; i++)
                rdeg[i] = shift[i];
            slong nz = nmod_poly_mat_kernel_via_approx(ker, pivind, rdeg, pmat);

            nmod_poly_mat_window_init(kernz, ker, 0, 0, nz, rdim);
            result = nmod_poly_mat_is_kernel(kernz, shift, pmat, ORD_WEAK_POPOV, ROW_LOWER);

            if (!result)
            {
                TEST_FUNCTION_FAIL("(via_approx, ker) -- rdim = %wu, cdim = %wu, length = %wu, p = %wu\n", \
                                   rdim, cdim, len, prime);
            }

            nmod_poly_mat_row_degree(rdeg_check, kernz, shift);
            nmod_poly_mat_pivot_index(pivind_check, kernz, shift, ROW_LOWER);
            result = 1;
            for (slong k = 0; k < nz; k++)
                if (rdeg_check[k] != rdeg[k] || pivind_check[k] != pivind[k])
                    result = 0;

            if (!result)
            {
                TEST_FUNCTION_FAIL("(via_approx, rdeg/pivind) -- rdim = %wu, cdim = %wu, length = %wu, p = %wu\n", \
                                   rdim, cdim, len, prime);
            }
            nmod_poly_mat_window_clear(kernz);
        }

        {
            for (slong i = 0; i < rdim; i++)
                rdeg[i] = shift[i];
            nmod_poly_mat_t pmat_c;
            nmod_poly_mat_init_set(pmat_c, pmat);
            slong nz = nmod_poly_mat_kernel_zls_approx(ker, pivind, rdeg, pmat_c);

            nmod_poly_mat_window_init(kernz, ker, 0, 0, nz, rdim);
            result = nmod_poly_mat_is_kernel(kernz, shift, pmat, ORD_WEAK_POPOV, ROW_LOWER);

            if (!result)
            {
                TEST_FUNCTION_FAIL("(zls_approx, ker) -- rdim = %wu, cdim = %wu, length = %wu, p = %wu\n", \
                                   rdim, cdim, len, prime);
            }

            nmod_poly_mat_row_degree(rdeg_check, kernz, shift);
            nmod_poly_mat_pivot_index(pivind_check, kernz, shift, ROW_LOWER);
            result = 1;
            for (slong k = 0; k < nz; k++)
                if (rdeg_check[k] != rdeg[k] || pivind_check[k] != pivind[k])
                    result = 0;

            if (!result)
            {
                TEST_FUNCTION_FAIL("(zls_approx, rdeg/pivind) -- rdim = %wu, cdim = %wu, length = %wu, p = %wu\n", \
                                   rdim, cdim, len, prime);
            }

            nmod_poly_mat_clear(pmat_c);
            nmod_poly_mat_window_clear(kernz);
        }

        nmod_poly_mat_window_clear(kernz);
        nmod_poly_mat_clear(pmat);
        nmod_poly_mat_clear(ker);
        flint_free(shift);
        flint_free(rdeg);
        flint_free(pivind);
        flint_free(rdeg_check);
        flint_free(pivind_check);
    }

    TEST_FUNCTION_END(state);
}
