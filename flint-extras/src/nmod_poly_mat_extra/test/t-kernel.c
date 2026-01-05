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
    int k, result;

    /* left kernel */
    for (k = 0; k < 50 * flint_test_multiplier(); k++)
    {
        ulong nbits = 2 + n_randint(state, 63);
        slong cdim = n_randint(state, 12);
        slong rdim = n_randint(state, 12);
        ulong len = 1 + n_randint(state, 75);

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

        /* kernel_via_approx */
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
            for (slong i = 0; i < nz; i++)
                if (rdeg_check[i] != rdeg[i] || pivind_check[i] != pivind[i])
                    result = 0;

            if (!result)
            {
                TEST_FUNCTION_FAIL("(via_approx, rdeg/pivind) -- rdim = %wu, cdim = %wu, length = %wu, p = %wu\n", \
                                   rdim, cdim, len, prime);
            }
            nmod_poly_mat_window_clear(kernz);
        }

        /* kernel_zls_approx */
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
            for (slong i = 0; i < nz; i++)
                if (rdeg_check[i] != rdeg[i] || pivind_check[i] != pivind[i])
                    result = 0;

            if (!result)
            {
                TEST_FUNCTION_FAIL("(zls_approx, rdeg/pivind) -- rdim = %wu, cdim = %wu, length = %wu, p = %wu\n", \
                                   rdim, cdim, len, prime);
            }

            nmod_poly_mat_clear(pmat_c);
            nmod_poly_mat_window_clear(kernz);
        }

        /* kernel interface */
        {
            poly_mat_form_t form = ORD_WEAK_POPOV;
            orientation_t orient = (n_randint(state, 2)) ? ROW_LOWER : ROW_UPPER;

            for (slong i = 0; i < rdim; i++)
                rdeg[i] = shift[i];
            slong nz = nmod_poly_mat_kernel(ker, pivind, rdeg, pmat, form, orient);

            nmod_poly_mat_window_init(kernz, ker, 0, 0, nz, rdim);
            result = nmod_poly_mat_is_kernel(kernz, shift, pmat, form, orient);

            if (!result)
            {
                TEST_FUNCTION_FAIL("(kernel, ker) -- rdim = %wu, cdim = %wu, length = %wu, p = %wu, orient = %wu\n", \
                                   rdim, cdim, len, prime, orient);
            }

            nmod_poly_mat_row_degree(rdeg_check, kernz, shift);
            nmod_poly_mat_pivot_index(pivind_check, kernz, shift, orient);
            result = 1;
            for (slong i = 0; i < nz; i++)
                if (rdeg_check[i] != rdeg[i] || pivind_check[i] != pivind[i])
                    result = 0;

            if (!result)
            {
                flint_printf("true rdeg == %{slong*}\n", rdeg, nz);
                flint_printf("got  rdeg == %{slong*}\n", rdeg_check, nz);
                flint_printf("true pivind == %{slong*}\n", pivind, nz);
                flint_printf("got  pivind == %{slong*}\n", pivind_check, nz);
                TEST_FUNCTION_FAIL("(kernel, rdeg/pivind) -- rdim = %wu, cdim = %wu, length = %wu, p = %wu, orient = %wu\n", \
                                   rdim, cdim, len, prime, orient);
            }

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

    /* right kernel */
    for (k = 0; k < 50 * flint_test_multiplier(); k++)
    {
        ulong nbits = 2 + n_randint(state, 63);
        slong cdim = n_randint(state, 12);
        slong rdim = n_randint(state, 12);
        ulong len = 1 + n_randint(state, 75);

        slong * shift = FLINT_ARRAY_ALLOC(cdim, slong);
        for (slong j = 0; j < cdim; j++)
            shift[j] = n_randint(state, len);

        slong * pivind = FLINT_ARRAY_ALLOC(cdim, slong);
        slong * pivind_check = FLINT_ARRAY_ALLOC(cdim, slong);
        slong * cdeg = FLINT_ARRAY_ALLOC(cdim, slong);
        slong * cdeg_check = FLINT_ARRAY_ALLOC(cdim, slong);

        ulong prime = n_randprime(state, nbits, 1);

        nmod_poly_mat_t pmat;
        nmod_poly_mat_init(pmat, rdim, cdim, prime);
        nmod_poly_mat_randtest(pmat, state, len);

        nmod_poly_mat_t ker;
        nmod_poly_mat_t kernz;
        nmod_poly_mat_init(ker, cdim, cdim, prime);

        /* kernel interface */
        {
            poly_mat_form_t form = ORD_WEAK_POPOV;
            /* orientation_t orient = (n_randint(state, 2)) ? COL_LOWER : COL_UPPER; */
            orientation_t orient = COL_UPPER;

            for (slong j = 0; j < cdim; j++)
                cdeg[j] = shift[j];
            slong nz = nmod_poly_mat_kernel(ker, pivind, cdeg, pmat, form, orient);

            nmod_poly_mat_window_init(kernz, ker, 0, 0, cdim, nz);
            result = nmod_poly_mat_is_kernel(kernz, shift, pmat, form, orient);

            if (!result)
            {
                TEST_FUNCTION_FAIL("(kernel, ker) -- rdim = %wu, cdim = %wu, length = %wu, p = %wu, orient = %wu\n", \
                                   rdim, cdim, len, prime, orient);
            }

            nmod_poly_mat_column_degree(cdeg_check, kernz, shift);
            nmod_poly_mat_pivot_index(pivind_check, kernz, shift, orient);
            result = 1;
            for (slong k = 0; k < nz; k++)
                if (cdeg_check[k] != cdeg[k] || pivind_check[k] != pivind[k])
                    result = 0;

            if (!result)
            {
                flint_printf("true cdeg == %{slong*}\n", cdeg, nz);
                flint_printf("got  cdeg == %{slong*}\n", cdeg_check, nz);
                flint_printf("true pivind == %{slong*}\n", pivind, nz);
                flint_printf("got  pivind == %{slong*}\n", pivind_check, nz);
                TEST_FUNCTION_FAIL("(kernel, rdeg/pivind) -- rdim = %wu, cdim = %wu, length = %wu, p = %wu, orient = %wu\n", \
                                   rdim, cdim, len, prime, orient);
            }

            nmod_poly_mat_window_clear(kernz);
        }

        nmod_poly_mat_window_clear(kernz);
        nmod_poly_mat_clear(pmat);
        nmod_poly_mat_clear(ker);
        flint_free(shift);
        flint_free(cdeg);
        flint_free(pivind);
        flint_free(cdeg_check);
        flint_free(pivind_check);
    }

    TEST_FUNCTION_END(state);
}
