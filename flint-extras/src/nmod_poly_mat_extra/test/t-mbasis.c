/*
    Copyright (C) 2025 Vincent Neiger, Kevin Tran

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/test_helpers.h>

#include "nmod_poly_mat_approximant.h"
#include "nmod_poly_mat_io.h"

#include "testing_collection.h"

/* test one given input for mbasis */
int core_test_mbasis(nmod_poly_mat_t mat, slong order, slong * shift)
{
    int res = 1;

    const slong rdim = mat->r;
    nmod_poly_mat_t appbas;

    slong cshift[rdim];
    for (slong i = 0; i < rdim; i++)
        cshift[i] = shift[i];

    nmod_poly_mat_init(appbas, rdim, rdim, mat->modulus);
    nmod_poly_mat_mbasis(appbas, cshift, mat, order);

    /* testing correctness of nmod_poly_mat_mbasis */
    if (!nmod_poly_mat_is_approximant_basis(appbas, mat, order, shift, ROW_LOWER))
    {
        res = 0;
        if (0)  /* kept for debug */
        {
            printf("nmod_poly_mat_mbasis output is not a minimal approximant basis\n");
            printf("Input matrix:\n");
            nmod_poly_mat_print_pretty(mat, "X");
            printf("Output matrix:\n");
            nmod_poly_mat_print(appbas,"X");
            printf("Residual matrix:\n");
            nmod_poly_mat_t res;
            nmod_poly_mat_init(res, appbas->r, mat->c, mat->modulus);
            nmod_poly_mat_mul(res, appbas, mat);
            nmod_poly_mat_print_pretty(res,"X");
            nmod_poly_mat_clear(res);
            printf("Degree matrix:\n");
            nmod_poly_mat_degree_matrix_print_pretty(appbas);
            printf("Input shift:\t");
            flint_printf("%{slong*}\n", shift, rdim);
            printf("Output shift:\t");
            flint_printf("%{slong*}\n", cshift, rdim);
            printf("\n");
        }
    }

    nmod_poly_mat_clear(appbas);

    return res;
}

/* Test with specified parameters, uniform shift */
int one_test_mbasis(slong prime, slong rdim, slong cdim, slong order, slong len, flint_rand_t state)
{
    int res;

    // random matrix
    nmod_poly_mat_t mat;
    nmod_poly_mat_init(mat, rdim, cdim, prime);
    nmod_poly_mat_randtest(mat, state, len);

    slong * shift = flint_malloc(rdim * sizeof(slong));

    _test_collection_shift_uniform(shift, rdim);
    res = core_test_mbasis(mat, order, shift);

    _test_collection_shift_increasing(shift, rdim);
    res = res && core_test_mbasis(mat, order, shift);

    _test_collection_shift_decreasing(shift, rdim);
    res = res && core_test_mbasis(mat, order, shift);

    _test_collection_shift_shuffle(shift, rdim, state);
    res = res && core_test_mbasis(mat, order, shift);

    _test_collection_shift_hermite(shift, rdim, order);
    res = res && core_test_mbasis(mat, order, shift);

    _test_collection_shift_rhermite(shift, rdim, order);
    res = res && core_test_mbasis(mat, order, shift);

    _test_collection_shift_plateau(shift, rdim, order);
    res = res && core_test_mbasis(mat, order, shift);

    _test_collection_shift_rplateau(shift, rdim, order);
    res = res && core_test_mbasis(mat, order, shift);

    nmod_poly_mat_clear(mat);
    flint_free(shift);

    return res;
}

/* Test against the whole testing collection */
int collection_test_mbasis(slong iter, flint_rand_t state)
{
    /* input matrix for approximation, input shift */
    nmod_poly_mat_t mat;
    slong * shift;

    int verbose = 0;
    if (verbose)
    {
        long total_nb_tests =
        iter
        * 40 /* number of mats (currently 5) x number of shifts (currently 8) */
        * _test_collection_nb_primes
        * _test_collection_nb_dims
        * _test_collection_nb_dims
        * _test_collection_nb_smalldegs;

        printf("Launching testing collection (%ld cases)\n", total_nb_tests);
    }

    for (slong i_primes = 0; i_primes < _test_collection_nb_primes; i_primes++)
        for (slong i_rdims = 0; i_rdims < _test_collection_nb_dims; i_rdims++)
            for (slong i_cdims = 0; i_cdims < _test_collection_nb_dims; i_cdims++)
                for (slong i_degs = 0; i_degs < _test_collection_nb_smalldegs; i_degs++)
                    for (slong it = 0; it < iter; it++)
                    {
                        const long prime = _test_collection_primes[i_primes];
                        const long rdim = _test_collection_rdims[i_rdims];
                        const long cdim = _test_collection_cdims[i_cdims];
                        const long order = _test_collection_smalldegs[i_degs];
                        if (verbose)
                            printf("prime %ld, rdim %ld, cdim %ld, order %ld.\n", prime, rdim, cdim, order);

                        shift = (slong *) flint_malloc(rdim * sizeof(slong));
                        nmod_poly_mat_init(mat, rdim, cdim, prime);

                        _test_collection_shift_uniform(shift, rdim);
                        _test_collection_mat_zero(mat);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "uniform", "zero"); return 0; }

                        _test_collection_shift_uniform(shift, rdim);
                        _test_collection_mat_uniform(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "uniform", "uniform"); return 0; }

                        _test_collection_shift_uniform(shift, rdim);
                        _test_collection_mat_test(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "uniform", "test"); return 0; }

                        _test_collection_shift_uniform(shift, rdim);
                        _test_collection_mat_sparse(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "uniform", "sparse"); return 0; }

                        _test_collection_shift_uniform(shift, rdim);
                        _test_collection_mat_rkdef(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "uniform", "rkdef"); return 0; }


                        _test_collection_shift_increasing(shift, rdim);
                        _test_collection_mat_zero(mat);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "increasing", "zero"); return 0; }

                        _test_collection_shift_increasing(shift, rdim);
                        _test_collection_mat_uniform(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "increasing", "uniform"); return 0; }

                        _test_collection_shift_increasing(shift, rdim);
                        _test_collection_mat_test(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "increasing", "test"); return 0; }

                        _test_collection_shift_increasing(shift, rdim);
                        _test_collection_mat_sparse(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "increasing", "sparse"); return 0; }

                        _test_collection_shift_increasing(shift, rdim);
                        _test_collection_mat_rkdef(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "increasing", "rkdef"); return 0; }


                        _test_collection_shift_decreasing(shift, rdim);
                        _test_collection_mat_zero(mat);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "decreasing", "zero"); return 0; }

                        _test_collection_shift_decreasing(shift, rdim);
                        _test_collection_mat_uniform(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "decreasing", "uniform"); return 0; }

                        _test_collection_shift_decreasing(shift, rdim);
                        _test_collection_mat_test(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "decreasing", "test"); return 0; }

                        _test_collection_shift_decreasing(shift, rdim);
                        _test_collection_mat_sparse(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "decreasing", "sparse"); return 0; }

                        _test_collection_shift_decreasing(shift, rdim);
                        _test_collection_mat_rkdef(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "decreasing", "rkdef"); return 0; }


                        _test_collection_shift_shuffle(shift, rdim, state);
                        _test_collection_mat_zero(mat);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "shuffle", "zero"); return 0; }

                        _test_collection_shift_shuffle(shift, rdim, state);
                        _test_collection_mat_uniform(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "shuffle", "uniform"); return 0; }

                        _test_collection_shift_shuffle(shift, rdim, state);
                        _test_collection_mat_test(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "shuffle", "test"); return 0; }

                        _test_collection_shift_shuffle(shift, rdim, state);
                        _test_collection_mat_sparse(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "shuffle", "sparse"); return 0; }

                        _test_collection_shift_shuffle(shift, rdim, state);
                        _test_collection_mat_rkdef(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "shuffle", "rkdef"); return 0; }


                        _test_collection_shift_hermite(shift, rdim, order);
                        _test_collection_mat_zero(mat);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "hermite", "zero"); return 0; }

                        _test_collection_shift_hermite(shift, rdim, order);
                        _test_collection_mat_uniform(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "hermite", "uniform"); return 0; }

                        _test_collection_shift_hermite(shift, rdim, order);
                        _test_collection_mat_test(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "hermite", "test"); return 0; }

                        _test_collection_shift_hermite(shift, rdim, order);
                        _test_collection_mat_sparse(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "hermite", "sparse"); return 0; }

                        _test_collection_shift_hermite(shift, rdim, order);
                        _test_collection_mat_rkdef(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "hermite", "rkdef"); return 0; }


                        _test_collection_shift_rhermite(shift, rdim, order);
                        _test_collection_mat_zero(mat);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "rhermite", "zero"); return 0; }

                        _test_collection_shift_rhermite(shift, rdim, order);
                        _test_collection_mat_uniform(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "rhermite", "uniform"); return 0; }

                        _test_collection_shift_rhermite(shift, rdim, order);
                        _test_collection_mat_test(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "rhermite", "test"); return 0; }

                        _test_collection_shift_rhermite(shift, rdim, order);
                        _test_collection_mat_sparse(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "rhermite", "sparse"); return 0; }

                        _test_collection_shift_rhermite(shift, rdim, order);
                        _test_collection_mat_rkdef(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "rhermite", "rkdef"); return 0; }


                        _test_collection_shift_plateau(shift, rdim, order);
                        _test_collection_mat_zero(mat);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "plateau", "zero"); return 0; }

                        _test_collection_shift_plateau(shift, rdim, order);
                        _test_collection_mat_uniform(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "plateau", "uniform"); return 0; }

                        _test_collection_shift_plateau(shift, rdim, order);
                        _test_collection_mat_test(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "plateau", "test"); return 0; }

                        _test_collection_shift_plateau(shift, rdim, order);
                        _test_collection_mat_sparse(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "plateau", "sparse"); return 0; }

                        _test_collection_shift_plateau(shift, rdim, order);
                        _test_collection_mat_rkdef(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "plateau", "rkdef"); return 0; }


                        _test_collection_shift_rplateau(shift, rdim, order);
                        _test_collection_mat_zero(mat);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "rplateau", "zero"); return 0; }

                        _test_collection_shift_rplateau(shift, rdim, order);
                        _test_collection_mat_uniform(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "rplateau", "uniform"); return 0; }

                        _test_collection_shift_rplateau(shift, rdim, order);
                        _test_collection_mat_test(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "rplateau", "test"); return 0; }

                        _test_collection_shift_rplateau(shift, rdim, order);
                        _test_collection_mat_sparse(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "rplateau", "sparse"); return 0; }

                        _test_collection_shift_rplateau(shift, rdim, order);
                        _test_collection_mat_rkdef(mat, order-1, state);
                        if (! core_test_mbasis(mat, order, shift))
                        { printf("failed %s -- %s,\n...exiting\n", "rplateau", "rkdef"); return 0; }

                        nmod_poly_mat_clear(mat);
                        flint_free(shift);
                    }

    return 1;
}

TEST_FUNCTION_START(nmod_poly_mat_mbasis, state)
{
    int i, result;

    /* TODO activate this for long test */
    if (0)
    {
        for (i = 0; i < 5 * flint_test_multiplier(); i++)
        {
            {
                result = collection_test_mbasis(5, state);

                if (!result)
                    TEST_FUNCTION_FAIL("Failed mbasis on testing collection\n");
            }
        }
    }
    else
    {
        for (i = 0; i < 100 * flint_test_multiplier(); i++)
        {
            ulong nbits = 2 + n_randint(state, 63);
            ulong rdim = 1 + n_randint(state, 10);
            ulong cdim = 1 + n_randint(state, 10);
            ulong order = n_randint(state, 40);
            ulong len = n_randint(state, 40);

            slong prime = n_randprime(state, nbits, 1);

            result = one_test_mbasis(prime, rdim, cdim, order, len, state);

            if (!result)
                TEST_FUNCTION_FAIL("prime = %wu, rdim = %wu, cdim = %wu, order = %wu, len = %wu\n",
                                   prime, rdim, cdim, order, len);
        }
    }

    TEST_FUNCTION_END(state);
}
