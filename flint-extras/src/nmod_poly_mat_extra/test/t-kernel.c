/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/flint.h>
#include <flint/nmod_poly_mat.h>
#include <flint/profiler.h>
#include <flint/test_helpers.h>

#include "nmod_poly_mat_forms.h"
#include "nmod_poly_mat_approximant.h"
#include "nmod_poly_mat_kernel.h"
#include "testing_collection.h"

#define _TEST_VERBOSE_ 0

/* test one given input for kernel */
/* note: length of shift should be at least max(mat->r, mat->c) */
int core_test_kernel(nmod_poly_mat_t mat, slong * shift)
{
    int res = 1;
    int local_res = 1;

    const slong prime = mat->modulus;
    const slong rdim = mat->r;
    const slong cdim = mat->c;

    /* left kernel */
    {
        nmod_poly_mat_t ker;
        nmod_poly_mat_init(ker, rdim, rdim, prime);

        slong * pivind = FLINT_ARRAY_ALLOC(rdim, slong);
        slong * pivind_check = FLINT_ARRAY_ALLOC(rdim, slong);
        slong * rdeg = FLINT_ARRAY_ALLOC(rdim, slong);
        slong * rdeg_check = FLINT_ARRAY_ALLOC(rdim, slong);

        nmod_poly_mat_t kernz;

        /* kernel_via_approx */
        {
            if (_TEST_VERBOSE_) flint_printf(". starting left kernel, via approx...\n");
            for (slong i = 0; i < rdim; i++)
                rdeg[i] = shift[i];
            slong nullity = nmod_poly_mat_kernel_via_approx(ker, pivind, rdeg, mat);

            nmod_poly_mat_window_init(kernz, ker, 0, 0, nullity, rdim);
            local_res = nmod_poly_mat_is_kernel(kernz, shift, mat, ORD_WEAK_POPOV, ROW_LOWER);
            res = res && local_res;

            if (!local_res)
                flint_printf("(via_approx, ker) -- rdim = %wu, cdim = %wu, length = %wu, p = %wu\n", \
                             rdim, mat->c, nmod_poly_mat_degree(mat), prime);

            nmod_poly_mat_row_degree(rdeg_check, kernz, shift);
            nmod_poly_mat_pivot_index(pivind_check, kernz, shift, ROW_LOWER);
            local_res = 1;
            for (slong i = 0; i < nullity; i++)
                if (rdeg_check[i] != rdeg[i] || pivind_check[i] != pivind[i])
                    local_res = 0;
            res = res && local_res;

            if (!local_res)
                flint_printf("(via_approx, rdeg/pivind) -- rdim = %wu, cdim = %wu, length = %wu, p = %wu\n", \
                             rdim, cdim, nmod_poly_mat_degree(mat), prime);
            nmod_poly_mat_window_clear(kernz);
        }

        /* kernel_zls_approx */
        {
            if (_TEST_VERBOSE_) flint_printf(". starting left kernel, ZLS-approx...\n");
            for (slong i = 0; i < rdim; i++)
                rdeg[i] = shift[i];
            nmod_poly_mat_t pmat_c;
            nmod_poly_mat_init_set(pmat_c, mat);
            slong nullity = nmod_poly_mat_kernel_zls_approx(ker, pivind, rdeg, pmat_c);

            nmod_poly_mat_window_init(kernz, ker, 0, 0, nullity, rdim);
            local_res = nmod_poly_mat_is_kernel(kernz, shift, mat, ORD_WEAK_POPOV, ROW_LOWER);
            res = res && local_res;

            if (!local_res)
                flint_printf("(zls_approx, ker) -- rdim = %wu, cdim = %wu, length = %wu, p = %wu\n", \
                             rdim, cdim, nmod_poly_mat_degree(mat), prime);

            nmod_poly_mat_row_degree(rdeg_check, kernz, shift);
            nmod_poly_mat_pivot_index(pivind_check, kernz, shift, ROW_LOWER);
            local_res = 1;
            for (slong i = 0; i < nullity; i++)
                if (rdeg_check[i] != rdeg[i] || pivind_check[i] != pivind[i])
                    local_res = 0;
            res = res && local_res;

            if (!local_res)
                flint_printf("(zls_approx, rdeg/pivind) -- rdim = %wu, cdim = %wu, length = %wu, p = %wu\n", \
                             rdim, cdim, nmod_poly_mat_degree(mat), prime);

            nmod_poly_mat_clear(pmat_c);
            nmod_poly_mat_window_clear(kernz);
        }

        /* kernel interface, ROW_LOWER */
        if (0)  /* TODO currently disabled: same as kernel_zls */
        {
            if (_TEST_VERBOSE_) flint_printf(". starting left kernel, interface ROW_LOWER...\n");
            poly_mat_form_t form = ORD_WEAK_POPOV;
            orientation_t orient = ROW_LOWER;

            for (slong i = 0; i < rdim; i++)
                rdeg[i] = shift[i];
            slong nullity = nmod_poly_mat_kernel(ker, pivind, rdeg, mat, form, orient);

            nmod_poly_mat_window_init(kernz, ker, 0, 0, nullity, rdim);
            local_res = nmod_poly_mat_is_kernel(kernz, shift, mat, form, orient);
            res = res && local_res;

            if (!local_res)
                flint_printf("(lkernel, row_lower, ker) -- rdim = %wu, cdim = %wu, length = %wu, p = %wu, orient = %wu\n", \
                             rdim, cdim, nmod_poly_mat_degree(mat), prime, orient);

            nmod_poly_mat_row_degree(rdeg_check, kernz, shift);
            nmod_poly_mat_pivot_index(pivind_check, kernz, shift, orient);
            local_res = 1;
            for (slong i = 0; i < nullity; i++)
                if (rdeg_check[i] != rdeg[i] || pivind_check[i] != pivind[i])
                    local_res = 0;
            res = res && local_res;

            if (!local_res)
                flint_printf("(lkernel, row_lower, rdeg/pivind) -- rdim = %wu, cdim = %wu, length = %wu, p = %wu, orient = %wu\n", \
                             rdim, cdim, nmod_poly_mat_degree(mat), prime, orient);

            nmod_poly_mat_window_clear(kernz);
        }

        /* kernel interface, ROW_UPPER */
        {
            if (_TEST_VERBOSE_) flint_printf(". starting left kernel, interface ROW_UPPER...\n");
            poly_mat_form_t form = ORD_WEAK_POPOV;
            orientation_t orient = ROW_UPPER;

            for (slong i = 0; i < rdim; i++)
                rdeg[i] = shift[i];
            slong nullity = nmod_poly_mat_kernel(ker, pivind, rdeg, mat, form, orient);

            nmod_poly_mat_window_init(kernz, ker, 0, 0, nullity, rdim);
            local_res = nmod_poly_mat_is_kernel(kernz, shift, mat, form, orient);
            res = res && local_res;

            if (!local_res)
                flint_printf("(lkernel, row_upper, ker) -- rdim = %wu, cdim = %wu, length = %wu, p = %wu, orient = %wu\n", \
                             rdim, cdim, nmod_poly_mat_degree(mat), prime, orient);

            nmod_poly_mat_row_degree(rdeg_check, kernz, shift);
            nmod_poly_mat_pivot_index(pivind_check, kernz, shift, orient);
            local_res = 1;
            for (slong i = 0; i < nullity; i++)
                if (rdeg_check[i] != rdeg[i] || pivind_check[i] != pivind[i])
                    local_res = 0;
            res = res && local_res;

            if (!local_res)
                flint_printf("(lkernel, row_upper, rdeg/pivind) -- rdim = %wu, cdim = %wu, length = %wu, p = %wu, orient = %wu\n", \
                             rdim, cdim, nmod_poly_mat_degree(mat), prime, orient);

            nmod_poly_mat_window_clear(kernz);
        }

        nmod_poly_mat_window_clear(kernz);
        nmod_poly_mat_clear(ker);
        flint_free(rdeg);
        flint_free(pivind);
        flint_free(rdeg_check);
        flint_free(pivind_check);
    }

    /* right kernel */
    {
        nmod_poly_mat_t ker;
        nmod_poly_mat_init(ker, cdim, cdim, prime);

        slong * pivind = FLINT_ARRAY_ALLOC(cdim, slong);
        slong * pivind_check = FLINT_ARRAY_ALLOC(cdim, slong);
        slong * cdeg = FLINT_ARRAY_ALLOC(cdim, slong);
        slong * cdeg_check = FLINT_ARRAY_ALLOC(cdim, slong);

        nmod_poly_mat_t kernz;

        /* kernel interface, COL_LOWER */
        {
            if (_TEST_VERBOSE_) flint_printf(". starting right kernel, interface COL_LOWER...\n");
            poly_mat_form_t form = ORD_WEAK_POPOV;
            orientation_t orient = COL_LOWER;

            for (slong j = 0; j < cdim; j++)
                cdeg[j] = shift[j];
            slong nullity = nmod_poly_mat_kernel(ker, pivind, cdeg, mat, form, orient);

            nmod_poly_mat_window_init(kernz, ker, 0, 0, cdim, nullity);
            local_res = nmod_poly_mat_is_kernel(kernz, shift, mat, form, orient);
            res = res && local_res;

            if (!local_res)
                flint_printf("(rkernel, col_lower, ker) -- rdim = %wu, cdim = %wu, length = %wu, p = %wu, orient = %wu\n", \
                             rdim, cdim, nmod_poly_mat_degree(mat), prime, orient);

            nmod_poly_mat_column_degree(cdeg_check, kernz, shift);
            nmod_poly_mat_pivot_index(pivind_check, kernz, shift, orient);
            local_res = 1;
            for (slong k = 0; k < nullity; k++)
                if (cdeg_check[k] != cdeg[k] || pivind_check[k] != pivind[k])
                    local_res = 0;
            res = res && local_res;

            if (!local_res)
                flint_printf("(rkernel, col_lower, cdeg/pivind) -- rdim = %wu, cdim = %wu, length = %wu, p = %wu, orient = %wu\n", \
                             rdim, cdim, nmod_poly_mat_degree(mat), prime, orient);

            nmod_poly_mat_window_clear(kernz);
        }

        /* kernel interface, COL_UPPER */
        {
            if (_TEST_VERBOSE_) flint_printf(". starting right kernel, interface COL_UPPER...\n");
            poly_mat_form_t form = ORD_WEAK_POPOV;
            orientation_t orient = COL_UPPER;

            for (slong j = 0; j < cdim; j++)
                cdeg[j] = shift[j];
            slong nullity = nmod_poly_mat_kernel(ker, pivind, cdeg, mat, form, orient);

            nmod_poly_mat_window_init(kernz, ker, 0, 0, cdim, nullity);
            local_res = nmod_poly_mat_is_kernel(kernz, shift, mat, form, orient);
            res = res && local_res;

            if (!local_res)
                flint_printf("(rkernel, col_upper, ker) -- rdim = %wu, cdim = %wu, length = %wu, p = %wu, orient = %wu\n", \
                             rdim, cdim, nmod_poly_mat_degree(mat), prime, orient);

            nmod_poly_mat_column_degree(cdeg_check, kernz, shift);
            nmod_poly_mat_pivot_index(pivind_check, kernz, shift, orient);
            local_res = 1;
            for (slong k = 0; k < nullity; k++)
                if (cdeg_check[k] != cdeg[k] || pivind_check[k] != pivind[k])
                    local_res = 0;
            res = res && local_res;

            if (!local_res)
                flint_printf("(rkernel, col_upper, cdeg/pivind) -- rdim = %wu, cdim = %wu, length = %wu, p = %wu, orient = %wu\n", \
                             rdim, cdim, nmod_poly_mat_degree(mat), prime, orient);

            nmod_poly_mat_window_clear(kernz);
        }

        nmod_poly_mat_window_clear(kernz);
        nmod_poly_mat_clear(ker);
        flint_free(cdeg);
        flint_free(pivind);
        flint_free(cdeg_check);
        flint_free(pivind_check);
    }

    return res;
}

/** Test with specified parameters, uniform shift */
int one_test_kernel(slong prime, slong rdim, slong cdim, slong len, flint_rand_t state)
{
    int res;

    nmod_poly_mat_t mat;
    nmod_poly_mat_init(mat, rdim, cdim, prime);
    nmod_poly_mat_randtest(mat, state, len);

    const slong maxdim = FLINT_MAX(rdim, cdim);
    slong * shift = flint_malloc(maxdim * sizeof(slong));

    if (_TEST_VERBOSE_) flint_printf("shift: uniform\n");
    _test_collection_shift_uniform(shift, maxdim);
    res = core_test_kernel(mat, shift);

    if (_TEST_VERBOSE_) flint_printf("shift: increasing\n");
    _test_collection_shift_increasing(shift, maxdim);
    res = res && core_test_kernel(mat, shift);

    if (_TEST_VERBOSE_) flint_printf("shift: decreasing\n");
    _test_collection_shift_decreasing(shift, maxdim);
    res = res && core_test_kernel(mat, shift);

    if (_TEST_VERBOSE_) flint_printf("shift: shuffled\n");
    _test_collection_shift_shuffle(shift, maxdim, state);
    res = res && core_test_kernel(mat, shift);

    /** disabling tests that take too long currently
    const slong degdet = FLINT_MAX(rdim, cdim) * len;

    _test_collection_shift_hermite(shift, rdim, degdet);
    res = res && core_test_kernel(mat, shift);

    _test_collection_shift_rhermite(shift, rdim, degdet);
    res = res && core_test_kernel(mat, shift);

    _test_collection_shift_plateau(shift, rdim, degdet);
    res = res && core_test_kernel(mat, shift);

    _test_collection_shift_rplateau(shift, rdim, degdet);
    res = res && core_test_kernel(mat, shift);
    */

    nmod_poly_mat_clear(mat);
    flint_free(shift);

    return res;
}

/** Test against the whole testing collection */
//int collection_test_pmbasis(slong iter, flint_rand_t state)
//{
//    // input matrix for approximation
//    nmod_poly_mat_t mat;
//
//    // input shift
//    slong * shift;
//
//    long total_nb_tests =
//            iter // number of iterations
//            * 40 // number of mats (currently 5) x number of shifts (currently 8)
//            * _test_collection_nb_primes
//            * _test_collection_nb_dims
//            * _test_collection_nb_dims
//            * _test_collection_nb_degs;
//    printf("Launching testing collection (%ld cases)\n", total_nb_tests);
//
//    for (slong i_primes = 0; i_primes < _test_collection_nb_primes; i_primes++)
//        for (slong i_rdims = 0; i_rdims < _test_collection_nb_dims; i_rdims++)
//            for (slong i_cdims = 0; i_cdims < _test_collection_nb_dims; i_cdims++)
//                for (slong i_degs = 0; i_degs < _test_collection_nb_degs; i_degs++)
//                    for (slong it = 0; it < iter; it++)
//                    {
//                        const long prime = _test_collection_primes[i_primes];
//                        const long rdim = _test_collection_rdims[i_rdims];
//                        const long cdim = _test_collection_cdims[i_cdims];
//                        const long order = _test_collection_degs[i_degs];
//                        printf("prime %ld, rdim %ld, cdim %ld, order %ld.\n", prime, rdim, cdim, order);
//
//                        shift = (slong *) flint_malloc(rdim * sizeof(slong));
//                        nmod_poly_mat_init(mat, rdim, cdim, prime);
//
//                        /* uniform shift */
//                        /* matrix: zero | uniform | rdeg | cdeg | test | sparse | rkdef */
//                        _test_collection_shift_uniform(shift, rdim);
//                        _test_collection_mat_zero(mat);
//                        if (! core_test_pmbasis(mat, order, shift))
//                        { printf("failed %s -- %s,\n...exiting\n", "uniform", "zero"); return 0; }
//
//                        _test_collection_shift_uniform(shift, rdim);
//                        _test_collection_mat_uniform(mat, order, state);
//                        if (! core_test_pmbasis(mat, order, shift))
//                        { printf("failed %s -- %s,\n...exiting\n", "uniform", "uniform"); return 0; }
//
//                        _test_collection_shift_uniform(shift, rdim);
//                        _test_collection_mat_unbalanced_rdeg(mat, order, state);
//                        if (! core_test_pmbasis(mat, order, shift))
//                        { printf("failed %s -- %s,\n...exiting\n", "uniform", "unbalanced rdeg"); return 0; }
//
//                        _test_collection_shift_uniform(shift, rdim);
//                        _test_collection_mat_unbalanced_cdeg(mat, order, state);
//                        if (! core_test_pmbasis(mat, order, shift))
//                        { printf("failed %s -- %s,\n...exiting\n", "uniform", "unbalanced cdeg"); return 0; }
//
//                        _test_collection_shift_uniform(shift, rdim);
//                        _test_collection_mat_test(mat, order, state);
//                        if (! core_test_pmbasis(mat, order, shift))
//                        { printf("failed %s -- %s,\n...exiting\n", "uniform", "test"); return 0; }
//
//                        _test_collection_shift_uniform(shift, rdim);
//                        _test_collection_mat_sparse(mat, order, state);
//                        if (! core_test_pmbasis(mat, order, shift))
//                        { printf("failed %s -- %s,\n...exiting\n", "uniform", "sparse"); return 0; }
//
//                        _test_collection_shift_uniform(shift, rdim);
//                        _test_collection_mat_rkdef(mat, order, state);
//                        if (! core_test_pmbasis(mat, order, shift))
//                        { printf("failed %s -- %s,\n...exiting\n", "uniform", "rkdef"); return 0; }
//
//                        /* increasing shift */
//                        /* matrix: zero | uniform | rdeg | cdeg | test | sparse | rkdef */
//                        _test_collection_shift_increasing(shift, rdim);
//                        _test_collection_mat_zero(mat);
//                        if (! core_test_pmbasis(mat, order, shift))
//                        { printf("failed %s -- %s,\n...exiting\n", "increasing", "zero"); return 0; }
//
//                        _test_collection_shift_increasing(shift, rdim);
//                        _test_collection_mat_uniform(mat, order, state);
//                        if (! core_test_pmbasis(mat, order, shift))
//                        { printf("failed %s -- %s,\n...exiting\n", "increasing", "uniform"); return 0; }
//
//                        _test_collection_shift_increasing(shift, rdim);
//                        _test_collection_mat_unbalanced_rdeg(mat, order, state);
//                        if (! core_test_pmbasis(mat, order, shift))
//                        { printf("failed %s -- %s,\n...exiting\n", "increasing", "unbalanced rdeg"); return 0; }
//
//                        _test_collection_shift_increasing(shift, rdim);
//                        _test_collection_mat_unbalanced_cdeg(mat, order, state);
//                        if (! core_test_pmbasis(mat, order, shift))
//                        { printf("failed %s -- %s,\n...exiting\n", "increasing", "unbalanced cdeg"); return 0; }
//
//                        _test_collection_shift_increasing(shift, rdim);
//                        _test_collection_mat_test(mat, order, state);
//                        if (! core_test_pmbasis(mat, order, shift))
//                        { printf("failed %s -- %s,\n...exiting\n", "increasing", "test"); return 0; }
//
//                        _test_collection_shift_increasing(shift, rdim);
//                        _test_collection_mat_sparse(mat, order, state);
//                        if (! core_test_pmbasis(mat, order, shift))
//                        { printf("failed %s -- %s,\n...exiting\n", "increasing", "sparse"); return 0; }
//
//                        _test_collection_shift_increasing(shift, rdim);
//                        _test_collection_mat_rkdef(mat, order, state);
//                        if (! core_test_pmbasis(mat, order, shift))
//                        { printf("failed %s -- %s,\n...exiting\n", "increasing", "rkdef"); return 0; }
//
//                        /* decreasing shift */
//                        /* matrix: zero | uniform | rdeg | cdeg | test | sparse | rkdef */
//                        _test_collection_shift_decreasing(shift, rdim);
//                        _test_collection_mat_zero(mat);
//                        if (! core_test_pmbasis(mat, order, shift))
//                        { printf("failed %s -- %s,\n...exiting\n", "decreasing", "zero"); return 0; }
//
//                        _test_collection_shift_decreasing(shift, rdim);
//                        _test_collection_mat_uniform(mat, order, state);
//                        if (! core_test_pmbasis(mat, order, shift))
//                        { printf("failed %s -- %s,\n...exiting\n", "decreasing", "uniform"); return 0; }
//
//                        _test_collection_shift_decreasing(shift, rdim);
//                        _test_collection_mat_unbalanced_rdeg(mat, order, state);
//                        if (! core_test_pmbasis(mat, order, shift))
//                        { printf("failed %s -- %s,\n...exiting\n", "decreasing", "unbalanced rdeg"); return 0; }
//
//                        _test_collection_shift_decreasing(shift, rdim);
//                        _test_collection_mat_unbalanced_cdeg(mat, order, state);
//                        if (! core_test_pmbasis(mat, order, shift))
//                        { printf("failed %s -- %s,\n...exiting\n", "decreasing", "unbalanced cdeg"); return 0; }
//
//                        _test_collection_shift_decreasing(shift, rdim);
//                        _test_collection_mat_test(mat, order, state);
//                        if (! core_test_pmbasis(mat, order, shift))
//                        { printf("failed %s -- %s,\n...exiting\n", "decreasing", "test"); return 0; }
//
//                        _test_collection_shift_decreasing(shift, rdim);
//                        _test_collection_mat_sparse(mat, order, state);
//                        if (! core_test_pmbasis(mat, order, shift))
//                        { printf("failed %s -- %s,\n...exiting\n", "decreasing", "sparse"); return 0; }
//
//                        _test_collection_shift_decreasing(shift, rdim);
//                        _test_collection_mat_rkdef(mat, order, state);
//                        if (! core_test_pmbasis(mat, order, shift))
//                        { printf("failed %s -- %s,\n...exiting\n", "decreasing", "rkdef"); return 0; }
//
//                        /* shuffled shift */
//                        /* matrix: zero | uniform | rdeg | cdeg | test | sparse | rkdef */
//                        _test_collection_shift_shuffle(shift, rdim, state);
//                        _test_collection_mat_zero(mat);
//                        if (! core_test_pmbasis(mat, order, shift))
//                        { printf("failed %s -- %s,\n...exiting\n", "shuffle", "zero"); return 0; }
//
//                        _test_collection_shift_shuffle(shift, rdim, state);
//                        _test_collection_mat_uniform(mat, order, state);
//                        if (! core_test_pmbasis(mat, order, shift))
//                        { printf("failed %s -- %s,\n...exiting\n", "shuffle", "uniform"); return 0; }
//
//                        _test_collection_shift_shuffle(shift, rdim, state);
//                        _test_collection_mat_unbalanced_rdeg(mat, order, state);
//                        if (! core_test_pmbasis(mat, order, shift))
//                        { printf("failed %s -- %s,\n...exiting\n", "shuffle", "unbalanced rdeg"); return 0; }
//
//                        _test_collection_shift_shuffle(shift, rdim, state);
//                        _test_collection_mat_unbalanced_cdeg(mat, order, state);
//                        if (! core_test_pmbasis(mat, order, shift))
//                        { printf("failed %s -- %s,\n...exiting\n", "shuffle", "unbalanced cdeg"); return 0; }
//
//                        _test_collection_shift_shuffle(shift, rdim, state);
//                        _test_collection_mat_test(mat, order, state);
//                        if (! core_test_pmbasis(mat, order, shift))
//                        { printf("failed %s -- %s,\n...exiting\n", "shuffle", "test"); return 0; }
//
//                        _test_collection_shift_shuffle(shift, rdim, state);
//                        _test_collection_mat_sparse(mat, order, state);
//                        if (! core_test_pmbasis(mat, order, shift))
//                        { printf("failed %s -- %s,\n...exiting\n", "shuffle", "sparse"); return 0; }
//
//                        _test_collection_shift_shuffle(shift, rdim, state);
//                        _test_collection_mat_rkdef(mat, order, state);
//                        if (! core_test_pmbasis(mat, order, shift))
//                        { printf("failed %s -- %s,\n...exiting\n", "shuffle", "rkdef"); return 0; }
//
//                        /* HNF shift */
//                        /* matrix: zero | uniform | rdeg | cdeg | test | sparse | rkdef */
//                        _test_collection_shift_hermite(shift, rdim, order);
//                        _test_collection_mat_zero(mat);
//                        if (! core_test_pmbasis(mat, order, shift))
//                        { printf("failed %s -- %s,\n...exiting\n", "hermite", "zero"); return 0; }
//
//                        _test_collection_shift_hermite(shift, rdim, order);
//                        _test_collection_mat_uniform(mat, order, state);
//                        if (! core_test_pmbasis(mat, order, shift))
//                        { printf("failed %s -- %s,\n...exiting\n", "hermite", "uniform"); return 0; }
//
//                        _test_collection_shift_hermite(shift, rdim, order);
//                        _test_collection_mat_unbalanced_rdeg(mat, order, state);
//                        if (! core_test_pmbasis(mat, order, shift))
//                        { printf("failed %s -- %s,\n...exiting\n", "hermite", "unbalanced rdeg"); return 0; }
//
//                        _test_collection_shift_hermite(shift, rdim, order);
//                        _test_collection_mat_unbalanced_cdeg(mat, order, state);
//                        if (! core_test_pmbasis(mat, order, shift))
//                        { printf("failed %s -- %s,\n...exiting\n", "hermite", "unbalanced cdeg"); return 0; }
//
//                        _test_collection_shift_hermite(shift, rdim, order);
//                        _test_collection_mat_test(mat, order, state);
//                        if (! core_test_pmbasis(mat, order, shift))
//                        { printf("failed %s -- %s,\n...exiting\n", "hermite", "test"); return 0; }
//
//                        _test_collection_shift_hermite(shift, rdim, order);
//                        _test_collection_mat_sparse(mat, order, state);
//                        if (! core_test_pmbasis(mat, order, shift))
//                        { printf("failed %s -- %s,\n...exiting\n", "hermite", "sparse"); return 0; }
//
//                        _test_collection_shift_hermite(shift, rdim, order);
//                        _test_collection_mat_rkdef(mat, order, state);
//                        if (! core_test_pmbasis(mat, order, shift))
//                        { printf("failed %s -- %s,\n...exiting\n", "hermite", "rkdef"); return 0; }
//
//                        /* HNF shift, bis */
//                        /* matrix: zero | uniform | rdeg | cdeg | test | sparse | rkdef */
//                        _test_collection_shift_rhermite(shift, rdim, order);
//                        _test_collection_mat_zero(mat);
//                        if (! core_test_pmbasis(mat, order, shift))
//                        { printf("failed %s -- %s,\n...exiting\n", "rhermite", "zero"); return 0; }
//
//                        _test_collection_shift_rhermite(shift, rdim, order);
//                        _test_collection_mat_uniform(mat, order, state);
//                        if (! core_test_pmbasis(mat, order, shift))
//                        { printf("failed %s -- %s,\n...exiting\n", "rhermite", "uniform"); return 0; }
//
//                        _test_collection_shift_rhermite(shift, rdim, order);
//                        _test_collection_mat_unbalanced_rdeg(mat, order, state);
//                        if (! core_test_pmbasis(mat, order, shift))
//                        { printf("failed %s -- %s,\n...exiting\n", "rhermite", "unbalanced rdeg"); return 0; }
//
//                        _test_collection_shift_rhermite(shift, rdim, order);
//                        _test_collection_mat_unbalanced_cdeg(mat, order, state);
//                        if (! core_test_pmbasis(mat, order, shift))
//                        { printf("failed %s -- %s,\n...exiting\n", "rhermite", "unbalanced cdeg"); return 0; }
//
//                        _test_collection_shift_rhermite(shift, rdim, order);
//                        _test_collection_mat_test(mat, order, state);
//                        if (! core_test_pmbasis(mat, order, shift))
//                        { printf("failed %s -- %s,\n...exiting\n", "rhermite", "test"); return 0; }
//
//                        _test_collection_shift_rhermite(shift, rdim, order);
//                        _test_collection_mat_sparse(mat, order, state);
//                        if (! core_test_pmbasis(mat, order, shift))
//                        { printf("failed %s -- %s,\n...exiting\n", "rhermite", "sparse"); return 0; }
//
//                        _test_collection_shift_rhermite(shift, rdim, order);
//                        _test_collection_mat_rkdef(mat, order, state);
//                        if (! core_test_pmbasis(mat, order, shift))
//                        { printf("failed %s -- %s,\n...exiting\n", "rhermite", "rkdef"); return 0; }
//
//                        /* plateau shift */
//                        /* matrix: zero | uniform | rdeg | cdeg | test | sparse | rkdef */
//                        _test_collection_shift_plateau(shift, rdim, order);
//                        _test_collection_mat_zero(mat);
//                        if (! core_test_pmbasis(mat, order, shift))
//                        { printf("failed %s -- %s,\n...exiting\n", "plateau", "zero"); return 0; }
//
//                        _test_collection_shift_plateau(shift, rdim, order);
//                        _test_collection_mat_uniform(mat, order, state);
//                        if (! core_test_pmbasis(mat, order, shift))
//                        { printf("failed %s -- %s,\n...exiting\n", "plateau", "uniform"); return 0; }
//
//                        _test_collection_shift_plateau(shift, rdim, order);
//                        _test_collection_mat_unbalanced_rdeg(mat, order, state);
//                        if (! core_test_pmbasis(mat, order, shift))
//                        { printf("failed %s -- %s,\n...exiting\n", "plateau", "unbalanced rdeg"); return 0; }
//
//                        _test_collection_shift_plateau(shift, rdim, order);
//                        _test_collection_mat_unbalanced_cdeg(mat, order, state);
//                        if (! core_test_pmbasis(mat, order, shift))
//                        { printf("failed %s -- %s,\n...exiting\n", "plateau", "unbalanced cdeg"); return 0; }
//
//                        _test_collection_shift_plateau(shift, rdim, order);
//                        _test_collection_mat_test(mat, order, state);
//                        if (! core_test_pmbasis(mat, order, shift))
//                        { printf("failed %s -- %s,\n...exiting\n", "plateau", "test"); return 0; }
//
//                        _test_collection_shift_plateau(shift, rdim, order);
//                        _test_collection_mat_sparse(mat, order, state);
//                        if (! core_test_pmbasis(mat, order, shift))
//                        { printf("failed %s -- %s,\n...exiting\n", "plateau", "sparse"); return 0; }
//
//                        _test_collection_shift_plateau(shift, rdim, order);
//                        _test_collection_mat_rkdef(mat, order, state);
//                        if (! core_test_pmbasis(mat, order, shift))
//                        { printf("failed %s -- %s,\n...exiting\n", "plateau", "rkdef"); return 0; }
//
//                        /* plateau shift, bis */
//                        /* matrix: zero | uniform | rdeg | cdeg | test | sparse | rkdef */
//                        _test_collection_shift_rplateau(shift, rdim, order);
//                        _test_collection_mat_zero(mat);
//                        if (! core_test_pmbasis(mat, order, shift))
//                        { printf("failed %s -- %s,\n...exiting\n", "rplateau", "zero"); return 0; }
//
//                        _test_collection_shift_rplateau(shift, rdim, order);
//                        _test_collection_mat_uniform(mat, order, state);
//                        if (! core_test_pmbasis(mat, order, shift))
//                        { printf("failed %s -- %s,\n...exiting\n", "rplateau", "uniform"); return 0; }
//
//                        _test_collection_shift_rplateau(shift, rdim, order);
//                        _test_collection_mat_unbalanced_rdeg(mat, order, state);
//                        if (! core_test_pmbasis(mat, order, shift))
//                        { printf("failed %s -- %s,\n...exiting\n", "rplateau", "unbalanced rdeg"); return 0; }
//
//                        _test_collection_shift_rplateau(shift, rdim, order);
//                        _test_collection_mat_unbalanced_cdeg(mat, order, state);
//                        if (! core_test_pmbasis(mat, order, shift))
//                        { printf("failed %s -- %s,\n...exiting\n", "rplateau", "unbalanced cdeg"); return 0; }
//
//                        _test_collection_shift_rplateau(shift, rdim, order);
//                        _test_collection_mat_test(mat, order, state);
//                        if (! core_test_pmbasis(mat, order, shift))
//                        { printf("failed %s -- %s,\n...exiting\n", "rplateau", "test"); return 0; }
//
//                        _test_collection_shift_rplateau(shift, rdim, order);
//                        _test_collection_mat_sparse(mat, order, state);
//                        if (! core_test_pmbasis(mat, order, shift))
//                        { printf("failed %s -- %s,\n...exiting\n", "rplateau", "sparse"); return 0; }
//
//                        _test_collection_shift_rplateau(shift, rdim, order);
//                        _test_collection_mat_rkdef(mat, order, state);
//                        if (! core_test_pmbasis(mat, order, shift))
//                        { printf("failed %s -- %s,\n...exiting\n", "rplateau", "rkdef"); return 0; }
//
//                        nmod_poly_mat_clear(mat);
//                        flint_free(shift);
//                    }
//
//    printf("--> Successful\n");
//    return 1;
//}

TEST_FUNCTION_START(nmod_poly_mat_kernel, state)
{
    int i, result;

    /* TODO activate this for long test */
    if (0)
    {
        for (i = 0; i < 5 * flint_test_multiplier(); i++)
        {
            {
                result = 1;
                /* result = collection_test_kernel(5, state); */

                if (!result)
                    TEST_FUNCTION_FAIL("Failed kernel on testing collection\n");
            }
        }
    }
    else
    {
        for (i = 0; i < 50 * flint_test_multiplier(); i++)
        {
            ulong nbits = 2 + n_randint(state, 63);
            ulong rdim = 1 + n_randint(state, 10);
            ulong cdim = 1 + n_randint(state, 10);
            ulong len = n_randint(state, 150);
            /* ulong rdim = 1 + n_randint(state, 3); */
            /* ulong cdim = 1 + n_randint(state, 3); */
            /* ulong len = n_randint(state, 15); */

            slong prime = n_randprime(state, nbits, 1);

            if (_TEST_VERBOSE_)
                flint_printf("prime = %wu, rdim = %wu, cdim = %wu, len = %wu\n",
                             prime, rdim, cdim, len);
            result = one_test_kernel(prime, rdim, cdim, len, state);

            if (!result)
                TEST_FUNCTION_FAIL("prime = %wu, rdim = %wu, cdim = %wu, len = %wu\n",
                                   prime, rdim, cdim, len);
        }
    }

    TEST_FUNCTION_END(state);
}

#undef _TEST_VERBOSE_
