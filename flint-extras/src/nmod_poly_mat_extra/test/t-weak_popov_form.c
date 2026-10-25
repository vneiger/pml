/*
    Copyright (C) 2025 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <flint/nmod_types.h>
#include <flint/test_helpers.h>
#include <flint/profiler.h>

#include "nmod_poly_mat_forms.h"
#include "nmod_poly_mat_utils.h"

#include "testing_collection.h"

/* TODO for long test, add unbalanced rdeg / unbalanced cdeg variants */

//#define NOTRANS

// stack row at the bottom of mat
// no dimension checking
static void stack_row(nmod_poly_mat_t mat, const nmod_poly_struct * row)
{
    nmod_poly_mat_t tmp;
    nmod_poly_mat_init(tmp, mat->r +1, mat->c, mat->modulus);
    for (slong i = 0; i < mat->r; i++)
        for (slong j = 0; j < mat->c; j++)
            nmod_poly_swap(nmod_poly_mat_entry(tmp,i,j), nmod_poly_mat_entry(mat,i,j));
    for (slong j = 0; j < mat->c; j++)
        nmod_poly_set(nmod_poly_mat_entry(tmp,mat->r,j), row+j);
    nmod_poly_mat_swap(mat, tmp);
    nmod_poly_mat_clear(tmp);
}

// replace last row of mat by ``row``
// no dimension checking
static void replace_last_row(nmod_poly_mat_t mat, const nmod_poly_struct * row)
{
    for (slong j = 0; j < mat->c; j++)
        nmod_poly_set(nmod_poly_mat_entry(mat,mat->r -1,j), row+j);
}

// helper for verifying rank profile
// rank must be <= the length of rrp
// assumption: rank is the rank of mat
int verify_row_rank_profile(const nmod_poly_mat_t mat, const slong * rrp, slong rank)
{
    if (mat->r == 0 || mat->c == 0)
        return 1;

    nmod_poly_mat_t submat;
    nmod_poly_mat_init(submat, 0, mat->c, mat->modulus);
    for (slong i = 0; i <= rank; i++)
    {
        // i is the number of verified rrp entries at this stage;
        // submat has full row rank and i rows

        // first verify that rows in rrp[i-1]+1:rrp[i] do not increase the rank
        slong start = (i==0) ? 0 : rrp[i-1]+1;
        slong end = (i==rank) ? mat->r : rrp[i];
        for (slong ii = start; ii < end; ii++)
        {
            if (ii == start) stack_row(submat, nmod_poly_mat_entry(mat, ii, 0));
            else replace_last_row(submat, nmod_poly_mat_entry(mat, ii, 0));
            if (nmod_poly_mat_rank(submat) != i)
            {
                nmod_poly_mat_clear(submat);
                printf("~~~ verify row rank profile ~~~ INCORRECT: rank increased between two rank profile rows\n");
                return 0;
            }
        }
        // now verify that the row rrp[i] does increase the rank
        if (i < rank)
        {
            if (submat->r == i) // happens when start == end; above for loop did nothing
                stack_row(submat, nmod_poly_mat_entry(mat, rrp[i], 0));
            else
                replace_last_row(submat, nmod_poly_mat_entry(mat, rrp[i], 0));
            if (nmod_poly_mat_rank(submat) != i+1)
            {
                nmod_poly_mat_clear(submat);
                printf("~~~ verify row rank profile ~~~ INCORRECT: rank profile row does not increase rank\n");
                return 0;
            }
        }
    }
    nmod_poly_mat_clear(submat);
    return 1;
}

// helper for verifying pivot index of first rank rows of mat
// rank must be the rank of mat and it is assumed that pivind has >= rank entries
int verify_pivot_index(const nmod_poly_mat_t mat, const slong * shift, const slong * pivind, slong rank)
{
    if (mat->r == 0 || mat->c == 0)
        return 1;

    slong * check_pivind = flint_malloc(mat->r * sizeof(slong));
    nmod_poly_mat_pivot_index(check_pivind, mat, shift, ROW_LOWER);
    int equal = 1;
    for (slong i = 0; i < rank; i++)
        if (pivind[i] != check_pivind[i])
            equal = 0;
    flint_free(check_pivind);
    if (!equal)
        printf("~~~ verify pivot index ~~~ INCORRECT: pivot index is not correct\n");
    return equal;
}

// verify weak Popov form
int verify_weak_popov_form(const nmod_poly_mat_t wpf, const slong * shift, const nmod_poly_mat_t tsf, slong rk, const nmod_poly_mat_t mat, flint_rand_t state)
{
    if (! tsf)
    {
        printf("~~~ verify weak Popov ~~~ transformation not provided, skipping verification\n");
        return 1;
    }

    // testing correctness, from fastest to slowest:
    // 0. dimensions
    if (mat->r != wpf->r || mat->c != wpf->c || mat->r != tsf->r || mat->r != tsf->c)
    {
            printf("~~~ verify weak Popov ~~~ INCORRECT: dimension mismatch\n");
            return 0;
    }

    // empty matrices are fine (checked by is_weak_popov, but done here otherwise
    // tsf does not pass test on tsf_deg >= 0)
    if (mat->r == 0 || mat->c == 0)
        return 1;

    // 1. unimodular transformation
    // upper bound on degree of determinant
    slong tsf_deg = nmod_poly_mat_degree(tsf);
    if (tsf_deg < 0)
    {
        printf("~~~ verify weak Popov ~~~ INCORRECT: transformation is zero\n");
        return 0;
    }
    ulong degdet = tsf->r * (ulong)tsf_deg;
    if (mat->modulus > 100*degdet)
    {
        if (!nmod_poly_mat_is_unimodular_randomized(tsf, state))
        {
            printf("~~~ verify weak Popov ~~~ INCORRECT: transformation is not unimodular (randomized test)\n");
            return 0;
        }
    }
    else
    {
        if (!nmod_poly_mat_is_unimodular(tsf))
        {
            printf("~~~ verify weak Popov ~~~ INCORRECT: transformation is not unimodular\n");
            return 0;
        }
    }

    // 2. tsf * mat == wpf
    nmod_poly_mat_t tmp;
    nmod_poly_mat_init(tmp, mat->r, mat->c, mat->modulus);
    nmod_poly_mat_mul(tmp, tsf, mat);
    if (! nmod_poly_mat_equal(tmp, wpf))
    {
        printf("~~~ verify weak Popov ~~~ INCORRECT: tsf * mat != wpf\n");
        nmod_poly_mat_clear(tmp);
        return 0;
    }

    // 3. wpf is in weak Popov form
    // first prune zero rows since this is a requirement of is_weak_popov
    // (note that nonzero rows are not necessarily the ones at rk ... rdim-1)
    nmod_poly_mat_clear(tmp);
    nmod_poly_mat_init(tmp, rk, wpf->c, wpf->modulus);
    slong current_row = 0;
    for (slong i = 0; i < wpf->r; i++)
    {
        slong jj = 0;
        while (jj < mat->c && nmod_poly_is_zero(nmod_poly_mat_entry(wpf, i, jj)))
            jj++;
        if (jj < mat->c) // row i is nonzero, jj is first nonzero
        {
            for (slong j = jj; j < wpf->c; j++)
                nmod_poly_set(nmod_poly_mat_entry(tmp, current_row, j), nmod_poly_mat_entry(wpf, i, j));
            current_row++;
        }
    }
    if (! nmod_poly_mat_is_weak_popov(tmp, shift, ROW_LOWER))
    {
        printf("~~~ verify weak Popov ~~~ INCORRECT: wpf not in shifted weak Popov form\n");
        nmod_poly_mat_clear(tmp);
        return 0;
    }

    nmod_poly_mat_clear(tmp);
    return 1;
}

// verify ordered weak Popov form
int verify_ordered_weak_popov_form(const nmod_poly_mat_t wpf, const slong * shift, const nmod_poly_mat_t tsf, slong rk, const nmod_poly_mat_t mat, flint_rand_t state)
{
    if (! tsf)
    {
        printf("~~~ verify weak Popov ~~~ transformation not provided, skipping verification\n");
        return 1;
    }

    // testing correctness, from fastest to slowest:
    // 0. dimensions
    if (mat->r != wpf->r || mat->c != wpf->c || mat->r != tsf->r || mat->r != tsf->c)
    {
            printf("~~~ verify weak Popov ~~~ INCORRECT: dimension mismatch\n");
            return 0;
    }

    // empty matrices are fine (checked by is_weak_popov, but done here otherwise
    // tsf does not pass test on tsf_deg >= 0)
    if (mat->r == 0 || mat->c == 0)
        return 1;

    // 1. unimodular transformation
    // upper bound on degree of determinant
    slong tsf_deg = nmod_poly_mat_degree(tsf);
    if (tsf_deg < 0)
    {
        printf("~~~ verify weak Popov ~~~ INCORRECT: transformation is zero\n");
        return 0;
    }
    ulong degdet = tsf->r * (ulong)tsf_deg;
    if (mat->modulus > 100*degdet)
    {
        if (!nmod_poly_mat_is_unimodular_randomized(tsf, state))
        {
            printf("~~~ verify weak Popov ~~~ INCORRECT: transformation is not unimodular (randomized test)\n");
            return 0;
        }
    }
    else
    {
        if (!nmod_poly_mat_is_unimodular(tsf))
        {
            printf("~~~ verify weak Popov ~~~ INCORRECT: transformation is not unimodular\n");
            return 0;
        }
    }

    // 2. tsf * mat == wpf
    nmod_poly_mat_t tmp;
    nmod_poly_mat_init(tmp, mat->r, mat->c, mat->modulus);
    nmod_poly_mat_mul(tmp, tsf, mat);
    if (! nmod_poly_mat_equal(tmp, wpf))
    {
        printf("~~~ verify weak Popov ~~~ INCORRECT: tsf * mat != wpf\n");
        nmod_poly_mat_clear(tmp);
        return 0;
    }

    // 3. wpf is in weak Popov form
    // first prune zero rows since this is a requirement of is_ordered_weak_popov
    // (note that zero rows are necessarily the ones at rk ... rdim-1)
    nmod_poly_mat_clear(tmp);
    nmod_poly_mat_init(tmp, rk, wpf->c, wpf->modulus);
    for (slong i = 0; i < rk; i++)
        for (slong j = 0; j < mat->c; j++)
            nmod_poly_set(nmod_poly_mat_entry(tmp, i, j), nmod_poly_mat_entry(wpf, i, j));

    if (! nmod_poly_mat_is_ordered_weak_popov(tmp, shift, ROW_LOWER))
    {
        printf("~~~ verify ord weak Popov ~~~ INCORRECT: wpf not in shifted ordered weak Popov form\n");
        nmod_poly_mat_clear(tmp);
        return 0;
    }

    nmod_poly_mat_clear(tmp);
    return 1;
}

// test one given input for weak Popov form
int core_test_weak_popov_form(const nmod_poly_mat_t mat, const slong * shift, int time, flint_rand_t state)
{
    // initialize pivind list and row rank profile
    slong max_rank = FLINT_MIN(mat->r, mat->c);
    slong * pivind = flint_malloc(FLINT_MIN(max_rank+1,mat->r) * sizeof(slong)); // see doc for the +1
    slong * rrp = flint_malloc(max_rank * sizeof(slong));
    // init copy of mat
    nmod_poly_mat_t wpf;
    nmod_poly_mat_init(wpf, mat->r, mat->c, mat->modulus);

    // init unimodular transformation tsf such that wpf = tsf * mat
    nmod_poly_mat_t tsf;
    nmod_poly_mat_init(tsf, mat->r, mat->r, mat->modulus);

    // verification of weak Popov form, verification of row rank profile
    int verif_form, verif_pivind, verif_rrp;

    { // Mulders and Storjohann's iterative algorithm, row by row
        nmod_poly_mat_set(wpf, mat);
        nmod_poly_mat_one(tsf);
        timeit_t timer;
        timeit_start(timer);
        // weak Popov form computation
#ifdef NOTRANS
        slong rk = _nmod_poly_mat_weak_popov_iter_submat_rowbyrow(wpf, shift, NULL, NULL, pivind, rrp, 0, 0, mat->r, mat->c, mat->r, ROW_LOWER);
#else
        slong rk = _nmod_poly_mat_weak_popov_iter_submat_rowbyrow(wpf, shift, tsf, NULL, pivind, rrp, 0, 0, mat->r, mat->c, mat->r, ROW_LOWER);
#endif /* ifdef NOTRANS */

        timeit_stop(timer);
        if (time)
            flint_printf("-- time (iter - rowbyrow - ur): %wd ms\n", timer->wall);

        // weak Popov form verification
        timeit_start(timer);
#ifdef NOTRANS
        verif_form = verify_weak_popov_form(wpf, shift, NULL, rk, mat, state);
#else
        verif_form = verify_weak_popov_form(wpf, shift, tsf, rk, mat, state);
#endif /* ifdef NOTRANS */
        verif_pivind = verify_pivot_index(wpf, shift, pivind, rk);
        verif_rrp = verify_row_rank_profile(mat, rrp, rk);
        if (! verif_form)
            printf("weak Popov form -- Mulders-Storjohann -- form failure.\n");
        if (! verif_pivind)
            printf("weak Popov form -- Mulders-Storjohann -- pivot index failure.\n");
        if (! verif_rrp)
            printf("weak Popov form -- Mulders-Storjohann -- rank profile failure.\n");
        timeit_stop(timer);
        if (time)
            flint_printf("-- time (verif-wpf+pivind+rrp): %wd ms\n", timer->wall);

        // ordered weak Popov form computation
        timeit_start(timer);
        slong * perm = flint_malloc(wpf->r * sizeof(slong));
        _nmod_poly_mat_permute_rows_by_sorting_vec(wpf, rk, pivind, perm);
#ifndef NOTRANS
        nmod_poly_mat_permute_rows(tsf, perm, NULL);
#endif
        timeit_stop(timer);
        if (time)
            flint_printf("-- time (wpf -> ordered wpf): %wd ms\n", timer->wall);

        // ordered weak Popov form verification
        timeit_start(timer);
#ifdef NOTRANS
        verif_form = verify_ordered_weak_popov_form(wpf, shift, NULL, rk, mat, state);
#else
        verif_form = verify_ordered_weak_popov_form(wpf, shift, tsf, rk, mat, state);
#endif /* ifdef NOTRANS */
        verif_pivind = verify_pivot_index(wpf, shift, pivind, rk);
        flint_free(perm);
        timeit_stop(timer);
        if (! verif_form)
            printf("weak Popov form -- Mulders-Storjohann -- form failure.\n");
        if (! verif_pivind)
            printf("weak Popov form -- Mulders-Storjohann -- pivot index failure.\n");
        if (time)
            flint_printf("-- time (verif-ord-wpf+pivind): %wd ms\n", timer->wall);
    }

    nmod_poly_mat_clear(wpf);
    nmod_poly_mat_clear(tsf);
    flint_free(pivind);
    flint_free(rrp);

    return verif_form && verif_rrp;
}

//** Test against the whole testing collection */
int collection_test_weak_popov_form(slong iter, flint_rand_t state)
{
    // input matrix
    nmod_poly_mat_t mat;

    // input shift
    slong * shift;

    long total_nb_tests =
            iter // number of iterations
            * 40 // number of mats (currently 5) x number of shifts (currently 8)
            * _test_collection_nb_primes
            * _test_collection_nb_smalldims
            * _test_collection_nb_smalldims
            * _test_collection_nb_minidegs;

    printf("Launching testing collection (%ld cases)\n", total_nb_tests);

    for (slong i_primes = 0; i_primes < _test_collection_nb_primes; i_primes++)
        for (slong i_rdims = 0; i_rdims < _test_collection_nb_minidims; i_rdims++)
            for (slong i_cdims = 0; i_cdims < _test_collection_nb_minidims; i_cdims++)
                for (slong i_degs = 0; i_degs < _test_collection_nb_minidegs; i_degs++)
                    for (slong it = 0; it < iter; it++)
                    {
                        const long prime = _test_collection_primes[i_primes];
                        const long rdim = _test_collection_minirdims[i_rdims];
                        const long cdim = _test_collection_minicdims[i_cdims];
                        const long len = _test_collection_minidegs[i_degs];
                        printf("prime %ld, rdim %ld, cdim %ld, length %ld.\n", prime, rdim, cdim, len);

                        nmod_poly_mat_init(mat, rdim, cdim, prime);
                        shift = (slong *) flint_malloc(cdim * sizeof(slong));

                        _test_collection_shift_uniform(shift, cdim);
                        _test_collection_mat_zero(mat);
                        if (! core_test_weak_popov_form(mat, shift, 0, state))
                        { printf("failed %s -- %s,\n...exiting\n", "uniform", "zero"); return 0; }

                        _test_collection_shift_uniform(shift, cdim);
                        _test_collection_mat_uniform(mat, len, state);
                        if (! core_test_weak_popov_form(mat, shift, 0, state))
                        { printf("failed %s -- %s,\n...exiting\n", "uniform", "uniform"); return 0; }

                        _test_collection_shift_uniform(shift, cdim);
                        _test_collection_mat_test(mat, len, state);
                        if (! core_test_weak_popov_form(mat, shift, 0, state))
                        { printf("failed %s -- %s,\n...exiting\n", "uniform", "test"); return 0; }

                        _test_collection_shift_uniform(shift, cdim);
                        _test_collection_mat_sparse(mat, len, state);
                        if (! core_test_weak_popov_form(mat, shift, 0, state))
                        { printf("failed %s -- %s,\n...exiting\n", "uniform", "sparse"); return 0; }

                        _test_collection_shift_uniform(shift, cdim);
                        _test_collection_mat_rkdef(mat, len, state);
                        if (! core_test_weak_popov_form(mat, shift, 0, state))
                        { printf("failed %s -- %s,\n...exiting\n", "uniform", "rkdef"); return 0; }


                        _test_collection_shift_increasing(shift, cdim);
                        _test_collection_mat_zero(mat);
                        if (! core_test_weak_popov_form(mat, shift, 0, state))
                        { printf("failed %s -- %s,\n...exiting\n", "increasing", "zero"); return 0; }

                        _test_collection_shift_increasing(shift, cdim);
                        _test_collection_mat_uniform(mat, len, state);
                        if (! core_test_weak_popov_form(mat, shift, 0, state))
                        { printf("failed %s -- %s,\n...exiting\n", "increasing", "uniform"); return 0; }

                        _test_collection_shift_increasing(shift, cdim);
                        _test_collection_mat_test(mat, len, state);
                        if (! core_test_weak_popov_form(mat, shift, 0, state))
                        { printf("failed %s -- %s,\n...exiting\n", "increasing", "test"); return 0; }

                        _test_collection_shift_increasing(shift, cdim);
                        _test_collection_mat_sparse(mat, len, state);
                        if (! core_test_weak_popov_form(mat, shift, 0, state))
                        { printf("failed %s -- %s,\n...exiting\n", "increasing", "sparse"); return 0; }

                        _test_collection_shift_increasing(shift, cdim);
                        _test_collection_mat_rkdef(mat, len, state);
                        if (! core_test_weak_popov_form(mat, shift, 0, state))
                        { printf("failed %s -- %s,\n...exiting\n", "increasing", "rkdef"); return 0; }


                        _test_collection_shift_decreasing(shift, cdim);
                        _test_collection_mat_zero(mat);
                        if (! core_test_weak_popov_form(mat, shift, 0, state))
                        { printf("failed %s -- %s,\n...exiting\n", "decreasing", "zero"); return 0; }

                        _test_collection_shift_decreasing(shift, cdim);
                        _test_collection_mat_uniform(mat, len, state);
                        if (! core_test_weak_popov_form(mat, shift, 0, state))
                        { printf("failed %s -- %s,\n...exiting\n", "decreasing", "uniform"); return 0; }

                        _test_collection_shift_decreasing(shift, cdim);
                        _test_collection_mat_test(mat, len, state);
                        if (! core_test_weak_popov_form(mat, shift, 0, state))
                        { printf("failed %s -- %s,\n...exiting\n", "decreasing", "test"); return 0; }

                        _test_collection_shift_decreasing(shift, cdim);
                        _test_collection_mat_sparse(mat, len, state);
                        if (! core_test_weak_popov_form(mat, shift, 0, state))
                        { printf("failed %s -- %s,\n...exiting\n", "decreasing", "sparse"); return 0; }

                        _test_collection_shift_decreasing(shift, cdim);
                        _test_collection_mat_rkdef(mat, len, state);
                        if (! core_test_weak_popov_form(mat, shift, 0, state))
                        { printf("failed %s -- %s,\n...exiting\n", "decreasing", "rkdef"); return 0; }


                        _test_collection_shift_shuffle(shift, cdim, state);
                        _test_collection_mat_zero(mat);
                        if (! core_test_weak_popov_form(mat, shift, 0, state))
                        { printf("failed %s -- %s,\n...exiting\n", "shuffle", "zero"); return 0; }

                        _test_collection_shift_shuffle(shift, cdim, state);
                        _test_collection_mat_uniform(mat, len, state);
                        if (! core_test_weak_popov_form(mat, shift, 0, state))
                        { printf("failed %s -- %s,\n...exiting\n", "shuffle", "uniform"); return 0; }

                        _test_collection_shift_shuffle(shift, cdim, state);
                        _test_collection_mat_test(mat, len, state);
                        if (! core_test_weak_popov_form(mat, shift, 0, state))
                        { printf("failed %s -- %s,\n...exiting\n", "shuffle", "test"); return 0; }

                        _test_collection_shift_shuffle(shift, cdim, state);
                        _test_collection_mat_sparse(mat, len, state);
                        if (! core_test_weak_popov_form(mat, shift, 0, state))
                        { printf("failed %s -- %s,\n...exiting\n", "shuffle", "sparse"); return 0; }

                        _test_collection_shift_shuffle(shift, cdim, state);
                        _test_collection_mat_rkdef(mat, len, state);
                        if (! core_test_weak_popov_form(mat, shift, 0, state))
                        { printf("failed %s -- %s,\n...exiting\n", "shuffle", "rkdef"); return 0; }


                        _test_collection_shift_hermite(shift, cdim, len);
                        _test_collection_mat_zero(mat);
                        if (! core_test_weak_popov_form(mat, shift, 0, state))
                        { printf("failed %s -- %s,\n...exiting\n", "hermite", "zero"); return 0; }

                        _test_collection_shift_hermite(shift, cdim, len);
                        _test_collection_mat_uniform(mat, len, state);
                        if (! core_test_weak_popov_form(mat, shift, 0, state))
                        { printf("failed %s -- %s,\n...exiting\n", "hermite", "uniform"); return 0; }

                        _test_collection_shift_hermite(shift, cdim, len);
                        _test_collection_mat_test(mat, len, state);
                        if (! core_test_weak_popov_form(mat, shift, 0, state))
                        { printf("failed %s -- %s,\n...exiting\n", "hermite", "test"); return 0; }

                        _test_collection_shift_hermite(shift, cdim, len);
                        _test_collection_mat_sparse(mat, len, state);
                        if (! core_test_weak_popov_form(mat, shift, 0, state))
                        { printf("failed %s -- %s,\n...exiting\n", "hermite", "sparse"); return 0; }

                        _test_collection_shift_hermite(shift, cdim, len);
                        _test_collection_mat_rkdef(mat, len, state);
                        if (! core_test_weak_popov_form(mat, shift, 0, state))
                        { printf("failed %s -- %s,\n...exiting\n", "hermite", "rkdef"); return 0; }


                        _test_collection_shift_rhermite(shift, cdim, len);
                        _test_collection_mat_zero(mat);
                        if (! core_test_weak_popov_form(mat, shift, 0, state))
                        { printf("failed %s -- %s,\n...exiting\n", "rhermite", "zero"); return 0; }

                        _test_collection_shift_rhermite(shift, cdim, len);
                        _test_collection_mat_uniform(mat, len, state);
                        if (! core_test_weak_popov_form(mat, shift, 0, state))
                        { printf("failed %s -- %s,\n...exiting\n", "rhermite", "uniform"); return 0; }

                        _test_collection_shift_rhermite(shift, cdim, len);
                        _test_collection_mat_test(mat, len, state);
                        if (! core_test_weak_popov_form(mat, shift, 0, state))
                        { printf("failed %s -- %s,\n...exiting\n", "rhermite", "test"); return 0; }

                        _test_collection_shift_rhermite(shift, cdim, len);
                        _test_collection_mat_sparse(mat, len, state);
                        if (! core_test_weak_popov_form(mat, shift, 0, state))
                        { printf("failed %s -- %s,\n...exiting\n", "rhermite", "sparse"); return 0; }

                        _test_collection_shift_rhermite(shift, cdim, len);
                        _test_collection_mat_rkdef(mat, len, state);
                        if (! core_test_weak_popov_form(mat, shift, 0, state))
                        { printf("failed %s -- %s,\n...exiting\n", "rhermite", "rkdef"); return 0; }


                        _test_collection_shift_plateau(shift, cdim, len);
                        _test_collection_mat_zero(mat);
                        if (! core_test_weak_popov_form(mat, shift, 0, state))
                        { printf("failed %s -- %s,\n...exiting\n", "plateau", "zero"); return 0; }

                        _test_collection_shift_plateau(shift, cdim, len);
                        _test_collection_mat_uniform(mat, len, state);
                        if (! core_test_weak_popov_form(mat, shift, 0, state))
                        { printf("failed %s -- %s,\n...exiting\n", "plateau", "uniform"); return 0; }

                        _test_collection_shift_plateau(shift, cdim, len);
                        _test_collection_mat_test(mat, len, state);
                        if (! core_test_weak_popov_form(mat, shift, 0, state))
                        { printf("failed %s -- %s,\n...exiting\n", "plateau", "test"); return 0; }

                        _test_collection_shift_plateau(shift, cdim, len);
                        _test_collection_mat_sparse(mat, len, state);
                        if (! core_test_weak_popov_form(mat, shift, 0, state))
                        { printf("failed %s -- %s,\n...exiting\n", "plateau", "sparse"); return 0; }

                        _test_collection_shift_plateau(shift, cdim, len);
                        _test_collection_mat_rkdef(mat, len, state);
                        if (! core_test_weak_popov_form(mat, shift, 0, state))
                        { printf("failed %s -- %s,\n...exiting\n", "plateau", "rkdef"); return 0; }


                        _test_collection_shift_rplateau(shift, cdim, len);
                        _test_collection_mat_zero(mat);
                        if (! core_test_weak_popov_form(mat, shift, 0, state))
                        { printf("failed %s -- %s,\n...exiting\n", "rplateau", "zero"); return 0; }

                        _test_collection_shift_rplateau(shift, cdim, len);
                        _test_collection_mat_uniform(mat, len, state);
                        if (! core_test_weak_popov_form(mat, shift, 0, state))
                        { printf("failed %s -- %s,\n...exiting\n", "rplateau", "uniform"); return 0; }

                        _test_collection_shift_rplateau(shift, cdim, len);
                        _test_collection_mat_test(mat, len, state);
                        if (! core_test_weak_popov_form(mat, shift, 0, state))
                        { printf("failed %s -- %s,\n...exiting\n", "rplateau", "test"); return 0; }

                        _test_collection_shift_rplateau(shift, cdim, len);
                        _test_collection_mat_sparse(mat, len, state);
                        if (! core_test_weak_popov_form(mat, shift, 0, state))
                        { printf("failed %s -- %s,\n...exiting\n", "rplateau", "sparse"); return 0; }

                        _test_collection_shift_rplateau(shift, cdim, len);
                        _test_collection_mat_rkdef(mat, len, state);
                        if (! core_test_weak_popov_form(mat, shift, 0, state))
                        { printf("failed %s -- %s,\n...exiting\n", "rplateau", "rkdef"); return 0; }

                        nmod_poly_mat_clear(mat);
                        flint_free(shift);
                    }

    printf("--> Successful\n");
    return 1;
}

/* pivot profile of zero vectors, for all orientations, shifted or not; and
 * row-wise upper shifted weak Popov forms of matrices with zero rows (which
 * relies on the pivot index convention for zero rows) */
static int test_upper_weak_popov_zero_rows(flint_rand_t state)
{
    const orientation_t orients[4] = {ROW_LOWER, ROW_UPPER, COL_LOWER, COL_UPPER};
    const slong len = 3;
    const slong vshift[3] = {2, -1, 5};
    nmod_poly_struct zvec[3];
    for (slong j = 0; j < len; j++)
        nmod_poly_init(zvec + j, 7);
    int res = 1;
    for (slong k = 0; k < 4; k++)
        for (slong sh = 0; sh < 2; sh++)
        {
            slong pivind = -2, pivdeg = -2;
            const slong zpiv = (orients[k] == ROW_LOWER || orients[k] == COL_UPPER) ? -1 : len;
            ulong nz = _nmod_poly_vec_pivot_profile(&pivind, &pivdeg, zvec, sh ? vshift : NULL,
                                                    len, orients[k]);
            if (nz != 0 || pivind != zpiv || pivdeg != -1)
            {
                printf("FAIL: pivot profile of zero vector, orientation %d, shifted %d\n",
                       (int) orients[k], (int) sh);
                res = 0;
            }
        }
    for (slong j = 0; j < len; j++)
        nmod_poly_clear(zvec + j);

    for (slong iter = 0; res && iter < 200; iter++)
    {
        const slong rdim = 1 + n_randint(state, 6);
        const slong cdim = 1 + n_randint(state, 6);
        const ulong prime = n_randprime(state, 2 + n_randint(state, 30), 1);
        nmod_poly_mat_t mat, wpf;
        nmod_poly_mat_init(mat, rdim, cdim, prime);
        nmod_poly_mat_init(wpf, rdim, cdim, prime);
        if (iter == 0)
        {
            /* [[1, 0], [0, 0]] */
            nmod_poly_mat_clear(mat);
            nmod_poly_mat_init(mat, 2, 2, 7);
            nmod_poly_one(nmod_poly_mat_entry(mat, 0, 0));
            nmod_poly_mat_clear(wpf);
            nmod_poly_mat_init(wpf, 2, 2, 7);
        }
        else
        {
            nmod_poly_mat_randtest_sparse(mat, state, 1 + n_randint(state, 5), 0.5);
            for (slong i = 0; i < rdim; i++)
                if (n_randint(state, 3) == 0)
                    for (slong j = 0; j < cdim; j++)
                        nmod_poly_zero(nmod_poly_mat_entry(mat, i, j));
        }
        slong * shift = flint_malloc(mat->c * sizeof(slong));
        for (slong j = 0; j < mat->c; j++)
            shift[j] = (iter == 0) ? 0 : (slong) n_randint(state, 7) - 3;
        slong * pivind = flint_malloc(mat->r * sizeof(slong));
        slong * check = flint_malloc(mat->r * sizeof(slong));
        nmod_poly_mat_set(wpf, mat);
        slong rk = nmod_poly_mat_weak_popov_iter(wpf, shift, NULL, pivind, NULL, ROW_UPPER);

        /* rank, zero rows at the bottom, weak Popov form of nonzero rows */
        if (rk != nmod_poly_mat_rank(mat))
            res = 0;
        for (slong i = rk; res && i < mat->r; i++)
            for (slong j = 0; j < mat->c; j++)
                if (!nmod_poly_is_zero(nmod_poly_mat_entry(wpf, i, j)))
                    res = 0;
        if (res)
        {
            nmod_poly_mat_t top;
            nmod_poly_mat_window_init(top, wpf, 0, 0, rk, mat->c);
            res = nmod_poly_mat_is_weak_popov(top, shift, ROW_UPPER);
            nmod_poly_mat_pivot_index(check, top, shift, ROW_UPPER);
            for (slong i = 0; res && i < rk; i++)
                if (check[i] != pivind[i])
                    res = 0;
            nmod_poly_mat_window_clear(top);
        }
        if (!res)
        {
            printf("FAIL: upper shifted weak Popov form with zero rows\n");
            nmod_poly_mat_print_pretty(mat, "x");
        }
        flint_free(check);
        flint_free(pivind);
        flint_free(shift);
        nmod_poly_mat_clear(wpf);
        nmod_poly_mat_clear(mat);
    }
    return res;
}

TEST_FUNCTION_START(nmod_poly_mat_weak_popov_form, state)
{
    int i, result;

    if (!test_upper_weak_popov_zero_rows(state))
        TEST_FUNCTION_FAIL("upper shifted weak Popov form, matrix with zero rows\n");

    /* TODO activate this for long test */
    if (0)
    {
        for (i = 0; i < 2 * flint_test_multiplier(); i++)
        {
            {
                result = collection_test_weak_popov_form(2, state);

                if (!result)
                    TEST_FUNCTION_FAIL("Failed weak Popov form on testing collection\n");
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
            ulong deg = n_randint(state, 40);

            slong prime = n_randprime(state, nbits, 1);

            nmod_poly_mat_t mat;
            nmod_poly_mat_init(mat, rdim, cdim, prime);
            slong * shift = (slong *) flint_malloc(cdim * sizeof(slong));

            {
                nmod_poly_mat_randtest(mat, state, deg);

                _test_collection_shift_uniform(shift, cdim);
                result = core_test_weak_popov_form(mat, shift, 0, state);

                _test_collection_shift_increasing(shift, cdim);
                result = result && core_test_weak_popov_form(mat, shift, 0, state);

                _test_collection_shift_decreasing(shift, cdim);
                result = result && core_test_weak_popov_form(mat, shift, 0, state);

                _test_collection_shift_shuffle(shift, cdim, state);
                result = result && core_test_weak_popov_form(mat, shift, 0, state);

                _test_collection_shift_hermite(shift, cdim, deg);
                result = result && core_test_weak_popov_form(mat, shift, 0, state);

                _test_collection_shift_rhermite(shift, cdim, deg);
                result = result && core_test_weak_popov_form(mat, shift, 0, state);

                _test_collection_shift_plateau(shift, cdim, deg);
                result = result && core_test_weak_popov_form(mat, shift, 0, state);

                _test_collection_shift_rplateau(shift, cdim, deg);
                result = result && core_test_weak_popov_form(mat, shift, 0, state);

                if (!result)
                    TEST_FUNCTION_FAIL("prime = %wu, rdim = %wu, cdim = %wu, deg = %wu\n",
                                       prime, rdim, cdim, deg);
            }

            {
                nmod_poly_mat_randtest_sparse(mat, state, deg, 0.2);

                _test_collection_shift_uniform(shift, cdim);
                result = core_test_weak_popov_form(mat, shift, 0, state);

                _test_collection_shift_increasing(shift, cdim);
                result = result && core_test_weak_popov_form(mat, shift, 0, state);

                _test_collection_shift_decreasing(shift, cdim);
                result = result && core_test_weak_popov_form(mat, shift, 0, state);

                _test_collection_shift_shuffle(shift, cdim, state);
                result = result && core_test_weak_popov_form(mat, shift, 0, state);

                _test_collection_shift_hermite(shift, cdim, deg);
                result = result && core_test_weak_popov_form(mat, shift, 0, state);

                _test_collection_shift_rhermite(shift, cdim, deg);
                result = result && core_test_weak_popov_form(mat, shift, 0, state);

                _test_collection_shift_plateau(shift, cdim, deg);
                result = result && core_test_weak_popov_form(mat, shift, 0, state);

                _test_collection_shift_rplateau(shift, cdim, deg);
                result = result && core_test_weak_popov_form(mat, shift, 0, state);

                if (!result)
                    TEST_FUNCTION_FAIL("prime = %wu, rdim = %wu, cdim = %wu, deg = %wu\n",
                                       prime, rdim, cdim, deg);
            }

            flint_free(shift);
            nmod_poly_mat_clear(mat);
        }
    }

    TEST_FUNCTION_END(state);
}
