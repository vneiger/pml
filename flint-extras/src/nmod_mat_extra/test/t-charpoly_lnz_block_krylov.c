/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/nmod.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_poly.h>
#include <flint/test_helpers.h>
#include <flint/ulong_extras.h>

#include "nmod_mat_extra.h"

/* matrix types: 0 random, 1 zero, 2 scalar, 3 diagonal, 4 upper triangular,
 * 5 lower triangular, 6 upper Hessenberg, 7 companion, 8 block diagonal
 * diag(B, C), 9 similar to diag(B, B) (derogatory), 10 sparse */
#define _NTYPES 11

static void _rand_structured(nmod_mat_t A, slong type, flint_rand_t state)
{
    const slong n = A->r;
    const ulong p = A->mod.n;
    nmod_mat_zero(A);

    if (type == 0)
        nmod_mat_rand(A, state);
    else if (type == 2)
    {
        const ulong c = n_randint(state, p);
        for (slong i = 0; i < n; i++)
            nmod_mat_entry(A, i, i) = c;
    }
    else if (type == 3 || type == 4 || type == 5 || type == 6 || type == 10)
    {
        for (slong i = 0; i < n; i++)
            for (slong j = 0; j < n; j++)
                if ((type == 3 && i == j) || (type == 4 && i <= j)
                    || (type == 5 && i >= j) || (type == 6 && i <= j + 1)
                    || (type == 10 && n_randint(state, 8) == 0))
                    nmod_mat_entry(A, i, j) = n_randint(state, p);
    }
    else if (type == 7)
    {
        for (slong i = 1; i < n; i++)
            nmod_mat_entry(A, i, i - 1) = 1;
        for (slong i = 0; i < n; i++)
            nmod_mat_entry(A, i, n - 1) = n_randint(state, p);
    }
    else if (type == 8 || type == 9)
    {
        const slong b = n / 2;
        for (slong i = 0; i < n; i++)
            for (slong j = 0; j < n; j++)
                if ((i < b) == (j < b))
                    nmod_mat_entry(A, i, j) = n_randint(state, p);
        if (type == 9)
        {
            for (slong i = 0; i < b; i++)
                for (slong j = 0; j < b; j++)
                    nmod_mat_entry(A, n - b + i, n - b + j) = nmod_mat_entry(A, i, j);
            /* random similarity, through random elementary operations */
            nmod_mat_t P, Pinv, T;
            nmod_mat_init(P, n, n, p);
            nmod_mat_init(Pinv, n, n, p);
            nmod_mat_init(T, n, n, p);
            nmod_mat_one(P);
            nmod_mat_randops(P, state, 2 * n);
            if (nmod_mat_inv(Pinv, P))
            {
                nmod_mat_mul(T, P, A);
                nmod_mat_mul(A, T, Pinv);
            }
            nmod_mat_clear(P);
            nmod_mat_clear(Pinv);
            nmod_mat_clear(T);
        }
    }
    /* type 1: zero */
}

TEST_FUNCTION_START(nmod_mat_charpoly_lnz_block_krylov, state)
{
    for (slong iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        /* modulus: prime of 2 to 64 bits, sometimes a large one, for which
         * the reduction must succeed on random matrices */
        const ulong p = (n_randint(state, 3) == 0)
                            ? n_randprime(state, 52 + n_randint(state, 13), 1)
                            : n_randprime(state, 2 + n_randint(state, 63), 1);

        const slong n = (iter % 10 == 0) ? n_randint(state, 120) : n_randint(state, 40);
        const slong type = (iter % 2 == 0) ? 0 : n_randint(state, _NTYPES);
        slong nlevels;
        switch (n_randint(state, 5))
        {
            case 0: nlevels = 1; break;
            case 1: nlevels = 2; break;
            case 2: nlevels = 4; break;
            case 3: nlevels = n; break;
            default: nlevels = n_randint(state, n + 3);
        }

        nmod_mat_t A;
        nmod_mat_init(A, n, n, p);
        _rand_structured(A, type, state);

        nmod_poly_t cp, cp_ref;
        nmod_poly_init(cp, p);
        nmod_poly_init(cp_ref, p);
        nmod_mat_charpoly(cp_ref, A);

        /* the reduction may fail, then cp is unchanged */
        nmod_poly_set(cp, cp_ref);
        nmod_poly_set_coeff_ui(cp, n + 1, 1);  /* marker, not a charpoly */
        const int ok = _nmod_mat_charpoly_lnz_block_krylov(cp, A, nlevels);
        if (ok && !nmod_poly_equal(cp, cp_ref))
            TEST_FUNCTION_FAIL("wrong result\n"
                               "n = %wd, type = %wd, p = %wu, nlevels = %wd\n",
                               n, type, p, nlevels);
        if (!ok && (nmod_poly_degree(cp) != n + 1))
            TEST_FUNCTION_FAIL("failed and modified cp\n"
                               "n = %wd, type = %wd, p = %wu, nlevels = %wd\n",
                               n, type, p, nlevels);

        /* no reduction (k == n): always succeeds; random matrices modulo a
         * prime of at least 52 bits: succeeds except with probability at
         * most about n^2 / 2^51 */
        const slong nl = FLINT_MAX(1, FLINT_MIN(nlevels, n));
        const slong k = (n + nl - 1) / nl;
        if (!ok && (k == n || (type == 0 && FLINT_BIT_COUNT(p) >= 52)))
            TEST_FUNCTION_FAIL("unexpected failure\n"
                               "n = %wd, type = %wd, p = %wu, nlevels = %wd\n",
                               n, type, p, nlevels);

        nmod_poly_clear(cp);
        nmod_poly_clear(cp_ref);
        nmod_mat_clear(A);
    }

    TEST_FUNCTION_END(state);
}

#undef _NTYPES
