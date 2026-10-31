/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>  /* qsort */

#include <flint/flint.h>
#include <flint/nmod_poly_mat.h>
#include <flint/test_helpers.h>

#include "nmod_poly_mat_forms.h"
#include "nmod_poly_mat_kernel.h"
#include "nmod_poly_mat_row_basis.h"
#include "testing_collection.h"
#include "hermite_test_utils.h"

static int _rbtest_slong_cmp(const void * a, const void * b)
{
    const slong x = *(const slong *) a;
    const slong y = *(const slong *) b;
    return (x < y) ? -1 : (x > y);
}

/* checks that the first r rows of rbas form a row basis of A, that the other
 * rows of rbas are zero, that lfac (if not NULL) is a left factor, and if
 * shift is not NULL, that the row degree of the row basis is bounded (once
 * sorted) by the r largest entries of shift (sorted) */
static int _verify_row_basis(const nmod_poly_mat_t A,
                             const nmod_poly_mat_t rbas, slong r,
                             const nmod_poly_mat_t lfac,
                             const slong * shift,
                             const char * label)
{
    const slong m = A->r;
    const slong n = A->c;

    const slong rk = (m == 0 || n == 0) ? 0 : nmod_poly_mat_rank(A);
    if (rk != r)
    {
        flint_printf("FAIL (%s): rank %wd, expected %wd\n", label, r, rk);
        return 0;
    }

    for (slong i = r; i < rbas->r; i++)
        for (slong j = 0; j < rbas->c; j++)
            if (!nmod_poly_is_zero(nmod_poly_mat_entry(rbas, i, j)))
            {
                flint_printf("FAIL (%s): nonzero row beyond rank\n", label);
                return 0;
            }

    nmod_poly_mat_t R;
    nmod_poly_mat_window_init(R, rbas, 0, 0, r, n);

    int res = 1;
    if (r > 0 && nmod_poly_mat_rank(R) != r)
    {
        flint_printf("FAIL (%s): row basis does not have full row rank\n", label);
        res = 0;
    }

    if (res && !_htest_same_row_space(A, m, rbas, r))
    {
        flint_printf("FAIL (%s): row spaces differ\n", label);
        res = 0;
    }

    if (res && lfac != NULL)
    {
        for (slong i = 0; res && i < m; i++)
            for (slong j = r; res && j < lfac->c; j++)
                if (!nmod_poly_is_zero(nmod_poly_mat_entry(lfac, i, j)))
                {
                    flint_printf("FAIL (%s): nonzero column of left factor beyond rank\n", label);
                    res = 0;
                }
        if (res)
        {
            nmod_poly_mat_t C, prod;
            nmod_poly_mat_window_init(C, lfac, 0, 0, m, r);
            nmod_poly_mat_init(prod, m, n, A->modulus);
            if (r > 0)
                nmod_poly_mat_mul(prod, C, R);
            if (!nmod_poly_mat_equal(prod, A))
            {
                flint_printf("FAIL (%s): pmat != lfac * rbas\n", label);
                res = 0;
            }
            nmod_poly_mat_clear(prod);
            nmod_poly_mat_window_clear(C);
        }
    }

    if (res && shift != NULL && r > 0)
    {
        slong * s = flint_malloc(m * sizeof(slong));
        slong * rd = flint_malloc(r * sizeof(slong));
        for (slong i = 0; i < m; i++)
            s[i] = shift[i];
        qsort(s, m, sizeof(slong), _rbtest_slong_cmp);
        nmod_poly_mat_row_degree(rd, R, NULL);
        qsort(rd, r, sizeof(slong), _rbtest_slong_cmp);
        for (slong i = 0; res && i < r; i++)
            if (rd[i] > s[m - r + i])
            {
                flint_printf("FAIL (%s): row degree of row basis not bounded as expected\n", label);
                res = 0;
            }
        flint_free(s);
        flint_free(rd);
    }

    nmod_poly_mat_window_clear(R);

    if (!res)
    {
        flint_printf("input matrix (%wd x %wd, modulus %wu):\n", m, n, A->modulus);
        nmod_poly_mat_print_pretty(A, "x");
        flint_printf("\n");
    }
    return res;
}

static int _one_test_row_basis(const nmod_poly_mat_t A, flint_rand_t state)
{
    const slong m = A->r;
    const slong n = A->c;
    const slong mn = FLINT_MIN(m, n);
    int res = 1;

    nmod_poly_mat_t rbas, lfac;
    nmod_poly_mat_init(rbas, mn, n, A->modulus);
    nmod_poly_mat_init(lfac, m, mn, A->modulus);
    slong * s = flint_malloc(FLINT_MAX(1, m) * sizeof(slong));
    nmod_poly_mat_row_degree_zero(s, A, NULL);

    /* main interface for Zhou-Labahn algorithm */
    slong r = nmod_poly_mat_row_basis_zl(rbas, lfac, A);
    res = res && _verify_row_basis(A, rbas, r, lfac, s, "row_basis_zl");

    /* without the unimodularity certificate */
    if (res)
    {
        r = _nmod_poly_mat_row_basis_zl(rbas, lfac, A, NULL, NULL, NULL, 0, 0);
        res = _verify_row_basis(A, rbas, r, lfac, s, "row_basis_zl, no certificate");
    }

    /* with a larger shift and a precomputed kernel basis; the shift is
     * translated by a random constant (only differences between its entries
     * matter), possibly making it smaller than the row degree of A */
    if (res && m > 0 && n > 0)
    {
        slong * shift = flint_malloc(m * sizeof(slong));
        slong * nshift = flint_malloc(m * sizeof(slong));
        slong * kdeg = flint_malloc(m * sizeof(slong));
        slong * pivind = flint_malloc(m * sizeof(slong));
        slong * rdeg = flint_malloc(m * sizeof(slong));
        nmod_poly_mat_row_degree(rdeg, A, NULL);
        const slong cst = (n_randint(state, 2)) ? 0 : (slong) n_randint(state, 41) - 20;
        for (slong i = 0; i < m; i++)
        {
            shift[i] = s[i] + cst + ((n_randint(state, 2)) ? n_randint(state, 4) : 0);
            if (rdeg[i] < 0)
                shift[i] += (slong) n_randint(state, 7) - 3;
            kdeg[i] = shift[i];
        }
        /* normalized shift, for the degree bound of the verification */
        slong c = WORD_MAX;
        for (slong i = 0; i < m; i++)
            if (rdeg[i] >= 0 && shift[i] - rdeg[i] < c)
                c = shift[i] - rdeg[i];
        for (slong i = 0; i < m; i++)
            nshift[i] = (c == WORD_MAX) ? shift[i] : shift[i] - c;
        nmod_poly_mat_t ker, copy;
        nmod_poly_mat_init(ker, m, m, A->modulus);
        nmod_poly_mat_init_set(copy, A);
        const slong nullity = nmod_poly_mat_kernel_zls_approx(ker, pivind, kdeg, copy);
        const int check = n_randint(state, 2);
        r = _nmod_poly_mat_row_basis_zl(rbas, lfac, A, shift, ker, kdeg, nullity, check);
        res = _verify_row_basis(A, rbas, r, lfac, nshift, "row_basis_zl, shift and kernel given");
        if (res)
        {
            r = _nmod_poly_mat_row_basis_zl(rbas, NULL, A, shift, NULL, NULL, 0, check);
            res = _verify_row_basis(A, rbas, r, NULL, nshift, "row_basis_zl, shift given");
        }
        if (!res)
        {
            flint_printf("shift:");
            for (slong i = 0; i < m; i++)
                flint_printf(" %wd", shift[i]);
            flint_printf("\n");
        }
        nmod_poly_mat_clear(copy);
        nmod_poly_mat_clear(ker);
        flint_free(rdeg);
        flint_free(pivind);
        flint_free(kdeg);
        flint_free(nshift);
        flint_free(shift);
    }

    /* general interface */
    if (res)
    {
        r = nmod_poly_mat_row_basis(rbas, A);
        res = _verify_row_basis(A, rbas, r, NULL, s, "row_basis");
    }

    flint_free(s);
    nmod_poly_mat_clear(rbas);
    nmod_poly_mat_clear(lfac);
    return res;
}

/* random input for row basis: favours matrices with more rows than
 * columns whose row space is not the whole of K[x]^{1 x n} */
static void _rbtest_rand_input(nmod_poly_mat_t A, slong len, flint_rand_t state)
{
    const slong m = A->r;
    const slong n = A->c;
    if (m == 0 || n == 0)
        return;

    const ulong c = n_randint(state, 20);
    if (c < 2)
        gen_E(A, len, state);
    else if (c < 6)
        _htest_rand_rect(A, len, state);
    else
    {
        /* A = W * B with W random m x k, and B the first k rows of an n x n
         * matrix with (often) nontrivial invariant factors; k = min(m,n) in
         * most cases, smaller otherwise (rank-deficient) */
        slong k = FLINT_MIN(m, n);
        if (c >= 17)
            k = 1 + n_randint(state, k);
        nmod_poly_mat_t W, Bsq, B;
        nmod_poly_mat_init(W, m, k, A->modulus);
        nmod_poly_mat_init(Bsq, n, n, A->modulus);
        nmod_poly_mat_init(B, k, n, A->modulus);
        if (c == 16)
            nmod_poly_mat_rand(W, state, 1);  /* constant left factor */
        else
            nmod_poly_mat_rand(W, state, FLINT_MAX(1, len/2));
        _htest_rand_square(Bsq, FLINT_MAX(1, len), state);
        for (slong i = 0; i < k; i++)
            for (slong j = 0; j < n; j++)
                nmod_poly_set(nmod_poly_mat_entry(B, i, j), nmod_poly_mat_entry(Bsq, i, j));
        nmod_poly_mat_mul(A, W, B);
        nmod_poly_mat_clear(W);
        nmod_poly_mat_clear(Bsq);
        nmod_poly_mat_clear(B);
    }
}

TEST_FUNCTION_START(nmod_poly_mat_row_basis, state)
{
    int i, result;

    /* general shapes */
    for (i = 0; i < 800 * flint_test_multiplier(); i++)
    {
        const ulong prime = _htest_rand_prime(state);
        const slong m = (n_randint(state, 20) == 0) ? 0 : 1 + n_randint(state, 12);
        const slong n = (n_randint(state, 20) == 0) ? 0 : 1 + n_randint(state, 8);
        const slong len = (n_randint(state, 20) == 0) ? 0 : 1 + n_randint(state, 6);

        nmod_poly_mat_t A;
        nmod_poly_mat_init(A, m, n, prime);
        _rbtest_rand_input(A, len, state);

        result = _one_test_row_basis(A, state);
        nmod_poly_mat_clear(A);

        if (!result)
            TEST_FUNCTION_FAIL("prime = %wu, m = %wd, n = %wd, len = %wd\n", prime, m, n, len);
    }

    /* tall matrices, processed by blocks in the general interface */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        const ulong prime = _htest_rand_prime(state);
        const slong n = 1 + n_randint(state, 4);
        const slong m = n * (4 + n_randint(state, 5)) + n_randint(state, n);
        const slong len = 1 + n_randint(state, 5);

        nmod_poly_mat_t A, rbas;
        nmod_poly_mat_init(A, m, n, prime);
        nmod_poly_mat_init(rbas, n, n, prime);
        _rbtest_rand_input(A, len, state);

        slong * s = flint_malloc(m * sizeof(slong));
        nmod_poly_mat_row_degree_zero(s, A, NULL);
        slong r = nmod_poly_mat_row_basis(rbas, A);
        result = _verify_row_basis(A, rbas, r, NULL, s, "row_basis, tall");
        if (result)
        {
            const slong bsize = 1 + n_randint(state, n + 1);
            r = _nmod_poly_mat_row_basis_zl_blocks(rbas, A, bsize);
            result = _verify_row_basis(A, rbas, r, NULL, s, "row_basis_zl_blocks");
        }
        flint_free(s);
        nmod_poly_mat_clear(A);
        nmod_poly_mat_clear(rbas);

        if (!result)
            TEST_FUNCTION_FAIL("(tall) prime = %wu, m = %wd, n = %wd, len = %wd\n", prime, m, n, len);
    }

    /* generic matrices with more rows than columns and arbitrary degree
     * patterns, over a large field: when the row space is the whole of
     * K[x]^{1 x n}, the unimodularity certificate is expected to succeed
     * (this fails with negligible probability), in which case the identity
     * is returned */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        const ulong prime = n_randprime(state, 60, 1);
        const slong n = 1 + n_randint(state, 6);
        const slong m = n + 1 + n_randint(state, 4);

        nmod_poly_mat_t A, rbas;
        nmod_poly_mat_init(A, m, n, prime);
        nmod_poly_mat_init(rbas, n, n, prime);
        for (slong ii = 0; ii < m; ii++)
            for (slong jj = 0; jj < n; jj++)
                if (n_randint(state, 3))
                    nmod_poly_rand(nmod_poly_mat_entry(A, ii, jj), state, 1 + n_randint(state, 8));

        slong * s = flint_malloc(m * sizeof(slong));
        nmod_poly_mat_row_degree_zero(s, A, NULL);
        slong r = nmod_poly_mat_row_basis_zl(rbas, NULL, A);
        result = _verify_row_basis(A, rbas, r, NULL, s, "row_basis_zl, degree patterns");
        if (result && r == n)
        {
            /* is the row space the whole of K[x]^{1 x n}? */
            nmod_poly_mat_t H, I;
            nmod_poly_mat_init(H, m, n, prime);
            nmod_poly_mat_init(I, n, n, prime);
            nmod_poly_mat_one(I);
            _htest_upper_hnf(H, A);
            if (_htest_same_row_space(H, n, I, n) && !nmod_poly_mat_is_one(rbas))
            {
                flint_printf("FAIL: unimodularity certificate did not succeed on a generic input\n");
                nmod_poly_mat_print_pretty(A, "x");
                flint_printf("\n");
                result = 0;
            }
            nmod_poly_mat_clear(H);
            nmod_poly_mat_clear(I);
        }
        flint_free(s);
        nmod_poly_mat_clear(A);
        nmod_poly_mat_clear(rbas);

        if (!result)
            TEST_FUNCTION_FAIL("(degree patterns) prime = %wu, m = %wd, n = %wd\n", prime, m, n);
    }

    /* regression over GF(2): this input used to make the Newton iteration
     * call the vandermonde2 middle product, which requires an odd modulus */
    {
        nmod_poly_mat_t A, rbas, lfac;
        nmod_poly_mat_init(A, 6, 5, 2);
        nmod_poly_mat_init(rbas, 5, 5, 2);
        nmod_poly_mat_init(lfac, 6, 5, 2);
        nmod_poly_set_coeff_ui(nmod_poly_mat_entry(A, 0, 0), 64, 1);
        nmod_poly_set_coeff_ui(nmod_poly_mat_entry(A, 1, 1), 0, 1);
        nmod_poly_set_coeff_ui(nmod_poly_mat_entry(A, 2, 2), 0, 1);
        nmod_poly_set_coeff_ui(nmod_poly_mat_entry(A, 2, 4), 32, 1);
        nmod_poly_set_coeff_ui(nmod_poly_mat_entry(A, 3, 3), 0, 1);
        nmod_poly_set_coeff_ui(nmod_poly_mat_entry(A, 4, 4), 0, 1);
        nmod_poly_set_coeff_ui(nmod_poly_mat_entry(A, 5, 2), 96, 1);
        slong * s = flint_malloc(6 * sizeof(slong));
        nmod_poly_mat_row_degree_zero(s, A, NULL);
        slong r = nmod_poly_mat_row_basis_zl(rbas, lfac, A);
        result = _verify_row_basis(A, rbas, r, lfac, s, "row_basis_zl, GF(2) regression");
        flint_free(s);
        nmod_poly_mat_clear(A);
        nmod_poly_mat_clear(rbas);
        nmod_poly_mat_clear(lfac);
        if (!result)
            TEST_FUNCTION_FAIL("GF(2) regression\n");
    }

    TEST_FUNCTION_END(state);
}
