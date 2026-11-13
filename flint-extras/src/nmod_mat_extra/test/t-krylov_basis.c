/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/test_helpers.h>
#include <flint/ulong_extras.h>
#include <flint/nmod.h>
#include <flint/nmod_mat.h>
#include "nmod_mat_extra/krylov_impl.h"

/* random Jordan matrix description with few distinct eigenvalues */
static slong
_rand_jordan(ulong * eig, slong * sizes, slong n, int kind,
        flint_rand_t state, ulong p)
{
    slong nb = 0, left = n;
    ulong vals[3];
    slong t;

    for (t = 0; t < 3; t++)
        vals[t] = n_randint(state, p);

    while (left > 0)
    {
        slong sz;

        if (kind == 2)  /* diagonal */
            sz = 1;
        else
        {
            slong lim = 1 + n_randint(state, 6);
            lim = FLINT_MIN(left, lim);
            sz = 1 + n_randint(state, lim);
        }

        sizes[nb] = sz;
        if (kind == 1)  /* nilpotent */
            eig[nb] = 0;
        else
            eig[nb] = vals[n_randint(state, 3)];
        left -= sz;
        nb++;
    }

    return nb;
}

/* dense matrix of the Jordan description */
static void
_jordan_to_dense(nmod_mat_t J, const ulong * eig, const slong * sizes, slong nb)
{
    slong b, c, o = 0;

    nmod_mat_zero(J);
    for (b = 0; b < nb; b++)
    {
        for (c = 0; c < sizes[b]; c++)
        {
            nmod_mat_entry(J, o + c, o + c) = (eig == NULL) ? 0 : eig[b];
            if (c + 1 < sizes[b])
                nmod_mat_entry(J, o + c, o + c + 1) = 1;
        }
        o += sizes[b];
    }
}

/* random square M with various structures */
static void
_rand_mult_matrix(nmod_mat_t M, flint_rand_t state)
{
    slong n = M->r;
    ulong p = M->mod.n;

    switch (n_randint(state, 6))
    {
        case 0:
            nmod_mat_randtest(M, state);
            break;
        case 1:
            nmod_mat_zero(M);
            break;
        case 2:
        case 3:
        {
            ulong * eig = flint_malloc((n + 1) * sizeof(ulong));
            slong * sizes = flint_malloc((n + 1) * sizeof(slong));
            slong nb = _rand_jordan(eig, sizes, n, n_randint(state, 3), state, p);
            _jordan_to_dense(M, eig, sizes, nb);
            flint_free(eig);
            flint_free(sizes);

            if (n_randint(state, 2) && n > 0)
            {
                /* similarity transform by a random invertible matrix */
                nmod_mat_t U, Ui, T;
                nmod_mat_init(U, n, n, p);
                nmod_mat_init(Ui, n, n, p);
                nmod_mat_init(T, n, n, p);
                nmod_mat_randrank(U, state, n);
                nmod_mat_randops(U, state, 2 * n);
                nmod_mat_inv(Ui, U);
                nmod_mat_mul(T, U, M);
                nmod_mat_mul(M, T, Ui);
                nmod_mat_clear(U);
                nmod_mat_clear(Ui);
                nmod_mat_clear(T);
            }
            break;
        }
        case 4:
            nmod_mat_one(M);
            break;
        default:
            nmod_mat_randfull(M, state);
            break;
    }
}

/* random E with various structures */
static void
_rand_input(nmod_mat_t E, flint_rand_t state)
{
    slong m = E->r, n = E->c, i;

    switch (n_randint(state, 4))
    {
        case 0:
            nmod_mat_randtest(E, state);
            break;
        case 1:
            if (m > 0 && n > 0)
            {
                nmod_mat_randrank(E, state, n_randint(state, FLINT_MIN(m, n) + 1));
                nmod_mat_randops(E, state, n_randint(state, 2 * m * n + 1));
            }
            break;
        case 2:
            nmod_mat_zero(E);
            for (i = 0; i < m && n > 0; i++)
                nmod_mat_entry(E, i, n_randint(state, n)) = 1;
            break;
        default:
            nmod_mat_randfull(E, state);
            if (m > 1)
                for (i = 0; i < n; i++)  /* duplicate a row */
                    nmod_mat_entry(E, m - 1, i) = nmod_mat_entry(E, 0, i);
            break;
    }
}

static void
_rand_shifts(slong * s, slong m, slong n, flint_rand_t state)
{
    slong i;
    int mode = n_randint(state, 5);

    for (i = 0; i < m; i++)
    {
        if (mode == 0)
            s[i] = 0;
        else if (mode == 1)
            s[i] = n_randint(state, 5) - 2;
        else if (mode == 2)
            s[i] = i * (n + 1);   /* "Hermite" shift */
        else if (mode == 3)
            s[i] = (m - i) * (n + 1);
        else
            s[i] = (slong) n_randint(state, 2000000) - 1000000;
    }
}

TEST_FUNCTION_START(nmod_mat_krylov_basis, state)
{
    slong iter;

    for (iter = 0; iter < 2000 * flint_test_multiplier(); iter++)
    {
        slong m, n, i, k, r, r2, N;
        ulong p;
        nmod_mat_t E, M, J, B, B2, K;
        slong * shifts, * degs, * rows, * exps, * rows2, * exps2, * deltas;
        slong * sizes;
        ulong * eig;
        slong nb;
        int ok, use_shifts, use_jordan;
        _nmod_mat_krylov_op_t op, op2;

        m = n_randint(state, 10);
        n = n_randint(state, 10);
        if (n_randint(state, 8) == 0)
        {
            m = n_randint(state, 30);
            n = n_randint(state, 30);
        }

        if (n_randint(state, 3) == 0)
            p = n_nth_prime(1 + n_randint(state, 4));   /* 2, 3, 5, 7 */
        else
            p = n_randtest_prime(state, 0);

        nmod_mat_init(E, m, n, p);
        nmod_mat_init(M, n, n, p);
        nmod_mat_init(J, n, n, p);
        nmod_mat_init(B, n, n, p);
        nmod_mat_init(B2, n, n, p);

        shifts = flint_malloc((m + 1) * sizeof(slong));
        degs = flint_malloc((m + 1) * sizeof(slong));
        deltas = flint_malloc((m + 1) * sizeof(slong));
        rows = flint_malloc((n + 1) * sizeof(slong));
        exps = flint_malloc((n + 1) * sizeof(slong));
        rows2 = flint_malloc((m * (n + 4) + 1) * sizeof(slong));
        exps2 = flint_malloc((m * (n + 4) + 1) * sizeof(slong));
        sizes = flint_malloc((n + 1) * sizeof(slong));
        eig = flint_malloc((n + 1) * sizeof(ulong));

        _rand_input(E, state);
        _rand_shifts(shifts, m, n, state);
        use_shifts = n_randint(state, 4) != 0;

        use_jordan = n_randint(state, 2);
        if (use_jordan)
        {
            nb = _rand_jordan(eig, sizes, n, n_randint(state, 3), state, p);
            _jordan_to_dense(M, eig, sizes, nb);
            _nmod_mat_krylov_op_init_jordan(op, eig, sizes, nb, M->mod);
        }
        else
        {
            nb = 0;
            _rand_mult_matrix(M, state);
            _nmod_mat_krylov_op_init_dense(op, M);
        }

        /* reference: row rank profile of the full Krylov matrix */
        r2 = _nmod_mat_krylov_basis_naive(B2, rows2, exps2, E, op,
                            use_shifts ? shifts : NULL, NULL);

        /* default degree bounds */
        for (k = 0; k < 5; k++)
        {
            slong chunk = (k == 0) ? 0 : (k == 1) ? 1 : (k == 2) ? 2 : 1 + n_randint(state, (k == 3) ? 10 : 100);

            if (k % 2 == 0)
            {
                r = _nmod_mat_krylov_basis_iter(B, rows, exps, deltas, NULL, NULL, NULL, &ok,
                        E, op, use_shifts ? shifts : NULL, NULL, chunk);
            }
            else
            {
                /* with the transformation G: G * B[:, piv] = I */
                nmod_mat_t G, C, P;
                slong * piv = flint_malloc((n + 1) * sizeof(slong));
                slong u, v;

                nmod_mat_t T, vv, ww, Mp;
                slong c;

                nmod_mat_init(G, n, n, p);
                nmod_mat_init(T, m, n, p);
                r = _nmod_mat_krylov_basis_iter(B, rows, exps, deltas, T, G, piv, &ok,
                        E, op, use_shifts ? shifts : NULL, NULL, chunk);

                /* target rows: T_i = E_i M^deltas[i] */
                nmod_mat_init(vv, 1, n, p);
                nmod_mat_init(ww, 1, n, p);
                nmod_mat_init(Mp, n, n, p);
                for (u = 0; u < m; u++)
                {
                    for (c = 0; c < n; c++)
                        nmod_mat_entry(vv, 0, c) = nmod_mat_entry(E, u, c);
                    nmod_mat_pow(Mp, M, deltas[u]);
                    nmod_mat_mul(ww, vv, Mp);
                    for (c = 0; c < n; c++)
                        if (nmod_mat_entry(ww, 0, c) != nmod_mat_entry(T, u, c))
                            TEST_FUNCTION_FAIL("target rows (chunk %wd)\n", chunk);
                }
                nmod_mat_clear(vv);
                nmod_mat_clear(ww);
                nmod_mat_clear(Mp);
                nmod_mat_clear(T);
                nmod_mat_init(C, r, r, p);
                nmod_mat_init(P, r, r, p);
                for (u = 0; u < r; u++)
                    for (v = 0; v < r; v++)
                        nmod_mat_entry(C, u, v) = nmod_mat_entry(B, u, piv[v]);
                {
                    nmod_mat_t Gw;
                    nmod_mat_window_init(Gw, G, 0, 0, r, r);
                    nmod_mat_mul(P, Gw, C);
                    nmod_mat_window_clear(Gw);
                }
                if (!nmod_mat_is_one(P))
                    TEST_FUNCTION_FAIL("transformation G (chunk %wd)\np = %wu, E = %{nmod_mat}\nM = %{nmod_mat}\n",
                            chunk, p, E, M);
                nmod_mat_clear(G);
                nmod_mat_clear(C);
                nmod_mat_clear(P);
                flint_free(piv);
            }

            if (r != r2 || !ok)
                TEST_FUNCTION_FAIL("rank (r = %wd, naive %wd, ok = %d, chunk %wd)\n"
                    "p = %wu, E = %{nmod_mat}\nM = %{nmod_mat}\n", r, r2, ok, chunk, p, E, M);

            for (i = 0; i < r; i++)
                if (rows[i] != rows2[i] || exps[i] != exps2[i])
                    TEST_FUNCTION_FAIL("profile (chunk %wd)\np = %wu, E = %{nmod_mat}\nM = %{nmod_mat}\n",
                            chunk, p, E, M);

            {
                nmod_mat_t W1, W2;
                nmod_mat_window_init(W1, B, 0, 0, r, n);
                nmod_mat_window_init(W2, B2, 0, 0, r, n);
                if (!nmod_mat_equal(W1, W2))
                    TEST_FUNCTION_FAIL("basis (chunk %wd)\n", chunk);
                nmod_mat_window_clear(W1);
                nmod_mat_window_clear(W2);
            }

            /* minimal degrees */
            for (i = 0; i < m; i++)
            {
                slong dl = 0, t;
                for (t = 0; t < r; t++)
                    if (rows[t] == i)
                        dl = FLINT_MAX(dl, exps[t] + 1);
                if (deltas[i] != dl)
                    TEST_FUNCTION_FAIL("deltas\n");
            }
        }

        /* Keller-Gehrig-type doubling algorithm */
        r = _nmod_mat_krylov_basis_doubling(B, rows, exps, deltas, &ok,
                E, op, use_shifts ? shifts : NULL, NULL);
        if (r != r2 || !ok)
            TEST_FUNCTION_FAIL("doubling: rank (r = %wd, naive %wd, ok = %d)\n"
                "p = %wu, E = %{nmod_mat}\nM = %{nmod_mat}\n", r, r2, ok, p, E, M);
        for (i = 0; i < r; i++)
            if (rows[i] != rows2[i] || exps[i] != exps2[i])
                TEST_FUNCTION_FAIL("doubling: profile\np = %wu, E = %{nmod_mat}\nM = %{nmod_mat}\n", p, E, M);
        {
            nmod_mat_t W1, W2;
            nmod_mat_window_init(W1, B, 0, 0, r, n);
            nmod_mat_window_init(W2, B2, 0, 0, r, n);
            if (!nmod_mat_equal(W1, W2))
                TEST_FUNCTION_FAIL("doubling: basis\n");
            nmod_mat_window_clear(W1);
            nmod_mat_window_clear(W2);
        }

        /* dense operator from the Jordan description gives the same */
        if (use_jordan)
        {
            _nmod_mat_krylov_op_init_dense(op2, M);
            r = _nmod_mat_krylov_basis_iter(B, rows, exps, NULL, NULL, NULL, NULL, &ok,
                    E, op2, use_shifts ? shifts : NULL, NULL, 0);
            if (r != r2 || !ok)
                TEST_FUNCTION_FAIL("jordan vs dense: rank\n");
            for (i = 0; i < r; i++)
                if (rows[i] != rows2[i] || exps[i] != exps2[i])
                    TEST_FUNCTION_FAIL("jordan vs dense: profile\n");
            _nmod_mat_krylov_op_clear(op2);
        }

        /* random degree bounds: when declared valid, the profile is that
           of the Krylov matrix with these bounds */
        for (i = 0; i < m; i++)
            degs[i] = n_randint(state, n + 3);

        r = _nmod_mat_krylov_basis_iter(B, rows, exps, deltas, NULL, NULL, NULL, &ok,
                E, op, use_shifts ? shifts : NULL, degs, 1 + n_randint(state, 8));

        {
            /* validity: for each i, the rows (i, j <= min(d_i, n)) of the
               full profile do not reach d_i */
            int valid = 1;
            for (i = 0; i < r2; i++)
                if (exps2[i] >= FLINT_MIN(degs[rows2[i]], n))
                    valid = 0;

            if (ok != valid)
                TEST_FUNCTION_FAIL("ok flag (ok = %d, valid = %d)\n", ok, valid);

            if (ok)
            {
                slong r3;
                nmod_mat_t B3;
                slong * rows3 = flint_malloc((n + 1) * sizeof(slong));
                slong * exps3 = flint_malloc((n + 1) * sizeof(slong));

                nmod_mat_init(B3, n, n, p);
                r3 = _nmod_mat_krylov_basis_naive(B3, rows3, exps3, E, op,
                            use_shifts ? shifts : NULL, degs);

                if (r3 != r || r != r2)
                    TEST_FUNCTION_FAIL("degree bounds: rank\n");
                for (i = 0; i < r; i++)
                    if (rows[i] != rows3[i] || exps[i] != exps3[i])
                        TEST_FUNCTION_FAIL("degree bounds: profile\n");

                /* the doubling algorithm with the same bounds */
                r3 = _nmod_mat_krylov_basis_doubling(B3, rows3, exps3, NULL, &ok,
                        E, op, use_shifts ? shifts : NULL, degs);
                if (r3 != r || !ok)
                    TEST_FUNCTION_FAIL("doubling, degree bounds: rank\n");
                for (i = 0; i < r; i++)
                    if (rows[i] != rows3[i] || exps[i] != exps3[i])
                        TEST_FUNCTION_FAIL("doubling, degree bounds: profile\n");

                nmod_mat_clear(B3);
                flint_free(rows3);
                flint_free(exps3);
            }
            else
            {
                r = _nmod_mat_krylov_basis_doubling(B, rows, exps, NULL, &ok,
                        E, op, use_shifts ? shifts : NULL, degs);
                if (ok)
                    TEST_FUNCTION_FAIL("doubling: invalid degree bounds not detected\n");
            }
        }

        /* extreme shifts give the same result as equivalent small ones:
           a low group with shifts near WORD_MIN and a high group with
           shifts near WORD_MAX, versus the same groups 1000 apart */
        {
            slong * sm = flint_malloc((m + 1) * sizeof(slong));
            slong * sx = flint_malloc((m + 1) * sizeof(slong));
            slong * rows3 = flint_malloc((n + 1) * sizeof(slong));
            slong * exps3 = flint_malloc((n + 1) * sizeof(slong));
            nmod_mat_t B3;
            slong r3, r4;

            for (i = 0; i < m; i++)
            {
                slong small = (slong) n_randint(state, 7) - 3;
                int high = n_randint(state, 2);
                sm[i] = small + (high ? 1000 : 0);
                sx[i] = high ? WORD_MAX - 3 + small : WORD_MIN + 3 + small;
            }

            nmod_mat_init(B3, n, n, p);
            r3 = _nmod_mat_krylov_basis_naive(B3, rows3, exps3, E, op, sm, NULL);

            r4 = _nmod_mat_krylov_basis_iter(B, rows, exps, NULL, NULL, NULL, NULL, &ok,
                    E, op, sx, NULL, 1 + n_randint(state, 20));
            if (r4 != r3 || !ok)
                TEST_FUNCTION_FAIL("extreme shifts: rank\n");
            for (i = 0; i < r3; i++)
                if (rows[i] != rows3[i] || exps[i] != exps3[i])
                    TEST_FUNCTION_FAIL("extreme shifts: profile\n");

            r4 = _nmod_mat_krylov_basis_doubling(B, rows, exps, NULL, &ok,
                    E, op, sx, NULL);
            if (r4 != r3 || !ok)
                TEST_FUNCTION_FAIL("extreme shifts, doubling: rank\n");
            for (i = 0; i < r3; i++)
                if (rows[i] != rows3[i] || exps[i] != exps3[i])
                    TEST_FUNCTION_FAIL("extreme shifts, doubling: profile\n");

            r4 = _nmod_mat_krylov_basis_naive(B, rows, exps, E, op, sx, NULL);
            if (r4 != r3)
                TEST_FUNCTION_FAIL("extreme shifts, naive: rank\n");
            for (i = 0; i < r3; i++)
                if (rows[i] != rows3[i] || exps[i] != exps3[i])
                    TEST_FUNCTION_FAIL("extreme shifts, naive: profile\n");

            nmod_mat_clear(B3);
            flint_free(sm);
            flint_free(sx);
            flint_free(rows3);
            flint_free(exps3);
        }

        /* public interfaces */
        if (use_jordan)
            r = nmod_mat_krylov_basis_jordan(B, rows, exps, E, eig, sizes, nb,
                    use_shifts ? shifts : NULL, NULL);
        else
            r = nmod_mat_krylov_basis(B, rows, exps, E, M,
                    use_shifts ? shifts : NULL, NULL);
        if (r != r2)
            TEST_FUNCTION_FAIL("public interface\n");
        for (i = 0; i < r; i++)
            if (rows[i] != rows2[i] || exps[i] != exps2[i])
                TEST_FUNCTION_FAIL("public interface: profile\n");

        /* Krylov matrix: rows are E_i M^j, sorted by priority */
        for (i = 0; i < m; i++)
            degs[i] = n_randint(state, 4);
        N = _nmod_mat_krylov_nrows(m, degs);
        nmod_mat_init(K, N, n, p);
        if (use_jordan)
            nmod_mat_krylov_matrix_jordan(K, E, eig, sizes, nb,
                    use_shifts ? shifts : NULL, degs);
        else
            nmod_mat_krylov_matrix(K, E, M, use_shifts ? shifts : NULL, degs);

        _nmod_mat_krylov_order(rows2, exps2, m, use_shifts ? shifts : NULL, degs);
        for (k = 0; k < N; k++)
        {
            nmod_mat_t v, w, Mp;
            slong s1, s2;

            if (k > 0)
            {
                s1 = (use_shifts ? shifts[rows2[k - 1]] : 0) + exps2[k - 1];
                s2 = (use_shifts ? shifts[rows2[k]] : 0) + exps2[k];
                if (s1 > s2 || (s1 == s2 && rows2[k - 1] >= rows2[k]))
                    TEST_FUNCTION_FAIL("order\n");
            }

            nmod_mat_init(v, 1, n, p);
            nmod_mat_init(w, 1, n, p);
            nmod_mat_init(Mp, n, n, p);
            for (i = 0; i < n; i++)
                nmod_mat_entry(v, 0, i) = nmod_mat_entry(E, rows2[k], i);
            nmod_mat_pow(Mp, M, exps2[k]);
            nmod_mat_mul(w, v, Mp);
            for (i = 0; i < n; i++)
                if (nmod_mat_entry(w, 0, i) != nmod_mat_entry(K, k, i))
                    TEST_FUNCTION_FAIL("krylov_matrix\n");
            nmod_mat_clear(v);
            nmod_mat_clear(w);
            nmod_mat_clear(Mp);
        }
        nmod_mat_clear(K);

        _nmod_mat_krylov_op_clear(op);
        nmod_mat_clear(E);
        nmod_mat_clear(M);
        nmod_mat_clear(J);
        nmod_mat_clear(B);
        nmod_mat_clear(B2);
        flint_free(shifts);
        flint_free(degs);
        flint_free(deltas);
        flint_free(rows);
        flint_free(exps);
        flint_free(rows2);
        flint_free(exps2);
        flint_free(sizes);
        flint_free(eig);
    }

    TEST_FUNCTION_END(state);
}
