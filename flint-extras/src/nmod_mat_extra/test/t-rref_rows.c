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
#include <flint/nmod_vec.h>
#include <flint/nmod_mat.h>
#include "nmod_mat_extra/krylov_impl.h"

TEST_FUNCTION_START(nmod_mat_rref_rows, state)
{
    slong iter;

    for (iter = 0; iter < 2000 * flint_test_multiplier(); iter++)
    {
        slong k, n, lim, r, a, t, t2, prev, rk, i;
        ulong p;
        nmod_mat_t X, X0, S, Big, Xw;
        int win;
        slong * Q;
        char * indep;

        k = n_randint(state, 40);
        n = n_randint(state, 40);
        lim = n_randint(state, n + 1);

        switch (n_randint(state, 4))
        {
            case 0: p = n_nth_prime(1 + n_randint(state, 4)); break;
            case 1: p = n_randprime(state, 2 + n_randint(state, 29), 1); break;
            case 2: p = n_randprime(state, 31 + n_randint(state, 30), 1); break;
            default: p = n_randtest_prime(state, 0); break;
        }

        nmod_mat_init(X, k, n, p);
        Q = flint_malloc((k + 1) * sizeof(slong));
        indep = flint_malloc(k + 1);

        /* random matrix of random rank, with repeated rows */
        if (k > 0 && n > 0 && n_randint(state, 2))
        {
            nmod_mat_randrank(X, state, n_randint(state, FLINT_MIN(k, n) + 1));
            nmod_mat_randops(X, state, n_randint(state, 2 * k * n + 1));
        }
        else
            nmod_mat_randtest(X, state);
        if (k > 2 && n_randint(state, 2))
            for (i = 0; i < n; i++)
                nmod_mat_entry(X, k - 1, i) = nmod_mat_entry(X, 1, i);

        nmod_mat_init_set(X0, X);

        /* in place, or in a window of a larger matrix (as in the Krylov
           engine), whose other entries must not change */
        win = n_randint(state, 2);
        if (win)
        {
            nmod_mat_t Big0;
            nmod_mat_init(Big, k + 3, n + 5, p);
            nmod_mat_randtest(Big, state);
            nmod_mat_window_init(Xw, Big, 1, 2, k + 1, n + 2);
            nmod_mat_set(Xw, X);
            nmod_mat_init_set(Big0, Big);
            r = _nmod_mat_rref_rows(Xw, lim, Q, indep);
            nmod_mat_set(X, Xw);
            nmod_mat_window_clear(Xw);
            for (a = 0; a < k + 3; a++)
                for (i = 0; i < n + 5; i++)
                    if ((a < 1 || a >= k + 1 || i < 2 || i >= n + 2)
                            && nmod_mat_entry(Big, a, i) != nmod_mat_entry(Big0, a, i))
                        TEST_FUNCTION_FAIL("entries outside the window changed\n");
            nmod_mat_clear(Big);
            nmod_mat_clear(Big0);
        }
        else
        {
            r = _nmod_mat_rref_rows(X, lim, Q, indep);
        }

        /* independent rows: those increasing the rank of the leading
           rows, restricted to the first lim columns */
        prev = 0;
        t = 0;
        for (a = 0; a < k; a++)
        {
            nmod_mat_t Wa;

            /* nmod_mat_rank works on a copy */
            nmod_mat_window_init(Wa, X0, 0, 0, a + 1, lim);
            rk = nmod_mat_rank(Wa);
            nmod_mat_window_clear(Wa);
            if (indep[a] != (rk > prev))
                TEST_FUNCTION_FAIL("independence of row %wd\nX = %{nmod_mat}\n", a, X0);
            t += (indep[a] != 0);
            prev = rk;
        }
        if (r != prev || r != t)
            TEST_FUNCTION_FAIL("rank\n");

        /* reduced echelon form with pivots Q[0] < Q[1] < ... < lim */
        for (t = 0; t < r; t++)
        {
            if (Q[t] < 0 || Q[t] >= lim || (t > 0 && Q[t] <= Q[t - 1]))
                TEST_FUNCTION_FAIL("pivot out of range or not increasing\n");
            for (t2 = 0; t2 < r; t2++)
                if (nmod_mat_entry(X, t2, Q[t]) != (t == t2))
                    TEST_FUNCTION_FAIL("not in reduced echelon form\n");
        }

        /* same span as the independent input rows (all the columns) */
        if (r > 0)
        {
            nmod_mat_t Sb;

            nmod_mat_init(S, 2 * r, n, p);
            t = 0;
            for (a = 0; a < k; a++)
                if (indep[a])
                    _nmod_vec_set(nmod_mat_entry_ptr(S, t++, 0),
                                  nmod_mat_entry_ptr(X0, a, 0), n);
            for (t = 0; t < r; t++)
                _nmod_vec_set(nmod_mat_entry_ptr(S, r + t, 0),
                              nmod_mat_entry_ptr(X, t, 0), n);
            nmod_mat_window_init(Sb, S, r, 0, 2 * r, n);
            if (nmod_mat_rank(S) != r || nmod_mat_rank(Sb) != r)
                TEST_FUNCTION_FAIL("span\n");
            nmod_mat_window_clear(Sb);
            nmod_mat_clear(S);
        }

        nmod_mat_clear(X);
        nmod_mat_clear(X0);
        flint_free(Q);
        flint_free(indep);
    }

    TEST_FUNCTION_END(state);
}
