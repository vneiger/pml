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
#include <flint/nmod_poly.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_poly_mat.h>
#include "nmod_mat_extra/krylov_impl.h"
#include "nmod_poly_mat_forms.h"
#include "nmod_poly_mat_kernel.h"

/* whether every row p of P satisfies sum_c p_c(M) applied to E_c = 0 */
static int
_is_kernel(const nmod_poly_mat_t P, const nmod_mat_t E, const nmod_mat_t M)
{
    slong m = E->r, n = E->c, i, c, t, j;
    nmod_mat_t acc, w, tmp;
    int res = 1;

    nmod_mat_init(acc, 1, n, E->mod.n);
    nmod_mat_init(w, 1, n, E->mod.n);
    nmod_mat_init(tmp, 1, n, E->mod.n);

    for (i = 0; i < m && res; i++)
    {
        nmod_mat_zero(acc);
        for (c = 0; c < m; c++)
        {
            const nmod_poly_struct * p = nmod_poly_mat_entry(P, i, c);
            nmod_mat_zero(w);
            for (t = p->length - 1; t >= 0; t--)
            {
                nmod_mat_mul(tmp, w, M);
                for (j = 0; j < n; j++)
                    nmod_mat_entry(w, 0, j) = nmod_add(nmod_mat_entry(tmp, 0, j),
                        nmod_mul(p->coeffs[t], nmod_mat_entry(E, c, j), E->mod), E->mod);
            }
            nmod_mat_add(acc, acc, w);
        }
        res = nmod_mat_is_zero(acc);
    }

    nmod_mat_clear(acc);
    nmod_mat_clear(w);
    nmod_mat_clear(tmp);

    return res;
}

static void
_sage_examples(void)
{
    ulong Ev[9] = {27, 49, 29, 50, 58, 0, 77, 10, 29};
    slong sh[3][3] = {{0, 0, 0}, {0, 3, 6}, {3, 0, 2}};
    /* expected coefficients, low degree first, rows of P, -1 terminates */
    slong expected[3][9][5] = {
        {{82, 40, 1, -1}, {76, -1}, {-1}, {13, 3, -1}, {57, 1, -1}, {-1}, {96, -1}, {96, -1}, {1, -1}},
        {{0, 0, 0, 1, -1}, {-1}, {-1}, {70, 72, 60, -1}, {1, -1}, {-1}, {69, 72, 60, -1}, {-1}, {1, -1}},
        {{1, -1}, {79, 49, 26, -1}, {-1}, {-1}, {0, 0, 0, 1, -1}, {-1}, {-1}, {78, 49, 26, -1}, {1, -1}}};
    slong s, i, j, t, sz = 3;
    nmod_mat_t E, M;
    nmod_poly_mat_t P;

    nmod_mat_init(E, 3, 3, 97);
    nmod_mat_init(M, 3, 3, 97);
    nmod_poly_mat_init(P, 3, 3, 97);

    for (i = 0; i < 9; i++)
        nmod_mat_entry(E, i / 3, i % 3) = Ev[i];
    nmod_mat_entry(M, 0, 1) = 1;
    nmod_mat_entry(M, 1, 2) = 1;

    for (s = 0; s < 3; s++)
    {
        for (int jordan = 0; jordan < 2; jordan++)
        {
            int ok = jordan ? nmod_poly_mat_krylov_kernel_basis_jordan(P, E, NULL, &sz, 1, sh[s], NULL)
                            : nmod_poly_mat_krylov_kernel_basis(P, E, M, sh[s], NULL);
            if (!ok)
                flint_throw(FLINT_TEST_FAIL, "Sage example %wd: failure\n", s);

            for (i = 0; i < 3; i++)
            {
                for (j = 0; j < 3; j++)
                {
                    const nmod_poly_struct * p = nmod_poly_mat_entry(P, i, j);
                    for (t = 0; expected[s][3 * i + j][t] != -1; t++)
                        if (nmod_poly_get_coeff_ui(p, t) != (ulong) expected[s][3 * i + j][t])
                            flint_throw(FLINT_TEST_FAIL, "Sage example %wd\n", s);
                    if (p->length != t)
                        flint_throw(FLINT_TEST_FAIL, "Sage example %wd: length\n", s);
                }
            }
        }
    }

    nmod_mat_clear(E);
    nmod_mat_clear(M);
    nmod_poly_mat_clear(P);
}

TEST_FUNCTION_START(nmod_poly_mat_krylov_kernel_basis, state)
{
    slong iter;

    _sage_examples();

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        slong m, n, i, r, rnaive, sumdeg, nb, t, k;
        ulong p;
        nmod_mat_t E, M, B, R, Rel;
        nmod_poly_mat_t P, P2;
        slong * shifts, * rows, * exps, * deltas, * degs, * sizes;
        ulong * eig;
        int use_shifts, use_jordan, res;
        _nmod_mat_krylov_op_t op;

        m = n_randint(state, 8);
        n = n_randint(state, 10);
        if (n_randint(state, 10) == 0)
        {
            m = n_randint(state, 40);
            n = n_randint(state, 40);
        }

        if (n_randint(state, 3) == 0)
            p = n_nth_prime(1 + n_randint(state, 4));
        else
            p = n_randtest_prime(state, 0);

        nmod_mat_init(E, m, n, p);
        nmod_mat_init(M, n, n, p);
        nmod_mat_init(B, n, n, p);
        nmod_mat_init(R, m, n, p);
        nmod_poly_mat_init(P, m, m, p);
        nmod_poly_mat_init(P2, m, m, p);

        shifts = flint_malloc((m + 1) * sizeof(slong));
        degs = flint_malloc((m + 1) * sizeof(slong));
        deltas = flint_malloc((m + 1) * sizeof(slong));
        rows = flint_malloc((n + 1) * sizeof(slong));
        exps = flint_malloc((n + 1) * sizeof(slong));
        sizes = flint_malloc((n + 1) * sizeof(slong));
        eig = flint_malloc((n + 1) * sizeof(ulong));

        /* input: random matrix with a random rank, sometimes zero */
        if (m > 0 && n > 0 && n_randint(state, 4) != 0)
        {
            nmod_mat_randrank(E, state, n_randint(state, FLINT_MIN(m, n) + 1));
            nmod_mat_randops(E, state, n_randint(state, 2 * m * n + 1));
        }
        else
            nmod_mat_randtest(E, state);

        for (i = 0; i < m; i++)
        {
            switch (iter % 4)
            {
                case 0: shifts[i] = 0; break;
                case 1: shifts[i] = n_randint(state, 7) - 3; break;
                case 2: shifts[i] = i * (n + 1); break;
                default: shifts[i] = (slong) n_randint(state, 100) - 50;
            }
        }
        use_shifts = n_randint(state, 5) != 0;

        /* multiplication matrix: Jordan, possibly with a similarity
           transform, or random */
        use_jordan = n_randint(state, 2);
        {
            slong left = n;
            ulong vals[2];
            int kind = n_randint(state, 3);
            vals[0] = n_randint(state, p);
            vals[1] = n_randint(state, p);
            nb = 0;
            while (left > 0)
            {
                slong lim = 1 + n_randint(state, 5);
                lim = FLINT_MIN(left, lim);
                sizes[nb] = (kind == 2) ? 1 : 1 + n_randint(state, lim);
                eig[nb] = (kind == 1) ? 0 : vals[n_randint(state, 2)];
                left -= sizes[nb];
                nb++;
            }
            nmod_mat_zero(M);
            t = 0;
            for (k = 0; k < nb; k++)
            {
                for (i = 0; i < sizes[k]; i++)
                {
                    nmod_mat_entry(M, t + i, t + i) = eig[k];
                    if (i + 1 < sizes[k])
                        nmod_mat_entry(M, t + i, t + i + 1) = 1;
                }
                t += sizes[k];
            }
        }

        if (!use_jordan && n_randint(state, 2))
            nmod_mat_randtest(M, state);

        if (use_jordan)
        {
            res = nmod_poly_mat_krylov_kernel_basis_jordan(P, E, eig, sizes, nb,
                    use_shifts ? shifts : NULL, NULL);
            _nmod_mat_krylov_op_init_jordan(op, eig, sizes, nb, E->mod);
        }
        else
        {
            res = nmod_poly_mat_krylov_kernel_basis(P, E, M,
                    use_shifts ? shifts : NULL, NULL);
            _nmod_mat_krylov_op_init_dense(op, M);
        }

        if (!res)
            TEST_FUNCTION_FAIL("failure with default degree bounds\n");

        if (!nmod_poly_mat_is_popov(P, use_shifts ? shifts : NULL, ROW_LOWER))
            TEST_FUNCTION_FAIL("not in shifted Popov form\nE = %{nmod_mat}\nM = %{nmod_mat}\n", E, M);

        if (!_is_kernel(P, E, M))
            TEST_FUNCTION_FAIL("not in the kernel\nE = %{nmod_mat}\nM = %{nmod_mat}\n", E, M);

        /* the sum of the diagonal degrees is the rank of the Krylov matrix:
           the rows generate the whole kernel (JNSV 2017, Lemma 7.5) */
        rnaive = _nmod_mat_krylov_basis_naive(B, rows, exps, E, op,
                    use_shifts ? shifts : NULL, NULL);
        sumdeg = 0;
        for (i = 0; i < m; i++)
            sumdeg += nmod_poly_degree(nmod_poly_mat_entry(P, i, i));
        if (sumdeg != rnaive)
            TEST_FUNCTION_FAIL("sum of degrees %wd != rank %wd\n", sumdeg, rnaive);

        /* dense version of the Jordan matrix gives the same basis */
        if (use_jordan)
        {
            res = nmod_poly_mat_krylov_kernel_basis(P2, E, M,
                    use_shifts ? shifts : NULL, NULL);
            if (!res || !nmod_poly_mat_equal(P, P2))
                TEST_FUNCTION_FAIL("jordan vs dense\n");
        }

        /* compact form: R * B = rows E_i M^{delta_i} */
        if (use_jordan)
            r = nmod_mat_krylov_kernel_jordan(R, rows, exps, deltas, E, eig, sizes, nb,
                    use_shifts ? shifts : NULL, NULL);
        else
            r = nmod_mat_krylov_kernel(R, rows, exps, deltas, E, M,
                    use_shifts ? shifts : NULL, NULL);
        if (r != rnaive)
            TEST_FUNCTION_FAIL("compact kernel: rank\n");
        if (r > 0)
        {
            nmod_mat_t Rw, Bw, Mp, row, rowM;
            r = _nmod_mat_krylov_basis_naive(B, rows, exps, E, op,
                    use_shifts ? shifts : NULL, NULL);
            nmod_mat_init(Rel, m, n, p);
            nmod_mat_window_init(Rw, R, 0, 0, m, r);
            nmod_mat_window_init(Bw, B, 0, 0, r, n);
            nmod_mat_mul(Rel, Rw, Bw);
            nmod_mat_init(Mp, n, n, p);
            nmod_mat_init(row, 1, n, p);
            nmod_mat_init(rowM, 1, n, p);
            for (i = 0; i < m; i++)
            {
                for (t = 0; t < n; t++)
                    nmod_mat_entry(row, 0, t) = nmod_mat_entry(E, i, t);
                nmod_mat_pow(Mp, M, deltas[i]);
                nmod_mat_mul(rowM, row, Mp);
                for (t = 0; t < n; t++)
                    if (nmod_mat_entry(rowM, 0, t) != nmod_mat_entry(Rel, i, t))
                        TEST_FUNCTION_FAIL("compact kernel: relation\n");
            }
            nmod_mat_clear(Mp);
            nmod_mat_clear(row);
            nmod_mat_clear(rowM);
            nmod_mat_window_clear(Rw);
            nmod_mat_window_clear(Bw);
            nmod_mat_clear(Rel);
        }

        /* compact kernel with small chunks gives the same relations */
        if (r >= 0)
        {
            nmod_mat_t R2;
            slong * rows2 = flint_malloc((n + 1) * sizeof(slong));
            slong * exps2 = flint_malloc((n + 1) * sizeof(slong));
            slong * deltas2 = flint_malloc((m + 1) * sizeof(slong));
            slong r2, chunk = 1 + n_randint(state, 6);

            nmod_mat_init(R2, m, n, p);
            r2 = _nmod_mat_krylov_kernel(R2, rows2, exps2, deltas2, E, op,
                    use_shifts ? shifts : NULL, NULL, chunk);
            if (r2 != r)
                TEST_FUNCTION_FAIL("compact kernel, chunk %wd: rank\n", chunk);
            for (i = 0; i < m; i++)
            {
                if (deltas2[i] != deltas[i])
                    TEST_FUNCTION_FAIL("compact kernel, chunk %wd: deltas\n", chunk);
                for (t = 0; t < r; t++)
                    if (nmod_mat_entry(R2, i, t) != nmod_mat_entry(R, i, t))
                        TEST_FUNCTION_FAIL("compact kernel, chunk %wd: relations\n", chunk);
            }
            nmod_mat_clear(R2);
            flint_free(rows2);
            flint_free(exps2);
            flint_free(deltas2);
        }

        /* degree bounds: valid ones give the same basis, the function
           fails exactly when they are not valid */
        for (i = 0; i < m; i++)
            degs[i] = n_randint(state, n + 2);
        {
            int valid = 1;
            for (i = 0; i < rnaive; i++)
                if (exps[i] >= FLINT_MIN(degs[rows[i]], n))
                    valid = 0;

            if (use_jordan)
                res = nmod_poly_mat_krylov_kernel_basis_jordan(P2, E, eig, sizes, nb,
                        use_shifts ? shifts : NULL, degs);
            else
                res = nmod_poly_mat_krylov_kernel_basis(P2, E, M,
                        use_shifts ? shifts : NULL, degs);

            if (res != valid)
                TEST_FUNCTION_FAIL("degree bounds: res = %d, valid = %d\n", res, valid);
            if (res && !nmod_poly_mat_equal(P, P2))
                TEST_FUNCTION_FAIL("degree bounds: different basis\n");
        }

        _nmod_mat_krylov_op_clear(op);
        nmod_mat_clear(E);
        nmod_mat_clear(M);
        nmod_mat_clear(B);
        nmod_mat_clear(R);
        nmod_poly_mat_clear(P);
        nmod_poly_mat_clear(P2);
        flint_free(shifts);
        flint_free(degs);
        flint_free(deltas);
        flint_free(rows);
        flint_free(exps);
        flint_free(sizes);
        flint_free(eig);
    }

    TEST_FUNCTION_END(state);
}
