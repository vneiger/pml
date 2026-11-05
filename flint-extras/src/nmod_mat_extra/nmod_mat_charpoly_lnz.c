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
#include <flint/nmod.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>

#include "nmod_mat_extra.h"
#include "nmod_poly_mat_determinant.h"

/* Block Krylov reduction of x I - A to a k x k polynomial matrix.
 *
 * Write n = q*k + r with 0 <= r < k, and delta_i = q+1 for i < r, delta_i =
 * q for r <= i < k (so that sum_i delta_i = n and max_i delta_i =
 * ceil(n/k)). The Krylov matrix of the first k unit vectors e_i with indices
 * delta_i is the n x n matrix
 *     K = [E | A E | ... | A^(q-1) E | A^q E_r]
 * where E = [I_k ; 0] is n x k and E_r is made of the first r columns of E;
 * the column A^j e_i is at position j*k + i. If K is invertible, the
 * vectors A^(delta_i) e_i are uniquely written in the basis K:
 *     Z = K^{-1} [A^(delta_0) e_0 | ... | A^(delta_(k-1)) e_(k-1)],
 * so that K^{-1} A K is block companion, and det(x I - A) = det(P) with
 *     P[l, i] = x^(delta_i) [l == i] - sum_j Z[j*k + l, i] x^j.
 * Indeed, x I - K^{-1} A K is equivalent to diag(I_(n-k), P) by unimodular
 * column operations, and both determinants are monic (P is column reduced
 * with column degrees delta_i and leading column matrix the identity).
 *
 * Since K = [[I_k, X], [0, Y]] with Y of dimension n-k, K is invertible if
 * and only if Y is, and then Z = [G_t - X Y^{-1} G_b ; Y^{-1} G_b] where G =
 * [G_t ; G_b] gathers the vectors A^(delta_i) e_i.
 *
 * Storage: V = [A E | A^2 E | ... | A^q E | A^(q+1) E_r] is n x n (blocks of
 * width k, except the last one of width r). Its first n-k columns are the
 * columns k to n-1 of K, and its last k columns are G, with the vectors
 * A^(delta_i) e_i ordered as i = r, ..., k-1, 0, ..., r-1.
 *
 * Cost: A E is free, the next blocks cost n^2 (n-k) operations in total;
 * then solving costs (n-k)^3/3 + (n-k)^2 k + (n-k) k^2 operations.
 *
 * Input: A square n x n, 1 <= k < n, P initialized k x k. Returns 1 and sets
 * P if Y is invertible, otherwise returns 0 (P undefined). */
static int _charpoly_krylov(nmod_poly_mat_t P, const nmod_mat_t A, slong k)
{
    const slong n = A->r;
    const slong q = n / k;
    const slong r = n % k;
    const slong m = n - k;

    nmod_mat_t V, Vin, Vout;
    nmod_mat_init(V, n, n, A->mod.n);

    /* first block: A E = A[:, 0:k] */
    nmod_mat_window_init(Vin, A, 0, 0, n, k);
    nmod_mat_window_init(Vout, V, 0, 0, n, k);
    nmod_mat_set(Vout, Vin);
    nmod_mat_window_clear(Vin);
    nmod_mat_window_clear(Vout);

    /* early exit: if K has at least two nontrivial blocks of columns, the
     * first one (rows k to n-1 of A E) must have full column rank k; this
     * quickly discards e.g. triangular, diagonal or companion matrices */
    int ok = 1;
    if (m > k)
    {
        nmod_mat_window_init(Vin, V, k, 0, n, k);
        ok = (nmod_mat_rank(Vin) == k);
        nmod_mat_window_clear(Vin);
    }

    /* blocks A^2 E, ..., A^q E, and A^(q+1) E_r */
    for (slong j = 1; ok && j <= q; j++)
    {
        const slong w = (j < q) ? k : r;
        if (w == 0)
            break;
        nmod_mat_window_init(Vin, V, 0, (j - 1) * k, n, (j - 1) * k + w);
        nmod_mat_window_init(Vout, V, 0, j * k, n, j * k + w);
        nmod_mat_mul(Vout, A, Vin);
        nmod_mat_window_clear(Vin);
        nmod_mat_window_clear(Vout);
    }

    if (ok)
    {
        /* Z_b = Y^{-1} G_b, Z_t = G_t - X Z_b */
        nmod_mat_t X, Y, Gt, Gb, Z, Zt, Zb;
        nmod_mat_window_init(X, V, 0, 0, k, m);
        nmod_mat_window_init(Y, V, k, 0, n, m);
        nmod_mat_window_init(Gt, V, 0, m, k, n);
        nmod_mat_window_init(Gb, V, k, m, n, n);
        nmod_mat_init(Z, n, k, A->mod.n);
        nmod_mat_window_init(Zt, Z, 0, 0, k, k);
        nmod_mat_window_init(Zb, Z, k, 0, n, k);

        ok = nmod_mat_solve(Zb, Y, Gb);
        if (ok)
        {
            nmod_mat_mul(Zt, X, Zb);
            nmod_mat_sub(Zt, Gt, Zt);

            /* column c of Z corresponds to i = (c + r) mod k */
            for (slong c = 0; c < k; c++)
            {
                const slong i = (c + r < k) ? c + r : c + r - k;
                for (slong l = 0; l < k; l++)
                {
                    /* entry P[l, i]: coefficients j with j*k + l < n, plus
                     * the leading coefficient x^(delta_i) if l == i */
                    nmod_poly_struct * f = nmod_poly_mat_entry(P, l, i);
                    const slong lenl = (l < r) ? q + 1 : q;
                    const slong len = (l == i) ? lenl + 1 : lenl;
                    nmod_poly_fit_length(f, len);
                    for (slong j = 0; j < lenl; j++)
                        f->coeffs[j] = nmod_neg(nmod_mat_entry(Z, j * k + l, c), A->mod);
                    if (l == i)
                        f->coeffs[lenl] = 1;
                    _nmod_poly_set_length(f, len);
                    _nmod_poly_normalise(f);
                }
            }
        }

        nmod_mat_window_clear(X);
        nmod_mat_window_clear(Y);
        nmod_mat_window_clear(Gt);
        nmod_mat_window_clear(Gb);
        nmod_mat_window_clear(Zt);
        nmod_mat_window_clear(Zb);
        nmod_mat_clear(Z);
    }

    nmod_mat_clear(V);
    return ok;
}

int _nmod_mat_charpoly_lnz_block_krylov(nmod_poly_t cp, const nmod_mat_t A, slong nlevels)
{
    const slong n = A->r;

    if (A->r != A->c)
        flint_throw(FLINT_ERROR, "Exception (_nmod_mat_charpoly_lnz_block_krylov). Non-square matrix.\n");

    if (n == 0)
    {
        nmod_poly_one(cp);
        return 1;
    }

    nlevels = FLINT_MAX(1, FLINT_MIN(nlevels, n));
    const slong k = (n + nlevels - 1) / nlevels;

    /* no reduction: Labahn-Neiger-Zhou triangularization of x I - A (not
     * through nmod_poly_mat_determinant, whose degree-1 case would compute
     * the characteristic polynomial of A by nmod_mat_charpoly) */
    if (k == n)
    {
        nmod_poly_mat_t F;
        nmod_poly_mat_init(F, n, n, A->mod.n);
        for (slong i = 0; i < n; i++)
        {
            for (slong j = 0; j < n; j++)
                nmod_poly_set_coeff_ui(nmod_poly_mat_entry(F, i, j), 0,
                                       nmod_neg(nmod_mat_entry(A, i, j), A->mod));
            nmod_poly_set_coeff_ui(nmod_poly_mat_entry(F, i, i), 1, 1);
        }
        _nmod_poly_mat_det_lnz(cp, F, 1);
        nmod_poly_mat_clear(F);
        return 1;
    }

    nmod_poly_mat_t P;
    nmod_poly_mat_init(P, k, k, A->mod.n);
    const int ok = _charpoly_krylov(P, A, k);
    if (ok)
        nmod_poly_mat_determinant(cp, P);
    nmod_poly_mat_clear(P);
    return ok;
}
