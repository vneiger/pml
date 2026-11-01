/*
    Copyright (C) 2025 Gilles Villard
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
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>
#include <flint/ulong_extras.h>

#include "nmod_poly_mat_forms.h"
#include "nmod_poly_mat_inverse.h"
#include "nmod_poly_mat_kernel.h"
#include "nmod_poly_mat_linsolve.h"
#include "nmod_poly_mat_multiply.h"
#include "nmod_poly_mat_utils.h"

/*------------------------------------------------------------*/
/* power series solution, Dixon's x-adic lifting              */
/*------------------------------------------------------------*/

/* With C = A^{-1} mod x^order, the chunks X_k of X (coefficients of degree k
 * to k+order-1, for k = 0, order, 2*order, ...) are obtained from the
 * residues R_k = (B - X_{<k} A) div x^k, where X_{<k} is the sum of the
 * chunks before k:
 *   - short product: X_k = R_k C mod x^order (one low product);
 *   - middle product: R_{k+order} = (R_k - X_k A) div x^order
 *   (note that X_k A = R_k mod x^order).
 * The residue is not stored as such, since it would have the length of B;
 * instead, R_k = (B div x^k) + E_k with E_k = -(X_{<k} A) div x^k, which has
 * length less than length(A) since X_{<k} has length at most k:
 *   - X_k = (W + E_k) C mod x^order, with W the window B[k, k+order);
 *   - E_{k+order} = (E_k div x^order) - (X_k A) div x^order.
 * Then X A = B mod x^prec for X = sum_k x^k X_k. At each iteration, only the
 * coefficients needed for the remaining precision rem = prec - k are
 * computed, and the last iteration computes X_k mod x^rem and no residue. */
int nmod_poly_mat_solve_series_dixon(nmod_poly_mat_t X,
                                     const nmod_poly_mat_t A,
                                     const nmod_poly_mat_t B,
                                     slong order,
                                     slong prec)
{
    const slong m = A->r;
    const slong n = B->r;

    if (prec <= 0 || m == 0)
    {
        nmod_poly_mat_zero(X);
        return 1;
    }

    /* B is read at each iteration, and A until the end */
    if (X == A || X == B)
    {
        nmod_poly_mat_t T;
        nmod_poly_mat_init(T, n, m, A->modulus);
        const int ok = nmod_poly_mat_solve_series_dixon(T, A, B, order, prec);
        if (ok)
            nmod_poly_mat_swap_entrywise(X, T);
        nmod_poly_mat_clear(T);
        return ok;
    }

    const slong lenA = nmod_poly_mat_max_length(A);
    if (order <= 0)
        order = FLINT_MAX(lenA, 1);
    if (order > prec)
        order = prec;

    nmod_poly_mat_t C;
    nmod_poly_mat_init(C, m, m, A->modulus);
    if (! nmod_poly_mat_inv_series(C, A, order))
    {
        nmod_poly_mat_clear(C);
        return 0;
    }

    nmod_poly_mat_t E, W, S, T;
    nmod_poly_mat_init(E, n, m, A->modulus);
    nmod_poly_mat_init(W, n, m, A->modulus);
    nmod_poly_mat_init(S, n, m, A->modulus);
    nmod_poly_mat_init(T, n, m, A->modulus);

    /* Ash = (A div x) mod x^(prec-1): since S has length at most order, the
     * constant coefficient of A does not contribute to (S A) div x^order,
     * which is (S Ash) div x^(order-1); precomputing Ash also avoids a
     * shifted copy of A in each middle product */
    nmod_poly_mat_t Ash;
    nmod_poly_mat_init(Ash, m, m, A->modulus);
    nmod_poly_mat_shift_right_trunc(Ash, A, 1, prec - 1);

    nmod_poly_mat_zero(X);

    /* k = offset of the current chunk;
     * with k + order < prec whenever k is increased (so k < prec always) */
    for (slong k = 0; ; k += order)
    {
        const slong rem = prec - k;
        const slong p = FLINT_MIN(order, rem);

        /* S = (B[k, k+p) + E) C mod x^p, and X = X + x^k S (taking fast
         * copy path without arithmetic, since X has length <= k) */
        nmod_poly_mat_shift_right_trunc(W, B, k, p);
        nmod_poly_mat_add(W, W, E);
        nmod_poly_mat_mulmid(S, W, C, 0, p);
        nmod_poly_mat_add_shift_left(X, X, S, k);

        if (rem <= order)
            break;

        /* E = (E div x^order) - (S A) div x^order, mod x^(rem - order);
         * length of S A is <= lenA + order - 1 */
        const slong hi = (lenA - 1 < rem - order) ? order + lenA - 1 : rem;
        nmod_poly_mat_mulmid(T, S, Ash, order - 1, hi - 1);
        nmod_poly_mat_shift_right_trunc(E, E, order, rem - order);
        nmod_poly_mat_sub(E, E, T);
    }

    nmod_poly_mat_clear(C);
    nmod_poly_mat_clear(Ash);
    nmod_poly_mat_clear(E);
    nmod_poly_mat_clear(W);
    nmod_poly_mat_clear(S);
    nmod_poly_mat_clear(T);
    return 1;
}

/*------------------------------------------------------------*/
/* solution over the polynomials, via kernel basis            */
/*------------------------------------------------------------*/

int nmod_poly_mat_linsolve_via_kernel(nmod_poly_mat_t u,
                                      nmod_poly_t den,
                                      const nmod_poly_mat_t A,
                                      const nmod_poly_mat_t b)
{
    const slong m = A->r;
    const slong n = A->c;

    if (n == 0)  /* no equation */
    {
        nmod_poly_mat_zero(u);
        nmod_poly_one(den);
        return 1;
    }

    /* pmat = [A ; b], of size (m+1) x n */
    nmod_poly_mat_t pmat;
    nmod_poly_mat_init(pmat, m + 1, n, A->modulus);
    for (slong j = 0; j < n; j++)
    {
        for (slong i = 0; i < m; i++)
            nmod_poly_set(nmod_poly_mat_entry(pmat, i, j), nmod_poly_mat_entry(A, i, j));
        nmod_poly_set(nmod_poly_mat_entry(pmat, m, j), nmod_poly_mat_entry(b, 0, j));
    }

    /* left kernel basis, shift = row degree of pmat; in output, shift
     * contains the shifted row degree of the kernel basis */
    slong * shift = FLINT_ARRAY_ALLOC(m + 1, slong);
    slong * pivind = FLINT_ARRAY_ALLOC(m + 1, slong);
    nmod_poly_mat_row_degree_zero(shift, pmat, NULL);

    nmod_poly_mat_t ker;
    nmod_poly_mat_init(ker, m + 1, m + 1, A->modulus);
    const slong nullity = nmod_poly_mat_kernel_zls_approx(ker, pivind, shift, pmat);

    /* among the kernel basis rows with nonzero last entry, find one with
     * minimal shifted row degree: by the predictable degree property, this
     * gives a solution of minimal shifted degree among all solutions */
    slong best = -1;
    for (slong i = 0; i < nullity; i++)
        if (! nmod_poly_is_zero(nmod_poly_mat_entry(ker, i, m))
            && (best == -1 || shift[i] < shift[best]))
            best = i;

    if (best >= 0)
    {
        /* den = ker[best, m] made monic, u = -ker[best, 0:m] / lc */
        const nmod_poly_struct * d = nmod_poly_mat_entry(ker, best, m);
        const nmod_t mod = d->mod;
        const ulong inv = n_invmod(d->coeffs[d->length - 1], mod.n);
        const ulong neginv = nmod_neg(inv, mod);
        for (slong i = 0; i < m; i++)
            nmod_poly_scalar_mul_nmod(nmod_poly_mat_entry(u, 0, i),
                                      nmod_poly_mat_entry(ker, best, i), neginv);
        nmod_poly_scalar_mul_nmod(den, d, inv);
    }

    nmod_poly_mat_clear(ker);
    nmod_poly_mat_clear(pmat);
    flint_free(shift);
    flint_free(pivind);
    return (best >= 0);
}
