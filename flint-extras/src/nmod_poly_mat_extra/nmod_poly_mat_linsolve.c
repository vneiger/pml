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
#include <flint/nmod_mat.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>
#include <flint/nmod_vec.h>
#include <flint/perm.h>
#include <flint/ulong_extras.h>

#include "nmod_extra.h"       /* for nmod_find_root */
#include "nmod_poly_extra.h"  /* for NMOD_POLY_CAN_USE_GEOMETRIC */
#include "nmod_poly_mat_determinant.h"
#include "nmod_poly_mat_forms.h"
#include "nmod_poly_mat_inverse.h"
#include "nmod_poly_mat_kernel.h"
#include "nmod_poly_mat_linsolve.h"
#include "nmod_poly_mat_multiply.h"
#include "nmod_poly_mat_utils.h"
#include "impl.h"  /* for _pml_fft_small_single_transform */

/*------------------------------------------------------------*/
/* power series solution, Dixon's x-adic lifting              */
/*------------------------------------------------------------*/

/* Whether the products of Dixon's lifting use stored transforms of the
 * truncated inverse C and of A (see nmod_poly_mat_mulmid_precomp_right_init), for
 * an m x m matrix A of degree d = lenA - 1 >= 1 and n rows in B, and the
 * length N of the transforms of the low product, which fixes the default
 * order (about N - d, see nmod_poly_mat_solve_series_dixon).
 *
 * With stored transforms, an iteration costs about n m^2 np N pointwise
 * products, np being the number of primes of the transforms, plus the
 * transforms of the n m entries of the left operands, for N - d new
 * coefficients; the truncated inverse costs about m^omega M(N - d), once. For
 * one row, the best N is the smallest power of two N >= max(16, 2 d), so that
 * the order is about d; with more rows the transforms of C are reused more
 * and longer orders pay: N >= 2 d f with f = 2 up to 4 rows, 4 up to 16, 8
 * above. Without stored transforms, each product transforms C again (or
 * multiplies short polynomials), and the default order is 8 lenA.
 *
 * Products of short polynomials (as the middle products use for short
 * operands) are faster for small degrees and dimensions.
 *
 * 2026-11 Measured on an Intel Xeon (AVX-512, session VM), against the
 * previous version (no stored transforms, order 8 lenA), for one row and
 * prec = 2 m d + 1, m = 2..256, d = 1..4096, primes of 20 to 60 bits: with
 * the rule below the stored transforms are 1.5-6x faster from d^2 m of a few
 * thousands on (e.g. m = 64, d = 32: 0.18 s against 1.24 s for a 60-bit
 * prime; geometric mean of the ratios over the grid 1.6-2.0), and lose at
 * most 10-20% near the thresholds; for 2 to 64 rows, 1.2-2.2x faster with the
 * factor f above, except for degrees below 4, where the products of
 * several rows by C and A without stored transforms use the constant matrix
 * products of the Vandermonde routines and are up to 2.6x faster (e.g. 128
 * rows, m = 128, d = 2, 30-bit prime): stored transforms are not used there
 * for more than 4 rows; for B of length prec, 2-6x faster. The thresholds depend
 * on the number of transforms np (1 for p of at most about 20 bits, or for p
 * supporting a direct transform; 2 for p of at most about 44 bits; 3 above):
 * - np = 1: d^2 m >= 64 and d m >= 32;
 * - np = 2: d^2 m >= 256 and d m >= 64;
 * - np = 3: d^2 m >= 256 and d m >= 32.
 * \todo thresholds and factor measured on one machine only (tuning PR)
 *
 * precomp: 0 or 1 to force the choice (1: when stored transforms are
 * available, i.e. with PML_HAVE_MACHINE_VECTORS), -1 for the rule above. */
static int _dixon_use_precomp(ulong modn, slong m, slong n, slong lenA,
                              int precomp, slong * N)
{
#if PML_HAVE_MACHINE_VECTORS
    const slong d = lenA - 1;
    if (precomp == 0 || d <= 0)
        return 0;

    const slong f = (n <= 1) ? 1 : (n <= 4 ? 2 : (n <= 16 ? 4 : 8));
    *N = FLINT_MAX(16, WORD(1) << n_clog2(2 * d * f));
    if (precomp == 1)
        return 1;

    if ((ulong) d * m < 32 || (n > 4 && d < 4))
        return 0;
    const int single = _pml_fft_small_single_transform(modn, n_clog2(*N),
                                                       (ulong) m * FLINT_MIN(lenA, *N));
    const slong np = single ? 1 : (FLINT_BIT_COUNT(modn) <= 44 ? 2 : 3);
    const ulong d2m = (ulong) d * d * m;
    const ulong dm = (ulong) d * m;
    if (np == 1)
        return (d2m >= 64);
    if (np == 2)
        return (d2m >= 256 && dm >= 64);
    return (d2m >= 256);
#else
    (void) modn; (void) m; (void) n; (void) lenA; (void) precomp; (void) N;
    return 0;
#endif
}

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
int _nmod_poly_mat_solve_series_dixon(nmod_poly_mat_t X,
                                      const nmod_poly_mat_t A,
                                      const nmod_poly_mat_t B,
                                      slong order,
                                      slong prec,
                                      int precomp)
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
        const int ok = _nmod_poly_mat_solve_series_dixon(T, A, B, order, prec, precomp);
        if (ok)
            nmod_poly_mat_swap_entrywise(X, T);
        nmod_poly_mat_clear(T);
        return ok;
    }

    const slong lenA = nmod_poly_mat_max_length(A);
    const slong lenB = nmod_poly_mat_max_length(B);

    /* products with stored transforms or not, and default order: with
     * stored transforms, the order is about deg(A), so that the transforms
     * of the low product have length N (the left operand W has length at
     * most w, see below), and the truncated inverse is cheap; without, each
     * product transforms C again (or uses products of short polynomials)
     * and longer orders amortize this, see nmod_poly_mat_solve_series_dixon
     * in nmod_poly_mat_linsolve.h */
    slong N = 0;
    const int use_precomp = _dixon_use_precomp(A->modulus, m, n, lenA, precomp, &N);
    if (order <= 0)
    {
        if (use_precomp)
        {
            const slong w = FLINT_MAX(lenA - 1, FLINT_MIN(lenB, (N + 1) / 2));
            order = N + 1 - w;
        }
        else
            order = FLINT_MAX(8 * lenA, 1);
    }
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
    const slong lenAsh = nmod_poly_mat_max_length(Ash);

    /* the right operands C and Ash are the same in all iterations: their
     * transforms are computed once (see nmod_poly_mat_mulmid_precomp_right_init).
     * The left operands are W, of length at most lenW (the window of B, and
     * E, of length less than lenA), truncated at order, and S, of length at
     * most order. */
    const slong lenW = FLINT_MIN(order, FLINT_MAX(FLINT_MIN(order, lenB), lenA - 1));
    nmod_poly_mat_mulmid_precomp_right_t preC, preA;
    _nmod_poly_mat_mulmid_precomp_right_init(preC, C, FLINT_MAX(lenW, 1), 0, order,
                                             use_precomp);
    if (prec > order)
        _nmod_poly_mat_mulmid_precomp_right_init(preA, Ash, order, order - 1, order + lenAsh - 1,
                                                 use_precomp);

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
        if (lenA - 1 > p)
            nmod_poly_mat_truncate(W, p);
        nmod_poly_mat_mulmid_precomp_right(S, W, preC, p);
        nmod_poly_mat_add_shift_left(X, X, S, k);

        if (rem <= order)
            break;

        /* E = (E div x^order) - (S A) div x^order, mod x^(rem - order);
         * length of S A is <= lenA + order - 1 */
        const slong hi = (lenA - 1 < rem - order) ? order + lenA - 1 : rem;
        nmod_poly_mat_mulmid_precomp_right(T, S, preA, hi - 1);
        nmod_poly_mat_shift_right_trunc(E, E, order, rem - order);
        nmod_poly_mat_sub(E, E, T);
    }

    nmod_poly_mat_mulmid_precomp_right_clear(preC);
    if (prec > order)
        nmod_poly_mat_mulmid_precomp_right_clear(preA);
    nmod_poly_mat_clear(C);
    nmod_poly_mat_clear(Ash);
    nmod_poly_mat_clear(E);
    nmod_poly_mat_clear(W);
    nmod_poly_mat_clear(S);
    nmod_poly_mat_clear(T);
    return 1;
}

int nmod_poly_mat_solve_series_dixon(nmod_poly_mat_t X,
                                     const nmod_poly_mat_t A,
                                     const nmod_poly_mat_t B,
                                     slong order,
                                     slong prec)
{
    return _nmod_poly_mat_solve_series_dixon(X, A, B, order, prec, -1);
}

/*------------------------------------------------------------*/
/* power series solution, high-order lifting                  */
/*------------------------------------------------------------*/

/* For m x m A of degree at most e >= 1 with A(0) invertible, and inv = A^{-1}
 * mod x^e: replaces the slice S_s of A^{-1} (its coefficients of degree s e -
 * e + 1 to s e + e - 1, for some s >= 1) by S_{2s}. This is the row-wise
 * version of high_order_lift_inverse_odd in ntl-extras (mat_lzz_pX_linsolve):
 *   T = (S_s mod x^e) A mod x^e
 *   N = ((T S_s) div x^(e-1)) mod x^e
 *   M = ((N A) div x^e) mod x^(e-1)
 *   S_{2s} = N - x^e (M inv mod x^(e-1)) */
static void _hol_lift_slice(nmod_poly_mat_t slice, const nmod_poly_mat_t A,
                            const nmod_poly_mat_t inv, slong e)
{
    const slong m = A->r;
    nmod_poly_mat_t T, N, M;
    nmod_poly_mat_init(T, m, m, A->modulus);
    nmod_poly_mat_init(N, m, m, A->modulus);
    nmod_poly_mat_init(M, m, m, A->modulus);

    nmod_poly_mat_set(T, slice);
    nmod_poly_mat_truncate(T, e);
    nmod_poly_mat_mulmid(T, T, A, 0, e);
    nmod_poly_mat_mulmid(N, T, slice, e - 1, 2 * e - 1);
    nmod_poly_mat_mulmid(M, N, A, e, 2 * e - 1);
    nmod_poly_mat_mulmid(M, M, inv, 0, e - 1);
    nmod_poly_mat_sub_shift_left(slice, N, M, e);

    nmod_poly_mat_clear(T);
    nmod_poly_mat_clear(N);
    nmod_poly_mat_clear(M);
}

/* With e >= deg(A), the residues R_k of Dixon's lifting at order e, R_0 = B
 * (of degree less than e) and R_{k+1} = (R_k - Y_k A) div x^e with Y_k = R_k
 * A^{-1} mod x^e, satisfy R_{k+s} = (((R_k A^{-1}) div x^(s e)) mod x^e) A mod
 * x^e: with R_k A^{-1} = P + x^(s e) Q, P = R_k A^{-1} mod x^(s e), x^(s e)
 * R_{k+s} = R_k - P A by definition, so that R_{k+s} = Q A, which has degree
 * less than e. Since R_k has degree less than e, ((R_k A^{-1}) div x^(s e))
 * mod x^e = ((R_k S_s) div x^(e-1)) mod x^e, where S_s is the slice of
 * coefficients of degree s e - e + 1 to s e + e - 1 of A^{-1}. Starting from
 * R_0 and S_1, each step computes R_s, ..., R_{2s-1} from R_0, ..., R_{s-1}
 * (all residues are stacked in one matrix) and S_{2s} from S_s; then X = sum_k
 * x^(k e) Y_k mod x^prec. */
int nmod_poly_mat_solve_series_high_order_lifting(nmod_poly_mat_t X,
                                                  const nmod_poly_mat_t A,
                                                  const nmod_poly_mat_t B,
                                                  slong prec)
{
    const slong m = A->r;
    const slong n = B->r;
    const ulong p = A->modulus;

    if (prec <= 0 || m == 0)
    {
        nmod_poly_mat_zero(X);
        return 1;
    }

    if (X == A || X == B)
    {
        nmod_poly_mat_t T;
        nmod_poly_mat_init(T, n, m, p);
        const int ok = nmod_poly_mat_solve_series_high_order_lifting(T, A, B, prec);
        if (ok)
            nmod_poly_mat_swap_entrywise(X, T);
        nmod_poly_mat_clear(T);
        return ok;
    }

    /* only A mod x^prec and B mod x^prec matter */
    if (nmod_poly_mat_max_length(A) > prec || nmod_poly_mat_max_length(B) > prec)
    {
        nmod_poly_mat_t Ap, Bp;
        nmod_poly_mat_init_set(Ap, A);
        nmod_poly_mat_init_set(Bp, B);
        nmod_poly_mat_truncate(Ap, prec);
        nmod_poly_mat_truncate(Bp, prec);
        const int ok = nmod_poly_mat_solve_series_high_order_lifting(X, Ap, Bp, prec);
        nmod_poly_mat_clear(Ap);
        nmod_poly_mat_clear(Bp);
        return ok;
    }

    /* constant A (or zero): a single product (Dixon at order prec) */
    const slong dA = nmod_poly_mat_degree(A);
    if (dA <= 0)
        return nmod_poly_mat_solve_series_dixon(X, A, B, prec, prec);

    /* block size e >= deg(A); the rows of B are split into nc chunks of e
     * coefficients, B = sum_t x^(t e) B_t, stacked into a matrix with R =
     * nc n rows; nb blocks of e coefficients of the solution */
    const slong lenB = nmod_poly_mat_max_length(B);
    const slong e = (lenB <= 2 * dA) ? FLINT_MAX(dA, lenB) : dA;
    const slong nc = FLINT_MAX(1, (lenB + e - 1) / e);
    const slong R = nc * n;
    const slong nb = (prec + e - 1) / e;

    /* inv = A^{-1} mod x^e, slice = S_1 (coefficients 1 to 2 e - 1) */
    nmod_poly_mat_t inv, slice;
    nmod_poly_mat_init(inv, m, m, p);
    if (! nmod_poly_mat_inv_series(inv, A, 2 * e))
    {
        nmod_poly_mat_clear(inv);
        return 0;
    }
    nmod_poly_mat_init(slice, m, m, p);
    nmod_poly_mat_shift_right_trunc(slice, inv, 1, 2 * e - 1);
    nmod_poly_mat_truncate(inv, e);

    /* res: the residues, R_k in rows [k R, (k+1) R) */
    nmod_poly_mat_t res, next, win;
    nmod_poly_mat_init(res, nb * R, m, p);
    for (slong t = 0; t < nc; t++)
        for (slong i = 0; i < n; i++)
            for (slong j = 0; j < m; j++)
            {
                nmod_poly_struct * r = nmod_poly_mat_entry(res, t * n + i, j);
                nmod_poly_shift_right(r, nmod_poly_mat_entry(B, i, j), t * e);
                nmod_poly_truncate(r, e);
            }

    for (slong s = 1; s < nb; s *= 2)
    {
        /* R_s, ..., R_{s+k-1} from R_0, ..., R_{k-1} */
        const slong k = FLINT_MIN(s, nb - s);
        nmod_poly_mat_window_init(win, res, 0, 0, k * R, m);
        nmod_poly_mat_init(next, k * R, m, p);
        nmod_poly_mat_mulmid(next, win, slice, e - 1, 2 * e - 1);
        nmod_poly_mat_window_clear(win);
        nmod_poly_mat_window_init(win, res, s * R, 0, (s + k) * R, m);
        nmod_poly_mat_mulmid(win, next, A, 0, e);
        nmod_poly_mat_window_clear(win);
        nmod_poly_mat_clear(next);
        if (2 * s < nb)
            _hol_lift_slice(slice, A, inv, e);
    }

    /* Y_k = R_k inv mod x^e, and X = sum_{t,k} x^((t+k) e) (rows of Y_k
     * for B_t) mod x^prec */
    nmod_poly_mat_mulmid(res, res, inv, 0, e);
    nmod_poly_mat_zero(X);
    for (slong i = 0; i < n; i++)
        for (slong j = 0; j < m; j++)
        {
            nmod_poly_struct * x = nmod_poly_mat_entry(X, i, j);
            nmod_poly_fit_length(x, prec);
            _nmod_vec_zero(x->coeffs, prec);
            for (slong t = 0; t < nc; t++)
                for (slong k = 0; k < nb && (t + k) * e < prec; k++)
                {
                    const nmod_poly_struct * y = nmod_poly_mat_entry(res, k * R + t * n + i, j);
                    const slong off = (t + k) * e;
                    _nmod_vec_add(x->coeffs + off, x->coeffs + off, y->coeffs,
                                  FLINT_MIN(y->length, prec - off), x->mod);
                }
            _nmod_poly_set_length(x, prec);
            _nmod_poly_normalise(x);
        }

    nmod_poly_mat_clear(res);
    nmod_poly_mat_clear(slice);
    nmod_poly_mat_clear(inv);
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
    nmod_poly_mat_concat_vertical(pmat, A, b);

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
        const ulong neginv = nmod_neg(nmod_inv(d->coeffs[d->length - 1], d->mod), d->mod);
        nmod_poly_mat_t kerrow;
        nmod_poly_mat_window_init(kerrow, ker, best, 0, best + 1, m);
        nmod_poly_mat_scalar_mul_nmod(u, kerrow, neginv);
        nmod_poly_mat_window_clear(kerrow);
        nmod_poly_make_monic(den, d);
    }

    nmod_poly_mat_clear(ker);
    nmod_poly_mat_clear(pmat);
    flint_free(shift);
    flint_free(pivind);
    return (best >= 0);
}

/*------------------------------------------------------------*/
/* square systems: common tools                               */
/*------------------------------------------------------------*/

/* Degree bounds for u A = den b, with A square n x n (n >= 1) and b nonzero
 * (see nmod_poly_mat_linsolve_via_series in nmod_poly_mat_linsolve.h):
 * deg(det(A)) <= D, and the entries of b adj(A) have degree at most N.
 * Returns 0 if A has a zero row or a zero column (then A is singular), and 1
 * otherwise. */
static int _linsolve_degree_bounds(slong * N, slong * D,
                                   const nmod_poly_mat_t A,
                                   const nmod_poly_mat_t b)
{
    const slong n = A->r;
    slong * rdeg = FLINT_ARRAY_ALLOC(n, slong);
    slong * cdeg = FLINT_ARRAY_ALLOC(n, slong);
    const ulong nzr = nmod_poly_mat_row_degree(rdeg, A, NULL);
    const ulong nzc = nmod_poly_mat_column_degree(cdeg, A, NULL);
    const int ok = (nzr == (ulong) n && nzc == (ulong) n);
    if (ok)
    {
        slong sum_r = 0, sum_c = 0, min_r = rdeg[0], degb = -1, max_bc = WORD_MIN;
        for (slong k = 0; k < n; k++)
        {
            const slong dbk = nmod_poly_degree(nmod_poly_mat_entry(b, 0, k));
            sum_r += rdeg[k];
            sum_c += cdeg[k];
            min_r = FLINT_MIN(min_r, rdeg[k]);
            degb = FLINT_MAX(degb, dbk);
            if (dbk >= 0)
                max_bc = FLINT_MAX(max_bc, dbk - cdeg[k]);
        }
        *D = FLINT_MIN(sum_r, sum_c);
        *N = FLINT_MIN(degb + sum_r - min_r, max_bc + sum_c);
    }
    flint_free(rdeg);
    flint_free(cdeg);
    return ok;
}

/* Finds a such that A(a) is invertible, for A square n x n with n >= 1: a =
 * 0 if A(0) is invertible, and otherwise a random value, among a few trials
 * (if A is nonsingular, at most deg(det(A)) values are unsuitable). Returns 1
 * if such an a was found, and 0 otherwise. */
static int _linsolve_expansion_point(ulong * a, const nmod_poly_mat_t A,
                                     flint_rand_t state)
{
    const slong n = A->r;
    nmod_mat_t M;
    nmod_mat_init(M, n, n, A->modulus);
    int found = 0;
    for (slong trial = 0; trial < 8 && !found; trial++)
    {
        *a = (trial == 0) ? 0 : n_randint(state, A->modulus);
        nmod_poly_mat_evaluate_nmod(M, A, *a);
        found = (nmod_mat_rank(M) == n);
    }
    nmod_mat_clear(M);
    return found;
}

/* B = A(x + a), entrywise (B may alias A) */
static void _linsolve_taylor_shift(nmod_poly_mat_t B, const nmod_poly_mat_t A, ulong a)
{
    for (slong i = 0; i < A->r; i++)
        for (slong j = 0; j < A->c; j++)
            nmod_poly_taylor_shift(nmod_poly_mat_entry(B, i, j),
                                   nmod_poly_mat_entry(A, i, j), a);
}

/* makes den (nonzero) monic, dividing the entries of u by the same constant */
static void _linsolve_make_monic(nmod_poly_mat_t u, nmod_poly_t den)
{
    const ulong lc = den->coeffs[den->length - 1];
    if (lc != 1)
    {
        const ulong inv = nmod_inv(lc, den->mod);
        nmod_poly_scalar_mul_nmod(den, den, inv);
        nmod_poly_mat_scalar_mul_nmod(u, u, inv);
    }
}

/* divides den (nonzero) and the entries of u by their gcd */
static void _linsolve_reduce(nmod_poly_mat_t u, nmod_poly_t den)
{
    nmod_poly_t g;
    nmod_poly_init_mod(g, den->mod);
    nmod_poly_set(g, den);
    for (slong j = 0; j < u->c && nmod_poly_degree(g) > 0; j++)
        nmod_poly_gcd(g, g, nmod_poly_mat_entry(u, 0, j));
    if (nmod_poly_degree(g) > 0)
    {
        nmod_poly_divexact(den, den, g);
        for (slong j = 0; j < u->c; j++)
            nmod_poly_divexact(nmod_poly_mat_entry(u, 0, j), nmod_poly_mat_entry(u, 0, j), g);
    }
    nmod_poly_clear(g);
}

/*------------------------------------------------------------*/
/* solution over the polynomials, via power series            */
/*------------------------------------------------------------*/

/* Rational reconstruction: for s = r / q mod x^(Np + D + 1), with deg(r) <=
 * Np, deg(q) <= D, q(0) != 0, and Np >= D - 1, computes q in lowest terms,
 * up to a constant factor (if s is not of this form, q is some polynomial of
 * degree at most D).
 *
 * The coefficients of degree Np + 1 to Np + D of q s are those of degree D to
 * 2 D - 1 of q t, where t = (s div x^(Np - D + 1)) mod x^(2 D): the low part
 * s mod x^(Np - D + 1), multiplied by q, has degree at most Np. Hence q is
 * the denominator of the fraction t' / q with deg(t') < D and t' = q t mod
 * x^(2 D): the cofactor of t in the first remainder of degree less than D in
 * the Euclidean algorithm on (x^(2 D), t), which also gives it in lowest
 * terms (any solution (t', q) of these constraints is a polynomial multiple
 * of this remainder and its cofactors, since deg(t') + deg(q) < 2 D). It is
 * m11 in the output of nmod_poly_hgcd: the remainders A, B with deg(A) >= D >
 * deg(B) satisfy B = +- (m11 t - m21 x^(2 D)). */
static void _linsolve_reconstruct_den(nmod_poly_t q, const nmod_poly_t s, slong Np, slong D)
{
    if (D == 0)
    {
        nmod_poly_one(q);
        return;
    }

    nmod_poly_t t;
    nmod_poly_init_mod(t, s->mod);
    nmod_poly_shift_right(t, s, Np - D + 1);
    nmod_poly_truncate(t, 2 * D);

    if (nmod_poly_is_zero(t))
        nmod_poly_one(q);
    else
    {
        nmod_poly_t a, m12, m21, m22, Ar, Br;
        nmod_poly_init_mod(a, s->mod);
        nmod_poly_init_mod(m12, s->mod);
        nmod_poly_init_mod(m21, s->mod);
        nmod_poly_init_mod(m22, s->mod);
        nmod_poly_init_mod(Ar, s->mod);
        nmod_poly_init_mod(Br, s->mod);
        nmod_poly_set_coeff_ui(a, 2 * D, 1);
        nmod_poly_hgcd(q, m12, m21, m22, Ar, Br, a, t);
        nmod_poly_clear(a);
        nmod_poly_clear(m12);
        nmod_poly_clear(m21);
        nmod_poly_clear(m22);
        nmod_poly_clear(Ar);
        nmod_poly_clear(Br);
    }
    nmod_poly_clear(t);
}

/* solver of X A = B mod x^prec, A(0) invertible (see @ref SeriesSolve) */
typedef int (*_linsolve_series_solver_t)(nmod_poly_mat_t X,
                                         const nmod_poly_mat_t A,
                                         const nmod_poly_mat_t B,
                                         slong prec);

static int _solve_series_dixon_default_order(nmod_poly_mat_t X,
                                             const nmod_poly_mat_t A,
                                             const nmod_poly_mat_t B,
                                             slong prec)
{
    return nmod_poly_mat_solve_series_dixon(X, A, B, 0, prec);
}

/* nmod_poly_mat_linsolve_via_series, with the power series solutions
 * computed by `series`; `name` is used in exceptions */
static int _linsolve_via_series(nmod_poly_mat_t u, nmod_poly_t den,
                                const nmod_poly_mat_t A, const nmod_poly_mat_t b,
                                linsolve_den_t den_method,
                                _linsolve_series_solver_t series,
                                const char * name)
{
    const slong n = A->r;

    if (A->c != n)
        flint_throw(FLINT_ERROR, "Exception (%s). Non-square matrix.\n", name);
    if (den_method != LINSOLVE_RECONSTRUCTION && den_method != LINSOLVE_DETERMINANT)
        flint_throw(FLINT_ERROR, "Exception (%s). Invalid den_method.\n", name);

    if (n == 0 || nmod_poly_mat_is_zero(b))
    {
        nmod_poly_mat_zero(u);
        nmod_poly_one(den);
        return 1;
    }

    /* degree bounds, and expansion point a: A(a) invertible */
    slong N, D;
    ulong a = 0;
    flint_rand_t state;
    flint_rand_init(state);
    const int ok = _linsolve_degree_bounds(&N, &D, A, b)
                   && _linsolve_expansion_point(&a, A, state);
    if (!ok)
    {
        flint_rand_clear(state);
        return nmod_poly_mat_linsolve_via_kernel(u, den, A, b);
    }

    const ulong p = A->modulus;

    /* solve with A1 = A(x + a), b1 = b(x + a), whose constant coefficient A1(0) = A(a) is invertible */
    nmod_poly_mat_t As, bs;
    const nmod_poly_mat_struct * A1 = A;
    const nmod_poly_mat_struct * b1 = b;
    if (a != 0)
    {
        nmod_poly_mat_init(As, n, n, p);
        nmod_poly_mat_init(bs, 1, n, p);
        _linsolve_taylor_shift(As, A, a);
        _linsolve_taylor_shift(bs, b, a);
        A1 = As;
        b1 = bs;
    }

    nmod_poly_mat_t X, w;
    nmod_poly_t q;
    nmod_poly_mat_init(X, 1, n, p);
    nmod_poly_mat_init(w, 1, n, p);
    nmod_poly_init(q, p);

    int done = 0;
    int have_X = 0;
    if (den_method == LINSOLVE_RECONSTRUCTION)
    {
        /* X = b1 A1^{-1} mod x^prec, prec = Np + D + 1, then q = denominator
         * of a random linear combination s of its entries, and w = q X mod
         * x^prec: if w has degree at most Np, then w / q is the solution in
         * lowest terms (see nmod_poly_mat_linsolve_via_series in
         * nmod_poly_mat_linsolve.h) */
        const slong Np = FLINT_MAX(N, D - 1);
        const slong prec = Np + D + 1;
        series(X, A1, b1, prec);
        have_X = 1;

        nmod_poly_t s;
        nmod_poly_init(s, p);
        for (slong trial = 0; trial < 3 && !done; trial++)
        {
            nmod_poly_zero(s);
            for (slong j = 0; j < n; j++)
                nmod_poly_scalar_addmul_nmod(s, nmod_poly_mat_entry(X, 0, j), n_randint(state, p));
            _linsolve_reconstruct_den(q, s, Np, D);

            done = (! nmod_poly_is_zero(q) && q->coeffs[0] != 0);
            for (slong j = 0; j < n && done; j++)
            {
                nmod_poly_struct * wj = nmod_poly_mat_entry(w, 0, j);
                nmod_poly_mullow(wj, q, nmod_poly_mat_entry(X, 0, j), prec);
                done = (wj->length <= Np + 1);
            }
        }
        nmod_poly_clear(s);
    }

    if (!done)
    {
        /* w = det(A1) X mod x^(N+1), which is b1 adj(A1), and then w and
         * det(A1) are divided by their gcd */
        if (!have_X)
            series(X, A1, b1, N + 1);
        nmod_poly_mat_determinant(q, A1);
        for (slong j = 0; j < n; j++)
            nmod_poly_mullow(nmod_poly_mat_entry(w, 0, j), q, nmod_poly_mat_entry(X, 0, j), N + 1);
        _linsolve_reduce(w, q);
    }

    /* back to x: w(x - a), q(x - a) */
    if (a != 0)
    {
        const ulong mina = p - a;
        _linsolve_taylor_shift(w, w, mina);
        nmod_poly_taylor_shift(q, q, mina);
        nmod_poly_mat_clear(As);
        nmod_poly_mat_clear(bs);
    }
    _linsolve_make_monic(w, q);

    nmod_poly_mat_swap_entrywise(u, w);
    nmod_poly_swap(den, q);

    nmod_poly_mat_clear(X);
    nmod_poly_mat_clear(w);
    nmod_poly_clear(q);
    flint_rand_clear(state);
    return 1;
}

int nmod_poly_mat_linsolve_via_series(nmod_poly_mat_t u,
                                      nmod_poly_t den,
                                      const nmod_poly_mat_t A,
                                      const nmod_poly_mat_t b,
                                      linsolve_den_t den_method)
{
    return _linsolve_via_series(u, den, A, b, den_method,
                                _solve_series_dixon_default_order,
                                "nmod_poly_mat_linsolve_via_series");
}

int nmod_poly_mat_linsolve_via_high_order_lifting(nmod_poly_mat_t u,
                                                  nmod_poly_t den,
                                                  const nmod_poly_mat_t A,
                                                  const nmod_poly_mat_t b,
                                                  linsolve_den_t den_method)
{
    return _linsolve_via_series(u, den, A, b, den_method,
                                nmod_poly_mat_solve_series_high_order_lifting,
                                "nmod_poly_mat_linsolve_via_high_order_lifting");
}

/*------------------------------------------------------------*/
/* solution over the polynomials, by evaluation/interpolation */
/*------------------------------------------------------------*/

/* u = b adj(M), for a constant n x n matrix M which is singular (n >= 1) and
 * a constant vector b of length n. For n >= 2: adj(M) is zero if rank(M) <
 * n - 1, and otherwise, since M adj(M) = adj(M) M = det(M) I = 0, adj(M) =
 * lambda v w for nonzero vectors v and w such that M v = 0 and w M = 0, so
 * that u = lambda (b v) w. The scalar lambda is found from one entry:
 * adj(M)[k,j] = (-1)^(j+k) det(M without row j and column k), for indices j,
 * k such that w[j] and v[k] are nonzero. */
static void _linsolve_adjugate_singular(nn_ptr u, const nmod_mat_t M, nn_srcptr b, nmod_t mod)
{
    const slong n = M->r;

    if (n == 1)  /* adj(M) = [1] */
    {
        u[0] = b[0];
        return;
    }

    _nmod_vec_zero(u, n);
    if (nmod_mat_rank(M) < n - 1)
        return;

    /* v: right kernel of M; w: right kernel of M^T */
    nmod_mat_t K, Mt, Kt;
    nmod_mat_init(K, n, n, mod.n);
    nmod_mat_init(Mt, n, n, mod.n);
    nmod_mat_init(Kt, n, n, mod.n);
    nmod_mat_nullspace(K, M);
    nmod_mat_transpose(Mt, M);
    nmod_mat_nullspace(Kt, Mt);

    slong j0 = 0, k0 = 0;
    while (nmod_mat_entry(Kt, j0, 0) == 0)
        j0++;
    while (nmod_mat_entry(K, k0, 0) == 0)
        k0++;

    /* c = (-1)^(j0+k0) det(M without row j0 and column k0) = adj(M)[k0,j0] */
    nmod_mat_t S;
    nmod_mat_init(S, n - 1, n - 1, mod.n);
    for (slong i = 0, ii = 0; i < n; i++)
    {
        if (i == j0)
            continue;
        for (slong k = 0, kk = 0; k < n; k++)
            if (k != k0)
                nmod_mat_entry(S, ii, kk++) = nmod_mat_entry(M, i, k);
        ii++;
    }
    ulong c = nmod_mat_det(S);
    if ((j0 + k0) & 1)
        c = nmod_neg(c, mod);

    /* lambda (b v) = c (b v) / (v[k0] w[j0]) */
    ulong bv = 0;
    for (slong k = 0; k < n; k++)
        bv = nmod_add(bv, nmod_mul(b[k], nmod_mat_entry(K, k, 0), mod), mod);
    const ulong vw = nmod_mul(nmod_mat_entry(K, k0, 0), nmod_mat_entry(Kt, j0, 0), mod);
    const ulong scale = nmod_mul(nmod_mul(c, bv, mod), nmod_inv(vw, mod), mod);
    for (slong j = 0; j < n; j++)
        u[j] = nmod_mul(scale, nmod_mat_entry(Kt, j, 0), mod);

    nmod_mat_clear(S);
    nmod_mat_clear(K);
    nmod_mat_clear(Mt);
    nmod_mat_clear(Kt);
}

int _nmod_poly_mat_linsolve_geometric(nmod_poly_mat_t u,
                                      nmod_poly_t den,
                                      const nmod_poly_mat_t A,
                                      const nmod_poly_mat_t b,
                                      slong chunk)
{
    const slong n = A->r;

    if (A->c != n)
        flint_throw(FLINT_ERROR, "Exception (nmod_poly_mat_linsolve_geometric). Non-square matrix.\n");

    if (n == 0 || nmod_poly_mat_is_zero(b))
    {
        nmod_poly_mat_zero(u);
        nmod_poly_one(den);
        return 1;
    }

    /* det(A) and b adj(A) have degree less than L; the L points are q^k,
     * k < L, with q = r^2 for r of order at least 2 L */
    slong N, D;
    if (! _linsolve_degree_bounds(&N, &D, A, b))
        return nmod_poly_mat_linsolve_via_kernel(u, den, A, b);
    const slong L = FLINT_MAX(N, D) + 1;
    const ulong p = A->modulus;
    if (! NMOD_POLY_CAN_USE_GEOMETRIC(p, L) || ! n_is_prime(p))
        return nmod_poly_mat_linsolve_via_kernel(u, den, A, b);
    nmod_t mod;
    nmod_init(&mod, p);
    const ulong r = nmod_find_root(2 * L, mod);
    if (r == 0)
        return nmod_poly_mat_linsolve_via_kernel(u, den, A, b);

    /* lengths: lenA <= D + 1 <= L (an entry has degree at most its row
     * degree) and lenb <= N + 1 <= L (N >= deg(b)) */
    const slong lenA = nmod_poly_mat_max_length(A);
    const slong lenb = nmod_poly_mat_max_length(b);

    /* G: the L points, for the evaluation of b and the interpolations;
     * Gc: the first c points, for the evaluation of A by chunks of c points:
     * f(q^(k0 + k)) = f_{k0}(q^k) with f_{k0}(x) = f(q^k0 x), whose
     * coefficients are those of f scaled by the powers of z = q^k0. By
     * default, the chunks hold at most max(n^2 lenA, 2^22) values: longer
     * chunks make fewer and longer products, which are more efficient */
    if (chunk <= 0)
        chunk = (WORD(1) << 22) / (n * n);
    const slong c = FLINT_MIN(L, FLINT_MAX(lenA, chunk));
    nmod_geometric_progression_t G, Gc;
    _nmod_geometric_progression_init_function(G, r, L, mod, UWORD(3));
    _nmod_geometric_progression_init_fft(G, lenb, L, 0);
    _nmod_geometric_progression_init_function(Gc, r, c, mod, UWORD(1));
    _nmod_geometric_progression_init_fft(Gc, lenA, 0, 0);

    /* evb[j*L + k] = b[j](q^k) */
    nn_ptr evb = _nmod_vec_init(n * L);
    for (slong j = 0; j < n; j++)
    {
        const nmod_poly_struct * f = nmod_poly_mat_entry(b, 0, j);
        _nmod_poly_evaluate_geometric_nmod_vec_fast_precomp(evb + j * L, f->coeffs, f->length, G, L, mod);
    }

    /* vals[k] = det(A(q^k)); vals[(1+j)*L + k] = (b adj(A))[j](q^k) */
    nn_ptr vals = _nmod_vec_init((n + 1) * L);

    /* chunk: ev[(i*n + j)*c + k] = A[i][j](q^(k0 + k)) */
    nn_ptr ev = _nmod_vec_init(n * n * c);
    nn_ptr zpow = _nmod_vec_init(lenA);
    nn_ptr scaled = _nmod_vec_init(lenA);
    const ulong q = nmod_mul(r, r, mod);
    const ulong qc = nmod_pow_ui(q, c, mod);
    ulong z = 1;

    nmod_mat_t M, PB, X;
    nmod_mat_init(M, n, n, p);
    nmod_mat_init(PB, n, 1, p);
    nmod_mat_init(X, n, 1, p);
    slong * perm = _perm_init(n);
    nn_ptr bpt = _nmod_vec_init(n);
    nn_ptr upt = _nmod_vec_init(n);

    for (slong k0 = 0; k0 < L; k0 += c)
    {
        const slong cl = FLINT_MIN(c, L - k0);

        zpow[0] = 1;
        for (slong t = 1; t < lenA; t++)
            zpow[t] = nmod_mul(zpow[t - 1], z, mod);
        for (slong i = 0; i < n; i++)
            for (slong j = 0; j < n; j++)
            {
                const nmod_poly_struct * f = nmod_poly_mat_entry(A, i, j);
                nn_ptr v = ev + (i * n + j) * c;
                if (f->length <= 1)
                {
                    const ulong f0 = (f->length == 1) ? f->coeffs[0] : 0;
                    for (slong k = 0; k < cl; k++)
                        v[k] = f0;
                }
                else
                {
                    for (slong t = 0; t < f->length; t++)
                        scaled[t] = nmod_mul(f->coeffs[t], zpow[t], mod);
                    _nmod_poly_evaluate_geometric_nmod_vec_fast_precomp(v, scaled, f->length, Gc, cl, mod);
                }
            }

        for (slong k = 0; k < cl; k++)
        {
            /* M = A(q^(k0+k))^T, LU decomposition P M = L U; then
             * det(A(pt)) = +- prod(U[i][i]) and, from M x = b(pt)^T,
             * b(pt) adj(A(pt)) = det(A(pt)) x^T */
            for (slong i = 0; i < n; i++)
                for (slong j = 0; j < n; j++)
                    nmod_mat_entry(M, j, i) = ev[(i * n + j) * c + k];
            for (slong j = 0; j < n; j++)
                bpt[j] = evb[j * L + k0 + k];

            if (nmod_mat_lu(perm, M, 1) == n)
            {
                ulong d = 1;
                for (slong i = 0; i < n; i++)
                    d = nmod_mul(d, nmod_mat_entry(M, i, i), mod);
                if (_perm_parity(perm, n))
                    d = nmod_neg(d, mod);
                for (slong i = 0; i < n; i++)
                    nmod_mat_entry(PB, i, 0) = bpt[perm[i]];
                nmod_mat_solve_tril(X, M, PB, 1);
                nmod_mat_solve_triu(X, M, X, 0);
                vals[k0 + k] = d;
                for (slong j = 0; j < n; j++)
                    vals[(1 + j) * L + k0 + k] = nmod_mul(d, nmod_mat_entry(X, j, 0), mod);
            }
            else
            {
                /* singular point (M was overwritten) */
                for (slong i = 0; i < n; i++)
                    for (slong j = 0; j < n; j++)
                        nmod_mat_entry(M, i, j) = ev[(i * n + j) * c + k];
                _linsolve_adjugate_singular(upt, M, bpt, mod);
                vals[k0 + k] = 0;
                for (slong j = 0; j < n; j++)
                    vals[(1 + j) * L + k0 + k] = upt[j];
            }
        }

        z = nmod_mul(z, qc, mod);
    }

    _perm_clear(perm);
    nmod_mat_clear(M);
    nmod_mat_clear(PB);
    nmod_mat_clear(X);
    _nmod_vec_clear(bpt);
    _nmod_vec_clear(upt);
    _nmod_vec_clear(ev);
    _nmod_vec_clear(zpow);
    _nmod_vec_clear(scaled);
    _nmod_vec_clear(evb);

    /* interpolation */
    nmod_poly_mat_t w;
    nmod_poly_t d;
    nmod_poly_mat_init(w, 1, n, p);
    nmod_poly_init(d, p);
    nmod_poly_interpolate_geometric_nmod_vec_fast_precomp(d, vals, G, L);
    int ok = ! nmod_poly_is_zero(d);
    for (slong j = 0; j < n && ok; j++)
        nmod_poly_interpolate_geometric_nmod_vec_fast_precomp(nmod_poly_mat_entry(w, 0, j),
                                                            vals + (1 + j) * L, G, L);
    _nmod_vec_clear(vals);
    nmod_geometric_progression_clear(G);
    nmod_geometric_progression_clear(Gc);

    if (ok)
    {
        _linsolve_reduce(w, d);
        _linsolve_make_monic(w, d);
        nmod_poly_mat_swap_entrywise(u, w);
        nmod_poly_swap(den, d);
    }

    nmod_poly_mat_clear(w);
    nmod_poly_clear(d);

    /* det(A) = 0: A is singular */
    if (!ok)
        return nmod_poly_mat_linsolve_via_kernel(u, den, A, b);
    return 1;
}

int nmod_poly_mat_linsolve_geometric(nmod_poly_mat_t u,
                                     nmod_poly_t den,
                                     const nmod_poly_mat_t A,
                                     const nmod_poly_mat_t b)
{
    return _nmod_poly_mat_linsolve_geometric(u, den, A, b, 0);
}
