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
