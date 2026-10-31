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
#include <flint/nmod_mat.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>

#include "nmod_poly_mat_forms.h"
#include "nmod_poly_mat_inverse.h"
#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_kernel.h"
#include "nmod_poly_mat_approximant.h"
#include "nmod_poly_mat_multiply.h"
#include "nmod_poly_mat_row_basis.h"

/* In the general interface, a matrix with more than ROW_BASIS_TALL_RATIO
 * times as many rows as columns is processed by blocks of rows */
#define ROW_BASIS_TALL_RATIO 3

/* In the general interface, the iterative weak Popov form algorithm of
 * Mulders and Storjohann is used when cdim * deg is at most
 * NMOD_POLY_MAT_ROW_BASIS_WEAK_POPOV_THRES (see the header), or
 * ROW_BASIS_WEAK_POPOV_THRES_TALL for a matrix processed by blocks: for
 * tall matrices with nontrivial row space, the weak Popov approach remained
 * faster in all tested cases (Intel Xeon, AVX-512, 60-bit prime, dimensions
 * up to 400 x 40, degree up to 1000) 
 * \todo thresholds remain to be studied more thoroughly
 * */
#define ROW_BASIS_WEAK_POPOV_THRES_TALL 20000

/*------------------------------------------------------------*/
/* helpers                                                    */
/*------------------------------------------------------------*/

/* zero the rows of mat starting from index i0 */
static void _zero_rows_from(nmod_poly_mat_t mat, slong i0)
{
    for (slong i = i0; i < mat->r; i++)
        for (slong j = 0; j < mat->c; j++)
            nmod_poly_zero(nmod_poly_mat_entry(mat, i, j));
}

/* zero the columns of mat starting from index j0 */
static void _zero_columns_from(nmod_poly_mat_t mat, slong j0)
{
    for (slong i = 0; i < mat->r; i++)
        for (slong j = j0; j < mat->c; j++)
            nmod_poly_zero(nmod_poly_mat_entry(mat, i, j));
}

/* the first src->r rows of mat become src (with mat->c == src->c), the other
 * rows of mat become zero; src is destroyed if swap is nonzero */
static void _set_top_rows(nmod_poly_mat_t mat, nmod_poly_mat_t src, int swap)
{
    if (swap)
    {
        for (slong i = 0; i < src->r; i++)
            for (slong j = 0; j < src->c; j++)
                nmod_poly_swap(nmod_poly_mat_entry(mat, i, j), nmod_poly_mat_entry(src, i, j));
    }
    else
    {
        for (slong i = 0; i < src->r; i++)
            for (slong j = 0; j < src->c; j++)
                nmod_poly_set(nmod_poly_mat_entry(mat, i, j), nmod_poly_mat_entry(src, i, j));
    }
    _zero_rows_from(mat, src->r);
}

/* the first src->c columns of mat become src (with mat->r == src->r), the
 * other columns of mat become zero; src is destroyed if swap is nonzero */
static void _set_left_columns(nmod_poly_mat_t mat, nmod_poly_mat_t src, int swap)
{
    if (swap)
    {
        for (slong i = 0; i < src->r; i++)
            for (slong j = 0; j < src->c; j++)
                nmod_poly_swap(nmod_poly_mat_entry(mat, i, j), nmod_poly_mat_entry(src, i, j));
    }
    else
    {
        for (slong i = 0; i < src->r; i++)
            for (slong j = 0; j < src->c; j++)
                nmod_poly_set(nmod_poly_mat_entry(mat, i, j), nmod_poly_mat_entry(src, i, j));
    }
    _zero_columns_from(mat, src->c);
}

/* the first n rows of mat become the n x n identity (n <= mat->r, n == mat->c),
 * the other rows become zero */
static void _set_top_identity(nmod_poly_mat_t mat, slong n)
{
    for (slong i = 0; i < n; i++)
        for (slong j = 0; j < mat->c; j++)
        {
            if (i == j)
                nmod_poly_one(nmod_poly_mat_entry(mat, i, j));
            else
                nmod_poly_zero(nmod_poly_mat_entry(mat, i, j));
        }
    _zero_rows_from(mat, n);
}

/* Let M be a matrix with r columns and rank r, and P M = L U its LU
 * decomposition as computed in place by nmod_mat_lu (stored in LU, with
 * permutation perm). The first r rows of P M form the invertible matrix
 * L1 U, where L1 is the unit lower triangular top r x r block of L and U is
 * upper triangular, both stored in the top r x r block of LU. This computes
 * inv = (L1 U)^{-1} = U^{-1} L1^{-1}, the inverse of the submatrix of M formed
 * by its rows perm[0], ..., perm[r-1]. */
static void _inverse_from_lu(nmod_mat_t inv, const nmod_mat_t LU, slong r)
{
    nmod_mat_t T;
    nmod_mat_window_init(T, LU, 0, 0, r, r);
    nmod_mat_one(inv);
    nmod_mat_solve_tril(inv, T, inv, 1);
    nmod_mat_solve_triu(inv, T, inv, 0);
    nmod_mat_window_clear(T);
}

/*------------------------------------------------------------*/
/* unimodularity certificate                                  */
/*------------------------------------------------------------*/

int _nmod_poly_mat_row_basis_certify_unimodular(const nmod_poly_mat_t pmat,
                                                const slong * shift,
                                                slong sum_t)
{
    const slong m = pmat->r;
    const slong n = pmat->c;

    slong sum_s = 0;
    for (slong i = 0; i < m; i++)
        sum_s += shift[i];

    /* first, bound w <= |cdeg_{-s}(pmat)| */
    slong sum_u = 0;
    for (slong j = 0; j < n; j++)
    {
        slong u = WORD_MIN;
        for (slong i = 0; i < m; i++)
        {
            const slong d = nmod_poly_degree(nmod_poly_mat_entry(pmat, i, j));
            if (d >= 0 && d - shift[i] > u)
                u = d - shift[i];
        }
        if (u == WORD_MIN)  /* zero column: not of full column rank */
            return 0;
        sum_u += u;
    }
    if (sum_t == sum_s + sum_u)
        return 1;

    /* then, the sharper bound given by the best assignment of the columns to
     * distinct rows (which exists since pmat has full column rank) */
    slong * ms = FLINT_ARRAY_ALLOC(m, slong);
    for (slong i = 0; i < m; i++)
        ms[i] = -shift[i];
    const slong w_bound = nmod_poly_mat_degree_det_bound(NULL, pmat, ms, COL_UPPER);
    flint_free(ms);
    return (sum_t == sum_s + w_bound);
}

/* Input: K, an s-reduced basis (first k rows of K) of the left kernel of
 * some m x n matrix of rank r = m - k >= 1, with t = rdeg_s(K); the shift s
 * is such that the left factor C of the matrix has (-s)-column degree <= 0
 * (this holds when s is at least the row degree of the matrix)
 *
 * Output: C (initialized as m x r), a (-s)-minimal right kernel basis of K,
 * and tt = -cdeg_{-s}(C) (length r, nonnegative entries).
 *
 * Following [ZL13, Lemma 4.3], C is formed by the approximants of
 * nonpositive (-s)-degree in a (-s)-minimal approximant basis of K^T at order
 * t+1, the latter being obtained from an approximant basis at the uniform
 * order max(t)+1 of K^T with its column l multiplied by x^(max(t) - t[l]).
 */
static void _right_kernel_nonpositive(nmod_poly_mat_t C, slong * tt,
                                      const nmod_poly_mat_t K, const slong * t, slong k,
                                      const slong * s)
{
    const slong m = K->c;
    const slong r = m - k;

    /* beta = max(t, 0): the order beta+1 of the approximant basis below is
     * positive, and the resulting order for column l of K^T is t[l]+1 (no
     * constraint if t[l] < 0, which may only happen for kernel vectors
     * supported on zero rows of the matrix) */
    slong beta = 0;
    for (slong l = 0; l < k; l++)
        if (t[l] > beta)
            beta = t[l];

    /* M = K^T * diag(x^(beta - t)) */
    nmod_poly_mat_t M;
    nmod_poly_mat_init(M, m, k, K->modulus);
    for (slong l = 0; l < k; l++)
        for (slong i = 0; i < m; i++)
            nmod_poly_shift_left(nmod_poly_mat_entry(M, i, l),
                                 nmod_poly_mat_entry(K, l, i),
                                 beta - t[l]);

    nmod_poly_mat_t P;
    nmod_poly_mat_init(P, m, m, K->modulus);
    slong * sh = FLINT_ARRAY_ALLOC(m, slong);
    for (slong i = 0; i < m; i++)
        sh[i] = -s[i];

    nmod_poly_mat_pmbasis(P, sh, M, beta + 1);

    /* rows of nonpositive (-s)-degree */
    slong nb = 0;
    for (slong i = 0; i < m; i++)
        if (sh[i] <= 0)
            nb++;
    if (nb != r)
        flint_throw(FLINT_ERROR, "Exception (_nmod_poly_mat_row_basis_zl). "
                    "Unexpected dimension of the right kernel (found %wd, expected %wd).\n", nb, r);

    nb = 0;
    for (slong i = 0; i < m; i++)
    {
        if (sh[i] <= 0)
        {
            for (slong j = 0; j < m; j++)
                nmod_poly_swap(nmod_poly_mat_entry(C, j, nb), nmod_poly_mat_entry(P, i, j));
            tt[nb] = -sh[i];
            nb++;
        }
    }

    flint_free(sh);
    nmod_poly_mat_clear(P);
    nmod_poly_mat_clear(M);
}

/*------------------------------------------------------------*/
/* Zhou-Labahn row basis, low level                           */
/*------------------------------------------------------------*/

slong _nmod_poly_mat_row_basis_zl(nmod_poly_mat_t rbas,
                                  nmod_poly_mat_t lfac,
                                  const nmod_poly_mat_t pmat,
                                  const slong * shift,
                                  const nmod_poly_mat_t ker,
                                  const slong * ker_rdeg,
                                  slong nullity,
                                  int check_unimodular)
{
    const slong m = pmat->r;
    const slong n = pmat->c;
    const ulong mod = pmat->modulus;

    /* empty matrix */
    if (m == 0 || n == 0)
    {
        _zero_rows_from(rbas, 0);
        if (lfac)
            _zero_columns_from(lfac, 0);
        return 0;
    }

    /* shift */
    slong * s = FLINT_ARRAY_ALLOC(m, slong);
    if (shift)
        for (slong i = 0; i < m; i++)
            s[i] = shift[i];
    else
        nmod_poly_mat_row_degree_zero(s, pmat, NULL);

    /* s-reduced left kernel basis: K (first k rows of Kmat), t = rdeg_s(K) */
    nmod_poly_mat_t kbuf;
    slong * tbuf = NULL;
    const slong * t;
    slong k;
    const nmod_poly_mat_struct * Kmat;
    if (ker == NULL)
    {
        nmod_poly_mat_t copy;
        nmod_poly_mat_init_set(copy, pmat);
        nmod_poly_mat_init(kbuf, m, m, mod);
        slong * pivind = FLINT_ARRAY_ALLOC(m, slong);
        tbuf = FLINT_ARRAY_ALLOC(m, slong);
        for (slong i = 0; i < m; i++)
            tbuf[i] = s[i];
        k = nmod_poly_mat_kernel_zls_approx(kbuf, pivind, tbuf, copy);
        flint_free(pivind);
        nmod_poly_mat_clear(copy);
        t = tbuf;
        Kmat = kbuf;
    }
    else
    {
        k = nullity;
        t = ker_rdeg;
        Kmat = ker;
    }

    const slong r = m - k;

    if (r == 0)
    {
        /* pmat is zero */
        _zero_rows_from(rbas, 0);
        if (lfac)
            _zero_columns_from(lfac, 0);
    }

    else if (k == 0)
    {
        /* full row rank: pmat is a row basis of itself (note m <= n) */
        for (slong i = 0; i < m; i++)
            for (slong j = 0; j < n; j++)
                nmod_poly_set(nmod_poly_mat_entry(rbas, i, j), nmod_poly_mat_entry(pmat, i, j));
        _zero_rows_from(rbas, m);
        if (lfac)
        {
            _zero_columns_from(lfac, 0);
            for (slong i = 0; i < m; i++)
                nmod_poly_one(nmod_poly_mat_entry(lfac, i, i));
        }
    }

    else
    {
        int done = 0;

        /* normalize the shift: only the differences between its entries
         * matter, and s - c with c = min(s[i] - rdeg(pmat)[i]) over the
         * nonzero rows of pmat is the smallest translate of s which bounds
         * the row degree of pmat; this is required in step 3 (see
         * documentation), and minimizes the order and degree bounds
         * involved; K is (s-c)-reduced with (s-c)-row degree t - c */
        slong * tn = FLINT_ARRAY_ALLOC(k, slong);
        {
            slong * rdeg = FLINT_ARRAY_ALLOC(m, slong);
            nmod_poly_mat_row_degree(rdeg, pmat, NULL);
            slong c = WORD_MAX;
            for (slong i = 0; i < m; i++)
                if (rdeg[i] >= 0 && s[i] - rdeg[i] < c)
                    c = s[i] - rdeg[i];
            flint_free(rdeg);
            /* here r >= 1, hence pmat is nonzero and c < WORD_MAX */
            for (slong i = 0; i < m; i++)
                s[i] -= c;
            for (slong l = 0; l < k; l++)
                tn[l] = t[l] - c;
            t = tn;
        }

        /* certificate of unimodularity of the row basis */
        if (check_unimodular && r == n)
        {
            slong sum_t = 0;
            for (slong l = 0; l < k; l++)
                sum_t += t[l];

            if (_nmod_poly_mat_row_basis_certify_unimodular(pmat, s, sum_t))
            {
                /* row space of pmat is K[x]^{1 x n} */
                _set_top_identity(rbas, n);
                if (lfac)
                {
                    for (slong i = 0; i < m; i++)
                        for (slong j = 0; j < n; j++)
                            nmod_poly_set(nmod_poly_mat_entry(lfac, i, j), nmod_poly_mat_entry(pmat, i, j));
                    _zero_columns_from(lfac, n);
                }
                done = 1;
            }
        }

        if (!done)
        {
            /* left factor C: (-s)-minimal right kernel basis of K */
            nmod_poly_mat_t C;
            nmod_poly_mat_init(C, m, r, mod);
            slong * tt = FLINT_ARRAY_ALLOC(r, slong);
            _right_kernel_nonpositive(C, tt, Kmat, t, k, s);

            /* J: r rows such that C(0)[J,:] is invertible, and the inverse
             * inv0 of C(0)[J,:] */
            nmod_mat_t C0, inv0;
            nmod_mat_init(C0, m, r, mod);
            nmod_poly_mat_get_coeff_mat(C0, C, 0);
            slong * perm = FLINT_ARRAY_ALLOC(m, slong);
            const slong rk0 = nmod_mat_lu(perm, C0, 0);
            if (rk0 != r)
                flint_throw(FLINT_ERROR, "Exception (_nmod_poly_mat_row_basis_zl). "
                            "Left factor is not of full rank at 0.\n");
            nmod_mat_init(inv0, r, r, mod);
            _inverse_from_lu(inv0, C0, r);
            nmod_mat_clear(C0);

            nmod_poly_mat_t CJ, AJ, X, R;
            nmod_poly_mat_init(CJ, r, r, mod);
            nmod_poly_mat_init(AJ, r, n, mod);
            for (slong i = 0; i < r; i++)
            {
                for (slong j = 0; j < r; j++)
                    nmod_poly_set(nmod_poly_mat_entry(CJ, i, j), nmod_poly_mat_entry(C, perm[i], j));
                for (slong j = 0; j < n; j++)
                    nmod_poly_set(nmod_poly_mat_entry(AJ, i, j), nmod_poly_mat_entry(pmat, perm[i], j));
            }
            flint_free(perm);

            /* R = CJ^{-1} * AJ, which has row degree <= tt */
            slong sigma = 0;
            for (slong i = 0; i < r; i++)
                if (tt[i] > sigma)
                    sigma = tt[i];
            sigma += 1;

            nmod_poly_mat_init(X, r, r, mod);
            _nmod_poly_mat_inv_series(X, CJ, inv0, sigma);
            nmod_mat_clear(inv0);
            nmod_poly_mat_init(R, r, n, mod);
            nmod_poly_mat_mulmid(R, X, AJ, 0, sigma);

            _set_top_rows(rbas, R, 1);
            if (lfac)
                _set_left_columns(lfac, C, 1);

            nmod_poly_mat_clear(R);
            nmod_poly_mat_clear(X);
            nmod_poly_mat_clear(AJ);
            nmod_poly_mat_clear(CJ);
            flint_free(tt);
            nmod_poly_mat_clear(C);
        }

        flint_free(tn);
    }

    if (ker == NULL)
    {
        nmod_poly_mat_clear(kbuf);
        flint_free(tbuf);
    }
    flint_free(s);

    return r;
}

/*------------------------------------------------------------*/
/* interfaces                                                 */
/*------------------------------------------------------------*/

slong nmod_poly_mat_row_basis_zl(nmod_poly_mat_t rbas,
                                 nmod_poly_mat_t lfac,
                                 const nmod_poly_mat_t pmat)
{
    return _nmod_poly_mat_row_basis_zl(rbas, lfac, pmat, NULL, NULL, NULL, 0, 1);
}

slong _nmod_poly_mat_row_basis_weak_popov(nmod_poly_mat_t rbas,
                                          const nmod_poly_mat_t pmat)
{
    const slong m = pmat->r;
    const slong n = pmat->c;

    if (m == 0 || n == 0)
    {
        _zero_rows_from(rbas, 0);
        return 0;
    }

    nmod_poly_mat_t copy;
    nmod_poly_mat_init_set(copy, pmat);
    slong * pivind = FLINT_ARRAY_ALLOC(m, slong);
    /* zero rows are moved to the bottom */
    const slong r = nmod_poly_mat_weak_popov_iter(copy, NULL, NULL, pivind, NULL, ROW_LOWER);
    for (slong i = 0; i < r; i++)
        for (slong j = 0; j < n; j++)
            nmod_poly_swap(nmod_poly_mat_entry(rbas, i, j), nmod_poly_mat_entry(copy, i, j));
    _zero_rows_from(rbas, r);
    flint_free(pivind);
    nmod_poly_mat_clear(copy);
    return r;
}

/* for sorting rows by degree */
typedef struct
{
    slong deg;
    slong idx;
}
_deg_idx_t;

static int _deg_idx_cmp(const void * a, const void * b)
{
    const _deg_idx_t * x = (const _deg_idx_t *) a;
    const _deg_idx_t * y = (const _deg_idx_t *) b;
    if (x->deg != y->deg)
        return (x->deg < y->deg) ? -1 : 1;
    return (x->idx < y->idx) ? -1 : (x->idx > y->idx);
}

slong _nmod_poly_mat_row_basis_zl_blocks(nmod_poly_mat_t rbas,
                                         const nmod_poly_mat_t pmat,
                                         slong bsize)
{
    const slong m = pmat->r;
    const slong n = pmat->c;

    if (m == 0 || n == 0)
    {
        _zero_rows_from(rbas, 0);
        return 0;
    }

    if (bsize < 1)
        bsize = 1;
    else if (bsize > m)
        bsize = m;

    /* tall matrix: nonzero rows sorted by increasing degree [ZL13, Section 6] */
    slong * rdeg = FLINT_ARRAY_ALLOC(m, slong);
    nmod_poly_mat_row_degree(rdeg, pmat, NULL);
    _deg_idx_t * order = FLINT_ARRAY_ALLOC(m, _deg_idx_t);
    slong nnz = 0;
    for (slong i = 0; i < m; i++)
    {
        if (rdeg[i] >= 0)
        {
            order[nnz].deg = rdeg[i];
            order[nnz].idx = i;
            nnz++;
        }
    }
    flint_free(rdeg);
    qsort(order, nnz, sizeof(_deg_idx_t), _deg_idx_cmp);

    if (nnz == 0)
    {
        flint_free(order);
        _zero_rows_from(rbas, 0);
        return 0;
    }

    /* first block: n + bsize rows (or all nonzero rows if fewer) */
    slong pos = FLINT_MIN(nnz, n + bsize);
    nmod_poly_mat_t W, Rcur;
    nmod_poly_mat_init(W, pos, n, pmat->modulus);
    for (slong i = 0; i < pos; i++)
        for (slong j = 0; j < n; j++)
            nmod_poly_set(nmod_poly_mat_entry(W, i, j), nmod_poly_mat_entry(pmat, order[i].idx, j));
    nmod_poly_mat_init(Rcur, n, n, pmat->modulus);
    slong r = _nmod_poly_mat_row_basis_zl(Rcur, NULL, W, NULL, NULL, NULL, 0, 1);
    nmod_poly_mat_clear(W);

    /* next blocks: bsize rows at a time, appended to the current row basis */
    while (pos < nnz)
    {
        /* if the row space is already K[x]^{1 x n}, we are done */
        if (r == n && nmod_poly_mat_is_one(Rcur))
            break;

        const slong nb = FLINT_MIN(bsize, nnz - pos);
        nmod_poly_mat_init(W, r + nb, n, pmat->modulus);
        for (slong i = 0; i < r; i++)
            for (slong j = 0; j < n; j++)
                nmod_poly_swap(nmod_poly_mat_entry(W, i, j), nmod_poly_mat_entry(Rcur, i, j));
        for (slong i = 0; i < nb; i++)
            for (slong j = 0; j < n; j++)
                nmod_poly_set(nmod_poly_mat_entry(W, r + i, j), nmod_poly_mat_entry(pmat, order[pos + i].idx, j));
        r = _nmod_poly_mat_row_basis_zl(Rcur, NULL, W, NULL, NULL, NULL, 0, 1);
        nmod_poly_mat_clear(W);
        pos += nb;
    }

    /* copy result */
    for (slong i = 0; i < r; i++)
        for (slong j = 0; j < n; j++)
            nmod_poly_swap(nmod_poly_mat_entry(rbas, i, j), nmod_poly_mat_entry(Rcur, i, j));
    _zero_rows_from(rbas, r);

    nmod_poly_mat_clear(Rcur);
    flint_free(order);
    return r;
}

/* For rank deficient inputs, the unimodularity certificate does not apply
 * and the weak Popov approach remains faster up to about
 * ROW_BASIS_RANK_DEFICIENT_FACTOR times the usual threshold (Intel Xeon,
 * AVX-512, 60-bit prime, 2n x n inputs of rank n/2, n = 4 ... 30) */
#define ROW_BASIS_RANK_DEFICIENT_FACTOR 3

slong nmod_poly_mat_row_basis(nmod_poly_mat_t rbas, const nmod_poly_mat_t pmat)
{
    const slong m = pmat->r;
    const slong n = pmat->c;
    const slong d = nmod_poly_mat_degree(pmat);

    if (m == 0 || n == 0 || d < 0)
    {
        _zero_rows_from(rbas, 0);
        return 0;
    }

    const int tall = (m > ROW_BASIS_TALL_RATIO * n);

    /* small degree: weak Popov form, iterative algorithm */
    if (n * d <= (tall ? ROW_BASIS_WEAK_POPOV_THRES_TALL : NMOD_POLY_MAT_ROW_BASIS_WEAK_POPOV_THRES))
        return _nmod_poly_mat_row_basis_weak_popov(rbas, pmat);

    /* moderate degree and (probably) rank deficient: the rank at a point is
     * a lower bound on the rank, this only affects the choice of algorithm */
    if (!tall && n * d <= ROW_BASIS_RANK_DEFICIENT_FACTOR * NMOD_POLY_MAT_ROW_BASIS_WEAK_POPOV_THRES
        && nmod_poly_mat_rank_at_some_point(pmat) < FLINT_MIN(m, n))
        return _nmod_poly_mat_row_basis_weak_popov(rbas, pmat);

    if (!tall)
        return _nmod_poly_mat_row_basis_zl(rbas, NULL, pmat, NULL, NULL, NULL, 0, 1);

    return _nmod_poly_mat_row_basis_zl_blocks(rbas, pmat, n);
}

/*------------------------------------------------------------*/
/* left factor from a reduced row basis                       */
/*------------------------------------------------------------*/

int nmod_poly_mat_row_basis_left_factor(nmod_poly_mat_t lfac,
                                        const nmod_poly_mat_t pmat,
                                        const nmod_poly_mat_t rbas,
                                        slong r,
                                        const slong * shift,
                                        int check)
{
    const slong m = pmat->r;
    const slong n = pmat->c;
    const ulong mod = pmat->modulus;

    nmod_poly_mat_zero(lfac);

    if (m == 0)
        return 1;

    /* r == 0: pmat must be zero */
    if (r == 0)
        return (!check || nmod_poly_mat_is_zero(pmat));

    if (r > n)
        return 0;

    slong * s = FLINT_ARRAY_ALLOC(n, slong);
    for (slong k = 0; k < n; k++)
        s[k] = shift ? shift[k] : 0;

    nmod_poly_mat_t R;
    nmod_poly_mat_window_init(R, rbas, 0, 0, r, n);

    int res = 1;

    /* d = rdeg_s(R); R must have no zero row */
    slong * d = FLINT_ARRAY_ALLOC(r, slong);
    if ((slong) nmod_poly_mat_row_degree(d, R, s) != r)
        res = 0;

    /* J: r columns such that the s-leading matrix L of R satisfies L[:,J]
     * invertible; L[j,k] = coefficient of degree d[j] - s[k] of R[j,k];
     * inv0 = L[:,J]^{-1} */
    slong * J = FLINT_ARRAY_ALLOC(n, slong);
    nmod_mat_t inv0;
    nmod_mat_init(inv0, r, r, mod);
    if (res)
    {
        nmod_mat_t LT;
        nmod_mat_init(LT, n, r, mod);
        for (slong j = 0; j < r; j++)
            for (slong k = 0; k < n; k++)
            {
                const nmod_poly_struct * f = nmod_poly_mat_entry(R, j, k);
                const slong e = d[j] - s[k];
                if (e >= 0 && e < f->length)
                    nmod_mat_entry(LT, k, j) = f->coeffs[e];
            }
        if (nmod_mat_lu(J, LT, 0) != r)
            res = 0;  /* R is not s-reduced */
        else
        {
            /* LT[J,:] = L[:,J]^T, whose inverse is the transpose of inv0 */
            nmod_mat_t tmp;
            nmod_mat_init(tmp, r, r, mod);
            _inverse_from_lu(tmp, LT, r);
            nmod_mat_transpose(inv0, tmp);
            nmod_mat_clear(tmp);
        }
        nmod_mat_clear(LT);
    }

    /* a = rdeg_s(pmat) on nonzero rows (0 on zero rows); sigma = max(a) -
     * min(d) + 1 */
    slong * a = flint_calloc(m, sizeof(slong));
    char * nz = flint_calloc(m, sizeof(char));
    slong sigma = WORD_MIN;
    slong mind = WORD_MAX;
    for (slong j = 0; res && j < r; j++)
        mind = FLINT_MIN(mind, d[j]);
    for (slong i = 0; res && i < m; i++)
    {
        for (slong k = 0; k < n; k++)
        {
            const slong dk = nmod_poly_degree(nmod_poly_mat_entry(pmat, i, k));
            if (dk >= 0 && (!nz[i] || dk + s[k] > a[i]))
            {
                a[i] = dk + s[k];
                nz[i] = 1;
            }
        }
        if (nz[i])
            sigma = FLINT_MAX(sigma, a[i] - mind + 1);
    }

    if (res && sigma != WORD_MIN)  /* otherwise pmat is zero, and so is C */
    {
        if (sigma <= 0)
            /* C must be zero, which is impossible since pmat is nonzero:
             * pmat is not in the row span of R */
            res = 0;
    }

    if (res && sigma > 0)
    {
        /* reversed matrices, as power series in y = 1/x (see
         * nmod_poly_mat_reverse_trunc), restricted to the columns J:
         *   Rt = diag(x^-d) * R[:,J] * diag(x^s_J)
         *   At = diag(x^-a) * pmat[:,J] * diag(x^s_J)
         * and Ct = diag(x^-a) * C * diag(x^d) satisfies Ct * Rt = At,
         * where Rt(0) = L[:,J] is invertible and Ct[i,j] has degree at most
         * a[i] - d[j] in y; the reversals are truncated at order sigma */
        nmod_poly_mat_t Rt, At, X, Ct;
        nmod_poly_mat_init(Rt, r, r, mod);
        nmod_poly_mat_init(At, m, r, mod);
        nmod_poly_mat_init(X, r, r, mod);
        nmod_poly_mat_init(Ct, m, r, mod);

        {
            slong * ms = FLINT_ARRAY_ALLOC(n, slong);
            for (slong k = 0; k < n; k++)
                ms[k] = -s[k];
            nmod_poly_mat_t Rrev, Arev;
            nmod_poly_mat_init(Rrev, r, n, mod);
            nmod_poly_mat_init(Arev, m, n, mod);
            nmod_poly_mat_reverse_trunc(Rrev, R, d, ms, sigma);
            nmod_poly_mat_reverse_trunc(Arev, pmat, a, ms, sigma);
            for (slong l = 0; l < r; l++)
            {
                for (slong j = 0; j < r; j++)
                    nmod_poly_swap(nmod_poly_mat_entry(Rt, j, l), nmod_poly_mat_entry(Rrev, j, J[l]));
                for (slong i = 0; i < m; i++)
                    nmod_poly_swap(nmod_poly_mat_entry(At, i, l), nmod_poly_mat_entry(Arev, i, J[l]));
            }
            nmod_poly_mat_clear(Rrev);
            nmod_poly_mat_clear(Arev);
            flint_free(ms);
        }

        /* X = Rt^{-1} mod y^sigma, with Rt(0) = L[:,J] */
        _nmod_poly_mat_inv_series(X, Rt, inv0, sigma);
        nmod_poly_mat_mulmid(Ct, At, X, 0, sigma);

        /* C = diag(x^a) * Ct(1/x) * diag(x^-d): C[i,j] is the reverse of
         * Ct[i,j] with respect to the bound a[i] - d[j] (zero if this bound
         * is negative) */
        slong * md = FLINT_ARRAY_ALLOC(r, slong);
        for (slong j = 0; j < r; j++)
            md[j] = -d[j];
        nmod_poly_mat_t C;
        nmod_poly_mat_window_init(C, lfac, 0, 0, m, r);
        nmod_poly_mat_reverse_trunc(C, Ct, a, md, -1);
        nmod_poly_mat_window_clear(C);
        flint_free(md);

        nmod_poly_mat_clear(Rt);
        nmod_poly_mat_clear(At);
        nmod_poly_mat_clear(X);
        nmod_poly_mat_clear(Ct);
    }

    if (res && check)
    {
        nmod_poly_mat_t C, prod;
        nmod_poly_mat_window_init(C, lfac, 0, 0, m, r);
        nmod_poly_mat_init(prod, m, n, mod);
        nmod_poly_mat_mul(prod, C, R);
        res = nmod_poly_mat_equal(prod, pmat);
        nmod_poly_mat_clear(prod);
        nmod_poly_mat_window_clear(C);
    }

    flint_free(nz);
    flint_free(a);
    nmod_mat_clear(inv0);
    flint_free(J);
    flint_free(d);
    nmod_poly_mat_window_clear(R);
    flint_free(s);
    return res;
}
