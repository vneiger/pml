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

/*------------------------------------------------------------*/
/* helpers                                                    */
/*------------------------------------------------------------*/

/* pmat (m x n) has full row rank, hence m <= n: it is its own row basis,
 * with left factor the identity; rbas (at least m rows) and lfac (if not
 * NULL, at least m columns) are assumed to be zero */
static void _set_full_row_rank(nmod_poly_mat_t rbas, nmod_poly_mat_t lfac,
                               const nmod_poly_mat_t pmat)
{
    nmod_poly_mat_t top;
    nmod_poly_mat_window_init(top, rbas, 0, 0, pmat->r, pmat->c);
    nmod_poly_mat_set(top, pmat);
    nmod_poly_mat_window_clear(top);
    if (lfac)
        nmod_poly_mat_one(lfac);
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

    /* first, bound w <= |cdeg_{-s}(pmat)|; a zero column means that pmat is
     * not of full column rank */
    slong * ms = FLINT_ARRAY_ALLOC(m, slong);
    for (slong i = 0; i < m; i++)
        ms[i] = -shift[i];
    slong sum_u;
    int cert = 0;
    if (nmod_poly_mat_column_degree_sum(&sum_u, pmat, ms) == (ulong) n)
    {
        cert = (sum_t == sum_s + sum_u);

        /* then, the sharper bound given by the best assignment of the columns
         * to distinct rows (which exists since pmat has full column rank) */
        if (!cert)
            cert = (sum_t == sum_s + nmod_poly_mat_degree_det_bound(NULL, pmat, ms, COL_UPPER));
    }
    flint_free(ms);
    return cert;
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

/* _nmod_poly_mat_row_basis_zl without the full row rank shortcut; rbas and
 * lfac (if not NULL) are assumed to be zero */
static slong _row_basis_zl(nmod_poly_mat_t rbas,
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
        return 0;

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
        /* pmat is zero: rbas and lfac stay zero */
    }

    else if (k == 0)
    {
        /* full row rank: pmat is a row basis of itself (note m <= n) */
        _set_full_row_rank(rbas, lfac, pmat);
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
                /* row space of pmat is K[x]^{1 x n}: R = I and C = pmat */
                nmod_poly_mat_one(rbas);
                if (lfac)
                {
                    nmod_poly_mat_t left;
                    nmod_poly_mat_window_init(left, lfac, 0, 0, m, n);
                    nmod_poly_mat_set(left, pmat);
                    nmod_poly_mat_window_clear(left);
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

            nmod_poly_mat_t win;
            nmod_poly_mat_window_init(win, rbas, 0, 0, r, n);
            nmod_poly_mat_swap_entrywise(win, R);
            nmod_poly_mat_window_clear(win);
            if (lfac)
            {
                nmod_poly_mat_window_init(win, lfac, 0, 0, m, r);
                nmod_poly_mat_swap_entrywise(win, C);
                nmod_poly_mat_window_clear(win);
            }

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

slong _nmod_poly_mat_row_basis_zl(nmod_poly_mat_t rbas,
                                  nmod_poly_mat_t lfac,
                                  const nmod_poly_mat_t pmat,
                                  const slong * shift,
                                  const nmod_poly_mat_t ker,
                                  const slong * ker_rdeg,
                                  slong nullity,
                                  int check_unimodular)
{
    nmod_poly_mat_zero(rbas);
    if (lfac)
        nmod_poly_mat_zero(lfac);

    /* no kernel basis given, m <= n, and full row rank at some point: then
     * pmat has full row rank (the rank at a point is a lower bound on the
     * rank), and computing a kernel basis is avoided */
    if (ker == NULL && pmat->r > 0 && pmat->r <= pmat->c
        && nmod_poly_mat_rank_at_some_point(pmat) == pmat->r)
    {
        _set_full_row_rank(rbas, lfac, pmat);
        return pmat->r;
    }

    return _row_basis_zl(rbas, lfac, pmat, shift, ker, ker_rdeg, nullity, check_unimodular);
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

    nmod_poly_mat_zero(rbas);
    if (m == 0 || n == 0)
        return 0;

    nmod_poly_mat_t copy;
    nmod_poly_mat_init_set(copy, pmat);
    slong * pivind = FLINT_ARRAY_ALLOC(m, slong);
    /* zero rows are moved to the bottom */
    const slong r = nmod_poly_mat_weak_popov_iter(copy, NULL, NULL, pivind, NULL, ROW_LOWER);
    nmod_poly_mat_t top, ctop;
    nmod_poly_mat_window_init(top, rbas, 0, 0, r, n);
    nmod_poly_mat_window_init(ctop, copy, 0, 0, r, n);
    nmod_poly_mat_swap_entrywise(top, ctop);
    nmod_poly_mat_window_clear(ctop);
    nmod_poly_mat_window_clear(top);
    flint_free(pivind);
    nmod_poly_mat_clear(copy);
    return r;
}

slong _nmod_poly_mat_row_basis_zl_blocks(nmod_poly_mat_t rbas,
                                         const nmod_poly_mat_t pmat,
                                         slong bsize)
{
    const slong m = pmat->r;
    const slong n = pmat->c;

    nmod_poly_mat_zero(rbas);
    if (m == 0 || n == 0)
        return 0;

    if (bsize < 1)
        bsize = 1;
    else if (bsize > m)
        bsize = m;

    /* tall matrix: rows sorted by increasing degree [ZL13, Section 6]; the
     * z = m - nnz zero rows (degree -1) come first, and are skipped */
    nmod_poly_mat_t P;
    nmod_poly_mat_init_set(P, pmat);
    slong * rdeg = FLINT_ARRAY_ALLOC(m, slong);
    slong * perm = FLINT_ARRAY_ALLOC(m, slong);
    const slong nnz = nmod_poly_mat_row_degree(rdeg, P, NULL);
    _nmod_poly_mat_permute_rows_by_sorting_vec(P, m, rdeg, perm);
    flint_free(perm);
    flint_free(rdeg);

    if (nnz == 0)
    {
        nmod_poly_mat_clear(P);
        return 0;
    }

    /* first block: n + bsize rows (or all nonzero rows if fewer) */
    const slong z = m - nnz;
    slong pos = FLINT_MIN(nnz, n + bsize);
    nmod_poly_mat_t W, Rcur, Rtop;
    nmod_poly_mat_window_init(W, P, z, 0, z + pos, n);
    nmod_poly_mat_init(Rcur, n, n, pmat->modulus);
    slong r = _nmod_poly_mat_row_basis_zl(Rcur, NULL, W, NULL, NULL, NULL, 0, 1);
    nmod_poly_mat_window_clear(W);

    /* next blocks: bsize rows at a time, appended to the current row basis */
    while (pos < nnz)
    {
        /* if the row space is already K[x]^{1 x n}, we are done */
        if (r == n && nmod_poly_mat_is_one(Rcur))
            break;

        const slong nb = FLINT_MIN(bsize, nnz - pos);
        nmod_poly_mat_t B;
        nmod_poly_mat_window_init(Rtop, Rcur, 0, 0, r, n);
        nmod_poly_mat_window_init(B, P, z + pos, 0, z + pos + nb, n);
        nmod_poly_mat_init(W, r + nb, n, pmat->modulus);
        nmod_poly_mat_concat_vertical(W, Rtop, B);
        nmod_poly_mat_window_clear(B);
        nmod_poly_mat_window_clear(Rtop);
        r = _nmod_poly_mat_row_basis_zl(Rcur, NULL, W, NULL, NULL, NULL, 0, 1);
        nmod_poly_mat_clear(W);
        pos += nb;
    }

    /* result: first r rows of Rcur (rbas is zero) */
    nmod_poly_mat_t top;
    nmod_poly_mat_window_init(top, rbas, 0, 0, r, n);
    nmod_poly_mat_window_init(Rtop, Rcur, 0, 0, r, n);
    nmod_poly_mat_swap_entrywise(top, Rtop);
    nmod_poly_mat_window_clear(Rtop);
    nmod_poly_mat_window_clear(top);

    nmod_poly_mat_clear(Rcur);
    nmod_poly_mat_clear(P);
    return r;
}

/* For rank deficient inputs, the unimodularity certificate does not apply
 * and the weak Popov approach remains faster up to about
 * ROW_BASIS_RANK_DEFICIENT_FACTOR times the usual threshold (Intel Xeon,
 * AVX-512, 60-bit prime, 2n x n inputs of rank n/2, n = 4 ... 30) */
#define ROW_BASIS_RANK_DEFICIENT_FACTOR 3

slong _nmod_poly_mat_row_basis(nmod_poly_mat_t rbas,
                               nmod_poly_mat_t lfac,
                               const nmod_poly_mat_t pmat,
                               slong thres,
                               slong thres_tall)
{
    const slong m = pmat->r;
    const slong n = pmat->c;
    const slong d = nmod_poly_mat_degree(pmat);

    nmod_poly_mat_zero(rbas);
    if (lfac)
        nmod_poly_mat_zero(lfac);

    if (m == 0 || n == 0 || d < 0)
        return 0;

    const int tall = (m > ROW_BASIS_TALL_RATIO * n);

    /* small degree: weak Popov form, iterative algorithm */
    int weak_popov = (n * d <= (tall ? thres_tall : thres));

    if (!weak_popov && !tall)
    {
        /* the rank at a point is a lower bound on the rank */
        const slong rk_pt = nmod_poly_mat_rank_at_some_point(pmat);

        /* full row rank (m <= n), certified: pmat is its own row basis */
        if (rk_pt == m)
        {
            _set_full_row_rank(rbas, lfac, pmat);
            return m;
        }

        /* moderate degree and (probably) rank deficient: weak Popov form
         * (this only affects the choice of algorithm); otherwise
         * Zhou-Labahn, which also computes the left factor */
        if (n * d <= ROW_BASIS_RANK_DEFICIENT_FACTOR * thres && rk_pt < FLINT_MIN(m, n))
            weak_popov = 1;
        else
            return _row_basis_zl(rbas, lfac, pmat, NULL, NULL, NULL, 0, 1);
    }

    slong r;
    if (weak_popov)
        r = _nmod_poly_mat_row_basis_weak_popov(rbas, pmat);
    else
    {
        /* tall matrix: Zhou-Labahn by blocks, which does not give the left
         * factor; for computing it, the row basis is first reduced (it has r
         * rows and degree at most deg(pmat), this is cheap compared to the
         * blocks) */
        r = _nmod_poly_mat_row_basis_zl_blocks(rbas, pmat, n);
        if (lfac && r > 0)
        {
            nmod_poly_mat_t top, R;
            nmod_poly_mat_window_init(top, rbas, 0, 0, r, n);
            nmod_poly_mat_init_set(R, top);
            nmod_poly_mat_window_clear(top);
            _nmod_poly_mat_row_basis_weak_popov(rbas, R);
            nmod_poly_mat_clear(R);
        }
    }

    /* left factor from the reduced row basis */
    if (lfac && r > 0)
        nmod_poly_mat_row_basis_left_factor(lfac, pmat, rbas, r, NULL, 0);

    return r;
}

slong nmod_poly_mat_row_basis(nmod_poly_mat_t rbas,
                              nmod_poly_mat_t lfac,
                              const nmod_poly_mat_t pmat)
{
    return _nmod_poly_mat_row_basis(rbas, lfac, pmat,
                                    NMOD_POLY_MAT_ROW_BASIS_WEAK_POPOV_THRES,
                                    NMOD_POLY_MAT_ROW_BASIS_WEAK_POPOV_THRES_TALL);
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

    /* a = rdeg_s(pmat), zero rows getting min(s) - 1 (they are zero in the
     * reversals below, whatever their bound); sigma = max(a) - min(d) + 1,
     * where the maximum of a is reached on a nonzero row if pmat is nonzero
     * (sigma stays WORD_MIN if pmat is zero) */
    slong * a = FLINT_ARRAY_ALLOC(m, slong);
    slong sigma = WORD_MIN;
    if (res && nmod_poly_mat_row_degree(a, pmat, s) > 0)
    {
        slong maxa = a[0];
        for (slong i = 1; i < m; i++)
            maxa = FLINT_MAX(maxa, a[i]);
        slong mind = d[0];
        for (slong j = 1; j < r; j++)
            mind = FLINT_MIN(mind, d[j]);
        sigma = maxa - mind + 1;
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

    flint_free(a);
    nmod_mat_clear(inv0);
    flint_free(J);
    flint_free(d);
    nmod_poly_mat_window_clear(R);
    flint_free(s);
    return res;
}
