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
#include <flint/perm.h>

#include "nmod_poly_mat_forms.h"
#include "nmod_poly_mat_kernel.h"
#include "nmod_poly_mat_multiply.h"
#include "nmod_poly_mat_row_basis.h"
#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_determinant.h"  /* for _nmod_poly_mat_upper_triangularize_iter */

/* This file implements the triangularization of
 *   [LNZ17] George Labahn, Vincent Neiger, Wei Zhou -- "Fast, deterministic
 *   computation of the Hermite normal form and determinant of a polynomial
 *   matrix", Journal of Complexity 42, 2017 -- arxiv.org/abs/1607.04176
 * in a row-wise version suited to the lower row-wise Hermite form
 * (ROW_LOWER), with a deterministic certificate which avoids row basis
 * computations in the generic case; this triangularization yields the
 * diagonal entries of the Hermite form [LNZ17, Algorithm 1] and the
 * determinant [LNZ17, Algorithm 2]. */

/*------------------------------------------------------------*/
/* one step of the triangularization                          */
/*------------------------------------------------------------*/

/* Input: F square n x n, n >= 2. Write F = [F_L | F_R], with F_L of
 * dimensions n x k with k = ceil(n/2), and F_R of dimensions n x m, m = n-k.
 *
 * This computes s = rdeg(F) (length n), an s-minimal left kernel basis K of
 * F_R in ker[0:k,:] (ker must be initialized n x n), its s-row degree t
 * (length n, only the first k entries are meaningful), and G1 = K * F_L
 * (G1 must be initialized k x k). It returns:
 *   - -1 if F is detected to be singular (zero row, or F_R not of full
 *   column rank); then ker, t, G1 are undefined
 *   - 0 if the row bases of F_R are certified to be unimodular (see
 *   _nmod_poly_mat_row_basis_certify_unimodular)
 *   - 1 otherwise
 *
 * Correctness: for any row basis R of F_R, [K ; W] * F = [[G1, 0], [*, R]]
 * for some W such that [K ; W] is unimodular [LNZ17, Lemma 3.1, transposed];
 * hence the lower row-wise Hermite form of F has the Hermite forms of G1 and R
 * as diagonal blocks.
 */
static int _lnz_kernel_step(nmod_poly_mat_t G1,
                            nmod_poly_mat_t ker,
                            slong * s,
                            slong * t,
                            const nmod_poly_mat_t F)
{
    const slong n = F->r;
    const slong k = (n + 1) / 2;
    const slong m = n - k;

    /* shift s = rdeg(F); a zero row means F singular */
    if (nmod_poly_mat_row_degree(s, F, NULL) < (ulong) n)
        return -1;

    /* s-minimal left kernel basis of F_R, in ker[0:k,:], t = rdeg_s(K)
     * (on a copy of F_R, which nmod_poly_mat_kernel_zls_approx modifies) */
    nmod_poly_mat_t FRw, FR;
    nmod_poly_mat_window_init(FRw, F, 0, k, n, n);
    nmod_poly_mat_init(FR, n, m, F->modulus);
    nmod_poly_mat_set(FR, FRw);
    slong * pivind = FLINT_ARRAY_ALLOC(n, slong);
    for (slong i = 0; i < n; i++)
        t[i] = s[i];
    const slong nullity = nmod_poly_mat_kernel_zls_approx(ker, pivind, t, FR);
    flint_free(pivind);
    nmod_poly_mat_clear(FR);

    /* F_R is rank deficient */
    if (nullity != k)
    {
        nmod_poly_mat_window_clear(FRw);
        return -1;
    }

    /* G1 = K * F_L */
    nmod_poly_mat_t K, FL;
    nmod_poly_mat_window_init(K, ker, 0, 0, k, n);
    nmod_poly_mat_window_init(FL, F, 0, 0, n, k);
    nmod_poly_mat_multiply(G1, K, FL);
    nmod_poly_mat_window_clear(FL);
    nmod_poly_mat_window_clear(K);

    /* unimodularity certificate */
    slong sum_t = 0;
    for (slong i = 0; i < k; i++)
        sum_t += t[i];
    const int cert = _nmod_poly_mat_row_basis_certify_unimodular(FRw, s, sum_t);
    nmod_poly_mat_window_clear(FRw);

    return cert ? 0 : 1;
}

/*------------------------------------------------------------*/
/* Hermite diagonal: recursion                                */
/*------------------------------------------------------------*/

/* Input: F square n x n, n >= 2, with notation as in _lnz_kernel_step.
 * Output: same return value as _lnz_kernel_step, G1 = K * F_L, and if the
 * return value is 1, R a row basis of F_R (R must be initialized m x m).
 * The row basis is computed by the iterative weak Popov form algorithm when
 * m * deg(F_R) <= wp_thres, and by the Zhou-Labahn approach reusing K
 * otherwise. */
static int _hnf_diagonal_lnz_step(nmod_poly_mat_t G1,
                                  nmod_poly_mat_t R,
                                  const nmod_poly_mat_t F,
                                  slong wp_thres)
{
    const slong n = F->r;
    const slong k = (n + 1) / 2;
    const slong m = n - k;

    nmod_poly_mat_t ker;
    nmod_poly_mat_init(ker, n, n, F->modulus);
    slong * s = FLINT_ARRAY_ALLOC(n, slong);
    slong * t = FLINT_ARRAY_ALLOC(n, slong);

    const int res = _lnz_kernel_step(G1, ker, s, t, F);

    if (res == 1)
    {
        nmod_poly_mat_t FRw;
        nmod_poly_mat_window_init(FRw, F, 0, k, n, n);
        if (wp_thres >= 0 && m * nmod_poly_mat_degree(FRw) <= wp_thres)
            _nmod_poly_mat_row_basis_weak_popov(R, FRw);
        else
            _nmod_poly_mat_row_basis_zl(R, NULL, FRw, s, ker, t, k, 0);
        nmod_poly_mat_window_clear(FRw);
    }

    nmod_poly_mat_clear(ker);
    flint_free(t);
    flint_free(s);
    return res;
}

/* base case n == 2: F = [[a, b], [c, d]], diag = ((ad-bc)/g, g) up to
 * normalization, where g = gcd(b, d) */
static int _hnf_diagonal_2x2(nmod_poly_struct * diag, const nmod_poly_mat_t F)
{
    nmod_poly_t det;
    nmod_poly_init(det, F->modulus);
    nmod_poly_mat_det(det, F);

    int res = 0;
    if (!nmod_poly_is_zero(det))
    {
        /* g == gcd(b, d) is nonzero since det is nonzero; it is monic */
        nmod_poly_gcd(diag + 1, nmod_poly_mat_entry(F, 0, 1), nmod_poly_mat_entry(F, 1, 1));
        nmod_poly_divexact(diag + 0, det, diag + 1);
        nmod_poly_make_monic(diag + 0, diag + 0);
        res = 1;
    }

    nmod_poly_clear(det);
    return res;
}

/* Base case for small dimension: Mulders-Storjohann triangularization (see
 * _nmod_poly_mat_upper_triangularize_iter) applied to J F J, where J is the
 * reversal permutation. All transformations are unimodular, so that the diagonal
 * entries of the obtained upper triangular matrix are, up to units, those of
 * the upper row-wise Hermite form of J F J, that is, those of the lower
 * row-wise Hermite form of F in reverse order. */
static int _hnf_diagonal_iter(nmod_poly_struct * diag, const nmod_poly_mat_t F)
{
    const slong n = F->r;

    nmod_poly_mat_t M;
    nmod_poly_mat_init_set(M, F);
    nmod_poly_mat_invert_rows(M, NULL);
    nmod_poly_mat_invert_columns(M, NULL);

    /* diagonal of the upper triangular form of J F J, in reverse order */
    const int res = (_nmod_poly_mat_upper_triangularize_iter(diag, M) != 0);
    if (res)
        for (slong i = 0; i < n / 2; i++)
            nmod_poly_swap(diag + i, diag + n - 1 - i);
    for (slong i = 0; res && i < n; i++)
        nmod_poly_make_monic(diag + i, diag + i);

    nmod_poly_mat_clear(M);
    return res;
}

static int _hnf_diagonal_lnz_rec(nmod_poly_struct * diag, const nmod_poly_mat_t F,
                                 slong wp_thres, slong iter_dim, slong iter_thres)
{
    const slong n = F->r;

    if (n == 0)
        return 1;

    if (n == 1)
    {
        if (nmod_poly_is_zero(nmod_poly_mat_entry(F, 0, 0)))
            return 0;
        nmod_poly_make_monic(diag, nmod_poly_mat_entry(F, 0, 0));
        return 1;
    }

    if (n == 2)
        return _hnf_diagonal_2x2(diag, F);

    /* TODO n * deg uses the maximum degree, which misjudges unbalanced inputs
     * (see the todo at NMOD_POLY_MAT_HNF_DIAGONAL_ITER_DIM): when retuning,
     * consider the sum of the row degrees and revisit the dimension bound */
    if (iter_thres >= 0 && n <= iter_dim
        && n * nmod_poly_mat_degree(F) <= iter_thres)
        return _hnf_diagonal_iter(diag, F);

    const slong k = (n + 1) / 2;
    const slong m = n - k;

    nmod_poly_mat_t G1, R;
    nmod_poly_mat_init(G1, k, k, F->modulus);
    nmod_poly_mat_init(R, m, m, F->modulus);

    int res = _hnf_diagonal_lnz_step(G1, R, F, wp_thres);
    if (res < 0)
        res = 0;
    else
    {
        const int unimodular = (res == 0);
        res = _hnf_diagonal_lnz_rec(diag, G1, wp_thres, iter_dim, iter_thres);
        if (res && unimodular)
            for (slong i = k; i < n; i++)
                nmod_poly_one(diag + i);
        else if (res)
            res = _hnf_diagonal_lnz_rec(diag + k, R, wp_thres, iter_dim, iter_thres);
    }

    nmod_poly_mat_clear(G1);
    nmod_poly_mat_clear(R);
    return res;
}

/*------------------------------------------------------------*/
/* determinant: recursion                                     */
/*------------------------------------------------------------*/

/* Constant factor of the determinant [LNZ17, Lemma 4.4, transposed].
 * With notation as in _lnz_kernel_step, let R be a row basis of F_R and C
 * the left factor such that F_R = C * R. Then det(F) = kappa * det(G1) *
 * det(R), where kappa is the nonzero constant
 *     kappa = det(C(0)[Jc,:]) / det([K(0) ; E_Jc])
 * for J any set of k columns such that K(0)[:,J] is invertible, Jc its
 * complement (increasing order), and E_Jc the m x n matrix whose rows are
 * the unit vectors indexed by Jc. Indeed [K ; W] * F = [[G1, 0], [*, R]]
 * with W such that W * C = I; then [K ; W] is unimodular and det([K ; W])
 * * det([Y | C]) = 1 for any Y such that K Y = I, so that det(F) = det(G1)
 * * det(R) * det([Y | C]), the latter being a constant which can be
 * evaluated at 0, and any Y0 such that K(0) Y0 = I gives det([Y0 | C(0)]) =
 * det(C(0)[Jc,:]) / det([K(0) ; E_Jc]). Both K(0) and C(0) have full rank
 * since K and C are kernel bases.
 *
 * Input: ker as in _lnz_kernel_step (first k rows), C0 = C(0) (n x m
 * constant matrix). Output: kappa. Raises an exception if a full rank
 * assumption is violated (which would indicate a bug). */
static ulong _det_lnz_constant(const nmod_poly_mat_t ker, slong k, const nmod_mat_t C0)
{
    const slong n = ker->c;
    const slong m = n - k;
    const ulong p = ker->modulus;

    nmod_t mod;
    nmod_init(&mod, p);

    /* K(0), and its pivot / nonpivot columns: P K(0) = L U with U in echelon
     * form, whose pivot columns are J = pivnp[0:k] */
    nmod_mat_t K0, CJc;
    nmod_mat_init(K0, k, n, p);
    for (slong i = 0; i < k; i++)
        for (slong j = 0; j < n; j++)
        {
            const nmod_poly_struct * f = nmod_poly_mat_entry(ker, i, j);
            nmod_mat_entry(K0, i, j) = (f->length > 0) ? f->coeffs[0] : 0;
        }
    slong * P = FLINT_ARRAY_ALLOC(FLINT_MAX(1, k), slong);
    slong * pivnp = FLINT_ARRAY_ALLOC(n, slong);
    const slong rk = nmod_mat_lu_with_pivots(P, pivnp, K0);
    if (rk != k)
        flint_throw(FLINT_ERROR, "Exception (nmod_poly_mat_det_lnz). Kernel basis not of full rank at 0.\n");

    /* det([K(0) ; E_Jc]) = sign(pivnp) * det(K(0)[:,J]) (permute the columns
     * as [J | Jc]: this gives [[K(0)[:,J], *], [0, I]]), and det(K(0)[:,J]) =
     * sign(P) * prod_i U[i, J[i]] */
    ulong dM = 1;
    for (slong i = 0; i < k; i++)
        dM = nmod_mul(dM, nmod_mat_entry(K0, i, pivnp[i]), mod);
    if (_perm_parity(P, k) != _perm_parity(pivnp, n))
        dM = nmod_neg(dM, mod);

    /* CJc = C(0)[Jc,:] */
    nmod_mat_init(CJc, m, m, p);
    for (slong i = 0; i < m; i++)
    {
        const slong jc = pivnp[k + i];
        for (slong j = 0; j < m; j++)
            nmod_mat_entry(CJc, i, j) = nmod_mat_entry(C0, jc, j);
    }
    const ulong dC = nmod_mat_det(CJc);
    if (dM == 0 || dC == 0)
        flint_throw(FLINT_ERROR, "Exception (nmod_poly_mat_det_lnz). Unexpected singular constant matrix.\n");

    const ulong kappa = nmod_mul(dC, nmod_inv(dM, mod), mod);

    flint_free(P);
    flint_free(pivnp);
    nmod_mat_clear(K0);
    nmod_mat_clear(CJc);
    return kappa;
}

/* hybrid == 0: recursive calls go to _det_lnz_rec, with expansion by minors
 * in dimension <= 4; hybrid != 0: recursive calls go to the general
 * interface nmod_poly_mat_determinant, which may itself call back this
 * function */
static void _det_lnz_rec(nmod_poly_t det, const nmod_poly_mat_t F, int hybrid)
{
    const slong n = F->r;

    if (n <= 4)
    {
        nmod_poly_mat_det_minors(det, F);
        return;
    }

    const slong k = (n + 1) / 2;
    const slong m = n - k;
    const ulong p = F->modulus;

    nmod_poly_mat_t G1, ker;
    nmod_poly_mat_init(G1, k, k, p);
    nmod_poly_mat_init(ker, n, n, p);
    slong * s = FLINT_ARRAY_ALLOC(n, slong);
    slong * t = FLINT_ARRAY_ALLOC(n, slong);

    const int res = _lnz_kernel_step(G1, ker, s, t, F);

    /* generic case: the constant factor only needs K(0) and F_R(0); compute
     * it now and release the kernel basis before the recursive call */
    ulong kappa0 = 0;
    if (res == 0)
    {
        nmod_poly_mat_t FRw;
        nmod_poly_mat_window_init(FRw, F, 0, k, n, n);
        nmod_mat_t C0;
        nmod_mat_init(C0, n, m, p);
        nmod_poly_mat_get_coeff_mat(C0, FRw, 0);
        kappa0 = _det_lnz_constant(ker, k, C0);
        nmod_mat_clear(C0);
        nmod_poly_mat_window_clear(FRw);
        nmod_poly_mat_clear(ker);
        nmod_poly_mat_init(ker, 0, 0, p);
    }

    if (res < 0)
        nmod_poly_zero(det);

    else
    {
        /* det(G1), first, since F is singular if it is zero */
        if (hybrid)
            nmod_poly_mat_determinant(det, G1);
        else
            _det_lnz_rec(det, G1, 0);
    }
    nmod_poly_mat_clear(G1);
    nmod_poly_mat_init(G1, 0, 0, p);

    if (res == 0)
        nmod_poly_scalar_mul_nmod(det, det, kappa0);

    else if (res > 0 && !nmod_poly_is_zero(det))
    {
        /* row basis R and left factor C, reusing the kernel basis */
        nmod_poly_mat_t FRw, R, C;
        nmod_poly_mat_window_init(FRw, F, 0, k, n, n);
        nmod_poly_mat_init(R, m, m, p);
        nmod_poly_mat_init(C, n, m, p);
        _nmod_poly_mat_row_basis_zl(R, C, FRw, s, ker, t, k, 0);
        nmod_poly_mat_window_clear(FRw);

        nmod_mat_t C0;
        nmod_mat_init(C0, n, m, p);
        nmod_poly_mat_get_coeff_mat(C0, C, 0);
        nmod_poly_mat_clear(C);
        const ulong kappa = _det_lnz_constant(ker, k, C0);
        nmod_mat_clear(C0);
        nmod_poly_mat_clear(ker);
        nmod_poly_mat_init(ker, 0, 0, p);

        nmod_poly_t detR;
        nmod_poly_init(detR, p);
        if (hybrid)
            nmod_poly_mat_determinant(detR, R);
        else
            _det_lnz_rec(detR, R, 0);
        nmod_poly_mat_clear(R);

        nmod_poly_mul(det, det, detR);
        nmod_poly_scalar_mul_nmod(det, det, kappa);
        nmod_poly_clear(detR);
    }

    nmod_poly_mat_clear(G1);
    nmod_poly_mat_clear(ker);
    flint_free(s);
    flint_free(t);
}

/*------------------------------------------------------------*/
/* interfaces                                                 */
/*------------------------------------------------------------*/

int _nmod_poly_mat_hnf_diagonal_lnz(nmod_poly_struct * diag,
                                    const nmod_poly_mat_t pmat,
                                    slong wp_thres,
                                    slong iter_dim,
                                    slong iter_thres)
{
    if (pmat->r != pmat->c)
        flint_throw(FLINT_ERROR, "Exception (nmod_poly_mat_hnf_diagonal_lnz). Non-square matrix.\n");

    return _hnf_diagonal_lnz_rec(diag, pmat, wp_thres, iter_dim, iter_thres);
}

int nmod_poly_mat_hnf_diagonal_lnz(nmod_poly_struct * diag,
                                   const nmod_poly_mat_t pmat)
{
    return _nmod_poly_mat_hnf_diagonal_lnz(diag, pmat,
                                           NMOD_POLY_MAT_ROW_BASIS_WEAK_POPOV_THRES,
                                           NMOD_POLY_MAT_HNF_DIAGONAL_ITER_DIM,
                                           NMOD_POLY_MAT_HNF_DIAGONAL_ITER_NDEG);
}

int nmod_poly_mat_hnf_diagonal(nmod_poly_struct * diag,
                               const nmod_poly_mat_t pmat,
                               orientation_t orient)
{
    const slong n = pmat->r;

    if (orient == ROW_LOWER)
        return nmod_poly_mat_hnf_diagonal_lnz(diag, pmat);

    if (pmat->r != pmat->c)
        flint_throw(FLINT_ERROR, "Exception (nmod_poly_mat_hnf_diagonal). Non-square matrix.\n");

    /* other orientations: transpose and/or reverse rows and columns
     *  ROW_UPPER(F) == reversed ROW_LOWER(J F J)
     *  COL_UPPER(F) == ROW_LOWER(F^T)
     *  COL_LOWER(F) == ROW_UPPER(F^T) == reversed ROW_LOWER(J F^T J)
     * where J is the reversal permutation */
    nmod_poly_mat_t F;
    nmod_poly_mat_init(F, n, n, pmat->modulus);
    if (orient == ROW_UPPER)
        nmod_poly_mat_set(F, pmat);
    else
        nmod_poly_mat_transpose(F, pmat);
    if (orient != COL_UPPER)
    {
        nmod_poly_mat_invert_rows(F, NULL);
        nmod_poly_mat_invert_columns(F, NULL);
    }

    int res = nmod_poly_mat_hnf_diagonal_lnz(diag, F);

    /* reverse the diagonal for ROW_UPPER and COL_LOWER */
    if (res && orient != COL_UPPER)
        for (slong i = 0; i < n/2; i++)
            nmod_poly_swap(diag + i, diag + n - 1 - i);

    nmod_poly_mat_clear(F);
    return res;
}

void _nmod_poly_mat_det_lnz(nmod_poly_t det, const nmod_poly_mat_t pmat, int hybrid)
{
    if (pmat->r != pmat->c)
        flint_throw(FLINT_ERROR, "Exception (nmod_poly_mat_det_lnz). Non-square matrix.\n");

    _det_lnz_rec(det, pmat, hybrid);
}

void nmod_poly_mat_det_lnz(nmod_poly_t det, const nmod_poly_mat_t pmat)
{
    _nmod_poly_mat_det_lnz(det, pmat, 0);
}
