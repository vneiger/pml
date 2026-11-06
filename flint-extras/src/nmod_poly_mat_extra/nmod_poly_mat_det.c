/*
    Copyright (C) 2025 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/nmod.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_vec.h>
#include <flint/perm.h>
#include <flint/ulong_extras.h>

#include "nmod_extra.h"  /* for nmod_find_root */
#include "nmod_poly_extra.h"  /* for NMOD_POLY_CAN_USE_GEOMETRIC */
#include "nmod_poly_mat_utils.h"  /* for permute_rows_by_sorting_vec */
#include "nmod_poly_mat_forms.h"
#include "nmod_poly_mat_determinant.h"
#include "impl.h"  /* for _nmod_mat_det_solve_lu, _nmod_poly_mat_evaluate_geometric */

// Mulders-Storjohann triangularization, see the documentation in
// nmod_poly_mat_determinant.h
slong _nmod_poly_mat_upper_triangularize_iter(nmod_poly_struct * diag, nmod_poly_mat_t mat)
{
    const slong n = mat->r;
    if (n == 0)
        return 1;

    // determinant (+1 or -1) of the unimodular transformation, 0 if singular
    slong udet = 1;
    slong * pivind = flint_malloc(n * sizeof(slong));
    // permutation for putting into ordered weak Popov
    slong * perm = flint_malloc(n * sizeof(slong));

    // window of original mat
    nmod_poly_mat_t view;
    nmod_poly_mat_window_init(view, mat, 0, 0, n, n);

    for (slong i = n - 1; i >= 1; i--)
    {
        //            [ V  * ]                 [ V ]
        // mat is now [ 0  D ] and view is now [ 0 ] of size n x (i+1),
        // with V of size (i+1) x (i+1)
        //                 [ V'  *  * ]
        // -> transform to [ 0   d  * ] via weak Popov form of V[:,:i]
        //                 [ 0   0  D ]
        // (weak Popov form of view[:i+1,:i], with transformations applied to
        // the whole rows view[:i+1,:] = V, early exit if rank < i)
        const slong rk = _nmod_poly_mat_weak_popov_iter_submat_rowbyrow(view, NULL, NULL, &udet, pivind, NULL, 0, 0, i+1, i, 2, ROW_UPPER);

        if (rk < i || nmod_poly_is_zero(nmod_poly_mat_entry(view, i, i)))
        {
            udet = 0;
            break;
        }

        // permute into ordered weak Popov form
        _nmod_poly_mat_permute_rows_by_sorting_vec(view, rk, pivind, perm);
        _nmod_poly_mat_window_resize_columns(view, -1);
        if (_perm_parity(perm, rk)) // odd permutation, negate udet
            udet = -udet;
    }

    // read the diagonal through view rather than mat, since with FLINT < 3.3
    // the rows of mat have not been permuted
    _nmod_poly_mat_window_resize_columns(view, n - view->c);
    if (udet != 0 && nmod_poly_is_zero(nmod_poly_mat_entry(view, 0, 0)))
        udet = 0;
    if (udet != 0)
        for (slong i = 0; i < n; i++)
            nmod_poly_swap(diag + i, nmod_poly_mat_entry(view, i, i));

    flint_free(pivind);
    flint_free(perm);
    nmod_poly_mat_window_clear(view);
    return udet;
}

// Mulders-Storjohann determinant algorithm
// matrix must be square, not tested; this modifies mat (see
// _nmod_poly_mat_upper_triangularize_iter)
void nmod_poly_mat_det_iter(nmod_poly_t det, nmod_poly_mat_t mat)
{
    const slong n = mat->r;
    nmod_poly_struct * diag = flint_malloc(FLINT_MAX(1, n) * sizeof(nmod_poly_struct));
    for (slong i = 0; i < n; i++)
        nmod_poly_init(diag + i, mat->modulus);

    const slong udet = _nmod_poly_mat_upper_triangularize_iter(diag, mat);
    if (udet == 0)
        nmod_poly_zero(det);
    else if (n == 0)
        nmod_poly_one(det);
    else
    {
        /* determinant = udet * product of the diagonal entries */
        nmod_poly_swap(det, diag + 0);
        for (slong i = 1; i < n; i++)
            nmod_poly_mul(det, det, diag + i);
        if (udet == -1)
            nmod_poly_neg(det, det);
    }

    for (slong i = 0; i < n; i++)
        nmod_poly_clear(diag + i);
    flint_free(diag);
}

/*------------------------------------------------------------*/
/* expansion by minors                                        */
/*------------------------------------------------------------*/

#define _ENTRY(i, j) (nmod_poly_mat_entry(pmat, (i), (j)))

/* d2 = det of the 2 x 2 submatrix of pmat at rows i0, i1 and columns j0, j1,
 * tmp is a temporary */
static void _det22(nmod_poly_t d2, const nmod_poly_mat_t pmat,
                   slong i0, slong i1, slong j0, slong j1, nmod_poly_t tmp)
{
    nmod_poly_mul(d2, _ENTRY(i0, j0), _ENTRY(i1, j1));
    nmod_poly_mul(tmp, _ENTRY(i0, j1), _ENTRY(i1, j0));
    nmod_poly_sub(d2, d2, tmp);
}

/* explicit formulas for dimension 3, as in the NTL-based part of PML:
 * expansion along the last row, with the 2 x 2 minors of the first two rows */
static void _det_minors_3(nmod_poly_t det, const nmod_poly_mat_t pmat)
{
    nmod_poly_t buf, det22;
    nmod_poly_init_mod(buf, det->mod);
    nmod_poly_init_mod(det22, det->mod);

    /* det = pmat[2][0] * det(pmat[0,1][1,2]) */
    _det22(det22, pmat, 0, 1, 1, 2, buf);
    nmod_poly_mul(det, _ENTRY(2, 0), det22);

    /* det -= pmat[2][1] * det(pmat[0,1][0,2]) */
    _det22(det22, pmat, 0, 1, 0, 2, buf);
    nmod_poly_mul(buf, _ENTRY(2, 1), det22);
    nmod_poly_sub(det, det, buf);

    /* det += pmat[2][2] * det(pmat[0,1][0,1]) */
    _det22(det22, pmat, 0, 1, 0, 1, buf);
    nmod_poly_mul(buf, _ENTRY(2, 2), det22);
    nmod_poly_add(det, det, buf);

    nmod_poly_clear(buf);
    nmod_poly_clear(det22);
}

/* explicit formulas for dimension 4, as in the NTL-based part of PML: the six
 * 2 x 2 minors of the first two rows are computed once, the four 3 x 3
 * minors of the first three rows are deduced, and then the expansion along
 * the last row */
static void _det_minors_4(nmod_poly_t det, const nmod_poly_mat_t pmat)
{
    nmod_poly_t buf, det33;
    nmod_poly_struct d22[6];
    nmod_poly_init_mod(buf, det->mod);
    nmod_poly_init_mod(det33, det->mod);
    for (slong i = 0; i < 6; i++)
        nmod_poly_init_mod(d22 + i, det->mod);

    /* 2 x 2 minors of rows 0,1, columns [0,1], [0,2], [0,3], [1,2], [1,3], [2,3] */
    _det22(d22 + 0, pmat, 0, 1, 0, 1, buf);
    _det22(d22 + 1, pmat, 0, 1, 0, 2, buf);
    _det22(d22 + 2, pmat, 0, 1, 0, 3, buf);
    _det22(d22 + 3, pmat, 0, 1, 1, 2, buf);
    _det22(d22 + 4, pmat, 0, 1, 1, 3, buf);
    _det22(d22 + 5, pmat, 0, 1, 2, 3, buf);

    /* det = - pmat[3][0] * det(pmat[0,1,2][1,2,3]) */
    nmod_poly_mul(det33, _ENTRY(2, 1), d22 + 5);
    nmod_poly_mul(buf, _ENTRY(2, 2), d22 + 4);
    nmod_poly_sub(det33, det33, buf);
    nmod_poly_mul(buf, _ENTRY(2, 3), d22 + 3);
    nmod_poly_add(det33, det33, buf);
    nmod_poly_mul(det, _ENTRY(3, 0), det33);
    nmod_poly_neg(det, det);

    /* det += pmat[3][1] * det(pmat[0,1,2][0,2,3]) */
    nmod_poly_mul(det33, _ENTRY(2, 0), d22 + 5);
    nmod_poly_mul(buf, _ENTRY(2, 2), d22 + 2);
    nmod_poly_sub(det33, det33, buf);
    nmod_poly_mul(buf, _ENTRY(2, 3), d22 + 1);
    nmod_poly_add(det33, det33, buf);
    nmod_poly_mul(buf, _ENTRY(3, 1), det33);
    nmod_poly_add(det, det, buf);

    /* det -= pmat[3][2] * det(pmat[0,1,2][0,1,3]) */
    nmod_poly_mul(det33, _ENTRY(2, 0), d22 + 4);
    nmod_poly_mul(buf, _ENTRY(2, 1), d22 + 2);
    nmod_poly_sub(det33, det33, buf);
    nmod_poly_mul(buf, _ENTRY(2, 3), d22 + 0);
    nmod_poly_add(det33, det33, buf);
    nmod_poly_mul(buf, _ENTRY(3, 2), det33);
    nmod_poly_sub(det, det, buf);

    /* det += pmat[3][3] * det(pmat[0,1,2][0,1,2]) */
    nmod_poly_mul(det33, _ENTRY(2, 0), d22 + 3);
    nmod_poly_mul(buf, _ENTRY(2, 1), d22 + 1);
    nmod_poly_sub(det33, det33, buf);
    nmod_poly_mul(buf, _ENTRY(2, 2), d22 + 0);
    nmod_poly_add(det33, det33, buf);
    nmod_poly_mul(buf, _ENTRY(3, 3), det33);
    nmod_poly_add(det, det, buf);

    nmod_poly_clear(buf);
    nmod_poly_clear(det33);
    for (slong i = 0; i < 6; i++)
        nmod_poly_clear(d22 + i);
}

/* general dimension n: minors[S] = det(pmat[0:k, S]) for all subsets S of
 * the columns with |S| = k, computed for k = 1, ..., n by expansion along
 * row k-1:
 *   det(pmat[0:k, S]) = sum_{j in S} (-1)^(k-1+pos(j,S)) pmat[k-1][j] det(pmat[0:k-1, S \ {j}])
 * where pos(j,S) is the number of elements of S smaller than j. Subsets are
 * represented by bitmasks; minors of size k-1 are cleared once those of size
 * k are known. */
static void _det_minors_subsets(nmod_poly_t det, const nmod_poly_mat_t pmat)
{
    const slong n = pmat->r;
    const ulong nsub = UWORD(1) << n;

    nmod_poly_struct * minors = flint_malloc(nsub * sizeof(nmod_poly_struct));
    nmod_poly_t tmp;
    nmod_poly_init_mod(tmp, det->mod);

    /* k == 1 */
    for (slong j = 0; j < n; j++)
    {
        const ulong S = UWORD(1) << j;
        nmod_poly_init_mod(minors + S, det->mod);
        nmod_poly_set(minors + S, _ENTRY(0, j));
    }

    for (slong k = 2; k <= n; k++)
    {
        /* subsets of size k: iterate over all masks with popcount k (Gosper's hack) */
        ulong S = (UWORD(1) << k) - 1;
        while (S < nsub)
        {
            nmod_poly_init_mod(minors + S, det->mod);
            slong pos = 0;
            for (ulong T = S; T != 0; T &= T - 1)
            {
                const slong j = flint_ctz(T);
                const ulong Sj = S & ~(UWORD(1) << j);
                if (!nmod_poly_is_zero(_ENTRY(k-1, j)) && !nmod_poly_is_zero(minors + Sj))
                {
                    nmod_poly_mul(tmp, _ENTRY(k-1, j), minors + Sj);
                    if ((k - 1 + pos) % 2 == 0)
                        nmod_poly_add(minors + S, minors + S, tmp);
                    else
                        nmod_poly_sub(minors + S, minors + S, tmp);
                }
                pos++;
            }
            /* next mask with the same popcount */
            const ulong c = S & (-S);
            const ulong r = S + c;
            S = (((r ^ S) >> 2) / c) | r;
        }

        /* clear minors of size k-1 (all of them have been initialized) */
        S = (UWORD(1) << (k-1)) - 1;
        while (S < nsub)
        {
            nmod_poly_clear(minors + S);
            const ulong c = S & (-S);
            const ulong r = S + c;
            S = (((r ^ S) >> 2) / c) | r;
        }
    }

    nmod_poly_swap(det, minors + (nsub - 1));
    nmod_poly_clear(minors + (nsub - 1));
    flint_free(minors);
    nmod_poly_clear(tmp);
}

void nmod_poly_mat_det_minors(nmod_poly_t det, const nmod_poly_mat_t pmat)
{
    const slong n = pmat->r;

    if (pmat->r != pmat->c)
        flint_throw(FLINT_ERROR, "Exception (nmod_poly_mat_det_minors). Non-square matrix.\n");

    if (n == 0)
        nmod_poly_one(det);
    else if (n == 1)
        nmod_poly_set(det, _ENTRY(0, 0));
    else if (n == 2)
    {
        nmod_poly_t tmp;
        nmod_poly_init_mod(tmp, det->mod);
        _det22(det, pmat, 0, 1, 0, 1, tmp);
        nmod_poly_clear(tmp);
    }
    else if (n == 3)
        _det_minors_3(det, pmat);
    else if (n == 4)
        _det_minors_4(det, pmat);
    else if (n <= 16)
        _det_minors_subsets(det, pmat);
    else
        flint_throw(FLINT_ERROR, "Exception (nmod_poly_mat_det_minors). Dimension %wd too large (max 16).\n", n);
}

#undef _ENTRY

/*------------------------------------------------------------*/
/* evaluation at a geometric progression                      */
/*------------------------------------------------------------*/

int nmod_poly_mat_det_geometric(nmod_poly_t det, const nmod_poly_mat_t pmat)
{
    const slong n = pmat->r;

    if (pmat->r != pmat->c)
        flint_throw(FLINT_ERROR, "Exception (nmod_poly_mat_det_geometric). Non-square matrix.\n");

    if (n == 0)
    {
        nmod_poly_one(det);
        return 1;
    }

    /* degree bound min(|rdeg|, |cdeg|); zero row or column means det == 0 */
    slong sum_rdeg, sum_cdeg;
    const ulong nzr = nmod_poly_mat_row_degree_sum(&sum_rdeg, pmat, NULL);
    const ulong nzc = nmod_poly_mat_column_degree_sum(&sum_cdeg, pmat, NULL);
    if (nzr < (ulong) n || nzc < (ulong) n)
    {
        nmod_poly_zero(det);
        return 1;
    }

    const slong len = FLINT_MIN(sum_rdeg, sum_cdeg) + 1;
    const ulong p = pmat->modulus;
    if (!NMOD_POLY_CAN_USE_GEOMETRIC(p, len) || !n_is_prime(p))
        return 0;

    nmod_t mod;
    nmod_init(&mod, p);
    const ulong w = nmod_find_root(2*len, mod);
    if (w == 0)
        return 0;

    /* note: each entry has length at most len, since its degree is at most
     * the row degree and column degree it belongs to. The evaluations are
     * done by chunks of c points (see _nmod_poly_mat_evaluate_geometric),
     * with Gc for the first c points: the n^2 evaluations of a chunk share
     * the stored transform for inputs of length at most maxlen (shorter than
     * for length c). G, for the len points, serves the single interpolation,
     * which does not need a stored transform. */
    const slong maxlen = nmod_poly_mat_max_length(pmat);
    const slong c = _nmod_poly_mat_evaluate_geometric_chunk_length(pmat, len, 0);
    nmod_geometric_progression_t G, Gc;
    _nmod_geometric_progression_init_function(G, w, len, mod, UWORD(2));
    _nmod_geometric_progression_init_function(Gc, w, c, mod, UWORD(1));
    _nmod_geometric_progression_init_fft(Gc, maxlen, 0, 0);

    /* determinants at the points q^k, k < len, with q = w^2, by chunks of
     * points q^(k0 + k), k < cl, with z = q^k0: evals[(i*n + j) * cl + k] =
     * pmat[i][j](q^(k0 + k)); for n > 8, in-place LU by _nmod_mat_det
     * (nmod_mat_det would copy the matrix and test the primality of p at
     * each point; for n <= 8 it uses explicit formulas or the
     * division-free Berkowitz algorithm, which are faster) */
    nn_ptr evals = _nmod_vec_init(n * n * c);
    nn_ptr dets = _nmod_vec_init(len);
    const ulong q = nmod_mul(w, w, mod);
    const ulong qc = nmod_pow_ui(q, c, mod);
    ulong z = 1;
    nmod_mat_t M;
    nmod_mat_init(M, n, n, p);
    for (slong k0 = 0; k0 < len; k0 += c)
    {
        const slong cl = FLINT_MIN(c, len - k0);
        _nmod_poly_mat_evaluate_geometric(evals, pmat, z, Gc, cl);
        for (slong k = 0; k < cl; k++)
        {
            for (slong i = 0; i < n; i++)
                for (slong j = 0; j < n; j++)
                    nmod_mat_entry(M, i, j) = evals[(i * n + j) * cl + k];
            dets[k0 + k] = (n <= 8) ? nmod_mat_det(M) : _nmod_mat_det(M);
        }
        z = nmod_mul(z, qc, mod);
    }
    nmod_mat_clear(M);
    _nmod_vec_clear(evals);

    /* interpolation */
    nmod_poly_interpolate_geometric_nmod_vec_fast_precomp(det, dets, G, len);

    _nmod_vec_clear(dets);
    nmod_geometric_progression_clear(G);
    nmod_geometric_progression_clear(Gc);
    return 1;
}

/*------------------------------------------------------------*/
/* general interface                                          */
/*------------------------------------------------------------*/

int _nmod_mat_det_solve_lu(ulong * c, nmod_mat_t X, nmod_mat_t A, const nmod_mat_t B,
                           slong * perm, nmod_mat_t PB)
{
    const slong n = A->r;
    if (nmod_mat_lu(perm, A, 1) != n)
        return 0;

    /* P A = L U: det(A) = +- prod(U[i][i]), and A^{-1} B = U^{-1} L^{-1} P B */
    ulong dA = 1;
    for (slong i = 0; i < n; i++)
        dA = nmod_mul(dA, nmod_mat_entry(A, i, i), A->mod);
    *c = _perm_parity(perm, n) ? nmod_neg(dA, A->mod) : dA;

    for (slong i = 0; i < n; i++)
        _nmod_vec_set(nmod_mat_entry_ptr(PB, i, 0), nmod_mat_entry_ptr(B, perm[i], 0), B->c);
    nmod_mat_solve_tril(X, A, PB, 1);
    nmod_mat_solve_triu(X, A, X, 0);
    return 1;
}

/* If the square matrix A is invertible, sets c = det(A) and X = A^{-1} B
 * and returns 1; otherwise returns 0 (c and X undefined). The modulus must
 * be prime. This is nmod_mat_det(A) followed by nmod_mat_solve(X, A, B),
 * with a single LU factorization. */
static int _det_and_solve(ulong * c, nmod_mat_t X, const nmod_mat_t A, const nmod_mat_t B)
{
    const slong n = A->r;
    nmod_mat_t LU, PB;
    nmod_mat_init_set(LU, A);
    nmod_mat_init(PB, n, B->c, A->mod.n);
    slong * perm = _perm_init(n);

    const int ok = _nmod_mat_det_solve_lu(c, X, LU, B, perm, PB);

    _perm_clear(perm);
    nmod_mat_clear(PB);
    nmod_mat_clear(LU);
    return ok;
}

/* Degree 1, pmat = x A + B: if A is invertible, det(pmat) = det(A) *
 * charpoly(-A^{-1} B); otherwise if B is invertible, det(pmat) = det(B) *
 * x^n * charpoly(-B^{-1} A)(1/x). Returns 0 (det unchanged) if both A and B
 * are singular or if the modulus is not prime. There is deliberately no
 * shortcut for A == I: callers who know that pmat is a characteristic matrix
 * call a characteristic polynomial function directly. */
/* TODO well-made degree-1 path, see the todo of nmod_poly_mat_determinant
 * about _nmod_mat_charpoly_lnz_block_krylov */
static int _det_degree1(nmod_poly_t det, const nmod_poly_mat_t pmat)
{
    const slong n = pmat->r;
    const ulong p = pmat->modulus;

    if (!n_is_prime(p))
        return 0;

    nmod_mat_t A, B, X;
    nmod_mat_init(A, n, n, p);
    nmod_mat_init(B, n, n, p);
    nmod_mat_init(X, n, n, p);
    nmod_poly_mat_get_coeff_mat(A, pmat, 1);
    nmod_poly_mat_get_coeff_mat(B, pmat, 0);

    int ok = 1;
    ulong c = 1;
    if (_det_and_solve(&c, X, A, B))
    {
        nmod_mat_neg(X, X);
        nmod_mat_charpoly(det, X);
    }
    else if (_det_and_solve(&c, X, B, A))
    {
        nmod_mat_neg(X, X);
        nmod_mat_charpoly(det, X);
        nmod_poly_reverse(det, det, n + 1);
    }
    else
        ok = 0;

    if (ok && c != 1)
        nmod_poly_scalar_mul_nmod(det, det, c);

    nmod_mat_clear(A);
    nmod_mat_clear(B);
    nmod_mat_clear(X);
    return ok;
}

void nmod_poly_mat_determinant(nmod_poly_t det, const nmod_poly_mat_t pmat)
{
    const slong n = pmat->r;

    if (pmat->r != pmat->c)
        flint_throw(FLINT_ERROR, "Exception (nmod_poly_mat_determinant). Non-square matrix.\n");

    if (n <= 5)
    {
        nmod_poly_mat_det_minors(det, pmat);
        return;
    }

    const slong d = nmod_poly_mat_degree(pmat);
    if (d < 0)
    {
        nmod_poly_zero(det);
        return;
    }

    /* constant matrix */
    if (d == 0)
    {
        nmod_mat_t M;
        nmod_mat_init(M, n, n, pmat->modulus);
        nmod_poly_mat_get_coeff_mat(M, pmat, 0);
        const ulong c = nmod_mat_det(M);
        nmod_mat_clear(M);
        nmod_poly_zero(det);
        nmod_poly_set_coeff_ui(det, 0, c);
        return;
    }

    /* zero row or column means det == 0; D = min(|rdeg|, |cdeg|) bounds the
     * degree of the determinant */
    slong sum_rdeg, sum_cdeg;
    const ulong nzr = nmod_poly_mat_row_degree_sum(&sum_rdeg, pmat, NULL);
    const ulong nzc = nmod_poly_mat_column_degree_sum(&sum_cdeg, pmat, NULL);
    if (nzr < (ulong) n || nzc < (ulong) n)
    {
        nmod_poly_zero(det);
        return;
    }
    const slong D = FLINT_MIN(sum_rdeg, sum_cdeg);

    /* degree 1, with invertible leading or constant coefficient:
     * characteristic polynomial of a constant matrix */
    if (d == 1 && _det_degree1(det, pmat))
        return;

    /* the algorithms below work on rows (LNZ uses rdeg(pmat) as a shift):
     * use the transpose if it has a smaller sum of row degrees */
    if (sum_cdeg < sum_rdeg)
    {
        nmod_poly_mat_t T;
        nmod_poly_mat_init(T, n, n, pmat->modulus);
        nmod_poly_mat_transpose(T, pmat);
        nmod_poly_mat_determinant(det, T);
        nmod_poly_mat_clear(T);
        return;
    }

    /* strongly unbalanced degrees (e.g. a few rows of much larger degree):
     * the cost of Mulders-Storjohann follows D, whereas that of LNZ follows
     * n*d since its approximation orders depend on the maximum row degree
     * (see the todo at ::nmod_poly_mat_det_lnz); the ratio 3 is the observed
     * crossover between the two */
    /* TODO revisit when LNZ handles unbalanced degrees (partial
     * linearization), and for dimensions beyond 256 (not measured) */
    const int unbalanced = (n * d > 3 * D);

    /* evaluation / interpolation: small dimension, large degree, balanced
     * degrees, provided the n^2 (D+1) stored evaluations fit in 256 MiB */
    if (!unbalanced
        && n <= 24 && (d >= 100 || (n <= 12 && d >= 30) || (n <= 8 && d >= 20))
        && (ulong) n * n * (D + 1) <= (UWORD(1) << 25)
        && nmod_poly_mat_det_geometric(det, pmat))
        return;

    /* Mulders-Storjohann: unbalanced degrees, or small dimension and degree */
    if (unbalanced || (n <= 16 && n * d <= 320) || (d <= 1 && n <= 32))
    {
        nmod_poly_mat_t copy;
        nmod_poly_mat_init_set(copy, pmat);
        nmod_poly_mat_det_iter(det, copy);
        nmod_poly_mat_clear(copy);
        return;
    }

    /* Labahn-Neiger-Zhou triangularization step, then recursively */
    _nmod_poly_mat_det_lnz(det, pmat, 1);
}
