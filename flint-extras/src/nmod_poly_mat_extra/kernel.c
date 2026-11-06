/*
    Copyright (C) 2025 Gilles Villard, Vincent Neiger
    Copyright (C) 2026 Gilles Villard, Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>  /* for qsort */
#include <flint/flint.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>
#include <flint/perm.h>

#include "nmod_poly_mat_extra.h"
#include "nmod_poly_mat_forms.h"
#include "nmod_poly_mat_kernel.h"
#include "nmod_poly_mat_multiply.h"
#include "nmod_poly_mat_utils.h"
#include "impl.h"    /* for _slong_comparator, views */

slong nmod_poly_mat_kernel(nmod_poly_mat_t ker,
                           slong * pivind,
                           slong * shift,
                           const nmod_poly_mat_t pmat,
                           poly_mat_form_t form,
                           orientation_t orient)
{
    if (form > ORD_WEAK_POPOV)
        flint_throw(FLINT_ERROR, "Exception (nmod_poly_mat_kernel). form > ORD_WEAK_POPOV not implemented.");

    slong nullity = -1;

    slong * _pivind = pivind;
    slong * _shift = shift;

    if (_pivind == NULL)
    {
        if (orient == ROW_LOWER || orient == ROW_UPPER)  /* left kernel */
            _pivind = FLINT_ARRAY_ALLOC(pmat->r, slong);
        else  /* right kernel */
            _pivind = FLINT_ARRAY_ALLOC(pmat->c, slong);
    }

    if (_shift == NULL)
    {
        if (orient == ROW_LOWER || orient == ROW_UPPER)  /* left kernel */
        {
            _shift = FLINT_ARRAY_ALLOC(pmat->r, slong);
            nmod_poly_mat_row_degree_zero(_shift, pmat, NULL);
        }
        else  /* right kernel */
        {
            _shift = FLINT_ARRAY_ALLOC(pmat->c, slong);
            nmod_poly_mat_column_degree_zero(_shift, pmat, NULL);
        }
    }

    if (orient == ROW_LOWER)
    {
        nmod_poly_mat_t mat;
        nmod_poly_mat_init_set(mat, pmat);
        nullity = nmod_poly_mat_kernel_zls_approx(ker, _pivind, _shift, mat);
        nmod_poly_mat_clear(mat);
    }

    else if (orient == COL_UPPER)
    {
        /* transpose input, call ROW_LOWER, transpose output */
        nmod_poly_mat_t mat_t;
        nmod_poly_mat_t ker_t;
        nmod_poly_mat_init(mat_t, pmat->c, pmat->r, pmat->modulus);
        nmod_poly_mat_init(ker_t, ker->c, ker->r, ker->modulus);
        nmod_poly_mat_transpose(mat_t, pmat);
        nullity = nmod_poly_mat_kernel_zls_approx(ker_t, _pivind, _shift, mat_t);
        nmod_poly_mat_transpose(ker, ker_t);
        nmod_poly_mat_clear(ker_t);
        nmod_poly_mat_clear(mat_t);
    }

    else if (orient == ROW_UPPER)
    {
        /* mirror rows of input, call ROW_LOWER, mirror columns+rows of output */
        nmod_poly_mat_t mat_i;
        nmod_poly_mat_init_set(mat_i, pmat);
        nmod_poly_mat_invert_rows(mat_i, _shift);
        nullity = nmod_poly_mat_kernel_zls_approx(ker, _pivind, _shift, mat_i);
        nmod_poly_mat_t kernz;
        nmod_poly_mat_window_init(kernz, ker, 0, 0, nullity, ker->c);
        nmod_poly_mat_invert_rows(kernz, _shift);
        nmod_poly_mat_invert_columns(kernz, NULL);
        for (slong i = 0; i < nullity/2; i++)
        {
            slong tmp = ker->c - 1 - _pivind[i];
            _pivind[i] = ker->c - 1 - _pivind[nullity - 1 - i];
            _pivind[nullity - 1 - i] = tmp;
        }
        if (nullity % 2)
            _pivind[nullity/2] = ker->c - 1 - _pivind[nullity/2];
        nmod_poly_mat_window_clear(kernz);
        nmod_poly_mat_clear(mat_i);
    }

    else if (orient == COL_LOWER)
    {
        /* transpose input, call ROW_UPPER, transpose output */
        nmod_poly_mat_t mat_it;
        nmod_poly_mat_t ker_it;
        nmod_poly_mat_init(mat_it, pmat->c, pmat->r, pmat->modulus);
        nmod_poly_mat_init(ker_it, ker->c, ker->r, ker->modulus);
        nmod_poly_mat_transpose(mat_it, pmat);
        nmod_poly_mat_invert_rows(mat_it, _shift);
        
        nullity = nmod_poly_mat_kernel_zls_approx(ker_it, _pivind, _shift, mat_it);
        nmod_poly_mat_t kernz;
        nmod_poly_mat_window_init(kernz, ker_it, 0, 0, nullity, ker_it->c);
        nmod_poly_mat_invert_rows(kernz, _shift);
        nmod_poly_mat_invert_columns(kernz, NULL);
        for (slong i = 0; i < nullity/2; i++)
        {
            slong tmp = ker_it->c - 1 - _pivind[i];
            _pivind[i] = ker_it->c - 1 - _pivind[nullity - 1 - i];
            _pivind[nullity - 1 - i] = tmp;
        }
        if (nullity % 2)
            _pivind[nullity/2] = ker_it->c - 1 - _pivind[nullity/2];
        nmod_poly_mat_window_clear(kernz);

        nmod_poly_mat_transpose(ker, ker_it);
        nmod_poly_mat_clear(ker_it);
        nmod_poly_mat_clear(mat_it);
    }

    if (pivind == NULL)
        flint_free(_pivind);
    if (shift == NULL)
        flint_free(_shift);

    return nullity;
}



slong nmod_poly_mat_kernel_via_approx(nmod_poly_mat_t ker,
                                      slong * pivind,
                                      slong * shift,
                                      const nmod_poly_mat_t pmat)
{
    const slong m = pmat->r;
    const slong n = pmat->c;
    const slong d = nmod_poly_mat_degree(pmat);

    /* pmat == 0 => kernel is identity*/
    if (d == -1)
    {
        nmod_poly_mat_one(ker);
        for (slong i = 0; i < m; i++)
            pivind[i] = i;
        return m;
    }

    /* FIXME may be reasonable to also have a case for d == 0 */

    /* pmat full row rank => empty kernel */
    /* full row rank implied by m == 1 or (m <= n && constant coefficient has full rank) */
    if (m == 1)
        return 0;

    if (m <= n)
    {
        nmod_mat_t coeff0;
        nmod_mat_init(coeff0, m, n, pmat->modulus);
        nmod_poly_mat_get_coeff_mat(coeff0, pmat, 0);
        /* could use: const slong r = nmod_mat_rank(coeff0);                */
        /* but let's call more directly its core computation                */
        /* (this uses pivind as temporary of length m for row permutation)  */
        const slong r = nmod_mat_lu(pivind, coeff0, 1);
        nmod_mat_clear(coeff0);
        if (r == m)
            return 0;
    }

    /* compute amplitude of `shift`, make `pivind` a temporary copy of `shift` */
    slong min = shift[0];
    slong amp = shift[0];
    for (slong i = 0; i < m; i++)
    {
        slong s = shift[i];
        pivind[i] = s;
        if (s < min)
            min = s;
        else if (s > amp)
            amp = s;
    }
    amp = amp - min;

    /* order for approximation: bound on rdeg_s(ker) + deg(pmat) + 1 */
    /* (see notes at end of file for rdeg_s(ker)) */
    const slong order = FLINT_MIN(m, n) * d + amp + d + 1;

    /* compute approximant basis and shifted row degree */
    nmod_poly_mat_pmbasis(ker, pivind, pmat, order);

    /* gather information for rows which belong to the kernel */
    /* note: at this stage, pivind == rdeg_s(appbas) */
    slong nullity = 0;
    for (slong i = 0; i < m; i++)
    {
        if (pivind[i] - shift[i] + amp < order - d)
        {
            shift[nullity] = pivind[i];
            pivind[nullity] = i;
            for (slong j = 0; j < m; j++)
                FLINT_SWAP(nmod_poly_struct,
                           *nmod_poly_mat_entry(ker, nullity, j),
                           *nmod_poly_mat_entry(ker, i, j));
            nullity += 1;
        }
    }

    return nullity;
}

/*------------------------------------------------------------*/
/* products with a few long rows                              */
/*------------------------------------------------------------*/
/* \todo this is temporary, until a more general mechanism is exposed,
 * for multiplication with unbalanced degrees */

/* Splits the rows of F into short and long ones: with med the median of the
 * row degrees of F (zero rows counting as degree 0), a row is long if its
 * degree exceeds 4 med + 4. Writes in idx the indices of the short rows, in
 * increasing order, followed by those of the long rows, and returns the
 * number of long rows. */
static slong _kernel_long_rows(slong * idx, const nmod_poly_mat_t F)
{
    const slong k = F->r;
    if (k == 0)
        return 0;

    slong * rdeg = FLINT_ARRAY_ALLOC(k, slong);
    slong * sorted = FLINT_ARRAY_ALLOC(k, slong);
    nmod_poly_mat_row_degree_zero(rdeg, F, NULL);
    for (slong i = 0; i < k; i++)
        sorted[i] = rdeg[i];
    qsort(sorted, k, sizeof(slong), _slong_comparator);
    const slong thres = 4 * sorted[(k - 1) / 2] + 4;

    slong ns = 0, nl = 0;
    for (slong i = 0; i < k; i++)
        if (rdeg[i] <= thres)
            idx[ns++] = i;
    for (slong i = 0; i < k; i++)
        if (rdeg[i] > thres)
            idx[ns + nl++] = i;

    flint_free(rdeg);
    flint_free(sorted);
    return nl;
}

/* res = P F if nlo < 0, otherwise res = (P F mod x^nhi) div x^nlo */
static void _kernel_mul(nmod_poly_mat_t res, const nmod_poly_mat_t P,
                        const nmod_poly_mat_t F, slong nlo, slong nhi)
{
    if (nlo < 0)
        nmod_poly_mat_multiply(res, P, F);
    else
        nmod_poly_mat_mulmid(res, P, F, nlo, nhi);
}

/* Same as _kernel_mul, with the rows of F split into short and long ones
 * (see _kernel_long_rows): res = Ps Fs + Pl Fl, where Fs, Fl are the short
 * and long rows of F and Ps, Pl the corresponding columns of P. */
static void _kernel_mul_split_right(nmod_poly_mat_t res, const nmod_poly_mat_t P,
                                    const nmod_poly_mat_t F, slong nlo, slong nhi)
{
    const slong k = F->r;
    slong * idx = FLINT_ARRAY_ALLOC(FLINT_MAX(k, 1), slong);
    const slong nl = _kernel_long_rows(idx, F);
    const slong ns = k - nl;

    if (nl == 0 || ns == 0)
        _kernel_mul(res, P, F, nlo, nhi);
    else
    {
        nmod_poly_mat_t Ps, Pl, Fs, Fl, T;
        _nmod_poly_mat_view_init_submatrix(Ps, P, NULL, P->r, idx, ns);
        _nmod_poly_mat_view_init_submatrix(Pl, P, NULL, P->r, idx + ns, nl);
        _nmod_poly_mat_view_init_submatrix(Fs, F, idx, ns, NULL, F->c);
        _nmod_poly_mat_view_init_submatrix(Fl, F, idx + ns, nl, NULL, F->c);
        nmod_poly_mat_init(T, P->r, F->c, P->modulus);
        _kernel_mul(res, Ps, Fs, nlo, nhi);
        _kernel_mul(T, Pl, Fl, nlo, nhi);
        nmod_poly_mat_add(res, res, T);
        nmod_poly_mat_clear(T);
        _nmod_poly_mat_view_clear(Ps);
        _nmod_poly_mat_view_clear(Pl);
        _nmod_poly_mat_view_clear(Fs);
        _nmod_poly_mat_view_clear(Fl);
    }
    flint_free(idx);
}

/* Product res = P F (if nlo < 0) or middle product res = (P F mod x^nhi) div
 * x^nlo (if nlo >= 0), for operands which may have a few rows of much larger
 * degree than the others, as happens in the ZLS algorithm when the input has
 * unbalanced row degrees (e.g. [A ; b] with deg(b) much larger than deg(A),
 * in linear system solving): the product routines work at the largest length
 * of the entries of their operands, so that a single long row would make the
 * whole product as expensive as if all entries were long. Here, the long rows
 * of P are multiplied by F apart from the short ones, and the short rows of P
 * are multiplied by F with the rows of F split in the same way (see
 * _kernel_mul_split_right). Without long rows, this is a single product.
 * res must not alias P or F (it may be a window). */
/* \todo this is temporary, until a more general mechanism is exposed,
 * for multiplication with unbalanced degrees */
static void _kernel_mul_split(nmod_poly_mat_t res, const nmod_poly_mat_t P,
                              const nmod_poly_mat_t F, slong nlo, slong nhi)
{
    const slong r = P->r;
    slong * idx = FLINT_ARRAY_ALLOC(FLINT_MAX(r, 1), slong);
    const slong nl = _kernel_long_rows(idx, P);
    const slong ns = r - nl;

    if (nl == 0 || ns == 0)
        _kernel_mul_split_right(res, P, F, nlo, nhi);
    else
    {
        nmod_poly_mat_t Ps, Pl, Rs, Rl;
        _nmod_poly_mat_view_init_submatrix(Ps, P, idx, ns, NULL, P->c);
        _nmod_poly_mat_view_init_submatrix(Pl, P, idx + ns, nl, NULL, P->c);
        nmod_poly_mat_init(Rs, ns, F->c, P->modulus);
        nmod_poly_mat_init(Rl, nl, F->c, P->modulus);
        _kernel_mul_split_right(Rs, Ps, F, nlo, nhi);
        _kernel_mul(Rl, Pl, F, nlo, nhi);
        for (slong i = 0; i < r; i++)
            for (slong j = 0; j < F->c; j++)
                nmod_poly_swap(nmod_poly_mat_entry(res, idx[i], j),
                               (i < ns) ? nmod_poly_mat_entry(Rs, i, j)
                                        : nmod_poly_mat_entry(Rl, i - ns, j));
        nmod_poly_mat_clear(Rs);
        nmod_poly_mat_clear(Rl);
        _nmod_poly_mat_view_clear(Ps);
        _nmod_poly_mat_view_clear(Pl);
    }
    flint_free(idx);
}

/*------------------------------------------------------------*/
/* zero rows                                                  */
/*------------------------------------------------------------*/

/* Kernel basis of pmat when pmat has at least one zero row and one nonzero
 * row, by removing the zero rows (see [note:zero_rows]); same specification
 * as nmod_poly_mat_kernel_zls_approx, with in addition rdeg = rdeg(pmat), -1
 * on zero rows. With Z and N the indices of the zero and nonzero rows of
 * pmat, the output rows are:
 * - the unit vectors e_i, i in Z, of shifted degree shift[i] and pivot i;
 * - the rows of a shift[N]-ordered weak Popov kernel basis of pmat[N,:]
 *   (recursive call), with zeros inserted at the indices in Z: their shifted
 *   degrees are unchanged, and the pivot index j becomes the j-th element
 *   of N.
 * The pivot indices are pairwise distinct; the rows are sorted by increasing
 * pivot index, which gives a shift-ordered weak Popov kernel basis of pmat.
 * pmat is overwritten: the entries of its nonzero rows are moved out. */
static slong _kernel_zls_zero_rows(nmod_poly_mat_t ker,
                                   slong * pivind,
                                   slong * shift,
                                   nmod_poly_mat_t pmat,
                                   const slong * rdeg)
{
    const slong m = pmat->r;
    const slong n = pmat->c;

    /* idx: the nonzero rows, then the zero rows, in increasing order */
    slong * idx = FLINT_ARRAY_ALLOC(m, slong);
    slong nnz = 0;
    for (slong i = 0; i < m; i++)
        if (rdeg[i] >= 0)
            idx[nnz++] = i;
    for (slong i = 0, k = nnz; i < m; i++)
        if (rdeg[i] < 0)
            idx[k++] = i;

    /* kernel of the nonzero rows */
    nmod_poly_mat_t sub, subker;
    nmod_poly_mat_init(sub, nnz, n, pmat->modulus);
    nmod_poly_mat_init(subker, nnz, nnz, pmat->modulus);
    slong * subshift = FLINT_ARRAY_ALLOC(nnz, slong);
    slong * subpiv = FLINT_ARRAY_ALLOC(nnz, slong);
    for (slong k = 0; k < nnz; k++)
    {
        subshift[k] = shift[idx[k]];
        for (slong j = 0; j < n; j++)
            nmod_poly_swap(nmod_poly_mat_entry(sub, k, j),
                           nmod_poly_mat_entry(pmat, idx[k], j));
    }
    const slong subnull = nmod_poly_mat_kernel_zls_approx(subker, subpiv, subshift, sub);

    /* merge, by increasing pivot index */
    const slong nullity = (m - nnz) + subnull;
    slong * zshift = FLINT_ARRAY_ALLOC(m - nnz, slong);
    for (slong k = 0; k < m - nnz; k++)
        zshift[k] = shift[idx[nnz + k]];
    for (slong r = 0, k1 = 0, k2 = 0; r < nullity; r++)
    {
        for (slong j = 0; j < m; j++)
            nmod_poly_zero(nmod_poly_mat_entry(ker, r, j));
        if (k2 == m - nnz || (k1 < subnull && idx[subpiv[k1]] < idx[nnz + k2]))
        {
            for (slong j = 0; j < nnz; j++)
                nmod_poly_swap(nmod_poly_mat_entry(ker, r, idx[j]),
                               nmod_poly_mat_entry(subker, k1, j));
            pivind[r] = idx[subpiv[k1]];
            shift[r] = subshift[k1];
            k1++;
        }
        else
        {
            nmod_poly_one(nmod_poly_mat_entry(ker, r, idx[nnz + k2]));
            pivind[r] = idx[nnz + k2];
            shift[r] = zshift[k2];
            k2++;
        }
    }

    nmod_poly_mat_clear(sub);
    nmod_poly_mat_clear(subker);
    flint_free(idx);
    flint_free(subshift);
    flint_free(subpiv);
    flint_free(zshift);
    return nullity;
}

/*------------------------------------------------------------*/
/* Zhou-Labahn-Storjohann                                     */
/*------------------------------------------------------------*/

/* Follows the description of Zhou-Labahn-Storjohann algorithm in [LNVZ22, Algo.1]  */
/* [LNVZ22] Labahn-Neiger-Vu-Zhou, Proceedings ISSAC 2022, arxiv.org/pdf/2202.09329 */
slong nmod_poly_mat_kernel_zls_approx(nmod_poly_mat_t ker,
                                      slong * pivind,
                                      slong * shift,
                                      nmod_poly_mat_t pmat)
{
    const slong m = pmat->r;
    const slong n = pmat->c;

    /* empty matrices, or single row/column vectors that are zero */
    if (m == 0 || n == 0
        || ((m == 1 || n == 1) && (nmod_poly_mat_max_length(pmat) == 0)))
    {
        nmod_poly_mat_one(ker);
        for (slong i = 0; i < m; i++)
            pivind[i] = i;
        return m;
    }

    /* one (nonzero) row vector */
    if (m == 1)
        return 0;

    /* if n > m/2 (with here m >= 2), straightforward recursion by splitting */
    /* the matrix into two submatrices with n/2 columns                      */
    /* TODO tune this threshold: for m slightly below 2n (e.g. after zero    */
    /* rows are removed, see [note:zero_rows]), a single approximant basis   */
    /* (case n <= m/2 below, order about 1 + sum(t) / (m - n)) can be faster */
    /* at small degrees. Measured, 60-bit prime, 60 x 32 with row degrees    */
    /* uniform in [0, 9], shift = row degrees: 9.5 ms with the split, 4.7 ms */
    /* when splitting only if n > 2m/3 (56 x 32: 9.1 ms -> 3.4 ms); with row */
    /* degrees up to 99, no difference.                                      */
    if (2*n > m)
    {
        /* we split pmat == [ n - n/2 cols | n/2 cols]                       */
        /* (we want more columns in left part, for the residual computation) */
        slong nullity;

        nmod_poly_mat_t submat;    /* window only */
        nmod_poly_mat_t residual;  /* window only */
        nmod_poly_mat_t ker1;      /* window only */
        nmod_poly_mat_t ker2;      /* actual init */
        nmod_poly_mat_t ker3;      /* actual init */

        /* first recursive call on left submatrix */
        nmod_poly_mat_window_init(submat, pmat, 0, 0, m, n - n/2);
        nullity = nmod_poly_mat_kernel_zls_approx(ker, pivind, shift, submat);
        nmod_poly_mat_window_clear(submat);

        /* if first kernel is empty, early exit */
        if (nullity == 0)
            return 0;

        /* actual kernel: `nullity` first rows */
        nmod_poly_mat_window_init(ker1, ker, 0, 0, nullity, m);

        /* residual: (first ker) * (right submatrix) */
        /* overwrite left columns of pmat */
        nmod_poly_mat_window_init(submat, pmat, 0, n - n/2, m, n);
        nmod_poly_mat_window_init(residual, pmat, 0, 0, nullity, n/2);
        _kernel_mul_split(residual, ker1, submat, -1, 0);
        nmod_poly_mat_window_clear(submat);

        /* recursive call 2, on residual */
        nmod_poly_mat_init(ker2, nullity, nullity, pmat->modulus);
        slong * pivind2 = FLINT_ARRAY_ALLOC(nullity, slong);
        nullity = nmod_poly_mat_kernel_zls_approx(ker2, pivind2, shift, residual);
        nmod_poly_mat_window_clear(residual);

        /* multiply bases and update pivind */
        nmod_poly_mat_init(ker3, nullity, m, pmat->modulus);
        nmod_poly_mat_window_init(submat, ker2, 0, 0, nullity, ker2->c);
        _kernel_mul_split(ker3, submat, ker1, -1, 0);
        nmod_poly_mat_window_clear(submat);
        nmod_poly_mat_window_clear(ker1);

        for (slong i = 0; i < nullity; i++)
            for (slong j = 0; j < m; j++)
                FLINT_SWAP(nmod_poly_struct,
                           *nmod_poly_mat_entry(ker, i, j),
                           *nmod_poly_mat_entry(ker3, i, j));

        for (slong i = 0; i < nullity; i++)
            pivind[i] = pivind[pivind2[i]];

        nmod_poly_mat_clear(ker2);
        nmod_poly_mat_clear(ker3);
        flint_free(pivind2);

        return nullity;
    }

    /* here, we are in the case 1 <= n <= m/2 */

    /* `nnz` <-- the number of nonzero rows of `pmat` */
    /* `buf` <-- `rdeg(pmat)`, -1 on zero rows        */
    slong * buf = FLINT_ARRAY_ALLOC(m, slong);
    ulong nnz = nmod_poly_mat_row_degree(buf, pmat, NULL);

    /* early exit: pmat == 0 => kernel is identity */
    if (nnz == 0)
    {
        nmod_poly_mat_one(ker);
        for (slong i = 0; i < m; i++)
            pivind[i] = i;
        flint_free(buf);
        return m;
    }

    /* zero rows: kernel of the nonzero rows, completed by unit vectors */
    /* (see [note:zero_rows])                                           */
    if (nnz < (ulong) m)
    {
        const slong nullity = _kernel_zls_zero_rows(ker, pivind, shift, pmat, buf);
        flint_free(buf);
        return nullity;
    }

    /* from here on, pmat has no zero row: buf = rdeg(pmat) >= 0 */

    /* diff_shift = min(shift - buf) (see [note:diff_shift]) */
    slong diff_shift = shift[0] - buf[0];
    for (slong i = 1; i < m; i++)
        if (diff_shift > shift[i] - buf[i])
            diff_shift = shift[i] - buf[i];

    /* approximation order (see [note:choice_of_order]), with                */
    /* t = shift - diff_shift >= rdeg(pmat):                                 */
    /*     order = 1 + dbound,                                               */
    /*     dbound = 1 + floor((sum(t) - 1) / (m - n))                        */
    /*            = ceil(sum(t) / (m - n))                                   */
    /* generically, every row of a t-ordered weak Popov kernel basis has     */
    /* t-degree <= dbound, so the approximant basis contains the kernel. The */
    /* order depends on sum(t), not on deg(pmat): a few rows of large degree */
    /* do not make it large. It is about half the order of                   */
    /* [LNVZ22, Algo. 1].                                                    */
    /* FIXME could be improved in case of unbalanced cdeg(pmat), e.g. based  */
    /* on [LNVZ22, Theorem 3.3, Proof of Item (iii)], see also               */
    /* [note:degree_bounds]                                                  */
    slong order = - m * diff_shift - 1;  /* order <- sum(t) - 1 */
    for (slong i = 0; i < m; i++)
        order += shift[i];
    /* C division rounds toward zero: the case sum(t) == 0, for which      */
    /* dbound = 0 (then pmat is constant), is treated apart                */
    order = (order < 0) ? 1 : 2 + order / (m - n);

    nmod_poly_mat_pmbasis(ker, shift, pmat, order);

    /* run an easy degree-based detection of rows in the kernel                   */
    /* (see [note:choice_of_order] and [note:diff_shift] for explanations)        */
    /* permute ker into [kernel rows \\ other rows], preserving increasing pivots */
    slong nullity = 0;
    slong sum_pmatdeg = 0;  /* sum of degrees of rows of pmat in complement of pivind */
    int no_constant_pivot_and_row = 1;  /* for early exit below, condition [3.] */
    for (slong i = 0; i < m; i++)
    {
        if (shift[i] < order + diff_shift)
        {
            /* this is a row in the kernel, see [note:row_in_kernel] */
            pivind[nullity] = i;
            nullity += 1;
        }
        else
        {
            if (buf[i] == 0 && nmod_poly_degree(nmod_poly_mat_entry(ker, i, i)) <= 0)
                no_constant_pivot_and_row = 0;
            sum_pmatdeg += buf[i];
            buf[i - nullity] = i;
        }
    }
    for (slong i = nullity; i < m; i++)
        pivind[i] = buf[i - nullity];

    nmod_poly_mat_permute_rows(ker, pivind, shift);

    /* early exit, subject to 3 conditions on the already identified `nullity`-many kernel rows  */
    /* [1.] nullity >= m - n  (otherwise, we are missing >= m - n - nullity > 0 kernel rows)     */
    /* [2.] sum_pmatdeg is equal to the sum of pivot degrees of already identified kernel rows   */
    /* [3.] among the rows of `ker` that are possibly not in the kernel (those at indices from   */
    /*    `nullity` to `m-1`), none has a pivot index `i` which is such that both the pivot      */
    /*     degree of that row is zero and the row `i` of the input matrix has degree 0           */
    /* (see [note:early_exit] for more details)                                                  */
    if (nullity >= m - n && no_constant_pivot_and_row)  /* condition [1.] and [3.] */
    {
        slong sum_pivdeg = 0;
        for (slong i = 0; i < nullity; i++)
            sum_pivdeg += nmod_poly_degree(nmod_poly_mat_entry(ker, i, pivind[i]));

        if (sum_pivdeg == sum_pmatdeg)  /* condition [2.] */
        {
            flint_free(buf);
            return nullity;
        }
    }

    /* proceed with remaining rows */
    nmod_poly_mat_t residual;  /* actual init */
    nmod_poly_mat_t ker2;      /* actual init */
    nmod_poly_mat_t prod;      /* actual init */
    nmod_poly_mat_t approx;    /* window only */
    nmod_poly_mat_t ker2nz;    /* window only */

    nmod_poly_mat_init(residual, m - nullity, n, pmat->modulus);
    nmod_poly_mat_window_init(approx, ker, nullity, 0, m, m);
    /* residual = (approx * pmat) div x^order; for each row p of approx,     */
    /* rdeg(p * pmat) <= rdeg_t(p) = rdeg_s(p) - diff_shift, where rdeg_s(p) */
    /* is the corresponding entry of shift (see [note:diff_shift])           */
    slong nhi = order + 1;
    for (slong i = nullity; i < m; i++)
        nhi = FLINT_MAX(nhi, shift[i] - diff_shift + 1);
    _kernel_mul_split(residual, approx, pmat, order, nhi);
    /* a zero row of residual comes from a row of approx which is already in */
    /* the kernel (with t-degree >= order); the recursive call handles it as */
    /* a unit vector (see [note:zero_rows]). [LNVZ22, Algo. 1, line 24]      */
    /* moves such rows to the kernel directly, which would save part of the  */
    /* final product.                                                        */

    /* kernel of residual */
    nmod_poly_mat_init(ker2, m - nullity, m - nullity, pmat->modulus);
    const slong nullity2 = nmod_poly_mat_kernel_zls_approx(ker2, buf, shift + nullity, residual);

    /* if nullity2 == 0, whole kernel already computed, just return */
    if (nullity2 == 0)
    {
        nmod_poly_mat_clear(residual);
        nmod_poly_mat_window_clear(approx);
        nmod_poly_mat_clear(ker2);
        flint_free(buf);
        return nullity;
    }

    /* multiply to get missing part of kernel */
    nmod_poly_mat_window_init(ker2nz, ker2, 0, 0, nullity2, ker2->c);
    nmod_poly_mat_init(prod, nullity2, m, pmat->modulus);
    _kernel_mul_split(prod, ker2nz, approx, -1, 0);
    for (slong i = 0; i < nullity2; i++)
        for (slong j = 0; j < m; j++)
            FLINT_SWAP(nmod_poly_struct,
                       *nmod_poly_mat_entry(prod, i, j),
                       *nmod_poly_mat_entry(approx, i, j));

    /* update pivind */
    for (slong i = 0; i < nullity2; i++)
        pivind[nullity + i] = pivind[nullity + buf[i]];

    /* permute rows to get increasing pivot indices */
    nmod_poly_mat_window_clear(approx);
    nmod_poly_mat_window_init(approx, ker, 0, 0, nullity + nullity2, m);
    slong i1 = 0;
    slong i2 = 0;
    for (slong i = 0; i < nullity + nullity2; i++)
    {
        if (i2 == nullity2 || (i1 < nullity && pivind[i1] < pivind[nullity + i2]))
        {
            buf[i] = i1;
            i1++;
        }
        else
        {
            buf[i] = nullity + i2;
            i2++;
        }
    }

    nmod_poly_mat_permute_rows(approx, buf, NULL);
    slong * tmp = FLINT_ARRAY_ALLOC(nullity + nullity2, slong);
    for (slong i = 0; i < nullity + nullity2; i++)
        tmp[i] = pivind[i];
    for (slong i = 0; i < nullity + nullity2; i++)
        pivind[i] = tmp[buf[i]];
    for (slong i = 0; i < nullity + nullity2; i++)
        tmp[i] = shift[i];
    for (slong i = 0; i < nullity + nullity2; i++)
        shift[i] = tmp[buf[i]];
    flint_free(tmp);

    nmod_poly_mat_window_clear(ker2nz);
    nmod_poly_mat_window_clear(approx);
    nmod_poly_mat_clear(residual);
    nmod_poly_mat_clear(ker2);
    nmod_poly_mat_clear(prod);
    flint_free(buf);

    return nullity + nullity2;
}


/*------------------------------------------------------------*/
/*                         NOTES                             */
/*------------------------------------------------------------*/

/** Degree bounds for kernel bases.
 *
 * [note:degree_bounds]
 * Assume pmat is m x n and has rank r. Let P be some shifted ordered weak
 * Popov kernel basis P of `pmat` (for a given, arbitrary shift).
 *
 * From Theorem 3.3 in
 *    Labahn, Neiger, Vu, Zhou
 *    Proceedings ISSAC 2022 (https://arxiv.org/pdf/2202.09329)
 * we get non-strict upper bounds on the sum of pivot degree of P:
 *   - it is <= sum of degrees of the r largest-degree rows of pmat
 *   - in particular, <= r * deg(pmat) <= min(m, n) * deg(pmat)
 *   - it is <= sum of degrees of the r largest-degree columns of pmat
 *   - in particular, <= sum of cdeg(pmat)
 * The third item above is not found directly in the mentioned Theorem 3.3,
 * but is deduced from its second item (part with deg(det(S)))
 *
 * The maximum of pivot degrees can be equal to their sum (although this is not
 * expected generically). This gives bounds on the non-pivot part of P: pick
 * one of the bounds above, and add the amplitude `max(shift) - min(shift)`.
 * Note that this might be pessimistic for shifts of large amplitude, since we
 * also have bounds coming from formula with adjugate of some submatrix of `pmat`.
 * For example, one can prove that the non-pivot part of the kernel basis in
 * shifted Popov form has degree at most n * deg(pmat)  (cf Lemma 12 of a research
 * internship report by Vu Thi Xuan).
 *
 * The maximum pivot degree also gives a bound on max(rdeg_s(ker)): pick one of
 * the bounds above for sum of pivot degrees, and add max(s) - min(s).
 */

/** Expected degree / order for approximation in ZLS algorithm.
 *
 * [note:choice_of_order]
 * In the branch n <= m/2 of zls_approx, pmat has no zero row (see
 * [note:zero_rows]) and the approximant basis is computed at order
 *     order = 1 + dbound,
 *     dbound = 1 + floor((sum(t) - 1) / (m - n)) = ceil(sum(t) / (m - n)),
 * where t = shift - diff_shift satisfies t >= rdeg(pmat) entry-wise, with
 * equality for at least one entry (see [note:diff_shift]). The two
 * expressions of dbound are equal for any integer sum(t) >= 0; the code uses
 * the first one.
 *
 * The choice of order only matters for efficiency. Any order > 0 gives a
 * correct result: the rows of the approximant basis of t-degree < order are
 * kernel rows (see [note:row_in_kernel]), and the missing kernel rows are
 * obtained from the recursive call on the residual. The aim is that, in the
 * generic case, the approximant basis already contains the whole kernel, so
 * that the recursive call is skipped (see [note:early_exit]).
 *
 * -> Why dbound. Let P be a t-ordered weak Popov kernel basis of pmat. It has
 * nullity >= m - n rows (the rank of pmat is at most n), and
 * sum(rdeg_t(P)) <= sum(t) since t >= rdeg(pmat) (Zhou-Labahn-Storjohann,
 * 2012; [LNVZ22, Thm. 3.3 (iii)]). Generically, the entries of rdeg_t(P) are
 * balanced: each is D = max(rdeg_t(P)) or D - 1. Then
 * D + (nullity - 1) (D - 1) <= sum(t), hence
 *     D <= ceil(sum(t) / nullity) <= ceil(sum(t) / (m - n)) = dbound < order,
 * so that the rows of the approximant basis of t-degree < order form a
 * kernel basis.
 *
 * -> Why not deg(pmat). The order can be <= deg(pmat). A kernel row p with a
 * nonzero entry at index i has rdeg_t(p) >= t[i] >= rdeg(pmat[i,:]): if
 * t[i] >= order, p cannot be detected at this step, and it is obtained from
 * the recursive call on the residual, as in the Zhou-Labahn-Storjohann
 * algorithm. For balanced row degrees this does not happen: if all rows of
 * pmat have degree d, then sum(t) >= m d, so dbound >= d and
 * order > deg(pmat). It does happen when a few rows have much larger degree,
 * e.g. [A ; b] in linear system solving with deg(b) about dim(A) deg(A):
 * then, requiring order > deg(pmat) would mean an approximant basis at order
 * about dim(A) deg(A), whereas dbound is a small multiple of deg(A). The
 * products of the residual stage then have operands with a few long rows,
 * which are handled by _kernel_mul_split.
 *
 * -> Difference with [LNVZ22, Algo. 1]. The paper uses the order
 * ceil(2 sum(t) / (m - n)), about twice ours: with it, at most (m + n) / 2
 * rows of the approximant basis are not kernel rows [LNVZ22, Lemma 4.3,
 * Cor. 4.4], which bounds the size of the residual in the cost analysis.
 * Here the target is the generic case above. The order 1 + dbound still
 * ensures that at least one kernel row is detected: if k rows of the
 * approximant basis have t-degree >= order, then, since all t-degrees are
 * nonnegative, k order <= sum(t) + n order (the sum of the t-degrees of
 * the approximant basis), so k <= sum(t) / order + n < (m - n) + n = m.
 *
 * [note:row_in_kernel]
 * Since t >= rdeg(pmat), we have rdeg(p*pmat) <= rdeg_t(p) for any row vector
 * p. So, an approximant p of pmat at order `order` with rdeg_t(p) < order
 * (in particular, with rdeg_t(p) <= dbound) is an actual kernel row:
 * p*pmat has degree < order and is zero modulo x^order, hence p*pmat == 0.
 *
 * [note:zero_rows]
 * Let Z and N be the indices of the zero rows and of the nonzero rows of
 * pmat. A row vector p is in the kernel of pmat if and only if p[N] is in the
 * kernel of pmat[N,:], whatever p[Z]. Hence a kernel basis of pmat is formed
 * by the unit vectors e_i for i in Z (s-degree s[i], s-pivot i), together
 * with a kernel basis of pmat[N,:] extended by zeros at the indices in Z
 * (which keeps its s-degrees and its s-pivots, renumbered). This is what
 * _kernel_zls_zero_rows computes.
 *
 * Removing the zero rows only matters for efficiency. The order of
 * [note:choice_of_order] assumes that the t-degrees of the m - n kernel rows
 * are balanced, which the e_i break:
 * - if t is small on Z (typically t[i] = 0 when shift = rdeg(pmat)), the e_i
 *   use almost none of sum(t), and the other kernel rows, about |N| - n of
 *   them, have t-degree about sum(t) / (|N| - n), which can exceed dbound =
 *   ceil(sum(t) / (m - n)): they are only found through the residual;
 * - if t is large on Z, t[Z] increases sum(t) and thus the order, although
 *   it only concerns the e_i; if t[i] >= order, e_i itself is not detected,
 *   and is recovered through the residual (typically from a zero row of it).
 * When there is no zero row, this costs nothing, since the number of nonzero
 * rows is computed anyway.
 *
 * Zero rows are only removed in the branch n <= m/2 (the column split passes
 * them to its recursive calls). If |N| < 2n, the kernel of pmat[N,:] is
 * computed by the column split, which at small degrees can be slower than a
 * single approximant basis (see the TODO on the threshold of the split).
 */

/** Minimizing the shift, and early exit.
 *
 * Here are some arguments, trying to make things clearer with respect to two
 * things: the definition of `diff_shift`, and the design of some correct early
 * exit.
 *
 * [note:diff_shift]
 * Let's call `s` the shift `shift` before the call to `pmbasis` in zls_approx
 * (at this point, `pmat` has no zero row, see [note:zero_rows]), and let `t` be
 * the same shift but to which we subtract `diff_shift` to all entries. The
 * choice of `diff_shift` is such that the entries of `t` are at least the row
 * degree of `pmat`, indeed,
 * ```t[i] = s[i] - diff_shift >= s[i] - (s[i] - rdeg(pmat[i,:])) = rdeg(pmat[i,:])```
 * and in fact, `diff_shift` is the largest value such that this is true, or in
 * other words, it minimizes `t` while ensuring `t >= rdeg(pmat)` (i.e., any
 * larger `diff_shift` would lead to some entry of `t` being strictly less than
 * the corresponding entry of `rdeg(pmat)`).
 *
 * This `pmbasis` call computes `ker` which is an `s`-ordered weak Popov
 * approximant basis of `pmat` at order `order`. It also in-place modifies `shift`
 * into the `s`-row degree of `ker`. Since `s` is `t` translated by `+diff_shift`,
 * then `ker` is also in `t`-ordered weak Popov form, and the modified `shift`
 * contains the `t`-row degree of `ker` translated by `+diff_shift`, that is,
 * after the `pmbasis` call,
 * ```shift[i] = rdeg_s(ker[i,:]) = rdeg_t(ker[i,:]) + diff_shift```
 * from which we can deduce that
 * ```shift[i] >= rdeg(ker[i,:] * pmat) + diff_shift```
 *
 * Then, what is said to be an "easy degree-based detection of rows in the kernel"
 * in the code is based on the fact that any `i` such that ```shift[i] < order +
 * diff_shift``` corresponds to a row in the kernel of `pmat`. Indeed this
 * inequality is equivalent to the fact that the `t`-row degree of the row
 * `ker[i,:]` is strictly less than `order`, and since `t >= rdeg(pmat)`
 * entry-wise, this means that the row degree of `ker[i,:] * pmat` is strictly
 * less than `order`. But at the same time this vector is zero modulo `x**order`,
 * so it must be actually zero, meaning that `ker[i,:]` is in the kernel of
 * `pmat`.
 *
 * This is in good part a mix of classical arguments concerning kernel bases via
 * approximant bases, but hopefully this (technically) highlights the interplay
 * between `diff_shift` and the row degree of `pmat`. Another classical thing is
 * that every such row is actually minimal, in the sense that there is no vector
 * in the kernel of `pmat` which would have `s`-pivot at index `i` and whose
 * `s`-pivot degree would be strictly less than the `s`-pivot degree
 * `deg(ker[i,i])` of the `i`-th row of the approximant basis `ker`.
 * 
 * [note:early_exit]
 * Concerning the early exit, the quantity `sum_pmatdeg` is the sum of degrees of
 * rows of `pmat` that are not indexed by some `i` that are already known to be an
 * `s`-pivot index of some kernel row (based on the above detection of kernel
 * rows). So this subset of rows includes all rows of `pmat` that are not indexed
 * by the `s`-pivot index of some row in the kernel, and this means that
 * `sum_pmatdeg` is a (non-strict) upper bound on the sum of `s`-pivot degrees of
 * an `s`-weak Popov kernel basis. In particular, if the sum of `s`-pivot degrees
 * of already found kernel rows is equal to `sum_pmatdeg`, then the only possibly
 * missing kernel rows must have `s`-pivot degree equal to zero.
 * 
 * Expanding upon this, we can say a bit more: "the corresponding row of `pmat`
 * must be constant as well".
 * 
 * Through `pmbasis`, we have found indices $i_0 < \cdots < i_{k-1}$ of rows that
 * are in the kernel of `pmat` (where $k$ is the partial nullity computed in the
 * code: number of already found kernel rows). Because they are extracted from an
 * `s`-ordered weak Popov basis, their `s`-pivots are on the diagonal, so at the
 * column positions $i_0 < \cdots < i_{k-1}$. Let $\delta_0,\ldots,\delta_{k-1}$
 * be their `s`-pivot degrees. Let $j_0 < \cdots < j_{m-k-1}$ be the $m-k$ indices
 * that form the complement, so that the union of the $i$'s and the $j$'s gives
 * all row indices $\{0,\ldots,m-1\}$ of the input (or column indices of the
 * kernel). Finally, let $r_0,\ldots,r_{m-k-1}$ be the degrees of the $m-k$ rows
 * indexed by these $j$'s in the input matrix `pmat` (which has no zero row at
 * this point, see [note:zero_rows]).
 * 
 * What we know is that $\delta_0 + \cdots + \delta_{k-1} \le r_0 + \cdots +
 * r_{m-k-1}$. This is possibly a strict inequality, even when we have found the
 * whole kernel. Now suppose this is an equality. If we are missing any row from
 * the kernel, say a row with `s`-pivot index $\pi$, adding this row to the
 * collection of "already known kernel rows" makes both the left-hand side
 * increase (or stay the same) and the right-hand side decrease (or stay the
 * same). But the inequality must remain valid, so, since we already had
 * equality, both sides must stay the same. What this means is that this kernel
 * row must have `s`-pivot degree `0`, and also, the corresponding row $\pi$ of
 * `pmat` must have degree `0`.
 *
 * This seems sufficient to make the existing early exit strategy correct without
 * impacting its efficiency: to make sure we can early exit, we just have to add a
 * check that, among the rows not detected to be in the kernel, none has pivot
 * degree $0$ with also the corresponding row of `pmat` having degree $0$.
 * 
 */
