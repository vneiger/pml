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
#include <flint/nmod_mat.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>
#include <flint/perm.h>

#include "nmod_poly_mat_extra.h"
#include "nmod_poly_mat_forms.h"
#include "nmod_poly_mat_kernel.h"
#include "nmod_poly_mat_multiply.h"
#include "nmod_poly_mat_utils.h"

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
        /* flint_printf("calling zls approx, left kernel, on:\n"); */
        /* nmod_poly_mat_print_pretty(mat_it, "x"); */
        /* flint_printf("with shift %{slong*}\n", _shift, pmat->c); */
        
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
        nmod_poly_mat_mul_classical(residual, ker1, submat);
        nmod_poly_mat_window_clear(submat);

        /* recursive call 2, on residual */
        nmod_poly_mat_init(ker2, nullity, nullity, pmat->modulus);
        slong * pivind2 = FLINT_ARRAY_ALLOC(nullity, slong);
        nullity = nmod_poly_mat_kernel_zls_approx(ker2, pivind2, shift, residual);
        nmod_poly_mat_window_clear(residual);

        /* multiply bases and update pivind */
        nmod_poly_mat_init(ker3, nullity, m, pmat->modulus);
        nmod_poly_mat_window_init(submat, ker2, 0, 0, nullity, ker2->c);
        nmod_poly_mat_mul(ker3, submat, ker1);
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
    /* `buf` <-- `max(0, rdeg(pmat))` */
    slong * buf = FLINT_ARRAY_ALLOC(m, slong);
    ulong nnz = nmod_poly_mat_row_degree_zero(buf, pmat, NULL);

    /* early exit: pmat == 0 => kernel is identity */
    if (nnz == 0)
    {
        nmod_poly_mat_one(ker);
        for (slong i = 0; i < m; i++)
            pivind[i] = i;
        flint_free(buf);
        return m;
    }

    /* find max(rdeg(pmat)), and diff_shift = min(shift - buf) */
    slong maxdeg = buf[0];
    slong diff_shift = shift[0] - maxdeg;
    for (slong i = 1; i < m; i++)
    {
        slong d = buf[i];
        slong diff = shift[i] - d;
        if (d > maxdeg)
            maxdeg = d;
        if (diff_shift > diff)
            diff_shift = diff;
    }

    /* build order for approximation (this choice diverges from [LNVZ22, Algo.1]): */
    /*     order = 1 + max(maxdeg, 1 + floor((sum(shift) - 1) / (m - n)))          */
    /* -> approx basis "generically" captures the whole kernel                     */
    /* (see [note:choice_of_order] for explanations)                               */
    /* FIXME could be improved in case of unbalanced cdeg(pmat), based on          */
    /* [LNVZ22, Theorem 3.3, Proof of Item (iii)], see also [note:degree_bounds]   */
    slong order = - m * diff_shift - 1;
    for (slong i = 0; i < m; i++)
        order += shift[i];
    order = 1 + order / (m - n);
    order = 1 + FLINT_MAX(maxdeg, order);

    /* flint_printf("order : %ld\n", order); */
    nmod_poly_mat_pmbasis(ker, shift, pmat, order);

    /* run an easy degree-based detection of rows in the kernel */
    /* permute ker into [kernel rows \\ other rows], preserving increasing pivots */
    slong nullity = 0;
    slong sum_pmatdeg = 0;  /* sum of degrees of rows of pmat in complement of pivind */
    for (slong i = 0; i < m; i++)
    {
        /* flint_printf("i : %ld, order+diff_shift : %ld, shift[i] : %ld --> ", i, order+diff_shift, shift[i]); */
        if (shift[i] < order + diff_shift)
        {
            /* FIXME remove */ /* flint_printf("yes!\n"); */
            /* this is a row in the kernel, see [note:row_in_kernel] */
            pivind[nullity] = i;
            nullity += 1;
        }
        else
        {
            /* FIXME remove */ /* flint_printf("no..\n"); */
            sum_pmatdeg += buf[i];
            buf[i - nullity] = i;
        }
    }
    for (slong i = nullity; i < m; i++)
        pivind[i] = buf[i - nullity];

    nmod_poly_mat_permute_rows(ker, pivind, shift);

    /* early exit */
    if (nullity >= m - n)  /* otherwise, there must be some kernel rows not in ker yet */
    {
        /* flint_printf("trying early exit...\n"); */
        slong sum_pivdeg = 0;
        for (slong i = 0; i < nullity; i++)
            sum_pivdeg += nmod_poly_degree(nmod_poly_mat_entry(ker, i, pivind[i]));

        if (sum_pivdeg == sum_pmatdeg)  /* whole kernel already found, see [note:degree_bounds] */
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
    nmod_poly_mat_middle_product_naive(residual, approx, pmat, order, order + maxdeg + 1);
    /* FIXME provide (and use) general middle_product interface */
    /* note: on non-generic input, one might want to check for zero rows in residual */

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
    nmod_poly_mat_mul(prod, ker2nz, approx);
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
 * For the moment, the code in zls_approx, when n <= m/2, uses the approximation order
 *     order = 1 + max(maxdeg, 1 + floor((sum(shift) - 1) / (m - n)))
 * where maxdeg is deg(pmat) and shift >= rdeg(pmat) entry-wise
 * (more precisely, the code implicitly ensures min(shift - rdeg(pmat)) == 0)
 *
 * Note that any order > 0 is enough to guarantee that the algorithm terminates,
 * but order too small may impact performance. We pick order that is large
 * enough to provide the whole kernel in some generic/frequent situations.
 *
 * -> `1 + maxdeg` is because any lower order is too low: it means ignoring top
 *  degree coefficients of `pmat`
 *
 * -> the other term, `1 + dbound` where
 *      dbound = 1 + floor((sum(shift) - 1) / (m - n))
 * is because this is the expected maximum entry of the shift-row degree of a
 * shift-weak Popov kernel basis P, if degrees are balanced in a generic way.
 *
 * Indeed: let dbound = max(rdeg_s(P))
 * . sum(rdeg_s(P)) <= sum(s) is known (Zhou-Labahn-Storjohann, 2012)
 * . if we assume degrees are well balanced (which happens generically):
 * then rdeg_s(P) consists of `nullity` values all equal to dbound or dbound - 1
 * -> dbound + (nullity - 1) * (dbound - 1) <= sum(s)
 * -> dbound <= floor((sum(s) + (nullity - 1)) / nullity)
 *            = 1 + floor((sum(s) - 1) / nullity)                                        
 *           <= 1 + floor((sum(s) - 1) / (m - n)).
 *
 * [note:row_in_kernel]
 * Note also that since s >= rdeg(pmat), we have rdeg(p*pmat) <= rdeg_s(p) < order
 * for any row vector p such that rdeg_s(p) <= dbound. So, approximants of pmat
 * whose s-degree is <= dbound are necessarily actual kernel rows p*pmat == 0.
 */
