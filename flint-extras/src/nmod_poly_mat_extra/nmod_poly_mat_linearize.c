/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/nmod_vec.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_poly_mat.h>

#include "nmod_poly_mat_utils.h"
#include "nmod_mat_poly_extra/impl.h"   /* _pml_transpose */

/*------------------------------------------------------------*/
/* by blocks of columns: one contiguous copy per entry        */
/*------------------------------------------------------------*/

void
nmod_poly_mat_linearize_columns(nmod_mat_t cmat,
                                const nmod_poly_mat_t pmat,
                                const slong * orders)
{
    for (slong i = 0; i < pmat->r; i++)
    {
        slong off = 0;
        for (slong b = 0; b < pmat->c; b++)
        {
            const slong d = orders[b];
            if (d == 0)
                continue;

            const nmod_poly_struct * f = nmod_poly_mat_entry(pmat, i, b);
            const slong len = FLINT_MIN(f->length, d);
            nn_ptr row = nmod_mat_entry_ptr(cmat, i, off);

            _nmod_vec_set(row, f->coeffs, len);
            _nmod_vec_zero(row + len, d - len);
            off += d;
        }
    }
}

void
nmod_poly_mat_polynomialize_columns(nmod_poly_mat_t pmat,
                                    const nmod_mat_t cmat,
                                    const slong * orders)
{
    for (slong i = 0; i < pmat->r; i++)
    {
        slong off = 0;
        for (slong b = 0; b < pmat->c; b++)
        {
            const slong d = orders[b];
            nmod_poly_struct * f = nmod_poly_mat_entry(pmat, i, b);

            if (d == 0)
            {
                nmod_poly_zero(f);
                continue;
            }

            nmod_poly_fit_length(f, d);
            _nmod_vec_set(f->coeffs, nmod_mat_entry_ptr(cmat, i, off), d);
            _nmod_poly_set_length(f, d);
            _nmod_poly_normalise(f);
            off += d;
        }
    }
}

/*------------------------------------------------------------*/
/* by blocks of rows: for each row a of pmat, one             */
/* transposition between the d_a x n block of cmat and the    */
/* coefficient arrays of the n entries of that row            */
/*------------------------------------------------------------*/

/* Kernels and schedules as in the conversions between nmod_poly_mat and
 * nmod_mat_poly, which do the same transpositions: the widest kernel that
 * fits, and the default schedule PML_CONV_DST_MAJOR towards the constant
 * matrix, the destination-major one towards the polynomials (see
 * PML_SETFROM_MAT_POLY_DST_MAJOR in nmod_poly_mat_utils.c). */

void
nmod_poly_mat_linearize_rows(nmod_mat_t cmat,
                             const nmod_poly_mat_t pmat,
                             const slong * orders)
{
    const slong m = pmat->r, n = pmat->c;
    slong dmax = 0;

    for (slong a = 0; a < m; a++)
        dmax = FLINT_MAX(dmax, orders[a]);

    if (n == 0 || dmax == 0)
        return;

    nn_ptr * dst = (nn_ptr *) flint_malloc(dmax * sizeof(nn_ptr));
    nn_srcptr * src = (nn_srcptr *) flint_malloc(n * sizeof(nn_srcptr));
    slong * slen = (slong *) flint_malloc(n * sizeof(slong));

    slong off = 0;
    for (slong a = 0; a < m; a++)
    {
        const slong d = orders[a];
        if (d == 0)
            continue;

        for (slong k = 0; k < d; k++)
            dst[k] = nmod_mat_entry_ptr(cmat, off + k, 0);

        for (slong j = 0; j < n; j++)
        {
            const nmod_poly_struct * f = nmod_poly_mat_entry(pmat, a, j);
            src[j] = f->coeffs;
            slen[j] = FLINT_MIN(f->length, d);
        }

        _pml_transpose(dst, d, src, slen, n,
                       _pml_transpose_narrow_kernel(_pml_transpose_widest_kernel(), d, n),
                       PML_CONV_DST_MAJOR);
        off += d;
    }

    flint_free(dst);
    flint_free(src);
    flint_free(slen);
}

void
nmod_poly_mat_polynomialize_rows(nmod_poly_mat_t pmat,
                                 const nmod_mat_t cmat,
                                 const slong * orders)
{
    const slong m = pmat->r, n = pmat->c;
    slong dmax = 0;

    for (slong a = 0; a < m; a++)
        dmax = FLINT_MAX(dmax, orders[a]);

    if (n == 0)
        return;

    nn_ptr * dst = (nn_ptr *) flint_malloc(n * sizeof(nn_ptr));
    nn_srcptr * src = (nn_srcptr *) flint_malloc(FLINT_MAX(dmax, 1) * sizeof(nn_srcptr));
    slong * slen = (slong *) flint_malloc(FLINT_MAX(dmax, 1) * sizeof(slong));

    /* the source rows, rows of cmat, are full */
    for (slong k = 0; k < dmax; k++)
        slen[k] = n;

    slong off = 0;
    for (slong a = 0; a < m; a++)
    {
        const slong d = orders[a];

        if (d == 0)
        {
            for (slong j = 0; j < n; j++)
                nmod_poly_zero(nmod_poly_mat_entry(pmat, a, j));
            continue;
        }

        for (slong j = 0; j < n; j++)
        {
            nmod_poly_struct * f = nmod_poly_mat_entry(pmat, a, j);
            nmod_poly_fit_length(f, d);
            dst[j] = f->coeffs;
        }

        for (slong k = 0; k < d; k++)
            src[k] = nmod_mat_entry_ptr(cmat, off + k, 0);

        _pml_transpose(dst, n, src, slen, d,
                       _pml_transpose_narrow_kernel(_pml_transpose_widest_kernel(), n, d),
                       1);

        for (slong j = 0; j < n; j++)
        {
            nmod_poly_struct * f = nmod_poly_mat_entry(pmat, a, j);
            _nmod_poly_set_length(f, d);
            _nmod_poly_normalise(f);
        }
        off += d;
    }

    flint_free(dst);
    flint_free(src);
    flint_free(slen);
}
