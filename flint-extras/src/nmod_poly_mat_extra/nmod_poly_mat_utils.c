/*
    Copyright (C) 2025 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h> // for qsort
#include <flint/nmod_poly.h>
#include "nmod_poly_mat_utils.h"
#include "nmod_mat_poly.h"
#include "nmod_mat_poly_extra/impl.h"   /* _pml_transpose */
#include "nmod_poly_mat_extra/impl.h"

/**********************************************************************
*                    ROW ROTATION DOWNWARD/UPWARD                    *
**********************************************************************/

void _nmod_poly_mat_rotate_rows_downward(nmod_poly_mat_t mat, slong * vec, slong i, slong j)
{
    if (i != j)
    {
        if (vec != NULL)
        {
            slong tmp_vec = vec[j];
            for (slong ii = j; ii > i; ii--)
                vec[ii] = vec[ii-1];
            vec[i] = tmp_vec;
        }

#if __FLINT_VERSION < 3 || (__FLINT_VERSION == 3 && __FLINT_VERSION_MINOR < 3)
        nmod_poly_struct * tmp_mat = mat->rows[j];
        for (slong ii = j; ii > i; ii--)
            mat->rows[ii] = mat->rows[ii-1];
        mat->rows[i] = tmp_mat;
#else
        nmod_poly_struct * tmp_row = flint_malloc(mat->c * sizeof(nmod_poly_struct));
        for (slong jj = 0; jj < mat->c; jj++)
            tmp_row[jj] = *nmod_poly_mat_entry(mat, j, jj);
        for (slong ii = j; ii > i; ii--)
            for (slong jj = 0; jj < mat->c; jj++)
                *nmod_poly_mat_entry(mat, ii, jj) = *nmod_poly_mat_entry(mat, ii-1, jj);
        for (slong jj = 0; jj < mat->c; jj++)
            *nmod_poly_mat_entry(mat, i, jj) = tmp_row[jj];
        flint_free(tmp_row);
#endif
    }
}

void _nmod_poly_mat_rotate_rows_upward(nmod_poly_mat_t mat, slong * vec, slong i, slong j)
{
    if (i != j)
    {
        if (vec != NULL)
        {
            slong tmp_vec = vec[i];
            for (slong ii = i; ii < j; ii++)
                vec[ii] = vec[ii+1];
            vec[j] = tmp_vec;
        }

#if __FLINT_VERSION < 3 || (__FLINT_VERSION == 3 && __FLINT_VERSION_MINOR < 3)
        nmod_poly_struct * tmp_mat = mat->rows[i];
        for (slong ii = i; ii < j; ii++)
            mat->rows[ii] = mat->rows[ii+1];
        mat->rows[j] = tmp_mat;
#else
        nmod_poly_struct * tmp_row = flint_malloc(mat->c * sizeof(nmod_poly_struct));
        for (slong jj = 0; jj < mat->c; jj++)
            tmp_row[jj] = *nmod_poly_mat_entry(mat, i, jj);
        for (slong ii = i; ii < j; ii++)
            for (slong jj = 0; jj < mat->c; jj++)
                *nmod_poly_mat_entry(mat, ii, jj) = *nmod_poly_mat_entry(mat, ii+1, jj);
        for (slong jj = 0; jj < mat->c; jj++)
            *nmod_poly_mat_entry(mat, j, jj) = tmp_row[jj];
        flint_free(tmp_row);
#endif
    }
}

/**********************************************************************
*                    PERMUTE ROWS BY SORTING VEC                     *
**********************************************************************/

/* type for stable sort while retaining the permutation */
typedef struct
{
    slong value;
    slong index;
} slong_pair;

/* comparator for quicksort, lexicographic (total) order to ensure stable sort */
static inline int _slong_pair_compare(const void * a, const void * b)
{
    slong_pair aa = * (const slong_pair *) a;
    slong_pair bb = * (const slong_pair *) b;
    if (aa.value == bb.value)
    {
        if (aa.index < bb.index)
            return -1;
        else if (aa.index > bb.index)
            return 1;
        else // aa.index == bb.index
            return 0;
    }
    else if (aa.value < bb.value)
        return -1;
    else // aa.value > bb.value
        return 1;
}

/** Creates a permutation from the sorting of a list of integers
 * After running this, perm is the unique list of integers which sorts
 * the pairs (vec,index) increasingly, i.e.
 * vec[perm[0]] <= vec[perm[1]] < ... < vec[perm[n-1]]
 * All inputs must be already initialized/allocated. sorted_vec can alias vec.
 * If sorted_vec is NULL, it is simply ignored, the permuted vec is not
 * returned.
 *
 * \param perm permutation (list of integers), length n
 * \param sorted_vec list of integer sorted nondecreasingly, length n
 * \param vec list of integer to be sorted nondecreasingly, length n
 * \param n length
 * \param pair_tmp temporary storage, length n
 *
 */
static inline void _vec_sort_permutation(slong * perm,
                                         slong * sorted_vec,
                                         const slong * vec,
                                         slong n,
                                         slong_pair * pair_tmp)
{
    for (slong i = 0; i < n; i++)
    {
        pair_tmp[i].value = vec[i];
        pair_tmp[i].index = i;
    }

    qsort(pair_tmp, n, sizeof(slong_pair), _slong_pair_compare);

    for (slong i = 0; i < n; i++)
        perm[i] = pair_tmp[i].index;
    if (sorted_vec)
        for (slong i = 0; i < n; i++)
            sorted_vec[i] = pair_tmp[i].value;
}

void _nmod_poly_mat_permute_rows_by_sorting_vec(nmod_poly_mat_t mat,
                                                slong r,
                                                slong * vec,
                                                slong * perm)
{
    slong_pair * tmp = flint_malloc(r * sizeof(slong_pair));
    _vec_sort_permutation(perm, vec, vec, r, tmp);
    for (slong i = r; i < mat->r; i++)
        perm[i] = i;
    flint_free(tmp);
    nmod_poly_mat_permute_rows(mat, perm, NULL);
}


void _nmod_poly_mat_permute_columns_by_sorting_vec(nmod_poly_mat_t mat,
                                                slong r,
                                                slong * vec,
                                                slong * perm)
{
    slong_pair * tmp = flint_malloc(r * sizeof(slong_pair));
    _vec_sort_permutation(perm, vec, vec, r, tmp);
    for (slong i = r; i < mat->c; i++)
        perm[i] = i;
    flint_free(tmp);
    nmod_poly_mat_permute_columns(mat, perm, NULL);
}


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* SET FROM MATRIX POLYNOMIAL                                 */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/* This is the transposition of the conversion in the other direction, so it
   runs on the same blocked kernels; see nmod_mat_poly_extra/impl.h.  The two
   sides are exchanged, and with them their alignments: here the source rows
   (the coefficients of `matp`) are the 64-byte aligned ones and the
   destination rows (the coefficient arrays of the entries of `pmat`) are
   whatever `flint_realloc` returned, which is why this direction does not
   pick the same default kernel as the other one. */

/* Below this many words the pointer tables and their allocation dominate;
   just run the naive loop. */
#define _PML_SETFROM_TINY 512

/* Default kernel.  Not the widest one, unlike the other direction: the
   destination rows here are the coefficient arrays of nmod_poly entries, 16
   byte aligned at best, so every 64-byte store straddles two cache lines
   whereas a 32-byte one does so only half the time.  Measured over matrices
   4x4 to 128x128 and lengths 32 to 8192, the 4-wide kernel and the widest
   available one are within 2% of each other while the data fits in cache,
   and the 4-wide one is 10% ahead once it does not. */
#ifndef PML_SETFROM_MAT_POLY_KERNEL
# define PML_SETFROM_MAT_POLY_KERNEL NMOD_MAT_POLY_CONV_VEC4
#endif

/* Default schedule: always sweep the destination rows, that is, write each
   output polynomial from its constant coefficient upwards.  This direction
   does not need the target test that PML_CONV_DST_MAJOR makes for the other
   one, and for the same reason: what one wants is the *aligned* side of the
   transposition to be the scattered one.  Here that side is the source (the
   coefficients of the nmod_mat_poly, 64-byte aligned), so scattering the
   loads and sweeping the stores is right everywhere -- and on Apple silicon
   it is what the 128-byte cache line asks for anyway.  Measured on Zen 4 the
   other schedule is 20% to 80% slower over the same grid. */
#ifndef PML_SETFROM_MAT_POLY_DST_MAJOR
# define PML_SETFROM_MAT_POLY_DST_MAJOR 1
#endif

void _nmod_poly_mat_set_trunc_from_mat_poly(nmod_poly_mat_t pmat,
                                            const nmod_mat_poly_t matp,
                                            slong order,
                                            int kern,
                                            int dmaj)
{
    if (order > matp->length)
        order = matp->length;

    const slong r = pmat->r;
    const slong c = pmat->c;

    // prepare memory
    for (slong i = 0; i < r; i++)
        for (slong j = 0; j < c; j++)
            nmod_poly_fit_length(nmod_poly_mat_entry(pmat, i, j), order);

    if (order == 0 || r == 0 || c == 0)
    {
        for (slong i = 0; i < r; i++)
            for (slong j = 0; j < c; j++)
                _nmod_poly_set_length(nmod_poly_mat_entry(pmat, i, j), 0);
        return;
    }

    if ((double) r * (double) c * (double) order < (double) _PML_SETFROM_TINY)
    {
        for (slong k = 0; k < order; k++)
            for (slong i = 0; i < r; i++)
                for (slong j = 0; j < c; j++)
                    nmod_poly_mat_entry(pmat, i, j)->coeffs[k] = nmod_mat_poly_entry(matp, k, i, j);
    }
    else
    {
        /* Effective number of rows of the destination: the whole matrix when
           the coefficients of `matp` are contiguous, one matrix row otherwise. */
        const slong ndst = (matp->stride == c) ? r * c : c;

        if (kern < 0 || kern > NMOD_MAT_POLY_CONV_VEC8)
            kern = PML_SETFROM_MAT_POLY_KERNEL;
        kern = _pml_transpose_narrow_kernel(kern, ndst, order);

        if (dmaj < 0 || dmaj > 1)
            dmaj = PML_SETFROM_MAT_POLY_DST_MAJOR;

        nn_ptr * dst = (nn_ptr *) flint_malloc(r * c * sizeof(nn_ptr));
        nn_srcptr * src = (nn_srcptr *) flint_malloc(order * sizeof(nn_srcptr));
        slong * slen = (slong *) flint_malloc(order * sizeof(slong));

        for (slong i = 0; i < r; i++)
            for (slong j = 0; j < c; j++)
                dst[i * c + j] = nmod_poly_mat_entry(pmat, i, j)->coeffs;

        /* every source row is full: the coefficients of `matp` hold all the
           entries of a matrix, with no truncation */
        for (slong k = 0; k < order; k++)
            slen[k] = ndst;

        if (ndst == r * c)
        {
            for (slong k = 0; k < order; k++)
                src[k] = matp->coeffs[k];
            _pml_transpose(dst, r * c, src, slen, order, kern, dmaj);
        }
        else
            for (slong i = 0; i < r; i++)
            {
                for (slong k = 0; k < order; k++)
                    src[k] = matp->coeffs[k] + i * matp->stride;
                _pml_transpose(dst + i * c, c, src, slen, order, kern, dmaj);
            }

        flint_free(dst);
        flint_free(src);
        flint_free(slen);
    }

    // normalize
    for (slong i = 0; i < r; i++)
        for (slong j = 0; j < c; j++)
        {
            _nmod_poly_set_length(nmod_poly_mat_entry(pmat, i, j), order);
            _nmod_poly_normalise(nmod_poly_mat_entry(pmat, i, j));
        }
}

void nmod_poly_mat_set_trunc_from_mat_poly(nmod_poly_mat_t pmat,
                                           const nmod_mat_poly_t matp,
                                           slong order)
{
    _nmod_poly_mat_set_trunc_from_mat_poly(pmat, matp, order, -1, -1);
}
