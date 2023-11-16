#include <stdlib.h> // for qsort
#include <flint/nmod_poly.h>
#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_forms.h"
#include "nmod_mat_poly.h"


/**********************************************************************
*                    ROW ROTATION DOWNWARD/UPWARD                    *
**********************************************************************/

void _nmod_poly_mat_rotate_rows_downward(nmod_poly_mat_t mat, slong * vec, slong i, slong j)
{
    if (i != j)
    {
        if (vec)
        {
            slong tmp_vec = vec[j];
            for (slong ii = j; ii > i; ii--)
                vec[ii] = vec[ii-1];
            vec[i] = tmp_vec;
        }

        nmod_poly_struct * tmp_mat = mat->rows[j];
        for (slong ii = j; ii > i; ii--)
            mat->rows[ii] = mat->rows[ii-1];
        mat->rows[i] = tmp_mat;
    }
}

void _nmod_poly_mat_rotate_rows_upward(nmod_poly_mat_t mat, slong * vec, slong i, slong j)
{
    if (i != j)
    {
        if (vec)
        {
            slong tmp_vec = vec[i];
            for (slong ii = i; ii < j; ii++)
                vec[ii] = vec[ii+1];
            vec[j] = tmp_vec;
        }

        nmod_poly_struct * tmp_mat = mat->rows[i];
        for (slong ii = i; ii < j; ii++)
            mat->rows[ii] = mat->rows[ii+1];
        mat->rows[j] = tmp_mat;
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


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* SET FROM MATRIX POLYNOMIAL                                 */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

void nmod_poly_mat_set_from_mat_poly0(nmod_poly_mat_t pmat, const nmod_mat_poly0_t matp)
{
	nmod_poly_mat_zero(pmat);
	for (slong k = 0; k <= matp->degree; k++)
        for (slong i = 0; i < matp->r; ++i)
            for (slong j = 0; j < matp->c; ++j)
                nmod_poly_set_coeff_ui(nmod_poly_mat_entry(pmat, i, j), k, nmod_mat_get_entry(matp->mat + k, i, j));
}

void nmod_poly_mat_set_trunc_from_mat_poly(nmod_poly_mat_t pmat,
                                           const nmod_mat_poly_t matp,
                                           slong order)
{
    if (order > matp->length)
        order = matp->length;

    // prepare memory
    for (slong i = 0; i < pmat->r; i++)
        for (slong j = 0; j < pmat->c; j++)
            nmod_poly_fit_length(nmod_poly_mat_entry(pmat, i, j), order);

    // fill data
    for (slong k = 0; k < order; k++)
        for (slong i = 0; i < pmat->r; i++)
            for (slong j = 0; j < pmat->c; j++)
                nmod_poly_mat_entry(pmat, i, j)->coeffs[k] = nmod_mat_poly_entry(matp, k, i, j);

    // normalize
    for (slong i = 0; i < pmat->r; i++)
        for (slong j = 0; j < pmat->c; j++)
        {
            _nmod_poly_set_length(nmod_poly_mat_entry(pmat, i, j), order);
            _nmod_poly_normalise(nmod_poly_mat_entry(pmat, i, j));
        }
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
