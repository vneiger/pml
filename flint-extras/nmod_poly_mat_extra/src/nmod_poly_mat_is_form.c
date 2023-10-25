#include <stdlib.h> // qsort
#include "nmod_poly_mat_forms.h"

// helper function for qsort
static int _slong_comparator(const void * a, const void * b)
{
    if ( *(slong *)a  <  *(slong *)b )
        return -1;
    if ( *(slong *)a  == *(slong *)b )
        return 0;
    else // if ( *(slong *)a  >  *(slong *)b )
        return 1;
}

/**********************************************************************
*                              reduced                               *
**********************************************************************/

int nmod_poly_mat_is_reduced(const nmod_poly_mat_t mat,
                             const slong * shift,
                             orientation_t orient)
{
    nmod_mat_t lmat;
    nmod_mat_init(lmat, mat->r, mat->c, mat->modulus);
    nmod_poly_mat_leading_matrix(lmat, mat, shift, orient);
    slong rank_lead = nmod_mat_rank(lmat);
    nmod_mat_clear(lmat);
    if (orient == ROW_LOWER || orient == ROW_UPPER)
        return mat->r == rank_lead;
    else // orient == COL_*
        return mat->c == rank_lead;
}

/**********************************************************************
*                             weak Popov                             *
**********************************************************************/

int nmod_poly_mat_is_weak_popov(const nmod_poly_mat_t mat,
                                const slong * shift,
                                orientation_t orient)
{
    slong dim, zpiv;
    if (orient == ROW_LOWER || orient == ROW_UPPER)
        { dim = mat->r; zpiv = mat->c; }
    else // orient == COL_*
        { dim = mat->c; zpiv = mat->r; }

    slong * pivind = flint_malloc(dim * sizeof(slong));
    nmod_poly_mat_pivot_index(pivind, mat, shift, orient);

    // sort pivot indices in nondecreasing order
    qsort(pivind, dim, sizeof(slong), _slong_comparator);

    // pivot index must not indicate zero row/col
    if (dim > 0 && pivind[0] == -1)
        return 0;

    // indices of pivots must increase strictly
    for (slong i = 0; i < dim - 1; i++)
        if (pivind[i] >= pivind[i+1])
            return 0;

    // pivot index must not indicate zero row/col
    if (dim > 0 && pivind[dim-1] == zpiv)
        return 0;

    flint_free(pivind);
    return 1;
}

/**********************************************************************
*                         ordered weak Popov                         *
**********************************************************************/

int nmod_poly_mat_is_ordered_weak_popov(const nmod_poly_mat_t mat,
                                        const slong * shift,
                                        orientation_t orient)
{
    slong dim, zpiv;
    if (orient == ROW_LOWER || orient == ROW_UPPER)
        { dim = mat->r; zpiv = mat->c; }
    else // orient == COL_*
        { dim = mat->c; zpiv = mat->r; }

	slong * pivind = flint_malloc(dim * sizeof(slong));
    nmod_poly_mat_pivot_index(pivind, mat, shift, orient);

    // first row/col must be nonzero
    if (dim > 0 && pivind[0] == -1)
        return 0;

    // indices of pivots must increase strictly
    for (slong i = 0; i < dim - 1; i++)
        if (pivind[i] >= pivind[i+1])
            return 0;

    // last row must be nonzero
    if (dim > 0 && pivind[dim-1] == zpiv)
        return 0;

    flint_free(pivind);
    return 1;
}

/**********************************************************************
*                               Popov                                *
**********************************************************************/

int nmod_poly_mat_is_popov(const nmod_poly_mat_t mat,
                           const slong * shift,
                           orientation_t orient)
{
    slong dim;
    if (orient == ROW_LOWER || orient == ROW_UPPER)
        dim = mat->r;
    else // orient == COL_*
        dim = mat->c;

    // matrix must be ordered weak Popov
    if (!nmod_poly_mat_is_ordered_weak_popov(mat, shift, orient))
        return 0;

    // retrieve pivot profile
    slong * pivind = flint_malloc(dim * sizeof(slong));
    slong * pivdeg = flint_malloc(dim * sizeof(slong));
    nmod_poly_mat_pivot_profile(pivind, pivdeg, mat, shift, orient);

    if (orient == ROW_LOWER || orient == ROW_LOWER)
    {
        // pivots must be monic
        for (slong i = 0; i < dim; i++)
            if (!nmod_poly_is_monic(nmod_poly_mat_entry(mat, i, pivind[i])))
                return 0;

        // pivots must have degree strictly larger than other entries in same column
        for (slong i = 0; i < dim; i++)
            for (slong ii = 0; ii < dim; ii++)
                if (i != ii && nmod_poly_degree(nmod_poly_mat_entry(mat, ii, pivind[i])) >= pivdeg[i])
                    return 0;
    }
    else // orient == COL_*
    {
        // pivots must be monic
        for (slong j = 0; j < dim; j++)
            if (!nmod_poly_is_monic(nmod_poly_mat_entry(mat, pivind[j], j)))
                return 0;

        // pivots must have degree strictly larger than other entries in same column
        for (slong j = 0; j < dim; j++)
            for (slong jj = 0; jj < dim; jj++)
                if (j != jj && nmod_poly_degree(nmod_poly_mat_entry(mat, pivind[j], jj)) >= pivdeg[j])
                    return 0;
    }

    flint_free(pivind);
    flint_free(pivdeg);
    return 1;
}

/**********************************************************************
*                              echelon                               *
**********************************************************************/

int nmod_poly_mat_is_echelon(const nmod_poly_mat_t mat,
                             orientation_t orient)
{
    slong dim, zpiv;
    if (orient == ROW_LOWER || orient == ROW_UPPER)
        { dim = mat->r; zpiv = mat->c; }
    else // orient == COL_*
        { dim = mat->c; zpiv = mat->r; }

    // retrieve pivot index
    slong * pivind = flint_malloc(dim * sizeof(slong));
    nmod_poly_mat_echelon_pivot_index(pivind, mat, orient);

    if (orient == ROW_LOWER || orient == COL_UPPER)
    {
        // scan rows/cols until first nonzero one
        slong i = 0;
        while (i < dim && pivind[i] == -1)
            i++;

        // check remaining pivot indices are strictly increasing
        // in particular, if we have already reached dim, the matrix is zero;
        // if we have already reached the last row/col, it has a single nonzero
        // row/col; in both cases it is echelon
        while (i+1 < dim)
        {
            if (pivind[i] >= pivind[i+1])
                return 0;
            // note: the above test includes the case where row i+1 is zero (has index -1),
            // i.e. zero row below a nonzero row
            i++;
        }
    }
    else // ROW_UPPER || COL_LOWER
    {
        // scan rows from last one until finding nonzero one
        slong i = dim - 1;
        while (i >= 0 && pivind[i] == zpiv)
            i--;

        // check remaining pivot indices are strictly increasing
        // in particular, if we have reached i == -1, the matrix is zero; if we
        // have reached the first row i==0, it has a single nonzero row; in both
        // cases it is echelon
        while (i > 0)
        {
            if (pivind[i-1] >= pivind[i])
                return 0;
            // note: the above test includes the case where row i-1 is zero (has index other_dim),
            // i.e. zero row above a nonzero row
            i--;
        }
    }

    flint_free(pivind);
    return 1;
}

/**********************************************************************
*                              Hermite                               *
**********************************************************************/

int nmod_poly_mat_is_hermite(const nmod_poly_mat_t mat, orientation_t orient)
{
    slong dim, zpiv;
    if (orient == ROW_LOWER || orient == ROW_UPPER)
        { dim = mat->r; zpiv = mat->c; }
    else // orient == COL_*
        { dim = mat->c; zpiv = mat->r; }

    if (dim == 0)
        return 1;

    // matrix must be in echelon form
    if (! nmod_poly_mat_is_echelon(mat, orient))
        return 0;

    // retrieve pivot index
    slong * pivind = flint_malloc(dim * sizeof(slong));
    slong * pivdeg = flint_malloc(dim * sizeof(slong));
    nmod_poly_mat_echelon_pivot_profile(pivind, pivdeg, mat, orient);

    // make sure no zero row/col  (note len(pivind) == dim > 0)
    if (pivind[0] == -1 || pivind[dim-1] == zpiv)
        return 0;
    // this is sufficient since is_echelon already forbids zero rows/cols
    // appearing after nonzero one (or before, depending on orientation)

    if (orient == ROW_LOWER || orient == ROW_UPPER)
    {
        // scan rows starting from topmost one
        for (slong i = 0; i < dim; i++)
        {
            // pivot entry must be monic
            if (! nmod_poly_is_monic(nmod_poly_mat_entry(mat, i, pivind[i])))
                return 0;
            if (orient == ROW_LOWER) // entries below pivot entry must have lower degree
            {
                for (slong ii = i+1; ii < dim; ii++)
                    if (nmod_poly_mat_entry(mat, ii, pivind[i])->length > pivdeg[i])
                        return 0;
            }
            else // ROW_UPPER, entries above pivot entry must have lower degree
            {
                for (slong ii = 0; ii < i; ii++)
                    if (nmod_poly_mat_entry(mat, ii, pivind[i])->length > pivdeg[i])
                        return 0;
            }
        }
    }
    else // orient == COL_*
    {
        for (slong j = 0; j < dim; j++)
        {
            // pivot entry must be monic
            if (! nmod_poly_is_monic(nmod_poly_mat_entry(mat, pivind[j], j)))
                return 0;

            if (orient == COL_LOWER) // entries to the left of pivot must have lower degree
            {
                for (slong jj = 0; jj < j; jj++)
                    if (nmod_poly_mat_entry(mat, pivind[j], jj)->length > pivdeg[j])
                        return 0;
            }
            else // COL_UPPER, entries to the right of pivot entry must have lower degree
            {
                for (slong jj = j+1; jj < dim; jj++)
                    if (nmod_poly_mat_entry(mat, pivind[j], jj)->length > pivdeg[j])
                        return 0;
            }
        }
    }

    flint_free(pivind);
    flint_free(pivdeg);
    return 1;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
