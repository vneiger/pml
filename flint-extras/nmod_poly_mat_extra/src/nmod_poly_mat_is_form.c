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

int nmod_poly_mat_is_reduced_rowwise(const nmod_poly_mat_t mat,
                                     const slong *shift)
{
    nmod_mat_t lmat;
    nmod_mat_init(lmat, mat->r, mat->c, mat->modulus);
    nmod_poly_mat_leading_matrix_rowwise(lmat, mat, shift);
    slong rank_lead = nmod_mat_rank(lmat);
    nmod_mat_clear(lmat);
    return mat->r == rank_lead;
}

int nmod_poly_mat_is_reduced_columnwise(const nmod_poly_mat_t mat,
                                        const slong *shift)
{
    nmod_mat_t lmat;
    nmod_mat_init(lmat, mat->r, mat->c, mat->modulus);
    nmod_poly_mat_leading_matrix_columnwise(lmat, mat, shift);
    slong rank_lead = nmod_mat_rank(lmat);
    nmod_mat_clear(lmat);
    return mat->c == rank_lead;
}

/**********************************************************************
*                             weak Popov                             *
**********************************************************************/

int nmod_poly_mat_is_weak_popov_rowwise(const nmod_poly_mat_t mat,
                                        const slong * shift)
{
	slong * pivind = flint_malloc(mat->r * sizeof(slong));
	nmod_poly_mat_pivot_index_rowwise(pivind, mat, shift);

	// sort pivot indices in nondecreasing order
	qsort(pivind, mat->r, sizeof(ulong), _slong_comparator);

	// smallest pivot index must not indicate zero row
	if (mat->r > 0 && pivind[0] == -1)
		return 0;

	// indices of pivots must increase strictly
	for (slong i = 0; i < mat->r - 1; i++)
		if (pivind[i] >= pivind[i+1])
			return 0;

    flint_free(pivind);
	return 1;
}

int nmod_poly_mat_is_weak_popov_columnwise(const nmod_poly_mat_t mat,
                                           const slong * shift)
{
	slong * pivind = flint_malloc(mat->c * sizeof(slong));
	nmod_poly_mat_pivot_index_columnwise(pivind, mat, shift);

	// sort pivot indices in nondecreasing order
	qsort(pivind, mat->c, sizeof(ulong), _slong_comparator);

	// smallest pivot index must not indicate zero column
	if (mat->c > 0 && pivind[0] == -1)
		return 0;

	// indices of pivots must increase strictly
	for (slong i = 0; i < mat->c - 1; i++)
		if (pivind[i] >= pivind[i+1])
			return 0;

    flint_free(pivind);
	return 1;
}

/**********************************************************************
*                         ordered weak Popov                         *
**********************************************************************/

int nmod_poly_mat_is_ordered_weak_popov_rowwise(const nmod_poly_mat_t mat,
                                                const slong * shift)
{
	slong * pivind = flint_malloc(mat->r * sizeof(slong));
    nmod_poly_mat_pivot_index_rowwise(pivind, mat, shift);

    // first row must be nonzero
    if (mat->r > 0 && pivind[0] == -1)
        return 0;

    // indices of pivots must increase strictly
    for (slong i = 0; i < mat->r - 1; i++)
        if (pivind[i] >= pivind[i+1])
            return 0;

    flint_free(pivind);
    return 1;
}

int nmod_poly_mat_is_ordered_weak_popov_columnwise(const nmod_poly_mat_t mat,
                                                   const slong * shift)
{
	slong * pivind = flint_malloc(mat->c * sizeof(slong));
    nmod_poly_mat_pivot_index_columnwise(pivind, mat, shift);

    // first column must be nonzero
    if (mat->c > 0 && pivind[0] == -1)
        return 0;

    // indices of pivots must increase strictly
    for (slong i = 0; i < mat->c - 1; i++)
        if (pivind[i] >= pivind[i+1])
            return 0;

    flint_free(pivind);
    return 1;
}

/**********************************************************************
*                               Popov                                *
**********************************************************************/

int nmod_poly_mat_is_popov_rowwise(const nmod_poly_mat_t mat,
                                   const slong * shift)
{
    // matrix must be ordered weak Popov
    if (!nmod_poly_mat_is_ordered_weak_popov_rowwise(mat, shift))
        return 0;

    // retrieve pivot profile
    slong * pivind = flint_malloc(mat->r * sizeof(slong));
    slong * pivdeg = flint_malloc(mat->r * sizeof(slong));
    nmod_poly_mat_pivot_profile_rowwise(pivind, pivdeg, mat, shift);

    // pivots must be monic
    for (slong i = 0; i < mat->r; i++)
        if (!nmod_poly_is_monic(nmod_poly_mat_entry(mat, i, pivind[i])))
            return 0;

    // pivots must have degree strictly larger than other entries in same column
    for (slong i = 0; i < mat->r; i++)
        for (slong ii = 0; ii < mat->r; ii++)
            if (i != ii && nmod_poly_degree(nmod_poly_mat_entry(mat, ii, pivind[i])) >= pivdeg[i])
                return 0;

    flint_free(pivind);
    flint_free(pivdeg);
    return 1;
}

int nmod_poly_mat_is_popov_columnwise(const nmod_poly_mat_t mat,
                                      const slong * shift)
{
    // matrix must be ordered weak Popov
    if (!nmod_poly_mat_is_ordered_weak_popov_columnwise(mat, shift))
        return 0;

    // retrieve pivot profile
    slong * pivind = flint_malloc(mat->c * sizeof(slong));
    slong * pivdeg = flint_malloc(mat->c * sizeof(slong));
    nmod_poly_mat_pivot_profile_columnwise(pivind, pivdeg, mat, shift);

    // pivots must be monic
    for (slong j = 0; j < mat->c; j++)
        if (!nmod_poly_is_monic(nmod_poly_mat_entry(mat, pivind[j], j)))
            return 0;

    // pivots must have degree strictly larger than other entries in same column
    for (slong j = 0; j < mat->c; j++)
        for (slong jj = 0; jj < mat->c; jj++)
            if (j != jj && nmod_poly_degree(nmod_poly_mat_entry(mat, pivind[j], jj)) >= pivdeg[j])
                return 0;

    flint_free(pivind);
    flint_free(pivdeg);
    return 1;
}



/**********************************************************************
*                              echelon                               *
**********************************************************************/

int nmod_poly_mat_is_lechelon_rowwise(const nmod_poly_mat_t mat)
{
    // retrieve pivot index
    slong * pivind = flint_malloc(mat->r * sizeof(slong));
    nmod_poly_mat_lechelon_pivot_index_rowwise(pivind, mat);

    // scan rows until first nonzero one
    slong i = 0;
    while (i < mat->r && pivind[i] == -1)
        i++;

    // check remaining pivot indices are strictly increasing
    // in particular, if we have already reached mat->r, the matrix is zero; if
    // we have already reached the last row, it has a single nonzero row; in
    // both cases it is echelon
    while (i+1 < mat->r)
    {
        if (pivind[i] >= pivind[i+1])
            return 0;
        // note: the above test includes the case where row i+1 is zero,
        // i.e. zero row below a nonzero row
        i++;
    }

    flint_free(pivind);

    return 1;
}

int nmod_poly_mat_is_uechelon_rowwise(const nmod_poly_mat_t mat)
{
    // retrieve pivot index
    slong * pivind = flint_malloc(mat->r * sizeof(slong));
    nmod_poly_mat_uechelon_pivot_index_rowwise(pivind, mat);

    // scan rows from last one until finding nonzero one
    slong i = mat->r - 1;
    while (i >= 0 && pivind[i] == mat->c)
        i--;

    // check remaining pivot indices are strictly increasing
    // in particular, if we have reached i == -1, the matrix is zero; if we
    // have reached the first row i==0, it has a single nonzero row; in both
    // cases it is echelon
    while (i > 0)
    {
        if (pivind[i-1] >= pivind[i])
            return 0;
        // note: the above test includes the case where row i-1 is zero (has index mat->c),
        // i.e. zero row above a nonzero row
        i--;
    }

    flint_free(pivind);

    return 1;
}

int nmod_poly_mat_is_lechelon_columnwise(const nmod_poly_mat_t mat)
{
    // retrieve pivot index
    slong * pivind = flint_malloc(mat->c * sizeof(slong));
    nmod_poly_mat_lechelon_pivot_index_columnwise(pivind, mat);

    // scan columns from last one until finding nonzero one
    slong j = mat->c - 1;
    while (j >= 0 && pivind[j] == mat->r)
        j--;

    // check remaining pivot indices are strictly increasing
    // in particular, if we have reached j == -1, the matrix is zero; if we
    // have reached the first column j==0, it has a single nonzero column; in both
    // cases it is echelon
    while (j > 0)
    {
        if (pivind[j-1] >= pivind[j])
            return 0;
        // note: the above test includes the case where column j-1 is zero (has index mat->r),
        // i.e. zero column to the left of a nonzero column
        j--;
    }

    flint_free(pivind);

    return 1;
}

int nmod_poly_mat_is_uechelon_columnwise(const nmod_poly_mat_t mat)
{
    // retrieve pivot index
    slong * pivind = flint_malloc(mat->c * sizeof(slong));
    nmod_poly_mat_lechelon_pivot_index_columnwise(pivind, mat);

    // scan columns until first nonzero one
    slong j = 0;
    while (j < mat->c && pivind[j] == -1)
        j++;

    // otherwise, check remaining pivot indices are strictly increasing
    // in particular: if we have already reached mat->c, the matrix is zero; if
    // we have already reached the last column, it has a single nonzero column;
    // in both cases it is echelon
    while (j+1 < mat->c)
    {
        if (pivind[j] >= pivind[j+1])
            return 0;
        // note: the above test includes the case where column j+1 is zero,
        // i.e. zero column to the right of a nonzero one
        j++;
    }

    flint_free(pivind);

    return 1;
}

/**********************************************************************
*                              Hermite                               *
**********************************************************************/

int nmod_poly_mat_is_lhermite_rowwise(const nmod_poly_mat_t mat)
{
    if (mat->r == 0)
        return 1;

    // matrix must be in lower echelon form row-wise, with no zero row
    // retrieve pivot index
    slong * pivind = flint_malloc(mat->r * sizeof(slong));
    nmod_poly_mat_lechelon_pivot_index_rowwise(pivind, mat);

    // first row: make sure it is nonzero  (note len(pivind) == mat->r > 0)
    if (pivind[0] == -1)
        return 0;

    nmod_poly_struct * pivot;
    // scan rows starting from topmost one
    for (slong i = 0; i < mat->r - 1; i++)
    {
        // row i+1 must be nonzero and pivot indices must be increasing
        if (pivind[i] >= pivind[i+1])
            return 0;
        // note: the above test includes the case where row i+1 is zero,
        // i.e. zero row below a nonzero row

        pivot = nmod_poly_mat_entry(mat, i, pivind[i]);
        // pivot entry must be monic
        if (! nmod_poly_is_monic(pivot))
            return 0;

        // entries below pivot entry must have lower degree
        for (slong ii = i+1; ii < mat->r; ii++)
            if (nmod_poly_mat_entry(mat, ii, pivind[i])->length >= pivot->length)
                return 0;
    }

    // last row: pivot must be monic, and that's it
    if (! nmod_poly_is_monic(nmod_poly_mat_entry(mat, mat->r -1, pivind[mat->r -1])))
        return 0;

    flint_free(pivind);

    return 1;
}

int nmod_poly_mat_is_uhermite_rowwise(const nmod_poly_mat_t mat)
{
    if (mat->r == 0)
        return 1;

    // matrix must be in upper echelon form row-wise, with no zero row
    // retrieve pivot index
    slong * pivind = flint_malloc(mat->r * sizeof(slong));
    nmod_poly_mat_uechelon_pivot_index_rowwise(pivind, mat);

    // last row: make sure it is nonzero  (note len(pivind) == mat->r > 0)
    if (pivind[mat->r - 1] == mat->c)
        return 0;

    nmod_poly_struct * pivot;
    // scan rows starting from bottommost one
    for (slong i = mat->r -1; i > 0; i--)
    {
        // row i-1 must be nonzero and pivot indices must be increasing
        if (pivind[i-1] >= pivind[i])
            return 0;
        // note: the above test includes the case where row i-1 is zero,
        // i.e. zero row above a nonzero row

        pivot = nmod_poly_mat_entry(mat, i, pivind[i]);
        // pivot entry must be monic
        if (! nmod_poly_is_monic(pivot))
            return 0;

        // entries above pivot entry must have lower degree
        for (slong ii = 0; ii < i; ii++)
            if (nmod_poly_mat_entry(mat, ii, pivind[i])->length >= pivot->length)
                return 0;
    }

    // first row: pivot must be monic, and that's it
    if (! nmod_poly_is_monic(nmod_poly_mat_entry(mat, 0, pivind[0])))
        return 0;

    flint_free(pivind);

    return 1;
}

int nmod_poly_mat_is_lhermite_columnwise(const nmod_poly_mat_t mat)
{
    if (mat->c == 0)
        return 1;

    // matrix must be in lower echelon form column-wise, with no zero column
    // retrieve pivot index
    slong * pivind = flint_malloc(mat->c * sizeof(slong));
    nmod_poly_mat_lechelon_pivot_index_columnwise(pivind, mat);

    // last column: make sure it is nonzero  (note len(pivind) == mat->c > 0)
    if (pivind[mat->c - 1] == mat->r)
        return 0;

    nmod_poly_struct * pivot;
    // scan columns starting from rightmost one
    for (slong j = mat->c -1; j > 0; j--)
    {
        // column j-1 must be nonzero and pivot indices must be increasing
        if (pivind[j-1] >= pivind[j])
            return 0;
        // note: the above test includes the case where column j-1 is zero,
        // i.e. zero column to the left of a nonzero column

        pivot = nmod_poly_mat_entry(mat, pivind[j], j);
        // pivot entry must be monic
        if (! nmod_poly_is_monic(pivot))
            return 0;

        // entries to the left of pivot entry must have lower degree
        for (slong jj = 0; jj < j; jj++)
            if (nmod_poly_mat_entry(mat, pivind[j], jj)->length >= pivot->length)
                return 0;
    }

    // first column: pivot must be monic, and that's it
    if (! nmod_poly_is_monic(nmod_poly_mat_entry(mat, pivind[0], 0)))
        return 0;

    flint_free(pivind);

    return 1;
}

int nmod_poly_mat_is_uhermite_columnwise(const nmod_poly_mat_t mat)
{
    if (mat->r == 0)
        return 1;

    // matrix must be in upper echelon form column-wise, with no zero column
    // retrieve pivot index
    slong * pivind = flint_malloc(mat->c * sizeof(slong));
    nmod_poly_mat_uechelon_pivot_index_columnwise(pivind, mat);

    // first column: make sure it is nonzero  (note len(pivind) == mat->c > 0)
    if (pivind[0] == -1)
        return 0;

    nmod_poly_struct * pivot;
    // scan columns starting from leftmost one
    for (slong j = 0; j < mat->c - 1; j++)
    {
        // column j+1 must be nonzero and pivot indices must be increasing
        if (pivind[j] >= pivind[j+1])
            return 0;
        // note: the above test includes the case where column j+1 is zero,
        // i.e. zero column to the right of a nonzero column

        pivot = nmod_poly_mat_entry(mat, pivind[j], j);
        // pivot entry must be monic
        if (! nmod_poly_is_monic(pivot))
            return 0;

        // entries to the right of pivot entry must have lower degree
        for (slong jj = j+1; jj < mat->c; jj++)
            if (nmod_poly_mat_entry(mat, pivind[j], jj)->length >= pivot->length)
                return 0;
    }

    // last column: pivot must be monic, and that's it
    if (! nmod_poly_is_monic(nmod_poly_mat_entry(mat, pivind[mat->c -1], mat->c -1)))
        return 0;

    flint_free(pivind);

    return 1;
}


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
