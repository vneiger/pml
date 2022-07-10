#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_forms.h"

// helper function for qsort
static int _int_comparator ( const void * first, const void * second ) {
    int firstInt = * (const int *) first;
    int secondInt = * (const int *) second;
    return firstInt - secondInt;
}

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* TESTING MATRIX FORMS - REDUCED                             */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

int is_row_reduced(const nmod_poly_mat_t mat)
{
    nmod_mat_t lmat;
    nmod_mat_init(lmat, mat->r, mat->c, mat->modulus);
    nmod_poly_mat_leading_matrix_rowwise(lmat, mat);
    slong rank_lead = nmod_mat_rank(lmat);
    nmod_mat_clear(lmat);
    return mat->r == rank_lead;
}

int is_row_reduced_shifted(const nmod_poly_mat_t mat, const slong *shift)
{
    nmod_mat_t lmat;
    nmod_mat_init(lmat, mat->r, mat->c, mat->modulus);
    nmod_poly_mat_leading_matrix_shifted_rowwise(lmat, mat, shift);
    slong rank_lead = nmod_mat_rank(lmat);
    nmod_mat_clear(lmat);
    return mat->r == rank_lead;
}

int is_column_reduced(const nmod_poly_mat_t mat)
{
    nmod_mat_t lmat;
    nmod_mat_init(lmat, mat->r, mat->c, mat->modulus);
    nmod_poly_mat_leading_matrix_columnwise(lmat, mat);
    slong rank_lead = nmod_mat_rank(lmat);
    nmod_mat_clear(lmat);
    return mat->c == rank_lead;
}

int is_column_reduced_shifted(const nmod_poly_mat_t mat, const slong *shift)
{
    nmod_mat_t lmat;
    nmod_mat_init(lmat, mat->r, mat->c, mat->modulus);
    nmod_poly_mat_leading_matrix_shifted_columnwise(lmat, mat, shift);
    slong rank_lead = nmod_mat_rank(lmat);
    nmod_mat_clear(lmat);
    return mat->c == rank_lead;
}

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* TESTING MATRIX FORMS - ORDERED WEAK POPOV                  */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

int nmod_poly_mat_is_ordered_weak_popov_rowwise(const nmod_poly_mat_t mat)
{
    slong pivind[mat->r];
    nmod_poly_mat_pivot_index_rowwise(pivind, mat);

    // first row must be nonzero
    if (mat->r > 0 && pivind[0] == -1)
        return 0;

    // indices of pivots must increase strictly
    for (slong i = 0; i < mat->r - 1; i++)
        if (pivind[i] >= pivind[i+1])
            return 0;

    return 1;
}

int nmod_poly_mat_is_ordered_weak_popov_shifted_rowwise(const nmod_poly_mat_t mat,
                                          const slong *shift)
{
    slong pivind[mat->r];
    nmod_poly_mat_pivot_index_shifted_rowwise(pivind, mat, shift);

    // first row must be nonzero
    if (mat->r > 0 && pivind[0] == -1)
        return 0;

    // indices of pivots must increase strictly
    for (slong i = 0; i < mat->r - 1; i++)
        if (pivind[i] >= pivind[i+1])
            return 0;

    return 1;
}

int nmod_poly_mat_is_ordered_weak_popov_columnwise(const nmod_poly_mat_t mat)
{
    slong pivind[mat->c];
    nmod_poly_mat_pivot_index_columnwise(pivind, mat);

    // first column must be nonzero
    if (mat->c > 0 && pivind[0] == -1)
        return 0;

    // indices of pivots must increase strictly
    for (slong i = 0; i < mat->c - 1; i++)
        if (pivind[i] >= pivind[i+1])
            return 0;

    return 1;
}

int nmod_poly_mat_is_ordered_weak_popov_shifted_columnwise(const nmod_poly_mat_t mat,
                                             const slong *shift)
{
    slong pivind[mat->c];
    nmod_poly_mat_pivot_index_shifted_columnwise(pivind, mat, shift);

    // first column must be nonzero
    if (mat->c > 0 && pivind[0] == -1)
        return 0;

    // indices of pivots must increase strictly
    for (slong i = 0; i < mat->c - 1; i++)
        if (pivind[i] >= pivind[i+1])
            return 0;

    return 1;
}


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* TESTING MATRIX FORMS - WEAK POPOV                          */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

int nmod_poly_mat_is_weak_popov_rowwise(const nmod_poly_mat_t mat)
{
    slong pivind[mat->r];
    nmod_poly_mat_pivot_index_rowwise(pivind, mat);

    // sort pivot indices in nondecreasing order
    qsort(pivind, mat->r, sizeof(ulong), _int_comparator);

    // smallest pivot index must not indicate zero row
    if (mat->r > 0 && pivind[0] == -1)
        return 0;

    // indices of pivots must increase strictly
    for (slong i = 0; i < mat->r - 1; i++)
        if (pivind[i] >= pivind[i+1])
            return 0;

    return 1;
}

int nmod_poly_mat_is_weak_popov_shifted_rowwise(const nmod_poly_mat_t mat,
                                  const slong *shift)
{
    slong pivind[mat->r];
    nmod_poly_mat_pivot_index_shifted_rowwise(pivind, mat, shift);

    // sort pivot indices in nondecreasing order
    qsort(pivind, mat->r, sizeof(ulong), _int_comparator);

    // smallest pivot index must not indicate zero row
    if (mat->r > 0 && pivind[0] == -1)
        return 0;

    // indices of pivots must increase strictly
    for (slong i = 0; i < mat->r - 1; i++)
        if (pivind[i] >= pivind[i+1])
            return 0;

    return 1;
}

int nmod_poly_mat_is_weak_popov_columnwise(const nmod_poly_mat_t mat)
{
    slong pivind[mat->c];
    nmod_poly_mat_pivot_index_columnwise(pivind, mat);

    // sort pivot indices in nondecreasing order
    qsort(pivind, mat->c, sizeof(ulong), _int_comparator);

    // smallest pivot index must not indicate zero column
    if (mat->c > 0 && pivind[0] == -1)
        return 0;

    // indices of pivots must increase strictly
    for (slong i = 0; i < mat->c - 1; i++)
        if (pivind[i] >= pivind[i+1])
            return 0;

    return 1;
}

int nmod_poly_mat_is_weak_popov_shifted_columnwise(const nmod_poly_mat_t mat,
                                     const slong *shift)
{
    slong pivind[mat->c];
    nmod_poly_mat_pivot_index_shifted_columnwise(pivind, mat, shift);

    // sort pivot indices in nondecreasing order
    qsort(pivind, mat->c, sizeof(ulong), _int_comparator);

    // smallest pivot index must not indicate zero column
    if (mat->c > 0 && pivind[0] == -1)
        return 0;

    // indices of pivots must increase strictly
    for (slong i = 0; i < mat->c - 1; i++)
        if (pivind[i] >= pivind[i+1])
            return 0;

    return 1;
}


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* TESTING MATRIX FORMS - POPOV                               */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

//int is_popov(const nmod_poly_mat_t mat, const slong *shifts, orientation_t row_wise, int ordered)
//{
//    if (!nmod_poly_mat_is_weak_popov(mat, shifts, row_wise, ordered))
//        return 0;
//    slong cdim = mat->c, rdim = mat->r, pivot_deg, d;
//    nmod_poly_struct *P, *pivot;
//
//    if (row_wise)
//    {
//        slong lead_pos[rdim];
//        leading_positions(lead_pos, mat, shifts, row_wise);
//
//        for(slong i = 0; i < rdim; i++)
//        {
//            pivot = nmod_poly_mat_entry(mat, i, lead_pos[i]);
//            pivot_deg = nmod_poly_degree(pivot);
//            if (nmod_poly_get_coeff_ui(pivot, pivot_deg) != 1)
//                return 0;
//            for(slong j = 0; j < rdim ; j++)
//            {
//                P = nmod_poly_mat_entry(mat, j, lead_pos[i]);
//                d = nmod_poly_degree(P);
//                if (d >= pivot_deg)
//                    return 0;
//            }
//        }
//        return 1;
//    }
//
//    slong lead_pos[cdim];
//    leading_positions(lead_pos, mat, shifts, row_wise);
//
//    for(slong i = 0; i < cdim; i++)
//    {
//        pivot = nmod_poly_mat_entry(mat, lead_pos[i], i);
//        pivot_deg = nmod_poly_degree(pivot);
//        if (nmod_poly_get_coeff_ui(pivot, pivot_deg) != 1)
//            return 0;
//        for(slong j = 0; j < cdim ; j++)
//        {
//            P = nmod_poly_mat_entry(mat, lead_pos[i], j);
//            d = nmod_poly_degree(P);
//            if (d >= pivot_deg)
//                return 0;
//        }
//    }
//    return 1;
//
//}

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* TESTING MATRIX FORMS - HERMITE                             */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

//int is_hermite(const nmod_poly_mat_t mat, orientation_t row_wise)
//{
//    slong rdim = mat->r, cdim = mat->c;
//    slong deg_mat = nmod_poly_mat_degree(mat);
//
//    if (row_wise)
//    {
//        slong shifts[cdim];
//        for (slong i = 0; i < cdim; i++)
//            shifts[i] = (cdim - i) * (deg_mat + 1);
//        return is_popov(mat, shifts, row_wise, 0);
//    }
//
//    slong shifts[rdim];
//    for (slong i = 0; i < rdim; i++)
//        shifts[i] = i * (deg_mat + 1);
//    return is_popov(mat, shifts, row_wise, 0);
//
//}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
