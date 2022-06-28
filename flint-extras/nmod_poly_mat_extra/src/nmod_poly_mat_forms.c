#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_forms.h"


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* TESTING MATRIX FORMS                                       */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

int is_hermite(const nmod_poly_mat_t mat, orientation_t row_wise)
{
    slong rdim = mat->r, cdim = mat->c;
    slong deg_mat = nmod_poly_mat_degree(mat);

    if (row_wise)
    {
        slong shifts[cdim];
        for (slong i = 0; i < cdim; i++)
            shifts[i] = (cdim - i) * (deg_mat + 1);
        return is_popov(mat, shifts, row_wise, 0);
    }

    slong shifts[rdim];
    for (slong i = 0; i < rdim; i++)
        shifts[i] = i * (deg_mat + 1);
    return is_popov(mat, shifts, row_wise, 0);

}

int is_popov(const nmod_poly_mat_t mat, const slong *shifts, orientation_t row_wise, int ordered)
{
    if (!is_weak_popov(mat, shifts, row_wise, ordered))
        return 0;
    slong cdim = mat->c, rdim = mat->r, pivot_deg, d;
    nmod_poly_struct *P, *pivot;

    if (row_wise)
    {
        slong lead_pos[rdim];
        leading_positions(lead_pos, mat, shifts, row_wise);

        for(slong i = 0; i < rdim; i++)
        {
            pivot = nmod_poly_mat_entry(mat, i, lead_pos[i]);
            pivot_deg = nmod_poly_degree(pivot);
            if (nmod_poly_get_coeff_ui(pivot, pivot_deg) != 1)
                return 0;
            for(slong j = 0; j < rdim ; j++)
            {
                P = nmod_poly_mat_entry(mat, j, lead_pos[i]);
                d = nmod_poly_degree(P);
                if (d >= pivot_deg)
                    return 0;
            }
        }
        return 1;
    }

    slong lead_pos[cdim];
    leading_positions(lead_pos, mat, shifts, row_wise);

    for(slong i = 0; i < cdim; i++)
    {
        pivot = nmod_poly_mat_entry(mat, lead_pos[i], i);
        pivot_deg = nmod_poly_degree(pivot);
        if (nmod_poly_get_coeff_ui(pivot, pivot_deg) != 1)
            return 0;
        for(slong j = 0; j < cdim ; j++)
        {
            P = nmod_poly_mat_entry(mat, lead_pos[i], j);
            d = nmod_poly_degree(P);
            if (d >= pivot_deg)
                return 0;
        }
    }
    return 1;

}

int is_reduced(const nmod_poly_mat_t mat, const slong *shifts, orientation_t row_wise)
{
    slong rdim = mat->r, cdim = mat->c;
    nmod_mat_t B;
    nmod_mat_init(B, rdim, cdim, nmod_poly_mat_modulus(mat));
    leading_matrix_shifted(B, mat, shifts, row_wise);
    slong rank_lead = nmod_mat_rank(B);
    nmod_mat_clear(B);
    return (int) (rdim == rank_lead);
}

static int intComparator ( const void * first, const void * second ) {
    int firstInt = * (const int *) first;
    int secondInt = * (const int *) second;
    return firstInt - secondInt;
}

int is_weak_popov(const nmod_poly_mat_t mat, const slong *shifts, orientation_t row_wise, int ordered)
{
    if (!is_reduced(mat, shifts, row_wise))
        return 0;

    slong rdim = mat->r, cdim = mat->c;
    if (row_wise)
    {
        slong lead_pos[rdim];
        leading_positions(lead_pos, mat, shifts, row_wise);

        if (!ordered)
            qsort(lead_pos, rdim, sizeof(ulong), intComparator);

        for (slong i = 0; i < rdim - 1; i++)
        {
            if (lead_pos[i] > lead_pos[i+1])
                return 0;
        }
        return 1;
    }

    slong lead_pos[cdim];
    leading_positions(lead_pos, mat, shifts, row_wise);

    if (!ordered)
    {
        qsort(lead_pos, cdim, sizeof(ulong), intComparator);

        for (slong i = 0; i < cdim; i++)
        {
            if (lead_pos[i] > lead_pos[i+1])
                return 0;
        }
    }
    return 1;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
