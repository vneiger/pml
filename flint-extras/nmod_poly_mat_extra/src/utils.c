#include "nmod_poly_mat_utils.h"

void coefficient_matrix(nmod_mat_t res, const nmod_poly_mat_t mat, slong degree)
{
    slong rdim = mat->r, cdim = mat->c;
    for(slong i = 0; i < rdim; i++)
        for(slong j = 0; j < cdim; j++)
            nmod_mat_set_entry(res, i, j,
                               nmod_poly_get_coeff_ui(nmod_poly_mat_entry(mat, i, j), degree));
}

void column_degrees(int64_t *res, const nmod_poly_mat_t mat, const int64_t *shifts)
{
    slong rdim = mat->r, cdim = mat->c;
    slong max, d;
    for (slong i = 0; i < cdim; i++)
    {
        max = -1;
        for (slong j = 0; j < rdim; j++)
        {
            d = nmod_poly_degree(nmod_poly_mat_entry(mat, j, i)) + shifts[j];
            if (max < d)
                max = d;
        }
        res[i] = max;
    }
}

void row_degrees(int64_t *res, const nmod_poly_mat_t mat, const int64_t *shifts)
{
    slong rdim = mat->r, cdim = mat->c;
    slong max, d;
    for (slong i = 0; i < rdim; i++)
    {
        max = -1;
        for (slong j = 0; j < cdim; j++)
        {
            d = nmod_poly_degree(nmod_poly_mat_entry(mat, i, j)) + shifts[j];
            if (max < d)
                max = d;
        }
        res[i] = max;
    }
}

slong nmod_poly_mat_degree(const nmod_poly_mat_t mat)
{
    slong rdim = mat->r, cdim = mat->c;
    slong d, max = 0;
    for(slong i = 0; i < rdim; i++)
        for(slong j = 0; j < cdim; j++)
        {
            d = nmod_poly_degree(nmod_poly_mat_entry(mat, i, j));
            if (max < d)
                max = d;

        }
    return max;
}

void int64_mat_print(const int64_t *mat, slong rdim, slong cdim)
{
    printf("[");
    for (slong i = 0; i < rdim; i++)
    {
        printf("[");
        for (slong j = 0; j < cdim; j++)
        {
            if (j != cdim - 1)
                printf("%ld, ", *(mat + (i * rdim) + j));
            else
                printf("%ld", *(mat + (i * rdim) + j));

        }
        if (i != rdim -1)
            printf("],\n");
        else
            printf("]");
    }
    printf("]\n");
}

void degree_matrix(int64_t *res, const nmod_poly_mat_t mat, const int64_t *shifts,
                   matrix_wise row_wise)
{

    slong rdim = mat->r, cdim = mat->c;
    slong d;

    if (row_wise)
    {
        for(slong i = 0; i < rdim; i++)
        {
            for(slong j = 0; j < cdim; j++)
            {

                d = nmod_poly_degree(nmod_poly_mat_entry(mat, i, j)) + shifts[j];
                *(res + (i * rdim) + j) = d;
            }
        }
        return;
    }
    for(slong i = 0; i < cdim; i++)
        for(slong j = 0; j < rdim; j++)
        {
            d = nmod_poly_degree(nmod_poly_mat_entry(mat, j, i)) + shifts[j];
            *(res + (i * cdim) + j) = d;
        }
}

void leading_matrix(nmod_mat_t res, const nmod_poly_mat_t mat,
                    const int64_t *shifts, matrix_wise row_wise)
{
    slong rdim = mat->r, cdim = mat->c;
    nmod_poly_struct *P;
    if (row_wise)
    {
        int64_t rdeg[rdim];
        row_degrees(rdeg, mat, shifts);
        for(slong i = 0; i < rdim; i++)
            for(slong j = 0; j < cdim; j++)
            {
                P = nmod_poly_mat_entry(mat, i, j);
                nmod_mat_set_entry(res, i, j, nmod_poly_get_coeff_ui(P, rdeg[i] - shifts[j]));
            }
        return;
    }
    {
        int64_t cdeg[cdim];
        column_degrees(cdeg, mat, shifts);
        for(slong i = 0; i < cdim; i++)
            for(slong j = 0; j < rdim; j++)
            {
                P = nmod_poly_mat_entry(mat, j, i);
                nmod_mat_set_entry(res, j, i, nmod_poly_get_coeff_ui(P, cdeg[i] - shifts[j]));
            }
    }
}

void leading_positions(int64_t *res, const nmod_poly_mat_t mat,
                       const int64_t *shifts, matrix_wise row_wise)
{
    slong rdim = mat->r, cdim = mat->c;
    slong max;
    slong d;
    int64_t ind; // TODO initial value
    if (row_wise)
    {
        for (slong i = 0; i < rdim; i++)
        {
            max = 0;
            for (slong j = 0; j < cdim; j++)
            {
                d = nmod_poly_degree(nmod_poly_mat_entry(mat, i, j)) + shifts[j];
                if (max <= d)
                {
                    max = d;
                    ind = j;
                }
            }
            res[i] = ind;
        }
        return;
    }
    {
        for (slong i = 0; i < cdim; i++)
        {
            max = 0;
            for (slong j = 0; j < rdim; j++)
            {
                d = nmod_poly_degree(nmod_poly_mat_entry(mat, j, i)) + shifts[j];
                if (max <= d)
                {
                    max = d;
                    ind = j;
                }
            }
            res[i] = ind;
        }
    }
}

int is_hermite(const nmod_poly_mat_t mat, matrix_wise row_wise)
{
    slong rdim = mat->r, cdim = mat->c;
    slong deg_mat = nmod_poly_mat_degree(mat);

    if (row_wise)
    {
        int64_t shifts[cdim];
        for (slong i = 0; i < cdim; i++)
            shifts[i] = (cdim - i) * (deg_mat + 1);
        return is_popov(mat, shifts, row_wise, 0);
    }

    int64_t shifts[rdim];
    for (slong i = 0; i < rdim; i++)
        shifts[i] = i * (deg_mat + 1);
    return is_popov(mat, shifts, row_wise, 0);

}

int is_popov(const nmod_poly_mat_t mat, const int64_t *shifts, matrix_wise row_wise, int ordered)
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

int is_reduced(const nmod_poly_mat_t mat, const int64_t *shifts, matrix_wise row_wise)
{
    slong rdim = mat->r, cdim = mat->c;
    nmod_mat_t B;
    nmod_mat_init(B, rdim, cdim, nmod_poly_mat_modulus(mat));
    leading_matrix(B, mat, shifts, row_wise);
    slong rank_lead = nmod_mat_rank(B);
    nmod_mat_clear(B);
    return (int) (rdim == rank_lead);
}

static int intComparator ( const void * first, const void * second ) {
    int firstInt = * (const int *) first;
    int secondInt = * (const int *) second;
    return firstInt - secondInt;
}

int is_weak_popov(const nmod_poly_mat_t mat, const int64_t *shifts, matrix_wise row_wise, int ordered)
{
    if (!is_reduced(mat, shifts, row_wise))
        return 0;

    slong rdim = mat->r, cdim = mat->c;
    if (row_wise)
    {
        int64_t lead_pos[rdim];
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

    int64_t lead_pos[cdim];
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

int is_zero_mod_xk(const nmod_poly_mat_t mat, int64_t k)
{
    nmod_poly_t P;
    nmod_poly_init(P, mat->modulus);
    for(slong i = 0; i < mat->r; i++)
        for(slong j = 0; j < mat->c; j++)
        {
            nmod_poly_set(P, nmod_poly_mat_entry(mat, i, j));
            nmod_poly_shift_right(P, P, k);
            if (!nmod_poly_is_zero(P))
                return 0;
        }
    return 1;
}


/** TO FIX **/
int is_minimal_approximant_basis(const nmod_poly_mat_t base,
                                 const nmod_mat_t mat, int64_t order,
                                 const int64_t *shifts)
{
    slong rdim = mat->r, cdim = mat->c;
    mp_limb_t prime = mat->mod.n;
    nmod_poly_t constant;
    nmod_poly_mat_t mat_poly, res_mul;

    if (base->c != base->r)
    {
        printf("not basis: wrong shape");
        return 0;
    }

    nmod_poly_mat_init(mat_poly, rdim, cdim, prime);
    slong alloc;
    nmod_poly_init(constant, prime);
    for (slong i = 0; i < rdim; i++)
        for (slong j = 0; j < cdim; j++)
        {
            alloc = (slong) nmod_mat_get_entry(mat, i, j);
            nmod_poly_set_coeff_ui(constant, 0, alloc);
            nmod_poly_set(nmod_poly_mat_entry(mat_poly, i, j), constant);
        }
    nmod_poly_mat_init(res_mul, rdim, cdim, prime);
    nmod_poly_mat_mul(res_mul, base, mat_poly);

    if (! is_zero_mod_xk(res_mul, order))
    {
        printf("not zero");
        return 0;
    }
    int64_t lead_pos[rdim];
    leading_positions(lead_pos, base, shifts, ROW_WISE);
    printf("\nleading positions\n");
    for (slong i = 0; i < rdim; i++)
        printf("%lu ", lead_pos[i]);
    return 1;
}



void nmod_poly_mat_shift(nmod_poly_mat_t res, slong k)
{
    nmod_poly_struct *P;
    if (k > 0)
    {
        for (slong i = 0; i < res->r; i++)
            for (slong j = 0; j < res->c; j++)
            {
                P = nmod_poly_mat_entry(res, i, j);
                if (!nmod_poly_is_zero(P))
                    nmod_poly_shift_left(P, P, k);
            }
        return;
    }


    if (k < 0)
    {
        for (slong i = 0; i < res->r; i++)
            for (slong j = 0; j < res->c; j++)
            {
                P = nmod_poly_mat_entry(res, i, j);
                if (!nmod_poly_is_zero(P))
                    nmod_poly_shift_right(P, P, k);
            }
        return;
    }
}

void int64_print_sage(const int64_t *shifts, slong length)
{
    printf("[");
    for (slong i = 0; i < length - 1; i++)
        printf("%ld,", shifts[i]);
    printf("%ld]\n", shifts[length - 1]);
}

void nmod_poly_mat_print_sage(const nmod_poly_mat_t mat)
{
    slong rdim = mat->r, cdim = mat->c;
    nmod_poly_struct *P;
    slong length;

    printf("[");
    for (slong i = 0; i < rdim; i++)
    {
        printf("[");
        for (slong j = 0; j < cdim; j++)
        {
            P = nmod_poly_mat_entry(mat, i, j);
            length = nmod_poly_length(P);
            if (length == 0)
            {
                if (j != cdim - 1)
                    printf("0, ");
                else
                    printf("0");
            }
            else
            {
                for (slong k = 0; k < length; k++)
                {
                    if (k != length - 1)
                    {
                        if (k == 0)
                            printf("%ld +", nmod_poly_get_coeff_ui(P, k));
                        else
                            printf("%ld*x^%ld + ", nmod_poly_get_coeff_ui(P, k), k);
                    }
                    else
                    {
                        if (j != cdim - 1)
                        {

                            if (k == 0)
                                printf("%ld,", nmod_poly_get_coeff_ui(P, k));
                            else
                                printf("%ld*x^%ld,", nmod_poly_get_coeff_ui(P, k), k);

                        }
                        else
                        {

                            if (k == 0)
                                printf("%ld", nmod_poly_get_coeff_ui(P, k));
                            else
                                printf("%ld*x^%ld", nmod_poly_get_coeff_ui(P, k), k);
                        }
                    }
                }
            }
        }
        if (i != rdim -1)
            printf("],\n");
        else
            printf("]");
    }
    printf("]\n");
}

void nmod_mat_print_sage(const nmod_mat_t mat)
{
    slong rdim = mat->r, cdim = mat->c;
    printf("[");
    for (slong i = 0; i < rdim; i++)
    {
        printf("[");
        for (slong j = 0; j < cdim; j++)
        {
            if (j != cdim - 1)
                printf("%ld, ",  nmod_mat_get_entry(mat, i, j));
            else
                printf("%ld",  nmod_mat_get_entry(mat, i, j));

        }
        if (i != rdim -1)
            printf("],\n");
        else
            printf("]");
    }
    printf("]\n");
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
