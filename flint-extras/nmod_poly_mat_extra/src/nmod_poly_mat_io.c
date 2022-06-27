#include "nmod_poly_mat_io.h"

void slong_print_sage(const slong *shifts, slong length)
{
    printf("[");
    for (slong i = 0; i < length - 1; i++)
        printf("%ld,", shifts[i]);
    printf("%ld]\n", shifts[length - 1]);
}

void nmod_poly_mat_print_sage(const nmod_poly_mat_t mat, const char * var)
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
                            printf("%ld*%s**%ld + ", nmod_poly_get_coeff_ui(P, k), var, k);
                    }
                    else
                    {
                        if (j != cdim - 1)
                        {

                            if (k == 0)
                                printf("%ld,", nmod_poly_get_coeff_ui(P, k));
                            else
                                printf("%ld*%s**%ld,", nmod_poly_get_coeff_ui(P, k), var, k);

                        }
                        else
                        {

                            if (k == 0)
                                printf("%ld", nmod_poly_get_coeff_ui(P, k));
                            else
                                printf("%ld*%s**%ld", nmod_poly_get_coeff_ui(P, k), var, k);
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

void slong_mat_print(const slong *mat, slong rdim, slong cdim)
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

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
