/*
    Copyright (C) 2025 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include "sagemath_extra.h"

void slongvec_print_sagemath(const slong * shift, slong length)
{
    printf("[");
    for (slong i = 0; i < length - 1; i++)
        printf("%ld,", shift[i]);
    printf("%ld]\n", shift[length - 1]);
}

void nmod_mat_print_sagemath(const nmod_mat_t mat)
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

void nmod_poly_mat_print_sagemath(const nmod_poly_mat_t mat, const char * var)
{
    slong rdim = mat->r, cdim = mat->c;

    printf("[");
    for (slong i = 0; i < rdim; i++)
    {
        printf("[");
        for (slong j = 0; j < cdim; j++)
        {
            nmod_poly_print_sagemath(nmod_poly_mat_entry(mat, i, j), var);
            if (j != cdim - 1)
                printf(", ");
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
