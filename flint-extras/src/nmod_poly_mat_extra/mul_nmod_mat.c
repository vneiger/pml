/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/nmod_mat.h>
#include <flint/nmod_poly_mat.h>

#include "nmod_poly_mat_arith.h"
#include "nmod_poly_mat_forms.h"   /* row and column degrees */
#include "nmod_poly_mat_utils.h"   /* linearization */

void
nmod_poly_mat_mul_left_nmod_mat(nmod_poly_mat_t c,
                                const nmod_mat_t a,
                                const nmod_poly_mat_t b)
{
    const slong n = b->c;
    slong * orders = (slong *) flint_malloc(FLINT_MAX(n, 1) * sizeof(slong));
    slong D = 0;
    nmod_mat_t bl, cl;

    /* orders: the lengths of the columns of b */
    nmod_poly_mat_column_degree(orders, b, NULL);
    for (slong j = 0; j < n; j++)
    {
        orders[j] += 1;
        D += orders[j];
    }

    nmod_mat_init(bl, b->r, D, a->mod.n);
    nmod_mat_init(cl, a->r, D, a->mod.n);

    nmod_poly_mat_linearize_columns(bl, b, orders);
    nmod_mat_mul(cl, a, bl);
    nmod_poly_mat_polynomialize_columns(c, cl, orders);

    nmod_mat_clear(cl);
    nmod_mat_clear(bl);
    flint_free(orders);
}

void
nmod_poly_mat_mul_right_nmod_mat(nmod_poly_mat_t c,
                                 const nmod_poly_mat_t a,
                                 const nmod_mat_t b)
{
    const slong m = a->r;
    slong * orders = (slong *) flint_malloc(FLINT_MAX(m, 1) * sizeof(slong));
    slong D = 0;
    nmod_mat_t al, cl;

    /* orders: the lengths of the rows of a */
    nmod_poly_mat_row_degree(orders, a, NULL);
    for (slong i = 0; i < m; i++)
    {
        orders[i] += 1;
        D += orders[i];
    }

    nmod_mat_init(al, D, a->c, b->mod.n);
    nmod_mat_init(cl, D, b->c, b->mod.n);

    nmod_poly_mat_linearize_rows(al, a, orders);
    nmod_mat_mul(cl, al, b);
    nmod_poly_mat_polynomialize_rows(c, cl, orders);

    nmod_mat_clear(cl);
    nmod_mat_clear(al);
    flint_free(orders);
}
