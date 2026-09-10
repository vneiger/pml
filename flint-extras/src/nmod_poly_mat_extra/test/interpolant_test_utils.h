/*
    Copyright (C) 2026 Gilles Villard

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

/* Shared test-only helpers, factored out of
   t-pmintbasis.c/t-pmintbasis_geometric.c/t-pmintbasis_geometric_auto.c.
   Same mechanism as testing_collection.h already which relies on for
   gen_shift/gen_E being shared the same way across these same files. Kept
   as its own small header rather than folded into testing_collection.h
   itself.  */

#ifndef INTERPOLANT_TEST_UTILS_H
#define INTERPOLANT_TEST_UTILS_H

#include <flint/nmod_mat.h>
#include <flint/nmod_poly_mat.h>

/* Converts a filled nmod_poly_mat_t (built via gen_E / nmod_poly_mat_rand,
   which operate on that entrywise representation) into the plain array
   nmod_poly_mat_pmintbasis/mintbasis/pmintbasis_geometric{,_auto} 
   expect. */
static nmod_mat_struct * poly_mat_to_array(const nmod_poly_mat_t E, slong d)
{
    nmod_mat_struct * arr = (nmod_mat_struct *) flint_malloc(FLINT_MAX(d, 1) * sizeof(nmod_mat_struct));
    for (slong k = 0; k < d; k++)
    {
        nmod_mat_init(arr + k, E->r, E->c, E->modulus);
        nmod_poly_mat_get_coeff_mat(arr + k, E, k);
    }
    return arr;
}

static void free_array(nmod_mat_struct * arr, slong d)
{
    for (slong k = 0; k < d; k++)
        nmod_mat_clear(arr + k);
    flint_free(arr);
}

#endif // INTERPOLANT_TEST_UTILS_H
