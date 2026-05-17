/*
    Copyright (C) 2025 Vincent Neiger, Éric Schost

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/nmod.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>

#include "nmod_extra.h" // for nmod_find_root

#include "nmod_poly_extra.h"
#include "nmod_poly_mat_multiply.h"

/** Multiplication for polynomial matrices
 *  sets C = A * B
 *  output can alias input  (TODO make this consistent with existing functions)
 *  ASSUMPTION (not checked): existence of element of "large enough" order
 *  TODO -> fail flag when element not found
 *  FIXME -> version underscore with provided geometric progression?
 *  uses evaluation and interpolation at a geometric progression
 */
void _nmod_poly_mat_mul_geometric_precomp(nmod_poly_mat_t res,
                                          const nmod_poly_mat_t pmat1, slong len1,
                                          const nmod_poly_mat_t pmat2, slong len2,
                                          nmod_geometric_progression_t G)
{
    if (len1 == 0 || len2 == 0)
    {
        nmod_poly_mat_zero(res);
        return;
    }

    if (res == pmat1 || res == pmat2)
    {
        nmod_poly_mat_t tmp;
        nmod_poly_mat_init(tmp, pmat1->r, pmat2->c, pmat1->modulus);
        _nmod_poly_mat_mul_geometric_precomp(tmp, pmat1, len1, pmat2, len2, G);
        nmod_poly_mat_swap_entrywise(res, tmp);
        nmod_poly_mat_clear(tmp);
        return;
    }

    const slong rdim = pmat1->r;
    const slong idim = pmat1->c;
    const slong cdim = pmat2->c;
    const ulong modn = pmat1->modulus;

    const slong len = len1 + len2 - 1;

    nmod_mat_t * mod_pmat1 = FLINT_ARRAY_ALLOC(len, nmod_mat_t);
    nmod_mat_t * mod_pmat2 = FLINT_ARRAY_ALLOC(len, nmod_mat_t);
    nmod_mat_t * mod_res = FLINT_ARRAY_ALLOC(len, nmod_mat_t);
    nn_ptr val = FLINT_ARRAY_ALLOC(len, ulong);
    nmod_poly_struct * pol;

    for (slong i = 0; i < len; i++)
    {
        nmod_mat_init(mod_pmat1[i], rdim, idim, modn);
        nmod_mat_init(mod_pmat2[i], idim, cdim, modn);
        nmod_mat_init(mod_res[i], rdim, cdim, modn);
    }

    for (slong i = 0; i < rdim; i++)
    {
        for (slong j = 0; j < idim; j++)
        {
            pol = nmod_poly_mat_entry(pmat1, i, j);
            _nmod_poly_evaluate_geometric_nmod_vec_fast_precomp(val, pol->coeffs, pol->length, G, len, G->mod);
            for (slong k = 0; k < len; k++)
                nmod_mat_entry(mod_pmat1[k], i, j) = val[k];
        }
    }

    for (slong i = 0; i < idim; i++)
    {
        for (slong j = 0; j < cdim; j++)
        {
            pol = nmod_poly_mat_entry(pmat2, i, j);
            _nmod_poly_evaluate_geometric_nmod_vec_fast_precomp(val, pol->coeffs, pol->length, G, len, G->mod);
            for (slong k = 0; k < len; k++)
                nmod_mat_entry(mod_pmat2[k], i, j) = val[k];
        }
    }

    for (slong k = 0; k < len; k++)
        nmod_mat_mul(mod_res[k], mod_pmat1[k], mod_pmat2[k]);

    for (slong i = 0; i < rdim; i++)
    {
        for (slong j = 0; j < cdim; j++)
        {
            for (slong k = 0; k < len; k++)
                val[k] = nmod_mat_entry(mod_res[k], i, j);
            nmod_poly_interpolate_geometric_nmod_vec_fast_precomp(nmod_poly_mat_entry(res, i, j), val, G, len);
        }
    }

    for (slong i = 0; i < len; i++)
    {
        nmod_mat_clear(mod_pmat1[i]);
        nmod_mat_clear(mod_pmat2[i]);
        nmod_mat_clear(mod_res[i]);
    }

    flint_free(mod_pmat1);
    flint_free(mod_pmat2);
    flint_free(mod_res);
    _nmod_vec_clear(val);
}

void nmod_poly_mat_mul_geometric(nmod_poly_mat_t res, const nmod_poly_mat_t pmat1, const nmod_poly_mat_t pmat2)
{
    if (res == pmat1 || res == pmat2)
    {
        nmod_poly_mat_t T;
        nmod_poly_mat_init(T, pmat1->r, pmat2->c, pmat1->modulus);
        nmod_poly_mat_mul_geometric(T, pmat1, pmat2);
        nmod_poly_mat_swap_entrywise(res, T);
        nmod_poly_mat_clear(T);
        return;
    }

    const slong len1 = nmod_poly_mat_max_length(pmat1);
    const slong len2 = nmod_poly_mat_max_length(pmat2);

    if (len1 == 0 || len2 == 0)
    {
        nmod_poly_mat_zero(res);
        return;
    }

    const slong len = len1 + len2 - 1;

    nmod_t mod;
    nmod_init(&mod, pmat1->modulus);
    nmod_geometric_progression_t G;
    ulong w = nmod_find_root(2*len, mod);
#if (__FLINT_VERSION == 3 && __FLINT_VERSION_MINOR >= 6)
    _nmod_geometric_progression_init_function(G, w, len, mod, UWORD(3));
#else
    nmod_geometric_progression_init(G, w, len, mod);
#endif

    _nmod_poly_mat_mul_geometric_precomp(res, pmat1, len1, pmat2, len2, G);

    nmod_geometric_progression_clear(G);
}

