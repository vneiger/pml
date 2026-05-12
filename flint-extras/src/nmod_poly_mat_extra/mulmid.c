/*
    Copyright (C) 2025 Vincent Neiger, Éric Schost

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/mpn_extras.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>

#include "nmod_extra.h"
#include "nmod_poly_mat_multiply.h"

/** sets C to the first nhi - nlo middle coefficients of the product of A of
 *  length len1 and B of length len2 starting at offset nlo
 */
void nmod_poly_mat_mulmid_naive(nmod_poly_mat_t C, const nmod_poly_mat_t A, const nmod_poly_mat_t B,
                                slong nlo, slong nhi)
{
    /* TODO should have branches depending on how deg(A) and deg(B) compare, etc. */
    /* below more or less assumes deg(A) <= deg(B) */
    slong lenA = nmod_poly_mat_max_length(A);
    if (lenA == 0)
    {
        nmod_poly_mat_zero(C);
        return;
    }
    if (lenA <= nlo)
    {
        nmod_poly_mat_t B_tmp;
        nmod_poly_mat_init(B_tmp, A->c, B->c, B->modulus);
        nmod_poly_mat_shift_right(B_tmp, B, nlo - lenA + 1);
        nmod_poly_mat_mul(C, A, B_tmp);
        nmod_poly_mat_clear(B_tmp);
        /* TODO add function for combined shift+truncate? */
        nmod_poly_mat_shift_right(C, C, lenA - 1);
        nmod_poly_mat_truncate(C, nhi - nlo);
        return;
    }

    nmod_poly_mat_mul(C, A, B);
    /* TODO add function for combined shift+truncate? */
    nmod_poly_mat_shift_right(C, C, nlo);
    nmod_poly_mat_truncate(C, nhi - nlo);
}


/** Middle product for polynomial matrices
 *  sets C = ((A * B) div x^dA) mod x^(dB+1)
 *  output can alias input
 *  ASSUME: len(A) <= nlo and len(B) <= nhi ( TODO really ? )
 *  ASSUME: existence of primitive root ( TODO replace by check!)
 *  uses geometric evaluation and interpolation
 */
void nmod_poly_mat_mulmid_geometric(nmod_poly_mat_t C, const nmod_poly_mat_t A, const nmod_poly_mat_t B,
                                    slong nlo, slong nhi)
{
    const slong m = A->r;
    const slong k = A->c;
    const slong n = B->c;

    if (m < 1 || n < 1 || k < 1)
    {
        nmod_poly_mat_zero(C);
        return;
    }

    if (C == A || C == B)
    {
        nmod_poly_mat_t T;
        nmod_poly_mat_init(T, m, n, C->modulus);
        nmod_poly_mat_mulmid_geometric(T, A, B, nlo, nhi);
        nmod_poly_mat_swap_entrywise(C, T);
        nmod_poly_mat_clear(T);
        return;
    }

    if (nmod_poly_mat_max_length(A) == 0 || nmod_poly_mat_max_length(B) == 0)
    {
        nmod_poly_mat_zero(C);
        return;
    }

    nmod_mat_t *mod_A, *mod_B, *mod_C;
    nmod_t mod;
    nmod_geometric_progression_t F;

    const ulong modn = A->modulus;
    nmod_init(&mod, modn);

    slong ellC = nhi;
    nmod_init(&mod, modn);
    ulong w = nmod_find_root(2*ellC, mod);  /* TODO check result for found! */
    nmod_geometric_progression_init(F, w, ellC, mod);

    mod_A = FLINT_ARRAY_ALLOC(ellC, nmod_mat_t);
    mod_B = FLINT_ARRAY_ALLOC(ellC, nmod_mat_t);
    mod_C = FLINT_ARRAY_ALLOC(ellC, nmod_mat_t);
    nn_ptr val = FLINT_ARRAY_ALLOC(ellC, ulong);
    nn_ptr poly = FLINT_ARRAY_ALLOC(ellC, ulong);

    for (slong i = 0; i < ellC; i++)
    {
        nmod_mat_init(mod_A[i], m, k, modn);
        nmod_mat_init(mod_B[i], k, n, modn);
        nmod_mat_init(mod_C[i], m, n, modn);
    }

    for (slong i = 0; i < m; i++)
    {
        for (slong j = 0; j < k; j++)
        {
            _nmod_poly_reverse(poly,
                               nmod_poly_mat_entry(A, i, j)->coeffs, 
                               nmod_poly_mat_entry(A, i, j)->length, 
                               nlo+1);

            _nmod_poly_evaluate_geometric_nmod_vec_fast_precomp(val, poly, nlo+1, F, ellC, F->mod);
            for (slong ell = 0; ell < ellC; ell++)
                nmod_mat_entry(mod_A[ell], i, j) = val[ell];
        }
    }

    for (slong i = 0; i < k; i++)
    {
        for (slong j = 0; j < n; j++)
        {
            slong len = nmod_poly_mat_entry(B, i, j)->length;
            _nmod_vec_set(val, nmod_poly_mat_entry(B, i, j)->coeffs, len);
            _nmod_vec_zero(val + len, ellC - len);

            _nmod_poly_interpolate_geometric_nmod_vec_fast_precomp(poly, val, F, ellC, mod);
            for (slong ell = 0; ell < ellC; ell++)
                nmod_mat_entry(mod_B[ell], i, j) = poly[ell];
        }
    }

    for (slong ell = 0; ell < ellC; ell++)
        nmod_mat_mul(mod_C[ell], mod_A[ell], mod_B[ell]);

    for (slong i = 0; i < m; i++)
    {
        for (slong j = 0; j < n; j++)
        {
            for (slong ell = 0; ell < ellC; ell++)
                poly[ell] = nmod_mat_entry(mod_C[ell], i, j);

            _nmod_poly_evaluate_geometric_nmod_vec_fast_precomp(val, poly, ellC, F, ellC, F->mod);

            nmod_poly_realloc(nmod_poly_mat_entry(C, i, j), nhi - nlo);
            nmod_poly_mat_entry(C, i, j)->length = nhi - nlo;
            nn_ptr dest = nmod_poly_mat_entry(C, i, j)->coeffs;
            for (slong u = 0; u < nhi - nlo; u++)
                dest[u] = val[u];
            _nmod_poly_normalise(nmod_poly_mat_entry(C, i, j));
        }
    }

    for (slong i = 0; i < ellC; i++)
    {
        nmod_mat_clear(mod_A[i]);
        nmod_mat_clear(mod_B[i]);
        nmod_mat_clear(mod_C[i]);
    }

    flint_free(mod_A);
    flint_free(mod_B);
    flint_free(mod_C);
    flint_free(poly);
    flint_free(val);
    nmod_geometric_progression_clear(F);
}
