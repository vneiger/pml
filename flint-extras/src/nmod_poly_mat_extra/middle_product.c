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
                                const slong nlo, const slong nhi)
{
    nmod_poly_mat_mul(C, A, B);
    nmod_poly_mat_truncate(C, nhi);
    nmod_poly_mat_shift_right(C, C, nlo);
}


/** Middle product for polynomial matrices
 *  sets C = ((A * B) div x^dA) mod x^(dB+1)
 *  output can alias input
 *  ASSUME: deg(A) <= dA and deg(B) <= dA + dB
 *  ASSUME: existence of primitive root
 *  uses geometric evaluation and interpolation
 */
void nmod_poly_mat_middle_product_geometric(nmod_poly_mat_t C, const nmod_poly_mat_t A, const nmod_poly_mat_t B,
                                            const ulong dA, const ulong dB)
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
        nmod_poly_mat_middle_product_geometric(T, A, B, dA, dB);
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

    slong ellC = dA + dB + 1;  // length(C) = length(A) + length(B) - 1
    nmod_init(&mod, modn);
    ulong w = nmod_find_root(2*ellC, mod);  /* TODO check necessary order */
    nmod_geometric_progression_init(F, w, ellC, mod);

    mod_A = FLINT_ARRAY_ALLOC(ellC, nmod_mat_t);
    mod_B = FLINT_ARRAY_ALLOC(ellC, nmod_mat_t);
    mod_C = FLINT_ARRAY_ALLOC(ellC, nmod_mat_t);
    nn_ptr val = FLINT_ARRAY_ALLOC(ellC, ulong);
    nn_ptr tmp_poly = FLINT_ARRAY_ALLOC(ellC, ulong);

#ifdef DIRTY_ALLOC_MATRIX
    // we alloc the memory for all matrices at once
    nn_ptr tmp = (nn_ptr) flint_malloc((m*k + k*n + m*n) * ellC * sizeof(ulong));
    nn_ptr bak;

    bak = tmp;
    for (slong i = 0; i < ellC; i++)
    {
        mod_A[i]->entries = tmp + i*m*k;
        mod_A[i]->stride = k;
        mod_A[i]->r = m;
        mod_A[i]->c = k;
        mod_A[i]->mod.n = mod.n;
        mod_A[i]->mod.norm = mod.norm;
        mod_A[i]->mod.ninv = mod.ninv;
    }
    tmp += ellC*m*k;

    for (slong i = 0; i < ellC; i++)
    {
        mod_B[i]->entries = tmp + i*k*n;
        mod_B[i]->stride = n;
        mod_B[i]->r = k;
        mod_B[i]->c = n;
        mod_B[i]->mod.n = mod.n;
        mod_B[i]->mod.norm = mod.norm;
        mod_B[i]->mod.ninv = mod.ninv;
    }
    tmp += ellC*k*n;

    for (slong i = 0; i < ellC; i++)
    {
        mod_C[i]->entries = tmp + i*m*n;
        mod_C[i]->stride = n;
        mod_C[i]->r = m;
        mod_C[i]->c = n;
        mod_C[i]->mod.n = mod.n;
        mod_C[i]->mod.norm = mod.norm;
        mod_C[i]->mod.ninv = mod.ninv;
    }
    tmp = bak;
#else
    for (slong i = 0; i < ellC; i++)
    {
        nmod_mat_init(mod_A[i], m, k, modn);
        nmod_mat_init(mod_B[i], k, n, modn);
        nmod_mat_init(mod_C[i], m, n, modn);
    }
#endif

    for (slong i = 0; i < m; i++)
    {
        for (slong j = 0; j < k; j++)
        {
            _nmod_poly_reverse(tmp_poly,
                               nmod_poly_mat_entry(A, i, j)->coeffs, 
                               nmod_poly_mat_entry(A, i, j)->length, 
                               dA+1);

            _nmod_poly_evaluate_geometric_nmod_vec_fast_precomp(val, tmp_poly, dA+1, F, ellC, F->mod);
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

            _nmod_poly_interpolate_geometric_nmod_vec_fast_precomp(tmp_poly, val, F, ellC, mod);
            for (slong ell = 0; ell < ellC; ell++)
                nmod_mat_entry(mod_B[ell], i, j) = tmp_poly[ell];
        }
    }

    for (slong ell = 0; ell < ellC; ell++)
        nmod_mat_mul(mod_C[ell], mod_A[ell], mod_B[ell]);

    for (slong i = 0; i < m; i++)
    {
        for (slong j = 0; j < n; j++)
        {
            for (slong ell = 0; ell < ellC; ell++)
                tmp_poly[ell] = nmod_mat_entry(mod_C[ell], i, j);

            _nmod_poly_evaluate_geometric_nmod_vec_fast_precomp(val, tmp_poly, ellC, F, ellC, F->mod);

            nmod_poly_realloc(nmod_poly_mat_entry(C, i, j), dB + 1);
            nmod_poly_mat_entry(C, i, j)->length = dB + 1;
            nn_ptr dest = nmod_poly_mat_entry(C, i, j)->coeffs;
            for (slong u = 0; u < (slong)dB+1; u++)
                dest[u] = val[u];
            _nmod_poly_normalise(nmod_poly_mat_entry(C, i, j));
        }
    }

#ifdef DIRTY_ALLOC_MATRIX
    flint_free(tmp);
#else
    for (slong i = 0; i < ellC; i++)
    {
        nmod_mat_clear(mod_A[i]);
        nmod_mat_clear(mod_B[i]);
        nmod_mat_clear(mod_C[i]);
    }
#endif

    flint_free(mod_A);
    flint_free(mod_B);
    flint_free(mod_C);
    flint_free(tmp_poly);
    flint_free(val);
    nmod_geometric_progression_clear(F);
}
