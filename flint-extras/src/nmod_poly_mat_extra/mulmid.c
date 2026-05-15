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
 *  assumes: 0 <= nlo, 0 <= nhi
 *  output can alias input
 */
void _nmod_poly_mat_mulmid_naive(nmod_poly_mat_t C,
                                 const nmod_poly_mat_t A, slong lenA,
                                 const nmod_poly_mat_t B, slong lenB,
                                 slong nlo, slong nhi)
{
    if (lenA == 0 || lenB == 0 || nlo >= nhi)
    {
        nmod_poly_mat_zero(C);
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
void _nmod_poly_mat_mulmid_geometric(nmod_poly_mat_t C,
                                    const nmod_poly_mat_t A, slong lenA,
                                    const nmod_poly_mat_t B, slong lenB,
                                    slong nlo, slong nhi)
{
    const slong m = A->r;
    const slong k = A->c;
    const slong n = B->c;

    if (lenA == 0 || lenB == 0)
    {
        nmod_poly_mat_zero(C);
        return;
    }

    if (C == A || C == B)
    {
        nmod_poly_mat_t temp;
        nmod_poly_mat_init(temp, A->r, B->c, A->modulus);
        _nmod_poly_mat_mulmid_geometric(temp, A, lenA, B, lenB, nlo, nhi);
        nmod_poly_mat_swap(C, temp);
        nmod_poly_mat_clear(temp);
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
    _nmod_geometric_progression_init_function(F, w, ellC, mod, UWORD(3));

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

            _nmod_poly_evaluate_geometric_nmod_vec_fast_precomp(val, poly, ellC, F, nhi - nlo, F->mod);

            nmod_poly_fit_length(nmod_poly_mat_entry(C, i, j), nhi - nlo);
            _nmod_poly_set_length(nmod_poly_mat_entry(C, i, j), nhi - nlo);
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

/** general interface
 * TODO description
 * requires: lenA <= lenB (required?)
 * 0 <= nlo < nhi <= lenA + lenB - 1 (required?)
 * lenA = max_length(A), lenB = max_length(B), max(lenA, lenB) <= nhi
 * aliasing is allowed
 */
void _nmod_poly_mat_mulmid(nmod_poly_mat_t C,
                           const nmod_poly_mat_t A, slong lenA,
                           const nmod_poly_mat_t B, slong lenB,
                           slong nlo, slong nhi)
{
    /* zero matrices or empty target coefficient indices */
    if (lenA == 0 || lenB == 0 || nlo >= nhi)
    {
        nmod_poly_mat_zero(C);
        return;
    }

    /* if (lenA == 1)  /1* TODO handle constant A or B properly *1/ */
    /*     _nmod_poly_mat_mulmid_naive(C, A, lenA, B, lenB, nlo, nhi); */
    /* else */
    if (lenA <= nlo+1 && lenB <= nhi && A->modulus > (ulong)10*nhi) /* FIXME temporary condition to ensure existence of geometric progression */
    /* if (lenA <= nlo+1)  /1* TODO version once truncation has been inserted in main function *1/ */
    {
        _nmod_poly_mat_mulmid_geometric(C, A, lenA, B, lenB, nlo, nhi);
    }
    /* else if (lenB <= nlo+1) /1* TODO version not implemented yet *1/ */
    else
        _nmod_poly_mat_mulmid_naive(C, A, lenA, B, lenB, nlo, nhi);
}

/** general interface
 * TODO description
 * aliasing is allowed
 */
void nmod_poly_mat_mulmid(nmod_poly_mat_t C, const nmod_poly_mat_t A, const nmod_poly_mat_t B,
                          slong nlo, slong nhi)
{
    PML_ASSERT(nlo >= 0);
    PML_ASSERT(nhi >= 0);

    slong lenA = nmod_poly_mat_max_length(A);
    slong lenB = nmod_poly_mat_max_length(B);
    nhi = FLINT_MIN(nhi, lenA + lenB - 1);

    /* zero matrices or empty target coefficient indices */
    if (lenA == 0 || lenB == 0 || nlo >= nhi)
    {
        nmod_poly_mat_zero(C);
        return;
    }

    /* TODO truncate A and B to not exceed length nhi */

    /* lenA <= lenB and len(A) <= nlo: shift away useless coefficients of B */
    if (lenA <= lenB && lenA <= nlo)
    {
        nmod_poly_mat_t B_tmp;
        nmod_poly_mat_init(B_tmp, B->r, B->c, B->modulus);
        nmod_poly_mat_shift_right(B_tmp, B, nlo - lenA + 1);
        _nmod_poly_mat_mulmid(C, A, lenA, B_tmp, lenB - nlo + lenA - 1, lenA - 1, nhi - nlo + lenA - 1);
        nmod_poly_mat_clear(B_tmp);
        return;
    }

    /* if len(B) <= nlo (implies lenA > lenB), shift away useless coefficients of A */
    if (lenB <= nlo)
    {
        nmod_poly_mat_t A_tmp;
        nmod_poly_mat_init(A_tmp, A->r, A->c, A->modulus);
        nmod_poly_mat_shift_right(A_tmp, A, nlo - lenB + 1);
        _nmod_poly_mat_mulmid(C, A_tmp, lenA - nlo + lenB - 1, B, lenB, lenB - 1, nhi - nlo + lenB - 1);
        nmod_poly_mat_clear(A_tmp);
        return;
    }

    _nmod_poly_mat_mulmid(C, A, lenA, B, lenB, nlo, nhi);
    return;
}
