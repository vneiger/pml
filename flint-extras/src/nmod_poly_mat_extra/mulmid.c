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
#include <flint/nmod_vec.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>

#include "nmod_poly_extra.h"  /* for NMOD_POLY_CAN_USE_GEOMETRIC */
#include "nmod_extra.h"  /* for nmod_find_root */
#include "nmod_poly_mat_multiply.h"
#include "impl.h"

/*------------------------------------------------------------*/
/* naive: multiply, truncate, shift                           */
/*------------------------------------------------------------*/

void _nmod_poly_mat_mulmid_naive(nmod_poly_mat_t res,
                                 const nmod_poly_mat_t pmat1, slong len1,
                                 const nmod_poly_mat_t pmat2, slong len2,
                                 slong nlo, slong nhi)
{
    if (len1 == 0 || len2 == 0 || nlo >= nhi)
    {
        nmod_poly_mat_zero(res);
        return;
    }

    nmod_poly_mat_multiply(res, pmat1, pmat2);
    nmod_poly_mat_shift_right(res, res, nlo);
    nmod_poly_mat_truncate(res, nhi - nlo);
}

/*------------------------------------------------------------*/
/* evaluation-interpolation at geometric progression          */
/*------------------------------------------------------------*/

#if (__FLINT_VERSION == 3 && __FLINT_VERSION_MINOR >= 6)
/** Middle product for polynomial matrices, via geometric evaluation-interpolation
 *  This is very close to geometric2: see differences there
 *  requires max_length(pmat1) <= nlo+1
 *  output cannot alias input
 */
void _nmod_poly_mat_mulmid_geometric1_precomp(nmod_poly_mat_t res,
                                              const nmod_poly_mat_t pmat1, const nmod_poly_mat_t pmat2,
                                              slong nlo, slong nhi, nmod_geometric_progression_t G)
{
    const slong rdim = pmat1->r;
    const slong idim = pmat1->c;
    const slong cdim = pmat2->c;
    const ulong modn = pmat1->modulus;

    nmod_mat_t * mod_pmat1 = FLINT_ARRAY_ALLOC(nhi, nmod_mat_t);
    nmod_mat_t * mod_pmat2 = FLINT_ARRAY_ALLOC(nhi, nmod_mat_t);
    nmod_mat_t * mod_res = FLINT_ARRAY_ALLOC(nhi, nmod_mat_t);
    nn_ptr val = FLINT_ARRAY_ALLOC(nhi, ulong);
    nn_ptr poly = FLINT_ARRAY_ALLOC(nhi, ulong);

    for (slong i = 0; i < nhi; i++)
    {
        nmod_mat_init(mod_pmat1[i], rdim, idim, modn);
        nmod_mat_init(mod_pmat2[i], idim, cdim, modn);
        nmod_mat_init(mod_res[i], rdim, cdim, modn);
    }

    for (slong i = 0; i < rdim; i++)
    {
        for (slong j = 0; j < idim; j++)
        {
            _nmod_poly_reverse(poly,
                               nmod_poly_mat_entry(pmat1, i, j)->coeffs,
                               nmod_poly_mat_entry(pmat1, i, j)->length,
                               nlo+1);

            _nmod_poly_evaluate_geometric_nmod_vec_fast_precomp(val, poly, nlo+1, G, nhi, G->mod);
            for (slong k = 0; k < nhi; k++)
                nmod_mat_entry(mod_pmat1[k], i, j) = val[k];
        }
    }

    for (slong i = 0; i < idim; i++)
    {
        for (slong j = 0; j < cdim; j++)
        {
            slong len = FLINT_MIN(nhi, nmod_poly_mat_entry(pmat2, i, j)->length);
            _nmod_vec_set(val, nmod_poly_mat_entry(pmat2, i, j)->coeffs, len);
            _nmod_vec_zero(val + len, nhi - len);

            _nmod_poly_interpolate_geometric_nmod_vec_fast_precomp(poly, val, G, nhi, G->mod);
            for (slong k = 0; k < nhi; k++)
                nmod_mat_entry(mod_pmat2[k], i, j) = poly[k];
        }
    }

    for (slong k = 0; k < nhi; k++)
        nmod_mat_mul(mod_res[k], mod_pmat1[k], mod_pmat2[k]);

    for (slong i = 0; i < rdim; i++)
    {
        for (slong j = 0; j < cdim; j++)
        {
            for (slong k = 0; k < nhi; k++)
                poly[k] = nmod_mat_entry(mod_res[k], i, j);

            _nmod_poly_evaluate_geometric_nmod_vec_fast_precomp(val, poly, nhi, G, nhi - nlo, G->mod);

            nmod_poly_fit_length(nmod_poly_mat_entry(res, i, j), nhi - nlo);
            _nmod_poly_set_length(nmod_poly_mat_entry(res, i, j), nhi - nlo);
            nn_ptr dest = nmod_poly_mat_entry(res, i, j)->coeffs;
            for (slong k = 0; k < nhi - nlo; k++)
                dest[k] = val[k];
            _nmod_poly_normalise(nmod_poly_mat_entry(res, i, j));
        }
    }

    for (slong i = 0; i < nhi; i++)
    {
        nmod_mat_clear(mod_pmat1[i]);
        nmod_mat_clear(mod_pmat2[i]);
        nmod_mat_clear(mod_res[i]);
    }

    flint_free(mod_pmat1);
    flint_free(mod_pmat2);
    flint_free(mod_res);
    flint_free(poly);
    flint_free(val);
}

/** Middle product for polynomial matrices, via geometric evaluation-interpolation
 *  This is very close to geometric1: differences are discussed in comments starting with [+-+-]
 *  requires max_length(pmat2) <= nlo+1   ([+-+-] differs from geometric1)
 *  output cannot alias input
 */
void _nmod_poly_mat_mulmid_geometric2_precomp(nmod_poly_mat_t res,
                                              const nmod_poly_mat_t pmat1, const nmod_poly_mat_t pmat2,
                                              slong nlo, slong nhi, nmod_geometric_progression_t G)
{
    const slong rdim = pmat1->r;
    const slong idim = pmat1->c;
    const slong cdim = pmat2->c;

    nmod_mat_t *mod_pmat1, *mod_pmat2, *mod_res;
    const ulong modn = pmat1->modulus;

    mod_pmat1 = FLINT_ARRAY_ALLOC(nhi, nmod_mat_t);
    mod_pmat2 = FLINT_ARRAY_ALLOC(nhi, nmod_mat_t);
    mod_res = FLINT_ARRAY_ALLOC(nhi, nmod_mat_t);
    nn_ptr val = FLINT_ARRAY_ALLOC(nhi, ulong);
    nn_ptr poly = FLINT_ARRAY_ALLOC(nhi, ulong);

    for (slong i = 0; i < nhi; i++)
    {
        nmod_mat_init(mod_pmat1[i], rdim, idim, modn);
        nmod_mat_init(mod_pmat2[i], idim, cdim, modn);
        nmod_mat_init(mod_res[i], rdim, cdim, modn);
    }

    /* [+-+-] : in this loop, pmat1/pmat2 swapped, dimensions adapted accordingly */
    for (slong i = 0; i < idim; i++)
    {
        for (slong j = 0; j < cdim; j++)
        {
            _nmod_poly_reverse(poly,
                               nmod_poly_mat_entry(pmat2, i, j)->coeffs,
                               nmod_poly_mat_entry(pmat2, i, j)->length,
                               nlo+1);

            _nmod_poly_evaluate_geometric_nmod_vec_fast_precomp(val, poly, nlo+1, G, nhi, G->mod);
            for (slong k = 0; k < nhi; k++)
                nmod_mat_entry(mod_pmat2[k], i, j) = val[k];
        }
    }

    /* [+-+-] : in this loop, pmat1/pmat2 swapped, dimensions adapted accordingly */
    for (slong i = 0; i < rdim; i++)
    {
        for (slong j = 0; j < idim; j++)
        {
            slong len = FLINT_MIN(nhi, nmod_poly_mat_entry(pmat1, i, j)->length);
            _nmod_vec_set(val, nmod_poly_mat_entry(pmat1, i, j)->coeffs, len);
            _nmod_vec_zero(val + len, nhi - len);

            _nmod_poly_interpolate_geometric_nmod_vec_fast_precomp(poly, val, G, nhi, G->mod);
            for (slong k = 0; k < nhi; k++)
                nmod_mat_entry(mod_pmat1[k], i, j) = poly[k];
        }
    }

    for (slong k = 0; k < nhi; k++)
        nmod_mat_mul(mod_res[k], mod_pmat1[k], mod_pmat2[k]);

    for (slong i = 0; i < rdim; i++)
    {
        for (slong j = 0; j < cdim; j++)
        {
            for (slong k = 0; k < nhi; k++)
                poly[k] = nmod_mat_entry(mod_res[k], i, j);

            _nmod_poly_evaluate_geometric_nmod_vec_fast_precomp(val, poly, nhi, G, nhi - nlo, G->mod);

            nmod_poly_fit_length(nmod_poly_mat_entry(res, i, j), nhi - nlo);
            _nmod_poly_set_length(nmod_poly_mat_entry(res, i, j), nhi - nlo);
            nn_ptr dest = nmod_poly_mat_entry(res, i, j)->coeffs;
            for (slong k = 0; k < nhi - nlo; k++)
                dest[k] = val[k];
            _nmod_poly_normalise(nmod_poly_mat_entry(res, i, j));
        }
    }

    for (slong i = 0; i < nhi; i++)
    {
        nmod_mat_clear(mod_pmat1[i]);
        nmod_mat_clear(mod_pmat2[i]);
        nmod_mat_clear(mod_res[i]);
    }

    flint_free(mod_pmat1);
    flint_free(mod_pmat2);
    flint_free(mod_res);
    flint_free(poly);
    flint_free(val);
}

void _nmod_poly_mat_mulmid_geometric_precomp(nmod_poly_mat_t res,
                                             const nmod_poly_mat_t pmat1, slong len1,
                                             const nmod_poly_mat_t pmat2, slong len2,
                                             slong nlo, slong nhi, nmod_geometric_progression_t G)
{
    PML_ASSERT(len1 <= nlo+1 || len2 <= nlo+1)

    if (len1 == 0 || len2 == 0)
    {
        nmod_poly_mat_zero(res);
        return;
    }

    if (res == pmat1 || res == pmat2)
    {
        nmod_poly_mat_t temp;
        nmod_poly_mat_init(temp, pmat1->r, pmat2->c, pmat1->modulus);
        _nmod_poly_mat_mulmid_geometric(temp, pmat1, len1, pmat2, len2, nlo, nhi);
        nmod_poly_mat_swap(res, temp);
        nmod_poly_mat_clear(temp);
        return;
    }

    if (len1 <= nlo+1)
        _nmod_poly_mat_mulmid_geometric1_precomp(res, pmat1, pmat2, nlo, nhi, G);
    else  /* len2 <= nlo+1 */
        _nmod_poly_mat_mulmid_geometric2_precomp(res, pmat1, pmat2, nlo, nhi, G);
}

void _nmod_poly_mat_mulmid_geometric(nmod_poly_mat_t res,
                                     const nmod_poly_mat_t pmat1, slong len1,
                                     const nmod_poly_mat_t pmat2, slong len2,
                                     slong nlo, slong nhi)
{
    PML_ASSERT(len1 <= nlo+1 || len2 <= nlo+1)

    if (len1 == 0 || len2 == 0)
    {
        nmod_poly_mat_zero(res);
        return;
    }

    if (res == pmat1 || res == pmat2)
    {
        nmod_poly_mat_t temp;
        nmod_poly_mat_init(temp, pmat1->r, pmat2->c, pmat1->modulus);
        _nmod_poly_mat_mulmid_geometric(temp, pmat1, len1, pmat2, len2, nlo, nhi);
        nmod_poly_mat_swap(res, temp);
        nmod_poly_mat_clear(temp);
        return;
    }

    nmod_t mod;
    nmod_geometric_progression_t G;

    nmod_init(&mod, pmat1->modulus);
    ulong w = nmod_find_root(2*nhi, mod);
    _nmod_geometric_progression_init_function(G, w, nhi, mod, UWORD(3));

    _nmod_poly_mat_mulmid_geometric_precomp(res, pmat1, len1, pmat2, len2, nlo, nhi, G);

    nmod_geometric_progression_clear(G);
}
#endif

/*------------------------------------------------------------*/
/* evaluation-interpolation using Vandermonde matrix          */
/*------------------------------------------------------------*/

/* TODO currently disabled: untested and unprofiled */
#if 0
/** Middle product for polynomial matrices
 *  sets C = ((A * B) div x^dA) mod x^(dB+1)
 *  output can alias input
 *  uses evaluation and interpolation at arithmetic points, done by matrix products
 *  ASSUME: large enough field 
 */
static void nmod_poly_mat_middle_product_vandermonde1(nmod_poly_mat_t c, const nmod_poly_mat_t a, const nmod_poly_mat_t b,
                                               const ulong dA, const ulong dB)
{

    ulong i, j, k, ellA, ellB, m, n, p, ell, u, v, nb_points;
    nmod_mat_t vA, vB, vB_t, iv, iv_t, tmp_mat, valA, valB, valAp, valBp, valC, valCp;
    nmod_t mod;
    
    ellA = nmod_poly_mat_max_length(a);
    ellB = nmod_poly_mat_max_length(b);

    if (ellA == 0 || ellB == 0)
    {
        nmod_poly_mat_zero(c);
        return;
    }

    m = a->r;
    n = a->c;
    p = b->c;

    if (c == a || c == b)
    {
        nmod_poly_mat_t T;
        nmod_poly_mat_init(T, m, n, p);
        nmod_poly_mat_middle_product_vandermonde1(T, a, b, dA, dB);
        nmod_poly_mat_swap_entrywise(c, T);
        nmod_poly_mat_clear(T);
        return;
    }
    
    nmod_init(&mod, a->modulus);
    vandermonde_init1(vA, vB, iv, dA, dB, mod);
    nb_points = dA + dB + 1;

    // evaluation of matrix a:
    // build tmp_mat, whose column ell = i*n + j is the coefficient vector 
    // of reverse(a[i][j]) 
    nmod_mat_init(tmp_mat, dA + 1, m * n, mod.n);
    ell = 0;
    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++, ell++)
        {
            ulong d;
            nn_ptr ptr;
            d = nmod_poly_degree(nmod_poly_mat_entry(a, i, j));
            ptr = nmod_poly_mat_entry(a, i, j)->coeffs;
            for (k = 0; k <= d; k++)
                nmod_mat_entry(tmp_mat, dA - k, ell) = ptr[k];
        }
    // note: d = deg(a[i][j]) is -1 if a[i][j] == 0
    // all non-touched entries already zero since tmp_mat was initialized as zero
   
    // valA: column ell = i*n + j contains the evaluations of a[i][j]
    nmod_mat_init(valA, vA->r, tmp_mat->c, mod.n);
    nmod_mat_mul_pml(valA, vA, tmp_mat);

    // transpose interpolation of matrix b:
    // build tmp_mat, whose column ell = i*n + j is the coefficient vector of
    // b[i][j] (padded with zeroes up to length dB+1 if necessary)
    nmod_mat_clear(tmp_mat);
    nmod_mat_init(tmp_mat, dA + dB + 1, n * p, mod.n);
    ell = 0;
    for (i = 0; i < n; i++)
        for ( j = 0; j < p; j++, ell++)
        {
            ulong d;
            nn_ptr ptr;
            d = nmod_poly_degree(nmod_poly_mat_entry(b, i, j));
            ptr = nmod_poly_mat_entry(b, i, j)->coeffs;
            for (k = 0; k <= d; k++)
                nmod_mat_entry(tmp_mat, k, ell) = ptr[k];
        }
    // mul by transpose(iv)
    nmod_mat_init(iv_t, iv->r, iv->c, mod.n);
    nmod_mat_transpose(iv_t, iv);
    nmod_mat_init(valB, iv_t->r, tmp_mat->c, mod.n);
    nmod_mat_mul_pml(valB, iv_t, tmp_mat);
 

    // perform the pointwise products
    nmod_mat_init(valAp, m, n, mod.n);
    nmod_mat_init(valBp, n, p, mod.n);
    nmod_mat_init(valC, nb_points, m * p, mod.n);
    nmod_mat_init(valCp, m, p, mod.n);
    
    for (i = 0; i < nb_points; i++)
    {
        // a evaluated at point i
        ell = 0;
        for (u = 0; u < m; u++)
            for (v = 0; v < n; v++, ell++)
                nmod_mat_entry(valAp, u, v) = nmod_mat_entry(valA, i, ell);

        ell = 0;
        for (u = 0; u < n; u++)
            for (v = 0; v < p; v++, ell++)
                nmod_mat_entry(valBp, u, v) = nmod_mat_entry(valB, i, ell);

        nmod_mat_mul_pml(valCp, valAp, valBp);

        // copy this into valC
        ell = 0;
        for (u = 0; u < m; u++)
            for (v = 0; v < p; v++, ell++)
                nmod_mat_entry(valC, i, ell) = nmod_mat_entry(valCp, u, v);
    }
    
    nmod_mat_init(vB_t, vB->c, vB->r, mod.n);
    nmod_mat_transpose(vB_t, vB);
    // transpose-evaluate to find the entries of c
    nmod_mat_clear(tmp_mat);
    nmod_mat_init(tmp_mat, vB_t->r, valC->c, mod.n);
    nmod_mat_mul_pml(tmp_mat, vB_t, valC);


    // copy to output (reorganize these entries into c)
    ell = 0;
    for (u = 0; u < m; u++)
        for (v = 0; v < p; v++, ell++)
        {
            nn_ptr coeffs;
            nmod_poly_realloc(nmod_poly_mat_entry(c, u, v), dB + 1);
            coeffs = nmod_poly_mat_entry(c, u, v)->coeffs;
            nmod_poly_mat_entry(c, u, v)->length = dB + 1; 
            for (i = 0; i <= dB; i++)
                coeffs[i] = nmod_mat_entry(tmp_mat, i, ell);
            _nmod_poly_normalise(nmod_poly_mat_entry(c, u, v));
        }
    
    nmod_mat_clear(vA);
    nmod_mat_clear(vB_t);
    nmod_mat_clear(vB);
    nmod_mat_clear(iv_t);
    nmod_mat_clear(iv);
    nmod_mat_clear(tmp_mat);
    nmod_mat_clear(valA);
    nmod_mat_clear(valB);
    nmod_mat_clear(valAp);
    nmod_mat_clear(valBp);
    nmod_mat_clear(valC);
    nmod_mat_clear(valCp);
}
#endif


/*------------------------------------------------------------*/
/* via 3-prime FFT                                            */
/*------------------------------------------------------------*/
/* TODO */




/*------------------------------------------------------------*/
/* general interfaces                                         */
/*------------------------------------------------------------*/

void _nmod_poly_mat_mulmid(nmod_poly_mat_t res,
                           const nmod_poly_mat_t pmat1, slong len1,
                           const nmod_poly_mat_t pmat2, slong len2,
                           slong nlo, slong nhi)
{
    /* zero matrices or empty target coefficient indices */
    if (len1 == 0 || len2 == 0 || nlo >= nhi)
    {
        nmod_poly_mat_zero(res);
        return;
    }

    /* TODO handle constant pmat1 or pmat2 properly */
    /* if (len1 == 1)*/
    /*     _nmod_poly_mat_mulmid_naive(res, pmat1, len1, pmat2, len2, nlo, nhi); */
    /* else */

    /* TODO rough thresholds, not finely tuned */
#if (__FLINT_VERSION == 3 && __FLINT_VERSION_MINOR >= 6)
    if (NMOD_POLY_CAN_USE_GEOMETRIC(pmat1->modulus, nhi)
        && ((pmat1->r >= 8 && pmat2->c >= 2) || (pmat1->r >= 2 && pmat2->c >= 8))
        && (len1 <= nlo+1 || len2 <= nlo+1))
        _nmod_poly_mat_mulmid_geometric(res, pmat1, len1, pmat2, len2, nlo, nhi);

    else
#endif
        _nmod_poly_mat_mulmid_naive(res, pmat1, len1, pmat2, len2, nlo, nhi);
}

void nmod_poly_mat_mulmid(nmod_poly_mat_t res, const nmod_poly_mat_t pmat1, const nmod_poly_mat_t pmat2,
                          slong nlo, slong nhi)
{
    PML_ASSERT(nlo >= 0);
    PML_ASSERT(nhi >= 0);

    slong len1 = nmod_poly_mat_max_length(pmat1);
    slong len2 = nmod_poly_mat_max_length(pmat2);
    nhi = FLINT_MIN(nhi, len1 + len2 - 1);

    /* zero matrices or empty target coefficient indices */
    if (len1 == 0 || len2 == 0 || nlo >= nhi)
    {
        nmod_poly_mat_zero(res);
        return;
    }

    if (len1 > nhi)
    {
        nmod_poly_mat_t pmat1_tmp;
        nmod_poly_mat_init(pmat1_tmp, pmat1->r, pmat1->c, pmat1->modulus);
        nmod_poly_mat_set_trunc(pmat1_tmp, pmat1, nhi);
        nmod_poly_mat_mulmid(res, pmat1_tmp, pmat2, nlo, nhi);
        nmod_poly_mat_clear(pmat1_tmp);
        return;
    }

    if (len2 > nhi)
    {
        nmod_poly_mat_t pmat2_tmp;
        nmod_poly_mat_init(pmat2_tmp, pmat2->r, pmat2->c, pmat2->modulus);
        nmod_poly_mat_set_trunc(pmat2_tmp, pmat2, nhi);
        nmod_poly_mat_mulmid(res, pmat1, pmat2_tmp, nlo, nhi);
        nmod_poly_mat_clear(pmat2_tmp);
        return;
    }

    /* len1 <= len2 and len(pmat1) <= nlo: shift away useless coefficients of pmat2 */
    if (len1 <= len2 && len1 <= nlo)
    {
        nmod_poly_mat_t pmat2_tmp;
        nmod_poly_mat_init(pmat2_tmp, pmat2->r, pmat2->c, pmat2->modulus);
        nmod_poly_mat_shift_right(pmat2_tmp, pmat2, nlo - len1 + 1);
        _nmod_poly_mat_mulmid(res, pmat1, len1, pmat2_tmp, len2 - nlo + len1 - 1, len1 - 1, nhi - nlo + len1 - 1);
        nmod_poly_mat_clear(pmat2_tmp);
        return;
    }

    /* if len(pmat2) <= nlo (implies len1 > len2), shift away useless coefficients of pmat1 */
    if (len2 <= nlo)
    {
        nmod_poly_mat_t pmat1_tmp;
        nmod_poly_mat_init(pmat1_tmp, pmat1->r, pmat1->c, pmat1->modulus);
        nmod_poly_mat_shift_right(pmat1_tmp, pmat1, nlo - len2 + 1);
        _nmod_poly_mat_mulmid(res, pmat1_tmp, len1 - nlo + len2 - 1, pmat2, len2, len2 - 1, nhi - nlo + len2 - 1);
        nmod_poly_mat_clear(pmat1_tmp);
        return;
    }

    _nmod_poly_mat_mulmid(res, pmat1, len1, pmat2, len2, nlo, nhi);
}
