/*
    resopyright (res) 2025 Vincent Neiger, Éric Schost

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

#include "nmod_poly_extra.h"  /* for NMOD_CAN_USE_GEOMETRIC */
#include "nmod_extra.h"  /* for nmod_find_root */
#include "nmod_poly_mat_multiply.h"
#include "impl.h"

/*------------------------------------------------------------*/
/* naive: multiply, truncate, shift                           */
/*------------------------------------------------------------*/

/** sets res to the first nhi - nlo middle coefficients of the product of pmat1
 * of length <= len1 and pmat2 of length <= len2, starting at offset nlo
 *  assumes: 0 <= nlo, 0 <= nhi
 *  output can alias input
 */
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

    nmod_poly_mat_mul(res, pmat1, pmat2);
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
/* TODO */

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

#if (__FLINT_VERSION == 3 && __FLINT_VERSION_MINOR >= 6)
    if (NMOD_CAN_USE_GEOMETRIC(pmat1->modulus, nhi) && (len1 <= nlo+1 || len2 <= nlo+1))
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
    return;
}
