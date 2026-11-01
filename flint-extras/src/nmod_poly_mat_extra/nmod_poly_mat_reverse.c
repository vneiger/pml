/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>

#include "nmod_poly_mat_utils.h"

/* rpol <- reverse of pol with respect to the bound hi, truncated at order
 * `order` (no truncation if order < 0): coefficient k of rpol is coefficient
 * hi - k of pol, for 0 <= k <= hi and k < order.
 * With len = min(hi+1, order), this is the plain reverse, as a polynomial of
 * length len, of the slice pol[hi-len+1 .. hi] (zero beyond pol->length). rpol
 * must not alias pol. */
static void _reverse_trunc(nmod_poly_t rpol, const nmod_poly_t pol, slong hi, slong order)
{
    slong len = hi + 1;  /* length of the reverse, before truncation */
    if (order >= 0 && order < len)
        len = order;
    const slong start = hi - len + 1;          /* first coefficient of the slice */
    const slong avail = pol->length - start;   /* its nonzero part has length <= avail */
    if (hi < 0 || len <= 0 || avail <= 0)
    {
        nmod_poly_zero(rpol);
        return;
    }
    nmod_poly_fit_length(rpol, len);
    _nmod_poly_reverse(rpol->coeffs, pol->coeffs + start, FLINT_MIN(avail, len), len);
    _nmod_poly_set_length(rpol, len);
    _nmod_poly_normalise(rpol);
}

void nmod_poly_mat_reverse_trunc(nmod_poly_mat_t rmat, const nmod_poly_mat_t pmat,
                                 const slong * rhi, const slong * chi, slong order)
{
    nmod_poly_t tmp;
    nmod_poly_init(tmp, pmat->modulus);
    for (slong i = 0; i < pmat->r; i++)
        for (slong j = 0; j < pmat->c; j++)
        {
            const slong hi = (rhi ? rhi[i] : 0) + (chi ? chi[j] : 0);
            /* through tmp, so that rmat may alias pmat */
            _reverse_trunc(tmp, nmod_poly_mat_entry(pmat, i, j), hi, order);
            nmod_poly_swap(nmod_poly_mat_entry(rmat, i, j), tmp);
        }
    nmod_poly_clear(tmp);
}

void nmod_poly_mat_reverse(nmod_poly_mat_t rmat, const nmod_poly_mat_t pmat, slong hi)
{
    nmod_poly_t tmp;
    nmod_poly_init(tmp, pmat->modulus);
    for (slong i = 0; i < pmat->r; i++)
        for (slong j = 0; j < pmat->c; j++)
        {
            _reverse_trunc(tmp, nmod_poly_mat_entry(pmat, i, j), hi, -1);
            nmod_poly_swap(nmod_poly_mat_entry(rmat, i, j), tmp);
        }
    nmod_poly_clear(tmp);
}

void nmod_poly_mat_reverse_rowwise(nmod_poly_mat_t rmat, const nmod_poly_mat_t pmat, const slong * hi)
{
    nmod_poly_mat_reverse_trunc(rmat, pmat, hi, NULL, -1);
}

void nmod_poly_mat_reverse_columnwise(nmod_poly_mat_t rmat, const nmod_poly_mat_t pmat, const slong * hi)
{
    nmod_poly_mat_reverse_trunc(rmat, pmat, NULL, hi, -1);
}
