/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/nmod_vec.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>

#include "nmod_poly_mat_utils.h"

/* The loops go along the rows of the matrices, and the compression allocates
 * each output entry once, to its final length. */

/* g <- (f div x^s) mod x^len, g not aliasing f */
static void
_expand_piece(nmod_poly_struct * g, const nmod_poly_struct * f,
              slong s, slong len)
{
    len = FLINT_MIN(len, f->length - s);
    if (len <= 0)
    {
        nmod_poly_zero(g);
        return;
    }
    nmod_poly_fit_length(g, len);
    _nmod_vec_set(g->coeffs, f->coeffs + s, len);
    _nmod_poly_set_length(g, len);
    _nmod_poly_normalise(g);
}

/* h[s, s + length(p)) += p, where h holds *filled coefficients and is
 * allocated for at least s + length(p): the coefficients from *filled are
 * copied, the zero gap below s is filled, the overlap with h is added */
static void
_compress_add_piece(nn_ptr h, slong * filled, const nmod_poly_struct * p,
                    slong s, nmod_t mod)
{
    const slong l = p->length;
    slong ov;

    if (l == 0)
        return;
    if (s > *filled)
    {
        _nmod_vec_zero(h + *filled, s - *filled);
        *filled = s;
    }
    ov = FLINT_MIN(*filled - s, l);
    _nmod_vec_add(h + s, h + s, p->coeffs, ov, mod);
    _nmod_vec_set(h + s + ov, p->coeffs + ov, l - ov);
    *filled = FLINT_MAX(*filled, s + l);
}

void
nmod_poly_mat_expand_columns(nmod_poly_mat_t epmat,
                             const nmod_poly_mat_t pmat,
                             const slong * orders,
                             slong d)
{
    for (slong i = 0; i < pmat->r; i++)
    {
        slong off = 0;
        for (slong b = 0; b < pmat->c; b++)
        {
            const nmod_poly_struct * f = nmod_poly_mat_entry(pmat, i, b);
            for (slong s = 0; s < orders[b]; s += d, off++)
                _expand_piece(nmod_poly_mat_entry(epmat, i, off), f, s,
                              FLINT_MIN(d, orders[b] - s));
        }
    }
}

void
nmod_poly_mat_compress_columns(nmod_poly_mat_t pmat,
                               const nmod_poly_mat_t epmat,
                               const slong * orders,
                               slong d)
{
    for (slong i = 0; i < pmat->r; i++)
    {
        slong off = 0;
        for (slong b = 0; b < pmat->c; b++)
        {
            nmod_poly_struct * h = nmod_poly_mat_entry(pmat, i, b);
            const slong e = (orders[b] + d - 1) / d;
            slong len = 0, filled = 0;

            for (slong k = 0; k < e; k++)
            {
                const slong l = nmod_poly_mat_entry(epmat, i, off + k)->length;
                if (l > 0)
                    len = FLINT_MAX(len, k * d + l);
            }

            nmod_poly_fit_length(h, len);
            for (slong k = 0; k < e; k++)
                _compress_add_piece(h->coeffs, &filled,
                                    nmod_poly_mat_entry(epmat, i, off + k),
                                    k * d, h->mod);
            _nmod_poly_set_length(h, len);
            _nmod_poly_normalise(h);
            off += e;
        }
    }
}

void
nmod_poly_mat_expand_rows(nmod_poly_mat_t epmat,
                          const nmod_poly_mat_t pmat,
                          const slong * orders,
                          slong d)
{
    slong off = 0;

    for (slong a = 0; a < pmat->r; a++)
        for (slong s = 0; s < orders[a]; s += d, off++)
            for (slong j = 0; j < pmat->c; j++)
                _expand_piece(nmod_poly_mat_entry(epmat, off, j),
                              nmod_poly_mat_entry(pmat, a, j), s,
                              FLINT_MIN(d, orders[a] - s));
}

/* row a of pmat from rows off, ..., off + e_a - 1 of epmat, one row of epmat
 * at a time; filled[j] = number of coefficients of pmat[a, j] written */
void
nmod_poly_mat_compress_rows(nmod_poly_mat_t pmat,
                            const nmod_poly_mat_t epmat,
                            const slong * orders,
                            slong d)
{
    const slong n = pmat->c;
    slong * filled;
    slong off = 0;

    if (n == 0)
        return;

    filled = (slong *) flint_malloc(n * sizeof(slong));

    for (slong a = 0; a < pmat->r; a++)
    {
        const slong e = (orders[a] + d - 1) / d;

        /* final lengths, then one allocation per entry */
        for (slong j = 0; j < n; j++)
            filled[j] = 0;
        for (slong k = 0; k < e; k++)
            for (slong j = 0; j < n; j++)
            {
                const slong l = nmod_poly_mat_entry(epmat, off + k, j)->length;
                if (l > 0)
                    filled[j] = FLINT_MAX(filled[j], k * d + l);
            }
        for (slong j = 0; j < n; j++)
        {
            nmod_poly_struct * h = nmod_poly_mat_entry(pmat, a, j);
            nmod_poly_fit_length(h, filled[j]);
            _nmod_poly_set_length(h, filled[j]);
            filled[j] = 0;
        }

        for (slong k = 0; k < e; k++)
            for (slong j = 0; j < n; j++)
            {
                nmod_poly_struct * h = nmod_poly_mat_entry(pmat, a, j);
                _compress_add_piece(h->coeffs, filled + j,
                                    nmod_poly_mat_entry(epmat, off + k, j),
                                    k * d, h->mod);
            }

        for (slong j = 0; j < n; j++)
            _nmod_poly_normalise(nmod_poly_mat_entry(pmat, a, j));
        off += e;
    }

    flint_free(filled);
}
