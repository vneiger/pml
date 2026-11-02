/*
    Copyright (C) 2025 Gilles Villard
    Copyright (C) 2026 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/flint.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>

#include "nmod_poly_mat_inverse.h"
#include "nmod_poly_mat_multiply.h"
#include "nmod_poly_mat_utils.h"

/*------------------------------------------------------------*/
/* thresholds                                                 */
/*------------------------------------------------------------*/

/* Provisional thresholds, depending on the dimension m (see the documentation
 * in nmod_poly_mat_inverse.h).
 * \todo tuning (also depending on the modulus?) */

slong _nmod_poly_mat_inv_series_cutoff(slong m)
{
    if (m <= 8)
        return 64;
    if (m <= 16)
        return 48;
    if (m <= 32)
        return 32;
    return 24;
}

/* \todo tuning */
slong _nmod_poly_mat_inv_series_short_length(slong m)
{
    return (m <= 32) ? 9 : 5;
}

/*------------------------------------------------------------*/
/* truncated inverse: basecase                                */
/*------------------------------------------------------------*/

/* Writing A = sum_i A_i x^i and X = A^{-1} mod x^order = sum_k X_k x^k, the
 * identity X A = I gives X_0 = A_0^{-1} and, for k >= 1,
 *     X_k = sum_{1 <= i <= min(k,d)} X_{k-i} * Ahat_i,
 * where Ahat_i = - A_i * X_0 and d = min(deg(A), order-1).
 *
 * The X_k are stored in a single m x (order*m) constant matrix
 * W = [X_0 | X_1 | ... | X_{order-1}], and the Ahat_i are stored in
 * reverse order in the (d*m) x m constant matrix Ahat = [Ahat_d; ...; Ahat_1].
 * With k' = min(k,d), the sum above is then the product of two windows:
 * [X_{k-k'} | ... | X_{k-1}] (columns (k-k')*m to k*m-1 of W) times
 * [Ahat_{k'}; ...; Ahat_1] (the last k'*m rows of Ahat). */
void _nmod_poly_mat_inv_series_basecase(nmod_poly_mat_t imat,
                                        const nmod_poly_mat_t pmat,
                                        const nmod_mat_t inv0,
                                        slong order)
{
    const slong m = pmat->r;

    if (order <= 0 || m == 0)
    {
        nmod_poly_mat_zero(imat);
        return;
    }

    /* number of useful nonconstant coefficients of pmat */
    const slong d = FLINT_MIN(nmod_poly_mat_max_length(pmat), order) - 1;

    if (d <= 0)
    {
        /* imat = inv0 (constant) */
        nmod_poly_mat_zero(imat);
        for (slong i = 0; i < m; i++)
            for (slong j = 0; j < m; j++)
                nmod_poly_set_coeff_ui(nmod_poly_mat_entry(imat, i, j), 0,
                                       nmod_mat_entry(inv0, i, j));
        return;
    }

    /* Ahat = - [A_d; ...; A_1] * inv0, of dimensions (d*m) x m */
    nmod_mat_t Ahat, tmp;
    nmod_mat_init(Ahat, d * m, m, pmat->modulus);
    nmod_mat_init(tmp, d * m, m, pmat->modulus);
    for (slong i = 0; i < m; i++)
        for (slong j = 0; j < m; j++)
        {
            const nmod_poly_struct * a = nmod_poly_mat_entry(pmat, i, j);
            const slong len = FLINT_MIN(a->length, d + 1);
            for (slong k = 1; k < len; k++)
                nmod_mat_entry(tmp, (d-k)*m + i, j) = a->coeffs[k];
        }
    nmod_mat_mul(Ahat, tmp, inv0);
    nmod_mat_neg(Ahat, Ahat);
    nmod_mat_clear(tmp);

    /* W = [X_0 | X_1 | ... | X_{order-1}] */
    nmod_mat_t W;
    nmod_mat_init(W, m, order * m, pmat->modulus);
    for (slong i = 0; i < m; i++)
        for (slong j = 0; j < m; j++)
            nmod_mat_entry(W, i, j) = nmod_mat_entry(inv0, i, j);

    nmod_mat_t Xwin, Awin, Ywin;
    for (slong k = 1; k < order; k++)
    {
        const slong kk = FLINT_MIN(k, d);
        nmod_mat_window_init(Xwin, W, 0, (k-kk)*m, m, k*m);
        nmod_mat_window_init(Awin, Ahat, (d-kk)*m, 0, d*m, m);
        nmod_mat_window_init(Ywin, W, 0, k*m, m, (k+1)*m);
        nmod_mat_mul(Ywin, Xwin, Awin);
        nmod_mat_window_clear(Ywin);
        nmod_mat_window_clear(Awin);
        nmod_mat_window_clear(Xwin);
    }

    /* copy W into imat */
    for (slong i = 0; i < m; i++)
        for (slong j = 0; j < m; j++)
        {
            nmod_poly_struct * x = nmod_poly_mat_entry(imat, i, j);
            nmod_poly_fit_length(x, order);
            for (slong k = 0; k < order; k++)
                x->coeffs[k] = nmod_mat_entry(W, i, k*m + j);
            _nmod_poly_set_length(x, order);
            _nmod_poly_normalise(x);
        }

    nmod_mat_clear(W);
    nmod_mat_clear(Ahat);
}

int nmod_poly_mat_inv_series_basecase(nmod_poly_mat_t imat,
                                      const nmod_poly_mat_t pmat,
                                      slong order)
{
    const slong m = pmat->r;

    if (order <= 0 || m == 0)
    {
        nmod_poly_mat_zero(imat);
        return 1;
    }

    nmod_mat_t inv0;
    nmod_mat_init(inv0, m, m, pmat->modulus);
    nmod_poly_mat_get_coeff_mat(inv0, pmat, 0);
    if (! nmod_mat_inv(inv0, inv0))
    {
        nmod_mat_clear(inv0);
        return 0;
    }

    if (imat == pmat)
    {
        nmod_poly_mat_t tmp;
        nmod_poly_mat_init(tmp, m, m, pmat->modulus);
        _nmod_poly_mat_inv_series_basecase(tmp, pmat, inv0, order);
        nmod_poly_mat_swap_entrywise(imat, tmp);
        nmod_poly_mat_clear(tmp);
    }
    else
        _nmod_poly_mat_inv_series_basecase(imat, pmat, inv0, order);

    nmod_mat_clear(inv0);
    return 1;
}

/*------------------------------------------------------------*/
/* truncated inverse: Newton iteration                        */
/*------------------------------------------------------------*/

/* Newton iteration X <- X - X (A X - I), which doubles the precision.
 *
 * Used truncation order is half of the requested one, rounding up
 * (until arriving below the basecase cutoff).
 *
 * Going from precision prec_src to precision prec_tgt <= 2 prec_src, with
 * X = A^{-1} mod x^prec_src of length at most prec_src:
 *  - E = (A X)[prec_src, prec_tgt), the coefficients of degree prec_src to
 *    prec_tgt-1 of A X, so that A X = I + x^prec_src E mod x^prec_tgt (one
 *    middle product);
 *  - F = X E mod x^{prec_tgt - prec_src} (one low product);
 *  - X <- X - x^prec_src F, that is, the coefficients of degree prec_src to
 *    prec_tgt-1 of X are set to those of -F (X has length at most prec_src
 *    before this update).
 */
void _nmod_poly_mat_inv_series_newton(nmod_poly_mat_t imat,
                                      const nmod_poly_mat_t pmat,
                                      const nmod_mat_t inv0,
                                      slong order,
                                      slong cutoff)
{
    const slong m = pmat->r;

    if (order <= 0 || m == 0)
    {
        nmod_poly_mat_zero(imat);
        return;
    }

    if (cutoff < 2)
        cutoff = 2;

    if (order < cutoff)
    {
        _nmod_poly_mat_inv_series_basecase(imat, pmat, inv0, order);
        return;
    }

    /* precisions: steps[0] = order > steps[1] > ... > steps[i] < cutoff */
    slong steps[FLINT_BITS];
    slong i = 0;
    steps[0] = order;
    while (steps[i] >= cutoff)
    {
        steps[i+1] = (steps[i] + 1) / 2;
        i++;
    }

    _nmod_poly_mat_inv_series_basecase(imat, pmat, inv0, steps[i]);

    nmod_poly_mat_t E, F;
    nmod_poly_mat_init(E, m, m, pmat->modulus);
    nmod_poly_mat_init(F, m, m, pmat->modulus);

    for (i--; i >= 0; i--)
    {
        const slong prec_src = steps[i+1];
        const slong prec_tgt = steps[i];

        /* A X = I + x^prec_src E mod x^prec_tgt */
        nmod_poly_mat_mulmid(E, pmat, imat, prec_src, prec_tgt);
        /* F = X E mod x^(prec_tgt - prec_src) */
        nmod_poly_mat_mulmid(F, imat, E, 0, prec_tgt - prec_src);

        /* X = X - x^prec_src F: coefficients prec_src .. prec_tgt-1 of X set
         * to those of -F, taking the fast path without arithmetic operations
         * since X has length <= prec_src */
        nmod_poly_mat_sub_shift_left(imat, imat, F, prec_src);
    }

    nmod_poly_mat_clear(E);
    nmod_poly_mat_clear(F);
}

int nmod_poly_mat_inv_series_newton(nmod_poly_mat_t imat,
                                    const nmod_poly_mat_t pmat,
                                    slong order,
                                    slong cutoff)
{
    const slong m = pmat->r;

    if (order <= 0 || m == 0)
    {
        nmod_poly_mat_zero(imat);
        return 1;
    }

    if (cutoff <= 0)
        cutoff = _nmod_poly_mat_inv_series_cutoff(m);

    nmod_mat_t inv0;
    nmod_mat_init(inv0, m, m, pmat->modulus);
    nmod_poly_mat_get_coeff_mat(inv0, pmat, 0);
    if (! nmod_mat_inv(inv0, inv0))
    {
        nmod_mat_clear(inv0);
        return 0;
    }

    if (imat == pmat)
    {
        nmod_poly_mat_t tmp;
        nmod_poly_mat_init(tmp, m, m, pmat->modulus);
        _nmod_poly_mat_inv_series_newton(tmp, pmat, inv0, order, cutoff);
        nmod_poly_mat_swap_entrywise(imat, tmp);
        nmod_poly_mat_clear(tmp);
    }
    else
        _nmod_poly_mat_inv_series_newton(imat, pmat, inv0, order, cutoff);

    nmod_mat_clear(inv0);
    return 1;
}

/*------------------------------------------------------------*/
/* truncated inverse: general interface                       */
/*------------------------------------------------------------*/

void _nmod_poly_mat_inv_series(nmod_poly_mat_t imat,
                               const nmod_poly_mat_t pmat,
                               const nmod_mat_t inv0,
                               slong order)
{
    const slong cutoff = _nmod_poly_mat_inv_series_cutoff(pmat->r);

    if (order < cutoff
        || nmod_poly_mat_max_length(pmat) <= _nmod_poly_mat_inv_series_short_length(pmat->r))
        _nmod_poly_mat_inv_series_basecase(imat, pmat, inv0, order);
    else
        _nmod_poly_mat_inv_series_newton(imat, pmat, inv0, order, cutoff);
}

int nmod_poly_mat_inv_series(nmod_poly_mat_t imat,
                             const nmod_poly_mat_t pmat,
                             slong order)
{
    const slong m = pmat->r;

    if (order <= 0 || m == 0)
    {
        nmod_poly_mat_zero(imat);
        return 1;
    }

    nmod_mat_t inv0;
    nmod_mat_init(inv0, m, m, pmat->modulus);
    nmod_poly_mat_get_coeff_mat(inv0, pmat, 0);
    if (! nmod_mat_inv(inv0, inv0))
    {
        nmod_mat_clear(inv0);
        return 0;
    }

    if (imat == pmat)
    {
        nmod_poly_mat_t tmp;
        nmod_poly_mat_init(tmp, m, m, pmat->modulus);
        _nmod_poly_mat_inv_series(tmp, pmat, inv0, order);
        nmod_poly_mat_swap_entrywise(imat, tmp);
        nmod_poly_mat_clear(tmp);
    }
    else
        _nmod_poly_mat_inv_series(imat, pmat, inv0, order);

    nmod_mat_clear(inv0);
    return 1;
}
