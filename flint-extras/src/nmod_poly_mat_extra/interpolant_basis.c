/*
    Copyright (C) 2026 Gilles Villard

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/


#include <flint/nmod_vec.h>

#include "nmod_poly_mat_multiply.h"
#include "nmod_poly_mat_interpolant.h"
#include "nmod_mat_poly.h"
#include "nmod_poly_mat_utils.h"
#include "impl.h"

/** nmod_mat_poly_mintbasis takes E as a plain, flat array of at 
 *     least `d` already-initialized nmod_mat_t's (nmod_mat_struct *)*/ 
void nmod_poly_mat_mintbasis(nmod_poly_mat_t intbas,
                             slong * shift,
                             const ulong * pts,
                             const nmod_mat_struct * E,
                             slong d)
{
    const slong m = intbas->r;
    nmod_mat_poly_t intbasmp;
    nmod_mat_poly_init(intbasmp, m, m, intbas->modulus);
    nmod_mat_poly_mintbasis(intbasmp, shift, pts, E, d);
    nmod_poly_mat_set_from_mat_poly(intbas, intbasmp);
    nmod_mat_poly_clear(intbasmp);
}



void _nmod_poly_mat_interpolant_basis_krylov_leaf(nmod_poly_mat_t P,
                                                  slong * shift,
                                                  const nmod_mat_t V,
                                                  nn_srcptr pts)
{
    nmod_poly_mat_interpolant_basis_krylov(P, V, pts, shift);

    /* P is shift-Popov: the pivot of row i is P[i,i] */
    for (slong i = 0; i < P->r; i++)
        shift[i] += nmod_poly_degree(nmod_poly_mat_entry(P, i, i));
}

void _nmod_poly_mat_evaluate_mul_columns(nmod_mat_t R, const nmod_poly_mat_t P,
                                         const nmod_mat_t V, nn_srcptr pts)
{
    const slong len = nmod_poly_mat_degree(P) + 1, c = V->c;
    const nmod_t mod = V->mod;
    nmod_mat_poly_t Pc;
    nmod_mat_t Pk, T;

    nmod_mat_zero(R);
    if (len <= 0 || c == 0)
        return;

    nmod_mat_poly_init(Pc, P->r, P->c, mod.n);
    nmod_mat_poly_set_from_poly_mat(Pc, P);
    nmod_mat_init(T, R->r, c, mod.n);
    for (slong k = len - 1; k >= 0; k--)
    {
        nmod_mat_poly_coeff_attach(Pk, Pc, k);
        nmod_mat_mul(T, Pk, V);
        for (slong i = 0; i < R->r; i++)
        {
            nn_ptr r = nmod_mat_entry_ptr(R, i, 0);
            nn_srcptr t = nmod_mat_entry_ptr(T, i, 0);
            for (slong b = 0; b < c; b++)
                r[b] = nmod_add(nmod_mul(r[b], pts[b], mod), t[b], mod);
        }
    }
    nmod_mat_clear(T);
    nmod_mat_poly_clear(Pc);
}

void _nmod_poly_mat_pmintbasis_krylov_leaf(nmod_poly_mat_t P,
                                           slong * shift,
                                           const ulong * pts,
                                           const nmod_mat_struct * E,
                                           slong d)
{
    const slong m = P->r;
    const slong n = (d > 0) ? E[0].c : 0;
    nmod_mat_t V;
    nn_ptr vpts;

    /* V = [E_0 | ... | E_{d-1}], column k n + j at the point pts[k] */
    nmod_mat_init(V, m, n * d, P->modulus);
    vpts = _nmod_vec_init(n * d + 1);
    for (slong k = 0; k < d; k++)
    {
        for (slong i = 0; i < m; i++)
            _nmod_vec_set(nmod_mat_entry_ptr(V, i, k * n),
                          nmod_mat_entry_ptr(E + k, i, 0), n);
        for (slong j = 0; j < n; j++)
            vpts[k * n + j] = pts[k];
    }

    _nmod_poly_mat_interpolant_basis_krylov_leaf(P, shift, V, vpts);

    _nmod_vec_clear(vpts);
    nmod_mat_clear(V);
}

void nmod_poly_mat_pmintbasis(nmod_poly_mat_t intbas,
                              slong * shift,
                              const ulong * pts,
                              const nmod_mat_struct * E,
                              slong d)
{
    _nmod_poly_mat_pmintbasis(intbas, shift, pts, E, d,
                              NMOD_POLY_MAT_PMBASIS_KRYLOV_CUTOFF);
}

void _nmod_poly_mat_pmintbasis(nmod_poly_mat_t intbas,
                               slong * shift,
                               const ulong * pts,
                               const nmod_mat_struct * E,
                               slong d,
                               double krylov_cutoff)
{
    /* Krylov leaf: replaces the subtree below this node */
    if (d > 0 && _nmod_poly_mat_krylov_wanted(intbas->r, E[0].c, E[0].c * d, krylov_cutoff))
    {
        _nmod_poly_mat_pmintbasis_krylov_leaf(intbas, shift, pts, E, d);
        return;
    }

    if (d <= PMINTBASIS_THRES)
    {
        nmod_poly_mat_mintbasis(intbas, shift, pts, E, d);
        return;
    }

    const slong m = intbas->r;

    const slong n = E[0].c;
    const slong d1 = (d + 1) / 2;
    const slong d2 = d - d1;

    /* No truncated copy of E needed for the first half -- E may be longer
       than d1, extra entries simply unused. */
    nmod_poly_mat_t P1;
    nmod_poly_mat_init(P1, m, m, intbas->modulus);
    _nmod_poly_mat_pmintbasis(P1, shift, pts, E, d1, krylov_cutoff);

    /* R_i = P1(pts[d1+i]) * E_{d1+i}, i = 0..d2-1 -- unlike pmbasis's
       residual (a middle product), these are d2 evaluate-then-multiply
       operations. With few columns (n^2 <= m and 8 n <= m), they are
       computed from the coefficients of P1, deg(P1) + 1 products by
       [E_{d1} | ... | E_{d-1}] (_nmod_poly_mat_evaluate_mul_columns),
       otherwise with one evaluation of P1 per point. Sandbox (30-bit
       prime, 128 points, P1 of degree 1 to 16): coefficients / evaluations
       0.04 to 0.4 for n = 1 and m = 16 to 256, about 1 at n^2 = m for
       m >= 64; for m = 16 and n = 4, up to 1.6 (pmintbasis 0.7x before
       8 n <= m). */
    nmod_mat_struct * R = (nmod_mat_struct *) flint_malloc(FLINT_MAX(d2, 1) * sizeof(nmod_mat_struct));
    nmod_mat_t Rcat;
    const int coeffs = (n * n <= m && 8 * n <= m);
    if (coeffs)
    {
        nmod_mat_t Ecat, W;
        nn_ptr cpts = _nmod_vec_init(FLINT_MAX(n * d2, 1));
        nmod_mat_init(Ecat, m, n * d2, intbas->modulus);
        nmod_mat_init(Rcat, m, n * d2, intbas->modulus);
        for (slong idx = 0; idx < d2; idx++)
        {
            nmod_mat_window_init(W, Ecat, 0, idx * n, m, (idx + 1) * n);
            nmod_mat_set(W, E + (d1 + idx));
            nmod_mat_window_clear(W);
            for (slong j = 0; j < n; j++)
                cpts[idx * n + j] = pts[d1 + idx];
        }
        _nmod_poly_mat_evaluate_mul_columns(Rcat, P1, Ecat, cpts);
        for (slong idx = 0; idx < d2; idx++)
            nmod_mat_window_init(R + idx, Rcat, 0, idx * n, m, (idx + 1) * n);
        nmod_mat_clear(Ecat);
        _nmod_vec_clear(cpts);
    }
    else
    {
        nmod_mat_t evalP1;
        nmod_mat_init(evalP1, m, m, intbas->modulus);
        for (slong idx = 0; idx < d2; idx++)
        {
            nmod_poly_mat_evaluate_nmod(evalP1, P1, pts[d1 + idx]);
            nmod_mat_init(R + idx, m, n, intbas->modulus);
            nmod_mat_mul(R + idx, evalP1, E + (d1 + idx));
        }
        nmod_mat_clear(evalP1);
    }

    nmod_poly_mat_t P2;
    nmod_poly_mat_init(P2, m, m, intbas->modulus);
    _nmod_poly_mat_pmintbasis(P2, shift, pts + d1, R, d2, krylov_cutoff);
    for (slong idx = 0; idx < d2; idx++)
    {
        if (coeffs)
            nmod_mat_window_clear(R + idx);
        else
            nmod_mat_clear(R + idx);
    }
    if (coeffs)
        nmod_mat_clear(Rcat);
    flint_free(R);

    nmod_poly_mat_multiply(intbas, P2, P1);

    nmod_poly_mat_clear(P1);
    nmod_poly_mat_clear(P2);
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

