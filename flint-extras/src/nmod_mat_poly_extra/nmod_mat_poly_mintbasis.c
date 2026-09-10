/*
    Copyright (C) 2026 Gilles Villard

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_mat.h>
#include <flint/perm.h>
#include "nmod_mat_poly.h"
#include "nmod_mat_extra.h"

/* pair used for a stable sort of row indices by shift value. */
typedef struct
{
    slong value;
    slong index;
} slong_pair;

static int _slong_pair_compare(const void * a, const void * b)
{
    slong_pair aa = *(const slong_pair *) a;
    slong_pair bb = *(const slong_pair *) b;
    if (aa.value == bb.value)
        return (aa.index < bb.index) ? -1 : (aa.index > bb.index ? 1 : 0);
    return (aa.value < bb.value) ? -1 : 1;
}

/* perm := stable permutation making shift[perm[0]] <= ... <= shift[perm[m-1]] */
static void _shift_sort_permutation(slong * perm,
                                    const slong * shift,
                                    slong m,
                                    slong_pair * tmp)
{
    for (slong i = 0; i < m; i++)
    {
        tmp[i].value = shift[i];
        tmp[i].index = i;
    }
    qsort(tmp, m, sizeof(slong_pair), _slong_pair_compare);
    for (slong i = 0; i < m; i++)
        perm[i] = tmp[i].index;
}

/* mintbasis, low-rank fold helper: bottom += nsbas * top, done as explicit
 * rank-1 updates instead of a generic nmod_mat_mul into a temporary buffer
 * followed by nmod_mat_add. Here
 * `rank = top->r = m-nullity`, which is generically `n` (the original 
 * problem's own column count, since nullity=m-n generically) -- so this
 * helps exactly when `n` is small relative to `m`.
 *
 *  THRES=4 is a crossover found by preliminary measurements
 *  TODO investigate this.  */
#ifndef NMOD_MAT_POLY_MINTBASIS_LOWRANK_THRES
#define NMOD_MAT_POLY_MINTBASIS_LOWRANK_THRES 4
#endif

static void _mintbasis_low_rank_addmul(nmod_mat_t bottom, const nmod_mat_t nsbas, const nmod_mat_t top)
{
    slong rank = top->r;
    slong nullity = bottom->r;
    slong c = bottom->c;
    for (slong k = 0; k < rank; k++)
        for (slong i = 0; i < nullity; i++)
        {
            ulong coeff = nmod_mat_entry(nsbas, i, k);
            if (coeff)
                _nmod_vec_scalar_addmul_nmod(nmod_mat_entry_ptr(bottom, i, 0),
                                             nmod_mat_entry_ptr(top, k, 0),
                                             c, coeff, bottom->mod);
        }
}

/* mintbasis, "rescomp" variant: at each iteration, recompute the residual
 * from scratch via a full Horner evaluation of the current `intbas` at
 * `pts[k]`, then multiply by E_k. See
 * nmod_mat_poly_mintbasis_resupdate below for the residual-maintenance
 * variant, and nmod_mat_poly_mintbasis (bottom of this file) for the
 * dispatcher between them. */
void nmod_mat_poly_mintbasis_rescomp(nmod_mat_poly_t intbas,
                                     slong * shift,
                                     const ulong * pts,
                                     const nmod_mat_struct * E,
                                     slong d)
{
    const slong m = intbas->r;

    nmod_mat_poly_one(intbas);

    if (d == 0)
        return;

    const slong n = E[0].c;

    nmod_mat_t evalP;
    nmod_mat_t res;
    nmod_mat_init(evalP, m, m, intbas->mod.n);
    nmod_mat_init(res, m, n, intbas->mod.n);

    slong * perm = _perm_init(m);
    slong_pair * tmp = (slong_pair *) flint_malloc(m * sizeof(slong_pair));
    slong * pivots = (slong *) flint_malloc(m * sizeof(slong));
    nmod_mat_t nsbas;
    /* trivial valid init before the loop, so the unconditional
       nmod_mat_clear at the end is always well-defined,
       d == 0 case, and  nullity==m branch. TO SEE.*/
    nmod_mat_init(nsbas, 0, 0, intbas->mod.n);

    for (slong k = 0; k < d; k++)
    {
        nmod_mat_poly_evaluate_nmod(evalP, intbas, pts[k]);
        nmod_mat_mul(res, evalP, E + k);

        _shift_sort_permutation(perm, shift, m, tmp);
        nmod_mat_permute_rows(res, perm, NULL);

        nmod_mat_clear(nsbas);
        slong nullity = nmod_mat_left_nullspace_compact(nsbas, pivots, res);

        if (nullity == m)
            /* res already zero (in the permuted world): intbas already
               satisfies this point's condition, nothing to update */
            continue;

        /* pivots[i] currently indexes rows of the permuted intbas;
           translate to original row indices */
        _perm_compose(pivots, perm, pivots, m);

        /* the first m-nullity (translated) indices are the rows that will
           be shifted by (X - pts[k]) below: bump their shifted degree now,
           before `pivots` is permuted again */
        for (slong i = 0; i < m - nullity; i++)
            shift[pivots[i]] += 1;

        nmod_mat_poly_permute_rows(intbas, pivots, NULL);

        /* fold: last `nullity` rows <- nsbas * (first m-nullity rows) + themselves.
           Dispatches to explicit rank-1 updates when the inner dimension
           (m-nullity, generically = n) is small -- see
           _mintbasis_low_rank_addmul's own doc comments */ 
        if (m - nullity <= NMOD_MAT_POLY_MINTBASIS_LOWRANK_THRES)
        {
            nmod_mat_t win_mul, win_add;
            for (slong deg = 0; deg < intbas->length; deg++)
            {
                nmod_mat_window_init(win_mul, intbas->coeffs + deg, 0, 0, m - nullity, m);
                nmod_mat_window_init(win_add, intbas->coeffs + deg, m - nullity, 0, m, m);
                _mintbasis_low_rank_addmul(win_add, nsbas, win_mul);
                nmod_mat_window_clear(win_mul);
                nmod_mat_window_clear(win_add);
            }
        }
        else
        {
            nmod_mat_t ns_app, win_mul, win_add;
            nmod_mat_init(ns_app, nullity, intbas->c, intbas->mod.n);
            for (slong deg = 0; deg < intbas->length; deg++)
            {
                nmod_mat_window_init(win_mul, intbas->coeffs + deg, 0, 0, m - nullity, m);
                nmod_mat_window_init(win_add, intbas->coeffs + deg, m - nullity, 0, m, m);
                nmod_mat_mul(ns_app, nsbas, win_mul);
                nmod_mat_add(win_add, win_add, ns_app);
                nmod_mat_window_clear(win_mul);
                nmod_mat_window_clear(win_add);
            }
            nmod_mat_clear(ns_app);
        }

       /* Multiply the first m-nullity rows by (X - pts[k]) (not plain X).
           Recentering Q via y = x - pts[k] (Q_x(x) = Q_y(x - pts[k])):
           a row with Q_y = X*e_i becomes Q_x = (X -
           pts[k])*e_i, i.e. new_row = (X-pts[k])*old_row = X*old_row -
           pts[k]*old_row, giving new_c[deg] = old_c[deg-1] -
           pts[k]*old_c[deg]. Grow length first if any of those rows
           currently has a nonzero top coefficient. */
        for (slong i = 0; i < m - nullity; i++)
            if (!_nmod_vec_is_zero(nmod_mat_poly_entry_ptr(intbas, intbas->length - 1, i, 0), m))
            {
                nmod_mat_poly_fit_length(intbas, intbas->length + 1);
                _nmod_mat_poly_set_length(intbas, intbas->length + 1);
                break;
            }
        /** Recentering uses _nmod_vec_scalar_mul_nmod (dispatches to a 
         * Shoup-precomputed multiply internally) + _nmod_vec_add instead 
         * of a manual per-entry nmod_mul/nmod_sub loop. */
        {
            ulong neg_pt = nmod_neg(pts[k], intbas->mod);
            for (slong deg = intbas->length - 1; deg > 0; deg--)
                for (slong i = 0; i < m - nullity; i++)
                {
                    ulong * dst = nmod_mat_poly_entry_ptr(intbas, deg, i, 0);
                    ulong * prev = nmod_mat_poly_entry_ptr(intbas, deg - 1, i, 0);
                    _nmod_vec_scalar_mul_nmod(dst, dst, m, neg_pt, intbas->mod);
                    _nmod_vec_add(dst, dst, prev, m, intbas->mod);
                }
            for (slong i = 0; i < m - nullity; i++)
            {
                ulong * dst = nmod_mat_poly_entry_ptr(intbas, 0, i, 0);
                _nmod_vec_scalar_mul_nmod(dst, dst, m, neg_pt, intbas->mod);
            }
        }

        _perm_inv(pivots, pivots, m);
        nmod_mat_poly_permute_rows(intbas, pivots, NULL);
    }

    nmod_mat_clear(evalP);
    nmod_mat_clear(res);
    _perm_clear(perm);
    flint_free(tmp);
    flint_free(pivots);
    nmod_mat_clear(nsbas);
}

/* mintbasis, "resupdate" variant: instead of recomputing the residual
 * `intbas(pts[k])*E_k` from scratch via a full Horner evaluation at every
 * iteration, maintain the vector of future residual values
 *   Res[j] = intbas(pts[j]) * E_j,   j = k..d-1
 * updating it via the same elementary row operations applied to `intbas`
 * at each iteration, instead of recomputing it.
 *
 * Validity: intbas_new(x) = M_k(x) * intbas_old(x), where M_k is the
 * (permuted-frame) elementary factor [[X-pts[k], 0], [nsbas, Id]] applied
 * this iteration. For any point y, intbas_new(y) = M_k(y) * intbas_old(y),
 * so Res[j]_new = intbas_new(pts[j])*E_j = M_k(pts[j]) * intbas_old(pts[j])
 * * E_j = M_k(pts[j]) * Res[j]_old -- i.e. Res[j] can be updated directly
 * by applying M_k AT pts[j], without ever re-evaluating intbas.
 *
 * Applying M_k(pts[j]) to Res[j] decomposes into exactly the same two
 * pieces as the polynomial-level update on `intbas` itself:
 *   - fold: bottom `nullity` rows of Res[j] += nsbas * top (m-nullity)
 *     rows of Res[j] -- identical fold applied to intbas's own
 *     coefficients, just on an m x n object (E's width) instead of m x m
 *     (intbas's own width) -- this narrower width is exactly where the
 *     saving comes from, mirroring mbasis_resupdate's own n-vs-m
 *     advantage.
 *   - "X-pts[k]" recentering: for intbas itself (a genuine polynomial)
 *     this needs the shift+rescale combine used in the rescomp version
 *     above. For Res[j] (a value at the fixed point pts[j], not a
 *     polynomial), the analogous transformation is simply a scalar
 *     multiply by (pts[j]-pts[k]) on the top m-nullity rows -- no
 *     coefficient-array bookkeeping needed at all, since evaluating
 *     "(X-pts[k])*old_row" AT pts[j] is just old_row(pts[j])*(pts[j]-
 *     pts[k]).
 *
 *  mintbasis's "X-pts[k]" recentering step on every live
 * Res[j] is a genuine scalar multiply on every entry of the top m-nullity
 * rows, for every remaining j, at every iteration -- an O(d^2*m*n)-ish
 * cost. That extra cost means this
 * variant is not generically cheaper the way mbasis's resupdate is.
 *
 * Rescomp pays one extra real matrix product
 * (`nmod_mat_mul(res, evalP, E->coeffs+k)`) every iteration that resupdate
 * never pays (it reads the already-maintained Res[k] directly); this
 * saving dominates when `d` is small relative to `rdim` and/or `nullity =
 * rdim-cdim` (generically `rdim-cdim`) is small -- the closer to square
 * (cdim near rdim) and the fewer points relative to rdim, the more this
 * saving outweighs resupdate's own extra recentering cost. The win is
 * strongest and scales with `rdim`.  See
 * nmod_mat_poly_mintbasis below for the dispatch condition measurements
 * led to. */
void nmod_mat_poly_mintbasis_resupdate(nmod_mat_poly_t intbas,
                                       slong * shift,
                                       const ulong * pts,
                                       const nmod_mat_struct * E,
                                       slong d)
{
    const slong m = intbas->r;

    nmod_mat_poly_one(intbas);

    if (d == 0)
        return;

    /* n dropped as a parameter -- see nmod_mat_poly_mintbasis_rescomp's
       own comment above for the rationale (identical here). */
    const slong n = E[0].c;

    /* Res[k] = intbas(pts[k])*E_k; initially intbas=Id so Res[k] = E_k. */
    nmod_mat_struct * Res = (nmod_mat_struct *) flint_malloc(d * sizeof(nmod_mat_struct));
    for (slong k = 0; k < d; k++)
    {
        nmod_mat_init(Res + k, m, n, intbas->mod.n);
        nmod_mat_set(Res + k, E + k);
    }

    nmod_mat_t res;
    nmod_mat_init(res, m, n, intbas->mod.n);

    slong * perm = _perm_init(m);
    slong_pair * tmp = (slong_pair *) flint_malloc(m * sizeof(slong_pair));
    slong * pivots = (slong *) flint_malloc(m * sizeof(slong));
    nmod_mat_t nsbas;
    nmod_mat_init(nsbas, 0, 0, intbas->mod.n);

    for (slong k = 0; k < d; k++)
    {
        nmod_mat_set(res, Res + k);
        _shift_sort_permutation(perm, shift, m, tmp);
        nmod_mat_permute_rows(res, perm, NULL);

        nmod_mat_clear(nsbas);
        slong nullity = nmod_mat_left_nullspace_compact(nsbas, pivots, res);

        if (nullity == m)
            continue;

        _perm_compose(pivots, perm, pivots, m);
        for (slong i = 0; i < m - nullity; i++)
            shift[pivots[i]] += 1;

        nmod_mat_poly_permute_rows(intbas, pivots, NULL);
        for (slong j = k + 1; j < d; j++)
            nmod_mat_permute_rows(Res + j, pivots, NULL);

        /* fold: bottom nullity rows += nsbas * top (m-nullity) rows --
           applied to intbas (all coefficient degrees) and to the LIVE
           residuals Res[k+1..d-1] (Res[k] is dead after this iteration).
           Dispatches to explicit rank-1 updates below the same threshold
           as the rescomp variant above. */
        nmod_mat_t win_mul, win_add;
        if (m - nullity <= NMOD_MAT_POLY_MINTBASIS_LOWRANK_THRES)
        {
            for (slong deg = 0; deg < intbas->length; deg++)
            {
                nmod_mat_window_init(win_mul, intbas->coeffs + deg, 0, 0, m - nullity, m);
                nmod_mat_window_init(win_add, intbas->coeffs + deg, m - nullity, 0, m, m);
                _mintbasis_low_rank_addmul(win_add, nsbas, win_mul);
                nmod_mat_window_clear(win_mul);
                nmod_mat_window_clear(win_add);
            }
            for (slong j = k + 1; j < d; j++)
            {
                nmod_mat_window_init(win_mul, Res + j, 0, 0, m - nullity, n);
                nmod_mat_window_init(win_add, Res + j, m - nullity, 0, m, n);
                _mintbasis_low_rank_addmul(win_add, nsbas, win_mul);
                nmod_mat_window_clear(win_mul);
                nmod_mat_window_clear(win_add);
            }
        }
        else
        {
            nmod_mat_t ns_app;
            nmod_mat_init(ns_app, nullity, m, intbas->mod.n);
            for (slong deg = 0; deg < intbas->length; deg++)
            {
                nmod_mat_window_init(win_mul, intbas->coeffs + deg, 0, 0, m - nullity, m);
                nmod_mat_window_init(win_add, intbas->coeffs + deg, m - nullity, 0, m, m);
                nmod_mat_mul(ns_app, nsbas, win_mul);
                nmod_mat_add(win_add, win_add, ns_app);
                nmod_mat_window_clear(win_mul);
                nmod_mat_window_clear(win_add);
            }
            nmod_mat_clear(ns_app);

            nmod_mat_t ns_res;
            nmod_mat_init(ns_res, nullity, n, intbas->mod.n);
            for (slong j = k + 1; j < d; j++)
            {
                nmod_mat_window_init(win_mul, Res + j, 0, 0, m - nullity, n);
                nmod_mat_window_init(win_add, Res + j, m - nullity, 0, m, n);
                nmod_mat_mul(ns_res, nsbas, win_mul);
                nmod_mat_add(win_add, win_add, ns_res);
                nmod_mat_window_clear(win_mul);
                nmod_mat_window_clear(win_add);
            }
            nmod_mat_clear(ns_res);
        }

        /* X-pts[k] recentering of intbas's own top m-nullity rows --
           identical to the rescomp version (genuine polynomial
           coefficients); see its own doc comment for the vectorized
           _nmod_vec_scalar_mul_nmod + _nmod_vec_add rewrite.  */
        for (slong i = 0; i < m - nullity; i++)
            if (!_nmod_vec_is_zero(nmod_mat_poly_entry_ptr(intbas, intbas->length - 1, i, 0), m))
            {
                nmod_mat_poly_fit_length(intbas, intbas->length + 1);
                _nmod_mat_poly_set_length(intbas, intbas->length + 1);
                break;
            }
        {
            ulong neg_pt = nmod_neg(pts[k], intbas->mod);
            for (slong deg = intbas->length - 1; deg > 0; deg--)
                for (slong i = 0; i < m - nullity; i++)
                {
                    ulong * dst = nmod_mat_poly_entry_ptr(intbas, deg, i, 0);
                    ulong * prev = nmod_mat_poly_entry_ptr(intbas, deg - 1, i, 0);
                    _nmod_vec_scalar_mul_nmod(dst, dst, m, neg_pt, intbas->mod);
                    _nmod_vec_add(dst, dst, prev, m, intbas->mod);
                }
            for (slong i = 0; i < m - nullity; i++)
            {
                ulong * dst = nmod_mat_poly_entry_ptr(intbas, 0, i, 0);
                _nmod_vec_scalar_mul_nmod(dst, dst, m, neg_pt, intbas->mod);
            }
        }

        /* top m-nullity rows of the LIVE residuals: scalar multiply by
           (pts[j]-pts[k]) -- see derivation above */
        for (slong j = k + 1; j < d; j++)
        {
            ulong factor = nmod_sub(pts[j], pts[k], intbas->mod);
            for (slong i = 0; i < m - nullity; i++)
                for (slong col = 0; col < n; col++)
                    nmod_mat_entry(Res + j, i, col) =
                        nmod_mul(factor, nmod_mat_entry(Res + j, i, col), intbas->mod);
        }

        _perm_inv(pivots, pivots, m);
        nmod_mat_poly_permute_rows(intbas, pivots, NULL);
        for (slong j = k + 1; j < d; j++)
            nmod_mat_permute_rows(Res + j, pivots, NULL);
    }

    for (slong k = 0; k < d; k++)
        nmod_mat_clear(Res + k);
    flint_free(Res);
    nmod_mat_clear(res);
    _perm_clear(perm);
    flint_free(tmp);
    flint_free(pivots);
    nmod_mat_clear(nsbas);
}

/** Main `mintbasis` function: chooses between
 * @ref nmod_mat_poly_mintbasis_rescomp and
 * @ref nmod_mat_poly_mintbasis_resupdate depending on the shape of `E`
 * and the number of points `d` s.t. `d <= (actual length of pts)` and 
 * `d <= (actual length of E)`.
 * 
 * Both variants always agree bit-for-bit, so this dispatch is a pure
 * timing decision, never a correctness one.
 *
 * The condition here, `d*(m-n+1) <= m`, was found by measurement
 * over a grid of `(m,n,d)`.
 * 
 * \todo investigate for a better dispatcher.  */
void nmod_mat_poly_mintbasis(nmod_mat_poly_t intbas,
                             slong * shift,
                             const ulong * pts,
                             const nmod_mat_struct * E,
                             slong d)
{
    /* d=0: both variants return the identity with no other side effect
       (checked first thing in each, before reading E), so
       there is nothing to dispatch on and E is never dereferenced --
       matches the established "E is never dereferenced when d=0"
       convention (see e.g. interpolant_basis_geometric_auto.c). */
    if (d == 0)
    {
        nmod_mat_poly_one(intbas);
        return;
    }

    const slong m = intbas->r;
    const slong n = E[0].c;

    if (d * (m - n + 1) <= m)
        nmod_mat_poly_mintbasis_resupdate(intbas, shift, pts, E, d);
    else
        nmod_mat_poly_mintbasis_rescomp(intbas, shift, pts, E, d);
}


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
