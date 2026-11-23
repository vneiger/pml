/*
    Copyright (C) 2026 Gilles Villard

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/nmod.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>
#include "nmod_mat_poly.h"
#include "nmod_poly_mat_interpolant.h"
#include "nmod_mat_poly_extra/impl.h"

/*
    Storage, as in mbasis (nmod_mat_poly_mbasis.c): the rows of the
    interpolant basis are stored as vectors of coefficients, row i of the
    matrix W being [coefficient 0 of row i of intbas | coefficient 1 | ...],
    with m entries per coefficient. W has room for cap coefficients, a bound
    on the length of intbas which is increased when needed (initially, about
    the length for generic input with uniform shift, n*d/m + 1). In the
    "resupdate" variant, the residuals R_j = intbas(pts[j]) E_j at the points
    j not processed yet are stored in the m x (n d) matrix R, whose row i is
    [row i of R_0 | row i of R_1 | ...]. Then, at the iteration for the point
    pts[k]:

    - in the "rescomp" variant, the residual intbas(pts[k]) E_k is computed
      from the evaluation of intbas at pts[k], one vector-matrix product per
      row (the powers of pts[k] times the coefficients of the row), and one
      product by E_k; for n = 1, it is directly one matrix-vector product, of
      the coefficients of intbas by the powers of pts[k] times E_k;

    - the update of the basis (and of the residuals at the points j > k) by
      the left nullspace of the residual is a single product X * (pivot rows),
      over all the coefficients at once, instead of one product per
      coefficient of intbas, plus a matrix addition each time;

    - the pivot rows are multiplied by x - pts[k] from a copy of them, and
      their residuals at the points j > k by pts[j] - pts[k], entrywise;

    - there is no row permutation of intbas (or of the residuals).

    The row operations are the same as those of the previous implementation
    (on nmod_mat_poly_t, with permutations of the rows), in the same order, so
    the output (intbas and shift) is the same, and the two variants give the
    same output.
*/

typedef struct
{
    slong m;            /* rows of the E_k */
    slong n;            /* columns of the E_k */
    slong d;            /* number of points */
    nmod_t mod;
    int resupd;         /* variant: residuals updated, or recomputed */
    nmod_mat_t W;       /* m x (cap m): the rows of intbas */
    slong cap;          /* capacity of W in coefficients, at most d + 1 */
    slong len;          /* length of intbas */
    nmod_mat_t R;       /* (resupd) m x (n d): the residuals at the points */
}
_mintbasis_struct;

typedef _mintbasis_struct _mintbasis_t[1];

/* W starts as the identity, R as [E_0 | ... | E_{d-1}] */
static void
_mintbasis_init(_mintbasis_t S, slong m, slong n, slong d,
                const nmod_mat_struct * E, nmod_t mod, int resupd)
{
    S->m = m;
    S->n = n;
    S->d = d;
    S->mod = mod;
    S->resupd = resupd;
    S->cap = FLINT_MIN(d + 1, (n * d + m - 1) / m + 2);
    nmod_mat_init(S->W, m, S->cap * m, mod.n);
    for (slong i = 0; i < m; i++)
        nmod_mat_entry(S->W, i, i) = 1;
    S->len = 1;
    nmod_mat_init(S->R, resupd ? m : 0, resupd ? n * d : 0, mod.n);
    if (resupd)
        for (slong k = 0; k < d; k++)
            for (slong i = 0; i < m; i++)
                _nmod_vec_set(nmod_mat_entry_ptr(S->R, i, k * n),
                              nmod_mat_entry_ptr(E + k, i, 0), n);
}

static void
_mintbasis_clear(_mintbasis_t S)
{
    nmod_mat_clear(S->W);
    nmod_mat_clear(S->R);
}

/* row i of intbas: coefficient e is at offset e*m */
static inline ulong *
_mintbasis_row(_mintbasis_t S, slong i)
{
    return nmod_mat_entry_ptr(S->W, i, 0);
}

/* The iterations of mintbasis: on input, S holds the identity (and the
   residuals); on output, intbas (and S->len); shift is updated. */
static void
_mintbasis_run(_mintbasis_t S, slong * shift, const ulong * pts,
               const nmod_mat_struct * E, slong lowrank)
{
    const slong m = S->m, n = S->n, d = S->d;
    const int resupd = S->resupd;
    const nmod_t mod = S->mod;
    slong len = S->len;

    /* T: the pivot rows (rank <= min(m, n)), U: the combinations of them
       added to the other rows (allocated when needed: not for rank <=
       lowrank); their columns: the coefficients of intbas, and the residuals
       (resupd); ev: (rescomp, n > 1) intbas(pts[k]) */
    nmod_mat_t res, K, X, Wv, ev, T, U, Tw, Uw;
    nmod_mat_init(res, m, n, mod.n);
    nmod_mat_init(K, m, m, mod.n);
    nmod_mat_init(ev, (resupd || n == 1) ? 0 : m, m, mod.n);
    nmod_mat_init(T, FLINT_MIN(m, n), S->cap * m + S->R->c, mod.n);
    nmod_mat_init(U, 0, 0, mod.n);

    /* perm: rows sorted by shift; rows: the pivot rows of the residual, then
       the others; P, Q: scratch for the PLE decomposition */
    slong * perm = (slong *) flint_malloc((3 * m + n) * sizeof(slong));
    slong * rows = perm + m;
    slong * P = rows + m;
    slong * Q = P + m;
    /* pw: (rescomp) powers of the point, and for n = 1 the column E_k times
       them (q); F: the factors pts[j] - pts[k] of the residuals of the pivot
       rows (resupdate) */
    nn_ptr pw = _nmod_vec_init(d + 1);
    nn_ptr q = (resupd || n > 1) ? NULL : _nmod_vec_init((d + 1) * m);
    nn_ptr v = _nmod_vec_init(m);
    nn_ptr F = _nmod_vec_init(FLINT_MAX(n * d, 1));
    for (slong i = 0; i < m; i++)
        perm[i] = i;

    for (slong k = 0; k < d; k++)
    {
        const ulong pt = pts[k];

        /* Here intbas is an s-ordered weak Popov interpolant basis for the
           points 0..k-1, where s is the input shift, and shift =
           rdeg_s(intbas). The residual: intbas(pts[k]) E_k, with its rows
           sorted by shift (stable sort, giving priority to the smaller
           shifts). */
        _nmod_mat_poly_shift_perm(perm, shift, m);
        if (resupd)
        {
            for (slong i = 0; i < m; i++)
                _nmod_vec_set(nmod_mat_entry_ptr(res, i, 0),
                              nmod_mat_entry_ptr(S->R, perm[i], k * n), n);
        }
        else
        {
            pw[0] = 1;
            for (slong t = 1; t < len; t++)
                pw[t] = nmod_mul(pw[t - 1], pt, mod);
            if (n == 1)
            {
                /* q = [pw[0] E_k | pw[1] E_k | ...], then (coefficients of
                   intbas) * q, with inner dimension len*m */
                for (slong i = 0; i < m; i++)
                    q[i] = nmod_mat_entry(E + k, i, 0);
                for (slong t = 1; t < len; t++)
                    _nmod_vec_scalar_mul_nmod(q + t * m, q, m, pw[t], mod);
                nmod_mat_window_init(Wv, S->W, 0, 0, m, len * m);
                nmod_mat_mul_nmod_vec(v, Wv, q, len * m);
                nmod_mat_window_clear(Wv);
                for (slong i = 0; i < m; i++)
                    nmod_mat_entry(res, i, 0) = v[perm[i]];
            }
            else
            {
                /* row i of ev: pw times the coefficients of the row perm[i]
                   of intbas, seen as a len x m matrix (a view, as in
                   nmod_mat_poly_coeff_attach) */
                nmod_mat_struct B[1];
                B->r = len;
                B->c = m;
                B->stride = m;
                B->mod = mod;
                for (slong i = 0; i < m; i++)
                {
                    B->entries = _mintbasis_row(S, perm[i]);
                    nmod_mat_nmod_vec_mul(nmod_mat_entry_ptr(ev, i, 0), pw, len, B);
                }
                nmod_mat_mul(res, ev, E + k);
            }
        }

        const slong rank = _nmod_mat_poly_left_nullspace_compact(X, rows, res, K, P, Q);
        const slong nullity = m - rank;

        if (nullity == m)
        {
            /* zero residual: nothing to do for this point */
            nmod_mat_window_clear(X);
            continue;
        }

        /* rows[0..rank): the pivot rows, rows[rank..m): the others, as row
           indices of intbas; the shifts of the pivot rows increase by 1 */
        for (slong i = 0; i < m; i++)
            rows[i] = perm[rows[i]];
        for (slong j = 0; j < rank; j++)
            shift[rows[j]] += 1;

        /* room for the multiplication of the pivot rows by x - pts[k] */
        _nmod_mat_poly_basis_reserve(S->W, &S->cap, len, len + 1, d + 1);
        _nmod_mat_poly_scratch_fit(T, T->r, S->cap * m + S->R->c);

        /* gather the pivot rows into T: intbas part (ca entries), and the
           residuals at the points k+1, ..., d-1 (cr entries) */
        const slong ca = len * m;
        const slong cr = resupd ? (d - 1 - k) * n : 0;
        const slong tc = ca + cr;
        for (slong j = 0; j < rank; j++)
        {
            ulong * t = nmod_mat_entry_ptr(T, j, 0);
            _nmod_vec_set(t, _mintbasis_row(S, rows[j]), ca);
            if (cr > 0)
                _nmod_vec_set(t + ca, nmod_mat_entry_ptr(S->R, rows[j], (k + 1) * n), cr);
        }

        /* the other rows: row rows[rank+k] += sum_j X[k][j] * row rows[j]
           (the residual at pts[k] becomes zero and is not used anymore) */
        if (nullity > 0 && rank <= lowrank)
        {
            for (slong i = 0; i < nullity; i++)
            {
                ulong * r = _mintbasis_row(S, rows[rank + i]);
                for (slong j = 0; j < rank; j++)
                {
                    const ulong c = nmod_mat_entry(X, i, j);
                    const ulong * t = nmod_mat_entry_ptr(T, j, 0);
                    if (c == 0)
                        continue;
                    _nmod_vec_scalar_addmul_nmod(r, t, ca, c, mod);
                    if (cr > 0)
                        _nmod_vec_scalar_addmul_nmod(nmod_mat_entry_ptr(S->R, rows[rank + i], (k + 1) * n),
                                                     t + ca, cr, c, mod);
                }
            }
        }
        else if (nullity > 0)
        {
            _nmod_mat_poly_scratch_fit(U, m - 1 - lowrank, S->cap * m + S->R->c);
            nmod_mat_window_init(Tw, T, 0, 0, rank, tc);
            nmod_mat_window_init(Uw, U, 0, 0, nullity, tc);
            nmod_mat_mul(Uw, X, Tw);
            for (slong i = 0; i < nullity; i++)
            {
                ulong * r = _mintbasis_row(S, rows[rank + i]);
                const ulong * u = nmod_mat_entry_ptr(U, i, 0);
                _nmod_vec_add(r, r, u, ca, mod);
                if (cr > 0)
                {
                    ulong * rr = nmod_mat_entry_ptr(S->R, rows[rank + i], (k + 1) * n);
                    _nmod_vec_add(rr, rr, u + ca, cr, mod);
                }
            }
            nmod_mat_window_clear(Tw);
            nmod_mat_window_clear(Uw);
        }
        nmod_mat_window_clear(X);

        /* the pivot rows are multiplied by x - pts[k]: the length of intbas
           increases if one of them has degree len-1; their residuals at the
           points pts[j], j > k, are multiplied by pts[j] - pts[k] */
        int grow = 0;
        for (slong j = 0; j < rank && !grow; j++)
            grow = !_nmod_vec_is_zero(nmod_mat_entry_ptr(T, j, (len - 1) * m), m);
        if (cr > 0)
        {
            for (slong j = k + 1; j < d; j++)
            {
                const ulong f = nmod_sub(pts[j], pt, mod);
                for (slong c = 0; c < n; c++)
                    F[(j - k - 1) * n + c] = f;
            }
        }
        const ulong neg_pt = nmod_neg(pt, mod);
        for (slong j = 0; j < rank; j++)
        {
            ulong * r = _mintbasis_row(S, rows[j]);
            const ulong * t = nmod_mat_entry_ptr(T, j, 0);
            _nmod_vec_zero(r, m);
            _nmod_vec_set(r + m, t, ca);
            _nmod_vec_scalar_addmul_nmod(r, t, ca, neg_pt, mod);
            if (cr > 0)
                _nmod_vec_mul(nmod_mat_entry_ptr(S->R, rows[j], (k + 1) * n), t + ca, F, cr, mod);
        }
        len += grow;
    }

    S->len = len;
    flint_free(perm);
    _nmod_vec_clear(pw);
    if (q != NULL)
        _nmod_vec_clear(q);
    _nmod_vec_clear(v);
    _nmod_vec_clear(F);
    nmod_mat_clear(res);
    nmod_mat_clear(K);
    nmod_mat_clear(ev);
    nmod_mat_clear(T);
    nmod_mat_clear(U);
}

void
_nmod_mat_poly_mintbasis(nmod_mat_poly_t intbas, slong * shift,
                         const ulong * pts, const nmod_mat_struct * E, slong d,
                         int resupd, slong lowrank)
{
    const slong m = intbas->r;
    _mintbasis_t S;

    nmod_mat_poly_one(intbas);
    if (d <= 0)
        return;

    _mintbasis_init(S, m, E[0].c, d, E, intbas->mod, resupd);
    _mintbasis_run(S, shift, pts, E, lowrank);

    nmod_mat_poly_fit_length(intbas, S->len);
    _nmod_mat_poly_set_length(intbas, S->len);
    for (slong e = 0; e < S->len; e++)
        for (slong i = 0; i < m; i++)
            _nmod_vec_set(nmod_mat_poly_entry_ptr(intbas, e, i, 0),
                          _mintbasis_row(S, i) + e * m, m);
    _mintbasis_clear(S);
}

void nmod_mat_poly_mintbasis_rescomp(nmod_mat_poly_t intbas,
                                     slong * shift,
                                     const ulong * pts,
                                     const nmod_mat_struct * E,
                                     slong d)
{
    _nmod_mat_poly_mintbasis(intbas, shift, pts, E, d, 0,
                             _nmod_mat_poly_basis_lowrank(intbas->mod));
}

void nmod_mat_poly_mintbasis_resupdate(nmod_mat_poly_t intbas,
                                       slong * shift,
                                       const ulong * pts,
                                       const nmod_mat_struct * E,
                                       slong d)
{
    _nmod_mat_poly_mintbasis(intbas, shift, pts, E, d, 1,
                             _nmod_mat_poly_basis_lowrank(intbas->mod));
}

void nmod_mat_poly_mintbasis(nmod_mat_poly_t intbas,
                             slong * shift,
                             const ulong * pts,
                             const nmod_mat_struct * E,
                             slong d)
{
    /* d = 0: E is never read */
    _nmod_mat_poly_mintbasis(intbas, shift, pts, E, d,
                             (d > 0) && _nmod_mat_poly_mintbasis_resupdate_wanted(intbas->r, E[0].c, d,
                                                                                  NMOD_BITS(intbas->mod)),
                             _nmod_mat_poly_basis_lowrank(intbas->mod));
}

/* Same computation, writing intbas directly, without conversion from
   nmod_mat_poly. */
void nmod_poly_mat_mintbasis(nmod_poly_mat_t intbas,
                             slong * shift,
                             const ulong * pts,
                             const nmod_mat_struct * E,
                             slong d)
{
    const slong m = intbas->r;
    nmod_t mod;
    _mintbasis_t S;

    if (d <= 0)
    {
        nmod_poly_mat_one(intbas);
        return;
    }

    nmod_init(&mod, intbas->modulus);
    _mintbasis_init(S, m, E[0].c, d, E, mod,
                    _nmod_mat_poly_mintbasis_resupdate_wanted(m, E[0].c, d, NMOD_BITS(mod)));
    _mintbasis_run(S, shift, pts, E, _nmod_mat_poly_basis_lowrank(mod));

    for (slong i = 0; i < m; i++)
    {
        const ulong * r = _mintbasis_row(S, i);
        for (slong j = 0; j < m; j++)
        {
            nmod_poly_struct * f = nmod_poly_mat_entry(intbas, i, j);
            nmod_poly_fit_length(f, S->len);
            for (slong e = 0; e < S->len; e++)
                f->coeffs[e] = r[e * m + j];
            _nmod_poly_set_length(f, S->len);
            _nmod_poly_normalise(f);
        }
    }
    _mintbasis_clear(S);
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
