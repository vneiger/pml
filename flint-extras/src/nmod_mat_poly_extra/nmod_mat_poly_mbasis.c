/*
    Copyright (C) 2025 Vincent Neiger, Éric Schost, Kevin Tran
    Copyright (C) 2026 Gilles Villard

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include <flint/nmod.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>
#include "nmod_mat_poly.h"
#include "nmod_poly_mat_approximant.h"
#include "nmod_mat_poly_extra/impl.h"

/*
    Storage. The rows of the approximant basis are stored as vectors of
    coefficients: row i of the matrix W is [coefficient 0 of row i of appbas |
    coefficient 1 | ...], with m entries per coefficient. W has room for cap
    coefficients, a bound on the length of appbas which is increased when
    needed (initially, about the length for generic input with uniform shift,
    n*order/m + 1). In the "resupdate" variant, row i of the matrix R is the
    row i of the residual series, that is, of the coefficients ord, ord+1,
    ..., order-1 of appbas*matp (n entries per coefficient, coefficient e at
    column e n). With this storage:

    - in the "rescomp" variant, the residual (coefficient ord of appbas*matp)
      is a single matrix product, with inner dimension m*len (len the length
      of appbas): the rows of appbas times the coefficients ord, ord-1, ...
      of matp stacked vertically, which is a window of the coefficients of
      matp stacked in reverse order (Mrev below);

    - the update of the basis (and of the residual series) by the left
      nullspace of the residual is a single product X * (pivot rows), over all
      the coefficients at once, instead of one product per coefficient of
      appbas, plus a matrix addition each time;

    - there is no row permutation of appbas (or of the residual series): the
      pivot rows are gathered into a contiguous matrix, which is also what
      their multiplication by x copies back.

    The row operations are the same as those of the previous implementation
    (on nmod_mat_poly_t, with permutations of the rows), in the same order, so
    the output (appbas and shift) is the same, and the two variants give the
    same output.
*/

typedef struct
{
    slong m;            /* rows of matp */
    slong n;            /* columns of matp */
    slong order;
    nmod_t mod;
    int resupd;         /* variant: residual series updated, or recomputed */
    nmod_mat_t W;       /* m x (cap m): the rows of appbas */
    slong cap;          /* capacity of W in coefficients, at most order + 1 */
    slong len;          /* length of appbas */
    nmod_mat_t R;       /* (resupd) m x (n order): the residual series */
    slong L;            /* number of coefficients of matp that are used */
    nmod_mat_t Mrev;    /* (rescomp) (m L) x n, block b is coefficient L-1-b */
}
_mbasis_struct;

typedef _mbasis_struct _mbasis_t[1];

/* W starts as the identity; the coefficients of the input are then written
   with _mbasis_input_row */
static void
_mbasis_init(_mbasis_t S, slong m, slong n, slong order, slong L, nmod_t mod,
             int resupd)
{
    S->m = m;
    S->n = n;
    S->order = order;
    S->mod = mod;
    S->resupd = resupd;
    S->cap = FLINT_MIN(order + 1, (n * order + m - 1) / m + 2);
    nmod_mat_init(S->W, m, S->cap * m, mod.n);
    for (slong i = 0; i < m; i++)
        nmod_mat_entry(S->W, i, i) = 1;
    S->len = 1;
    nmod_mat_init(S->R, resupd ? m : 0, resupd ? n * order : 0, mod.n);
    S->L = L;
    nmod_mat_init(S->Mrev, resupd ? 0 : m * L, n, mod.n);
}

static void
_mbasis_clear(_mbasis_t S)
{
    nmod_mat_clear(S->W);
    nmod_mat_clear(S->R);
    nmod_mat_clear(S->Mrev);
}

/* where row i of the coefficient d < L of the input goes (n entries) */
static inline ulong *
_mbasis_input_row(_mbasis_t S, slong d, slong i)
{
    if (S->resupd)
        return nmod_mat_entry_ptr(S->R, i, d * S->n);
    else
        return nmod_mat_entry_ptr(S->Mrev, (S->L - 1 - d) * S->m + i, 0);
}

/* row i of the coefficient d < len of appbas (m entries) */
static inline const ulong *
_mbasis_output_row(const _mbasis_t S, slong d, slong i)
{
    return nmod_mat_entry_ptr(S->W, i, d * S->m);
}

void
_nmod_mat_poly_basis_reserve(nmod_mat_t W, slong * cap, slong len,
                             slong need, slong maxcap)
{
    if (need <= *cap)
        return;

    const slong m = W->r;
    const slong newcap = FLINT_MIN(maxcap, FLINT_MAX(need, 2 * *cap));
    nmod_mat_t W2;
    nmod_mat_init(W2, m, newcap * m, W->mod.n);
    for (slong i = 0; i < m; i++)
        _nmod_vec_set(nmod_mat_entry_ptr(W2, i, 0), nmod_mat_entry_ptr(W, i, 0), len * m);
    nmod_mat_swap(W, W2);
    nmod_mat_clear(W2);
    *cap = newcap;
}

void
_nmod_mat_poly_shift_perm(slong * perm, const slong * shift, slong m)
{
    for (slong i = 1; i < m; i++)
    {
        const slong p = perm[i], v = shift[p];
        slong j = i;
        while (j > 0 && (shift[perm[j - 1]] > v
                         || (shift[perm[j - 1]] == v && perm[j - 1] > p)))
        {
            perm[j] = perm[j - 1];
            j--;
        }
        perm[j] = p;
    }
}

slong
_nmod_mat_poly_left_nullspace_compact(nmod_mat_t X, slong * rows, nmod_mat_t res,
                                      const nmod_mat_t K, slong * P, slong * Q)
{
    const slong m = res->r;
    const nmod_t mod = res->mod;
    slong rank;

    if (res->c == 1)
    {
        /* one column: the pivot is the first nonzero entry */
        slong i0 = 0;
        while (i0 < m && nmod_mat_entry(res, i0, 0) == 0)
            i0++;
        rank = (i0 < m);
        nmod_mat_window_init(X, K, 0, 0, m - rank, rank);
        if (rank == 1 && m > 1)
        {
            const ulong c = nmod_neg(nmod_inv(nmod_mat_entry(res, i0, 0), mod), mod);
            rows[0] = i0;
            for (slong i = 0; i < i0; i++)
            {
                rows[1 + i] = i;
                nmod_mat_entry(X, i, 0) = 0;
            }
            for (slong i = i0 + 1; i < m; i++)
            {
                rows[i] = i;
                nmod_mat_entry(X, i - 1, 0) = nmod_mul(nmod_mat_entry(res, i, 0), c, mod);
            }
        }
    }
    else
    {
        rank = nmod_mat_ple(P, Q, res, 0);
        nmod_mat_window_init(X, K, 0, 0, m - rank, rank);
        if (rank > 0 && rank < m)
            nmod_mat_ple_left_nullspace_compact(X, rows, res, P, rank);
    }

    if (rank == m)
        for (slong i = 0; i < m; i++)
            rows[i] = i;

    return rank;
}

/* The iterations of mbasis: on input, S holds the identity and the input;
   on output, appbas (and S->len); shift is updated. */
static void
_mbasis_run(_mbasis_t S, slong * shift, slong lowrank)
{
    const slong m = S->m, n = S->n, order = S->order, L = S->L;
    const int resupd = S->resupd;
    const nmod_t mod = S->mod;
    slong len = S->len;

    /* T: the pivot rows (rank <= min(m, n)), U: the combinations of them
       added to the other rows (allocated when needed: not for rank <=
       lowrank); their columns: the coefficients of appbas, and the residual
       series (resupd) */
    nmod_mat_t res, K, X, Aw, Mw, T, U, Tw, Uw;
    nmod_mat_init(res, m, n, mod.n);
    nmod_mat_init(K, m, m, mod.n);
    nmod_mat_init(T, FLINT_MIN(m, n), S->cap * m + S->R->c, mod.n);
    nmod_mat_init(U, 0, 0, mod.n);

    /* perm: rows sorted by shift; rows: the pivot rows of the residual, then
       the others; P, Q: scratch for the PLE decomposition; v: for n = 1 */
    slong * perm = (slong *) flint_malloc((3 * m + n) * sizeof(slong));
    slong * rows = perm + m;
    slong * P = rows + m;
    slong * Q = P + m;
    nn_ptr v = _nmod_vec_init(m);
    for (slong i = 0; i < m; i++)
        perm[i] = i;

    for (slong ord = 0; ord < order; ord++)
    {
        /* Here appbas is an s-ordered weak Popov approximant basis at order
           ord, where s is the input shift, and shift = rdeg_s(appbas). The
           residual: coefficient ord of appbas*matp, with its rows sorted by
           shift (stable sort, giving priority to the smaller shifts). */
        _nmod_mat_poly_shift_perm(perm, shift, m);
        if (resupd)
        {
            for (slong i = 0; i < m; i++)
                _nmod_vec_set(nmod_mat_entry_ptr(res, i, 0),
                              nmod_mat_entry_ptr(S->R, perm[i], ord * n), n);
        }
        else
        {
            /* sum of appbas[d] * matp[ord-d] for lb <= d < ub */
            const slong lb = FLINT_MAX(0, ord + 1 - L);
            const slong ub = FLINT_MIN(ord + 1, len);
            if (lb >= ub)
                nmod_mat_zero(res);
            else
            {
                nmod_mat_window_init(Aw, S->W, 0, lb * m, m, ub * m);
                nmod_mat_window_init(Mw, S->Mrev, (L - 1 - ord + lb) * m, 0,
                                     (L - 1 - ord + ub) * m, n);
                if (n == 1)
                {
                    /* matrix-vector product (for few rows, nmod_mat_mul
                       does not use it); the column Mw is contiguous */
                    nmod_mat_mul_nmod_vec(v, Aw, nmod_mat_entry_ptr(Mw, 0, 0), Aw->c);
                    for (slong i = 0; i < m; i++)
                        nmod_mat_entry(res, i, 0) = v[perm[i]];
                }
                else
                {
                    nmod_mat_mul(res, Aw, Mw);
                    nmod_mat_permute_rows(res, perm, NULL);
                }
                nmod_mat_window_clear(Aw);
                nmod_mat_window_clear(Mw);
            }
        }

        const slong rank = _nmod_mat_poly_left_nullspace_compact(X, rows, res, K, P, Q);
        const slong nullity = m - rank;

        if (nullity == 0)
        {
            /* the residual has full rank: the basis at order `order` is
               x^(order-ord) times the current one */
            const slong k = order - ord;
            _nmod_mat_poly_basis_reserve(S->W, &S->cap, len, len + k, order + 1);
            for (slong i = 0; i < m; i++)
            {
                ulong * r = nmod_mat_entry_ptr(S->W, i, 0);
                memmove(r + k * m, r, len * m * sizeof(ulong));
                _nmod_vec_zero(r, k * m);
                shift[i] += k;
            }
            len += k;
            nmod_mat_window_clear(X);
            break;
        }

        if (nullity == m)
        {
            /* zero residual: nothing to do at this order */
            nmod_mat_window_clear(X);
            continue;
        }

        /* rows[0..rank): the pivot rows, rows[rank..m): the others, as row
           indices of appbas; the shifts of the pivot rows increase by 1 */
        for (slong i = 0; i < m; i++)
            rows[i] = perm[rows[i]];
        for (slong j = 0; j < rank; j++)
            shift[rows[j]] += 1;

        /* room for the multiplication of the pivot rows by x */
        _nmod_mat_poly_basis_reserve(S->W, &S->cap, len, len + 1, order + 1);
        _nmod_mat_poly_scratch_fit(T, T->r, S->cap * m + S->R->c);

        /* gather the pivot rows into T: appbas part (ca entries), and
           residual part from the coefficient ord on (cr entries) */
        const slong ca = len * m;
        const slong cr = resupd ? (order - ord) * n : 0;
        const slong tc = ca + cr;
        for (slong j = 0; j < rank; j++)
        {
            ulong * t = nmod_mat_entry_ptr(T, j, 0);
            _nmod_vec_set(t, nmod_mat_entry_ptr(S->W, rows[j], 0), ca);
            if (resupd)
                _nmod_vec_set(t + ca, nmod_mat_entry_ptr(S->R, rows[j], ord * n), cr);
        }

        /* the other rows: row rows[rank+k] += sum_j X[k][j] * row rows[j];
           in the residual part, the coefficient ord becomes zero and is not
           used anymore: it is skipped */
        if (rank <= lowrank)
        {
            for (slong k = 0; k < nullity; k++)
            {
                ulong * r = nmod_mat_entry_ptr(S->W, rows[rank + k], 0);
                for (slong j = 0; j < rank; j++)
                {
                    const ulong c = nmod_mat_entry(X, k, j);
                    const ulong * t = nmod_mat_entry_ptr(T, j, 0);
                    if (c == 0)
                        continue;
                    _nmod_vec_scalar_addmul_nmod(r, t, ca, c, mod);
                    if (cr > n)
                        _nmod_vec_scalar_addmul_nmod(nmod_mat_entry_ptr(S->R, rows[rank + k], (ord + 1) * n),
                                                     t + ca + n, cr - n, c, mod);
                }
            }
        }
        else
        {
            _nmod_mat_poly_scratch_fit(U, m - 1 - lowrank, S->cap * m + S->R->c);
            nmod_mat_window_init(Tw, T, 0, 0, rank, tc);
            nmod_mat_window_init(Uw, U, 0, 0, nullity, tc);
            nmod_mat_mul(Uw, X, Tw);
            for (slong k = 0; k < nullity; k++)
            {
                ulong * r = nmod_mat_entry_ptr(S->W, rows[rank + k], 0);
                const ulong * u = nmod_mat_entry_ptr(U, k, 0);
                _nmod_vec_add(r, r, u, ca, mod);
                if (cr > n)
                {
                    ulong * rr = nmod_mat_entry_ptr(S->R, rows[rank + k], (ord + 1) * n);
                    _nmod_vec_add(rr, rr, u + ca + n, cr - n, mod);
                }
            }
            nmod_mat_window_clear(Tw);
            nmod_mat_window_clear(Uw);
        }
        nmod_mat_window_clear(X);

        /* the pivot rows are multiplied by x: the length of appbas increases
           if one of them has degree len-1; in the residual part, the
           coefficient that would reach degree `order` is dropped */
        int grow = 0;
        for (slong j = 0; j < rank && !grow; j++)
            grow = !_nmod_vec_is_zero(nmod_mat_entry_ptr(T, j, (len - 1) * m), m);
        for (slong j = 0; j < rank; j++)
        {
            ulong * r = nmod_mat_entry_ptr(S->W, rows[j], 0);
            const ulong * t = nmod_mat_entry_ptr(T, j, 0);
            _nmod_vec_zero(r, m);
            _nmod_vec_set(r + m, t, ca);
            if (cr > n)
                _nmod_vec_set(nmod_mat_entry_ptr(S->R, rows[j], (ord + 1) * n), t + ca, cr - n);
        }
        len += grow;
    }

    S->len = len;
    flint_free(perm);
    _nmod_vec_clear(v);
    nmod_mat_clear(res);
    nmod_mat_clear(K);
    nmod_mat_clear(T);
    nmod_mat_clear(U);
}

void
_nmod_mat_poly_mbasis(nmod_mat_poly_t appbas, slong * shift,
                      const nmod_mat_poly_t matp, slong order,
                      int resupd, slong lowrank)
{
    const slong m = matp->r, n = matp->c;
    _mbasis_t S;

    nmod_mat_poly_one(appbas);
    if (nmod_mat_poly_is_zero(matp) || order <= 0)
        return;

    const slong L = FLINT_MIN(matp->length, order);
    _mbasis_init(S, m, n, order, L, matp->mod, resupd);
    for (slong d = 0; d < L; d++)
        for (slong i = 0; i < m; i++)
            _nmod_vec_set(_mbasis_input_row(S, d, i),
                          nmod_mat_poly_entry_ptr(matp, d, i, 0), n);

    _mbasis_run(S, shift, lowrank);

    nmod_mat_poly_fit_length(appbas, S->len);
    _nmod_mat_poly_set_length(appbas, S->len);
    for (slong d = 0; d < S->len; d++)
        for (slong i = 0; i < m; i++)
            _nmod_vec_set(nmod_mat_poly_entry_ptr(appbas, d, i, 0),
                          _mbasis_output_row(S, d, i), m);
    _mbasis_clear(S);
}

void nmod_mat_poly_mbasis_rescomp(nmod_mat_poly_t appbas,
                                  slong * shift,
                                  const nmod_mat_poly_t matp,
                                  slong order)
{
    _nmod_mat_poly_mbasis(appbas, shift, matp, order, 0,
                          _nmod_mat_poly_basis_lowrank(matp->mod));
}

void nmod_mat_poly_mbasis_resupdate(nmod_mat_poly_t appbas,
                                    slong * shift,
                                    const nmod_mat_poly_t matp,
                                    slong order)
{
    _nmod_mat_poly_mbasis(appbas, shift, matp, order, 1,
                          _nmod_mat_poly_basis_lowrank(matp->mod));
}

void nmod_mat_poly_mbasis(nmod_mat_poly_t appbas,
                          slong * shift,
                          const nmod_mat_poly_t matp,
                          slong order)
{
    _nmod_mat_poly_mbasis(appbas, shift, matp, order,
                          _nmod_mat_poly_mbasis_resupdate_wanted(matp->r, matp->c, order,
                                                                 NMOD_BITS(matp->mod)),
                          _nmod_mat_poly_basis_lowrank(matp->mod));
}

/* Same computation, reading pmat and writing appbas directly, without
   conversions to nmod_mat_poly. */
void nmod_poly_mat_mbasis(nmod_poly_mat_t appbas,
                          slong * shift,
                          const nmod_poly_mat_t pmat,
                          ulong order)
{
    const slong m = pmat->r, n = pmat->c;
    slong L = 0;
    nmod_t mod;
    _mbasis_t S;

    nmod_init(&mod, pmat->modulus);

    for (slong i = 0; i < m; i++)
        for (slong j = 0; j < n; j++)
            L = FLINT_MAX(L, nmod_poly_length(nmod_poly_mat_entry(pmat, i, j)));
    L = FLINT_MIN(L, (slong) order);

    if (L <= 0)
    {
        nmod_poly_mat_one(appbas);
        return;
    }

    _mbasis_init(S, m, n, order, L, mod,
                 _nmod_mat_poly_mbasis_resupdate_wanted(m, n, order, NMOD_BITS(mod)));
    for (slong i = 0; i < m; i++)
        for (slong j = 0; j < n; j++)
        {
            const nmod_poly_struct * f = nmod_poly_mat_entry(pmat, i, j);
            const slong l = FLINT_MIN(f->length, L);
            for (slong d = 0; d < l; d++)
                _mbasis_input_row(S, d, i)[j] = f->coeffs[d];
        }

    _mbasis_run(S, shift, _nmod_mat_poly_basis_lowrank(mod));

    for (slong i = 0; i < m; i++)
        for (slong j = 0; j < m; j++)
        {
            nmod_poly_struct * f = nmod_poly_mat_entry(appbas, i, j);
            nmod_poly_fit_length(f, S->len);
            for (slong d = 0; d < S->len; d++)
                f->coeffs[d] = _mbasis_output_row(S, d, i)[j];
            _nmod_poly_set_length(f, S->len);
            _nmod_poly_normalise(f);
        }
    _mbasis_clear(S);
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
