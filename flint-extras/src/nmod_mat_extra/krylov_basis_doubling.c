/*
    Copyright (C) 2026 Vincent Neiger

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
#include <flint/nmod_poly.h>
#include <flint/nmod_mat.h>
#include "nmod_mat_extra/krylov_impl.h"

/*
    Powers M^(2^l) of the operator, applied to rows. Dense kind: the
    matrix itself, by repeated squaring. Jordan kinds: on a block of size t
    with eigenvalue x, the row v is multiplied by the polynomial
    h = (x + X)^(2^l) mod X^t, in the sense v <- coefficients of
    v(X) h(X) mod X^t; h is computed by repeated squaring.
*/
typedef struct
{
    const _nmod_mat_krylov_op_struct * op;
    slong e;                /* current exponent 2^l */
    nmod_mat_t Mp;          /* dense: M^e */
    nn_ptr * h;             /* Jordan: per block, h of length size */
}
_krylov_pow_struct;

static void
_krylov_pow_init(_krylov_pow_struct * P, const _nmod_mat_krylov_op_struct * op)
{
    slong b;

    P->op = op;
    P->e = 1;
    P->h = NULL;

    if (op->kind == _NMOD_MAT_KRYLOV_DENSE)
    {
        nmod_mat_init_set(P->Mp, op->M);
    }
    else if (op->kind != _NMOD_MAT_KRYLOV_SHIFT)
    {
        P->h = flint_malloc(op->nblocks * sizeof(nn_ptr));
        for (b = 0; b < op->nblocks; b++)
        {
            slong t = op->off[b + 1] - op->off[b];
            P->h[b] = _nmod_vec_init(t);
            _nmod_vec_zero(P->h[b], t);
            P->h[b][0] = op->eig[b];
            if (t > 1)
                P->h[b][1] = 1;
        }
    }
}

static void
_krylov_pow_clear(_krylov_pow_struct * P)
{
    slong b;

    if (P->op->kind == _NMOD_MAT_KRYLOV_DENSE)
        nmod_mat_clear(P->Mp);
    else if (P->h != NULL)
    {
        for (b = 0; b < P->op->nblocks; b++)
            _nmod_vec_clear(P->h[b]);
        flint_free(P->h);
    }
}

/* e <- 2 e */
static void
_krylov_pow_square(_krylov_pow_struct * P)
{
    const _nmod_mat_krylov_op_struct * op = P->op;
    slong b;

    if (op->kind == _NMOD_MAT_KRYLOV_DENSE)
    {
        nmod_mat_t T;
        nmod_mat_init(T, op->n, op->n, op->mod.n);
        nmod_mat_mul(T, P->Mp, P->Mp);
        nmod_mat_swap(T, P->Mp);
        nmod_mat_clear(T);
    }
    else if (P->h != NULL)
    {
        for (b = 0; b < op->nblocks; b++)
        {
            slong t = op->off[b + 1] - op->off[b];
            nn_ptr g = _nmod_vec_init(t);
            _nmod_poly_mullow(g, P->h[b], t, P->h[b], t, t, op->mod);
            _nmod_vec_set(P->h[b], g, t);
            _nmod_vec_clear(g);
        }
    }

    P->e *= 2;
}

/* Out = In * M^e (Out and In distinct) */
static void
_krylov_pow_apply(nmod_mat_t Out, const nmod_mat_t In, const _krylov_pow_struct * P)
{
    const _nmod_mat_krylov_op_struct * op = P->op;
    slong a, b, c;

    if (In->r == 0)
        return;

    if (op->kind == _NMOD_MAT_KRYLOV_DENSE)
    {
        nmod_mat_mul(Out, In, P->Mp);
        return;
    }

    for (a = 0; a < In->r; a++)
    {
        nn_srcptr v = nmod_mat_entry_ptr(In, a, 0);
        nn_ptr w = nmod_mat_entry_ptr(Out, a, 0);

        for (b = 0; b < op->nblocks; b++)
        {
            slong o = op->off[b], t = op->off[b + 1] - op->off[b];

            if (op->kind == _NMOD_MAT_KRYLOV_SHIFT)
            {
                slong e = FLINT_MIN(P->e, t);
                for (c = 0; c < e; c++)
                    w[o + c] = 0;
                for (c = e; c < t; c++)
                    w[o + c] = v[o + c - e];
            }
            else if (t == 1)
            {
                w[o] = nmod_mul(v[o], P->h[b][0], op->mod);
            }
            else
            {
                _nmod_poly_mullow(w + o, v + o, t, P->h[b], t, t, op->mod);
            }
        }
    }
}

/* pairs (i, j) sorted by priority (s_i + j, i), for normalized shifts */
static int
_krylov_pair_less(slong i1, slong j1, slong i2, slong j2, const slong * shifts)
{
    slong k1 = shifts[i1] + j1;
    slong k2 = shifts[i2] + j2;
    return k1 < k2 || (k1 == k2 && i1 < i2);
}

slong
_nmod_mat_krylov_basis_doubling(nmod_mat_t B, slong * rows, slong * exps,
        slong * deltas, int * ok, const nmod_mat_t E,
        const _nmod_mat_krylov_op_t op, const slong * shifts,
        const slong * degrees)
{
    slong m = E->r, n = E->c, i, a, t, r, maxd, cnt, nnew, nmerged;
    slong * d, * order, * nsh, * ri, * rj, * Q, * mi, * mj, * src;
    char * indep;
    nmod_mat_t X, Cur, New, Mrg;
    _krylov_pow_struct P;

    if (op->n != n)
        flint_throw(FLINT_ERROR, "nmod_mat_krylov_basis: incompatible dimensions\n");

    if (ok != NULL)
        *ok = 1;

    if (n == 0 || m == 0)
    {
        if (deltas != NULL)
            for (i = 0; i < m; i++)
                deltas[i] = 0;
        return 0;
    }

    d = flint_malloc((3 * m + 1) * sizeof(slong));
    order = d + m;
    nsh = order + m;
    _nmod_mat_krylov_degrees(d, m, degrees, n, 1);
    _nmod_mat_krylov_normalize_shifts(nsh, shifts, m, n);

    maxd = 0;
    for (i = 0; i < m; i++)
        maxd = FLINT_MAX(maxd, d[i]);

    ri = flint_malloc((7 * (m + 2 * n) + 1) * sizeof(slong));
    rj = ri + (m + 2 * n);
    Q = rj + (m + 2 * n);
    mi = Q + (m + 2 * n);
    mj = mi + (m + 2 * n);
    src = mj + (m + 2 * n);
    indep = flint_malloc(m + 2 * n + 1);

    /* step 1: row rank profile of E, rows sorted by priority (s_i, i) */
    {
        slong * zeros = flint_calloc(m + 1, sizeof(slong));
        _nmod_mat_krylov_order(order, ri, m, nsh, zeros);
        flint_free(zeros);
    }

    nmod_mat_init(X, m, n, E->mod.n);
    for (a = 0; a < m; a++)
        memcpy(nmod_mat_entry_ptr(X, a, 0), nmod_mat_entry_ptr(E, order[a], 0),
               n * sizeof(ulong));
    r = _nmod_mat_row_rank_profile(indep, X);
    nmod_mat_clear(X);

    nmod_mat_init(Cur, FLINT_MAX(r, 1), n, E->mod.n);
    t = 0;
    for (a = 0; a < m; a++)
    {
        if (indep[a])
        {
            ri[t] = order[a];
            rj[t] = 0;
            memcpy(nmod_mat_entry_ptr(Cur, t, 0),
                   nmod_mat_entry_ptr(E, order[a], 0), n * sizeof(ulong));
            t++;
        }
    }

    /* step 2: doubling */
    _krylov_pow_init(&P, op);

    while (r > 0 && P.e <= maxd)
    {
        nmod_mat_t Sel;

        /* new rows: (i, j + e) for the profile rows (i, j) */
        nnew = 0;
        for (t = 0; t < r; t++)
            if (rj[t] + P.e <= d[ri[t]])
                nnew++;

        if (nnew == 0)
            break;

        nmod_mat_init(Sel, nnew, n, E->mod.n);
        nmod_mat_init(New, nnew, n, E->mod.n);
        cnt = 0;
        for (t = 0; t < r; t++)
        {
            if (rj[t] + P.e <= d[ri[t]])
            {
                memcpy(nmod_mat_entry_ptr(Sel, cnt, 0),
                       nmod_mat_entry_ptr(Cur, t, 0), n * sizeof(ulong));
                mi[r + cnt] = ri[t];
                mj[r + cnt] = rj[t] + P.e;
                cnt++;
            }
        }
        _krylov_pow_apply(New, Sel, &P);
        nmod_mat_clear(Sel);

        /* merge the two sorted lists (both are sorted by priority) */
        {
            slong u = 0, v = 0;
            nmerged = r + nnew;
            for (t = 0; t < nmerged; t++)
            {
                if (v >= nnew || (u < r && _krylov_pair_less(ri[u], rj[u],
                                            mi[r + v], mj[r + v], nsh)))
                    src[t] = u++;
                else
                    src[t] = r + v++;
            }
        }

        nmod_mat_init(Mrg, nmerged, n, E->mod.n);
        for (t = 0; t < nmerged; t++)
        {
            nn_srcptr row = (src[t] < r) ? nmod_mat_entry_ptr(Cur, src[t], 0)
                                         : nmod_mat_entry_ptr(New, src[t] - r, 0);
            memcpy(nmod_mat_entry_ptr(Mrg, t, 0), row, n * sizeof(ulong));
        }

        /* row rank profile of the merged rows */
        {
            slong r_old = r;

            nmod_mat_init_set(X, Mrg);
            r = _nmod_mat_row_rank_profile(indep, X);
            nmod_mat_clear(X);

            nmod_mat_clear(Cur);
            nmod_mat_init(Cur, FLINT_MAX(r, 1), n, E->mod.n);
            cnt = 0;
            for (t = 0; t < nmerged; t++)
            {
                if (indep[t])
                {
                    memcpy(nmod_mat_entry_ptr(Cur, cnt, 0),
                           nmod_mat_entry_ptr(Mrg, t, 0), n * sizeof(ulong));
                    cnt++;
                }
            }

            /* pairs of the kept rows, written to temporary arrays */
            {
                slong * ni = flint_malloc((2 * r + 1) * sizeof(slong));
                slong * nj = ni + r;
                cnt = 0;
                for (t = 0; t < nmerged; t++)
                {
                    if (indep[t])
                    {
                        slong sidx = src[t];
                        ni[cnt] = (sidx < r_old) ? ri[sidx] : mi[sidx];
                        nj[cnt] = (sidx < r_old) ? rj[sidx] : mj[sidx];
                        cnt++;
                    }
                }
                for (t = 0; t < r; t++)
                {
                    ri[t] = ni[t];
                    rj[t] = nj[t];
                }
                flint_free(ni);
            }
        }

        nmod_mat_clear(Mrg);
        nmod_mat_clear(New);

        if (2 * P.e > maxd)
            break;

        _krylov_pow_square(&P);
    }

    _krylov_pow_clear(&P);

    /* output */
    for (t = 0; t < r; t++)
    {
        memcpy(nmod_mat_entry_ptr(B, t, 0), nmod_mat_entry_ptr(Cur, t, 0),
               n * sizeof(ulong));
        rows[t] = ri[t];
        exps[t] = rj[t];
    }

    if (ok != NULL)
        *ok = 1;

    for (i = 0; i < m; i++)
        order[i] = 0;
    for (t = 0; t < r; t++)
        order[ri[t]] = FLINT_MAX(order[ri[t]], rj[t] + 1);
    for (i = 0; i < m; i++)
    {
        if (deltas != NULL)
            deltas[i] = order[i];
        if (order[i] > d[i] && ok != NULL)
            *ok = 0;
    }

    nmod_mat_clear(Cur);
    flint_free(d);
    flint_free(ri);
    flint_free(indep);

    return r;
}
