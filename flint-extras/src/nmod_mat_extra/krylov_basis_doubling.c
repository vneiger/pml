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
    matrix itself, by repeated squaring. Jordan kinds: on the block b, of
    size D_b and eigenvalue x_b, the row v is multiplied by the polynomial
    h = (x_b + X)^(2^l) mod X^(D_b), in the sense v <- coefficients of
    v(X) h(X) mod X^(D_b); h is computed by repeated squaring.
*/
typedef struct
{
    const _nmod_mat_krylov_op_struct * op;
    slong q;                /* current exponent 2^l */
    nmod_mat_t Mp;          /* dense: M^q */
    nn_ptr * h;             /* Jordan: per block, h of length size */
}
_krylov_pow_struct;

static void
_krylov_pow_init(_krylov_pow_struct * P, const _nmod_mat_krylov_op_struct * op)
{
    slong b;

    P->op = op;
    P->q = 1;
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

/* q <- 2 q */
static void
_krylov_pow_square(_krylov_pow_struct * P)
{
    const _nmod_mat_krylov_op_struct * op = P->op;
    slong b;

    if (op->kind == _NMOD_MAT_KRYLOV_DENSE)
    {
        nmod_mat_t T;
        nmod_mat_init(T, op->D, op->D, op->mod.n);
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

    P->q *= 2;
}

/* Out = In * M^q (Out and In distinct) */
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
                slong qe = FLINT_MIN(P->q, t);
                for (c = 0; c < qe; c++)
                    w[o + c] = 0;
                for (c = qe; c < t; c++)
                    w[o + c] = v[o + c - qe];
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

/* pairs (i, e) sorted by priority (s_i + e, i), for normalized shift */
static int
_krylov_pair_less(slong i1, slong e1, slong i2, slong e2, const slong * shift)
{
    slong k1 = shift[i1] + e1;
    slong k2 = shift[i2] + e2;
    return k1 < k2 || (k1 == k2 && i1 < i2);
}

slong
_nmod_mat_krylov_basis_doubling(nmod_mat_t B, slong * rows, slong * exps,
        slong * deltas, int * ok, const nmod_mat_t V,
        const _nmod_mat_krylov_op_t op, const slong * shift,
        const slong * degrees)
{
    slong m = V->r, D = V->c, i, a, t, r, maxbeta, cnt, nnew, nmerged;
    slong * beta, * order, * nsh, * ri, * re, * Q, * mi, * me, * src;
    char * indep;
    nmod_mat_t X, Cur, New, Mrg;
    _krylov_pow_struct P;

    if (op->D != D)
        flint_throw(FLINT_ERROR, "nmod_mat_krylov_basis: incompatible dimensions\n");

    if (ok != NULL)
        *ok = 1;

    if (D == 0 || m == 0)
    {
        if (deltas != NULL)
            for (i = 0; i < m; i++)
                deltas[i] = 0;
        return 0;
    }

    beta = flint_malloc((3 * m + 1) * sizeof(slong));
    order = beta + m;
    nsh = order + m;
    _nmod_mat_krylov_degrees(beta, m, degrees, D, 1);
    _nmod_mat_krylov_normalize_shift(nsh, shift, m, D);

    maxbeta = 0;
    for (i = 0; i < m; i++)
        maxbeta = FLINT_MAX(maxbeta, beta[i]);

    ri = flint_malloc((7 * (m + 2 * D) + 1) * sizeof(slong));
    re = ri + (m + 2 * D);
    Q = re + (m + 2 * D);
    mi = Q + (m + 2 * D);
    me = mi + (m + 2 * D);
    src = me + (m + 2 * D);
    indep = flint_malloc(m + 2 * D + 1);

    /* step 1: row rank profile of V, rows sorted by priority (s_i, i) */
    {
        slong * zeros = flint_calloc(m + 1, sizeof(slong));
        _nmod_mat_krylov_order(order, ri, m, nsh, zeros);
        flint_free(zeros);
    }

    nmod_mat_init(X, m, D, V->mod.n);
    for (a = 0; a < m; a++)
        memcpy(nmod_mat_entry_ptr(X, a, 0), nmod_mat_entry_ptr(V, order[a], 0),
               D * sizeof(ulong));
    r = _nmod_mat_row_rank_profile(indep, X);
    nmod_mat_clear(X);

    nmod_mat_init(Cur, FLINT_MAX(r, 1), D, V->mod.n);
    t = 0;
    for (a = 0; a < m; a++)
    {
        if (indep[a])
        {
            ri[t] = order[a];
            re[t] = 0;
            memcpy(nmod_mat_entry_ptr(Cur, t, 0),
                   nmod_mat_entry_ptr(V, order[a], 0), D * sizeof(ulong));
            t++;
        }
    }

    /* step 2: doubling */
    _krylov_pow_init(&P, op);

    while (r > 0 && P.q <= maxbeta)
    {
        nmod_mat_t Sel;

        /* new rows: (i, e + q) for the profile rows (i, e) */
        nnew = 0;
        for (t = 0; t < r; t++)
            if (re[t] + P.q <= beta[ri[t]])
                nnew++;

        if (nnew == 0)
            break;

        nmod_mat_init(Sel, nnew, D, V->mod.n);
        nmod_mat_init(New, nnew, D, V->mod.n);
        cnt = 0;
        for (t = 0; t < r; t++)
        {
            if (re[t] + P.q <= beta[ri[t]])
            {
                memcpy(nmod_mat_entry_ptr(Sel, cnt, 0),
                       nmod_mat_entry_ptr(Cur, t, 0), D * sizeof(ulong));
                mi[r + cnt] = ri[t];
                me[r + cnt] = re[t] + P.q;
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
                if (v >= nnew || (u < r && _krylov_pair_less(ri[u], re[u],
                                            mi[r + v], me[r + v], nsh)))
                    src[t] = u++;
                else
                    src[t] = r + v++;
            }
        }

        nmod_mat_init(Mrg, nmerged, D, V->mod.n);
        for (t = 0; t < nmerged; t++)
        {
            nn_srcptr row = (src[t] < r) ? nmod_mat_entry_ptr(Cur, src[t], 0)
                                         : nmod_mat_entry_ptr(New, src[t] - r, 0);
            memcpy(nmod_mat_entry_ptr(Mrg, t, 0), row, D * sizeof(ulong));
        }

        /* row rank profile of the merged rows */
        {
            slong r_old = r;

            nmod_mat_init_set(X, Mrg);
            r = _nmod_mat_row_rank_profile(indep, X);
            nmod_mat_clear(X);

            nmod_mat_clear(Cur);
            nmod_mat_init(Cur, FLINT_MAX(r, 1), D, V->mod.n);
            cnt = 0;
            for (t = 0; t < nmerged; t++)
            {
                if (indep[t])
                {
                    memcpy(nmod_mat_entry_ptr(Cur, cnt, 0),
                           nmod_mat_entry_ptr(Mrg, t, 0), D * sizeof(ulong));
                    cnt++;
                }
            }

            /* pairs of the kept rows, written to temporary arrays */
            {
                slong * ni = flint_malloc((2 * r + 1) * sizeof(slong));
                slong * ne = ni + r;
                cnt = 0;
                for (t = 0; t < nmerged; t++)
                {
                    if (indep[t])
                    {
                        slong sidx = src[t];
                        ni[cnt] = (sidx < r_old) ? ri[sidx] : mi[sidx];
                        ne[cnt] = (sidx < r_old) ? re[sidx] : me[sidx];
                        cnt++;
                    }
                }
                for (t = 0; t < r; t++)
                {
                    ri[t] = ni[t];
                    re[t] = ne[t];
                }
                flint_free(ni);
            }
        }

        nmod_mat_clear(Mrg);
        nmod_mat_clear(New);

        if (2 * P.q > maxbeta)
            break;

        _krylov_pow_square(&P);
    }

    _krylov_pow_clear(&P);

    /* output */
    for (t = 0; t < r; t++)
    {
        memcpy(nmod_mat_entry_ptr(B, t, 0), nmod_mat_entry_ptr(Cur, t, 0),
               D * sizeof(ulong));
        rows[t] = ri[t];
        exps[t] = re[t];
    }

    if (ok != NULL)
        *ok = 1;

    for (i = 0; i < m; i++)
        order[i] = 0;
    for (t = 0; t < r; t++)
        order[ri[t]] = FLINT_MAX(order[ri[t]], re[t] + 1);
    for (i = 0; i < m; i++)
    {
        if (deltas != NULL)
            deltas[i] = order[i];
        if (order[i] > beta[i] && ok != NULL)
            *ok = 0;
    }

    nmod_mat_clear(Cur);
    flint_free(beta);
    flint_free(ri);
    flint_free(indep);

    return r;
}
