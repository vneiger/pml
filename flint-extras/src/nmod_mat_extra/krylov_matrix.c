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
#include <flint/nmod_mat.h>
#include "nmod_mat_extra/krylov_impl.h"

void
_nmod_mat_krylov_matrix(nmod_mat_t K, slong * rows, slong * exps,
        const nmod_mat_t V, const _nmod_mat_krylov_op_t op,
        const slong * shift, const slong * degrees)
{
    slong m = V->r, D = V->c, N, i, e, k, maxbeta;
    slong * ri, * re, * off, * pos;

    if (op->D != D)
        flint_throw(FLINT_ERROR, "nmod_mat_krylov_matrix: incompatible dimensions\n");

    N = _nmod_mat_krylov_nrows(m, degrees);

    if (K->r != N || K->c != D)
        flint_throw(FLINT_ERROR, "nmod_mat_krylov_matrix: K must have sum (degrees[i] + 1) rows and V->c columns\n");

    if (N == 0)
        return;

    ri = flint_malloc((3 * N + m + 1) * sizeof(slong));
    re = ri + N;
    pos = re + N;
    off = pos + N;

    _nmod_mat_krylov_order(ri, re, m, shift, degrees);

    /* pos[off[i] + e] = index in K of the row (i, e) */
    off[0] = 0;
    maxbeta = 0;
    for (i = 0; i < m; i++)
    {
        off[i + 1] = off[i] + degrees[i] + 1;
        maxbeta = FLINT_MAX(maxbeta, degrees[i]);
    }
    for (k = 0; k < N; k++)
        pos[off[ri[k]] + re[k]] = k;

    if (D > 0)
    {
        if (op->kind == _NMOD_MAT_KRYLOV_DENSE)
        {
            nmod_mat_t S, O;
            slong c, u, * act;

            for (i = 0; i < m; i++)
                memcpy(nmod_mat_entry_ptr(K, pos[off[i]], 0),
                       nmod_mat_entry_ptr(V, i, 0), D * sizeof(ulong));

            /* indices i with degrees[i] >= e, compacted at each step */
            act = flint_malloc((m + 1) * sizeof(slong));
            c = 0;
            for (i = 0; i < m; i++)
                if (degrees[i] >= 1)
                    act[c++] = i;

            /* all rows of exponent e at once, from those of exponent e-1 */
            for (e = 1; e <= maxbeta && c > 0; e++)
            {
                slong c2 = 0;
                for (u = 0; u < c; u++)
                    if (degrees[act[u]] >= e)
                        act[c2++] = act[u];
                c = c2;

                nmod_mat_init(S, c, D, V->mod.n);
                nmod_mat_init(O, c, D, V->mod.n);

                for (u = 0; u < c; u++)
                    memcpy(nmod_mat_entry_ptr(S, u, 0),
                           nmod_mat_entry_ptr(K, pos[off[act[u]] + e - 1], 0),
                           D * sizeof(ulong));

                nmod_mat_mul(O, S, op->M);

                for (u = 0; u < c; u++)
                    memcpy(nmod_mat_entry_ptr(K, pos[off[act[u]] + e], 0),
                           nmod_mat_entry_ptr(O, u, 0), D * sizeof(ulong));

                nmod_mat_clear(S);
                nmod_mat_clear(O);
            }

            flint_free(act);
        }
        else
        {
            /* the row (i, e-1) precedes the row (i, e) */
            for (k = 0; k < N; k++)
            {
                i = ri[k];
                e = re[k];
                if (e == 0)
                    memcpy(nmod_mat_entry_ptr(K, k, 0),
                           nmod_mat_entry_ptr(V, i, 0), D * sizeof(ulong));
                else
                    _nmod_mat_krylov_op_apply_vec(nmod_mat_entry_ptr(K, k, 0),
                        nmod_mat_entry_ptr(K, pos[off[i] + e - 1], 0), op);
            }
        }
    }

    if (rows != NULL)
        memcpy(rows, ri, N * sizeof(slong));
    if (exps != NULL)
        memcpy(exps, re, N * sizeof(slong));

    flint_free(ri);
}

void
nmod_mat_krylov_matrix(nmod_mat_t K, const nmod_mat_t V, const nmod_mat_t M,
        const slong * shift, const slong * degrees)
{
    _nmod_mat_krylov_op_t op;
    slong * beta;

    if (M->r != V->c || M->c != V->c)
        flint_throw(FLINT_ERROR, "nmod_mat_krylov_matrix: incompatible dimensions\n");

    beta = flint_malloc((V->r + 1) * sizeof(slong));
    _nmod_mat_krylov_degrees(beta, V->r, degrees, V->c, 0);
    _nmod_mat_krylov_op_init_dense(op, M);
    _nmod_mat_krylov_matrix(K, NULL, NULL, V, op, shift, beta);
    _nmod_mat_krylov_op_clear(op);
    flint_free(beta);
}

void
nmod_mat_krylov_matrix_jordan(nmod_mat_t K, const nmod_mat_t V,
        nn_srcptr eig, const slong * sizes, slong nblocks,
        const slong * shift, const slong * degrees)
{
    _nmod_mat_krylov_op_t op;
    slong * beta;

    _nmod_mat_krylov_op_init_jordan(op, eig, sizes, nblocks, V->mod);

    if (op->D != V->c)
        flint_throw(FLINT_ERROR, "nmod_mat_krylov_matrix_jordan: incompatible dimensions\n");

    beta = flint_malloc((V->r + 1) * sizeof(slong));
    _nmod_mat_krylov_degrees(beta, V->r, degrees, V->c, 0);
    _nmod_mat_krylov_matrix(K, NULL, NULL, V, op, shift, beta);
    _nmod_mat_krylov_op_clear(op);
    flint_free(beta);
}
