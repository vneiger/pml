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
        const nmod_mat_t E, const _nmod_mat_krylov_op_t op,
        const slong * shift, const slong * degrees)
{
    slong m = E->r, n = E->c, N, i, j, k, maxd;
    slong * ri, * rj, * off, * pos;

    if (op->n != n)
        flint_throw(FLINT_ERROR, "nmod_mat_krylov_matrix: incompatible dimensions\n");

    N = _nmod_mat_krylov_nrows(m, degrees);

    if (K->r != N || K->c != n)
        flint_throw(FLINT_ERROR, "nmod_mat_krylov_matrix: K must have sum (degrees[i] + 1) rows and E->c columns\n");

    if (N == 0)
        return;

    ri = flint_malloc((3 * N + m + 1) * sizeof(slong));
    rj = ri + N;
    pos = rj + N;
    off = pos + N;

    _nmod_mat_krylov_order(ri, rj, m, shift, degrees);

    /* pos[off[i] + j] = index in K of the row (i, j) */
    off[0] = 0;
    maxd = 0;
    for (i = 0; i < m; i++)
    {
        off[i + 1] = off[i] + degrees[i] + 1;
        maxd = FLINT_MAX(maxd, degrees[i]);
    }
    for (k = 0; k < N; k++)
        pos[off[ri[k]] + rj[k]] = k;

    if (n > 0)
    {
        if (op->kind == _NMOD_MAT_KRYLOV_DENSE)
        {
            nmod_mat_t S, O;
            slong c, u, * act;

            for (i = 0; i < m; i++)
                memcpy(nmod_mat_entry_ptr(K, pos[off[i]], 0),
                       nmod_mat_entry_ptr(E, i, 0), n * sizeof(ulong));

            /* indices i with degrees[i] >= j, compacted at each step */
            act = flint_malloc((m + 1) * sizeof(slong));
            c = 0;
            for (i = 0; i < m; i++)
                if (degrees[i] >= 1)
                    act[c++] = i;

            /* all rows of exponent j at once, from those of exponent j-1 */
            for (j = 1; j <= maxd && c > 0; j++)
            {
                slong c2 = 0;
                for (u = 0; u < c; u++)
                    if (degrees[act[u]] >= j)
                        act[c2++] = act[u];
                c = c2;

                nmod_mat_init(S, c, n, E->mod.n);
                nmod_mat_init(O, c, n, E->mod.n);

                for (u = 0; u < c; u++)
                    memcpy(nmod_mat_entry_ptr(S, u, 0),
                           nmod_mat_entry_ptr(K, pos[off[act[u]] + j - 1], 0),
                           n * sizeof(ulong));

                nmod_mat_mul(O, S, op->M);

                for (u = 0; u < c; u++)
                    memcpy(nmod_mat_entry_ptr(K, pos[off[act[u]] + j], 0),
                           nmod_mat_entry_ptr(O, u, 0), n * sizeof(ulong));

                nmod_mat_clear(S);
                nmod_mat_clear(O);
            }

            flint_free(act);
        }
        else
        {
            /* the row (i, j-1) precedes the row (i, j) */
            for (k = 0; k < N; k++)
            {
                i = ri[k];
                j = rj[k];
                if (j == 0)
                    memcpy(nmod_mat_entry_ptr(K, k, 0),
                           nmod_mat_entry_ptr(E, i, 0), n * sizeof(ulong));
                else
                    _nmod_mat_krylov_op_apply_vec(nmod_mat_entry_ptr(K, k, 0),
                        nmod_mat_entry_ptr(K, pos[off[i] + j - 1], 0), op);
            }
        }
    }

    if (rows != NULL)
        memcpy(rows, ri, N * sizeof(slong));
    if (exps != NULL)
        memcpy(exps, rj, N * sizeof(slong));

    flint_free(ri);
}

void
nmod_mat_krylov_matrix(nmod_mat_t K, const nmod_mat_t E, const nmod_mat_t M,
        const slong * shift, const slong * degrees)
{
    _nmod_mat_krylov_op_t op;
    slong * d;

    if (M->r != E->c || M->c != E->c)
        flint_throw(FLINT_ERROR, "nmod_mat_krylov_matrix: incompatible dimensions\n");

    d = flint_malloc((E->r + 1) * sizeof(slong));
    _nmod_mat_krylov_degrees(d, E->r, degrees, E->c, 0);
    _nmod_mat_krylov_op_init_dense(op, M);
    _nmod_mat_krylov_matrix(K, NULL, NULL, E, op, shift, d);
    _nmod_mat_krylov_op_clear(op);
    flint_free(d);
}

void
nmod_mat_krylov_matrix_jordan(nmod_mat_t K, const nmod_mat_t E,
        nn_srcptr eig, const slong * sizes, slong nblocks,
        const slong * shift, const slong * degrees)
{
    _nmod_mat_krylov_op_t op;
    slong * d;

    _nmod_mat_krylov_op_init_jordan(op, eig, sizes, nblocks, E->mod);

    if (op->n != E->c)
        flint_throw(FLINT_ERROR, "nmod_mat_krylov_matrix_jordan: incompatible dimensions\n");

    d = flint_malloc((E->r + 1) * sizeof(slong));
    _nmod_mat_krylov_degrees(d, E->r, degrees, E->c, 0);
    _nmod_mat_krylov_matrix(K, NULL, NULL, E, op, shift, d);
    _nmod_mat_krylov_op_clear(op);
    flint_free(d);
}
