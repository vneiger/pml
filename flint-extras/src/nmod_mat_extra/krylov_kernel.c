/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/nmod_mat.h>
#include "nmod_mat_extra/krylov_impl.h"

/*
    With B the Krylov basis (r x n, full rank), W = G B its reduced row
    echelon form with pivot columns piv, and T the target rows (m x n),
    which are in the span of B: T = T[:, piv] W = T[:, piv] G B, so that
    the relations are R = T[:, piv] G (this is JNSV 2017, Lemma 7.11, with
    G = B[:, piv]^(-1) obtained along the elimination).
*/
slong
_nmod_mat_krylov_kernel(nmod_mat_t R, slong * rows, slong * exps,
        slong * deltas, const nmod_mat_t E, const _nmod_mat_krylov_op_t op,
        const slong * shifts, const slong * degrees, slong chunk)
{
    slong m = E->r, n = E->c, r, i, t;
    slong * piv;
    int ok;
    nmod_mat_t B, T, G, Tp, Gw, Rw;

    if (R->r != m)
        flint_throw(FLINT_ERROR, "nmod_mat_krylov_kernel: R must have E->r rows\n");

    nmod_mat_init(B, n, n, E->mod.n);
    nmod_mat_init(T, m, n, E->mod.n);
    nmod_mat_init(G, n, n, E->mod.n);
    piv = flint_malloc((n + 1) * sizeof(slong));

    r = _nmod_mat_krylov_basis_iter(B, rows, exps, deltas, T, G, piv, &ok,
                                    E, op, shifts, degrees, chunk);

    if (!ok)
    {
        r = -1;
    }
    else if (r > 0)
    {
        if (R->c < r)
            flint_throw(FLINT_ERROR, "nmod_mat_krylov_kernel: R has too few columns\n");

        nmod_mat_init(Tp, m, r, E->mod.n);
        for (i = 0; i < m; i++)
        {
            nn_srcptr y = nmod_mat_entry_ptr(T, i, 0);
            nn_ptr z = nmod_mat_entry_ptr(Tp, i, 0);
            for (t = 0; t < r; t++)
                z[t] = y[piv[t]];
        }

        nmod_mat_window_init(Gw, G, 0, 0, r, r);
        nmod_mat_window_init(Rw, R, 0, 0, m, r);
        nmod_mat_mul(Rw, Tp, Gw);
        nmod_mat_window_clear(Gw);
        nmod_mat_window_clear(Rw);
        nmod_mat_clear(Tp);
    }

    nmod_mat_clear(B);
    nmod_mat_clear(T);
    nmod_mat_clear(G);
    flint_free(piv);

    return r;
}

slong
nmod_mat_krylov_kernel(nmod_mat_t R, slong * rows, slong * exps,
        slong * deltas, const nmod_mat_t E, const nmod_mat_t M,
        const slong * shifts, const slong * degrees)
{
    _nmod_mat_krylov_op_t op;
    slong r;

    if (M->r != E->c || M->c != E->c)
        flint_throw(FLINT_ERROR, "nmod_mat_krylov_kernel: incompatible dimensions\n");

    _nmod_mat_krylov_op_init_dense(op, M);
    r = _nmod_mat_krylov_kernel(R, rows, exps, deltas, E, op, shifts, degrees, 0);
    _nmod_mat_krylov_op_clear(op);

    return r;
}

slong
nmod_mat_krylov_kernel_jordan(nmod_mat_t R, slong * rows, slong * exps,
        slong * deltas, const nmod_mat_t E,
        nn_srcptr eig, const slong * sizes, slong nblocks,
        const slong * shifts, const slong * degrees)
{
    _nmod_mat_krylov_op_t op;
    slong r;

    _nmod_mat_krylov_op_init_jordan(op, eig, sizes, nblocks, E->mod);

    if (op->n != E->c)
        flint_throw(FLINT_ERROR, "nmod_mat_krylov_kernel_jordan: incompatible dimensions\n");

    r = _nmod_mat_krylov_kernel(R, rows, exps, deltas, E, op, shifts, degrees, 0);
    _nmod_mat_krylov_op_clear(op);

    return r;
}
