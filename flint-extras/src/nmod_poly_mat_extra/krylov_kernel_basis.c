/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/nmod.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_poly_mat.h>
#include "nmod_mat_extra/krylov_impl.h"
#include "nmod_poly_mat_relation_basis.h"

/*
    P = compression of [-R | I] (JNSV 2017, Lemma 7.11): row i of P is
    X^{delta_i} e_i - sum_k R[i, k] X^{exps[k]} e_{rows[k]}. The rows
    (c, e) of the profile are exactly those with e < delta_c, so that
    column c of P has degree less than delta_c except for its diagonal
    entry, which is monic of degree delta_c.
*/
static void
_nmod_poly_mat_set_krylov_kernel(nmod_poly_mat_t P, const nmod_mat_t R,
        const slong * rows, const slong * exps, const slong * deltas, slong r)
{
    slong m = P->r, i, c, e, k, len;
    slong * off, * kidx;
    nmod_t mod = R->mod;

    off = flint_malloc((m + 1 + r + 1) * sizeof(slong));
    kidx = off + m + 1;

    off[0] = 0;
    for (c = 0; c < m; c++)
        off[c + 1] = off[c] + deltas[c];

    for (k = 0; k < r; k++)
        kidx[off[rows[k]] + exps[k]] = k;

    for (i = 0; i < m; i++)
    {
        nn_srcptr Ri = (r > 0) ? nmod_mat_entry_ptr(R, i, 0) : NULL;

        for (c = 0; c < m; c++)
        {
            nmod_poly_struct * p = nmod_poly_mat_entry(P, i, c);

            len = deltas[c] + (i == c);
            nmod_poly_fit_length(p, len);

            for (e = 0; e < deltas[c]; e++)
                p->coeffs[e] = nmod_neg(Ri[kidx[off[c] + e]], mod);

            if (i == c)
                p->coeffs[deltas[c]] = 1;

            _nmod_poly_set_length(p, len);
            _nmod_poly_normalise(p);
        }
    }

    flint_free(off);
}

static int
_nmod_poly_mat_krylov_kernel_basis(nmod_poly_mat_t P, const nmod_mat_t E,
        const _nmod_mat_krylov_op_t op, const slong * shift,
        const slong * degrees)
{
    slong m = E->r, n = E->c, r;
    slong * rows, * exps, * deltas;
    nmod_mat_t R;

    if (P->r != m || P->c != m)
        flint_throw(FLINT_ERROR, "nmod_poly_mat_krylov_kernel_basis: P must be E->r x E->r\n");

    rows = flint_malloc((2 * n + m + 1) * sizeof(slong));
    exps = rows + n;
    deltas = exps + n;

    nmod_mat_init(R, m, n, E->mod.n);

    r = _nmod_mat_krylov_kernel(R, rows, exps, deltas, E, op, shift, degrees, 0,
                                _NMOD_MAT_KRYLOV_KERNEL_AUTO);

    if (r >= 0)
        _nmod_poly_mat_set_krylov_kernel(P, R, rows, exps, deltas, r);

    nmod_mat_clear(R);
    flint_free(rows);

    return (r >= 0);
}

int
nmod_poly_mat_krylov_kernel_basis(nmod_poly_mat_t P, const nmod_mat_t E,
        const nmod_mat_t M, const slong * shift, const slong * degrees)
{
    _nmod_mat_krylov_op_t op;
    int res;

    if (M->r != E->c || M->c != E->c)
        flint_throw(FLINT_ERROR, "nmod_poly_mat_krylov_kernel_basis: incompatible dimensions\n");

    _nmod_mat_krylov_op_init_dense(op, M);
    res = _nmod_poly_mat_krylov_kernel_basis(P, E, op, shift, degrees);
    _nmod_mat_krylov_op_clear(op);

    return res;
}

int
nmod_poly_mat_krylov_kernel_basis_jordan(nmod_poly_mat_t P,
        const nmod_mat_t E, nn_srcptr eig, const slong * sizes, slong nblocks,
        const slong * shift, const slong * degrees)
{
    _nmod_mat_krylov_op_t op;
    int res;

    _nmod_mat_krylov_op_init_jordan(op, eig, sizes, nblocks, E->mod);

    if (op->n != E->c)
        flint_throw(FLINT_ERROR, "nmod_poly_mat_krylov_kernel_basis_jordan: incompatible dimensions\n");

    res = _nmod_poly_mat_krylov_kernel_basis(P, E, op, shift, degrees);
    _nmod_mat_krylov_op_clear(op);

    return res;
}
