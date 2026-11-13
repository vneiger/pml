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
#include <flint/nmod_mat.h>
#include "nmod_mat_extra/krylov_impl.h"

/*
    Automatic choice of the method (see krylov_impl.h). The row echelon
    basis (LU) saves the transformation G, whose cost is about r^3
    whatever m is, but its final step is two triangular solves with m
    right-hand sides instead of one product: it pays off unless m is large
    compared with n, all the more as the modulus is large (FLINT's
    triangular solves are closest to the speed of nmod_mat_mul above
    2^52). Measured with uniform shifts, Hermite-Pade and dense M, n = 64
    to 1024, x86_64 AVX-512 without IFMA (relative to TRACK): 0.31 to 0.98
    for m <= n / 2 at 30 to 60 bits; for m = n, 0.79 to 0.83 at n = 1024
    and 0.94 to 1.19 below at 30 to 50 bits, 0.71 to 0.82 at 60 bits; for
    m = 2 n, 1.02 to 1.38 at 30 to 50 bits and 0.84 to 0.98 at 60 bits;
    for m = 4 n, 1.11 to 1.94 at 30 to 50 bits and 0.94 to 1.24 at 60 bits.
*/
static int
_krylov_kernel_lu_wanted(slong m, slong n, nmod_t mod)
{
    if (2 * m <= n)
        return 1;
    if (NMOD_BITS(mod) > 52)
        return m <= 2 * n;
    return m <= n && n >= 512;
}

/*
    With B the Krylov basis (r x n, full rank), piv r columns such that
    C = B[:, piv] is invertible, and T the target rows (m x n), which are
    in the span of B: T = T[:, piv] C^-1 B, so that the relations are
    R = T[:, piv] C^-1 (JNSV 2017, Lemma 7.11). Either G = C^-1 is obtained
    along the elimination (TRACK, reduced row echelon basis), and
    R = T[:, piv] G, or the elimination gives a factorization C = Pi^T L U
    (LU, row echelon basis), and R Pi^T = T[:, piv] U^-1 L^-1 is obtained by
    two triangular solves with the matrix on the right.
*/
slong
_nmod_mat_krylov_kernel(nmod_mat_t R, slong * rows, slong * exps,
        slong * deltas, const nmod_mat_t E, const _nmod_mat_krylov_op_t op,
        const slong * shifts, const slong * degrees, slong chunk, int method)
{
    slong m = E->r, n = E->c, r, i, t;
    slong * piv, * bpos = NULL;
    int ok;
    nmod_mat_t B, T, G;

    if (R->r != m)
        flint_throw(FLINT_ERROR, "nmod_mat_krylov_kernel: R must have E->r rows\n");

    if (method == _NMOD_MAT_KRYLOV_KERNEL_AUTO)
        method = _krylov_kernel_lu_wanted(m, n, E->mod)
                    ? _NMOD_MAT_KRYLOV_KERNEL_LU : _NMOD_MAT_KRYLOV_KERNEL_TRACK;

    nmod_mat_init(B, n, n, E->mod.n);
    nmod_mat_init(T, m, n, E->mod.n);
    nmod_mat_init(G, n, n, E->mod.n);   /* G, or L and U */
    piv = flint_malloc((n + 1) * sizeof(slong));

    if (method == _NMOD_MAT_KRYLOV_KERNEL_LU)
    {
        bpos = flint_malloc((n + 1) * sizeof(slong));
        r = _nmod_mat_krylov_basis_iter_lu(B, rows, exps, deltas, T, G, bpos,
                                           piv, &ok, E, op, shifts, degrees, chunk);
    }
    else
    {
        r = _nmod_mat_krylov_basis_iter(B, rows, exps, deltas, T, G, piv, &ok,
                                        E, op, shifts, degrees, chunk);
    }

    if (!ok)
    {
        r = -1;
    }
    else if (r > 0)
    {
        nmod_mat_t Tp, Gw;

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

        if (method == _NMOD_MAT_KRYLOV_KERNEL_LU)
        {
            /* x[bpos[t]] = (T[:, piv] U^-1 L^-1)[t] */
            nmod_mat_solve_triu_right(Tp, Gw, Tp, 0);
            nmod_mat_solve_tril_right(Tp, Gw, Tp, 1);

            for (i = 0; i < m; i++)
            {
                nn_srcptr y = nmod_mat_entry_ptr(Tp, i, 0);
                nn_ptr z = nmod_mat_entry_ptr(R, i, 0);
                for (t = 0; t < r; t++)
                    z[bpos[t]] = y[t];
            }
        }
        else
        {
            nmod_mat_t Rw;
            nmod_mat_window_init(Rw, R, 0, 0, m, r);
            nmod_mat_mul(Rw, Tp, Gw);
            nmod_mat_window_clear(Rw);
        }

        nmod_mat_window_clear(Gw);
        nmod_mat_clear(Tp);
    }

    nmod_mat_clear(G);
    nmod_mat_clear(B);
    nmod_mat_clear(T);
    flint_free(piv);
    if (bpos != NULL)
        flint_free(bpos);

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
    r = _nmod_mat_krylov_kernel(R, rows, exps, deltas, E, op, shifts, degrees, 0,
                                _NMOD_MAT_KRYLOV_KERNEL_AUTO);
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

    r = _nmod_mat_krylov_kernel(R, rows, exps, deltas, E, op, shifts, degrees, 0,
                                _NMOD_MAT_KRYLOV_KERNEL_AUTO);
    _nmod_mat_krylov_op_clear(op);

    return r;
}
