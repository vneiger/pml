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
    compared with D. Measured with uniform shift, Hermite-Pade (one
    nilpotent block) and dense M, D = 64 to 1024 (Krylov tuning kit, Nov
    2026, Zen 4, Ice Lake, Cascade Lake, Apple M4), LU relative to TRACK:
    0.32 to 0.97 for m <= D / 2; for m = D, 0.77 to 1.04 below 2^32, 0.60
    to 0.85 at 60 bits, and at 50 bits 0.91 to 1.21 for D <= 256 and 0.77
    to 1.00 for D >= 512 on x86, 0.69 to 0.83 on the M4; for m = 2 D, 0.94
    to 1.59 below 2^52 except 0.88 to 0.99 on the M4 at 50 bits, and 0.70
    to 1.03 at 60 bits; for m = 4 D, 1.06 to 2.13 below 2^52 and 0.79 to
    1.31 at 60 bits. The aarch64 rule between 2^32 and 2^52 comes from the
    M4 alone.
*/
static int
_krylov_kernel_lu_wanted(slong m, slong D, nmod_t mod)
{
    flint_bitcnt_t bits = NMOD_BITS(mod);

    if (2 * m <= D)
        return 1;
    if (bits > 52)
        return m <= 2 * D;
    if (bits <= 32)
        return m <= D;
#if defined(__aarch64__)
    return m <= 2 * D;
#else
    return m <= D && D >= 512;
#endif
}

/*
    With B the Krylov basis (r x D, full rank), piv r columns such that
    C = B[:, piv] is invertible, and T the target rows (m x D), which are
    in the span of B: T = T[:, piv] C^-1 B, so that the relations are
    R = T[:, piv] C^-1 (JNSV 2017, Lemma 7.11). Either G = C^-1 is obtained
    along the elimination (TRACK, reduced row echelon basis), and
    R = T[:, piv] G, or the elimination gives a factorization C = Pi^T L U
    (LU, row echelon basis), and R Pi^T = T[:, piv] U^-1 L^-1 is obtained by
    two triangular solves with the matrix on the right.
*/
slong
_nmod_mat_krylov_kernel(nmod_mat_t R, slong * rows, slong * exps,
        slong * deltas, const nmod_mat_t V, const _nmod_mat_krylov_op_t op,
        const slong * shift, const slong * degrees, slong chunk, int method)
{
    slong m = V->r, D = V->c, r, i, t;
    slong * piv, * bpos = NULL;
    int ok;
    nmod_mat_t B, T, G;

    if (R->r != m)
        flint_throw(FLINT_ERROR, "nmod_mat_krylov_kernel: R must have V->r rows\n");

    if (method == _NMOD_MAT_KRYLOV_KERNEL_AUTO)
        method = _krylov_kernel_lu_wanted(m, D, V->mod)
                    ? _NMOD_MAT_KRYLOV_KERNEL_LU : _NMOD_MAT_KRYLOV_KERNEL_TRACK;

    nmod_mat_init(B, D, D, V->mod.n);
    nmod_mat_init(T, m, D, V->mod.n);
    nmod_mat_init(G, D, D, V->mod.n);   /* G, or L and U */
    piv = flint_malloc((D + 1) * sizeof(slong));

    if (method == _NMOD_MAT_KRYLOV_KERNEL_LU)
    {
        bpos = flint_malloc((D + 1) * sizeof(slong));
        r = _nmod_mat_krylov_basis_iter_lu(B, rows, exps, deltas, T, G, bpos,
                                           piv, &ok, V, op, shift, degrees, chunk);
    }
    else
    {
        r = _nmod_mat_krylov_basis_iter(B, rows, exps, deltas, T, G, piv, &ok,
                                        V, op, shift, degrees, chunk);
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

        nmod_mat_init(Tp, m, r, V->mod.n);
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
        slong * deltas, const nmod_mat_t V, const nmod_mat_t M,
        const slong * shift, const slong * degrees)
{
    _nmod_mat_krylov_op_t op;
    slong r;

    if (M->r != V->c || M->c != V->c)
        flint_throw(FLINT_ERROR, "nmod_mat_krylov_kernel: incompatible dimensions\n");

    _nmod_mat_krylov_op_init_dense(op, M);
    r = _nmod_mat_krylov_kernel(R, rows, exps, deltas, V, op, shift, degrees, 0,
                                _NMOD_MAT_KRYLOV_KERNEL_AUTO);
    _nmod_mat_krylov_op_clear(op);

    return r;
}

slong
nmod_mat_krylov_kernel_jordan(nmod_mat_t R, slong * rows, slong * exps,
        slong * deltas, const nmod_mat_t V,
        nn_srcptr eig, const slong * sizes, slong nblocks,
        const slong * shift, const slong * degrees)
{
    _nmod_mat_krylov_op_t op;
    slong r;

    _nmod_mat_krylov_op_init_jordan(op, eig, sizes, nblocks, V->mod);

    if (op->D != V->c)
        flint_throw(FLINT_ERROR, "nmod_mat_krylov_kernel_jordan: incompatible dimensions\n");

    r = _nmod_mat_krylov_kernel(R, rows, exps, deltas, V, op, shift, degrees, 0,
                                _NMOD_MAT_KRYLOV_KERNEL_AUTO);
    _nmod_mat_krylov_op_clear(op);

    return r;
}
