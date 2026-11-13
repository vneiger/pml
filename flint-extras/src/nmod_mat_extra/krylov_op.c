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
#include <flint/ulong_extras.h>
#include <flint/nmod.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_mat.h>
#include "nmod_mat_extra/krylov_impl.h"

void
_nmod_mat_krylov_op_init_dense(_nmod_mat_krylov_op_t op, const nmod_mat_t M)
{
    if (M->r != M->c)
        flint_throw(FLINT_ERROR, "nmod_mat_krylov: M must be square\n");

    op->kind = _NMOD_MAT_KRYLOV_DENSE;
    op->n = M->r;
    op->mod = M->mod;
    op->M = M;
    op->nblocks = 0;
    op->eig = NULL;
    op->eig_pre = NULL;
    op->off = NULL;
}

void
_nmod_mat_krylov_op_init_jordan(_nmod_mat_krylov_op_t op, nn_srcptr eig,
        const slong * sizes, slong nblocks, nmod_t mod)
{
    slong b;
    int all_zero = 1, all_one = 1;

    if (nblocks < 0)
        flint_throw(FLINT_ERROR,
            "nmod_mat_krylov: the number of Jordan blocks must be nonnegative\n");

    op->mod = mod;
    op->M = NULL;
    op->nblocks = nblocks;
    op->off = flint_malloc((nblocks + 1) * sizeof(slong));
    op->off[0] = 0;

    for (b = 0; b < nblocks; b++)
    {
        slong t = (sizes == NULL) ? 1 : sizes[b];

        if (t <= 0)
            flint_throw(FLINT_ERROR,
                "nmod_mat_krylov: Jordan block sizes must be positive\n");

        op->off[b + 1] = op->off[b] + t;
        if (t != 1)
            all_one = 0;
    }

    op->n = op->off[nblocks];

    op->eig = NULL;
    op->eig_pre = NULL;

    if (eig != NULL && nblocks > 0)
    {
        op->eig = _nmod_vec_init(nblocks);
        for (b = 0; b < nblocks; b++)
        {
            op->eig[b] = (eig[b] < mod.n) ? eig[b] : n_mod2_preinv(eig[b], mod.n, mod.ninv);
            if (op->eig[b] != 0)
                all_zero = 0;
        }
    }

    if (all_zero)
    {
        op->kind = _NMOD_MAT_KRYLOV_SHIFT;
        if (op->eig != NULL)
        {
            _nmod_vec_clear(op->eig);
            op->eig = NULL;
        }
        return;
    }

    op->kind = all_one ? _NMOD_MAT_KRYLOV_DIAGONAL : _NMOD_MAT_KRYLOV_JORDAN;

    /* Shoup's multiplication needs n < 2^(FLINT_BITS - 1) */
    if (mod.n < (UWORD(1) << (FLINT_BITS - 1)))
    {
        op->eig_pre = _nmod_vec_init(nblocks);
        for (b = 0; b < nblocks; b++)
            op->eig_pre[b] = n_mulmod_precomp_shoup(op->eig[b], mod.n);
    }
}

void
_nmod_mat_krylov_op_clear(_nmod_mat_krylov_op_t op)
{
    if (op->eig != NULL)
        _nmod_vec_clear(op->eig);
    if (op->eig_pre != NULL)
        _nmod_vec_clear(op->eig_pre);
    if (op->off != NULL)
        flint_free(op->off);
}

void
_nmod_mat_krylov_op_apply_vec(nn_ptr out, nn_srcptr in,
        const _nmod_mat_krylov_op_t op)
{
    slong b, c, nb = op->nblocks;
    const slong * off = op->off;
    nmod_t mod = op->mod;

    switch (op->kind)
    {
        case _NMOD_MAT_KRYLOV_SHIFT:
            /* within each block, shift the entries by one position */
            for (b = 0; b < nb; b++)
            {
                slong o = off[b], t = off[b + 1] - off[b];
                out[o] = 0;
                if (t > 1)
                    memcpy(out + o + 1, in + o, (t - 1) * sizeof(ulong));
            }
            break;

        case _NMOD_MAT_KRYLOV_DIAGONAL:
            if (op->eig_pre != NULL)
            {
                for (c = 0; c < nb; c++)
                    out[c] = n_mulmod_shoup(op->eig[c], in[c], op->eig_pre[c], mod.n);
            }
            else
            {
                for (c = 0; c < nb; c++)
                    out[c] = nmod_mul(op->eig[c], in[c], mod);
            }
            break;

        case _NMOD_MAT_KRYLOV_JORDAN:
            for (b = 0; b < nb; b++)
            {
                slong o = off[b], t = off[b + 1] - off[b];
                ulong x = op->eig[b];

                if (x == 0)
                {
                    out[o] = 0;
                    if (t > 1)
                        memcpy(out + o + 1, in + o, (t - 1) * sizeof(ulong));
                }
                else if (op->eig_pre != NULL)
                {
                    ulong xp = op->eig_pre[b];
                    out[o] = n_mulmod_shoup(x, in[o], xp, mod.n);
                    for (c = 1; c < t; c++)
                        out[o + c] = nmod_add(n_mulmod_shoup(x, in[o + c], xp, mod.n),
                                              in[o + c - 1], mod);
                }
                else
                {
                    out[o] = nmod_mul(x, in[o], mod);
                    for (c = 1; c < t; c++)
                        out[o + c] = nmod_add(nmod_mul(x, in[o + c], mod),
                                              in[o + c - 1], mod);
                }
            }
            break;

        default:
            flint_throw(FLINT_ERROR, "_nmod_mat_krylov_op_apply_vec: dense operator\n");
    }
}
