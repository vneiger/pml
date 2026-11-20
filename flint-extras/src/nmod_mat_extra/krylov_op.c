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
    op->D = M->r;
    op->mod = M->mod;
    op->M = M;
    op->logpow = -1;
    op->powruns = 0;
    op->mpdeg = M->r;
    op->predict = 1;
    op->nblocks = 0;
    op->eig = NULL;
    op->off = NULL;
}

/* (eigenvalue, block size) */
typedef struct
{
    ulong e;
    slong t;
}
_krylov_block_struct;

static int
_krylov_block_cmp(const void * a, const void * b, void * data)
{
    ulong x = ((const _krylov_block_struct *) a)->e;
    ulong y = ((const _krylov_block_struct *) b)->e;

    (void) data;
    return (x > y) - (x < y);
}

/* degree of the minimal polynomial of the Jordan matrix: the sum, over
   the distinct eigenvalues, of their largest block size (eig == NULL: all
   eigenvalues zero) */
static slong
_krylov_jordan_minpoly_degree(nn_srcptr eig, const slong * off, slong nblocks)
{
    slong b, deg = 0, tmax;
    _krylov_block_struct * blk;

    if (nblocks == 0)
        return 0;

    if (eig == NULL)
    {
        for (b = 0; b < nblocks; b++)
            deg = FLINT_MAX(deg, off[b + 1] - off[b]);
        return deg;
    }

    blk = flint_malloc(nblocks * sizeof(_krylov_block_struct));
    for (b = 0; b < nblocks; b++)
    {
        blk[b].e = eig[b];
        blk[b].t = off[b + 1] - off[b];
    }
    flint_sort(blk, nblocks, sizeof(_krylov_block_struct),
               _krylov_block_cmp, NULL);

    tmax = blk[0].t;
    for (b = 1; b < nblocks; b++)
    {
        if (blk[b].e != blk[b - 1].e)
        {
            deg += tmax;
            tmax = 0;
        }
        tmax = FLINT_MAX(tmax, blk[b].t);
    }
    deg += tmax;

    flint_free(blk);
    return deg;
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
    op->logpow = 0;
    op->powruns = 0;
    op->predict = 1;
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

    op->D = op->off[nblocks];

    op->eig = NULL;

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

    op->mpdeg = _krylov_jordan_minpoly_degree(op->eig, op->off, nblocks);

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
}

void
_nmod_mat_krylov_op_clear(_nmod_mat_krylov_op_t op)
{
    if (op->eig != NULL)
        _nmod_vec_clear(op->eig);
    if (op->off != NULL)
        flint_free(op->off);
}

void
_nmod_mat_krylov_op_apply_vec(nn_ptr out, nn_srcptr in,
        const _nmod_mat_krylov_op_t op)
{
    slong b, nb = op->nblocks;
    const slong * off = op->off;

    switch (op->kind)
    {
        case _NMOD_MAT_KRYLOV_SHIFT:
        case _NMOD_MAT_KRYLOV_JORDAN:
            /* within each block, shift the entries by one position:
               -> one shifted copy of the whole row, then the first entry of each
               block is set to zero (measured faster than one copy per block
               when the blocks are small) */
            if (op->D > 1)
                memcpy(out + 1, in, (op->D - 1) * sizeof(ulong));
            for (b = 0; b < nb; b++)
                out[off[b]] = 0;
            /* Jordan: the block times its eigenvalue added */
            if (op->kind == _NMOD_MAT_KRYLOV_JORDAN)
                for (b = 0; b < nb; b++)
                    _nmod_vec_scalar_addmul_nmod(out + off[b], in + off[b],
                            off[b + 1] - off[b], op->eig[b], op->mod);
            break;

        case _NMOD_MAT_KRYLOV_DIAGONAL:
            _nmod_vec_mul(out, op->eig, in, nb, op->mod);
            break;

        default:
            flint_throw(FLINT_ERROR, "_nmod_mat_krylov_op_apply_vec: dense operator\n");
    }
}
