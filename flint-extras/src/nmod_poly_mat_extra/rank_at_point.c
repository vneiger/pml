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
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>

#include "nmod_poly_mat_utils.h"

/* rank of the constant matrix M, which is overwritten (LU in place, without
 * the copy made by nmod_mat_rank) */
static slong _nmod_mat_rank_inplace(nmod_mat_t M)
{
    if (M->r == 0 || M->c == 0)
        return 0;
    slong * perm = FLINT_ARRAY_ALLOC(M->r, slong);
    const slong rk = nmod_mat_lu(perm, M, 0);
    flint_free(perm);
    return rk;
}

slong nmod_poly_mat_rank_at_point(const nmod_poly_mat_t pmat, ulong pt)
{
    nmod_mat_t M;
    nmod_mat_init(M, pmat->r, pmat->c, pmat->modulus);
    nmod_poly_mat_evaluate_nmod(M, pmat, pt);
    const slong rk = _nmod_mat_rank_inplace(M);
    nmod_mat_clear(M);
    return rk;
}

slong nmod_poly_mat_rank_at_zero(const nmod_poly_mat_t pmat)
{
    nmod_mat_t M;
    nmod_mat_init(M, pmat->r, pmat->c, pmat->modulus);
    nmod_poly_mat_get_coeff_mat(M, pmat, 0);
    const slong rk = _nmod_mat_rank_inplace(M);
    nmod_mat_clear(M);
    return rk;
}
