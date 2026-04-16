/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#ifndef NMOD_POLY_MAT_EXTRA_IMPL_H
#define NMOD_POLY_MAT_EXTRA_IMPL_H

#include <flint/nmod_poly_mat.h>

/* Hermite form helpers */
void _atomic_solve_pivot_collision_uechelon_rowwise(nmod_poly_mat_t mat, nmod_poly_mat_t other,
                                                    slong pi1, slong pi2, slong j);
void _pivot_collision_xgcd_uref(nmod_poly_mat_t mat, nmod_poly_mat_t other,
                                                      slong pi, slong ii, slong j,
                                                      nmod_poly_t g, nmod_poly_t u, nmod_poly_t v,
                                                      nmod_poly_t pivg, nmod_poly_t nonzg);
void _reduce_against_pivot_uref(nmod_poly_mat_t mat, nmod_poly_mat_t other,
                                            slong i, slong j, slong ii,
                                            nmod_poly_t u, nmod_poly_t v);
ulong _normalize_pivot_uref(nmod_poly_mat_t mat, nmod_poly_mat_t other, slong i, slong j);
void _normalize_uref(nmod_poly_mat_t mat, nmod_poly_mat_t other, slong * pivind, slong rk);

/* weak Popov form helpers */
void _atomic_solve_pivot_collision_rowwise(nmod_poly_mat_t mat, nmod_poly_mat_t other,
                                           slong pi1, slong pi2, slong j);
void _reduce_against_pivot_general_rowwise(nmod_poly_mat_t mat, nmod_poly_mat_t other,
                                   slong i, slong j, slong ii,
                                   nmod_poly_t u, nmod_poly_t v);
ulong _normalize_pivot_general_rowwise(nmod_poly_mat_t mat, nmod_poly_mat_t other, slong i, slong j);

#endif  /* NMOD_POLY_MAT_EXTRA_IMPL_H */

