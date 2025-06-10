/*
    Copyright (C) 2025 Vincent Neiger, Éric Schost

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#ifndef __FMPZ_POLY_MAT_EXTRA_H
#define __FMPZ_POLY_MAT_EXTRA_H

#include <flint/fmpz.h>
#include <flint/fmpz_poly_mat.h>

#include "pml.h"

/** ------------------------------------------------------------ */
/** Waksman's algorithm for matrix multiplication                */
/** does n^3/2+O(n^2) products, but many additions               */
/** good for small matrices with large entries                   */
/** ------------------------------------------------------------ */
void fmpz_poly_mat_mul_waksman(fmpz_poly_mat_t C, const fmpz_poly_mat_t A, const fmpz_poly_mat_t B);


/** ------------------------------------------------------------ */
/** a multimodular algorithm for matrix multiplication           */
/** ------------------------------------------------------------ */
void fmpz_poly_mat_mul_multimod(fmpz_poly_mat_t C, const fmpz_poly_mat_t A, const fmpz_poly_mat_t B);

#endif
