/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#ifndef NMOD_MAT_POLY_EXTRA_IMPL_H
#define NMOD_MAT_POLY_EXTRA_IMPL_H

#include <flint/nmod_types.h>

#include "nmod_mat_poly.h"

/* Block kernels for the nmod_poly_mat to nmod_mat_poly conversion: values for
 * the `kern` parameter of _nmod_mat_poly_set_trunc_from_poly_mat. A kernel
 * that the build does not provide is silently replaced by the widest
 * available one, so passing any of these is always safe. */
#define NMOD_MAT_POLY_CONV_SCALAR 0  /* 8x8 blocks, no vector instructions */
#define NMOD_MAT_POLY_CONV_VEC4   1  /* 4x4 blocks, FLINT machine vectors (AVX2, NEON) */
#define NMOD_MAT_POLY_CONV_VEC8   2  /* 8x8 blocks, AVX-512 */

/* Same as nmod_mat_poly_set_trunc_from_poly_mat, with explicit control over
 * the block kernel, the loop schedule and the prefetching
 *
 * `kern` is one of NMOD_MAT_POLY_CONV_{SCALAR,VEC4,VEC8};
 * any other value selects the widest kernel the build provides, narrowed if
 * its block does not fit inside the problem.
 *
 * `cmaj` is 1 for the coefficient-major schedule (the blocks are visited so
 * that the stores into the output matrices are long sequential streams), 0
 * for the entry-major schedule (the loads from the input polynomials are long
 * sequential streams); any other value (e.g. -1) selects the default.
 *
 * `pf` is 1 to software-prefetch the side that is accessed one cache line per
 * row (the input polynomials under the coefficient-major schedule, the output
 * matrices under the entry-major one), 0 to disable it; any other value
 * (e.g. -1) enables it only when the data is large enough, see
 * PML_CONV_PREFETCH_BYTES. */
void _nmod_mat_poly_set_trunc_from_poly_mat(nmod_mat_poly_t matp,
                                            const nmod_poly_mat_t pmat,
                                            slong order,
                                            int kern,
                                            int cmaj,
                                            int pf);

#endif /* ifndef NMOD_MAT_POLY_EXTRA_IMPL_H */
