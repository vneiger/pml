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

/* Block kernels for the transposition that underlies both conversions
 * between nmod_poly_mat and nmod_mat_poly: values for the `kern` parameter
 * below.  A kernel that the build does not provide is silently replaced by
 * the widest available one, so passing any of these is always safe. */
#define NMOD_MAT_POLY_CONV_SCALAR 0  /* 8x8 blocks, no vector instructions */
#define NMOD_MAT_POLY_CONV_VEC4   1  /* 4x4 blocks, FLINT machine vectors (AVX2, NEON) */
#define NMOD_MAT_POLY_CONV_VEC8   2  /* 8x8 blocks, AVX-512 */

/* Default loop schedule: 1 to sweep the destination rows sequentially and
 * scatter the loads, 0 to sweep the source rows sequentially and scatter the
 * stores.  (For the nmod_poly_mat -> nmod_mat_poly direction the first is the
 * coefficient-major schedule and the second the entry-major one; for the
 * other direction the names swap, which is why they are not used here.)
 *
 * Which one wins depends on how much of a cache line one scattered store
 * covers.  On x86, where a line is 64 bytes, a block is a whole line (8-wide
 * kernel) or half of one (4-wide), the store needs no line to be kept around,
 * and scattering the stores is nearly free: on Zen 4 that schedule is 1.6x to
 * 2.4x faster from 16 MB of data upwards, and within noise below.
 *
 * On Apple silicon a line is 128 bytes and the widest kernel is the 4-wide
 * NEON one, so a scattered store covers a quarter of a line: the line has to
 * be fetched for ownership and then kept until the remaining quarters are
 * written, which only happens after a full sweep of the inner loop.  That
 * schedule then loses badly -- on an M4 it is 3x to 6x slower as soon as the
 * matrix reaches 16x16, and slower than the naive loop itself from 32x32 on.
 * With the stores sequential instead, consecutive blocks fill each line back
 * to back and the line size stops mattering.
 *
 * Overridable at build time. */
#ifndef PML_CONV_DST_MAJOR
# if defined(__APPLE__) && defined(__aarch64__)
#  define PML_CONV_DST_MAJOR 1
# else
#  define PML_CONV_DST_MAJOR 0
# endif
#endif

/* Transposition between two collections of separately allocated rows:
 *
 *     dst[i][j] = (i < slen[j]) ? src[j][i] : 0,   i < ndst,  j < nsrc
 *
 * so `dst` has `ndst` rows of at least `nsrc` words each, and `src` has
 * `nsrc` rows, row `j` holding `slen[j] <= ndst` words (the rest reads as
 * zero and is never touched).  Pass `slen[j] == ndst` for every `j` when the
 * source rows are all full.
 *
 * `kern` is one of NMOD_MAT_POLY_CONV_{SCALAR,VEC4,VEC8}, and must be a
 * kernel whose block fits: see _pml_transpose_narrow_kernel.  `dmaj` is 1 for
 * the destination-major schedule, 0 for the source-major one. */
void _pml_transpose(nn_ptr * dst, slong ndst,
                    nn_srcptr * src, const slong * slen, slong nsrc,
                    int kern, int dmaj);

/* Steps `kern` down to a kernel whose block fits inside a `ndst x nsrc`
 * transposition, and returns it. */
int _pml_transpose_narrow_kernel(int kern, slong ndst, slong nsrc);

/* The widest kernel this build provides. */
int _pml_transpose_widest_kernel(void);

/* Same as nmod_mat_poly_set_trunc_from_poly_mat, with explicit control over
 * the block kernel and the loop schedule.
 *
 * `kern` is one of NMOD_MAT_POLY_CONV_{SCALAR,VEC4,VEC8};
 * any other value selects the widest kernel the build provides, narrowed if
 * its block does not fit inside the problem.
 *
 * `dmaj` is 1 for the coefficient-major schedule (the stores into the output
 * matrices are long sequential streams), 0 for the entry-major schedule (the
 * loads from the input polynomials are long sequential streams); any other
 * value (e.g. -1) selects the default, PML_CONV_DST_MAJOR. */
void _nmod_mat_poly_set_trunc_from_poly_mat(nmod_mat_poly_t matp,
                                            const nmod_poly_mat_t pmat,
                                            slong order,
                                            int kern,
                                            int dmaj);

#endif /* ifndef NMOD_MAT_POLY_EXTRA_IMPL_H */
