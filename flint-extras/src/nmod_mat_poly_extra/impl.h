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
#include <flint/nmod.h>
#include <flint/nmod_mat.h>

#include "nmod_mat_poly.h"

/* ------------------------------------------------------------------------ */
/* mbasis and mintbasis (nmod_mat_poly_mbasis.c, nmod_mat_poly_mintbasis.c) */
/* ------------------------------------------------------------------------ */

/* Sorts perm (a permutation of 0..m-1) so that the pairs (shift[perm[i]],
 * perm[i]) increase: the rows sorted by shift, with priority to the smaller
 * index for equal shifts. Insertion sort: from one iteration of mbasis to
 * the next, perm is nearly sorted. */
void _nmod_mat_poly_shift_perm(slong * perm, const slong * shift, slong m);

/* Left nullspace of the m x n matrix A, in the compact form computed by
 * FLINT's nmod_mat_ple_left_nullspace_compact (directly from the first
 * nonzero entry for n = 1); used at each step of mbasis and mintbasis. A is
 * overwritten (except for n = 1). Returns the rank r of A. On return, X is a window (cleared by
 * the caller with nmod_mat_window_clear) of the m x m matrix K, of size
 * (m - r) x r: the basis without its identity part. If 0 < r < m, `rows`
 * (length m) holds the row rank profile of A followed by the other rows, both
 * in increasing order: row k of the basis has a 1 in position rows[r+k],
 * zeros at the other rows[r+j], and X[k][j] in position rows[j]; if r = m,
 * rows is the identity permutation (all rows are pivot rows); if r = 0, rows
 * is not set. P and Q are scratch arrays of length m and n. */
slong _nmod_mat_poly_left_nullspace_compact(nmod_mat_t X, slong * rows, nmod_mat_t A,
                                            const nmod_mat_t K, slong * P, slong * Q);

/* W is m x (cap m), its row i holding the coefficients 0..len-1 of the row i
 * of an m x m polynomial matrix (m entries each), followed by zeros. If need >
 * cap, reallocates W with room for max(need, min(2*cap, maxcap)) coefficients,
 * keeping the coefficients, and updates cap. Requires need <= maxcap. */
void _nmod_mat_poly_basis_reserve(nmod_mat_t W, slong * cap, slong len,
                                  slong need, slong maxcap);

/* Makes sure that A has at least r rows and c columns, reallocating it (its
 * entries are not kept) if needed: for scratch matrices whose size grows. */
static inline void
_nmod_mat_poly_scratch_fit(nmod_mat_t A, slong r, slong c)
{
    if (A->r < r || A->c < c)
    {
        const ulong p = A->mod.n;
        r = FLINT_MAX(r, A->r);
        c = FLINT_MAX(c, A->c);
        nmod_mat_clear(A);
        nmod_mat_init(A, r, c, p);
    }
}

/* In mbasis and mintbasis, the other rows receive a combination of the rank
 * pivot rows of the residual: up to the rank returned here, by scalar times
 * row operations, and otherwise by one matrix product, whose inner dimension
 * is this rank. First crossovers from measurements in 2026-11, balanced
 * uniform random input, m = 2 to 64, orders 8, 32, 128, Zen 4, Comet Lake,
 * Apple M4, Emerald Rapids VM: up to 2^52, the row operations win up to rank 3
 * (Zen 4, M4) to 6 (Comet Lake, VM); above 2^52, where nmod_mat_mul has no
 * SIMD kernel, up to rank 4 (M4) to 8 (Zen 4, VM). Hence 3 and 4 on aarch64
 * (measured on the M4 only), 4 and 8 elsewhere.
 * TODO rough values, to be refined in a tuning stage (per machine?) */
#ifndef NMOD_MAT_POLY_BASIS_LOWRANK
# if defined(__aarch64__)
#  define NMOD_MAT_POLY_BASIS_LOWRANK 3
# else
#  define NMOD_MAT_POLY_BASIS_LOWRANK 4
# endif
#endif
#ifndef NMOD_MAT_POLY_BASIS_LOWRANK_LARGE
# if defined(__aarch64__)
#  define NMOD_MAT_POLY_BASIS_LOWRANK_LARGE 4
# else
#  define NMOD_MAT_POLY_BASIS_LOWRANK_LARGE 8
# endif
#endif

static inline slong
_nmod_mat_poly_basis_lowrank(nmod_t mod)
{
    return (NMOD_BITS(mod) > 52) ? NMOD_MAT_POLY_BASIS_LOWRANK_LARGE
                                 : NMOD_MAT_POLY_BASIS_LOWRANK;
}

/* The choice of the variant of mbasis, for an m x n input at order `order`
 * modulo a prime of `bits` bits. The update of the residual series
 * (resupdate) costs about (m-n) n^2 order^2 / 2 operations in total, its
 * recomputation (rescomp) about m n^2 order^2 / 2; but the products of
 * rescomp have a larger inner dimension (m len instead of the rank,
 * generically n), hence fewer reductions per operation. Rescomp from order 64
 * on for n <= m/2 and m >= 3, resupdate otherwise.
 *
 * FIXME TEMPORARY NOTES: 2026-11 (grid, m = 2 to 64, orders 8, 32, 128, 31 and 60
 * bits; Zen 4, Comet Lake, Apple M4), 0.9% to 2.6% slower than the better
 * variant on average on each machine (at most 1.42x), against 1.7% to 6.5% for
 * the previous rule by modulus; for m = 2 the machines disagree, and resupdate
 * is up to 1.15x faster at order 64 to 128 on Zen 4 and the Emerald Rapids VM.
 * */
static inline int
_nmod_mat_poly_mbasis_resupdate_wanted(slong m, slong n, slong order, flint_bitcnt_t bits)
{
    (void) bits;
    return order < 64 || 2 * n > m || m <= 2;
}

/* mbasis, as nmod_mat_poly_mbasis, with the variant (resupdate: 1, rescomp:
 * 0) and the rank up to which the update is done by scalar times row
 * operations (see _nmod_mat_poly_basis_lowrank) given */
void _nmod_mat_poly_mbasis(nmod_mat_poly_t appbas, slong * shift,
                           const nmod_mat_poly_t matp, slong order,
                           int resupd, slong lowrank);

/* The choice of the variant of mintbasis, for d points and m x n matrices
 * E_k, modulo a prime of `bits` bits: resupdate when n d <= 6 m, or for
 * n = 1 and m <= 4, rescomp otherwise.
 *
 * FIXME TEMPORARY NOTES 2026-11 (grid, m = 2 to 64,
 * d = 8, 32, 128, 31 and 60 bits; Zen 4, Comet Lake, Apple M4): 0.9% to 2.3%
 * slower than the better variant on average on each machine (at most 1.37x),
 * against 6.8% to 17% for the previous rule d (m-n+1) <= m, which nearly
 * always chose rescomp. */
static inline int
_nmod_mat_poly_mintbasis_resupdate_wanted(slong m, slong n, slong d, flint_bitcnt_t bits)
{
    (void) bits;
    return n * d <= 6 * m || (n == 1 && m <= 4);
}

/* mintbasis, as nmod_mat_poly_mintbasis, with the variant (resupdate: 1,
 * rescomp: 0) and the rank up to which the update is done by scalar times
 * row operations (see _nmod_mat_poly_basis_lowrank) given */
void _nmod_mat_poly_mintbasis(nmod_mat_poly_t intbas, slong * shift,
                              const ulong * pts, const nmod_mat_struct * E, slong d,
                              int resupd, slong lowrank);

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

/* Below this many words, the pointer tables of _pml_transpose and their
 * allocation dominate: the conversions use a plain loop instead. */
#define _PML_CONV_TINY 512

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
