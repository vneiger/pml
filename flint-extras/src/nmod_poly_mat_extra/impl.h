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

#include <flint/nmod_types.h>
#include <flint/nmod_poly.h>  /* for geometric_progression_t */
#include <flint/thread_pool.h>
#include <flint/ulong_extras.h>  /* for n_sqrt */

#include "nmod_mat_poly.h"     /* for nmod_mat_poly_t */

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


/* ------------------------------------------------------------------------ */
/* a few functions to help tests/profiles be more complete                  */
/* ------------------------------------------------------------------------ */

/* Same as nmod_poly_mat_set_trunc_from_mat_poly, with explicit control over
 * the block kernel and the loop schedule of the underlying transposition.
 *
 * `kern` is one of NMOD_MAT_POLY_CONV_{SCALAR,VEC4,VEC8} (see
 * nmod_mat_poly_extra/impl.h); any other value selects the default, narrowed
 * if its block does not fit inside the problem.
 *
 * `dmaj` is 1 for the entry-major schedule (the stores into the output
 * polynomials are long sequential streams) and 0 for the coefficient-major
 * one (the loads from the input matrices are long sequential streams); any
 * other value (e.g. -1) selects the default.  Note that the two schedules
 * carry the opposite names in the other direction: `dmaj` always means
 * "sweep the destination rows sequentially". */
void _nmod_poly_mat_set_trunc_from_mat_poly(nmod_poly_mat_t pmat,
                                            const nmod_mat_poly_t matp,
                                            slong order,
                                            int kern,
                                            int dmaj);

/* Same as _nmod_poly_mat_mul_geometric_precomp, with the soft bound (in
 * bytes) on the memory used for the constant matrices given explicitly;
 * 0 selects the default. Below the bound the product is performed in one
 * go, above it the rows of A and then the columns of B are processed by
 * groups. Exposed so that the tests can exercise the grouping at sizes
 * that do not need hundreds of megabytes. */
void _nmod_poly_mat_mul_geometric_precomp_bounded(nmod_poly_mat_t res,
                                          const nmod_poly_mat_t pmat1, slong len1,
                                          const nmod_poly_mat_t pmat2, slong len2,
                                          nmod_geometric_progression_t G,
                                          ulong membytes);

/* Same for the middle product _nmod_poly_mat_mulmid_geometric_precomp. */
void _nmod_poly_mat_mulmid_geometric_precomp_bounded(nmod_poly_mat_t res,
                                             const nmod_poly_mat_t pmat1, slong len1,
                                             const nmod_poly_mat_t pmat2, slong len2,
                                             slong nlo, slong nhi,
                                             nmod_geometric_progression_t G,
                                             ulong membytes);

/* ------------------------------------------------------------------------ */
/* shared machinery of the products at a geometric progression              */
/* ------------------------------------------------------------------------ */

/*
    What follows is shared by nmod_poly_mat_mul_geometric and
    nmod_poly_mat_mulmid_geometric : the layout of the constant matrices, the
    retained buffer that holds them, the splitting of the operands into groups
    of rows and columns that fit a memory bound, and the dispatch of a phase
    over the thread pool.
*/

/* Soft bound, in bytes, on the memory used for the constant matrices.
   This is not a constant but a floor: the bound used is the larger of it
   and twice the size of the operands and the result. */
#define NMOD_POLY_MAT_GEOMETRIC_MEM_FLOOR (UWORD(1) << 28)

/* a rounded up to the next multiple of the power of two b (n_round_up
   itself lives in fft_small.h, which these files must not depend on) */
static inline ulong _round_up2(ulong a, ulong b)
{
    return (a + b - 1) & ~(b - 1);
}

/*
    Retained scratch space for the constant matrices, thread local and only
    ever grown; flint_cleanup() releases it. Defined in mul_geometric.c.
*/
ulong * _nmod_poly_mat_geometric_fit_buffer(ulong nbytes);

/*
    The evaluations of one region (the block of B, of A or of C currently
    held) at one point form one matrix, laid out contiguously and row major
    so that an nmod_mat can borrow it:

        region + t * mstride,   mstride = _mats_stride(r, c),

    where r and c are the dimensions of the region and the row stride
    inside a matrix is c. The distance between consecutive points is a
    whole number of cache lines, on a 64-byte aligned base, which keeps
    each group of 8 consecutive entries of one point inside a single
    line: the scatter below writes the region one such group at a time,
    and two threads never write the same line. It is an odd number of
    them, see _odd_lines.
*/
typedef struct
{
    ulong * data;
    slong r;
    slong c;
    slong mstride;
    slong len;
}
_mats_struct;

static inline ulong * _mats_at(const _mats_struct * M, slong t)
{
    return M->data + t * M->mstride;
}

/*
    words rounded up to an odd number of cache lines (of 8 words). The
    distance between the matrices of consecutive points, and between the
    three regions, is made such a number: with a power of two of entries
    in a matrix the natural distances are large powers of two, and every
    matrix then lands on the same few cache sets -- measured at 10-20% of
    the whole product for 128 x 128 matrices, in the pointwise stage and
    in the strided scatter. One line per point is the most this costs.
*/
static inline ulong _odd_lines(ulong words)
{
    words = _round_up2(words, 8);
    return ((words / 8) % 2 == 0) ? words + 8 : words;
}

/* distance in words between the matrices of consecutive points */
static inline slong _mats_stride(slong r, slong c)
{
    return (slong) _odd_lines((ulong) r * c);
}

/* words occupied by a region of r x c matrices at len points */
static inline ulong _mats_words(slong r, slong c, slong len)
{
    return _odd_lines((ulong) len * (ulong) _mats_stride(r, c));
}

/* the len values of entry e, scattered over the matrices */
static inline void _mats_scatter(const _mats_struct * M, slong e, nn_srcptr v)
{
    const slong len = M->len, ms = M->mstride;
    ulong * dst = M->data + e;
    slong t;

    for (t = 0; t < len; t++)
        dst[t * ms] = v[t];
}

/* the same for an entry which is the zero polynomial */
static inline void _mats_scatter_zero(const _mats_struct * M, slong e)
{
    const slong len = M->len, ms = M->mstride;
    ulong * dst = M->data + e;
    slong t;

    for (t = 0; t < len; t++)
        dst[t * ms] = 0;
}

/* the first zl values of entry e, gathered from the matrices */
static inline void _mats_gather(nn_ptr v, const _mats_struct * M, slong e, slong zl)
{
    const slong ms = M->mstride;
    const ulong * src = M->data + e;
    slong t;

    for (t = 0; t < zl; t++)
        v[t] = src[t * ms];
}

/*
    The groups: NRG rows of A and of C, NCG columns of B and of C are held
    at a time, so that k*NCG + NRG*k + NRG*NCG entries fit in `avail`
    (counted in entries, i.e. in sets of values at all the points). The
    two are balanced -- they are the outer dimensions of the matrices
    handed to nmod_mat_mul, which changes algorithm on the smallest of the
    three -- rather than the rows of A alone being cut.
*/
static inline void _geometric_groups(slong m, slong k, slong n, ulong avail,
                                     slong * NRG, slong * NCG)
{
    *NRG = m;
    *NCG = n;
    if ((ulong) k * n + (ulong) m * k + (ulong) m * n > avail)
    {
        /* NRG = NCG = g uses 2*k*g + g^2 entries */
        slong g = (slong) (n_sqrt((ulong) k * k + avail) - (ulong) k);
        g = FLINT_MAX(g, 1);
        *NRG = FLINT_MIN(m, g);
        /* the largest NCG for that NRG, and then the largest NRG for
           that NCG: this gives back to one side what the other could
           not use because m or n was the binding cap */
        *NCG = (avail > (ulong) *NRG * k)
               ? (slong) ((avail - (ulong) *NRG * k) / ((ulong) k + *NRG)) : 1;
        *NCG = FLINT_MAX(*NCG, 1);
        *NCG = FLINT_MIN(*NCG, n);
        *NRG = (avail > (ulong) k * *NCG)
               ? (slong) ((avail - (ulong) k * *NCG) / ((ulong) k + *NCG)) : 1;
        *NRG = FLINT_MAX(*NRG, 1);
        *NRG = FLINT_MIN(*NRG, m);
    }
}

/*
    A phase is a loop over independent tasks, split evenly over the worker
    descriptors W[0], ..., W[nthreads-1] and run over the thread pool, the
    calling thread taking W[0]. Every worker descriptor type begins with
    this range, which is how the split is written into it.
*/
typedef struct
{
    slong start;
    slong stop;
}
_geometric_range_struct;

static inline void _geometric_run(void (* func)(void *),
                                  void * W, size_t wsize, slong nthreads,
                                  thread_pool_handle * handles, slong ntasks)
{
    char * base = (char *) W;
    slong i;

    if (ntasks < 1)
        return;

    nthreads = FLINT_MIN(nthreads, ntasks);

    for (i = 0; i < nthreads; i++)
    {
        _geometric_range_struct * R = (_geometric_range_struct *) (base + i * wsize);
        R->start = (i + 0) * ntasks / nthreads;
        R->stop  = (i + 1) * ntasks / nthreads;
    }

    for (i = nthreads - 1; i > 0; i--)
        thread_pool_wake(global_thread_pool, handles[i - 1], 0, func, base + i * wsize);
    func(base);
    for (i = nthreads - 1; i > 0; i--)
        thread_pool_wait(global_thread_pool, handles[i - 1]);
}

#endif  /* NMOD_POLY_MAT_EXTRA_IMPL_H */
