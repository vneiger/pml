/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <string.h>  /* memcpy */

#include <flint/nmod.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>
#include <flint/ulong_extras.h>

#include "nmod_poly_mat_multiply.h"

#if PML_HAVE_MACHINE_VECTORS

#include <flint/fft_small.h>
#include <flint/machine_vectors.h>
#include <flint/thread_pool.h>
#include <flint/thread_support.h>

/* timing breakdown of the phases, printed on stderr (development aid) */
#ifdef PML_MUL_SD_FFT_DIRECT_TIMING
# include <stdio.h>
# include <stdlib.h>
# include <time.h>
static double _now(void) { struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return ts.tv_sec + 1e-9 * ts.tv_nsec; }
# define TIMING_DECL double _t0 = _now(), _t1; double _tB = 0, _tA = 0, _tP = 0, _tC = 0;
# define TIMING_MARK(acc) do { _t1 = _now(); acc += _t1 - _t0; _t0 = _t1; } while (0)
# define TIMING_PRINT flint_fprintf(stderr, "  [sd_fft_direct] np=%wu depth=%wu ztrunc=%wu T=%wu NRG=%wd NCG=%wd | B fft %.3f, A fft %.3f, pointwise %.3f, ifft+crt %.3f\n", np, P->depth, ztrunc, T, NRG, NCG, _tB, _tA, _tP, _tC)
#else
# define TIMING_DECL
# define TIMING_MARK(acc)
# define TIMING_PRINT
#endif

/*
    Polynomial matrix multiplication C = A * B by evaluation-interpolation
    with FLINT's fft_small transforms, using the naive (cubic) matrix
    multiplication algorithm directly on the transforms:

        1. plan: fix the FFT depth from lenA + lenB - 1 and the set of
           primes from the bound on the output coefficients (before
           reduction modulo p they are sums of at most
           k * min(lenA, lenB) products of two residues, k the inner
           dimension), so that np primes are used with
           prod(primes) >= k * min(lenA, lenB) * p^2. For a 64-bit p and
           the 50-bit FFT primes this means 3 primes unless
           k * min(lenA, lenB) exceeds 2^22, in which case 4; smaller p
           may need only 1 or 2. When p is itself an FFT prime of at most
           50 bits with enough 2-adicity, a single transform modulo p is
           used instead (no CRT at all).
        2. transform every entry of B, and every entry of A (all at once
           when memory allows, otherwise by groups of rows of A and, if
           needed, of columns of B);
        3. for each prime and each tile of the transform, compute
           sum_l A[i][l] * B[l][j] pointwise for every (i, j), with the
           normalization factor folded once into the transforms of B,
           and with the reductions of the sums batched (see _dot_tile);
        4. inverse transform and chinese remainder each entry of C.

    Step 3 is organized tile by tile -- a tile is a range of at most
    BLK_SZ points of the transforms -- with the transforms stored tile
    major (see _tiles_struct), so that one pass over the tiles of one
    prime reads each transform of A and of B exactly once, in order,
    while the innermost work reuses the tiles of a block of 8 columns of
    B out of L1 across all the rows of A. The pointwise stage is then
    close to the arithmetic peak rather than bandwidth bound.

    Each of the four steps is a loop over independent tasks and is run
    over the FLINT thread pool (see _sd_fft_direct_run).

    Memory: (m*k + k*n + m*n) transforms of np * ztrunc doubles when it
    fits the budget (see NMOD_POLY_MAT_MUL_SD_FFT_DIRECT_MEM_FLOOR), otherwise
    the rows of A -- and, if that is still not enough, the columns of B
    -- are processed by groups until it does (see the choice of NRG and
    NCG below).
*/

/*
   Soft bound, in bytes, on the memory used for the transforms.
   Grouping the rows of A within it is free in transform count -- each
   entry of A and of B is still transformed exactly once -- and only
   trades memory for bandwidth, since every group of rows streams the
   transforms of B once. Grouping the columns of B, which only happens
   when the k*n transforms of B alone exceed the budget, does cost
   transforms: A is then transformed once per group of columns.

   The bound is not a constant but scales with the problem: the
   transforms of a product whose operands are themselves a gigabyte have
   no reason to be capped at a few hundred megabytes, and the bandwidth
   the grouping trades away is not free. Measured at dimension 256,
   length 2*256, modulo a generic 60-bit prime (three primes): 26.2 s
   within the constant below, against 20.9 s undivided. The constant is
   therefore only a floor, for the small products where a fixed working
   set is what one wants to bound, and the factor is about what an
   undivided square product needs relative to its operands when a single
   prime suffices.

   The transforms are taken from the retained fft_small scratch buffer
   for any request within that bound (see _sd_fft_direct_alloc), so a
   thread that has performed one large product keeps a working set of
   the order of a small multiple of the operands it was given.
*/
#define NMOD_POLY_MAT_MUL_SD_FFT_DIRECT_MEM_FLOOR (UWORD(1) << 28)

/* Soft bound, in bytes, on the tiles of B read by one call of the
   register-blocked kernel (8 columns of B, i.e. 8 * k tiles): kept within the
   L1 cache, these tiles are reused from L1 for all the rows of A. */
#define NMOD_POLY_MAT_MUL_SD_FFT_DIRECT_TILE_BUDGET (UWORD(1) << 15)

/*
    Allocation of the transform storage, up to the budget above. Freshly
    mapped memory of that size page-faults on first touch, which costs a
    sizeable fraction of a single product, so the storage is taken from
    the scratch buffer of the fft_small context, as FLINT's fused
    multiplication drivers do. That buffer is thread local and retained
    across calls, so a sequence of products of similar sizes -- what the
    polynomial matrix algorithms of PML do -- pays the page faults once
    rather than once per product. It is exclusive to us here: none of the
    fft_small entry points called below uses it.

    The grouping keeps the request within the budget except when even one
    row of A against one column of B exceeds it; that degenerate case is
    served by an ordinary allocation rather than by growing the retained
    buffer beyond the budget for the rest of the thread's life.

    (Asking for transparent huge pages on the buffer was measured slower:
    the madvise call, once per product, costs more than the page walks it
    saves.)
*/
static double * _sd_fft_direct_alloc(mpn_ctx_struct * R, ulong nbytes, ulong budget)
{
    if (nbytes <= budget)
        return (double *) mpn_ctx_fit_buffer(R, nbytes);
    return flint_aligned_alloc(FLINT_FFT_SMALL_ALIGNMENT,
                               n_round_up(nbytes, FLINT_FFT_SMALL_ALIGNMENT));
}

static void _sd_fft_direct_free(double * buf, ulong nbytes, ulong budget)
{
    if (nbytes > budget)
        flint_aligned_free(buf);
}

/* transform-space truncation of an operand of length an (see
   _len_trunc/_op_trunc in FLINT's fft_small/nmod_poly_mul.c): the full
   transform below one block, whole blocks otherwise */
static inline ulong _op_trunc(ulong an, const fft_small_plan_t P)
{
    return (P->depth < LG_BLK_SZ) ? n_pow2(P->depth)
                                  : n_round_up(an, BLK_SZ);
}

/*
    Multiply each of the np transforms of X by the plan's normalization
    factor, in place, restricted to the ztrunc points read by the inverse
    transform. Forward transforms have entries in (-4n, 4n) and m is
    reduced to (-n/2, n/2], so the products are in (-2n^2, 2n^2) and the
    results in (-n, n) (see the range documentation of mulmod in
    machine_vectors.h).
*/
static void _op_scale_by_m(fft_small_op_t X, const fft_small_plan_t P)
{
    ulong pi, t;
    for (pi = 0; pi < P->np; pi++)
    {
        const sd_fft_ctx_struct * Q = P->ffts + P->offset + pi;
        vec8d n = vec8d_set_d(Q->p);
        vec8d ninv = vec8d_set_d(Q->pinv);
        vec8d m = vec8d_set_d(vec1d_reduce_0n_to_pmhn((slong) P->m[pi], Q->p));
        double * x = X->data + pi * P->stride;
        for (t = 0; t < P->ztrunc; t += 8)
            vec8d_store(x + t, vec8d_mulmod(vec8d_load(x + t), m, n, ninv));
    }
}

/*
    Pointwise dot products on one tile of one prime: for 0 <= c < NC,

        z[c][t] = sum_{l < nj} a[l][t] * b[c][l][t],     0 <= t < T,

    reduced into [-n, n]. Requires entries of a in (-4n, 4n) (raw forward
    transforms) and entries of b in (-n, n) (forward transforms scaled by
    the normalization factor), T a multiple of 8, n < 2^50.

    Each product a*b is in (-4n^2, 4n^2) hence mulmod returns it in
    (-3n/2, 3n/2). Sums are reduced only every four products: the running
    total is then bounded by n + 4 * 3n/2 = 7n < 2^53, so it is exact in
    double precision, and reduce_to_pm1n brings it back to [-n, n].

    The NC outputs share the loads of a (register blocking); NC = 8 was
    measured fastest, about twice as fast as NC = 1 with L2-resident
    tiles, on AVX2 as well as AVX-512 (there the NC = 8 accumulators do
    not fit the registers, and the compiler spills, without much harm).
*/

#define DOT_TILE_GROUP(NC, Q) \
        for (t = 0; t < T; t += 8) \
        { \
            vec8d acc[NC], x; \
            int c, u; \
            for (c = 0; c < NC; c++) \
                acc[c] = (l == 0) ? vec8d_set_d(0.0) : vec8d_load(z[c] + t); \
            for (u = 0; u < (Q); u++) \
            { \
                x = vec8d_load(a[l + u] + t); \
                for (c = 0; c < NC; c++) \
                    acc[c] = vec8d_add(acc[c], \
                        vec8d_mulmod(x, vec8d_load(b[c][l + u] + t), n, ninv)); \
            } \
            for (c = 0; c < NC; c++) \
                vec8d_store(z[c] + t, vec8d_reduce_to_pm1n(acc[c], n, ninv)); \
        }

#define DEFINE_DOT_TILE(NC) \
static void _dot_tile_##NC(double ** z, const double ** a, \
                           const double * const * const * b, slong nj, \
                           ulong T, vec8d n, vec8d ninv) \
{ \
    slong l; \
    ulong t; \
    for (l = 0; l + 4 <= nj; l += 4) \
    { \
        DOT_TILE_GROUP(NC, 4) \
    } \
    if (l < nj) \
    { \
        const int q = (int) (nj - l); \
        DOT_TILE_GROUP(NC, q) \
    } \
}

DEFINE_DOT_TILE(1)
DEFINE_DOT_TILE(2)
DEFINE_DOT_TILE(4)
DEFINE_DOT_TILE(8)

#undef DEFINE_DOT_TILE
#undef DOT_TILE_GROUP

/* z[c] = sum_l a[l] * b[c][l] for 0 <= c < nc, any nc > 0, on one tile */
static void _dot_tile(double ** z, const double ** a,
                      const double * const * const * b, slong nc, slong nj,
                      ulong T, vec8d n, vec8d ninv)
{
    while (nc >= 8)
    {
        _dot_tile_8(z, a, b, nj, T, n, ninv);
        z += 8; b += 8; nc -= 8;
    }
    if (nc >= 4)
    {
        _dot_tile_4(z, a, b, nj, T, n, ninv);
        z += 4; b += 4; nc -= 4;
    }
    if (nc >= 2)
    {
        _dot_tile_2(z, a, b, nj, T, n, ninv);
        z += 2; b += 2; nc -= 2;
    }
    if (nc >= 1)
        _dot_tile_1(z, a, b, nj, T, n, ninv);
}

/*
    Tile-major storage of the transforms. The pointwise stage works on
    one tile (T consecutive points) of one prime of many transforms at a
    time; storing the transforms as fft_small does (each transform
    contiguous, np * 2^depth doubles) would make it read k*n + ... short
    chunks scattered over as many pages, which TLB misses and poor
    prefetching make about twice as slow. Instead, for each prime pi and
    tile ti, the tiles of all the transforms of a region (B, or the current
    rows of A, or of C) are stored consecutively:

        region + ((pi * ntiles + ti) * nops + o) * T,   0 <= o < nops.

    Transforms are computed in a contiguous scratch op and scattered into
    this layout (and gathered back for the inverse transforms); the copies
    cost one extra pass over the data. Only the ztrunc points read by the
    inverse transform are stored, not the 2^depth of a full op.
*/
typedef struct
{
    double * data;
    ulong nops;
    ulong ntiles;
    ulong T;
} _tiles_struct;

static inline double * _tile(const _tiles_struct * R, ulong pi, ulong ti, ulong o)
{
    return R->data + ((pi * R->ntiles + ti) * R->nops + o) * R->T;
}

/* scatter the ztrunc points of each prime of X into the tiles of op o */
static void _tiles_scatter(const _tiles_struct * R, ulong o, const fft_small_op_t X,
                           const fft_small_plan_t P)
{
    ulong pi, ti;
    const ulong T = R->T;
    for (pi = 0; pi < P->np; pi++)
    {
        const double * src = X->data + pi * P->stride;
        for (ti = 0; ti < R->ntiles; ti++)
            memcpy(_tile(R, pi, ti, o), src + ti * T, T * sizeof(double));
    }
}

/* gather the tiles of op o into the first ztrunc points of each prime of X */
static void _tiles_gather(fft_small_op_t X, const _tiles_struct * R, ulong o,
                          const fft_small_plan_t P)
{
    ulong pi, ti;
    const ulong T = R->T;
    for (pi = 0; pi < P->np; pi++)
    {
        double * dst = X->data + pi * P->stride;
        for (ti = 0; ti < R->ntiles; ti++)
            memcpy(dst + ti * T, _tile(R, pi, ti, o), T * sizeof(double));
    }
}

/* forward transform of the entry pol (assumed nonzero) into the scratch
   op X, optionally scaled by the normalization factors, then scattered
   into the tiles of op o */
static void _transform_entry(const _tiles_struct * R, ulong o, fft_small_op_t X,
                             const nmod_poly_struct * pol, int scale,
                             nmod_t mod, const fft_small_plan_t P)
{
    fft_small_fft_nmod(X, pol->coeffs, pol->length,
                       _op_trunc(pol->length, P), mod, P);
    if (scale)
        _op_scale_by_m(X, P);
    _tiles_scatter(R, o, X, P);
}

/* ------------------------------------------------------------------------ */
/* parallel phases                                                          */
/* ------------------------------------------------------------------------ */

/*
    Each of the four phases -- transform B, transform a group of rows of
    A, pointwise products, inverse transforms with chinese remaindering
    -- is a loop over independent tasks: one entry of B, one entry of A,
    one (prime, tile) pair, one entry of C. They run over the FLINT
    thread pool with a shared descriptor holding the read-only context
    plus, per worker, a scratch transform and the pointer tables the
    kernels need.

    Threads are requested once for the whole product, so the fft_small
    entry points called inside a worker find none left and run serially.
    That is the intent: one transform per thread is coarser, and scales
    further, than their own threading over the (at most four) primes.
    When the matrix is too small for this to pay, no threads are
    requested and their per-prime threading applies as usual.
*/
typedef struct
{
    /* shared, read only */
    const fft_small_plan_struct * P;
    const _tiles_struct * At;
    const _tiles_struct * Bt;
    const _tiles_struct * Ct;
    const nmod_poly_mat_struct * A;
    const nmod_poly_mat_struct * B;
    nmod_poly_mat_struct * C;
    nmod_t mod;
    slong k;
    slong n;
    slong g;        /* first row of the current group of rows of A */
    slong nrows;    /* its size */
    slong h;        /* first column of the current group of columns of B */
    slong ncols;    /* its size */
    slong cstride;  /* NCG: column stride of the C and B blocks */
    ulong T;
    ulong ntiles;
    const slong * jlist;
    const slong * jcount;
    const slong * zlen;
    const slong * nzA;
    const slong * cntA;
    int Bdense;

    /* per worker */
    fft_small_op_t X;
    const double ** aptr;
    const double ** bptr;
    const double * const ** bcol;
    const double ** bfull;
    const double * const ** bfullcol;
    double ** zptr;
    slong start;
    slong stop;
}
_sd_fft_direct_worker_struct;

/* tasks o = jj*k + l: the transform of B[l][h+jj], scaled by the
   normalization factors (folded in here, once per product, rather than
   in each of the m pointwise passes that read it) */
static void _worker_fft_B(void * varg)
{
    _sd_fft_direct_worker_struct * W = (_sd_fft_direct_worker_struct *) varg;
    slong o;

    for (o = W->start; o < W->stop; o++)
    {
        const nmod_poly_struct * b =
            nmod_poly_mat_entry(W->B, o % W->k, W->h + o / W->k);
        if (b->length > 0)
            _transform_entry(W->Bt, o, W->X, b, 1, W->mod, W->P);
    }
}

/* tasks o = r*k + l: the transform of A[g+r][l] */
static void _worker_fft_A(void * varg)
{
    _sd_fft_direct_worker_struct * W = (_sd_fft_direct_worker_struct *) varg;
    slong o;

    for (o = W->start; o < W->stop; o++)
    {
        const nmod_poly_struct * a = nmod_poly_mat_entry(W->A, W->g + o / W->k, o % W->k);
        if (a->length > 0)
            _transform_entry(W->At, o, W->X, a, 0, W->mod, W->P);
    }
}

/* tasks o = pi*ntiles + ti: all the products of one tile of one prime.
   Distinct tasks write distinct tiles of C, so no synchronization is
   needed. */
static void _worker_pointwise(void * varg)
{
    _sd_fft_direct_worker_struct * W = (_sd_fft_direct_worker_struct *) varg;
    const fft_small_plan_struct * P = W->P;
    const slong k = W->k, nrows = W->nrows;
    const slong n = W->ncols;       /* columns in the current block */
    const slong cs = W->cstride;    /* their stride in the C and B blocks */
    const ulong T = W->T, ntiles = W->ntiles;
    slong o;

    for (o = W->start; o < W->stop; o++)
    {
        const ulong pi = (ulong) o / ntiles;
        const ulong ti = (ulong) o % ntiles;
        const sd_fft_ctx_struct * Q = P->ffts + P->offset + pi;
        const vec8d nn = vec8d_set_d(Q->p);
        const vec8d ninv = vec8d_set_d(Q->pinv);
        const double * Ablk = _tile(W->At, pi, ti, 0);
        const double * Bblk = _tile(W->Bt, pi, ti, 0);
        double * Cblk = _tile(W->Ct, pi, ti, 0);
        slong r, i, j, l;

        if (W->Bdense)
        {
            /* B without zero entries: every column is full for every
               row. Blocks of 8 columns of B outermost, so that their
               8*k tiles are reused from L1 across the rows of A. */
            for (j = 0; j < n; j += 8)
            {
                const slong nc = FLINT_MIN(8, n - j);
                slong c;

                for (c = 0; c < nc; c++)
                    for (l = 0; l < k; l++)
                        W->bfull[c * k + l] = Bblk + ((j + c) * k + l) * T;

                for (r = 0; r < nrows; r++)
                {
                    const slong ca = W->cntA[r];
                    const slong * lstA = W->nzA + r * k;
                    const double * const ** btab = W->bfullcol;

                    if (ca == 0)
                        continue;

                    for (i = 0; i < ca; i++)
                        W->aptr[i] = Ablk + (r * k + lstA[i]) * T;
                    if (ca < k)   /* A[g+r][l] == 0 for some l */
                    {
                        for (c = 0; c < nc; c++)
                            for (i = 0; i < ca; i++)
                                W->bptr[c * k + i] = Bblk + ((j + c) * k + lstA[i]) * T;
                        btab = W->bcol;
                    }
                    for (c = 0; c < nc; c++)
                        W->zptr[c] = Cblk + (r * cs + j + c) * T;

                    _dot_tile(W->zptr, W->aptr, btab, nc, ca, T, nn, ninv);
                }
            }
            continue;
        }

        for (r = 0; r < nrows; r++)
        {
            const slong ca = W->cntA[r];
            const slong * lstA = W->nzA + r * k;
            slong nfull = 0;

            if (ca == 0)
                continue;

            for (i = 0; i < ca; i++)
                W->aptr[i] = Ablk + (r * k + lstA[i]) * T;

            /* the columns whose list of nonzero products is the whole
               list of nonzero entries of the row: blocked together */
            for (j = 0; j < n; j++)
                if (W->jcount[r * cs + j] == ca)
                {
                    for (i = 0; i < ca; i++)
                        W->bptr[nfull * k + i] = Bblk + (j * k + lstA[i]) * T;
                    W->zptr[nfull++] = Cblk + (r * cs + j) * T;
                }
            if (nfull > 0)
                _dot_tile(W->zptr, W->aptr, W->bcol, nfull, ca, T, nn, ninv);

            /* the others, one at a time */
            for (j = 0; j < n; j++)
            {
                const slong cnt = W->jcount[r * cs + j];
                const slong * lst = W->jlist + ((ulong) r * cs + j) * k;
                if (cnt == 0 || cnt == ca)
                    continue;
                for (i = 0; i < cnt; i++)
                {
                    W->aptr[i] = Ablk + (r * k + lst[i]) * T;
                    W->bptr[i] = Bblk + (j * k + lst[i]) * T;
                }
                W->zptr[0] = Cblk + (r * cs + j) * T;
                _dot_tile(W->zptr, W->aptr, W->bcol, 1, cnt, T, nn, ninv);
            }
        }
    }
}

/* tasks o = r*ncols + jj: the entry C[g+r][h+jj] */
static void _worker_ifft(void * varg)
{
    _sd_fft_direct_worker_struct * W = (_sd_fft_direct_worker_struct *) varg;
    const slong cs = W->cstride, nc = W->ncols;
    slong o;

    for (o = W->start; o < W->stop; o++)
    {
        /* the block is ncols wide but strided by cstride */
        const slong oo = (o / nc) * cs + o % nc;
        nmod_poly_struct * c =
            nmod_poly_mat_entry(W->C, W->g + o / nc, W->h + o % nc);
        const slong zl = W->zlen[oo];

        if (W->jcount[oo] == 0)
        {
            nmod_poly_zero(c);
            continue;
        }

        _tiles_gather(W->X, W->Ct, oo, W->P);
        W->X->domain = FFT_SMALL_OP_PRODUCT;
        fft_small_ifft(W->X, W->P);
        nmod_poly_fit_length(c, zl);
        fft_small_export_nmod_range(c->coeffs, W->X, 0, zl, W->mod, W->P);
        _nmod_poly_set_length(c, zl);
        _nmod_poly_normalise(c);
    }
}

/* split [0, ntasks) over the first `nthreads` worker descriptors and run */
static void _sd_fft_direct_run(void (* func)(void *),
                            _sd_fft_direct_worker_struct * W, slong nthreads,
                            thread_pool_handle * handles, slong ntasks)
{
    slong i;

    if (ntasks < 1)
        return;

    nthreads = FLINT_MIN(nthreads, ntasks);

    for (i = 0; i < nthreads; i++)
    {
        W[i].start = (i + 0) * ntasks / nthreads;
        W[i].stop  = (i + 1) * ntasks / nthreads;
    }

    for (i = nthreads - 1; i > 0; i--)
        thread_pool_wake(global_thread_pool, handles[i - 1], 0, func, W + i);
    func(W + 0);
    for (i = nthreads - 1; i > 0; i--)
        thread_pool_wait(global_thread_pool, handles[i - 1]);
}

void nmod_poly_mat_mul_sd_fft_direct(nmod_poly_mat_t C,
                                  const nmod_poly_mat_t A,
                                  const nmod_poly_mat_t B)
{
    const slong m = A->r;
    const slong k = A->c;
    const slong n = B->c;
    const slong lenA = nmod_poly_mat_max_length(A);
    const slong lenB = nmod_poly_mat_max_length(B);

    if (m == 0 || n == 0)
        return;

    if (k == 0 || lenA == 0 || lenB == 0)
    {
        nmod_poly_mat_zero(C);
        return;
    }

    if (C == A || C == B)
    {
        nmod_poly_mat_t T;
        nmod_poly_mat_init(T, m, n, A->modulus);
        nmod_poly_mat_mul_sd_fft_direct(T, A, B);
        nmod_poly_mat_swap_entrywise(C, T);
        nmod_poly_mat_clear(T);
        return;
    }

    nmod_t mod;
    nmod_init(&mod, A->modulus);

    /* plan: output window is the full product; at most
       k * min(lenA, lenB) products of two residues accumulate onto one
       output coefficient. The direct single-prime transform modulo p is
       forced whenever p allows it (direct_len = UWORD_MAX): its setup
       cost, proportional to the transform length, is amortized over the
       m*k + k*n + m*n transforms performed here, unlike in a single
       polynomial product. */
    const ulong zn = (ulong) lenA + lenB - 1;
    const ulong xtrunc_max = n_max(n_round_up(lenA, BLK_SZ), n_round_up(lenB, BLK_SZ));
    const ulong len_bound = (ulong) k * n_min(lenA, lenB);
    mpn_ctx_struct * R = get_default_mpn_ctx();
    fft_small_plan_t P;

    if (!fft_small_plan_init_nmod(P, R, 0, zn, zn, xtrunc_max, len_bound,
                                  2 * NMOD_BITS(mod), mod, UWORD_MAX))
    {
        /* bound beyond the capacity of the eight primes: needs
           k * min(lenA, lenB) > 2^(400 - 2*64), which no representable
           matrix reaches */
        nmod_poly_mat_mul(C, A, B);
        return;
    }

    TIMING_DECL

    const ulong np = P->np;
    const ulong ztrunc = P->ztrunc;
    const ulong scratchdbls = fft_small_op_sizeof_data(P) / sizeof(double);
    const ulong opdbls = np * ztrunc;   /* tile-major footprint of one transform */

    /* tile: T points, T | ztrunc, 16 <= T <= BLK_SZ, 8*k*T*8 within budget */
    ulong T = n_min(ztrunc, BLK_SZ);
#ifdef PML_MUL_SD_FFT_DIRECT_TIMING
    ulong tile_budget = getenv("PML_TILE_BUDGET")
        ? strtoul(getenv("PML_TILE_BUDGET"), NULL, 10)
        : NMOD_POLY_MAT_MUL_SD_FFT_DIRECT_TILE_BUDGET;
    if (getenv("PML_T")) { T = strtoul(getenv("PML_T"), NULL, 10); tile_budget = UWORD_MAX; }
#else
    const ulong tile_budget = NMOD_POLY_MAT_MUL_SD_FFT_DIRECT_TILE_BUDGET;
#endif
    while (T > 16 && 8 * (ulong) k * T * sizeof(double) > tile_budget)
        T /= 2;
    const ulong ntiles = ztrunc / T;

    /* threads: one task is one transform or one (prime, tile) pair, so
       there is little to gain below a few entries; there, letting the
       fft_small calls thread over the primes themselves is better */
    thread_pool_handle * handles = NULL;
    slong nworkers = 0;
    slong nthreads = 1;

    if ((ulong) m * k * n >= 4
        && (ulong) m * k * n * np * ztrunc >= (UWORD(1) << 20))
    {
        nworkers = flint_request_threads(&handles, flint_get_num_threads());
        nthreads = nworkers + 1;
    }

    /* Grouping. NRG rows of A and C, NCG columns of B and C are held at
       a time, so that the storage
           k*NCG (B) + NRG*k (A) + NRG*NCG (C) + nthreads scratch
       transforms stays within the budget. Grouping rows is free, so it
       is used first and as much as needed; grouping columns costs one
       transform of A per group of columns, so it is used only when the
       transforms of B alone (plus one row) do not fit. */
    /* words of the operands and of the result, as the scale of the bound
       on the transforms (see NMOD_POLY_MAT_MUL_SD_FFT_DIRECT_MEM_FLOOR) */
    const ulong opwords = (ulong) m * k * lenA + (ulong) k * n * lenB
                          + (ulong) m * n * zn;
    ulong bbytes = n_max(NMOD_POLY_MAT_MUL_SD_FFT_DIRECT_MEM_FLOOR,
                         2 * opwords * sizeof(ulong));
#ifdef PML_MUL_SD_FFT_DIRECT_TIMING
    if (getenv("PML_MEM_BUDGET"))
        bbytes = strtoul(getenv("PML_MEM_BUDGET"), NULL, 10);
#endif

    slong NRG, NCG;
    {
        const ulong budget = bbytes / (opdbls * sizeof(double));
        const ulong scratch = ((ulong) nthreads * scratchdbls + opdbls - 1) / opdbls;
        const ulong avail = (budget > scratch + 2) ? budget - scratch : 2;

        NCG = n;
        if ((ulong) k * n + (ulong) (k + n) > avail)
        {
            /* give about half of what is left to the transforms of B */
            NCG = (slong) (avail / (2 * (ulong) k + 1));
            NCG = FLINT_MAX(NCG, 1);
            NCG = FLINT_MIN(NCG, n);
        }

        {
            const ulong fixed = (ulong) k * NCG;
            const ulong per_row = (ulong) k + NCG;
            NRG = (avail > fixed) ? (slong) ((avail - fixed) / per_row) : 1;
            NRG = FLINT_MAX(NRG, 1);
            NRG = FLINT_MIN(NRG, m);
        }
    }

    /* storage: the per-thread scratch ops first, since they are the only
       part with an alignment requirement (fft_small_op_init_borrowed
       wants 4096 bytes, and fft_small_op_sizeof_data is a multiple of
       that, so the tile regions that follow are aligned as well), then
       the tiles of B (k*NCG transforms), of the current rows of A
       (NRG*k) and of C (NRG*NCG) */
    const ulong nbytes = ((ulong) nthreads * scratchdbls
                          + ((ulong) k * NCG + (ulong) NRG * (k + NCG)) * opdbls)
                         * sizeof(double);
    double * buf = _sd_fft_direct_alloc(R, nbytes, bbytes);
    double * tilebuf = buf + (ulong) nthreads * scratchdbls;
    _tiles_struct Bt = { tilebuf, (ulong) k * NCG, ntiles, T };
    _tiles_struct At = { Bt.data + Bt.nops * opdbls, (ulong) NRG * k, ntiles, T };
    _tiles_struct Ct = { At.data + At.nops * opdbls, (ulong) NRG * NCG, ntiles, T };

    /* structure of the products, for the current block of rows and
       columns: nzA[r*k ..] lists the l with A[g+r][l] != 0 (cntA[r] of
       them); for each (r, jj), jlist[(r*NCG+jj)*k ..] lists the l with
       A[g+r][l] and B[l][h+jj] both nonzero (jcount[r*NCG+jj] of them),
       and zlen[r*NCG+jj] is the length of C[g+r][h+jj] */
    slong * jlist = FLINT_ARRAY_ALLOC((ulong) NRG * NCG * k, slong);
    slong * jcount = FLINT_ARRAY_ALLOC((ulong) NRG * NCG, slong);
    slong * zlen = FLINT_ARRAY_ALLOC((ulong) NRG * NCG, slong);
    slong * nzA = FLINT_ARRAY_ALLOC((ulong) NRG * k, slong);
    slong * cntA = FLINT_ARRAY_ALLOC(NRG, slong);

    _sd_fft_direct_worker_struct * W =
        FLINT_ARRAY_ALLOC(nthreads, _sd_fft_direct_worker_struct);
    {
        slong w, j;
        for (w = 0; w < nthreads; w++)
        {
            W[w].P = P;
            W[w].At = &At; W[w].Bt = &Bt; W[w].Ct = &Ct;
            W[w].A = A; W[w].B = B; W[w].C = C;
            W[w].mod = mod;
            W[w].k = k; W[w].n = n;
            W[w].g = 0; W[w].nrows = 0;
            W[w].h = 0; W[w].ncols = 0; W[w].cstride = NCG;
            W[w].T = T; W[w].ntiles = ntiles;
            W[w].jlist = jlist; W[w].jcount = jcount; W[w].zlen = zlen;
            W[w].nzA = nzA; W[w].cntA = cntA;
            W[w].Bdense = 0;

            fft_small_op_init_borrowed(W[w].X, P, buf + (ulong) w * scratchdbls);
            /* the transforms only write the ztrunc points they use, but
               sd_ifft_trunc works on the whole 2^depth array; zeroing
               once keeps every value it ever sees a finite residue
               (after this, each pass leaves bounded values behind) */
            {
                ulong pi;
                for (pi = 0; pi < np; pi++)
                    memset(W[w].X->data + pi * P->stride + ztrunc, 0,
                           (sd_fft_ctx_data_size(P->depth) - ztrunc) * sizeof(double));
            }

            W[w].aptr = FLINT_ARRAY_ALLOC(k, const double *);
            W[w].bptr = FLINT_ARRAY_ALLOC((ulong) NCG * k, const double *);
            W[w].bcol = FLINT_ARRAY_ALLOC(NCG, const double * const *);
            W[w].bfull = FLINT_ARRAY_ALLOC(8 * (ulong) k, const double *);
            W[w].bfullcol = FLINT_ARRAY_ALLOC(8, const double * const *);
            W[w].zptr = FLINT_ARRAY_ALLOC(NCG, double *);
            for (j = 0; j < NCG; j++)
                W[w].bcol[j] = W[w].bptr + j * k;
            for (j = 0; j < 8; j++)
                W[w].bfullcol[j] = W[w].bfull + j * k;
        }
    }

    slong g, h;
    for (h = 0; h < n; h += NCG)
    {
        const slong ncols = FLINT_MIN(NCG, n - h);
        slong r, j, l, w;
        int Bdense = 1;

        /* a block of B without zero entries has every column full for
           every row of A, which lets the pointwise stage put the
           columns of B outermost (see _worker_pointwise) */
        for (l = 0; l < k; l++)
            for (j = 0; j < ncols; j++)
                if (nmod_poly_mat_entry(B, l, h + j)->length == 0)
                    Bdense = 0;

        for (w = 0; w < nthreads; w++)
        {
            W[w].h = h;
            W[w].ncols = ncols;
            W[w].Bdense = Bdense;
        }

        /* B: op jj*k + l <- transform of B[l][h+jj] (column major, so
           that the transforms of a column of B are consecutive) */
        _sd_fft_direct_run(_worker_fft_B, W, nthreads, handles, ncols * k);

        TIMING_MARK(_tB);

        for (g = 0; g < m; g += NRG)
        {
            const slong nrows = FLINT_MIN(NRG, m - g);

            for (w = 0; w < nthreads; w++)
            {
                W[w].g = g;
                W[w].nrows = nrows;
            }

            _sd_fft_direct_run(_worker_fft_A, W, nthreads, handles, nrows * k);

            for (r = 0; r < nrows; r++)
            {
                slong cnt = 0;
                for (l = 0; l < k; l++)
                    if (nmod_poly_mat_entry(A, g + r, l)->length > 0)
                        nzA[r * k + cnt++] = l;
                cntA[r] = cnt;

                for (j = 0; j < ncols; j++)
                {
                    slong zl = 0;
                    slong * lst = jlist + ((ulong) r * NCG + j) * k;
                    cnt = 0;
                    for (l = 0; l < k; l++)
                    {
                        const slong la = nmod_poly_mat_entry(A, g + r, l)->length;
                        const slong lb = nmod_poly_mat_entry(B, l, h + j)->length;
                        if (la > 0 && lb > 0)
                        {
                            lst[cnt++] = l;
                            zl = FLINT_MAX(zl, la + lb - 1);
                        }
                    }
                    jcount[r * NCG + j] = cnt;
                    zlen[r * NCG + j] = zl;
                }
            }
            TIMING_MARK(_tA);

            _sd_fft_direct_run(_worker_pointwise, W, nthreads, handles,
                            (slong) (np * ntiles));

            TIMING_MARK(_tP);

            _sd_fft_direct_run(_worker_ifft, W, nthreads, handles, nrows * ncols);

            TIMING_MARK(_tC);
        }
    }

    TIMING_PRINT;

    {
        slong w;
        for (w = 0; w < nthreads; w++)
        {
            flint_free(W[w].aptr);
            flint_free(W[w].bptr);
            flint_free(W[w].bcol);
            flint_free(W[w].bfull);
            flint_free(W[w].bfullcol);
            flint_free(W[w].zptr);
        }
    }
    flint_free(W);
    flint_free(jlist);
    flint_free(jcount);
    flint_free(zlen);
    flint_free(nzA);
    flint_free(cntA);
    _sd_fft_direct_free(buf, nbytes, bbytes);
    flint_give_back_threads(handles, nworkers);
    fft_small_plan_clear(P);
}

#else  /* PML_HAVE_MACHINE_VECTORS */

void nmod_poly_mat_mul_sd_fft_direct(nmod_poly_mat_t C,
                                  const nmod_poly_mat_t A,
                                  const nmod_poly_mat_t B)
{
    nmod_poly_mat_mul(C, A, B);
}

#endif  /* PML_HAVE_MACHINE_VECTORS */
