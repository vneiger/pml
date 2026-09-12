/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <string.h>  /* memset */

#include <flint/nmod.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>
#include <flint/ulong_extras.h>

#include "nmod_poly_mat_multiply.h"

#if PML_HAVE_MACHINE_VECTORS

#include <flint/fft_small.h>
#include <flint/thread_pool.h>
#include <flint/thread_support.h>

#include "machine_vectors.h"

/* timing breakdown of the phases, printed on stderr (development aid) */
#ifdef PML_MUL_SD_FFT_MATMUL_TIMING
# include <stdio.h>
# include <stdlib.h>
# include <time.h>
static double _now(void) { struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return ts.tv_sec + 1e-9 * ts.tv_nsec; }
# define TIMING_DECL double _t0 = _now(), _t1; double _tB = 0, _tA = 0, _tP = 0, _tC = 0;
# define TIMING_MARK(acc) do { _t1 = _now(); acc += _t1 - _t0; _t0 = _t1; } while (0)
# define TIMING_PRINT flint_fprintf(stderr, "  [sd_fft_matmul] np=%wu depth=%wu ztrunc=%wu NRG=%wd NCG=%wd | B fft %.3f, A fft %.3f, matmul %.3f, ifft+crt %.3f\n", np, P->depth, ztrunc, NRG, NCG, _tB, _tA, _tP, _tC)
#else
# define TIMING_DECL
# define TIMING_MARK(acc)
# define TIMING_PRINT
#endif

/*
    Polynomial matrix multiplication C = A * B by evaluation-interpolation
    with FLINT's fft_small transforms, calling nmod_mat_mul on the
    evaluations -- the counterpart of nmod_poly_mat_mul_sd_fft_direct,
    which instead runs a cubic kernel directly on the transforms, and the
    fft_small analogue of what nmod_poly_mat_mul_geometric does with a
    geometric progression of evaluation points.

        1. plan: exactly as in mul_sd_fft_direct. The FFT depth comes from
           lenA + lenB - 1 and the set of primes from the bound on the
           output coefficients: before reduction modulo p these are sums
           of at most k * min(lenA, lenB) products of two residues, so np
           primes with prod(primes) >= k * min(lenA, lenB) * p^2 are
           needed (1 to 4 of the 50-bit fft_small primes; a single direct
           transform modulo p when p is itself a suitable FFT prime).
        2. transform every entry of A and of B, and store the residues
           point major (see _pts_struct): for each prime and each
           evaluation point, the m*k values of A form one contiguous
           m x k matrix over Z/p_i, and likewise for B and C.
        3. for each prime and each point, one call to nmod_mat_mul on
           borrowed storage -- no copy in or out.
        4. inverse transform and chinese remainder each entry of C.

    Sparsity is not exploited, unlike in mul_sd_fft_direct: zero entries
    of A and B are materialised as zero residues, since nmod_mat_mul
    reads whole matrices.
*/

/*
    Soft bound, in bytes, on the memory used for the residues. Same role
    as in mul_sd_fft_direct, chosen differently in two ways.

    The shape of a group matters here: it is the dimension of the
    matrices handed to nmod_mat_mul, and that routine changes algorithm
    on the smallest of the three (Strassen above 200, BLAS above 450 in
    FLINT's current tuning). So the two group sizes are balanced rather
    than the rows of A alone being cut, and a group that is large enough
    to reach those algorithms is worth a good deal: measured at dimension
    512, length 2*128, modulo a 50-bit FFT prime, the pointwise stage
    takes 25.2 s in groups of 115 (a 256MB bound), 23.1 s in groups of
    374 (1GB) and 14.5 s undivided (3GB).

    The bound is therefore not a constant but scales with the problem:
    the transforms of a product whose operands are themselves a gigabyte
    have no business being capped at a few hundred megabytes. The
    constant below is only a floor, for the small products where a fixed
    working set is what one wants to bound. The factor 2 is about what
    an undivided square product of a generic modulus needs relative to
    its operands when a single prime suffices; with 3 or 4 primes the
    grouping still applies, and the balance above still puts the largest
    possible square blocks in front of nmod_mat_mul.

    The residues are taken from the retained fft_small scratch buffer
    for any request within that bound, as in mul_sd_fft_direct: the
    first-touch faults of a fresh mapping of this size are not lost in
    the arithmetic (at the dimension 512 point above, a freshly mapped
    1.6GB costs 5 s of the 20, since the faults land inside the product
    loops), and PML's algorithms perform sequences of products of
    similar shape, which pay them once. The buffer is thread local and
    is not given back, so a thread that has done one such product keeps
    a working set of the order of twice the operands it was given.
*/
#define NMOD_POLY_MAT_MUL_SD_FFT_MATMUL_MEM_FLOOR (UWORD(1) << 28)

/* see the same allocator in mul_sd_fft_direct.c */
static ulong * _matmul_alloc(mpn_ctx_struct * R, ulong nbytes, ulong budget)
{
    if (nbytes <= budget)
        return (ulong *) mpn_ctx_fit_buffer(R, nbytes);
    return (ulong *) flint_aligned_alloc(FLINT_FFT_SMALL_ALIGNMENT,
                               n_round_up(nbytes, FLINT_FFT_SMALL_ALIGNMENT));
}

static void _matmul_free(ulong * buf, ulong nbytes, ulong budget)
{
    if (nbytes > budget)
        flint_aligned_free(buf);
}

/* transform-space truncation of an operand of length an */
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
    results in (-n, n).
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

/* ------------------------------------------------------------------------ */
/* point-major residues                                                     */
/* ------------------------------------------------------------------------ */

/*
    The residues of one region (the block of B, of A or of C currently
    held) modulo one prime at one evaluation point form one matrix, laid
    out contiguously and row major so that an nmod_mat can borrow it:

        region + (pi * ztrunc + t) * stride,   stride = nops rounded up to 8,

    where nops is the number of entries of the region (m*k, k*n or m*n for
    the current groups) and the row stride inside the matrix is its number
    of columns. Rounding the distance between consecutive points up to 8
    words, on a 64-byte aligned base, keeps each group of 8 consecutive
    entries of one point inside a single cache line: the transposition
    below writes the array one such group at a time, and two threads never
    share a line.
*/
typedef struct
{
    ulong * data;
    ulong nops;
    ulong stride;
    ulong ztrunc;
}
_pts_struct;

static inline ulong * _point(const _pts_struct * R, ulong pi, ulong t)
{
    return R->data + (pi * R->ztrunc + t) * R->stride;
}

/* number of ulongs of a region of nops entries */
static inline ulong _pts_words(ulong nops, ulong np, ulong ztrunc)
{
    return np * ztrunc * n_round_up(nops, 8);
}

/*
    Scatter the ztrunc points of each prime of X into entry o of the
    point-major region, reduced into [0, p).

    The stores are one word per point, i.e. one cache line per point,
    which is why the entries are visited in increasing o: the ztrunc
    lines written for entry o are the ones written for entries o+1, ...,
    o+7, so a group of eight entries fills them and (for the sizes where
    this routine matters) they stay in L2 meanwhile. The loads and the
    arithmetic are vectorised; the stores cannot be, short of buffering
    eight whole transforms.
*/
static void _pts_scatter(const _pts_struct * R, ulong o, const fft_small_op_t X,
                         const fft_small_plan_t P)
{
    const ulong ztrunc = R->ztrunc, str = R->stride;
    ulong pi, t;

    for (pi = 0; pi < P->np; pi++)
    {
        const sd_fft_ctx_struct * Q = P->ffts + P->offset + pi;
        const vec4d n = vec4d_set_d(Q->p);
        const vec4d ninv = vec4d_set_d(Q->pinv);
        const double * src = X->data + pi * P->stride;
        ulong * dst = _point(R, pi, 0) + o;

        for (t = 0; t < ztrunc; t += 4)
        {
            ulong buf[4];
            vec4d_store_unaligned_nn_ptr(buf,
                vec4d_reduce_to_0n(vec4d_load(src + t), n, ninv));
            dst[(t + 0) * str] = buf[0];
            dst[(t + 1) * str] = buf[1];
            dst[(t + 2) * str] = buf[2];
            dst[(t + 3) * str] = buf[3];
        }
    }
}

/* the same for an entry which is the zero polynomial */
static void _pts_scatter_zero(const _pts_struct * R, ulong o, const fft_small_plan_t P)
{
    const ulong ztrunc = R->ztrunc, str = R->stride;
    ulong pi, t;

    for (pi = 0; pi < P->np; pi++)
    {
        ulong * dst = _point(R, pi, 0) + o;
        for (t = 0; t < ztrunc; t++)
            dst[t * str] = 0;
    }
}

/* gather entry o of the point-major region into the first ztrunc points
   of each prime of X, as doubles */
static void _pts_gather(fft_small_op_t X, const _pts_struct * R, ulong o,
                        const fft_small_plan_t P)
{
    const ulong ztrunc = R->ztrunc, str = R->stride;
    ulong pi, t;

    for (pi = 0; pi < P->np; pi++)
    {
        double * dst = X->data + pi * P->stride;
        const ulong * src = _point(R, pi, 0) + o;

        for (t = 0; t < ztrunc; t += 4)
        {
            ulong buf[4];
            buf[0] = src[(t + 0) * str];
            buf[1] = src[(t + 1) * str];
            buf[2] = src[(t + 2) * str];
            buf[3] = src[(t + 3) * str];
            vec4d_store(dst + t, vec4d_load_unaligned_nn_ptr(buf));
        }
    }
}

/* forward transform of the entry pol (assumed nonzero) into the scratch
   op X, optionally scaled by the normalization factors, then scattered */
static void _transform_entry(const _pts_struct * R, ulong o, fft_small_op_t X,
                             const nmod_poly_struct * pol, int scale,
                             nmod_t mod, const fft_small_plan_t P)
{
    fft_small_fft_nmod(X, pol->coeffs, pol->length,
                       _op_trunc(pol->length, P), mod, P);
    if (scale)
        _op_scale_by_m(X, P);
    _pts_scatter(R, o, X, P);
}

/* ------------------------------------------------------------------------ */
/* parallel phases                                                          */
/* ------------------------------------------------------------------------ */

/*
    Four phases -- transform B, transform a group of rows of A, the matrix
    products, inverse transforms with chinese remaindering -- each a loop
    over independent tasks run over the FLINT thread pool. Threads are
    requested once for the whole product, so the fft_small entry points
    and nmod_mat_mul called inside a worker find none left and run
    serially: one task per thread is coarser, and scales further, than
    their own threading.

    Tasks of the two transform phases are groups of eight entries, so that
    the cache lines of the point-major arrays are not shared between
    threads (see _pts_struct).
*/
typedef struct
{
    /* shared, read only */
    const fft_small_plan_struct * P;
    const _pts_struct * At;
    const _pts_struct * Bt;
    const _pts_struct * Ct;
    const nmod_poly_mat_struct * A;
    const nmod_poly_mat_struct * B;
    nmod_poly_mat_struct * C;
    nmod_t mod;
    slong k;
    slong g;        /* first row of the current group of rows of A */
    slong nrows;    /* its size */
    slong h;        /* first column of the current group of columns of B */
    slong ncols;    /* its size */
    slong cstride;  /* NCG: column stride of the C and B blocks */
    ulong ztrunc;
    const slong * zlen;

    /* per worker */
    fft_small_op_t X;
    slong start;
    slong stop;
}
_matmul_worker_struct;

/* tasks are groups of 8 of the entries o = l*NCG + jj of the block of B,
   holding the transform of B[l][h+jj] scaled by the normalization factors
   (folded in here, once per product) */
static void _worker_fft_B(void * varg)
{
    _matmul_worker_struct * W = (_matmul_worker_struct *) varg;
    const slong cs = W->cstride, nc = W->ncols;
    const ulong nops = W->Bt->nops;
    slong o;

    for (o = 8 * W->start; o < 8 * W->stop && (ulong) o < nops; o++)
    {
        const slong l = o / cs, jj = o % cs;
        const nmod_poly_struct * b;

        if (jj >= nc)   /* padding column of the block: never read */
            continue;

        b = nmod_poly_mat_entry(W->B, l, W->h + jj);
        if (b->length > 0)
            _transform_entry(W->Bt, (ulong) o, W->X, b, 1, W->mod, W->P);
        else
            _pts_scatter_zero(W->Bt, (ulong) o, W->P);
    }
}

/* tasks are groups of 8 of the entries o = r*k + l, holding the transform
   of A[g+r][l] */
static void _worker_fft_A(void * varg)
{
    _matmul_worker_struct * W = (_matmul_worker_struct *) varg;
    const slong k = W->k;
    const ulong nops = (ulong) W->nrows * k;
    slong o;

    for (o = 8 * W->start; o < 8 * W->stop && (ulong) o < nops; o++)
    {
        const nmod_poly_struct * a =
            nmod_poly_mat_entry(W->A, W->g + o / k, o % k);
        if (a->length > 0)
            _transform_entry(W->At, (ulong) o, W->X, a, 0, W->mod, W->P);
        else
            _pts_scatter_zero(W->At, (ulong) o, W->P);
    }
}

/* tasks o = pi*ztrunc + t: the matrix product at one point of one prime.
   Distinct tasks write distinct matrices of C, so no synchronization is
   needed. */
static void _worker_matmul(void * varg)
{
    _matmul_worker_struct * W = (_matmul_worker_struct *) varg;
    const fft_small_plan_struct * P = W->P;
    const ulong ztrunc = W->ztrunc;
    slong o;

    for (o = W->start; o < W->stop; o++)
    {
        const ulong pi = (ulong) o / ztrunc;
        const ulong t = (ulong) o % ztrunc;
        const nmod_t pmod = P->ffts[P->offset + pi].mod;
        nmod_mat_struct MA, MB, MC;

        MA.entries = _point(W->At, pi, t);
        MA.r = W->nrows; MA.c = W->k; MA.stride = W->k; MA.mod = pmod;

        MB.entries = _point(W->Bt, pi, t);
        MB.r = W->k; MB.c = W->ncols; MB.stride = W->cstride; MB.mod = pmod;

        MC.entries = _point(W->Ct, pi, t);
        MC.r = W->nrows; MC.c = W->ncols; MC.stride = W->cstride; MC.mod = pmod;

        nmod_mat_mul(&MC, &MA, &MB);
    }
}

/* tasks o = r*ncols + jj: the entry C[g+r][h+jj] */
static void _worker_ifft(void * varg)
{
    _matmul_worker_struct * W = (_matmul_worker_struct *) varg;
    const slong cs = W->cstride, nc = W->ncols;
    slong o;

    for (o = W->start; o < W->stop; o++)
    {
        /* the block is ncols wide but strided by cstride */
        const slong oo = (o / nc) * cs + o % nc;
        nmod_poly_struct * c =
            nmod_poly_mat_entry(W->C, W->g + o / nc, W->h + o % nc);
        const slong zl = W->zlen[oo];

        if (zl == 0)
        {
            nmod_poly_zero(c);
            continue;
        }

        _pts_gather(W->X, W->Ct, (ulong) oo, W->P);
        W->X->domain = FFT_SMALL_OP_PRODUCT;
        fft_small_ifft(W->X, W->P);
        nmod_poly_fit_length(c, zl);
        fft_small_export_nmod_range(c->coeffs, W->X, 0, zl, W->mod, W->P);
        _nmod_poly_set_length(c, zl);
        _nmod_poly_normalise(c);
    }
}

/* split [0, ntasks) over the first `nthreads` worker descriptors and run */
static void _matmul_run(void (* func)(void *),
                        _matmul_worker_struct * W, slong nthreads,
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

void nmod_poly_mat_mul_sd_fft_matmul(nmod_poly_mat_t C,
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
        nmod_poly_mat_mul_sd_fft_matmul(T, A, B);
        nmod_poly_mat_swap_entrywise(C, T);
        nmod_poly_mat_clear(T);
        return;
    }

    nmod_t mod;
    nmod_init(&mod, A->modulus);

    const ulong zn = (ulong) lenA + lenB - 1;
    const ulong xtrunc_max = n_max(n_round_up(lenA, BLK_SZ), n_round_up(lenB, BLK_SZ));
    const ulong len_bound = (ulong) k * n_min(lenA, lenB);
    mpn_ctx_struct * R = get_default_mpn_ctx();
    fft_small_plan_t P;

    if (!fft_small_plan_init_nmod(P, R, 0, zn, zn, xtrunc_max, len_bound,
                                  2 * NMOD_BITS(mod), mod, UWORD_MAX))
    {
        nmod_poly_mat_mul(C, A, B);
        return;
    }

    TIMING_DECL

    const ulong np = P->np;
    const ulong ztrunc = P->ztrunc;
    const ulong scratchdbls = fft_small_op_sizeof_data(P) / sizeof(double);

    /* threads: one task is a group of transforms or one matrix product */
    thread_pool_handle * handles = NULL;
    slong nworkers = 0;
    slong nthreads = 1;

    if ((ulong) m * k * n >= 4
        && (ulong) m * k * n * np * ztrunc >= (UWORD(1) << 20))
    {
        nworkers = flint_request_threads(&handles, flint_get_num_threads());
        nthreads = nworkers + 1;
    }

    /*
        Grouping. NRG rows of A and C, NCG columns of B and C are held at
        a time, within the memory budget. Unlike in mul_sd_fft_direct the
        two are balanced: they are the outer dimensions of the matrices
        handed to nmod_mat_mul, and a flat block is worth less to it than
        a square one of the same area (it reaches neither Strassen nor
        BLAS, whose thresholds are on the smallest dimension).
    */
    /* words of the operands and of the result, as a scale for the bound
       on the residues (see NMOD_POLY_MAT_MUL_SD_FFT_MATMUL_MEM_FLOOR) */
    const ulong opwords = (ulong) m * k * lenA + (ulong) k * n * lenB
                          + (ulong) m * n * zn;
    ulong bbytes = n_max(NMOD_POLY_MAT_MUL_SD_FFT_MATMUL_MEM_FLOOR,
                         2 * opwords * sizeof(ulong));
#ifdef PML_MUL_SD_FFT_MATMUL_TIMING
    if (getenv("PML_MEM_BUDGET"))
        bbytes = strtoul(getenv("PML_MEM_BUDGET"), NULL, 10);
#endif

    slong NRG, NCG;
    {
        /* budget and requirements counted in transforms of np*ztrunc words */
        const ulong budget = bbytes / (np * ztrunc * sizeof(ulong));
        const ulong scratch = ((ulong) nthreads * scratchdbls + np * ztrunc - 1)
                              / (np * ztrunc);
        const ulong avail = (budget > scratch + 3) ? budget - scratch : 3;

        NRG = m;
        NCG = n;
        if ((ulong) k * n + (ulong) m * k + (ulong) m * n > avail)
        {
            /* NRG = NCG = g uses 2*k*g + g^2 transforms */
            slong g = (slong) (n_sqrt((ulong) k * k + avail) - (ulong) k);
            g = FLINT_MAX(g, 1);
            NRG = FLINT_MIN(m, g);
            /* the largest NCG for that NRG, and then the largest NRG for
               that NCG: this gives back to one side what the other could
               not use because m or n was the binding cap */
            NCG = (avail > (ulong) NRG * k)
                  ? (slong) ((avail - (ulong) NRG * k) / ((ulong) k + NRG)) : 1;
            NCG = FLINT_MAX(NCG, 1);
            NCG = FLINT_MIN(NCG, n);
            NRG = (avail > (ulong) k * NCG)
                  ? (slong) ((avail - (ulong) k * NCG) / ((ulong) k + NCG)) : 1;
            NRG = FLINT_MAX(NRG, 1);
            NRG = FLINT_MIN(NRG, m);
        }
    }

    /* storage: the per-thread scratch ops first (they are the only part
       with an alignment requirement, and fft_small_op_sizeof_data is a
       multiple of it, so the regions that follow are aligned too), then
       the residues of B (k x NCG), of the current rows of A (NRG x k) and
       of C (NRG x NCG) */
    const ulong wordsB = _pts_words((ulong) k * NCG, np, ztrunc);
    const ulong wordsA = _pts_words((ulong) NRG * k, np, ztrunc);
    const ulong wordsC = _pts_words((ulong) NRG * NCG, np, ztrunc);
    const ulong nbytes = ((ulong) nthreads * scratchdbls
                          + wordsB + wordsA + wordsC) * sizeof(ulong);
    ulong * buf = _matmul_alloc(R, nbytes, bbytes);
    ulong * ptsbuf = buf + (ulong) nthreads * scratchdbls;

    _pts_struct Bt = { ptsbuf, (ulong) k * NCG,
                       n_round_up((ulong) k * NCG, 8), ztrunc };
    _pts_struct At = { ptsbuf + wordsB, (ulong) NRG * k,
                       n_round_up((ulong) NRG * k, 8), ztrunc };
    _pts_struct Ct = { ptsbuf + wordsB + wordsA, (ulong) NRG * NCG,
                       n_round_up((ulong) NRG * NCG, 8), ztrunc };

    /* zlen[r*NCG+jj] is the length of C[g+r][h+jj], zero when that entry
       has no nonzero product contributing to it */
    slong * zlen = FLINT_ARRAY_ALLOC((ulong) NRG * NCG, slong);

    _matmul_worker_struct * W =
        FLINT_ARRAY_ALLOC(nthreads, _matmul_worker_struct);
    {
        slong w;
        for (w = 0; w < nthreads; w++)
        {
            W[w].P = P;
            W[w].At = &At; W[w].Bt = &Bt; W[w].Ct = &Ct;
            W[w].A = A; W[w].B = B; W[w].C = C;
            W[w].mod = mod;
            W[w].k = k;
            W[w].g = 0; W[w].nrows = 0;
            W[w].h = 0; W[w].ncols = 0; W[w].cstride = NCG;
            W[w].ztrunc = ztrunc;
            W[w].zlen = zlen;

            fft_small_op_init_borrowed(W[w].X, P,
                                       (double *) (buf + (ulong) w * scratchdbls));
            /* the transforms only write the ztrunc points they use, but
               sd_ifft_trunc works on the whole 2^depth array; zeroing
               once keeps every value it ever sees a finite residue */
            {
                ulong pi;
                for (pi = 0; pi < np; pi++)
                    memset(W[w].X->data + pi * P->stride + ztrunc, 0,
                           (sd_fft_ctx_data_size(P->depth) - ztrunc) * sizeof(double));
            }
        }
    }

    slong g, h;
    for (h = 0; h < n; h += NCG)
    {
        const slong ncols = FLINT_MIN(NCG, n - h);
        slong r, j, l, w;

        for (w = 0; w < nthreads; w++)
        {
            W[w].h = h;
            W[w].ncols = ncols;
        }

        /* B: entry l*NCG + jj <- transform of B[l][h+jj] */
        _matmul_run(_worker_fft_B, W, nthreads, handles,
                    (slong) ((Bt.nops + 7) / 8));

        TIMING_MARK(_tB);

        for (g = 0; g < m; g += NRG)
        {
            const slong nrows = FLINT_MIN(NRG, m - g);

            for (w = 0; w < nthreads; w++)
            {
                W[w].g = g;
                W[w].nrows = nrows;
            }

            _matmul_run(_worker_fft_A, W, nthreads, handles,
                        (slong) (((ulong) nrows * k + 7) / 8));

            for (r = 0; r < nrows; r++)
                for (j = 0; j < ncols; j++)
                {
                    slong zl = 0;
                    for (l = 0; l < k; l++)
                    {
                        const slong la = nmod_poly_mat_entry(A, g + r, l)->length;
                        const slong lb = nmod_poly_mat_entry(B, l, h + j)->length;
                        if (la > 0 && lb > 0)
                            zl = FLINT_MAX(zl, la + lb - 1);
                    }
                    zlen[r * NCG + j] = zl;
                }

            TIMING_MARK(_tA);

            _matmul_run(_worker_matmul, W, nthreads, handles,
                        (slong) (np * ztrunc));

            TIMING_MARK(_tP);

            _matmul_run(_worker_ifft, W, nthreads, handles, nrows * ncols);

            TIMING_MARK(_tC);
        }
    }

    TIMING_PRINT;

    flint_free(W);
    flint_free(zlen);
    _matmul_free(buf, nbytes, bbytes);
    flint_give_back_threads(handles, nworkers);
    fft_small_plan_clear(P);
}

#else  /* PML_HAVE_MACHINE_VECTORS */

void nmod_poly_mat_mul_sd_fft_matmul(nmod_poly_mat_t C,
                                     const nmod_poly_mat_t A,
                                     const nmod_poly_mat_t B)
{
    nmod_poly_mat_mul(C, A, B);
}

#endif  /* PML_HAVE_MACHINE_VECTORS */
