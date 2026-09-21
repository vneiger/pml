/*
    Copyright (C) 2025 Vincent Neiger, Éric Schost

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/nmod.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>
#include <flint/thread_pool.h>
#include <flint/thread_support.h>
#include <flint/ulong_extras.h>

#include "nmod_extra.h" /* for nmod_find_root */

#include "nmod_poly_extra.h"
#include "nmod_poly_mat_forms.h" /* for row/column degrees */
#include "nmod_poly_mat_multiply.h"
#include "impl.h"

/* timing breakdown of the phases, printed on stderr (development aid) */
#ifdef PML_MUL_GEOMETRIC_TIMING
# include <stdio.h>
# include <stdlib.h>
# include <time.h>
static double _now(void) { struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return ts.tv_sec + 1e-9 * ts.tv_nsec; }
# define TIMING_DECL double _t0 = _now(), _t1; double _tB = 0, _tA = 0, _tP = 0, _tC = 0;
# define TIMING_MARK(acc) do { _t1 = _now(); acc += _t1 - _t0; _t0 = _t1; } while (0)
# define TIMING_PRINT flint_fprintf(stderr, "  [geometric] len=%wd NRG=%wd NCG=%wd | eval B %.3f, eval A %.3f, matmul %.3f, interp %.3f\n", len, NRG, NCG, _tB, _tA, _tP, _tC)
#else
# define TIMING_DECL
# define TIMING_MARK(acc)
# define TIMING_PRINT
#endif

/*
    Multiplication of polynomial matrices by evaluation and interpolation
    at a geometric progression: C = A * B is obtained by evaluating every
    entry of A and of B at len = lenA + lenB - 1 points 1, q, q^2, ...,
    multiplying the resulting constant matrices pointwise, and
    interpolating every entry of C.

    The organisation mirrors nmod_poly_mat_mul_sd_fft_direct and
    nmod_poly_mat_mul_sd_fft_matmul, for the same reasons:

      * the 3*len constant matrices are slices of one buffer, so that
        nmod_mat_mul reads and writes borrowed storage and the scattered
        side of the two transpositions has a constant, cache-line aligned
        stride (see _mats_struct);
      * that buffer is retained across calls (see _geometric_fit_buffer):
        the first-touch page faults of a freshly mapped working set are
        a sizeable fraction of a large product, and the polynomial matrix
        algorithms of PML perform sequences of products of similar shape;
      * the four phases -- evaluate B, evaluate a group of rows of A,
        the products at each point, interpolate C -- are loops over
        independent tasks run over the FLINT thread pool;
      * the memory is bounded, with the rows of A and the columns of B
        processed by groups when it does not fit;
      * an entry of C with no nonzero product contributing to it is not
        interpolated at all, and the others are interpolated from as many
        points as their true length rather than from all len of them.
*/

/* ------------------------------------------------------------------------ */
/* shared machinery of the product and of the middle product                */
/* ------------------------------------------------------------------------ */

/*
    What follows is shared by nmod_poly_mat_mul_geometric and
    nmod_poly_mat_mulmid_geometric (both in this file): the memory bound,
    the retained buffer, the layout of the constant matrices, the
    splitting of the operands into groups of rows and columns that fit the
    bound, and the dispatch of a phase over the thread pool.
*/

/* Soft bound, in bytes, on the memory used for the constant matrices.
   This is not a constant but a floor: the bound used is the larger of it
   and twice the size of the operands and the result. The shape of a group
   matters -- the group sizes are the outer dimensions of the matrices
   handed to nmod_mat_mul, and that routine changes algorithm on the
   smallest of the three -- so the two are balanced rather than the rows
   of A alone being cut (see _geometric_groups). */
#define NMOD_POLY_MAT_GEOMETRIC_MEM_FLOOR (UWORD(1) << 28)

/* a rounded up to the next multiple of the power of two b (n_round_up
   itself lives in fft_small.h, which this file must not depend on) */
static ulong _round_up2(ulong a, ulong b)
{
    return (a + b - 1) & ~(b - 1);
}

/*
    Retained scratch space for the constant matrices. The buffer is thread
    local and is not given back between products, only grown;
    flint_cleanup() releases it (and forgets the cleanup functions it has
    run, so the flag below is reset for the next product to register
    again). A thread that has performed one large product therefore keeps
    a working set of the order of the bound above.
*/
static FLINT_TLS_PREFIX ulong * _geometric_buf = NULL;
static FLINT_TLS_PREFIX ulong _geometric_buf_alloc = 0;   /* bytes */
static FLINT_TLS_PREFIX int _geometric_buf_registered = 0;

static void _geometric_cleanup(void)
{
    if (_geometric_buf != NULL)
        flint_aligned_free(_geometric_buf);
    _geometric_buf = NULL;
    _geometric_buf_alloc = 0;
    _geometric_buf_registered = 0;
}

static ulong * _geometric_fit_buffer(ulong nbytes)
{
    if (nbytes > _geometric_buf_alloc)
    {
        nbytes = _round_up2(nbytes, UWORD(4096));
        if (_geometric_buf != NULL)
            flint_aligned_free(_geometric_buf);
        _geometric_buf = (ulong *) flint_aligned_alloc(64, nbytes);
        _geometric_buf_alloc = nbytes;
        if (!_geometric_buf_registered)
        {
            flint_register_cleanup_function(_geometric_cleanup);
            _geometric_buf_registered = 1;
        }
    }
    return _geometric_buf;
}

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

static ulong * _mats_at(const _mats_struct * M, slong t)
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
static ulong _odd_lines(ulong words)
{
    words = _round_up2(words, 8);
    return ((words / 8) % 2 == 0) ? words + 8 : words;
}

/* distance in words between the matrices of consecutive points */
static slong _mats_stride(slong r, slong c)
{
    return (slong) _odd_lines((ulong) r * c);
}

/* words occupied by a region of r x c matrices at len points */
static ulong _mats_words(slong r, slong c, slong len)
{
    return _odd_lines((ulong) len * (ulong) _mats_stride(r, c));
}

/* the len values of entry e, scattered over the matrices */
static void _mats_scatter(const _mats_struct * M, slong e, nn_srcptr v)
{
    const slong len = M->len, ms = M->mstride;
    ulong * dst = M->data + e;
    slong t;

    for (t = 0; t < len; t++)
        dst[t * ms] = v[t];
}

/* the same for an entry which is the zero polynomial */
static void _mats_scatter_zero(const _mats_struct * M, slong e)
{
    const slong len = M->len, ms = M->mstride;
    ulong * dst = M->data + e;
    slong t;

    for (t = 0; t < len; t++)
        dst[t * ms] = 0;
}

/* the first zl values of entry e, gathered from the matrices */
static void _mats_gather(nn_ptr v, const _mats_struct * M, slong e, slong zl)
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
static void _geometric_groups(slong m, slong k, slong n, ulong avail,
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

static void _geometric_run(void (* func)(void *),
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

/* ------------------------------------------------------------------------ */
/* parallel phases                                                          */
/* ------------------------------------------------------------------------ */

/*
    Threads are requested once for the whole product, so the evaluation,
    the interpolation and nmod_mat_mul called inside a worker find none
    left and run serially: one task per thread is coarser, and scales
    further, than their own threading. Tasks of the two evaluation
    phases are groups of eight entries, so that the cache lines of the
    matrices are not shared between threads (see _mats_struct).
*/
typedef struct
{
    _geometric_range_struct range;   /* the tasks of this worker; first, see _geometric_run */

    /* shared, read only */
    const nmod_geometric_progression_struct * G;
    const _mats_struct * Am;
    const _mats_struct * Bm;
    const _mats_struct * Cm;
    const nmod_poly_mat_struct * A;
    const nmod_poly_mat_struct * B;
    nmod_poly_mat_struct * C;
    nmod_t mod;
    slong k;
    slong len;
    slong g;        /* first row of the current group of rows of A */
    slong nrows;    /* its size */
    slong h;        /* first column of the current group of columns of B */
    slong ncols;    /* its size */
    slong cstride;  /* NCG: column stride of the C and B blocks */
    const slong * zlen;

    /* per worker */
    nn_ptr val;
}
_geometric_worker_struct;

/* tasks are groups of 8 of the entries e = l*NCG + jj of the block of B,
   holding the evaluations of B[l][h+jj] */
static void _worker_eval_B(void * varg)
{
    _geometric_worker_struct * W = (_geometric_worker_struct *) varg;
    const slong cs = W->cstride, nc = W->ncols;
    const slong nops = W->Bm->r * W->Bm->c;
    slong e;

    for (e = 8 * W->range.start; e < 8 * W->range.stop && e < nops; e++)
    {
        const slong l = e / cs, jj = e % cs;
        const nmod_poly_struct * b;

        if (jj >= nc)   /* padding column of the block: never read */
            continue;

        b = nmod_poly_mat_entry(W->B, l, W->h + jj);
        if (b->length > 0)
        {
            _nmod_poly_evaluate_geometric_nmod_vec_fast_precomp(W->val, b->coeffs,
                                            b->length, W->G, W->len, W->mod);
            _mats_scatter(W->Bm, e, W->val);
        }
        else
            _mats_scatter_zero(W->Bm, e);
    }
}

/* tasks are groups of 8 of the entries e = r*k + l, holding the
   evaluations of A[g+r][l] */
static void _worker_eval_A(void * varg)
{
    _geometric_worker_struct * W = (_geometric_worker_struct *) varg;
    const slong k = W->k;
    const slong nops = W->nrows * k;
    slong e;

    for (e = 8 * W->range.start; e < 8 * W->range.stop && e < nops; e++)
    {
        const nmod_poly_struct * a =
            nmod_poly_mat_entry(W->A, W->g + e / k, e % k);

        if (a->length > 0)
        {
            _nmod_poly_evaluate_geometric_nmod_vec_fast_precomp(W->val, a->coeffs,
                                            a->length, W->G, W->len, W->mod);
            _mats_scatter(W->Am, e, W->val);
        }
        else
            _mats_scatter_zero(W->Am, e);
    }
}

/* tasks are the evaluation points; distinct tasks write distinct
   matrices of C, so no synchronization is needed */
static void _worker_matmul(void * varg)
{
    _geometric_worker_struct * W = (_geometric_worker_struct *) varg;
    slong t;

    for (t = W->range.start; t < W->range.stop; t++)
    {
        nmod_mat_struct MA, MB, MC;

        MA.entries = _mats_at(W->Am, t);
        MA.r = W->nrows; MA.c = W->k; MA.stride = W->k; MA.mod = W->mod;

        MB.entries = _mats_at(W->Bm, t);
        MB.r = W->k; MB.c = W->ncols; MB.stride = W->cstride; MB.mod = W->mod;

        MC.entries = _mats_at(W->Cm, t);
        MC.r = W->nrows; MC.c = W->ncols; MC.stride = W->cstride; MC.mod = W->mod;

        nmod_mat_mul(&MC, &MA, &MB);
    }
}

/* tasks e = r*ncols + jj: the entry C[g+r][h+jj] */
static void _worker_interp(void * varg)
{
    _geometric_worker_struct * W = (_geometric_worker_struct *) varg;
    const slong cs = W->cstride, nc = W->ncols;
    slong e;

    for (e = W->range.start; e < W->range.stop; e++)
    {
        /* the block is ncols wide but strided by cstride */
        const slong ee = (e / nc) * cs + e % nc;
        nmod_poly_struct * c =
            nmod_poly_mat_entry(W->C, W->g + e / nc, W->h + e % nc);
        const slong zl = W->zlen[ee];

        if (zl == 0)
        {
            nmod_poly_zero(c);
            continue;
        }

        _mats_gather(W->val, W->Cm, ee, zl);
        nmod_poly_interpolate_geometric_nmod_vec_fast_precomp(c, W->val, W->G, zl);
    }
}

/** Multiplication for polynomial matrices
 *  sets C = A * B
 *  output can alias input  (TODO make this consistent with existing functions)
 *  ASSUMPTION (not checked): existence of element of "large enough" order
 *  TODO -> fail flag when element not found
 *  uses evaluation and interpolation at a geometric progression
 */
void _nmod_poly_mat_mul_geometric_precomp(nmod_poly_mat_t res,
                                          const nmod_poly_mat_t pmat1, slong len1,
                                          const nmod_poly_mat_t pmat2, slong len2,
                                          nmod_geometric_progression_t G)
{
    _nmod_poly_mat_mul_geometric_precomp_bounded(res, pmat1, len1, pmat2, len2, G, 0);
}

void _nmod_poly_mat_mul_geometric_precomp_bounded(nmod_poly_mat_t res,
                                          const nmod_poly_mat_t pmat1, slong len1,
                                          const nmod_poly_mat_t pmat2, slong len2,
                                          nmod_geometric_progression_t G,
                                          ulong membytes)
{
    if (len1 == 0 || len2 == 0)
    {
        nmod_poly_mat_zero(res);
        return;
    }

    if (res == pmat1 || res == pmat2)
    {
        nmod_poly_mat_t tmp;
        nmod_poly_mat_init(tmp, pmat1->r, pmat2->c, pmat1->modulus);
        _nmod_poly_mat_mul_geometric_precomp_bounded(tmp, pmat1, len1, pmat2, len2, G, membytes);
        nmod_poly_mat_swap_entrywise(res, tmp);
        nmod_poly_mat_clear(tmp);
        return;
    }

    const slong m = pmat1->r;
    const slong k = pmat1->c;
    const slong n = pmat2->c;
    const slong len = len1 + len2 - 1;

    if (m == 0 || n == 0)
        return;

    if (k == 0)
    {
        nmod_poly_mat_zero(res);
        return;
    }

    TIMING_DECL

    /* threads: one task is a group of evaluations, or the product at one
       point; there is little to gain below a few entries */
    thread_pool_handle * handles = NULL;
    slong nworkers = 0;
    slong nthreads = 1;

    if ((ulong) m * k * n >= 4 && (ulong) m * k * n * len >= (UWORD(1) << 18))
    {
        nworkers = flint_request_threads(&handles, flint_get_num_threads());
        nthreads = nworkers + 1;
    }

    /* words of the operands and of the result, as the scale of the bound
       on the constant matrices (see the MEM_FLOOR constant) */
    const ulong opwords = (ulong) m * k * len1 + (ulong) k * n * len2
                          + (ulong) m * n * len;
    ulong bbytes = (membytes != 0)
                   ? membytes
                   : FLINT_MAX(NMOD_POLY_MAT_GEOMETRIC_MEM_FLOOR,
                               2 * opwords * sizeof(ulong));
#ifdef PML_MUL_GEOMETRIC_TIMING
    if (membytes == 0 && getenv("PML_MEM_BUDGET"))
        bbytes = strtoul(getenv("PML_MEM_BUDGET"), NULL, 10);
#endif

    /* NRG rows of A and C, NCG columns of B and C are held at a time */
    slong NRG, NCG;
    {
        /* budget and requirements counted in entries, i.e. len words */
        const ulong budget = bbytes / ((ulong) len * sizeof(ulong));
        const ulong scratch = (ulong) nthreads;   /* the val buffers */
        const ulong avail = (budget > scratch + 3) ? budget - scratch : 3;

        _geometric_groups(m, k, n, avail, &NRG, &NCG);
    }

    /* storage: the per-thread value buffers, then the evaluations of B
       (k x NCG), of the current rows of A (NRG x k) and of C (NRG x NCG) */
    const ulong wordsB = _mats_words(k, NCG, len);
    const ulong wordsA = _mats_words(NRG, k, len);
    const ulong wordsC = _mats_words(NRG, NCG, len);
    const ulong valwords = _round_up2((ulong) nthreads * len, 8);
    ulong * buf = _geometric_fit_buffer((valwords + wordsB + wordsA + wordsC)
                                        * sizeof(ulong));
    ulong * matbuf = buf + valwords;

    _mats_struct Bm = { matbuf, k, NCG, _mats_stride(k, NCG), len };
    _mats_struct Am = { matbuf + wordsB, NRG, k,
                        _mats_stride(NRG, k), len };
    _mats_struct Cm = { matbuf + wordsB + wordsA, NRG, NCG,
                        _mats_stride(NRG, NCG), len };

    /* zlen[r*NCG+jj] is the length of C[g+r][h+jj], zero when that entry
       has no nonzero product contributing to it */
    slong * zlen = FLINT_ARRAY_ALLOC((ulong) NRG * NCG, slong);

    /* The length of C[i][j] is bounded by rdegA[i] + cdegB[j] + 1. A zero row
       or column is reported as degree -1, which gives zlen = 0.
       Note that evaluating the more precise max_l (len(A[i][l]) +
       len(B[l][j]) - 1) costs O(m*k*n) and has been observed too costly for
       large dimensions and small lengths (for example, slowdown by a factor
       of 2 for 1024x1024x1024 in input length 10). */
    slong * rdegA = FLINT_ARRAY_ALLOC(m, slong);
    slong * cdegB = FLINT_ARRAY_ALLOC(n, slong);
    nmod_poly_mat_row_degree(rdegA, pmat1, NULL);
    nmod_poly_mat_column_degree(cdegB, pmat2, NULL);

    _geometric_worker_struct * W =
        FLINT_ARRAY_ALLOC(nthreads, _geometric_worker_struct);
    {
        slong w;
        for (w = 0; w < nthreads; w++)
        {
            W[w].G = G;
            W[w].Am = &Am; W[w].Bm = &Bm; W[w].Cm = &Cm;
            W[w].A = pmat1; W[w].B = pmat2; W[w].C = res;
            W[w].mod = G->mod;
            W[w].k = k;
            W[w].len = len;
            W[w].g = 0; W[w].nrows = 0;
            W[w].h = 0; W[w].ncols = 0; W[w].cstride = NCG;
            W[w].zlen = zlen;
            W[w].val = buf + (ulong) w * len;
        }
    }

    slong g, h;
    for (h = 0; h < n; h += NCG)
    {
        const slong ncols = FLINT_MIN(NCG, n - h);
        slong r, j, w;

        for (w = 0; w < nthreads; w++)
        {
            W[w].h = h;
            W[w].ncols = ncols;
        }

        _geometric_run(_worker_eval_B, W, sizeof(*W), nthreads, handles,
                       ((slong) k * NCG + 7) / 8);

        TIMING_MARK(_tB);

        for (g = 0; g < m; g += NRG)
        {
            const slong nrows = FLINT_MIN(NRG, m - g);

            for (w = 0; w < nthreads; w++)
            {
                W[w].g = g;
                W[w].nrows = nrows;
            }

            _geometric_run(_worker_eval_A, W, sizeof(*W), nthreads, handles,
                           (nrows * k + 7) / 8);

            for (r = 0; r < nrows; r++)
                for (j = 0; j < ncols; j++)
                    zlen[r * NCG + j] =
                        (rdegA[g + r] >= 0 && cdegB[h + j] >= 0)
                        ? rdegA[g + r] + cdegB[h + j] + 1 : 0;

            TIMING_MARK(_tA);

            _geometric_run(_worker_matmul, W, sizeof(*W), nthreads, handles, len);

            TIMING_MARK(_tP);

            _geometric_run(_worker_interp, W, sizeof(*W), nthreads, handles, nrows * ncols);

            TIMING_MARK(_tC);
        }
    }

    TIMING_PRINT;

    flint_free(W);
    flint_free(zlen);
    flint_free(rdegA);
    flint_free(cdegB);
    flint_give_back_threads(handles, nworkers);
}

void nmod_poly_mat_mul_geometric(nmod_poly_mat_t res, const nmod_poly_mat_t pmat1, const nmod_poly_mat_t pmat2)
{
    if (res == pmat1 || res == pmat2)
    {
        nmod_poly_mat_t T;
        nmod_poly_mat_init(T, pmat1->r, pmat2->c, pmat1->modulus);
        nmod_poly_mat_mul_geometric(T, pmat1, pmat2);
        nmod_poly_mat_swap_entrywise(res, T);
        nmod_poly_mat_clear(T);
        return;
    }

    const slong len1 = nmod_poly_mat_max_length(pmat1);
    const slong len2 = nmod_poly_mat_max_length(pmat2);

    if (len1 == 0 || len2 == 0)
    {
        nmod_poly_mat_zero(res);
        return;
    }

    const slong len = len1 + len2 - 1;

    nmod_t mod;
    nmod_init(&mod, pmat1->modulus);
    nmod_geometric_progression_t G;
    ulong w = nmod_find_root(2*len, mod);
#if (__FLINT_VERSION == 3 && __FLINT_VERSION_MINOR >= 6)
    _nmod_geometric_progression_init_function(G, w, len, mod, UWORD(3));
#else
    nmod_geometric_progression_init(G, w, len, mod);
#endif

    _nmod_poly_mat_mul_geometric_precomp(res, pmat1, len1, pmat2, len2, G);

    nmod_geometric_progression_clear(G);
}

/* ------------------------------------------------------------------------ */
/* the middle product                                                       */
/* ------------------------------------------------------------------------ */

/*
    Middle product by the transposition principle, at a geometric
    progression 1, q, q^2, ... of nhi points. The coefficients nlo, ...,
    nhi-1 of A*B are wanted, and one of the two operands -- call it X,
    the other Y -- has entries of length at most nlo+1. Then

      * every entry of X is reversed as a polynomial of length nlo+1 and
        evaluated at the nhi points;
      * every entry of Y, padded to nhi coefficients, is read as a vector
        of values and interpolated (the transpose of evaluation);
      * the constant matrices are multiplied pointwise, in the order of
        A and B;
      * the nhi values of every entry of the product are read as the
        coefficients of a polynomial and evaluated at the first nhi-nlo
        points: those are the wanted coefficients, in order.

    The organisation is that of nmod_poly_mat_mul_geometric above, whose
    machinery this shares: the 3*nhi constant matrices are slices of the
    retained buffer, the four phases are loops over independent tasks run
    over the FLINT thread pool, the memory is bounded with the rows of A
    and the columns of B processed by groups when it does not fit, and an
    entry of the result is evaluated at only as many points as it has
    nonzero coefficients -- none at all when no product contributes to the
    wanted range.
*/

typedef struct
{
    _geometric_range_struct range;   /* the tasks of this worker; first, see _geometric_run */

    /* shared, read only */
    const nmod_geometric_progression_struct * G;
    const _mats_struct * Am;
    const _mats_struct * Bm;
    const _mats_struct * Cm;
    const nmod_poly_mat_struct * A;
    const nmod_poly_mat_struct * B;
    nmod_poly_mat_struct * C;
    nmod_t mod;
    slong k;
    slong nlo;
    slong nhi;
    int revA;       /* A is the reversed operand X (else B is) */
    slong g;        /* first row of the current group of rows of A */
    slong nrows;    /* its size */
    slong h;        /* first column of the current group of columns of B */
    slong ncols;    /* its size */
    slong cstride;  /* NCG: column stride of the C and B blocks */
    const slong * zlen;

    /* per worker */
    nn_ptr val;     /* nhi words */
    nn_ptr poly;    /* nhi words */
}
_mulmid_worker_struct;

/* the nhi values of entry p, transformed as the operand X or Y */
static void _mulmid_transform(const _mulmid_worker_struct * W,
                              const nmod_poly_struct * p, int reversed)
{
    if (reversed)
    {
        /* X: reverse as a polynomial of length nlo+1, evaluate */
        _nmod_poly_reverse(W->poly, p->coeffs, p->length, W->nlo + 1);
        _nmod_poly_evaluate_geometric_nmod_vec_fast_precomp(W->val, W->poly,
                                    W->nlo + 1, W->G, W->nhi, W->mod);
    }
    else
    {
        /* Y: coefficients as values, interpolate */
        const slong len = FLINT_MIN(W->nhi, p->length);
        _nmod_vec_set(W->poly, p->coeffs, len);
        _nmod_vec_zero(W->poly + len, W->nhi - len);
        _nmod_poly_interpolate_geometric_nmod_vec_fast_precomp(W->val, W->poly,
                                    W->G, W->nhi, W->mod);
    }
}

/* tasks are groups of 8 of the entries e = l*NCG + jj of the block of B,
   holding the transform of B[l][h+jj] */
static void _mulmid_worker_B(void * varg)
{
    _mulmid_worker_struct * W = (_mulmid_worker_struct *) varg;
    const slong cs = W->cstride, nc = W->ncols;
    const slong nops = W->Bm->r * W->Bm->c;
    slong e;

    for (e = 8 * W->range.start; e < 8 * W->range.stop && e < nops; e++)
    {
        const slong l = e / cs, jj = e % cs;
        const nmod_poly_struct * b;

        if (jj >= nc)   /* padding column of the block: never read */
            continue;

        b = nmod_poly_mat_entry(W->B, l, W->h + jj);
        if (b->length > 0)
        {
            _mulmid_transform(W, b, !W->revA);
            _mats_scatter(W->Bm, e, W->val);
        }
        else
            _mats_scatter_zero(W->Bm, e);
    }
}

/* tasks are groups of 8 of the entries e = r*k + l, holding the
   transform of A[g+r][l] */
static void _mulmid_worker_A(void * varg)
{
    _mulmid_worker_struct * W = (_mulmid_worker_struct *) varg;
    const slong k = W->k;
    const slong nops = W->nrows * k;
    slong e;

    for (e = 8 * W->range.start; e < 8 * W->range.stop && e < nops; e++)
    {
        const nmod_poly_struct * a =
            nmod_poly_mat_entry(W->A, W->g + e / k, e % k);

        if (a->length > 0)
        {
            _mulmid_transform(W, a, W->revA);
            _mats_scatter(W->Am, e, W->val);
        }
        else
            _mats_scatter_zero(W->Am, e);
    }
}

/* tasks are the points; distinct tasks write distinct matrices of C */
static void _mulmid_worker_matmul(void * varg)
{
    _mulmid_worker_struct * W = (_mulmid_worker_struct *) varg;
    slong t;

    for (t = W->range.start; t < W->range.stop; t++)
    {
        nmod_mat_struct MA, MB, MC;

        MA.entries = _mats_at(W->Am, t);
        MA.r = W->nrows; MA.c = W->k; MA.stride = W->k; MA.mod = W->mod;

        MB.entries = _mats_at(W->Bm, t);
        MB.r = W->k; MB.c = W->ncols; MB.stride = W->cstride; MB.mod = W->mod;

        MC.entries = _mats_at(W->Cm, t);
        MC.r = W->nrows; MC.c = W->ncols; MC.stride = W->cstride; MC.mod = W->mod;

        nmod_mat_mul(&MC, &MA, &MB);
    }
}

/* tasks e = r*ncols + jj: the entry C[g+r][h+jj], evaluated at as many
   points as it has coefficients in the wanted range */
static void _mulmid_worker_C(void * varg)
{
    _mulmid_worker_struct * W = (_mulmid_worker_struct *) varg;
    const slong cs = W->cstride, nc = W->ncols;
    slong e;

    for (e = W->range.start; e < W->range.stop; e++)
    {
        /* the block is ncols wide but strided by cstride */
        const slong ee = (e / nc) * cs + e % nc;
        nmod_poly_struct * c =
            nmod_poly_mat_entry(W->C, W->g + e / nc, W->h + e % nc);
        const slong zl = W->zlen[ee];

        if (zl == 0)
        {
            nmod_poly_zero(c);
            continue;
        }

        _mats_gather(W->poly, W->Cm, ee, W->nhi);
        _nmod_poly_evaluate_geometric_nmod_vec_fast_precomp(W->val, W->poly,
                                    W->nhi, W->G, zl, W->mod);
        nmod_poly_fit_length(c, zl);
        _nmod_vec_set(c->coeffs, W->val, zl);
        _nmod_poly_set_length(c, zl);
        _nmod_poly_normalise(c);
    }
}

void _nmod_poly_mat_mulmid_geometric_precomp_bounded(nmod_poly_mat_t res,
                                             const nmod_poly_mat_t pmat1, slong len1,
                                             const nmod_poly_mat_t pmat2, slong len2,
                                             slong nlo, slong nhi,
                                             nmod_geometric_progression_t G,
                                             ulong membytes)
{
    PML_ASSERT(len1 <= nlo+1 || len2 <= nlo+1)

    if (len1 == 0 || len2 == 0 || nlo >= nhi)
    {
        nmod_poly_mat_zero(res);
        return;
    }

    if (res == pmat1 || res == pmat2)
    {
        nmod_poly_mat_t temp;
        nmod_poly_mat_init(temp, pmat1->r, pmat2->c, pmat1->modulus);
        _nmod_poly_mat_mulmid_geometric_precomp_bounded(temp, pmat1, len1, pmat2, len2,
                                                        nlo, nhi, G, membytes);
        nmod_poly_mat_swap(res, temp);
        nmod_poly_mat_clear(temp);
        return;
    }

    const slong m = pmat1->r;
    const slong k = pmat1->c;
    const slong n = pmat2->c;
    const int revA = (len1 <= nlo + 1);

    if (m == 0 || n == 0)
        return;

    if (k == 0)
    {
        nmod_poly_mat_zero(res);
        return;
    }

    /* threads: one task is a group of transforms, or the product at one
       point; there is little to gain below a few entries */
    thread_pool_handle * handles = NULL;
    slong nworkers = 0;
    slong nthreads = 1;

    if ((ulong) m * k * n >= 4 && (ulong) m * k * n * nhi >= (UWORD(1) << 18))
    {
        nworkers = flint_request_threads(&handles, flint_get_num_threads());
        nthreads = nworkers + 1;
    }

    /* words of the operands and of the result, as the scale of the bound
       on the constant matrices (see NMOD_POLY_MAT_GEOMETRIC_MEM_FLOOR) */
    const ulong opwords = (ulong) m * k * len1 + (ulong) k * n * len2
                          + (ulong) m * n * (nhi - nlo);
    const ulong bbytes = (membytes != 0)
                         ? membytes
                         : FLINT_MAX(NMOD_POLY_MAT_GEOMETRIC_MEM_FLOOR,
                                     2 * opwords * sizeof(ulong));

    /* NRG rows of A and C, NCG columns of B and C are held at a time */
    slong NRG, NCG;
    {
        /* budget and requirements counted in entries, i.e. nhi words */
        const ulong budget = bbytes / ((ulong) nhi * sizeof(ulong));
        const ulong scratch = 2 * (ulong) nthreads;   /* val and poly */
        const ulong avail = (budget > scratch + 3) ? budget - scratch : 3;

        _geometric_groups(m, k, n, avail, &NRG, &NCG);
    }

    /* storage: the per-thread scratch, then the transforms of B (k x NCG),
       of the current rows of A (NRG x k) and of C (NRG x NCG) */
    const ulong wordsB = _mats_words(k, NCG, nhi);
    const ulong wordsA = _mats_words(NRG, k, nhi);
    const ulong wordsC = _mats_words(NRG, NCG, nhi);
    const ulong scratchwords = _round_up2(2 * (ulong) nthreads * nhi, 8);
    ulong * buf = _geometric_fit_buffer(
                      (scratchwords + wordsB + wordsA + wordsC) * sizeof(ulong));
    ulong * matbuf = buf + scratchwords;

    _mats_struct Bm = { matbuf, k, NCG, _mats_stride(k, NCG), nhi };
    _mats_struct Am = { matbuf + wordsB, NRG, k,
                        _mats_stride(NRG, k), nhi };
    _mats_struct Cm = { matbuf + wordsB + wordsA, NRG, NCG,
                        _mats_stride(NRG, NCG), nhi };

    /* zlen[r*NCG+jj] is the number of coefficients of C[g+r][h+jj] in the
       range nlo..nhi-1, zero when no product reaches that range */
    slong * zlen = FLINT_ARRAY_ALLOC((ulong) NRG * NCG, slong);

    /* The length of C[i][j] is bounded by rdegA[i] + cdegB[j] + 1. A zero row
       or column is reported as degree -1, which gives zlen = 0.
       Note that evaluating the more precise max_l (len(A[i][l]) +
       len(B[l][j]) - 1) costs O(m*k*n) and has been observed too costly for
       large dimensions and small lengths (for example, slowdown by a factor
       of 2 for 1024x1024x1024 in input length 10). */
    slong * rdegA = FLINT_ARRAY_ALLOC(m, slong);
    slong * cdegB = FLINT_ARRAY_ALLOC(n, slong);
    nmod_poly_mat_row_degree(rdegA, pmat1, NULL);
    nmod_poly_mat_column_degree(cdegB, pmat2, NULL);

    _mulmid_worker_struct * W = FLINT_ARRAY_ALLOC(nthreads, _mulmid_worker_struct);
    {
        slong w;
        for (w = 0; w < nthreads; w++)
        {
            W[w].G = G;
            W[w].Am = &Am; W[w].Bm = &Bm; W[w].Cm = &Cm;
            W[w].A = pmat1; W[w].B = pmat2; W[w].C = res;
            W[w].mod = G->mod;
            W[w].k = k;
            W[w].nlo = nlo; W[w].nhi = nhi;
            W[w].revA = revA;
            W[w].g = 0; W[w].nrows = 0;
            W[w].h = 0; W[w].ncols = 0; W[w].cstride = NCG;
            W[w].zlen = zlen;
            W[w].val = buf + (2 * (ulong) w) * nhi;
            W[w].poly = buf + (2 * (ulong) w + 1) * nhi;
        }
    }

    slong g, h;
    for (h = 0; h < n; h += NCG)
    {
        const slong ncols = FLINT_MIN(NCG, n - h);
        slong r, j, w;

        for (w = 0; w < nthreads; w++)
        {
            W[w].h = h;
            W[w].ncols = ncols;
        }

        _geometric_run(_mulmid_worker_B, W, sizeof(*W), nthreads, handles,
                       ((slong) k * NCG + 7) / 8);

        for (g = 0; g < m; g += NRG)
        {
            const slong nrows = FLINT_MIN(NRG, m - g);

            for (w = 0; w < nthreads; w++)
            {
                W[w].g = g;
                W[w].nrows = nrows;
            }

            _geometric_run(_mulmid_worker_A, W, sizeof(*W), nthreads, handles,
                           (nrows * k + 7) / 8);

            /* the product A[g+r][l] * B[l][h+j] has coefficients up to
               index la+lb-2, of which those from nlo on are wanted;
               min(nhi, .) - nlo is nondecreasing in the product length,
               so the bound above applies here too */
            for (r = 0; r < nrows; r++)
                for (j = 0; j < ncols; j++)
                {
                    const slong zl =
                        (rdegA[g + r] >= 0 && cdegB[h + j] >= 0)
                        ? FLINT_MIN(nhi, rdegA[g + r] + cdegB[h + j] + 1) - nlo
                        : 0;
                    zlen[r * NCG + j] = FLINT_MAX(zl, 0);
                }

            _geometric_run(_mulmid_worker_matmul, W, sizeof(*W), nthreads, handles, nhi);

            _geometric_run(_mulmid_worker_C, W, sizeof(*W), nthreads, handles,
                           nrows * ncols);
        }
    }

    flint_free(W);
    flint_free(zlen);
    flint_free(rdegA);
    flint_free(cdegB);
    flint_give_back_threads(handles, nworkers);
}

void _nmod_poly_mat_mulmid_geometric_precomp(nmod_poly_mat_t res,
                                             const nmod_poly_mat_t pmat1, slong len1,
                                             const nmod_poly_mat_t pmat2, slong len2,
                                             slong nlo, slong nhi, nmod_geometric_progression_t G)
{
    _nmod_poly_mat_mulmid_geometric_precomp_bounded(res, pmat1, len1, pmat2, len2,
                                                    nlo, nhi, G, 0);
}

void _nmod_poly_mat_mulmid_geometric(nmod_poly_mat_t res,
                                     const nmod_poly_mat_t pmat1, slong len1,
                                     const nmod_poly_mat_t pmat2, slong len2,
                                     slong nlo, slong nhi)
{
    PML_ASSERT(len1 <= nlo+1 || len2 <= nlo+1)

    if (len1 == 0 || len2 == 0 || nlo >= nhi)
    {
        nmod_poly_mat_zero(res);
        return;
    }

    nmod_t mod;
    nmod_geometric_progression_t G;

    nmod_init(&mod, pmat1->modulus);
    ulong w = nmod_find_root(2*nhi, mod);
    _nmod_geometric_progression_init_function(G, w, nhi, mod, UWORD(3));

    _nmod_poly_mat_mulmid_geometric_precomp(res, pmat1, len1, pmat2, len2, nlo, nhi, G);

    nmod_geometric_progression_clear(G);
}
