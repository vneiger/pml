/*
    Copyright (C) 2025 Vincent Neiger, Éric Schost

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/mpn_extras.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>
#include <flint/thread_pool.h>
#include <flint/thread_support.h>

#include "nmod_poly_extra.h"  /* for NMOD_POLY_CAN_USE_GEOMETRIC */
#include "nmod_extra.h"  /* for nmod_find_root */
#include "nmod_poly_mat_multiply.h"
#include "impl.h"

/*------------------------------------------------------------*/
/* naive: multiply, truncate, shift                           */
/*------------------------------------------------------------*/

void _nmod_poly_mat_mulmid_naive(nmod_poly_mat_t res,
                                 const nmod_poly_mat_t pmat1, slong len1,
                                 const nmod_poly_mat_t pmat2, slong len2,
                                 slong nlo, slong nhi)
{
    if (len1 == 0 || len2 == 0 || nlo >= nhi)
    {
        nmod_poly_mat_zero(res);
        return;
    }

    nmod_poly_mat_multiply(res, pmat1, pmat2);
    nmod_poly_mat_shift_right(res, res, nlo);
    nmod_poly_mat_truncate(res, nhi - nlo);
}

/*------------------------------------------------------------*/
/* evaluation-interpolation at geometric progression          */
/*------------------------------------------------------------*/

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

    The organisation is that of nmod_poly_mat_mul_geometric, whose
    machinery this shares (see impl.h): the 3*nhi constant matrices are
    slices of one retained buffer, the four phases are loops over
    independent tasks run over the FLINT thread pool, the memory is
    bounded with the rows of A and the columns of B processed by groups
    when it does not fit, and an entry of the result is evaluated at only
    as many points as it has nonzero coefficients -- none at all when no
    product contributes to the wanted range.
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
    ulong * buf = _nmod_poly_mat_geometric_fit_buffer(
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
        slong r, j, l, w;

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
               index la+lb-2, of which those from nlo on are wanted */
            for (r = 0; r < nrows; r++)
                for (j = 0; j < ncols; j++)
                {
                    slong zl = 0;
                    for (l = 0; l < k; l++)
                    {
                        const slong la = nmod_poly_mat_entry(pmat1, g + r, l)->length;
                        const slong lb = nmod_poly_mat_entry(pmat2, l, h + j)->length;
                        if (la > 0 && lb > 0)
                            zl = FLINT_MAX(zl, FLINT_MIN(nhi, la + lb - 1) - nlo);
                    }
                    zlen[r * NCG + j] = FLINT_MAX(zl, 0);
                }

            _geometric_run(_mulmid_worker_matmul, W, sizeof(*W), nthreads, handles, nhi);

            _geometric_run(_mulmid_worker_C, W, sizeof(*W), nthreads, handles,
                           nrows * ncols);
        }
    }

    flint_free(W);
    flint_free(zlen);
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

/*------------------------------------------------------------*/
/* evaluation-interpolation using Vandermonde matrix          */
/*------------------------------------------------------------*/

/* TODO currently disabled: untested and unprofiled */
#if 0
/** Middle product for polynomial matrices
 *  sets C = ((A * B) div x^dA) mod x^(dB+1)
 *  output can alias input
 *  uses evaluation and interpolation at arithmetic points, done by matrix products
 *  ASSUME: large enough field 
 */
static void nmod_poly_mat_middle_product_vandermonde1(nmod_poly_mat_t c, const nmod_poly_mat_t a, const nmod_poly_mat_t b,
                                               const ulong dA, const ulong dB)
{

    ulong i, j, k, ellA, ellB, m, n, p, ell, u, v, nb_points;
    nmod_mat_t vA, vB, vB_t, iv, iv_t, tmp_mat, valA, valB, valAp, valBp, valC, valCp;
    nmod_t mod;
    
    ellA = nmod_poly_mat_max_length(a);
    ellB = nmod_poly_mat_max_length(b);

    if (ellA == 0 || ellB == 0)
    {
        nmod_poly_mat_zero(c);
        return;
    }

    m = a->r;
    n = a->c;
    p = b->c;

    if (c == a || c == b)
    {
        nmod_poly_mat_t T;
        nmod_poly_mat_init(T, m, n, p);
        nmod_poly_mat_middle_product_vandermonde1(T, a, b, dA, dB);
        nmod_poly_mat_swap_entrywise(c, T);
        nmod_poly_mat_clear(T);
        return;
    }
    
    nmod_init(&mod, a->modulus);
    vandermonde_init1(vA, vB, iv, dA, dB, mod);
    nb_points = dA + dB + 1;

    // evaluation of matrix a:
    // build tmp_mat, whose column ell = i*n + j is the coefficient vector 
    // of reverse(a[i][j]) 
    nmod_mat_init(tmp_mat, dA + 1, m * n, mod.n);
    ell = 0;
    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++, ell++)
        {
            ulong d;
            nn_ptr ptr;
            d = nmod_poly_degree(nmod_poly_mat_entry(a, i, j));
            ptr = nmod_poly_mat_entry(a, i, j)->coeffs;
            for (k = 0; k <= d; k++)
                nmod_mat_entry(tmp_mat, dA - k, ell) = ptr[k];
        }
    // note: d = deg(a[i][j]) is -1 if a[i][j] == 0
    // all non-touched entries already zero since tmp_mat was initialized as zero
   
    // valA: column ell = i*n + j contains the evaluations of a[i][j]
    nmod_mat_init(valA, vA->r, tmp_mat->c, mod.n);
    nmod_mat_mul_pml(valA, vA, tmp_mat);

    // transpose interpolation of matrix b:
    // build tmp_mat, whose column ell = i*n + j is the coefficient vector of
    // b[i][j] (padded with zeroes up to length dB+1 if necessary)
    nmod_mat_clear(tmp_mat);
    nmod_mat_init(tmp_mat, dA + dB + 1, n * p, mod.n);
    ell = 0;
    for (i = 0; i < n; i++)
        for ( j = 0; j < p; j++, ell++)
        {
            ulong d;
            nn_ptr ptr;
            d = nmod_poly_degree(nmod_poly_mat_entry(b, i, j));
            ptr = nmod_poly_mat_entry(b, i, j)->coeffs;
            for (k = 0; k <= d; k++)
                nmod_mat_entry(tmp_mat, k, ell) = ptr[k];
        }
    // mul by transpose(iv)
    nmod_mat_init(iv_t, iv->r, iv->c, mod.n);
    nmod_mat_transpose(iv_t, iv);
    nmod_mat_init(valB, iv_t->r, tmp_mat->c, mod.n);
    nmod_mat_mul_pml(valB, iv_t, tmp_mat);
 

    // perform the pointwise products
    nmod_mat_init(valAp, m, n, mod.n);
    nmod_mat_init(valBp, n, p, mod.n);
    nmod_mat_init(valC, nb_points, m * p, mod.n);
    nmod_mat_init(valCp, m, p, mod.n);
    
    for (i = 0; i < nb_points; i++)
    {
        // a evaluated at point i
        ell = 0;
        for (u = 0; u < m; u++)
            for (v = 0; v < n; v++, ell++)
                nmod_mat_entry(valAp, u, v) = nmod_mat_entry(valA, i, ell);

        ell = 0;
        for (u = 0; u < n; u++)
            for (v = 0; v < p; v++, ell++)
                nmod_mat_entry(valBp, u, v) = nmod_mat_entry(valB, i, ell);

        nmod_mat_mul_pml(valCp, valAp, valBp);

        // copy this into valC
        ell = 0;
        for (u = 0; u < m; u++)
            for (v = 0; v < p; v++, ell++)
                nmod_mat_entry(valC, i, ell) = nmod_mat_entry(valCp, u, v);
    }
    
    nmod_mat_init(vB_t, vB->c, vB->r, mod.n);
    nmod_mat_transpose(vB_t, vB);
    // transpose-evaluate to find the entries of c
    nmod_mat_clear(tmp_mat);
    nmod_mat_init(tmp_mat, vB_t->r, valC->c, mod.n);
    nmod_mat_mul_pml(tmp_mat, vB_t, valC);


    // copy to output (reorganize these entries into c)
    ell = 0;
    for (u = 0; u < m; u++)
        for (v = 0; v < p; v++, ell++)
        {
            nn_ptr coeffs;
            nmod_poly_realloc(nmod_poly_mat_entry(c, u, v), dB + 1);
            coeffs = nmod_poly_mat_entry(c, u, v)->coeffs;
            nmod_poly_mat_entry(c, u, v)->length = dB + 1; 
            for (i = 0; i <= dB; i++)
                coeffs[i] = nmod_mat_entry(tmp_mat, i, ell);
            _nmod_poly_normalise(nmod_poly_mat_entry(c, u, v));
        }
    
    nmod_mat_clear(vA);
    nmod_mat_clear(vB_t);
    nmod_mat_clear(vB);
    nmod_mat_clear(iv_t);
    nmod_mat_clear(iv);
    nmod_mat_clear(tmp_mat);
    nmod_mat_clear(valA);
    nmod_mat_clear(valB);
    nmod_mat_clear(valAp);
    nmod_mat_clear(valBp);
    nmod_mat_clear(valC);
    nmod_mat_clear(valCp);
}
#endif


/*------------------------------------------------------------*/
/* via 3-prime FFT                                            */
/*------------------------------------------------------------*/
/* TODO */




/*------------------------------------------------------------*/
/* general interfaces                                         */
/*------------------------------------------------------------*/

void _nmod_poly_mat_mulmid(nmod_poly_mat_t res,
                           const nmod_poly_mat_t pmat1, slong len1,
                           const nmod_poly_mat_t pmat2, slong len2,
                           slong nlo, slong nhi)
{
    /* zero matrices or empty target coefficient indices */
    if (len1 == 0 || len2 == 0 || nlo >= nhi)
    {
        nmod_poly_mat_zero(res);
        return;
    }

    /* TODO handle constant pmat1 or pmat2 properly */
    /* if (len1 == 1)*/
    /*     _nmod_poly_mat_mulmid_naive(res, pmat1, len1, pmat2, len2, nlo, nhi); */
    /* else */

    /* TODO rough thresholds, not finely tuned */
#if (__FLINT_VERSION == 3 && __FLINT_VERSION_MINOR >= 6)
    if (NMOD_POLY_CAN_USE_GEOMETRIC(pmat1->modulus, nhi)
        && ((pmat1->r >= 8 && pmat2->c >= 2) || (pmat1->r >= 2 && pmat2->c >= 8))
        && (len1 <= nlo+1 || len2 <= nlo+1))
        _nmod_poly_mat_mulmid_geometric(res, pmat1, len1, pmat2, len2, nlo, nhi);

    else
#endif
        _nmod_poly_mat_mulmid_naive(res, pmat1, len1, pmat2, len2, nlo, nhi);
}

void nmod_poly_mat_mulmid(nmod_poly_mat_t res, const nmod_poly_mat_t pmat1, const nmod_poly_mat_t pmat2,
                          slong nlo, slong nhi)
{
    PML_ASSERT(nlo >= 0);
    PML_ASSERT(nhi >= 0);

    slong len1 = nmod_poly_mat_max_length(pmat1);
    slong len2 = nmod_poly_mat_max_length(pmat2);
    nhi = FLINT_MIN(nhi, len1 + len2 - 1);

    /* zero matrices or empty target coefficient indices */
    if (len1 == 0 || len2 == 0 || nlo >= nhi)
    {
        nmod_poly_mat_zero(res);
        return;
    }

    if (len1 > nhi)
    {
        nmod_poly_mat_t pmat1_tmp;
        nmod_poly_mat_init(pmat1_tmp, pmat1->r, pmat1->c, pmat1->modulus);
        nmod_poly_mat_set_trunc(pmat1_tmp, pmat1, nhi);
        nmod_poly_mat_mulmid(res, pmat1_tmp, pmat2, nlo, nhi);
        nmod_poly_mat_clear(pmat1_tmp);
        return;
    }

    if (len2 > nhi)
    {
        nmod_poly_mat_t pmat2_tmp;
        nmod_poly_mat_init(pmat2_tmp, pmat2->r, pmat2->c, pmat2->modulus);
        nmod_poly_mat_set_trunc(pmat2_tmp, pmat2, nhi);
        nmod_poly_mat_mulmid(res, pmat1, pmat2_tmp, nlo, nhi);
        nmod_poly_mat_clear(pmat2_tmp);
        return;
    }

    /* len1 <= len2 and len(pmat1) <= nlo: shift away useless coefficients of pmat2 */
    if (len1 <= len2 && len1 <= nlo)
    {
        nmod_poly_mat_t pmat2_tmp;
        nmod_poly_mat_init(pmat2_tmp, pmat2->r, pmat2->c, pmat2->modulus);
        nmod_poly_mat_shift_right(pmat2_tmp, pmat2, nlo - len1 + 1);
        _nmod_poly_mat_mulmid(res, pmat1, len1, pmat2_tmp, len2 - nlo + len1 - 1, len1 - 1, nhi - nlo + len1 - 1);
        nmod_poly_mat_clear(pmat2_tmp);
        return;
    }

    /* if len(pmat2) <= nlo (implies len1 > len2), shift away useless coefficients of pmat1 */
    if (len2 <= nlo)
    {
        nmod_poly_mat_t pmat1_tmp;
        nmod_poly_mat_init(pmat1_tmp, pmat1->r, pmat1->c, pmat1->modulus);
        nmod_poly_mat_shift_right(pmat1_tmp, pmat1, nlo - len2 + 1);
        _nmod_poly_mat_mulmid(res, pmat1_tmp, len1 - nlo + len2 - 1, pmat2, len2, len2 - 1, nhi - nlo + len2 - 1);
        nmod_poly_mat_clear(pmat1_tmp);
        return;
    }

    _nmod_poly_mat_mulmid(res, pmat1, len1, pmat2, len2, nlo, nhi);
}
