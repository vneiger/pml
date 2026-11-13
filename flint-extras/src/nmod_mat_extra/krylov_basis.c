/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include <flint/nmod.h>
#include <flint/nmod_mat.h>
#include "nmod_mat_extra/krylov_impl.h"

/* default maximal number of rows processed at once: n / 8 (n = number of
   columns), clamped to [16, 128] (profile p-krylov chunk; on Zen 4, Ice
   Lake, Cascade Lake and Apple M4, Nov 2026, within 1.05 times the best
   size except at 60 bits for sigma = 128: 1.08 to 1.17) */
#define KRYLOV_CHUNK_DEFAULT(n) FLINT_MAX(16, FLINT_MIN(128, (n) / 8))

/*
    Automatic choice of the powers of M used to generate the rows (dense
    M, see krylov_impl.h): for n >= KRYLOV_POW_MIN_DIM, where the products
    of a few rows by M are limited by the memory bandwidth, and only in the
    chunks where KRYLOV_POW_MAX_RUNS = 1 index generates rows (m = 1, or a
    Hermite-type shift; with more indices, the rows are already computed by
    products with several rows). The squarings and the products of a few
    rows are done by nmod_mat_mul, so the gain depends on its speed for the
    modulus: M^2 and M^4 (logpow 2) for moduli below 2^32, and below 2^52
    with AVX-512 IFMA (FLINT's IFMA products).

    Measured at n = 1600, dense M, relative to products by M only (profile
    p-krylov pow; Nov 2026, Zen 4, Ice Lake, Cascade Lake, Apple M4): with
    one index, 0.66 to 0.95 below 2^32 (0.88 to 1.08 on Cascade Lake), where
    logpow 2 is within 1.17 times the best logpow and the former logpow 5 up
    to 1.37; 0.68 to 0.80 at 50 bits with IFMA, against 1.04 to 1.22
    without (logpow 1); with two indices, up to 1.88 times slower on x86
    (logpow 5). Above 2^52 there is no gain (1.25 to 1.72 with logpow 1)
    since FLINT's SIMD vector-matrix products for these moduli (FLINT
    #2907).
*/
#define KRYLOV_POW_MIN_DIM 1200
#define KRYLOV_POW_MAX_RUNS 1

static slong
_krylov_auto_logpow(slong n, ulong p)
{
    flint_bitcnt_t bits = FLINT_BIT_COUNT(p);

    if (n < KRYLOV_POW_MIN_DIM)
        return 0;
#if defined(__AVX512IFMA__)
    if (bits <= 52)
        return 2;
#else
    if (bits <= 32)
        return 2;
#endif
    return 0;
}

/* row states */
#define ACTIVE 0
#define DEAD 1       /* a dependent row (i, delta_i) was found */
#define EXHAUSTED 2  /* the row (i, d_i) is independent */

/* gather the columns cols[0..c) of A (rows 0..r) into the r x c matrix G */
static void
_gather_cols(nmod_mat_t G, const nmod_mat_t A, const slong * cols, slong r, slong c)
{
    slong u, t;

    for (u = 0; u < r; u++)
    {
        nn_srcptr a = nmod_mat_entry_ptr(A, u, 0);
        nn_ptr g = nmod_mat_entry_ptr(G, u, 0);
        for (t = 0; t < c; t++)
            g[t] = a[cols[t]];
    }
}

/* M and its highest power computed so far, P = M^(2^len) (len > 0); only
   these are kept, to limit the memory to one more n x n matrix */
typedef struct
{
    const nmod_mat_struct * M;
    nmod_mat_t P;
    slong len;
}
_krylov_pow_cache_struct;

/* raises P to M^(2^t), for t >= len */
static void
_krylov_pow_cache_raise(_krylov_pow_cache_struct * C, slong t)
{
    while (C->len < t)
    {
        nmod_mat_t T;

        nmod_mat_init(T, C->M->r, C->M->c, C->M->mod.n);
        if (C->len == 0)
        {
            nmod_mat_mul(T, C->M, C->M);
            nmod_mat_init(C->P, C->M->r, C->M->c, C->M->mod.n);
        }
        else
        {
            nmod_mat_mul(T, C->P, C->P);
        }
        nmod_mat_swap(C->P, T);
        nmod_mat_clear(T);
        C->len++;
    }
}

/* M^(2^t), for t = 0 or t = len */
static const nmod_mat_struct *
_krylov_pow_cache_get(_krylov_pow_cache_struct * C, slong t)
{
    FLINT_ASSERT(t == 0 || t == C->len);
    return (t == 0) ? C->M : C->P;
}

static void
_krylov_pow_cache_clear(_krylov_pow_cache_struct * C)
{
    if (C->len > 0)
        nmod_mat_clear(C->P);
}

/*
    Blocked iterative computation of the row rank profile, see
    krylov_impl.h for the specification.

    The rows (i, j) are processed in increasing priority (s_i + j, i), by
    chunks of consecutive rows. As soon as a row (i, j) is found to depend
    on the previous ones, all the rows (i, j') with j' > j are known to
    depend on the rows preceding them (JNSV 2017, Lemma 7.6), so they are
    never considered. Thus only r + m rows are processed, up to the rows of
    a chunk that follow a dependent row of the same index i.

    A chunk is a run of consecutive rows in the priority order, of at most
    `chunk` rows, with at most max(1, f_i) rows for each index i, where f_i
    is the number of rows of i found independent before the chunk: the
    rows of an index i which appears several times in a chunk are computed
    from each other before knowing whether they are independent, and this
    bounds the number of such rows which may be wasted.

    The span of the processed rows is kept as a basis W in reduced row
    echelon form, stored with permuted columns: the first r columns of the
    permutation are the pivot columns, so that W = [I_r | Wf]. A chunk V of
    rows (in the permuted column order) is reduced as Vf - Vp * Wf, where
    Vp and Vf are the pivot and free columns of V, then echelonized by
    _nmod_mat_rref_rows, which gives which rows are independent and the
    reduced row echelon form Y of the new rows, with pivots Q among the free
    columns; then Wf is updated as Wf - Wf[:, Q] * Y, Y is appended, and the
    columns Q are moved to the front of the free columns.

    If G is requested, the transformation W = G B is maintained (B being
    the basis rows, that is, the independent rows of the Krylov matrix):
    the chunk is echelonized with the identity matrix appended to it, which
    gives Y = H V - H Vp W = H V - H Vp G B for some H, supported on the
    independent rows of the chunk.

    In the row echelon mode (lumode, _nmod_mat_krylov_basis_iter_lu), the
    basis is W = [U | Wf] with U upper triangular, and its rows are never
    updated: a chunk is reduced as Vf - (Vp U^-1) Wf, by a triangular solve
    with the matrix on the right and one product, then decomposed by
    nmod_mat_ple, whose echelon rows are appended to W. The multipliers
    Vp U^-1 of the independent rows, followed by the unit lower triangular
    factor of the chunk, form the new rows of L, so that B[:, piv] = Pi^T L U
    for a row permutation Pi (given by bpos).
*/
static slong
_krylov_basis_iter(nmod_mat_t B, slong * rows, slong * exps,
        slong * deltas, nmod_mat_t T, nmod_mat_t G, int lumode, nmod_mat_t LU,
        slong * bpos, slong * piv, int * ok,
        const nmod_mat_t E, const _nmod_mat_krylov_op_t op,
        const slong * shift, const slong * degrees, slong chunk)
{
    slong m = E->r, n = E->c;
    ulong p = E->mod.n;
    slong i, j, a, t, u, k, kk, r, rb, nf, hlen, budget, maxstep, aug, logpow, powruns;
    slong * d, * key, * nextj, * last, * cnt, * prevpos, * heap;
    slong * ci, * cj, * cdepth, * csrc, * Q, * L, * colperm, * orig_at, * cur_of;
    slong * crun, * roff, * rlist, * bsrc, * clog, * cstep, * bprev;
    slong * Pc = NULL, * Qc = NULL, * bidx = NULL;
    char * state, * indep;
    int dense = (op->kind == _NMOD_MAT_KRYLOV_DENSE);
    nmod_mat_t V, Vp, W, Xf, Xa, Vpiv, Wf, Y, Wq, S, O;
    _krylov_pow_cache_struct C;

    if (op->n != n)
        flint_throw(FLINT_ERROR, "nmod_mat_krylov_basis: incompatible dimensions\n");

    if (ok != NULL)
        *ok = 1;

    if (chunk <= 0)
        chunk = KRYLOV_CHUNK_DEFAULT(n);

    /* a chunk has at most max(1, f_i) rows of index i, where the f_i sum
       to at most n */
    chunk = FLINT_MAX(1, FLINT_MIN(chunk, m + n));

    /* number of columns appended to the chunk to record row operations */
    aug = (G != NULL) ? chunk : 0;

    d = flint_malloc((7 * m + 1) * sizeof(slong));
    key = d + m;
    nextj = key + m;
    last = nextj + m;
    cnt = last + m;
    prevpos = cnt + m;
    heap = prevpos + m;
    state = flint_malloc(m + 1);

    ci = flint_malloc((12 * chunk + 4 * n + 1) * sizeof(slong));
    cj = ci + chunk;
    cdepth = cj + chunk;
    csrc = cdepth + chunk;
    Q = csrc + chunk;
    L = Q + chunk;
    crun = L + chunk;
    roff = crun + chunk;
    rlist = roff + chunk;
    bsrc = rlist + chunk;
    clog = bsrc + chunk;
    cstep = clog + chunk;
    colperm = cstep + chunk;
    orig_at = colperm + n;
    cur_of = orig_at + n;
    bprev = cur_of + n;
    indep = flint_malloc(chunk + 1);

    /* LU mode: PLE of the chunks, the pivot rows and the multipliers */
    if (lumode)
    {
        Pc = flint_malloc((2 * chunk + n + 1) * sizeof(slong));
        Qc = Pc + chunk;
        bidx = Qc + n;
    }

    /* powers of M for generating rows, at most 2^logpow <= chunk, in the
       chunks with at most powruns indices (0: no limit) */
    logpow = 0;
    powruns = 0;
    if (dense)
    {
        if (op->logpow < 0)
        {
            logpow = _krylov_auto_logpow(n, p);
            powruns = KRYLOV_POW_MAX_RUNS;
        }
        else
        {
            logpow = op->logpow;
            powruns = op->powruns;
        }

        logpow = FLINT_MIN(logpow, FLINT_BITS - 2);
        logpow = FLINT_MIN(logpow, (slong) FLINT_BIT_COUNT(chunk) - 1);
    }
    C.M = op->M;
    C.len = 0;

    _nmod_mat_krylov_degrees(d, m, degrees, n, 1);

    /* at most min(n, sum (d_i + 1)) rows are written in B */
    {
        slong N = 0;

        for (i = 0; i < m; i++)
            N += d[i] + 1;

        if (B->c != n || B->r < FLINT_MIN(n, N))
            flint_throw(FLINT_ERROR, "nmod_mat_krylov_basis: B must have E->c "
                    "columns and at least min(E->c, sum (d_i + 1)) rows\n");
    }

    /* priorities: keys s_i + j with equivalent shift in a small range */
    _nmod_mat_krylov_normalize_shift(key, shift, m, n);

    hlen = 0;
    for (i = 0; i < m; i++)
    {
        nextj[i] = 0;
        last[i] = -1;
        cnt[i] = 0;
        state[i] = ACTIVE;
        _nmod_mat_krylov_heap_push(heap, &hlen, i, key);
    }

    for (t = 0; t < n; t++)
        colperm[t] = t;

    nmod_mat_init(V, chunk, n, p);
    nmod_mat_init(Vp, chunk, n + aug, p);
    nmod_mat_init(W, n, n, p);

    r = 0;

    while (hlen > 0 && r < n)
    {
        /* ------------------------------------------------------------ */
        /* form the chunk                                               */
        /* ------------------------------------------------------------ */
        k = 0;

        while (k < chunk && hlen > 0)
        {
            i = heap[0];

            if (state[i] != ACTIVE)
            {
                _nmod_mat_krylov_heap_pop(heap, &hlen, key);
                continue;
            }

            /* nextj[i] - cnt[i] rows of i are in the profile */
            budget = FLINT_MAX(1, nextj[i] - cnt[i]);
            if (cnt[i] >= budget)
                break;

            _nmod_mat_krylov_heap_pop(heap, &hlen, key);

            ci[k] = i;
            cj[k] = nextj[i];
            cdepth[k] = cnt[i];
            csrc[k] = (cnt[i] == 0) ? -1 : prevpos[i];
            prevpos[i] = k;
            cnt[i]++;
            k++;

            nextj[i]++;
            if (nextj[i] <= d[i])
            {
                key[i]++;
                _nmod_mat_krylov_heap_push(heap, &hlen, i, key);
            }
        }

        /* run lengths: cnt[i] is the number of rows of i in the chunk */
        for (a = 0; a < k; a++)
            crun[a] = cnt[ci[a]];
        for (a = 0; a < k; a++)
            cnt[ci[a]] = 0;

        if (k == 0)
            break;

        /* ------------------------------------------------------------ */
        /* compute the rows of the chunk                                */
        /* ------------------------------------------------------------ */
        for (a = 0; a < k; a++)
            if (cj[a] == 0)
                memcpy(nmod_mat_entry_ptr(V, a, 0),
                       nmod_mat_entry_ptr(E, ci[a], 0), n * sizeof(ulong));

        if (dense)
        {
            /* row a, at position cdepth[a] of the run of its index, is
               computed from the row at distance e = 2^clog[a] in the run
               (or in the basis, for the first e rows of the run), at step
               cdepth[a] / e; rlist[roff[a] + u] is the chunk position of
               the row at position u of the run, and bsrc[roff[a] + u], for
               u < e, the basis position of the row at distance e before.
               The exponent of a run of c rows is the largest power of two
               at most min(2^lp, c) if it is the highest power of M
               available after raising it for the chunk, and 1 otherwise */
            slong nxt = 0, s, tt, c, e, lp = logpow;

            /* number of runs to generate */
            if (powruns > 0)
            {
                c = 0;
                for (a = 0; a < k; a++)
                    c += (cdepth[a] == 0 && cj[a] > 0);
                if (c > powruns)
                    lp = 0;
            }

            if (lp > 0)
            {
                slong tmax = 0;

                for (a = 0; a < k; a++)
                    if (cdepth[a] == 0 && cj[a] > 0)
                        tmax = FLINT_MAX(tmax, FLINT_MIN(lp,
                                    (slong) FLINT_BIT_COUNT(crun[a]) - 1));

                _krylov_pow_cache_raise(&C, tmax);
            }

            maxstep = 0;
            for (a = 0; a < k; a++)
            {
                if (cdepth[a] == 0)
                {
                    roff[a] = nxt;
                    nxt += crun[a];
                }
                else
                {
                    roff[a] = roff[csrc[a]];
                }

                rlist[roff[a] + cdepth[a]] = a;
                clog[a] = FLINT_MIN(lp, (slong) FLINT_BIT_COUNT(crun[a]) - 1);
                if (clog[a] != C.len)
                    clog[a] = 0;
                cstep[a] = cdepth[a] >> clog[a];
                if (cj[a] > 0)
                    maxstep = FLINT_MAX(maxstep, cstep[a]);

                /* the basis rows (i, j0 - e), ..., (i, j0 - 1) before the
                   first row (i, j0) of the run (e <= crun[a] <= j0) */
                if (cdepth[a] == 0 && cj[a] > 0)
                {
                    slong pos = last[ci[a]];

                    e = WORD(1) << clog[a];
                    for (u = e - 1; u >= 0; u--)
                    {
                        bsrc[roff[a] + u] = pos;
                        if (u > 0)
                            pos = bprev[pos];
                    }
                }
            }

            /* one product per step and exponent */
            for (s = 0; s <= maxstep; s++)
            {
                for (tt = 0; tt <= logpow; tt++)
                {
                    c = 0;
                    for (a = 0; a < k; a++)
                        if (cj[a] > 0 && cstep[a] == s && clog[a] == tt)
                            L[c++] = a;

                    if (c == 0)
                        continue;

                    e = WORD(1) << tt;
                    nmod_mat_init(S, c, n, p);
                    nmod_mat_init(O, c, n, p);

                    for (u = 0; u < c; u++)
                    {
                        a = L[u];
                        memcpy(nmod_mat_entry_ptr(S, u, 0), (cdepth[a] < e)
                            ? nmod_mat_entry_ptr(B, bsrc[roff[a] + cdepth[a]], 0)
                            : nmod_mat_entry_ptr(V, rlist[roff[a] + cdepth[a] - e], 0),
                            n * sizeof(ulong));
                    }

                    nmod_mat_mul(O, S, _krylov_pow_cache_get(&C, tt));

                    for (u = 0; u < c; u++)
                        memcpy(nmod_mat_entry_ptr(V, L[u], 0),
                               nmod_mat_entry_ptr(O, u, 0), n * sizeof(ulong));

                    nmod_mat_clear(S);
                    nmod_mat_clear(O);
                }
            }
        }
        else
        {
            for (a = 0; a < k; a++)
            {
                if (cj[a] == 0)
                    continue;
                _nmod_mat_krylov_op_apply_vec(nmod_mat_entry_ptr(V, a, 0),
                    (cdepth[a] == 0) ? nmod_mat_entry_ptr(B, last[ci[a]], 0)
                                     : nmod_mat_entry_ptr(V, csrc[a], 0), op);
            }
        }

        /* ------------------------------------------------------------ */
        /* reduce by the basis, in the permuted column order            */
        /* ------------------------------------------------------------ */
        nf = n - r;

        for (a = 0; a < k; a++)
        {
            nn_srcptr v = nmod_mat_entry_ptr(V, a, 0);
            nn_ptr w = nmod_mat_entry_ptr(Vp, a, 0);
            for (t = 0; t < n; t++)
                w[t] = v[colperm[t]];
            if (aug)
            {
                for (t = 0; t < k; t++)
                    w[n + t] = 0;
                w[n + a] = 1;
            }
        }

        nmod_mat_window_init(Vpiv, Vp, 0, 0, k, r);

        if (r > 0)
        {
            nmod_mat_window_init(Xf, Vp, 0, r, k, n);
            nmod_mat_window_init(Wf, W, 0, r, r, n);
            if (lumode)
            {
                /* W = [U | Wf], U upper triangular: Vpiv <- Vpiv U^-1,
                   the coefficients of the chunk rows on the rows of W */
                nmod_mat_t Uw;
                nmod_mat_window_init(Uw, W, 0, 0, r, r);
                nmod_mat_solve_triu_right(Vpiv, Uw, Vpiv, 0);
                nmod_mat_window_clear(Uw);
            }
            nmod_mat_submul(Xf, Xf, Vpiv, Wf);
            nmod_mat_window_clear(Xf);
            nmod_mat_window_clear(Wf);
        }

        /* ------------------------------------------------------------ */
        /* echelonize the chunk                                         */
        /* ------------------------------------------------------------ */
        if (lumode)
        {
            nmod_mat_window_init(Xa, Vp, 0, r, k, n);
            kk = nmod_mat_ple(Pc, Qc, Xa, 0);
            nmod_mat_window_clear(Xa);
            memset(indep, 0, k);
            for (t = 0; t < kk; t++)
            {
                indep[Pc[t]] = 1;
                Q[t] = Qc[t];
            }
        }
        else
        {
            nmod_mat_window_init(Xa, Vp, 0, r, k, n + (aug ? k : 0));
            kk = _nmod_mat_rref_rows(Xa, nf, Q, indep);
            nmod_mat_window_clear(Xa);
        }

        /* ------------------------------------------------------------ */
        /* record independent and dependent rows                        */
        /* ------------------------------------------------------------ */
        rb = r;
        for (a = 0; a < k; a++)
        {
            i = ci[a];

            if (indep[a])
            {
                memcpy(nmod_mat_entry_ptr(B, rb, 0),
                       nmod_mat_entry_ptr(V, a, 0), n * sizeof(ulong));
                rows[rb] = i;
                exps[rb] = cj[a];
                L[rb - r] = a;   /* chunk position of the new basis row */
                if (bidx != NULL)
                    bidx[a] = rb;
                bprev[rb] = -1;

                if (state[i] == ACTIVE)
                {
                    bprev[rb] = last[i];
                    last[i] = rb;

                    if (cj[a] == d[i])
                    {
                        state[i] = EXHAUSTED;
                        if (ok != NULL)
                            *ok = 0;
                        if (deltas != NULL)
                            deltas[i] = d[i] + 1;
                    }
                }
                else if (ok != NULL)
                {
                    /* a row following a dependent row of the same index
                       can only be independent if the degree bounds are
                       not valid, which has been detected before */
                    *ok = 0;
                }

                rb++;
            }
            else if (state[i] == ACTIVE)
            {
                state[i] = DEAD;
                if (deltas != NULL)
                    deltas[i] = cj[a];
                if (T != NULL)
                    memcpy(nmod_mat_entry_ptr(T, i, 0),
                           nmod_mat_entry_ptr(V, a, 0), n * sizeof(ulong));
            }
        }

        /* ------------------------------------------------------------ */
        /* update the basis                                             */
        /* ------------------------------------------------------------ */
        if (kk > 0 && lumode)
        {
            /* new rows of W: the echelon form of the chunk (zero before the
               pivot, where the PLE output holds multipliers); new rows of
               the multipliers: those on the old rows of W (Vpiv, whose rows
               the PLE did not move) and the unit lower triangular factor
               of the chunk */
            for (t = 0; t < kk; t++)
            {
                nn_ptr w = nmod_mat_entry_ptr(W, r + t, 0);
                nn_srcptr x = nmod_mat_entry_ptr(Vp, t, 0);

                for (u = 0; u < r; u++)
                    w[u] = 0;
                for (u = 0; u < nf; u++)
                    w[r + u] = (u >= Qc[t]) ? x[r + u] : 0;

                if (LU != NULL)
                {
                    nn_ptr l = nmod_mat_entry_ptr(LU, r + t, 0);
                    memcpy(l, nmod_mat_entry_ptr(Vp, Pc[t], 0), r * sizeof(ulong));
                    for (u = 0; u < t; u++)
                        l[r + u] = x[r + u];
                    bpos[r + t] = bidx[Pc[t]];
                }
            }
        }

        if (kk > 0 && !lumode)
        {
            nmod_mat_window_init(Y, Vp, 0, r, kk, n);
            nmod_mat_init(Wq, r, kk, p);

            if (r > 0)
            {
                nmod_mat_window_init(Wf, W, 0, r, r, n);
                _gather_cols(Wq, Wf, Q, r, kk);
            }

            if (G != NULL)
            {
                /* coefficients of the new rows Y: on the chunk rows,
                   H = Vp[0:kk, n:n+k], restricted to the independent rows
                   L[0..kk), and on the old basis rows, -H Vpiv G */
                nmod_mat_t Hc, Hv, Yg, Gw, Gw2;
                slong t2;

                nmod_mat_window_init(Gw, G, r, r, r + kk, r + kk);
                for (t = 0; t < kk; t++)
                    for (t2 = 0; t2 < kk; t2++)
                        nmod_mat_entry(Gw, t, t2) = nmod_mat_entry(Vp, t, n + L[t2]);

                if (r > 0)
                {
                    nmod_mat_window_init(Hc, Vp, 0, n, kk, n + k);
                    nmod_mat_init(Hv, kk, r, p);
                    nmod_mat_mul(Hv, Hc, Vpiv);
                    nmod_mat_window_clear(Hc);

                    nmod_mat_init(Yg, kk, r, p);
                    nmod_mat_window_init(Gw2, G, 0, 0, r, r);
                    nmod_mat_mul(Yg, Hv, Gw2);
                    nmod_mat_neg(Yg, Yg);
                    nmod_mat_clear(Hv);

                    /* old rows: G[0:r, 0:r] -= Wq Yg,
                       G[0:r, r:r+kk] = -Wq Gw */
                    nmod_mat_submul(Gw2, Gw2, Wq, Yg);
                    nmod_mat_window_clear(Gw2);

                    nmod_mat_window_init(Gw2, G, 0, r, r, r + kk);
                    nmod_mat_mul(Gw2, Wq, Gw);
                    nmod_mat_neg(Gw2, Gw2);
                    nmod_mat_window_clear(Gw2);

                    /* new rows: G[r:r+kk, 0:r] = Yg */
                    nmod_mat_window_init(Gw2, G, r, 0, r + kk, r);
                    nmod_mat_set(Gw2, Yg);
                    nmod_mat_window_clear(Gw2);
                    nmod_mat_clear(Yg);
                }

                nmod_mat_window_clear(Gw);
            }

            if (r > 0)
            {
                /* Wf = Wf - Wf[:, Q] * Y */
                nmod_mat_submul(Wf, Wf, Wq, Y);
                nmod_mat_window_clear(Wf);
            }

            nmod_mat_clear(Wq);

            for (t = 0; t < kk; t++)
                memcpy(nmod_mat_entry_ptr(W, r + t, r),
                       nmod_mat_entry_ptr(Y, t, 0), nf * sizeof(ulong));

            nmod_mat_window_clear(Y);
        }

        if (kk > 0)
        {
            /* move the free columns Q[0], ..., Q[kk-1] to the front */
            for (t = 0; t < nf; t++)
            {
                orig_at[t] = t;
                cur_of[t] = t;
            }

            for (t = 0; t < kk; t++)
            {
                slong pt = cur_of[Q[t]], o;

                if (pt != t)
                {
                    for (u = 0; u < r + kk; u++)
                    {
                        nn_ptr w = nmod_mat_entry_ptr(W, u, r);
                        ulong tmp = w[t];
                        w[t] = w[pt];
                        w[pt] = tmp;
                    }

                    o = colperm[r + t];
                    colperm[r + t] = colperm[r + pt];
                    colperm[r + pt] = o;

                    o = orig_at[t];
                    orig_at[t] = Q[t];
                    orig_at[pt] = o;
                    cur_of[o] = pt;
                    cur_of[Q[t]] = t;
                }
            }

            r += kk;
        }

        nmod_mat_window_clear(Vpiv);
    }

    /* ---------------------------------------------------------------- */
    /* full rank: the next row of each active index is dependent        */
    /* ---------------------------------------------------------------- */
    if (r == n)
    {
        slong c = 0;

        for (i = 0; i < m; i++)
        {
            if (state[i] != ACTIVE)
                continue;

            j = nextj[i];
            state[i] = DEAD;
            if (deltas != NULL)
                deltas[i] = j;

            if (T != NULL && n > 0)
            {
                if (j == 0)
                    memcpy(nmod_mat_entry_ptr(T, i, 0),
                           nmod_mat_entry_ptr(E, i, 0), n * sizeof(ulong));
                else if (!dense)
                    _nmod_mat_krylov_op_apply_vec(nmod_mat_entry_ptr(T, i, 0),
                            nmod_mat_entry_ptr(B, last[i], 0), op);
                else
                    heap[c++] = i;  /* reuse heap as a list */
            }
        }

        if (c > 0)
        {
            nmod_mat_init(S, c, n, p);
            nmod_mat_init(O, c, n, p);
            for (u = 0; u < c; u++)
                memcpy(nmod_mat_entry_ptr(S, u, 0),
                       nmod_mat_entry_ptr(B, last[heap[u]], 0), n * sizeof(ulong));
            nmod_mat_mul(O, S, op->M);
            for (u = 0; u < c; u++)
                memcpy(nmod_mat_entry_ptr(T, heap[u], 0),
                       nmod_mat_entry_ptr(O, u, 0), n * sizeof(ulong));
            nmod_mat_clear(S);
            nmod_mat_clear(O);
        }
    }

    if (piv != NULL)
        for (t = 0; t < r; t++)
            piv[t] = colperm[t];

    /* U: the leading r x r block of W, upper triangular */
    if (LU != NULL)
        for (t = 0; t < r; t++)
            memcpy(nmod_mat_entry_ptr(LU, t, t), nmod_mat_entry_ptr(W, t, t),
                   (r - t) * sizeof(ulong));

    _krylov_pow_cache_clear(&C);
    nmod_mat_clear(V);
    nmod_mat_clear(Vp);
    nmod_mat_clear(W);
    flint_free(d);
    flint_free(state);
    flint_free(ci);
    flint_free(indep);
    if (Pc != NULL)
        flint_free(Pc);

    return r;
}

slong
_nmod_mat_krylov_basis_iter(nmod_mat_t B, slong * rows, slong * exps,
        slong * deltas, nmod_mat_t T, nmod_mat_t G, slong * piv, int * ok,
        const nmod_mat_t E, const _nmod_mat_krylov_op_t op,
        const slong * shift, const slong * degrees, slong chunk)
{
    return _krylov_basis_iter(B, rows, exps, deltas, T, G, 0, NULL, NULL, piv, ok,
                              E, op, shift, degrees, chunk);
}

slong
_nmod_mat_krylov_basis_iter_lu(nmod_mat_t B, slong * rows, slong * exps,
        slong * deltas, nmod_mat_t T, nmod_mat_t LU, slong * bpos, slong * piv,
        int * ok, const nmod_mat_t E, const _nmod_mat_krylov_op_t op,
        const slong * shift, const slong * degrees, slong chunk)
{
    return _krylov_basis_iter(B, rows, exps, deltas, T, NULL, 1, LU, bpos, piv, ok,
                              E, op, shift, degrees, chunk);
}

slong
_nmod_mat_krylov_basis_naive(nmod_mat_t B, slong * rows, slong * exps,
        const nmod_mat_t E, const _nmod_mat_krylov_op_t op,
        const slong * shift, const slong * degrees)
{
    slong m = E->r, n = E->c, N, a, r, t;
    slong * d, * ri, * rj;
    char * indep;
    nmod_mat_t K, X;

    d = flint_malloc((m + 1) * sizeof(slong));
    _nmod_mat_krylov_degrees(d, m, degrees, n, 0);
    N = _nmod_mat_krylov_nrows(m, d);

    ri = flint_malloc((2 * N + 1) * sizeof(slong));
    rj = ri + N;
    indep = flint_malloc(N + 1);

    nmod_mat_init(K, N, n, E->mod.n);
    _nmod_mat_krylov_matrix(K, ri, rj, E, op, shift, d);

    nmod_mat_init_set(X, K);
    r = _nmod_mat_row_rank_profile(indep, X);

    t = 0;
    for (a = 0; a < N; a++)
    {
        if (indep[a])
        {
            memcpy(nmod_mat_entry_ptr(B, t, 0), nmod_mat_entry_ptr(K, a, 0),
                   n * sizeof(ulong));
            rows[t] = ri[a];
            exps[t] = rj[a];
            t++;
        }
    }

    nmod_mat_clear(K);
    nmod_mat_clear(X);
    flint_free(d);
    flint_free(ri);
    flint_free(indep);

    return r;
}

/*
    The basis alone is computed in row echelon form (LU mode, see
    _nmod_mat_krylov_basis_iter_lu) rather than in reduced form: measured
    with uniform and Hermite shifts, Hermite-Pade and dense M, n = 64 to
    1024 (Nov 2026, Zen 4, Ice Lake, Cascade Lake, Apple M4), it takes 0.54
    to 1.03 times the time, except at n = 128 on x86 below 2^52 (0.78 to
    1.19 at 30 bits, 1.01 to 1.06 at 50 bits).
*/
slong
nmod_mat_krylov_basis(nmod_mat_t B, slong * rows, slong * exps,
        const nmod_mat_t E, const nmod_mat_t M,
        const slong * shift, const slong * degrees)
{
    _nmod_mat_krylov_op_t op;
    slong r;

    if (M->r != E->c || M->c != E->c)
        flint_throw(FLINT_ERROR, "nmod_mat_krylov_basis: incompatible dimensions\n");

    _nmod_mat_krylov_op_init_dense(op, M);
    r = _nmod_mat_krylov_basis_iter_lu(B, rows, exps, NULL, NULL, NULL, NULL,
                                       NULL, NULL, E, op, shift, degrees, 0);
    _nmod_mat_krylov_op_clear(op);

    return r;
}

slong
nmod_mat_krylov_basis_jordan(nmod_mat_t B, slong * rows, slong * exps,
        const nmod_mat_t E, nn_srcptr eig, const slong * sizes, slong nblocks,
        const slong * shift, const slong * degrees)
{
    _nmod_mat_krylov_op_t op;
    slong r;

    _nmod_mat_krylov_op_init_jordan(op, eig, sizes, nblocks, E->mod);

    if (op->n != E->c)
        flint_throw(FLINT_ERROR, "nmod_mat_krylov_basis_jordan: incompatible dimensions\n");

    r = _nmod_mat_krylov_basis_iter_lu(B, rows, exps, NULL, NULL, NULL, NULL,
                                       NULL, NULL, E, op, shift, degrees, 0);
    _nmod_mat_krylov_op_clear(op);

    return r;
}
