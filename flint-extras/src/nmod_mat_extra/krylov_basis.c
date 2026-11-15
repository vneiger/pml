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

/* default maximal number of rows processed at once: D / 8 (D = number of
   columns), clamped to [16, 128] (profile p-krylov chunk; on Zen 4, Ice
   Lake, Cascade Lake and Apple M4, Nov 2026, within 1.05 times the best
   size except at 60 bits for D = 128: 1.08 to 1.17) */
#define KRYLOV_CHUNK_DEFAULT(D) FLINT_MAX(16, FLINT_MIN(128, (D) / 8))

/*
    Automatic choice of the powers of M used to generate the rows (dense
    M, see krylov_impl.h): for D >= KRYLOV_POW_MIN_DIM, where the products
    of a few rows by M are limited by the memory bandwidth, and only in the
    chunks where KRYLOV_POW_MAX_RUNS = 1 index generates rows (m = 1, or a
    Hermite-type shift; with more indices, the rows are already computed by
    products with several rows). The squarings and the products of a few
    rows are done by nmod_mat_mul, so the gain depends on its speed for the
    modulus: M^2 and M^4 (logpow 2) for moduli below 2^32, and below 2^52
    with AVX-512 IFMA (FLINT's IFMA products).

    Measured at D = 1600, dense M, relative to products by M only (profile
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
_krylov_auto_logpow(slong D, ulong p)
{
    flint_bitcnt_t bits = FLINT_BIT_COUNT(p);

    if (D < KRYLOV_POW_MIN_DIM)
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
#define EXHAUSTED 2  /* the row (i, beta_i) is independent */

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
   these are kept, to limit the memory to one more D x D matrix */
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

    The rows (i, e) are processed in increasing priority (s_i + e, i), by
    chunks of consecutive rows. As soon as a row (i, e) is found to depend
    on the previous ones, all the rows (i, e') with e' > e are known to
    depend on the rows preceding them (JNSV 2017, Lemma 7.6), so they are
    never considered. Thus only r + m rows are processed, up to the rows of
    a chunk that follow a dependent row of the same index i.

    A chunk is a run of consecutive rows in the priority order, of at most
    `chunk` rows, with a budget of rows for each index i: the rows of an
    index which appears several times in a chunk are computed from each
    other before knowing whether they are independent, and the budget
    bounds the number of such rows which may be wasted. Let f_i be the
    number of rows of i found independent before the chunk, and p_i the
    number of rows (i, e), e < beta_i, among the first D rows in priority
    order. For generic inputs (M cyclic and V random), these D rows are
    independent whatever the shift (their determinant is a nonzero
    polynomial in the entries of V), so they form the profile and
    p_i = delta_i. As long as no index has stopped before its p_i rows, the
    budget of i is max(1, f_i, p_i - f_i) (once f_i >= 2^logpow when rows
    come from powers of M), so that the chunks are full and generically no
    row is wasted; afterwards, it is max(1, f_i). Hence fewer than
    r + chunk rows are wasted in total: the budget max(1, f_i) wastes fewer
    than delta_i rows of i, the rows taken beyond it have e < p_i and are
    not wasted as long as no index stops early, and fewer than `chunk` rows
    are wasted in the chunk where one does.

    Other strategy measured and not kept: taking the first D rows in priority
    order as one speculative chunk, eliminated by one nmod_mat_ple. Correct,
    and no work is lost when the guess fails (the chunk is a prefix of the
    priority order, so its independent rows and the first dependent row of each
    index are final, and the iteration resumes from there), but 0.91 to 1.03
    times as fast as the chunks for the kernel: one PLE of D rows is about 1.2
    times slower than the elimination by chunks (0.6 to 1.0 against 0.5 to 0.9
    products at 30 bits, D = 256 to 1024), and the D x D chunk costs copies. It
    only gains in corner cases, e.g. for one nilpotent block driven by one
    index at 60 bits (1.7 to 5 times, as Crout's PLE skips the zeros of the
    triangular Krylov matrix).

    The span of the processed rows is kept as a basis W in reduced row
    echelon form, stored with permuted columns: the first r columns of the
    permutation are the pivot columns, so that W = [I_r | Wf]. A chunk X of
    rows (in the permuted column order) is reduced as Xf - Xpiv * Wf, where
    Xpiv and Xf are the pivot and free columns of X, then echelonized by
    _nmod_mat_rref_rows, which gives which rows are independent and the
    reduced row echelon form Y of the new rows, with pivots Q among the free
    columns; then Wf is updated as Wf - Wf[:, Q] * Y, Y is appended, and the
    columns Q are moved to the front of the free columns.

    If G is requested, the transformation W = G B is maintained (B being
    the basis rows, that is, the independent rows of the Krylov matrix):
    the chunk is echelonized with the identity matrix appended to it, which
    gives Y = H X - H Xpiv W = H X - H Xpiv G B for some H, supported on the
    independent rows of the chunk.

    In the row echelon mode (lumode, _nmod_mat_krylov_basis_iter_lu), the
    basis is W = [U | Wf] with U upper triangular, and its rows are never
    updated: a chunk is reduced as Xf - (Xpiv U^-1) Wf, by a triangular solve
    with the matrix on the right and one product, then decomposed by
    nmod_mat_ple, whose echelon rows are appended to W. The multipliers
    Xpiv U^-1 of the independent rows, followed by the unit lower triangular
    factor of the chunk, form the new rows of L, so that B[:, piv] = Pi^T L U
    for a row permutation Pi (given by bpos).
*/
static slong
_krylov_basis_iter(nmod_mat_t B, slong * rows, slong * exps,
        slong * deltas, nmod_mat_t T, nmod_mat_t G, int lumode, nmod_mat_t LU,
        slong * bpos, slong * piv, int * ok,
        const nmod_mat_t V, const _nmod_mat_krylov_op_t op,
        const slong * shift, const slong * degrees, slong chunk)
{
    slong m = V->r, D = V->c;
    ulong p = V->mod.n;
    slong i, e, a, t, u, k, kk, r, rb, nf, hlen, budget, maxstep, aug, logpow, powruns;
    slong * beta, * key, * nexte, * last, * cnt, * prevpos, * heap, * pcnt;
    slong pmin;
    int predict = op->predict;
    slong * ci, * ce, * cdepth, * csrc, * Q, * L, * colperm, * orig_at, * cur_of;
    slong * crun, * roff, * rlist, * bsrc, * clog, * cstep, * bprev;
    slong * Pc = NULL, * Qc = NULL, * bidx = NULL;
    char * state, * indep;
    int dense = (op->kind == _NMOD_MAT_KRYLOV_DENSE);
    nmod_mat_t X, Xp, W, Xf, Xa, Xpiv, Wf, Y, Wq, S, O;
    _krylov_pow_cache_struct C;

    if (op->D != D)
        flint_throw(FLINT_ERROR, "nmod_mat_krylov_basis: incompatible dimensions\n");

    if (ok != NULL)
        *ok = 1;

    if (chunk <= 0)
        chunk = KRYLOV_CHUNK_DEFAULT(D);

    /* larger chunks are not useful: at most r + m <= D + m rows are
       needed */
    chunk = FLINT_MAX(1, FLINT_MIN(chunk, m + D));

    /* number of columns appended to the chunk to record row operations */
    aug = (G != NULL) ? chunk : 0;

    beta = flint_malloc((8 * m + 1) * sizeof(slong));
    key = beta + m;
    nexte = key + m;
    last = nexte + m;
    cnt = last + m;
    prevpos = cnt + m;
    heap = prevpos + m;
    pcnt = heap + m;
    state = flint_malloc(m + 1);

    ci = flint_malloc((12 * chunk + 4 * D + 1) * sizeof(slong));
    ce = ci + chunk;
    cdepth = ce + chunk;
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
    orig_at = colperm + D;
    cur_of = orig_at + D;
    bprev = cur_of + D;
    indep = flint_malloc(chunk + 1);

    /* LU mode: PLE of the chunks, the pivot rows and the multipliers */
    if (lumode)
    {
        Pc = flint_malloc((2 * chunk + D + 1) * sizeof(slong));
        Qc = Pc + chunk;
        bidx = Qc + D;
    }

    /* powers of M for generating rows, at most 2^logpow <= chunk, in the
       chunks with at most powruns indices (0: no limit) */
    logpow = 0;
    powruns = 0;
    if (dense)
    {
        if (op->logpow < 0)
        {
            logpow = _krylov_auto_logpow(D, p);
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

    /* with powers of M, a run of rows (i, e0), ... is generated from the
       rows (i, e0 - q), ... with q <= e0: the larger budget only starts
       once f_i >= 2^logpow, so that the run can use the power (measured:
       otherwise 2 to 6 percent slower with one index at D = 1600) */
    pmin = (logpow > 0) ? (WORD(1) << logpow) : 0;

    _nmod_mat_krylov_degrees(beta, m, degrees, D, 1);

    /* at most min(D, sum (beta_i + 1)) rows are written in B */
    {
        slong N = 0;

        for (i = 0; i < m; i++)
            N += beta[i] + 1;

        if (B->c != D || B->r < FLINT_MIN(D, N))
            flint_throw(FLINT_ERROR, "nmod_mat_krylov_basis: B must have V->c "
                    "columns and at least min(V->c, sum (degrees[i] + 1)) rows\n");
    }

    /* V_i M^e depends on the rows V_i M^e' with e' < e as soon as e is
       the degree of the minimal polynomial of M (op->mpdeg, which is D for
       a dense operator): the bounds min(beta_i, mpdeg) give the same
       profile, and are valid if the beta_i are. For a Jordan matrix which
       is not cyclic, this saves the rows of an index beyond mpdeg
       (measured: 4 nilpotent blocks of size D / 4 with a Hermite shift,
       1.16 to 1.21 times faster for the kernel, 1.30 to 1.34 for the
       basis) */
    for (i = 0; i < m; i++)
        beta[i] = FLINT_MIN(beta[i], op->mpdeg);

    /* priorities: keys s_i + e with equivalent shift in a small range */
    _nmod_mat_krylov_normalize_shift(key, shift, m, D);

    /* predicted number of rows p_i of each index: its rows (i, e), e < beta_i,
       among the first D rows in priority order (the heap is used as
       workspace, with keys in nexte) */
    if (predict)
    {
        slong tot = 0;

        hlen = 0;
        for (i = 0; i < m; i++)
        {
            pcnt[i] = 0;
            nexte[i] = key[i];
            if (beta[i] > 0)
                _nmod_mat_krylov_heap_push(heap, &hlen, i, nexte);
        }

        while (tot < D && hlen > 0)
        {
            i = _nmod_mat_krylov_heap_pop(heap, &hlen, nexte);
            pcnt[i]++;
            tot++;
            if (pcnt[i] < beta[i])
            {
                nexte[i]++;
                _nmod_mat_krylov_heap_push(heap, &hlen, i, nexte);
            }
        }
    }

    hlen = 0;
    for (i = 0; i < m; i++)
    {
        nexte[i] = 0;
        last[i] = -1;
        cnt[i] = 0;
        state[i] = ACTIVE;
        _nmod_mat_krylov_heap_push(heap, &hlen, i, key);
    }

    for (t = 0; t < D; t++)
        colperm[t] = t;

    nmod_mat_init(X, chunk, D, p);
    nmod_mat_init(Xp, chunk, D + aug, p);
    nmod_mat_init(W, D, D, p);

    r = 0;

    while (hlen > 0 && r < D)
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

            /* f_i = nexte[i] - cnt[i] rows of i are in the profile; the
               budget is max(1, f_i), or max(1, f_i, p_i - f_i) while the
               prediction holds */
            budget = FLINT_MAX(1, nexte[i] - cnt[i]);
            if (predict && nexte[i] - cnt[i] >= pmin)
                budget = FLINT_MAX(budget, pcnt[i] - (nexte[i] - cnt[i]));
            if (cnt[i] >= budget)
                break;

            _nmod_mat_krylov_heap_pop(heap, &hlen, key);

            ci[k] = i;
            ce[k] = nexte[i];
            cdepth[k] = cnt[i];
            csrc[k] = (cnt[i] == 0) ? -1 : prevpos[i];
            prevpos[i] = k;
            cnt[i]++;
            k++;

            nexte[i]++;
            if (nexte[i] <= beta[i])
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
            if (ce[a] == 0)
                memcpy(nmod_mat_entry_ptr(X, a, 0),
                       nmod_mat_entry_ptr(V, ci[a], 0), D * sizeof(ulong));

        if (dense)
        {
            /* row a, at position cdepth[a] of the run of its index, is
               computed from the row at distance q = 2^clog[a] in the run
               (or in the basis, for the first q rows of the run), at step
               cdepth[a] / q; rlist[roff[a] + u] is the chunk position of
               the row at position u of the run, and bsrc[roff[a] + u], for
               u < q, the basis position of the row at distance q before.
               The exponent of a run of c rows starting at (i, e0) is the
               largest power of two at most min(2^lp, c, e0) if it is the
               highest power of M available after raising it for the
               chunk, and 1 otherwise (e0 < c is possible with the
               predicted budget) */
            slong nxt = 0, s, tt, c, q, lp = logpow;

            /* number of runs to generate */
            if (powruns > 0)
            {
                c = 0;
                for (a = 0; a < k; a++)
                    c += (cdepth[a] == 0 && ce[a] > 0);
                if (c > powruns)
                    lp = 0;
            }

            if (lp > 0)
            {
                slong tmax = 0;

                for (a = 0; a < k; a++)
                    if (cdepth[a] == 0 && ce[a] > 0)
                        tmax = FLINT_MAX(tmax, FLINT_MIN(lp, (slong)
                                    FLINT_BIT_COUNT(FLINT_MIN(crun[a], ce[a])) - 1));

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
                clog[a] = FLINT_MAX(0, FLINT_MIN(lp, (slong) FLINT_BIT_COUNT(
                            FLINT_MIN(crun[a], ce[a] - cdepth[a])) - 1));
                if (clog[a] != C.len)
                    clog[a] = 0;
                cstep[a] = cdepth[a] >> clog[a];
                if (ce[a] > 0)
                    maxstep = FLINT_MAX(maxstep, cstep[a]);

                /* the basis rows (i, e0 - q), ..., (i, e0 - 1) before the
                   first row (i, e0) of the run (q <= min(crun[a], e0)) */
                if (cdepth[a] == 0 && ce[a] > 0)
                {
                    slong pos = last[ci[a]];

                    q = WORD(1) << clog[a];
                    for (u = q - 1; u >= 0; u--)
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
                        if (ce[a] > 0 && cstep[a] == s && clog[a] == tt)
                            L[c++] = a;

                    if (c == 0)
                        continue;

                    q = WORD(1) << tt;
                    nmod_mat_init(S, c, D, p);
                    nmod_mat_init(O, c, D, p);

                    for (u = 0; u < c; u++)
                    {
                        a = L[u];
                        memcpy(nmod_mat_entry_ptr(S, u, 0), (cdepth[a] < q)
                            ? nmod_mat_entry_ptr(B, bsrc[roff[a] + cdepth[a]], 0)
                            : nmod_mat_entry_ptr(X, rlist[roff[a] + cdepth[a] - q], 0),
                            D * sizeof(ulong));
                    }

                    nmod_mat_mul(O, S, _krylov_pow_cache_get(&C, tt));

                    for (u = 0; u < c; u++)
                        memcpy(nmod_mat_entry_ptr(X, L[u], 0),
                               nmod_mat_entry_ptr(O, u, 0), D * sizeof(ulong));

                    nmod_mat_clear(S);
                    nmod_mat_clear(O);
                }
            }
        }
        else
        {
            for (a = 0; a < k; a++)
            {
                if (ce[a] == 0)
                    continue;
                _nmod_mat_krylov_op_apply_vec(nmod_mat_entry_ptr(X, a, 0),
                    (cdepth[a] == 0) ? nmod_mat_entry_ptr(B, last[ci[a]], 0)
                                     : nmod_mat_entry_ptr(X, csrc[a], 0), op);
            }
        }

        /* ------------------------------------------------------------ */
        /* reduce by the basis, in the permuted column order            */
        /* ------------------------------------------------------------ */
        nf = D - r;

        for (a = 0; a < k; a++)
        {
            nn_srcptr v = nmod_mat_entry_ptr(X, a, 0);
            nn_ptr w = nmod_mat_entry_ptr(Xp, a, 0);
            for (t = 0; t < D; t++)
                w[t] = v[colperm[t]];
            if (aug)
            {
                for (t = 0; t < k; t++)
                    w[D + t] = 0;
                w[D + a] = 1;
            }
        }

        nmod_mat_window_init(Xpiv, Xp, 0, 0, k, r);

        if (r > 0)
        {
            nmod_mat_window_init(Xf, Xp, 0, r, k, D);
            nmod_mat_window_init(Wf, W, 0, r, r, D);
            if (lumode)
            {
                /* W = [U | Wf], U upper triangular: Xpiv <- Xpiv U^-1,
                   the coefficients of the chunk rows on the rows of W */
                nmod_mat_t Uw;
                nmod_mat_window_init(Uw, W, 0, 0, r, r);
                nmod_mat_solve_triu_right(Xpiv, Uw, Xpiv, 0);
                nmod_mat_window_clear(Uw);
            }
            nmod_mat_submul(Xf, Xf, Xpiv, Wf);
            nmod_mat_window_clear(Xf);
            nmod_mat_window_clear(Wf);
        }

        /* ------------------------------------------------------------ */
        /* echelonize the chunk                                         */
        /* ------------------------------------------------------------ */
        if (lumode)
        {
            nmod_mat_window_init(Xa, Xp, 0, r, k, D);
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
            nmod_mat_window_init(Xa, Xp, 0, r, k, D + (aug ? k : 0));
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
                       nmod_mat_entry_ptr(X, a, 0), D * sizeof(ulong));
                rows[rb] = i;
                exps[rb] = ce[a];
                L[rb - r] = a;   /* chunk position of the new basis row */
                if (bidx != NULL)
                    bidx[a] = rb;
                bprev[rb] = -1;

                if (state[i] == ACTIVE)
                {
                    bprev[rb] = last[i];
                    last[i] = rb;

                    if (ce[a] == beta[i])
                    {
                        state[i] = EXHAUSTED;
                        if (ok != NULL)
                            *ok = 0;
                        if (deltas != NULL)
                            deltas[i] = beta[i] + 1;
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
                /* the prediction fails: the budgets go back to max(1, f_i) */
                if (predict == 1 && ce[a] < pcnt[i])
                    predict = 0;
                if (deltas != NULL)
                    deltas[i] = ce[a];
                if (T != NULL)
                    memcpy(nmod_mat_entry_ptr(T, i, 0),
                           nmod_mat_entry_ptr(X, a, 0), D * sizeof(ulong));
            }
        }

        /* ------------------------------------------------------------ */
        /* update the basis                                             */
        /* ------------------------------------------------------------ */
        if (kk > 0 && lumode)
        {
            /* new rows of W: the echelon form of the chunk (zero before the
               pivot, where the PLE output holds multipliers); new rows of
               the multipliers: those on the old rows of W (Xpiv, whose rows
               the PLE did not move) and the unit lower triangular factor
               of the chunk */
            for (t = 0; t < kk; t++)
            {
                nn_ptr w = nmod_mat_entry_ptr(W, r + t, 0);
                nn_srcptr x = nmod_mat_entry_ptr(Xp, t, 0);

                for (u = 0; u < r; u++)
                    w[u] = 0;
                for (u = 0; u < nf; u++)
                    w[r + u] = (u >= Qc[t]) ? x[r + u] : 0;

                if (LU != NULL)
                {
                    nn_ptr l = nmod_mat_entry_ptr(LU, r + t, 0);
                    memcpy(l, nmod_mat_entry_ptr(Xp, Pc[t], 0), r * sizeof(ulong));
                    for (u = 0; u < t; u++)
                        l[r + u] = x[r + u];
                    bpos[r + t] = bidx[Pc[t]];
                }
            }
        }

        if (kk > 0 && !lumode)
        {
            nmod_mat_window_init(Y, Xp, 0, r, kk, D);
            nmod_mat_init(Wq, r, kk, p);

            if (r > 0)
            {
                nmod_mat_window_init(Wf, W, 0, r, r, D);
                _gather_cols(Wq, Wf, Q, r, kk);
            }

            if (G != NULL)
            {
                /* coefficients of the new rows Y: on the chunk rows,
                   H = Xp[0:kk, D:D+k], restricted to the independent rows
                   L[0..kk), and on the old basis rows, -H Xpiv G */
                nmod_mat_t Hc, Hv, Yg, Gw, Gw2;
                slong t2;

                nmod_mat_window_init(Gw, G, r, r, r + kk, r + kk);
                for (t = 0; t < kk; t++)
                    for (t2 = 0; t2 < kk; t2++)
                        nmod_mat_entry(Gw, t, t2) = nmod_mat_entry(Xp, t, D + L[t2]);

                if (r > 0)
                {
                    nmod_mat_window_init(Hc, Xp, 0, D, kk, D + k);
                    nmod_mat_init(Hv, kk, r, p);
                    nmod_mat_mul(Hv, Hc, Xpiv);
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

        nmod_mat_window_clear(Xpiv);
    }

    /* ---------------------------------------------------------------- */
    /* full rank: the next row of each active index is dependent        */
    /* ---------------------------------------------------------------- */
    if (r == D)
    {
        slong c = 0;

        for (i = 0; i < m; i++)
        {
            if (state[i] != ACTIVE)
                continue;

            e = nexte[i];
            state[i] = DEAD;
            if (deltas != NULL)
                deltas[i] = e;

            if (T != NULL && D > 0)
            {
                if (e == 0)
                    memcpy(nmod_mat_entry_ptr(T, i, 0),
                           nmod_mat_entry_ptr(V, i, 0), D * sizeof(ulong));
                else if (!dense)
                    _nmod_mat_krylov_op_apply_vec(nmod_mat_entry_ptr(T, i, 0),
                            nmod_mat_entry_ptr(B, last[i], 0), op);
                else
                    heap[c++] = i;  /* reuse heap as a list */
            }
        }

        if (c > 0)
        {
            nmod_mat_init(S, c, D, p);
            nmod_mat_init(O, c, D, p);
            for (u = 0; u < c; u++)
                memcpy(nmod_mat_entry_ptr(S, u, 0),
                       nmod_mat_entry_ptr(B, last[heap[u]], 0), D * sizeof(ulong));
            nmod_mat_mul(O, S, op->M);
            for (u = 0; u < c; u++)
                memcpy(nmod_mat_entry_ptr(T, heap[u], 0),
                       nmod_mat_entry_ptr(O, u, 0), D * sizeof(ulong));
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
    nmod_mat_clear(X);
    nmod_mat_clear(Xp);
    nmod_mat_clear(W);
    flint_free(beta);
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
        const nmod_mat_t V, const _nmod_mat_krylov_op_t op,
        const slong * shift, const slong * degrees, slong chunk)
{
    return _krylov_basis_iter(B, rows, exps, deltas, T, G, 0, NULL, NULL, piv, ok,
                              V, op, shift, degrees, chunk);
}

slong
_nmod_mat_krylov_basis_iter_lu(nmod_mat_t B, slong * rows, slong * exps,
        slong * deltas, nmod_mat_t T, nmod_mat_t LU, slong * bpos, slong * piv,
        int * ok, const nmod_mat_t V, const _nmod_mat_krylov_op_t op,
        const slong * shift, const slong * degrees, slong chunk)
{
    return _krylov_basis_iter(B, rows, exps, deltas, T, NULL, 1, LU, bpos, piv, ok,
                              V, op, shift, degrees, chunk);
}

slong
_nmod_mat_krylov_basis_naive(nmod_mat_t B, slong * rows, slong * exps,
        const nmod_mat_t V, const _nmod_mat_krylov_op_t op,
        const slong * shift, const slong * degrees)
{
    slong m = V->r, D = V->c, N, a, r, t;
    slong * beta, * ri, * re;
    char * indep;
    nmod_mat_t K, X;

    beta = flint_malloc((m + 1) * sizeof(slong));
    _nmod_mat_krylov_degrees(beta, m, degrees, D, 0);
    N = _nmod_mat_krylov_nrows(m, beta);

    ri = flint_malloc((2 * N + 1) * sizeof(slong));
    re = ri + N;
    indep = flint_malloc(N + 1);

    nmod_mat_init(K, N, D, V->mod.n);
    _nmod_mat_krylov_matrix(K, ri, re, V, op, shift, beta);

    nmod_mat_init_set(X, K);
    r = _nmod_mat_row_rank_profile(indep, X);

    t = 0;
    for (a = 0; a < N; a++)
    {
        if (indep[a])
        {
            memcpy(nmod_mat_entry_ptr(B, t, 0), nmod_mat_entry_ptr(K, a, 0),
                   D * sizeof(ulong));
            rows[t] = ri[a];
            exps[t] = re[a];
            t++;
        }
    }

    nmod_mat_clear(K);
    nmod_mat_clear(X);
    flint_free(beta);
    flint_free(ri);
    flint_free(indep);

    return r;
}

/*
    The basis alone is computed in row echelon form (LU mode, see
    _nmod_mat_krylov_basis_iter_lu) rather than in reduced form: measured
    with uniform and Hermite shifts, Hermite-Pade and dense M, D = 64 to
    1024 (Nov 2026, Zen 4, Ice Lake, Cascade Lake, Apple M4), it takes 0.54
    to 1.03 times the time, except at D = 128 on x86 below 2^52 (0.78 to
    1.19 at 30 bits, 1.01 to 1.06 at 50 bits).
*/
slong
nmod_mat_krylov_basis(nmod_mat_t B, slong * rows, slong * exps,
        const nmod_mat_t V, const nmod_mat_t M,
        const slong * shift, const slong * degrees)
{
    _nmod_mat_krylov_op_t op;
    slong r;

    if (M->r != V->c || M->c != V->c)
        flint_throw(FLINT_ERROR, "nmod_mat_krylov_basis: incompatible dimensions\n");

    _nmod_mat_krylov_op_init_dense(op, M);
    r = _nmod_mat_krylov_basis_iter_lu(B, rows, exps, NULL, NULL, NULL, NULL,
                                       NULL, NULL, V, op, shift, degrees, 0);
    _nmod_mat_krylov_op_clear(op);

    return r;
}

slong
nmod_mat_krylov_basis_jordan(nmod_mat_t B, slong * rows, slong * exps,
        const nmod_mat_t V, nn_srcptr eig, const slong * sizes, slong nblocks,
        const slong * shift, const slong * degrees)
{
    _nmod_mat_krylov_op_t op;
    slong r;

    _nmod_mat_krylov_op_init_jordan(op, eig, sizes, nblocks, V->mod);

    if (op->D != V->c)
        flint_throw(FLINT_ERROR, "nmod_mat_krylov_basis_jordan: incompatible dimensions\n");

    r = _nmod_mat_krylov_basis_iter_lu(B, rows, exps, NULL, NULL, NULL, NULL,
                                       NULL, NULL, V, op, shift, degrees, 0);
    _nmod_mat_krylov_op_clear(op);

    return r;
}
