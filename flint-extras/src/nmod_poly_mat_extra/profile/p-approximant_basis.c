/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

/*
    Krylov in the approximant and interpolant bases: pmbasis and pmintbasis
    with and without Krylov leaves, and the divide and conquer behind
    nmod_poly_mat_approximant_basis and nmod_poly_mat_interpolant_basis with
    several values of its parameters. Times in ms (best of up to 3 runs, one
    run beyond 1 s); each line ends with "ok", or with "MISMATCH" if two
    variants give different shifted row degrees.

        p-approximant_basis mulmid BITS [MMAX]
            middle products (m x m, length len1) times (m x n, window of L
            coefficients) as in the residual of the approximant bases: the
            dispatcher, the products of coefficient matrices, and the other
            routines (naive, sd_fft_direct, sd_fft_matmul, vandermonde1,
            vandermonde2; -1 if not applicable); m = 2 .. MMAX, n = 1 .. 16
            (at most m/2), len1 = 2 .. 64, L = len1, 8 len1, 64 len1
        p-approximant_basis leaf BITS [MMAX] [RMAX]
            pmbasis on random m x n matrices at order d, D = n d, with the
            Krylov cutoff c = 0 (no leaf), 4, 6 (default), 8, 12, 16;
            m = 16 .. MMAX (default 256), n = 1, 2, 4, 8, m/8, m/4, m/2,
            D = m/4 .. RMAX m (default 32, at most 16 for m = 128, 8 for
            m = 256)
        p-approximant_basis ileaf BITS [MMAX] [RMAX]
            the same with pmintbasis, at d distinct points
        p-approximant_basis orders BITS [MMAX] [RMAX]
            nmod_poly_mat_approximant_basis on m x n matrices with an order
            per column (profiles below), against pmbasis on the input raised
            to the largest order and the direct Krylov function (D <= 8m);
            divide and conquer with the defaults (direct cutoff 8,
            quasi-uniform cut, eps = 0.25), then with direct cutoff 4, 16;
            with the split in half of each order; with eps = -1 (never the
            emulation), 0.5, 1
        p-approximant_basis interp BITS [MMAX] [RMAX]
            nmod_poly_mat_interpolant_basis with one point per column, against
            the direct Krylov function (D <= 8m); the default (direct
            cutoff 8), then direct cutoff 4, 16
*/

#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <flint/flint.h>
#include <flint/ulong_extras.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>

#include "nmod_poly_mat_multiply.h"
#include "nmod_poly_mat_approximant.h"
#include "nmod_poly_mat_interpolant.h"
#include "nmod_poly_mat_extra/impl.h"

static double
_now_ms(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec * 1e3 + ts.tv_nsec * 1e-6;
}

/* best of up to 3 runs, one run if it takes more than 1 s; CODE is run with
 * the shift reset to SH0 (m entries) in SH */
#define BEST(res, sh, sh0, m, code)                                   \
    do {                                                              \
        double _b = 1e300, _tot = 0;                                  \
        for (int _r = 0; _r < 3 && _tot < 1000; _r++)                 \
        {                                                             \
            memcpy(sh, sh0, (m) * sizeof(slong));                     \
            double _t0 = _now_ms();                                   \
            code;                                                     \
            double _t = _now_ms() - _t0;                              \
            _tot += _t;                                               \
            if (_t < _b) _b = _t;                                     \
        }                                                             \
        res = _b;                                                     \
    } while (0)

static int
_same(const slong * a, const slong * b, slong m)
{
    for (slong i = 0; i < m; i++)
        if (a[i] != b[i])
            return 0;
    return 1;
}

static slong
_rmax_for(slong m, slong rmax)
{
    if (m >= 256)
        return FLINT_MIN(rmax, 8);
    if (m >= 128)
        return FLINT_MIN(rmax, 16);
    return rmax;
}

/*------------------------------------------------------------*/
/* mulmid                                                     */
/*------------------------------------------------------------*/

/* best of up to 3 runs of CODE, one run if it takes more than 0.5 s */
#define BEST0(res, code)                                              \
    do {                                                              \
        double _b = 1e300, _tot = 0;                                  \
        for (int _r = 0; _r < 3 && _tot < 500; _r++)                  \
        {                                                             \
            double _t0 = _now_ms();                                   \
            code;                                                     \
            double _t = _now_ms() - _t0;                              \
            _tot += _t;                                               \
            if (_t < _b) _b = _t;                                     \
        }                                                             \
        res = _b;                                                     \
    } while (0)

static void
profile_mulmid(ulong p, slong mmax, flint_rand_t state)
{
    const slong ms[5] = {2, 4, 16, 64, 256};
    const slong ns[5] = {1, 2, 4, 8, 16};
    const slong lens[6] = {2, 4, 8, 16, 32, 64};
    const slong facs[3] = {1, 8, 64};

    flint_printf("# mulmid; m n len1 L | ms: dispatch coeffs | naive sd_fft_direct sd_fft_matmul vandermonde1 vandermonde2\n");

    for (int a = 0; a < 5 && ms[a] <= mmax; a++)
    for (int b = 0; b < 5; b++)
    for (int c = 0; c < 6; c++)
    for (int e = 0; e < 3; e++)
    {
        const slong m = ms[a], n = ns[b], len1 = lens[c], L = len1 * facs[e];
        const slong nlo = len1 - 1 + L / 2, nhi = nlo + L, len2 = nhi;
        const slong npts = len1 + L - 1;
        double t[7];
        nmod_poly_mat_t A, B, R, R2;
        int ok = 1;

        if (n > FLINT_MAX(1, m / 2) || L > 4096
            || (double) m * m * n * L * len1 > 2e10)
            continue;

        nmod_poly_mat_init(A, m, m, p);
        nmod_poly_mat_init(B, m, n, p);
        nmod_poly_mat_init(R, m, n, p);
        nmod_poly_mat_init(R2, m, n, p);
        nmod_poly_mat_rand(A, state, len1);
        nmod_poly_mat_rand(B, state, len2);

        BEST0(t[0], nmod_poly_mat_mulmid(R, A, B, nlo, nhi));
        BEST0(t[1], _nmod_poly_mat_mulmid_coeffs(R2, A, len1, B, len2, nlo, nhi));
        ok = ok && nmod_poly_mat_equal(R, R2);
        BEST0(t[2], _nmod_poly_mat_mulmid_naive(R2, A, len1, B, len2, nlo, nhi));
        ok = ok && nmod_poly_mat_equal(R, R2);
        BEST0(t[3], _nmod_poly_mat_mulmid_sd_fft_direct(R2, A, len1, B, len2, nlo, nhi));
        ok = ok && nmod_poly_mat_equal(R, R2);
        BEST0(t[4], _nmod_poly_mat_mulmid_sd_fft_matmul(R2, A, len1, B, len2, nlo, nhi));
        ok = ok && nmod_poly_mat_equal(R, R2);
        t[5] = t[6] = -1;
        if (NMOD_POLY_CAN_USE_VANDERMONDE1(p, npts))
        {
            BEST0(t[5], _nmod_poly_mat_mulmid_vandermonde1(R2, A, len1, B, len2, nlo, nhi));
            ok = ok && nmod_poly_mat_equal(R, R2);
        }
        if (NMOD_POLY_CAN_USE_VANDERMONDE2(p, npts))
        {
            BEST0(t[6], _nmod_poly_mat_mulmid_vandermonde2(R2, A, len1, B, len2, nlo, nhi));
            ok = ok && nmod_poly_mat_equal(R, R2);
        }

        flint_printf("mulmid %wd %wd %wd %wd | %.4f %.4f | %.4f %.4f %.4f %.4f %.4f %s\n",
                     m, n, len1, L, t[0], t[1], t[2], t[3], t[4], t[5], t[6],
                     ok ? "ok" : "MISMATCH");
        fflush(stdout);

        nmod_poly_mat_clear(A);
        nmod_poly_mat_clear(B);
        nmod_poly_mat_clear(R);
        nmod_poly_mat_clear(R2);
    }
}

/*------------------------------------------------------------*/
/* leaf, ileaf                                                */
/*------------------------------------------------------------*/

static void
profile_leaf(int interp, ulong p, slong mmax, slong rmax, flint_rand_t state)
{
    const double cutoffs[6] = {0.0, 4.0, 6.0, 8.0, 12.0, 16.0};
    const double fr[8] = {0.25, 0.5, 1, 2, 4, 8, 16, 32};

    flint_printf("# %s; m n d D/m | ms: cutoff 0 4 6 8 12 16 | cutoff0/cutoff6\n",
                 interp ? "pmintbasis (ileaf)" : "pmbasis (leaf)");

    for (slong m = 16; m <= mmax; m *= 2)
    {
        slong ns[7] = {1, 2, 4, 8, m / 8, m / 4, m / 2};
        for (int a = 0; a < 7; a++)
        {
            const slong n = ns[a];
            int dup = 0;
            for (int b = 0; b < a; b++)
                if (ns[b] == n)
                    dup = 1;
            if (dup || n < 1 || n >= m)
                continue;

            for (int c = 0; c < 8 && fr[c] <= _rmax_for(m, rmax); c++)
            {
                const slong d = (slong) (fr[c] * m / n + 0.5);
                double t[6];
                slong * sh0, * sh, * ref;
                int ok = 1;

                if (d < 1 || (fr[c] < 1 && n * d > 1.5 * fr[c] * m))
                    continue;
                if (interp && (ulong) d >= p)
                    continue;

                sh0 = flint_calloc(m, sizeof(slong));
                sh = flint_calloc(m, sizeof(slong));
                ref = flint_calloc(m, sizeof(slong));

                if (!interp)
                {
                    nmod_poly_mat_t F, P;
                    nmod_poly_mat_init(F, m, n, p);
                    nmod_poly_mat_init(P, m, m, p);
                    nmod_poly_mat_rand(F, state, d);
                    for (int k = 0; k < 6; k++)
                    {
                        BEST(t[k], sh, sh0, m, _nmod_poly_mat_pmbasis(P, sh, F, d, cutoffs[k]));
                        if (k == 0)
                            memcpy(ref, sh, m * sizeof(slong));
                        else if (!_same(ref, sh, m))
                            ok = 0;
                    }
                    nmod_poly_mat_clear(F);
                    nmod_poly_mat_clear(P);
                }
                else
                {
                    nmod_mat_struct * E = flint_malloc(d * sizeof(nmod_mat_struct));
                    ulong * pts = flint_malloc(d * sizeof(ulong));
                    nmod_poly_mat_t P;
                    nmod_poly_mat_init(P, m, m, p);
                    for (slong k = 0; k < d; k++)
                    {
                        nmod_mat_init(E + k, m, n, p);
                        nmod_mat_rand(E + k, state);
                        pts[k] = k + 1;  /* distinct since d < p */
                    }
                    for (int k = 0; k < 6; k++)
                    {
                        BEST(t[k], sh, sh0, m, _nmod_poly_mat_pmintbasis(P, sh, pts, E, d, cutoffs[k]));
                        if (k == 0)
                            memcpy(ref, sh, m * sizeof(slong));
                        else if (!_same(ref, sh, m))
                            ok = 0;
                    }
                    for (slong k = 0; k < d; k++)
                        nmod_mat_clear(E + k);
                    flint_free(E);
                    flint_free(pts);
                    nmod_poly_mat_clear(P);
                }

                flint_printf("%s %wd %wd %wd %g | %.3f %.3f %.3f %.3f %.3f %.3f | %.2f %s\n",
                             interp ? "ileaf" : "leaf", m, n, d, (double) n * d / m,
                             t[0], t[1], t[2], t[3], t[4], t[5], t[0] / t[2], ok ? "ok" : "MISMATCH");
                fflush(stdout);
                flint_free(sh0);
                flint_free(sh);
                flint_free(ref);
            }
        }
    }
}

/*------------------------------------------------------------*/
/* orders                                                     */
/*------------------------------------------------------------*/

#define NPROF 6
static const char * prof_name[NPROF] = {"geo", "one", "hp", "steps", "close", "wide"};

/* orders of profile k for about D = target conditions; returns n */
static slong
_profile_orders(slong * orders, int k, slong m, double target, flint_rand_t state)
{
    double w[64];
    slong n = 0;
    double sw = 0;

    if (k == 0)        /* geometric: d, d/4, d/16, d/64 */
    {
        n = 4; w[0] = 1; w[1] = 0.25; w[2] = 1.0 / 16; w[3] = 1.0 / 64;
    }
    else if (k == 1)   /* one column at d, seven at d/8 */
    {
        n = 8; w[0] = 1;
        for (slong b = 1; b < n; b++) w[b] = 0.125;
    }
    else if (k == 2)   /* one column at d, fifteen at d/32 */
    {
        n = 16; w[0] = 1;
        for (slong b = 1; b < n; b++) w[b] = 1.0 / 32;
    }
    else if (k == 3)   /* d, d, d/2, d/2 */
    {
        n = 4; w[0] = w[1] = 1; w[2] = w[3] = 0.5;
    }
    else if (k == 4)   /* close: 8 columns in [0.85 d, d] */
    {
        n = 8;
        for (slong b = 0; b < n; b++) w[b] = 0.85 + 0.15 * b / (n - 1);
    }
    else               /* wide: m/2 columns in [d/4, d] (many columns) */
    {
        n = FLINT_MIN(m / 2, 64);
        for (slong b = 0; b < n; b++) w[b] = 0.25 + 0.75 * b / (n - 1);
    }

    for (slong b = 0; b < n; b++)
        sw += w[b];
    for (slong b = 0; b < n; b++)
        orders[b] = FLINT_MAX(1, (slong) (w[b] * target / sw + 0.5));
    (void) state;
    return n;
}

static void
profile_orders(ulong p, slong mmax, slong rmax, flint_rand_t state)
{
    const double fr[6] = {1, 2, 4, 8, 16, 32};
    /* variants of the divide and conquer: direct cutoff, split, eps (the
     * default first: NMOD_POLY_MAT_APPROXIMANT_BASIS_KRYLOV_CUTOFF,
     * quasi-uniform cut, NMOD_POLY_MAT_APPROXIMANT_BASIS_RAISE_SLACK) */
    const double va[7] = {-1, 4, 16, -1, -1, -1, -1};
    const int vs[7] = {0, 0, 0, 1, 0, 0, 0};
    const double ve[7] = {-2, -2, -2, -2, -1.0, 0.5, 1.0};

    flint_printf("# orders; prof m n dmax D/m D Demul | ms: emul krylov | dac default, direct cutoff 4 16, half, eps -1 0.5 1 | emul/default\n");

    for (int k = 0; k < NPROF; k++)
    for (slong m = 16; m <= mmax; m *= 4)
    for (int c = 0; c < 6 && fr[c] <= _rmax_for(m, rmax); c++)
    {
        slong orders[64], n, D = 0, dmax = 0;
        double te, tk = -1, t[7];
        slong * sh0, * sh, * ref;
        nmod_poly_mat_t F, Fe, P;
        _nmod_poly_mat_basis_dac_param_t prm;
        int ok = 1;

        n = _profile_orders(orders, k, m, fr[c] * m, state);
        if (k == 5 && n < 2)
            continue;
        for (slong b = 0; b < n; b++)
        {
            D += orders[b];
            dmax = FLINT_MAX(dmax, orders[b]);
        }

        sh0 = flint_calloc(m, sizeof(slong));
        sh = flint_calloc(m, sizeof(slong));
        ref = flint_calloc(m, sizeof(slong));
        nmod_poly_mat_init(F, m, n, p);
        nmod_poly_mat_init(Fe, m, n, p);
        nmod_poly_mat_init(P, m, m, p);
        for (slong i = 0; i < m; i++)
            for (slong b = 0; b < n; b++)
            {
                nmod_poly_struct * f = nmod_poly_mat_entry(F, i, b);
                nmod_poly_fit_length(f, orders[b]);
                for (slong l = 0; l < orders[b]; l++)
                    f->coeffs[l] = n_randint(state, p);
                _nmod_poly_set_length(f, orders[b]);
                _nmod_poly_normalise(f);
                nmod_poly_shift_left(nmod_poly_mat_entry(Fe, i, b), f, dmax - orders[b]);
            }

        BEST(te, sh, sh0, m, nmod_poly_mat_pmbasis(P, sh, Fe, dmax));
        memcpy(ref, sh, m * sizeof(slong));
        if (D <= 8 * m)
        {
            BEST(tk, sh, sh0, m, _nmod_poly_mat_approximant_basis_krylov_leaf(P, sh, F, orders));
            if (!_same(ref, sh, m))
                ok = 0;
        }
        for (int v = 0; v < 7; v++)
        {
            _nmod_poly_mat_basis_dac_param_init(prm);
            if (va[v] > 0)
                prm->direct_cutoff = va[v];
            prm->split = vs[v];
            if (ve[v] > -2)
                prm->raise_slack = ve[v];
            if (v == 0)  /* warm-up run, not timed */
            {
                memcpy(sh, sh0, m * sizeof(slong));
                _nmod_poly_mat_approximant_basis(P, sh, F, orders, prm);
            }
            BEST(t[v], sh, sh0, m, _nmod_poly_mat_approximant_basis(P, sh, F, orders, prm));
            if (!_same(ref, sh, m))
                ok = 0;
        }

        flint_printf("orders %s %wd %wd %wd %g %wd %wd | %.3f %.3f | %.3f %.3f %.3f %.3f %.3f %.3f %.3f | %.2f %s\n",
                     prof_name[k], m, n, dmax, fr[c], D, n * dmax, te, tk,
                     t[0], t[1], t[2], t[3], t[4], t[5], t[6], te / t[0], ok ? "ok" : "MISMATCH");
        fflush(stdout);

        nmod_poly_mat_clear(F);
        nmod_poly_mat_clear(Fe);
        nmod_poly_mat_clear(P);
        flint_free(sh0);
        flint_free(sh);
        flint_free(ref);
    }
}

/*------------------------------------------------------------*/
/* interp                                                     */
/*------------------------------------------------------------*/

#define NIPROF 4
static const char * iprof_name[NIPROF] = {"twice", "pool", "hp", "pairs"};

/* points of profile k for D columns, sorted: twice (D/2 points, each twice,
 * except one once), pool (D/4 points, random multiplicities), hp (one point
 * with D/2 columns, the others distinct), pairs (D/2 points twice, D/4 once,
 * mixed) */
static void
_profile_points(ulong * pts, int k, slong D, ulong p, flint_rand_t state)
{
    slong b = 0;
    if (k == 0)
    {
        ulong x = 1;
        for (b = 0; b + 1 < D; b += 2, x++)
            pts[b] = pts[b + 1] = x;
        if (b < D)
            pts[b] = x;
    }
    else if (k == 1)
    {
        const slong np = FLINT_MAX(1, D / 4);
        for (b = 0; b < D; b++)
            pts[b] = 1 + n_randint(state, np);
        /* sort */
        for (slong i = 1; i < D; i++)
        {
            ulong x = pts[i];
            slong j = i;
            while (j > 0 && pts[j - 1] > x)
            {
                pts[j] = pts[j - 1];
                j--;
            }
            pts[j] = x;
        }
    }
    else if (k == 2)
    {
        for (b = 0; b < D / 2; b++)
            pts[b] = 1;
        for (ulong x = 2; b < D; b++, x++)
            pts[b] = x;
    }
    else
    {
        ulong x = 1;
        while (b < D)
        {
            const slong c = (x % 3 == 0) ? 1 : 2;
            for (slong l = 0; l < c && b < D; l++)
                pts[b++] = x;
            x++;
        }
    }
    (void) p;
}

static void
profile_interp(ulong p, slong mmax, slong rmax, flint_rand_t state)
{
    const double fr[6] = {1, 2, 4, 8, 16, 32};
    /* variants: the direct cutoff (the default first) */
    const double va[3] = {-1, 4, 16};

    flint_printf("# interp; prof m D/m D | ms: krylov | dac default, direct cutoff 4 16\n");

    for (int k = 0; k < NIPROF; k++)
    for (slong m = 16; m <= mmax; m *= 4)
    for (int c = 0; c < 6 && fr[c] <= _rmax_for(m, rmax); c++)
    {
        const slong D = (slong) (fr[c] * m);
        double tk = -1, t[3];
        slong * sh0, * sh, * ref;
        ulong * pts;
        nmod_mat_t V;
        nmod_poly_mat_t P;
        _nmod_poly_mat_basis_dac_param_t prm;
        int ok = 1;

        if ((ulong) D + 2 >= p)
            continue;

        sh0 = flint_calloc(m, sizeof(slong));
        sh = flint_calloc(m, sizeof(slong));
        ref = flint_calloc(m, sizeof(slong));
        pts = flint_malloc(D * sizeof(ulong));
        nmod_mat_init(V, m, D, p);
        nmod_poly_mat_init(P, m, m, p);
        nmod_mat_rand(V, state);
        _profile_points(pts, k, D, p, state);

        for (int v = 0; v < 3; v++)
        {
            _nmod_poly_mat_basis_dac_param_init(prm);
            if (va[v] > 0)
                prm->direct_cutoff = va[v];
            if (v == 0)  /* warm-up run, not timed */
            {
                memcpy(sh, sh0, m * sizeof(slong));
                _nmod_poly_mat_interpolant_basis(P, sh, V, pts, prm);
            }
            BEST(t[v], sh, sh0, m, _nmod_poly_mat_interpolant_basis(P, sh, V, pts, prm));
            if (v == 0)
                memcpy(ref, sh, m * sizeof(slong));
            else if (!_same(ref, sh, m))
                ok = 0;
        }
        if (D <= 8 * m)
        {
            BEST(tk, sh, sh0, m, _nmod_poly_mat_interpolant_basis_krylov_leaf(P, sh, V, pts));
            if (!_same(ref, sh, m))
                ok = 0;
        }

        flint_printf("interp %s %wd %g %wd | %.3f | %.3f %.3f %.3f %s\n",
                     iprof_name[k], m, fr[c], D, tk, t[0], t[1], t[2],
                     ok ? "ok" : "MISMATCH");
        fflush(stdout);

        nmod_mat_clear(V);
        nmod_poly_mat_clear(P);
        flint_free(pts);
        flint_free(sh0);
        flint_free(sh);
        flint_free(ref);
    }
}

int main(int argc, char ** argv)
{
    flint_rand_t state;
    const char * mode;
    slong bits, mmax, rmax;
    ulong p;

    if (argc < 3)
    {
        flint_printf("usage: %s mulmid|leaf|ileaf|orders|interp BITS [MMAX] [RMAX]\n", argv[0]);
        return 1;
    }
    mode = argv[1];
    bits = atol(argv[2]);
    mmax = (argc > 3) ? atol(argv[3]) : 256;
    rmax = (argc > 4) ? atol(argv[4]) : 32;

    flint_rand_init(state);
    p = n_nextprime(UWORD(1) << (bits - 1), 1);
    flint_printf("# %s, %wd-bit prime %wu, MMAX %wd, RMAX %wd\n", mode, bits, p, mmax, rmax);

    if (!strcmp(mode, "mulmid"))
        profile_mulmid(p, mmax, state);
    else if (!strcmp(mode, "leaf"))
        profile_leaf(0, p, mmax, rmax, state);
    else if (!strcmp(mode, "ileaf"))
        profile_leaf(1, p, mmax, rmax, state);
    else if (!strcmp(mode, "orders"))
        profile_orders(p, mmax, rmax, state);
    else if (!strcmp(mode, "interp"))
        profile_interp(p, mmax, rmax, state);
    else
    {
        flint_printf("unknown mode %s\n", mode);
        return 1;
    }

    flint_rand_clear(state);
    return 0;
}
