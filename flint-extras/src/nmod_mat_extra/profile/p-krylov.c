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
    Timings of the Krylov functions.

    Usage: p-krylov [mode] [bits]

      mode = dense   (default) multiplication matrix dense random: Krylov
                     basis (iterative, doubling), compact kernel, with
                     uniform and "Hermite" shifts, for 1, D/8 and D rows
                     (V is m x D, M is D x D)
      mode = hp      Hermite-Pade: approximant basis of an m x n matrix at
                     order d in each column (Jordan matrix with n nilpotent
                     blocks of size d, D = n d), in shifted Popov form
      mode = interp  interpolation at D points (diagonal matrix)
      mode = chunk   effect of the chunk size (Hermite-Pade, m = D)
      mode = kernel  the two methods for the relations of the compact kernel
                     (reduced basis and G, or row echelon basis and
                     triangular solves) for m / D from 1/16 to 4,
                     Hermite-Pade (one nilpotent block) and dense M
      mode = pow     dense random M, generation of the rows from powers
                     M^(2^t) (Keller-Gehrig-type doubling), D up to Dmax:
                     p-krylov pow [bits] [Dmax]
      mode = budget  budget of the indices in a chunk (op->predict): from
                     the generic profile (default) against max(1, f_i)
                     alone, uniform, staggered, random and Hermite shifts,
                     D up to Dmax: p-krylov budget [bits] [Dmax]

    Times in seconds; "mul" is the time of nmod_mat_mul for D x D
    matrices.
*/

#include <stdlib.h>
#include <string.h>
#include <flint/profiler.h>
#include <flint/ulong_extras.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_poly_mat.h>
#include "nmod_mat_extra/krylov_impl.h"
#include "nmod_poly_mat_approximant.h"
#include "nmod_poly_mat_interpolant.h"

#define TIME_IT(res, code)                                  \
    do {                                                    \
        timeit_t _t;                                        \
        slong _reps = 1;                                    \
        while (1)                                           \
        {                                                   \
            slong _i;                                       \
            timeit_start(_t);                               \
            for (_i = 0; _i < _reps; _i++) { code; }        \
            timeit_stop(_t);                                \
            if (_t->wall >= 100) break;                     \
            _reps *= (_t->wall < 10) ? 10 : 2;              \
        }                                                   \
        (res) = 0.001 * _t->wall / _reps;                   \
    } while (0)

static double
_time_mul(slong D, ulong p, flint_rand_t state)
{
    nmod_mat_t A, B, C;
    double t;

    nmod_mat_init(A, D, D, p);
    nmod_mat_init(B, D, D, p);
    nmod_mat_init(C, D, D, p);
    nmod_mat_rand(A, state);
    nmod_mat_rand(B, state);
    TIME_IT(t, nmod_mat_mul(C, A, B));
    nmod_mat_clear(A);
    nmod_mat_clear(B);
    nmod_mat_clear(C);

    return t;
}

/* time of c x D times D x D */
static double
_time_mul_rows(slong c, slong D, ulong p, flint_rand_t state)
{
    nmod_mat_t A, B, C;
    double t;

    nmod_mat_init(A, c, D, p);
    nmod_mat_init(B, D, D, p);
    nmod_mat_init(C, c, D, p);
    nmod_mat_rand(A, state);
    nmod_mat_rand(B, state);
    TIME_IT(t, nmod_mat_mul(C, A, B));
    nmod_mat_clear(A);
    nmod_mat_clear(B);
    nmod_mat_clear(C);

    return t;
}

static void
profile_dense(flint_bitcnt_t bits, flint_rand_t state)
{
    slong Ds[] = {50, 100, 200, 400, 800, 0};
    slong k, j, s;
    ulong p = n_randprime(state, bits, 1);

    flint_printf("dense random M, %wu-bit prime\n", bits);
    flint_printf("    D     m  shift     mul    basis  doubling  kernel   (basis/mul  kernel/mul  doubling/basis)\n");

    for (k = 0; Ds[k] != 0; k++)
    {
        slong D = Ds[k];
        slong ms[3];
        double tmul = _time_mul(D, p, state);

        ms[0] = 1;
        ms[1] = FLINT_MAX(1, D / 8);
        ms[2] = D;

        for (j = 0; j < 3; j++)
        {
            for (s = 0; s < 2; s++)
            {
                slong m = ms[j], i, r = 0;
                slong * shift = flint_malloc(m * sizeof(slong));
                slong * rows = flint_malloc(D * sizeof(slong));
                slong * exps = flint_malloc(D * sizeof(slong));
                slong * deltas = flint_malloc(m * sizeof(slong));
                nmod_mat_t V, M, B, R;
                _nmod_mat_krylov_op_t op;
                double tb, td, tk;
                int ok;

                if (s == 1 && m == 1)
                {
                    flint_free(shift); flint_free(rows); flint_free(exps); flint_free(deltas);
                    continue;
                }

                for (i = 0; i < m; i++)
                    shift[i] = (s == 0) ? 0 : i * (D + 1);

                nmod_mat_init(V, m, D, p);
                nmod_mat_init(M, D, D, p);
                nmod_mat_init(B, D, D, p);
                nmod_mat_init(R, m, D, p);
                nmod_mat_rand(V, state);
                nmod_mat_rand(M, state);
                _nmod_mat_krylov_op_init_dense(op, M);

                TIME_IT(tb, r = _nmod_mat_krylov_basis_iter(B, rows, exps, NULL,
                            NULL, NULL, NULL, &ok, V, op, shift, NULL, 0));
                TIME_IT(td, _nmod_mat_krylov_basis_doubling(B, rows, exps, NULL,
                            &ok, V, op, shift, NULL));
                TIME_IT(tk, _nmod_mat_krylov_kernel(R, rows, exps, deltas, V, op,
                            shift, NULL, 0, _NMOD_MAT_KRYLOV_KERNEL_AUTO));

                flint_printf("%5wd %5wd  %s  %.2e %.2e %.2e %.2e   (%6.2f %6.2f %6.2f)   rank %wd\n",
                    D, m, s == 0 ? "unif" : "herm", tmul, tb, td, tk,
                    tb / tmul, tk / tmul, td / tb, r);

                _nmod_mat_krylov_op_clear(op);
                nmod_mat_clear(V);
                nmod_mat_clear(M);
                nmod_mat_clear(B);
                nmod_mat_clear(R);
                flint_free(shift);
                flint_free(rows);
                flint_free(exps);
                flint_free(deltas);
            }
        }
    }
}

/*
    Generation of the rows from powers M^(2^t) (op->logpow), dense random
    M: Hermite shift (the iterates of each row of V in turn), or uniform
    shift with few rows. Times relative to logpow = 0; "KG" uses chunks of
    up to m + D rows and all the powers needed (Keller-Gehrig doubling with
    incremental elimination); "auto" is the default (automatic choice);
    "paper" is _nmod_mat_krylov_basis_doubling.
*/
static void
profile_pow(flint_bitcnt_t bits, slong Dmax, flint_rand_t state)
{
    slong Ds[] = {200, 400, 800, 1600, 3200, 6400, 0};
    slong cfg_m[] = {1, 4, 2, 4, 16};
    int cfg_herm[] = {0, 1, 0, 0, 0};
    slong lps[] = {1, 2, 3, 4, 5, 7};
    slong nlp = 6, k, c, l, i;
    ulong p = n_randprime(state, bits, 1);

    flint_printf("powers of M for generating rows, dense random M, %wu-bit prime\n", bits);
    flint_printf("rho1, rho8: time of D rows as products of 1 or 8 rows by M, relative to mul\n");
    flint_printf("    D     m shift      mul  logpow=0  | ratio for logpow =");
    for (l = 0; l < nlp; l++)
        flint_printf(" %4wd", lps[l]);
    flint_printf("    KG  auto  paper\n");

    for (k = 0; Ds[k] != 0 && Ds[k] <= Dmax; k++)
    {
        slong D = Ds[k];
        double tmul = _time_mul(D, p, state);
        double rho1 = D * _time_mul_rows(1, D, p, state) / tmul;
        double rho8 = (D / 8) * _time_mul_rows(8, D, p, state) / tmul;

        flint_printf("D = %wd: rho1 = %.2f, rho8 = %.2f\n", D, rho1, rho8);

        for (c = 0; c < 5; c++)
        {
            slong m = cfg_m[c], r0, r;
            slong * shift = flint_malloc(m * sizeof(slong));
            slong * rows0 = flint_malloc(D * sizeof(slong));
            slong * exps0 = flint_malloc(D * sizeof(slong));
            slong * rows = flint_malloc(D * sizeof(slong));
            slong * exps = flint_malloc(D * sizeof(slong));
            nmod_mat_t V, M, B;
            _nmod_mat_krylov_op_t op;
            double t0, t;
            int ok, bad = 0;

            for (i = 0; i < m; i++)
                shift[i] = cfg_herm[c] ? i * (D + 1) : 0;

            nmod_mat_init(V, m, D, p);
            nmod_mat_init(M, D, D, p);
            nmod_mat_init(B, D, D, p);
            nmod_mat_rand(V, state);
            nmod_mat_rand(M, state);
            _nmod_mat_krylov_op_init_dense(op, M);

            op->logpow = 0;
            TIME_IT(t0, r0 = _nmod_mat_krylov_basis_iter(B, rows0, exps0, NULL,
                        NULL, NULL, NULL, &ok, V, op, shift, NULL, 0));

            flint_printf("%5wd %5wd  %s  %.2e  %.2e  |                  ",
                    D, m, cfg_herm[c] ? "herm" : "unif", tmul, t0);

            for (l = 0; l <= nlp + 1; l++)
            {
                slong chunk = (l == nlp) ? m + D : 0;

                op->logpow = (l < nlp) ? lps[l] : (l == nlp) ? FLINT_BITS : -1;
                TIME_IT(t, r = _nmod_mat_krylov_basis_iter(B, rows, exps, NULL,
                            NULL, NULL, NULL, &ok, V, op, shift, NULL, chunk));
                flint_printf(" %5.2f", t / t0);
                fflush(stdout);

                if (r != r0)
                    bad = 1;
                for (i = 0; i < r && !bad; i++)
                    if (rows[i] != rows0[i] || exps[i] != exps0[i])
                        bad = 1;
            }

            if (D <= 1600)
            {
                TIME_IT(t, r = _nmod_mat_krylov_basis_doubling(B, rows, exps, NULL,
                            &ok, V, op, shift, NULL));
                flint_printf("  %5.2f", t / t0);
                if (r != r0)
                    bad = 1;
            }
            else
                flint_printf("      -");

            flint_printf("%s\n", bad ? "  MISMATCH" : "");

            _nmod_mat_krylov_op_clear(op);
            nmod_mat_clear(V);
            nmod_mat_clear(M);
            nmod_mat_clear(B);
            flint_free(shift);
            flint_free(rows0);
            flint_free(exps0);
            flint_free(rows);
            flint_free(exps);
        }
    }
}

/*
    Methods for the relations (_NMOD_MAT_KRYLOV_KERNEL_TRACK: reduced row
    echelon basis and transformation G; _LU: row echelon basis and
    triangular solves) and the automatic choice, for m / D from 1/16 to 4:
    Hermite-Pade (one nilpotent block of size D, uniform shift) and dense
    random M. The methods are run in turn, and the minimum of three runs of
    each is kept.
*/
static void
profile_kernel(flint_bitcnt_t bits, flint_rand_t state)
{
    slong Ds[] = {64, 256, 1024, 0};
    slong num[] = {1, 1, 1, 1, 1, 2, 4, 0};
    slong den[] = {16, 8, 4, 2, 1, 1, 1, 0};
    slong k, c, kind;
    ulong p = n_randprime(state, bits, 1);

    flint_printf("compact kernel: reduced (G) vs row echelon basis (LU), %wu-bit prime\n", bits);
    flint_printf(" kind     D     m        track           lu     lu/track  auto/track\n");

    for (kind = 0; kind < 2; kind++)
    {
        for (k = 0; Ds[k] != 0; k++)
        {
            slong D = Ds[k];

            if (kind == 1 && D > 256)
                continue;

            for (c = 0; num[c] != 0; c++)
            {
                slong m = FLINT_MAX(1, (D * num[c]) / den[c]);
                slong * rows = flint_malloc(D * sizeof(slong));
                slong * exps = flint_malloc(D * sizeof(slong));
                slong * deltas = flint_malloc(m * sizeof(slong));
                nmod_mat_t V, M, R;
                _nmod_mat_krylov_op_t op;
                double t1, t2, t0;

                nmod_mat_init(V, m, D, p);
                nmod_mat_init(M, D, D, p);
                nmod_mat_init(R, m, D, p);
                nmod_mat_rand(V, state);
                if (kind == 0)
                {
                    slong sz = D;
                    _nmod_mat_krylov_op_init_jordan(op, NULL, &sz, 1, V->mod);
                }
                else
                {
                    nmod_mat_rand(M, state);
                    _nmod_mat_krylov_op_init_dense(op, M);
                }

                {
                    double u0, u1, u2;
                    slong rep;
                    t0 = t1 = t2 = 1e300;
                    for (rep = 0; rep < 3; rep++)
                    {
                        TIME_IT(u1, _nmod_mat_krylov_kernel(R, rows, exps, deltas, V, op,
                                    NULL, NULL, 0, _NMOD_MAT_KRYLOV_KERNEL_TRACK));
                        TIME_IT(u2, _nmod_mat_krylov_kernel(R, rows, exps, deltas, V, op,
                                    NULL, NULL, 0, _NMOD_MAT_KRYLOV_KERNEL_LU));
                        TIME_IT(u0, _nmod_mat_krylov_kernel(R, rows, exps, deltas, V, op,
                                    NULL, NULL, 0, _NMOD_MAT_KRYLOV_KERNEL_AUTO));
                        t0 = FLINT_MIN(t0, u0);
                        t1 = FLINT_MIN(t1, u1);
                        t2 = FLINT_MIN(t2, u2);
                    }
                }

                flint_printf("%5s %5wd %5wd  %.3e  %.3e  %11.2f  %10.2f\n",
                    kind == 0 ? "hp" : "dense", D, m, t1, t2, t2 / t1, t0 / t1);

                _nmod_mat_krylov_op_clear(op);
                nmod_mat_clear(V);
                nmod_mat_clear(M);
                nmod_mat_clear(R);
                flint_free(rows);
                flint_free(exps);
                flint_free(deltas);
            }
        }
    }
}

/*
    Budget of the indices in a chunk: from the generic profile until an
    index stops early (op->predict = 1, the default) against max(1, f_i)
    alone (op->predict = 0), for the compact kernel (automatic method) and
    the basis (row echelon form). Shifts: uniform, staggered s_i = i
    ceil(D / 2m), random in [0, 2D / m], Hermite s_i = i (D + 1); M dense
    random, one nilpotent Jordan block of size D (Hermite-Pade with one
    column), or 4 nilpotent blocks of size D / 4 (four columns, not
    cyclic). Ratios time(predict = 0) / time(predict = 1) (above 1: the
    default is faster), best of 3 interleaved runs; the outputs of both are
    compared.
*/
static void
profile_budget(flint_bitcnt_t bits, slong Dmax, flint_rand_t state)
{
    slong Ds[] = {256, 512, 1024, 2048, 0};
    const char * kinds[] = {"dense", "hp", "hp4"};
    const char * shname[] = {"unif", "stag", "rand", "herm"};
    slong k, kind, c, sh, i;
    ulong p = n_randprime(state, bits, 1);

    flint_printf("budget from the generic profile (default) against max(1, f_i), %wu-bit prime\n", bits);
    flint_printf(" kind     D     m shift   kernel     basis   (time predict=0 / predict=1)\n");

    for (k = 0; Ds[k] != 0 && Ds[k] <= Dmax; k++)
    for (kind = 0; kind < 3; kind++)
    {
        slong D = Ds[k];
        slong ms[4];
        nmod_mat_t M;
        _nmod_mat_krylov_op_t op;

        ms[0] = 1; ms[1] = D / 16; ms[2] = D / 4; ms[3] = D;

        nmod_mat_init(M, D, D, p);
        if (kind == 0)
        {
            nmod_mat_rand(M, state);
            _nmod_mat_krylov_op_init_dense(op, M);
        }
        else
        {
            slong sz[4];
            slong nb = (kind == 1) ? 1 : 4;
            for (i = 0; i < nb; i++)
                sz[i] = D / nb;
            _nmod_mat_krylov_op_init_jordan(op, NULL, sz, nb, M->mod);
        }

        for (c = 0; c < 4; c++)
        for (sh = 0; sh < 4; sh++)
        {
            slong m = ms[c], q, rep, r0 = 0, r1 = 0;
            slong * shift, * rows0, * exps0, * deltas0, * rows1, * exps1, * deltas1;
            nmod_mat_t V, R0, R1, B;
            double tk[2] = {1e300, 1e300}, tb[2] = {1e300, 1e300}, t;
            int bad = 0;

            if (m == 1 && sh != 0)
                continue;

            shift = flint_malloc(m * sizeof(slong));
            deltas0 = flint_malloc(m * sizeof(slong));
            deltas1 = flint_malloc(m * sizeof(slong));
            rows0 = flint_malloc(D * sizeof(slong));
            exps0 = flint_malloc(D * sizeof(slong));
            rows1 = flint_malloc(D * sizeof(slong));
            exps1 = flint_malloc(D * sizeof(slong));
            q = FLINT_MAX(1, (D + 2 * m - 1) / (2 * m));
            for (i = 0; i < m; i++)
                shift[i] = (sh == 0) ? 0 : (sh == 1) ? i * q
                         : (sh == 2) ? (slong) n_randint(state, 2 * D / m + 1) : i * (D + 1);

            nmod_mat_init(V, m, D, p);
            nmod_mat_init(R0, m, D, p);
            nmod_mat_init(R1, m, D, p);
            nmod_mat_init(B, D, D, p);
            nmod_mat_rand(V, state);

            for (rep = 0; rep < 3; rep++)
            {
                op->predict = 0;
                TIME_IT(t, r0 = _nmod_mat_krylov_kernel(R0, rows0, exps0, deltas0, V, op,
                            shift, NULL, 0, _NMOD_MAT_KRYLOV_KERNEL_AUTO));
                tk[0] = FLINT_MIN(tk[0], t);
                op->predict = 1;
                TIME_IT(t, r1 = _nmod_mat_krylov_kernel(R1, rows1, exps1, deltas1, V, op,
                            shift, NULL, 0, _NMOD_MAT_KRYLOV_KERNEL_AUTO));
                tk[1] = FLINT_MIN(tk[1], t);
                op->predict = 0;
                TIME_IT(t, _nmod_mat_krylov_basis_iter_lu(B, rows0, exps0, NULL, NULL,
                            NULL, NULL, NULL, NULL, V, op, shift, NULL, 0));
                tb[0] = FLINT_MIN(tb[0], t);
                op->predict = 1;
                TIME_IT(t, _nmod_mat_krylov_basis_iter_lu(B, rows1, exps1, NULL, NULL,
                            NULL, NULL, NULL, NULL, V, op, shift, NULL, 0));
                tb[1] = FLINT_MIN(tb[1], t);
            }

            if (r0 != r1 || (r0 > 0 && !nmod_mat_equal(R0, R1)))
                bad = 1;
            for (i = 0; i < r0 && !bad; i++)
                if (rows0[i] != rows1[i] || exps0[i] != exps1[i])
                    bad = 1;
            for (i = 0; i < m && !bad; i++)
                if (deltas0[i] != deltas1[i])
                    bad = 1;

            flint_printf("%5s %5wd %5wd  %s  %8.2f  %8.2f%s\n", kinds[kind], D, m,
                    shname[sh], tk[0] / tk[1], tb[0] / tb[1], bad ? "  MISMATCH" : "");
            fflush(stdout);

            nmod_mat_clear(V);
            nmod_mat_clear(R0);
            nmod_mat_clear(R1);
            nmod_mat_clear(B);
            flint_free(shift);
            flint_free(deltas0);
            flint_free(deltas1);
            flint_free(rows0);
            flint_free(exps0);
            flint_free(rows1);
            flint_free(exps1);
        }

        _nmod_mat_krylov_op_clear(op);
        nmod_mat_clear(M);
    }
}

static void
profile_hp(flint_bitcnt_t bits, flint_rand_t state)
{
    /* m, n, d */
    slong cases[][3] = {
        {16, 1, 16}, {32, 1, 32}, {64, 1, 64}, {128, 1, 128}, {256, 1, 256},
        {512, 1, 512}, {64, 1, 256}, {256, 1, 64}, {1024, 1, 256},
        {64, 2, 64}, {64, 4, 32}, {64, 8, 16}, {64, 32, 2}, {64, 64, 1},
        {0, 0, 0}};
    slong k, i, j, c;
    ulong p = n_randprime(state, bits, 1);

    flint_printf("Hermite-Pade approximation (shifted Popov basis), %wu-bit prime\n", bits);
    flint_printf("    m     n     d      D      mul   approximant_basis_krylov  (/mul)\n");

    for (k = 0; cases[k][0] != 0; k++)
    {
        slong m = cases[k][0], n = cases[k][1], d = cases[k][2];
        slong * orders = flint_malloc(n * sizeof(slong));
        nmod_poly_mat_t F, P;
        double t, tmul;

        nmod_poly_mat_init(F, m, n, p);
        nmod_poly_mat_init(P, m, m, p);
        for (i = 0; i < m; i++)
            for (j = 0; j < n; j++)
            {
                nmod_poly_struct * f = nmod_poly_mat_entry(F, i, j);
                nmod_poly_fit_length(f, d);
                for (c = 0; c < d; c++)
                    f->coeffs[c] = n_randint(state, p);
                _nmod_poly_set_length(f, d);
                _nmod_poly_normalise(f);
            }
        for (j = 0; j < n; j++)
            orders[j] = d;

        tmul = _time_mul(n * d, p, state);
        TIME_IT(t, nmod_poly_mat_approximant_basis_krylov(P, F, orders, NULL));

        flint_printf("%5wd %5wd %5wd %6wd  %.2e  %.2e  (%6.2f)\n",
            m, n, d, n * d, tmul, t, t / tmul);

        nmod_poly_mat_clear(F);
        nmod_poly_mat_clear(P);
        flint_free(orders);
    }
}

static void
profile_interp(flint_bitcnt_t bits, flint_rand_t state)
{
    slong cases[][2] = {{16, 16}, {64, 64}, {256, 256}, {16, 256}, {64, 512},
                        {512, 128}, {0, 0}};
    slong k, j;
    ulong p = n_randprime(state, bits, 1);

    flint_printf("interpolation at D points (shifted Popov basis), %wu-bit prime\n", bits);
    flint_printf("    m      D      mul   interpolant_basis_krylov  (/mul)\n");

    for (k = 0; cases[k][0] != 0; k++)
    {
        slong m = cases[k][0], D = cases[k][1];
        nn_ptr pts = _nmod_vec_init(D);
        nmod_mat_t V;
        nmod_poly_mat_t P;
        double t, tmul;

        nmod_mat_init(V, m, D, p);
        nmod_mat_rand(V, state);
        nmod_poly_mat_init(P, m, m, p);
        for (j = 0; j < D; j++)
            pts[j] = n_randint(state, p);

        tmul = _time_mul(D, p, state);
        TIME_IT(t, nmod_poly_mat_interpolant_basis_krylov(P, V, pts, NULL));

        flint_printf("%5wd %6wd  %.2e  %.2e  (%6.2f)\n", m, D, tmul, t, t / tmul);

        nmod_mat_clear(V);
        nmod_poly_mat_clear(P);
        _nmod_vec_clear(pts);
    }
}

static void
profile_chunk(flint_bitcnt_t bits, flint_rand_t state)
{
    slong Ds[] = {64, 128, 256, 512, 1024, 0};
    slong chunks[] = {16, 32, 64, 128, 256, 0};
    slong k, c, i;
    ulong p = n_randprime(state, bits, 1);

    flint_printf("chunk size, Hermite-Pade with m = D (compact kernel), %wu-bit prime\n", bits);
    flint_printf("    D");
    for (c = 0; chunks[c] != 0; c++)
        flint_printf("   chunk %3wd", chunks[c]);
    flint_printf("\n");

    for (k = 0; Ds[k] != 0; k++)
    {
        slong D = Ds[k], m = Ds[k];
        slong * rows = flint_malloc(D * sizeof(slong));
        slong * exps = flint_malloc(D * sizeof(slong));
        slong * deltas = flint_malloc(m * sizeof(slong));
        nmod_mat_t V, R;
        _nmod_mat_krylov_op_t op;

        nmod_mat_init(V, m, D, p);
        nmod_mat_init(R, m, D, p);
        nmod_mat_rand(V, state);
        i = D;
        _nmod_mat_krylov_op_init_jordan(op, NULL, &i, 1, V->mod);

        flint_printf("%5wd", D);
        for (c = 0; chunks[c] != 0; c++)
        {
            double t;
            TIME_IT(t, _nmod_mat_krylov_kernel(R, rows, exps, deltas, V, op,
                        NULL, NULL, chunks[c], _NMOD_MAT_KRYLOV_KERNEL_AUTO));
            flint_printf("   %.3e", t);
        }
        flint_printf("\n");

        _nmod_mat_krylov_op_clear(op);
        nmod_mat_clear(V);
        nmod_mat_clear(R);
        flint_free(rows);
        flint_free(exps);
        flint_free(deltas);
    }
}

int
main(int argc, char ** argv)
{
    const char * mode = (argc > 1) ? argv[1] : "dense";
    flint_bitcnt_t bits = (argc > 2) ? atoi(argv[2]) : 60;
    flint_rand_t state;

    flint_rand_init(state);

    if (!strcmp(mode, "dense"))
        profile_dense(bits, state);
    else if (!strcmp(mode, "hp"))
        profile_hp(bits, state);
    else if (!strcmp(mode, "interp"))
        profile_interp(bits, state);
    else if (!strcmp(mode, "chunk"))
        profile_chunk(bits, state);
    else if (!strcmp(mode, "kernel"))
        profile_kernel(bits, state);
    else if (!strcmp(mode, "pow"))
        profile_pow(bits, (argc > 3) ? atol(argv[3]) : 1600, state);
    else if (!strcmp(mode, "budget"))
        profile_budget(bits, (argc > 3) ? atol(argv[3]) : 1024, state);
    else
        flint_printf("unknown mode %s (dense, hp, interp, chunk, kernel, pow, budget)\n", mode);

    flint_rand_clear(state);
    flint_cleanup_master();
    return 0;
}
