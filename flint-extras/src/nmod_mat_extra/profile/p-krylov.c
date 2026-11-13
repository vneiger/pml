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
                     uniform and "Hermite" shifts, for 1, n/8 and n rows
      mode = hp      Hermite-Pade: approximant basis of an m x n matrix at
                     order d in each column (Jordan matrix with n nilpotent
                     blocks of size d), in shifted Popov form
      mode = interp  interpolation at sigma points (diagonal matrix)
      mode = chunk   effect of the chunk size (Hermite-Pade, m = sigma)
      mode = pow     dense random M, generation of the rows from powers
                     M^(2^t) (Keller-Gehrig-type doubling), n up to nmax:
                     p-krylov pow [bits] [nmax]

    Times in seconds; "mul" is the time of nmod_mat_mul for n x n (dense)
    or sigma x sigma (hp, interp) matrices.
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
_time_mul(slong n, ulong p, flint_rand_t state)
{
    nmod_mat_t A, B, C;
    double t;

    nmod_mat_init(A, n, n, p);
    nmod_mat_init(B, n, n, p);
    nmod_mat_init(C, n, n, p);
    nmod_mat_rand(A, state);
    nmod_mat_rand(B, state);
    TIME_IT(t, nmod_mat_mul(C, A, B));
    nmod_mat_clear(A);
    nmod_mat_clear(B);
    nmod_mat_clear(C);

    return t;
}

/* time of c x n times n x n */
static double
_time_mul_rows(slong c, slong n, ulong p, flint_rand_t state)
{
    nmod_mat_t A, B, C;
    double t;

    nmod_mat_init(A, c, n, p);
    nmod_mat_init(B, n, n, p);
    nmod_mat_init(C, c, n, p);
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
    slong ns[] = {50, 100, 200, 400, 800, 0};
    slong k, j, s;
    ulong p = n_randprime(state, bits, 1);

    flint_printf("dense random M, %wu-bit prime\n", bits);
    flint_printf("    n     m  shift     mul    basis  doubling  kernel   (basis/mul  kernel/mul  doubling/basis)\n");

    for (k = 0; ns[k] != 0; k++)
    {
        slong n = ns[k];
        slong ms[3];
        double tmul = _time_mul(n, p, state);

        ms[0] = 1;
        ms[1] = FLINT_MAX(1, n / 8);
        ms[2] = n;

        for (j = 0; j < 3; j++)
        {
            for (s = 0; s < 2; s++)
            {
                slong m = ms[j], i, r = 0;
                slong * shifts = flint_malloc(m * sizeof(slong));
                slong * rows = flint_malloc(n * sizeof(slong));
                slong * exps = flint_malloc(n * sizeof(slong));
                slong * deltas = flint_malloc(m * sizeof(slong));
                nmod_mat_t E, M, B, R;
                _nmod_mat_krylov_op_t op;
                double tb, td, tk;
                int ok;

                if (s == 1 && m == 1)
                {
                    flint_free(shifts); flint_free(rows); flint_free(exps); flint_free(deltas);
                    continue;
                }

                for (i = 0; i < m; i++)
                    shifts[i] = (s == 0) ? 0 : i * (n + 1);

                nmod_mat_init(E, m, n, p);
                nmod_mat_init(M, n, n, p);
                nmod_mat_init(B, n, n, p);
                nmod_mat_init(R, m, n, p);
                nmod_mat_rand(E, state);
                nmod_mat_rand(M, state);
                _nmod_mat_krylov_op_init_dense(op, M);

                TIME_IT(tb, r = _nmod_mat_krylov_basis_iter(B, rows, exps, NULL,
                            NULL, NULL, NULL, &ok, E, op, shifts, NULL, 0));
                TIME_IT(td, _nmod_mat_krylov_basis_doubling(B, rows, exps, NULL,
                            &ok, E, op, shifts, NULL));
                TIME_IT(tk, _nmod_mat_krylov_kernel(R, rows, exps, deltas, E, op,
                            shifts, NULL, 0));

                flint_printf("%5wd %5wd  %s  %.2e %.2e %.2e %.2e   (%6.2f %6.2f %6.2f)   rank %wd\n",
                    n, m, s == 0 ? "unif" : "herm", tmul, tb, td, tk,
                    tb / tmul, tk / tmul, td / tb, r);

                _nmod_mat_krylov_op_clear(op);
                nmod_mat_clear(E);
                nmod_mat_clear(M);
                nmod_mat_clear(B);
                nmod_mat_clear(R);
                flint_free(shifts);
                flint_free(rows);
                flint_free(exps);
                flint_free(deltas);
            }
        }
    }
}

/*
    Generation of the rows from powers M^(2^t) (op->logpow), dense random
    M: Hermite shift (the iterates of each row of E in turn), or uniform
    shift with few rows. Times relative to logpow = 0; "KG" uses chunks of
    up to m + n rows and all the powers needed (Keller-Gehrig doubling with
    incremental elimination); "auto" is the default (automatic choice);
    "paper" is _nmod_mat_krylov_basis_doubling.
*/
static void
profile_pow(flint_bitcnt_t bits, slong nmax, flint_rand_t state)
{
    slong ns[] = {200, 400, 800, 1600, 3200, 6400, 0};
    slong cfg_m[] = {1, 4, 2, 4, 16};
    int cfg_herm[] = {0, 1, 0, 0, 0};
    slong lps[] = {1, 2, 3, 4, 5, 7};
    slong nlp = 6, k, c, l, i;
    ulong p = n_randprime(state, bits, 1);

    flint_printf("powers of M for generating rows, dense random M, %wu-bit prime\n", bits);
    flint_printf("rho1, rho8: time of n rows as products of 1 or 8 rows by M, relative to mul\n");
    flint_printf("    n     m shift      mul  logpow=0  | ratio for logpow =");
    for (l = 0; l < nlp; l++)
        flint_printf(" %4wd", lps[l]);
    flint_printf("    KG  auto  paper\n");

    for (k = 0; ns[k] != 0 && ns[k] <= nmax; k++)
    {
        slong n = ns[k];
        double tmul = _time_mul(n, p, state);
        double rho1 = n * _time_mul_rows(1, n, p, state) / tmul;
        double rho8 = (n / 8) * _time_mul_rows(8, n, p, state) / tmul;

        flint_printf("n = %wd: rho1 = %.2f, rho8 = %.2f\n", n, rho1, rho8);

        for (c = 0; c < 5; c++)
        {
            slong m = cfg_m[c], r0, r;
            slong * shifts = flint_malloc(m * sizeof(slong));
            slong * rows0 = flint_malloc(n * sizeof(slong));
            slong * exps0 = flint_malloc(n * sizeof(slong));
            slong * rows = flint_malloc(n * sizeof(slong));
            slong * exps = flint_malloc(n * sizeof(slong));
            nmod_mat_t E, M, B;
            _nmod_mat_krylov_op_t op;
            double t0, t;
            int ok, bad = 0;

            for (i = 0; i < m; i++)
                shifts[i] = cfg_herm[c] ? i * (n + 1) : 0;

            nmod_mat_init(E, m, n, p);
            nmod_mat_init(M, n, n, p);
            nmod_mat_init(B, n, n, p);
            nmod_mat_rand(E, state);
            nmod_mat_rand(M, state);
            _nmod_mat_krylov_op_init_dense(op, M);

            op->logpow = 0;
            TIME_IT(t0, r0 = _nmod_mat_krylov_basis_iter(B, rows0, exps0, NULL,
                        NULL, NULL, NULL, &ok, E, op, shifts, NULL, 0));

            flint_printf("%5wd %5wd  %s  %.2e  %.2e  |                  ",
                    n, m, cfg_herm[c] ? "herm" : "unif", tmul, t0);

            for (l = 0; l <= nlp + 1; l++)
            {
                slong chunk = (l == nlp) ? m + n : 0;

                op->logpow = (l < nlp) ? lps[l] : (l == nlp) ? FLINT_BITS : -1;
                TIME_IT(t, r = _nmod_mat_krylov_basis_iter(B, rows, exps, NULL,
                            NULL, NULL, NULL, &ok, E, op, shifts, NULL, chunk));
                flint_printf(" %5.2f", t / t0);
                fflush(stdout);

                if (r != r0)
                    bad = 1;
                for (i = 0; i < r && !bad; i++)
                    if (rows[i] != rows0[i] || exps[i] != exps0[i])
                        bad = 1;
            }

            if (n <= 1600)
            {
                TIME_IT(t, r = _nmod_mat_krylov_basis_doubling(B, rows, exps, NULL,
                            &ok, E, op, shifts, NULL));
                flint_printf("  %5.2f", t / t0);
                if (r != r0)
                    bad = 1;
            }
            else
                flint_printf("      -");

            flint_printf("%s\n", bad ? "  MISMATCH" : "");

            _nmod_mat_krylov_op_clear(op);
            nmod_mat_clear(E);
            nmod_mat_clear(M);
            nmod_mat_clear(B);
            flint_free(shifts);
            flint_free(rows0);
            flint_free(exps0);
            flint_free(rows);
            flint_free(exps);
        }
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
    flint_printf("    m     n     d  sigma      mul   approximant_basis_krylov  (/mul)\n");

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

    flint_printf("interpolation at sigma points (shifted Popov basis), %wu-bit prime\n", bits);
    flint_printf("    m  sigma      mul   interpolant_basis_krylov  (/mul)\n");

    for (k = 0; cases[k][0] != 0; k++)
    {
        slong m = cases[k][0], sigma = cases[k][1];
        nn_ptr pts = _nmod_vec_init(sigma);
        nmod_mat_t V;
        nmod_poly_mat_t P;
        double t, tmul;

        nmod_mat_init(V, m, sigma, p);
        nmod_mat_rand(V, state);
        nmod_poly_mat_init(P, m, m, p);
        for (j = 0; j < sigma; j++)
            pts[j] = n_randint(state, p);

        tmul = _time_mul(sigma, p, state);
        TIME_IT(t, nmod_poly_mat_interpolant_basis_krylov(P, V, pts, NULL));

        flint_printf("%5wd %6wd  %.2e  %.2e  (%6.2f)\n", m, sigma, tmul, t, t / tmul);

        nmod_mat_clear(V);
        nmod_poly_mat_clear(P);
        _nmod_vec_clear(pts);
    }
}

static void
profile_chunk(flint_bitcnt_t bits, flint_rand_t state)
{
    slong sigmas[] = {64, 128, 256, 512, 1024, 0};
    slong chunks[] = {16, 32, 64, 128, 256, 0};
    slong k, c, i;
    ulong p = n_randprime(state, bits, 1);

    flint_printf("chunk size, Hermite-Pade with m = sigma (compact kernel), %wu-bit prime\n", bits);
    flint_printf("sigma");
    for (c = 0; chunks[c] != 0; c++)
        flint_printf("   chunk %3wd", chunks[c]);
    flint_printf("\n");

    for (k = 0; sigmas[k] != 0; k++)
    {
        slong n = sigmas[k], m = sigmas[k];
        slong * rows = flint_malloc(n * sizeof(slong));
        slong * exps = flint_malloc(n * sizeof(slong));
        slong * deltas = flint_malloc(m * sizeof(slong));
        nmod_mat_t E, R;
        _nmod_mat_krylov_op_t op;

        nmod_mat_init(E, m, n, p);
        nmod_mat_init(R, m, n, p);
        nmod_mat_rand(E, state);
        i = n;
        _nmod_mat_krylov_op_init_jordan(op, NULL, &i, 1, E->mod);

        flint_printf("%5wd", n);
        for (c = 0; chunks[c] != 0; c++)
        {
            double t;
            TIME_IT(t, _nmod_mat_krylov_kernel(R, rows, exps, deltas, E, op,
                        NULL, NULL, chunks[c]));
            flint_printf("   %.3e", t);
        }
        flint_printf("\n");

        _nmod_mat_krylov_op_clear(op);
        nmod_mat_clear(E);
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
    else if (!strcmp(mode, "pow"))
        profile_pow(bits, (argc > 3) ? atol(argv[3]) : 1600, state);
    else
        flint_printf("unknown mode %s (dense, hp, interp, chunk, pow)\n", mode);

    flint_rand_clear(state);
    flint_cleanup_master();
    return 0;
}
