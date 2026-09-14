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
   Shared by the tuning profiles p-mul_tune and p-mulmid_tune: the timing
   of one parameter point within a time budget, and the sweep of a grid of
   points, pruned so that no point which the points already timed predict
   to exceed the budget is even attempted. The profiles differ in what a
   point is (which lengths, which window of coefficients), how a cell of
   the grid maps to one, and which routine runs it -- that is what the
   callbacks of tune_grid_struct provide.

   Two ways of reading the results. With a single algorithm, a table with
   the rows and columns of the grid, whose shape shows where the algorithm
   gives up. With a list of algorithms, one line per point with the
   parameters followed by one column per algorithm, so that the fastest one
   at each point is read off directly; each algorithm is still pruned on
   its own. The former is the cleaner measurement -- these variants keep
   large scratch buffers, and several of them in one process share the
   cache and the memory bandwidth -- the latter the easier to read.

   The budget is enforced before the fact as well as after it, which
   matters for the quadratic variants: vandermonde1, vandermonde2 and
   waksman can be slower than the fft variants by two orders of magnitude
   at the top of the grid, and a single call at such a point can run for
   an hour. A point is attempted only when the points already timed at a
   smaller dimension or a shorter length predict that it fits the budget,
   the prediction extrapolating each of them with the growth rate
   measured from its own two predecessors. Nothing is ever repeated once
   the budget is spent, so the worst case is one call at a point the
   prediction underestimated.
*/

#ifndef NMOD_POLY_MAT_EXTRA_TUNE_IMPL_H
#define NMOD_POLY_MAT_EXTRA_TUNE_IMPL_H

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <flint/longlong.h>
#include <flint/nmod_poly_mat.h>
#include <flint/profiler.h>
#include <flint/ulong_extras.h>

#include "nmod_poly_mat_utils.h"

/* one of the primes of FLINT's default fft_small context: 50 bits, and
   enough 2-adicity for the transform lengths reachable here */
#define TUNE_FFT_PRIME UWORD(1108307720798209)

/* the dimensions of the square grid */
static const slong TUNE_SQ_DIMS[] = { 2, 4, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256, 384, 512, 1024 };
#define TUNE_NSQ_DIMS (slong)(sizeof(TUNE_SQ_DIMS) / sizeof(TUNE_SQ_DIMS[0]))

/* what a cell holds when it was not timed; either prints as "-", but
   only TUNE_SKIPPED propagates to the points beyond it */
#define TUNE_NOT_APPLICABLE (-1.0)   /* the algorithm needs a larger field */
#define TUNE_SKIPPED        (-2.0)   /* too large, or predicted too slow */

/* a parameter point: the product of a d1 x d2 matrix of length len1 by a
   d2 x d3 matrix of length len2, the coefficients nlo, ..., nhi-1 of it
   being wanted (a plain product has nlo = 0, nhi = len1 + len2 - 1) */
typedef struct
{
    slong d1, d2, d3;
    slong len1, len2;
    slong nlo, nhi;
}
tune_point_struct;

/* runs the routine under test at a point, C = the wanted part of A * B */
typedef void (* tune_fun)(nmod_poly_mat_t C, const nmod_poly_mat_t A,
                          const nmod_poly_mat_t B, const tune_point_struct * P);

/* bytes of the operands and of the result: their coefficients, and the
   nmod_poly_struct and its allocation header for each entry -- at a
   large dimension and a short length the latter is the larger of the two
   and leaving it out makes the memory limit meaningless there */
static double tune_opbytes(const tune_point_struct * P)
{
    const double entries = (double) P->d1 * P->d2 + (double) P->d2 * P->d3
                           + (double) P->d1 * P->d3;
    const double coeffs = (double) P->d1 * P->d2 * P->len1
                          + (double) P->d2 * P->d3 * P->len2
                          + (double) P->d1 * P->d3 * (P->nhi - P->nlo);

    return 8.0 * coeffs + (double) (sizeof(nmod_poly_struct) + 32) * entries;
}

/*
    Minimum wall time over a few repetitions. The first call is timed as
    well -- it is the only guard against a point the prediction
    underestimated -- and no call is started once the budget is spent, so
    a point costs at most about `budget` seconds whatever it turns out to
    be.
*/
static double tune_time_one(tune_fun fun, const tune_point_struct * P,
                            ulong modn, double budget)
{
    flint_rand_t state;
    nmod_poly_mat_t A, B, C;
    timeit_t timer;
    double best, total;
    slong reps, i;

    flint_rand_init(state);
    flint_rand_set_seed(state, 1234, 5678);

    nmod_poly_mat_init(A, P->d1, P->d2, modn);
    nmod_poly_mat_init(B, P->d2, P->d3, modn);
    nmod_poly_mat_init(C, P->d1, P->d3, modn);
    nmod_poly_mat_rand(A, state, P->len1);
    nmod_poly_mat_rand(B, state, P->len2);

    timeit_start_us(timer);
    fun(C, A, B, P);
    timeit_stop_us(timer);
    best = 1e-6 * timer->wall;
    total = best;

    reps = (slong) (0.5 / (best > 1e-6 ? best : 1e-6));
    reps = FLINT_MAX(reps, 1);
    reps = FLINT_MIN(reps, 20);

    for (i = 0; i < reps && total + best <= budget; i++)
    {
        double t;
        timeit_start_us(timer);
        fun(C, A, B, P);
        timeit_stop_us(timer);
        t = 1e-6 * timer->wall;
        best = FLINT_MIN(best, t);
        total += t;
    }

    nmod_poly_mat_clear(A);
    nmod_poly_mat_clear(B);
    nmod_poly_mat_clear(C);
    flint_rand_clear(state);

    return best;
}

/*
    t1 is the time at x1 and t0, if positive, the time at x0 < x1: the
    time at x > x1, extrapolated with the growth rate between those two,
    and with defexp when there is only one of them. The rate is clamped:
    at the cheap end of the grid the timer resolution is a microsecond
    and a ratio of two such readings says nothing, and a measured rate
    below the floor is a cache effect of one particular step rather than
    the growth of the algorithm -- no variant here is subcubic in the
    dimension by more than the Strassen exponent. *fitted is set when the
    rate was measured rather than assumed.
*/
static double tune_extrapolate(double t0, double t1, slong x0, slong x1, slong x,
                               double defexp, double minexp, double maxexp,
                               int * fitted)
{
    double e = defexp;

    *fitted = 0;
    if (t0 > 0.0 && t1 > 0.0 && x1 > x0 && t1 > 4e-6)
    {
        e = log(t1 / t0) / log((double) x1 / (double) x0);
        e = FLINT_MAX(e, minexp);
        e = FLINT_MIN(e, maxexp);
        *fitted = 1;
    }

    return t1 * pow((double) x / (double) x1, e);
}

/* the sentinels above are the only negative values a cell can hold, so
   a reading of zero -- a product that took less than the microsecond the
   timer resolves -- still prints as a number */
static void tune_print_cell(double t)
{
    if (t >= 0.0)
        flint_printf(" %10.3e", t);
    else
        flint_printf(" %10s", "-");
    fflush(stdout);
}

/*
    A grid of points. Cell (i, j) is the point point(i, j); select(alg,
    modn, P) returns the routine to run algorithm alg there, or NULL when
    it does not apply, and may adjust the point to what that routine
    actually computes (the memory estimate is for the adjusted point, the
    printed parameters are the original ones).

    The pruning needs, for a cell, the cells it is a monotone step away
    from: row i - rowstep (a larger dimension, coordinate rowx[i]) and
    column j - 1 (a longer length, coordinate colx[j]); a rowstep of one
    is a plain table, while the rectangular grid of p-mul_tune has its
    shapes interleaved and steps over all of them. The labels are for the
    single-algorithm table.
*/
typedef struct
{
    slong nrows;
    slong ncols;
    slong rowstep;
    const slong * rowx;
    const slong * colx;
    const char * corner;
    void (* rowlabel)(slong i, char * buf);
    void (* collabel)(slong j, char * buf);
    void (* point)(slong i, slong j, tune_point_struct * P);
    tune_fun (* select)(const char * alg, ulong modn, tune_point_struct * P);
}
tune_grid_struct;

static double tune_predict(const tune_grid_struct * G, const double * T,
                           slong i, slong j)
{
    const slong nc = G->ncols, rs = G->rowstep;
    double pred = 0.0, guess = 0.0, p;
    int fitted;

#define _T(a, b) T[(a) * nc + (b)]

    if (i >= rs && _T(i - rs, j) > 0.0)
    {
        p = tune_extrapolate(i >= 2 * rs ? _T(i - 2 * rs, j) : -1.0, _T(i - rs, j),
                             i >= 2 * rs ? G->rowx[i - 2 * rs] : 0,
                             G->rowx[i - rs], G->rowx[i], 3.0, 2.5, 3.5, &fitted);
        if (fitted)
            pred = FLINT_MAX(pred, p);
        else
            guess = FLINT_MAX(guess, p);
    }

    if (j > 0 && _T(i, j - 1) > 0.0)
    {
        p = tune_extrapolate(j > 1 ? _T(i, j - 2) : -1.0, _T(i, j - 1),
                             j > 1 ? G->colx[j - 2] : 0, G->colx[j - 1], G->colx[j],
                             1.5, 0.5, 2.5, &fitted);
        if (fitted)
            pred = FLINT_MAX(pred, p);
        else
            guess = FLINT_MAX(guess, p);
    }

#undef _T

    return (pred > 0.0) ? pred : guess;
}

/* the cell (i, j) of algorithm alg, given its table T so far */
static double tune_cell(const tune_grid_struct * G, const char * alg, double * T,
                        slong i, slong j, ulong modn, double budget, double membytes)
{
    const slong nc = G->ncols, rs = G->rowstep;
    tune_point_struct P;
    tune_fun fun;

    G->point(i, j, &P);
    fun = G->select(alg, modn, &P);
    if (fun == NULL)
        return TUNE_NOT_APPLICABLE;

    /* neither parameter makes the product cheaper, so a point that was
       skipped, or that spent its whole budget, stops the column and the
       row at once */
    if ((i >= rs && (T[(i - rs) * nc + j] == TUNE_SKIPPED || T[(i - rs) * nc + j] > budget))
        || (j > 0 && (T[i * nc + j - 1] == TUNE_SKIPPED || T[i * nc + j - 1] > budget))
        || tune_opbytes(&P) > membytes
        || tune_predict(G, T, i, j) > budget)
        return TUNE_SKIPPED;

    return tune_time_one(fun, &P, modn, budget);
}

/*
    Runs the nalgs algorithms over the grid. Not `lines`: the table of the
    grid, for a single algorithm. `lines`: one line per point, the
    parameters (with nlo and nhi when lohi is set) followed by one column
    per algorithm.
*/
static void tune_run(const tune_grid_struct * G, const char ** algs, slong nalgs,
                     ulong modn, slong nthreads, double budget, double membytes,
                     int lines, int lohi)
{
    const slong nr = G->nrows, nc = G->ncols;
    double * T = flint_malloc(nalgs * nr * nc * sizeof(double));
    char buf[64];
    slong a, i, j;

    if (!lines)
    {
        flint_printf("%-*s", (int) strlen(G->corner), G->corner);
        for (j = 0; j < nc; j++)
        {
            G->collabel(j, buf);
            flint_printf(" %10s", buf);
        }
        flint_printf("\n");

        for (i = 0; i < nr; i++)
        {
            G->rowlabel(i, buf);
            flint_printf("%-*s", (int) strlen(G->corner), buf);
            for (j = 0; j < nc; j++)
            {
                T[i * nc + j] = tune_cell(G, algs[0], T, i, j, modn, budget, membytes);
                tune_print_cell(T[i * nc + j]);
            }
            flint_printf("\n");
        }
    }
    else
    {
        flint_printf("%20s %8s %6s %6s %6s %6s %6s", "prime", "nthreads",
                     "dim1", "dim2", "dim3", "len1", "len2");
        if (lohi)
            flint_printf(" %6s %6s", "nlo", "nhi");
        for (a = 0; a < nalgs; a++)
            flint_printf(" %14s", algs[a]);
        flint_printf("\n");

        for (i = 0; i < nr; i++)
            for (j = 0; j < nc; j++)
            {
                tune_point_struct P;

                G->point(i, j, &P);
                flint_printf("%20wu %8wd %6wd %6wd %6wd %6wd %6wd", modn, nthreads,
                             P.d1, P.d2, P.d3, P.len1, P.len2);
                if (lohi)
                    flint_printf(" %6wd %6wd", P.nlo, P.nhi);
                for (a = 0; a < nalgs; a++)
                {
                    double * Ta = T + a * nr * nc;
                    Ta[i * nc + j] = tune_cell(G, algs[a], Ta, i, j, modn, budget, membytes);
                    flint_printf("    ");
                    tune_print_cell(Ta[i * nc + j]);
                }
                flint_printf("\n");
            }
    }

    flint_free(T);
}

/*
    The algorithm argument: a name, or a bracketed list of names --
    "[a,b]", "[a, b]", "[ a b ]" are all accepted, and "[]" stands for
    every algorithm the profile knows. Fills algs (at most maxalgs of
    them, from `all` for the empty list) and returns the index of the
    first argument after it, or -1 when a name is not one of `all`.
    *listed is set when the argument was a list, even of one name: that
    asks for the line-per-point output.
*/
static int tune_parse_algs(int argc, char ** argv, int first,
                           const char ** all, slong nall,
                           const char ** algs, slong * nalgs, int * listed)
{
    static char text[4096];
    slong n = 0, k;
    int i = first, done = 0;
    char * tok;

    *listed = (argv[first][0] == '[');
    if (!*listed)
    {
        algs[0] = argv[first];
        *nalgs = 1;
        for (k = 0; k < nall; k++)
            if (!strcmp(algs[0], all[k]))
                return first + 1;
        return -1;
    }

    /* gather the list, which may span several arguments */
    text[0] = '\0';
    while (i < argc && !done)
    {
        if (strlen(text) + strlen(argv[i]) + 2 > sizeof(text))
            return -1;
        strcat(text, argv[i]);
        strcat(text, " ");
        done = (strchr(argv[i], ']') != NULL);
        i++;
    }
    if (!done)
        return -1;

    for (tok = strtok(text, "[], \t"); tok != NULL; tok = strtok(NULL, "[], \t"))
    {
        for (k = 0; k < nall; k++)
            if (!strcmp(tok, all[k]))
                break;
        if (k == nall || n == nall)
            return -1;
        algs[n++] = all[k];
    }

    if (n == 0)   /* "[]": all of them */
        for (n = 0; n < nall; n++)
            algs[n] = all[n];

    *nalgs = n;
    return i;
}

/* the header lines common to the two profiles */
static void tune_print_header(const char ** algs, slong nalgs, ulong modn,
                              slong nthreads, double budget, double membytes)
{
    slong a;

    flint_printf("#");
    for (a = 0; a < nalgs; a++)
        flint_printf(" %s", algs[a]);
    flint_printf(", modulus %wu (%wu bits, 2-adicity %wu), %wd thread%s, "
                 "budget %gs, mem %.1fGB\n",
                 modn, FLINT_BIT_COUNT(modn),
                 modn > 1 ? (ulong) flint_ctz(modn - 1) : UWORD(0),
                 nthreads, nthreads == 1 ? "" : "s",
                 budget, membytes / (double) (1 << 30));
    flint_printf("# a modulus of at most 50 bits whose 2-adicity reaches the "
                 "transform depth is\n# transformed directly by fft_small, "
                 "with no chinese remaindering\n");
}

#endif  /* NMOD_POLY_MAT_EXTRA_TUNE_IMPL_H */
