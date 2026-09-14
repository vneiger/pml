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
   Timings of the polynomial matrix multiplication variants that
   nmod_poly_mat_multiply chooses between, laid out so that the
   crossovers between them can be read off.

   One algorithm per run: these variants allocate the evaluations of
   their operands, up to hundreds of megabytes, and running several of
   them in a single process measures the memory subsystem as much as the
   algorithms. Run the same command once per algorithm and compare the
   tables.

   Usage:
     p-mul_tune prime nthreads fun [opts]
         table of timings over a grid of parameters
     p-mul_tune prime nthreads fun dim1 dim2 dim3 len1 len2
         one timing, for the product of a dim1 x dim2 matrix of length
         len1 by a dim2 x dim3 matrix of length len2

     prime      the modulus; 0 selects one of FLINT's 50-bit FFT primes,
                for which the fft_small variants use a single transform
                modulo the prime itself and no chinese remaindering
     nthreads   number of threads
     fun        sd_fft_direct | sd_fft_matmul | geometric
                | vandermonde1 | vandermonde2 | waksman | multiply

   Options of the table mode, in any order:
     rect         sweep rectangular shapes and unbalanced lengths
                  instead of the square grid
     budget=SECS  spend at most about SECS seconds on any one parameter
                  point (default 20)
     mem=GB       skip a parameter point whose operands and result would
                  exceed about GB gigabytes (default 8)

   A cell reads "-" when the algorithm does not apply to that modulus and
   length, or when the point was skipped by one of the two limits above.

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

/* bash script for an overall tune (possibly adapt budget=10) */

/*
# ╰─ ./build/nmod_poly_mat_extra/profile/p-mul_tune
# Usage: ./build/nmod_poly_mat_extra/profile/p-mul_tune prime nthreads fun [opts]
#        ./build/nmod_poly_mat_extra/profile/p-mul_tune prime nthreads fun dim1 dim2 dim3 len1 len2
#    prime: the modulus; 0 selects a 50-bit FFT prime
#    fun:   sd_fft_direct | sd_fft_matmul | geometric | multiply
#    opts:  rect | budget=SECS | mem=GB

# 21-bit FFT prime + 50 bit FFT prime + 30 bit prime + 60 bit prime
for nthreads in 1 2 4 8
do
	for fun in "sd_fft_direct" "sd_fft_matmul" "geometric" "multiply"
	do
		for prime in 1179649 0 1073741789 1152921504606846883 
		do
			echo "--------------------------------------------------"
			echo "nthreads = $nthreads | fun = $fun | prime = $prime"
			./build/nmod_poly_mat_extra/profile/p-mul_tune $prime $nthreads $fun budget=10
			echo "--------------------------------------------------"
		done;
	done;
done
*/


#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <flint/longlong.h>
#include <flint/nmod_poly_mat.h>
#include <flint/profiler.h>
#include <flint/ulong_extras.h>

#include "nmod_poly_extra.h"
#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_multiply.h"

typedef void (*mulfun)(nmod_poly_mat_t, const nmod_poly_mat_t, const nmod_poly_mat_t);

/* one of the primes of FLINT's default fft_small context: 50 bits, and
   enough 2-adicity for the transform lengths reachable here */
#define FFT_PRIME UWORD(1108307720798209)

static const slong SQ_DIMS[] = { 2, 4, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256, 384, 512, 1024 };
static const slong SQ_LENS[] = { 2, 5, 8, 16, 32, 64, 128, 256, 512, 1024, 2048 };

/* rectangular shapes, as multiples of a base dimension d: the fat and
   thin products that the thresholds, tuned on square ones, extrapolate
   to */
static const slong RECT_SHAPES[][3] = {
    {1, 1, 1}, {1, 4, 1}, {4, 1, 4}, {1, 1, 4}, {4, 1, 1}, {2, 1, 1}
};
static const slong RECT_BASES[] = { 8, 32, 128 };
/* (len1, len2) in eighths of a base length */
static const slong RECT_LENS[][2] = { {8, 8}, {8, 1}, {1, 8} };
static const slong RECT_BASELENS[] = { 64, 512 };

#define NSQ_DIMS (slong)(sizeof(SQ_DIMS) / sizeof(SQ_DIMS[0]))
#define NSQ_LENS (slong)(sizeof(SQ_LENS) / sizeof(SQ_LENS[0]))
#define NRECT_SHAPES (slong)(sizeof(RECT_SHAPES) / sizeof(RECT_SHAPES[0]))
#define NRECT_BASES (slong)(sizeof(RECT_BASES) / sizeof(RECT_BASES[0]))
#define NRECT_LENS (slong)(sizeof(RECT_LENS) / sizeof(RECT_LENS[0]))
#define NRECT_BASELENS (slong)(sizeof(RECT_BASELENS) / sizeof(RECT_BASELENS[0]))

/* what a cell holds when it was not timed; any of these prints as "-",
   but only _SKIPPED propagates to the points beyond it */
#define _NOT_APPLICABLE (-1.0)   /* the algorithm needs a larger field */
#define _SKIPPED        (-2.0)   /* too large, or predicted too slow */

static mulfun _select(const char * alg, ulong modn, slong len)
{
    if (!strcmp(alg, "sd_fft_direct"))  return nmod_poly_mat_mul_sd_fft_direct;
    if (!strcmp(alg, "sd_fft_matmul"))  return nmod_poly_mat_mul_sd_fft_matmul;
    if (!strcmp(alg, "multiply"))       return nmod_poly_mat_multiply;
    if (!strcmp(alg, "geometric"))
        return NMOD_POLY_CAN_USE_GEOMETRIC(modn, len)
               ? nmod_poly_mat_mul_geometric : NULL;
    if (!strcmp(alg, "vandermonde1"))
        return NMOD_POLY_CAN_USE_VANDERMONDE1(modn, len)
               ? nmod_poly_mat_mul_vandermonde1 : NULL;
    if (!strcmp(alg, "vandermonde2"))
        return NMOD_POLY_CAN_USE_VANDERMONDE2(modn, len)
               ? nmod_poly_mat_mul_vandermonde2 : NULL;
    if (!strcmp(alg, "waksman"))
        return NMOD_POLY_MAT_CAN_USE_WAKSMAN(modn)
               ? nmod_poly_mat_mul_waksman : NULL;
    return NULL;
}

static int _known(const char * alg)
{
    return !strcmp(alg, "sd_fft_direct") || !strcmp(alg, "sd_fft_matmul")
        || !strcmp(alg, "multiply") || !strcmp(alg, "geometric")
        || !strcmp(alg, "vandermonde1") || !strcmp(alg, "vandermonde2")
        || !strcmp(alg, "waksman");
}

/* bytes of the operands and of the result: their coefficients, and the
   nmod_poly_struct and its allocation header for each entry -- at a
   large dimension and a short length the latter is the larger of the two
   and leaving it out makes the memory limit meaningless there */
static double _opbytes(slong d1, slong d2, slong d3, slong l1, slong l2)
{
    const double entries = (double) d1 * d2 + (double) d2 * d3
                           + (double) d1 * d3;
    const double coeffs = (double) d1 * d2 * l1 + (double) d2 * d3 * l2
                          + (double) d1 * d3 * (l1 + l2 - 1);

    return 8.0 * coeffs + (double) (sizeof(nmod_poly_struct) + 32) * entries;
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
static double _extrapolate(double t0, double t1, slong x0, slong x1, slong x,
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

/*
    An estimate of the time at (i, j) of the square grid from the points
    already timed in the same column and in the same row. The larger of
    the two is returned, and 0 when neither is available: the profile
    must not launch a product it cannot afford, and the cost of erring on
    that side is a missing cell rather than an hour of wall time.

    An estimate whose growth rate was measured is preferred to one that
    had to assume it, however: the grid jumps from length 3 to length 11
    and from dimension 512 to 1024, and an assumed rate over a step that
    wide would prune half the table.
*/
static double _predict_square(const double * T, slong i, slong j)
{
    const slong * D = SQ_DIMS;
    const slong * L = SQ_LENS;
    double pred = 0.0, guess = 0.0, p;
    int fitted;

#define _T(a, b) T[(a) * NSQ_LENS + (b)]

    if (i > 0 && _T(i - 1, j) > 0.0)
    {
        p = _extrapolate(i > 1 ? _T(i - 2, j) : -1.0, _T(i - 1, j),
                         i > 1 ? D[i - 2] : 0, D[i - 1], D[i],
                         3.0, 2.5, 3.5, &fitted);
        if (fitted)
            pred = FLINT_MAX(pred, p);
        else
            guess = FLINT_MAX(guess, p);
    }

    if (j > 0 && _T(i, j - 1) > 0.0)
    {
        p = _extrapolate(j > 1 ? _T(i, j - 2) : -1.0, _T(i, j - 1),
                         j > 1 ? 2 * L[j - 2] - 1 : 0, 2 * L[j - 1] - 1,
                         2 * L[j] - 1, 1.5, 0.5, 2.5, &fitted);
        if (fitted)
            pred = FLINT_MAX(pred, p);
        else
            guess = FLINT_MAX(guess, p);
    }

#undef _T

    return (pred > 0.0) ? pred : guess;
}

/*
    Minimum wall time over a few repetitions. The first call is timed as
    well -- it is the only guard against a point the prediction
    underestimated -- and no call is started once the budget is spent, so
    a point costs at most about `budget` seconds whatever it turns out to
    be.
*/
static double _time_one(mulfun fun, slong d1, slong d2, slong d3,
                        slong len1, slong len2, ulong modn, double budget)
{
    flint_rand_t state;
    nmod_poly_mat_t A, B, C;
    timeit_t timer;
    double best, total;
    slong reps, i;

    flint_rand_init(state);
    flint_rand_set_seed(state, 1234, 5678);

    nmod_poly_mat_init(A, d1, d2, modn);
    nmod_poly_mat_init(B, d2, d3, modn);
    nmod_poly_mat_init(C, d1, d3, modn);
    nmod_poly_mat_rand(A, state, len1);
    nmod_poly_mat_rand(B, state, len2);

    timeit_start_us(timer);
    fun(C, A, B);
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
        fun(C, A, B);
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

/* the sentinels above are the only negative values a cell can hold, so
   a reading of zero -- a product that took less than the microsecond the
   timer resolves -- still prints as a number */
static void _print_cell(double t)
{
    if (t >= 0.0)
        flint_printf(" %10.3e", t);
    else
        flint_printf(" %10s", "-");
    fflush(stdout);
}

static void _table_square(const char * alg, ulong modn, double budget, double membytes)
{
    double * T = flint_malloc(NSQ_DIMS * NSQ_LENS * sizeof(double));
    slong i, j;

#define _T(a, b) T[(a) * NSQ_LENS + (b)]

    flint_printf("%-6s", "dim\\len");
    for (j = 0; j < NSQ_LENS; j++)
        flint_printf(" %10wd", 2 * SQ_LENS[j] - 1);
    flint_printf("\n");

    for (i = 0; i < NSQ_DIMS; i++)
    {
        const slong d = SQ_DIMS[i];

        flint_printf("%-6wd", d);
        for (j = 0; j < NSQ_LENS; j++)
        {
            const slong l = SQ_LENS[j];
            mulfun fun = _select(alg, modn, 2 * l - 1);
            double pred;

            if (fun == NULL)
            {
                _T(i, j) = _NOT_APPLICABLE;
                _print_cell(_NOT_APPLICABLE);
                continue;
            }

            /* neither parameter makes the product cheaper, so a point
               that was skipped, or that spent its whole budget, stops
               the column and the row at once */
            if ((i > 0 && (_T(i - 1, j) == _SKIPPED || _T(i - 1, j) > budget))
                || (j > 0 && (_T(i, j - 1) == _SKIPPED || _T(i, j - 1) > budget))
                || _opbytes(d, d, d, l, l) > membytes)
            {
                _T(i, j) = _SKIPPED;
                _print_cell(_SKIPPED);
                continue;
            }

            pred = _predict_square(T, i, j);
            if (pred > budget)
            {
                _T(i, j) = _SKIPPED;
                _print_cell(_SKIPPED);
                continue;
            }

            _T(i, j) = _time_one(fun, d, d, d, l, l, modn, budget);
            _print_cell(_T(i, j));
        }
        flint_printf("\n");
    }

#undef _T

    flint_free(T);
}

static void _rect_len(slong c, slong * len1, slong * len2)
{
    const slong base = RECT_BASELENS[c / NRECT_LENS];
    const slong l = c % NRECT_LENS;

    *len1 = FLINT_MAX(base * RECT_LENS[l][0] / 8, 1);
    *len2 = FLINT_MAX(base * RECT_LENS[l][1] / 8, 1);
}

static void _table_rect(const char * alg, ulong modn, double budget, double membytes)
{
    const slong ncols = NRECT_BASELENS * NRECT_LENS;
    /*
        P[s][c] is what the same shape and column did at the previous
        base: its time, _SKIPPED if it was skipped there, or
        _NOT_APPLICABLE if the algorithm does not apply to that column.
        The shape multipliers are fixed, so the product grows as the cube
        of the base and a point skipped at one base is skipped at every
        larger one -- without this the first column of a row has nothing
        to predict from and the row runs unguarded.
    */
    double * P = flint_malloc(NRECT_SHAPES * ncols * sizeof(double));
    slong b, s, c;

    for (s = 0; s < NRECT_SHAPES * ncols; s++)
        P[s] = _NOT_APPLICABLE;

    flint_printf("%-18s", "m x k x n \\ len");
    for (c = 0; c < ncols; c++)
    {
        slong len1, len2;
        char lab[32];

        _rect_len(c, &len1, &len2);
        flint_sprintf(lab, "%wd+%wd", len1, len2);
        flint_printf(" %10s", lab);
    }
    flint_printf("\n");

    for (b = 0; b < NRECT_BASES; b++)
        for (s = 0; s < NRECT_SHAPES; s++)
        {
            const slong d = RECT_BASES[b];
            const slong d1 = d * RECT_SHAPES[s][0];
            const slong d2 = d * RECT_SHAPES[s][1];
            const slong d3 = d * RECT_SHAPES[s][2];
            double prev = -1.0, prevprev = -1.0;
            slong prevlen = 0, prevprevlen = 0;
            int stop = 0;
            char shape[32];

            flint_sprintf(shape, "%wd x %wd x %wd", d1, d2, d3);
            flint_printf("%-18s", shape);

            for (c = 0; c < ncols; c++)
            {
                const double base_prev = (b > 0) ? P[s * ncols + c]
                                                 : _NOT_APPLICABLE;
                slong len1, len2;
                mulfun fun;
                double pred = 0.0, t;

                _rect_len(c, &len1, &len2);
                fun = _select(alg, modn, len1 + len2 - 1);

                if (fun == NULL)
                {
                    P[s * ncols + c] = _NOT_APPLICABLE;
                    _print_cell(_NOT_APPLICABLE);
                    continue;
                }

                if (stop || base_prev == _SKIPPED || base_prev > budget
                    || _opbytes(d1, d2, d3, len1, len2) > membytes)
                {
                    stop = 1;
                    P[s * ncols + c] = _SKIPPED;
                    _print_cell(_SKIPPED);
                    continue;
                }

                /* from the previous column of this row, and from this
                   column at the previous base */
                if (prev > 0.0)
                {
                    int fitted;
                    pred = _extrapolate(prevprev, prev, prevprevlen, prevlen,
                                        len1 + len2 - 1, 1.5, 0.5, 2.5, &fitted);
                }
                if (base_prev > 0.0)
                {
                    const double q = base_prev
                        * pow((double) d / (double) RECT_BASES[b - 1], 3.0);
                    pred = FLINT_MAX(pred, q);
                }

                if (pred > budget)
                {
                    stop = 1;
                    P[s * ncols + c] = _SKIPPED;
                    _print_cell(_SKIPPED);
                    continue;
                }

                t = _time_one(fun, d1, d2, d3, len1, len2, modn, budget);
                P[s * ncols + c] = t;
                _print_cell(t);

                prevprev = prev; prevprevlen = prevlen;
                prev = t; prevlen = len1 + len2 - 1;
                if (t > budget)
                    stop = 1;
            }
            flint_printf("\n");
        }

    flint_free(P);
}

int main(int argc, char ** argv)
{
    ulong modn;
    slong nthreads;
    const char * alg;
    double budget = 20.0, membytes = 8.0 * (1 << 30);
    int rect = 0, i;

    if (argc < 4)
    {
        flint_printf("Usage: %s prime nthreads fun [opts]\n", argv[0]);
        flint_printf("       %s prime nthreads fun dim1 dim2 dim3 len1 len2\n", argv[0]);
        flint_printf("   prime: the modulus; 0 selects a 50-bit FFT prime\n");
        flint_printf("   fun:   sd_fft_direct | sd_fft_matmul | geometric\n");
        flint_printf("          | vandermonde1 | vandermonde2 | waksman | multiply\n");
        flint_printf("   opts:  rect | budget=SECS | mem=GB\n");
        return 0;
    }

    modn = strtoul(argv[1], NULL, 10);
    if (modn == 0)
        modn = FFT_PRIME;
    nthreads = atol(argv[2]);
    alg = argv[3];
    flint_set_num_threads(nthreads);

    if (!_known(alg))
    {
        flint_printf("unknown algorithm %s\n", alg);
        return 1;
    }

    if (argc == 9)   /* one timing */
    {
        const slong d1 = atol(argv[4]), d2 = atol(argv[5]), d3 = atol(argv[6]);
        const slong len1 = atol(argv[7]), len2 = atol(argv[8]);
        mulfun fun = _select(alg, modn, len1 + len2 - 1);

        if (fun == NULL)
            flint_printf("-\n");
        else
            flint_printf("%.3e\n", _time_one(fun, d1, d2, d3, len1, len2, modn,
                                             DBL_MAX));
        return 0;
    }

    for (i = 4; i < argc; i++)
    {
        if (!strcmp(argv[i], "rect"))
            rect = 1;
        else if (!strncmp(argv[i], "budget=", 7))
            budget = atof(argv[i] + 7);
        else if (!strncmp(argv[i], "mem=", 4))
            membytes = atof(argv[i] + 4) * (double) (1 << 30);
        else
        {
            flint_printf("unknown option %s\n", argv[i]);
            return 1;
        }
    }

    flint_printf("# %s, modulus %wu (%wu bits, 2-adicity %wu), %wd thread%s, "
                 "budget %gs, mem %.1fGB\n",
                 alg, modn, FLINT_BIT_COUNT(modn),
                 modn > 1 ? (ulong) flint_ctz(modn - 1) : UWORD(0),
                 nthreads, nthreads == 1 ? "" : "s",
                 budget, membytes / (double) (1 << 30));
    flint_printf("# a modulus of at most 50 bits whose 2-adicity reaches the "
                 "transform depth is\n# transformed directly by fft_small, "
                 "with no chinese remaindering\n");
    flint_printf("# columns are len1 + len2 - 1, the length of the result\n");

    if (rect)
        _table_rect(alg, modn, budget, membytes);
    else
        _table_square(alg, modn, budget, membytes);

    return 0;
}
