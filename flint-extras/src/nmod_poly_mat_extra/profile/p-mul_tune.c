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
     fun        sd_fft_direct | sd_fft_matmul | geometric | multiply

   Options of the table mode, in any order:
     rect         sweep rectangular shapes and unbalanced lengths
                  instead of the square grid
     budget=SECS  give up on a parameter point after about SECS seconds
                  (default 20); the grid is pruned monotonically from
                  each point that exceeds it
     mem=GB       skip a parameter point whose operands and result would
                  exceed about GB gigabytes (default 8)

   A cell reads "-" when the point was pruned by one of the two limits
   above, or when the algorithm does not apply to that modulus.
*/

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

static const slong SQ_DIMS[] = { 2, 4, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256, 384, 512 };
static const slong SQ_LENS[] = { 8, 16, 32, 64, 128, 256, 512, 1024, 2048 };

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

static mulfun _select(const char * alg, ulong modn, slong len)
{
    if (!strcmp(alg, "sd_fft_direct"))  return nmod_poly_mat_mul_sd_fft_direct;
    if (!strcmp(alg, "sd_fft_matmul"))  return nmod_poly_mat_mul_sd_fft_matmul;
    if (!strcmp(alg, "multiply"))       return nmod_poly_mat_multiply;
    if (!strcmp(alg, "geometric"))
        return NMOD_POLY_CAN_USE_GEOMETRIC(modn, len)
               ? nmod_poly_mat_mul_geometric : NULL;
    return NULL;
}

/* bytes of the operands and of the result */
static double _opbytes(slong d1, slong d2, slong d3, slong l1, slong l2)
{
    return 8.0 * ((double) d1 * d2 * l1 + (double) d2 * d3 * l2
                  + (double) d1 * d3 * (l1 + l2 - 1));
}

/* minimum wall time over a few repetitions, after one untimed call which
   sizes the run and warms the scratch buffers that a sequence of
   products would find warm anyway */
static double _time_one(mulfun fun, slong d1, slong d2, slong d3,
                        slong len1, slong len2, ulong modn)
{
    flint_rand_t state;
    nmod_poly_mat_t A, B, C;
    timeit_t timer;
    double best;
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

    reps = (slong) (0.5 / (best > 1e-6 ? best : 1e-6));
    reps = FLINT_MAX(reps, 1);
    reps = FLINT_MIN(reps, 20);

    for (i = 0; i < reps; i++)
    {
        double t;
        timeit_start_us(timer);
        fun(C, A, B);
        timeit_stop_us(timer);
        t = 1e-6 * timer->wall;
        best = FLINT_MIN(best, t);
    }

    nmod_poly_mat_clear(A);
    nmod_poly_mat_clear(B);
    nmod_poly_mat_clear(C);
    flint_rand_clear(state);

    return best;
}

static void _table_square(const char * alg, ulong modn, double budget, double membytes)
{
    /* blown[j] is the smallest row index at which column j exceeded a
       limit; a point is skipped as soon as one at least as small in both
       parameters has, since neither is cheaper */
    slong blown[NSQ_LENS];
    slong i, j;

    for (j = 0; j < NSQ_LENS; j++)
        blown[j] = NSQ_DIMS;

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
            double t;

            if (fun == NULL || i >= blown[j]
                || _opbytes(d, d, d, l, l) > membytes)
            {
                flint_printf(" %10s", "-");
                continue;
            }

            t = _time_one(fun, d, d, d, l, l, modn);
            flint_printf(" %10.3e", t);
            fflush(stdout);

            if (t > budget)
            {
                slong jj;
                for (jj = j; jj < NSQ_LENS; jj++)
                    blown[jj] = FLINT_MIN(blown[jj], i);
            }
        }
        flint_printf("\n");
    }
}

static void _table_rect(const char * alg, ulong modn, double budget, double membytes)
{
    const slong ncols = NRECT_BASELENS * NRECT_LENS;
    slong b, s, c;

    flint_printf("%-18s", "m x k x n \\ len");
    for (c = 0; c < ncols; c++)
    {
        const slong base = RECT_BASELENS[c / NRECT_LENS];
        const slong l = c % NRECT_LENS;
        const slong len1 = FLINT_MAX(base * RECT_LENS[l][0] / 8, 1);
        const slong len2 = FLINT_MAX(base * RECT_LENS[l][1] / 8, 1);
        char lab[32];

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
            int stop = 0;
            char shape[32];

            flint_sprintf(shape, "%wd x %wd x %wd", d1, d2, d3);
            flint_printf("%-18s", shape);

            for (c = 0; c < ncols; c++)
            {
                const slong base = RECT_BASELENS[c / NRECT_LENS];
                const slong l = c % NRECT_LENS;
                const slong len1 = FLINT_MAX(base * RECT_LENS[l][0] / 8, 1);
                const slong len2 = FLINT_MAX(base * RECT_LENS[l][1] / 8, 1);
                mulfun fun = _select(alg, modn, len1 + len2 - 1);
                double t;

                if (stop || fun == NULL
                    || _opbytes(d1, d2, d3, len1, len2) > membytes)
                {
                    flint_printf(" %10s", "-");
                    continue;
                }

                t = _time_one(fun, d1, d2, d3, len1, len2, modn);
                flint_printf(" %10.3e", t);
                fflush(stdout);
                if (t > budget)
                    stop = 1;
            }
            flint_printf("\n");
        }
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
        flint_printf("   fun:   sd_fft_direct | sd_fft_matmul | geometric | multiply\n");
        flint_printf("   opts:  rect | budget=SECS | mem=GB\n");
        return 0;
    }

    modn = strtoul(argv[1], NULL, 10);
    if (modn == 0)
        modn = FFT_PRIME;
    nthreads = atol(argv[2]);
    alg = argv[3];
    flint_set_num_threads(nthreads);

    if (_select(alg, modn, 2) == NULL
        && strcmp(alg, "geometric"))   /* geometric may just not apply */
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
            flint_printf("%.3e\n", _time_one(fun, d1, d2, d3, len1, len2, modn));
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
