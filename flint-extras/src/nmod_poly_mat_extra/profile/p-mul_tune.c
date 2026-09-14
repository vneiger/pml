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

   Usage:
     p-mul_tune prime nthreads fun [opts]
         table of timings of one algorithm over a grid of parameters
     p-mul_tune prime nthreads [fun1,fun2,...] [opts]
         one line per point of the grid, the parameters followed by the
         timings of the listed algorithms side by side; "[]" lists them
         all, and a list of one name still gives this output
     p-mul_tune prime nthreads fun dim1 dim2 dim3 len1 len2
         one timing, for the product of a dim1 x dim2 matrix of length
         len1 by a dim2 x dim3 matrix of length len2

     prime      the modulus; 0 selects one of FLINT's 50-bit FFT primes,
                for which the fft_small variants use a single transform
                modulo the prime itself and no chinese remaindering
     nthreads   number of threads
     fun        sd_fft_direct | sd_fft_matmul | geometric
                | vandermonde1 | vandermonde2 | waksman | multiply

   Options of the table modes, in any order:
     rect         sweep rectangular shapes and unbalanced lengths
                  instead of the square grid
     budget=SECS  spend at most about SECS seconds on any one parameter
                  point (default 20)
     mem=GB       skip a parameter point whose operands and result would
                  exceed about GB gigabytes (default 8)

   A cell reads "-" when the algorithm does not apply to that modulus and
   length, or when the point was skipped by one of the two limits above;
   see tune_impl.h for how the budget is enforced, and for why a table
   per algorithm is the cleaner measurement and the list the easier read.

   The companion p-mulmid_tune does the same for the middle product.
*/

/* bash script for an overall tune (possibly adapt budget=10) */

/*
# 21-bit FFT prime + 50 bit FFT prime + 30 bit prime + 60 bit prime
for nthreads in 1 2 4 8
do
	for prime in 1179649 0 1073741789 1152921504606846883
	do
		echo "--------------------------------------------------"
		echo "nthreads = $nthreads | prime = $prime"
		./build/nmod_poly_mat_extra/profile/p-mul_tune $prime $nthreads [] budget=10
		echo "--------------------------------------------------"
	done;
done
*/

#include "tune_impl.h"

#include "nmod_poly_extra.h"
#include "nmod_poly_mat_multiply.h"

/* operand lengths of the square grid; the columns are the result lengths
   2*len - 1 */
static const slong SQ_LENS[] = { 2, 6, 8, 16, 32, 64, 128, 256, 512, 1024, 2048 };
#define NSQ_LENS (slong)(sizeof(SQ_LENS) / sizeof(SQ_LENS[0]))

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

#define NRECT_SHAPES (slong)(sizeof(RECT_SHAPES) / sizeof(RECT_SHAPES[0]))
#define NRECT_BASES (slong)(sizeof(RECT_BASES) / sizeof(RECT_BASES[0]))
#define NRECT_LENS (slong)(sizeof(RECT_LENS) / sizeof(RECT_LENS[0]))
#define NRECT_BASELENS (slong)(sizeof(RECT_BASELENS) / sizeof(RECT_BASELENS[0]))

/* the routines under test, in the calling convention of tune_impl.h */
#define WRAP(name, call)                                                  \
    static void name(nmod_poly_mat_t C, const nmod_poly_mat_t A,          \
                     const nmod_poly_mat_t B, const tune_point_struct * P) \
    { (void) P; call(C, A, B); }
WRAP(_run_sd_fft_direct, nmod_poly_mat_mul_sd_fft_direct)
WRAP(_run_sd_fft_matmul, nmod_poly_mat_mul_sd_fft_matmul)
WRAP(_run_geometric,     nmod_poly_mat_mul_geometric)
WRAP(_run_vandermonde1,  nmod_poly_mat_mul_vandermonde1)
WRAP(_run_vandermonde2,  nmod_poly_mat_mul_vandermonde2)
WRAP(_run_waksman,       nmod_poly_mat_mul_waksman)
WRAP(_run_multiply,      nmod_poly_mat_multiply)
#undef WRAP

static const char * ALGS[] = { "sd_fft_direct", "sd_fft_matmul", "geometric",
                               "vandermonde1", "vandermonde2", "waksman", "multiply" };
#define NALGS (slong)(sizeof(ALGS) / sizeof(ALGS[0]))

/* the routine for alg at the point P, NULL when it does not apply */
static tune_fun _select(const char * alg, ulong modn, tune_point_struct * P)
{
    const slong len = P->len1 + P->len2 - 1;

    if (!strcmp(alg, "sd_fft_direct"))  return _run_sd_fft_direct;
    if (!strcmp(alg, "sd_fft_matmul"))  return _run_sd_fft_matmul;
    if (!strcmp(alg, "multiply"))       return _run_multiply;
    if (!strcmp(alg, "geometric"))
        return NMOD_POLY_CAN_USE_GEOMETRIC(modn, len) ? _run_geometric : NULL;
    if (!strcmp(alg, "vandermonde1"))
        return NMOD_POLY_CAN_USE_VANDERMONDE1(modn, len) ? _run_vandermonde1 : NULL;
    if (!strcmp(alg, "vandermonde2"))
        return NMOD_POLY_CAN_USE_VANDERMONDE2(modn, len) ? _run_vandermonde2 : NULL;
    if (!strcmp(alg, "waksman"))
        return NMOD_POLY_MAT_CAN_USE_WAKSMAN(modn) ? _run_waksman : NULL;
    return NULL;
}

static void _set_point(tune_point_struct * P, slong d1, slong d2, slong d3,
                       slong len1, slong len2)
{
    P->d1 = d1; P->d2 = d2; P->d3 = d3;
    P->len1 = len1; P->len2 = len2;
    P->nlo = 0; P->nhi = len1 + len2 - 1;
}

/* the square grid: dimension i, both operands of length SQ_LENS[j] */
static slong SQ_X[NSQ_LENS];

static void _sq_point(slong i, slong j, tune_point_struct * P)
{
    _set_point(P, TUNE_SQ_DIMS[i], TUNE_SQ_DIMS[i], TUNE_SQ_DIMS[i],
               SQ_LENS[j], SQ_LENS[j]);
}
static void _sq_rowlabel(slong i, char * buf) { flint_sprintf(buf, "%wd", TUNE_SQ_DIMS[i]); }
static void _sq_collabel(slong j, char * buf) { flint_sprintf(buf, "%wd", SQ_X[j]); }

/* the rectangular grid: row i is shape i % NRECT_SHAPES at base
   RECT_BASES[i / NRECT_SHAPES], column j the lengths RECT_LENS[j %
   NRECT_LENS] in eighths of RECT_BASELENS[j / NRECT_LENS]; the rows
   step over the shapes to reach the same shape at the previous base */
#define NRECT_ROWS (NRECT_BASES * NRECT_SHAPES)
#define NRECT_COLS (NRECT_BASELENS * NRECT_LENS)
static slong RECT_ROWX[NRECT_ROWS], RECT_COLX[NRECT_COLS];

static void _rect_len(slong j, slong * len1, slong * len2)
{
    const slong base = RECT_BASELENS[j / NRECT_LENS];
    const slong l = j % NRECT_LENS;

    *len1 = FLINT_MAX(base * RECT_LENS[l][0] / 8, 1);
    *len2 = FLINT_MAX(base * RECT_LENS[l][1] / 8, 1);
}

static void _rect_point(slong i, slong j, tune_point_struct * P)
{
    const slong d = RECT_BASES[i / NRECT_SHAPES];
    const slong * sh = RECT_SHAPES[i % NRECT_SHAPES];
    slong len1, len2;

    _rect_len(j, &len1, &len2);
    _set_point(P, d * sh[0], d * sh[1], d * sh[2], len1, len2);
}

static void _rect_rowlabel(slong i, char * buf)
{
    tune_point_struct P;
    _rect_point(i, 0, &P);
    flint_sprintf(buf, "%wd x %wd x %wd", P.d1, P.d2, P.d3);
}

static void _rect_collabel(slong j, char * buf)
{
    slong len1, len2;
    _rect_len(j, &len1, &len2);
    flint_sprintf(buf, "%wd+%wd", len1, len2);
}

int main(int argc, char ** argv)
{
    ulong modn;
    slong nthreads, nalgs;
    const char * algs[NALGS];
    double budget = 20.0, membytes = 8.0 * (1 << 30);
    int rect = 0, listed, i, j;

    if (argc < 4)
    {
        flint_printf("Usage: %s prime nthreads fun [opts]\n", argv[0]);
        flint_printf("       %s prime nthreads [fun1,fun2,...] [opts]\n", argv[0]);
        flint_printf("       %s prime nthreads fun dim1 dim2 dim3 len1 len2\n", argv[0]);
        flint_printf("   prime: the modulus; 0 selects a 50-bit FFT prime\n");
        flint_printf("   fun:   sd_fft_direct | sd_fft_matmul | geometric\n");
        flint_printf("          | vandermonde1 | vandermonde2 | waksman | multiply\n");
        flint_printf("          a list gives one line per point with the algorithms\n");
        flint_printf("          side by side; [] lists them all\n");
        flint_printf("   opts:  rect | budget=SECS | mem=GB\n");
        return 0;
    }

    modn = strtoul(argv[1], NULL, 10);
    if (modn == 0)
        modn = TUNE_FFT_PRIME;
    nthreads = atol(argv[2]);
    flint_set_num_threads(nthreads);

    i = tune_parse_algs(argc, argv, 3, ALGS, NALGS, algs, &nalgs, &listed);
    if (i < 0)
    {
        flint_printf("unknown algorithm in %s\n", argv[3]);
        return 1;
    }

    if (!listed && argc == 9)   /* one timing */
    {
        tune_point_struct P;
        tune_fun fun;

        _set_point(&P, atol(argv[4]), atol(argv[5]), atol(argv[6]),
                   atol(argv[7]), atol(argv[8]));
        fun = _select(algs[0], modn, &P);

        if (fun == NULL)
            flint_printf("-\n");
        else
            flint_printf("%.3e\n", tune_time_one(fun, &P, modn, DBL_MAX));
        return 0;
    }

    for ( ; i < argc; i++)
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

    tune_print_header(algs, nalgs, modn, nthreads, budget, membytes);
    flint_printf("# columns are len1 + len2 - 1, the length of the result\n");

    if (rect)
    {
        tune_grid_struct G = { NRECT_ROWS, NRECT_COLS, NRECT_SHAPES,
                               RECT_ROWX, RECT_COLX, "m x k x n \\ len ",
                               _rect_rowlabel, _rect_collabel, _rect_point, _select };
        for (i = 0; i < NRECT_ROWS; i++)
            RECT_ROWX[i] = RECT_BASES[i / NRECT_SHAPES];
        for (j = 0; j < NRECT_COLS; j++)
        {
            slong len1, len2;
            _rect_len(j, &len1, &len2);
            RECT_COLX[j] = len1 + len2 - 1;
        }
        tune_run(&G, algs, nalgs, modn, nthreads, budget, membytes, listed, 0);
    }
    else
    {
        tune_grid_struct G = { TUNE_NSQ_DIMS, NSQ_LENS, 1,
                               TUNE_SQ_DIMS, SQ_X, "dim\\len",
                               _sq_rowlabel, _sq_collabel, _sq_point, _select };
        for (j = 0; j < NSQ_LENS; j++)
            SQ_X[j] = 2 * SQ_LENS[j] - 1;
        tune_run(&G, algs, nalgs, modn, nthreads, budget, membytes, listed, 0);
    }

    return 0;
}
