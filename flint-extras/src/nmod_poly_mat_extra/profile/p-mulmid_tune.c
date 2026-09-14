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
   Timings of the polynomial matrix middle product variants that
   nmod_poly_mat_mulmid chooses between, laid out so that the crossovers
   between them can be read off. The companion of p-mul_tune, with the
   same conventions (one algorithm per run, a time budget enforced before
   the fact, a memory limit; see there and tune_impl.h).

   A point is the product of a dim1 x dim2 matrix A of length len1 by a
   dim2 x dim3 matrix B of length len2, of which the coefficients nlo,
   ..., nhi-1 are wanted.

   Usage:
     p-mulmid_tune prime nthreads fun [opts]
         table of timings of one algorithm over a grid of parameters
     p-mulmid_tune prime nthreads [fun1,fun2,...] [opts]
         one line per point of the grid, the parameters followed by the
         timings of the listed algorithms side by side; "[]" lists them
         all, and a list of one name still gives this output
     p-mulmid_tune prime nthreads fun dim1 dim2 dim3 len1 len2 nlo nhi
         one timing

     prime      the modulus; 0 selects one of FLINT's 50-bit FFT primes
     nthreads   number of threads
     fun        geometric | naive | mulmid | multiply

                geometric   _nmod_poly_mat_mulmid_geometric, which needs
                            len1 <= nlo+1 or len2 <= nlo+1 and a field
                            large enough for a progression of nhi points
                naive       the full product, shifted and truncated
                mulmid      nmod_poly_mat_mulmid, the dispatcher
                multiply    not a middle product: nmod_poly_mat_multiply on
                            the product that the middle product is the
                            transpose of, i.e. A (length nlo+1) times a
                            dim2 x dim3 matrix of length nhi-nlo, whose
                            result has length nhi. By the transposition
                            principle the two have the same cost, up to
                            the constant of the algorithm; this column
                            shows how far the middle product is from it.

   The table has one row per dimension of a square product and one column
   per base length L; the point of a column is

     nlo = lo*L/8,  nhi = hi*L/8,  len1 = nlo + 1,  len2 = nhi,

   with lo = 8 and hi = 16 by default -- the balanced middle product, A of
   length L+1 against B of length 2L for L coefficients -- and the column
   is labelled by nhi.

   Options of the table modes, in any order:
     lo=N hi=N    nlo and nhi in eighths of the base length, as above
     budget=SECS  spend at most about SECS seconds on any one parameter
                  point (default 20)
     mem=GB       skip a parameter point whose operands and result would
                  exceed about GB gigabytes (default 8)

   A cell reads "-" when the algorithm does not apply to that modulus and
   point, or when the point was skipped by one of the two limits above.
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
		./build/nmod_poly_mat_extra/profile/p-mulmid_tune $prime $nthreads [] budget=10
		echo "--------------------------------------------------"
	done;
done
*/

#include "tune_impl.h"

#include "nmod_poly_extra.h"
#include "nmod_poly_mat_multiply.h"

/* base lengths of the square grid */
static const slong SQ_LENS[] = { 2, 6, 8, 16, 32, 64, 128, 256, 512, 1024, 2048 };
#define NSQ_LENS (slong)(sizeof(SQ_LENS) / sizeof(SQ_LENS[0]))

/* nlo and nhi of a column, in eighths of its base length */
static slong LO8 = 8, HI8 = 16;

/* the routines under test, in the calling convention of tune_impl.h */
static void _run_geometric(nmod_poly_mat_t C, const nmod_poly_mat_t A,
                           const nmod_poly_mat_t B, const tune_point_struct * P)
{
    _nmod_poly_mat_mulmid_geometric(C, A, P->len1, B, P->len2, P->nlo, P->nhi);
}

static void _run_naive(nmod_poly_mat_t C, const nmod_poly_mat_t A,
                       const nmod_poly_mat_t B, const tune_point_struct * P)
{
    _nmod_poly_mat_mulmid_naive(C, A, P->len1, B, P->len2, P->nlo, P->nhi);
}

static void _run_mulmid(nmod_poly_mat_t C, const nmod_poly_mat_t A,
                        const nmod_poly_mat_t B, const tune_point_struct * P)
{
    nmod_poly_mat_mulmid(C, A, B, P->nlo, P->nhi);
}

static void _run_multiply(nmod_poly_mat_t C, const nmod_poly_mat_t A,
                          const nmod_poly_mat_t B, const tune_point_struct * P)
{
    (void) P;
    nmod_poly_mat_multiply(C, A, B);
}

static const char * ALGS[] = { "geometric", "naive", "mulmid", "multiply" };
#define NALGS (slong)(sizeof(ALGS) / sizeof(ALGS[0]))

/*
    The routine for alg at the middle product P, NULL when it does not
    apply. For multiply, P is replaced by the transposed product: the
    same A, a right operand of length nhi-nlo, and the whole result.
*/
static tune_fun _select(const char * alg, ulong modn, tune_point_struct * P)
{
    if (P->nlo >= P->nhi || P->len1 <= 0 || P->len2 <= 0)
        return NULL;

    if (!strcmp(alg, "geometric"))
        return (NMOD_POLY_CAN_USE_GEOMETRIC(modn, P->nhi)
                && (P->len1 <= P->nlo + 1 || P->len2 <= P->nlo + 1))
               ? _run_geometric : NULL;
    if (!strcmp(alg, "naive"))   return _run_naive;
    if (!strcmp(alg, "mulmid"))  return _run_mulmid;
    if (!strcmp(alg, "multiply"))
    {
        P->len2 = P->nhi - P->nlo;
        P->nlo = 0;
        P->nhi = P->len1 + P->len2 - 1;
        return _run_multiply;
    }
    return NULL;
}

/* the square grid: dimension i, base length SQ_LENS[j] */
static slong SQ_X[NSQ_LENS];

static void _sq_point(slong i, slong j, tune_point_struct * P)
{
    const slong d = TUNE_SQ_DIMS[i], L = SQ_LENS[j];

    P->d1 = d; P->d2 = d; P->d3 = d;
    P->nlo = LO8 * L / 8;
    P->nhi = HI8 * L / 8;
    P->len1 = P->nlo + 1;
    P->len2 = P->nhi;
}
static void _sq_rowlabel(slong i, char * buf) { flint_sprintf(buf, "%wd", TUNE_SQ_DIMS[i]); }
static void _sq_collabel(slong j, char * buf) { flint_sprintf(buf, "%wd", SQ_X[j]); }

int main(int argc, char ** argv)
{
    ulong modn;
    slong nthreads, nalgs;
    const char * algs[NALGS];
    double budget = 20.0, membytes = 8.0 * (1 << 30);
    int listed, i, j;

    if (argc < 4)
    {
        flint_printf("Usage: %s prime nthreads fun [opts]\n", argv[0]);
        flint_printf("       %s prime nthreads [fun1,fun2,...] [opts]\n", argv[0]);
        flint_printf("       %s prime nthreads fun dim1 dim2 dim3 len1 len2 nlo nhi\n", argv[0]);
        flint_printf("   prime: the modulus; 0 selects a 50-bit FFT prime\n");
        flint_printf("   fun:   geometric | naive | mulmid | multiply\n");
        flint_printf("          a list gives one line per point with the algorithms\n");
        flint_printf("          side by side; [] lists them all\n");
        flint_printf("   opts:  lo=N | hi=N | budget=SECS | mem=GB\n");
        flint_printf("          (nlo = lo*L/8 and nhi = hi*L/8 for a column of base length L;\n");
        flint_printf("           default lo=8 hi=16, the balanced middle product)\n");
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

    if (!listed && argc == 11)   /* one timing */
    {
        tune_point_struct P;
        tune_fun fun;

        P.d1 = atol(argv[4]); P.d2 = atol(argv[5]); P.d3 = atol(argv[6]);
        P.len1 = atol(argv[7]); P.len2 = atol(argv[8]);
        P.nlo = atol(argv[9]); P.nhi = atol(argv[10]);
        fun = _select(algs[0], modn, &P);

        if (fun == NULL)
            flint_printf("-\n");
        else
            flint_printf("%.3e\n", tune_time_one(fun, &P, modn, DBL_MAX));
        return 0;
    }

    for ( ; i < argc; i++)
    {
        if (!strncmp(argv[i], "lo=", 3))
            LO8 = atol(argv[i] + 3);
        else if (!strncmp(argv[i], "hi=", 3))
            HI8 = atol(argv[i] + 3);
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

    if (LO8 < 0 || HI8 <= LO8)
    {
        flint_printf("need 0 <= lo < hi\n");
        return 1;
    }

    tune_print_header(algs, nalgs, modn, nthreads, budget, membytes);
    flint_printf("# nlo = nhi * %wd / %wd, len1 = nlo + 1, len2 = nhi, "
                 "result length nhi - nlo; columns of a table are nhi\n", LO8, HI8);
    flint_printf("# multiply is the transposed product: A of length nlo + 1 "
                 "times a matrix of length nhi - nlo\n");

    {
        tune_grid_struct G = { TUNE_NSQ_DIMS, NSQ_LENS, 1,
                               TUNE_SQ_DIMS, SQ_X, "dim\\nhi",
                               _sq_rowlabel, _sq_collabel, _sq_point, _select };
        for (j = 0; j < NSQ_LENS; j++)
            SQ_X[j] = HI8 * SQ_LENS[j] / 8;
        tune_run(&G, algs, nalgs, modn, nthreads, budget, membytes, listed, 1);
    }

    return 0;
}
