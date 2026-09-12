/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

/* Timings for the nmod_mat_poly -> nmod_poly_mat conversion.
 *
 * Unit: nanoseconds per converted word
 *
 * The conversion is a transposition of r*c*order words and moves 8 bytes in
 * and 8 bytes out per word, so a value close to the cost of a plain memcpy of
 * the same volume is the target.
 *
 * Unlike the other direction, the destination rows here are the coefficient
 * arrays of the entries of an nmod_poly_mat, whose alignment is whatever
 * `flint_realloc` returned; the source rows, the coefficients of the
 * nmod_mat_poly, are the 64-byte aligned ones. */

#include <stdlib.h>  // for atol, atoi
#include <time.h>    // for time

#include <flint/flint.h>
#include <flint/ulong_extras.h>
#include <flint/profiler.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>

#include "pml.h"
#include "nmod_mat_poly.h"
#include "nmod_poly_mat_utils.h"
#include "nmod_mat_poly_extra/impl.h"
#include "nmod_poly_mat_extra/impl.h"

/* Straightforward entry-by-entry conversion: the reference that the blocked
   implementation is checked against. */
static void _set_trunc_from_mat_poly_naive(nmod_poly_mat_t pmat,
                                           const nmod_mat_poly_t matp,
                                           slong order)
{
    if (order > matp->length)
        order = matp->length;

    for (slong i = 0; i < pmat->r; i++)
        for (slong j = 0; j < pmat->c; j++)
            nmod_poly_fit_length(nmod_poly_mat_entry(pmat, i, j), order);

    for (slong k = 0; k < order; k++)
        for (slong i = 0; i < pmat->r; i++)
            for (slong j = 0; j < pmat->c; j++)
                nmod_poly_mat_entry(pmat, i, j)->coeffs[k] = nmod_mat_poly_entry(matp, k, i, j);

    for (slong i = 0; i < pmat->r; i++)
        for (slong j = 0; j < pmat->c; j++)
        {
            _nmod_poly_set_length(nmod_poly_mat_entry(pmat, i, j), order);
            _nmod_poly_normalise(nmod_poly_mat_entry(pmat, i, j));
        }
}

#define NFUNS 8

static const char * description[NFUNS] = {
    "#0  --> naive (entry by entry, local reference)           ",
    "#1  --> blocked 8x8, scalar,          entry-major (st seq)",
    "#2  --> blocked 8x8, scalar,          coeff-major (ld seq)",
    "#3  --> blocked 4x4, machine vectors, entry-major (st seq)",
    "#4  --> blocked 4x4, machine vectors, coeff-major (ld seq)",
    "#5  --> blocked 8x8, AVX-512,         entry-major (st seq)",
    "#6  --> blocked 8x8, AVX-512,         coeff-major (ld seq)",
    "#7  --> default kernel and schedule                       ",
};

/* (kernel, schedule) of function number f >= 1; -1 means "default".
   `dmaj == 1` sweeps the destination rows, i.e. the output polynomials. */
static const int fun_kern[NFUNS] = {0, NMOD_MAT_POLY_CONV_SCALAR, NMOD_MAT_POLY_CONV_SCALAR,
                                       NMOD_MAT_POLY_CONV_VEC4, NMOD_MAT_POLY_CONV_VEC4,
                                       NMOD_MAT_POLY_CONV_VEC8, NMOD_MAT_POLY_CONV_VEC8,
                                       -1};
static const int fun_dmaj[NFUNS] = {0, 1, 0, 1, 0, 1, 0, -1};

/* random matrix polynomial; if `ragged`, some coefficients are zeroed, so the
   output entries come out with unequal lengths */
static void _rand_matp(nmod_mat_poly_t matp, flint_rand_t state, slong len, int ragged)
{
    nmod_mat_poly_rand(matp, state, len);
    if (ragged)
        for (slong k = 0; k < matp->length; k++)
            for (slong i = 0; i < matp->r; i++)
                for (slong j = 0; j < matp->c; j++)
                    if (n_randint(state, 3) == 0)
                        nmod_mat_poly_entry(matp, k, i, j) = UWORD(0);
}

static double time_fun(slong fun_nb, slong dim1, slong dim2, slong len,
                       ulong modn, int ragged, flint_rand_t state)
{
    nmod_mat_poly_t matp;
    nmod_mat_poly_init(matp, dim1, dim2, modn);
    _rand_matp(matp, state, len, ragged);

    nmod_poly_mat_t pmat;
    nmod_poly_mat_init(pmat, dim1, dim2, modn);

    double FLINT_SET_BUT_UNUSED(tcpu), twall;

    if (fun_nb == 0)
    {
        TIMEIT_START;
        _set_trunc_from_mat_poly_naive(pmat, matp, len);
        TIMEIT_STOP_VALUES(tcpu, twall);
    }
    else
    {
        TIMEIT_START;
        _nmod_poly_mat_set_trunc_from_mat_poly(pmat, matp, len,
                                               fun_kern[fun_nb], fun_dmaj[fun_nb]);
        TIMEIT_STOP_VALUES(tcpu, twall);
    }

    nmod_poly_mat_clear(pmat);
    nmod_mat_poly_clear(matp);

    /* nanoseconds per converted word */
    return 1e9 * twall / ((double) dim1 * (double) dim2 * (double) len);
}

int main(int argc, char ** argv)
{
    flint_rand_t state;
    flint_rand_init(state);
    flint_rand_set_seed(state, time(NULL), time(NULL) + 129384125L);

    flint_printf("build: 4x4 machine-vectors kernel: %d, 8x8 AVX-512 kernel: %d,"
                 " default schedule: %s\n",
#if PML_HAVE_MACHINE_VECTORS
                 1,
#else
                 0,
#endif
#if PML_HAVE_AVX512
                 1,
#else
                 0,
#endif
                 PML_CONV_DST_MAJOR ? "entry-major" : "coeff-major"
                );
    flint_printf("(unavailable kernels silently fall back, columns then repeat)\n");

    if (argc == 1)
    {
        flint_printf("Usage: `%s [nbits] [fun] [dim1] [dim2] [len] [opt:ragged]`\n", argv[0]);
        flint_printf("   No argument runs a default grid with all functions.\n");
        flint_printf("   - nbits: number of bits in (1..64] for the modulus, nextprime(2**(nbits-1))\n");
        flint_printf("   - fun: id of the timed function, -1 for all (see below)\n");
        flint_printf("   - dim1, dim2: the input matrix polynomial is dim1 x dim2\n");
        flint_printf("   - len: the input has length len (order of the conversion)\n");
        flint_printf("   - ragged: optional, if nonzero some input coefficients are zeroed\n");
        flint_printf("\nAvailable functions:\n");
        for (slong j = 0; j < NFUNS; j++)
            flint_printf("   %s\n", description[j]);
        flint_printf("\nRunning default grid (ns per converted word, 60-bit modulus):\n");

        const ulong modn = n_nextprime(UWORD(1) << 59, 0);
        const slong dims[] = {2, 4, 8, 16, 32, 64, 128};
        const slong lens[] = {8, 32, 128, 512, 2048};

        flint_printf("%5s %6s %8s", "dim", "len", "MB");
        for (slong j = 0; j < NFUNS; j++)
            flint_printf(" %8s%-2wd", "fun#", j);
        flint_printf("   speedup\n");

        for (slong di = 0; di < 7; di++)
            for (slong li = 0; li < 5; li++)
            {
                const slong d = dims[di], len = lens[li];
                if ((double) d * d * len > 2.0e7)
                    continue;
                flint_printf("%5wd %6wd %8.2f", d, len, (double) d * d * len * 8 / 1048576.0);
                double t[NFUNS], best = 1e30;
                for (slong j = 0; j < NFUNS; j++)
                {
                    t[j] = time_fun(j, d, d, len, modn, 0, state);
                    if (j > 0 && t[j] < best)
                        best = t[j];
                    flint_printf(" %10.3f", t[j]);
                }
                flint_printf("   %5.2fx\n", t[0] / best);
                fflush(stdout);
            }

        flint_rand_clear(state);
        return 0;
    }

    if (argc >= 6)
    {
        const slong b = atol(argv[1]);
        const slong ifun = atol(argv[2]);
        const slong dim1 = atol(argv[3]);
        const slong dim2 = atol(argv[4]);
        const slong len = atol(argv[5]);
        const int ragged = (argc >= 7) ? atoi(argv[6]) : 0;

        const ulong modn = n_nextprime(UWORD(1) << (b - 1), 0);

        flint_printf("Available functions:\n");
        for (slong j = 0; j < NFUNS; j++)
            flint_printf("   %s\n", description[j]);

        flint_printf("%-5s%-6s%-6s%-8s", "bits", "dim1", "dim2", "len");
        if (ifun >= 0)
            flint_printf("fun#%-6wd\n", ifun);
        else
        {
            for (slong j = 0; j < NFUNS; j++)
                flint_printf("fun#%-6wd", j);
            flint_printf("\n");
        }

        flint_printf("%-5wd%-6wd%-6wd%-8wd", b, dim1, dim2, len);
        if (ifun >= 0)
            flint_printf("%-10.3f", time_fun(ifun, dim1, dim2, len, modn, ragged, state));
        else
            for (slong j = 0; j < NFUNS; j++)
                flint_printf("%-10.3f", time_fun(j, dim1, dim2, len, modn, ragged, state));
        flint_printf("\n");
    }

    flint_rand_clear(state);
    return 0;
}
