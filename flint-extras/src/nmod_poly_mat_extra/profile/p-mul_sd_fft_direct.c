/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

/* One timing per invocation, so that the measurement is not perturbed by
   the allocations of the other algorithms: the transforms of this one
   reach hundreds of megabytes, and running several variants in a single
   process measures the memory subsystem as much as the algorithms.
   Drive it from a shell loop to build a table.

   Unlike p-mul, the modulus is given explicitly, which is what allows
   the single-prime path (a modulus that is itself an FFT prime, i.e. of
   at most 50 bits with 2-adicity at least the transform depth) to be
   compared against the 3-prime CRT path.

   Usage: p-mul_sd_fft_direct alg dim1 dim2 dim3 len1 len2 modulus [threads]
     alg: sd_fft_direct | geometric | waksman | mul | multiply
   Prints the minimum wall time over a few repetitions, or "-" when the
   algorithm does not apply to that modulus.
*/

#include <stdlib.h>
#include <string.h>

#include <flint/nmod_poly_mat.h>
#include <flint/profiler.h>
#include <flint/ulong_extras.h>

#include "nmod_poly_extra.h"
#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_multiply.h"

typedef void (*mulfun)(nmod_poly_mat_t, const nmod_poly_mat_t, const nmod_poly_mat_t);

int main(int argc, char ** argv)
{
    if (argc < 8)
    {
        flint_printf("Usage: %s alg dim1 dim2 dim3 len1 len2 modulus [threads]\n", argv[0]);
        flint_printf("   alg: sd_fft_direct | geometric | waksman | mul | multiply\n");
        return 0;
    }

    const char * alg = argv[1];
    const slong d1 = atol(argv[2]), d2 = atol(argv[3]), d3 = atol(argv[4]);
    const slong len1 = atol(argv[5]), len2 = atol(argv[6]);
    const ulong modn = strtoul(argv[7], NULL, 10);
    const slong len = len1 + len2 - 1;

    if (argc > 8)
        flint_set_num_threads(atol(argv[8]));

    mulfun fun;
    if (!strcmp(alg, "sd_fft_direct"))     fun = nmod_poly_mat_mul_sd_fft_direct;
    else if (!strcmp(alg, "mul"))       fun = nmod_poly_mat_mul;
    else if (!strcmp(alg, "multiply"))  fun = nmod_poly_mat_multiply;
    else if (!strcmp(alg, "geometric"))
    {
        if (!NMOD_POLY_CAN_USE_GEOMETRIC(modn, len)) { flint_printf("-\n"); return 0; }
        fun = nmod_poly_mat_mul_geometric;
    }
    else if (!strcmp(alg, "waksman"))
    {
        if (!NMOD_POLY_MAT_CAN_USE_WAKSMAN(modn)) { flint_printf("-\n"); return 0; }
        fun = nmod_poly_mat_mul_waksman;
    }
    else { flint_printf("unknown algorithm %s\n", alg); return 1; }

    flint_rand_t state;
    flint_rand_init(state);
    flint_rand_set_seed(state, 1234, 5678);

    nmod_poly_mat_t A, B, C;
    nmod_poly_mat_init(A, d1, d2, modn);
    nmod_poly_mat_init(B, d2, d3, modn);
    nmod_poly_mat_init(C, d1, d3, modn);
    nmod_poly_mat_rand(A, state, len1);
    nmod_poly_mat_rand(B, state, len2);

    /* one untimed call: it sizes the run, and warms the retained
       scratch buffer that repeated products would find warm anyway */
    timeit_t timer;
    double best;

    timeit_start_us(timer);
    fun(C, A, B);
    timeit_stop_us(timer);
    best = 1e-6 * timer->wall;

    slong reps = (slong) (0.5 / (best > 1e-6 ? best : 1e-6));
    reps = FLINT_MAX(reps, 1);
    reps = FLINT_MIN(reps, 50);

    for (slong i = 0; i < reps; i++)
    {
        double t;
        timeit_start_us(timer);
        fun(C, A, B);
        timeit_stop_us(timer);
        t = 1e-6 * timer->wall;
        best = FLINT_MIN(best, t);
    }

    flint_printf("%.3e\n", best);

    nmod_poly_mat_clear(A);
    nmod_poly_mat_clear(B);
    nmod_poly_mat_clear(C);
    flint_rand_clear(state);
    return 0;
}
