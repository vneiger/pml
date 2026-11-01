/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>  // for atoi, strtoul

#include <flint/flint.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_poly_mat.h>
#include <flint/profiler.h>
#include <flint/ulong_extras.h>

#include "nmod_poly_mat_inverse.h"

/* times the truncated inverse of a random dim x dim matrix of degree < deg,
 * at order `order`, with the general interface, the basecase algorithm,
 * and Newton iteration (with the given cutoff, or the default one); prints
 * the times in seconds */

#define TIME_INV(call)                                          \
    do {                                                        \
        double best = 0.0;                                      \
        for (slong rep = 0; rep < 3; rep++)                     \
        {                                                       \
            slong reps = 1;                                     \
            double el;                                          \
            timeit_t t;                                         \
            while (1)                                           \
            {                                                   \
                timeit_start(t);                                \
                for (slong r = 0; r < reps; r++)                \
                    call;                                       \
                timeit_stop(t);                                 \
                if (t->wall >= 50)                              \
                {                                               \
                    el = t->wall / 1000.0 / reps;               \
                    break;                                      \
                }                                               \
                reps *= 2;                                      \
            }                                                   \
            if (rep == 0 || el < best)                          \
                best = el;                                      \
        }                                                       \
        flint_printf("  %.3e", best);                           \
    } while (0)

int main(int argc, char ** argv)
{
    if (argc < 5 || argc > 6)
    {
        flint_printf("Usage: `%s [nbits] [dim] [deg] [order] [cutoff]`\n", argv[0]);
        flint_printf("   - nbits: number of bits in [2..64] for the modulus, chosen as nextprime(2**(nbits-1))\n");
        flint_printf("   - dim: the input matrix is dim x dim\n");
        flint_printf("   - deg: the input matrix is random of degree < deg\n");
        flint_printf("   - order: truncation order\n");
        flint_printf("   - cutoff: optional, cutoff for Newton iteration (default: 0, the default cutoff)\n");
        flint_printf("Prints the times (in seconds) of: general interface, basecase, Newton iteration\n");
        return 0;
    }

    const ulong b = strtoul(argv[1], NULL, 10);
    const slong dim = atoi(argv[2]);
    const slong deg = atoi(argv[3]);
    const slong order = atoi(argv[4]);
    const slong cutoff = (argc >= 6) ? atoi(argv[5]) : 0;
    const ulong n = n_nextprime(UWORD(1) << (b-1), 0);

    flint_rand_t state;
    flint_rand_init(state);
    flint_rand_set_seed(state, 1234, 5678);

    nmod_poly_mat_t pmat, imat;
    nmod_poly_mat_init(pmat, dim, dim, n);
    nmod_poly_mat_init(imat, dim, dim, n);
    nmod_mat_t inv0;
    nmod_mat_init(inv0, dim, dim, n);
    do
    {
        nmod_poly_mat_rand(pmat, state, deg);
        nmod_poly_mat_get_coeff_mat(inv0, pmat, 0);
    } while (! nmod_mat_inv(inv0, inv0));

    flint_printf("nbits %wu, dim %wd, deg %wd, order %wd, cutoff %wd (default %wd):",
                 b, dim, deg, order, cutoff, _nmod_poly_mat_inv_series_cutoff(dim));
    TIME_INV(nmod_poly_mat_inv_series(imat, pmat, order));
    TIME_INV(nmod_poly_mat_inv_series_basecase(imat, pmat, order));
    TIME_INV(nmod_poly_mat_inv_series_newton(imat, pmat, order, cutoff));
    flint_printf("\n");

    nmod_mat_clear(inv0);
    nmod_poly_mat_clear(pmat);
    nmod_poly_mat_clear(imat);
    flint_rand_clear(state);
    return 0;
}
