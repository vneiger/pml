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
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>
#include <flint/profiler.h>
#include <flint/ulong_extras.h>

#include "nmod_poly_mat_linsolve.h"

/* times the solution of u A = den b, for a random dim x dim matrix A of
 * degree < deg (with A(0) invertible) and a random 1 x dim vector b of
 * degree < deg:
 *  - by nmod_poly_mat_linsolve_via_kernel;
 *  - for reference, the power series solution b A^{-1} mod x^prec by
 *    nmod_poly_mat_solve_series_dixon (default order), at the precision
 *    prec = 2 dim (deg-1) + 1 which suffices to recover u / den by rational
 *    reconstruction (not included);
 *  - optionally, by FLINT's nmod_poly_mat_solve on the transposed system
 *    (fraction-free elimination, den = +-det(A), not reduced).
 * prints the times in seconds, and the degree of den */

#define TIME_SOLVE(call)                                        \
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
    if (argc < 4 || argc > 5)
    {
        flint_printf("Usage: `%s [nbits] [dim] [deg] [flint]`\n", argv[0]);
        flint_printf("   - nbits: number of bits in [2..64] for the modulus, chosen as nextprime(2**(nbits-1))\n");
        flint_printf("   - dim: A is dim x dim, b is 1 x dim\n");
        flint_printf("   - deg: A and b are random of degree < deg\n");
        flint_printf("   - flint: optional, 1 to also time FLINT's nmod_poly_mat_solve (default 0)\n");
        flint_printf("Prints the times (in seconds) of: nmod_poly_mat_linsolve_via_kernel,\n");
        flint_printf("    nmod_poly_mat_solve_series_dixon at precision 2*dim*(deg-1)+1, [nmod_poly_mat_solve]; and deg(den)\n");
        return 0;
    }

    const ulong b = strtoul(argv[1], NULL, 10);
    const slong dim = atoi(argv[2]);
    const slong deg = atoi(argv[3]);
    const int with_flint = (argc == 5) ? atoi(argv[4]) : 0;
    const ulong n = n_nextprime(UWORD(1) << (b-1), 0);
    const slong prec = 2 * dim * FLINT_MAX(deg - 1, 0) + 1;

    flint_rand_t state;
    flint_rand_init(state);
    flint_rand_set_seed(state, 1234, 5678);

    nmod_poly_mat_t A, vb, u, X;
    nmod_poly_t den;
    nmod_poly_mat_init(A, dim, dim, n);
    nmod_poly_mat_init(vb, 1, dim, n);
    nmod_poly_mat_init(u, 1, dim, n);
    nmod_poly_mat_init(X, 1, dim, n);
    nmod_poly_init(den, n);
    nmod_mat_t A0;
    nmod_mat_init(A0, dim, dim, n);
    do
    {
        nmod_poly_mat_rand(A, state, deg);
        nmod_poly_mat_get_coeff_mat(A0, A, 0);
    } while (nmod_mat_rank(A0) < dim);
    nmod_poly_mat_rand(vb, state, deg);

    flint_printf("nbits %wu, dim %wd, deg %wd:", b, dim, deg);
    TIME_SOLVE(nmod_poly_mat_linsolve_via_kernel(u, den, A, vb));
    TIME_SOLVE(nmod_poly_mat_solve_series_dixon(X, A, vb, 0, prec));
    if (with_flint)
    {
        nmod_poly_mat_t At, bt, w;
        nmod_poly_t d;
        nmod_poly_mat_init(At, dim, dim, n);
        nmod_poly_mat_init(bt, dim, 1, n);
        nmod_poly_mat_init(w, dim, 1, n);
        nmod_poly_init(d, n);
        nmod_poly_mat_transpose(At, A);
        nmod_poly_mat_transpose(bt, vb);
        TIME_SOLVE(nmod_poly_mat_solve(w, d, At, bt));
        nmod_poly_mat_clear(At);
        nmod_poly_mat_clear(bt);
        nmod_poly_mat_clear(w);
        nmod_poly_clear(d);
    }
    else
        flint_printf("  %9s", "-");
    flint_printf("  %wd\n", nmod_poly_degree(den));

    nmod_mat_clear(A0);
    nmod_poly_mat_clear(A);
    nmod_poly_mat_clear(vb);
    nmod_poly_mat_clear(u);
    nmod_poly_mat_clear(X);
    nmod_poly_clear(den);
    flint_rand_clear(state);
    return 0;
}
