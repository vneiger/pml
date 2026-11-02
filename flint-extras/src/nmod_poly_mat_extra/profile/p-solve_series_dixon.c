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
#include "nmod_poly_mat_linsolve.h"
#include "nmod_poly_mat_multiply.h"

/* times the power series solution X = B A^{-1} mod x^prec, for a random
 * dim x dim matrix A of degree < deg (with A(0) invertible) and a random
 * rdim x dim matrix B of degree < degB (default: deg):
 *  - by nmod_poly_mat_solve_series_dixon, with the given order (truncated
 *    inverse of A at this order, then about prec/order lifting steps; order
 *    0 means the default, 8 (deg(A) + 1));
 *  - optionally, directly, as the truncated inverse of A at order prec
 *    followed by the product with B modulo x^prec;
 * prints the times in seconds */

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

/* X = B A^{-1} mod x^prec, via the truncated inverse at order prec */
static void _solve_via_inverse(nmod_poly_mat_t X, nmod_poly_mat_t Ainv,
                               const nmod_poly_mat_t A, const nmod_poly_mat_t B,
                               slong prec)
{
    nmod_poly_mat_inv_series(Ainv, A, prec);
    nmod_poly_mat_mulmid(X, B, Ainv, 0, prec);
}

int main(int argc, char ** argv)
{
    if (argc < 7 || argc > 9)
    {
        flint_printf("Usage: `%s [nbits] [dim] [rdim] [deg] [order] [prec] [direct] [degB]`\n", argv[0]);
        flint_printf("   - nbits: number of bits in [2..64] for the modulus, chosen as nextprime(2**(nbits-1))\n");
        flint_printf("   - dim: A is dim x dim\n");
        flint_printf("   - rdim: B is rdim x dim\n");
        flint_printf("   - deg: A is random of degree < deg\n");
        flint_printf("   - order: order of the truncated inverse used in Dixon's iteration (0: default, 8 (deg(A) + 1))\n");
        flint_printf("   - prec: the solution is computed modulo x^prec\n");
        flint_printf("   - direct: optional, 0 to skip the timing of the direct approach (default 1)\n");
        flint_printf("   - degB: optional, B is random of degree < degB (default: deg)\n");
        flint_printf("Prints the times (in seconds) of: nmod_poly_mat_solve_series_dixon, [truncated inverse at order prec + product]\n");
        return 0;
    }

    const ulong b = strtoul(argv[1], NULL, 10);
    const slong dim = atoi(argv[2]);
    const slong rdim = atoi(argv[3]);
    const slong deg = atoi(argv[4]);
    const slong order = atoi(argv[5]);
    const slong prec = atoi(argv[6]);
    const int direct = (argc >= 8) ? atoi(argv[7]) : 1;
    const slong degB = (argc >= 9) ? atoi(argv[8]) : deg;
    const ulong n = n_nextprime(UWORD(1) << (b-1), 0);

    flint_rand_t state;
    flint_rand_init(state);
    flint_rand_set_seed(state, 1234, 5678);

    nmod_poly_mat_t A, B, X, Ainv;
    nmod_poly_mat_init(A, dim, dim, n);
    nmod_poly_mat_init(B, rdim, dim, n);
    nmod_poly_mat_init(X, rdim, dim, n);
    nmod_poly_mat_init(Ainv, dim, dim, n);
    nmod_mat_t A0;
    nmod_mat_init(A0, dim, dim, n);
    do
    {
        nmod_poly_mat_rand(A, state, deg);
        nmod_poly_mat_get_coeff_mat(A0, A, 0);
    } while (nmod_mat_rank(A0) < dim);
    nmod_poly_mat_rand(B, state, degB);

    flint_printf("nbits %wu, dim %wd, rdim %wd, deg %wd, order %wd, prec %wd, degB %wd:",
                 b, dim, rdim, deg, order, prec, degB);
    TIME_SOLVE(nmod_poly_mat_solve_series_dixon(X, A, B, order, prec));
    if (direct)
        TIME_SOLVE(_solve_via_inverse(X, Ainv, A, B, prec));
    else
        flint_printf("  %9s", "-");
    flint_printf("\n");

    nmod_mat_clear(A0);
    nmod_poly_mat_clear(A);
    nmod_poly_mat_clear(B);
    nmod_poly_mat_clear(X);
    nmod_poly_mat_clear(Ainv);
    flint_rand_clear(state);
    return 0;
}
