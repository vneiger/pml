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
#include <string.h>  // for strstr

#include <flint/flint.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>
#include <flint/profiler.h>
#include <flint/ulong_extras.h>

#include "nmod_poly_mat_linsolve.h"

/* times the solution of u A = den b over the polynomials, for a random
 * dim x dim matrix A of degree < deg (with A(0) invertible) and a random
 * 1 x dim vector b of degree < degb (by default degb = deg), by the
 * functions below; prints the times in seconds, the degree of den, and
 * "MISMATCH" if the outputs differ (they should not: A is nonsingular, and
 * all functions give the solution in lowest terms with den monic) */

typedef int (*linsolve_fun_t)(nmod_poly_mat_t u, nmod_poly_t den,
                              const nmod_poly_mat_t A, const nmod_poly_mat_t b);

static int _ls_kernel(nmod_poly_mat_t u, nmod_poly_t den,
                      const nmod_poly_mat_t A, const nmod_poly_mat_t b)
{
    return nmod_poly_mat_linsolve_via_kernel(u, den, A, b);
}

static int _ls_series_rr(nmod_poly_mat_t u, nmod_poly_t den,
                         const nmod_poly_mat_t A, const nmod_poly_mat_t b)
{
    return nmod_poly_mat_linsolve_via_series(u, den, A, b, LINSOLVE_RECONSTRUCTION);
}

static int _ls_series_det(nmod_poly_mat_t u, nmod_poly_t den,
                          const nmod_poly_mat_t A, const nmod_poly_mat_t b)
{
    return nmod_poly_mat_linsolve_via_series(u, den, A, b, LINSOLVE_DETERMINANT);
}

typedef struct
{
    const char * name;
    linsolve_fun_t fun;
}
linsolve_method_t;

static const linsolve_method_t methods[] =
{
    {"kernel", _ls_kernel},
    {"series_rr", _ls_series_rr},
    {"series_det", _ls_series_det},
};

#define NMETHODS ((slong) (sizeof(methods) / sizeof(methods[0])))

/* best of 3, each repeated until it takes at least 50 ms */
static double _time_method(linsolve_fun_t fun, nmod_poly_mat_t u, nmod_poly_t den,
                           const nmod_poly_mat_t A, const nmod_poly_mat_t b)
{
    double best = 0.0;
    for (slong rep = 0; rep < 3; rep++)
    {
        slong reps = 1;
        double el;
        timeit_t t;
        while (1)
        {
            timeit_start(t);
            for (slong r = 0; r < reps; r++)
                fun(u, den, A, b);
            timeit_stop(t);
            if (t->wall >= 50)
            {
                el = t->wall / 1000.0 / reps;
                break;
            }
            reps *= 2;
        }
        if (rep == 0 || el < best)
            best = el;
    }
    return best;
}

int main(int argc, char ** argv)
{
    if (argc < 4 || argc > 6)
    {
        flint_printf("Usage: `%s [nbits] [dim] [deg] [degb] [methods]`\n", argv[0]);
        flint_printf("   - nbits: number of bits in [2..64] for the modulus, chosen as nextprime(2**(nbits-1))\n");
        flint_printf("   - dim: A is dim x dim, b is 1 x dim\n");
        flint_printf("   - deg: A is random of degree < deg\n");
        flint_printf("   - degb: optional, b is random of degree < degb (default deg)\n");
        flint_printf("   - methods: optional, comma-separated names of the methods to time (default: all)\n");
        flint_printf("Methods:");
        for (slong k = 0; k < NMETHODS; k++)
            flint_printf(" %s", methods[k].name);
        flint_printf("\nPrints the times in seconds, deg(den), and MISMATCH if the outputs differ\n");
        return 0;
    }

    const ulong bits = strtoul(argv[1], NULL, 10);
    const slong dim = atoi(argv[2]);
    const slong deg = atoi(argv[3]);
    const slong degb = (argc >= 5) ? atoi(argv[4]) : deg;
    const char * which = (argc >= 6) ? argv[5] : NULL;
    const ulong n = n_nextprime(UWORD(1) << (bits - 1), 0);

    flint_rand_t state;
    flint_rand_init(state);
    flint_rand_set_seed(state, 1234, 5678);

    nmod_poly_mat_t A, b, u, u1;
    nmod_poly_t den, den1;
    nmod_poly_mat_init(A, dim, dim, n);
    nmod_poly_mat_init(b, 1, dim, n);
    nmod_poly_mat_init(u, 1, dim, n);
    nmod_poly_mat_init(u1, 1, dim, n);
    nmod_poly_init(den, n);
    nmod_poly_init(den1, n);
    nmod_mat_t A0;
    nmod_mat_init(A0, dim, dim, n);
    do
    {
        nmod_poly_mat_rand(A, state, deg);
        nmod_poly_mat_get_coeff_mat(A0, A, 0);
    } while (nmod_mat_rank(A0) < dim);
    nmod_poly_mat_rand(b, state, degb);

    flint_printf("nbits %wu, dim %wd, deg %wd, degb %wd:", bits, dim, deg, degb);
    int first = 1, mismatch = 0;
    for (slong k = 0; k < NMETHODS; k++)
    {
        if (which != NULL && strstr(which, methods[k].name) == NULL)
            continue;
        const double t = _time_method(methods[k].fun, u, den, A, b);
        flint_printf("  %s %.3e", methods[k].name, t);
        if (first)
        {
            nmod_poly_mat_set(u1, u);
            nmod_poly_set(den1, den);
            first = 0;
        }
        else if (!nmod_poly_mat_equal(u, u1) || !nmod_poly_equal(den, den1))
            mismatch = 1;
    }
    flint_printf("  deg(den) %wd%s\n", nmod_poly_degree(den1), mismatch ? "  MISMATCH" : "");

    nmod_mat_clear(A0);
    nmod_poly_mat_clear(A);
    nmod_poly_mat_clear(b);
    nmod_poly_mat_clear(u);
    nmod_poly_mat_clear(u1);
    nmod_poly_clear(den);
    nmod_poly_clear(den1);
    flint_rand_clear(state);
    return 0;
}
