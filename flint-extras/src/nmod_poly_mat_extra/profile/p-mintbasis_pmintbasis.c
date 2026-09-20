/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

/* Note:
   All timed functions are run on the same instance (E,pts): the points used
   throughout are a geometric progression pts_k = r^{2k},
   with r = nmod_find_root(2d) */

#include <stdlib.h>  // for atoi, strtoul
#include <string.h>  // for strtok

#include <flint/ulong_extras.h>
#include <flint/profiler.h>
#include <flint/nmod.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_poly_mat.h>

#include "nmod_extra.h"  // for nmod_find_root
#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_interpolant.h"

/* Shift types are:
 * 0 -> uniform, shift == (0,...,0)
 * */

/* one of the primes of FLINT's default fft_small context: 50 bits, and
   enough 2-adicity for the transform lengths reachable here; the
   polynomial matrix products inside the interpolant basis algorithms then
   run on a single transform modulo the prime itself */
#define FFT_PRIME UWORD(1108307720798209)

typedef struct
{
    slong rdim;  /* row dimension */
    slong cdim;  /* column dimension */
    slong npts;  /* number of points d */
    slong rank;  /* rank */
    slong stype; /* shift type */
    ulong root;  /* r, such that the points are r^{2k}, k = 0..d-1 */
    ulong modn;  /* modulus */
}
time_args;

/* shift is modified by the functions,
 * so it must be zero-ed before the call */

#define TIME_INTERP(name, call)                                     \
void time_##name(time_args targs, const nmod_mat_struct * E,        \
                 const ulong * pts, double * tcpu, double * twall)  \
{                                                                   \
    const slong rdim = targs.rdim;                                  \
    const slong d = targs.npts;                                     \
    /* const slong rank = targs.rank; */ /* TODO */                 \
    /* const slong stype = targs.stype; */ /* TODO */               \
    const ulong r = targs.root;                                     \
    const ulong n = targs.modn;                                     \
    (void) pts; (void) r;                                           \
                                                                    \
    slong * shift = FLINT_ARRAY_ALLOC(rdim, slong);                 \
    nmod_poly_mat_t P;                                              \
    nmod_poly_mat_init(P, rdim, rdim, n);                           \
                                                                    \
    TIMEIT_START;                                                   \
    for (slong i = 0; i < rdim; i++)                                \
        shift[i] = 0;                                               \
    call;                                                           \
    TIMEIT_STOP_VALUES(*tcpu, *twall);                              \
                                                                    \
    nmod_poly_mat_clear(P);                                         \
    flint_free(shift);                                              \
}

TIME_INTERP(pmintbasis,
            nmod_poly_mat_pmintbasis(P, shift, pts, E, d))
TIME_INTERP(mintbasis,
            nmod_poly_mat_mintbasis(P, shift, pts, E, d))
TIME_INTERP(pmintbasis_geometric,
            nmod_poly_mat_pmintbasis_geometric(P, shift, NULL, E, r, d))
TIME_INTERP(pmintbasis_geometric_auto,
            nmod_poly_mat_pmintbasis_geometric_auto(P, shift, NULL, E, d))

/*-------------------------*/
/*  main                   */
/*-------------------------*/

int main(int argc, char ** argv)
{
    flint_rand_t state;
    flint_rand_init(state);
    /* a fixed seed: the same input across thread counts and across runs */
    flint_rand_set_seed(state, 1234, 5678);

    /* TODO rank */
    /* TODO shift type */

    // bench functions
    const slong nfuns = 4;
    typedef void (*timefun) (time_args, const nmod_mat_struct *, const ulong *, double *, double *);
    const timefun funs[] = {
        time_pmintbasis,                   // 0
        time_mintbasis,                    // 1
        time_pmintbasis_geometric,         // 2
        time_pmintbasis_geometric_auto,    // 3
    };

    const char * description[] = {
        "#0  --> pmintbasis (general points)              ",
        "#1  --> mintbasis (general points)               ",
        "#2  --> pmintbasis_geometric (r given)           ",
        "#3  --> pmintbasis_geometric_auto (r found here) ",
    };

    if (argc < 6 || argc > 7)  // show usage
    {
        flint_printf("Usage: `%s [nbits] [fun] [rdim] [cdim] [npts] [nthreads]`\n", argv[0]);
        flint_printf("   No argument shows this help.\n");
        flint_printf("   - nbits: number of bits in [2..64] for the modulus, chosen as nextprime(2**(nbits-1));\n");
        flint_printf("            0 selects a 50-bit FFT prime; a value above 64 is taken as the modulus itself\n");
        flint_printf("   - fun: id number of the timed function (see below),\n");
        flint_printf("   - rdim, cdim: input matrices E_k are rdim x cdim\n");
        flint_printf("   - npts: number d of points; the input is d random matrices E_k, at the\n");
        flint_printf("            geometric points pts_k = r**(2k), k = 0..d-1 (see this file's header)\n");
        flint_printf("   - nthreads: optional, default 1; a comma-separated list, e.g. 1,2,4,8, times\n");
        flint_printf("            the same input at each count and prints the speedups against the first\n");
        flint_printf("   - rank, shift: [unsupported yet]\n");
        flint_printf("\nAvailable functions:\n");
        for (slong j = 0; j < nfuns; j++)
            flint_printf("   %s\n", description[j]);

        return 0;
    }

    // nbits, fun, rdim, cdim, npts, [nthreads]
    {
        const ulong b     = strtoul(argv[1], NULL, 10);
        const slong ifun  = atoi(argv[2]);
        const slong rdim  = atoi(argv[3]);
        const slong cdim  = atoi(argv[4]);
        const slong d     = atoi(argv[5]);
        const ulong n = (b == 0) ? FFT_PRIME
                      : (b <= 64) ? n_nextprime(UWORD(1) << (b-1), 0) : b;
        const timefun tfun = funs[ifun];

        nmod_t mod;
        nmod_init(&mod, n);

        /* the same r as nmod_poly_mat_pmintbasis_geometric_auto's own; it
           exists as soon as n > 2*d+1, which is precisely that function's
           documented requirement */
        const ulong r = (d > 0) ? nmod_find_root(2 * d, mod) : 1;
        if (r == 0)
        {
            flint_printf("modulus %wu too small for npts %wd (requires modulus > 2*npts+1 = %wd)\n",
                         n, d, 2 * d + 1);
            flint_rand_clear(state);
            return 1;
        }

        /* rdim; cdim; npts; rank; stype; root; modn; */
        time_args targs = {rdim, cdim, d, FLINT_MIN(rdim, cdim), 0, r, n};

        /* the thread counts */
        slong nthreads[64], ncounts = 0;
        if (argc == 7)
        {
            char * tok;
            for (tok = strtok(argv[6], ", "); tok != NULL && ncounts < 64; tok = strtok(NULL, ", "))
                nthreads[ncounts++] = atol(tok);
        }
        if (ncounts == 0)
            nthreads[ncounts++] = 1;

        /* E: a flat array of d random rdim x cdim constant matrices, and
           pts: the d geometric points r**(2k) they are the values at */
        nmod_mat_struct * E = FLINT_ARRAY_ALLOC(FLINT_MAX(d, 1), nmod_mat_struct);
        for (slong k = 0; k < d; k++)
        {
            nmod_mat_init(E + k, rdim, cdim, n);
            nmod_mat_rand(E + k, state);
        }

        ulong * pts = FLINT_ARRAY_ALLOC(FLINT_MAX(d, 1), ulong);
        {
            const ulong rho = nmod_mul(r, r, mod);
            ulong cur = 1;
            for (slong k = 0; k < d; k++)
            {
                pts[k] = cur;
                cur = nmod_mul(cur, rho, mod);
            }
        }

        flint_printf("modulus %wu (%wu bits), fun #%wd, rdim %wd, cdim %wd, npts %wd, r %wu\n",
                     n, FLINT_BIT_COUNT(n), ifun, rdim, cdim, d, r);
        flint_printf("%8s %10s %9s %8s\n", "nthreads", "wall", "cpu/wall", "speedup");

        double twall1 = 0.0;
        for (slong c = 0; c < ncounts; c++)
        {
            double tcpu, twall;

            flint_set_num_threads(nthreads[c]);
            tfun(targs, E, pts, &tcpu, &twall);
            if (c == 0)
                twall1 = twall;

            flint_printf("%8wd %10.2e %9.2f %8.2f\n", nthreads[c], twall,
                         twall > 0 ? tcpu / twall : 0.0,
                         twall > 0 ? twall1 / twall : 0.0);
        }

        for (slong k = 0; k < d; k++)
            nmod_mat_clear(E + k);
        flint_free(E);
        flint_free(pts);
    }

    flint_rand_clear(state);
    return 0;
}
