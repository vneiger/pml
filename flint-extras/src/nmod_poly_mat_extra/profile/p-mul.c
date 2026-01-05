/*
    Copyright (C) 2025 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>  // for atoi

#include <flint/ulong_extras.h>
#include <flint/profiler.h>
#include <flint/nmod.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_poly_mat.h>

#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_multiply.h"

#define MEASURE_SAMPLE 0

typedef struct
{
    slong rdim;  /* row outer dimension */
    slong idim;  /* inner dimension */
    slong cdim;  /* column outer dimension */
    slong deg;   /* degree */
    slong modn;  /* modulus */
}
time_args;

#define TIME_MUL(fun)                                   \
void time_##fun(time_args targs, flint_rand_t state)    \
{                                                       \
    const slong rdim = targs.rdim;                      \
    const slong idim = targs.idim;                      \
    const slong cdim = targs.cdim;                      \
    const slong deg = targs.deg;                        \
    const slong n = targs.modn;                         \
                                                        \
    nmod_t mod;                                         \
    nmod_init(&mod, n);                                 \
                                                        \
    nmod_poly_mat_t A;                                  \
    nmod_poly_mat_init(A, rdim, idim, n);               \
    nmod_poly_mat_rand(A, state, deg);                  \
    nmod_poly_mat_t B;                                  \
    nmod_poly_mat_init(B, idim, cdim, n);               \
    nmod_poly_mat_rand(B, state, deg);                  \
    nmod_poly_mat_t C;                                  \
    nmod_poly_mat_init(C, rdim, cdim, n);               \
                                                        \
    double FLINT_SET_BUT_UNUSED(tcpu), twall;           \
                                                        \
    TIMEIT_START;                                       \
    nmod_poly_mat_##fun(C, A, B);                       \
    TIMEIT_STOP_VALUES(tcpu, twall);                    \
                                                        \
    printf("%.2e", twall);                              \
                                                        \
    nmod_poly_mat_clear(A);                             \
    nmod_poly_mat_clear(B);                             \
    nmod_poly_mat_clear(C);                             \
}

TIME_MUL(mul)
TIME_MUL(mul_geometric)

/*-------------------------*/
/*  main                   */
/*-------------------------*/

int main(int argc, char ** argv)
{
    flint_rand_t state;
    flint_rand_init(state);
    flint_rand_set_seed(state, time(NULL), time(NULL)+129384125L);

    // modulus bitsize
    const slong nbits = 7;
    const ulong bits[] = {12, 24, 30, 40, 50, 60, 63};

    // matrix dimensions (all square for the moment)
    const slong ndims = 10;
    const ulong dims[] = {2, 4, 6, 8, 11, 15, 20, 30, 50, 100};

    // matrix degrees
    const slong ndegs = 12;
    const ulong degs[] = {5, 10, 20, 40, 80, 160, 320, 640, 1280, 2560, 5120, 10240};

    // bench functions
    const slong nfuns = 2;
    typedef void (*timefun) (time_args, flint_rand_t);
    const timefun funs[] = {
        time_mul,                      // 0
        time_mul_geometric,            // 1
    };

    // TODO
    //typedef void (*samplefun) (void*, ulong);
    //const samplefun sfuns[] = {
    //    sample_mul,                      // 0
    //    sample_mul_geometric,            // 1
    //};

    const char * description[] = {
        "#0  --> mul                          ",
        "#1  --> mul_geometric                ",
    };

    if (argc == 1)  // show usage
    {
        printf("Usage: `%s [nbits] [dim] [deg] [fun]`\n", argv[0]);
        printf("   Each argument is optional; no argument shows this help.\n");
        printf("   - nbits: number of bits (in (1..64]) for the modulus, chosen as nextprime(2**(nbits-1))\n");
        printf("        (nbits == -1 launches full suite)\n");
        printf("   - dim: matrices are dim x dim\n");
        printf("   - deg: matrices are random of degree < deg\n");
        printf("   - fun: id number of the timed function (see below),\n");
        printf("\nAvailable functions:\n");
        for (slong j = 0; j < nfuns; j++)
            printf("   %s\n", description[j]);

        return 0;
    }

    printf("#warmup...\n");
    for (slong i = 0; i < 3; i++)
    {
        time_args targs = {4, 4, 4, 1000, UWORD(1) << 20};
        time_mul(targs, state);
        printf(" ");
    }
    printf("\n\n");

    if (argc == 2 && atoi(argv[1]) == -1)  // launching full suite
    {
        printf("           dim");
        for (slong i = 0; i < ndims; i++)
            printf("%17ld", dims[i]);
        printf("\n");
        printf("bits fun deg\n");
        for (slong j = 0; j < nbits; j++)
        {
            const slong b = bits[j];
            ulong n;
            n = n_nextprime(UWORD(1) << (b-1), 0);
            for (slong ifun = 0; ifun < nfuns; ifun++)
            {
                for (slong d = 0; d < ndegs; d++)
                {
                    printf("%-5ld#%-3ld%-8ld", b, ifun, degs[d]);
                    for (slong i = 0; i < ndims; i++)
                    {
                        time_args targs = {dims[i], dims[i], dims[i], degs[d], n};

#if MEASURE_SAMPLE
                        const samplefun sfun = sfuns[ifun];
                        double min, max;
                        prof_repeat(&min, &max, sfun, (void*) &targs);
                        printf("%.2e", min/1000000);
#else
                        const timefun tfun = funs[ifun];
                        tfun(targs, state);
#endif
                        printf(" ");
                    }
                }
                printf("\n");
            }
        }
    }
    else if (argc == 2)  // nbits is given
    {
        printf("       dim");
        for (slong i = 0; i < ndims; i++)
            printf("%17ld", dims[i]);
        printf("\n");
        printf("bits fun deg\n");
        const slong b = atoi(argv[1]);
        ulong n;
        n = n_nextprime(UWORD(1) << (b-1), 0);
        for (slong ifun = 0; ifun < nfuns; ifun++)
        {
            const timefun tfun = funs[ifun];
            for (slong d = 0; d < ndegs; d++)
            {
                printf("%-5ld#%-3ld%-8ld", b, ifun, degs[d]);
                for (slong i = 0; i < ndims; i++)
                {
                    time_args targs = {dims[i], dims[i], dims[i], degs[d], n};
                    tfun(targs, state);
                    printf(" ");
                }
            }
            printf("\n");
        }
    }
    else if (argc == 3)  // nbits + dim given
    {
        const slong dim = atoi(argv[2]);
        printf("       dim");
        printf("%17ld", dim);
        printf("\n");
        printf("bits fun deg\n");
        const slong b = atoi(argv[1]);
        ulong n;
        n = n_nextprime(UWORD(1) << (b-1), 0);
        for (slong ifun = 0; ifun < nfuns; ifun++)
        {
            const timefun tfun = funs[ifun];
            for (slong d = 0; d < ndegs; d++)
            {
                printf("%-5ld#%-3ld%-8ld", b, ifun, degs[d]);
                time_args targs = {dim, dim, dim, degs[d], n};
                tfun(targs, state);
                printf(" ");
                printf("\n");
            }
        }
    }
    else if (argc == 4)  // nbits + dim + deg given
    {
        const slong dim = atoi(argv[2]);
        const slong deg = atoi(argv[3]);
        printf("       dim");
        printf("%17ld", dim);
        printf("\n");
        printf("bits fun deg\n");
        const slong b = atoi(argv[1]);
        ulong n;
        n = n_nextprime(UWORD(1) << (b-1), 0);
        for (slong ifun = 0; ifun < nfuns; ifun++)
        {
            const timefun tfun = funs[ifun];
            printf("%-5ld#%-3ld%-8ld", b, ifun, deg);
            time_args targs = {dim, dim, dim, deg, n};
            tfun(targs, state);
            printf(" ");
            printf("\n");
        }
    }
    else if (argc == 5)  // nbits + dim + deg + fun given
    {
        const slong b = atoi(argv[1]);
        const slong dim = atoi(argv[2]);
        const slong deg = atoi(argv[3]);
        const slong ifun = atoi(argv[4]);
        const timefun tfun = funs[ifun];
        printf("       dim");
        printf("%17ld", dim);
        printf("\n");
        printf("bits fun deg\n");
        ulong n;
        n = n_nextprime(UWORD(1) << (b-1), 0);

        printf("%-5ld#%-3ld%-8ld", b, ifun, deg);
        time_args targs = {dim, dim, dim, deg, n};
        tfun(targs, state);

        printf(" ");
        printf("\n");
    }

    flint_rand_clear(state);
    return 0;
}
