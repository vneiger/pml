/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/flint.h>
#include <flint/nmod_poly.h>
#include <stdlib.h>  // for atoi

#include <flint/ulong_extras.h>
#include <flint/profiler.h>
#include <flint/nmod.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_poly_mat.h>

#include "nmod_poly_mat_forms.h"
#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_kernel.h"

#define MEASURE_SAMPLE 0

/* Shift types are:
 * -1 -> [TODO not implemented yet] negative shift that represents known degree bound on a kernel basis
 * 0 -> input degrees: shift == row degrees of input matrix
 * 1 -> uniform: shift == (0,...,0)
 * 100 -> [TODO not implemented yet] shift that yields the kernel basis Hermite form
 * */

typedef struct
{
    slong rdim;  /* row dimension */
    slong cdim;  /* column dimension */
    slong deg;   /* degree */
    slong rank;  /* rank */
    slong stype; /* shift type */
    slong modn;  /* modulus */
}
time_args;

#define TIME_KER(fun)                                   \
void time_##fun(time_args targs, flint_rand_t state)    \
{                                                       \
    const slong rdim = targs.rdim;                      \
    const slong cdim = targs.cdim;                      \
    const slong deg = targs.deg;                        \
    /* const slong rank = targs.rank; */ /* TODO */     \
    const slong stype = targs.stype;                    \
    const slong n = targs.modn;                         \
                                                        \
    nmod_t mod;                                         \
    nmod_init(&mod, n);                                 \
                                                        \
    nmod_poly_mat_t pmat;                               \
    nmod_poly_mat_init(pmat, rdim, cdim, n);            \
    nmod_poly_mat_rand(pmat, state, deg);               \
                                                        \
    slong * pivind = FLINT_ARRAY_ALLOC(rdim, slong);    \
    slong * shift = FLINT_ARRAY_ALLOC(rdim, slong);     \
                                                        \
    if (stype != 0 && stype != 1)                       \
        flint_printf("~warning~ unsupported stype\n");  \
                                                        \
    nmod_poly_mat_t ker;                                \
    nmod_poly_mat_init(ker, rdim, rdim, n);             \
                                                        \
    double FLINT_SET_BUT_UNUSED(tcpu), twall;           \
                                                        \
    TIMEIT_START;                                       \
    nmod_poly_mat_t copy_pmat;                          \
    nmod_poly_mat_init_set(copy_pmat, pmat);            \
    if (stype == 1)                                     \
        for (slong i = 0; i < rdim; i++)                \
            shift[i] = 0;                               \
    else if (stype == 0)                                \
        nmod_poly_mat_row_degree_zero(shift,            \
                                      pmat, NULL);      \
                                                        \
    nmod_poly_mat_##fun(ker, pivind, shift, copy_pmat); \
    TIMEIT_STOP_VALUES(tcpu, twall);                    \
                                                        \
    flint_printf("%.2e", twall);                        \
                                                        \
    flint_free(pivind);                                 \
    flint_free(shift);                                  \
    nmod_poly_mat_clear(pmat);                          \
    nmod_poly_mat_clear(ker);                           \
}

/* similar function, but creating matrix with unbalanced row degrees */
#define TIME_KER_UNBALANCED(fun)                                    \
void time_##fun##_unbalanced(time_args targs, flint_rand_t state)     \
{                                                                   \
    const slong rdim = targs.rdim;                                  \
    const slong cdim = targs.cdim;                                  \
    const slong deg = targs.deg;                                    \
    /* const slong rank = targs.rank; */ /* TODO */                 \
    const slong stype = targs.stype;                                \
    const slong n = targs.modn;                                     \
                                                                    \
    nmod_t mod;                                                     \
    nmod_init(&mod, n);                                             \
                                                                    \
    nmod_poly_mat_t pmat;                                           \
    nmod_poly_mat_init(pmat, rdim, cdim, n);                        \
    for (slong i = 0; i < rdim; i++)                                \
    {                                                               \
        slong d = n_randint(state, deg);                            \
        for (slong j = 0; j < cdim; j++)                            \
            nmod_poly_randtest(nmod_poly_mat_entry(pmat, i, j),     \
                               state, d);                           \
    }                                                               \
                                                                    \
    slong * pivind = FLINT_ARRAY_ALLOC(rdim, slong);                \
    slong * shift = FLINT_ARRAY_ALLOC(rdim, slong);                 \
                                                                    \
    if (stype != 0 && stype != 1)                                   \
        flint_printf("~warning~ unsupported stype\n");              \
                                                                    \
    nmod_poly_mat_t ker;                                            \
    nmod_poly_mat_init(ker, rdim, rdim, n);                         \
                                                                    \
    double FLINT_SET_BUT_UNUSED(tcpu), twall;                       \
                                                                    \
    TIMEIT_START;                                                   \
    nmod_poly_mat_t copy_pmat;                                      \
    nmod_poly_mat_init_set(copy_pmat, pmat);                        \
    if (stype == 1)                                                 \
        for (slong i = 0; i < rdim; i++)                            \
            shift[i] = 0;                                           \
    else if (stype == 0)                                            \
        nmod_poly_mat_row_degree_zero(shift, pmat, NULL);           \
    nmod_poly_mat_##fun(ker, pivind, shift, copy_pmat);             \
    TIMEIT_STOP_VALUES(tcpu, twall);                                \
                                                                    \
    flint_printf("%.2e", twall);                                    \
                                                                    \
    flint_free(pivind);                                             \
    flint_free(shift);                                              \
    nmod_poly_mat_clear(pmat);                                      \
    nmod_poly_mat_clear(ker);                                       \
}

void time_nullspace(time_args targs, flint_rand_t state)
{
    const slong rdim = targs.rdim;
    const slong cdim = targs.cdim;
    const slong deg = targs.deg;
    /* const slong rank = targs.rank; */ /* TODO */
    const slong n = targs.modn;

    nmod_t mod;
    nmod_init(&mod, n);

    nmod_poly_mat_t pmat;
    nmod_poly_mat_init(pmat, cdim, rdim, n);
    nmod_poly_mat_rand(pmat, state, deg);

    nmod_poly_mat_t ker;
    nmod_poly_mat_init(ker, rdim, rdim, n);

    double FLINT_SET_BUT_UNUSED(tcpu), twall;

    TIMEIT_START;
    nmod_poly_mat_nullspace(ker, pmat);
    TIMEIT_STOP_VALUES(tcpu, twall);

    flint_printf("%.2e", twall);

    nmod_poly_mat_clear(pmat);
    nmod_poly_mat_clear(ker);
}

void time_nullspace_unbalanced(time_args targs, flint_rand_t state)
{
    const slong rdim = targs.rdim;
    const slong cdim = targs.cdim;
    const slong deg = targs.deg;
    /* const slong rank = targs.rank; */ /* TODO */
    const slong n = targs.modn;

    nmod_t mod;
    nmod_init(&mod, n);

    nmod_poly_mat_t pmat;
    nmod_poly_mat_init(pmat, cdim, rdim, n);
    for (slong j = 0; j < rdim; j++)
    {
        slong d = n_randint(state, deg);
        for (slong i = 0; i < cdim; i++)
            nmod_poly_randtest(nmod_poly_mat_entry(pmat, i, j),
                               state, d);
    }

    nmod_poly_mat_t ker;
    nmod_poly_mat_init(ker, rdim, rdim, n);

    double FLINT_SET_BUT_UNUSED(tcpu), twall;

    TIMEIT_START;
    nmod_poly_mat_nullspace(ker, pmat);
    TIMEIT_STOP_VALUES(tcpu, twall);

    flint_printf("%.2e", twall);

    nmod_poly_mat_clear(pmat);
    nmod_poly_mat_clear(ker);
}


TIME_KER(kernel_via_approx)
TIME_KER(kernel_zls_approx)
TIME_KER_UNBALANCED(kernel_via_approx)
TIME_KER_UNBALANCED(kernel_zls_approx)

/*-------------------------*/
/*  main                   */
/*-------------------------*/

int main(int argc, char ** argv)
{
    flint_rand_t state;
    flint_rand_init(state);
    flint_rand_set_seed(state, time(NULL), time(NULL)+129384125L);

    /* TODO rank */

    // bench functions
    const slong nfuns = 6;
    typedef void (*timefun) (time_args, flint_rand_t);
    const timefun funs[] = {
        time_nullspace,                       // 0
        time_kernel_via_approx,               // 1
        time_kernel_zls_approx,               // 2
        time_nullspace_unbalanced,            // 3
        time_kernel_via_approx_unbalanced,    // 4
        time_kernel_zls_approx_unbalanced,    // 5
    };

    // TODO
    //typedef void (*samplefun) (void*, ulong);
    //const samplefun sfuns[] = {
    //    sample_mul,                      // 0
    //    sample_mul_geometric,            // 1
    //};
    //#if MEASURE_SAMPLE
    //                        const samplefun sfun = sfuns[ifun];
    //                        double min, max;
    //                        prof_repeat(&min, &max, sfun, (void*) &targs);
    //                        flint_printf("%.2e", min/1000000);
    //#else
    //                        const timefun tfun = funs[ifun];
    //                        tfun(targs, state);
    //#endif

    const char * description[] = {
        "#0  --> FLINT's nullspace                            ",
        "#1  --> via approx                                   ",
        "#2  --> ZLS via approx                               ",
        "#3  --> FLINT's nullspace, unbalanced row degrees    ",
        "#4  --> via approx, unbalanced row degrees           ",
        "#5  --> ZLS via approx, unbalanced row degrees       ",
    };

    if (argc != 7)  // show usage
    {
        flint_printf("Usage: `%s [fun] [nbits] [rdim] [cdim] [deg] [shift] [rank]`\n", argv[0]);
        flint_printf("   No argument shows this help.\n");
        flint_printf("   All arguments are mandatory (except the unsupported rank).\n");
        flint_printf("   [rank] not supported yet.\n");
        flint_printf("   - fun: id number of the timed function (see list below, -1 times all)\n");
        flint_printf("   - nbits: number of bits in [2..64] for the modulus, chosen as nextprime(2**(nbits-1))\n");
        flint_printf("   - rdim, cdim: input matrix is rdim x cdim\n");
        flint_printf("   - deg: matrix is random of degree < deg\n");
        flint_printf("   - shift: type of shift (0 : row degree | 1 : (0,..,0))\n");
        flint_printf("   - rank: [unsupported] matrix is random of this rank\n");
        flint_printf("\nAvailable functions:\n");
        for (slong j = 0; j < nfuns; j++)
            flint_printf("   %s\n", description[j]);

        return 0;
    }

    flint_printf("#warmup...\n");
    for (slong i = 0; i < 3; i++)
    {
        /* rdim; cdim; deg; rank; stype; modn; */
        time_args targs = {8, 4, 1000, 4, 0, n_nextprime(UWORD(1) << 20, 0)};
        time_kernel_via_approx(targs, state);
        flint_printf(" ");
    }
    flint_printf("\n\n");

    if (argc == 7)  // fun + nbits + rdim + cdim + deg + shift given
    {
        const slong ifun  = atoi(argv[1]);
        const slong bits  = atoi(argv[2]);
        const slong rdim  = atoi(argv[3]);
        const slong cdim  = atoi(argv[4]);
        const slong deg   = atoi(argv[5]);
        const slong stype = atoi(argv[6]);
        const ulong n = n_nextprime(UWORD(1) << (bits-1), 0);
        if (ifun >= 0)
        {
            flint_printf("bits fun rdim cdim deg stype\n");
            flint_printf("%-5ld#%-3ld%-5ld%-5ld%-8ld%-2ld", bits, ifun, rdim, cdim, deg, stype);
            const timefun tfun = funs[ifun];
            /* rdim; cdim; deg; rank; stype; modn; */
            time_args targs = {rdim, cdim, deg, FLINT_MIN(rdim, cdim), stype, n};
            tfun(targs, state);
            flint_printf("\n");
        }
        else if (ifun == -1)
        {
            flint_printf("bits fun rdim cdim deg stype    flint   approx  zls-app  u_flint    u_app  u_zls-a   \n");
            flint_printf("%-5ld#%-3ld%-5ld%-5ld%-8ld%-2ld", bits, ifun, rdim, cdim, deg, stype);
            for (slong i = 0; i < nfuns; i++)
            {
                const timefun tfun = funs[i];
                time_args targs = {rdim, cdim, deg, FLINT_MIN(rdim, cdim), stype, n};
                tfun(targs, state);
                flint_printf(" ");
            }
            flint_printf("\n");
        }
    }

    flint_rand_clear(state);
    return 0;
}
