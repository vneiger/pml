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
    slong len1;   /* length of left operand */
    slong len2;   /* length of right operand */
    slong modn;  /* modulus */
}
time_args;

#define TIME_MUL(fun)                                   \
void time_##fun(time_args targs, flint_rand_t state)    \
{                                                       \
    const slong rdim = targs.rdim;                      \
    const slong idim = targs.idim;                      \
    const slong cdim = targs.cdim;                      \
    const slong len1 = targs.len1;                      \
    const slong len2 = targs.len2;                      \
    const slong modn = targs.modn;                      \
                                                        \
    nmod_t mod;                                         \
    nmod_init(&mod, modn);                              \
                                                        \
    nmod_poly_mat_t A;                                  \
    nmod_poly_mat_init(A, rdim, idim, modn);            \
    nmod_poly_mat_rand(A, state, len1);                 \
    nmod_poly_mat_t B;                                  \
    nmod_poly_mat_init(B, idim, cdim, modn);            \
    nmod_poly_mat_rand(B, state, len2);                 \
    nmod_poly_mat_t C;                                  \
    nmod_poly_mat_init(C, rdim, cdim, modn);            \
                                                        \
    double FLINT_SET_BUT_UNUSED(tcpu), twall;           \
                                                        \
    TIMEIT_START;                                       \
    nmod_poly_mat_##fun(C, A, B);                       \
    TIMEIT_STOP_VALUES(tcpu, twall);                    \
                                                        \
    flint_printf("%.2e", twall);                        \
                                                        \
    nmod_poly_mat_clear(A);                             \
    nmod_poly_mat_clear(B);                             \
    nmod_poly_mat_clear(C);                             \
}

TIME_MUL(multiply)
TIME_MUL(mul)
TIME_MUL(mul_geometric)
TIME_MUL(mul_waksman)
TIME_MUL(mul_vandermonde1)
TIME_MUL(mul_vandermonde2)

/*-------------------------*/
/*  main                   */
/*-------------------------*/

int main(int argc, char ** argv)
{
    flint_rand_t state;
    flint_rand_init(state);
    flint_rand_set_seed(state, time(NULL), time(NULL)+129384125L);

    // bench functions
    const slong nfuns = 3;
    typedef void (*timefun) (time_args, flint_rand_t);
    const timefun funs[] = {
        time_multiply,                      // 0
        time_mul,                           // 1
        time_mul_waksman,                   // 2
        time_mul_geometric,                 // 3
        time_mul_vandermonde1,              // 4
        time_mul_vandermonde2,              // 5
    };

    const char * description[] = {
        "#0  --> multiply                     ",
        "#1  --> mul                          ",
        "#2  --> mul_waksman                  ",
        "#3  --> mul_geometric                ",
        "#4  --> mul_vandermonde1             ",
        "#5  --> mul_vandermonde2             ",
    };

    if (argc == 1)  // show usage
    {
        flint_printf("Usage: `%s [nbits] [fun] [dim1] [dim2] [dim3] [len1] [len2] [opt:hide_desc]`\n", argv[0]);
        flint_printf("   No argument shows this help.\n");
        flint_printf("   - nbits: number of bits (in (1..64]) for the modulus, chosen as nextprime(2**(nbits-1))\n");
        flint_printf("   - fun: id number of the timed function (see below),\n");
        flint_printf("      (fun == -1 launches all available functions)\n");
        flint_printf("   - dim1, dim2, dim3: matrices are dim1 x dim2 and dim2 x dim3\n");
        flint_printf("   - len1, len2: matrices are random of length <= len1 and <= len2\n");
        flint_printf("   - hide_desc: optional, if nonzero, will not show the first line that describes parameters name\n");
        flint_printf("\nAvailable functions:\n");
        for (slong j = 0; j < nfuns; j++)
            flint_printf("   %s\n", description[j]);

        return 0;
    }

    if (argc >= 8)  // nbits + fun + dim1 + dim2 + dim3 + len1 + len2
    {
        const slong b = atol(argv[1]);
        const slong ifun = atol(argv[2]);
        const slong dim1 = atol(argv[3]);
        const slong dim2 = atol(argv[4]);
        const slong dim3 = atol(argv[5]);
        const slong len1 = atol(argv[6]);
        const slong len2 = atol(argv[7]);

        const int hide_desc = (argc == 9) ? atoi(argv[8]) : 0;
        if (!hide_desc)
        {
            flint_printf("Available functions:\n");
            for (slong j = 0; j < nfuns; j++)
                flint_printf("   %s\n", description[j]);
            if (ifun >= 0)
                flint_printf("bits dim1 dim2 dim3 len1    len2    fun#%-5lu\n", ifun);
            else
            {
                flint_printf("bits dim1 dim2 dim3 len1    len2    ");
                for (slong fun_nb = 0; fun_nb < nfuns; fun_nb++)
                    flint_printf("fun#%-5lu", fun_nb);
                flint_printf("\n");
            }
        }

        ulong modn = n_nextprime(UWORD(1) << (b-1), 0);

        if (ifun == -1)
        {
            flint_printf("%-5ld%-5ld%-5ld%-5ld%-8ld%-8ld", b, dim1, dim2, dim3, len1, len2);
            for (slong fun_nb = 0; fun_nb < nfuns; fun_nb++)
            {
                const timefun tfun = funs[fun_nb];
                time_args targs = {dim1, dim2, dim3, len1, len2, modn};
                tfun(targs, state);
                flint_printf(" ");
            }
        }
        else
        {
            const timefun tfun = funs[ifun];
            flint_printf("%-5ld%-5ld%-5ld%-5ld%-8ld%-8ld", b, ifun, dim1, dim2, dim3, len1, len2);
            time_args targs = {dim1, dim2, dim3, len1, len2, modn};
            tfun(targs, state);
        }
        flint_printf("\n");
    }

    flint_rand_clear(state);
    return 0;
}
