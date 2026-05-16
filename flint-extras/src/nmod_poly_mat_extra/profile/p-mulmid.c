/*
    Copyright (C) 2025 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/nmod_poly.h>
#include <stdlib.h>  // for atoi

#include <flint/ulong_extras.h>
#include <flint/profiler.h>
#include <flint/nmod.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_poly_mat.h>

#include "nmod_extra.h"
#include "nmod_poly_extra.h"
#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_multiply.h"

typedef struct
{
    slong rdim;  /* row outer dimension */
    slong idim;  /* inner dimension */
    slong cdim;  /* column outer dimension */
    slong len1;  /* length of left operand */
    slong len2;  /* length of right operand */
    slong nlo;   /* nlo of mulmid */
    slong nhi;   /* nhi of mulmid */
    slong modn;  /* modulus */
}
time_args;

void time_mulmid(time_args targs, flint_rand_t state)
{
    const slong rdim = targs.rdim;
    const slong idim = targs.idim;
    const slong cdim = targs.cdim;
    const slong len1 = targs.len1;
    const slong len2 = targs.len2;
    const slong nlo  = targs.nlo;
    const slong nhi  = targs.nhi;
    const slong modn = targs.modn;

    nmod_t mod;
    nmod_init(&mod, modn);

    nmod_poly_mat_t pmat1;
    nmod_poly_mat_init(pmat1, rdim, idim, modn);
    nmod_poly_mat_rand(pmat1, state, len1);
    nmod_poly_mat_t pmat2;
    nmod_poly_mat_init(pmat2, idim, cdim, modn);
    nmod_poly_mat_rand(pmat2, state, len2);
    nmod_poly_mat_t res;
    nmod_poly_mat_init(res, rdim, cdim, modn);

    double FLINT_SET_BUT_UNUSED(tcpu), twall;

    TIMEIT_START;
    nmod_poly_mat_mulmid(res, pmat1, pmat2, nlo, nhi);
    TIMEIT_STOP_VALUES(tcpu, twall);

    flint_printf("%.2e", twall);

    nmod_poly_mat_clear(pmat1);
    nmod_poly_mat_clear(pmat2);
    nmod_poly_mat_clear(res);
}


#define TIME_MULMID_METHOD(fun)                          \
void time_##fun(time_args targs, flint_rand_t state)     \
{                                                        \
    const slong rdim = targs.rdim;                       \
    const slong idim = targs.idim;                       \
    const slong cdim = targs.cdim;                       \
    const slong len1 = targs.len1;                       \
    const slong len2 = targs.len2;                       \
    const slong nlo  = targs.nlo;                        \
    const slong nhi  = targs.nhi;                        \
    const slong modn = targs.modn;                       \
                                                         \
    nmod_t mod;                                          \
    nmod_init(&mod, modn);                               \
                                                         \
    nmod_poly_mat_t pmat1;                               \
    nmod_poly_mat_init(pmat1, rdim, idim, modn);         \
    nmod_poly_mat_rand(pmat1, state, len1);              \
    nmod_poly_mat_t pmat2;                               \
    nmod_poly_mat_init(pmat2, idim, cdim, modn);         \
    nmod_poly_mat_rand(pmat2, state, len2);              \
    nmod_poly_mat_t res;                                 \
    nmod_poly_mat_init(res, rdim, cdim, modn);           \
                                                         \
    double FLINT_SET_BUT_UNUSED(tcpu), twall;            \
                                                         \
    TIMEIT_START;                                        \
    _nmod_poly_mat_##fun(res, pmat1, len1, pmat2, len2, nlo, nhi); \
    TIMEIT_STOP_VALUES(tcpu, twall);                     \
                                                         \
    flint_printf("%.2e", twall);                         \
                                                         \
    nmod_poly_mat_clear(pmat1);                          \
    nmod_poly_mat_clear(pmat2);                          \
    nmod_poly_mat_clear(res);                            \
}

TIME_MULMID_METHOD(mulmid_naive)
#if (__FLINT_VERSION == 3 && __FLINT_VERSION_MINOR >= 6)
TIME_MULMID_METHOD(mulmid_geometric)
#endif

#if (__FLINT_VERSION == 3 && __FLINT_VERSION_MINOR >= 6)
void time_mulmid_geometric_precomp(time_args targs, flint_rand_t state)
{
    const slong rdim = targs.rdim;
    const slong idim = targs.idim;
    const slong cdim = targs.cdim;
    const slong len1 = targs.len1;
    const slong len2 = targs.len2;
    const slong nlo  = targs.nlo;
    const slong nhi  = targs.nhi;
    const slong modn = targs.modn;

    nmod_t mod;
    nmod_init(&mod, modn);

    nmod_poly_mat_t pmat1;
    nmod_poly_mat_init(pmat1, rdim, idim, modn);
    nmod_poly_mat_rand(pmat1, state, len1);
    nmod_poly_mat_t pmat2;
    nmod_poly_mat_init(pmat2, idim, cdim, modn);
    nmod_poly_mat_rand(pmat2, state, len2);
    nmod_poly_mat_t res;
    nmod_poly_mat_init(res, rdim, cdim, modn);

    double FLINT_SET_BUT_UNUSED(tcpu), twall;

    nmod_geometric_progression_t G;
    ulong w = nmod_find_root(2*nhi, mod);
    _nmod_geometric_progression_init_function(G, w, nhi, mod, UWORD(3));

    TIMEIT_START;
    _nmod_poly_mat_mulmid_geometric_precomp(res, pmat1, len1, pmat2, len2, nlo, nhi, G);
    TIMEIT_STOP_VALUES(tcpu, twall);

    nmod_geometric_progression_clear(G);

    flint_printf("%.2e", twall);

    nmod_poly_mat_clear(pmat1);
    nmod_poly_mat_clear(pmat2);
    nmod_poly_mat_clear(res);
}
#endif

/*-------------------------*/
/*  main                   */
/*-------------------------*/

int main(int argc, char ** argv)
{
    flint_rand_t state;
    flint_rand_init(state);
    flint_rand_set_seed(state, time(NULL), time(NULL)+129384125L);

    // bench functions
#if (__FLINT_VERSION == 3 && __FLINT_VERSION_MINOR >= 6)
    const slong nfuns = 4;
    typedef void (*timefun) (time_args, flint_rand_t);
    const timefun funs[] = {
        time_mulmid,                              // 0
        time_mulmid_naive,                        // 1
        time_mulmid_geometric,                    // 2
        time_mulmid_geometric_precomp,            // 3
#else
    const slong nfuns = 2;
    typedef void (*timefun) (time_args, flint_rand_t);
    const timefun funs[] = {
        time_mulmid,                              // 0
        time_mulmid_naive,                        // 1
#endif
    };

    const char * description[] = {
        "#0  --> general interface                       ",
        "#1  --> mulmid_naive                            ",
        "#2  --> mulmid_geometric                        ",
        "#3  --> mulmid_geometric_precomp                ",
    };

    if (argc == 1)  // show usage
    {
        flint_printf("Usage: `%s [nbits] [fun] [dim1] [dim2] [dim3] [len1] [len2] [nlo] [nhi] [opt:hide_desc]`\n", argv[0]);
        flint_printf("   No argument shows this help.\n");
        flint_printf("   - nbits: number of bits (in (1..64]) for the modulus, chosen as nextprime(2**(nbits-1))\n");
        flint_printf("   - fun: id number of the timed function (see below),\n");
        flint_printf("      (fun == -1 launches all available functions)\n");
        flint_printf("   - dim1, dim2, dim3: matrices are dim1 x dim2 and dim2 x dim3\n");
        flint_printf("   - len1, len2: matrices are random of length <= len1 and <= len2\n");
        flint_printf("   - nlo, nhi: compute middle coefficients nlo...nhi-1 of product\n");
        flint_printf("   - hide_desc: optional, if nonzero, will not show the first line that describes parameters name\n");
        flint_printf("\nAvailable functions:\n");
        for (slong j = 0; j < nfuns; j++)
            flint_printf("   %s\n", description[j]);

        return 0;
    }

    if (argc >= 8)  // nbits + fun + dim1 + dim2 + dim3 + len1 + len2 + nlo + nhi
    {
        const slong b = atol(argv[1]);
        const slong ifun = atol(argv[2]);
        const slong dim1 = atol(argv[3]);
        const slong dim2 = atol(argv[4]);
        const slong dim3 = atol(argv[5]);
        const slong len1 = atol(argv[6]);
        const slong len2 = atol(argv[7]);
        const slong nlo = atol(argv[8]);
        const slong nhi = atol(argv[9]);

        const int hide_desc = (argc == 11) ? atoi(argv[10]) : 0;
        if (!hide_desc)
        {
            flint_printf("Available functions:\n");
            for (slong j = 0; j < nfuns; j++)
                flint_printf("   %s\n", description[j]);
            if (ifun >= 0)
                flint_printf("bits dim1 dim2 dim3 len1    len2    nlo     nhi     fun#%-5lu\n", ifun);
            else
            {
                flint_printf("bits dim1 dim2 dim3 len1    len2    nlo     nhi     ");
                for (slong fun_nb = 0; fun_nb < nfuns; fun_nb++)
                    flint_printf("fun#%-5lu", fun_nb);
                flint_printf("\n");
            }
        }

        ulong modn = n_nextprime(UWORD(1) << (b-1), 0);
        const int can_use_geom = NMOD_CAN_USE_GEOMETRIC(modn, FLINT_MAX(len1, len2));

        flint_printf("%-5ld%-5ld%-5ld%-5ld%-8ld%-8ld%-8ld%-8ld", b, dim1, dim2, dim3, len1, len2, nlo, nhi);
        if (ifun == -1)
        {
            for (slong fun_nb = 0; fun_nb < nfuns; fun_nb++)
            {
                if (fun_nb >= 2 && fun_nb <= 3 && !can_use_geom)  /* geometric not feasible */
                    flint_printf("n/a      ");
                else
                {
                    const timefun tfun = funs[fun_nb];
                    time_args targs = {dim1, dim2, dim3, len1, len2, nlo, nhi, modn};
                    tfun(targs, state);
                    flint_printf(" ");
                }
            }
        }
        else
        {
            if (ifun >= 2 && ifun <= 3 && !can_use_geom)  /* geometric not feasible */
                flint_printf("n/a      ");
            else
            {
                const timefun tfun = funs[ifun];
                time_args targs = {dim1, dim2, dim3, len1, len2, nlo, nhi, modn};
                tfun(targs, state);
            }
        }
        flint_printf("\n");
    }

    flint_rand_clear(state);
    return 0;
}
