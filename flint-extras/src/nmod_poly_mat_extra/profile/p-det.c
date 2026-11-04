/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>  // for atoi

#include <flint/flint.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>
#include <flint/profiler.h>
#include <flint/ulong_extras.h>

#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_determinant.h"

/* Input types are:
 * 0 -> random n x n matrix of degree < deg
 * 1 -> U * D * V with U, V random of degree < deg/2 and D diagonal with
 *      D[i,i] = 1 for i < n/2 and D[i,i] = x^2+x otherwise
 * 2 -> x*I - M with M random constant (characteristic matrix; deg ignored)
 * 3 -> x*I - M with M similar to a block diagonal matrix diag(B, B) with B
 *      random (derogatory: two nontrivial invariant factors; deg ignored)
 * 4 -> unbalanced row degrees: row 0 has degree < deg, the other rows have
 *      degree < 2
 * 5 -> unbalanced column degrees: the last column has degree < deg, the other
 *      columns have degree < 2
 */

typedef struct
{
    slong dim;   /* dimension */
    slong deg;   /* degree */
    slong itype; /* input type */
    slong modn;  /* modulus */
}
time_args;

static void _build_input(nmod_poly_mat_t F, time_args targs, flint_rand_t state)
{
    const slong n = targs.dim;
    const slong deg = targs.deg;
    if (targs.itype == 0)
        nmod_poly_mat_rand(F, state, deg);
    else if (targs.itype == 4 || targs.itype == 5)
    {
        slong * len = flint_malloc(n * sizeof(slong));
        for (slong i = 0; i < n; i++)
            len[i] = FLINT_MIN(2, deg);
        if (targs.itype == 4)
        {
            len[0] = deg;
            nmod_poly_mat_rand_row_degree(F, state, len);
        }
        else
        {
            len[n-1] = deg;
            nmod_poly_mat_rand_column_degree(F, state, len);
        }
        flint_free(len);
    }
    else if (targs.itype == 1)
    {
        nmod_poly_mat_t U, V;
        nmod_poly_mat_init(U, n, n, targs.modn);
        nmod_poly_mat_init(V, n, n, targs.modn);
        nmod_poly_mat_rand(U, state, FLINT_MAX(1, deg/2));
        nmod_poly_mat_rand(V, state, FLINT_MAX(1, deg/2));
        nmod_poly_t f;
        nmod_poly_init(f, targs.modn);
        nmod_poly_set_coeff_ui(f, 2, 1);
        nmod_poly_set_coeff_ui(f, 1, 1);
        for (slong i = 0; i < n; i++)
            for (slong j = n/2; j < n; j++)
                nmod_poly_mul(nmod_poly_mat_entry(U, i, j), nmod_poly_mat_entry(U, i, j), f);
        nmod_poly_clear(f);
        nmod_poly_mat_mul(F, U, V);
        nmod_poly_mat_clear(U);
        nmod_poly_mat_clear(V);
    }
    else
    {
        nmod_mat_t M;
        nmod_mat_init(M, n, n, targs.modn);
        if (targs.itype == 2)
            nmod_mat_rand(M, state);
        else
        {
            /* random similarity of diag(B, B), as in p-hermite_diagonal.c */
            const slong b = n / 2;
            nmod_mat_t B, P, Pinv, T;
            nmod_mat_init(B, b, b, targs.modn);
            nmod_mat_rand(B, state);
            nmod_mat_zero(M);
            for (slong i = 0; i < b; i++)
                for (slong j = 0; j < b; j++)
                {
                    nmod_mat_entry(M, i, j) = nmod_mat_entry(B, i, j);
                    nmod_mat_entry(M, b + i, b + j) = nmod_mat_entry(B, i, j);
                }
            if (n % 2)
                nmod_mat_entry(M, n-1, n-1) = n_randint(state, targs.modn);
            nmod_mat_clear(B);
            nmod_mat_init(P, n, n, targs.modn);
            nmod_mat_init(Pinv, n, n, targs.modn);
            do
                nmod_mat_rand(P, state);
            while (!nmod_mat_inv(Pinv, P));
            nmod_mat_init(T, n, n, targs.modn);
            nmod_mat_mul(T, P, M);
            nmod_mat_mul(M, T, Pinv);
            nmod_mat_clear(T);
            nmod_mat_clear(P);
            nmod_mat_clear(Pinv);
        }
        nmod_poly_mat_zero(F);
        for (slong i = 0; i < n; i++)
        {
            for (slong j = 0; j < n; j++)
                nmod_poly_set_coeff_ui(nmod_poly_mat_entry(F, i, j), 0,
                                       nmod_neg(nmod_mat_entry(M, i, j), M->mod));
            nmod_poly_set_coeff_ui(nmod_poly_mat_entry(F, i, i), 1, 1);
        }
        nmod_mat_clear(M);
    }
}

/* applicable: whether the method applies (evaluated once, before timing: an
 * empty timed loop may be optimized away, and the timing loop of
 * TIMEIT_START, which multiplies the number of repetitions by 10 until 0.1 s
 * of cpu time is reached, would then never stop) */
#define TIME_DET(name, applicable, call)                         \
static void time_##name(time_args targs, flint_rand_t state)     \
{                                                                \
    const slong n = targs.dim;                                   \
    nmod_poly_mat_t F, copy;                                     \
    nmod_poly_mat_init(F, n, n, targs.modn);                     \
    nmod_poly_mat_init(copy, n, n, targs.modn);                  \
    _build_input(F, targs, state);                               \
    nmod_poly_t det;                                             \
    nmod_poly_init(det, targs.modn);                             \
                                                                 \
    if (applicable)                                              \
    {                                                            \
        double FLINT_SET_BUT_UNUSED(tcpu), twall;                \
        TIMEIT_START;                                            \
        call;                                                    \
        TIMEIT_STOP_VALUES(tcpu, twall);                         \
        flint_printf("%.2e", twall);                             \
    }                                                            \
    else                                                         \
        flint_printf("    n/a ");                                \
                                                                 \
    nmod_poly_clear(det);                                        \
    nmod_poly_mat_clear(F);                                      \
    nmod_poly_mat_clear(copy);                                   \
}

TIME_DET(det_iter, 1, nmod_poly_mat_set(copy, F); nmod_poly_mat_det_iter(det, copy))
TIME_DET(det_minors, n <= 12, nmod_poly_mat_det_minors(det, F))
TIME_DET(det_geometric, nmod_poly_mat_det_geometric(det, F), nmod_poly_mat_det_geometric(det, F))
TIME_DET(det_flint, 1, nmod_poly_mat_det(det, F))
TIME_DET(det_lnz, 1, nmod_poly_mat_det_lnz(det, F))
TIME_DET(determinant, 1, nmod_poly_mat_determinant(det, F))

/* FLINT's charpoly of the constant matrix M, for the characteristic matrix
 * F = x I - M (types 2 and 3), since det(F) = charpoly(M); n/a for other
 * types */
static void time_charpoly(time_args targs, flint_rand_t state)
{
    if (targs.itype != 2 && targs.itype != 3)
    {
        flint_printf("    n/a ");
        return;
    }
    const slong n = targs.dim;
    nmod_poly_mat_t F;
    nmod_poly_mat_init(F, n, n, targs.modn);
    _build_input(F, targs, state);
    /* M = - F(0) */
    nmod_mat_t M;
    nmod_mat_init(M, n, n, targs.modn);
    nmod_poly_mat_get_coeff_mat(M, F, 0);
    nmod_mat_neg(M, M);
    nmod_poly_t cp;
    nmod_poly_init(cp, targs.modn);

    double FLINT_SET_BUT_UNUSED(tcpu), twall;
    TIMEIT_START;
    nmod_mat_charpoly(cp, M);
    TIMEIT_STOP_VALUES(tcpu, twall);
    flint_printf("%.2e", twall);

    nmod_poly_clear(cp);
    nmod_mat_clear(M);
    nmod_poly_mat_clear(F);
}

/*-------------------------*/
/*  main                   */
/*-------------------------*/

int main(int argc, char ** argv)
{
    flint_rand_t state;
    flint_rand_init(state);
    flint_rand_set_seed(state, 1234, 5678);

    const slong nfuns = 7;
    typedef void (*timefun) (time_args, flint_rand_t);
    const timefun funs[] = {
        time_det_iter,               // 0
        time_det_minors,             // 1
        time_det_geometric,          // 2
        time_det_flint,              // 3
        time_det_lnz,                // 4
        time_determinant,            // 5
        time_charpoly,               // 6
    };

    const char * description[] = {
        "#0  --> det_iter (Mulders-Storjohann)                      ",
        "#1  --> det_minors (expansion by minors, dim <= 12 only)   ",
        "#2  --> det_geometric (evaluation/interpolation)           ",
        "#3  --> FLINT's nmod_poly_mat_det                          ",
        "#4  --> det_lnz (Labahn-Neiger-Zhou)                       ",
        "#5  --> determinant (general interface)                    ",
        "#6  --> FLINT's nmod_mat_charpoly of M (types 2, 3 only)   ",
    };

    if (argc != 6)  // show usage
    {
        flint_printf("Usage: `%s [fun] [nbits] [dim] [deg] [type]`\n", argv[0]);
        flint_printf("   No argument shows this help.\n");
        flint_printf("   - fun: id number of the timed function (see list below), -1 times all\n");
        flint_printf("   - nbits: number of bits in [2..64] for the modulus, chosen as nextprime(2**(nbits-1))\n");
        flint_printf("   - dim: input matrix is dim x dim\n");
        flint_printf("   - deg: input matrix has degree < deg (for types 0, 1, 4, 5)\n");
        flint_printf("   - type: 0 random, 1 product U * diag(1, .., 1, x^2+x, .., x^2+x) * V,\n");
        flint_printf("           2 characteristic matrix x I - M, 3 same with M derogatory,\n");
        flint_printf("           4 row 0 of degree < deg and other rows of degree < 2,\n");
        flint_printf("           5 last column of degree < deg and other columns of degree < 2\n");
        flint_printf("\nAvailable functions:\n");
        for (slong j = 0; j < nfuns; j++)
            flint_printf("   %s\n", description[j]);
        return 0;
    }

    const slong ifun  = atoi(argv[1]);
    const slong bits  = atoi(argv[2]);
    const slong dim   = atoi(argv[3]);
    const slong deg   = atoi(argv[4]);
    const slong itype = atoi(argv[5]);
    const ulong n = n_nextprime(UWORD(1) << (bits-1), 0);
    time_args targs = {dim, deg, itype, n};

    if (itype < 0 || itype > 5)
    {
        flint_printf("Unknown input type %wd (see the usage, without arguments)\n", itype);
        return 1;
    }

    if (ifun >= 0 && ifun < nfuns)
    {
        flint_printf("bits fun dim  deg     type\n");
        flint_printf("%-5ld#%-3ld%-5ld%-8ld%-6ld", bits, ifun, dim, deg, itype);
        funs[ifun](targs, state);
        flint_printf("\n");
    }
    else if (ifun == -1)
    {
        flint_printf("bits dim  deg     type     iter   minors     geom    flint      lnz  general charpoly\n");
        flint_printf("%-5ld%-5ld%-8ld%-6ld", bits, dim, deg, itype);
        for (slong i = 0; i < nfuns; i++)
        {
            funs[i](targs, state);
            flint_printf(" ");
        }
        flint_printf("\n");
    }

    flint_rand_clear(state);
    return 0;
}
