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
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>
#include <flint/profiler.h>
#include <flint/ulong_extras.h>

#include "nmod_poly_mat_forms.h"
#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_row_basis.h"

/* Input types are:
 * 0 -> random m x n matrix of degree < deg (generic: when m > n, its row
 *      space is the whole of K[x]^{1 x n})
 * 1 -> W * B with W random m x n of degree < deg/2 and B random n x n of
 *      degree < deg - deg/2 (when m >= n, B is generically a row basis)
 * 2 -> W * B as above, but with B of rank n/2 (the input has rank n/2)
 */

typedef struct
{
    slong rdim;  /* row dimension */
    slong cdim;  /* column dimension */
    slong deg;   /* degree */
    slong itype; /* input type */
    slong modn;  /* modulus */
}
time_args;

static void _build_input(nmod_poly_mat_t pmat, time_args targs, flint_rand_t state)
{
    const slong m = targs.rdim;
    const slong n = targs.cdim;
    const slong deg = targs.deg;
    if (targs.itype == 0)
        nmod_poly_mat_rand(pmat, state, deg);
    else
    {
        const slong k = (targs.itype == 1) ? n : FLINT_MAX(1, n/2);
        nmod_poly_mat_t W, B;
        nmod_poly_mat_init(W, m, k, targs.modn);
        nmod_poly_mat_init(B, k, n, targs.modn);
        nmod_poly_mat_rand(W, state, FLINT_MAX(1, deg/2));
        nmod_poly_mat_rand(B, state, FLINT_MAX(1, deg - deg/2));
        nmod_poly_mat_mul(pmat, W, B);
        nmod_poly_mat_clear(W);
        nmod_poly_mat_clear(B);
    }
}

void time_row_basis(time_args targs, flint_rand_t state)
{
    nmod_poly_mat_t pmat, rbas;
    nmod_poly_mat_init(pmat, targs.rdim, targs.cdim, targs.modn);
    nmod_poly_mat_init(rbas, FLINT_MIN(targs.rdim, targs.cdim), targs.cdim, targs.modn);
    _build_input(pmat, targs, state);

    double FLINT_SET_BUT_UNUSED(tcpu), twall;
    TIMEIT_START;
    nmod_poly_mat_row_basis(rbas, pmat);
    TIMEIT_STOP_VALUES(tcpu, twall);
    flint_printf("%.2e", twall);

    nmod_poly_mat_clear(pmat);
    nmod_poly_mat_clear(rbas);
}

void time_row_basis_zl(time_args targs, flint_rand_t state)
{
    nmod_poly_mat_t pmat, rbas;
    nmod_poly_mat_init(pmat, targs.rdim, targs.cdim, targs.modn);
    nmod_poly_mat_init(rbas, FLINT_MIN(targs.rdim, targs.cdim), targs.cdim, targs.modn);
    _build_input(pmat, targs, state);

    double FLINT_SET_BUT_UNUSED(tcpu), twall;
    TIMEIT_START;
    nmod_poly_mat_row_basis_zl(rbas, NULL, pmat);
    TIMEIT_STOP_VALUES(tcpu, twall);
    flint_printf("%.2e", twall);

    nmod_poly_mat_clear(pmat);
    nmod_poly_mat_clear(rbas);
}

void time_row_basis_zl_nocert(time_args targs, flint_rand_t state)
{
    nmod_poly_mat_t pmat, rbas;
    nmod_poly_mat_init(pmat, targs.rdim, targs.cdim, targs.modn);
    nmod_poly_mat_init(rbas, FLINT_MIN(targs.rdim, targs.cdim), targs.cdim, targs.modn);
    _build_input(pmat, targs, state);

    double FLINT_SET_BUT_UNUSED(tcpu), twall;
    TIMEIT_START;
    _nmod_poly_mat_row_basis_zl(rbas, NULL, pmat, NULL, NULL, NULL, 0, 0);
    TIMEIT_STOP_VALUES(tcpu, twall);
    flint_printf("%.2e", twall);

    nmod_poly_mat_clear(pmat);
    nmod_poly_mat_clear(rbas);
}

void time_row_basis_zl_lfac(time_args targs, flint_rand_t state)
{
    nmod_poly_mat_t pmat, rbas, lfac;
    nmod_poly_mat_init(pmat, targs.rdim, targs.cdim, targs.modn);
    nmod_poly_mat_init(rbas, FLINT_MIN(targs.rdim, targs.cdim), targs.cdim, targs.modn);
    nmod_poly_mat_init(lfac, targs.rdim, FLINT_MIN(targs.rdim, targs.cdim), targs.modn);
    _build_input(pmat, targs, state);

    double FLINT_SET_BUT_UNUSED(tcpu), twall;
    TIMEIT_START;
    nmod_poly_mat_row_basis_zl(rbas, lfac, pmat);
    TIMEIT_STOP_VALUES(tcpu, twall);
    flint_printf("%.2e", twall);

    nmod_poly_mat_clear(pmat);
    nmod_poly_mat_clear(rbas);
    nmod_poly_mat_clear(lfac);
}

void time_row_basis_zl_blocks(time_args targs, flint_rand_t state)
{
    nmod_poly_mat_t pmat, rbas;
    nmod_poly_mat_init(pmat, targs.rdim, targs.cdim, targs.modn);
    nmod_poly_mat_init(rbas, FLINT_MIN(targs.rdim, targs.cdim), targs.cdim, targs.modn);
    _build_input(pmat, targs, state);

    double FLINT_SET_BUT_UNUSED(tcpu), twall;
    TIMEIT_START;
    _nmod_poly_mat_row_basis_zl_blocks(rbas, pmat, targs.cdim);
    TIMEIT_STOP_VALUES(tcpu, twall);
    flint_printf("%.2e", twall);

    nmod_poly_mat_clear(pmat);
    nmod_poly_mat_clear(rbas);
}

/* row basis as the nonzero rows of a weak Popov form (Mulders-Storjohann) */
void time_weak_popov(time_args targs, flint_rand_t state)
{
    nmod_poly_mat_t pmat, copy;
    nmod_poly_mat_init(pmat, targs.rdim, targs.cdim, targs.modn);
    nmod_poly_mat_init(copy, targs.rdim, targs.cdim, targs.modn);
    _build_input(pmat, targs, state);
    slong * pivind = flint_malloc(FLINT_MAX(1, targs.rdim) * sizeof(slong));

    double FLINT_SET_BUT_UNUSED(tcpu), twall;
    TIMEIT_START;
    nmod_poly_mat_set(copy, pmat);
    nmod_poly_mat_weak_popov_iter(copy, NULL, NULL, pivind, NULL, ROW_LOWER);
    TIMEIT_STOP_VALUES(tcpu, twall);
    flint_printf("%.2e", twall);

    flint_free(pivind);
    nmod_poly_mat_clear(pmat);
    nmod_poly_mat_clear(copy);
}

/* row basis as the nonzero rows of the upper Hermite normal form */
void time_hnf(time_args targs, flint_rand_t state)
{
    nmod_poly_mat_t pmat, copy;
    nmod_poly_mat_init(pmat, targs.rdim, targs.cdim, targs.modn);
    nmod_poly_mat_init(copy, targs.rdim, targs.cdim, targs.modn);
    _build_input(pmat, targs, state);
    slong * pivind = flint_malloc(FLINT_MAX(1, FLINT_MIN(targs.rdim, targs.cdim)) * sizeof(slong));
    slong * mrp = flint_malloc(FLINT_MAX(1, targs.rdim) * sizeof(slong));

    double FLINT_SET_BUT_UNUSED(tcpu), twall;
    TIMEIT_START;
    nmod_poly_mat_set(copy, pmat);
    nmod_poly_mat_hnf_ur_revlex_xgcd(copy, NULL, pivind, mrp);
    TIMEIT_STOP_VALUES(tcpu, twall);
    flint_printf("%.2e", twall);

    flint_free(pivind);
    flint_free(mrp);
    nmod_poly_mat_clear(pmat);
    nmod_poly_mat_clear(copy);
}

/*-------------------------*/
/*  main                   */
/*-------------------------*/

int main(int argc, char ** argv)
{
    flint_rand_t state;
    flint_rand_init(state);
    flint_rand_set_seed(state, time(NULL), time(NULL)+129384125L);

    const slong nfuns = 7;
    typedef void (*timefun) (time_args, flint_rand_t);
    const timefun funs[] = {
        time_row_basis,              // 0
        time_row_basis_zl,           // 1
        time_row_basis_zl_nocert,    // 2
        time_row_basis_zl_lfac,      // 3
        time_row_basis_zl_blocks,    // 4
        time_weak_popov,             // 5
        time_hnf,                    // 6
    };

    const char * description[] = {
        "#0  --> row_basis (general interface)                ",
        "#1  --> row_basis_zl                                 ",
        "#2  --> row_basis_zl, without unimodularity check    ",
        "#3  --> row_basis_zl, with left factor               ",
        "#4  --> row_basis_zl_blocks, blocks of cdim rows     ",
        "#5  --> weak Popov form (Mulders-Storjohann)         ",
        "#6  --> upper Hermite form (revlex, xgcd)            ",
    };

    if (argc != 7)  // show usage
    {
        flint_printf("Usage: `%s [fun] [nbits] [rdim] [cdim] [deg] [type]`\n", argv[0]);
        flint_printf("   No argument shows this help.\n");
        flint_printf("   - fun: id number of the timed function (see list below), -1 times all\n");
        flint_printf("          but the last one (Hermite form), -2 times all\n");
        flint_printf("   - nbits: number of bits in [2..64] for the modulus, chosen as nextprime(2**(nbits-1))\n");
        flint_printf("   - rdim, cdim: input matrix is rdim x cdim\n");
        flint_printf("   - deg: input matrix has degree < deg\n");
        flint_printf("   - type: 0 for random input (generic), 1 for product W * B with B\n");
        flint_printf("           random square (nontrivial row basis), 2 for such a product\n");
        flint_printf("           with B of rank cdim/2\n");
        flint_printf("\nAvailable functions:\n");
        for (slong j = 0; j < nfuns; j++)
            flint_printf("   %s\n", description[j]);
        return 0;
    }

    const slong ifun  = atoi(argv[1]);
    const slong bits  = atoi(argv[2]);
    const slong rdim  = atoi(argv[3]);
    const slong cdim  = atoi(argv[4]);
    const slong deg   = atoi(argv[5]);
    const slong itype = atoi(argv[6]);
    const ulong n = n_nextprime(UWORD(1) << (bits-1), 0);
    time_args targs = {rdim, cdim, deg, itype, n};

    if (ifun >= 0 && ifun < nfuns)
    {
        flint_printf("bits fun rdim cdim deg type\n");
        flint_printf("%-5ld#%-3ld%-5ld%-5ld%-8ld%-2ld", bits, ifun, rdim, cdim, deg, itype);
        funs[ifun](targs, state);
        flint_printf("\n");
    }
    else if (ifun == -1 || ifun == -2)
    {
        const slong nf = (ifun == -1) ? nfuns - 1 : nfuns;
        flint_printf("bits rdim cdim deg type      rbas       zl   zl-ncrt  zl-lfac zl-block  wkpopov      hnf\n");
        flint_printf("%-5ld%-5ld%-5ld%-8ld%-6ld", bits, rdim, cdim, deg, itype);
        for (slong i = 0; i < nf; i++)
        {
            funs[i](targs, state);
            flint_printf(" ");
        }
        flint_printf("\n");
    }

    flint_rand_clear(state);
    return 0;
}
