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

#include "nmod_poly_mat_forms.h"
#include "nmod_poly_mat_extra.h"  /* for det_iter */

/* Timings of the Mulders-Storjohann iteration (weak Popov form) and of the
 * functions built on it: determinant (det_iter) and Hermite form
 * (hnf_ur_matrixgcd_iter).
 *
 * Input types are:
 * 0 -> random rdim x cdim matrix of degree < deg
 * 1 -> W * B with W random rdim x cdim of degree < deg/2 and B random
 *      cdim x cdim of degree < deg - deg/2
 * 2 -> W * B as above, but with B of rank cdim/2 (the input has rank cdim/2)
 * 3 -> x I - M with M a random constant rdim x rdim matrix (requires
 *      rdim == cdim, deg is ignored)
 *
 * Each timing includes a copy of the input, since the functions work in
 * place. Functions that require a square nonsingular input (det_iter,
 * hnf_ur_matrixgcd_iter) print "-" on other inputs.
 */

typedef struct
{
    slong rdim;  /* row dimension */
    slong cdim;  /* column dimension */
    slong deg;   /* degree */
    slong itype; /* input type */
    ulong modn;  /* modulus */
}
time_args;

static void _build_input(nmod_poly_mat_t pmat, time_args targs, flint_rand_t state)
{
    const slong m = targs.rdim;
    const slong n = targs.cdim;
    const slong deg = targs.deg;
    if (targs.itype == 0)
        nmod_poly_mat_rand(pmat, state, deg);
    else if (targs.itype == 3)
    {
        nmod_mat_t M;
        nmod_mat_init(M, m, m, targs.modn);
        nmod_mat_rand(M, state);
        nmod_mat_neg(M, M);
        nmod_poly_mat_set_nmod_mat(pmat, M);
        for (slong i = 0; i < m; i++)
            nmod_poly_set_coeff_ui(nmod_poly_mat_entry(pmat, i, i), 1, 1);
        nmod_mat_clear(M);
    }
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

/* whether the input is square and (generically) nonsingular */
static int _square_nonsingular(time_args targs)
{
    return targs.rdim == targs.cdim && targs.itype != 2;
}

/* weak Popov form, without transformation */
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

/* weak Popov form, with the unimodular transformation */
void time_weak_popov_tsf(time_args targs, flint_rand_t state)
{
    nmod_poly_mat_t pmat, copy, tsf;
    nmod_poly_mat_init(pmat, targs.rdim, targs.cdim, targs.modn);
    nmod_poly_mat_init(copy, targs.rdim, targs.cdim, targs.modn);
    nmod_poly_mat_init(tsf, targs.rdim, targs.rdim, targs.modn);
    _build_input(pmat, targs, state);
    slong * pivind = flint_malloc(FLINT_MAX(1, targs.rdim) * sizeof(slong));

    double FLINT_SET_BUT_UNUSED(tcpu), twall;
    TIMEIT_START;
    nmod_poly_mat_set(copy, pmat);
    nmod_poly_mat_one(tsf);
    nmod_poly_mat_weak_popov_iter(copy, NULL, tsf, pivind, NULL, ROW_LOWER);
    TIMEIT_STOP_VALUES(tcpu, twall);
    flint_printf("%.2e", twall);

    flint_free(pivind);
    nmod_poly_mat_clear(pmat);
    nmod_poly_mat_clear(copy);
    nmod_poly_mat_clear(tsf);
}

/* determinant (Mulders-Storjohann) */
void time_det_iter(time_args targs, flint_rand_t state)
{
    if (! _square_nonsingular(targs))
    {
        flint_printf("%8s", "-");
        return;
    }

    nmod_poly_mat_t pmat, copy;
    nmod_poly_t det;
    nmod_poly_mat_init(pmat, targs.rdim, targs.cdim, targs.modn);
    nmod_poly_mat_init(copy, targs.rdim, targs.cdim, targs.modn);
    nmod_poly_init(det, targs.modn);
    _build_input(pmat, targs, state);

    double FLINT_SET_BUT_UNUSED(tcpu), twall;
    TIMEIT_START;
    nmod_poly_mat_set(copy, pmat);
    nmod_poly_mat_det_iter(det, copy);
    TIMEIT_STOP_VALUES(tcpu, twall);
    flint_printf("%.2e", twall);

    nmod_poly_mat_clear(pmat);
    nmod_poly_mat_clear(copy);
    nmod_poly_clear(det);
}

/* upper Hermite form (Mulders-Storjohann) */
void time_hnf_matrixgcd(time_args targs, flint_rand_t state)
{
    if (! _square_nonsingular(targs))
    {
        flint_printf("%8s", "-");
        return;
    }

    nmod_poly_mat_t pmat, copy;
    nmod_poly_mat_init(pmat, targs.rdim, targs.cdim, targs.modn);
    nmod_poly_mat_init(copy, targs.rdim, targs.cdim, targs.modn);
    _build_input(pmat, targs, state);
    slong * pivind = flint_malloc(FLINT_MAX(1, targs.rdim) * sizeof(slong));

    double FLINT_SET_BUT_UNUSED(tcpu), twall;
    TIMEIT_START;
    nmod_poly_mat_set(copy, pmat);
    nmod_poly_mat_hnf_ur_matrixgcd_iter(copy, NULL, pivind, NULL, NULL);
    TIMEIT_STOP_VALUES(tcpu, twall);
    flint_printf("%.2e", twall);

    flint_free(pivind);
    nmod_poly_mat_clear(pmat);
    nmod_poly_mat_clear(copy);
}

/*-------------------------*/
/*  main                   */
/*-------------------------*/

int main(int argc, char ** argv)
{
    /* fixed seed: two builds (e.g. before and after a change) are timed on
     * the same inputs, as long as they use the same FLINT */
    flint_rand_t state;
    flint_rand_init(state);

    const slong nfuns = 4;
    typedef void (*timefun) (time_args, flint_rand_t);
    const timefun funs[] = {
        time_weak_popov,             // 0
        time_weak_popov_tsf,         // 1
        time_det_iter,               // 2
        time_hnf_matrixgcd,          // 3
    };

    const char * description[] = {
        "#0  --> weak Popov form                              ",
        "#1  --> weak Popov form, with transformation         ",
        "#2  --> determinant, det_iter                        ",
        "#3  --> upper Hermite form, hnf_ur_matrixgcd_iter    ",
    };

    if (argc != 7)  // show usage
    {
        flint_printf("Usage: `%s [fun] [nbits] [rdim] [cdim] [deg] [type]`\n", argv[0]);
        flint_printf("   No argument shows this help.\n");
        flint_printf("   - fun: id number of the timed function (see list below), -1 times all\n");
        flint_printf("   - nbits: number of bits in [2..64] for the modulus, chosen as nextprime(2**(nbits-1))\n");
        flint_printf("   - rdim, cdim: input matrix is rdim x cdim\n");
        flint_printf("   - deg: input matrix has degree < deg\n");
        flint_printf("   - type: 0 for random input, 1 for product W * B with B random square,\n");
        flint_printf("           2 for such a product with B of rank cdim/2, 3 for x I - M\n");
        flint_printf("           with M random constant (requires rdim == cdim, deg is ignored)\n");
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

    if (itype == 3 && rdim != cdim)
    {
        flint_printf("type 3 (x I - M) requires rdim == cdim\n");
        return 1;
    }

    if (ifun >= 0 && ifun < nfuns)
    {
        flint_printf("bits fun rdim cdim deg type\n");
        flint_printf("%-5ld#%-3ld%-5ld%-5ld%-8ld%-2ld", bits, ifun, rdim, cdim, deg, itype);
        funs[ifun](targs, state);
        flint_printf("\n");
    }
    else if (ifun == -1)
    {
        flint_printf("bits rdim cdim deg type   wkpopov  wkp-tsf  det_iter  hnf_ms\n");
        flint_printf("%-5ld%-5ld%-5ld%-8ld%-6ld", bits, rdim, cdim, deg, itype);
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
