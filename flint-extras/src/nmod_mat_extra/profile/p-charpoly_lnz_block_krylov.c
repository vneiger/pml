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
#include <flint/nmod.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_poly.h>
#include <flint/profiler.h>
#include <flint/ulong_extras.h>

#include "nmod_mat_extra.h"

/* Input types are:
 * 0 -> random n x n matrix
 * 1 -> similar to a block diagonal matrix diag(B, B) with B random
 *      (derogatory: two nontrivial invariant factors)
 * 2 -> upper triangular (the Krylov reduction fails at once)
 * 3 -> block diagonal diag(B, C) with B, C random (the Krylov reduction
 *      fails after computing the Krylov matrix, when n/2 > k)
 */

#define NLEVELS_MAX 9  /* nlevels = 1, 2, 4, ..., 256 */

static void _build_input(nmod_mat_t A, slong itype, flint_rand_t state)
{
    const slong n = A->r;
    const ulong p = A->mod.n;
    nmod_mat_zero(A);
    if (itype == 0)
        nmod_mat_rand(A, state);
    else if (itype == 2)
    {
        for (slong i = 0; i < n; i++)
            for (slong j = i; j < n; j++)
                nmod_mat_entry(A, i, j) = n_randint(state, p);
    }
    else
    {
        const slong b = n / 2;
        for (slong i = 0; i < n; i++)
            for (slong j = 0; j < n; j++)
                if ((i < b) == (j < b))
                    nmod_mat_entry(A, i, j) = n_randint(state, p);
        if (itype == 1)
        {
            for (slong i = 0; i < b; i++)
                for (slong j = 0; j < b; j++)
                    nmod_mat_entry(A, n - b + i, n - b + j) = nmod_mat_entry(A, i, j);
            nmod_mat_t P, Pinv, T;
            nmod_mat_init(P, n, n, p);
            nmod_mat_init(Pinv, n, n, p);
            nmod_mat_init(T, n, n, p);
            do
                nmod_mat_rand(P, state);
            while (!nmod_mat_inv(Pinv, P));
            nmod_mat_mul(T, P, A);
            nmod_mat_mul(A, T, Pinv);
            nmod_mat_clear(P);
            nmod_mat_clear(Pinv);
            nmod_mat_clear(T);
        }
    }
}

/* ifun: 0 -> nmod_mat_charpoly (FLINT), 1 + j ->
 * _nmod_mat_charpoly_lnz_block_krylov with nlevels = 2^j */
static void time_fun(slong ifun, const nmod_mat_t A)
{
    nmod_poly_t cp;
    nmod_poly_init(cp, A->mod.n);
    int ok = 1;
    const slong nlevels = (ifun >= 1) ? (WORD(1) << (ifun - 1)) : 0;

    if (ifun >= 1 && nlevels > A->r && (nlevels / 2) >= A->r)
    {
        /* same as a smaller number of levels */
        flint_printf("      - ");
        nmod_poly_clear(cp);
        return;
    }

    double FLINT_SET_BUT_UNUSED(tcpu), twall;
    TIMEIT_START;
    if (ifun == 0)
        nmod_mat_charpoly(cp, A);
    else
        ok = _nmod_mat_charpoly_lnz_block_krylov(cp, A, nlevels);
    TIMEIT_STOP_VALUES(tcpu, twall);

    if (ok)
        flint_printf("%.2e", twall);
    else
        flint_printf("fail%.0e", twall);
    nmod_poly_clear(cp);
}

int main(int argc, char ** argv)
{
    flint_rand_t state;
    flint_rand_init(state);
    flint_rand_set_seed(state, 1234, 5678);

    const slong nfuns = 1 + NLEVELS_MAX;

    if (argc != 5)  // show usage
    {
        flint_printf("Usage: `%s [fun] [nbits] [dim] [type]`\n", argv[0]);
        flint_printf("   No argument shows this help.\n");
        flint_printf("   - fun: id number of the timed function (see list below), -1 times all\n");
        flint_printf("   - nbits: number of bits in [2..64] for the modulus, chosen as nextprime(2**(nbits-1))\n");
        flint_printf("   - dim: input matrix is dim x dim\n");
        flint_printf("   - type: 0 random, 1 derogatory (similar to diag(B, B)), 2 upper triangular,\n");
        flint_printf("           3 block diagonal diag(B, C)\n");
        flint_printf("\nAvailable functions:\n");
        flint_printf("   #0  --> FLINT's nmod_mat_charpoly\n");
        for (slong j = 0; j < NLEVELS_MAX; j++)
            flint_printf("   #%-2wd --> _nmod_mat_charpoly_lnz_block_krylov with nlevels = %wd\n", 1 + j, WORD(1) << j);
        flint_printf("\n'fail': the Krylov reduction failed (time to detect it); '-': same as fewer levels\n");
        return 0;
    }

    const slong ifun  = atoi(argv[1]);
    const slong bits  = atoi(argv[2]);
    const slong dim   = atoi(argv[3]);
    const slong itype = atoi(argv[4]);
    const ulong p = n_nextprime(UWORD(1) << (bits-1), 0);

    if (itype < 0 || itype > 3 || ifun < -1 || ifun >= nfuns)
    {
        flint_printf("Unknown function or input type (see the usage, without arguments)\n");
        return 1;
    }

    nmod_mat_t A;
    nmod_mat_init(A, dim, dim, p);
    _build_input(A, itype, state);

    if (ifun >= 0)
    {
        flint_printf("bits fun dim  type\n");
        flint_printf("%-5wd#%-3wd%-5wd%-5wd", bits, ifun, dim, itype);
        time_fun(ifun, A);
        flint_printf("\n");
    }
    else
    {
        flint_printf("bits dim  type    flint |  nlevels=1        2        4        8       16       32       64      128      256\n");
        flint_printf("%-5wd%-5wd%-5wd", bits, dim, itype);
        for (slong i = 0; i < nfuns; i++)
        {
            flint_printf(" ");
            time_fun(i, A);
            if (i == 0)
                flint_printf(" |");
        }
        flint_printf("\n");
    }

    nmod_mat_clear(A);
    flint_rand_clear(state);
    return 0;
}

#undef NLEVELS_MAX
