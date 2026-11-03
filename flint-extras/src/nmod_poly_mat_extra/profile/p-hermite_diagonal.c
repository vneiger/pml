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
#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_extra.h"  /* for det_iter */

/* Input types are:
 * 0 -> random n x n matrix of degree < deg (generic: the diagonal is
 *      (det, 1, ..., 1))
 * 1 -> U * D * V with U, V random of degree < deg/2 and D diagonal with
 *      D[i,i] = 1 for i < n/2 and D[i,i] = x^2+x otherwise (nontrivial
 *      Hermite diagonal, about half of the diagonal entries are nontrivial)
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

static void _charmat(nmod_poly_mat_t F, const nmod_mat_t M)
{
    const slong n = F->r;
    nmod_poly_mat_zero(F);
    for (slong i = 0; i < n; i++)
    {
        for (slong j = 0; j < n; j++)
            nmod_poly_set_coeff_ui(nmod_poly_mat_entry(F, i, j), 0,
                                   nmod_neg(nmod_mat_entry(M, i, j), M->mod));
        nmod_poly_set_coeff_ui(nmod_poly_mat_entry(F, i, i), 1, 1);
    }
}

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
        /* U <- U * D */
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
        nmod_mat_t M, P, Pinv;
        nmod_mat_init(M, n, n, targs.modn);
        if (targs.itype == 2)
            nmod_mat_rand(M, state);
        else
        {
            const slong b = n / 2;
            nmod_mat_t B;
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
            /* random similarity */
            nmod_mat_init(P, n, n, targs.modn);
            nmod_mat_init(Pinv, n, n, targs.modn);
            do
                nmod_mat_rand(P, state);
            while (!nmod_mat_inv(Pinv, P));
            nmod_mat_t T;
            nmod_mat_init(T, n, n, targs.modn);
            nmod_mat_mul(T, P, M);
            nmod_mat_mul(M, T, Pinv);
            nmod_mat_clear(T);
            nmod_mat_clear(P);
            nmod_mat_clear(Pinv);
        }
        _charmat(F, M);
        nmod_mat_clear(M);
    }
}

/* hnf_diagonal_lnz with given thresholds */
static void _time_hnf_diagonal_lnz(time_args targs, flint_rand_t state,
                                   slong wp_thres, slong iter_dim, slong iter_thres)
{
    const slong n = targs.dim;
    nmod_poly_mat_t F;
    nmod_poly_mat_init(F, n, n, targs.modn);
    _build_input(F, targs, state);
    nmod_poly_struct * diag = flint_malloc(FLINT_MAX(1, n) * sizeof(nmod_poly_struct));
    for (slong i = 0; i < n; i++)
        nmod_poly_init(diag + i, targs.modn);

    double FLINT_SET_BUT_UNUSED(tcpu), twall;
    TIMEIT_START;
    _nmod_poly_mat_hnf_diagonal_lnz(diag, F, wp_thres, iter_dim, iter_thres);
    TIMEIT_STOP_VALUES(tcpu, twall);
    flint_printf("%.2e", twall);

    for (slong i = 0; i < n; i++)
        nmod_poly_clear(diag + i);
    flint_free(diag);
    nmod_poly_mat_clear(F);
}

/* default thresholds */
static void time_hnf_diagonal_lnz(time_args targs, flint_rand_t state)
{
    _time_hnf_diagonal_lnz(targs, state, NMOD_POLY_MAT_ROW_BASIS_WEAK_POPOV_THRES,
                           NMOD_POLY_MAT_HNF_DIAGONAL_ITER_DIM,
                           NMOD_POLY_MAT_HNF_DIAGONAL_ITER_NDEG);
}

/* recursion only: Mulders-Storjohann base case disabled */
static void time_hnf_diagonal_lnz_rec(time_args targs, flint_rand_t state)
{
    _time_hnf_diagonal_lnz(targs, state, NMOD_POLY_MAT_ROW_BASIS_WEAK_POPOV_THRES,
                           NMOD_POLY_MAT_HNF_DIAGONAL_ITER_DIM, -1);
}

/* Zhou-Labahn row bases only (default base case) */
static void time_hnf_diagonal_lnz_zl(time_args targs, flint_rand_t state)
{
    _time_hnf_diagonal_lnz(targs, state, -1, NMOD_POLY_MAT_HNF_DIAGONAL_ITER_DIM,
                           NMOD_POLY_MAT_HNF_DIAGONAL_ITER_NDEG);
}

/* Mulders-Storjohann base case whenever dim <= NMOD_POLY_MAT_HNF_DIAGONAL_ITER_DIM */
static void time_hnf_diagonal_lnz_iter(time_args targs, flint_rand_t state)
{
    _time_hnf_diagonal_lnz(targs, state, NMOD_POLY_MAT_ROW_BASIS_WEAK_POPOV_THRES,
                           NMOD_POLY_MAT_HNF_DIAGONAL_ITER_DIM, WORD_MAX);
}

/* one recursion step at the top level, default thresholds below: compared
 * with the default where the latter uses the base case, this tells whether
 * the base case thresholds should be lowered */
static void time_hnf_diagonal_lnz_step(time_args targs, flint_rand_t state)
{
    _time_hnf_diagonal_lnz(targs, state, NMOD_POLY_MAT_ROW_BASIS_WEAK_POPOV_THRES,
                           FLINT_MIN(targs.dim - 1, NMOD_POLY_MAT_HNF_DIAGONAL_ITER_DIM),
                           NMOD_POLY_MAT_HNF_DIAGONAL_ITER_NDEG);
}

/* existing upper row-wise HNF algorithms, applied to the matrix with reversed
 * rows and columns (whose diagonal gives the lower HNF diagonal) */
#define TIME_HNF(name, call)                                         \
static void time_##name(time_args targs, flint_rand_t state)              \
{                                                                    \
    const slong n = targs.dim;                                       \
    nmod_poly_mat_t F, H;                                            \
    nmod_poly_mat_init(F, n, n, targs.modn);                         \
    nmod_poly_mat_init(H, n, n, targs.modn);                         \
    _build_input(F, targs, state);                                   \
    nmod_poly_mat_invert_rows(F, NULL);                              \
    nmod_poly_mat_invert_columns(F, NULL);                           \
    slong * pivind = flint_malloc(n * sizeof(slong));                \
    slong * mrp = flint_malloc(n * sizeof(slong));                   \
                                                                     \
    double FLINT_SET_BUT_UNUSED(tcpu), twall;                        \
    TIMEIT_START;                                                    \
    nmod_poly_mat_set(H, F);                                         \
    call;                                                            \
    TIMEIT_STOP_VALUES(tcpu, twall);                                 \
    flint_printf("%.2e", twall);                                     \
                                                                     \
    flint_free(pivind);                                              \
    flint_free(mrp);                                                 \
    nmod_poly_mat_clear(F);                                          \
    nmod_poly_mat_clear(H);                                          \
}

TIME_HNF(hnf_revlex_xgcd, nmod_poly_mat_hnf_ur_revlex_xgcd(H, NULL, pivind, mrp))
TIME_HNF(hnf_lex_xgcd, nmod_poly_mat_hnf_ur_lex_xgcd(H, NULL, pivind, mrp))
TIME_HNF(hnf_maxdeg_atomic, nmod_poly_mat_hnf_ur_maxdeg_atomic(H, NULL, pivind))
TIME_HNF(hnf_matrixgcd_iter, nmod_poly_mat_hnf_ur_matrixgcd_iter(H, NULL, pivind, NULL, NULL))
TIME_HNF(det_iter, nmod_poly_t det; nmod_poly_init(det, targs.modn); nmod_poly_mat_det_iter(det, H); nmod_poly_clear(det))

/*-------------------------*/
/*  main                   */
/*-------------------------*/

int main(int argc, char ** argv)
{
    flint_rand_t state;
    flint_rand_init(state);
    flint_rand_set_seed(state, 1234, 5678);

    const slong nfuns = 10;
    typedef void (*timefun) (time_args, flint_rand_t);
    const timefun funs[] = {
        time_hnf_diagonal_lnz,       // 0
        time_hnf_diagonal_lnz_rec,   // 1
        time_hnf_diagonal_lnz_zl,    // 2
        time_hnf_diagonal_lnz_iter,  // 3
        time_det_iter,               // 4
        time_hnf_matrixgcd_iter,     // 5
        time_hnf_revlex_xgcd,        // 6
        time_hnf_lex_xgcd,           // 7
        time_hnf_maxdeg_atomic,      // 8
        time_hnf_diagonal_lnz_step,  // 9
    };

    const char * description[] = {
        "#0  --> hnf_diagonal_lnz (default thresholds)              ",
        "#1  --> hnf_diagonal_lnz, recursion only (no base case)    ",
        "#2  --> hnf_diagonal_lnz, Zhou-Labahn row bases only       ",
        "#3  --> hnf_diagonal_lnz, base case whenever dim <= ITER_DIM",
        "#4  --> det_iter (Mulders-Storjohann determinant)          ",
        "#5  --> HNF, matrixgcd_iter (Mulders-Storjohann)           ",
        "#6  --> HNF, revlex_xgcd (Bradley-like)                    ",
        "#7  --> HNF, lex_xgcd                                      ",
        "#8  --> HNF, maxdeg_atomic (Rosser-like)                   ",
        "#9  --> hnf_diagonal_lnz, one recursion step at the top level",
    };

    if (argc != 6)  // show usage
    {
        flint_printf("Usage: `%s [fun] [nbits] [dim] [deg] [type]`\n", argv[0]);
        flint_printf("   No argument shows this help.\n");
        flint_printf("   - fun: id number of the timed function (see list below),\n");
        flint_printf("          -1 times #0 to #5, -2 times all\n");
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

    if (ifun >= 0 && ifun < nfuns)
    {
        flint_printf("bits fun dim  deg     type\n");
        flint_printf("%-5ld#%-3ld%-5ld%-8ld%-6ld", bits, ifun, dim, deg, itype);
        funs[ifun](targs, state);
        flint_printf("\n");
    }
    else if (ifun == -1 || ifun == -2)
    {
        const slong nf = (ifun == -1) ? 6 : nfuns;
        flint_printf("bits dim  deg     type   hd_lnz   hd_rec    hd_zl  hd_iter det_iter hnf_mgcd  hnf_rlx  hnf_lex hnf_rssr  hd_step\n");
        flint_printf("%-5ld%-5ld%-8ld%-6ld", bits, dim, deg, itype);
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
