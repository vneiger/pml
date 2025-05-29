/*
    Copyright 2024 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>  // for atoi

#include "flint/ulong_extras.h"
#include "flint/profiler.h"
#include "flint/nmod.h"
#include "flint/nmod_vec.h"

#include "nmod32_vec.h"

// utility (nmod vec uniform random)
static inline
void _nmod32_vec_rand(n32_ptr vec, flint_rand_t state, slong len, nmod_t mod)
{
    for (slong i = 0; i < len; i++)
        vec[i] = n_randint(state, mod.n);
}

/*------------------------------------*/
/* direct: dot / dot_rev / dot expr   */
/*------------------------------------*/

void time_dot_split(ulong len, ulong n, flint_rand_t state)
{
    nmod_t mod;
    nmod_init(&mod, n);
    ulong pow2_precomp;
    NMOD_RED(pow2_precomp, (UWORD(1) << DOT_SPLIT_BITS), mod);

    n32_ptr v1 = _nmod32_vec_init(len);
    _nmod32_vec_rand(v1, state, len, mod);
    n32_ptr v2 = _nmod32_vec_init(len);
    _nmod32_vec_rand(v2, state, len, mod);

    // store results in volatile variable to avoid that they
    // are "optimized away" (especially for inlined part)
    volatile uint FLINT_SET_BUT_UNUSED(res);

    double FLINT_SET_BUT_UNUSED(tcpu), twall;

    TIMEIT_START
    res = _nmod32_vec_dot_split(v1, v2, len, mod, pow2_precomp);
    TIMEIT_STOP_VALUES(tcpu, twall)

    printf("%.2e", twall);

    _nmod32_vec_clear(v1);
    _nmod32_vec_clear(v2);
}

void time_dot_split_avx2(ulong len, ulong n, flint_rand_t state)
{
    nmod_t mod;
    nmod_init(&mod, n);
    ulong pow2_precomp;
    NMOD_RED(pow2_precomp, (UWORD(1) << DOT_SPLIT_BITS), mod);

    n32_ptr v1 = _nmod32_vec_init(len);
    _nmod32_vec_rand(v1, state, len, mod);
    n32_ptr v2 = _nmod32_vec_init(len);
    _nmod32_vec_rand(v2, state, len, mod);

    // store results in volatile variable to avoid that they
    // are "optimized away" (especially for inlined part)
    volatile uint FLINT_SET_BUT_UNUSED(res);

    double FLINT_SET_BUT_UNUSED(tcpu), twall;

    TIMEIT_START
    res = _nmod32_vec_dot_split_avx2(v1, v2, len, mod, pow2_precomp);
    TIMEIT_STOP_VALUES(tcpu, twall)

    printf("%.2e", twall);

    _nmod32_vec_clear(v1);
    _nmod32_vec_clear(v2);
}

void time_dot_split_avx512(ulong len, ulong n, flint_rand_t state)
{
    nmod_t mod;
    nmod_init(&mod, n);
    ulong pow2_precomp;
    NMOD_RED(pow2_precomp, (UWORD(1) << DOT_SPLIT_BITS), mod);

    n32_ptr v1 = _nmod32_vec_init(len);
    _nmod32_vec_rand(v1, state, len, mod);
    n32_ptr v2 = _nmod32_vec_init(len);
    _nmod32_vec_rand(v2, state, len, mod);

    // store results in volatile variable to avoid that they
    // are "optimized away" (especially for inlined part)
    volatile uint FLINT_SET_BUT_UNUSED(res);

    double FLINT_SET_BUT_UNUSED(tcpu), twall;

    TIMEIT_START
    res = _nmod32_vec_dot_split_avx512(v1, v2, len, mod, pow2_precomp);
    TIMEIT_STOP_VALUES(tcpu, twall)

    printf("%.2e", twall);

    _nmod32_vec_clear(v1);
    _nmod32_vec_clear(v2);
}

void time_dot_msolve_avx2(ulong len, ulong n, flint_rand_t state)
{
    nmod_t mod;
    nmod_init(&mod, n);

    n32_ptr v1 = _nmod32_vec_init(len);
    _nmod32_vec_rand(v1, state, len, mod);
    n32_ptr v2 = _nmod32_vec_init(len);
    _nmod32_vec_rand(v2, state, len, mod);

    // store results in volatile variable to avoid that they
    // are "optimized away" (especially for inlined part)
    volatile uint FLINT_SET_BUT_UNUSED(res);

    double FLINT_SET_BUT_UNUSED(tcpu), twall;

    TIMEIT_START
    res = _nmod32_vec_dot_msolve_avx2(v1, v2, len, mod.n);
    TIMEIT_STOP_VALUES(tcpu, twall)

    printf("%.2e", twall);

    _nmod32_vec_clear(v1);
    _nmod32_vec_clear(v2);
}


/*-------------------------*/
/*  main                   */
/*-------------------------*/

int main(int argc, char ** argv)
{
    flint_rand_t state;
    flint_rand_init(state);
    flint_rand_set_seed(state, time(NULL), time(NULL)+129384125L);

    // modulus bitsize
    //const slong nbits = 5;
    //const ulong bits[] = {12, 26, 27, 28, 31, 32};
    const slong nbits = 1;
    const ulong bits[] = {31};

    // vector lengths
    //const slong nlens = 20;
    //const ulong lens[] = {1, 2, 3, 4, 5, 7, 10, 15, 25, 35, 50, 100, 250, 500, 1000, 2500, 5000, 10000, 100000, 1000000};
    const slong nlens = 10;
    const ulong lens[] = {50, 100, 250, 500, 1000, 2500, 5000, 10000, 100000, 1000000};

    // bench functions
    const slong nfuns = 4;
    typedef void (*timefun) (ulong, ulong, flint_rand_t);
    const timefun funs[] = {
        time_dot_split,             // 0
        time_dot_split_avx2,        // 1
        time_dot_split_avx512,      // 2
        time_dot_msolve_avx2,       // 3
    };

    const char * description[] = {
        "#0  --> dot_split            ",
        "#1  --> dot_split_avx2       ",
        "#2  --> dot_split_avx512     ",
        "#3  --> dot_msolve_avx2      ",
    };

    if (argc == 1)  // show usage
    {
        printf("Usage: `%s [fun] [nbits] [len]`\n", argv[0]);
        printf("   Each argument is optional; no argument shows this help.\n");
        printf("   - fun: id number of the timed function (see below),\n");
        printf("          exception: fun == -1 times all available functions successively\n");
        printf("   - nbits: number of bits for the modulus, chosen as nextprime(2**(nbits-1))\n");
        printf("   - len: length of the vectors\n");
        printf("\nAvailable functions:\n");
        for (slong j = 0; j < nfuns; j++)
            printf("   %s\n", description[j]);

        return 0;
    }

    printf("#warmup... ");
    for (slong i = 0; i < 10; i++)
    {
        time_dot_split_avx2(1000000, (UWORD(1)<<20)+5, state);
        printf(" ");
    }
    printf("\n");

    if (argc == 2 && atoi(argv[1]) == -1)  // launching full suite
    {
        for (slong ifun = 0; ifun < nfuns; ifun++)
        {
            const timefun tfun = funs[ifun];

            printf("\n%s\n", description[ifun]);
            printf("#bits\\len");
            for (slong i = 0; i < nlens; i++)
                printf("%9ld", lens[i]);
            printf("\n");

            for (slong j = 0; j < nbits; j++)
            {
                const slong b = bits[j];

                printf("%-10ld", b);
                ulong n;
                if (b == 232)
                    n = UWORD(1) << 32;
                else if (b == 263)
                    n = UWORD(1) << 63;
                else
                    n = n_nextprime(UWORD(1) << (b-1), 0);
                for (slong i = 0; i < nlens; i++)
                {
                    tfun(lens[i], n, state);
                    printf(" ");
                }
                printf("\n");
            }
        }
    }
    else if (argc == 2)  // function is given
    {
        const timefun tfun = funs[atoi(argv[1])];

        printf("\n%s\n", description[atoi(argv[1])]);
        printf("#bits\\len");
        for (slong i = 0; i < nlens; i++)
            printf("%9ld", lens[i]);
        printf("\n");

        for (slong j = 0; j < nbits; j++)
        {
            const slong b = bits[j];

            printf("%-10ld", b);
            ulong n;
            if (b == 232)
                n = UWORD(1) << 32;
            else if (b == 263)
                n = UWORD(1) << 63;
            else
                n = n_nextprime(UWORD(1) << (b-1), 0);
            for (slong i = 0; i < nlens; i++)
            {
                tfun(lens[i], n, state);
                printf(" ");
            }
            printf("\n");
        }
    }
    else if (argc == 3)  // function + nbits given
    {
        const timefun tfun = funs[atoi(argv[1])];
        const slong b = atoi(argv[2]);

        printf("\n%s\n", description[atoi(argv[1])]);
        printf("#bits\\len");
        for (slong i = 0; i < nlens; i++)
            printf("%9ld", lens[i]);
        printf("\n");

        printf("%-10ld", b);
        ulong n;
        if (b == 232)
            n = UWORD(1) << 32;
        else if (b == 263)
            n = UWORD(1) << 63;
        else
            n = n_nextprime(UWORD(1) << (b-1), 0);
        for (slong i = 0; i < nlens; i++)
        {
            tfun(lens[i], n, state);
            printf(" ");
        }
        printf("\n");
    }
    else if (argc == 4)  // function + nbits + len given
    {
        const timefun tfun = funs[atoi(argv[1])];
        const slong b = atoi(argv[2]);
        const slong len = atoi(argv[3]);

        printf("\n%s\n", description[atoi(argv[1])]);
        printf("#bits\\len");
        for (slong i = 0; i < nlens; i++)
            printf("%9ld", lens[i]);
        printf("\n");

        printf("%-10ld", b);
        ulong n;
        if (b == 232)
            n = UWORD(1) << 32;
        else if (b == 263)
            n = UWORD(1) << 63;
        else
            n = n_nextprime(UWORD(1) << (b-1), 0);

        tfun(len, n, state);
        printf("\n");
    }

    flint_rand_clear(state);
    return 0;
}
