#include <stdlib.h>  // for atoi

#include <flint/ulong_extras.h>
#include <flint/profiler.h>
#include <flint/nmod.h>
#include <flint/nmod_vec.h>

#include "nmod32_vec.h"

typedef struct {slong nrows; slong len; slong modn;} time_args;

// utility (nmod vec uniform random)
static inline
void _nmod32_vec_rand(n32_ptr vec, flint_rand_t state, slong len, nmod_t mod)
{
    for (slong i = 0; i < len; i++)
        vec[i] = n_randint(state, mod.n);
}

/*---------------*/
/* direct: dot   */
/*---------------*/

#define TIME_DOT(fun, arg) \
void time_##fun(time_args targs, flint_rand_t state) \
{ \
    const slong len = targs.len; \
    const slong n = targs.modn; \
    \
    nmod_t mod; \
    nmod_init(&mod, n); \
    ulong pow2_precomp; \
    NMOD_RED(pow2_precomp, (UWORD(1) << DOT_SPLIT_BITS), mod); \
    \
    n32_ptr v1 = _nmod32_vec_init(len); \
    _nmod32_vec_rand(v1, state, len, mod); \
    n32_ptr v2 = _nmod32_vec_init(len); \
    _nmod32_vec_rand(v2, state, len, mod); \
    \
    volatile uint FLINT_SET_BUT_UNUSED(res); \
    \
    double FLINT_SET_BUT_UNUSED(tcpu), twall; \
    \
    TIMEIT_START \
    res = _nmod32_vec_##fun(v1, v2, len, mod, arg); \
    TIMEIT_STOP_VALUES(tcpu, twall) \
    \
    printf("%.2e", twall); \
    \
    _nmod32_vec_clear(v1); \
    _nmod32_vec_clear(v2); \
}

TIME_DOT(dot_split, pow2_precomp);
TIME_DOT(dot_split_avx2, pow2_precomp);
TIME_DOT(dot_split_avx512, pow2_precomp);
TIME_DOT(dot_ifma_avx2, pow2_precomp);
TIME_DOT(dot_ifma_avx512, pow2_precomp);

void time_dot_msolve_avx2(time_args targs, flint_rand_t state)
{
    const slong len = targs.len;
    const slong n = targs.modn;

    nmod_t mod;
    nmod_init(&mod, n);

    n32_ptr v1 = _nmod32_vec_init(len);
    _nmod32_vec_rand(v1, state, len, mod);
    n32_ptr v2 = _nmod32_vec_init(len);
    _nmod32_vec_rand(v2, state, len, mod);

    volatile uint FLINT_SET_BUT_UNUSED(res);

    double FLINT_SET_BUT_UNUSED(tcpu), twall;

    TIMEIT_START
    res = _nmod32_vec_dot_msolve_avx2(v1, v2, len, mod.n);
    TIMEIT_STOP_VALUES(tcpu, twall)

    printf("%.2e", twall);

    _nmod32_vec_clear(v1);
    _nmod32_vec_clear(v2);
}

/*-----------------------*/
/* indirect: multi-dot   */
/*-----------------------*/

#define TIME_MDOT(fun)                                  \
void time_##fun(time_args targs, flint_rand_t state)    \
{                                                       \
    const slong nrows = targs.nrows;                    \
    const slong len = targs.len;                        \
    const slong n = targs.modn;                         \
                                                        \
    nmod_t mod;                                         \
    nmod_init(&mod, n);                                 \
                                                        \
    n32_ptr mat = _nmod32_vec_init(nrows*len);          \
    _nmod32_vec_rand(mat, state, nrows*len, mod);       \
    n32_ptr vec = _nmod32_vec_init(len);                \
    _nmod32_vec_rand(vec, state, len, mod);             \
                                                        \
    n32_ptr res = _nmod32_vec_init(nrows);              \
                                                        \
    double FLINT_SET_BUT_UNUSED(tcpu), twall;           \
                                                        \
    TIMEIT_START                                        \
    _nmod32_vec_##fun(res, mat, vec,                    \
                      nrows, len, len, mod);            \
    TIMEIT_STOP_VALUES(tcpu, twall)                     \
                                                        \
    printf("%.2e", twall);                              \
                                                        \
    _nmod32_vec_clear(mat);                             \
    _nmod32_vec_clear(vec);                             \
    _nmod32_vec_clear(res);                             \
}

TIME_MDOT(mdot_split);
TIME_MDOT(mdot_split_avx2);
TIME_MDOT(mdot_split_avx512);
TIME_MDOT(mdot_msolve_native_avx2);
TIME_MDOT(mdot_msolve_via_dot_avx2);

TIME_MDOT(mdot2_split);
TIME_MDOT(mdot2_split_avx2);
TIME_MDOT(mdot2_split_avx512);

TIME_MDOT(mdot3_split_avx2);
TIME_MDOT(mdot4_split_avx512);

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
    const slong nfuns = 15;
    typedef void (*timefun) (time_args, flint_rand_t);
    const timefun funs[] = {
        time_dot_split,                 // 0
        time_dot_split_avx2,            // 1
        time_dot_split_avx512,          // 2
        time_dot_msolve_avx2,           // 3
        time_dot_ifma_avx2,             // 4
        time_dot_ifma_avx512,           // 5
        time_mdot_split,                // 6
        time_mdot_split_avx2,           // 7
        time_mdot_split_avx512,         // 8
        time_mdot_msolve_via_dot_avx2,  // 9
        time_mdot_msolve_native_avx2,   // 10
        time_mdot2_split,               // 11
        time_mdot2_split_avx2,          // 12
        time_mdot2_split_avx512,        // 13
        time_mdot3_split_avx2,          // 14
        time_mdot4_split_avx512,        // 15
    };

    const char * description[] = {
        "#0  --> dot_split                     ",
        "#1  --> dot_split_avx2                ",
        "#2  --> dot_split_avx512              ",
        "#3  --> dot_msolve_avx2               ",
        "#4  --> dot_ifma_avx2                 ",
        "#5  --> dot_ifma_avx512               ",
        "#6  --> mdot_split                    ",
        "#7  --> mdot_split_avx2               ",
        "#8  --> mdot_split_avx512             ",
        "#9  --> mdot_msolve_via_dot_avx2      ",
        "#10 --> mdot_msolve_native_avx2       ",
        "#11 --> mdot2_split                   ",
        "#12 --> mdot2_split_avx2              ",
        "#13 --> mdot2_split_avx512            ",
        "#14 --> mdot3_split_avx2              ",
        "#15 --> mdot4_split_avx512            ",
    };

    if (argc == 1)  // show usage
    {
        printf("Usage: `%s [nbits] [len] [fun]`\n", argv[0]);
        printf("   Each argument is optional; no argument shows this help.\n");
        printf("   - nbits: number of bits (in (0..31]) for the modulus, chosen as nextprime(2**(nbits-1))\n");
        printf("        (nbits == -1 launches full suite)\n");
        printf("   - len: length of the vectors\n");
        printf("   - fun: id number of the timed function (see below),\n");
        printf("\nAvailable functions:\n");
        for (slong j = 0; j < nfuns; j++)
            printf("   %s\n", description[j]);

        return 0;
    }
    else
    {
        printf("\nAvailable functions:\n");
        for (slong j = 0; j < nfuns; j++)
            printf("   %s\n", description[j]);
        printf("\nnrows types (indicated in column nrows):\n");
        printf("#0  --> 1 row (== standard dot product) \n");
        printf("#1  --> 4 rows                          \n");
        printf("#2  --> max(1, len/10) rows             \n");
        printf("#3  --> max(1, len/2) rows              \n");
        printf("#4  --> len rows                        \n");
        printf("\n");
    }

    printf("#warmup...\n");
    for (slong i = 0; i < 3; i++)
    {
        time_args targs = {1, 10000, UWORD(1) << 20};
        time_dot_split_avx2(targs, state);
        printf(" ");
    }
    printf("\n\n");

    if (argc == 2 && atoi(argv[1]) == -1)  // launching full suite
    {
        printf("           len");
        for (slong i = 0; i < nlens; i++)
            printf("%9ld", lens[i]);
        printf("\n");
        printf("bits nrows fun\n");
        for (slong j = 0; j < nbits; j++)
        {
            const slong b = bits[j];
            ulong n;
            n = n_nextprime(UWORD(1) << (b-1), 0);
            for (slong nrows = 2; nrows < 5; nrows++)
            {
                for (slong ifun = 0; ifun < nfuns; ifun++)
                {
                    const timefun tfun = funs[ifun];
                    printf("%-5ld#%-5ld#%-3ld", b, nrows, ifun);
                    for (slong i = 0; i < nlens; i++)
                    {
                        time_args targs = {1, lens[i], n};
                        if (nrows == 1)
                            targs.nrows = 4;
                        if (nrows == 2)
                            targs.nrows = FLINT_MAX(1, lens[i]/10);
                        if (nrows == 3)
                            targs.nrows = FLINT_MAX(1, lens[i]/2);
                        if (nrows == 4)
                            targs.nrows = lens[i];

                        if (nrows == 0 ||
                            (ifun >= 4 && lens[i] * targs.nrows < 10000000000))
                            tfun(targs, state);
                        printf(" ");
                    }
                    printf("\n");
                }
            }
        }
    }
    else if (argc == 2)  // nbits is given
    {
        printf("       len");
        for (slong i = 0; i < nlens; i++)
            printf("%9ld", lens[i]);
        printf("\n");
        printf("bits fun\n");
        const slong b = atoi(argv[1]);
        ulong n;
        n = n_nextprime(UWORD(1) << (b-1), 0);
        for (slong ifun = 0; ifun < nfuns; ifun++)
        {
            const timefun tfun = funs[ifun];
            printf("%-5ld#%-4ld", b, ifun);
            for (slong i = 0; i < nlens; i++)
            {
                time_args targs = {1, lens[i], n};
                tfun(targs, state);
                printf(" ");
            }
            printf("\n");
        }
    }
    else if (argc == 3)  // nbits + len given
    {
        const slong len = atoi(argv[2]);
        printf("       len");
        printf("%9ld", len);
        printf("\n");
        printf("bits fun\n");
        const slong b = atoi(argv[1]);
        ulong n;
        n = n_nextprime(UWORD(1) << (b-1), 0);
        for (slong ifun = 0; ifun < nfuns; ifun++)
        {
            const timefun tfun = funs[ifun];
            printf("%-5ld#%-4ld", b, ifun);
            time_args targs = {1, len, n};
            tfun(targs, state);
            printf(" ");
            printf("\n");
        }
    }
    else if (argc == 4)  // function + nbits + len given
    {
        const slong b = atoi(argv[1]);
        const slong len = atoi(argv[2]);
        const slong ifun = atoi(argv[3]);
        const timefun tfun = funs[ifun];
        printf("       len");
        printf("%9ld", len);
        printf("\n");
        printf("bits nrows fun\n");
        ulong n;
        n = n_nextprime(UWORD(1) << (b-1), 0);

        for (slong nrows = 0; nrows < 5; nrows++)
        {
            printf("%-5ld#%-5ld#%-3ld", b, nrows, ifun);
            time_args targs = {1, len, n};
            if (nrows == 1)
                targs.nrows = 4;
            if (nrows == 2)
                targs.nrows = FLINT_MAX(1, len/10);
            if (nrows == 3)
                targs.nrows = FLINT_MAX(1, len/2);
            if (nrows == 4)
                targs.nrows = len;
            if (nrows == 0 ||
                (ifun >= 4 && len * targs.nrows < WORD(10000000000)))
                tfun(targs, state);
            printf(" ");
            printf("\n");
        }
    }

    flint_rand_clear(state);
    return 0;
}
