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

#define TIME_MDOT(fun) \
void time_##fun(time_args targs, flint_rand_t state) \
{ \
    const slong nrows = targs.nrows; \
    const slong len = targs.len; \
    const slong n = targs.modn; \
 \
    nmod_t mod; \
    nmod_init(&mod, n); \
 \
    n32_ptr mat = _nmod32_vec_init(nrows*len); \
    _nmod32_vec_rand(mat, state, nrows*len, mod); \
    n32_ptr vec = _nmod32_vec_init(len); \
    _nmod32_vec_rand(vec, state, len, mod); \
 \
    n32_ptr res = _nmod32_vec_init(nrows); \
 \
    double FLINT_SET_BUT_UNUSED(tcpu), twall; \
 \
    TIMEIT_START \
    _nmod32_vec_##fun(res, mat, vec, nrows, len, len, mod); \
    TIMEIT_STOP_VALUES(tcpu, twall) \
 \
    printf("%.2e", twall); \
 \
    _nmod32_vec_clear(mat); \
    _nmod32_vec_clear(vec); \
    _nmod32_vec_clear(res); \
}

TIME_MDOT(mdot_split);
TIME_MDOT(mdot_split_avx2);
TIME_MDOT(mdot_split_avx512);
TIME_MDOT(mdot_msolve_native_avx2);
TIME_MDOT(mdot_msolve_via_dot_avx2);

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
    const slong nfuns = 9;
    typedef void (*timefun) (time_args, flint_rand_t);
    const timefun funs[] = {
        time_dot_split,                 // 0
        time_dot_split_avx2,            // 1
        time_dot_split_avx512,          // 2
        time_dot_msolve_avx2,           // 3
        time_mdot_split,                // 4
        time_mdot_split_avx2,           // 5
        time_mdot_split_avx512,         // 6
        time_mdot_msolve_via_dot_avx2,  // 7
        time_mdot_msolve_native_avx2,   // 8
    };

    const char * description[] = {
        "#0  --> dot_split                     ",
        "#1  --> dot_split_avx2                ",
        "#2  --> dot_split_avx512              ",
        "#3  --> dot_msolve_avx2               ",
        "#4  --> mdot_split                    ",
        "#5  --> mdot_split_avx2               ",
        "#6  --> mdot_split_avx2               ",
        "#7  --> mdot_msolve_via_dot_avx2      ",
        "#8  --> mdot_msolve_native_avx2       ",
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
        time_args targs = {1, 10000, UWORD(1) << 20};
        time_dot_split_avx2(targs, state);
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
                n = n_nextprime(UWORD(1) << (b-1), 0);
                for (slong i = 0; i < nlens; i++)
                {
                    time_args targs = {1, lens[i], n};
                    tfun(targs, state);
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
                time_args targs = {1, lens[i], n};
                tfun(targs, state);
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
            time_args targs = {1, lens[i], n};
            tfun(targs, state);
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

        time_args targs = {1, len, n};
        tfun(targs, state);
        printf("\n");
    }

    flint_rand_clear(state);
    return 0;
}
