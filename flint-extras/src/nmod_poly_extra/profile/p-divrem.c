#include <flint/flint.h>
#include <stdlib.h>  // for atoi

#include <flint/ulong_extras.h>
#include <flint/profiler.h>
#include <flint/nmod.h>
#include <flint/nmod_poly.h>
#include "nmod_poly_extra.h"

#define MEASURE_SAMPLE 0

/* for divrem for lengths len1 and len2 */
typedef struct {
    slong len1;
    slong len2;
    slong modn;
} time_args;

#define TIME_DIV(fun) \
void time_##fun(time_args targs, flint_rand_t state) \
{ \
    const slong len1 = targs.len1; \
    const slong len2 = targs.len2; \
    const slong modn = targs.modn; \
    \
    nmod_poly_t pol1; \
    nmod_poly_t pol2; \
    nmod_poly_init(pol1, modn); \
    nmod_poly_init(pol2, modn); \
    nmod_poly_rand(pol1, state, len1); \
    nmod_poly_rand(pol2, state, len2); \
    nmod_poly_t quo; \
    nmod_poly_t rem; \
    nmod_poly_init(quo, modn); \
    nmod_poly_init(rem, modn); \
    \
    double FLINT_SET_BUT_UNUSED(tcpu), twall; \
    \
    TIMEIT_START; \
    fun(quo, rem, pol1, pol2); \
    TIMEIT_STOP_VALUES(tcpu, twall); \
    \
    flint_printf("%.2e", twall); \
    \
    nmod_poly_clear(pol1); \
    nmod_poly_clear(pol2); \
    nmod_poly_clear(quo); \
    nmod_poly_clear(rem); \
}

#define SAMPLE_DIV(fun)                                         \
void sample_##fun(void * arg, ulong count)                      \
{                                                               \
    time_args * targs = (time_args *) arg;                      \
    const slong len1 = targs->len1;                             \
    const slong len2 = targs->len2;                             \
    const slong modn = targs->modn;                             \
                                                                \
    FLINT_TEST_INIT(state);                                     \
                                                                \
    nmod_poly_t pol1;                                           \
    nmod_poly_t pol2;                                           \
    nmod_poly_init(pol1, modn);                                 \
    nmod_poly_init(pol2, modn);                                 \
    nmod_poly_rand(pol1, state, len1);                          \
    nmod_poly_rand(pol2, state, len2);                          \
    nmod_poly_t quo;                                            \
    nmod_poly_t rem;                                            \
    nmod_poly_init(quo, modn);                                  \
    nmod_poly_init(rem, modn);                                  \
	for (ulong i = 0; i < count; i++)                           \
    {                                                           \
        nmod_poly_t tmp1;                                       \
        nmod_poly_t tmp2;                                       \
        nmod_poly_init(tmp1, modn);                             \
        nmod_poly_init(tmp2, modn);                             \
        nmod_poly_set(tmp1, pol1);                              \
        nmod_poly_set(tmp2, pol2);                              \
        prof_start();                                           \
        fun(quo, rem, tmp1, tmp2);                              \
        prof_stop();                                            \
        nmod_poly_clear(tmp1);                                  \
        nmod_poly_clear(tmp2);                                  \
    }                                                           \
                                                                \
    nmod_poly_clear(pol1);                                      \
    nmod_poly_clear(pol2);                                      \
    nmod_poly_clear(quo);                                      \
    nmod_poly_clear(rem);                                      \
    FLINT_TEST_CLEAR(state);                                    \
}

TIME_DIV(nmod_poly_divrem)
TIME_DIV(nmod_poly_divrem_basecase)

SAMPLE_DIV(nmod_poly_divrem)
SAMPLE_DIV(nmod_poly_divrem_basecase)

/*-------------------------*/
/*  main                   */
/*-------------------------*/
int main(int argc, char ** argv)
{
    flint_rand_t state;
    flint_rand_init(state);
    flint_rand_set_seed(state, time(NULL), time(NULL)+129384125L);

    // bench functions
    const slong nfuns = 2;
    typedef void (*timefun) (time_args, flint_rand_t);
    const timefun funs[] = {
        time_nmod_poly_divrem,               // 0
        time_nmod_poly_divrem_basecase,      // 1
    };

    typedef void (*samplefun) (void*, ulong);
    const samplefun sfuns[] = {
        sample_nmod_poly_divrem,               // 0
        sample_nmod_poly_divrem_basecase,      // 1
    };

    const char * description[] = {
        "#0  --> nmod_poly_divrem           ", 
        "#1  --> nmod_poly_divrem_basecase  ", 
    };

    if (argc == 1)  // show usage
    {
        flint_printf("Usage: `%s [nbits] [fun] [len1] [len2]`\n", argv[0]);
        flint_printf("   No argument shows this help.\n");
        flint_printf("   - nbits: number of bits, in (0..63], for the randomly-chosen modulus\n");
        flint_printf("   - fun: id number of the timed function (see below),\n");
        flint_printf("   - len1, len2: length of polynomials\n");
        flint_printf("\nAvailable functions:\n");
        for (slong j = 0; j < nfuns; j++)
            flint_printf("   %s\n", description[j]);

        return 0;
    }

    /* flint_printf("#warmup...\n"); */
    for (slong i = 0; i < 3; i++)
    {
        time_args targs = {
            (i+1)*100,
            (i+1)*100,
            n_nextprime(17 + (UWORD(1) << (i*17)), 0)
        };
        double min, max;
        prof_repeat(&min, &max, sfuns[0], &targs);
    }
    /* flint_printf("\n\n"); */

    if (argc == 5)
    {
        const slong b = atoi(argv[1]);
        const slong ifun = atoi(argv[2]);
        const slong len1 = atoi(argv[3]);
        const slong len2 = atoi(argv[4]);
        flint_printf("bits fun len1      len2      \n");

        ulong n = n_nextprime(UWORD(1) << (b-1), 1);

        flint_printf("%-4ld %-3ld %-10ld%-10ld", b, ifun, len1, len2);
        time_args targs = {len1, len2, n};

#if MEASURE_SAMPLE
        double min, max;
        prof_repeat(&min, &max, sfuns[ifun], &targs);
        flint_printf("%.2e", min/1000000);
#else
        funs[ifun](targs, state);
#endif
        flint_printf("\n");
    }

    flint_rand_clear(state);
    return 0;
}
