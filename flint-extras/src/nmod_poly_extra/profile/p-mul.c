#include <flint/flint.h>
#include <stdlib.h>  // for atoi

#include <flint/ulong_extras.h>
#include <flint/profiler.h>
#include <flint/nmod.h>
#include <flint/nmod_poly.h>
#include "nmod_poly_extra.h"
#include "nmod_vec_extra.h"


#define MEASURE_SAMPLE 0

/* input polynomials of length len1 x len2: multiply, or middle product with parameters nlo, nhi */
typedef struct {
    slong len1;
    slong len2;
    slong nlo;
    slong nhi;
    slong modn;
} time_args;

#define TIME_MUL(fun)                                   \
void time_nmod_poly_##fun(time_args targs,              \
                          flint_rand_t state)           \
{                                                       \
    const slong len1 = targs.len1;                      \
    const slong len2 = targs.len2;                      \
    const slong modn = targs.modn;                      \
    nmod_t mod;                                         \
    nmod_init(&mod, modn);                              \
                                                        \
    nn_ptr pol1 = FLINT_ARRAY_ALLOC(len1, ulong);       \
    nn_ptr pol2 = FLINT_ARRAY_ALLOC(len2, ulong);       \
    nn_ptr res = FLINT_ARRAY_ALLOC(len1+len2-1, ulong); \
    _nmod_vec_rand(pol1, state, len1, mod);             \
    _nmod_vec_rand(pol2, state, len2, mod);             \
                                                        \
    double FLINT_SET_BUT_UNUSED(tcpu), twall;           \
                                                        \
    TIMEIT_START;                                       \
    _nmod_poly_##fun(res, pol1, len1, pol2, len2, mod); \
    TIMEIT_STOP_VALUES(tcpu, twall);                    \
                                                        \
    flint_printf("%.2e", twall);                        \
                                                        \
    flint_free(pol1);                                   \
    flint_free(pol2);                                   \
    flint_free(res);                                    \
}

#define TIME_MULMID(fun)                                \
void time_nmod_poly_##fun(time_args targs,              \
                          flint_rand_t state)           \
{                                                       \
    const slong len1 = targs.len1;                      \
    const slong len2 = targs.len2;                      \
    const slong nlo = targs.nlo;                        \
    const slong nhi = targs.nhi;                        \
    const slong modn = targs.modn;                      \
    nmod_t mod;                                         \
    nmod_init(&mod, modn);                              \
                                                        \
    nn_ptr pol1 = FLINT_ARRAY_ALLOC(len1, ulong);       \
    nn_ptr pol2 = FLINT_ARRAY_ALLOC(len2, ulong);       \
    nn_ptr res = FLINT_ARRAY_ALLOC(nhi - nlo, ulong);   \
    _nmod_vec_rand(pol1, state, len1, mod);             \
    _nmod_vec_rand(pol2, state, len2, mod);             \
                                                        \
    double FLINT_SET_BUT_UNUSED(tcpu), twall;           \
                                                        \
    TIMEIT_START;                                       \
    _nmod_poly_##fun(res, pol1, len1, pol2, len2,       \
                      nlo, nhi, mod);                   \
    TIMEIT_STOP_VALUES(tcpu, twall);                    \
                                                        \
    flint_printf("%.2e", twall);                        \
                                                        \
    flint_free(pol1);                                   \
    flint_free(pol2);                                   \
    flint_free(res);                                    \
}

#define SAMPLE_MUL(fun)                                         \
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
    nmod_poly_t prod;                                           \
    nmod_poly_init(pol1, modn);                                 \
    nmod_poly_init(pol2, modn);                                 \
    nmod_poly_init(prod, modn);                                 \
    nmod_poly_rand(pol1, state, len1);                          \
    nmod_poly_rand(pol2, state, len2);                          \
	for (ulong i = 0; i < count; i++)                           \
    {                                                           \
        nmod_poly_t tmp1;                                       \
        nmod_poly_t tmp2;                                       \
        nmod_poly_init(tmp1, modn);                             \
        nmod_poly_init(tmp2, modn);                             \
        nmod_poly_set(tmp1, pol1);                              \
        nmod_poly_set(tmp2, pol2);                              \
        prof_start();                                           \
        fun(prod, tmp1, tmp2);                                  \
        prof_stop();                                            \
        nmod_poly_clear(tmp1);                                  \
        nmod_poly_clear(tmp2);                                  \
    }                                                           \
                                                                \
    nmod_poly_clear(pol1);                                      \
    nmod_poly_clear(pol2);                                      \
    nmod_poly_clear(prod);                                      \
    FLINT_TEST_CLEAR(state);                                    \
}

TIME_MUL(mul)
TIME_MUL(mul_classical)
TIME_MUL(mul_KS)
TIME_MUL(mul_KS2)
TIME_MUL(mul_KS4)

TIME_MULMID(mulmid)
TIME_MULMID(mulmid_classical)
TIME_MULMID(mulmid_KS)
TIME_MULMID(mulmid_fft_small)

SAMPLE_MUL(nmod_poly_mul)
SAMPLE_MUL(nmod_poly_mul_classical)
SAMPLE_MUL(nmod_poly_mul_KS)
SAMPLE_MUL(nmod_poly_mul_KS2)
SAMPLE_MUL(nmod_poly_mul_KS4)

/*-------------------------*/
/*  main                   */
/*-------------------------*/
int main(int argc, char ** argv)
{
    flint_rand_t state;
    flint_rand_init(state);
    flint_rand_set_seed(state, time(NULL), time(NULL)+129384125L);

    // bench functions
    const slong nfuns = 9;
    typedef void (*timefun) (time_args, flint_rand_t);
    const timefun funs[] = {
        time_nmod_poly_mul,               // 0
        time_nmod_poly_mul_classical,     // 1
        time_nmod_poly_mul_KS,            // 2
        time_nmod_poly_mul_KS2,           // 3
        time_nmod_poly_mul_KS4,           // 4
        time_nmod_poly_mulmid,            // 5
        time_nmod_poly_mulmid_classical,  // 6
        time_nmod_poly_mulmid_KS,         // 7
        time_nmod_poly_mulmid_fft_small,  // 8
    };

    typedef void (*samplefun) (void*, ulong);
    const samplefun sfuns[] = {
        sample_nmod_poly_mul,               // 0
        sample_nmod_poly_mul_classical,     // 1
        sample_nmod_poly_mul_KS,            // 2
        sample_nmod_poly_mul_KS2,           // 3
        sample_nmod_poly_mul_KS4,           // 4
    };

    const char * description[] = {
        "#0  --> nmod_poly_mul              ", 
        "#1  --> nmod_poly_mul_classical    ", 
        "#2  --> nmod_poly_mul_KS           ", 
        "#3  --> nmod_poly_mul_KS2          ", 
        "#4  --> nmod_poly_mul_KS4          ", 
        "#5  --> nmod_poly_mulmid           ", 
        "#6  --> nmod_poly_mulmid_classical ", 
        "#7  --> nmod_poly_mulmid_KS        ", 
        "#8  --> nmod_poly_mulmid_fft_small ", 
    };

    if (argc == 1)  // show usage
    {
        flint_printf("Usage: `%s [nbits] [fun] [len1] [len2] [nlo] [nhi]`\n", argv[0]);
        flint_printf("   No argument shows this help.\n");
        flint_printf("   - nbits: number of bits, in (0..63], for the randomly-chosen modulus\n");
        flint_printf("   - fun: id number of the timed function (see below),\n");
        flint_printf("   - len1, len2: length of polynomials to be multiplied\n");
        flint_printf("   - nlo, nhi: parameters for mulmid (can be ignored for 0 <= fun < 5)\n");
        flint_printf("     (requires 0 <= nlo < nhi <= len1 + len2 - 1)\n");
        flint_printf("\nAvailable functions:\n");
        for (slong j = 0; j < nfuns; j++)
            flint_printf("   %s\n", description[j]);

        return 0;
    }

    /* flint_printf("#warmup...\n"); */
    for (slong i = 0; i < 3; i++)
    {
        time_args targs = {
            i*100,
            i*100,
            i*100,
            2*i*100 - 1,
            n_nextprime(17 + (UWORD(1) << (i*17)), 0)
        };
        double min, max;
        prof_repeat(&min, &max, sfuns[0], &targs);
    }
    /* flint_printf("\n\n"); */

    if (argc == 5 && atoi(argv[2]) < 5)
    {
        const slong b = atoi(argv[1]);
        const slong ifun = atoi(argv[2]);
        const slong len1 = atoi(argv[3]);
        const slong len2 = atoi(argv[4]);
        flint_printf("bits fun len1      len2      \n");

        ulong n = n_nextprime(UWORD(1) << (b-1), 1);

        flint_printf("%-4ld %-3ld %-10ld%-10ld", b, ifun, len1, len2);
        time_args targs = {len1, len2, 0, 0, n};

#if MEASURE_SAMPLE
        double min, max;
        prof_repeat(&min, &max, sfuns[ifun], &targs);
        flint_printf("%.2e", min/1000000);
#else
        funs[ifun](targs, state);
#endif
        flint_printf("\n");
    }

    if (argc == 7)
    {
        const slong b = atoi(argv[1]);
        const slong ifun = atoi(argv[2]);
        const slong len1 = atoi(argv[3]);
        const slong len2 = atoi(argv[4]);
        const slong nlo = atoi(argv[5]);
        const slong nhi = atoi(argv[6]);
        flint_printf("bits fun len1      len2      nlo       nhi       \n");

        ulong n = n_nextprime(UWORD(1) << (b-1), 1);

        flint_printf("%-4ld %-3ld %-10ld%-10ld%-10ld%-10ld", b, ifun, len1, len2, nlo, nhi);
        time_args targs = {len1, len2, nlo, nhi, n};

#if MEASURE_SAMPLE
        double min, max;
        prof_repeat(&min, &max, sfuns[ifun], &targs);
        flint_printf("%.2e", min/1000000);
#else
        funs[ifun](targs, state);
#endif
        flint_printf("\n");
    }

    /* Special: bench all mul */
    if (argc == 4)
    {
        const slong b    = atoi(argv[1]);
        const slong len1 = atoi(argv[2]);
        const slong len2 = atoi(argv[3]);
        ulong n = n_nextprime(UWORD(1) << (b-1), 1);

        flint_printf("%-4ld %-5ld%-5ld", b, len1, len2);
        for (slong ifun = 0; ifun < 5; ifun++)
        {
            time_args targs = {len1, len2, 0, 0, n};
            funs[ifun](targs, state);
            flint_printf(" ");
        }
        flint_printf("\n");
    }

    /* Special: bench all mulmid */
    if (argc == 6)
    {
        const slong b    = atoi(argv[1]);
        const slong len1 = atoi(argv[2]);
        const slong len2 = atoi(argv[3]);
        const slong nlo  = atoi(argv[4]);
        const slong nhi  = atoi(argv[5]);
        ulong n = n_nextprime(UWORD(1) << (b-1), 1);

        flint_printf("%-4ld %-5ld%-5ld%-5ld%-5ld", b, len1, len2, nlo, nhi);
        for (slong ifun = 5; ifun < 9; ifun++)
        {
            time_args targs = {len1, len2, nlo, nhi, n};
            funs[ifun](targs, state);
            flint_printf(" ");
        }
        flint_printf("\n");
    }


    flint_rand_clear(state);
    return 0;
}
