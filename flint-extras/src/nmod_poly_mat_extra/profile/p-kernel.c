#include <stdlib.h>  // for atoi

#include <flint/ulong_extras.h>
#include <flint/profiler.h>
#include <flint/nmod.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_poly_mat.h>

#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_kernel.h"

#define MEASURE_SAMPLE 0

/* Shift types are:
 * -1 -> [TODO not implemented yet] negative shift that represents known degree bound on a kernel basis
 * 0 -> uniform, shift == (0,...,0)
 * 1 -> input degrees: shift == row degrees of input matrix
 * 100 -> shift that yields the kernel basis Hermite form
 * */

typedef struct
{
    slong rdim;  /* row dimension */
    slong cdim;  /* column dimension */
    slong deg;   /* degree */
    slong rank;  /* rank */
    slong stype; /* shift type */
    slong modn;  /* modulus */
}
time_args;

#define TIME_KER_ZLS(fun)                               \
void time_##fun(time_args targs, flint_rand_t state)    \
{                                                       \
    const slong rdim = targs.rdim;                      \
    const slong cdim = targs.cdim;                      \
    const slong deg = targs.deg;                        \
    /* const slong rank = targs.rank; */ /* TODO */           \
    /* const slong stype = targs.stype; */ /* TODO */         \
    const slong n = targs.modn;                         \
                                                        \
    nmod_t mod;                                         \
    nmod_init(&mod, n);                                 \
                                                        \
    nmod_poly_mat_t F;                                  \
    nmod_poly_mat_init(F, rdim, cdim, n);               \
    nmod_poly_mat_rand(F, state, deg);                  \
                                                        \
    /* TMP */                                           \
    nmod_poly_mat_t Ft;                                 \
    nmod_poly_mat_init(Ft, cdim, rdim, n);              \
    nmod_poly_mat_transpose(Ft, F);                     \
    /* END TMP */                                       \
                                                        \
    slong degN[rdim];                                   \
    nmod_poly_mat_t K;                                  \
    nmod_poly_mat_init(K, rdim, rdim, n);               \
                                                        \
    double FLINT_SET_BUT_UNUSED(tcpu), twall;           \
                                                        \
    TIMEIT_START;                                       \
    nmod_poly_mat_##fun(K, degN, Ft, NULL, 3.);         \
    TIMEIT_STOP_VALUES(tcpu, twall);                    \
                                                        \
    printf("%.2e", twall);                              \
                                                        \
    nmod_poly_mat_clear(F);                             \
    nmod_poly_mat_clear(Ft); /* TMP */                  \
    nmod_poly_mat_clear(K);                             \
}

#define TIME_KER(fun)                                   \
void time_##fun(time_args targs, flint_rand_t state)    \
{                                                       \
    const slong rdim = targs.rdim;                      \
    const slong cdim = targs.cdim;                      \
    const slong deg = targs.deg;                        \
    /* const slong rank = targs.rank; */ /* TODO */           \
    /* const slong stype = targs.stype; */ /* TODO */         \
    const slong n = targs.modn;                         \
                                                        \
    nmod_t mod;                                         \
    nmod_init(&mod, n);                                 \
                                                        \
    nmod_poly_mat_t pmat;                               \
    nmod_poly_mat_init(pmat, rdim, cdim, n);            \
    nmod_poly_mat_rand(pmat, state, deg);               \
                                                        \
    slong pivind[rdim];                                 \
    slong shift[rdim];                                  \
    for (slong i = 0; i < rdim; i++)                    \
        shift[i] = 0;                                   \
    nmod_poly_mat_t ker;                                \
    nmod_poly_mat_init(ker, rdim, rdim, n);             \
                                                        \
    double FLINT_SET_BUT_UNUSED(tcpu), twall;           \
                                                        \
    TIMEIT_START;                                       \
    nmod_poly_mat_##fun(ker, pivind, shift, pmat);      \
    TIMEIT_STOP_VALUES(tcpu, twall);                    \
                                                        \
    printf("%.2e", twall);                              \
                                                        \
    nmod_poly_mat_clear(pmat);                          \
    nmod_poly_mat_clear(ker);                           \
}

TIME_KER_ZLS(kernel_zls)
TIME_KER(kernel_via_approx)

/*-------------------------*/
/*  main                   */
/*-------------------------*/

int main(int argc, char ** argv)
{
    flint_rand_t state;
    flint_rand_init(state);
    flint_rand_set_seed(state, time(NULL), time(NULL)+129384125L);

    // modulus bitsize
    /* const slong nbits = 7; */
    /* const ulong bits[] = {12, 24, 30, 40, 50, 60, 63}; */

    // matrix dimensions
    /* const slong ndims = 10; */
    /* const ulong rdims[] = {2, 4, 6, 8, 11, 15, 20, 30, 50, 100}; */
    /* const ulong cdims[] = {1, 3, 5, 7,  9, 13, 17, 25, 40, 75}; */

    // matrix degrees
    const slong ndegs = 13;
    const ulong degs[] = {2, 5, 10, 20, 40, 80, 160, 320, 640, 1280, 2560, 5120, 10240};

    /* TODO rank */
    /* TODO shift type */

    // bench functions
    const slong nfuns = 2;
    typedef void (*timefun) (time_args, flint_rand_t);
    const timefun funs[] = {
        time_kernel_via_approx,               // 0
        time_kernel_zls,                      // 1
    };

    // TODO
    //typedef void (*samplefun) (void*, ulong);
    //const samplefun sfuns[] = {
    //    sample_mul,                      // 0
    //    sample_mul_geometric,            // 1
    //};
    //#if MEASURE_SAMPLE
    //                        const samplefun sfun = sfuns[ifun];
    //                        double min, max;
    //                        prof_repeat(&min, &max, sfun, (void*) &targs);
    //                        printf("%.2e", min/1000000);
    //#else
    //                        const timefun tfun = funs[ifun];
    //                        tfun(targs, state);
    //#endif

    const char * description[] = {
        "#0  --> kernel (general interface)   ",
        "#1  --> .... TODO                ",
    };

    if (argc < 4 || argc > 6)  // show usage
    {
        printf("Usage: `%s [nbits] [rdim] [cdim] [deg] [rank] [shift] [fun]`\n", argv[0]);
        printf("   No argument shows this help.\n");
        printf("   First 3 arguments are mandatory.\n");
        printf("   [rank] and [shift] not supported yet.\n");
        printf("   - nbits: number of bits in [2..64] for the modulus, chosen as nextprime(2**(nbits-1))\n");
        printf("   - rdim, cdim: input matrix is rdim x cdim\n");
        printf("   - deg: matrix is random of degree < deg (default: predefined list)\n");
        printf("   - rank: [unsupported] matrix is random of this rank\n");
        printf("   - shift: [unsupported] type of shift\n");
        printf("   - fun: id number of the timed function (see below),\n");
        printf("\nAvailable functions:\n");
        for (slong j = 0; j < nfuns; j++)
            printf("   %s\n", description[j]);

        return 0;
    }

    printf("#warmup...\n");
    for (slong i = 0; i < 3; i++)
    {
        /* rdim; cdim; deg; rank; stype; modn; */
        time_args targs = {8, 4, 1000, 4, 0, n_nextprime(UWORD(1) << 20, 0)};
        time_kernel_via_approx(targs, state);
        printf(" ");
    }
    printf("\n\n");

    if (argc == 4)  // nbits + rdim + cdim given
    {
        const slong rdim = atoi(argv[2]);
        const slong cdim = atoi(argv[3]);
        const slong b = atoi(argv[1]);
        printf("bits fun rdim cdim deg\n");
        ulong n = n_nextprime(UWORD(1) << (b-1), 0);
        for (slong ifun = 0; ifun < nfuns; ifun++)
        {
            const timefun tfun = funs[ifun];
            for (slong d = 0; d < ndegs; d++)
            {
                printf("%-5ld#%-3ld%-5ld%-5ld%-8ld", b, ifun, rdim, cdim, degs[d]);
                time_args targs = {rdim, cdim, degs[d], FLINT_MIN(rdim, cdim), 0, n};
                tfun(targs, state);
                printf(" ");
                printf("\n");
            }
        }
    }
    else if (argc == 5)  // nbits + rdim + cdim + deg given
    {
        const slong rdim = atoi(argv[2]);
        const slong cdim = atoi(argv[3]);
        const slong deg  = atoi(argv[4]);
        const slong b = atoi(argv[1]);
        printf("bits fun rdim cdim deg\n");
        ulong n = n_nextprime(UWORD(1) << (b-1), 0);
        for (slong ifun = 0; ifun < nfuns; ifun++)
        {
            const timefun tfun = funs[ifun];
            printf("%-5ld#%-3ld%-5ld%-5ld%-8ld", b, ifun, rdim, cdim, deg);
            /* rdim; cdim; deg; rank; stype; modn; */
            time_args targs = {rdim, cdim, deg, FLINT_MIN(rdim, cdim), 0, n};
            tfun(targs, state);
            printf(" ");
            printf("\n");
        }
    }

    flint_rand_clear(state);
    return 0;
}
