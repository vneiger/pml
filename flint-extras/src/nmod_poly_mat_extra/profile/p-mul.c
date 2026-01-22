#include <stdlib.h>  // for atoi

#include <flint/ulong_extras.h>
#include <flint/profiler.h>
#include <flint/nmod.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_poly_mat.h>

#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_multiply.h"

#define MEASURE_SAMPLE 0

typedef struct
{
    slong rdim;  /* row outer dimension */
    slong idim;  /* inner dimension */
    slong cdim;  /* column outer dimension */
    slong degl;   /* degree of left operand */
    slong degr;   /* degree of right operand */
    slong modn;  /* modulus */
}
time_args;

#define TIME_MUL(fun)                                   \
void time_##fun(time_args targs, flint_rand_t state)    \
{                                                       \
    const slong rdim = targs.rdim;                      \
    const slong idim = targs.idim;                      \
    const slong cdim = targs.cdim;                      \
    const slong degl = targs.degl;                      \
    const slong degr = targs.degr;                      \
    const slong n = targs.modn;                         \
                                                        \
    nmod_t mod;                                         \
    nmod_init(&mod, n);                                 \
                                                        \
    nmod_poly_mat_t A;                                  \
    nmod_poly_mat_init(A, rdim, idim, n);               \
    nmod_poly_mat_rand(A, state, degl);                  \
    nmod_poly_mat_t B;                                  \
    nmod_poly_mat_init(B, idim, cdim, n);               \
    nmod_poly_mat_rand(B, state, degr);                  \
    nmod_poly_mat_t C;                                  \
    nmod_poly_mat_init(C, rdim, cdim, n);               \
                                                        \
    double FLINT_SET_BUT_UNUSED(tcpu), twall;           \
                                                        \
    TIMEIT_START;                                       \
    nmod_poly_mat_##fun(C, A, B);                       \
    TIMEIT_STOP_VALUES(tcpu, twall);                    \
                                                        \
    flint_printf("%.2e", twall);                        \
                                                        \
    nmod_poly_mat_clear(A);                             \
    nmod_poly_mat_clear(B);                             \
    nmod_poly_mat_clear(C);                             \
}

TIME_MUL(mul)
TIME_MUL(mul_geometric)
TIME_MUL(mul_waksman)

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

    // matrix dimensions (all square for the moment)
    /* const slong ndims = 10; */
    /* const ulong dims[] = {2, 4, 6, 8, 11, 15, 20, 30, 50, 100}; */

    // matrix degrees
    const slong ndegs = 12;
    const ulong degs[] = {5, 10, 20, 40, 80, 160, 320, 640, 1280, 2560, 5120, 10240};

    // bench functions
    const slong nfuns = 3;
    typedef void (*timefun) (time_args, flint_rand_t);
    const timefun funs[] = {
        time_mul,                      // 0
        time_mul_waksman,              // 1
        time_mul_geometric,            // 2
    };

    /* typedef void (*samplefun) (void*, ulong); */
    /* const samplefun sfuns[] = { */
    /*     sample_mul,                      // 0 */
    /*     sample_mul_waksman,              // 1 */
    /*     sample_mul_geometric,            // 2 */
    /* }; */

    const char * description[] = {
        "#0  --> mul                          ",
        "#1  --> mul_waksman                  ",
        "#2  --> mul_geometric                ",
    };

    if (argc == 1)  // show usage
    {
        flint_printf("Usage: `%s [nbits] [fun] [dim1] [dim2] [dim3] [degl] [degr]`\n", argv[0]);
        flint_printf("   No argument shows this help.\n");
        flint_printf("   6 arguments [nbits] [fun] [dim1] [dim2] [dim3] runs for several degrees.\n");
        /* flint_printf("   4 arguments [nbits] [fun] [dim] [deg] is the square dim x dim case.\n"); */
        /* flint_printf("   3 arguments [nbits] [fun] [dim] does several degrees in the square dim x dim case.\n"); */
        flint_printf("   - nbits: number of bits (in (1..64]) for the modulus, chosen as nextprime(2**(nbits-1))\n");
        flint_printf("   - fun: id number of the timed function (see below),\n");
        flint_printf("      (fun == -1 launches all available functions)\n");
        flint_printf("   - dim1, dim2, dim3: matrices are dim1 x dim2 and dim2 x dim3\n");
        flint_printf("   - deg: matrices are random of degree < deg\n");
        flint_printf("\nAvailable functions:\n");
        for (slong j = 0; j < nfuns; j++)
            flint_printf("   %s\n", description[j]);

        return 0;
    }

    flint_printf("#warmup...\n");
    for (slong i = 0; i < 3; i++)
    {
        time_args targs = {4, 4, 4, 1000, 1000, UWORD(1) << 20};
        time_mul(targs, state);
        flint_printf(" ");
    }
    flint_printf("\n\n");

    if (argc == 8)  // nbits + fun + dim1 + dim2 + dim3 + degl + degr
    {
        const slong b = atoi(argv[1]);
        slong ifun = atoi(argv[2]);
        const slong dim1 = atoi(argv[3]);
        const slong dim2 = atoi(argv[4]);
        const slong dim3 = atoi(argv[5]);
        const slong degl = atoi(argv[6]);
        const slong degr = atoi(argv[7]);
        flint_printf("bits fun dim1 dim2 dim3 degl    degr\n");

        ulong n;
        n = n_nextprime(UWORD(1) << (b-1), 0);

        if (ifun == -1)
        {
            for (ifun = 0; ifun < nfuns; ifun++)
            {
                const timefun tfun = funs[ifun];
                flint_printf("%-5ld#%-3ld%-5ld%-5ld%-5ld%-8ld%-8ld", b, ifun, dim1, dim2, dim3, degl, degr);
                time_args targs = {dim1, dim2, dim3, degl, degr, n};
                tfun(targs, state);
                flint_printf("\n");
            }
        }
        else
        {
            const timefun tfun = funs[ifun];
            flint_printf("%-5ld#%-3ld%-5ld%-5ld%-5ld%-8ld%-8ld", b, ifun, dim1, dim2, dim3, degl, degr);
            time_args targs = {dim1, dim2, dim3, degl, degr, n};
            tfun(targs, state);
        }

        flint_printf(" ");
        flint_printf("\n");
    }
    else if (argc == 6)  // nbits + fun + dim1 + dim2 + dim3
    {
        const slong b = atoi(argv[1]);
        slong ifun = atoi(argv[2]);
        const slong dim1 = atoi(argv[3]);
        const slong dim2 = atoi(argv[4]);
        const slong dim3 = atoi(argv[5]);
        flint_printf("bits fun dim1 dim2 dim3 degl    degr\n");

        ulong n;
        n = n_nextprime(UWORD(1) << (b-1), 0);

        if (ifun == -1)
        {
            for (ifun = 0; ifun < nfuns; ifun++)
            {
                const timefun tfun = funs[ifun];
                for (slong ideg = 0; ideg < ndegs; ideg++)
                {
                    const slong deg = degs[ideg];
                    flint_printf("%-5ld#%-3ld%-5ld%-5ld%-5ld%-8ld%-8ld", b, ifun, dim1, dim2, dim3, deg, deg);
                    time_args targs = {dim1, dim2, dim3, deg, deg, n};
                    tfun(targs, state);
                    flint_printf("\n");
                }
            }
        }
        else
        {
            const timefun tfun = funs[ifun];
            for (slong ideg = 0; ideg < ndegs; ideg++)
            {
                const slong deg = degs[ideg];
                flint_printf("%-5ld#%-3ld%-5ld%-5ld%-5ld%-8ld%-8ld", b, ifun, dim1, dim2, dim3, deg, deg);
                time_args targs = {dim1, dim2, dim3, deg, deg, n};
                tfun(targs, state);
            }
        }

        flint_printf(" ");
        flint_printf("\n");
    }


    flint_rand_clear(state);
    return 0;
}
