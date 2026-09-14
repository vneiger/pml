#include <stdlib.h>  // for atoi, strtoul
#include <string.h>  // for strtok

#include <flint/ulong_extras.h>
#include <flint/profiler.h>
#include <flint/nmod.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_poly_mat.h>

#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_approximant.h"

/* Shift types are:
 * 0 -> uniform, shift == (0,...,0)
 * */

/* one of the primes of FLINT's default fft_small context: 50 bits, and
   enough 2-adicity for the transform lengths reachable here; the
   polynomial matrix products inside the approximant basis algorithms then
   run on a single transform modulo the prime itself */
#define FFT_PRIME UWORD(1108307720798209)

typedef struct
{
    slong rdim;  /* row dimension */
    slong cdim;  /* column dimension */
    slong deg;   /* degree */
    slong order; /* approx order */
    slong rank;  /* rank */
    slong stype; /* shift type */
    slong modn;  /* modulus */
}
time_args;

/* wall and cpu time of one call on the given input; cpu is that of the
   whole process, so cpu/wall is the number of cores effectively busy */
#define TIME_APPROX(fun)                                            \
void time_##fun(time_args targs, const nmod_poly_mat_t F,           \
                double * tcpu, double * twall)                      \
{                                                                   \
    const slong rdim = targs.rdim;                                  \
    const slong order = targs.order;                                \
    /* const slong rank = targs.rank; */ /* TODO */                 \
    /* const slong stype = targs.stype; */ /* TODO */               \
    const slong n = targs.modn;                                     \
                                                                    \
    slong * shift = FLINT_ARRAY_ALLOC(rdim, slong);                 \
    for (slong i = 0; i < rdim; i++)                                \
        shift[i] = 0;                                               \
    nmod_poly_mat_t P;                                              \
    nmod_poly_mat_init(P, rdim, rdim, n);                           \
                                                                    \
    TIMEIT_START;                                                   \
    nmod_poly_mat_##fun(P, shift, F, order);                        \
    TIMEIT_STOP_VALUES(*tcpu, *twall);                              \
                                                                    \
    nmod_poly_mat_clear(P);                                         \
    flint_free(shift);                                              \
}

TIME_APPROX(pmbasis)
TIME_APPROX(mbasis)

/*-------------------------*/
/*  main                   */
/*-------------------------*/

int main(int argc, char ** argv)
{
    flint_rand_t state;
    flint_rand_init(state);
    /* a fixed seed: the same input across thread counts and across runs */
    flint_rand_set_seed(state, 1234, 5678);

    /* TODO rank */
    /* TODO shift type */

    // bench functions
    const slong nfuns = 2;
    typedef void (*timefun) (time_args, const nmod_poly_mat_t, double *, double *);
    const timefun funs[] = {
        time_pmbasis,                      // 0
        time_mbasis,                       // 1
    };

    const char * description[] = {
        "#0  --> pmbasis (general interface)   ",
        "#1  --> mbasis (general interface)    ",
    };

    if (argc < 7 || argc > 8)  // show usage
    {
        flint_printf("Usage: `%s [nbits] [fun] [rdim] [cdim] [deg] [order] [nthreads]`\n", argv[0]);
        flint_printf("   No argument shows this help.\n");
        flint_printf("   - nbits: number of bits in [2..64] for the modulus, chosen as nextprime(2**(nbits-1));\n");
        flint_printf("            0 selects a 50-bit FFT prime; a value above 64 is taken as the modulus itself\n");
        flint_printf("   - fun: id number of the timed function (see below),\n");
        flint_printf("   - rdim, cdim: input matrix is rdim x cdim\n");
        flint_printf("   - deg: matrix is random of degree < deg\n");
        flint_printf("   - order: order of approximation\n");
        flint_printf("   - nthreads: optional, default 1; a comma-separated list, e.g. 1,2,4,8, times\n");
        flint_printf("            the same input at each count and prints the speedups against the first\n");
        flint_printf("   - rank, shift: [unsupported yet]\n");
        flint_printf("\nAvailable functions:\n");
        for (slong j = 0; j < nfuns; j++)
            flint_printf("   %s\n", description[j]);

        return 0;
    }

    // nbits, fun, rdim, cdim, deg, order, [nthreads]
    {
        const ulong b     = strtoul(argv[1], NULL, 10);
        const slong ifun  = atoi(argv[2]);
        const slong rdim  = atoi(argv[3]);
        const slong cdim  = atoi(argv[4]);
        const slong deg   = atoi(argv[5]);
        const slong order = atoi(argv[6]);
        const ulong n = (b == 0) ? FFT_PRIME
                      : (b <= 64) ? n_nextprime(UWORD(1) << (b-1), 0) : b;
        const timefun tfun = funs[ifun];
        /* rdim; cdim; deg; order; rank; stype; modn; */
        time_args targs = {rdim, cdim, deg, order, FLINT_MIN(rdim, cdim), 0, n};

        /* the thread counts */
        slong nthreads[64], ncounts = 0;
        if (argc == 8)
        {
            char * tok;
            for (tok = strtok(argv[7], ", "); tok != NULL && ncounts < 64; tok = strtok(NULL, ", "))
                nthreads[ncounts++] = atol(tok);
        }
        if (ncounts == 0)
            nthreads[ncounts++] = 1;

        nmod_poly_mat_t F;
        nmod_poly_mat_init(F, rdim, cdim, n);
        nmod_poly_mat_rand(F, state, deg);

        flint_printf("modulus %wu (%wu bits), fun #%wd, rdim %wd, cdim %wd, deg %wd, order %wd\n",
                     n, FLINT_BIT_COUNT(n), ifun, rdim, cdim, deg, order);
        flint_printf("%8s %10s %9s %8s\n", "nthreads", "wall", "cpu/wall", "speedup");

        double twall1 = 0.0;
        for (slong c = 0; c < ncounts; c++)
        {
            double tcpu, twall;

            flint_set_num_threads(nthreads[c]);
            tfun(targs, F, &tcpu, &twall);
            if (c == 0)
                twall1 = twall;

            flint_printf("%8wd %10.2e %9.2f %8.2f\n", nthreads[c], twall,
                         twall > 0 ? tcpu / twall : 0.0,
                         twall > 0 ? twall1 / twall : 0.0);
        }

        nmod_poly_mat_clear(F);
    }

    flint_rand_clear(state);
    return 0;
}
