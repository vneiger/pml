#include <stdlib.h>  // for atoi

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

#define TIME_APPROX(fun)                                \
void time_##fun(time_args targs, flint_rand_t state)    \
{                                                       \
    const slong rdim = targs.rdim;                      \
    const slong cdim = targs.cdim;                      \
    const slong deg = targs.deg;                        \
    const slong order = targs.order;                    \
    /* const slong rank = targs.rank; */ /* TODO */     \
    /* const slong stype = targs.stype; */ /* TODO */   \
    const slong n = targs.modn;                         \
                                                        \
    nmod_t mod;                                         \
    nmod_init(&mod, n);                                 \
                                                        \
    nmod_poly_mat_t F;                                  \
    nmod_poly_mat_init(F, rdim, cdim, n);               \
    nmod_poly_mat_rand(F, state, deg);                  \
                                                        \
    slong * shift = FLINT_ARRAY_ALLOC(rdim, slong);     \
    for (slong i = 0; i < rdim; i++)                    \
        shift[i] = 0;                                   \
    nmod_poly_mat_t P;                                  \
    nmod_poly_mat_init(P, rdim, rdim, n);               \
                                                        \
    double FLINT_SET_BUT_UNUSED(tcpu), twall;           \
                                                        \
    TIMEIT_START;                                       \
    nmod_poly_mat_##fun(P, shift, F, order);            \
    TIMEIT_STOP_VALUES(tcpu, twall);                    \
                                                        \
    flint_printf("%.2e", twall);                        \
                                                        \
    nmod_poly_mat_clear(F);                             \
    nmod_poly_mat_clear(P);                             \
    flint_free(shift);                                  \
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
    flint_rand_set_seed(state, time(NULL), time(NULL)+129384125L);

    /* TODO rank */
    /* TODO shift type */

    // bench functions
    const slong nfuns = 2;
    typedef void (*timefun) (time_args, flint_rand_t);
    const timefun funs[] = {
        time_pmbasis,                      // 0
        time_mbasis,                       // 1
    };

    const char * description[] = {
        "#0  --> pmbasis (general interface)   ",
        "#1  --> mbasis (general interface)    ",
    };

    if (argc == 1 || argc != 7)  // show usage
    {
        flint_printf("Usage: `%s [nbits] [fun] [rdim] [cdim] [deg] [order] [shift] [rank]`\n", argv[0]);
        flint_printf("   No argument shows this help.\n");
        flint_printf("   [rank] and [shift] not supported yet.\n");
        flint_printf("   - nbits: number of bits in [2..64] for the modulus, chosen as nextprime(2**(nbits-1))\n");
        flint_printf("   - fun: id number of the timed function (see below),\n");
        flint_printf("   - rdim, cdim: input matrix is rdim x cdim\n");
        flint_printf("   - deg: matrix is random of degree < deg\n");
        flint_printf("   - rank: [unsupported] matrix is random of this rank\n");
        flint_printf("   - shift: [unsupported] type of shift\n");
        flint_printf("\nAvailable functions:\n");
        for (slong j = 0; j < nfuns; j++)
            flint_printf("   %s\n", description[j]);

        return 0;
    }

    flint_printf("#warmup...\n");
    for (slong i = 0; i < 3; i++)
    {
        /* rdim; cdim; deg; rank; stype; modn; */
        time_args targs = {8, 4, 1000, 1000, 4, 0, n_nextprime(UWORD(1) << 20, 0)};
        time_pmbasis(targs, state);
        flint_printf(" ");
    }
    flint_printf("\n\n");

    if (argc == 7)  // nbits, fun, rdim, cdim, deg, order
    {
        const slong b     = atoi(argv[1]);
        const slong ifun  = atoi(argv[2]);
        const slong rdim  = atoi(argv[3]);
        const slong cdim  = atoi(argv[4]);
        const slong deg   = atoi(argv[5]);
        const slong order = atoi(argv[6]);
        flint_printf("bits fun rdim cdim deg     order\n");
        ulong n = n_nextprime(UWORD(1) << (b-1), 0);
        const timefun tfun = funs[ifun];
        flint_printf("%-5ld#%-3ld%-5ld%-5ld%-8ld%-8ld", b, ifun, rdim, cdim, deg, order);
        /* rdim; cdim; deg; order; rank; stype; modn; */
        time_args targs = {rdim, cdim, deg, order, FLINT_MIN(rdim, cdim), 0, n};
        tfun(targs, state);
        flint_printf(" ");
        flint_printf("\n");
    }

    flint_rand_clear(state);
    return 0;
}
