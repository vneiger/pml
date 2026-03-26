#include <flint/flint.h>
#include <stdlib.h>  // for atoi

#include <flint/ulong_extras.h>
#include <flint/profiler.h>
#include <flint/nmod.h>
#include <flint/nmod_mat.h>

#define MEASURE_SAMPLE 0

/* for LU of dim1 x dim2 */
typedef struct {
    slong dim1;
    slong dim2;
    slong modn;
} time_args;

#define TIME_LU(fun) \
void time_##fun(time_args targs, flint_rand_t state) \
{ \
    const slong dim1 = targs.dim1; \
    const slong dim2 = targs.dim2; \
    const slong modn = targs.modn; \
    \
    nmod_mat_t mat; \
    nmod_mat_init(mat, dim1, dim2, modn); \
    nmod_mat_rand(mat, state); \
    slong * P = FLINT_ARRAY_ALLOC(dim1, slong); \
    \
    double FLINT_SET_BUT_UNUSED(tcpu), twall; \
    \
    TIMEIT_START; \
    fun(P, mat, 0); \
    TIMEIT_STOP_VALUES(tcpu, twall); \
    \
    flint_printf("%.2e", twall); \
    \
    flint_free(P); \
    nmod_mat_clear(mat); \
}

#define SAMPLE_LU(fun)                                          \
void sample_##fun(void * arg, ulong count)                      \
{                                                               \
    time_args * targs = (time_args *) arg;                      \
    const slong dim1 = targs->dim1;                             \
    const slong dim2 = targs->dim2;                             \
    const slong modn = targs->modn;                             \
                                                                \
    FLINT_TEST_INIT(state);                                     \
                                                                \
    nmod_mat_t mat;                                             \
    nmod_mat_init(mat, dim1, dim2, modn);                       \
    nmod_mat_rand(mat, state);                                  \
    slong * P = FLINT_ARRAY_ALLOC(dim1, slong);                 \
	for (ulong i = 0; i < count; i++)                           \
    {                                                           \
        nmod_mat_t tmp;                                         \
        nmod_mat_init(tmp, dim1, dim2, modn);                   \
        for (slong k = 0; k < dim1*dim2; k++)                   \
            tmp->entries[k] = mat->entries[k];                  \
        prof_start();                                           \
        fun(P, mat, 0);                                         \
        prof_stop();                                            \
        nmod_mat_clear(tmp);                                    \
    }                                                           \
                                                                \
    flint_free(P);                                              \
    nmod_mat_clear(mat);                                        \
    FLINT_TEST_CLEAR(state);                                    \
}

TIME_LU(nmod_mat_lu)
TIME_LU(nmod_mat_lu_classical)
TIME_LU(nmod_mat_lu_classical_delayed)
TIME_LU(nmod_mat_lu_recursive)

SAMPLE_LU(nmod_mat_lu)
SAMPLE_LU(nmod_mat_lu_classical)
SAMPLE_LU(nmod_mat_lu_classical_delayed)
SAMPLE_LU(nmod_mat_lu_recursive)

/*-------------------------*/
/*  main                   */
/*-------------------------*/
int main(int argc, char ** argv)
{
    flint_rand_t state;
    flint_rand_init(state);
    flint_rand_set_seed(state, time(NULL), time(NULL)+129384125L);

    // bench functions
    const slong nfuns = 4;
    typedef void (*timefun) (time_args, flint_rand_t);
    const timefun funs[] = {
        time_nmod_mat_lu,                     // 0
        time_nmod_mat_lu_classical,           // 1
        time_nmod_mat_lu_classical_delayed,   // 2
        time_nmod_mat_lu_recursive,           // 3
    };

    typedef void (*samplefun) (void*, ulong);
    const samplefun sfuns[] = {
        sample_nmod_mat_lu,                     // 0
        sample_nmod_mat_lu_classical,           // 1
        sample_nmod_mat_lu_classical_delayed,   // 2
        sample_nmod_mat_lu_recursive,           // 3
    };

    const char * description[] = {
        "#0  --> nmod_mat_lu                     ",
        "#1  --> nmod_mat_lu_classical           ",
        "#2  --> nmod_mat_lu_classical_delayed   ",
        "#3  --> nmod_mat_lu_recursive           ",
    };

    if (argc == 1)  // show usage
    {
        flint_printf("Usage: `%s [nbits] [fun] [dim1] [dim2]`\n", argv[0]);
        flint_printf("   No argument shows this help.\n");
        flint_printf("   - nbits: number of bits, in (0..63], for the randomly-chosen modulus\n");
        flint_printf("   - fun: id number of the timed function (see below),\n");
        flint_printf("   - dim1, dim2: matrix dimensions\n");
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
        const slong dim1 = atoi(argv[3]);
        const slong dim2 = atoi(argv[4]);
        flint_printf("bits fun dim1  dim2  \n");

        ulong n = n_nextprime(UWORD(1) << (b-1), 1);

        flint_printf("%-4ld %-3ld %-6ld%-6ld", b, ifun, dim1, dim2);
        time_args targs = {dim1, dim2, n};

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
