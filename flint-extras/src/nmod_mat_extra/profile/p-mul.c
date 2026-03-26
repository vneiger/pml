#include <stdlib.h>  // for atoi

#include <flint/ulong_extras.h>
#include <flint/profiler.h>
#include <flint/nmod.h>
#include <flint/nmod_mat.h>

#define MEASURE_SAMPLE 0

/* for multiplying dim1 x dim2 x dim3 */
typedef struct {
    slong dim1;
    slong dim2;
    slong dim3;
    slong modn;
} time_args;

#define TIME_MUL(fun) \
void time_##fun(time_args targs, flint_rand_t state) \
{ \
    const slong dim1 = targs.dim1; \
    const slong dim2 = targs.dim2; \
    const slong dim3 = targs.dim3; \
    const slong modn = targs.modn; \
    \
    nmod_mat_t mat1; \
    nmod_mat_t mat2; \
    nmod_mat_t prod; \
    nmod_mat_init(mat1, dim1, dim2, modn); \
    nmod_mat_init(mat2, dim2, dim3, modn); \
    nmod_mat_init(prod, dim1, dim3, modn); \
    nmod_mat_rand(mat1, state); \
    nmod_mat_rand(mat2, state); \
    \
    double FLINT_SET_BUT_UNUSED(tcpu), twall; \
    \
    TIMEIT_START; \
    fun(prod, mat1, mat2); \
    TIMEIT_STOP_VALUES(tcpu, twall); \
    \
    flint_printf("%.2e", twall); \
    \
    nmod_mat_clear(mat1); \
    nmod_mat_clear(mat2); \
    nmod_mat_clear(prod); \
}

#define SAMPLE_MUL(fun)                                         \
void sample_##fun(void * arg, ulong count)                      \
{                                                               \
    time_args * targs = (time_args *) arg;                      \
    const slong dim1 = targs->dim1;                             \
    const slong dim2 = targs->dim2;                             \
    const slong dim3 = targs->dim3;                             \
    const slong modn = targs->modn;                             \
                                                                \
    FLINT_TEST_INIT(state);                                     \
                                                                \
    nmod_mat_t mat1;                                            \
    nmod_mat_t mat2;                                            \
    nmod_mat_t prod;                                            \
    nmod_mat_init(mat1, dim1, dim2, modn);                      \
    nmod_mat_init(mat2, dim2, dim3, modn);                      \
    nmod_mat_rand(mat1, state);                                 \
    nmod_mat_rand(mat2, state);                                 \
	for (ulong i = 0; i < count; i++)                           \
    {                                                           \
        nmod_mat_t tmp1;                                        \
        nmod_mat_t tmp2;                                        \
        nmod_mat_init(tmp1, dim1, dim2, modn);                  \
        nmod_mat_init(tmp2, dim2, dim3, modn);                  \
        nmod_mat_init(prod, dim1, dim3, modn);                  \
        for (slong k = 0; k < dim1*dim2; k++)                   \
            tmp1->entries[k] = mat1->entries[k];                \
        for (slong k = 0; k < dim2*dim3; k++)                   \
            tmp2->entries[k] = mat2->entries[k];                \
        prof_start();                                           \
        fun(prod, tmp1, tmp2);                                  \
        prof_stop();                                            \
        nmod_mat_clear(prod);                                   \
        nmod_mat_clear(tmp1);                                   \
        nmod_mat_clear(tmp2);                                   \
    }                                                           \
                                                                \
    nmod_mat_clear(mat1);                                       \
    nmod_mat_clear(mat2);                                       \
    FLINT_TEST_CLEAR(state);                                    \
}

TIME_MUL(nmod_mat_mul)
TIME_MUL(nmod_mat_mul_classical)
TIME_MUL(nmod_mat_mul_classical_threaded)
TIME_MUL(nmod_mat_mul_blas)
TIME_MUL(nmod_mat_mul_strassen)

SAMPLE_MUL(nmod_mat_mul)
SAMPLE_MUL(nmod_mat_mul_classical)
SAMPLE_MUL(nmod_mat_mul_classical_threaded)
SAMPLE_MUL(nmod_mat_mul_blas)
SAMPLE_MUL(nmod_mat_mul_strassen)

/*-------------------------*/
/*  main                   */
/*-------------------------*/
int main(int argc, char ** argv)
{
    flint_rand_t state;
    flint_rand_init(state);
    flint_rand_set_seed(state, time(NULL), time(NULL)+129384125L);

    // bench functions
    const slong nfuns = 5;
    typedef void (*timefun) (time_args, flint_rand_t);
    const timefun funs[] = {
        time_nmod_mat_mul,                     // 0
        time_nmod_mat_mul_classical,           // 1
        time_nmod_mat_mul_classical_threaded,  // 2
        time_nmod_mat_mul_blas,                // 3
        time_nmod_mat_mul_strassen,            // 4
    };

    typedef void (*samplefun) (void*, ulong);
    const samplefun sfuns[] = {
        sample_nmod_mat_mul,                     // 0
        sample_nmod_mat_mul_classical,           // 1
        sample_nmod_mat_mul_classical_threaded,  // 2
        sample_nmod_mat_mul_blas,                // 3
        sample_nmod_mat_mul_strassen,            // 4
    };

    const char * description[] = {
        "#0  --> nmod_mat_mul                     ",
        "#1  --> nmod_mat_mul_classical           ",
        "#2  --> nmod_mat_mul_classical_threaded  ",
        "#3  --> nmod_mat_mul_blas                ",
        "#4  --> nmod_mat_mul_strassen            ",
    };

    if (argc == 1)  // show usage
    {
        flint_printf("Usage: `%s [nbits] [fun] [dim1] [dim2] [dim3]`\n", argv[0]);
        flint_printf("   No argument shows this help.\n");
        flint_printf("   - nbits: number of bits, in (0..63], for the randomly-chosen modulus\n");
        flint_printf("   - fun: id number of the timed function (see below),\n");
        flint_printf("   - dim1, dim2, dim3: matrix dimensions\n");
        flint_printf("\nAvailable functions:\n");
        for (slong j = 0; j < nfuns; j++)
            flint_printf("   %s\n", description[j]);

        return 0;
    }

    /* flint_printf("#warmup...\n"); */
    for (slong i = 0; i < 3; i++)
    {
        time_args targs = {i*100, i*100, i*100, 17 + (UWORD(1) << (i*17))};
        double min, max;
        prof_repeat(&min, &max, sfuns[0], &targs);
    }
    /* flint_printf("\n\n"); */

    if (argc == 6)
    {
        const slong b = atoi(argv[1]);
        const slong ifun = atoi(argv[2]);
        const slong dim1 = atoi(argv[3]);
        const slong dim2 = atoi(argv[4]);
        const slong dim3 = atoi(argv[5]);
        flint_printf("bits fun dim1  dim2  dim3  \n");

        ulong n = n_nextprime(UWORD(1) << (b-1), 1);

        flint_printf("%-4ld %-3ld %-6ld%-6ld%-6ld", b, ifun, dim1, dim2, dim3);
        time_args targs = {dim1, dim2, dim3, n};

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
