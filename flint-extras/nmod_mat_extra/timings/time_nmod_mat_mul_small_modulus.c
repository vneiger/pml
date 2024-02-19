#include <time.h>
#include <flint/nmod.h>

#include "nmod_mat_extra.h"

/*--------------------------------------------------------------*/
/* computes a square matrix product in size len modulo n        */
/*--------------------------------------------------------------*/
void time_nmod_mat_mul_small_modulus(ulong len, ulong n)
{
    flint_rand_t state;
    nmod_mat_t a, b, c1, c2;
    nmod_t mod;
    double t;
    clock_t tt;
    long nb_iter;

    flint_randinit(state);
    nmod_init(&mod, n);

    nmod_mat_init(a, len, len, mod.n);
    nmod_mat_init(b, len, len, mod.n);
    nmod_mat_init(c1, len, len, mod.n);
    nmod_mat_init(c2, len, len, mod.n);

    printf("%lu\t", len);

    t = 0.0;
    nb_iter = 0;
    while (t < 0.5)
    {
        nmod_mat_rand(a, state);
        nmod_mat_rand(b, state);
        nmod_mat_rand(c1, state);
        tt = clock();
        nmod_mat_mul(c1, a, b);
        nmod_mat_mul(c1, a, b);
        nmod_mat_mul(c1, a, b);
        nmod_mat_mul(c1, a, b);
        nmod_mat_mul(c1, a, b);
        t += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += 5;
    }
    t = 1000 * t;
    t /= nb_iter;
    printf("%lf\t", t);

    t = 0.0;
    nb_iter = 0;
    while (t < 0.5)
    {
        nmod_mat_rand(a, state);
        nmod_mat_rand(b, state);
        nmod_mat_rand(c1, state);
        tt = clock();
        nmod_mat_mul_blas(c1, a, b);
        nmod_mat_mul_blas(c1, a, b);
        nmod_mat_mul_blas(c1, a, b);
        nmod_mat_mul_blas(c1, a, b);
        nmod_mat_mul_blas(c1, a, b);
        t += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += 5;
    }
    t = 1000 * t;
    t /= nb_iter;
    printf("%lf\t", t);

    t = 0.0;
    nb_iter = 0;
    while (t < 0.5)
    {
        nmod_mat_rand(a, state);
        nmod_mat_rand(b, state);
        nmod_mat_rand(c2, state);
        tt = clock();
        nmod_mat_mul_small_modulus(c2, a, b);
        nmod_mat_mul_small_modulus(c2, a, b);
        nmod_mat_mul_small_modulus(c2, a, b);
        nmod_mat_mul_small_modulus(c2, a, b);
        nmod_mat_mul_small_modulus(c2, a, b);
        t += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += 5;
    }
    t = 1000 * t;
    t /= nb_iter;
    printf("%lf", t);

    printf("\n");

    nmod_mat_clear(a);
    nmod_mat_clear(b);
    nmod_mat_clear(c1);
    nmod_mat_clear(c2);
    flint_randclear(state);
}

/*--------------------------------------------------------------*/
/* main calls time                                              */
/*--------------------------------------------------------------*/
int main()
{
    flint_set_num_threads(1);
    printf("Modulus %ld, times in ms\n", 1L << 29);
    printf("#size\tt_flint\t\tt_blas\t\tt_here\n");
    ulong i;
    for (i = 1; i < 1000; i += 10)
        time_nmod_mat_mul_small_modulus(i, (1L << 29) + 1);

    return 0;
}
