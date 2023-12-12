#include <assert.h>
#include <math.h>
#include <gmp.h>
#include <flint/flint.h>
#include <flint/nmod.h>
#include <flint/nmod_vec.h>

#include "nmod_mat_extra.h"

/*--------------------------------------------------------------*/
/* computes a square matrix product in size len modulo n        */
/*--------------------------------------------------------------*/
void time_nmod_mat_mul(ulong len, ulong n)
{
    flint_rand_t state;
    nmod_mat_t a, b, c1, c2;
    nmod_t mod;
    double t, rate;
    clock_t tt;
    long nb_iter;

    flint_randinit(state);
    nmod_init(&mod, n);

    nmod_mat_init(a, len, len, mod.n);
    nmod_mat_init(b, len, len, mod.n);
    nmod_mat_init(c1, len, len, mod.n);
    nmod_mat_init(c2, len, len, mod.n);

    nmod_mat_rand(a, state);
    nmod_mat_rand(b, state);
    nmod_mat_rand(c1, state);
    nmod_mat_rand(c2, state);

    rate = log2(n);
    printf("%lu %lu ", len, (ulong)rate);
    
    t = 0.0;
    nb_iter = 0;
    while (t < 0.5)
    {
        tt = clock();
        nmod_mat_mul(c1, a, b);
        nmod_mat_mul(c1, a, b);
        nmod_mat_mul(c1, a, b);
        nmod_mat_mul(c1, a, b);
        nmod_mat_mul(c1, a, b);
        t += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += 5;
    }
    t /= nb_iter;
    printf("%4g ", t/rate);

    
    t = 0.0;
    nb_iter = 0;
    while (t < 0.5)
    {
        tt = clock();
        nmod_mat_mul_classical(c1, a, b);
        nmod_mat_mul_classical(c1, a, b);
        nmod_mat_mul_classical(c1, a, b);
        nmod_mat_mul_classical(c1, a, b);
        nmod_mat_mul_classical(c1, a, b);
        t += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += 5;
    }
    t /= nb_iter;
    printf("%4g ", t/rate);

    t = 0.0;
    nb_iter = 0;
    while (t < 0.5)
    {
        tt = clock();
        nmod_mat_mul_blas(c1, a, b);
        nmod_mat_mul_blas(c1, a, b);
        nmod_mat_mul_blas(c1, a, b);
        nmod_mat_mul_blas(c1, a, b);
        nmod_mat_mul_blas(c1, a, b);
        t += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += 5;
    }
    t /= nb_iter;
    printf("%4g ", t/rate);

    if (n < (1L << 30))
    {
        t = 0.0;
        nb_iter = 0;
        while (t < 0.5)
        {
            tt = clock();
            nmod_mat_mul_small_modulus(c2, a, b);
            nmod_mat_mul_small_modulus(c2, a, b);
            nmod_mat_mul_small_modulus(c2, a, b);
            nmod_mat_mul_small_modulus(c2, a, b);
            nmod_mat_mul_small_modulus(c2, a, b);
            t += (double)(clock()-tt) / CLOCKS_PER_SEC;
            nb_iter += 5;
        }
        t /= nb_iter;
        printf("%4g ", t/rate);
    }

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
int main(int argc, char **argv)
{
    printf("#size n mul classical blas (small_modulus)\t\t (time per bit)\n");
    ulong i;

    i = 1;
    for (; i < 100; i += 5)
    {
        time_nmod_mat_mul(i, (1L << 29) + 1);
        time_nmod_mat_mul(i, (1L << 60) + 1);
        printf("\n");
    }
    /* for (; i < 200; i += 10) */
    /* { */
    /*     time_nmod_mat_mul(i, (1L << 29) + 1); */
    /*     time_nmod_mat_mul(i, (1L << 60) + 1); */
    /* } */

    return 0;
}
