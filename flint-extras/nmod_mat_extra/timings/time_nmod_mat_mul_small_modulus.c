#include <assert.h>
#include <gmp.h>
#include <flint/flint.h>
#include <flint/nmod.h>
#include <flint/nmod_vec.h>

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

    nmod_mat_rand(a, state);
    nmod_mat_rand(b, state);
    nmod_mat_rand(c1, state);
    nmod_mat_rand(c2, state);

    printf("%lu ", len);
    
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
    t = 1000 * t;
    t /= nb_iter;
    printf("%lf", t/1000);

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
    t = 1000 * t;
    t /= nb_iter;
    printf(" %lf", t/1000);

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
    printf("#size t_old t_new\n");
    ulong i;
    for (i = 1; i < 110; i += 10)
        time_nmod_mat_mul_small_modulus(i, (1L << 29) + 1);

    return 0;
}
