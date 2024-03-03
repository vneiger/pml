#include <time.h>
#include <gmp.h>
#include <flint/flint.h>
#include <flint/nmod.h>
#include <flint/nmod_vec.h>

#include "nmod_mat_extra.h"

/*--------------------------------------------------------------*/
/* computes a square matrix product in size len modulo n        */
/*--------------------------------------------------------------*/
void time_nmod_mat_mul(ulong len, ulong nbits, ulong n, flint_rand_t state)
{
    printf("%lu\t%lu\t", nbits, len);

    nmod_mat_t a, b, c;
    nmod_t mod;

    nmod_init(&mod, n);

    nmod_mat_init(a, len, len, mod.n);
    nmod_mat_init(b, len, len, mod.n);
    nmod_mat_init(c, len, len, mod.n);

    nmod_mat_rand(a, state);
    nmod_mat_rand(b, state);
    nmod_mat_rand(c, state);

    double t;
    clock_t tt;
    long nb_iter;

    t = 0.0;
    nb_iter = 0;
    while (t < 0.5)
    {
        tt = clock();
        nmod_mat_mul(c, a, b);
        nmod_mat_mul(c, a, b);
        nmod_mat_mul(c, a, b);
        nmod_mat_mul(c, a, b);
        nmod_mat_mul(c, a, b);
        t += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += 5;
    }
    t = 1000*t;
    t /= nb_iter;
    printf("%.1e\t", t);

    t = 0.0;
    nb_iter = 0;
    while (t < 0.5)
    {
        tt = clock();
        nmod_mat_mul_classical(c, a, b);
        nmod_mat_mul_classical(c, a, b);
        nmod_mat_mul_classical(c, a, b);
        nmod_mat_mul_classical(c, a, b);
        nmod_mat_mul_classical(c, a, b);
        t += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += 5;
    }
    t = 1000*t;
    t /= nb_iter;
    printf("%.1e\t", t);

    t = 0.0;
    nb_iter = 0;
    while (t < 0.5)
    {
        tt = clock();
        nmod_mat_mul_blas(c, a, b);
        nmod_mat_mul_blas(c, a, b);
        nmod_mat_mul_blas(c, a, b);
        nmod_mat_mul_blas(c, a, b);
        nmod_mat_mul_blas(c, a, b);
        t += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += 5;
    }
    t = 1000*t;
    t /= nb_iter;
    printf("%.1e\t", t);

    t = 0.0;
    nb_iter = 0;
    while (t < 0.5)
    {
        tt = clock();
        nmod_mat_mul_newdot(c, a, b);
        nmod_mat_mul_newdot(c, a, b);
        nmod_mat_mul_newdot(c, a, b);
        nmod_mat_mul_newdot(c, a, b);
        nmod_mat_mul_newdot(c, a, b);
        t += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += 5;
    }
    t = 1000*t;
    t /= nb_iter;
    printf("%.1e\t", t);

    if (n < (UWORD(1) << 30))
    {
        t = 0.0;
        nb_iter = 0;
        while (t < 0.5)
        {
            tt = clock();
            nmod_mat_mul_small_modulus(c, a, b);
            nmod_mat_mul_small_modulus(c, a, b);
            nmod_mat_mul_small_modulus(c, a, b);
            nmod_mat_mul_small_modulus(c, a, b);
            nmod_mat_mul_small_modulus(c, a, b);
            t += (double)(clock()-tt) / CLOCKS_PER_SEC;
            nb_iter += 5;
        }
        t = 1000*t;
        t /= nb_iter;
        printf("%.1e\t", t);
    }
    else
        printf("xxxxx\t");

    printf("\n");
    fflush(stdout);

    nmod_mat_clear(a);
    nmod_mat_clear(b);
    nmod_mat_clear(c);
}

/*--------------------------------------------------------------*/
/* main calls time                                              */
/*--------------------------------------------------------------*/
int main()
{
    flint_rand_t state;
    flint_randinit(state);

    printf("Square matrix multiplication over nmod:\n");
    printf("  - mul: flint's nmod_mat_mul\n");
    printf("  - classic: flint's classical cubic mul\n");
    printf("  - blas: flint's BLAS-based mul\n");
    printf("  - newdot: pml's mul using new nmod_vec_dot_product\n");
    printf("  - small_mod: pml's mul for small moduli using AVX-based nmod_vec_dot_product\n\n");
    const ulong mods[11] =
    {
             (1L << 3) + 1,
             (1L << 10) + 1,
             (1L << 20) + 1,
             (1L << 25) + 1,
             (1L << 29) + 1,
             (1L << 30) + 1,
             (1L << 31) + 1,
             (1L << 40) + 1,
             (1L << 50) + 1,
             (1L << 60) + 1,
             UWORD_MAX
    };
    const ulong nbits[11] = { 4, 11, 21, 26, 30, 31, 32, 41, 51, 61, 64 };

    for (slong i = 0; i < 11; i++)
    {
        printf("nbits\tsize\tmul\tclassic\tblas\tnewdot\tsmall_mod\n");
        //ulong len = 1;
        //for (; len < 100; len += 5)
        //    time_nmod_mat_mul(len, nbits[i], mods[i], state);
        //for (; len < 200; len += 20)
        //    time_nmod_mat_mul(len, nbits[i], mods[i], state);
        for (ulong len = 200; len < 2000; len += 50)
            time_nmod_mat_mul(len, nbits[i], mods[i], state);
        printf("\n");
    }

    flint_randclear(state);
    return 0;
}
