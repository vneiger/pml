#include <time.h>
#include <flint/nmod.h>
#include "nmod_poly_fft.h"

#define num_primes 5

/*------------------------------------------------------------*/
/* computes init for FFT for several bit lengths and orders   */
/*------------------------------------------------------------*/
void time_init_set()
{
    flint_rand_t state;
    flint_randinit(state);

    printf("- order is log(fft length)\n");
    printf("- timing init of FFt tables for several bit lengths and orders\n");

    ulong primes[num_primes] = {
        1108307720798209, // 20 bits, 1 + 2**44 * 3**2 * 7
        1108307720798209, // 30 bits, 1 + 2**44 * 3**2 * 7
        1108307720798209, // 40 bits, 1 + 2**44 * 3**2 * 7
        1108307720798209, // 50 bits, 1 + 2**44 * 3**2 * 7
        1108307720798209, // 60 bits, 1 + 2**44 * 3**2 * 7
    };
    ulong max_orders[num_primes] = { 26, 26, 26, 26, 26 };

    for (ulong k = 0; k < num_primes; k++)
    {
        // prime, modulus
        nmod_t mod;
        ulong p = primes[k];
        nmod_init(&mod, p);

        // find root of unity of specified maximum order
        mp_limb_t w0 = nmod_pow_ui(n_primitive_root_prime(p), (p - 1) >> max_orders[k], mod);

        for (ulong order = 3; order <= max_orders[k]; order++)
        {
            printf("%ld\t", order);

            // root of unity of order 2**order
            mp_limb_t w = nmod_pow_ui(w0, 1UL<<(max_orders[k]-order), mod);

            double t;
            clock_t tt;
            long nb_iter;

            { // current version
                t = 0.0;
                nb_iter = 0;
                while (t < 0.5)
                {
                    tt = clock();
                    nmod_integer_fft_t F;
                    nmod_integer_fft_init_set(F, w, order, mod);
                    nmod_integer_fft_clear(F);
                    nmod_integer_fft_init_set(F, w, order, mod);
                    nmod_integer_fft_clear(F);
                    nmod_integer_fft_init_set(F, w, order, mod);
                    nmod_integer_fft_clear(F);
                    nmod_integer_fft_init_set(F, w, order, mod);
                    nmod_integer_fft_clear(F);
                    nmod_integer_fft_init_set(F, w, order, mod);
                    nmod_integer_fft_clear(F);
                    nmod_integer_fft_init_set(F, w, order, mod);
                    nmod_integer_fft_clear(F);
                    nmod_integer_fft_init_set(F, w, order, mod);
                    nmod_integer_fft_clear(F);
                    nmod_integer_fft_init_set(F, w, order, mod);
                    nmod_integer_fft_clear(F);
                    nmod_integer_fft_init_set(F, w, order, mod);
                    nmod_integer_fft_clear(F);
                    nmod_integer_fft_init_set(F, w, order, mod);
                    nmod_integer_fft_clear(F);
                    t += (double)(clock()-tt) / CLOCKS_PER_SEC;
                    nb_iter+=10;
                }
                t /= nb_iter;
                printf("%.1e\t", t);
            }

            { // v2
                t = 0.0;
                nb_iter = 0;
                while (t < 0.5)
                {
                    tt = clock();
                    nmod_integer_fft_t F;
                    nmod_integer_fft_init_set2(F, w, order, mod);
                    nmod_integer_fft_clear(F);
                    nmod_integer_fft_init_set2(F, w, order, mod);
                    nmod_integer_fft_clear(F);
                    nmod_integer_fft_init_set2(F, w, order, mod);
                    nmod_integer_fft_clear(F);
                    nmod_integer_fft_init_set2(F, w, order, mod);
                    nmod_integer_fft_clear(F);
                    nmod_integer_fft_init_set2(F, w, order, mod);
                    nmod_integer_fft_clear(F);
                    nmod_integer_fft_init_set2(F, w, order, mod);
                    nmod_integer_fft_clear(F);
                    nmod_integer_fft_init_set2(F, w, order, mod);
                    nmod_integer_fft_clear(F);
                    nmod_integer_fft_init_set2(F, w, order, mod);
                    nmod_integer_fft_clear(F);
                    nmod_integer_fft_init_set2(F, w, order, mod);
                    nmod_integer_fft_clear(F);
                    nmod_integer_fft_init_set2(F, w, order, mod);
                    nmod_integer_fft_clear(F);
                    t += (double)(clock()-tt) / CLOCKS_PER_SEC;
                    nb_iter+=10;
                }
                t /= nb_iter;
                printf("%.1e\t", t);
            }

            { // v3
                t = 0.0;
                nb_iter = 0;
                while (t < 0.5)
                {
                    tt = clock();
                    nmod_integer_fft_t F;
                    nmod_integer_fft_init_set3(F, w, order, mod);
                    nmod_integer_fft_clear(F);
                    nmod_integer_fft_init_set3(F, w, order, mod);
                    nmod_integer_fft_clear(F);
                    nmod_integer_fft_init_set3(F, w, order, mod);
                    nmod_integer_fft_clear(F);
                    nmod_integer_fft_init_set3(F, w, order, mod);
                    nmod_integer_fft_clear(F);
                    nmod_integer_fft_init_set3(F, w, order, mod);
                    nmod_integer_fft_clear(F);
                    nmod_integer_fft_init_set3(F, w, order, mod);
                    nmod_integer_fft_clear(F);
                    nmod_integer_fft_init_set3(F, w, order, mod);
                    nmod_integer_fft_clear(F);
                    nmod_integer_fft_init_set3(F, w, order, mod);
                    nmod_integer_fft_clear(F);
                    nmod_integer_fft_init_set3(F, w, order, mod);
                    nmod_integer_fft_clear(F);
                    nmod_integer_fft_init_set3(F, w, order, mod);
                    nmod_integer_fft_clear(F);
                    t += (double)(clock()-tt) / CLOCKS_PER_SEC;
                    nb_iter+=10;
                }
                t /= nb_iter;
                printf("%.1e\t", t);
            }
            printf("\n");
        }
    }

    flint_randclear(state);
}


/*------------------------------------------------------------*/
/* main just calls time_init_set()                            */
/*------------------------------------------------------------*/
int main()
{
    time_init_set();
    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
