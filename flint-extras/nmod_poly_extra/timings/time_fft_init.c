#include <time.h>
#include <flint/nmod.h>
#include <flint/nmod_poly.h>
#include "nmod_poly_extra.h"
#include "nmod_poly_fft.h"

#define num_primes 5

/*-----------------------------------------------------------------*/
/* initialize context for FFT for several bit lengths and orders   */
/*-----------------------------------------------------------------*/
void time_fft_init(ulong * primes, ulong * max_orders, flint_rand_t state)
{
    for (ulong k = 4; k < num_primes; k++)
    {
        // prime, modulus
        nmod_t mod;
        ulong p = primes[k];
        nmod_init(&mod, p);

        // find root of unity of specified maximum order
        ulong w0 = nmod_pow_ui(n_primitive_root_prime(p), (p - 1) >> max_orders[k], mod);

        for (ulong order = 3; order <= max_orders[k]; order++)
        {
            printf("%ld\t", order);

            // root of unity of order 2**order
            ulong w = nmod_pow_ui(w0, 1UL<<(max_orders[k]-order), mod);

            double t;
            clock_t tt;
            long nb_iter;

            { // fft_init
                nmod_fft_ctx_t F;
                t = 0.0;
                nb_iter = 0;
                while (t < 0.5)
                {
                    tt = clock();
                    nmod_fft_ctx_init_set(F, w, order, mod);
                    nmod_fft_ctx_clear(F);
                    nmod_fft_ctx_init_set(F, w, order, mod);
                    nmod_fft_ctx_clear(F);
                    nmod_fft_ctx_init_set(F, w, order, mod);
                    nmod_fft_ctx_clear(F);
                    nmod_fft_ctx_init_set(F, w, order, mod);
                    nmod_fft_ctx_clear(F);
                    nmod_fft_ctx_init_set(F, w, order, mod);
                    nmod_fft_ctx_clear(F);
                    nmod_fft_ctx_init_set(F, w, order, mod);
                    nmod_fft_ctx_clear(F);
                    nmod_fft_ctx_init_set(F, w, order, mod);
                    nmod_fft_ctx_clear(F);
                    nmod_fft_ctx_init_set(F, w, order, mod);
                    nmod_fft_ctx_clear(F);
                    nmod_fft_ctx_init_set(F, w, order, mod);
                    nmod_fft_ctx_clear(F);
                    nmod_fft_ctx_init_set(F, w, order, mod);
                    nmod_fft_ctx_clear(F);
                    t += (double)(clock()-tt) / CLOCKS_PER_SEC;
                    nb_iter+=10;
                }
                t /= nb_iter;
                printf("%.1e\t", t);
            }

            { // fft_init_red
                nmod_fft_ctx_t Fred;
                t = 0.0;
                nb_iter = 0;
                while (t < 0.5)
                {
                    tt = clock();
                    nmod_fft_ctx_init_set_red(Fred, w, order, mod);
                    nmod_fft_ctx_clear_red(Fred);
                    nmod_fft_ctx_init_set_red(Fred, w, order, mod);
                    nmod_fft_ctx_clear_red(Fred);
                    nmod_fft_ctx_init_set_red(Fred, w, order, mod);
                    nmod_fft_ctx_clear_red(Fred);
                    nmod_fft_ctx_init_set_red(Fred, w, order, mod);
                    nmod_fft_ctx_clear_red(Fred);
                    nmod_fft_ctx_init_set_red(Fred, w, order, mod);
                    nmod_fft_ctx_clear_red(Fred);
                    nmod_fft_ctx_init_set_red(Fred, w, order, mod);
                    nmod_fft_ctx_clear_red(Fred);
                    nmod_fft_ctx_init_set_red(Fred, w, order, mod);
                    nmod_fft_ctx_clear_red(Fred);
                    nmod_fft_ctx_init_set_red(Fred, w, order, mod);
                    nmod_fft_ctx_clear_red(Fred);
                    nmod_fft_ctx_init_set_red(Fred, w, order, mod);
                    nmod_fft_ctx_clear_red(Fred);
                    nmod_fft_ctx_init_set_red(Fred, w, order, mod);
                    nmod_fft_ctx_clear_red(Fred);
                    t += (double)(clock()-tt) / CLOCKS_PER_SEC;
                    nb_iter+=10;
                }
                t /= nb_iter;
                printf("%.1e\t", t);
            }

            { // fft_init_new
                nmod_fft_ctx_t F;
                t = 0.0;
                nb_iter = 0;
                while (t < 0.5)
                {
                    tt = clock();
                    nmod_fft_ctx_init_set_new(F, w, order, mod);
                    nmod_fft_ctx_clear(F);
                    nmod_fft_ctx_init_set_new(F, w, order, mod);
                    nmod_fft_ctx_clear(F);
                    nmod_fft_ctx_init_set_new(F, w, order, mod);
                    nmod_fft_ctx_clear(F);
                    nmod_fft_ctx_init_set_new(F, w, order, mod);
                    nmod_fft_ctx_clear(F);
                    nmod_fft_ctx_init_set_new(F, w, order, mod);
                    nmod_fft_ctx_clear(F);
                    nmod_fft_ctx_init_set_new(F, w, order, mod);
                    nmod_fft_ctx_clear(F);
                    nmod_fft_ctx_init_set_new(F, w, order, mod);
                    nmod_fft_ctx_clear(F);
                    nmod_fft_ctx_init_set_new(F, w, order, mod);
                    nmod_fft_ctx_clear(F);
                    nmod_fft_ctx_init_set_new(F, w, order, mod);
                    nmod_fft_ctx_clear(F);
                    nmod_fft_ctx_init_set_new(F, w, order, mod);
                    nmod_fft_ctx_clear(F);
                    t += (double)(clock()-tt) / CLOCKS_PER_SEC;
                    nb_iter+=10;
                }
                t /= nb_iter;
                printf("%.1e\t", t);
            }

            printf("\n");
        }
    }

}

/*------------------------------------------------------------*/
/* main just calls time_init_set()                            */
/*------------------------------------------------------------*/
int main()
{
    flint_rand_t state;
    flint_rand_init(state);

    printf("- order is log(fft length)\n");
    printf("- timing init FFT context at this order\n");
    printf("order\tdif\tred\n");

    ulong primes[num_primes] = {
        786433,              // 20 bits, 1 + 2**18 * 3
        2013265921,          // 31 bits, 1 + 2**27 * 3 * 5
        2748779069441,       // 42 bits, 1 + 2**39 * 5
        1108307720798209,    // 50 bits, 1 + 2**44 * 3**2 * 7
        1139410705724735489, // 60 bits, 1 + 2**52 * 11 * 23
    };
    ulong max_orders[num_primes] = { 18, 25, 25, 25, 25 };

    time_fft_init(primes, max_orders, state);

    flint_rand_clear(state);

    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
