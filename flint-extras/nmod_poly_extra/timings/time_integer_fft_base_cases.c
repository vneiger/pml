#include <time.h>
#include <flint/nmod.h>
#include <flint/nmod_poly.h>
#include "nmod_poly_extra.h"
#include "nmod_poly_fft.h"
#include "../src/nmod_integer_fft_base_cases.c"

#define VERSIONS 1

#define num_primes 5

/*------------------------------------------------------------*/
/* computes init for FFT for several bit lengths and orders   */
/*------------------------------------------------------------*/
void time_evaluate()
{
    flint_rand_t state;
    flint_randinit(state);

    printf("- order is log(fft length)\n");
    printf("- timing init FFT tables + DIF evaluate for several bit lengths and orders\n");
    printf("order\tdif\tred\n");

    ulong primes[num_primes] = {
        786433,              // 20 bits, 1 + 2**18 * 3
        2013265921,          // 31 bits, 1 + 2**27 * 3 * 5
        2748779069441,       // 42 bits, 1 + 2**39 * 5
        1108307720798209,    // 50 bits, 1 + 2**44 * 3**2 * 7
        1139410705724735489, // 60 bits, 1 + 2**52 * 11 * 23
    };
    ulong max_orders[num_primes] = { 3, 3, 3, 3, 3 };

    for (ulong k = 0; k < 5; k++)
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

            const ulong len = (1<<order);

            // root of unity of order 2**order
            mp_limb_t w = nmod_pow_ui(w0, 1UL<<(max_orders[k]-order), mod);

            double t;
            clock_t tt;
            long nb_iter;

            if (VERSIONS >= 1)
            { // dft8_dif_lazy
                nmod_integer_fft_t Fpre;
                nmod_integer_fft_init_set_pre(Fpre, w, order, mod);
                t = 0.0;
                nb_iter = 0;
                while (t < 500)
                {
                    nmod_poly_t pol;
                    nmod_poly_init(pol, mod.n);
                    nmod_poly_rand(pol, state, len);
                    tt = clock();
                    dft8_dif_lazy(pol->coeffs, Fpre);
                    dft8_dif_lazy(pol->coeffs, Fpre);
                    dft8_dif_lazy(pol->coeffs, Fpre);
                    dft8_dif_lazy(pol->coeffs, Fpre);
                    dft8_dif_lazy(pol->coeffs, Fpre);
                    dft8_dif_lazy(pol->coeffs, Fpre);
                    dft8_dif_lazy(pol->coeffs, Fpre);
                    dft8_dif_lazy(pol->coeffs, Fpre);
                    dft8_dif_lazy(pol->coeffs, Fpre);
                    dft8_dif_lazy(pol->coeffs, Fpre);
                    t += (double)(clock()-tt) / 1000;
                    nb_iter+=10;
                    nmod_poly_clear(pol);
                }
                t /= nb_iter;
                nmod_integer_fft_clear_pre(Fpre);
                printf("%.1e\t", t);
            }

            if (VERSIONS >= 1)
            { // dft8_red_lazy
                nmod_integer_fft_t Fpre;
                nmod_integer_fft_init_set_pre(Fpre, w, order, mod);
                t = 0.0;
                nb_iter = 0;
                while (t < 500)
                {
                    nmod_poly_t pol;
                    nmod_poly_init(pol, mod.n);
                    nmod_poly_rand(pol, state, len);
                    tt = clock();
                    dft8_red_lazy(pol->coeffs, Fpre);
                    dft8_red_lazy(pol->coeffs, Fpre);
                    dft8_red_lazy(pol->coeffs, Fpre);
                    dft8_red_lazy(pol->coeffs, Fpre);
                    dft8_red_lazy(pol->coeffs, Fpre);
                    dft8_red_lazy(pol->coeffs, Fpre);
                    dft8_red_lazy(pol->coeffs, Fpre);
                    dft8_red_lazy(pol->coeffs, Fpre);
                    dft8_red_lazy(pol->coeffs, Fpre);
                    dft8_red_lazy(pol->coeffs, Fpre);
                    t += (double)(clock()-tt) / 1000;
                    nb_iter+=10;
                    nmod_poly_clear(pol);
                }
                t /= nb_iter;
                nmod_integer_fft_clear_pre(Fpre);
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
    time_evaluate();
    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
