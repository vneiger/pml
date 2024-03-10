#include <time.h>
#include "nmod_poly_extra.h"

/*------------------------------------------------------------*/
/* computes FFTs/TFTs up to order 2^25                        */
/*------------------------------------------------------------*/
void get_time()
{
    flint_rand_t state;
    flint_randinit(state);

    // prime, modulus
    nmod_t mod;
    ulong p = 1108307720798209; // 50 bits, 1 + 2**44 * 3**2 * 7
    nmod_init(&mod, p);

    // fft context
    sd_fft_ctx_t Q;
    sd_fft_lctx_t QL;
    sd_fft_ctx_init_prime(Q, p);
    sd_fft_lctx_init(QL, Q, 27);

    // nmod_sd_fft prepare
    mp_limb_t w0 = nmod_pow_ui(n_primitive_root_prime(p), (p - 1) >> 27, mod);
    mp_limb_t w = nmod_pow_ui(w0, 1L<<(27-27), mod);
    nmod_sd_fft_t F;
    nmod_sd_fft_init_set(F, w, 27, mod);

    printf("- order is log(fft length)\n");
    printf("- timing various tft lengths, e.g. tft1.25 means length 1.25 * 2**order\n");
    printf("order\tfft1\ttft1\ttft1+\ttft1.25\ttft1.5\ttft1.75\ttft2-\n");

    for (ulong order = 2; order < 26; order++)
    {
        double t;
        clock_t tt;
        long nb_iter;

        ulong Ns[6] = {
            1L<<order,
            (1L<<order) + 1,
            (1L<<order) + (1L<<(order-2)),
            (1L<<order) + (1L<<(order-1)),
            (1L<<order) + 3*(1L<<(order-2)),
            (1L<<(order+1)) - 1,
        };

        printf("%ld\t", order);
        for (ulong k = 0; k < 6; k++)
        {
            ulong N = Ns[k];

            nmod_poly_t P;
            nmod_poly_init2(P, p, N);
            mp_ptr val = _nmod_vec_init(N);

            for (slong i = 0; i < (slong)N; i++)
                nmod_poly_set_coeff_ui(P, i, n_randtest(state) % p);

            if (k == 0)
            {
                t = 0.0;
                nb_iter = 0;
                while (t < 0.5)
                {
                    tt = clock();
                    nmod_sd_fft_evaluate(val, P, QL, order);
                    t += (double)(clock()-tt) / CLOCKS_PER_SEC;
                    ++nb_iter;
                }
                t /= nb_iter;
                printf("%.1e\t", t);
            }

            t = 0.0;
            nb_iter = 0;
            while (t < 0.5)
            {
                tt = clock();
                nmod_sd_tft_evaluate(val, P, QL, F, N);
                t += (double)(clock()-tt) / CLOCKS_PER_SEC;
                ++nb_iter;
            }
            t /= nb_iter;
            printf("%.1e\t", t);
            _nmod_vec_clear(val);
            nmod_poly_clear(P);
        }
        printf("\n");
    }

    sd_fft_lctx_clear(QL, Q);
    sd_fft_ctx_clear(Q);
    nmod_sd_fft_clear(F);
    flint_randclear(state);
}


/*------------------------------------------------------------*/
/* main just calls get_time()                                 */
/*------------------------------------------------------------*/
int main()
{
    get_time();
    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
