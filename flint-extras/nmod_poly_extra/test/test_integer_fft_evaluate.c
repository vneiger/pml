#include <flint/nmod_vec.h>
#include <time.h>
#include <flint/nmod.h>
#include "nmod_poly_extra.h"
#include "nmod_poly_fft.h"

#define num_primes 5

/***********************
*  bit reversed copy  *
***********************/

long RevInc(long a, long k)
{
    long j, m;

    j = k;
    m = 1L << (k-1);

    while (j && (m & a)) {
        a ^= m;
        m >>= 1;
        j--;
    }
    if (j) a ^= m;
    return a;
}

// indices initialized with length >= k
void brc_indices(mp_limb_t * indices, long k)
{
    const long n = (1L << k);
    for (long i = 0, j = 0; i < n; i++, j = RevInc(j, k))
        indices[i] = j;
}


/*------------------------------------------------------------*/
/* computes init for FFT for several bit lengths and orders   */
/*------------------------------------------------------------*/
void test_fft_eval()
{
    flint_rand_t state;
    flint_randinit(state);

    printf("- order is log(fft length)\n");
    printf("- testing fft-eval for several bit lengths and orders\n");

    ulong primes[num_primes] = {
        786433,              // 20 bits, 1 + 2**18 * 3
        2013265921,          // 31 bits, 1 + 2**27 * 3 * 5
        2748779069441,       // 42 bits, 1 + 2**39 * 5
        1108307720798209,    // 50 bits, 1 + 2**44 * 3**2 * 7
        1139410705724735489, // 60 bits, 1 + 2**52 * 11 * 23
    };
    ulong max_orders[num_primes] = { 13, 13, 13, 13, 13 };

    for (ulong nb_prime = 0; nb_prime < num_primes; nb_prime++)
    //for (ulong nb_prime = 0; nb_prime < 1; nb_prime++)
    {
        // prime, modulus
        nmod_t mod;
        ulong p = primes[nb_prime];
        nmod_init(&mod, p);

        // find root of unity of specified maximum order
        mp_limb_t w0 = nmod_pow_ui(n_primitive_root_prime(p), (p - 1) >> max_orders[nb_prime], mod);

        printf("prime %ld, orders: ", nb_prime);

        for (ulong order = 3; order <= max_orders[nb_prime]; order++)
        {
            const ulong len = (1<<order);

            // root of unity of order 2**order
            mp_limb_t w = nmod_pow_ui(w0, 1UL<<(max_orders[nb_prime]-order), mod);

            // build FFT tables
            nmod_integer_fft_t F;
            nmod_integer_fft_init_set(F, w, order, mod);
            nmod_integer_fft_t Fpre;
            nmod_integer_fft_init_set_pre(Fpre, w, order, mod);

            // choose random poly
            nmod_poly_t pol;
            nmod_poly_init(pol, mod.n);
            nmod_poly_rand(pol, state, len);

            // naive evals by Horner
            mp_ptr evals = _nmod_vec_init(len);
            for (ulong k = 0; k < len/2; k++)
            {
                mp_limb_t point = F->tab_w[order-2][k];
                evals[k] = nmod_poly_evaluate_nmod(pol, point);
                evals[k+len/2] = nmod_poly_evaluate_nmod(pol, nmod_neg(point, F->mod));
            }
            // put in bit reversed order
            mp_ptr br_ind = flint_malloc(len * sizeof(mp_limb_t));
            brc_indices(br_ind, order);
            mp_ptr evals_br = _nmod_vec_init(len);
            for (ulong k = 0; k < len; k++)
                evals_br[k] = evals[br_ind[k]];

            // FFT evals, inplace radix 2 recursive
            nmod_poly_t pol2;
            nmod_poly_init(pol2, mod.n);
            nmod_poly_set(pol2, pol);
            nmod_poly_t pol3;
            nmod_poly_init(pol3, mod.n);
            nmod_poly_set(pol3, pol);
            nmod_poly_t pol4;
            nmod_poly_init(pol4, mod.n);
            nmod_poly_set(pol4, pol);
            nmod_poly_t pol5;
            nmod_poly_init(pol5, mod.n);
            nmod_poly_set(pol5, pol);
            nmod_poly_t pol6;
            nmod_poly_init(pol6, mod.n);
            nmod_poly_set(pol6, pol);
            nmod_poly_t pol7;
            nmod_poly_init(pol7, mod.n);
            nmod_poly_set(pol7, pol);
            nmod_poly_t pol8;
            nmod_poly_init(pol8, mod.n);
            nmod_poly_set(pol8, pol);

            _nmod_poly_dif_inplace_radix2_rec(pol->coeffs, len, order, F);
            _nmod_poly_dif_inplace_radix2_rec_v2(pol2->coeffs, len, order, F);
            _nmod_poly_dif_inplace_radix2_rec_v3(pol3->coeffs, len, order, F);
            _nmod_poly_dif_inplace_radix2_rec_v4(pol4->coeffs, len, order, Fpre);
            _nmod_poly_dif_inplace_radix2_iter(pol5->coeffs, len, order, F);
            _nmod_poly_dif_inplace_radix2_iter_v2(pol6->coeffs, len, order, Fpre);
            _nmod_poly_dif_inplace_radix4_rec(pol7->coeffs, len, order, F);
            _nmod_poly_dif_inplace_radix4_iter(pol8->coeffs, len, order, F);

            if (! _nmod_vec_equal(evals_br, pol->coeffs, len))
            {
                printf("\n\nERROR! in _nmod_poly_dif_inplace_radix2_rec\n\n");
                return;
            }
            else if (! _nmod_vec_equal(evals_br, pol2->coeffs, len))
            {
                printf("\n\nERROR! in _nmod_poly_dif_inplace_radix2_rec_v2\n\n");
                return;
            }
            else if (! _nmod_vec_equal(evals_br, pol3->coeffs, len))
            {
                printf("\n\nERROR! in _nmod_poly_dif_inplace_radix2_rec_v3\n\n");
                return;
            }
            else if (! _nmod_vec_equal(evals_br, pol4->coeffs, len))
            {
                printf("\n\nERROR! in _nmod_poly_dif_inplace_radix2_rec_v4\n\n");
                return;
            }
            else if (! _nmod_vec_equal(evals_br, pol5->coeffs, len))
            {
                printf("\n\nERROR! in _nmod_poly_dif_inplace_radix2_iter\n\n");
                return;
            }
            else if (! _nmod_vec_equal(evals_br, pol6->coeffs, len))
            {
                printf("\n\nERROR! in _nmod_poly_dif_inplace_radix2_iter_v2\n\n");
                return;
            }
            else if (! _nmod_vec_equal(evals_br, pol7->coeffs, len))
            {
                printf("\n\nERROR! in _nmod_poly_dif_inplace_radix4_rec\n\n");
                _nmod_vec_print(pol7->coeffs, len, mod);
                _nmod_vec_print(evals_br, len, mod);
                return;
            }
            else if (! _nmod_vec_equal(evals_br, pol8->coeffs, len))
            {
                printf("\n\nERROR! in _nmod_poly_dif_inplace_radix4_iter\n\n");
                _nmod_vec_print(pol8->coeffs, len, mod);
                _nmod_vec_print(evals_br, len, mod);
                return;
            }
            else
                printf("%ld ", order);

            nmod_poly_clear(pol);
            nmod_poly_clear(pol2);
            nmod_poly_clear(pol3);
            nmod_poly_clear(pol4);
            nmod_poly_clear(pol5);
            nmod_poly_clear(pol6);
            nmod_poly_clear(pol7);
            nmod_poly_clear(pol8);
            _nmod_vec_clear(evals);
            _nmod_vec_clear(evals_br);
            flint_free(br_ind);
            nmod_integer_fft_clear(F);
            nmod_integer_fft_clear_pre(Fpre);
        }
        printf("\n");
    }

    flint_randclear(state);
}


/*------------------------------------------------------------*/
/* main just calls time_init_set()                            */
/*------------------------------------------------------------*/
int main()
{
    setbuf(stdout, NULL);
    test_fft_eval();
    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
