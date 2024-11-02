#include <flint/nmod_vec.h>
#include <flint/ulong_extras.h>
#include <time.h>
#include <flint/nmod.h>
#include "nmod_poly_extra.h"
#include "nmod_poly_integer_fft.h"

#define num_primes 5

/***********************
*  bit reversed copy  *
***********************/

//long RevInc(long a, long k)
//{
//    long j, m;
//
//    j = k;
//    m = 1L << (k-1);
//
//    while (j && (m & a)) {
//        a ^= m;
//        m >>= 1;
//        j--;
//    }
//    if (j) a ^= m;
//    return a;
//}
//
//// indices initialized with length >= k
//void brc_indices(ulong * indices, long k)
//{
//    const long n = (1L << k);
//    for (long i = 0, j = 0; i < n; i++, j = RevInc(j, k))
//        indices[i] = j;
//}


// vector equality up to reduction mod
int nmod_vec_red_equal(nn_srcptr vec1, nn_srcptr vec2, ulong len, nmod_t mod)
{
    for (ulong k = 0; k < len; k++)
    {
        ulong v1;
        ulong v2;
        NMOD_RED(v1, vec1[k], mod);
        NMOD_RED(v2, vec2[k], mod);
        if (v1 != v2)
            return 0;
    }

    return 1;
}

int nmod_vec_range(nn_srcptr vec, ulong len, ulong bound)
{
    for (ulong k = 0; k < len; k++)
        if (vec[k] >= bound)
            return 0;

    return 1;
}

/*------------------------------------------------------------*/
/* computes init for FFT for several bit lengths and orders   */
/*------------------------------------------------------------*/
void test_fft_eval()
{
    flint_rand_t state;
    flint_rand_init(state);
    flint_rand_set_seed(state, time(NULL), time(NULL)+57);

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
    //ulong max_orders[num_primes] = { 4, 4, 4, 4, 4 };

    for (ulong nb_prime = 0; nb_prime < num_primes; nb_prime++)
    //for (ulong nb_prime = 0; nb_prime < 1; nb_prime++)
    {
        // prime, modulus
        nmod_t mod;
        ulong p = primes[nb_prime];
        nmod_init(&mod, p);

        // find root of unity of specified maximum order
        const ulong prt = n_primitive_root_prime(p);
        ulong w0 = nmod_pow_ui(prt, (p - 1) >> max_orders[nb_prime], mod);

        printf("prime %ld, orders: ", nb_prime);

        for (ulong order = 3; order <= max_orders[nb_prime]; order++)
        {
            const ulong len = (1UL<<order);

            // root of unity of order 2**order
            ulong w = nmod_pow_ui(w0, 1UL<<(max_orders[nb_prime]-order), mod);

            // build FFT tables
            nmod_integer_fft_t F;
            nmod_integer_fft_init_set(F, w, order, mod);
            nmod_integer_fft_t Fpre;
            nmod_integer_fft_init_set_pre(Fpre, w, order, mod);
            nmod_integer_fft_t Fred;
            nmod_integer_fft_init_set_red(Fred, w, order, mod);
            nmod_integer_fft_t Fredpre;
            nmod_integer_fft_init_set_red_pre(Fredpre, w, order, mod);

            // choose random poly
            nmod_poly_t pol;
            nmod_poly_init(pol, mod.n);
            nmod_poly_rand(pol, state, len);

            //// naive evals by Horner
            //nn_ptr evals = _nmod_vec_init(len);
            //for (ulong k = 0; k < len/2; k++)
            //{
            //    ulong point = F->tab_w[order-2][k];
            //    evals[k] = nmod_poly_evaluate_nmod(pol, point);
            //    evals[k+len/2] = nmod_poly_evaluate_nmod(pol, nmod_neg(point, F->mod));
            //}

            // naive evals by Horner, in bit reversed order
            nn_ptr evals_br = _nmod_vec_init(len);
            for (ulong k = 0; k < len/2; k++)
            {
                ulong point = Fred->tab_w[1][k];
                evals_br[2*k] = nmod_poly_evaluate_nmod(pol, point);
                evals_br[2*k+1] = nmod_poly_evaluate_nmod(pol, nmod_neg(point,Fred->mod));
            }

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
            nmod_poly_t pol9;
            nmod_poly_init(pol9, mod.n);
            nmod_poly_set(pol9, pol);
            nmod_poly_t pol10;
            nmod_poly_init(pol10, mod.n);
            nmod_poly_set(pol10, pol);
            nmod_poly_t pol11;
            nmod_poly_init(pol11, mod.n);
            nmod_poly_set(pol11, pol);
            nmod_poly_t pol12;
            nmod_poly_init(pol12, mod.n);
            nmod_poly_set(pol12, pol);
            nmod_poly_t pol13;
            nmod_poly_init(pol13, mod.n);
            nmod_poly_set(pol13, pol);
            nmod_poly_t pol14;
            nmod_poly_init(pol14, mod.n);
            nmod_poly_set(pol14, pol);

            _nmod_poly_dif_inplace_radix2_rec_prenorm(pol->coeffs, len, order, F);
            _nmod_poly_dif_inplace_radix2_rec_prenorm_unroll4(pol2->coeffs, len, order, F);
            _nmod_poly_dif_inplace_radix2_rec_shoup(pol3->coeffs, len, order, Fpre);
            _nmod_poly_dif_inplace_radix2_rec_shoup_unroll4(pol4->coeffs, len, order, Fpre);
            _nmod_poly_dif_inplace_radix2_iter_prenorm(pol5->coeffs, len, order, F);
            _nmod_poly_dif_inplace_radix2_iter_shoup(pol6->coeffs, len, order, Fpre);
            _nmod_poly_dif_inplace_radix4_rec(pol7->coeffs, len, order, F);
            _nmod_poly_dif_inplace_radix4_iter(pol8->coeffs, len, order, F);
            _nmod_poly_red_inplace_radix2_rec_prenorm(pol9->coeffs, len, order, 0, Fred);
            _nmod_poly_red_inplace_radix2_rec_shoup(pol10->coeffs, len, order, 0, Fredpre);

            _nmod_poly_dif_inplace_radix2_rec_shoup_lazy(pol11->coeffs, len, order, Fpre);
            _nmod_poly_dif_inplace_radix2_iter_shoup_lazy(pol12->coeffs, len, order, Fpre);
            _nmod_poly_red_inplace_radix2_rec_shoup_lazy(pol13->coeffs, len, order, 0, Fredpre);
            _nmod_poly_dif_inplace_radix4_rec_shoup_lazy(pol14->coeffs, len, order, Fpre);

            if (! _nmod_vec_equal(evals_br, pol->coeffs, len))
            {
                printf("\n\nERROR! in _nmod_poly_dif_inplace_radix2_rec_prenorm\n\n");
                if (len < 33)
                {
                    _nmod_vec_print(pol->coeffs, len, mod);
                    _nmod_vec_print(evals_br, len, mod);
                }
                return;
            }
            else if (! _nmod_vec_equal(evals_br, pol2->coeffs, len))
            {
                printf("\n\nERROR! in _nmod_poly_dif_inplace_radix2_rec_prenorm_unroll4\n\n");
                if (len < 33)
                {
                    _nmod_vec_print(pol2->coeffs, len, mod);
                    _nmod_vec_print(evals_br, len, mod);
                }
                return;
            }
            else if (! _nmod_vec_equal(evals_br, pol3->coeffs, len))
            {
                printf("\n\nERROR! in _nmod_poly_dif_inplace_radix2_rec_shoup\n\n");
                if (len < 33)
                {
                    _nmod_vec_print(pol3->coeffs, len, mod);
                    _nmod_vec_print(evals_br, len, mod);
                }
                return;
            }
            else if (! _nmod_vec_equal(evals_br, pol4->coeffs, len))
            {
                printf("\n\nERROR! in _nmod_poly_dif_inplace_radix2_rec_shoup_unroll4\n\n");
                if (len < 33)
                {
                    _nmod_vec_print(pol4->coeffs, len, mod);
                    _nmod_vec_print(evals_br, len, mod);
                }
                return;
            }
            else if (! _nmod_vec_equal(evals_br, pol5->coeffs, len))
            {
                printf("\n\nERROR! in _nmod_poly_dif_inplace_radix2_iter_prenorm\n\n");
                if (len < 33)
                {
                    _nmod_vec_print(pol5->coeffs, len, mod);
                    _nmod_vec_print(evals_br, len, mod);
                }
                return;
            }
            else if (! _nmod_vec_equal(evals_br, pol6->coeffs, len))
            {
                printf("\n\nERROR! in _nmod_poly_dif_inplace_radix2_iter_shoup\n\n");
                if (len < 33)
                {
                    _nmod_vec_print(pol6->coeffs, len, mod);
                    _nmod_vec_print(evals_br, len, mod);
                }
                return;
            }
            else if (! _nmod_vec_equal(evals_br, pol7->coeffs, len))
            {
                printf("\n\nERROR! in _nmod_poly_dif_inplace_radix4_rec\n\n");
                if (len < 33)
                {
                    _nmod_vec_print(pol7->coeffs, len, mod);
                    _nmod_vec_print(evals_br, len, mod);
                }
                return;
            }
            else if (! _nmod_vec_equal(evals_br, pol8->coeffs, len))
            {
                printf("\n\nERROR! in _nmod_poly_dif_inplace_radix4_iter\n\n");
                if (len < 33)
                {
                    _nmod_vec_print(pol8->coeffs, len, mod);
                    _nmod_vec_print(evals_br, len, mod);
                }
                return;
            }
            else if (! _nmod_vec_equal(evals_br, pol9->coeffs, len))
            {
                printf("\n\nERROR! in _nmod_poly_red_inplace_radix2_rec_prenorm\n\n");
                if (len < 33)
                {
                    _nmod_vec_print(pol9->coeffs, len, mod);
                    _nmod_vec_print(evals_br, len, mod);
                }
                return;
            }
            else if (! _nmod_vec_equal(evals_br, pol10->coeffs, len))
            {
                printf("\n\nERROR! in _nmod_poly_red_inplace_radix2_rec_shoup\n\n");
                if (len < 33)
                {
                    _nmod_vec_print(pol10->coeffs, len, mod);
                    _nmod_vec_print(evals_br, len, mod);
                }
                return;
            }
            else if (! nmod_vec_red_equal(evals_br, pol11->coeffs, len, mod)
                     || !nmod_vec_range(pol11->coeffs, len, 4*mod.n))
            {
                printf("\n\nERROR! in _nmod_poly_dif_inplace_radix2_rec_shoup_lazy\n\n");
                if (len < 33)
                {
                    _nmod_vec_print(pol11->coeffs, len, mod);
                    _nmod_vec_print(evals_br, len, mod);
                }
                return;
            }
            else if (! _nmod_vec_equal(pol11->coeffs, pol12->coeffs, len)
                     || !nmod_vec_range(pol12->coeffs, len, 4*mod.n))
            {
                printf("\n\nERROR! in _nmod_poly_dif_inplace_radix2_iter_shoup_lazy\n\n");
                if (len < 33)
                {
                    _nmod_vec_print(pol11->coeffs, len, mod);
                    _nmod_vec_print(pol12->coeffs, len, mod);
                }
                return;
            }
            else if (! nmod_vec_red_equal(evals_br, pol13->coeffs, len, mod)
                     || !nmod_vec_range(pol13->coeffs, len, 8*mod.n))
            {
                printf("\n\nERROR! in _nmod_poly_red_inplace_radix2_rec_shoup_lazy\n\n");
                if (len < 33)
                {
                    _nmod_vec_print(pol13->coeffs, len, mod);
                    _nmod_vec_print(evals_br, len, mod);
                }
                return;
            }
            else if (! nmod_vec_red_equal(evals_br, pol14->coeffs, len, mod)
                     || !nmod_vec_range(pol14->coeffs, len, 8*mod.n))
            {
                printf("\n\nERROR! in _nmod_poly_dif_inplace_radix4_rec_shoup_lazy\n\n");
                if (len < 33)
                {
                    _nmod_vec_print(pol14->coeffs, len, mod);
                    _nmod_vec_print(evals_br, len, mod);
                }
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
            nmod_poly_clear(pol9);
            nmod_poly_clear(pol10);
            nmod_poly_clear(pol11);
            nmod_poly_clear(pol12);
            nmod_poly_clear(pol13);
            nmod_poly_clear(pol14);
            _nmod_vec_clear(evals_br);
            nmod_integer_fft_clear(F);
            nmod_integer_fft_clear_pre(Fpre);
            nmod_integer_fft_clear_red(Fred);
            nmod_integer_fft_clear_red_pre(Fredpre);
        }
        printf("\n");
    }

    flint_rand_clear(state);
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
