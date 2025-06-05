#include <flint/nmod_vec.h>
#include <flint/ulong_extras.h>
#include <time.h>
#include <flint/nmod.h>
#include "nmod_poly_extra.h"
#include "nmod_poly_fft.h"

#define num_primes 5

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

// testing that all elements of "vec" are less than "bound"
int nmod_vec_range(nn_srcptr vec, ulong len, ulong bound)
{
    for (ulong k = 0; k < len; k++)
        if (vec[k] >= bound)
            return 0;

    return 1;
}

/*------------------------------------------------------------------*/
/* computes and check FFT eval for several bit lengths and depths   */
/*------------------------------------------------------------------*/
void test_fft_eval()
{
    flint_rand_t state;
    flint_rand_init(state);
    flint_rand_set_seed(state, time(NULL), time(NULL)+57);

    printf("- depth is log(fft length)\n");
    printf("- testing fft-eval for several bit lengths and depths\n");

    ulong primes[num_primes] = {
        786433,              // 20 bits, 1 + 2**18 * 3
        2013265921,          // 31 bits, 1 + 2**27 * 3 * 5
        2748779069441,       // 42 bits, 1 + 2**39 * 5
        1108307720798209,    // 50 bits, 1 + 2**44 * 3**2 * 7
        1139410705724735489, // 60 bits, 1 + 2**52 * 11 * 23
    };
    ulong max_depths[num_primes] = { 13, 13, 13, 13, 13 };
    //ulong max_depths[num_primes] = { 4, 4, 4, 4, 4 };

    for (ulong nb_prime = 0; nb_prime < num_primes; nb_prime++)
    {
        // prime, modulus
        nmod_t mod;
        ulong p = primes[nb_prime];
        nmod_init(&mod, p);

        // find root of unity of specified maximum depth
        const ulong prt = n_primitive_root_prime(p);
        ulong w0 = nmod_pow_ui(prt, (p - 1) >> max_depths[nb_prime], mod);

        printf("prime %ld, depths: ", nb_prime);

        for (ulong depth = 3; depth <= max_depths[nb_prime]; depth++)
        {
            const ulong len = (1UL<<depth);

            // build FFT tables
            n_fft_old_ctx_t F_old;
            ulong w = nmod_pow_ui(w0, 1UL<<(max_depths[nb_prime]-depth), mod);
            n_fft_old_ctx_init_set(F_old, w, depth, p);

            n_fft_ctx_t F;
            n_fft_ctx_init2(F, depth, p);
//            printf("prime: %lu\n", p);
//            printf("tab_w:\n");
//            for (ulong kk = 0; kk < 8; kk++)
//                printf("%lu, ", F->tab_w[kk]);
//            printf("\n");
//            printf("tab_w2:\n");
//            for (ulong kk = 0; kk < 8; kk++)
//                printf("%lu, ", F->tab_w2[kk]);
//            printf("\n");


            // choose random poly
            nmod_poly_t pol;
            nmod_poly_init(pol, mod.n);
            nmod_poly_rand(pol, state, len);

            // naive evals by Horner, in bit reversed order
            nn_ptr evals_br = _nmod_vec_init(len);
            for (ulong k = 0; k < len/2; k++)
            {
                ulong point = F->tab_w[2*k];
                evals_br[2*k] = nmod_poly_evaluate_nmod(pol, point);
                evals_br[2*k+1] = nmod_poly_evaluate_nmod(pol, nmod_neg(point, mod));
            }

            {  // dif red_rec2 lazy
                ulong * coeffs = _nmod_vec_init(len);
                _nmod_vec_set(coeffs, pol->coeffs, len);

                _n_fft_red_rec2_lazy(coeffs, len, depth, F);

                if (! nmod_vec_red_equal(evals_br, coeffs, len, mod)
                        || !nmod_vec_range(coeffs, len, 4*mod.n))
                {
                    printf("\n\nERROR! in _n_fft_red_rec2_lazy\n\n");
                    if (len < 33)
                    {
                        _nmod_vec_print(coeffs, len, mod);
                        _nmod_vec_print(evals_br, len, mod);
                    }
                    return;
                }

                _nmod_vec_clear(coeffs);
            }

            {  // dif red_rec4 lazy
                ulong * coeffs = _nmod_vec_init(len);
                _nmod_vec_set(coeffs, pol->coeffs, len);

                _n_fft_red_rec4_lazy(coeffs, len, depth, F);

                if (! nmod_vec_red_equal(evals_br, coeffs, len, mod)
                        || !nmod_vec_range(coeffs, len, 4*mod.n))
                {
                    printf("\n\nERROR! in _n_fft_red_rec2_lazy\n\n");
                    if (len < 33)
                    {
                        _nmod_vec_print(coeffs, len, mod);
                        _nmod_vec_print(evals_br, len, mod);
                    }
                    return;
                }

                _nmod_vec_clear(coeffs);
            }

            {  // red iter2 lazy
                ulong * coeffs = _nmod_vec_init(len);
                _nmod_vec_set(coeffs, pol->coeffs, len);

                _n_fft_red_iter2_lazy(coeffs, len, depth, F);

                if (! nmod_vec_red_equal(evals_br, coeffs, len, mod)
                        || !nmod_vec_range(coeffs, len, 8*mod.n))
                {
                    printf("\n\nERROR! in _n_fft_red_iter2_lazy\n\n");
                    if (len < 33)
                    {
                        _nmod_vec_print(coeffs, len, mod);
                        _nmod_vec_print(evals_br, len, mod);
                    }
                    return;
                }

                _nmod_vec_clear(coeffs);
            }

            {  // dif rec2 lazy
                ulong * coeffs = _nmod_vec_init(len);
                _nmod_vec_set(coeffs, pol->coeffs, len);

                _n_fft_old_dif_rec2_lazy(coeffs, len, depth, F_old);

                if (! nmod_vec_red_equal(evals_br, coeffs, len, mod)
                        || !nmod_vec_range(coeffs, len, 4*mod.n))
                {
                    printf("\n\nERROR! in _n_fft_old_dif_rec2_lazy\n\n");
                    if (len < 33)
                    {
                        _nmod_vec_print(coeffs, len, mod);
                        _nmod_vec_print(evals_br, len, mod);
                    }
                    return;
                }

                _nmod_vec_clear(coeffs);
            }

            {  // dif iter2 lazy
                ulong * coeffs = _nmod_vec_init(len);
                _nmod_vec_set(coeffs, pol->coeffs, len);

                _n_fft_old_dif_iter2_lazy(coeffs, len, depth, F_old);

                if (! nmod_vec_red_equal(evals_br, coeffs, len, mod)
                        || !nmod_vec_range(coeffs, len, 4*mod.n))
                {
                    printf("\n\nERROR! in _n_fft_old_dif_iter2_lazy\n\n");
                    if (len < 33)
                    {
                        _nmod_vec_print(coeffs, len, mod);
                        _nmod_vec_print(evals_br, len, mod);
                    }
                    return;
                }

                _nmod_vec_clear(coeffs);
            }

            {  // dif rec4 lazy
                ulong * coeffs = _nmod_vec_init(len);
                _nmod_vec_set(coeffs, pol->coeffs, len);

                _n_fft_old_dif_rec4_lazy(coeffs, len, depth, F_old);

                if (! nmod_vec_red_equal(evals_br, coeffs, len, mod)
                        || !nmod_vec_range(coeffs, len, 8*mod.n))
                {
                    printf("\n\nERROR! in _n_fft_old_dif_rec4_lazy\n\n");
                    if (len < 33)
                    {
                        _nmod_vec_print(coeffs, len, mod);
                        _nmod_vec_print(evals_br, len, mod);
                    }
                    return;
                }

                _nmod_vec_clear(coeffs);
            }

            {  // dif rec8 lazy
                ulong * coeffs = _nmod_vec_init(len);
                _nmod_vec_set(coeffs, pol->coeffs, len);

                _n_fft_old_dif_rec8_lazy(coeffs, len, depth, F_old);

                if (! nmod_vec_red_equal(evals_br, coeffs, len, mod)
                        || !nmod_vec_range(coeffs, len, 8*mod.n))
                {
                    printf("\n\nERROR! in _n_fft_old_dif_rec8_lazy\n\n");
                    if (len < 33)
                    {
                        _nmod_vec_print(coeffs, len, mod);
                        _nmod_vec_print(evals_br, len, mod);
                    }
                    return;
                }

                _nmod_vec_clear(coeffs);
            }

            printf("%ld ", depth);

            nmod_poly_clear(pol);
            _nmod_vec_clear(evals_br);
            n_fft_ctx_clear(F);
            n_fft_old_ctx_clear(F_old);
        }
        printf("\n");
    }

    flint_rand_clear(state);
}

/*------------------------------------------------------------*/
/* main just calls test_fft_eval()                            */
/*------------------------------------------------------------*/
int main()
{
    setbuf(stdout, NULL);
    test_fft_eval();
    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
