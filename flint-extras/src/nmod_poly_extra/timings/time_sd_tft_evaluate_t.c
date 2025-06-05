#include <assert.h>
#include <flint/flint.h>
#include <flint/nmod_mat.h>

#include "nmod_extra.h"
#include "nmod_poly_extra.h"

/*------------------------------------------------------------*/
/* computes TFTs and transpose TFTs                           */
/* uses fft_small-based TFT and transpose TFT                 */
/*------------------------------------------------------------*/
void get_time()
{
    ulong i, nmin, nmax, nb_iter;
    flint_rand_t state;
    ulong p, w0, w;
    nmod_t mod;
    nmod_sd_fft_t F;
    sd_fft_ctx_t Q;
    sd_fft_lctx_t QL;
    nn_ptr val, val2;
    nmod_poly_t P;
    clock_t tt;
    double t;

        
    flint_rand_init(state);

    p = 1108307720798209;
    sd_fft_ctx_init_prime(Q, p);
    sd_fft_lctx_init(QL, Q, 16);
    
    nmod_init(&mod, p);
    w0 = nmod_pow_ui(n_primitive_root_prime(p), (p - 1) >> 16, mod);
    w = nmod_pow_ui(w0, 1L<<(16-16), mod);
    nmod_sd_fft_init_set(F, w, 16, mod);

    nmin = 120;
    nmax = 2000;

    for (long n = nmin; n < nmax+1; n++)
    {
        nmod_poly_init2(P, p, n);
        for (i = 0; i < n; i++)
            nmod_poly_set_coeff_ui(P, i, n_randtest(state) % p);

        val = _nmod_vec_init(n);
        val2 = _nmod_vec_init(n);

        printf("%lu ", n);
        
        t = 0.0;
        nb_iter = 0;
        while (t < 0.5)
        {
            tt = clock();
            nmod_sd_tft_evaluate(val, P, QL, F, n);
            t += (double)(clock()-tt) / CLOCKS_PER_SEC;
            ++nb_iter;
        }
        t = 1000 * t;
        t /= nb_iter;
        printf("%f ", t);

        t = 0.0;
        nb_iter = 0;
        while (t < 0.5)
        {
            tt = clock();
            nmod_sd_tft_evaluate_t(val2, val, QL, F, n);
            t += (double)(clock()-tt) / CLOCKS_PER_SEC;
            ++nb_iter;
        }
        t = 1000 * t;
        t /= nb_iter;
        printf("%f ", t);

        printf("\n");
        
        nmod_poly_clear(P);
        _nmod_vec_clear(val);
        _nmod_vec_clear(val2);
    }

    sd_fft_ctx_clear(Q);
    nmod_sd_fft_clear(F);
    flint_rand_clear(state);
}



/*------------------------------------------------------------*/
/* main just calls get_time()                                 */
/*------------------------------------------------------------*/
int main(int argc, char **argv){
    get_time();
    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
