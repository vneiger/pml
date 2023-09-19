#include <assert.h>
#include <flint/flint.h>
#include <flint/nmod_mat.h>

#include "nmod_extra.h"
#include "nmod_poly_extra.h"

/*------------------------------------------------------------*/
/* computes TFTs                                              */
/* uses fft_small FFT, my old TFT, the fft_small-based one    */
/*------------------------------------------------------------*/
void get_time()
{
    ulong i, nmin, nmax, nb_iter;
    flint_rand_t state;
    mp_limb_t p, w0, w;
    nmod_t mod;
    nmod_double_fft_t Fd;
    nmod_sd_fft_t F;
    sd_fft_ctx_t Q;
    sd_fft_lctx_t QL;
    mp_ptr val;
    nmod_poly_t P;
    clock_t tt;
    double t;

        
    flint_randinit(state);
    p = 1108307720798209;
    nmod_init(&mod, p);
// sd initialization
    sd_fft_ctx_init_prime(Q, p);
    sd_fft_lctx_init(QL, Q, 16);
    w0 = nmod_pow_ui(n_primitive_root_prime(p), (p - 1) >> 16, mod);
    w = nmod_pow_ui(w0, 1L<<(16-16), mod);
    nmod_sd_fft_init_set(F, w, 16, mod);
// double initialization
    nmod_double_fft_init_set(Fd, w, 16, mod);
    
    nmin = 1;
    nmax = 2000;

    printf("n double_tft sd_tft\n");
    
    for (long n = nmin; n < nmax+1; n++)
    {
        nmod_poly_init2(P, p, n);
        for (i = 0; i < n; i++)
        {
            nmod_poly_set_coeff_ui(P, i, n_randtest(state) % p);
        }
        val = _nmod_vec_init(1L << n_clog2(n));

        printf("%lu ", n);
        t = 0.0;
        nb_iter = 0;
        while (t < 0.5)
        {
            tt = clock();
            nmod_double_tft_evaluate(val, P, Fd, n);
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
            nmod_sd_tft_evaluate(val, P, QL, F, n);
            t += (double)(clock()-tt) / CLOCKS_PER_SEC;
            ++nb_iter;
        }
        t = 1000 * t;
        t /= nb_iter;
        printf("%f ", t);

        printf("\n");
        
        nmod_poly_clear(P);
        _nmod_vec_clear(val);
            
    }

    sd_fft_ctx_clear(Q);
    nmod_double_fft_clear(Fd);
    nmod_sd_fft_clear(F);
    flint_randclear(state);
}



/*------------------------------------------------------------*/
/* main just calls get_time()                                 */
/*------------------------------------------------------------*/
int main(int argc, char **argv){
    get_time();
    return 0;
}
