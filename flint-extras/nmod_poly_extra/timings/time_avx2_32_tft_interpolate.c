#include "nmod_poly_extra.h"

/*------------------------------------------------------------*/
/* computes FFTs up to order 2^14                             */
/*------------------------------------------------------------*/
void get_time()
{
#ifdef HAS_AVX2
    ulong order_w0, order_max, nmin, nmax;
    slong i;
    flint_rand_t state;
    mp_limb_t p, w0, w;
    nmod_t mod;
    nmod_32_fft_t F;
    mp_ptr val;
    nmod_poly_t P, P2;
    
    flint_randinit(state);

    p = 537133057;
    nmod_init(&mod, p);
    w0 = 144173337;
    order_w0 = 16; 
    order_max = 10;
    w = nmod_pow_ui(w0, 1L<<(order_w0-order_max), mod);
    nmod_32_fft_init_set(F, w, order_max, mod);

    nmin = 500;
    nmax = 500;
    
    for (long n = nmin; n < nmax+1; n++)
    {
        double t;
        clock_t tt;
        long nb_iter;
        
        nmod_poly_init2(P, p, n);
        nmod_poly_init(P2, p);
        for (i = 0; i < n; i++)
        {
            nmod_poly_set_coeff_ui(P, i, n_randtest(state) % p);
        }
        
        val = _nmod_vec_init(n);
        nmod_avx2_32_tft_evaluate(val, P, F, n);
                
        t = 0.0;
        nb_iter = 0;
        while (t < 0.5)
        {
            tt = clock();
            nmod_avx2_32_tft_interpolate(P2, val, F, n);
            t += (double)(clock()-tt) / CLOCKS_PER_SEC;
            ++nb_iter;
        }
        t = 1000 * t;
        t /= nb_iter;
        printf("%ld %f\n", n, t);

        _nmod_vec_clear(val);
        nmod_poly_clear(P);
        nmod_poly_clear(P2);
    }
    
    nmod_32_fft_clear(F);
    flint_randclear(state);
#endif
}

/*------------------------------------------------------------*/
/* main just calls get_time()                                 */
/*------------------------------------------------------------*/
int main(int argc, char **argv){
    get_time();
    return 0;
}
