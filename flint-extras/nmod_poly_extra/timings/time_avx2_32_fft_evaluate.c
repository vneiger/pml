#include "nmod_poly_extra.h"

/*------------------------------------------------------------*/
/* computes FFTs up to order 2^14                             */
/*------------------------------------------------------------*/
void get_time()
{
#ifdef HAS_AVX2
    ulong order, order_max, N;
    slong i;
    flint_rand_t state;
    mp_limb_t p, w0, w;
    nmod_t mod;
    nmod_32_fft_t F;
    mp_ptr val;
    nmod_poly_t P;
    
    flint_randinit(state);

    p = 7340033;
    nmod_init(&mod, p);
    w0 = 3308891;
    order_max = 14;
    w = nmod_pow_ui(w0, 1L<<(16-order_max), mod);
    nmod_32_fft_init_set(F, w, order_max, mod);
    
    for (order = 1; order <= order_max; order++)
    {
        double t;
        clock_t tt;
        long nb_iter;

        w = nmod_pow_ui(w0, 1L<<(16-order), mod);
        N = 1L << order;
        nmod_poly_init2(P, p, N);
        val = _nmod_vec_init(N);
        
        for (i = 0; i < N; i++)
        {
            nmod_poly_set_coeff_ui(P, i, n_randtest(state) % p);
        }
        
        t = 0.0;
        nb_iter = 0;
        while (t < 0.5)
        {
            tt = clock();
            nmod_avx2_32_fft_evaluate(val, P, F, order);
            t += (double)(clock()-tt) / CLOCKS_PER_SEC;
            ++nb_iter;
        }
        t = 1000 * t;
        t /= nb_iter;
        printf("%ld %f\n", N, t);

        _nmod_vec_clear(val);
        nmod_poly_clear(P);
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
