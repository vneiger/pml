#include "nmod_poly_extra.h"

/*------------------------------------------------------------*/
/* computes FFTs up to order 2^14                             */
/*------------------------------------------------------------*/
void get_time()
{
    ulong order, order_min, order_max, N;
    slong i;
    flint_rand_t state;
    mp_limb_t p, w0, w;
    sd_fft_ctx_t Q;
    sd_fft_lctx_t QL;
    nmod_t mod;
    nmod_double_fft_t F;
    mp_ptr val;
    nmod_poly_t P;
    
    flint_randinit(state);

    p = 1108307720798209;
    nmod_init(&mod, p);
    sd_fft_ctx_init_prime(Q, p);
    sd_fft_lctx_init(QL, Q, 16);

    order_min = 6;
    order_max = 11;

    w0 = nmod_pow_ui(n_primitive_root_prime(p), (p - 1) >> 16, mod);
    w = nmod_pow_ui(w0, 1L<<(16-order_max), mod);
    nmod_double_fft_init_set(F, w, order_max, mod);
    
    for (order = order_min; order < order_max; order++)
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
        
        printf("%ld ", N);

        t = 0.0;
        nb_iter = 0;
        while (t < 0.5)
        {
            tt = clock();
            nmod_double_fft_evaluate(val, P, F, order);
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
            nmod_double_tft_evaluate(val, P, F, N);
            t += (double)(clock()-tt) / CLOCKS_PER_SEC;
            ++nb_iter;
        }
        t = 1000 * t;
        t /= nb_iter;
        printf("%f | ", t);

        t = 0.0;
        nb_iter = 0;
        while (t < 0.5)
        {
            tt = clock();
            nmod_sd_fft_evaluate(val, P, QL, order);
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
            nmod_sd_tft_evaluate(val, P, QL, N);
            t += (double)(clock()-tt) / CLOCKS_PER_SEC;
            ++nb_iter;
        }
        t = 1000 * t;
        t /= nb_iter;
        printf("%f", t);

        printf("\n");
        _nmod_vec_clear(val);
        nmod_poly_clear(P);
    }
    
    sd_fft_ctx_clear(Q);
    nmod_double_fft_clear(F);
    flint_randclear(state);
}


/*------------------------------------------------------------*/
/* main just calls get_time()                                 */
/*------------------------------------------------------------*/
int main(int argc, char **argv){
    get_time();
    return 0;
}
