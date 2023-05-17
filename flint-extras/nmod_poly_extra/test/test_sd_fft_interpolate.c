#include <assert.h>
#include <flint/flint.h>

#include "nmod_extra.h"
#include "nmod_poly_extra.h"

/*------------------------------------------------------------*/
/* FFT evaluation and interpolation                           */
/*------------------------------------------------------------*/
void check()
{
    ulong order, order_min, order_max, N;
    slong i;
    flint_rand_t state;
    mp_limb_t p;
    nmod_t mod;
    mp_ptr val;
    sd_fft_ctx_t Q;
    nmod_poly_t P, P2;
    
    flint_randinit(state);
    
    p = 1108307720798209;
    sd_fft_ctx_init_prime(Q, p);
    nmod_init(&mod, p);

    order_min = 0;
    order_max = 10;
    
    for (order = order_min; order <= order_max; order++)
    {
        sd_fft_lctx_t QL;
        
        N = 1L << order;
        nmod_poly_init2(P, p, N);
        nmod_poly_init(P2, p);
        val = _nmod_vec_init(N);
        
        sd_fft_lctx_init(QL, Q, order);
        
        for (i = 0; i < N; i++)
        {
            mp_limb_t t = rand() % p;
            t = (t * t) % p;
            t = (t * t) % p;
            t = (t * t) % p;
            nmod_poly_set_coeff_ui(P, i, t);
        }

        nmod_sd_fft_evaluate(val, P, QL, order);
        nmod_sd_fft_interpolate(P2, val, QL, order);

        for (i = 0; i < N; i++)
        {
            /* printf("%lu %lu\n", nmod_mul(P->coeffs[i], 1L << order, mod), P2->coeffs[i]); */
            assert(nmod_mul(P->coeffs[i], 1L << order, mod) == P2->coeffs[i]);
        }
        /* printf("%lu\n", order); */
        
        _nmod_vec_clear(val);
        nmod_poly_clear(P);
        nmod_poly_clear(P2);
    }
    
    sd_fft_ctx_clear(Q);
    flint_randclear(state);
}



/*------------------------------------------------------------*/
/* main just calls check()                                    */
/*------------------------------------------------------------*/
int main(int argc, char **argv){
    check();
    return 0;
}
