#include <assert.h>
#include <flint/flint.h>

#include "nmod_extra.h"
#include "nmod_poly_extra.h"

/*------------------------------------------------------------*/
/* FFT evaluation/interpolation                               */
/*------------------------------------------------------------*/
void check()
{
#ifdef HAS_AVX2
    ulong order, order_max, N;
    slong i;
    flint_rand_t state;
    mp_limb_t p, w0, w;
    nmod_t mod;
    nmod_avx2_32_fft_t F;
    mp_ptr val;
    nmod_poly_t P, P2;
    
    flint_randinit(state);

    p = 7340033;
    nmod_init(&mod, p);
    w0 = 3308891;
    order_max = 10;
    w = nmod_pow_ui(w0, 1L<<(16-order_max), mod);
    nmod_avx2_32_fft_init_set(F, w, order_max, mod);
    //
    for (order = 1; order <= order_max; order++)
    {
        w = nmod_pow_ui(w0, 1L<<(16-order), mod);
        N = 1L << order;
        nmod_poly_init2(P, p, N);
        nmod_poly_init(P2, p);
        val = _nmod_vec_init(N);
        
        for (i = 0; i < N; i++)
        {
            nmod_poly_set_coeff_ui(P, i, n_randtest(state) % p);
        }
        nmod_avx2_32_fft_evaluate(val, P, F, order);
        nmod_avx2_32_fft_interpolate(P2, val, F, order);

        // check if interpolate(evaluate(P)) == P
        assert (nmod_poly_equal(P, P2));
        
        _nmod_vec_clear(val);
        nmod_poly_clear(P);
        nmod_poly_clear(P2);
    }
    
    nmod_avx2_32_fft_clear(F);
    flint_randclear(state);
#endif
}



/*------------------------------------------------------------*/
/* main just calls check()                                    */
/*------------------------------------------------------------*/
int main(int argc, char **argv){
    check();
    return 0;
}
