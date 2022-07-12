#include <assert.h>
#include <flint/flint.h>

#include "nmod_extra.h"
#include "nmod_poly_extra.h"

/*------------------------------------------------------------*/
/* TFT evaluation_t, then interpolation_t                     */
/* check if matches input vector                              */
/*------------------------------------------------------------*/
void check()
{
#ifdef HAS_AVX2
    ulong order_w0, order_max, nmin, nmax, n, i;
    flint_rand_t state;
    mp_limb_t p, w0, w;
    nmod_t mod;
    nmod_32_fft_t F;
    mp_ptr val, val2, val3;
    
    flint_randinit(state);

    p = 537133057;
    nmod_init(&mod, p);
    w0 = 144173337;
    order_w0 = 16; 
    order_max = 10;
    
    w = nmod_pow_ui(w0, 1L<<(order_w0-order_max), mod);
    nmod_32_fft_init_set(F, w, order_max, mod);

    nmin = 1;
    nmax = 1000;

    for (n = nmin; n < nmax+1; n++)
    {
        val = _nmod_vec_init(n);
        val2 = _nmod_vec_init(n);
        val3 = _nmod_vec_init(n);

        for (i = 0; i < n; i++)
        {
            val[i] = n_randtest(state) % p;
        }
        nmod_avx2_32_tft_evaluate_t(val2, val, F, n);
        nmod_avx2_32_tft_interpolate_t(val3, val2, F, n);

        for (i = 0; i < n; i++)
        {
            assert(val[i] == val3[i]);
        }
        _nmod_vec_clear(val);
        _nmod_vec_clear(val2);
        _nmod_vec_clear(val3);
    }

        
    nmod_32_fft_clear(F);
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
