#include <assert.h>
#include <flint/flint.h>

#include "nmod_extra.h"
#include "nmod_poly_extra.h"

/*------------------------------------------------------------*/
/* TFT evaluation, then interpolation                         */
/* checks if matches input poly                               */
/*------------------------------------------------------------*/
void check()
{
#ifdef HAS_INT128
    ulong order_w0, order_max, nmin, nmax, n, i;
    flint_rand_t state;
    mp_limb_t p, w0, w;
    nmod_t mod;
    nmod_64_fft_t F;
    mp_ptr val;
    nmod_poly_t P, P2;
    
    flint_randinit(state);

    p = 1073872897;
    nmod_init(&mod, p);
    w0 = 1056585332;
    order_w0 = 16; 
    order_max = 10;
    
    w = nmod_pow_ui(w0, 1L<<(order_w0-order_max), mod);
    nmod_64_fft_init_set(F, w, order_max, mod);

    nmin = 1;
    nmax = 1000;

    for (n = nmin; n < nmax+1; n++)
    {

        nmod_poly_init2(P, p, n);
        nmod_poly_init(P2, p);
        for (i = 0; i < n; i++)
        {
            nmod_poly_set_coeff_ui(P, i, n_randtest(state) % p);
        }
        
        val = _nmod_vec_init(n);
                
        nmod_64_tft_evaluate(val, P, F, n);
        nmod_64_tft_interpolate(P2, val, F, n);
        
        assert (nmod_poly_equal(P, P2));
        
        nmod_poly_clear(P2);
        nmod_poly_clear(P);
        _nmod_vec_clear(val);
    }
        
    nmod_64_fft_clear(F);
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
