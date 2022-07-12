#include <assert.h>
#include <flint/flint.h>

#include "nmod_extra.h"
#include "nmod_poly_extra.h"

/*------------------------------------------------------------*/
/* FFT evaluation, compared to naive evaluation               */
/*------------------------------------------------------------*/
void check()
{
#ifdef HAS_AVX2
    ulong order_w0, order_max, nmin, nmax;
    flint_rand_t state;
    mp_limb_t p, w0, w;
    nmod_t mod;
    nmod_32_fft_t F;
    mp_ptr val;
    nmod_poly_t P;
    
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
    
    for (long n = nmin; n < nmax+1; n++)
    {
        slong i, j;

        len = find_length(n);

        
        nmod_poly_init2(P, p, n);
        for (i = 0; i < n; i++)
        {
            nmod_poly_set_coeff_ui(P, i, n_randtest(state) % p);
        }
        
        val = _nmod_vec_init(n);
                
        nmod_avx2_32_tft_evaluate(val, P, F, n);

        i = 0;
        shift = 0;
        do
        {
            slong e, j;
            mp_limb_t rho, g;
            
            e = exponents[i];

            rho = nmod_pow_ui(w, 1L << (order_max-e-1), F->mod);
            g = nmod_pow_ui(w, 1L << (order_max-e), F->mod);

            for (j = 0; j < (1L << e); j++)
            {
                assert (val[j + shift] == nmod_poly_evaluate_nmod(P, nmod_mul(rho, nmod_pow_ui(g, bit_reverse(j, e), F->mod), F->mod)));
            }
            shift += (1L << e);
            i++;
        } while (exponents[i] != -1);

        nmod_poly_clear(P);
        _nmod_vec_clear(val);
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
