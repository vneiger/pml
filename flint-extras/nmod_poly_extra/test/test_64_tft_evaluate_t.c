#include <assert.h>
#include <flint/flint.h>
#include <flint/nmod_mat.h>

#include "nmod_extra.h"
#include "nmod_poly_extra.h"

/*------------------------------------------------------------*/
/* TFT evaluation_t                                           */
/* builds the matrix of TFT evaluation to compare results     */
/*------------------------------------------------------------*/
void check()
{
    ulong order_w0, order_max, nmin, nmax;
    flint_rand_t state;
    mp_limb_t p, w0, w;
    nmod_t mod;
    nmod_64_fft_t F;
    mp_ptr val, val2, val3;
    nmod_poly_t P;
    
    flint_randinit(state);

    p = 537133057;
    nmod_init(&mod, p);
    w0 = 144173337;
    order_w0 = 16; 
    order_max = 10;
    
    w = nmod_pow_ui(w0, 1L<<(order_w0-order_max), mod);
    nmod_64_fft_init_set(F, w, order_max, mod);

    nmin = 1;
    nmax = 1000;
    
    for (long n = nmin; n < nmax+1; n++)
    {
        slong i, j;
        nmod_mat_t M;
        val = _nmod_vec_init(n);
        val2 = _nmod_vec_init(n);
        val3 = _nmod_vec_init(n);
        nmod_poly_init2(P, p, n);
        nmod_mat_init(M, n, n, p);

        
        for (i = 0; i < n; i++)
        {
            nmod_poly_zero(P);
            nmod_poly_set_coeff_ui(P, i, 1);
            nmod_64_tft_evaluate(val, P, F, n);
            for (j = 0; j < n; j++)
            {
                nmod_mat_set_entry(M, j, i, val[j]);
            }
        }
        
        for (i = 0; i < n; i++)
        {
            val[i] = n_randtest(state) % p;
        }

        nmod_mat_nmod_vec_mul(val2, val, n, M);
        nmod_64_tft_evaluate_t(val3, val, F, n);

        for (i = 0; i < n; i++)
        {
            assert (val2[i]==val3[i]);
        }

        
        nmod_mat_clear(M);
        nmod_poly_clear(P);
        _nmod_vec_clear(val);
        _nmod_vec_clear(val2);
        _nmod_vec_clear(val3);
    }
        
    nmod_64_fft_clear(F);
    flint_randclear(state);
}



/*------------------------------------------------------------*/
/* main just calls check()                                    */
/*------------------------------------------------------------*/
int main(int argc, char **argv){
    check();
    return 0;
}
