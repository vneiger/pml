#include <assert.h>
#include <flint/flint.h>

#include "nmod_extra.h"
#include "nmod_poly_extra.h"

/*------------------------------------------------------------*/
/* TFT evaluation and interpolation                           */
/*------------------------------------------------------------*/
void check()
{
    ulong nmin, nmax;
    flint_rand_t state;
    mp_limb_t p;
    nmod_t mod;
    sd_fft_ctx_t Q;
    mp_ptr val;
    nmod_poly_t P, P2;
    
    flint_randinit(state);

    p = 1108307720798209;
    sd_fft_ctx_init_prime(Q, p);
    nmod_init(&mod, p);

    nmin = 1;
    nmax = 2000;

    for (long n = nmin; n < nmax+1; n++)
    {
        sd_fft_lctx_t QL;
        ulong i;
        
        sd_fft_lctx_init(QL, Q, 16);

        nmod_poly_init2(P, p, n);
        nmod_poly_init(P2, p);
        for (i = 0; i < n; i++)
        {
            nmod_poly_set_coeff_ui(P, i, n_randtest(state) % p);
        }
        
        val = _nmod_vec_init(n);
        nmod_sd_tft_evaluate(val, P, QL, n);
        nmod_sd_tft_interpolate(P2, val, QL, n);

        assert(nmod_poly_equal(P, P2));
        
        nmod_poly_clear(P2);
        nmod_poly_clear(P);
        _nmod_vec_clear(val);
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
