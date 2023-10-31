#include <assert.h>
#include <flint/flint.h>
#include <flint/nmod_mat.h>

#include "nmod_extra.h"
#include "nmod_poly_extra.h"

/*------------------------------------------------------------*/
/* FFT interpolation and transposed interpolation             */
/*------------------------------------------------------------*/
void check()
{
    ulong order, order_min, order_max, N, i, j;
    nmod_mat_t M, Mt;
    flint_rand_t state;
    mp_limb_t p;
    nmod_t mod;
    mp_ptr val;
    sd_fft_ctx_t Q, Qt;
    nmod_poly_t P;
    
    flint_randinit(state);
    
    p = 1108307720798209;
    sd_fft_ctx_init_prime(Q, p);
    sd_fft_ctx_init_inverse(Qt, Q);

    nmod_init(&mod, p);

    order_min = 1;
    order_max = 10;
    
    for (order = order_min; order <= order_max; order++)
    {
        sd_fft_lctx_t QL, QtL;
    
        N = 1L << order;
        nmod_poly_init2(P, p, N);
        nmod_mat_init(M, N, N, p);
        nmod_mat_init(Mt, N, N, p);
        val = _nmod_vec_init(N);
        
        sd_fft_lctx_init(QL, Q, order);
        sd_fft_lctx_init(QtL, Qt, order);
        
        for (i = 0; i < N; i++)
        {
            for (j = 0; j < N; j++)
            {
                val[j] = 0;
            }
            val[i] = 1;
            nmod_sd_fft_interpolate(P, val, QL, order);
            for (j = 0; j < N; j++)
            {
                nmod_mat_set_entry(M, j, i, nmod_poly_get_coeff_ui(P, j));
            }
        }

        for (i = 0; i < N; i++)
        {
            for (j = 0; j < N; j++)
            {
                val[j] = 0;
            }
            val[i] = 1;
            nmod_sd_fft_interpolate_t(P, val, QL, order);
            for (j = 0; j < N; j++)
            {
                nmod_mat_set_entry(Mt, i, j, nmod_poly_get_coeff_ui(P, j));
            }
        }

        if (!nmod_mat_equal(M, Mt))
        {
            printf("%lu \n", order);
            nmod_mat_print(M);
            nmod_mat_print(Mt);
            continue;
        }

        _nmod_vec_clear(val);
        nmod_poly_clear(P);
        nmod_mat_clear(M);
        nmod_mat_clear(Mt);
    }
    
    sd_fft_ctx_clear(Q);
    sd_fft_ctx_clear(Qt);

    flint_randclear(state);
}    




/*------------------------------------------------------------*/
/* main just calls check()                                    */
/*------------------------------------------------------------*/
int main(int argc, char **argv){
    check();
    return 0;
}
