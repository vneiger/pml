#include <assert.h>
#include <flint/flint.h>
#include <flint/nmod_mat.h>

#include "nmod_extra.h"
#include "nmod_poly_extra.h"

/*------------------------------------------------------------*/
/* TFT evaluation and transposed evaluation                   */
/*------------------------------------------------------------*/
void check()
{
    ulong N, i, j;
    nmod_mat_t M, Mt;
    flint_rand_t state;
    mp_limb_t w0, w, p;
    nmod_t mod;
    mp_ptr val, val2;
    sd_fft_ctx_t Q, Qt;
    nmod_sd_fft_t F;
    nmod_poly_t P;
    
    flint_randinit(state);
    
    p = 1108307720798209;
    sd_fft_ctx_init_prime(Q, p);
    sd_fft_ctx_init_inverse(Qt, Q);

    nmod_init(&mod, p);

    w0 = nmod_pow_ui(n_primitive_root_prime(p), (p - 1) >> 16, mod);
    w = nmod_pow_ui(w0, 1L<<(16-16), mod);
    nmod_sd_fft_init_set(F, w, 16, mod);

    for (N = 1; N < 400; N++)
    {
        sd_fft_lctx_t QL, QtL;

        nmod_poly_init(P, p);
        nmod_mat_init(M, N, N, p);
        nmod_mat_init(Mt, N, N, p);
        val = _nmod_vec_init(N);
        val2 = _nmod_vec_init(N);
        
        sd_fft_lctx_init(QL, Q, SD_FFT_CTX_INIT_DEPTH);
        sd_fft_lctx_init(QtL, Qt, SD_FFT_CTX_INIT_DEPTH);
        
        for (i = 0; i < N; i++)
        {
            nmod_poly_zero(P);
            nmod_poly_set_coeff_ui(P, i, 1);
            nmod_sd_tft_evaluate(val, P, QL, F, N);
            for (j = 0; j < N; j++)
            {
                nmod_mat_set_entry(M, j, i, val[j]);
            }
        }

        for (i = 0; i < N; i++)
        {
            for (j = 0; j < N; j++)
            {
                val2[j] = 0;
            }
            val2[i] = 1;
            nmod_sd_tft_evaluate_t(val, val2, QtL, F, N);
            for (j = 0; j < N; j++)
            {
                nmod_mat_set_entry(Mt, i, j, val[j]);
            }
        }

        if (!nmod_mat_equal(M, Mt))
        {
            printf("%lu \n", N);
            continue;
            _nmod_vec_clear(val);
            _nmod_vec_clear(val2);
            nmod_poly_clear(P);
            nmod_mat_clear(M);
            nmod_mat_clear(Mt);
        }

        _nmod_vec_clear(val);
        _nmod_vec_clear(val2);
        nmod_poly_clear(P);
        nmod_mat_clear(M);
        nmod_mat_clear(Mt);
    }
    
    sd_fft_ctx_clear(Q);
    sd_fft_ctx_clear(Qt);
    nmod_sd_fft_clear(F);
    
    flint_randclear(state);
}    




/*------------------------------------------------------------*/
/* main just calls check()                                    */
/*------------------------------------------------------------*/
int main(int argc, char **argv){
    check();
    return 0;
}
