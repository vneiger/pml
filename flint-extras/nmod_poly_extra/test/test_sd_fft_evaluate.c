#include <assert.h>
#include <flint/flint.h>
#include <flint/nmod_mat.h>

#include "nmod_extra.h"
#include "nmod_poly_extra.h"

/*------------------------------------------------------------*/
/* FFT evaluation, compared to naive evaluation               */
/*------------------------------------------------------------*/
void check()
{
    ulong order, order_min, order_max, N;
    slong i;
    flint_rand_t state;
    mp_limb_t p, w0, w, w2;
    nmod_t mod;
    mp_ptr val;
    sd_fft_ctx_t Q;
    nmod_poly_t P;
    
    flint_randinit(state);
    
    p = 1108307720798209;
    sd_fft_ctx_init_prime(Q, p);
    nmod_init(&mod, p);

    order_min = 0;
    order_max = 9;
    w0 = nmod_pow_ui(n_primitive_root_prime(p), (p - 1) >> 16, mod);
    
    for (order = order_min; order <= order_max; order++)
    {
        sd_fft_lctx_t QL;
        
        w = nmod_pow_ui(w0, 1L << (16-order), mod);
        N = 1L << order;
        nmod_poly_init2(P, p, N);
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

        printf("S:=[");
        for (i = 0; i < N; i++)
        {
            printf("%ld", val[i]);
            if (i < N-1)
            {
                printf(", ");
            }
        }
        printf("];");
        printf("\n");
        printf("U:=[");
        w2 = 1;
        for (i = 0; i < N; i++)
        {
            printf("%lu", nmod_poly_evaluate_nmod(P, w2));
            if (i < N-1)
            {
                printf(", ");
            }
            w2 = nmod_mul(w2, w, mod);
        }
        printf("];\n");

        printf("SequenceToSet(S) eq SequenceToSet(U);\n");
        
        _nmod_vec_clear(val);
        nmod_poly_clear(P);
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
