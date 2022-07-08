#include <assert.h>
#include <flint/flint.h>

#include "nmod_extra.h"
#include "nmod_poly_extra.h"

/*------------------------------------------------------------*/
/* returns reverse_nb(u)                                      */
/*------------------------------------------------------------*/
ulong bit_reverse(ulong u, ulong nb)
{
    ulong res, i;
    res = 0;
    for (i = 0; i < nb; i++)
    {
        res <<= 1;
        res |= (u & 1);
        u >>= 1;
    }
    return res;
}

/*------------------------------------------------------------*/
/* FFT evaluation, compared to naive evaluation               */
/*------------------------------------------------------------*/
void check()
{
#ifdef HAS_AVX512
    ulong order, order_max, N;
    slong i;
    flint_rand_t state;
    mp_limb_t p, w0, w, w2;
    nmod_t mod;
    nmod_32_fft_t F;
    mp_ptr val;
    nmod_poly_t P;
    
    flint_randinit(state);

    p = 7340033;
    nmod_init(&mod, p);
    w0 = 3308891;
    order_max = 10;
    w = nmod_pow_ui(w0, 1L<<(16-order_max), mod);
    nmod_32_fft_init_set(F, w, order_max, mod);
    
    for (order = 0; order <= order_max; order++)
    {
        w = nmod_pow_ui(w0, 1L<<(16-order), mod);
        N = 1L << order;
        nmod_poly_init2(P, p, N);
        val = _nmod_vec_init(N);
        
        for (i = 0; i < N; i++)
        {
            nmod_poly_set_coeff_ui(P, i, n_randtest(state) % p);
        }
        
        nmod_avx512_32_fft_evaluate(val, P, F, order);

        w2 = 1;
        for (i = 0; i < N; i++)
        {
            assert (val[bit_reverse(i, order)] ==  nmod_poly_evaluate_nmod(P, w2));
            w2 = nmod_mul(w2, w, mod);
        }
        _nmod_vec_clear(val);
        nmod_poly_clear(P);
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

