#include <flint/nmod_poly.h>

#include "nmod_poly_extra.h"

#ifdef HAS_INT128

/*------------------------------------------------------------*/
/* initializes all entries of F                               */
/* w primitive and w^(2^order))=1                             */
/* DFTs of size up to 2^order are supported                   */ 
/*------------------------------------------------------------*/
void nmod_64_fft_init_set(nmod_64_fft_t F, mp_limb_t w, ulong order, nmod_t mod)
{
    mp_limb_t inv_2, inv;
    slong k;
    
    F->order = order;
    F->w = w;
    F->inv_w = nmod_inv(w, mod);
    F->mod = mod;

    /* prepare the positive powers for inverse FFT                              */
    F->powers_w = flint_malloc((order+1) * sizeof(mp_ptr));
    F->i_powers_w = flint_malloc((order+1) * sizeof(mp_ptr));
    
    for (k = 0; k <= order; k++)
    {
        mp_limb_t wr;
        slong j;
        ulong K;
        
        wr = w;
        for (j = 0; j < order-k; j++)
        {
            wr = nmod_mul(wr, wr, mod); 
        }
        // wr = w^{2^(order-k)}
        
        K = 1L << k;
        F->powers_w[k] = _nmod_vec_init(K);
        F->i_powers_w[k] = _nmod_vec_init(K);
        
        if (K == 1)
        {
            F->powers_w[k][0] = 1;
            F->i_powers_w[k][0] = prep_mul_mod_precon(1, mod.n);
        }
        else
        {
            slong nb = 0;
            while (K >= 2)
            {
                mp_limb_t wi;
                slong i;

                wi = 1;     // w^i
                for (i = 0; i < K/2; i++)
                {
                    F->powers_w[k][nb] = wi;
                    F->i_powers_w[k][nb] = prep_mul_mod_precon(wi, mod.n);
                    wi = nmod_mul(wi, wr, mod);
                    nb++;
                }
                wr = nmod_mul(wr, wr, mod);
                K /= 2;
            }
        }
    }
    

    /* prepare the negative powers for inverse FFT                              */
    /* the order of the levels is reversed                                      */
    F->powers_inv_w = flint_malloc((order+1) * sizeof(mp_ptr));
    F->i_powers_inv_w = flint_malloc((order+1) * sizeof(mp_ptr));
    
    for (k = 0; k <= order; k++)
    {
        mp_limb_t iwr;
        slong j;
        ulong K;
        
        iwr = F->inv_w;
        for (j = 0; j < order-k; j++)
        {
            iwr = nmod_mul(iwr, iwr, mod); 
        }
        
        K = 1L << k;
        F->powers_inv_w[k] = _nmod_vec_init(K);
        F->i_powers_inv_w[k] = _nmod_vec_init(K);
        
        if (K == 1)
        {
            F->powers_inv_w[k][0] = 1;
            F->i_powers_inv_w[k][0] = prep_mul_mod_precon(1, mod.n);
        }
        else
        {
            long KK = K;
            while (K >= 2)
            {
                mp_limb_t iwi;
                slong i;

                iwi = 1;     // 1/w^i
                for (i = 0; i < K/2; i++)
                {
                    F->powers_inv_w[k][KK-K/2+i] = iwi;
                    F->i_powers_inv_w[k][KK-K/2+i] = prep_mul_mod_precon(iwi, mod.n);
                    iwi = nmod_mul(iwi, iwr, mod);
                }
                iwr = nmod_mul(iwr, iwr, mod);
                K /= 2;
                KK -= K;
            }
        }
    }

    F->powers_inv_2 = _nmod_vec_init(order + 1);
    F->i_powers_inv_2 = _nmod_vec_init(order + 1);
    inv_2 = nmod_inv(2, mod);
    inv = 1;
    for (k = 0; k <= order; k++)
    {
        F->powers_inv_2[k] = inv;
        F->i_powers_inv_2[k] = prep_mul_mod_precon(inv, mod.n);
        inv = nmod_mul(inv, inv_2, mod);
    }
    
    F->powers_inv_w_over_2 = flint_malloc(sizeof(mp_ptr) * order);
    F->i_powers_inv_w_over_2 = flint_malloc(sizeof(mp_ptr) * order);
    for (k = 0; k < order; k++)
    {
        ulong i, K;
        mp_ptr src;
        
        K = 1L << k;
        F->powers_inv_w_over_2[k] = _nmod_vec_init(K);
        F->i_powers_inv_w_over_2[k] = _nmod_vec_init(K);
        src = F->powers_inv_w[k+1] + K;
        
        for (i = 0; i < K; i++)
        {
            F->powers_inv_w_over_2[k][i] = nmod_mul(src[i], F->powers_inv_2[k], mod);
            F->i_powers_inv_w_over_2[k][i] = prep_mul_mod_precon(F->powers_inv_w_over_2[k][i], mod.n);
        }
    }

}

#endif
