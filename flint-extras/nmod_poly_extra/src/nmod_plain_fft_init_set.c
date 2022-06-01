#include <flint/nmod_poly.h>

#include "nmod_poly_extra.h"

/*------------------------------------------------------------*/
/* initializes all entries of F                               */
/* w primitive and w^(2^order))=1                             */
/* DFTs of size up to 2^order are supported                   */ 
/*------------------------------------------------------------*/
void nmod_plain_fft_init_set(nmod_plain_fft_t F, mp_limb_t w, ulong order, nmod_t mod)
{
    mp_limb_t inv_2, inv;
    slong k;

    F->order = order;
    F->w = w;
    F->inv_w = nmod_inv(w, mod);
    F->mod = mod;

    F->powers_w = flint_malloc(sizeof(mp_ptr) * (order + 1));
    F->powers_inv_w = flint_malloc(sizeof(mp_ptr) * (order + 1));

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
        
        if (K == 1)
        {
            F->powers_w[k][0] = 1;
        }
        else
        {
            slong nb = 0;
            while (K >= 2)
            {
                mp_limb_t wri;
                wri = 1;     // wr^i
                slong i;
                for (i = 0; i < K/2; i++)
                {
                    F->powers_w[k][nb] = wri;
                    wri = nmod_mul(wri, wr, mod);
                    nb++;
                }
                wr = nmod_mul(wr, wr, mod);
                K /= 2;
            }
        }
    }


    for (k = 0; k <= order; k++)
    {
        mp_limb_t inv_wr;
        slong j;

        inv_wr = F->inv_w;
        for (j = 0; j < order-k; j++)
        {
            inv_wr = nmod_mul(inv_wr, inv_wr, mod);
        } 
        // inv_wr = 1/w^{2^(order-k)}
        
        long K = 1L << k;
        F->powers_inv_w[k] = _nmod_vec_init(K);
        
        if (K == 1)
        {
            F->powers_inv_w[k][0] = 1;
        }
        else
        {
            slong KK;
            KK = K;
            while (K >= 2)
            {
                mp_limb_t inv_wri;
                slong i;
                inv_wri = 1; // 1/wr^i
                for (i = 0; i < K/2; i++)
                {
                    F->powers_inv_w[k][KK - K/2 + i] = inv_wri;
                    inv_wri = nmod_mul(inv_wri, inv_wr, mod);
                }
                inv_wr = nmod_mul(inv_wr, inv_wr, mod);
                K /= 2;
                KK -= K;
            }
        }
    }


    
    F->powers_inv_2 = _nmod_vec_init(order + 1);
    inv_2 = nmod_inv(2, mod);
    inv = 1;
    for (k = 0; k <= order; k++)
    {
        F->powers_inv_2[k] = inv;
        inv = nmod_mul(inv, inv_2, mod);
    }
}
