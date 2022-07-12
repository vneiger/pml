#include <flint/nmod_poly.h>

#include "nmod_poly_extra.h"

#ifdef HAS_AVX2

/*------------------------------------------------------------*/
/* initializes all entries of F                               */
/* w primitive and w^(2^order))=1                             */
/* DFTs of size up to 2^order are supported                   */ 
/*------------------------------------------------------------*/
void nmod_32_fft_init_set(nmod_32_fft_t F, mp_limb_t w, ulong order, nmod_t mod)
{
    mp_hlimb_t inv_2, inv;
    ulong k;
    
    F->order = order;
    F->w = (mp_hlimb_t) w;
    F->inv_w = (mp_hlimb_t) nmod_inv(w, mod);
    F->mod = mod;

    /* prepare the positive powers for inverse FFT                              */
    F->powers_w = flint_malloc((order+1) * sizeof(mp_hlimb_t *));
    F->powers_inv_w = flint_malloc((order+1) * sizeof(mp_hlimb_t *));
    F->powers_w_t = flint_malloc((order+1) * sizeof(mp_hlimb_t *));
    F->powers_inv_w_t = flint_malloc((order+1) * sizeof(mp_hlimb_t *));
    
    F->i_powers_w = flint_malloc((order+1) * sizeof(mp_hlimb_t *));
    F->i_powers_inv_w = flint_malloc((order+1) * sizeof(mp_hlimb_t *));
    F->i_powers_w_t = flint_malloc((order+1) * sizeof(mp_hlimb_t *));
    F->i_powers_inv_w_t = flint_malloc((order+1) * sizeof(mp_hlimb_t *));
    
    for (k = 0; k <= order; k++)
    {
        mp_hlimb_t wr, inv_wr;
        ulong j, K;
        
        wr = w;
        inv_wr = F->inv_w;
        
        for (j = 0; j < order-k; j++)
        {
            wr = nmod_mul(wr, wr, mod); 
            inv_wr = nmod_mul(inv_wr, inv_wr, mod); 
        }
        // wr = w^{2^(order-k)}
        
        K = 1L << k;
        if (k <= 2)
        {
            F->powers_w[k] = malloc(K * sizeof(mp_hlimb_t));
            F->powers_inv_w[k] = malloc(K * sizeof(mp_hlimb_t));
            F->powers_w_t[k] = malloc(K * sizeof(mp_hlimb_t));
            F->powers_inv_w_t[k] = malloc(K * sizeof(mp_hlimb_t));

            F->i_powers_w[k] = malloc(K * sizeof(mp_hlimb_t));
            F->i_powers_inv_w[k] = malloc(K * sizeof(mp_hlimb_t));
            F->i_powers_w_t[k] = malloc(K * sizeof(mp_hlimb_t));
            F->i_powers_inv_w_t[k] = malloc(K * sizeof(mp_hlimb_t));
        }
        if (k == 3)
        {
            F->powers_w[k] = aligned_alloc(32, 4 * K);
            F->powers_inv_w[k] = aligned_alloc(32, 4 * K);
            F->powers_w_t[k] = aligned_alloc(32, 4 * K);
            F->powers_inv_w_t[k] = aligned_alloc(32, 4 * K);

            F->i_powers_w[k] = aligned_alloc(32, 4 * K);
            F->i_powers_inv_w[k] = aligned_alloc(32, 4 * K);
            F->i_powers_w_t[k] = aligned_alloc(32, 4 * K);
            F->i_powers_inv_w_t[k] = aligned_alloc(32, 4 * K);
        }
        if (k >= 4)
        {
            F->powers_w[k] = aligned_alloc(64, 4 * K);
            F->powers_inv_w[k] = aligned_alloc(64, 4 * K);
            F->powers_w_t[k] = aligned_alloc(64, 4 * K);
            F->powers_inv_w_t[k] = aligned_alloc(64, 4 * K);

            F->i_powers_w[k] = aligned_alloc(64, 4 * K);
            F->i_powers_inv_w[k] = aligned_alloc(64, 4 * K);
            F->i_powers_w_t[k] = aligned_alloc(64, 4 * K);
            F->i_powers_inv_w_t[k] = aligned_alloc(64, 4 * K);
        }

        
        if (K == 1)
        {
            mp_hlimb_t t;
            F->powers_w[k][0] = 1;
            F->powers_inv_w[k][0] = 1;
            F->powers_w_t[k][0] = 1;
            F->powers_inv_w_t[k][0] = 1;

            t = prep_mul_mod_precon_32(1, mod.n);
            F->i_powers_w[k][0] = t;
            F->i_powers_inv_w[k][0] = t;
            F->i_powers_w_t[k][0] = t;
            F->i_powers_inv_w_t[k][0] = t;
        }
        else
        {
            ulong nb;
            nb = 0;
            K = K/2;
            
            while (K >= 1)
            {
                mp_hlimb_t wri, inv_wri;
                ulong i;

                wri = 1;     // wr^i
                inv_wri = 1;

                for (i = 0; i < K; i++)
                {
                    mp_hlimb_t p_wri, p_inv_wri;
                    F->powers_w[k][nb] = wri;
                    F->powers_w_t[k][K + i] = wri;
                    F->powers_inv_w_t[k][nb] = inv_wri;
                    F->powers_inv_w[k][K + i] = inv_wri;

                    p_wri = prep_mul_mod_precon_32(wri, mod.n);
                    p_inv_wri = prep_mul_mod_precon_32(inv_wri, mod.n);

                    F->i_powers_w[k][nb] = p_wri;
                    F->i_powers_w_t[k][K + i] = p_wri;
                    F->i_powers_inv_w_t[k][nb] = p_inv_wri;
                    F->i_powers_inv_w[k][K + i] = p_inv_wri;
                    
                    wri = nmod_mul(wri, wr, mod);
                    inv_wri = nmod_mul(inv_wri, inv_wr, mod);
                    
                    nb++;
                }

                wr = nmod_mul(wr, wr, mod);
                inv_wr = nmod_mul(inv_wr, inv_wr, mod);
                K /= 2;
            }
        }
        /* { */
        /*     ulong nb = 0; */
        /*     while (K >= 2) */
        /*     { */
        /*         mp_limb_t wi; */
        /*         ulong i; */

        /*         wi = 1;     // w^i */
        /*         for (i = 0; i < K/2; i++) */
        /*         { */
        /*             F->powers_w[k][nb] = wi; */
        /*             F->i_powers_w[k][nb] = prep_mul_mod_precon_32(wi, mod.n); */
        /*             wi = nmod_mul(wi, wr, mod); */
        /*             nb++; */
        /*         } */
        /*         wr = nmod_mul(wr, wr, mod); */
        /*         K /= 2; */
        /*     } */
        /* } */
    }
    

    
    /*     else */
    /*     { */
    /*         ulong KK = K; */
    /*         while (K >= 2) */
    /*         { */
    /*             mp_limb_t iwi; */
    /*             ulong i; */

    /*             iwi = 1;     // 1/w^i */
    /*             for (i = 0; i < K/2; i++) */
    /*             { */
    /*                 F->powers_inv_w[k][KK-K/2+i] = iwi; */
    /*                 F->i_powers_inv_w[k][KK-K/2+i] = prep_mul_mod_precon_32(iwi, mod.n); */
    /*                 iwi = nmod_mul(iwi, inv_wr, mod); */
    /*             } */
    /*             inv_wr = nmod_mul(inv_wr, inv_wr, mod); */
    /*             K /= 2; */
    /*             KK -= K; */
    /*         } */
    /*     } */
    /* } */

    F->powers_inv_2 = flint_malloc((order+1) * sizeof(mp_hlimb_t));
    F->i_powers_inv_2 = flint_malloc((order+1) * sizeof(mp_hlimb_t));
    inv_2 = (mp_hlimb_t) nmod_inv(2, mod);
    inv = 1;
    for (k = 0; k <= order; k++)
    {
        F->powers_inv_2[k] = inv;
        F->i_powers_inv_2[k] = prep_mul_mod_precon_32(inv, mod.n);
        inv = nmod_mul(inv, inv_2, mod);
    }


    F->powers_inv_w_over_2 = flint_malloc(sizeof(mp_hlimb_t *) * order);
    F->i_powers_inv_w_over_2 = flint_malloc(sizeof(mp_hlimb_t *) * order);
    for (k = 0; k < order; k++)
    {
        ulong i, K;
        mp_hlimb_t * src;
        
        K = 1L << k;
        F->powers_inv_w_over_2[k] = (mp_hlimb_t *) aligned_alloc(32, 4*K > 32 ? 4*K : 32);
        F->i_powers_inv_w_over_2[k] = (mp_hlimb_t *) aligned_alloc(32, 4*K > 32 ? 4*K : 32);
        src = F->powers_inv_w[k+1] + K;
        
        for (i = 0; i < K; i++)
        {
            F->powers_inv_w_over_2[k][i] = nmod_mul(src[i], F->powers_inv_2[k], mod);
            F->i_powers_inv_w_over_2[k][i] = prep_mul_mod_precon_32(F->powers_inv_w_over_2[k][i], mod.n);
        }
    }
}

#endif
