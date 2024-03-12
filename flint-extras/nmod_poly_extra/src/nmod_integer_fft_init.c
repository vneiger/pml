#include <flint/nmod.h>
#include <flint/nmod_vec.h>
#include "nmod_poly_fft.h"

/*------------------------------------------------------------*/
/* initializes all entries of F                               */
/* w primitive and w^(2^order))=1                             */
/* DFTs of size up to 2^order are supported                   */
/* order >= 3 required                                        */
/*------------------------------------------------------------*/
void nmod_integer_fft_init_set(nmod_integer_fft_t F, mp_limb_t w, ulong order, nmod_t mod)
{
    // basic attributes
    F->mod = mod;
    F->order = order;
    F->w = w;
    F->inv_w = nmod_inv(w, mod);

    // 1. fill tables of powers of w
    F->tab_w = (mp_limb_t **) flint_malloc(sizeof(mp_limb_t *) * order);

    // fill largest array of powers of w
    slong ell = order-2;  // >= 1
    ulong len = (1 << (order-1));  // len == 2**(ell+1) >= 4
    F->tab_w[ell] = _nmod_vec_init(len);
    F->tab_w[ell][0] = UWORD(1);
    F->tab_w[ell][1] = w;
    F->tab_w[ell][2] = n_mulmod2_preinv(w, w, mod.n, mod.ninv);
    F->tab_w[ell][3] = n_mulmod2_preinv(F->tab_w[ell][2], w, mod.n, mod.ninv);
    if (order > 3)
    {
        mp_limb_t w4 = n_mulmod2_preinv(F->tab_w[ell][2], F->tab_w[ell][2], mod.n, mod.ninv);
        mp_limb_t w4_pr = n_mulmod_precomp_shoup(w4, mod.n);
        for (ulong k = 0; k+7 < len; k+=4)
        {
            F->tab_w[ell][k+4] = n_mulmod_shoup(w4, F->tab_w[ell][k+0], w4_pr, mod.n);
            F->tab_w[ell][k+5] = n_mulmod_shoup(w4, F->tab_w[ell][k+1], w4_pr, mod.n);
            F->tab_w[ell][k+6] = n_mulmod_shoup(w4, F->tab_w[ell][k+2], w4_pr, mod.n);
            F->tab_w[ell][k+7] = n_mulmod_shoup(w4, F->tab_w[ell][k+3], w4_pr, mod.n);
        }
        // finished here, k reached exactly len since len is a power of 2
    }

    // copy into other arrays
    // NOTES: did not try other methods
    //   -> if order up to ~10, this is negligible (<~10%) compared to the above
    //   -> beyond, its impact grows and reaches a factor around 2 up to order 30 at least
    //       (meaning the next for-ell loop takes about as much time as the above for-k loop)
    //ell--;
    //for (; ell >= 0; ell--)
    //{
    //    len = len >> 1;  // len == 2**(ell+1)
    //    F->tab_w[ell] = _nmod_vec_init(len);
    //    for (ulong k = 0; k < len; k++)
    //        F->tab_w[ell][k] = F->tab_w[ell+1][2*k];
    //}

    //// VERSION2:
    // ell--;
    // for (; ell >= 3; ell-=4)
    // {
    //     len = len >> 4;  // len == 2**(ell-3+1)
    //     F->tab_w[ell] = _nmod_vec_init(8*len);
    //     F->tab_w[ell-1] = _nmod_vec_init(4*len);
    //     F->tab_w[ell-2] = _nmod_vec_init(2*len);
    //     F->tab_w[ell-3] = _nmod_vec_init(len);
    //     for (ulong k = 0; k < len; k++)
    //     {
    //         F->tab_w[ell][8*k+0] = F->tab_w[ell+1][16*k+0];
    //         F->tab_w[ell][8*k+1] = F->tab_w[ell+1][16*k+2];
    //         F->tab_w[ell][8*k+2] = F->tab_w[ell+1][16*k+4];
    //         F->tab_w[ell][8*k+3] = F->tab_w[ell+1][16*k+6];
    //         F->tab_w[ell][8*k+4] = F->tab_w[ell+1][16*k+8];
    //         F->tab_w[ell][8*k+5] = F->tab_w[ell+1][16*k+10];
    //         F->tab_w[ell][8*k+6] = F->tab_w[ell+1][16*k+12];
    //         F->tab_w[ell][8*k+7] = F->tab_w[ell+1][16*k+14];
    //         F->tab_w[ell-1][4*k+0] = F->tab_w[ell+1][16*k+0];
    //         F->tab_w[ell-1][4*k+1] = F->tab_w[ell+1][16*k+4];
    //         F->tab_w[ell-1][4*k+2] = F->tab_w[ell+1][16*k+8];
    //         F->tab_w[ell-1][4*k+3] = F->tab_w[ell+1][16*k+12];
    //         F->tab_w[ell-2][2*k+0] = F->tab_w[ell+1][16*k+0];
    //         F->tab_w[ell-2][2*k+1] = F->tab_w[ell+1][16*k+8];
    //         F->tab_w[ell-3][k] = F->tab_w[ell+1][16*k+0];
    //     }
    // }
    // 
    // for (; ell >= 0; ell--)
    // {
    //     len = len >> 1;  // len == 2**(ell+1)
    //     F->tab_w[ell] = _nmod_vec_init(len);
    //     for (ulong k = 0; k < len; k++)
    //         F->tab_w[ell][k] = F->tab_w[ell+1][2*k];
    // }

    ell--;
    for (; ell >= 3; ell-=4)
    {
        len = len >> 4;  // len == 2**(ell-3+1)
        F->tab_w[ell] = _nmod_vec_init(8*len);
        F->tab_w[ell-1] = _nmod_vec_init(4*len);
        F->tab_w[ell-2] = _nmod_vec_init(2*len);
        F->tab_w[ell-3] = _nmod_vec_init(len);
        for (ulong k = 0; k < len; k++)
        {
            F->tab_w[ell+0][8*k+0] = F->tab_w[ell+1][16*k+0];
            F->tab_w[ell+0][8*k+1] = F->tab_w[ell+1][16*k+2];
            F->tab_w[ell+0][8*k+2] = F->tab_w[ell+1][16*k+4];
            F->tab_w[ell+0][8*k+3] = F->tab_w[ell+1][16*k+6];
            F->tab_w[ell+0][8*k+4] = F->tab_w[ell+1][16*k+8];
            F->tab_w[ell+0][8*k+5] = F->tab_w[ell+1][16*k+10];
            F->tab_w[ell+0][8*k+6] = F->tab_w[ell+1][16*k+12];
            F->tab_w[ell+0][8*k+7] = F->tab_w[ell+1][16*k+14];
            F->tab_w[ell-1][4*k+0] = F->tab_w[ell+0][8*k+0];
            F->tab_w[ell-1][4*k+1] = F->tab_w[ell+0][8*k+2];
            F->tab_w[ell-1][4*k+2] = F->tab_w[ell+0][8*k+4];
            F->tab_w[ell-1][4*k+3] = F->tab_w[ell+0][8*k+6];
            F->tab_w[ell-2][2*k+0] = F->tab_w[ell-1][4*k+0];
            F->tab_w[ell-2][2*k+1] = F->tab_w[ell-1][4*k+2];
            F->tab_w[ell-3][k] = F->tab_w[ell-2][2*k];
        }
    }

    for (; ell >= 0; ell--)
    {
        len = len >> 1;  // len == 2**(ell+1)
        F->tab_w[ell] = _nmod_vec_init(len);
        for (ulong k = 0; k < len; k++)
            F->tab_w[ell][k] = F->tab_w[ell+1][2*k];
    }

    //F->tab_inv_w = (mp_limb_t **) flint_malloc(sizeof(mp_limb_t *) * order);

    //mp_limb_t inv_2, inv;

    //F->powers_inv_2 = (mp_ptr) flint_malloc(sizeof(vec1d) * (order + 1));
    //inv_2 = nmod_inv(2, mod);
    //inv = 1;
    //for (k = 0; k <= order; k++)
    //{
    //    F->powers_inv_2[k] = inv;
    //    inv = nmod_mul(inv, inv_2, mod);
    //}

    //if (order > 0)
    //{
    //    F->powers_inv_w_over_2 = (vec1d **) flint_malloc(sizeof(vec1d *) * order);
    //    for (k = 0; k < order; k++)
    //    {
    //        ulong i, K;
    //        mp_limb_t *src;

    //        K = 1L << k;

    //        if ((sizeof(vec1d) * K) >= 32)
    //        {
    //            F->powers_inv_w_over_2[k] = (vec1d *) aligned_alloc(32, sizeof(vec1d) * K);
    //        }
    //        else
    //        {
    //            F->powers_inv_w_over_2[k] = (vec1d *) flint_malloc(sizeof(vec1d) * K);
    //        }
    //        src = F->powers_inv_w_t[k+1];

    //        for (i = 0; i < K; i++)
    //        {
    //            F->powers_inv_w_over_2[k][i] = vec1d_reduce_0n_to_pmhn((vec1d) nmod_mul(src[i], F->powers_inv_2[k], mod), dp);
    //        }
    //    }
    //}
}

void nmod_integer_fft_clear(nmod_integer_fft_t F)
{
    for (ulong ell = 0; ell <= F->order-2; ell++)
        _nmod_vec_clear(F->tab_w[ell]);
    flint_free(F->tab_w);
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
