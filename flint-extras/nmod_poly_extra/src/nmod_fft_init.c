#include <flint/nmod.h>
#include <flint/nmod_vec.h>

#include "nmod_poly_fft.h"

/***********************
*  bit reversed copy  *
***********************/

inline long RevInc(long a, long k)
{
    long j, m;

    j = k;
    m = 1L << (k-1);

    while (j && (m & a)) {
        a ^= m;
        m >>= 1;
        j--;
    }
    if (j) a ^= m;
    return a;
}

// indices initialized with length >= k
//static inline void brc_indices(ulong * indices, long k)
//{
//    const long n = (1L << k);
//    for (long i = 0, j = 0; i < n; i++, j = RevInc(j, k))
//        indices[i] = j;
//}

//// counts in bit reversed order, in C
//// or faster to build list as in dft.sage??
//void iter_reversed(ulong bits) {
//    ulong n = 1 << bits;
//
//    for (ulong i = 0, j = 0; i < n; i++) {
//        printf("%ld\n", j);
//
//        // Compute a mask of LSBs.
//        ulong mask = i ^ (i + 1);
//        // Length of the mask.
//        ulong len = __builtin_ctz(~mask);
//        // Align the mask to MSB of n.
//        mask <<= bits - len;
//        // XOR with mask.
//        j ^= mask;
//    }
//}




/*------------------------------------------------------------*/
/* initializes all entries of F                               */
/* w primitive and w^(2^order))=1                             */
/* DFTs of size up to 2^order are supported                   */
/* order >= 3 required                                        */
/*------------------------------------------------------------*/
void nmod_fft_ctx_init_set(nmod_fft_ctx_t F, ulong w, ulong order, nmod_t mod)
{
    // basic attributes
    F->mod = mod;
    F->modn2 = 2*mod.n;
    F->modn4 = 4*mod.n;
    F->order = order;
    F->w = w;
    F->inv_w = nmod_inv(w, mod);

    // 1. fill tables of powers of w

    // fill largest array of powers of w
    slong ell = order-2;  // >= 1
    ulong len = (1 << (order-1));  // len == 2**(ell+1) >= 4
    ulong w_pr_quo, w_pr_rem;
    n_mulmod_precomp_shoup_quo_rem(&w_pr_quo, &w_pr_rem, w, mod.n);

    F->tab_w[ell] = _nmod_vec_init(len);
    F->tab_w_pre[ell] = _nmod_vec_init(len);
    F->tab_w[ell][0] = UWORD(1);
    F->tab_w_pre[ell][0] = n_mulmod_precomp_shoup(UWORD(1), mod.n);
    F->tab_w[ell][1] = w;
    F->tab_w_pre[ell][1] = w_pr_quo;
    n_mulmod_and_precomp_shoup(F->tab_w[ell] + 2, F->tab_w_pre[ell] + 2, w, w, w_pr_quo, w_pr_rem, w_pr_quo, mod.n);
    n_mulmod_and_precomp_shoup(F->tab_w[ell] + 3, F->tab_w_pre[ell] + 3, w, F->tab_w[ell][2], w_pr_quo, w_pr_rem, F->tab_w_pre[ell][2], mod.n);
    if (order > 3)
    {
        n_mulmod_and_precomp_shoup(&w, &w_pr_quo, w, F->tab_w[ell][3], w_pr_quo, w_pr_rem, F->tab_w_pre[ell][3], mod.n);
        w_pr_rem = n_mulmod_precomp_shoup_rem_from_quo(w_pr_quo, mod.n);
        for (ulong k = 0; k+7 < len; k+=4)
        {
            n_mulmod_and_precomp_shoup(F->tab_w[ell] + k+4, F->tab_w_pre[ell] + k+4, w, F->tab_w[ell][k+0], w_pr_quo, w_pr_rem, F->tab_w_pre[ell][k+0], mod.n);
            n_mulmod_and_precomp_shoup(F->tab_w[ell] + k+5, F->tab_w_pre[ell] + k+5, w, F->tab_w[ell][k+1], w_pr_quo, w_pr_rem, F->tab_w_pre[ell][k+1], mod.n);
            n_mulmod_and_precomp_shoup(F->tab_w[ell] + k+6, F->tab_w_pre[ell] + k+6, w, F->tab_w[ell][k+2], w_pr_quo, w_pr_rem, F->tab_w_pre[ell][k+2], mod.n);
            n_mulmod_and_precomp_shoup(F->tab_w[ell] + k+7, F->tab_w_pre[ell] + k+7, w, F->tab_w[ell][k+3], w_pr_quo, w_pr_rem, F->tab_w_pre[ell][k+3], mod.n);
        }
        // finished here, k reached exactly len since len is a power of 2
    }

    // copy into other arrays
    // NAIVE VERSION:
    //   -> if order up to ~10, this is negligible (<~10%) compared to the above
    //   -> beyond, its impact grows and reaches a factor around 2 up to order 30 at least
    //       (meaning the next for-ell loop takes about as much time as the above for-k loop)
    ell--;
    for (; ell >= 0; ell--)
    {
        len = len >> 1;  // len == 2**(ell+1)
        F->tab_w[ell] = _nmod_vec_init(len);
        F->tab_w_pre[ell] = _nmod_vec_init(len);
        for (ulong k = 0; k < len; k++)
        {
            F->tab_w[ell][k] = F->tab_w[ell+1][2*k];
            F->tab_w_pre[ell][k] = F->tab_w_pre[ell+1][2*k];
        }
    }

    F->J  = F->tab_w[1][1];
    F->I  = F->tab_w[1][2];
    F->IJ = F->tab_w[1][3];
    F->Jpre  = F->tab_w_pre[1][1];
    F->Ipre  = F->tab_w_pre[1][2];
    F->IJpre = F->tab_w_pre[1][3];
}

void nmod_fft_ctx_init_set_new(nmod_fft_ctx_t F, ulong w, ulong order, nmod_t mod)
{
    // basic attributes
    F->mod = mod;
    F->modn2 = 2*mod.n;
    F->modn4 = 4*mod.n;
    F->order = order;
    F->w = w;
    F->inv_w = nmod_inv(w, mod);

    // 1. fill tables of powers of w

    // fill largest array of powers of w
    slong ell = order-2;  // >= 1
    ulong len = (1 << (order-1));  // len == 2**(ell+1) >= 4
    ulong w_pr_quo, w_pr_rem;
    n_mulmod_precomp_shoup_quo_rem(&w_pr_quo, &w_pr_rem, w, mod.n);

    F->tab_w[ell] = _nmod_vec_init(2*len);
    F->tab_w_pre[ell] = _nmod_vec_init(0);
    F->tab_w[ell][0] = UWORD(1);
    F->tab_w[ell][1] = n_mulmod_precomp_shoup(UWORD(1), mod.n);
    F->tab_w[ell][2] = w;
    F->tab_w[ell][3] = w_pr_quo;
    n_mulmod_and_precomp_shoup(F->tab_w[ell] + 4, F->tab_w[ell] + 5, w, w, w_pr_quo, w_pr_rem, w_pr_quo, mod.n);
    n_mulmod_and_precomp_shoup(F->tab_w[ell] + 6, F->tab_w[ell] + 7, w, F->tab_w[ell][4], w_pr_quo, w_pr_rem, F->tab_w[ell][5], mod.n);
    if (order > 3)
    {
        n_mulmod_and_precomp_shoup(&w, &w_pr_quo, w, F->tab_w[ell][6], w_pr_quo, w_pr_rem, F->tab_w[ell][7], mod.n);
        w_pr_rem = n_mulmod_precomp_shoup_rem_from_quo(w_pr_quo, mod.n);
        for (ulong k = 0; 2*(k+7)+1 < 2*len; k+=4)
        {
            n_mulmod_and_precomp_shoup(F->tab_w[ell] + 2*(k+4), F->tab_w[ell] + 2*(k+4)+1, w, F->tab_w[ell][2*(k+0)], w_pr_quo, w_pr_rem, F->tab_w[ell][2*(k+0)+1], mod.n);
            n_mulmod_and_precomp_shoup(F->tab_w[ell] + 2*(k+5), F->tab_w[ell] + 2*(k+5)+1, w, F->tab_w[ell][2*(k+1)], w_pr_quo, w_pr_rem, F->tab_w[ell][2*(k+1)+1], mod.n);
            n_mulmod_and_precomp_shoup(F->tab_w[ell] + 2*(k+6), F->tab_w[ell] + 2*(k+6)+1, w, F->tab_w[ell][2*(k+2)], w_pr_quo, w_pr_rem, F->tab_w[ell][2*(k+2)+1], mod.n);
            n_mulmod_and_precomp_shoup(F->tab_w[ell] + 2*(k+7), F->tab_w[ell] + 2*(k+7)+1, w, F->tab_w[ell][2*(k+3)], w_pr_quo, w_pr_rem, F->tab_w[ell][2*(k+3)+1], mod.n);
        }
        // finished here, k reached exactly len since len is a power of 2
    }

    // copy into other arrays
    // NAIVE VERSION:
    //   -> if order up to ~10, this is negligible (<~10%) compared to the above
    //   -> beyond, its impact grows and reaches a factor around 2 up to order 30 at least
    //       (meaning the next for-ell loop takes about as much time as the above for-k loop)
    ell--;
    for (; ell >= 0; ell--)
    {
        len = len >> 1;  // len == 2**(ell+1)
        F->tab_w[ell] = _nmod_vec_init(2*len);
        F->tab_w_pre[ell] = _nmod_vec_init(1);
        for (ulong k = 0; k < len; k++)
        {
            F->tab_w[ell][2*k] = F->tab_w[ell+1][4*k];
            F->tab_w[ell][2*k+1] = F->tab_w[ell+1][4*k+1];
        }
    }

    F->J  = F->tab_w[1][2];
    F->I  = F->tab_w[1][4];
    F->IJ = F->tab_w[1][6];
    F->Jpre  = F->tab_w[1][3];
    F->Ipre  = F->tab_w[1][5];
    F->IJpre = F->tab_w[1][7];
}


void nmod_fft_ctx_init_set_red(nmod_fft_ctx_t F, ulong w, ulong order, nmod_t mod)
{
    // basic attributes
    F->mod = mod;
    F->modn2 = 2*mod.n;
    F->modn4 = 4*mod.n;
    F->order = order;
    F->w = w;
    F->inv_w = nmod_inv(w, mod);

    // 1. fill tables of powers of w
    // F->tab_w[0][i] == w**i, i = 0 ... len-1   where len = 2**(order-1)
    // F->tab_w[1] same in bit-reversed: 1, w**(len/2), w**(len/4), w**(3*len/4), ...
    ulong len = (1UL << (order-1));  // len == 2**(ell+1) >= 4
    F->tab_w[0] = _nmod_vec_init(len);
    F->tab_w_pre[0] = _nmod_vec_init(len);
    F->tab_w[0][0] = UWORD(1);
    F->tab_w[0][1] = w;
    F->tab_w[0][2] = n_mulmod2_preinv(w, w, mod.n, mod.ninv);
    F->tab_w[0][3] = n_mulmod2_preinv(F->tab_w[0][2], w, mod.n, mod.ninv);
    F->tab_w_pre[0][0] = n_mulmod_precomp_shoup(F->tab_w[0][0], mod.n);
    F->tab_w_pre[0][1] = n_mulmod_precomp_shoup(F->tab_w[0][1], mod.n);
    F->tab_w_pre[0][2] = n_mulmod_precomp_shoup(F->tab_w[0][2], mod.n);
    F->tab_w_pre[0][3] = n_mulmod_precomp_shoup(F->tab_w[0][3], mod.n);
    if (order > 3)
    {
        ulong w4 = n_mulmod2_preinv(F->tab_w[0][2], F->tab_w[0][2], mod.n, mod.ninv);
        ulong w4_pr = n_mulmod_precomp_shoup(w4, mod.n);
        for (ulong k = 0; k+7 < len; k+=4)
        {
            F->tab_w[0][k+4] = n_mulmod_shoup(w4, F->tab_w[0][k+0], w4_pr, mod.n);
            F->tab_w[0][k+5] = n_mulmod_shoup(w4, F->tab_w[0][k+1], w4_pr, mod.n);
            F->tab_w[0][k+6] = n_mulmod_shoup(w4, F->tab_w[0][k+2], w4_pr, mod.n);
            F->tab_w[0][k+7] = n_mulmod_shoup(w4, F->tab_w[0][k+3], w4_pr, mod.n);
            F->tab_w_pre[0][k+4] = n_mulmod_precomp_shoup(F->tab_w[0][k+4], mod.n);
            F->tab_w_pre[0][k+5] = n_mulmod_precomp_shoup(F->tab_w[0][k+5], mod.n);
            F->tab_w_pre[0][k+6] = n_mulmod_precomp_shoup(F->tab_w[0][k+6], mod.n);
            F->tab_w_pre[0][k+7] = n_mulmod_precomp_shoup(F->tab_w[0][k+7], mod.n);
        }
        // finished here, k reached exactly len since len is a power of 2
    }

    // put in bit reversed order for tab_w[1]
    F->tab_w[1] = _nmod_vec_init(len);
    F->tab_w_pre[1] = _nmod_vec_init(len);
    for (ulong k = 0, j = 0; k < len; k++, j = RevInc(j, order-1))
    {
        F->tab_w[1][k] = F->tab_w[0][j];
        F->tab_w_pre[1][k] = F->tab_w_pre[0][j];
    }

    F->J  = F->tab_w[1][2];
    F->I  = F->tab_w[1][1];
    F->IJ = F->tab_w[1][3];
    F->Jpre  = F->tab_w_pre[1][2];
    F->Ipre  = F->tab_w_pre[1][1];
    F->IJpre = F->tab_w_pre[1][3];
}

void nmod_fft_ctx_clear(nmod_fft_ctx_t F)
{
    for (ulong ell = 0; ell <= F->order-2; ell++)
    {
        _nmod_vec_clear(F->tab_w[ell]);
        _nmod_vec_clear(F->tab_w_pre[ell]);
    }
}

void nmod_fft_ctx_clear_red(nmod_fft_ctx_t F)
{
    for (ulong ell = 0; ell < 2; ell++)
    {
        _nmod_vec_clear(F->tab_w[ell]);
        _nmod_vec_clear(F->tab_w_pre[ell]);
    }
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
