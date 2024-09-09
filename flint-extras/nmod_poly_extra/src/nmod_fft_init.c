#include <flint/nmod.h>
#include <flint/nmod_vec.h>
#include <flint/ulong_extras.h>

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
void nmod_fft_ctx_init_set(nmod_fft_ctx_t F, ulong w, ulong order, ulong mod)
{
    // basic attributes
    F->mod = mod;
    F->mod2 = 2*mod;
    F->mod4 = 4*mod;
    F->order = order;
    F->w = w;
    //F->inv_w = nmod_inv(w, mod);  // TODO

    // 1. fill tables of powers of w

    // fill largest array of powers of w
    slong ell = order-2;  // >= 1
    ulong len = (1 << (order-1));  // len == 2**(ell+1) >= 4

    F->tab_w[ell] = _nmod_vec_init(2*len);
    _n_geometric_sequence_with_precomp(F->tab_w[ell], w, len, mod);

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
        for (ulong k = 0; k < len; k++)
        {
            F->tab_w[ell][2*k] = F->tab_w[ell+1][4*k];
            F->tab_w[ell][2*k+1] = F->tab_w[ell+1][4*k+1];
        }
    }

    F->J     = F->tab_w[1][2];
    F->Jpre  = F->tab_w[1][3];
    F->I     = F->tab_w[1][4];
    F->Ipre  = F->tab_w[1][5];
    F->IJ    = F->tab_w[1][6];
    F->IJpre = F->tab_w[1][7];
}

void nmod_fft_ctx_init_set_new(nmod_fft_ctx_t F, ulong w, ulong order, ulong mod)
{
    // basic attributes
    F->mod = mod;
    F->mod2 = 2*mod;
    F->mod4 = 4*mod;
    F->order = order;
    F->w = w;
    //F->inv_w = nmod_inv(w, mod);  // TODO

    // 1. fill tables of powers of w

    // fill largest array of powers of w
    slong ell = order-2;  // >= 1
    ulong len = (1 << (order-1));  // len == 2**(ell+1) >= 4

    F->tab_w[ell] = _nmod_vec_init(2*len);
    _n_geometric_sequence_with_precomp(F->tab_w[ell], w, len, mod);

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


void nmod_fft_ctx_init_set_red(nmod_fft_ctx_t F, ulong w, ulong order, ulong mod)
{
    // basic attributes
    F->mod = mod;
    F->mod2 = 2*mod;
    F->mod4 = 4*mod;
    F->order = order;
    F->w = w;
    //F->inv_w = nmod_inv(w, mod);  // TODO

    // fill table of powers of w:
    // F->tab_w[0][2*i] == w**i and F->tab_w[0][2*i+1] its precomp, i = 0 ... len-1   where len = 2**(order-1)
    // F->tab_w[1] same in bit-reversed: 1, w**(len/2), w**(len/4), w**(3*len/4), ...
    ulong len = (1UL << (order-1));  // len == 2**(ell+1) >= 4
    F->tab_w[0] = _nmod_vec_init(2*len);
    _n_geometric_sequence_with_precomp(F->tab_w[0], w, len, mod);

    // put in bit reversed order for tab_w[1]
    F->tab_w[1] = _nmod_vec_init(2*len);
    for (ulong k = 0, j = 0; k < len; k++, j = RevInc(j, order-1))
    {
        F->tab_w[1][2*k] = F->tab_w[0][2*j];
        F->tab_w[1][2*k+1] = F->tab_w[0][2*j+1];
    }

    F->I     = F->tab_w[1][2];
    F->Ipre  = F->tab_w[1][3];
    F->J     = F->tab_w[1][4];
    F->Jpre  = F->tab_w[1][5];
    F->IJ    = F->tab_w[1][6];
    F->IJpre = F->tab_w[1][7];
}

void nmod_fft_ctx_clear(nmod_fft_ctx_t F)
{
    for (ulong ell = 0; ell <= F->order-2; ell++)
        _nmod_vec_clear(F->tab_w[ell]);
}

void nmod_fft_ctx_clear_new(nmod_fft_ctx_t F)
{
    for (ulong ell = 0; ell <= F->order-2; ell++)
        _nmod_vec_clear(F->tab_w[ell]);
}


void nmod_fft_ctx_clear_red(nmod_fft_ctx_t F)
{
    for (ulong ell = 0; ell < 2; ell++)
        _nmod_vec_clear(F->tab_w[ell]);
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
