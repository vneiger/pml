#include <flint/flint.h>
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




void n_fft_ctx_init(n_fft_ctx_t F, ulong p)
{
    // modulus p should be prime, with 2 < p < 2**61 (TODO check the latter suffices)
    // we do not check primality
    if (p <= 2 || flint_clz(p) <= 3)
        flint_throw(FLINT_ERROR, "Provided prime p = %wu for n_fft does not satisfy bounds: 2 < p < 2**61", p);

    // find the constant and exponent such that p == c * 2**max_order + 1
    const ulong max_order = flint_ctz(p - UWORD(1));
    const ulong c = (p - UWORD(1)) >> max_order;

    // find primitive root w of order 2**max_order
    const ulong prim_root = n_primitive_root_prime(p);
    const ulong w = n_powmod2(prim_root, c, p);

    // build powers w**2**k
}

void n_fft_ctx_init_root(n_fft_ctx_t F, ulong root, ulong p)
{
}

/*------------------------------------------------------------*/
/* initializes all entries of F                               */
/* w primitive and w^(2^order))=1                             */
/* DFTs of size up to 2^order are supported                   */
/* order >= 3 required                                        */
/*------------------------------------------------------------*/

void n_fft_ctx_init_set(n_fft_ctx_t F, ulong w, ulong order, ulong mod)
{
    // basic attributes
    F->mod = mod;
    F->mod2 = 2*mod;
    F->mod4 = 4*mod;
    F->order = order;
    //F->w = w;
    //F->inv_w = nmod_inv(w, mod);  // TODO

    // fill table of powers of w:
    // buf[2*i] == w**i and buf[2*i+1] its precomp, i = 0 ... len-1   where len = 2**(order-1)
    // F->tab_w same in bit-reversed: 1, w**(len/2), w**(len/4), w**(3*len/4), ...
    ulong len = (UWORD(1) << (order-1));  // len == 2**(ell+1) >= 4
    nn_ptr buf = _nmod_vec_init(2*len);
    n_geometric_sequence_with_precomp(buf, w, len, mod);

    // put in bit reversed order for tab_w[1]
    F->tab_w = _nmod_vec_init(2*len);
    for (ulong k = 0, j = 0; k < len; k++, j = RevInc(j, order-1))
    {
        F->tab_w[2*k]   = buf[2*j];
        F->tab_w[2*k+1] = buf[2*j+1];
    }

    F->I     = F->tab_w[2];
    F->Ipre  = F->tab_w[3];
    F->J     = F->tab_w[4];
    F->Jpre  = F->tab_w[5];
    F->IJ    = F->tab_w[6];
    F->IJpre = F->tab_w[7];

    _nmod_vec_clear(buf);
}

void n_fft_ctx_clear(n_fft_ctx_t F)
{
    _nmod_vec_clear(F->tab_w);
}




/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
