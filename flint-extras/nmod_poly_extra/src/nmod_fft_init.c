#include <flint/flint.h>
#include <flint/longlong.h>
#include <flint/nmod.h>
#include <flint/nmod_vec.h>
#include <flint/ulong_extras.h>

#include "nmod_poly_fft.h"

/***********************
*  bit reversed copy  *
***********************/

//inline long RevInc(long a, long k)
//{
//    long j, m;
//
//    j = k;
//    m = 1L << (k-1);
//
//    while (j && (m & a)) {
//        a ^= m;
//        m >>= 1;
//        j--;
//    }
//    if (j) a ^= m;
//    return a;
//}

// indices initialized with length >= k
//static inline void brc_indices(ulong * indices, long k)
//{
//    const long n = (1L << k);
//    for (long i = 0, j = 0; i < n; i++, j = RevInc(j, k))
//        indices[i] = j;
//}

//// counts in bit reversed depth, in C
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
/* main initialization function                               */
/* initializes all entries of F                               */
/* w primitive root of 1, of order 2**max_depth               */
/* max_depth must be >= 3 (so, 8 must divide p-1)             */
/* p is a prime with p < 2**61 TODO confirm                   */
/* FIXME currently large depth can lead to heavy memory usage */
/*                                                            */
/* none of the requirements is checked                        */
/*                                                            */
/* after this, DFTs of length up to 2**depth are supported    */
/*------------------------------------------------------------*/

void n_fft_ctx_init2_root(n_fft_ctx_t F, ulong w, ulong max_depth, ulong depth, ulong p)
{
    if (depth < 3)
        depth = 3;
    if (max_depth < depth)
        depth = max_depth;

    // fill basic attributes
    F->mod = p;
    F->mod2 = 2*p;
    F->mod4 = 4*p;
    F->max_depth = max_depth;
    F->depth = depth;

    // fill tab_w2
    ulong pr_quo, pr_rem, ww;
    ww = w;
    n_mulmod_precomp_shoup_quo_rem(&pr_quo, &pr_rem, ww, p);
    F->tab_w2[2*(max_depth-2)] = ww;
    F->tab_w2[2*(max_depth-2)+1] = pr_quo;
    for (slong k = max_depth-3; k >= 0; k--)
    {
        // ww <- ww**2 and its precomputed quotient
        n_mulmod_and_precomp_shoup(&ww, &pr_quo, ww, ww, pr_quo, pr_rem, pr_quo, p);
        pr_rem = n_mulmod_precomp_shoup_rem_from_quo(pr_quo, p);
        F->tab_w2[2*k] = ww;
        F->tab_w2[2*k+1] = pr_quo;
    }
    // at this stage, pr_quo and pr_rem are for k == 0 i.e. for I == tab_w2[0]

    // fill I, J, IJ
    F->I    = F->tab_w2[0];
    F->I_pr = F->tab_w2[1];
    F->J    = F->tab_w2[2];
    F->J_pr = F->tab_w2[3];
    n_mulmod_and_precomp_shoup(&F->IJ, &F->IJ_pr, F->I, F->J, pr_quo, pr_rem, F->J_pr, p);

    // fill tab_w
    ulong len = (UWORD(1) << (depth-1));  // len >= 4
    F->tab_w = _nmod_vec_init(2*len);

    F->tab_w[0] = UWORD(1);
    F->tab_w[1] = n_mulmod_precomp_shoup(UWORD(1), p);
    F->tab_w[2] = F->I;
    F->tab_w[3] = F->I_pr;
    F->tab_w[4] = F->J;
    F->tab_w[5] = F->J_pr;
    F->tab_w[6] = F->IJ;
    F->tab_w[7] = F->IJ_pr;

    // tab_w[2*4:2*8] is w**(L/16) * tab_w[2*0:2*4] where L = 2**max_depth,
    // tab_w[2*8:2*16] is w**(L/32) * tab_w[2*0:2*8], etc.
    // recall tab_w2[2*d] == w**(L / 2**(d+2))
    ulong d = 2;  // we start with w**(L/16) == sqrt(J)
    for (ulong llen = 4; llen < len; llen <<= 1, d += 1)
    {
        ww = F->tab_w2[2*d];
        pr_quo = F->tab_w2[2*d+1];
        pr_rem = n_mulmod_precomp_shoup_rem_from_quo(pr_quo, p);
        // for each k, tab_w[2*(k+llen)] <- ww * tab_w[2*k], and deduce precomputation
        for (ulong k = 0; k < llen; k++)
            n_mulmod_and_precomp_shoup(F->tab_w + 2*llen + 2*k, F->tab_w + 2*llen + 2*k+1,
                                        ww, F->tab_w[2*k],
                                        pr_quo, pr_rem, F->tab_w[2*k+1], p);
    }
}

void n_fft_ctx_init2(n_fft_ctx_t F, ulong depth, ulong p)
{
    FLINT_ASSERT(p > 2 && flint_clz(p) >= 3);    // 2 < p < 2**61
    FLINT_ASSERT(flint_ctz(p - UWORD(1)) >= 3);  // p-1 divisible by 8

    // find the constant and exponent such that p == c * 2**max_depth + 1
    const ulong max_depth = flint_ctz(p - UWORD(1));
    const ulong c = (p - UWORD(1)) >> max_depth;

    // find primitive root w of order 2**max_depth
    const ulong prim_root = n_primitive_root_prime(p);
    const ulong w = n_powmod2(prim_root, c, p);

    // fill all attributes and tables
    n_fft_ctx_init2_root(F, w, max_depth, depth, p);
}

void n_fft_ctx_clear(n_fft_ctx_t F)
{
    _nmod_vec_clear(F->tab_w);
}




/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
