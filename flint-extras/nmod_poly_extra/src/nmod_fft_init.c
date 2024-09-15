#include <flint/flint.h>
#include <flint/longlong.h>
#include <flint/nmod.h>
#include <flint/nmod_vec.h>
#include <flint/ulong_extras.h>

#include "nmod_poly_fft.h"

#define N_FFT_CTX_DEFAULT_DEPTH 12

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
/* initializes all entries of F                               */
/* w primitive root of 1, of order 2**depth                   */
/* DFTs of size up to 2^depth are supported                   */
/* depth >= 3 required                                        */
/*------------------------------------------------------------*/

void n_fft_ctx_init2_root(n_fft_ctx_t F, ulong w, ulong w_depth, ulong depth, ulong p)
{
    // fill basic attributes
    F->mod = p;
    F->mod2 = 2*p;
    F->mod4 = 4*p;
    F->depth = w_depth;
    F->alloc = depth;

    // fill tab_ww
    ulong pr_quo, pr_rem, ww;
    ww = w;
    n_mulmod_precomp_shoup_quo_rem(&pr_quo, &pr_rem, ww, p);
    F->tab_ww[2*(depth-2)] = ww;
    F->tab_ww[2*(depth-2)+1] = pr_quo;
    for (slong k = depth-3; k >= 0; k--)
    {
        // ww <- ww**2 and its precomputed quotient
        n_mulmod_and_precomp_shoup(&ww, &pr_quo, ww, ww, pr_quo, pr_rem, pr_quo, p);
        pr_rem = n_mulmod_precomp_shoup_rem_from_quo(pr_quo, p);
        F->tab_ww[2*k] = ww;
        F->tab_ww[2*k+1] = pr_quo;
    }
    // at this stage, pr_quo and pr_rem are for k == 0 i.e. for I == tab_ww[0]

    // fill I, J, IJ
    F->I    = F->tab_ww[0];
    F->I_pr = F->tab_ww[1];
    F->J    = F->tab_ww[2];
    F->J_pr = F->tab_ww[3];
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

    ulong d = 2;  // tab_ww[2*d] == sqrt(J)
    for (ulong llen = 4; llen < len; llen <<= 1, d += 1)
    {
        ww = F->tab_ww[2*d];
        pr_quo = F->tab_ww[2*d+1];
        pr_rem = n_mulmod_precomp_shoup_rem_from_quo(pr_quo, p);
        // for each k, tab_w[2*(k+llen)] <- ww * tab_w[2*k], and deduce precomputation
        for (ulong k = 0; k < llen; k++)
            n_mulmod_and_precomp_shoup(F->tab_w + 2*llen + 2*k, F->tab_w + 2*llen + 2*k+1,
                                       ww, F->tab_w[2*k],
                                       pr_quo, pr_rem, F->tab_w[2*k+1], p);
    }
}


void n_fft_ctx_init_root(n_fft_ctx_t F, ulong w, ulong depth, ulong p)
{
    n_fft_ctx_init2_root(F, w, depth, FLINT_MIN(depth, N_FFT_CTX_DEFAULT_DEPTH), p);
}

void n_fft_ctx_init(n_fft_ctx_t F, ulong p)
{
    // modulus p should be prime, with 2 < p < 2**61 (TODO check the latter suffices)
    // we do not check primality
    if (p <= 2 || flint_clz(p) <= 3)
        flint_throw(FLINT_ERROR, "Provided prime p = %wu for n_fft does not satisfy bounds: 2 < p < 2**61", p);

    // find the constant and exponent such that p == c * 2**depth + 1
    const ulong depth = flint_ctz(p - UWORD(1));
    const ulong c = (p - UWORD(1)) >> depth;
    if (depth <= 3)
        flint_throw(FLINT_ERROR, "Provided prime p = %wu for n_fft does not satisfy `8 divides p-1`", p);

    // find primitive root w of depth 2**depth
    const ulong prim_root = n_primitive_root_prime(p);
    const ulong w = n_powmod2(prim_root, c, p);

    // fill all attributes and tables, with default depth
    n_fft_ctx_init_root(F, w, depth, p);
}

void n_fft_ctx_init2(n_fft_ctx_t F, ulong len, ulong p)
{
    // modulus p should be prime, with 2 < p < 2**61 (TODO check the latter suffices)
    // we do not check primality
    if (p <= 2 || flint_clz(p) <= 3)
        flint_throw(FLINT_ERROR, "Provided prime p = %wu for n_fft does not satisfy bounds: 2 < p < 2**61", p);

    // find the constant and exponent such that p == c * 2**depth + 1
    const ulong depth = flint_ctz(p - UWORD(1));
    const ulong c = (p - UWORD(1)) >> depth;
    if (depth <= 3)
        flint_throw(FLINT_ERROR, "Provided prime p = %wu for n_fft does not satisfy `8 divides p-1`", p);

    // verify required length is not too large
    const ulong req_depth = FLINT_BIT_COUNT(len);
    if (depth < req_depth)
        flint_throw(FLINT_ERROR, "Required length %wu exceeds the limit for the modulus %wu", len, p);

    // find primitive root w of depth 2**depth
    const ulong prim_root = n_primitive_root_prime(p);
    const ulong w = n_powmod2(prim_root, c, p);

    // fill all attributes and tables
    n_fft_ctx_init2_root(F, w, depth, req_depth, p);
}

void n_fft_ctx_clear(n_fft_ctx_t F)
{
    _nmod_vec_clear(F->tab_w);
}




/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
