#ifndef __NMOD_POLY_FFT__H
#define __NMOD_POLY_FFT__H

#include <flint/nmod_types.h>

#ifdef __cplusplus
extern "C" {
#endif

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* struct for basic info + precomputed root tables            */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

typedef struct
{
    nmod_t mod;
    ulong modn2;     // 2*mod.n  FIXME try not storing this
    ulong modn4;     // 4*mod.n  FIXME try not storing this
    ulong I;                   // sqrt(-1)
    ulong Ipre;                // precomp on I
    ulong J;                   // sqrt(I)
    ulong Jpre;                // precomp on J
    ulong IJ;                  // sqrt(-I) == I*J
    ulong IJpre;               // precomp on IJ
    ulong order;               // maximum supported order (currently: order of w + max order in precomputed tables)
    ulong w;                   // primitive (2**order)th root of 1
    ulong inv_w;               // inverse of w
    ulong * tab_w[64];            // tabulated powers of w
    ulong * tab_w_pre[64];        // tabulated powers of precomputations for multiplication by w mod mod.n
} nmod_fft_ctx_struct;
typedef nmod_fft_ctx_struct nmod_fft_ctx_t[1];

// FFT tables of powers / twiddle factors, say for w:
// for 0 <= ell < order-1, "powers_w[ell]" has length 2**(ell+1)
// and contains [w**(2**(order-2-ell)*i), 0 <= i < 2**(ell+1)]
// -> powers_w[0] = [1, I] where I is w**(2**(order-2))
//         (note I**2 = -1, I**3 = -I)
// -> powers_w[1] = [1, J, J**2, J**3] where J is w**(2**(order-3))
//         (note [J**4, J**5, J**6, J**7] = [-1, -J, -J**2, -J**3])
// -> etc..
// -> powers_w[order-2] = [1, w, w**2, ..., w**(2**(order-1)-1)]
//         (note next powers until 2**order-1 would be -1, -w, -w**2, ..)

/*------------------------------------------------------------*/
/* initializes all entries of F                               */
/* w primitive and w^(2^order))=1                             */
/* DFTs of size up to 2^order are supported                   */ 
/* order >= 3 and order < FLINT_BITS (-sth?) required         */
/*------------------------------------------------------------*/

// TODO :
// - optimize speed (see integer_fft_init)
// - allow initialization with NULL tables
// - allow fit_depth to precompute more tables when wanted/needed
// - separate computation of the tables from basic init
// - is storing I, J, IJ useful for efficiency?
// - allocate first tables on stack --> no noticeable effect
void nmod_fft_ctx_init_set(nmod_fft_ctx_t F, ulong w, ulong order, nmod_t mod);
void nmod_fft_ctx_init_set_new(nmod_fft_ctx_t F, ulong w, ulong order, nmod_t mod);

// version with just a list of roots in bit reversed order
void nmod_fft_ctx_init_set_red(nmod_fft_ctx_t F, ulong w, ulong order, nmod_t mod);




/*------------------------------------------------------------*/
/* clears all memory assigned to F                            */
/*------------------------------------------------------------*/
void nmod_fft_ctx_clear(nmod_fft_ctx_t F);
void nmod_fft_ctx_clear_red(nmod_fft_ctx_t F);



/*------------------------------------------------------------*/
/* fft evaluation, in place                                   */
/* x[i] = poly(w^i), len=2^k, in bit reverse order            */
/* x must have length >= len                                  */
/*------------------------------------------------------------*/

// TODO
//   - clarify input/output bounds below
//   - what constraints on modulus to avoid overflow?
//   - truncated: handle itrunc / otrunc
//   - try DIT variant (and Harvey's DIT tft)
//   - 

// recursive, decimation in frequency, radix 2
// input coefficients in [0..??)
// output coefficients in [0..??)
void _nmod_fft_dif_rec2_lazy(nn_ptr p, ulong len, ulong order, nmod_fft_ctx_t F);

// iterative, decimation in frequency, radix 2
// input coefficients in [0..??)
// output coefficients in [0..??)
void _nmod_fft_dif_iter2_lazy(nn_ptr p, ulong len, ulong order, nmod_fft_ctx_t F);

// recursive, reduction-tree approach, radix 2
// input coefficients in [0..??)
// output coefficients in [0..??)

// general: modulus is x**len - w**2 = (x**(len/2) - w) (x**(len/2) + w),
// where w = F->tab_w[0][node]
// e.g. for node == 1 this is w = I and x**len + 1
void _nmod_fft_red_rec2_lazy_general(nn_ptr p, ulong len, ulong order, ulong node, nmod_fft_ctx_t F);
// entry point: case where node == 0, modulus is x**len - 1 = (x**(len/2) - 1) (x**(len/2) + 1)
void _nmod_fft_red_rec2_lazy(nn_ptr p, ulong len, ulong order, nmod_fft_ctx_t F);

void _nmod_fft_red_iter2_lazy(nn_ptr p, ulong len, ulong order, nmod_fft_ctx_t F);

// recursive, decimation in frequency, radix 4
// input coefficients in [0..??)
// output coefficients in [0..??)
void _nmod_fft_dif_rec4_lazy(nn_ptr p, ulong len, ulong order, nmod_fft_ctx_t F);

#ifdef __cplusplus
}
#endif

#endif // __NMOD_POLY_FFT__H

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
