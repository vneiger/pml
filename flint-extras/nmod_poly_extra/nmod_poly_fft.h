#ifndef __NMOD_POLY_FFT__H
#define __NMOD_POLY_FFT__H

//#include "nmod_mat_poly.h"
#include <flint/nmod_types.h>

#ifdef __cplusplus
extern "C" {
#endif

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* integer FFT                                                */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

typedef struct
{
    nmod_t mod;
    ulong order;
    mp_limb_t w;                   // primitive (2**order)th root of 1
    mp_limb_t inv_w;               // inverse of w
    ulong ** tab_w;              // tabulated powers of w
    ulong ** tab_w_pre;              // tabulated powers of precomputations for multiplication by w mod mod.n
    //mp_limb_t ** tab_inv_w;      // tabulated powers of 1/w
    //mp_limb_t ** tab_inv_w_over_2; // length order, level k is [1/w{k+1}^i/2^k 
    //mp_limb_t * powers_inv_2;    // length order+1, with powers_inv_2[i] = 1/2^i 
} nmod_integer_fft_struct;
typedef nmod_integer_fft_struct nmod_integer_fft_t[1];

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
void nmod_integer_fft_init_set(nmod_integer_fft_t F, mp_limb_t w, ulong order, nmod_t mod);
void nmod_integer_fft_init_set2(nmod_integer_fft_t F, mp_limb_t w, ulong order, nmod_t mod);
void nmod_integer_fft_init_set_pre(nmod_integer_fft_t F, mp_limb_t w, ulong order, nmod_t mod);

// allow initialization with NULL tables
// allow fit_depth to precompute more tables when wanted/needed
// separate computation of the tables from basic init

/*------------------------------------------------------------*/
/* clears all memory assigned to F                            */
/*------------------------------------------------------------*/
void nmod_integer_fft_clear(nmod_integer_fft_t F);
void nmod_integer_fft_clear_pre(nmod_integer_fft_t F);



/*------------------------------------------------------------*/
/* fft evaluation, in place                                   */
/* returns x[i] = poly(w^i), len=2^k, in bit reverse order    */
/* x must have length >= len                                  */
/*------------------------------------------------------------*/
void _nmod_poly_dif_inplace_radix2_rec(mp_ptr p, ulong len, ulong order, nmod_integer_fft_t F);
void _nmod_poly_dif_inplace_radix2_rec_v2(mp_ptr p, ulong len, ulong order, nmod_integer_fft_t F);
void _nmod_poly_dif_inplace_radix2_rec_v3(mp_ptr p, ulong len, ulong order, nmod_integer_fft_t F);
void _nmod_poly_dif_inplace_radix2_rec_v4(mp_ptr p, ulong len, ulong order, nmod_integer_fft_t F);

void _nmod_poly_dif_inplace_radix4_rec(mp_ptr p, ulong len, ulong order, nmod_integer_fft_t F);
// TODO :
// - radix 4
// - inverse fft
// - not in place (e.g. if evaluating polynomial for fft mul, but polynomial needs to be kept)



#ifdef __cplusplus
}
#endif

#endif // __NMOD_POLY_FFT__H

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
