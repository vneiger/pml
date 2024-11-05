#ifndef __NMOD_POLY_INTEGER_FFT__H
#define __NMOD_POLY_INTEGER_FFT__H

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
    ulong modn2;     // 2*mod.n
    ulong modn4;     // 4*mod.n
    ulong I;                   // sqrt(-1)
    ulong Ipre;                // precomp on I
    ulong J;                   // sqrt(I)
    ulong Jpre;                // precomp on J
    ulong IJ;                  // sqrt(-I) == I*J
    ulong IJpre;               // precomp on IJ
    ulong order;               // maximum supported order (currently: order of w + max order in precomputed tables)
    ulong w;                   // primitive (2**order)th root of 1
    ulong inv_w;               // inverse of w
    ulong ** tab_w;            // tabulated powers of w
    ulong ** tab_w_pre;        // tabulated powers of precomputations for multiplication by w mod mod.n
    //ulong ** tab_inv_w;      // tabulated powers of 1/w
    //ulong ** tab_inv_w_over_2; // length order, level k is [1/w{k+1}^i/2^k 
    //ulong * powers_inv_2;    // length order+1, with powers_inv_2[i] = 1/2^i 
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
void nmod_integer_fft_init_set(nmod_integer_fft_t F, ulong w, ulong order, nmod_t mod);
void nmod_integer_fft_init_set2(nmod_integer_fft_t F, ulong w, ulong order, nmod_t mod); // faster

// not tried to optimize:
// version with table of precomputed things for Shoup multiplication
void nmod_integer_fft_init_set_pre(nmod_integer_fft_t F, ulong w, ulong order, nmod_t mod);

// version with just a list of roots in bit reversed order
void nmod_integer_fft_init_set_red(nmod_integer_fft_t F, ulong w, ulong order, nmod_t mod);
void nmod_integer_fft_init_set_red_pre(nmod_integer_fft_t F, ulong w, ulong order, nmod_t mod);

// allow initialization with NULL tables
// allow fit_depth to precompute more tables when wanted/needed
// separate computation of the tables from basic init
// allocate first tables on stack



/*------------------------------------------------------------*/
/* clears all memory assigned to F                            */
/*------------------------------------------------------------*/
void nmod_integer_fft_clear(nmod_integer_fft_t F);
void nmod_integer_fft_clear_pre(nmod_integer_fft_t F);
void nmod_integer_fft_clear_red(nmod_integer_fft_t F);
void nmod_integer_fft_clear_red_pre(nmod_integer_fft_t F);



/*------------------------------------------------------------*/
/* fft evaluation, in place                                   */
/* returns x[i] = poly(w^i), len=2^k, in bit reverse order    */
/* x must have length >= len                                  */
/*------------------------------------------------------------*/
void _nmod_poly_dif_inplace_radix2_rec_prenorm(nn_ptr p, ulong len, ulong order, nmod_integer_fft_t F);
void _nmod_poly_dif_inplace_radix2_rec_prenorm_unroll4(nn_ptr p, ulong len, ulong order, nmod_integer_fft_t F);
void _nmod_poly_dif_inplace_radix2_rec_shoup(nn_ptr p, ulong len, ulong order, nmod_integer_fft_t F);
void _nmod_poly_dif_inplace_radix2_rec_shoup_unroll4(nn_ptr p, ulong len, ulong order, nmod_integer_fft_t F);
void _nmod_poly_dif_inplace_radix2_rec_bench(nn_ptr p, ulong len, ulong order, nmod_integer_fft_t F);

void _nmod_poly_dif_inplace_radix2_iter_prenorm(nn_ptr p, ulong len, ulong order, nmod_integer_fft_t F);
void _nmod_poly_dif_inplace_radix2_iter_shoup(nn_ptr p, ulong len, ulong order, nmod_integer_fft_t F);
void _nmod_poly_dif_inplace_radix2_iter_bench(nn_ptr p, ulong len, ulong order, nmod_integer_fft_t F);

// not very clean
void _nmod_poly_dif_inplace_radix4_rec(nn_ptr p, ulong len, ulong order, nmod_integer_fft_t F);
void _nmod_poly_dif_inplace_radix4_rec_bench(nn_ptr p, ulong len, ulong order, nmod_integer_fft_t F);
void _nmod_poly_dif_inplace_radix4_iter(nn_ptr p, ulong len, ulong order, nmod_integer_fft_t F);


// reduction tree attempt
void _nmod_poly_red_inplace_radix2_rec_prenorm(nn_ptr p, ulong len, ulong order, ulong node, nmod_integer_fft_t F);
void _nmod_poly_red_inplace_radix2_rec_shoup(nn_ptr p, ulong len, ulong order, ulong node, nmod_integer_fft_t F);


// lazy attempts
void _nmod_poly_dif_inplace_radix2_rec_shoup_lazy(nn_ptr p, ulong len, ulong order, nmod_integer_fft_t F);
void _nmod_poly_dif_inplace_radix2_iter_shoup_lazy(nn_ptr p, ulong len, ulong order, nmod_integer_fft_t F);
void _nmod_poly_red_inplace_radix2_rec_shoup_lazy(nn_ptr p, ulong len, ulong order, ulong node, nmod_integer_fft_t F);
void _nmod_poly_dif_inplace_radix4_rec_shoup_lazy(nn_ptr p, ulong len, ulong order, nmod_integer_fft_t F);



#ifdef __cplusplus
}
#endif

#endif // __NMOD_POLY_INTEGER_FFT__H

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
