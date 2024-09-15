#ifndef __NMOD_POLY_FFT__H
#define __NMOD_POLY_FFT__H

#include "flint/flint.h"
#include "flint/long_extras.h"

#define N_FFT_OLD_MAX_TAB_SIZE 32

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
    ulong mod;
    ulong mod2;                // 2*mod  (storing this does help for speed, in _dif_rec2 at least)
    ulong mod4;                // 4*mod  (storing this does help for speed, in _dif_rec2 at least)
    ulong I;                   // sqrt(-1)
    ulong Ipre;                // precomp on I
    ulong J;                   // sqrt(I)
    ulong Jpre;                // precomp on J
    ulong IJ;                  // sqrt(-I) == I*J
    ulong IJpre;               // precomp on IJ
    ulong order;               // maximum supported order (currently: order of w)
    ulong * tab_w;             // tabulated powers of w
} n_fft_ctx_struct;
typedef n_fft_ctx_struct n_fft_ctx_t[1];

// TODO describe fields

void n_fft_ctx_init_set(n_fft_ctx_t F, ulong w, ulong order, ulong mod);
void n_fft_ctx_clear(n_fft_ctx_t F);



/*------------------------------------------------------------*/
/* fft evaluation, in place                                   */
/* x[i] = poly(w^i), len=2^k, in bit reverse order            */
/* x must have length >= len                                  */
/*------------------------------------------------------------*/

// TODO
//   - clarify input/output bounds below
//   - what constraints on modulus to avoid overflow?
//   - truncated: handle itrunc / otrunc

// recursive, reduction-tree approach, radix 2
// input coefficients in [0..??)
// output coefficients in [0..??)

// general: modulus is x**len - w**2 = (x**(len/2) - w) (x**(len/2) + w),
// where w = F->tab_w[node]
// e.g. for node == 1 this is w = I and x**len + 1
void _n_fft_red_rec2_lazy_general(nn_ptr p, ulong len, ulong order, ulong node, n_fft_ctx_t F);
// entry point: case where node == 0, modulus is x**len - 1 = (x**(len/2) - 1) (x**(len/2) + 1)
void _n_fft_red_rec2_lazy(nn_ptr p, ulong len, ulong order, n_fft_ctx_t F);

void _n_fft_red_rec4_lazy_general(nn_ptr p, ulong len, ulong order, ulong node, n_fft_ctx_t F);
void _n_fft_red_rec4_lazy(nn_ptr p, ulong len, ulong order, n_fft_ctx_t F);

void _n_fft_red_iter2_lazy(nn_ptr p, ulong len, ulong order, n_fft_ctx_t F);




/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* CURRENTLY UNUSED                                           */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/


typedef struct
{
    ulong mod;
    ulong mod2;                // 2*mod  (storing this does help for speed, in _dif_rec2 at least)
    ulong mod4;                // 4*mod  (storing this does help for speed, in _dif_rec2 at least)
    ulong I;                   // sqrt(-1)
    ulong Ipre;                // precomp on I
    ulong J;                   // sqrt(I)
    ulong Jpre;                // precomp on J
    ulong IJ;                  // sqrt(-I) == I*J
    ulong IJpre;               // precomp on IJ
    ulong order;               // maximum supported order (currently: order of w)
    ulong w;                   // primitive (2**order)th root of 1
    ulong winv;               // inverse of w
    ulong * tab_w[N_FFT_OLD_MAX_TAB_SIZE]; // tabulated powers of w and winv, see below
    //ulong * tab_winv[n_fft_MAX_TAB_SIZE]; // tabulated powers of winv
} n_fft_old_ctx_struct;
typedef n_fft_old_ctx_struct n_fft_old_ctx_t[1];

// FFT tables of powers / twiddle factors, say for w of order == `2**order`:
// for 0 <= ell <= order-4, "powers_w[ell]" has length 2**(ell+5)
// and contains [wi, wi_pre for 0 <= i < 2**(ell+4)]
//    where wi = w**(2**(order-4-ell)*i) are the 2**(ell+4)-th roots of unity
//      and wi_pre is the corresponding precomputation for mulmod
// -> powers_w[-2] would have been [1, 1_pre, I, I_pre, -1, -1_pre, -I, -I_pre]
//          where I is w**(2**(order-2))
// -> powers_w[-1] would have been similar with [1, J, J**2, J**3, -1, -J, -J**2, -J**3]
//          where J is w**(2**(order-3))
//         (note J**2 == I, and J**3 == I*J)
// -> etc..
// -> powers_w[order-4] = [1, 1_pre, w, w_pre, w**2, w**2_pre, ..., w**(2**order -1), w**(2**order -1)_pre]
//
// Observe that the second part of each table gives the inverses:
//    -I is the inverse of I; 
//    -J**3 is the inverse of J, -J**2 is the inverse of J**2, -J is the inverse of J**3;
//    etc.


/*------------------------------------------------------------*/
/* initializes all entries of F                               */
/* w primitive and w^(2^order))=1                             */
/* DFTs of size up to 2^order are supported                   */ 
/* order >= 3 and order < FLINT_BITS (-sth?) required         */
/*------------------------------------------------------------*/

// fills seq with, for i = 0...d-1:
// seq[2*i] = a**i
// seq[2*i+1] = precomp(a**i)   (precomputed quotient for n_mulmod_shoup)
// seq already allocated with >= 2*d entries; seq[i] not accessed for i >= 2*d

// this version assumes d >= 5
void _n_geometric_sequence_with_precomp(ulong * seq, ulong a, ulong d, ulong n);
// general, any d
void n_geometric_sequence_with_precomp(ulong * seq, ulong a, ulong d, ulong n);

// rather fft specific:
// as above (and restricted to d >= 5),
// for i < d, seq[2*i] and seq[2*i+1] are as above
// but also adds similar stored data for -a**i:
// for d <= i < 2*d, seq[2*i] = -a**(i-d) and seq[2*i+1] = precomp(-a**(i-d))
// (if w is a d-th root of unity, these are inverses of the first part of the array)
// seq already allocated with >= 4*d entries; seq[i] not accessed for i >= 4*d
void _n_geometric_sequence_and_opposites_with_precomp(ulong * seq, ulong a, ulong d, ulong n);

// NOTES:
// - not thoroughly optimized (especially copying part; and try building by increasing ell?)
// - separate computation of the tables from basic init,
//   and add fit_depth to precompute more tables when wanted/needed
void n_fft_old_ctx_init_set(n_fft_old_ctx_t F, ulong w, ulong order, ulong mod);
void n_fft_old_ctx_clear(n_fft_old_ctx_t F);



// recursive, Gentleman-Sande butterflies, radix 2
// input coefficients in [0..??)
// output coefficients in [0..??)
void _n_fft_old_dif_rec2_lazy(nn_ptr p, ulong len, ulong order, n_fft_old_ctx_t F);

// recursive, Gentleman-Sande butterflies, radix 4
// input coefficients in [0..??)
// output coefficients in [0..??)
void _n_fft_old_dif_rec4_lazy(nn_ptr p, ulong len, ulong order, n_fft_old_ctx_t F);

// recursive, Gentleman-Sande butterflies, radix 8
// input coefficients in [0..??)
// output coefficients in [0..??)
void _n_fft_old_dif_rec8_lazy(nn_ptr p, ulong len, ulong order, n_fft_old_ctx_t F);

// iterative, Gentleman-Sande butterflies, radix 2
// input coefficients in [0..??)
// output coefficients in [0..??)
void _n_fft_old_dif_iter2_lazy(nn_ptr p, ulong len, ulong order, n_fft_old_ctx_t F);


#ifdef __cplusplus
}
#endif

#endif // __NMOD_POLY_FFT__H

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
