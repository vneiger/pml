#ifndef __NMOD_POLY_EXTRA__H
#define __NMOD_POLY_EXTRA__H

#include <flint/nmod_poly.h>
#include "nmod_extra.h"

#ifdef __cplusplus
extern "C" {
#endif


/** Generates random polynomial `pol` of length up to `len` with uniformly
 * random coefficients. */
void nmod_poly_rand(nmod_poly_t pol,
                    flint_rand_t state,
                    slong len);



/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* a structure for FFT modulo p                               */
/* only uses flint arithmetic                                 */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
typedef struct
{
    nmod_t mod;
    mp_limb_t w, inv_w;       // root of 1 and its inverse
    ulong order;       // its order
    
    /*  table of roots of unity, with (order+1) levels                              */
    /*  at level k > 0, set wk = w^(2^(order-k)). this is a K-th root, with K=2^k   */
    /*  the entries of powers_w[k] are                                              */
    /*        1, wk, wk^2, ..., wk^{K/2-1},           (len K/2)                     */
    /*        1, wk^2, wk^4, ..., wk^{K/2-2},         (len K/4)                     */
    /*        1, wk^4, wk^8, ..., wk^{K/2-4}          (len K/8)                     */
    /*        ...                                                                   */
    /*        1, wk^{K/8}, wk^{K/4}, wk^{3K/8},       (len (K/2)/(K/8)=4)           */
    /*        1, wk^{K/4},                            (len 2)                       */
    /*        1.                                      (len 1)                       */
    /*  total length K-1 (we allocate K)                                            */ 
    mp_ptr * powers_w;
    mp_ptr * powers_inv_w_t;  // same table for 1/w as for w

    // length order+1, with powers_inv_2[i] = 1/2^i 
    mp_ptr powers_inv_2;
    // length order, level k has length 2^k with entries 1/w{k+1}^i/2^k 
    mp_ptr * powers_inv_w_over_2;
    
} nmod_fft_struct;
typedef nmod_fft_struct nmod_fft_t[1];

/*------------------------------------------------------------*/
/* a butterfly with reduction mod p                           */
/*------------------------------------------------------------*/
#define BUTTERFLY(o1, o2, i1, i2, mod)                 \
    do {                                               \
        o1 = i1 + i2;                                  \
        o2 = i1 - i2;                                  \
        CORRECT_0_2P(o1, mod.n);                       \
        CORRECT_MINUSP_P(o2, mod.n);                   \
    } while(0)

/*------------------------------------------------------------*/
/* initializes all entries of F                               */
/* w primitive and w^(2^order))=1                             */
/* DFTs of size up to 2^order are supported                   */ 
/*------------------------------------------------------------*/
void nmod_fft_init_set(nmod_fft_t F, mp_limb_t w, ulong order, nmod_t mod);

/*------------------------------------------------------------*/
/* clears all memory assigned to F                            */
/*------------------------------------------------------------*/
void nmod_fft_clear(nmod_fft_t F);

/*------------------------------------------------------------*/
/* fft evaluation                                             */
/* returns x[i] = poly(w^bitreverse_n(i)), n=2^k              */
/* x must have length >= n                                    */
/*------------------------------------------------------------*/
void nmod_fft_evaluate(mp_ptr x, const nmod_poly_t poly, const nmod_fft_t F, const ulong k);

/*------------------------------------------------------------*/
/* inverse fft                                                */
/* given x[i] = poly(w^bitreverse(i,2^k)), returns poly       */
/*------------------------------------------------------------*/
void nmod_fft_interpolate(nmod_poly_t poly, mp_srcptr x, const nmod_fft_t F, const ulong k);

/*------------------------------------------------------------*/
/* tft evaluation                                             */
/* returns a cyclotomic tft evaluation of x at n points       */
/* x must have length >= N                                    */
/*------------------------------------------------------------*/
void nmod_tft_evaluate(mp_ptr x, const nmod_poly_t poly, const nmod_fft_t F, const ulong N);

/*------------------------------------------------------------*/
/* inverse tft                                                */
/*------------------------------------------------------------*/
void nmod_tft_interpolate(nmod_poly_t poly, mp_srcptr x, const nmod_fft_t F, const ulong N);

/*-----------------------------------------------------------*/
/* transpose tft in length N                                 */
/*-----------------------------------------------------------*/
void nmod_tft_evaluate_t(mp_ptr x, mp_srcptr A, const nmod_fft_t F, const ulong N);

/*-----------------------------------------------------------*/
/* transpose inverse tft in length N                         */
/*-----------------------------------------------------------*/
void nmod_tft_interpolate_t(mp_ptr x, mp_srcptr A, const nmod_fft_t F, const ulong N);


    
#ifdef HAS_INT128

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* a structure for FFT modulo p                               */
/* assumes C compiler has __int128 built-in type              */
/* __int128 used for mod p reductions                         */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
typedef struct
{
    nmod_t mod;
    mp_limb_t w, inv_w;       // root of 1 and its inverse
    ulong order;       // its order
    
    mp_ptr * powers_w; // the powers of the successive roots, as in nmod_fft_t
    mp_ptr * i_powers_w; // the same, after precomputation

    mp_ptr * powers_inv_w_t;  // same table for 1/w as for w
    mp_ptr * i_powers_inv_w_t;  // same table for 1/w as for w
    
    mp_ptr * powers_r; // the powers of rho such that rho^2 = w
    mp_ptr * i_powers_r; // the same, after precomputation
    mp_ptr * powers_inv_r; // the powers of rho such that rho^2 = w
    mp_ptr * i_powers_inv_r; // the same, after precomputation
    
    // inverses of the powers of 2
    mp_ptr powers_inv_2;
    mp_ptr i_powers_inv_2;

    // length order, level k has length 2^k with entries 1/w{k+1}^i/2^k 
    mp_ptr * powers_inv_w_over_2;
    mp_ptr * i_powers_inv_w_over_2;
} nmod_64_fft_struct;
typedef nmod_64_fft_struct nmod_64_fft_t[1];


/*------------------------------------------------------------*/
/* initializes all entries of F                               */
/* w primitive and w^(2^order))=1                             */
/* DFTs of size up to 2^order are supported                   */ 
/*------------------------------------------------------------*/
void nmod_64_fft_init_set(nmod_64_fft_t F, mp_limb_t w, ulong order, nmod_t mod);

/*------------------------------------------------------------*/
/* clears all memory assigned to F                            */
/*------------------------------------------------------------*/
void nmod_64_fft_clear(nmod_64_fft_t F);

/*------------------------------------------------------------*/
/* fft evaluation                                             */
/* returns x[i] = poly(w^bitreverse_n(i)), n=2^k              */
/* x must have length >= n                                    */
/*------------------------------------------------------------*/
void nmod_64_fft_evaluate(mp_ptr x, const nmod_poly_t poly, const nmod_64_fft_t F, const ulong k);

/*------------------------------------------------------------*/
/* inverse fft                                                */
/* given x[i] = poly(w^bitreverse(i,2^k)), returns poly       */
/*------------------------------------------------------------*/
void nmod_64_fft_interpolate(nmod_poly_t poly, mp_srcptr x, const nmod_64_fft_t F, const ulong k);

/*------------------------------------------------------------*/
/* tft evaluation                                             */
/* returns a cyclotomic tft evaluation of x at n points       */
/* x must have length >= N                                    */
/*------------------------------------------------------------*/
void nmod_64_tft_evaluate(mp_ptr x, const nmod_poly_t poly, const nmod_64_fft_t F, const ulong N);

/*------------------------------------------------------------*/
/* tft interpolation                                          */
/*------------------------------------------------------------*/
void nmod_64_tft_interpolate(nmod_poly_t poly, mp_srcptr x, const nmod_64_fft_t F, const ulong N);

/*-----------------------------------------------------------*/
/* transpose tft in length N                                 */
/*-----------------------------------------------------------*/
void nmod_64_tft_evaluate_t(mp_ptr x, mp_srcptr A, const nmod_64_fft_t F, const ulong N);

/*-----------------------------------------------------------*/
/* transpose inverse tft in length N                         */
/*-----------------------------------------------------------*/
void nmod_64_tft_interpolate_t(mp_ptr x, mp_srcptr A, const nmod_64_fft_t F, const ulong N);

#endif




#ifdef HAS_AVX2

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* a structure for FFT modulo p                               */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
typedef struct
{
    nmod_t mod;
    mp_hlimb_t w, inv_w;       // root of 1 and its inverse
    ulong order;       // its order
    
    mp_hlimb_t ** powers_w; // the powers of the successive roots, ordered for fft_k
    mp_hlimb_t ** i_powers_w; // the same, after precomputation

    mp_hlimb_t ** powers_inv_w_t; // the inverse powers of the successive roots, ordered for fft_k
    mp_hlimb_t ** i_powers_inv_w_t; // the same, after precomputation

    // inverses of the powers of 2
    mp_hlimb_t * powers_inv_2;
    mp_hlimb_t * i_powers_inv_2;

    // length order, level k has length 2^k with entries 1/w{k+1}^i/2^k 
    mp_hlimb_t ** powers_inv_w_over_2;
    mp_hlimb_t ** i_powers_inv_w_over_2;

} nmod_32_fft_struct;
typedef nmod_32_fft_struct nmod_32_fft_t[1];

/*------------------------------------------------------------*/
/* a butterfly with reduction mod p                           */
/*------------------------------------------------------------*/
#define BUTTERFLY32(o1, o2, i1, i2, mod)                    \
    do {                                                    \
        o1 = i1 + i2;                                       \
        o2 = i1 - i2;                                       \
        CORRECT_32_0_2P(o1, (mp_hlimb_t) mod.n);            \
        CORRECT_32_MINUSP_P(o2, (mp_hlimb_t) mod.n);        \
    } while(0)

/*------------------------------------------------------------*/
/* initializes all entries of F                               */
/* w primitive and w^(2^order))=1                             */
/* DFTs of size up to 2^order are supported                   */ 
/*------------------------------------------------------------*/
void nmod_32_fft_init_set(nmod_32_fft_t F, mp_limb_t w, ulong order, nmod_t mod);

/*------------------------------------------------------------*/
/* clears all memory assigned to F                            */
/*------------------------------------------------------------*/
void nmod_32_fft_clear(nmod_32_fft_t F);

/*------------------------------------------------------------*/
/* fft evaluation                                             */
/* returns x[i] = poly(w^bitreverse_n(i)), n=2^k              */
/* x must have length >= n                                    */
/*------------------------------------------------------------*/
void nmod_avx2_32_fft_evaluate(mp_ptr x, const nmod_poly_t poly, const nmod_32_fft_t F, const ulong k);

/*------------------------------------------------------------*/
/* inverse fft                                                */
/* given x[i] = poly(w^bitreverse(i,2^k)), returns poly       */
/*------------------------------------------------------------*/
void nmod_avx2_32_fft_interpolate(nmod_poly_t poly, mp_srcptr x, const nmod_32_fft_t F, const ulong k);

/*------------------------------------------------------------*/
/* tft evaluation                                             */
/* returns a cyclotomic tft evaluation of x at n points       */
/* x must have length >= N                                    */
/*------------------------------------------------------------*/
void nmod_avx2_32_tft_evaluate(mp_ptr x, const nmod_poly_t poly, const nmod_32_fft_t F, const ulong N);

/*------------------------------------------------------------*/
/* inverse tft                                                */
/*------------------------------------------------------------*/
void nmod_avx2_32_tft_interpolate(nmod_poly_t poly, mp_srcptr x, const nmod_32_fft_t F, const ulong N);

/*-----------------------------------------------------------*/
/* transpose tft in length N                                 */
/*-----------------------------------------------------------*/
void nmod_avx2_32_tft_evaluate_t(mp_ptr x, mp_srcptr A, const nmod_32_fft_t F, const ulong N);

/*-----------------------------------------------------------*/
/* transpose inverse tft in length N                         */
/*-----------------------------------------------------------*/
void nmod_avx2_32_tft_interpolate_t(mp_ptr x, mp_srcptr A, const nmod_32_fft_t F, const ulong N);

#endif

 

#ifdef HAS_AVX512

/*------------------------------------------------------------*/
/* fft evaluation                                             */
/* returns x[i] = poly(w^bitreverse_n(i)), n=2^k              */
/* x must have length >= n                                    */
/*------------------------------------------------------------*/
void nmod_avx512_32_fft_evaluate(mp_ptr x, const nmod_poly_t poly, const nmod_32_fft_t F, const ulong k);

/*------------------------------------------------------------*/
/* inverse fft                                                */
/* given x[i] = poly(w^bitreverse(i,2^k)), returns poly       */
/*------------------------------------------------------------*/
void nmod_avx512_32_fft_interpolate(nmod_poly_t poly, mp_srcptr x, const nmod_32_fft_t F, const ulong k);


#endif

/*------------------------------------------------------------*/
/* a structure for geometric evaluation / interpolation       */
/* TODO: handle non-existence                                 */
/*------------------------------------------------------------*/
typedef struct
{
    mp_ptr x, t, w, y, z;       // five vectors of precomputed constants
    nmod_poly_t f, g1, g2;      // three precomputed polys
    nmod_t mod;
    slong d;                    // number of points

} nmod_geometric_progression_struct;

typedef nmod_geometric_progression_struct nmod_geometric_progression_t[1];

/*------------------------------------------------------------*/
/* initializes all quantities attached to G                   */
/* evaluates/interpolates at powers of q = r^2                */
/*------------------------------------------------------------*/
void nmod_geometric_progression_init_set(nmod_geometric_progression_t G, mp_limb_t r, slong n, nmod_t mod);
 
/*------------------------------------------------------------*/
/* frees all memory attached to G                             */
/*------------------------------------------------------------*/
void nmod_geometric_progression_clear(nmod_geometric_progression_t G);

/*------------------------------------------------------------*/
/* in: polynomial to evaluate, deg(poly) < d                  */
/* out: v[i] = poly(q^i), i = 0 .. d-1                        */
/*------------------------------------------------------------*/
void nmod_geometric_progression_evaluate(mp_ptr v, const nmod_poly_t poly, const nmod_geometric_progression_t G);

/*------------------------------------------------------------*/
/* in: coeffs must have size at least d                       */
/* out: interpolating polynomial C s.t. C(q^i) = v[i], i<d    */
/*------------------------------------------------------------*/
void nmod_geometric_progression_interpolate(nmod_poly_t poly, mp_srcptr v, const nmod_geometric_progression_t G);

#ifdef __cplusplus
}
#endif

#endif // __NMOD_POLY_EXTRA__H

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
