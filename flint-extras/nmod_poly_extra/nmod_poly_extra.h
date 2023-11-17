#ifndef __NMOD_POLY_EXTRA__H
#define __NMOD_POLY_EXTRA__H

#include <flint/nmod_poly.h>
#include <flint/fft_small.h>
#include <flint/machine_vectors.h>

#include <math.h>
#include "nmod_extra.h"

#ifdef __cplusplus
extern "C" {
#endif


/** Generates random polynomial `pol` of length up to `len` with uniformly
 * random coefficients. If `len` is nonpositive, `pol` is set to zero. */
void nmod_poly_rand(nmod_poly_t pol,
                    flint_rand_t state,
                    slong len);


/** Generates random monic polynomial `pol` of length exactly `len` with
 * uniformly random coefficients. If `len` is nonpositive, `pol` is set to
 * zero. */
void nmod_poly_rand_monic(nmod_poly_t pol,
                          flint_rand_t state,
                          slong len);


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* fft_small FFT, using doubles                               */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

typedef struct
{
    nmod_t mod;
    vec1d p, pinv;

    mp_limb_t w, inv_w;       // root of 1 and its inverse
    ulong order;       // its order
    
    vec1d **powers_w;
    mp_limb_t **powers_inv_w_t;  // same table for 1/w as for w

    // length order+1, with powers_inv_2[i] = 1/2^i 
    mp_limb_t *powers_inv_2;
    
    // length order, level k has length 2^k with entries 1/w{k+1}^i/2^k 
    vec1d **powers_inv_w_over_2;
} nmod_sd_fft_struct;
typedef nmod_sd_fft_struct nmod_sd_fft_t[1];


/*------------------------------------------------------------*/
/* initializes all entries of F                               */
/* w primitive and w^(2^order))=1                             */
/* DFTs of size up to 2^order are supported                   */ 
/*------------------------------------------------------------*/
void nmod_sd_fft_init_set(nmod_sd_fft_t F, mp_limb_t w, ulong order, nmod_t mod);
    
/*------------------------------------------------------------*/
/* clears all memory assigned to F                            */
/*------------------------------------------------------------*/
void nmod_sd_fft_clear(nmod_sd_fft_t F);


void sd_fft_ctx_init_inverse(sd_fft_ctx_t Qt, sd_fft_ctx_t Q);

/*------------------------------------------------------------*/
/* fft evaluation                                             */
/* returns x[i] = poly(w^i), n=2^k, up to permutation         */
/* x must have length >= n                                    */
/*------------------------------------------------------------*/
void nmod_sd_fft_evaluate(mp_ptr x, const nmod_poly_t poly, sd_fft_lctx_t Q, const ulong k);
void nmod_sd_fft_evaluate_t(mp_ptr x, const nmod_poly_t poly, sd_fft_lctx_t Q, const ulong k);

/*------------------------------------------------------------*/
/* tft evaluation and its transpose                           */
/* x must have length >= N                                    */
/*------------------------------------------------------------*/
void nmod_sd_tft_evaluate(mp_ptr x, const nmod_poly_t poly, sd_fft_lctx_t Q, nmod_sd_fft_t F, const ulong N);
void nmod_sd_tft_evaluate_t(mp_ptr x, mp_srcptr A, sd_fft_lctx_t Q, nmod_sd_fft_t F, ulong N);

/*------------------------------------------------------------*/
/* tft interpolation and its transpose                        */
/* inverts nmod_sd_tft_evaluate                               */
/*------------------------------------------------------------*/
void nmod_sd_tft_interpolate(nmod_poly_t poly, mp_ptr x, sd_fft_lctx_t Q, nmod_sd_fft_t F, const ulong N);
void nmod_sd_tft_interpolate_t(mp_ptr a, mp_srcptr A, sd_fft_lctx_t Q, nmod_sd_fft_t F, const ulong N);

/*------------------------------------------------------------*/
/* fft interpolation                                          */
/* inverts nmod_sd_fft_evaluate                               */
/*------------------------------------------------------------*/
void nmod_sd_fft_interpolate(nmod_poly_t poly, mp_ptr x, sd_fft_lctx_t Q, const ulong k); 
void nmod_sd_fft_interpolate_t(nmod_poly_t poly, mp_ptr x, sd_fft_lctx_t Q, const ulong k); 


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
