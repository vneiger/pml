#ifndef __NMOD_POLY_EXTRA__H
#define __NMOD_POLY_EXTRA__H

#include <flint/nmod_poly.h>
#include <flint/fft_small.h>
#include <flint/machine_vectors.h>

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
/* a structure for geometric evaluation / interpolation       */
/* TODO: handle non-existence                                 */
/*------------------------------------------------------------*/
typedef struct
{
    nn_ptr x, t, w, y, z;       // five vectors of precomputed constants
    nmod_poly_t f, g1, g2;      // three precomputed polys
    nmod_t mod;
    slong d;                    // number of points

} nmod_geometric_progression_struct;

typedef nmod_geometric_progression_struct nmod_geometric_progression_t[1];

/*------------------------------------------------------------*/
/* initializes all quantities attached to G                   */
/* evaluates/interpolates at powers of q = r^2                */
/*------------------------------------------------------------*/
void nmod_geometric_progression_init_set(nmod_geometric_progression_t G, ulong r, slong n, nmod_t mod);
 
/*------------------------------------------------------------*/
/* frees all memory attached to G                             */
/*------------------------------------------------------------*/
void nmod_geometric_progression_clear(nmod_geometric_progression_t G);

/*------------------------------------------------------------*/
/* in: polynomial to evaluate, deg(poly) < d                  */
/* out: v[i] = poly(q^i), i = 0 .. d-1                        */
/*------------------------------------------------------------*/
void nmod_geometric_progression_evaluate(nn_ptr v, const nmod_poly_t poly, const nmod_geometric_progression_t G);

/*------------------------------------------------------------*/
/* in: coeffs must have size at least d                       */
/* out: interpolating polynomial C s.t. C(q^i) = v[i], i<d    */
/*------------------------------------------------------------*/
void nmod_geometric_progression_interpolate(nmod_poly_t poly, nn_srcptr v, const nmod_geometric_progression_t G);




#ifdef __cplusplus
}
#endif

#endif // __NMOD_POLY_EXTRA__H

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
