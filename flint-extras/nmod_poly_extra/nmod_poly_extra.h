#ifndef __NMOD_POLY_EXTRA__H
#define __NMOD_POLY_EXTRA__H

#include <flint/nmod_poly.h>
#include "nmod_extra.h"

#ifdef __cplusplus
extern "C" {
#endif


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

    // same table for 1/w, but the order of the rows is reversed
    mp_ptr * powers_inv_w;

    // length order+1, with powers_inv_2[i] = 1/2^i 
    mp_ptr powers_inv_2;
    
} nmod_plain_fft_struct;
typedef nmod_plain_fft_struct nmod_plain_fft_t[1];

/*------------------------------------------------------------*/
/* initializes all entries of F                               */
/* w primitive and w^(2^order))=1                             */
/* DFTs of size up to 2^order are supported                   */ 
/*------------------------------------------------------------*/
void nmod_plain_fft_init_set(nmod_plain_fft_t F, mp_limb_t w, ulong order, nmod_t mod);

/*------------------------------------------------------------*/
/* clears all memory assigned to F                            */
/*------------------------------------------------------------*/
void nmod_plain_fft_clear(nmod_plain_fft_t F);

#if 1 
#define BUTTERFLY(o1, o2, i1, i2, mod)                 \
    do {                                               \
        o1 = i1 + i2;                                  \
        o2 = i1 - i2;                                  \
        CORRECT_0_2P(o1, mod.n);                       \
        CORRECT_MINUSP_P(o2, mod.n);                   \
    } while(0)
#else
#define BUTTERFLY(o1, o2, i1, i2, mod)               \
    do {                                             \
        o1 = _nmod_add(i1, i2, mod);                 \
        o2 = _nmod_sub(i1, i2, mod);                 \
    } while(0)
#endif


/*------------------------------------------------------------*/
/*  in-place size 1 butterfly                                 */
/*  o1, o2 = i1+i2, i1-i2 mod p                               */
/*------------------------------------------------------------*/
static inline void _fft_1(mp_ptr x, const nmod_t mod)
{
    mp_limb_t u0, u1, t0, t1;
    
    u0 = x[0];
    u1 = x[1];
    
    t0 = u0 + u1;
    t1 = u0 - u1;
    CORRECT_0_2P(t0, mod.n);
    CORRECT_MINUSP_P(t1, mod.n);
    
    x[0] = t0;
    x[1] = t1;
}

/*------------------------------------------------------------*/
/*  in-place size 2 butterfly                                 */
/*  t1, t2, t3, t4 = i1+i3, i2+i4, i1-i3, w(i2-i4) mod p      */
/*  o1, o2, o3, o4 = t1+t2, t1-t2, t3+t4, t3-t4 mod p         */
/*------------------------------------------------------------*/
static inline void _fft_2(mp_ptr x, const mp_limb_t w, const nmod_t mod)
{
    mp_limb_t u0, u1, u2, u3, v0, v1, v2, v3, z0, z1, z2, z3;

    u0 = x[0];
    u1 = x[1];
    u2 = x[2];
    u3 = x[3];
    
    v0 = u0 + u2;
    CORRECT_0_2P(v0, mod.n);
    v2 = u0 - u2;
    CORRECT_MINUSP_P(v2, mod.n);
    v1 = u1 + u3;
    CORRECT_0_2P(v1, mod.n);
    v3 = u1 - u3;
    CORRECT_MINUSP_P(v3, mod.n);
    v3 = nmod_mul(v3, w, mod);
    z0 = v0 + v1;
    CORRECT_0_2P(z0, mod.n);
    z1 = v0 - v1;
    CORRECT_MINUSP_P(z1, mod.n);
    z2 = v2 + v3;
    CORRECT_0_2P(z2, mod.n);
    z3 = v2 - v3;
    CORRECT_MINUSP_P(z3, mod.n);
    
    x[0] = z0;
    x[1] = z1;
    x[2] = z2;
    x[3] = z3;
}

/*------------------------------------------------------------*/
/* fft evaluation                                             */
/* returns x[i] = poly(w^bitreverse_n(i)), n=2^k              */
/* x must have length >= n                                    */
/*------------------------------------------------------------*/
void nmod_plain_fft_evaluate(mp_ptr x, const nmod_poly_t poly, const nmod_plain_fft_t F, const ulong k);

/*------------------------------------------------------------*/
/*  in-place size 1 inverse butterfly                         */
/*------------------------------------------------------------*/
static inline void _inv_fft_1(mp_ptr x, const nmod_t mod)
{
    mp_limb_t u0, u1, t0, t1;
    
    u0 = x[0];
    u1 = x[1];
    BUTTERFLY(t0, t1, u0, u1, mod);
    x[0] = t0;
    x[1] = t1;
}

/*------------------------------------------------------------*/
/*  in-place size 2 inverse butterfly                         */
/*------------------------------------------------------------*/
static inline void _inv_fft_2(mp_ptr x, const mp_limb_t w, const nmod_t mod)
{
    mp_limb_t u0, u1, u2, u3, v0, v1, v2, v3, z0, z1, z2, z3;
    
    u0 = x[0];
    u1 = x[1];
    u2 = x[2];
    u3 = x[3];

    BUTTERFLY(v0, v1, u0, u1, mod);
    BUTTERFLY(v2, v3, u2, u3, mod);
    BUTTERFLY(z0, z1, v0, v2, mod);
    v3 = nmod_mul(v3, w, mod);
    BUTTERFLY(z2, z3, v1, v3, mod);

    x[0] = z0;
    x[1] = z2;
    x[2] = z1;
    x[3] = z3;
}

/*------------------------------------------------------------*/
/* inverse fft                                                */
/* given x[i] = poly(w^bitreverse(i,2^k)), returns poly       */
/*------------------------------------------------------------*/
void nmod_plain_fft_interpolate(nmod_poly_t poly, mp_srcptr x, const nmod_plain_fft_t F, const ulong k);





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
    
    // inverses of the powers of 2
    mp_ptr powers_inv_2;
    mp_ptr i_powers_inv_2;
    
    mp_ptr * powers_w; // the powers of the successive roots, as in nmod_plain_fft_t
    mp_ptr * i_powers_w; // the same, after precomputation
    
    mp_ptr * powers_inv_w; // the inverse powers of the successive roots 
    mp_ptr * i_powers_inv_w; // the same, after precomputation

    /* mp_ptr * powersR; // the powers of rho such that rho^2 = w */
    /* mp_ptr * i_powersR; // the same, after precomputation */
    
    /* mp_ptr * powers_invR; // the powers of rho such that rho^2 = w */
    /* mp_ptr * i_powers_invR; // the same, after precomputation */
    
} nmod_int128_fft_struct;
typedef nmod_int128_fft_struct nmod_int128_fft_t[1];


/*------------------------------------------------------------*/
/* initializes all entries of F                               */
/* w primitive and w^(2^order))=1                             */
/* DFTs of size up to 2^order are supported                   */ 
/*------------------------------------------------------------*/
void nmod_int128_fft_init_set(nmod_int128_fft_t F, mp_limb_t w, ulong order, nmod_t mod);

/*------------------------------------------------------------*/
/* clears all memory assigned to F                            */
/*------------------------------------------------------------*/
void nmod_int128_fft_clear(nmod_int128_fft_t F);

/*------------------------------------------------------------*/
/* fft evaluation                                             */
/* returns x[i] = poly(w^bitreverse_n(i)), n=2^k              */
/* x must have length >= n                                    */
/*------------------------------------------------------------*/
void nmod_int128_fft_evaluate(mp_ptr x, const nmod_poly_t poly, const nmod_int128_fft_t F, const ulong k);

/*------------------------------------------------------------*/
/* inverse fft                                                */
/* given x[i] = poly(w^bitreverse(i,2^k)), returns poly       */
/*------------------------------------------------------------*/
void nmod_int128_fft_interpolate(nmod_poly_t poly, mp_srcptr x, const nmod_int128_fft_t F, const ulong k);

#endif






#ifdef HAS_AVX2

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* a structure for FFT modulo p                               */
/* assumes that avx2 is supported                             */
/* __int128 used for mod p reductions                         */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
typedef struct
{
    nmod_t mod;
    mp_hlimb_t w, inv_w;       // root of 1 and its inverse
    ulong order;       // its order
    
    // inverses of the powers of 2
    mp_hlimb_t * powers_inv_2;
    mp_hlimb_t * i_powers_inv_2;
    
    mp_hlimb_t ** powers_w; // the powers of the successive roots, as in nmod_plain_fft_t
    mp_hlimb_t ** i_powers_w; // the same, after precomputation
    
    mp_hlimb_t ** powers_inv_w; // the inverse powers of the successive roots 
    mp_hlimb_t ** i_powers_inv_w; // the same, after precomputation

    /* mp_ptr * powersR; // the powers of rho such that rho^2 = w */
    /* mp_ptr * i_powersR; // the same, after precomputation */
    
    /* mp_ptr * powers_invR; // the powers of rho such that rho^2 = w */
    /* mp_ptr * i_powers_invR; // the same, after precomputation */
    
} nmod_avx2_32_fft_struct;
typedef nmod_avx2_32_fft_struct nmod_avx2_32_fft_t[1];

/*------------------------------------------------------------*/
/* initializes all entries of F                               */
/* w primitive and w^(2^order))=1                             */
/* DFTs of size up to 2^order are supported                   */ 
/*------------------------------------------------------------*/
void nmod_avx2_32_fft_init_set(nmod_avx2_32_fft_t F, mp_limb_t w, ulong order, nmod_t mod);

/*------------------------------------------------------------*/
/* clears all memory assigned to F                            */
/*------------------------------------------------------------*/
void nmod_avx2_32_fft_clear(nmod_avx2_32_fft_t F);

/*------------------------------------------------------------*/
/* fft evaluation                                             */
/* returns x[i] = poly(w^bitreverse_n(i)), n=2^k              */
/* x must have length >= n                                    */
/*------------------------------------------------------------*/
void nmod_avx2_32_fft_evaluate(mp_ptr x, const nmod_poly_t poly, const nmod_avx2_32_fft_t F, const ulong k);

/*------------------------------------------------------------*/
/* inverse fft                                                */
/* given x[i] = poly(w^bitreverse(i,2^k)), returns poly       */
/*------------------------------------------------------------*/
void nmod_avx2_32_fft_interpolate(nmod_poly_t poly, mp_srcptr x, const nmod_avx2_32_fft_t F, const ulong k);

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

#endif
