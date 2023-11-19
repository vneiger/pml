#ifndef __LZZ_PX_MIDDLE_PRODUCT__H
#define __LZZ_PX_MIDDLE_PRODUCT__H

/** Middle product for univariate polynomials over `zz_p`
 *
 * \file lzz_pX_middle_product.h
 * \author Seung Gyu Hyun, Vincent Neiger, Eric Schost
 * \version 0.1
 * \date 2018-12-19
 *
 */

#define KARX (32)

#include <NTL/lzz_pX.h>

NTL_CLIENT

/** Computes the transposed product `xp` of `a` by `b` (naive algorithm) */
void tPlainMul(zz_p *xp, const zz_p *ap, const zz_p *bp, long N);

/*------------------------------------------------------------*/
/* size of the extra storage for Karatsuba                    */
/*------------------------------------------------------------*/
long Kar_stk_size(long N);

/*------------------------------------------------------------*/
/* Karatsuba transposed product, switches to naive            */
/*------------------------------------------------------------*/
void tKarMul_aux(zz_p *b, const long sb, const zz_p *a, const long sa, const zz_p *c, const long sc, zz_p *stk);

/*------------------------------------------------------------*/
/* middle product via FFT                                     */
/* returns trunc( trunc(a, dA+1)*c div x^dA, dB+1 )           */
/* output can alias input                                     */
/*------------------------------------------------------------*/
void middle_FFT(zz_pX& x, const zz_pX& a, const zz_pX& c, long dA, long dB);

/*------------------------------------------------------------*/
/* returns x=(a*b div x^(N-1)) mod x^N                        */
/*------------------------------------------------------------*/
inline void middle_FFT(zz_pX& x, const zz_pX& a, const zz_pX& b, long N)
{ middle_FFT(x, a, b, N-1, N-1); }

/*------------------------------------------------------------*/
/* middle product of (a,c)                                    */
/* returns trunc( trunc(a, dA+1)*c div x^dA, dB+1 )           */
/* output can alias input                                     */
/*------------------------------------------------------------*/
void middle_product(zz_pX& b, const zz_pX& a, const zz_pX& c, long dA, long dB);

inline zz_pX middle_product(const zz_pX& a, const zz_pX& c, long dA, long dB)
{
    zz_pX b;
    middle_product(b, a, c, dA, dB);
    return b;
}

/*------------------------------------------------------------*/
/* middle product of (a,c)                                    */
/*   a has length <= N                                        */
/*   c has length <= 2*N-1                                    */
/* returns b=(a*c div x^(N-1)) mod x^N                        */
/*------------------------------------------------------------*/
inline zz_pX middle_product(const zz_pX& a, const zz_pX& c, long N)
{
    zz_pX b;
    middle_product(b, a, c, N-1, N-1);
    return b;
}

#endif

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
