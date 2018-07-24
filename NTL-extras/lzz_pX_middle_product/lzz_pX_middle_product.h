#ifndef __LZZ_PX_MIDDLE_PRODUCT__H
#define __LZZ_PX_MIDDLE_PRODUCT__H

#define KARX (32)

#include <NTL/lzz_pX.h>

NTL_CLIENT

/*------------------------------------------------------------*/
/* naive transposed product of (a,b)                          */
/*------------------------------------------------------------*/
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
/* returns x=(a*b div x^(N-1)) mod x^N                        */
/*------------------------------------------------------------*/
void middle_FFT(zz_pX& x, const zz_pX& a, const zz_pX& b, long N);

/*------------------------------------------------------------*/
/* middle product of (a,b)                                    */
/*   c has length <= 2*N-1                                    */
/*   a has length <= N                                        */
/*   x has length <= N                                        */
/* returns x=(a*b div x^(N-1)) mod x^N                        */
/*------------------------------------------------------------*/
zz_pX middle_product(const zz_pX& c, const zz_pX& a, long N);


#endif
