#ifndef __LZZ_P_EXTRA__H
#define __LZZ_P_EXTRA__H

NTL_CLIENT

/*------------------------------------------------------------*/
/* multiplicative order of a                                  */
/* -1 if a = 0                                                */
/*------------------------------------------------------------*/
long order(const zz_p& a);

/*------------------------------------------------------------*/
/* finds an element of order at least ord                     */
/* assumes it exists, does not verify                         */
/*------------------------------------------------------------*/
void element_of_order(zz_p& a, long ord);

/*------------------------------------------------------------*/
/* 1 if the current prime can be used as an FFT prime         */
/*------------------------------------------------------------*/
long is_FFT_ready();

#endif
