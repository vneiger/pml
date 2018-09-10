#ifndef __THRESHOLDS_MATRIX_MIDDLE_PRODUCT__H
#define __THRESHOLDS_MATRIX_MIDDLE_PRODUCT__H

#include <limits.h>
#include "thresholds_mp_naive_evaluate.h"

/*------------------------------------------------------------*/
/* max degree for evaluate                                    */
/*------------------------------------------------------------*/
inline long max_degree_mp_evaluate(long sz)
{
    long t = type_of_prime();
    if (t == TYPE_FFT_PRIME)
        return LONG_MAX;

    long i;
    for (i = 0; i < MATRIX_MP_THRESHOLDS_LEN - 1; )
        if (sz > MATRIX_MP_THRESHOLDS_SIZES[i])
            i++;
        else 
            break;

    if (t == TYPE_SMALL_PRIME)
        return MATRIX_MP_DEGREE_THRESHOLDS_SMALL[i];
    else
        return MATRIX_MP_DEGREE_THRESHOLDS_LARGE[i];
}

/*------------------------------------------------------------*/
/* max degree for naive                                       */
/*------------------------------------------------------------*/
inline long max_degree_mp_naive(long sz)
{
    long t = type_of_prime();
    long i;
    for (i = 0; i < MATRIX_MP_THRESHOLDS_LEN - 1; )
        if (sz > MATRIX_MP_THRESHOLDS_SIZES[i])
            i++;
        else 
            break;

    if (t == TYPE_FFT_PRIME)
        return MATRIX_MP_NAIVE_THRESHOLDS_FFT[i];
    else 
        if (t == TYPE_SMALL_PRIME)
            return MATRIX_MP_NAIVE_THRESHOLDS_SMALL[i];
        else
            return MATRIX_MP_NAIVE_THRESHOLDS_LARGE[i];
}

#endif

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
