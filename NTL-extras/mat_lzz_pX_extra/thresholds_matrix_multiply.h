#ifndef __THRESHOLDS_MATRIX_MULTIPLY__H
#define __THRESHOLDS_MATRIX_MULTIPLY__H

#include <limits.h>
#include "thresholds_waksman_evaluate.h"

#define MAX_DEGREE_TRANSFORM_FFT 3
#define MAX_DEGREE_TRANSFORM_SMALL 4 
#define MAX_DEGREE_TRANSFORM_LARGE 7

/*------------------------------------------------------------*/
/* max degree for using transforms                            */
/*------------------------------------------------------------*/
inline long max_degree_transform()
{
    long t = type_of_prime();

    if (t == TYPE_FFT_PRIME)
        return MAX_DEGREE_TRANSFORM_FFT;
    else 
        if (t == TYPE_SMALL_PRIME)
            return MAX_DEGREE_TRANSFORM_SMALL;
        else
            return MAX_DEGREE_TRANSFORM_LARGE;
}

/*------------------------------------------------------------*/
/* max degree for evaluate                                    */
/*------------------------------------------------------------*/
inline long max_degree_evaluate(long sz)
{
    long t = type_of_prime();
    if (t == TYPE_FFT_PRIME)
        return LONG_MAX;

    long i;
    for (i = 0; i < MATRIX_THRESHOLDS_LEN - 1; )
        if (sz > MATRIX_THRESHOLDS_SIZES[i])
            i++;
        else 
            break;

    if (t == TYPE_SMALL_PRIME)
        return MATRIX_DEGREE_THRESHOLDS_SMALL[i];
    else
        return MATRIX_DEGREE_THRESHOLDS_LARGE[i];
}

/*------------------------------------------------------------*/
/* max degree for waskman                                     */
/*------------------------------------------------------------*/
inline long max_degree_waksman(long sz)
{
    long t = type_of_prime();
    long i;
    for (i = 0; i < MATRIX_THRESHOLDS_LEN - 1; )
        if (sz > MATRIX_THRESHOLDS_SIZES[i])
            i++;
        else 
            break;

    if (t == TYPE_FFT_PRIME)
        return MATRIX_WAKSMAN_THRESHOLDS_FFT[i];
    else 
        if (t == TYPE_SMALL_PRIME)
            return MATRIX_WAKSMAN_THRESHOLDS_SMALL[i];
        else
            return MATRIX_WAKSMAN_THRESHOLDS_LARGE[i];
}

#endif

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
