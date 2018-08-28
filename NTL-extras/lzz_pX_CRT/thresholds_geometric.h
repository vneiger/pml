#ifndef __THRESHOLDS_GEOMETRIC_FFT_H
#define __THRESHOLDS_GEOMETRIC_FFT_H

#include <NTL/ZZ.h>

#include "thresholds_geometric_evaluate_FFT.h"
#include "thresholds_geometric_interpolate_FFT.h"
#include "lzz_p_extra.h"

NTL_CLIENT

#define MIN_GEOMETRIC_FFT_SMALL ((MIN_GEOMETRIC_FFT_EVALUATE_SMALL+MIN_GEOMETRIC_FFT_INTERPOLATE_SMALL)/2)
#define MIN_GEOMETRIC_FFT_LARGE ((MIN_GEOMETRIC_FFT_EVALUATE_LARGE+MIN_GEOMETRIC_FFT_INTERPOLATE_LARGE)/2)
#define MIN_GEOMETRIC_FFT_FFT ((MIN_GEOMETRIC_FFT_EVALUATE_FFT+MIN_GEOMETRIC_FFT_INTERPOLATE_FFT)/2)

#define MAX_GEOMETRIC_FFT_SMALL ((MAX_GEOMETRIC_FFT_EVALUATE_SMALL+MAX_GEOMETRIC_FFT_INTERPOLATE_SMALL)/2)
#define MAX_GEOMETRIC_FFT_LARGE ((MAX_GEOMETRIC_FFT_EVALUATE_LARGE+MAX_GEOMETRIC_FFT_INTERPOLATE_LARGE)/2)
#define MAX_GEOMETRIC_FFT_FFT ((MAX_GEOMETRIC_FFT_EVALUATE_FFT+MAX_GEOMETRIC_FFT_INTERPOLATE_FFT)/2)

/*------------------------------------------------------------*/
/* min degree for using no-FFT operations                     */
/*------------------------------------------------------------*/
inline long min_geometric_FFT()
{
    long t = type_of_prime();
    switch (t){
    case TYPE_FFT_PRIME:
        return MIN_GEOMETRIC_FFT_FFT;
    case TYPE_SMALL_PRIME:
        return MIN_GEOMETRIC_FFT_SMALL;
    case TYPE_LARGE_PRIME:
        return MIN_GEOMETRIC_FFT_LARGE;
    }
    LogicError("Type of prime unknown");
    return -1;
}

/*------------------------------------------------------------*/
/* max degree for using no-FFT operations                     */
/*------------------------------------------------------------*/
inline long max_geometric_FFT()
{
    long t = type_of_prime();
    switch (t){
    case TYPE_FFT_PRIME:
        return MAX_GEOMETRIC_FFT_FFT;
    case TYPE_SMALL_PRIME:
        return MAX_GEOMETRIC_FFT_SMALL;
    case TYPE_LARGE_PRIME:
        return MAX_GEOMETRIC_FFT_LARGE;
    }
    LogicError("Type of prime unknown");
    return -1;
}


#endif

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
