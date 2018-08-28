#ifndef __UTIL_H
#define __UTIL_H

#include <NTL/version.h>

/*------------------------------------------------------------*/
/* bug in size-2 lzz_pX FFT in version 11.1.0                 */
/*------------------------------------------------------------*/
#if ((NTL_MAJOR_VERSION == 11) && (NTL_MINOR_VERSION == 1) && (NTL_REVISION == 0))
#define __NTL_FIX_SIZE_2_FFT
#endif

/*------------------------------------------------------------*/
/* wraps either GetTime (old NTL's) or GetWallTime (>= 11)    */
/*------------------------------------------------------------*/
#if (NTL_MAJOR_VERSION >= 11)
#define get_time GetWallTime
#else
#define get_time GetTime
#endif

/*------------------------------------------------------------*/
/* warm-up the CPU                                            */
/*------------------------------------------------------------*/
void warmup();

#endif

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
