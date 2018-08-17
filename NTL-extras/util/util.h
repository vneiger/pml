#ifndef __UTIL_H
#define __UTIL_H

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
