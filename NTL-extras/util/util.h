#ifndef __UTIL_H
#define __UTIL_H

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
