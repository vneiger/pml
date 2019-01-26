#ifndef __UTIL_H
#define __UTIL_H

/** \brief Some useful functions and macros
 *
 * \file util.h
 * \author Seung Gyu Hyun, Vincent Neiger, Eric Schost
 * \version 0.1
 * \date 2018-12-19
 *
 */

#include <NTL/version.h>
#include <NTL/tools.h>

NTL_CLIENT

#if ((NTL_MAJOR_VERSION == 11) && (NTL_MINOR_VERSION == 1) && (NTL_REVISION == 0))
/** If NTL's version is 11.1.0, define the macro `__NTL_FIX_SIZE_2_FFT` to fix
 * a bug in size-2 FFT in `lzz_pX` */
#define __NTL_FIX_SIZE_2_FFT
#endif

#if (NTL_MAJOR_VERSION >= 11)
/** Macro `get_time` which wraps either GetTime (NTL prior to v11) or
 * GetWallTime (NTL from v11) */
#define get_time GetWallTime
#else
/** Macro `get_time` which wraps either GetTime (NTL prior to v11) or
 * GetWallTime (NTL from v11) */
#define get_time GetTime
#endif

/** Warms the CPU up (currently naive: while loop with empty body, lasting one
 * second) */
void warmup();

#endif

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
