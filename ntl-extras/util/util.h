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

#include <NTL/tools.h>

//#include <NTL/version.h>
// #if ((NTL_MAJOR_VERSION == 11) && (NTL_MINOR_VERSION == 1) && (NTL_REVISION == 0))

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
