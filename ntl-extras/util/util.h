#ifndef __UTIL_H
#define __UTIL_H

#include <NTL/tools.h>

/** \brief Some useful functions and macros
 *
 * \file util.h
 * \author Seung Gyu Hyun, Vincent Neiger, Eric Schost
 * \date 2025-05-26
 *
 */

//#include <NTL/version.h>
// #if ((NTL_MAJOR_VERSION == 11) && (NTL_MINOR_VERSION == 1) && (NTL_REVISION == 0))

#define PML_OPEN_NNS namespace PML {
#define PML_CLOSE_NNS  }
#define PML_USE_PNS using namespace PML;
#define PML_START_IMPL PML_OPEN_NNS NTL_USE_NNS NTL_IMPORT_FROM_STD
#define PML_END_IMPL PML_CLOSE_NNS
#define PML_START_IMPL_IO \
            PML_OPEN_NNS \
            NTL_USE_NNS \
            NTL_IMPORT_FROM_STD \
            using std::string; \
            using std::cout; \
            using std::endl;
#define PML_END_IMPL_IO PML_CLOSE_NNS
#define PML_CLIENT NTL_USE_SNS NTL_USE_NNS PML_USE_PNS

PML_OPEN_NNS

/** Warms the CPU up (currently naive: while loop with empty body, lasting one
 * second) */
void warmup();

PML_CLOSE_NNS

#endif

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
