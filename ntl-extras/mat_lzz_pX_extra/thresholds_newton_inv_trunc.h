#ifndef __THRESHOLDS_NEWTON_INV_TRUNC__H
#define __THRESHOLDS_NEWTON_INV_TRUNC__H

#include <limits.h>
#include "thresholds_plain_newton_geometric_mp_inv_trunc.h"

// The behavior middle product vs geometric is too irregular.
// The geometric option is disabled for the moment.
inline long max_degree_middle_product_inv_trunc()
{
    return LONG_MAX;
}

#endif

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
