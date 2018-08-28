#include <NTL/lzz_p.h>

#include "util.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* warm-up the CPU                                            */
/*------------------------------------------------------------*/
void warmup()
{
    double t = get_time();
    do
    {
        ;
    } 
    while (get_time() - t < 1);
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
