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
