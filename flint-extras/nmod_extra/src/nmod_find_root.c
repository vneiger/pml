#include <flint/flint.h>

#include "nmod_extra.h"

/*------------------------------------------------------------*/
/* finds an element of order at least n                       */
/* returns 0 if not found                                     */
/*------------------------------------------------------------*/
mp_limb_t nmod_find_root(long n, nmod_t mod)
{
    long q;
    for(q = 2; q < mod.n; q++)
    {
        slong k = 1;
        slong qk = q;
        while (qk != 1 && k < n)
        {
            qk = nmod_mul(qk, q, mod);
            k++;
        }
        if (qk != 1)
        {
            return q;
        }
    }
    return 0; 
}

