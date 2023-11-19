#include "nmod_extra.h"

/*------------------------------------------------------------*/
/* clears all data in C                                       */
/*------------------------------------------------------------*/
void nmod_multimod_CRT_clear(nmod_multimod_CRT_t C)
{
    if (C->p >= (1L << 50)) // large modulus
        _nmod_vec_clear(C->data);
}
