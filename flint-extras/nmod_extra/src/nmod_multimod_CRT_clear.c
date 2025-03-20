#include "nmod_extra.h"

/*------------------------------------------------------------*/
/* clears all data in C                                       */
/*------------------------------------------------------------*/
void nmod_multimod_CRT_clear(nmod_multimod_CRT_t C)
{
#if FLINT_HAVE_FFT_SMALL
    if (C->p >= (1L << 50)) // large modulus
#endif  // FLINT_HAVE_FFT_SMALL
        if (C->num_primes > 1)
            _nmod_vec_clear(C->data);
}
