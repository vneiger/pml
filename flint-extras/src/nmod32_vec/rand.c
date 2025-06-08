#include "nmod32_vec.h"

void _nmod32_vec_rand(n32_ptr vec, flint_rand_t state, slong len, nmod_t mod)
{
    for (slong i = 0; i < len; i++)
        vec[i] = (uint)n_randint(state, mod.n);
}
