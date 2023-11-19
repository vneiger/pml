#include <assert.h>
#include <gmp.h>
#include <flint/nmod_vec.h>

#include "nmod_vec_extra.h"
 

/*--------------------------------------------------------------*/
/* computes a dot modulus in size len modulo n                  */
/*--------------------------------------------------------------*/
void check_nmod_vec_dot_small_modulus(ulong len, ulong n)
{
    flint_rand_t state;
    mp_limb_t res1, res2;
    mp_ptr v1, v2;
    nmod_t mod;

    flint_randinit(state);
    nmod_init(&mod, n);

    v1 = _nmod_vec_init(len);
    v2 = _nmod_vec_init(len);

    _nmod_vec_rand(v1, state, len, mod);
    _nmod_vec_rand(v2, state, len, mod);
    
    res1 = _nmod_vec_dot(v1, v2, len, mod, _nmod_vec_dot_bound_limbs(len, mod));
    res2 = _nmod_vec_dot_small_modulus(v1, v2, len, mod);
    assert (res1 == res2);
    
    _nmod_vec_clear(v1);
    _nmod_vec_clear(v2);
    flint_randclear(state);
}

/*--------------------------------------------------------------*/
/* main calls check                                             */
/*--------------------------------------------------------------*/
int main(int argc, char **argv)
{
    slong i;
    for (i = 1; i < 1000; i += 21)
	check_nmod_vec_dot_small_modulus(i, (1L << 29) + 1);

    return 0;
}
