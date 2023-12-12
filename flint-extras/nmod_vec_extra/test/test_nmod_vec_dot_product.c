#include <assert.h>
#include <gmp.h>
#include <flint/nmod_vec.h>

#include "nmod_vec_extra.h"
 

/*--------------------------------------------------------------*/
/* computes a dot product in size len modulo n                  */
/*--------------------------------------------------------------*/
void check_nmod_vec_dot_product(ulong len, ulong bits1, ulong bits2, ulong n)
{
    flint_rand_t state;
    mp_limb_t res1, res2;
    mp_ptr v1, v2, v1r, v2r;
    nmod_t mod, mod1, mod2;
    ulong i;
    
    flint_randinit(state);
    nmod_init(&mod, n);
    nmod_init(&mod1, 1L << bits1);
    nmod_init(&mod2, 1L << bits2);

    v1 = _nmod_vec_init(len);
    v2 = _nmod_vec_init(len);
    v1r = _nmod_vec_init(len);
    v2r = _nmod_vec_init(len);

    _nmod_vec_rand(v1, state, len, mod1);
    _nmod_vec_rand(v2, state, len, mod2);

    for (i = 0; i < len; i++)
    {
        v1r[i] = v1[i] % n;
        v2r[i] = v2[i] % n;
    }
    
    res1 = _nmod_vec_dot(v1r, v2r, len, mod, _nmod_vec_dot_bound_limbs(len, mod));
    res2 = nmod_vec_dot_product(v1, v2, len, bits1, bits2, mod);
    assert (res1 == res2);

    _nmod_vec_clear(v1r);
    _nmod_vec_clear(v2r);
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

    for (i = 1; i < 1000; i += 1)
	check_nmod_vec_dot_product(i, 29, 63, (1L << 29) + 1);

    for (i = 1; i < 1000; i += 1)
	check_nmod_vec_dot_product(i, 63, 63, (1L << 63) + 1);
    

    return 0;
}
