#include <assert.h>
#include <gmp.h>
#include <flint/nmod_vec.h>

#include "nmod_vec_extra.h"
 

/*--------------------------------------------------------------*/
/* computes an integer dot-product, length len, bitsize bit_len */
/* checks modulo a prime                                        */
/*--------------------------------------------------------------*/
void check_nmod_vec_integer_dot_product(slong len, mp_bitcnt_t bit_len)
{
    flint_rand_t state;
    mp_limb_t p, res1, res2;
    mp_ptr v1, v2, res;
    nmod_t mod;

    flint_randinit(state);
    p = n_randprime(state, bit_len, 0);
    nmod_init(&mod, p);

    v1 = _nmod_vec_init(len);
    v2 = _nmod_vec_init(len);

    _nmod_vec_randtest(v1, state, len, mod);
    _nmod_vec_randtest(v2, state, len, mod);

    res = _nmod_vec_init(3);
    nmod_vec_integer_dot_product(res, v1, v2, len, bit_len, bit_len);
    NMOD_RED3(res1, res[2], res[1], res[0], mod);
    res2 = _nmod_vec_dot(v1, v2, len, mod, 3);

    assert (res1 == res2);
    
    _nmod_vec_clear(res);
    _nmod_vec_clear(v1);
    _nmod_vec_clear(v2);
    flint_randclear(state);
}

/*--------------------------------------------------------------*/
/* main calls check                                             */
/*--------------------------------------------------------------*/
int main(int argc, char**argv)
{
    slong i;
    for (i = 1; i < 1000; i += 20)
    {
	check_nmod_vec_integer_dot_product(i, 59);
    }
    return 0;
}
