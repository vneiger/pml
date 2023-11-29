#include <assert.h>
#include <gmp.h>
#include <flint/nmod_vec.h>
#include <flint/machine_vectors.h>

#include "nmod_vec_extra.h"


/*--------------------------------------------------------------*/
/* computes a double dot product in size len modulo n           */
/*--------------------------------------------------------------*/
void check_nmod_vec_dot2_small_modulus(ulong len, ulong n)
{
    flint_rand_t state;
    mp_limb_t res1[2], res2[2];
    mp_ptr v11, v12, v2;
    nmod_t mod;
    vec2d p2, pinv2;

    flint_randinit(state);
    nmod_init(&mod, n);

    v11 = aligned_alloc(32, (4 + ((len >> 2) << 2)) * sizeof(mp_limb_t));
    v12 = aligned_alloc(32, (4 + ((len >> 2) << 2)) * sizeof(mp_limb_t));
    v2 = _nmod_vec_init(len);

    _nmod_vec_rand(v11, state, len, mod);
    _nmod_vec_rand(v12, state, len, mod);
    _nmod_vec_rand(v2, state, len, mod);
    p2[0] = n;
    p2[1] = n;
    pinv2[0] = 1 / p2[0];
    pinv2[1] = 1 / p2[1];

    res1[0] = _nmod_vec_dot(v11, v2, len, mod, _nmod_vec_dot_bound_limbs(len, mod));
    res1[1] = _nmod_vec_dot(v12, v2, len, mod, _nmod_vec_dot_bound_limbs(len, mod));
    _nmod_vec_dot2_small_modulus(res2, v11, v12, v2, len, (1L << 45) % n, p2, pinv2);
    assert (res1[0] == res2[0]);
    assert (res1[1] == res2[1]);

    _nmod_vec_clear(v11);
    _nmod_vec_clear(v12);
    _nmod_vec_clear(v2);
    flint_randclear(state);
}

/*--------------------------------------------------------------*/
/* main calls check                                             */
/*--------------------------------------------------------------*/
int main(int argc, char **argv)
{
    slong i;
    for (i = 1; i < 320; i += 1)
	check_nmod_vec_dot2_small_modulus(i, (1L << 29) + 1);

    return 0;
}
