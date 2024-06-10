#include <assert.h>
#include <flint/nmod.h>  // for nmod_init
#include <flint/nmod_vec.h>

#include "nmod_vec_extra.h"


/*--------------------------------------------------------------*/
/* computes a double dot product in size len modulo n           */
/*--------------------------------------------------------------*/
void check_nmod_vec_dot2_small_modulus(ulong len, ulong n, flint_rand_t state)
{
    nmod_t mod;
    nmod_init(&mod, n);

    nn_ptr v11 = aligned_alloc(32, (4 + ((len >> 2) << 2)) * sizeof(ulong));
    nn_ptr v12 = aligned_alloc(32, (4 + ((len >> 2) << 2)) * sizeof(ulong));
    nn_ptr v2 = _nmod_vec_init(len);

    _nmod_vec_rand(v11, state, len, mod);
    _nmod_vec_rand(v12, state, len, mod);
    _nmod_vec_rand(v2, state, len, mod);

    vec2d p2, pinv2;
    p2[0] = n;
    p2[1] = n;
    pinv2[0] = 1 / p2[0];
    pinv2[1] = 1 / p2[1];

    ulong res1[2];
    res1[0] = _nmod_vec_dot(v11, v2, len, mod, _nmod_vec_dot_bound_limbs(len, mod));
    res1[1] = _nmod_vec_dot(v12, v2, len, mod, _nmod_vec_dot_bound_limbs(len, mod));

    ulong res2[2];
    _nmod_vec_dot2_small_modulus(res2, v11, v12, v2, len, (UWORD(1) << 45) % n, p2, pinv2);

    assert (res1[0] == res2[0]);
    assert (res1[1] == res2[1]);

    _nmod_vec_clear(v11);
    _nmod_vec_clear(v12);
    _nmod_vec_clear(v2);
}

/*--------------------------------------------------------------*/
/* main calls check                                             */
/*--------------------------------------------------------------*/
int main()
{
    flint_rand_t state;
    flint_rand_init(state);

    printf("test running, various bitlengths, len from 1 to 999 (no error message means success)...\n");
    for (slong len = 1; len < 1000; len += 1)
    {
        if (len % 20 == 0)
        {
            printf("%ld..", len);
            fflush(stdout);
        }
        for (slong repeat = 0; repeat < 30; repeat++)
        {
            check_nmod_vec_dot2_small_modulus(len, (UWORD(1) << 2) + 1, state);
            check_nmod_vec_dot2_small_modulus(len, (UWORD(1) << 6) + 1, state);
            check_nmod_vec_dot2_small_modulus(len, (UWORD(1) << 10) + 1, state);
            check_nmod_vec_dot2_small_modulus(len, (UWORD(1) << 14) + 1, state);
            check_nmod_vec_dot2_small_modulus(len, (UWORD(1) << 18) + 1, state);
            check_nmod_vec_dot2_small_modulus(len, (UWORD(1) << 22) + 1, state);
            check_nmod_vec_dot2_small_modulus(len, (UWORD(1) << 26) + 1, state);
            check_nmod_vec_dot2_small_modulus(len, (UWORD(1) << 29) + 1, state);
        }
    }
    printf("\n");

    flint_rand_clear(state);

    return 0;
}
