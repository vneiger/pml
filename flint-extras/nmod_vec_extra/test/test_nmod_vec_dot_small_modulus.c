#include <assert.h>
#include <flint/nmod.h>
#include <flint/nmod_vec.h>

#include "nmod_vec_extra.h"


/*--------------------------------------------------------------*/
/* computes a dot product in size len modulo n                  */
/*--------------------------------------------------------------*/
void check_nmod_vec_dot_small_modulus(ulong len, ulong n, flint_rand_t state)
{
    ulong res1, res2, res3;
    nn_ptr v1, v2;
    nmod_t mod;

    nmod_init(&mod, n);

    v1 = _nmod_vec_init(len);
    v2 = _nmod_vec_init(len);

    _nmod_vec_rand(v1, state, len, mod);
    _nmod_vec_rand(v2, state, len, mod);

    res1 = _nmod_vec_dot(v1, v2, len, mod, _nmod_vec_dot_bound_limbs(len, mod));
    res2 = _nmod_vec_dot_mod32(v1, v2, len, mod, (uint) ((1L << DOT_SP_NB) % n));
    res3 = _nmod_vec_dot_mod32_avx2(v1, v2, len, mod, (uint) ((1L << DOT_SP_NB) % n));
    assert (res1 == res2);
    assert (res1 == res3);

    _nmod_vec_clear(v1);
    _nmod_vec_clear(v2);
}

/*--------------------------------------------------------------*/
/* main calls check                                             */
/*--------------------------------------------------------------*/
int main()
{
    flint_rand_t state;
    flint_rand_init(state);

    printf("test running, various bitlengths, various lengths (no error message means success)...\n");
    printf("lens: ");
    for (slong len = 1; len < (slong)100000000; len *= 3)
    {
        printf("%ld..", len);
        fflush(stdout);
        for (slong repeat = 0; repeat < 10; repeat++)
        {
            check_nmod_vec_dot_small_modulus(len, (UWORD(1) << 2) + 1, state);
            check_nmod_vec_dot_small_modulus(len, (UWORD(1) << 6) + 1, state);
            check_nmod_vec_dot_small_modulus(len, (UWORD(1) << 10) + 1, state);
            check_nmod_vec_dot_small_modulus(len, (UWORD(1) << 14) + 1, state);
            check_nmod_vec_dot_small_modulus(len, (UWORD(1) << 18) + 1, state);
            check_nmod_vec_dot_small_modulus(len, (UWORD(1) << 22) + 1, state);
            check_nmod_vec_dot_small_modulus(len, (UWORD(1) << 26) + 1, state);
            check_nmod_vec_dot_small_modulus(len, (UWORD(1) << 30) + 1, state);
        }
    }
    // test max claimed modulus+len
    printf("%ld..", 134744072L);
    fflush(stdout);
    for (slong repeat = 0; repeat < 10; repeat++)
        check_nmod_vec_dot_small_modulus(134744072L, 1515531528L, state);
    printf("\n");

    flint_rand_clear(state);
    return 0;
}
