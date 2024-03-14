#include <assert.h>
#include <flint/nmod.h>
#include <flint/nmod_vec.h>

#include "nmod_vec_extra.h"
 

/*--------------------------------------------------------------*/
/* computes a dot product in size len modulo n                  */
/*--------------------------------------------------------------*/
void check_nmod_vec_dot_small_modulus(ulong len, ulong n, flint_rand_t state)
{
    mp_limb_t res1, res2;
    mp_ptr v1, v2;
    nmod_t mod;

    nmod_init(&mod, n);

    v1 = _nmod_vec_init(len);
    v2 = _nmod_vec_init(len);

    _nmod_vec_rand(v1, state, len, mod);
    _nmod_vec_rand(v2, state, len, mod);
    
    res1 = _nmod_vec_dot(v1, v2, len, mod, _nmod_vec_dot_bound_limbs(len, mod));
    res2 = _nmod_vec_dot_small_modulus(v1, v2, len, (1L << 45) % n, n, 1 / (double) n);
    assert (res1 == res2);
    
    _nmod_vec_clear(v1);
    _nmod_vec_clear(v2);
}

/*--------------------------------------------------------------*/
/* main calls check                                             */
/*--------------------------------------------------------------*/
int main()
{
    flint_rand_t state;
    flint_randinit(state);

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
            check_nmod_vec_dot_small_modulus(len, (UWORD(1) << 2) + 1, state);
            check_nmod_vec_dot_small_modulus(len, (UWORD(1) << 6) + 1, state);
            check_nmod_vec_dot_small_modulus(len, (UWORD(1) << 10) + 1, state);
            check_nmod_vec_dot_small_modulus(len, (UWORD(1) << 14) + 1, state);
            check_nmod_vec_dot_small_modulus(len, (UWORD(1) << 18) + 1, state);
            check_nmod_vec_dot_small_modulus(len, (UWORD(1) << 22) + 1, state);
            check_nmod_vec_dot_small_modulus(len, (UWORD(1) << 26) + 1, state);
            check_nmod_vec_dot_small_modulus(len, (UWORD(1) << 29) + 1, state);
        }
    }
    printf("\n");

    flint_randclear(state);
    return 0;
}
