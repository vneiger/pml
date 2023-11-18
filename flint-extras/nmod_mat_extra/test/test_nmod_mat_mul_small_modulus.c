#include <assert.h>
#include <gmp.h>
#include <flint/flint.h>
#include <flint/nmod.h>
#include <flint/nmod_vec.h>

#include "nmod_mat_extra.h"
 

/*--------------------------------------------------------------*/
/* computes a square matrix product in size len modulo n        */
/*--------------------------------------------------------------*/
void check_nmod_mat_mul_small_modulus(ulong len, ulong n)
{
    flint_rand_t state;
    nmod_mat_t a, b, c1, c2;
    nmod_t mod;

    flint_randinit(state);
    nmod_init(&mod, n);

    nmod_mat_init(a, len, len, mod.n);
    nmod_mat_init(b, len, len, mod.n);
    nmod_mat_init(c1, len, len, mod.n);
    nmod_mat_init(c2, len, len, mod.n);

    nmod_mat_rand(a, state);
    nmod_mat_rand(b, state);
    nmod_mat_rand(c1, state);
    nmod_mat_rand(c2, state);
    
    nmod_mat_mul(c1, a, b);
    nmod_mat_mul_small_modulus(c2, a, b);

    assert (nmod_mat_equal(c1, c2));
    
    nmod_mat_clear(a);
    nmod_mat_clear(b);
    nmod_mat_clear(c1);
    nmod_mat_clear(c2);
    flint_randclear(state);
}

/*--------------------------------------------------------------*/
/* main calls check                                             */
/*--------------------------------------------------------------*/
int main(int argc, char **argv)
{
    ulong i;
    for (i = 1; i < 1000; i += 11)
	check_nmod_mat_mul_small_modulus(i, (1L << 29) + 1);

    return 0;
}
