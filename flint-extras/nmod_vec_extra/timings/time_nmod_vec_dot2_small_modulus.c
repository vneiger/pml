#include <assert.h>
#include <gmp.h>
#include <flint/nmod_vec.h>
#include <flint/fmpz.h>
#include <flint/machine_vectors.h>

#include "nmod_vec_extra.h"
 

/*--------------------------------------------------------------*/
/* computes a dot modulus in size len modulo n                  */
/*--------------------------------------------------------------*/
void time_nmod_vec_dot2_small_modulus(ulong len, ulong n)
{
    flint_rand_t state;
    mp_limb_t res2, res[2], two_32;
    mp_ptr v11, v12, v2;
    nmod_t mod;
    double t;
    clock_t tt;
    long nb_iter;
    vec2d p2, pinv2;
    
    flint_randinit(state);
    nmod_init(&mod, n);

    v11 = aligned_alloc(32, (4 + ((len >> 2) << 2)) * sizeof(mp_limb_t));
    v12 = aligned_alloc(32, (4 + ((len >> 2) << 2)) * sizeof(mp_limb_t));
    v2 = _nmod_vec_init(len);
    p2[0] = n;
    p2[1] = n;
    pinv2[0] = 1 / (double) n;
    pinv2[1] = 1 / (double) n;
    
    _nmod_vec_rand(v11, state, len, mod);
    _nmod_vec_rand(v12, state, len, mod);
    _nmod_vec_rand(v2, state, len, mod);
    two_32 = (1L << 32) % n;
    
    printf("%lu ", len);
    
    t = 0.0;
    nb_iter = 0;
    while (t < 0.5)
    {
        tt = clock();
        _nmod_vec_dot2_small_modulus(res, v11, v12, v2, len, two_32, p2, pinv2, mod);
        _nmod_vec_dot2_small_modulus(res, v11, v12, v2, len, two_32, p2, pinv2, mod);
        _nmod_vec_dot2_small_modulus(res, v11, v12, v2, len, two_32, p2, pinv2, mod);
        _nmod_vec_dot2_small_modulus(res, v11, v12, v2, len, two_32, p2, pinv2, mod);
        _nmod_vec_dot2_small_modulus(res, v11, v12, v2, len, two_32, p2, pinv2, mod);
        t += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += 5;
    }
    t /= nb_iter;
    printf("%4g", t);

    t = 0.0;
    nb_iter = 0;
    while (t < 0.5)
    {
        tt = clock();
        res2 = _nmod_vec_dot_small_modulus(v11, v2, len, mod);
        res2 = _nmod_vec_dot_small_modulus(v12, v2, len, mod);
        res2 = _nmod_vec_dot_small_modulus(v11, v2, len, mod);
        res2 = _nmod_vec_dot_small_modulus(v12, v2, len, mod);
        res2 = _nmod_vec_dot_small_modulus(v11, v2, len, mod);
        res2 = _nmod_vec_dot_small_modulus(v12, v2, len, mod);
        res2 = _nmod_vec_dot_small_modulus(v11, v2, len, mod);
        res2 = _nmod_vec_dot_small_modulus(v12, v2, len, mod);
        res2 = _nmod_vec_dot_small_modulus(v11, v2, len, mod);
        res2 = _nmod_vec_dot_small_modulus(v12, v2, len, mod);
        t += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += 5;
    }
    t /= nb_iter;
    printf(" %4g", t);

    printf("\n");

    
    _nmod_vec_clear(v11);
    _nmod_vec_clear(v12);
    _nmod_vec_clear(v2);
    flint_randclear(state);
}

/*--------------------------------------------------------------*/
/* main calls time                                              */
/*--------------------------------------------------------------*/
int main(int argc, char **argv)
{
    slong i;
    for (i = 1; i < 1000; i += 5)
	time_nmod_vec_dot2_small_modulus(i, (1L << 29) + 1);

    return 0;
}
