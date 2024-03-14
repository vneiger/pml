#include <time.h>
#include <gmp.h>
#include <flint/nmod.h>
#include <flint/nmod_vec.h>

#include "nmod_vec_extra.h"


/*--------------------------------------------------------------*/
/* computes a dot modulus in size len modulo n                  */
/*--------------------------------------------------------------*/
void time_nmod_vec_dot2_small_modulus(ulong len, ulong n, flint_rand_t state)
{
    nmod_t mod;
    nmod_init(&mod, n);

    mp_ptr v11 = aligned_alloc(32, (4 + ((len >> 2) << 2)) * sizeof(mp_limb_t));
    mp_ptr v12 = aligned_alloc(32, (4 + ((len >> 2) << 2)) * sizeof(mp_limb_t));
    mp_ptr v2 = _nmod_vec_init(len);

    const vec1d p = n;
    const vec1d pinv = 1 / (double) n;

    vec2d p2, pinv2;
    p2[0] = n;
    p2[1] = n;
    pinv2[0] = 1 / p2[0];
    pinv2[1] = 1 / p2[1];

    const mp_limb_t two_32 = (UWORD(1) << 45) % n;

    double t;
    clock_t tt;
    long nb_iter;

    t = 0.0;
    nb_iter = 0;
    while (t < 0.5)
    {
        mp_limb_t res[2];
    _nmod_vec_rand(v11, state, len, mod);
    _nmod_vec_rand(v12, state, len, mod);

        tt = clock();
        _nmod_vec_dot2_small_modulus(res, v11, v12, v2, len, two_32, p2, pinv2);
        _nmod_vec_dot2_small_modulus(res, v11, v12, v2, len, two_32, p2, pinv2);
        _nmod_vec_dot2_small_modulus(res, v11, v12, v2, len, two_32, p2, pinv2);
        _nmod_vec_dot2_small_modulus(res, v11, v12, v2, len, two_32, p2, pinv2);
        _nmod_vec_dot2_small_modulus(res, v11, v12, v2, len, two_32, p2, pinv2);
        t += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += 5;
    }
    t = 1000*t;
    t /= nb_iter;
    printf("%.1e\t", t);

    t = 0.0;
    nb_iter = 0;
    while (t < 0.5)
    {
    _nmod_vec_rand(v2, state, len, mod);
        tt = clock();
        _nmod_vec_dot_small_modulus(v11, v2, len, two_32, p, pinv);
        _nmod_vec_dot_small_modulus(v12, v2, len, two_32, p, pinv);
        _nmod_vec_dot_small_modulus(v11, v2, len, two_32, p, pinv);
        _nmod_vec_dot_small_modulus(v12, v2, len, two_32, p, pinv);
        _nmod_vec_dot_small_modulus(v11, v2, len, two_32, p, pinv);
        _nmod_vec_dot_small_modulus(v12, v2, len, two_32, p, pinv);
        _nmod_vec_dot_small_modulus(v11, v2, len, two_32, p, pinv);
        _nmod_vec_dot_small_modulus(v12, v2, len, two_32, p, pinv);
        _nmod_vec_dot_small_modulus(v11, v2, len, two_32, p, pinv);
        _nmod_vec_dot_small_modulus(v12, v2, len, two_32, p, pinv);
        t += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += 5;
    }
    t = 1000*t;
    t /= nb_iter;
    printf("%.1e\t", t);

    _nmod_vec_clear(v11);
    _nmod_vec_clear(v12);
    _nmod_vec_clear(v2);
}

/*--------------------------------------------------------------*/
/* main calls time                                              */
/*--------------------------------------------------------------*/
int main()
{
    flint_rand_t state;
    flint_randinit(state);

    printf("len\t4\t\t10\t\t20\t\t25\t\t29\t\t30\n");
    for (slong len = 1; len < 1000; len += 21)
    {
        printf("%ld\t", len);
        time_nmod_vec_dot2_small_modulus(len, (UWORD(1) << 3) + 1, state);
        time_nmod_vec_dot2_small_modulus(len, (UWORD(1) << 9) + 1, state);
        time_nmod_vec_dot2_small_modulus(len, (UWORD(1) << 19) + 1, state);
        time_nmod_vec_dot2_small_modulus(len, (UWORD(1) << 24) + 1, state);
        time_nmod_vec_dot2_small_modulus(len, (UWORD(1) << 28) + 1, state);
        time_nmod_vec_dot2_small_modulus(len, (UWORD(1) << 29) + 1, state);
        printf("\n");
    }

    flint_randclear(state);

    return 0;
}
