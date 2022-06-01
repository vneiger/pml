#include <assert.h>
#include <flint/flint.h>

#include "nmod_vec_extra.h"

/*------------------------------------------------------------*/
/* multiplies a vector by a scalar                            */
/*------------------------------------------------------------*/
void check()
{
#ifdef HAS_AVX2
    slong i;
    flint_rand_t state;
    mp_limb_t nn;
    nmod_t mod;
    
    flint_randinit(state);
    nn = 65537;
    nmod_init(&mod, nn);

    for (i = 1; i < 100; i+=1)
    {
        slong ia, j;
        mp_hlimb_t r, i_r;
        mp_hlimb_t * val, * val2;

        ia = 8 * ((i+7)/8); // next multiple of 8
        val = (mp_hlimb_t *) aligned_alloc(32, 4 * ia);
        val2 = (mp_hlimb_t *) aligned_alloc(32, 4 * ia);

        for (j = 0; j < i; j++)
        {
            val[j] = n_randtest(state) % nn;
            val2[j] = val[j];
        }

        r = nn/2;
        i_r = prep_mul_mod_precon_32(r, mod.n);
        _nmod_vec_avx2_32_scalar_mul(val, val, i, r, i_r, mod);

        for (j = 0; j < i; j++)
        {
            assert (val[j] == nmod_mul(val2[j], r, mod));
        }
        
        free(val);
        free(val2);
    }
    flint_randclear(state);
#endif
}

/*------------------------------------------------------------*/
/* main just calls check()                                    */
/*------------------------------------------------------------*/
int main(int argc, char **argv){
    check();
    return 0;
}
