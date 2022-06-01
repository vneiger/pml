#include <flint/nmod_poly.h>
#include "nmod_vec_extra.h"

/*------------------------------------------------------------*/
/* multiplies a vector by a scalar                            */
/*------------------------------------------------------------*/
void get_time()
{
#ifdef HAS_AVX2
    slong i;
    flint_rand_t state;
    mp_limb_t nn;
    nmod_t mod;
    
    flint_randinit(state);
    nn = 65537;
    nmod_init(&mod, nn);

    for (i = 100; i < 10000; i+=1000)
    {
        slong ia, j;
        mp_hlimb_t r, i_r;
        mp_hlimb_t * val;
        double t;
        clock_t tt;
        long nb_iter;

        ia = 8 * ((i+7)/8); // next multiple of 8
        val = (mp_hlimb_t *) aligned_alloc(32, 4 * ia);
        for (j = 0; j < i; j++)
        {
            val[j] = n_randtest(state) % nn;

        }

        r = nn/2;
        i_r = prep_mul_mod_precon_32(r, mod.n);
        
        t = 0.0;
        nb_iter = 0;
        while (t < 0.5)
        {
            tt = clock();
            _nmod_vec_avx2_32_scalar_mul(val, val, i, r, i_r, mod);
            t += (double)(clock()-tt) / CLOCKS_PER_SEC;
            ++nb_iter;
        }
        t = 1000 * t;
        t /= nb_iter;
        printf("%ld %f\n", i, t);

        free(val);
    }
    flint_randclear(state);
#endif
}

/*------------------------------------------------------------*/
/* main just calls get_time()                                 */
/*------------------------------------------------------------*/
int main(int argc, char **argv){
    get_time();
    return 0;
}
