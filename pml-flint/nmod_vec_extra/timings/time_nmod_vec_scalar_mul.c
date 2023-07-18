#include <flint/nmod_poly.h>
#include "nmod_vec_extra.h"

/*------------------------------------------------------------*/
/* multiplies a vector by a scalar                            */
/*------------------------------------------------------------*/
void get_time()
{
    slong i;
    flint_rand_t state;
    mp_limb_t nn;
    nmod_t mod;
    
    flint_randinit(state);
    nn = 65537;
    nmod_init(&mod, nn);

    for (i = 100; i < 10000; i+=1000)
    {
        mp_limb_t r;
        mp_ptr val;
        double t;
        clock_t tt;
        long nb_iter;
        
        val = _nmod_vec_init(i);
        _nmod_vec_randtest(val, state, i, mod);

        r = nn/2;
        t = 0.0;
        nb_iter = 0;
        while (t < 0.5)
        {
            tt = clock();
            _nmod_vec_scalar_mul(val, val, i, r, mod);
            t += (double)(clock()-tt) / CLOCKS_PER_SEC;
            ++nb_iter;
        }
        t = 1000 * t;
        t /= nb_iter;
        printf("%ld %f\n", i, t);

        _nmod_vec_clear(val);
    }
    flint_randclear(state);
}

/*------------------------------------------------------------*/
/* main just calls get_time()                                 */
/*------------------------------------------------------------*/
int main(int argc, char **argv){
    get_time();
    return 0;
}
