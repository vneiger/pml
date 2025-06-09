#include <assert.h>
#include <flint/flint.h>

#include "nmod_extra.h"
#include "nmod_poly_extra.h"

/*------------------------------------------------------------*/
/* multi-point evaluation/interpolation                       */
/*------------------------------------------------------------*/
void check()
{
    slong i;
    flint_rand_t state;
    ulong nn;
    nmod_t mod;
    
    flint_rand_init(state);
    nn = 65537;
    nmod_init(&mod, nn);

    for (i = 1000; i < 10000; i+=1000)
    {
        ulong r;
        nmod_geometric_progression_t G;
        nmod_poly_t P, P2;
        nn_ptr val;
        
        r = nmod_find_root(2*i, mod);
        nmod_geometric_progression_init_set(G, r, i, mod);

        nmod_poly_init(P, nn);
        nmod_poly_init(P2, nn);
        nmod_poly_randtest(P, state, i);
        
        val = _nmod_vec_init(i);
        nmod_geometric_progression_evaluate(val, P, G);
        nmod_geometric_progression_interpolate(P2, val, G);
        
        assert (nmod_poly_equal(P, P2));
        
        _nmod_vec_clear(val);
        nmod_poly_clear(P2);
        nmod_poly_clear(P);
        nmod_geometric_progression_clear(G);
    }
    flint_rand_clear(state);
}

/*------------------------------------------------------------*/
/* main just calls check()                                    */
/*------------------------------------------------------------*/
int main(int argc, char **argv){
    check();
    return 0;
}
