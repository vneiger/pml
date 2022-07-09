#include <assert.h>
#include <flint/flint.h>

#include "nmod_extra.h"

/*------------------------------------------------------------*/
/* testing high product of packed 32-bit integers             */
/*------------------------------------------------------------*/
void check()
{
#ifdef HAS_AVX512
    slong i;
    mp_hlimb_t p;
    mp_hlimb_t * ar, * br, * cr;
    __m512i a, b, c;
    flint_rand_t state;

    flint_randinit(state);

    ar = (mp_hlimb_t *) aligned_alloc(64, 64);
    br = (mp_hlimb_t *) aligned_alloc(64, 64);
    cr = (mp_hlimb_t *) aligned_alloc(64, 64);

    p = 1 << 30;
    
    for (i = 0; i < 16; i++)
    {
        ar[i] = n_randtest(state) % p;
        br[i] = n_randtest(state) % p;
        cr[i] = n_randtest(state) % p;
    }

    a = _mm512_load_si512(ar);
    b = _mm512_load_si512(br);
    c = mm512_mulhi_epi32(a, b);
    
    _mm512_store_si512((__m512i*) cr, c);

    for (i = 0; i < 16; i++)
    {
        assert(((((mp_limb_t) ar[i]) * ((mp_limb_t) br[i])) >> 32) == (mp_limb_t) cr[i]);
    }
    
    free(ar);
    free(br);
    free(cr);
    
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
