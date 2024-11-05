#include <assert.h>
#include <flint/nmod.h>

#include "nmod_mat_extra.h"
 

/*--------------------------------------------------------------*/
/* computes a square matrix product in size len modulo n        */
/*--------------------------------------------------------------*/
void check_nmod_mat_mul_newdot(ulong len, ulong n, flint_rand_t state)
{
    nmod_mat_t a, b, c1, c2;
    nmod_t mod;

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
    nmod_mat_mul_newdot(c2, a, b);

    assert(nmod_mat_equal(c1, c2));
    
    nmod_mat_clear(a);
    nmod_mat_clear(b);
    nmod_mat_clear(c1);
    nmod_mat_clear(c2);
}

/*--------------------------------------------------------------*/
/* main calls check                                             */
/*--------------------------------------------------------------*/
int main()
{
    flint_rand_t state;
    flint_rand_init(state);

    printf("test running, various bitlengths, n x n product from 1 to 400 (no error message means success)...\n");
    for (ulong n = 1; n < 400; n += 10)
    {
        if ((n-1) % 20 == 0)
        {
            printf("%ld..", n);
            fflush(stdout);
        }
        check_nmod_mat_mul_newdot(n, (UWORD(1) << 3) + 1, state);
        check_nmod_mat_mul_newdot(n, (UWORD(1) << 10) + 1, state);
        check_nmod_mat_mul_newdot(n, (UWORD(1) << 20) + 1, state);
        check_nmod_mat_mul_newdot(n, (UWORD(1) << 29) + 1, state);
        check_nmod_mat_mul_newdot(n, (UWORD(1) << 30) + 1, state);
        check_nmod_mat_mul_newdot(n, (UWORD(1) << 31) + 1, state);
        check_nmod_mat_mul_newdot(n, (UWORD(1) << 32) + 1, state);
        check_nmod_mat_mul_newdot(n, (UWORD(1) << 40) + 1, state);
        check_nmod_mat_mul_newdot(n, (UWORD(1) << 50) + 1, state);
        check_nmod_mat_mul_newdot(n, (UWORD(1) << 60) + 1, state);
        check_nmod_mat_mul_newdot(n, UWORD_MAX, state);
    }
    printf("\n");

    flint_rand_clear(state);
    return 0;
}
