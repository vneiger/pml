#include <assert.h>
#include <gmp.h>
#include <flint/nmod_vec.h>
#include <flint/fmpz.h>

#include "nmod_extra.h"
#include "nmod_vec_extra.h"


/*----------------------------------------------------------------*/
/* reduces 1002 elements of n_bits size modulo num_primes primes  */
/*----------------------------------------------------------------*/
void check_nmod_multimod_CRT_reduce(ulong num_primes, ulong n_bits)
{
    flint_rand_t state;
    nmod_multimod_CRT_t CRT; 
    mp_ptr *vec_residues;
    mp_ptr input;
    ulong i, j, N;
    mp_limb_t n;
    nmod_t mod;
    
    N=1002;
    
    flint_randinit(state);

    n = n_urandint(state, 1L << n_bits);
    nmod_init(&mod, n);
    input = _nmod_vec_init(N);

    _nmod_vec_rand(input, state, N, mod);
    
    nmod_multimod_CRT_init(CRT, n, num_primes);

    vec_residues = (mp_ptr *) malloc(num_primes * sizeof(mp_ptr));
    for (i = 0; i < num_primes; i++)
        vec_residues[i] = _nmod_vec_init(N);

    nmod_multimod_CRT_reduce(vec_residues, input, N, CRT);
   
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < num_primes; j++)
            if (vec_residues[j][i] != input[i] % CRT->mod_primes[j].n)
            {
                printf("error with i=%lu, j=%lu\n", i, j);
                exit(-1);
            }
    }
        
    _nmod_vec_clear(input);
    
    for (i = 0; i < num_primes; i++)
        _nmod_vec_clear(vec_residues[i]);
    free(vec_residues);

    nmod_multimod_CRT_clear(CRT);
    flint_randclear(state);
}

/*--------------------------------------------------------------*/
/* main calls check                                             */
/*--------------------------------------------------------------*/
int main(int argc, char**argv)
{
    check_nmod_multimod_CRT_reduce(1, 30);
    check_nmod_multimod_CRT_reduce(1, 60);
    check_nmod_multimod_CRT_reduce(2, 30);
    check_nmod_multimod_CRT_reduce(2, 60);
    check_nmod_multimod_CRT_reduce(3, 30);
    check_nmod_multimod_CRT_reduce(3, 60);
    check_nmod_multimod_CRT_reduce(4, 30);
    check_nmod_multimod_CRT_reduce(4, 60);
    return 0;
}
