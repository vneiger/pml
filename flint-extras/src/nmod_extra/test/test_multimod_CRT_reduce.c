#include <assert.h>
#include <stdlib.h>
#include <gmp.h>
#include <flint/nmod_vec.h>
#include <flint/fmpz.h>

#include "nmod_extra.h"
#include "nmod_vec_extra.h"


/*----------------------------------------------------------------*/
/* reduces N elements of n_bits size modulo num_primes primes     */
/*----------------------------------------------------------------*/
void check_nmod_multimod_CRT_reduce(ulong N, ulong num_primes, ulong n_bits)
{
    flint_rand_t state;
    nmod_multimod_CRT_t CRT; 
    nn_ptr *vec_residues;
    nn_ptr input;
    ulong i, j;
    ulong n;
    nmod_t mod;
    
    flint_rand_init(state);
    n = n_urandint(state, 1L << n_bits);
    nmod_init(&mod, n);
    
    input = _nmod_vec_init(N);
    _nmod_vec_rand(input, state, N, mod);
    
    nmod_multimod_CRT_init(CRT, n, num_primes);

    vec_residues = (nn_ptr *) malloc(num_primes * sizeof(nn_ptr *));
    for (i = 0; i < num_primes; i++)
        vec_residues[i] = _nmod_vec_init(N);

    nmod_multimod_CRT_reduce(vec_residues, input, N, CRT);
    for (i = 0; i < N; i++)
        for (j = 0; j < num_primes; j++)
            if (vec_residues[j][i] != input[i] % CRT->mod_primes[j].n)
            {
                printf("error with i=%lu, j=%lu\n", i, j);
                exit(-1);
            }
        
    _nmod_vec_clear(input);
    for (i = 0; i < num_primes; i++)
        _nmod_vec_clear(vec_residues[i]);
    free(vec_residues);
    nmod_multimod_CRT_clear(CRT);
    flint_rand_clear(state);
}

/*--------------------------------------------------------------*/
/* main calls check                                             */
/*--------------------------------------------------------------*/
int main(int argc, char**argv)
{
    ulong i;

    for (i = 1; i < 500; i++)
    {
        check_nmod_multimod_CRT_reduce(i, 1, 30);
        check_nmod_multimod_CRT_reduce(i, 1, 60);
        check_nmod_multimod_CRT_reduce(i, 2, 30);
        check_nmod_multimod_CRT_reduce(i, 2, 60);
        check_nmod_multimod_CRT_reduce(i, 3, 30);
        check_nmod_multimod_CRT_reduce(i, 3, 60);
        check_nmod_multimod_CRT_reduce(i, 4, 30);
        check_nmod_multimod_CRT_reduce(i, 4, 60);
    }
    return 0;
}
