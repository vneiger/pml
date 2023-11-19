#include <assert.h>
#include <gmp.h>
#include <flint/nmod_vec.h>
#include <flint/fmpz.h>

#include "nmod_vec_extra.h"


/*----------------------------------------------------------------*/
/* CRT in size num_bits, reduced mod n of size p_bits             */
/* does 1002 reductions                                           */
/*----------------------------------------------------------------*/
void check_nmod_multimod_CRT_CRT(ulong num_bits, ulong p_bits)
{
    flint_rand_t state;
    nmod_multimod_CRT_t CRT; 
    mp_ptr *vec_residues;
    mp_ptr output;
    ulong i, j, num_primes, N;
    mp_limb_t n;
    fmpz * input;
    
    N=1006;
    
    flint_randinit(state);

    n = n_urandint(state, 1L << p_bits);
    input = (fmpz *) malloc(N * sizeof(fmpz));

    for (i = 0; i < N; i++)
    {
        fmpz_init(input + i);
        fmpz_randbits(input + i, state, num_bits);
        fmpz_abs(input + i, input + i);
    }
    
    if (num_bits > 4*49)
    {
        printf("not enough primes available\n");
        exit(-1);
    }
    num_primes = 1;
    if (num_bits > 49)
        num_primes = 2;
    if (num_bits > 2*49)
        num_primes = 3;
    if (num_bits > 3*49)
        num_primes = 4;

    nmod_multimod_CRT_init(CRT, n, num_primes);

    vec_residues = (mp_ptr *) malloc(num_primes * sizeof(mp_ptr));
    for (i = 0; i < num_primes; i++)
    {
        vec_residues[i] = _nmod_vec_init(N);
        for (j = 0; j < N; j++)
            vec_residues[i][j] = fmpz_fdiv_ui(input + j, CRT->mod_primes[i].n);
    }

    output = _nmod_vec_init(N);
    nmod_multimod_CRT_CRT(output, vec_residues, N, CRT);

    for (i = 0; i < N; i++)
    {
        if (output[i] != fmpz_fdiv_ui(input + i, n))
        {
            printf("error with i=%lu for %lu bits\n", i, num_bits);
            exit(-1);
        }
    }
    
    _nmod_vec_clear(output);
    
    for (i = 0; i < num_primes; i++)
        _nmod_vec_clear(vec_residues[i]);
    free(vec_residues);

    nmod_multimod_CRT_clear(CRT);

    for (i = 0; i < N; i++)
        fmpz_clear(input + i);
    free(input);
    
    flint_randclear(state);
}

/*--------------------------------------------------------------*/
/* main calls check                                             */
/*--------------------------------------------------------------*/
int main(int argc, char**argv)
{
    check_nmod_multimod_CRT_CRT(30, 25);
    check_nmod_multimod_CRT_CRT(50, 25);
    check_nmod_multimod_CRT_CRT(100, 25);
    check_nmod_multimod_CRT_CRT(150, 25);

    check_nmod_multimod_CRT_CRT(30, 49);
    check_nmod_multimod_CRT_CRT(50, 49);
    check_nmod_multimod_CRT_CRT(100, 49);
    check_nmod_multimod_CRT_CRT(150, 49);

    check_nmod_multimod_CRT_CRT(30, 60);
    check_nmod_multimod_CRT_CRT(50, 60);
    check_nmod_multimod_CRT_CRT(100, 60);
    check_nmod_multimod_CRT_CRT(150, 60);

    return 0;
}
