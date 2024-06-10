#include <assert.h>
#include <gmp.h>
#include <flint/nmod_vec.h>
#include <flint/fmpz.h>

#include "nmod_vec_extra.h"


/*----------------------------------------------------------------*/
/* CRT in size num_bits, reduced mod n of size n_bits             */
/* does 1000 reductions                                           */
/*----------------------------------------------------------------*/
void check_nmod_multimod_CRT_CRT(ulong num_bits, ulong n_bits)
{
    flint_rand_t state;
    nmod_multimod_CRT_t CRT;
    nn_ptr *vec_residues;
    nn_ptr output;
    ulong i, j, num_primes, N, nb_iter;
    ulong n;
    fmpz * input;
    clock_t tt;
    double t;
    
    N=1000;
    
    flint_rand_init(state);

    n = n_urandint(state, 1L << n_bits);
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
        
    vec_residues = (nn_ptr *) malloc(num_primes * sizeof(nn_ptr));
    for (i = 0; i < num_primes; i++)
    {
        vec_residues[i] = _nmod_vec_init(N);
        for (j = 0; j < N; j++)
            vec_residues[i][j] = fmpz_fdiv_ui(input + j, CRT->primes[i]);
    }
    
    output = _nmod_vec_init(N);

    printf("#%lu primes, %lu bits\n", num_primes, n_bits);
    
    t = 0.0;
    nb_iter = 0;
    while (t < 0.5)
    {
        tt = clock();
        nmod_multimod_CRT_CRT(output, vec_residues, N, CRT);    
        t += (double)(clock()-tt) / CLOCKS_PER_SEC;
        ++nb_iter;
    }
    t = 1000 * t;
    t /= nb_iter;
    printf("%f ", t);

    printf("\n");
        
    _nmod_vec_clear(output);
    for (i = 0; i < num_primes; i++)
        _nmod_vec_clear(vec_residues[i]);
    free(vec_residues);
    nmod_multimod_CRT_clear(CRT);
    for (i = 0; i < N; i++)
        fmpz_clear(input + i);
    free(input);
    flint_rand_clear(state);
}

/*--------------------------------------------------------------*/
/* main calls check                                             */
/*--------------------------------------------------------------*/
int main(int argc, char**argv)
{
    check_nmod_multimod_CRT_CRT(30, 30);
    check_nmod_multimod_CRT_CRT(30, 60);
    check_nmod_multimod_CRT_CRT(50, 30);
    check_nmod_multimod_CRT_CRT(50, 60);
    check_nmod_multimod_CRT_CRT(100, 30);
    check_nmod_multimod_CRT_CRT(100, 60);
    check_nmod_multimod_CRT_CRT(150, 30);
    check_nmod_multimod_CRT_CRT(150, 60);

    return 0;
}
