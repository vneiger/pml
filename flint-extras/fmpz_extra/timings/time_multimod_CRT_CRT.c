#include <assert.h>
#include <gmp.h>
#include <flint/nmod_vec.h>

#include "fmpz_extra.h"
#include "nmod_vec_extra.h"


/*--------------------------------------------------------------*/
/* creates and deletes a multimod, does CRT                     */
/* 50 bit primes, total bit length is max_bit_length            */
/*--------------------------------------------------------------*/
void check_fmpz_multimod_CRT(ulong max_bit_length)
{
    flint_rand_t state;
    flint_bitcnt_t prime_length;
    fmpz_CRT_naive_t mmod; 
    fmpz_multimod_CRT_t mmod_C; 
    mp_ptr primes, residues;
    ulong num_primes;
    fmpz_comb_t C;
    fmpz_comb_temp_t Ct;
    fmpz_t A;
    double t, tp;
    clock_t tt;
    long nb_iter;

    flint_randinit(state);
    prime_length = 50;
    num_primes = 1 + (max_bit_length / prime_length);
    
    primes = _nmod_vec_init(num_primes);
    residues = _nmod_vec_init(num_primes);
    nmod_vec_primes(primes, num_primes, prime_length);

    printf("%ld ", num_primes);
    
    fmpz_CRT_naive_init(mmod, primes, num_primes);
    fmpz_multimod_CRT_init(mmod_C, primes, num_primes);
    fmpz_comb_init(C, primes, num_primes);
    fmpz_comb_temp_init(Ct, C);

    fmpz_init(A);
    fmpz_randbits(A, state, max_bit_length);

    t = 0.0;
    nb_iter = 0;
    while (t < 0.5)
    {
        tt = clock();
        fmpz_CRT_naive_CRT(A, residues, mmod);
        fmpz_CRT_naive_CRT(A, residues, mmod);
        fmpz_CRT_naive_CRT(A, residues, mmod);
        fmpz_CRT_naive_CRT(A, residues, mmod);
        fmpz_CRT_naive_CRT(A, residues, mmod);
        fmpz_CRT_naive_CRT(A, residues, mmod);
        fmpz_CRT_naive_CRT(A, residues, mmod);
        fmpz_CRT_naive_CRT(A, residues, mmod);
        fmpz_CRT_naive_CRT(A, residues, mmod);
        fmpz_CRT_naive_CRT(A, residues, mmod);
        t += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += 10;
    }
    t = 1000 * t;
    t /= nb_iter;
    printf("%lf ", t);

    
    t = 0.0;
    nb_iter = 0;
    while (t < 0.5)
    {
        tt = clock();
        fmpz_multimod_CRT_CRT(A, residues, mmod_C);
        fmpz_multimod_CRT_CRT(A, residues, mmod_C);
        fmpz_multimod_CRT_CRT(A, residues, mmod_C);
        fmpz_multimod_CRT_CRT(A, residues, mmod_C);
        fmpz_multimod_CRT_CRT(A, residues, mmod_C);
        fmpz_multimod_CRT_CRT(A, residues, mmod_C);
        fmpz_multimod_CRT_CRT(A, residues, mmod_C);
        fmpz_multimod_CRT_CRT(A, residues, mmod_C);
        fmpz_multimod_CRT_CRT(A, residues, mmod_C);
        fmpz_multimod_CRT_CRT(A, residues, mmod_C);
        t += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += 10;
    }
    t = 1000 * t;
    t /= nb_iter;
    printf("%lf ", t);

    
    tp = 0.0;
    nb_iter = 0;
    while (tp < 0.5)
    {
        tt = clock();
        fmpz_multi_CRT_ui(A, residues, C, Ct, 0);
        fmpz_multi_CRT_ui(A, residues, C, Ct, 0);
        fmpz_multi_CRT_ui(A, residues, C, Ct, 0);
        fmpz_multi_CRT_ui(A, residues, C, Ct, 0);
        fmpz_multi_CRT_ui(A, residues, C, Ct, 0);
        fmpz_multi_CRT_ui(A, residues, C, Ct, 0);
        fmpz_multi_CRT_ui(A, residues, C, Ct, 0);
        fmpz_multi_CRT_ui(A, residues, C, Ct, 0);
        fmpz_multi_CRT_ui(A, residues, C, Ct, 0);
        fmpz_multi_CRT_ui(A, residues, C, Ct, 0);
        tp += (double)(clock()-tt) / CLOCKS_PER_SEC;
        nb_iter += 10;
    }
    tp = 1000 * tp;
    tp /= nb_iter;
    printf("%lf ", tp);

    printf("\n");
    
    fmpz_clear(A);
    fmpz_comb_temp_clear(Ct);
    fmpz_comb_clear(C);
    fmpz_multimod_CRT_clear(mmod_C);
    fmpz_CRT_naive_clear(mmod);
    _nmod_vec_clear(residues);
    _nmod_vec_clear(primes);
    flint_randclear(state);
}

/*--------------------------------------------------------------*/
/* main calls check                                             */
/*--------------------------------------------------------------*/
int main(int argc, char **argv)
{
    ulong i;

    printf("# num_primes t_naive t_DAC t_old\n");
    for (i = 2000; i < 100000; i += 2000000)
	check_fmpz_multimod_CRT(i);

    return 0;
}
