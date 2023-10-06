#include <assert.h>
#include <gmp.h>
#include <flint/nmod_vec.h>

#include "fmpz_extra.h"
#include "nmod_vec_extra.h"


/*--------------------------------------------------------------*/
/* creates and deletes a multimod                               */
/* total bit length is max_bit_length                           */
/*--------------------------------------------------------------*/
void check_fmpz_multimod_reduce(ulong max_bit_length)
{
    flint_rand_t state;
    flint_bitcnt_t prime_length;
    fmpz_multimod_t mmod; 
    mp_ptr primes, residues;
    ulong num_primes;
    fmpz_comb_t C;
    fmpz_comb_temp_t Ct;
    fmpz_t A;
    double t, tp;
    clock_t tt;
    long nb_iter;

    flint_randinit(state);
    prime_length = 59;
    num_primes = 1 + (max_bit_length / prime_length);
    
    primes = _nmod_vec_init(num_primes);
    residues = _nmod_vec_init(num_primes);
    nmod_vec_primes(primes, num_primes, prime_length);
    fmpz_multimod_init(mmod, primes, num_primes, max_bit_length);

    printf("%ld ", num_primes);

    fmpz_comb_init(C, primes, num_primes);
    fmpz_comb_temp_init(Ct, C);

    fmpz_init(A);
    fmpz_randtest(A, state, max_bit_length);

    t = 0.0;
    nb_iter = 0;
    while (t < 0.5)
    {
        tt = clock();
        fmpz_multimod_reduce(residues, A, mmod);
        t += (double)(clock()-tt) / CLOCKS_PER_SEC;
        ++nb_iter;
    }
    t = 1000 * t;
    t /= nb_iter;
    printf("%lf ", t);

    tp = 0.0;
    nb_iter = 0;
    while (tp < 0.5)
    {
        tt = clock();
        fmpz_multi_mod_ui(residues, A, C, Ct);
        tp += (double)(clock()-tt) / CLOCKS_PER_SEC;
        ++nb_iter;
    }
    tp = 1000 * tp;
    tp /= nb_iter;
    printf("%lf %lf ", tp, tp / t);

    printf("\n");
    
    fmpz_clear(A);
    fmpz_comb_temp_clear(Ct);
    fmpz_comb_clear(C);
    fmpz_multimod_clear(mmod);
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

    printf("# num_primes t_new t_old t_old/t_new\n");
    for (i = 10; i < 100000; i += 200)
	check_fmpz_multimod_reduce(i);

    return 0;
}
