#include <assert.h>
#include <gmp.h>
#include <flint/nmod_vec.h>

#include "fmpz_extra.h"
#include "nmod_vec_extra.h"


/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
void check_fmpz_CRT_naive_CRT(ulong num_primes)
{
    flint_rand_t state;
    fmpz_t comb;
    fmpz_multimod_naive_t mmod; 
    fmpz_CRT_naive_t mCRT; 
    fmpz_comb_t C;
    fmpz_comb_temp_t Ct;
    nn_ptr primes, residues;
    ulong i;
    double t, tp;
    clock_t tt;
    long nb_iter;

    flint_rand_init(state);
    primes = _nmod_vec_init(num_primes);
    residues = _nmod_vec_init(num_primes);
    fmpz_init(comb);

    nmod_vec_primes(primes, num_primes, FLINT_BITS-8);

    // used to time remainder
    fmpz_multimod_naive_init(mmod, primes, num_primes);
    
    for (i = 0; i < num_primes; i++)
	residues[i] = n_randlimb(state) % primes[i];

    fmpz_CRT_naive_init(mCRT, primes, num_primes);

    printf("%ld ", num_primes);

    t = 0.0;
    nb_iter = 0;
    while (t < 0.5)
    {
        tt = clock();
        fmpz_CRT_naive_CRT(comb, residues, mCRT);
        t += (double)(clock()-tt) / CLOCKS_PER_SEC;
        ++nb_iter;
    }
    t = 1000 * t;
    t /= nb_iter;
    printf("%lf ", t);

    fmpz_comb_init(C, primes, num_primes);
    fmpz_comb_temp_init(Ct, C);

    tp = 0.0;
    nb_iter = 0;
    while (tp < 0.5)
    {
        tt = clock();
        fmpz_multi_CRT_ui(comb, residues, C, Ct, 0);
        tp += (double)(clock()-tt) / CLOCKS_PER_SEC;
        ++nb_iter;
    }
    tp = 1000 * tp;
    tp /= nb_iter;
    printf("%lf %lf ", tp, tp / t);

    t = 0.0;
    nb_iter = 0;
    while (t < 0.5)
    {
        tt = clock();
        fmpz_multimod_naive_reduce(residues, comb, mmod);
        t += (double)(clock()-tt) / CLOCKS_PER_SEC;
        ++nb_iter;
    }
    t = 1000 * t;
    t /= nb_iter;
    printf("%lf ", t);


    printf("\n");
    
    fmpz_CRT_naive_clear(mCRT);
    fmpz_clear(comb);
    fmpz_multimod_naive_clear(mmod);
    _nmod_vec_clear(residues);
    _nmod_vec_clear(primes);
    flint_rand_clear(state);
}

/*--------------------------------------------------------------*/
/* main calls check                                             */
/*--------------------------------------------------------------*/
int main(int argc, char**argv)
{
    ulong i;

    printf("# num_primes t_new t_old t_old/t_new t_multimod\n");
    for (i = 1; i < 1000; i += 10)
	check_fmpz_CRT_naive_CRT(i);

    return 0;
}
