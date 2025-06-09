#include <time.h>
#include <gmp.h>
#include <flint/nmod_vec.h>

#include "fmpz_extra.h"
#include "nmod_vec_extra.h"


/*--------------------------------------------------------------*/
/* creates and deletes a multimod                               */
/* uses num_primes of size n_bits                               */
/* reduces a random integer A modulo all primes                 */
/* A initially reduced modulo the producto of primes            */
/*--------------------------------------------------------------*/
void check_fmpz_multimod_naive_reduce(ulong num_primes, ulong n_bits)
{
    flint_rand_t state;
    fmpz_multimod_naive_t mmod; 
    nn_ptr primes, residues;
    fmpz_comb_t C;
    fmpz_comb_temp_t Ct;
    fmpz_t A;
    double t, tp;
    clock_t tt;
    long nb_iter;

    flint_rand_init(state);
    
    primes = _nmod_vec_init(num_primes);
    residues = _nmod_vec_init(num_primes);
    nmod_vec_primes(primes, num_primes, n_bits);
    fmpz_multimod_naive_init(mmod, primes, num_primes);

    printf("%lu %lu ", num_primes, n_bits);

    fmpz_comb_init(C, primes, num_primes);
    fmpz_comb_temp_init(Ct, C);
    fmpz_init(A);
    fmpz_randbits(A, state, fmpz_bits(mmod->prod) - 1);

    t = 0.0;
    nb_iter = 0;
    while (t < 0.5)
    {
        tt = clock();
        fmpz_multimod_naive_reduce(residues, A, mmod);
        fmpz_multimod_naive_reduce(residues, A, mmod);
        fmpz_multimod_naive_reduce(residues, A, mmod);
        fmpz_multimod_naive_reduce(residues, A, mmod);
        fmpz_multimod_naive_reduce(residues, A, mmod);
        fmpz_multimod_naive_reduce(residues, A, mmod);
        fmpz_multimod_naive_reduce(residues, A, mmod);
        fmpz_multimod_naive_reduce(residues, A, mmod);
        fmpz_multimod_naive_reduce(residues, A, mmod);
        fmpz_multimod_naive_reduce(residues, A, mmod);
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
        fmpz_multi_mod_ui(residues, A, C, Ct);
        fmpz_multi_mod_ui(residues, A, C, Ct);
        fmpz_multi_mod_ui(residues, A, C, Ct);
        fmpz_multi_mod_ui(residues, A, C, Ct);
        fmpz_multi_mod_ui(residues, A, C, Ct);
        fmpz_multi_mod_ui(residues, A, C, Ct);
        fmpz_multi_mod_ui(residues, A, C, Ct);
        fmpz_multi_mod_ui(residues, A, C, Ct);
        fmpz_multi_mod_ui(residues, A, C, Ct);
        fmpz_multi_mod_ui(residues, A, C, Ct);
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
    fmpz_multimod_naive_clear(mmod);
    _nmod_vec_clear(residues);
    _nmod_vec_clear(primes);
    flint_rand_clear(state);
}

/*--------------------------------------------------------------*/
/* main calls check                                             */
/*--------------------------------------------------------------*/
int main(int argc, char **argv)
{
    ulong i;

    printf("# num_primes n_bits t_new t_old\n");
    for (i = 1; i < 100; i += 5)
        check_fmpz_multimod_naive_reduce(i, 50);

    for (i = 1; i < 100; i += 5)
        check_fmpz_multimod_naive_reduce(i, 29);

    return 0;
}
