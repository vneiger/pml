#ifndef __FMPZ_EXTRA__H
#define __FMPZ_EXTRA__H

#include <gmp.h> // TODO for mpz_t temp; could probably use fmpz instead?
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/nmod.h>
#include <flint/nmod_vec.h>
#include <flint/machine_vectors.h>

#ifdef __cplusplus
extern "C" {
#endif


/* ------------------------------------------------------------ */
/* ------------------------------------------------------------ */
/* data structure for multi-modular reduction                   */
/* naive algorithms, word-size moduli                           */
/* ------------------------------------------------------------ */
/* ------------------------------------------------------------ */
typedef struct
{
    flint_bitcnt_t prime_bit_length; 
    nn_ptr primes;
    nn_ptr *powers_of_two;
    nmod_t *mod;
    ulong num_primes;
    ulong num_limbs;
    fmpz_t prod;
    int small_moduli;
}
fmpz_multimod_naive_struct;

typedef fmpz_multimod_naive_struct fmpz_multimod_naive_t[1];

/* ------------------------------------------------------------ */
/* copies primes in mmod-> primes (length is num_primes)        */
/* ------------------------------------------------------------ */
void fmpz_multimod_naive_init(fmpz_multimod_naive_t mmod, nn_srcptr primes, ulong num_primes);

/* ------------------------------------------------------------ */
/* frees all memory in mmod                                     */
/* ------------------------------------------------------------ */
void fmpz_multimod_naive_clear(fmpz_multimod_naive_t mmod);

/* ------------------------------------------------------------ */
/* computes A mod mmod[i] for all i                             */
/* ------------------------------------------------------------ */
void fmpz_multimod_naive_reduce(ulong * out, const fmpz_t A, const fmpz_multimod_naive_t mmod);

/* ------------------------------------------------------------ */
/* ------------------------------------------------------------ */
/* data structure for Chinese Remaindering                      */
/* naive algorithms, word-size moduli                           */
/* ------------------------------------------------------------ */
/* ------------------------------------------------------------ */
typedef struct
{
    fmpz_t prod; /* product of all primes */
    flint_bitcnt_t prime_bit_length; /* common bit-length of all primes */
    nn_ptr primes; 
    nn_ptr inverses;  /* 1/prod_{j ne i} primes[j] mod primes[i] */
    nn_ptr * coefficients; /* coefficients[i][j] = i-th limb of prod_{k ne j} primes[k] */
    mpz_t temp;  /* tmp space for reconstruction */
    nmod_t * mod; /* mod[i] = Z/primes[i]Z */
    ulong num_primes; /* number of primes */
    ulong num_limbs; /* number of limbs of prod */
}
fmpz_CRT_naive_struct;

typedef fmpz_CRT_naive_struct fmpz_CRT_naive_t[1];

/* ------------------------------------------------------------ */
/* prepares the vectors of coefficients and inverses            */
/* ------------------------------------------------------------ */
void fmpz_CRT_naive_init(fmpz_CRT_naive_t mmod, nn_srcptr primes, ulong num_primes);

/* ------------------------------------------------------------ */
/* clears all memory                                            */
/* ------------------------------------------------------------ */
void fmpz_CRT_naive_clear(fmpz_CRT_naive_t mCRT);

/* ------------------------------------------------------------ */
/* compute A = sum_i m[i] \prod_{j \ne i} primes[j]             */
/* ------------------------------------------------------------ */
void fmpz_CRT_naive_combine(fmpz_t A, nn_srcptr m, const fmpz_CRT_naive_t mCRT);

/* ------------------------------------------------------------ */
/* compute A s.t. A mod primes[j] = m[j] for all j              */
/* ------------------------------------------------------------ */
void fmpz_CRT_naive_CRT(fmpz_t A, nn_srcptr m, const fmpz_CRT_naive_t mCRT);


/* ------------------------------------------------------------ */
/* ------------------------------------------------------------ */
/* data structure for multimod and Chinese Remaindering         */
/* subproduct trees + naive algorithms, word-size moduli        */
/* ------------------------------------------------------------ */
/* ------------------------------------------------------------ */
typedef struct
{
    ulong num_primes;
    ulong size_leaves;
    ulong num_leaves;
    nn_ptr inverse_cofactors;
    fmpz_multimod_naive_t * leaves_mod;
    fmpz_CRT_naive_t * leaves_CRT;
    fmpz * products_leaves;
    fmpz_multi_mod_t top_mod;
    fmpz_multi_CRT_t top_CRT;
    fmpz_t product_primes;
}
fmpz_multimod_CRT_struct;

typedef fmpz_multimod_CRT_struct fmpz_multimod_CRT_t[1];

#define MULTIMOD_CRT_LEAF_SIZE 200

/* ------------------------------------------------------------ */
/* prepares the vectors of coefficients and inverses            */
/* ------------------------------------------------------------ */
void fmpz_multimod_CRT_init(fmpz_multimod_CRT_t mmod, nn_srcptr primes, ulong num_primes);

/* ------------------------------------------------------------ */
/*  clears all memory used by mmod                              */
/* ------------------------------------------------------------ */
void fmpz_multimod_CRT_clear(fmpz_multimod_CRT_t mmod);

/* ------------------------------------------------------------ */
/* computes A mod mmod[i] for all i                             */
/* ------------------------------------------------------------ */
void fmpz_multimod_CRT_reduce(ulong * out, const fmpz_t A, const fmpz_multimod_CRT_t mmod);

/* ------------------------------------------------------------ */
/* compute A s.t. A mod primes[j] = m[j] for all j              */
/* ------------------------------------------------------------ */
void fmpz_multimod_CRT_CRT(fmpz_t A, nn_srcptr m, const fmpz_multimod_CRT_t mmod);


/*-------------------------------------------------------------------*/
/* helper function                                                   */
/* similar to fmpz_multi_CRT_precomp, except that we do not multiply */
/* by cofactors at the leaves                                        */
/* result is still reduced modulo the product of moduli              */
/*-------------------------------------------------------------------*/
void fmpz_multi_CRT_combine(fmpz_t output, const fmpz_multi_CRT_t P, fmpz * inputs);



#ifdef __cplusplus
}
#endif


#endif

