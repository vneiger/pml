#ifndef __FMPZ_EXTRA__H
#define __FMPZ_EXTRA__H

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
/* ------------------------------------------------------------ */
/* ------------------------------------------------------------ */
typedef struct
{
    mp_bitcnt_t prime_bit_length; 
    mp_ptr primes;
    mp_ptr *powers_of_two;
    nmod_t *mod;
    ulong num_primes;
    ulong num_limbs;
}
fmpz_multimod_struct;

typedef fmpz_multimod_struct fmpz_multimod_t[1];

/* ------------------------------------------------------------ */
/* copies primes in mmod-> primes (length is num_primes)        */
/* max_bit_length = maximum bit-length of inputs we will reduce */
/* ------------------------------------------------------------ */
void fmpz_multimod_init(fmpz_multimod_t mmod, mp_srcptr primes, ulong num_primes, slong max_bit_length);

/* ------------------------------------------------------------ */
/* frees all memory in mmod                                     */
/* ------------------------------------------------------------ */
void fmpz_multimod_clear(fmpz_multimod_t mmod);

/* ------------------------------------------------------------ */
/* computes A mod mmod[i] for all i                             */
/* ------------------------------------------------------------ */
void fmpz_multimod_reduce(mp_limb_t * out, const fmpz_t A, const fmpz_multimod_t mmod);


/* ------------------------------------------------------------ */
/* ------------------------------------------------------------ */
/* data structure for Chinese Remaindering                      */
/* ------------------------------------------------------------ */
/* ------------------------------------------------------------ */
typedef struct
{
    fmpz_t prod; /* product of all primes */
    mp_bitcnt_t prime_bit_length; /* common bit-length of all primes */
    mp_ptr primes; 
    mp_ptr inverses;  /* 1/prod_{j ne i} primes[j] mod primes[i] */
    mp_ptr * coefficients; /* coefficients[i][j] = i-th limb of prod_{k ne j} primes[k] */
    mpz_t temp;  /* tmp space for reconstruction */
    nmod_t * mod; /* mod[i] = Z/primes[i]Z */
    ulong num_primes; /* number of primes */
    ulong num_limbs; /* number of limbs of prod */
}
fmpz_CRT_struct;

typedef fmpz_CRT_struct fmpz_CRT_t[1];

/* ------------------------------------------------------------ */
/* prepares the vectors of coefficients and inverses            */
/* ------------------------------------------------------------ */
void fmpz_CRT_init(fmpz_CRT_t mmod, mp_srcptr primes, ulong num_primes);

/* ------------------------------------------------------------ */
/* clears all memory                                            */
/* ------------------------------------------------------------ */
void fmpz_CRT_clear(fmpz_CRT_t mCRT);

/* ------------------------------------------------------------ */
/* compute A = sum_i m[i] \prod_{j \ne i} primes[j]             */
/* ------------------------------------------------------------ */
void fmpz_CRT_combine(fmpz_t A, mp_srcptr m, const fmpz_CRT_t mCRT);

/* ------------------------------------------------------------ */
/* compute A s.t. A mod primes[j] = m[j] for all j              */
/* ------------------------------------------------------------ */
void fmpz_CRT_CRT(fmpz_t A, mp_srcptr m, const fmpz_CRT_t mCRT);


#ifdef __cplusplus
}
#endif


#endif

