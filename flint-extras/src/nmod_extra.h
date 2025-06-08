#ifndef __NMOD_EXTRA__H
#define __NMOD_EXTRA__H

#include <flint/flint.h>
#include <flint/nmod.h>
#include <flint/nmod_vec.h>
#include <flint/mpn_extras.h>
#include <flint/longlong.h>
#include <flint/ulong_extras.h>

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*        A few extra functionalities for Fp                  */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* returns the smallest i such that 2^i >= x                  */
/*------------------------------------------------------------*/
int next_power_of_two(ulong x);

/*------------------------------------------------------------*/
/* returns 1/p mod 2^k, assuming p is odd                     */
/* ill-defined when p is even                                 */
/*------------------------------------------------------------*/
ulong inverse_mod_power_of_two(ulong p, int k);

/*------------------------------------------------------------*/
/* finds an element of order at least n                       */
/* returns 0 if not found                                     */
/*------------------------------------------------------------*/
ulong nmod_find_root(long n, nmod_t mod);



/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* 64 bit modular multiplication using a preconditionner      */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/* TODO unused for the moment:
 * - what does this bring compared to n_mulmod_shoup in FLINT?
 * - the lazy/unreduced variant could be defined locally in relevant files like FFT ones (?))
 */

/*------------------------------------------------------------*/
/* high word of a*b, assuming words are 64 bits               */
/*------------------------------------------------------------*/
FLINT_FORCE_INLINE ulong mul_hi(ulong a, ulong b)
{
    ulong p_hi, p_lo;
    umul_ppmm(p_hi, p_lo, a, b);
    return p_hi;
}

/*------------------------------------------------------------*/
/* returns floor (2^64*b)/p                                   */
/*------------------------------------------------------------*/
FLINT_FORCE_INLINE ulong prep_mul_mod_precon(ulong b, ulong p)
{
    return n_mulmod_precomp_shoup(b, p);
}

/*------------------------------------------------------------*/
/* preconditioned product                                     */
/* returns ab mod p                                           */
/* a is in [0..2^64), b is in [0..p)                          */
/*------------------------------------------------------------*/
FLINT_FORCE_INLINE ulong mul_mod_precon(ulong a, ulong b, ulong p, ulong i_b)
{
    return n_mulmod_shoup(b, a, i_b, p);
}

/*------------------------------------------------------------*/
/* returns ab mod p in [0..2p)                                */
/* a is in [0..2^64), b is in [0..p)                          */
/*------------------------------------------------------------*/
FLINT_FORCE_INLINE ulong mul_mod_precon_unreduced(ulong a, ulong b, ulong p, ulong i_b)
{
    ulong q;

    q = mul_hi(i_b, a);
    return a*b - q*p;
}

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*                    CRT and multimod                        */
/* -we only handle 1 to 4 hard-coded primes p0..pk            */
/* -we take as extra parameter another modulus p              */           
/* -multimod reduces vectors of ulong mod all pi's        */
/* -CRT does Chinese Remainders, and reduces the result mod p */
/* -use AVX2 for p < 2^50, long arithmetic otherwise          */
/* TODO: set a macro for testing large/small modulus          */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

// 4 primes should be enough for all practical purposes involving nmods
// TODO: make sure we test if so, just in case?
//
// PRIME0 < PRIME1 < PRIME2 < PRIME3 needed for the double implementation
#define PRIME0 659706976665601
#define PRIME1 910395627798529
#define PRIME2 1086317488242689
#define PRIME3 1108307720798209

typedef struct {
    ulong num_primes;
    ulong p;
    nmod_t mod;

    nn_ptr data;
    double pinv;
    nmod_t mod_primes[4];
    ulong primes[4];
    double primes_inv[4];
    ulong inverse_cofactors[4];
    double p0_red, p1_red, p0p1_red, p0p1p2_red, p0_red2, p0p1red_2, invp0_p1, invp0p1_p2, p0p1_red3, invp0p1p2_p3;
} nmod_multimod_CRT_struct;

typedef nmod_multimod_CRT_struct nmod_multimod_CRT_t[1];


/*------------------------------------------------------------*/
/* initializes all data in C                                  */
/*------------------------------------------------------------*/
void nmod_multimod_CRT_init(nmod_multimod_CRT_t C, ulong p, ulong num_primes);

/*------------------------------------------------------------*/
/* clears all data in C                                       */
/*------------------------------------------------------------*/
void nmod_multimod_CRT_clear(nmod_multimod_CRT_t C);

/*------------------------------------------------------------*/
/* out[i] = CRT(residues[j][i], j < num_primes) mod p, i < nb */
/*------------------------------------------------------------*/
void nmod_multimod_CRT_CRT(nn_ptr out, nn_ptr *residues, ulong nb, nmod_multimod_CRT_t C);

/*------------------------------------------------------------*/
/* residues[j][i] = input[i] mod prime[j]                     */
/* for i < nb, j < num_primes                                 */
/*------------------------------------------------------------*/
void nmod_multimod_CRT_reduce(nn_ptr *residues, nn_ptr input, ulong nb, nmod_multimod_CRT_t C);

#endif
