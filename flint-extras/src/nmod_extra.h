/*
    Copyright (C) 2025 Vincent Neiger, Éric Schost

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#ifndef __NMOD_EXTRA__H
#define __NMOD_EXTRA__H

#include <flint/flint.h>
#include <flint/nmod.h>
#include <flint/nmod_vec.h>
#include <flint/mpn_extras.h>
#include <flint/longlong.h>
#include <flint/ulong_extras.h>
#include <stdint.h>

#include "pml.h"

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
ulong nmod_find_root(ulong n, nmod_t mod);


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
// TODO make sure we test if so, just in case?
//
// PRIME0 < PRIME1 < PRIME2 < PRIME3 needed for the double implementation
#define PRIME0 UWORD(659706976665601)
#define PRIME1 UWORD(910395627798529)
#define PRIME2 UWORD(1086317488242689)
#define PRIME3 UWORD(1108307720798209)

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


/*******************************************
*  modular reduction with precomputation  *
*******************************************/

/* general */
typedef struct
{
   ulong n;
   ulong redp;
   ulong shift;
}
nmod_redp_t;

/* FLINT_BITS == 64, nbits < 32, storage uint64 */
typedef struct
{
   ulong n;
   ulong redp;
   ulong shift;
}
nmod32_redp_t;

/* nbits < 32, storage uint32 */
typedef struct
{
   uint32_t n;
   uint32_t redp;
   uint32_t shift;
}
nmod32s_redp_t;

/**************************
*  ulong_extras, general  *
**************************/

/** Precomputation for modular reduction, general
 * (supports up to FLINT_BITS in the number to reduce)
 * Requirements: nbits == nbits(n) > 0
 * Guarantee: reduced number is in [0,n+1]; for some n it is in [0,n]
 * (TODO say more precisely which n have this)
 * NOTE first row 64|nu|64 of table
 */

ULONG_EXTRAS_INLINE void n_red_precomp(ulong * redp, ulong * shift, ulong n, ulong nbits)
{
    FLINT_ASSERT(nbits > 0 && nbits <= FLINT_BITS);
    *shift = nbits + FLINT_BITS - 1;
    ulong FLINT_SET_BUT_UNUSED(rem);
    udiv_qrnnd(*redp, rem, (UWORD(1) << (nbits - 1)), UWORD(0), n);
}

/** Modular reduction, given precomputation
 * Requirements: those of n_red_precomp
 * Output: in [0,n] or [0,n+1] (depends on n)
 */
ULONG_EXTRAS_INLINE ulong n_mod_redp_fast(ulong a, ulong n, ulong redp, ulong shift)
{
    ulong q, tmp;
    umul_ppmm(q, tmp, a, redp);
    return a - (q >> shift) * n;
}

ULONG_EXTRAS_INLINE ulong n_mod_redp(ulong a, ulong n, ulong redp, ulong shift)
{
    /* TODO call n_mod_redp_fast */
    ulong q, tmp;
    umul_ppmm(q, tmp, a, redp);
    tmp = a - (q >> shift) * n;
    if (tmp >= n)
        tmp -= n;
    return tmp;
}

/**********************************
*  ulong_extras, 64/32, uint64_t  *
***********************************/

/** Precomputation for modular reduction, general
 * (supports up to 32 bits in the number to reduce)
 * Requirements: nbits == nbits(n) > 0, FLINT_BITS == 64
 * Guarantee: reduced number is in [0,n+1]; for some n it is in [0,n]
 * (TODO say more precisely which n have this)
 * NOTE first row 64|nu|32 of table
 */

ULONG_EXTRAS_INLINE void n32_red_precomp(ulong * redp, ulong * shift, ulong n, ulong nbits)
{
    FLINT_ASSERT(nbits > 0 && nbits <= FLINT_BITS/2);
    *shift = nbits + 31;
    *redp = (UWORD(1) << (nbits + 31)) / n;
}

/** Modular reduction, given precomputation
 * Requirements: those of n_red_precomp
 * Output: in [0,n] or [0,n+1] (depends on n)
 */
ULONG_EXTRAS_INLINE ulong n32_mod_redp_fast(ulong a, ulong n, ulong redp, ulong shift)
{
    ulong q = (a * redp) >> shift;
    return a - q * n;
}

ULONG_EXTRAS_INLINE ulong n32_mod_redp(ulong a, ulong n, ulong redp, ulong shift)
{
    ulong a_red = n32_mod_redp_fast(a, n, redp, shift);
    if (a_red >= n)
        a_red -= n;
    return a_red;
}


/**********************************
*  ulong_extras, 64/32, uint32_t  *
***********************************/

/** Precomputation for modular reduction, general
 * (supports up to 32 bits in the number to reduce)
 * Requirements: nbits == nbits(n) > 0, FLINT_BITS == 64
 * Guarantee: reduced number is in [0,n+1]; for some n it is in [0,n]
 * (TODO say more precisely which n have this)
 * NOTE first row 64|nu|32 of table
 * TODO is the (u32) cast necessary?
 */

ULONG_EXTRAS_INLINE void n32s_red_precomp(uint32_t * redp, uint32_t * shift, uint32_t n, uint32_t nbits)
{
    FLINT_ASSERT(nbits > 0 && nbits <= FLINT_BITS/2);
    *shift = nbits + 31;
    *redp = (UWORD(1) << (nbits + 31)) / n;
}

/** Modular reduction, given precomputation
 * Requirements: those of n_red_precomp
 * Output: in [0,n] or [0,n+1] (depends on n)
 */
ULONG_EXTRAS_INLINE uint32_t n32s_mod_redp_fast(uint32_t a, uint32_t n, uint32_t redp, uint32_t shift)
{
    uint32_t q = ((ulong)a * redp) >> shift;
    return a - q * n;
}

ULONG_EXTRAS_INLINE uint32_t n32s_mod_redp(uint32_t a, uint32_t n, uint32_t redp, uint32_t shift)
{
    uint32_t a_red = n32s_mod_redp_fast(a, n, redp, shift);
    if (a_red >= n)
        a_red -= n;
    return a_red;
}




/**********
*  nmod  *
**********/

NMOD_INLINE void nmod_redp_init(nmod_redp_t * mod, ulong n)
{
    mod->n = n;
    ulong nbits = FLINT_BITS - flint_clz(n);
    n_red_precomp(&mod->redp, &mod->shift, n, nbits);  /* selection to be done here? */
}

NMOD_INLINE void nmod32_redp_init(nmod32_redp_t * mod, ulong n)
{
    mod->n = n;
    ulong nbits = FLINT_BITS - flint_clz(n);
    n32_red_precomp(&mod->redp, &mod->shift, n, nbits);  /* selection to be done here? */
}

NMOD_INLINE void nmod32s_redp_init(nmod32s_redp_t * mod, uint32_t n)
{
    mod->n = n;
    uint32_t nbits = FLINT_BITS - flint_clz(n);
    n32s_red_precomp(&mod->redp, &mod->shift, n, nbits);  /* selection to be done here? */
}




/**********************
*  !!experimental!!  *
**********************/

/** Modular reduction, given precomputation
 * Requirements: those of n_red_precomp, and nbits >= 33
 * Output: ??
 * TODO experimental!!
 * based on 64 | nu | 64 of table, when nu >= 33
 */
ULONG_EXTRAS_INLINE ulong n_mod_redp_lazy(ulong a, ulong n, ulong redp, ulong shift)
{
    ulong q, a_red;
    q = ((a >> 32) * (redp >> 32)) >> (shift - 64);
    a_red = a - q * n;
    return a_red;
}

ULONG_EXTRAS_INLINE ulong n_mod_redp_lazy_correct(ulong a, ulong n, ulong redp, ulong shift)
{
    ulong q, a_red;
    q = ((a >> 32) * (redp >> 32)) >> (shift - 64);
    a_red = a - q * n;
    if (a_red >= n)
        a_red -= n;
    return a_red;
}

/** Modular multiplication, given precomputation */
/* assumes a<n, b<n (or close to this) */
/* nbits == number of bits in n */
/* ULONG_EXTRAS_INLINE ulong n_mod_mulmod_redp_lazy(ulong a, ulong b, ulong n, ulong redp, ulong nbits) */
/* { */
/*     ulong ab_hi, ab_lo, q, r; */

/*     /1* ab_lo = (a * b) % 2**64  |  ab_hi = (a * b) >> nbits *1/ */
/*     n_mul2(&ab_hi, &ab_lo, a, b); */
/*     ab_hi = (ab_lo >> nbits) + (ab_hi << (64 - nbits)); */

/*     q = n_mulhi(ab_hi<<1, redp);  /1* (ab_hi * redp) >> 63 *1/ */
/*     r = ab_lo - q * n; */
/*     return r; */
/* } */


/** Modular multiplication, given precomputation */
/* assumes a<n, b<n (or close to this) */
/* nbits == number of bits in n */
ULONG_EXTRAS_INLINE ulong n_mod_mulmod_redp_lazy(ulong a, ulong b, ulong n, ulong redp, ulong nbits)
{
    ulong ab_hi, ab_lo, q, r;

    /* ab_lo = (a * b) % 2**64  |  ab_hi = (a * b) >> nbits */
    umul_ppmm(ab_hi, ab_lo, a, b);
    ab_hi = (ab_lo >> nbits) + (ab_hi << (64 - nbits));

    umul_ppmm(q, r, ab_hi<<1, redp);  /* (ab_hi * redp) >> 63 */
    r = ab_lo - q * n;
    return r;
}

/* assumes a<n, b<n (or close to this) */
ULONG_EXTRAS_INLINE ulong n_mod_mulmod_redp(ulong a, ulong b, ulong n, ulong redp, ulong nbits)
{
    ulong r = n_mod_mulmod_redp_lazy(a, b, n, redp, nbits);
    if (r >= 2*n)
        return r - 2*n;
    if (r >= n)
        return r - n;
    return r;
}

#endif
