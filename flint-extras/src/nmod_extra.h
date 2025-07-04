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
ulong nmod_find_root(long n, nmod_t mod);


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

#endif
