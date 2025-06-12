/*
    Copyright (C) 2025 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#ifndef __NMOD_POLY_FFT__H
#define __NMOD_POLY_FFT__H


/*-----------*/
/*-----------*/
/* TEMPORARY */
/*-----------*/
/*-----------*/


/***********************
*  bit reversed copy  *
***********************/

// #if defined(__GNUC__)
// # define FLINT_NO_VECTORIZE __attribute__((optimize("no-tree-vectorize")))
// #else
// # define FLINT_NO_VECTORIZE
// #endif

//inline long RevInc(long a, long k)
//{
//    long j, m;
//
//    j = k;
//    m = 1L << (k-1);
//
//    while (j && (m & a)) {
//        a ^= m;
//        m >>= 1;
//        j--;
//    }
//    if (j) a ^= m;
//    return a;
//}

// indices initialized with length >= k
//static inline void brc_indices(ulong * indices, long k)
//{
//    const long n = (1L << k);
//    for (long i = 0, j = 0; i < n; i++, j = RevInc(j, k))
//        indices[i] = j;
//}

//// counts in bit reversed depth, in C
//// or faster to build list as in dft.sage??
//void iter_reversed(ulong bits) {
//    ulong n = 1 << bits;
//
//    for (ulong i = 0, j = 0; i < n; i++) {
//        printf("%ld\n", j);
//
//        // Compute a mask of LSBs.
//        ulong mask = i ^ (i + 1);
//        // Length of the mask.
//        ulong len = __builtin_ctz(~mask);
//        // Align the mask to MSB of n.
//        mask <<= bits - len;
//        // XOR with mask.
//        j ^= mask;
//    }
//}

#endif
