#ifndef __LZZ_PX_SEQ__H
#define __LZZ_PX_SEQ__H

#include <NTL/vector.h>
#include <NTL/lzz_pX.h>
#include "lzz_pXY.h"
NTL_CLIENT

/* utility functions */
// seq times univariate poly
Vec<zz_pX> mul(const Vec<zz_pX> &S, const zz_pX &a); 
// seq times bivariate poly
Vec<zz_pX> mul(const Vec<zz_pX> &S, const zz_pXY &f);
// truncated multiplication
Vec<zz_pX> mulTrunc(const Vec<zz_pX> &S, const zz_pX &a, const long d);

// add two sequences together
Vec<zz_pX> add(const Vec<zz_pX> &S, const Vec<zz_pX> &T); 

// returns index of left-most non-zero term
long index_non_zero(const Vec<zz_pX> &S);

// checks if the sequence is all zeros
bool is_zero_seq(const Vec<zz_pX> &S);

// implementations of Kurakin's algorithm for compute the
// generators of ann(S) for some sequence S over k[x]/x^d
void kurakin(const long d, const Vec<zz_pX> &S, Vec<zz_pXY> &gens);

#endif
