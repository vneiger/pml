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

Vec<Vec<zz_pX>> mul(const Vec<Vec<zz_pX>> &S, const zz_pX &a);
Vec<Vec<zz_pX>> mul(const Vec<Vec<zz_pX>> &S, const zz_pXY &f);

// truncated multiplication
Vec<zz_pX> mulTrunc(const Vec<zz_pX> &S, const zz_pX &a, const long d);
Vec<zz_pX> mulTrunc(const Vec<zz_pX> &S, const zz_pXY &f, const long d);

Vec<Vec<zz_pX>> mulTrunc(const Vec<Vec<zz_pX>> &S, const zz_pX &a, 
		         const long d);
Vec<Vec<zz_pX>> mulTrunc(const Vec<Vec<zz_pX>> &S, const zz_pXY &f, 
		         const long d);

// add two sequences together
Vec<zz_pX> add(const Vec<zz_pX> &S, const Vec<zz_pX> &T);
Vec<Vec<zz_pX>> add(const Vec<Vec<zz_pX>> &S, const Vec<Vec<zz_pX>> &T);

// returns index of left-most non-zero term
long index_non_zero(const Vec<zz_pX> &S);

// checks if the sequence is all zeros
bool is_zero_seq(const Vec<zz_pX> &S);

/* single sequence */
// implementations of Kurakin's algorithm for computing the
// generators of ann(S) for some sequence S over k[x]/x^d
void kurakin(const long d, const Vec<zz_pX> &S, Vec<zz_pXY> &gens);

// the output will have only generators that are potentially useful;
// if the full list of generators are needed, use the function fill_in
void modified_kurakin(const long d, const Vec<zz_pX> &S, Vec<zz_pXY> &gens);

/* modules */
// implementation of Kurakin's algorithm for computing generators of
// ann(S) for sequence S over (k[x]/x^d)^tau
void kurakin(const long d, const Vec<Mat<zz_pX>> &S, Vec<zz_pXY> &gens);

// requires that gens is the output of modified_kurakin
void fill_in(const long d, Vec<zz_pXY> &gens);

/* for verifying */
// checks if each polynomial in gens cancels S
bool check_cancel(const Vec<zz_pX> &S, const Vec<zz_pXY> &gens, 
		  const long d=-1);

/* Wiedemann algorithms */

// Computes the minpoly using structured
// lifting from a (precomputed) sequence S, assuming S is non-degenerate
void minpoly_nondegenerate(const long d, const Vec<zz_pX> &S, zz_pXY &P);













#endif
