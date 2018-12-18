#ifndef MAT_LZZ_P_EXTRA__H
#define MAT_LZZ_P_EXTRA__H

/** Additional functions for matrices over `zz_p`.
 *
 * \file mat_lzz_p_extra.h
 * \author Seung Gyu Hyun, Vincent Neiger, Eric Schost
 * \version 0.1
 * \date 2018-12-18
 *
 */

#include <NTL/mat_lzz_p.h>

NTL_CLIENT

/* [0 0 0 c] */
/* [1 0 0 0] */
/* [0 1 0 0] */
/* [0 0 1 0] */
/* in size n */
/** Builds and returns the `n x n` matrix with `1` on the subdiagonal and `c`
 * in the top-right position; precisely, the entries `(i+1,i)` are `1` for all
 * `i`, and the entry `(0,n-1)` is `c`, and the other entries are zero. */
Mat<zz_p> Z_lzz_p(long n, const zz_p & c);

/* [0 0 0 1] */
/* [0 0 1 0] */
/* [0 1 0 0] */
/* [1 0 0 0] */
/* in size n */
/** Builds and returns the `n x n` matrix with `1` on the antidiagonal;
 * precisely, the entries `(n-1-i,i)` are `1` for all `i`, and the other
 * entries are zero.  */
Mat<zz_p> J_lzz_p(long n);

/** Computes and returns the `n x n` diagonal matrix with entry at `(i,i)`
 * given by `d[i]` for all `i`, where the input vector is `d = (d[0], ...,
 * d[n-1])` */
Mat<zz_p> diagonal_matrix(const Vec<zz_p> & d);

#endif

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
