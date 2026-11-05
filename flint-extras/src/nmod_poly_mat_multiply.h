/*
    Copyright (C) 2025 Vincent Neiger, Éric Schost

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#ifndef NMOD_POLY_MAT_MULTIPLY_H
#define NMOD_POLY_MAT_MULTIPLY_H

#include <flint/nmod_types.h>
#include <flint/nmod_poly.h>

#include "pml.h"

/* for Vandermonde1, we need 0, ..., len-1 to be distinct points in Z/modn Z */
#define NMOD_POLY_CAN_USE_VANDERMONDE1(modn, len) ((modn) >= (ulong)(len))
/* for Vandermonde2, we need 1**2, 2**2, ..., len**2 to be distinct points in
   Z/modn Z, and 2 to be invertible; for modn prime, this is modn odd and
   modn >= 2*len (the parity test only excludes modn == 2, len <= 1) */
#define NMOD_POLY_CAN_USE_VANDERMONDE2(modn, len) \
    (((modn) & UWORD(1)) && (modn) >= UWORD(2)*(len))
/* for Waksman, we need 2 to be invertible in Z/modn Z, i.e. modn odd */
#define NMOD_POLY_MAT_CAN_USE_WAKSMAN(modn) (((modn) & UWORD(1)) != 0)

/** Multiplication for polynomial matrices
 *  sets C = A * B
 *  output can alias input
 *  uses evaluation and interpolation at arithmetic points, done by matrix products
 *  v1 requires that the field contains distinct points 0, 1, 2, ..., len-1,
 *  where len == max_length(A) + max_length(B) - 1
 *     (for Z/pZ, equivalent to cardinality >= len, see NMOD_POLY_CAN_USE_VANDERMONDE1)
 *  v2 requires that the field contains distinct points 1**2, 2**2, ..., len**2
 *     and has characteristic != 2
 *     (for Z/pZ, equivalent to p odd and p >= 2 * len, see NMOD_POLY_CAN_USE_VANDERMONDE2)
 */
void nmod_poly_mat_mul_vandermonde1(nmod_poly_mat_t C, const nmod_poly_mat_t A, const nmod_poly_mat_t B);
void nmod_poly_mat_mul_vandermonde2(nmod_poly_mat_t C, const nmod_poly_mat_t A, const nmod_poly_mat_t B);

/** Multiplication for polynomial matrices
 *  sets C = A * B
 *  output can alias input
 *  uses Waksman's algorithm
 *  requires an odd modulus, see NMOD_POLY_MAT_CAN_USE_WAKSMAN
 */
void nmod_poly_mat_mul_waksman(nmod_poly_mat_t C, const nmod_poly_mat_t A,  const nmod_poly_mat_t B);

/** Multiplication for polynomial matrices, using evaluation-interpolation at geometric progression
 *  output can alias input
 *  in _precomp, G is precomputed to support evaluation and interpolation with G->len >= len1+len2-1
 *  the other variant computes G, assuming this is feasible (not checked)
 *  see NMOD_POLY_CAN_USE_GEOMETRIC
 */
void nmod_poly_mat_mul_geometric(nmod_poly_mat_t res, const nmod_poly_mat_t pmat1, const nmod_poly_mat_t pmat2);
void _nmod_poly_mat_mul_geometric_precomp(nmod_poly_mat_t res,
                                          const nmod_poly_mat_t pmat1, slong len1,
                                          const nmod_poly_mat_t pmat2, slong len2,
                                          nmod_geometric_progression_t G);


/** Multiplication for polynomial matrices, using FLINT's fft_small
 * evaluation-interpolation with the naive matrix multiplication algorithm
 * performed directly on the transforms
 *  sets C = A * B
 *  output can alias input
 *  Works for any modulus via a multi-modular strategy. The main computation
 *  is done modulo 1 to 4 50-bit FFT primes, the number of which depends on the
 *  modulus, on the inner dimension and on the lengths of the operands; this is
 *  followed by CRT. In terms of memory, it takes about (rdim*idim + idim*cdim
 *  + rdim*cdim) * nprimes * 2^depth doubles, where 2^depth is the transform
 *  length, up to a soft budget beyond which the rows of A are processed by
 *  groups; that budget is the larger of a fixed floor and twice the size of
 *  the operands and the result, and the storage is retained across calls.
 *  This falls back on nmod_poly_mat_mul if FLINT was built without fft_small
 */
void nmod_poly_mat_mul_sd_fft_direct(nmod_poly_mat_t C, const nmod_poly_mat_t A, const nmod_poly_mat_t B);

/** Multiplication for polynomial matrices, using FLINT's fft_small
 * evaluation-interpolation and calling nmod_mat_mul on the evaluations
 *  sets C = A * B
 *  output can alias input
 *  Same evaluation-interpolation scheme as nmod_poly_mat_mul_sd_fft_direct,
 *  and the same multi-modular strategy (1 to 4 50-bit FFT primes, then CRT),
 *  but the products in the transformed domain are performed by nmod_mat_mul
 *  on the matrix of evaluations at each point, rather than by a cubic kernel
 *  running along the transforms.
 *
 *  In terms of memory, this takes about (rdim*idim + idim*cdim + rdim*cdim) *
 *  nprimes * 2^depth words, where 2^depth is the transform length, up to a
 *  soft budget beyond which the rows of A and the columns of B are processed
 *  by groups; that budget is the larger of a fixed floor and twice the size
 *  of the operands and the result, and the storage is retained across calls.
 *  This falls back on nmod_poly_mat_mul if FLINT was built without fft_small
 */
void nmod_poly_mat_mul_sd_fft_matmul(nmod_poly_mat_t C, const nmod_poly_mat_t A, const nmod_poly_mat_t B);

/** general interface, picks an algorithm depending on parameters
 * TODO thresholds to be tuned
 * TODO add function to multiply with constant
 * TODO naming "multiply" because FLINT already has "mul"
 */
void nmod_poly_mat_multiply(nmod_poly_mat_t res, const nmod_poly_mat_t pmat1, const nmod_poly_mat_t pmat2);


/** Middle product for polynomial matrices
 * For all functions below:
 * They compute
 *       res  <-  (pmat1 * pmat2 mod x^nhi) div x^nlo
 * i.e., they set res to the first nhi - nlo middle coefficients [nlo, nhi) of
 * the product pmat1 * pmat2
 * Input len1 and len2 is such that pmat1 has length <= len1 and pmat2 has
 * length <= len2
 * Output can alias input
 *
 * TODO computation flow currently not ideal: reductions (e.g., if pmat1 has length << nlo)
 * are done in the general interface, and so will not be done in case one wants to directly
 * call the geometric variant
 *
 * TODO for evaluation-interpolation (currently, geometric), would be good to
 * have separate functions for matrix evaluation and matrix interpolation, and just call them?
 * (but then we also need it for "matrix evaluation of reverse"...)
 *
 * TODO faster matrix eval/interp by better handling the memory aspects

 * TODO improve the case of constant pmat1 or pmat2, currently no specific function
 */

/** general interface
 * truncates the operands at nhi, drops the useless low coefficients of one
 * of them (see the code for the rule), and calls _nmod_poly_mat_mulmid; the
 * reduced operands are views of the entries, no coefficient is copied (when
 * res aliases pmat1 or pmat2 and an operand is reduced, the result is
 * computed in a temporary) */
void nmod_poly_mat_mulmid(nmod_poly_mat_t res,
                          const nmod_poly_mat_t pmat1, const nmod_poly_mat_t pmat2,
                          slong nlo, slong nhi);

/** main interface
 * picks the fastest variant depending on parameters */
void _nmod_poly_mat_mulmid(nmod_poly_mat_t res,
                           const nmod_poly_mat_t pmat1, slong len1,
                           const nmod_poly_mat_t pmat2, slong len2,
                           slong nlo, slong nhi);

/** Middle products by a fixed right operand.
 *
 * (The suffix `_right` leaves room for the same with a fixed left operand.)
 *
 * Precomputation for repeated middle products
 *       res  <-  (pmat1 * pmat2 mod x^nhi') div x^nlo
 * with pmat2 fixed, pmat1 varying with entries of length at most `len1`, nlo
 * fixed, and nhi' at most the `nhi` given at initialization (e.g. the low
 * products by a truncated inverse in Dixon's lifting, where pmat1 is a row
 * vector). When the product goes through the fft_small transforms
 * (::_nmod_poly_mat_mulmid_sd_fft_direct), the plan and the transforms of the
 * entries of pmat2 are computed once and stored, and each product only
 * transforms pmat1, multiplies pointwise and transforms back the result.
 * Otherwise nothing is precomputed and each product calls
 * ::nmod_poly_mat_mulmid.
 *
 * - `pmat2` must not be modified or cleared while the precomputation is in
 *   use (it is read again by the products that do not use stored transforms)
 * - in ::nmod_poly_mat_mulmid_precomp_right, `res` may alias `pmat1` but not
 *   `pmat2`; if `pmat1` has an entry longer than `len1`, or if `nhi` is
 *   larger than at initialization (it is then reduced to that value), the
 *   result is still correct (::nmod_poly_mat_mulmid is called in the first
 *   case)
 * - memory: the stored transforms take about np * N doubles per entry of
 *   pmat2, where N is the transform length (about len1 + length(pmat2), or
 *   less for a window, see ::_nmod_poly_mat_mulmid_sd_fft_direct) and np the
 *   number of primes (1 to 4); see the code for the bound used by
 *   ::nmod_poly_mat_mulmid_precomp_right_init
 */
typedef struct
{
    const nmod_poly_mat_struct * pmat2;
    slong len1;     /* bound on the lengths of the entries of pmat1 */
    slong nlo;
    slong nhi;
    slong shift;    /* coefficients of pmat2 dropped below shift */
    void * fft;     /* stored transforms, or NULL */
} nmod_poly_mat_mulmid_precomp_right_struct;

typedef nmod_poly_mat_mulmid_precomp_right_struct nmod_poly_mat_mulmid_precomp_right_t[1];

/** Initializes the precomputation; stores transforms whenever they are
 * available (fft_small transforms, and a plan exists), unless they would take
 * more than the larger of 1 GiB and 16 times the coefficients of pmat2
 *
 * \todo for short operands and small dimensions, products of short
 * polynomials or the Vandermonde routines are faster than stored transforms
 * (see the rule of ::nmod_poly_mat_solve_series_dixon, which makes this
 * choice itself) */
void nmod_poly_mat_mulmid_precomp_right_init(nmod_poly_mat_mulmid_precomp_right_t pre,
                                             const nmod_poly_mat_t pmat2,
                                             slong len1, slong nlo, slong nhi);

/** Same as ::nmod_poly_mat_mulmid_precomp_right_init, with the choice given by
 * `method`: 0 for no precomputation, 1 for stored transforms when available
 * (without memory bound), -1 for the automatic choice */
void _nmod_poly_mat_mulmid_precomp_right_init(nmod_poly_mat_mulmid_precomp_right_t pre,
                                              const nmod_poly_mat_t pmat2,
                                              slong len1, slong nlo, slong nhi,
                                              int method);

/** Middle product `res = (pmat1 * pmat2 mod x^nhi) div x^nlo`, with pmat2
 * and nlo those of the precomputation */
void nmod_poly_mat_mulmid_precomp_right(nmod_poly_mat_t res,
                                        const nmod_poly_mat_t pmat1,
                                        const nmod_poly_mat_mulmid_precomp_right_t pre,
                                        slong nhi);

/** Whether the precomputation stores transforms */
static inline int
nmod_poly_mat_mulmid_precomp_right_is_fft(const nmod_poly_mat_mulmid_precomp_right_t pre)
{
    return pre->fft != NULL;
}

void nmod_poly_mat_mulmid_precomp_right_clear(nmod_poly_mat_mulmid_precomp_right_t pre);

/** actual worker: naive method, multiply + shift+truncate */
void _nmod_poly_mat_mulmid_naive(nmod_poly_mat_t res,
                                 const nmod_poly_mat_t pmat1, slong len1,
                                 const nmod_poly_mat_t pmat2, slong len2,
                                 slong nlo, slong nhi);

/** actual workers: middle product via the fft_small transforms, as the
 * multiplications nmod_poly_mat_mul_sd_fft_direct and
 * nmod_poly_mat_mul_sd_fft_matmul but with the plan built for the window
 * [nlo, nhi) of the convolution and only that window reconstructed
 * - no constraint on len1, len2 (entries of pmat1, pmat2 have length at
 *   most len1, len2); the transforms have the cyclic length N, a power of
 *   two, when max(len1, len2, nhi) <= N and len1 + len2 - 1 - N <= nlo,
 *   e.g. N = 2^ceil(log2(nhi)) when len1 <= nlo + 1 and len2 <= nhi, and
 *   the length of the full product otherwise
 * - output may alias input; multithreaded; memory bounded as in the
 *   multiplications, the result counted as its nhi - nlo coefficients
 */
void _nmod_poly_mat_mulmid_sd_fft_direct(nmod_poly_mat_t res,
                                         const nmod_poly_mat_t pmat1, slong len1,
                                         const nmod_poly_mat_t pmat2, slong len2,
                                         slong nlo, slong nhi);

void _nmod_poly_mat_mulmid_sd_fft_matmul(nmod_poly_mat_t res,
                                         const nmod_poly_mat_t pmat1, slong len1,
                                         const nmod_poly_mat_t pmat2, slong len2,
                                         slong nlo, slong nhi);

/** actual worker: transposed multiplication, via evaluation-interpolation at geometric progression
 * - requires len1 <= nlo+1 or len2 <= nlo+1
 * - requires the existence of a geometric progression of order (at least) nhi
 * - finds such a progression, builds temporary precomputation data,
 *   and calls `_nmod_poly_mat_mulmid_geometric_precomp`
 * - output may alias input; multithreaded (FLINT thread pool); the memory
 *   used for the constant matrices is bounded by the larger of a fixed floor
 *   and twice the size of the operands and the result, the rows of pmat1
 *   and the columns of pmat2 being processed by groups beyond it
 */
void _nmod_poly_mat_mulmid_geometric(nmod_poly_mat_t res,
                                     const nmod_poly_mat_t pmat1, slong len1,
                                     const nmod_poly_mat_t pmat2, slong len2,
                                     slong nlo, slong nhi);

/** same as above, with additional G, precomputed data for evaluation and
 * interpolation at a geometric progression of order at least nhi
 * requires len1 <= nlo+1 or len2 <= nlo+1
 */
void _nmod_poly_mat_mulmid_geometric_precomp(nmod_poly_mat_t res,
                                             const nmod_poly_mat_t pmat1, slong len1,
                                             const nmod_poly_mat_t pmat2, slong len2,
                                             slong nlo, slong nhi, nmod_geometric_progression_t G);


/** actual workers: transposed multiplication, via evaluation-interpolation
 * at small points with the evaluations and interpolations done by products
 * of constant matrices (the transposes of nmod_poly_mat_mul_vandermonde1
 * and nmod_poly_mat_mul_vandermonde2)
 * - no constraint on len1, len2, nlo, nhi (entries of pmat1, pmat2 have
 *   length at most len1, len2); the cost is that of the product of the
 *   shorter operand, of length ls = min(len1, len2), by a matrix of length
 *   nhi - nlo, i.e. P = ls + nhi - nlo - 1 points (after nhi is clamped to
 *   len1 + len2 - 1), against len1 + len2 - 1 for the naive approach
 * - v1 requires that 0, 1, ..., P-1 be distinct with invertible differences
 *   (prime modulus, NMOD_POLY_CAN_USE_VANDERMONDE1(modn, P))
 * - v2 requires that 1**2, 2**2, ..., P**2 be distinct with invertible
 *   differences, and 2 be invertible (odd prime modulus,
 *   NMOD_POLY_CAN_USE_VANDERMONDE2(modn, P))
 * - output may alias input
 */
void _nmod_poly_mat_mulmid_vandermonde1(nmod_poly_mat_t res,
                                        const nmod_poly_mat_t pmat1, slong len1,
                                        const nmod_poly_mat_t pmat2, slong len2,
                                        slong nlo, slong nhi);

void _nmod_poly_mat_mulmid_vandermonde2(nmod_poly_mat_t res,
                                        const nmod_poly_mat_t pmat1, slong len1,
                                        const nmod_poly_mat_t pmat2, slong len2,
                                        slong nlo, slong nhi);

#endif
