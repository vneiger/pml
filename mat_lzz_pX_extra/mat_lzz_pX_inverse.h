#ifndef MAT_LZZ_PX_INVERSE__H
#define MAT_LZZ_PX_INVERSE__H

/** \brief Inverse of a univariate polynomial matrix over `zz_p`
 *
 * \file mat_lzz_pX_inverse.h
 * \author Seung Gyu Hyun, Vincent Neiger, Eric Schost
 * \version 0.1
 * \date 2018-12-10
 *
 * Functions for computing the inverse of a polynomial matrix and related
 * operations such as the truncated inverse.
 *
 * \todo implementation of inversion algorithms:
 *   - evaluation-interpolation approach
 *   - Storjohann 2015 (fast Las Vegas)
 *   - Jeannerod-Villard J.Complexity 2005 (fast, correct for generic matrix)
 *   - Zhou-Labahn-Storjohann 2014 (fast deterministic, generalizes the above
 *   to any matrix, lower cost bound by compact representation, and also gives
 *   the largest Smith factor)
 *
 */

#include "mat_lzz_pX_multiply.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*                INVERSE TRUNCATED EXPANSION                 */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/** @name Truncated inverse expansion
 * \anchor TruncatedInverse
 * For a given square polynomial matrix `pmat` and a positive integer `d`, the
 * following functions compute the truncated expansion at order `d` of the
 * inverse `pmat^(-1)`, that is, the unique polynomial matrix `imat` of degree
 * less than `d` such that `imat*pmat = pmat*imat` is the identity matrix
 * modulo `x^d`. This expansion exists if and only if the constant coefficient
 * of `pmat` is invertible.
 * 
 * In all these functions, the OUT parameter `imat` may alias the IN parameter
 * `pmat`. The truncation order `d` must be a positive integer. These functions
 * throw an error if `pmat` is not square or if the constant coefficient
 * `pmat(0)` is not invertible.
 *  
 */
//@{

/** Computes the truncated inverse `imat` of `pmat` at order `d` (see @ref
 * TruncatedInverse) by Newton iteration, whose cost is quasi-linear in `d`.
 * This function is the main interface, which chooses the fastest between
 * the available implementations. Currently, these are:
 *   - if the modulus is a sufficiently large FFT prime, it relies on FFT
 *   evaluation-interpolation
 *   - otherwise, depending on how `d` compares to the relevant threshold
 *   value, it chooses between the middle product approach or the
 *   evaluation-interpolation at geometric points.
 */
void inv_trunc(Mat<zz_pX> & imat, const Mat<zz_pX> & pmat, long d);

/** Computes and returns the truncated inverse of `pmat` at order `d` (see @ref
 * TruncatedInverse) by Newton iteration, whose cost is quasi-linear in `d`.
 * This function is the main interface, which chooses the fastest between the
 * available implementations; see #inv_trunc for more details.
 */
inline Mat<zz_pX> inv_trunc(const Mat<zz_pX> & pmat, long d)
{ Mat<zz_pX> imat; inv_trunc(imat, pmat, d); return imat; }


/** Computes the truncated inverse `imat` of `pmat` at order `d` (see @ref
 * TruncatedInverse) by the naive algorithm, whose cost is quadratic in `d`
 */
void plain_inv_trunc(Mat<zz_pX> & imat, const Mat<zz_pX> & pmat, long d);

/** Computes the truncated inverse `imat` of `pmat` at order `d` (see @ref
 * TruncatedInverse) by Newton iteration, whose cost is quasi-linear in `d`.
 * This function is specifically for the case where the modulus is an FFT
 * prime; Newton iteration is implemented using FFT evaluation and
 * interpolation. This function calls the naive algo #plain_inv_trunc for
 * starting the iteration or if the order `d` is small; the crossover point is
 * for `d = thresh`; predetermined values are used if thresh = -1, which is the
 * default value.
 */
void newton_inv_trunc_FFT(
                          Mat<zz_pX> & imat,
                          const Mat<zz_pX> & pmat,
                          long d,
                          long thresh = -1
                         );

/** Computes the truncated inverse `imat` of `pmat` at order `d` (see @ref
 * TruncatedInverse) by Newton iteration, whose cost is quasi-linear in `d`,
 * implemented using the #middle_product function. This function calls the naive
 * algo #plain_inv_trunc for starting the iteration or if the order `d` is
 * small; the crossover point is for `d = thresh`; predetermined values are
 * used if thresh = -1, which is the default value.
 */
void newton_inv_trunc_middle_product(
                                     Mat<zz_pX> & imat,
                                     const Mat<zz_pX> & pmat,
                                     long d,
                                     long thresh = -1
                                    );

/** Computes the truncated inverse `imat` of `pmat` at order `d` (see @ref
 * TruncatedInverse) by Newton iteration, whose cost is quasi-linear in `d`.
 * The iteration is implemented using evaluation and interpolation at geometric
 * points. This function calls the naive algo #plain_inv_trunc for starting the
 * iteration or if the order `d` is small; the crossover point is for `d =
 * thresh`; predetermined values are used if thresh = -1, which is the
 * default value.
 */
void newton_inv_trunc_geometric(
                                Mat<zz_pX> & imat,
                                const Mat<zz_pX> & pmat,
                                long d,
                                long thresh = -1
                               );

//@} // doxygen group: Truncated inverse expansion


/** @name High-order lifting
 *
 *  Computes a high-degree slice of the inverse expansion of a matrix.
 */
//@{

/** Helper for high order lifting. As input, we have a polynomial matrix `pmat`
 * of degree `d`, and its truncated inverse `inv` at order `d`, that is, `inv =
 * pmat^{-1} mod x^d` (see #TruncatedInverse). In particular, the constant
 * coefficient `pmat(0)` must be invertible.
 *
 * For a given nonnegative integer `i`, we define `Si` to be the slice of the
 * expansion of `pmat^{-1}` of its terms of degree between `i-(2d-1)` and `i-1`
 * (both included).  Given `src = Si`, this function computes `next =
 * S_{2i-d}`. Note that deg(Si) < 2d-1.
 *
 * The OUT parameter `next` can alias the IN parameter `src`.
 */
void high_order_lift_inverse_odd(
                                 Mat<zz_pX> & next,
                                 const Mat<zz_pX> & src, 
                                 const std::unique_ptr<mat_lzz_pX_lmultiplier> & pmat, 
                                 const std::unique_ptr<mat_lzz_pX_lmultiplier> & inv,
                                 long d
                                );

//@} // doxygen group: High-order lifting

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*                         INVERSION                          */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/



#endif /* end of include guard: MAT_LZZ_PX_INVERSE__H */

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
