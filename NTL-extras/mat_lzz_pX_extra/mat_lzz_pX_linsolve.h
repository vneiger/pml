#ifndef MAT_LZZ_PX_LINSOLVE__H
#define MAT_LZZ_PX_LINSOLVE__H

/** Functions for linear system solving.
 *
 * \file mat_lzz_pX_linsolve.h
 * \author Seung Gyu Hyun, Vincent Neiger, Eric Schost
 * \version 0.1
 * \date 2019-01-01
 *
 */

#include <NTL/matrix.h>
#include <NTL/lzz_pX.h>

NTL_CLIENT

/** @name Solve linear system over the power series
 * \anchor SeriesSolve
 *
 * These functions take as input a matrix `A`, a matrix or vector `b`, a
 * precision `prec`, and compute the power series linear system solution `u =
 * A^{-1} b mod X^prec`. The OUT parameter `u` can alias the IN parameters `A`
 * or `b`. The matrix `A` must be square, and most of these functions also
 * require that the constant coefficient of `A` is invertible.
 */
//@{

/** Computes `u = A^{-1} b mod X^prec` (see @ref SeriesSolve). The constant
 * coefficient `A(0)` must be invertible. To be used when the precision `prec`
 * is close to the degree of `A`. This starts by computing the truncated
 * inverse A^-1 mod x^thresh; thresh=-1 is the default value, uses lookup table
 * (see @ref TruncatedInverse). */
void solve_series_low_precision(
                                Mat<zz_pX> & u,
                                const Mat<zz_pX> & A,
                                const Mat<zz_pX> & b,
                                long prec,
                                long thresh = -1
                               );

/** Computes and returns `A^{-1} b mod X^prec` (see @ref SeriesSolve). The
 * constant coefficient `A(0)` must be invertible. To be used when the
 * precision `prec` is close to the degree of `A`. This starts by computing the
 * truncated inverse A^-1 mod x^thresh; thresh=-1 is the default value, uses
 * lookup table (see @ref TruncatedInverse). */
inline Mat<zz_pX> solve_series_low_precision(
                                             const Mat<zz_pX> & A,
                                             const Mat<zz_pX> & b,
                                             long prec,
                                             long thresh = -1
                                            )
{ Mat<zz_pX> u; solve_series_low_precision(u, A, b, prec, thresh); return u; }

/** Computes `u = A^{-1} b mod X^prec` (see @ref SeriesSolve). The constant
 * coefficient `A(0)` must be invertible. To be used when the precision `prec`
 * is significantly larger than the degree of `A`. */
void solve_series_high_precision(
                                 Mat<zz_pX> &u,
                                 const Mat<zz_pX>& A,
                                 const Mat<zz_pX>& b,
                                 long prec
                                );

/** Computes and returns `A^{-1} b mod X^prec` (see @ref SeriesSolve). The
 * constant coefficient `A(0)` must be invertible. To be used when the
 * precision `prec` is significantly larger than the degree of `A`. */
inline Mat<zz_pX> solve_series_high_precision(
                                              const Mat<zz_pX>& A,
                                              const Mat<zz_pX>& b,
                                              long prec
                                             )
{
    Mat<zz_pX> u;
    solve_series_high_precision(u, A, b, prec);
    return u;
}

/** Computes `u = A^{-1} b mod X^prec` (see @ref SeriesSolve). The constant
 * coefficient `A(0)` must be invertible, and `deg(b)` must be less than
 * `deg(A)`. To be used when the precision `prec` is significantly larger than
 * the degree of `A`. Implements a minor variation of Storjohann's algorithm
 *
 * Note: for the moment, this is slower than #solve_series_high_precision for
 * most input (if not all).
 */
void solve_series_high_order_lifting(
                                     Mat<zz_pX> & u,
                                     const Mat<zz_pX> & A,
                                     const Mat<zz_pX> & b,
                                     long prec
                                    );

/** Computes and returns `u = A^{-1} b mod X^prec` (see @ref SeriesSolve). The
 * constant coefficient `A(0)` must be invertible, and `deg(b)` must be less
 * than `deg(A)`. To be used when the precision `prec` is significantly larger
 * than the degree of `A`.
 *
 * Note: for the moment, this is slower than #solve_series_high_precision for
 * most input (if not all).
 */
inline Mat<zz_pX> solve_series_high_order_lifting(
                                                  const Mat<zz_pX> & A,
                                                  const Mat<zz_pX> & b,
                                                  long prec
                                                 )
{ Mat<zz_pX> u; solve_series_high_order_lifting(u, A, b, prec); return u; }



/** Computes `u = A^{-1} b mod X^prec` (see @ref SeriesSolve). The constant
 * coefficient `A(0)` must be invertible. Interface which tries to choose the
 * fastest available algorithm.
 *
 * \todo not properly tuned yet: relies on high precision variant when `prec`
 * exceeds `4 deg(A)`
 */
void solve_series(
                  Mat<zz_pX> & u, 
                  const Mat<zz_pX> & A,
                  const Mat<zz_pX> & b,
                  long prec
                 );

/** Computes and returns `A^{-1} b mod X^prec` (see @ref SeriesSolve). The
 * constant coefficient `A(0)` must be invertible. Interface which tries to
 * choose the fastest available algorithm.
 *
 * \todo not properly tuned yet: relies on high precision variant when `prec`
 * exceeds `4 deg(A)`
 */
inline Mat<zz_pX> solve_series(
                               const Mat<zz_pX> & A,
                               const Mat<zz_pX> & b,
                               long prec
                              )
{ Mat<zz_pX> u; solve_series(u, A, b, prec); return u; }

/** Computes `u = A^{-1} b mod X^prec` (see @ref SeriesSolve). The constant
 * coefficient `A(0)` must be invertible. Interface which tries to choose the
 * fastest available algorithm.
 *
 * \todo not properly tuned yet: relies on high precision variant when `prec`
 * exceeds `4 deg(A)`
 */
void solve_series(
                  Vec<zz_pX> & u,
                  const Mat<zz_pX> & A,
                  const Vec<zz_pX> & b,
                  long prec
                 );

/** Computes and returns `A^{-1} b mod X^prec` (see @ref SeriesSolve). The
 * constant coefficient `A(0)` must be invertible. Interface which tries to
 * choose the fastest available algorithm.
 *
 * \todo not properly tuned yet: relies on high precision variant when `prec`
 * exceeds `4 deg(A)`
 */
inline Vec<zz_pX> solve_series(
                               const Mat<zz_pX> & A,
                               const Vec<zz_pX> & b,
                               long prec
                              )
{ Vec<zz_pX> u; solve_series(u, A, b, prec); return u; }

//@} // doxygen group: Solve linear system over the power series

/*------------------------------------------------------------*/
/* solve A (u/den) = b                                        */
/* A must be square, A(0) invertible                          */
/* output can alias input                                     */
/* uses lifting and rational reconstruction                   */
/* nb_max is the max. number of vectors for use in            */
/* vector rational reconstruction (we use min(size, nb_max))  */
/* nb_max = -1 means a look-up table value is used            */
/*------------------------------------------------------------*/
long linsolve_via_series(
                         Vec<zz_pX> & u,
                         zz_pX & den,
                         const Mat<zz_pX> & A,
                         const Vec<zz_pX> & b,
                         long nb_max = -1
                        );

// solve aM = b via kernel basis
// return a and denominator d
// assumes M is invertible
// TODO not well tested yet
void linsolve_via_kernel(
                         Vec<zz_pX> & a,
                         zz_pX & d,
                         const Mat<zz_pX> & pmat,
                         const Vec<zz_pX> & b
                        );


#endif /* ifndef MAT_LZZ_PX_LINSOLVE__H */

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
