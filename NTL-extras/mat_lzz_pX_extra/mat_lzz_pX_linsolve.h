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



/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* x-adic algorithms for solving systems                      */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* solve A u = b mod x^prec                                   */
/* A square, A(0) invertible                                  */
/* use when deg(A) close to prec                              */
/* computes A^-1 mod x^{2^thresh}                             */
/* thresh=-1 is the default value, uses lookup table          */
/* output can alias input                                     */
/*------------------------------------------------------------*/
void solve_series_low_precision(
                                Mat<zz_pX> & u,
                                const Mat<zz_pX> & A,
                                const Mat<zz_pX> & b,
                                long prec,
                                long thresh = -1
                               );

inline Mat<zz_pX> solve_series_low_precision(
                                             const Mat<zz_pX> & A,
                                             const Mat<zz_pX> & b,
                                             long prec,
                                             long thresh = -1
                                            )
{ Mat<zz_pX> u; solve_series_low_precision(u, A, b, prec, thresh); return u; }

/*------------------------------------------------------------*/
/* solve A u = b mod x^prec                                   */
/* A square, A(0) invertible, deg(A), deg(b) < prec           */
/* use when deg(A) << prec                                    */
/*------------------------------------------------------------*/
void solve_series_high_precision(
                                 Mat<zz_pX> &u,
                                 const Mat<zz_pX>& A,
                                 const Mat<zz_pX>& b,
                                 long prec
                                );

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

/*------------------------------------------------------------*/
/* solve A u = b mod x^prec                                   */
/* A must be square, A(0) invertible                          */
/* output can alias input                                     */
/*------------------------------------------------------------*/
void solve_series(
                  Mat<zz_pX> & u, 
                  const Mat<zz_pX> & A,
                  const Mat<zz_pX> & b,
                  long prec
                 );

inline Mat<zz_pX> solve_series(
                               const Mat<zz_pX> & A,
                               const Mat<zz_pX> & b,
                               long prec
                              )
{ Mat<zz_pX> u; solve_series(u, A, b, prec); return u; }

/*------------------------------------------------------------*/
/* solve A u = b mod x^prec                                   */
/* A must be square, A(0) invertible                          */
/* output can alias input                                     */
/*------------------------------------------------------------*/
void solve_series(
                  Vec<zz_pX> & u,
                  const Mat<zz_pX> & A,
                  const Vec<zz_pX> & b,
                  long prec
                 );

inline Vec<zz_pX> solve_series(
                               const Mat<zz_pX> & A,
                               const Vec<zz_pX> & b,
                               long prec
                              )
{ Vec<zz_pX> u; solve_series(u, A, b, prec); return u; }


/*------------------------------------------------------------*/
/* Implements a minor variation of Storjohann's algorithm     */
/* A must be square, A(0) invertible, deg(b) < deg(A)         */
/* output can alias input                                     */
/*------------------------------------------------------------*/
void solve_series_high_order_lifting(
                                     Mat<zz_pX> & u,
                                     const Mat<zz_pX> & A,
                                     const Mat<zz_pX> & b,
                                     long prec
                                    );

inline Mat<zz_pX> solve_series_high_order_lifting(
                                                  const Mat<zz_pX> & A,
                                                  const Mat<zz_pX> & b,
                                                  long prec
                                                 )
{ Mat<zz_pX> u; solve_series_high_order_lifting(u, A, b, prec); return u; }


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
