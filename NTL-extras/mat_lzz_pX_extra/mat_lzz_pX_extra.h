#ifndef MAT_LZZ_PX_EXTRA__H
#define MAT_LZZ_PX_EXTRA__H

#include <NTL/matrix.h>
#include <NTL/lzz_pX.h>
#include <iostream>
#include <vector> // std vector, for shifts, degrees, pivot indices
#include <memory>
#include <algorithm>

#include "lzz_p_extra.h"
#include "lzz_pX_CRT.h"
#include "thresholds_matrix_multiply.h"
#include "thresholds_matrix_middle_product.h"
#include "thresholds_newton_inv_trunc.h"
#include "thresholds_solve_lift.h"

#include "mat_lzz_pX_utils.h"
#include "mat_lzz_pX_forms.h"
#include "mat_lzz_pX_arith.h"
#include "mat_lzz_pX_partial_linearization.h"
#include "mat_lzz_pX_multiply.h"
#include "mat_lzz_pX_approximant.h"
#include "mat_lzz_pX_kernel.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* UTILS                                                      */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* in-place reduction modulo the current prime                */
/*------------------------------------------------------------*/
void reduce_mod_p(Mat<zz_pX> & a);


/*------------------------------------------------------------*/
/* horizontal join                                            */
/* requires a.NumRows() == b.NumRows()                        */
/*------------------------------------------------------------*/
void horizontal_join(Mat<zz_pX>& c, const Mat<zz_pX>& a, const Mat<zz_pX>& b);

inline Mat<zz_pX> horizontal_join(const Mat<zz_pX>& a, const Mat<zz_pX>& b)
{
    Mat<zz_pX> c;
    horizontal_join(c, a, b);
    return c;
}

// TODO vertical join
// TODO vertical/horizonal splits (then update kernel basis)

/*------------------------------------------------------------*/
/* collapses s consecutive columns of a into one column of c  */
/* let t=a.NumCols(). For i=0..t/s-1, the i-th column of c is */
/* a[i*s] + x^d a[i*s+1] + ... + x^{(s-1)*d} a[i*s+s-1)]      */
/* requires that s divides t exactly                          */
/*------------------------------------------------------------*/
void collapse_consecutive_columns(Mat<zz_pX>& c, const Mat<zz_pX>& a, long d, long s);

inline Mat<zz_pX> collapse_consecutive_columns(const Mat<zz_pX>& a, long d, long s)
{
    Mat<zz_pX> c;
    collapse_consecutive_columns(c, a, d, s);
    return c;
}

/*------------------------------------------------------------*/
/* collapses columns with stepsize s of a into a column of c  */
/* let t=a.NumCols(). For i=0..s-1, the i-th column of c is   */
/* a[i] + x^d a[i+s] + ... + x^{(t/s-1)*d} a[i+(t/s-1)*s)]    */
/* requires that s divides t exactly                          */
/*------------------------------------------------------------*/
void collapse_nonconsecutive_columns(Mat<zz_pX>& c, const Mat<zz_pX>& a, long d, long s);

inline Mat<zz_pX> collapse_nonconsecutive_columns(const Mat<zz_pX>& a, long d, long s)
{
    Mat<zz_pX> c;
    collapse_nonconsecutive_columns(c, a, d, s);
    return c;
}



/*------------------------------------------------------------*/
/* multipoint evaluation for matrices                         */
/*------------------------------------------------------------*/
void matrix_evaluate (Vec<Mat<zz_p>> &evals,
                      const Mat<zz_pX> &pmat,
                      const zz_pX_Multipoint &ev);


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* Division with remainder                                    */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

// computes Q,R such that A = BQ + R
void quo_rem(Mat<zz_pX> &Q, 
             Mat<zz_pX> &R, 
             const Mat<zz_pX> &A,
             const Mat<zz_pX> &B);




/**********************************************************************
 *                     MINIMAL INTERPOLANT BASES                      *
 **********************************************************************/

// list of points: for each column, we have a list of points (elt from zz_p and
// multiplicity)
// FIXME not thought thorougly yet, subject to change
typedef std::vector<std::vector<std::pair<zz_p,long>>> Points;

////Definition (interpolant basis)
// Given:
//   * m x n matrix of univariate polynomials 'pmat',
//   * list 'points' of n pairs [root,multiplicity], which define n products of linear factors M_0,...,M_{n-1},
// An interpolant basis for (pmat,points) is a matrix over K[X]
// whose rows form a basis for the K[X]-module
// { 'int' in K[X]^{1 x m}  |  the column j of 'app' 'pmat' is 0 modulo M_j }

/*------------------------------------------------------------*/
/* general user-friendly interface                            */
/*------------------------------------------------------------*/

// TODO
DegVec interpolant_basis(
                         Mat<zz_pX> & intbas,
                         const Mat<zz_pX> & pmat,
                         const Points & pts,
                         const Shift & shift = Shift(),
                         const PolMatForm form = ORD_WEAK_POPOV,
                         const bool row_wise = true,
                         const bool generic = false
                        );

// TODO uncomment below when above is ready
// below, the following is called "uniform interpolation case" (or case with uniform points)
// this means that we have the same points on all columns, all with multiplicity one)
// (FIXME could be easily generalized to any constant multiplicity for all...?)
DegVec interpolant_basis(
                         Mat<zz_pX> & intbas,
                         const Mat<zz_pX> & pmat,
                         const Vec<zz_p> & pts,
                         const Shift & shift = Shift(),
                         const PolMatForm form = ORD_WEAK_POPOV,
                         const bool row_wise = true,
                         const bool generic = false
                        );
//{
//  std::vector<std::pair<zz_p,long>> list_pts(pts.size());
//  for ( long i=0; i<list_pts.size(); ++i )
//  {
//    list_pts[i] = std::pair<zz_p,long>(pts[i],1);
//  }
//  Points points(mat.NumCols(),list_pts);
//  return interpolant_basis(appbas,mat,points,shift,canonical,row_wise,generic);
//}


/*------------------------------------------------------------*/
/* Verifying that intbas is a shift-minimal interpolant       */
/* basis for input matrix 'pmat' and points 'points'          */
/* 'form' gives the minimal requirement to check (matrix must */
/* be at least in the form 'form')                            */
/* 'randomized' says whether using a Monte Carlo or Las Vegas */
/* verification algorithm is acceptable                       */
/* Note: currently, deterministic verification is, for most   */
/* instances, as long as re-computing the basis               */
/*------------------------------------------------------------*/

// TODO not implemented yet
bool is_interpolant_basis(
                          const Mat<zz_pX> & intbas,
                          const Mat<zz_pX> & pmat,
                          const Points & pts,
                          const Shift & shift = Shift(),
                          const PolMatForm & form = ORD_WEAK_POPOV,
                          const bool row_wise = true,
                          const bool randomized = false
                         );

// TODO (naive version written)
bool is_interpolant_basis(
                          const Mat<zz_pX> & intbas,
                          const Vec<Mat<zz_p>> & pmat, // vector of evaluations
                          const Vec<zz_p> & pts, // "uniform" case
                          const Shift & shift = Shift(),
                          const PolMatForm & form = ORD_WEAK_POPOV,
                          const bool row_wise = true,
                          const bool randomized = false
                         );

/*------------------------------------------------------------*/
/* Iterative algorithm for arbitrary points and shift         */
/* References:                                                */
/*   - Beckermann 1992                                        */
/*   - Van Barel-Bultheel 1991+1992                           */
/*   - Beckermann-Labahn 2000 (ensuring s-Popov)              */
/*------------------------------------------------------------*/
DegVec intbas_iterative(
                        Mat<zz_pX> & intbas,
                        const Mat<zz_pX> & pmat,
                        const Points & pts,
                        const Shift & shift,
                        bool point_wise=true // TODO to be thought about
                       );

DegVec popov_intbas_iterative(
                              Mat<zz_pX> & intbas,
                              const Mat<zz_pX> & pmat,
                              const Points & pts,
                              const Shift & shift,
                              bool point_wise=true // TODO to be thought about
                             );


/*------------------------------------------------------------*/
/* Adaptation of M-Basis for uniform interpolation points     */
/*------------------------------------------------------------*/

// --> popov_mbasis1 can be called as such (with, as input, pmat evaluated at a
// point)

// TODO input pmat = polynomial matrix, not implemented yet
DegVec mbasis(
              Mat<zz_pX> & intbas,
              const Mat<zz_pX> & pmat,
              const Vec<zz_p> & pts,
              const Shift & shift
             );

// input pmat = list of evaluations, implemented
DegVec mbasis(
              Mat<zz_pX> & intbas,
              const Vec<Mat<zz_p>> & evals,
              const Vec<zz_p> & pts,
              const Shift & shift
             );

DegVec popov_mbasis(
                    Mat<zz_pX> &intbas,
                    const Mat<zz_pX> & pmat,
                    const Vec<zz_p> & pts,
                    const Shift & shift
                   );

/*------------------------------------------------------------*/
/* PM-Basis algorithm for uniform interpolation points        */
/*------------------------------------------------------------*/

// TODO there are two variants, test them to be sure if they are similar / which is faster
//   either compute more in the evaluated world and interpolate intbas at the end,
//   or compute in the polynomial world and evaluate to obtain the residuals
// (in any case, there will still be interpolation/evaluation in the middle)
// TODO input pmat = polynomial matrix, not implemented yet
DegVec pmbasis(
               Mat<zz_pX> & intbas,
               const Mat<zz_pX> & pmat,
               const Vec<zz_p> & pts,
               const Shift & shift
              );

// returns the points and matrix evaluations used       
DegVec pmbasis_geometric(
               Mat<zz_pX> & intbas,
               const Mat<zz_pX> & pmat,
               const zz_p & r,
               const long order,
               const Shift & shift,
               Vec<Mat<zz_p>> &evals,
               Vec<zz_p> &pts
              );

// requires that pts contain powers of r
// with entries of evals evaluated at pts
DegVec pmbasis_geometric(
                         Mat<zz_pX> & intbas,
                         const Vec<Mat<zz_p>> & evals,
                         const Vec<zz_p> & pts,
                         const zz_p & r,
                         const Shift & shift
                        );

// input pmat = list of evaluations, implemented
DegVec pmbasis(
              Mat<zz_pX> & intbas,
              const Vec<Mat<zz_p>> & evals,
              const Vec<zz_p> & pts,
              const Shift & shift
             );

DegVec popov_pmbasis(
                     Mat<zz_pX> &intbas,
                     const Mat<zz_pX> & pmat,
                     const Vec<zz_p> & pts,
                     const Shift & shift
                    );


/**********************************************************************
 *                             ROW BASIS                              *
 **********************************************************************/

/*************************************
 *  Li and Storjohann's compression  *
 *************************************/

// TODO, cf Chao Li's 2007 Master's thesis
// does not really give a row basis, but allows to reduce to almost
// square case by an efficient Las Vegas algorithm

/***************************
 *  Zhou-Labahn algorithm  *
 ***************************/

// TODO

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*                      INVERSE EXPANSION                     */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* returns x = 1/a mod z^m, quadratic algorithm               */
/* throws an error if a(0) not invertible                     */
/* x can alias a                                              */
/*------------------------------------------------------------*/
void plain_inv_trunc(Mat<zz_pX>& x, const Mat<zz_pX>& a, long m);

/*------------------------------------------------------------*/
/* returns x = 1/a mod z^m, Newton iteration                  */
/* throws an error if a(0) not invertible                     */
/* x can alias a                                              */
/* crossover point with naive algo is m = 2^thresh            */
/* thresh = -1 means we use predetermined values              */
/*------------------------------------------------------------*/
void newton_inv_trunc_FFT(Mat<zz_pX>& x, const Mat<zz_pX>& a, long m, long thresh = -1);
void newton_inv_trunc_middle_product(Mat<zz_pX>& x, const Mat<zz_pX>& a, long m, long thresh = -1);
void newton_inv_trunc_geometric(Mat<zz_pX>& x, const Mat<zz_pX>& a, long m, long thresh = -1);

/*------------------------------------------------------------*/
/* returns x = 1/a mod z^m                                    */
/* throws an error if a(0) not invertible                     */
/* x can alias a                                              */
/*------------------------------------------------------------*/
void inv_trunc(Mat<zz_pX>& x, const Mat<zz_pX>& a, long m);

inline Mat<zz_pX> inv_trunc(const Mat<zz_pX>& a, long m)
{
    Mat<zz_pX> y;
    inv_trunc(y, a, m);
    return y;
}

/*------------------------------------------------------------*/
/* for i >= 0, define Si = coefficients of A^{-1} of degrees  */
/*             i-(2d-1) .. i-1, with d=deg(A)                 */
/* given src = Si, this computes S_{2i-d}                     */
/* invA = A^{-1} mod x^d                                      */
/* note: deg(Si) < 2d-1                                       */
/* output can alias input                                     */
/*------------------------------------------------------------*/
void high_order_lift_inverse_odd(Mat<zz_pX> & next, const Mat<zz_pX>& src, 
                                 std::unique_ptr<mat_lzz_pX_lmultiplier> & A, 
                                 std::unique_ptr<mat_lzz_pX_lmultiplier> & invA, long d);

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
void solve_series_low_precision(Mat<zz_pX> &u, const Mat<zz_pX>& A, const Mat<zz_pX>& b, long prec, long thresh = -1);

inline Mat<zz_pX> solve_series_low_precision(const Mat<zz_pX>& A, const Mat<zz_pX>& b, long prec, long thresh = -1)
{
    Mat<zz_pX> u;
    solve_series_low_precision(u, A, b, prec, thresh);
    return u;
}

/*------------------------------------------------------------*/
/* solve A u = b mod x^prec                                   */
/* A square, A(0) invertible, deg(A), deg(b) < prec           */
/* use when deg(A) << prec                                    */
/*------------------------------------------------------------*/
void solve_series_high_precision(Mat<zz_pX> &u, const Mat<zz_pX>& A, const Mat<zz_pX>& b, long prec);

inline Mat<zz_pX> solve_series_high_precision(const Mat<zz_pX>& A, const Mat<zz_pX>& b, long prec)
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
void solve_series(Mat<zz_pX> &u, const Mat<zz_pX>& A, const Mat<zz_pX>& b, long prec);

inline Mat<zz_pX> solve_series(const Mat<zz_pX>& A, const Mat<zz_pX>& b, long prec)
{
    Mat<zz_pX> u;
    solve_series(u, A, b, prec);
    return u;
}

/*------------------------------------------------------------*/
/* solve A u = b mod x^prec                                   */
/* A must be square, A(0) invertible                          */
/* output can alias input                                     */
/*------------------------------------------------------------*/
void solve_series(Vec<zz_pX> &u, const Mat<zz_pX>& A, const Vec<zz_pX>& b, long prec);

inline Vec<zz_pX> solve_series(const Mat<zz_pX>& A, const Vec<zz_pX>& b, long prec)
{
    Vec<zz_pX> u;
    solve_series(u, A, b, prec);
    return u;
}


/*------------------------------------------------------------*/
/* Implements a minor variation of Storjohann's algorithm     */
/* A must be square, A(0) invertible, deg(b) < deg(A)         */
/* output can alias input                                     */
/*------------------------------------------------------------*/
void solve_series_high_order_lifting(Mat<zz_pX> &u, const Mat<zz_pX>& A, const Mat<zz_pX>& b, long prec);

inline Mat<zz_pX> solve_series_high_order_lifting(const Mat<zz_pX>& A, const Mat<zz_pX>& b, long prec)
{
    Mat<zz_pX> u;
    solve_series_high_order_lifting(u, A, b, prec);
    return u;
}


/*------------------------------------------------------------*/
/* solve A (u/den) = b                                        */
/* A must be square, A(0) invertible                          */
/* output can alias input                                     */
/* uses lifting and rational reconstruction                   */
/* nb_max is the max. number of vectors for use in            */
/* vector rational reconstruction (we use min(size, nb_max))  */
/* nb_max = -1 means a look-up table value is used            */
/*------------------------------------------------------------*/
long linsolve_via_series(Vec<zz_pX> &u, zz_pX& den, const Mat<zz_pX>& A, const Vec<zz_pX>& b, long nb_max = -1);


// TODO: polynomial matrix division with remainder (cf. e.g. Neiger-Vu 2017)
// TODO see Storjohann 2003 for high-order lifting
// TODO application to linear system solving
// --> see which other consequences of high-order lifting may be worth implementing


/**********************************************************************
 *                             INVERSION                              *
 **********************************************************************/

// TODO evaluation-interpolation approach

// TODO Jeannerod-Villard: worth implementing?

// TODO Storjohann 2015 fast Las Vegas

// TODO Zhou-Labahn-Storjohann 2014 fast deterministic

/**********************************************************************
 *                          BASIS REDUCTION                           *
 *            (shifted reduced form and shifted normal forms)         *
 **********************************************************************/

// TODO general reduction to uniform shift via pre-multiplication
// worthwile at least when shift close to uniform

// TODO naive algorithms (see Mulders-Storjohann for good reference)

// TODO general shifted Popov form via kernel (itself via approximant basis)

// TODO understand if there is any chance Alekhnovich improves over the
// kernel approach

// TODO nonsingular: Giorgi-Jeannerod-Villard's Las Vegas reduction
// (worth implementing for shifts other than uniform?)


/**********************************************************************
 *                    TRIANGULARIZATION ALGORITHMS                    *
 **********************************************************************/

/**************************************************
 *  Labahn-Neiger-Zhou partial triangularization  *
 **************************************************/

// TODO one step (uses kernel + row basis)

// TODO diagonal entries of Hermite form, not implemented yet
void diagonal_of_hermite(Vec<zz_pX> & diag, const Mat<zz_pX> & pmat);

/***********************************************
 *  Labahn-Neiger-Zhou Hermite form algorithm  *
 ***********************************************/

// TODO (requires partial linearization + basis reduction)


/**********************************************************************
 *                       DETERMINANT ALGORITHMS                       *
 **********************************************************************/

// general user interface
// TODO (not implemented yet)
void determinant(zz_pX & det, const Mat<zz_pX> & pmat);

inline zz_pX determinant(const Mat<zz_pX> & pmat)
{
    zz_pX det;
    determinant(det, pmat);
    return det;
}

// verifies that det = c det(pmat),
// for some c a nonzero field element if up_to_constant==false; and c=1 otherwise
// if randomized==true, it is allowed to use a Monte Carlo randomized approach
// TODO: only randomized implemented for now
bool verify_determinant(const zz_pX & det, const Mat<zz_pX> & pmat, bool up_to_constant, bool randomized);

/*******************************************************************
 *  Labahn-Neiger-Zhou: via diagonal entries of triangularization  *
 *******************************************************************/

void determinant_via_diagonal_of_hermite(zz_pX & det, const Mat<zz_pX> & pmat);

// TODO version assuming diagonal of (lower triangular) Hermite form is
// [det 1 .. 1] and that we have generic degree profiles in the partial
// triangularization of Labahn-Neiger-Zhou

// Version 1 (deterministic algorithm): degree of determinant is known (e.g. if
// matrix is reduced or if the determinant corresponds to some invariant of an
// object with known degree). Runs the partial triangularization and returns
// false if the computed diagonal entry of Hermite form does not have the right
// degree. True is returned iff determinant is correct. 
// TODO currently returns determinant up to constant factor!!
bool determinant_generic_knowing_degree(zz_pX & det, const Mat<zz_pX> & pmat, long degree);

// Version 2 (Las Vegas randomized algorithm): runs the partial
// triangularization and checks determinant by Zippel-Schwartz; returns false
// if determinant is wrong or if field size is too small for allowing the
// Zippel-Schwartz check. If true is returned, then determinant is correct.
// TODO: not implemented yet
bool determinant_generic_las_vegas(zz_pX & det, const Mat<zz_pX> & pmat);

// Version 3 (randomized; via random linear system solving)
// TODO first version, should be improved. Make Las Vegas.
void determinant_via_linsolve(zz_pX & det, const Mat<zz_pX> & pmat);

// TODO other determinant algorithms??
// --> could rely on x-Smith decomposition of Gupta et al (worth
// implementing??), cf Appendix of LaNeZh17 

#endif // MAT_LZZ_PX_EXTRA__H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
