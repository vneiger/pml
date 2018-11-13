#ifndef MAT_LZZ_PX_INTERPOLANT__H
#define MAT_LZZ_PX_INTERPOLANT__H

#include <NTL/matrix.h>
#include <NTL/lzz_pX.h>
#include "mat_lzz_pX_forms.h" // for DegVec, Shift, PolMatForm

NTL_CLIENT

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
                          const Shift & shift,
                          const PolMatForm & form = ORD_WEAK_POPOV,
                          const bool randomized = false
                         );

// TODO (naive version written)
bool is_interpolant_basis(
                          const Mat<zz_pX> & intbas,
                          const Vec<Mat<zz_p>> & pmat, // vector of evaluations
                          const Vec<zz_p> & pts, // "uniform" case
                          const Shift & shift,
                          const PolMatForm & form = ORD_WEAK_POPOV,
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

#endif /* ifndef MAT_LZZ_PX_INTERPOLANT__H */

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
