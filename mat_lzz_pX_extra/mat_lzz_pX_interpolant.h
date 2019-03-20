#ifndef MAT_LZZ_PX_INTERPOLANT__H
#define MAT_LZZ_PX_INTERPOLANT__H

/** \brief Minimal interpolant bases.
 *
 * \file mat_lzz_pX_interpolant.h
 * \author Seung Gyu Hyun, Vincent Neiger, Eric Schost
 * \version 0.1
 * \date 2018-12-22
 *
 */

#include "mat_lzz_pX_forms.h" // for VecLong, VecLong, PolMatForm

NTL_CLIENT

/** \file mat_lzz_pX_interpolant.h
 * Definition (interpolant basis).
 * -------------------------------
 * Consider
 *   - an m x n matrix of univariate polynomials F,
 *   - a list of d elements pts = (x_0,...,x_{d-1}) from the base field.
 *
 * Then an interpolant basis for (F,pts) is a matrix over the univariate
 * polynomials whose rows form a basis for the following module:
 * { p in K[X]^{1 x m}  |  p(x_i) F(x_i) = 0 for 0 <= i < d }.
 * Note that such a matrix is square, m x m, and nonsingular.
 */

/** \file mat_lzz_pX_interpolant.h
 * Definition (shifted minimal approximant basis).
 * -----------------------------------------------
 * Considering furthermore:
 *   - a degree shift s (a list of m integers).
 *
 * Then an interpolant basis for (F,d) is said to be <em>a shift-minimal</em>
 * (resp. <em>a shift-ordered weak Popov</em>, resp. <em>the shift-Popov</em>)
 * interpolant basis if it is in shift-reduced form (resp. in shift-ordered
 * weak Popov form, resp. in shift-Popov form). See mat_lzz_pX_forms.h
 * for definitions of these forms.
 */

/** \file mat_lzz_pX_interpolant.h
 * Conventions.
 * ------------
 * Apart from the general interfaces which offer the choice between left or
 * right interpolants, all other functions compute left interpolant bases
 * (interpolants operate on the left of the matrix F; the basis elements are
 * the rows of the matrix).
 *
 * Most functions below use the following parameters.
 *
 * \param[out] intbas the output interpolant basis (cannot alias `pmat`)
 * \param[in] pmat the input polynomial matrix (no restriction)
 * \param[in] pts the input points (list of field elements)
 * \param[in,out] shift the input shift and output shifted row degree of
 * `intbas` (list of integers, length must be `pmat.NumRows()`)
 *
 * Note that the latter restriction on the length of the list is assuming left
 * interpolants (for right interpolants, it would be `pmat.NumCols()`).
 *
 */

/** \file mat_lzz_pX_interpolant.h
 * Not implemented yet -- possible future work.
 * --------------------------------------------
 * Currently, only interpolants with matrix-wise points and no multiplicity are
 * considered; future work may extend functionalities to more general cases,
 * namely:
 *  - column-wise points: instead of requiring that `p(x_i) F(x_i) = 0`, one
 *  may more generally ask that each column of `p F` vanishes at a given point
 *  (which can differ from column to column)
 *  - multiplicities: one may further ask that each column vanishes with a given
 *  multiplicity.
 *
 *  Put together, these two generalizations would allow one to consider
 *  interpolant equations where each column of `p F` vanishes modulo some
 *  polynomial which splits, with known roots and multiplicities.
 *
 */


/** Verifies that `intbas` is a `shift`-minimal interpolant basis for `(pmat,pts)`
 * with the given form `form`.
 *
 * \param[in] intbas approximant basis
 * \param[in] pmat polynomial matrix
 * \param[in] pts interpolation points
 * \param[in] shift shift
 * \param[in] form required form for `intbas` (see #PolMatForm)
 * \param[in] row_wise indicates whether we consider left interpolants (working row-wise) or right interpolants (working column-wise)
 * \param[in] randomized if `true`, the algorithm may use a Monte Carlo or Las Vegas verification algorithm
 *
 * \return boolean, result of the verification
 *
 * \todo add parameter row_wise
 * \todo support all options, make doc more clear concerning Las Vegas / Monte Carlo
 *
 **/
bool is_interpolant_basis(
                          const Mat<zz_pX> & intbas,
                          const Mat<zz_pX> & pmat,
                          const Vec<zz_p> & pts, // "uniform" case
                          const VecLong & shift,
                          const PolMatForm & form = ORD_WEAK_POPOV,
                          const bool randomized = false
                         );


/** Verifies that `intbas` is a `shift`-minimal interpolant basis for `(pmat,pts)`
 * with the given form `form`.
 *
 * \param[in] intbas approximant basis
 * \param[in] evals evaluation matrices
 * \param[in] pts interpolation points
 * \param[in] shift shift
 * \param[in] form required form for `intbas` (see #PolMatForm)
 * \param[in] row_wise indicates whether we consider left interpolants (working row-wise) or right interpolants (working column-wise)
 * \param[in] randomized if `true`, the algorithm may use a Monte Carlo or Las Vegas verification algorithm
 *
 * \return boolean, result of the verification
 *
 * \todo add parameter row_wise
 * \todo support all options, make doc more clear concerning Las Vegas / Monte Carlo
 *
 **/
bool is_interpolant_basis(
                          const Mat<zz_pX> & intbas,
                          const Vec<Mat<zz_p>> & evals, // vector of evaluations
                          const Vec<zz_p> & pts, // "uniform" case
                          const VecLong & shift,
                          const PolMatForm & form = ORD_WEAK_POPOV,
                          const bool randomized = false
                         );

/** Verifies that `intbas` is a `shift`-minimal interpolant basis for `(pmat,pts)`
 * with the given form `form`, where `pts` is the geometric sequence defined by
 * `pt` and `order`.
 *
 * \param[in] intbas approximant basis
 * \param[in] pmat polynomial matrix
 * \param[in] pt interpolation initial point
 * \param[in] order length of sequence
 * \param[in] shift shift
 * \param[in] form required form for `intbas` (see #PolMatForm)
 * \param[in] row_wise indicates whether we consider left interpolants (working row-wise) or right interpolants (working column-wise)
 * \param[in] randomized if `true`, the algorithm may use a Monte Carlo or Las Vegas verification algorithm
 *
 * \return boolean, result of the verification
 *
 * \todo add parameter row_wise
 * \todo support all options, make doc more clear concerning Las Vegas / Monte Carlo
 *
 **/
bool is_interpolant_basis_geometric(
                                    const Mat<zz_pX> & intbas,
                                    const Mat<zz_pX> & pmat,
                                    const zz_p & pt, // geometric case
                                    const long order,
                                    const VecLong & shift,
                                    const PolMatForm & form = ORD_WEAK_POPOV,
                                    const bool randomized = false
                                   );

/** Verifies that `intbas` is a `shift`-minimal interpolant basis for `(pmat,pts)`
 * with the given form `form`, where `pts` is the geometric sequence defined by
 * `pt` and `order`.
 *
 * \param[in] intbas approximant basis
 * \param[in] evals evaluation matrices
 * \param[in] pt interpolation initial point
 * \param[in] order length of sequence
 * \param[in] shift shift
 * \param[in] form required form for `intbas` (see #PolMatForm)
 * \param[in] row_wise indicates whether we consider left interpolants (working row-wise) or right interpolants (working column-wise)
 * \param[in] randomized if `true`, the algorithm may use a Monte Carlo or Las Vegas verification algorithm
 *
 * \return boolean, result of the verification
 *
 * \todo add parameter row_wise
 * \todo support all options, make doc more clear concerning Las Vegas / Monte Carlo
 *
 **/
bool is_interpolant_basis_geometric(
                                    const Mat<zz_pX> & intbas,
                                    const Vec<Mat<zz_p>> & evals,
                                    const zz_p & pt, // geometric case
                                    const long order,
                                    const VecLong & shift,
                                    const PolMatForm & form = ORD_WEAK_POPOV,
                                    const bool randomized = false
                                   );



/** @name M-Basis algorithm
 * \anchor MBasisInt
 *
 * These functions compute a `shift`-minimal ordered weak Popov approximant
 * basis for `(pmat,pts)`. They use an iteration on the points, computing at
 * each step a basis for a single point (using #mbasis1 with input matrix an
 * evaluation of `pmat` at this point), and using it to update the output
 * `intbas`, the so-called _residual matrix_, and the considered shift. After
 * `d` iterations, `intbas*pmat` is zero at the first `d` points.
 *
 * At the end of the computation, the vector `shift` contains the shifted row
 * degree of `intbas`, for the input shift. 
 *
 * In this context, the residual matrix is a constant matrix with the same
 * dimensions as `pmat` which, at the iteration `d`, is equal to the evaluation
 * of `intbas*pmat` at the `d`-th point.
 *
 */
//@{

/** Computes a `shift`-ordered weak Popov interpolant basis `intbas` for
 * `(pmat,pts)`, with the matrix `pmat` given by its list of evaluations
 * `evals` at all the points in `pts` (see @ref MBasisInt). Variant where the
 * residual is computed at each iteration, by evaluating `intbas` and
 * multiplying by the corresponding evaluation in `evals`. The positive
 * integers `offset` and `order` indicate that we consider the entries `offset,
 * offset+1, .., offset+order-1` of `evals` and `pts` (no check is performed to
 * verify that these indices stay within the allowed bounds).
 *
 * Requirements: this function assumes that there are no repeated points in
 * `pts` (undefined behaviour otherwise), that `evals` and `pts` have the
 * same strictly positive length, and that all matrices in `evals` have the
 * same dimensions.
 * 
 **/
void mbasis_rescomp(
                    Mat<zz_pX> & intbas,
                    const Vec<Mat<zz_p>> & evals,
                    const Vec<zz_p> & pts,
                    VecLong & shift,
                    long offset,
                    long order
                   );

/** Computes a `shift`-ordered weak Popov interpolant basis `intbas` for
 * `(pmat,pts)`, with the matrix `pmat` given by its list of evaluations
 * `evals` at all the points in `pts` (see @ref MBasisInt). Variant where the
 * evaluations in `evals` are continuously updated along the iterations. The
 * positive integers `offset` and `order` indicate that we consider the entries
 * `offset, offset+1, .., offset+order-1` of `evals` and `pts` (no check is
 * performed to verify that these indices stay within the allowed bounds).
 *
 * Requirements: this function assumes that there are no repeated points in
 * `pts` (undefined behaviour otherwise), and that `evals` and `pts` have the
 * same strictly positive length, and that all matrices in `evals` have the
 * same dimensions.
 * 
 **/
void mbasis_resupdate(
                      Mat<zz_pX> & intbas,
                      Vec<Mat<zz_p>> & evals,
                      const Vec<zz_p> & pts,
                      VecLong & shift,
                      long offset,
                      long order
                     );


// input pmat = list of evaluations, implemented
// REQUIREMENT : len(evals) == len(pts) > 0
// assumes no repeated points (will not fail but undefined behaviour)
// (one could e.g. do a cleaning of pts beforehand)
// 
//** `intbas` represented as evaluations
// *
// * \todo currently experimental, not properly tested
// * \todo deal with case where intbas reaches degree = nbpoints
// **/
//void mbasis_rescomp_eval(
//                         Vec<Mat<zz_p>> & intbas,
//                         const Vec<Mat<zz_p>> & evals,
//                         const Vec<zz_p> & pts,
//                         VecLong & shift,
//                         long offset,
//                         long order
//                        );


/** Computes a `shift`-ordered weak Popov interpolant basis `intbas` for
 * `(pmat,pts)`, with the matrix `pmat` given by its list of evaluations
 * `evals` at all the points in `pts` (see @ref MBasisInt). Tries to choose the
 * fastest of the available `mbasis` variants. The positive integers `offset`
 * and `order` indicate that we consider the entries `offset, offset+1, ..,
 * offset+order-1` of `evals` and `pts` (no check is performed to verify that
 * these indices stay within the allowed bounds).
 *
 * Requirements: this function assumes that there are no repeated points in
 * `pts` (undefined behaviour otherwise), and that `evals` and `pts` have the
 * same strictly positive length, and that all matrices in `evals` have the
 * same dimensions.
 *
 * \todo tune threshold (currently only ensures to choose the residual update
 * variant in the "almost square" case)
 *
 * \todo add checks and allow repeated points
 * 
 **/
inline void mbasis(
                   Mat<zz_pX> & intbas,
                   Vec<Mat<zz_p>> & evals,
                   const Vec<zz_p> & pts,
                   VecLong & shift,
                   long offset,
                   long order
                  )
{
    if (order <= 5 || (evals[0].NumCols() == evals[0].NumRows()-1 && evals[0].NumRows()>20))
        mbasis_resupdate(intbas,evals,pts,shift,offset,order);
    else
        mbasis_rescomp(intbas,evals,pts,shift,offset,order);
}

/** Computes a `shift`-ordered weak Popov interpolant basis `intbas` for
 * `(pmat,pts)`, with the matrix `pmat` given by its list of evaluations
 * `evals` at all the points in `pts`. Tries to choose the fastest of the
 * available `mbasis` variants (see @ref MBasisInt).
 *
 * Requirements: this function assumes that there are no repeated points in
 * `pts` (undefined behaviour otherwise), and that `evals` and `pts` have the
 * same strictly positive length, and that all matrices in `evals` have the
 * same dimensions.
 *
 * \todo tune threshold (currently only ensures to choose the residual update
 * variant in the "almost square" case)
 *
 * \todo add checks and allow repeated points
 * 
 **/
inline void mbasis(
                   Mat<zz_pX> & intbas,
                   Vec<Mat<zz_p>> & evals,
                   const Vec<zz_p> & pts,
                   VecLong & shift
                  )
{
    if (pts.length() <= 5 || (evals[0].NumCols() == evals[0].NumRows()-1 && evals[0].NumRows()>20))
        mbasis_resupdate(intbas,evals,pts,shift,0,pts.length());
    else
        mbasis_rescomp(intbas,evals,pts,shift,0,pts.length());
}


/** Computes a `shift`-Popov interpolant basis `intbas` for `(pmat,pts)`, with
 * the matrix `pmat` given by its list of evaluations `evals` at all the points
 * in `pts`.
 *
 * Requirements: this function assumes that there are no repeated points in
 * `pts` (undefined behaviour otherwise), and that `evals` and `pts` have the
 * same strictly positive length, and that all matrices in `evals` have the
 * same dimensions.
 *
 * \todo tune threshold (currently only ensures to choose the residual update
 * variant in the "almost square" case)
 *
 * \todo add checks and allow repeated points
 * 
 **/
void popov_mbasis(
                  Mat<zz_pX> & intbas,
                  Vec<Mat<zz_p>> & evals,
                  const Vec<zz_p> & pts,
                  VecLong & shift
                 );

//@} // doxygen group: M-Basis algorithm


/** @name PM-Basis algorithm
 * \anchor PMBasisInt
 *
 * These functions compute a `shift`-minimal ordered weak Popov interpolant
 * basis for `(pmat,pts)`. They use a divide and conquer approach , computing a
 * first basis for the first half of the points, finding the so-called
 * _residual matrix_, computing a second basis for the remaining half of the
 * points, and deducing the sought basis by multiplying the two obtained bases.
 *
 * At the end of the computation, the vector `shift` contains the shifted row
 * degree of `intbas`, for the input shift. 
 *
 * The first recursive call returns an interpolant basis `intbas1` such that
 * `intbas1*pmat` vanishes at the first half of the points, and the residual
 * matrix has the same dimensions as `pmat` and is a matrix whose evaluations
 * at the second half of the points are the same as the evaluations of
 * `intbas1*pmat` at these points.
 */
//@{

/** Computes a `shift`-ordered weak Popov interpolant basis `intbas` for
 * `(pmat,pts)` (see @ref PMBasisInt).
 *
 * Requirements: this function assumes that there are no repeated points in
 * `pts` (undefined behaviour otherwise).
 *
 * \todo tune threshold for call to `mbasis`
 * \todo add checks and allow repeated points
 * 
 **/
void pmbasis(
             Mat<zz_pX> & intbas,
             const Mat<zz_pX> & pmat,
             const Vec<zz_p> & pts,
             VecLong & shift
            );

/** Computes a `shift`-ordered weak Popov interpolant basis `intbas` for
 * `(pmat,pts)`, where the points are geometric sequence defined by `r` and
 * length `order` (see @ref PMBasisInt). This fills the vector `pts` with the
 * actual list of points.
 *
 * Requirement: the degree of `pmat` must be less than `order` (if that is not
 * the case, one may reduce `pmat` modulo the `order` interpolation points).
 *
 * \todo tune threshold for call to `mbasis`
 * 
 **/
void pmbasis_geometric(
                       Mat<zz_pX> & intbas,
                       const Mat<zz_pX> & pmat,
                       const zz_p & r,
                       const long order,
                       VecLong & shift,
                       Vec<zz_p> & pts
                      );

/** Computes a `shift`-ordered weak Popov interpolant basis `intbas` for
 * `(pmat,pts)` (see @ref PMBasisInt), where the points `pts` are geometric
 * sequence defined by `r` and length `order`, and `evals` is the list of
 * evaluations of `pmat` at these points. The positive integers `offset` and
 * `order` indicate that we consider the entries `offset, offset+1, ..,
 * offset+order-1` of `evals` and `pts` (no check is performed to verify that
 * these indices stay within the allowed bounds).
 *
 * Note that `evals` is not `const`.
 *
 * \todo tune threshold for call to `mbasis`
 * 
 **/
void pmbasis_geometric(
                       Mat<zz_pX> & intbas,
                       Vec<Mat<zz_p>> & evals,
                       const Vec<zz_p> & pts,
                       const zz_p & r,
                       VecLong & shift,
                       long offset,
                       long order
                      );


/** Computes a `shift`-ordered weak Popov interpolant basis `intbas` for
 * `(pmat,pts)`, where `evals` is the list of evaluations of `pmat` at these
 * points (see @ref PMBasisInt). The positive integers `offset` and `order`
 * indicate that we consider the entries `offset, offset+1, .., offset+order-1`
 * of `evals` and `pts` (no check is performed to verify that these indices
 * stay within the allowed bounds).
 *
 * Note that `evals` is not `const`.
 *
 * \todo tune threshold for call to `mbasis`
 * 
 **/
void pmbasis(
             Mat<zz_pX> & intbas,
             Vec<Mat<zz_p>> & evals,
             const Vec<zz_p> & pts,
             VecLong & shift,
             long offset,
             long order
            );

/** Computes a `shift`-Popov interpolant basis `intbas` for `(pmat,pts)`.
 *
 * \todo tune threshold for call to `mbasis`
 **/
void popov_pmbasis(
                   Mat<zz_pX> & intbas,
                   const Mat<zz_pX> & pmat,
                   const Vec<zz_p> & pts,
                   VecLong & shift
                  );

//@} // doxygen group: PM-Basis algorithm


// list of points: for each column, we have a list of points (elt from zz_p and
// multiplicity)
// FIXME not thought thorougly yet, subject to change
//typedef std::vector<std::vector<std::pair<zz_p,long>>> Points;

/*------------------------------------------------------------*/
/* general user-friendly interface                            */
/*------------------------------------------------------------*/

// TODO
//VecLong interpolant_basis(
//                         Mat<zz_pX> & intbas,
//                         const Mat<zz_pX> & pmat,
//                         const Points & pts,
//                         const VecLong & shift = VecLong(),
//                         const PolMatForm form = ORD_WEAK_POPOV,
//                         const bool row_wise = true,
//                         const bool generic = false
//                        );

// below, the following is called "uniform interpolation case" (or case with uniform points)
// this means that we have the same points on all columns, all with multiplicity one)
// (FIXME could be easily generalized to any constant multiplicity for all...?)
//VecLong interpolant_basis(
//                         Mat<zz_pX> & intbas,
//                         const Mat<zz_pX> & pmat,
//                         const Vec<zz_p> & pts,
//                         const VecLong & shift = VecLong(),
//                         const PolMatForm form = ORD_WEAK_POPOV,
//                         const bool row_wise = true,
//                         const bool generic = false
//                        );
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
/* Iterative algorithm for arbitrary points and shift         */
/* References:                                                */
/*   - Beckermann 1992                                        */
/*   - Van Barel-Bultheel 1991+1992                           */
/*   - Beckermann-Labahn 2000 (ensuring s-Popov)              */
/*------------------------------------------------------------*/
//VecLong intbas_iterative(
//                        Mat<zz_pX> & intbas,
//                        const Mat<zz_pX> & pmat,
//                        const Points & pts,
//                        const VecLong & shift,
//                        bool point_wise=true // TODO to be thought about
//                       );
//
//VecLong popov_intbas_iterative(
//                              Mat<zz_pX> & intbas,
//                              const Mat<zz_pX> & pmat,
//                              const Points & pts,
//                              const VecLong & shift,
//                              bool point_wise=true // TODO to be thought about
//                             );



#endif /* ifndef MAT_LZZ_PX_INTERPOLANT__H */

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
