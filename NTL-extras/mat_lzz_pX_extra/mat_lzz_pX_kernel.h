#ifndef MAT_LZZ_PX_KERNEL__H
#define MAT_LZZ_PX_KERNEL__H

#include <NTL/matrix.h>
#include <NTL/lzz_pX.h>
#include "mat_lzz_pX_forms.h" // for VecLong, VecLong, PolMatForm

NTL_CLIENT

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*                     KERNEL BASIS                           */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

////Definition (kernel basis)
// Given an m x n matrix of univariate polynomials 'pmat',
// a (left) kernel basis for pmat is a matrix over K[X]
// whose rows form a basis for the (left) kernel of pmat, that is,
// the K[X]-module { v in K[X]^{1 x m}  |  v * pmat = 0 }.

// Degree bounds for kernel bases: assuming pmat has full column rank n, an
// unpublished result by Vu Thi Xuan (Master's research report, Lemma 10) shows
// that the sum of the pivot degrees of the shifted Popov kernel basis (for any
// shift) is at most the degree of the determinant of the complement part in
// 'pmat' (complement: rows not in the shifted pivot support of the kernel).
// This can probably be derived from knowledge about irreducible fractions.
//
// As a result, we have the following non-strict upper bounds on the sum of
// pivot degree of the shifted Popov kernel basis of pmat:
//   * pmat.NumCols() * degree(pmat)
//   * sum of column degrees of pmat
//   * sum of degrees of the rows in the complement of pmat
//          (useful if pivot support is known)
//   * if pivot support is unknown, since the complement has <= n rows, this
//   can be relaxed as: sum of degrees of the n largest-degree rows of pmat.
//
// The max of the pivot degrees, which is also the degree of the pivot part of
// the shifted Popov kernel, can be equal to the sum of pivot degrees (although
// this is not expected generically). Hence we have the same bounds on degree
// of pivot part. This directly gives bounds on the non-pivot part: pick one of
// the bounds above, and add max(shift) - min(shift). Note that this might be
// pessimistic for shifts of large amplitude, since we also have bounds on the
// non-pivot part involving the adjugate of the pivot support complement part
// of pmat, ensuring in particular that the non-pivot part of the kernel has
// degree at most n * deg(pmat)   (cf Lemma 12 of Vu Thi Xuan's research
// report).
//
// The max of pivot degrees of the kernel also gives a bound on
// max(rdeg_s(kernel)): pick one of the bounds above for sum of pivot degrees,
// and add max(s).
//
// Case where pmat does not have full column rank: the left kernel of pmat is
// equal to the left kernel of any column basis of pmat; taking a column
// reduced form of pmat allows us to reduce to the full column rank case. This
// shows that the bounds above still hold without full column rank assumption.

/*------------------------------------------------------------*/
/* general user-friendly interface                            */
/*------------------------------------------------------------*/
// TODO options for row-wise, normal form, etc
// TODO thresholds, ..
VecLong kernel_basis(
                    Mat<zz_pX> & kerbas,
                    const Mat<zz_pX> & pmat,
                    const VecLong & shift
                   );

/*------------------------------------------------------------*/
/* Verification that a matrix is a shifted-minimal kernel     */
/* basis                                                      */
/*------------------------------------------------------------*/
bool is_kernel_basis(
                     Mat<zz_pX> & kerbas,
                     const Mat<zz_pX> & pmat,
                     const VecLong & shift,
                     const PolMatForm & form = ORD_WEAK_POPOV,
                     const bool randomized = false
                    );

/*------------------------------------------------------------*/
/* Kernel basis: naive via large order approximant basis      */
/*------------------------------------------------------------*/
// TODO describe input-output?
// can this be faster than ZLS for some cases? (e.g. if shifts does not
// correspond at all to row degrees of pmat, or generally if shifts are "bad")

// The order of approximation is designed as follows.
VecLong kernel_basis_via_approximation(
                                      Mat<zz_pX> & kerbas,
                                      const Mat<zz_pX> & pmat,
                                      const VecLong & shift
                                     );
// TODO same via interpolant?

/*------------------------------------------------------------*/
/* Kernel basis: Zhou-Labahn Storjohann algorithm,            */
/* original version via approximant bases                     */
/*------------------------------------------------------------*/
// TODO describe input-output?
VecLong kernel_basis_zls_via_approximation(
                                          Mat<zz_pX> & kerbas,
                                          const Mat<zz_pX> & pmat,
                                          const VecLong & shift
                                         );

/*------------------------------------------------------------*/
/* Kernel basis: Zhou-Labahn Storjohann algorithm,            */
/* modified version via interpolation bases                   */
/*------------------------------------------------------------*/
// TODO describe input-output?
VecLong kernel_basis_zls_via_interpolation(
                                          Mat<zz_pX> & kerbas,
                                          const Mat<zz_pX> & pmat,
                                          const VecLong & shift
                                         );

// TODO generic case

// TODO general case via fast s-Popov appbas (fastest known approach for bad shifts)


/*************************************************
 *  Fast left kernel for small column dimension  *
 *************************************************/

// TODO-long shot via relation basis (itself todo), cf work with Xuan

#endif /* ifndef MAT_LZZ_PX_KERNEL__H */

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
