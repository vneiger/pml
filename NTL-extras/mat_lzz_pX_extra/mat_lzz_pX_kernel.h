#ifndef MAT_LZZ_PX_KERNEL__H
#define MAT_LZZ_PX_KERNEL__H

/** Minimal kernel basis.
 *
 * \file mat_lzz_pX_kernel.h
 * \author Seung Gyu Hyun, Vincent Neiger, Eric Schost
 * \version 0.1
 * \date 2019-01-01
 *
 */

#include "mat_lzz_pX_forms.h" // for VecLong, VecLong, PolMatForm

NTL_CLIENT

/** \file mat_lzz_pX_kernel.h
 *
 * Definition (kernel basis).
 * --------------------------
 * For an m x n matrix of univariate polynomials A, a (left) kernel basis for A
 * is a matrix of univariate polynomials whose rows form a basis for the (left)
 * kernel of pmat, that is, of the module { v in K[X]^{1 x m}  |  v * A = 0 }.
 *
 */

/** \file mat_lzz_pX_kernel.h
 * Definition (shifted minimal kernel basis).
 * ------------------------------------------
 * Considering furthermore a degree shift s (a list of m integers), a kernel
 * basis for F is said to be <em>a shift-minimal</em> (resp.  <em>a
 * shift-ordered weak Popov</em>, resp. <em>the shift-Popov</em>) kernel basis
 * if it is in shift-reduced form (resp. in shift-ordered weak Popov form,
 * resp. in shift-Popov form). See mat_lzz_pX_forms.h for definitions of these
 * forms.
 */

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

/** Computes a `shift`-minimal kernel basis `kerbas` for `pmat`. General
 * interface, tries to choose the fastest available method.
 *
 * \todo options for row-wise, required form, etc
 * \todo thresholds (currently uses direct method via approximation)
 **/
VecLong kernel_basis(
                     Mat<zz_pX> & kerbas,
                     const Mat<zz_pX> & pmat,
                     const VecLong & shift
                    );

/** Verifies that `kerbas` is a `shift`-minimal kernel basis for `pmat` with
 * form at least `form`. Uses a Monte-Carlo randomized algorithm if
 * `randomized` is `true` (left and right projections, Freivalds-like); default
 * is `false`, that is, the check is deterministic.
 */
bool is_kernel_basis(
                     Mat<zz_pX> & kerbas,
                     const Mat<zz_pX> & pmat,
                     const VecLong & shift,
                     const PolMatForm & form = ORD_WEAK_POPOV,
                     const bool randomized = false
                    );

/** Computes a `shift`-minimal kernel basis `kerbas` for `pmat`, using a single
 * call to minimal appoximant basis at sufficiently large order.
 *
 * \todo same approach via interpolant basis may bring some speed up in some
 * cases.
 */
VecLong kernel_basis_via_approximation(
                                       Mat<zz_pX> & kerbas,
                                       const Mat<zz_pX> & pmat,
                                       const VecLong & shift
                                      );

/** Computes a `shift`-minimal kernel basis `kerbas` for `pmat` using the
 * Zhou-Labahn-Storjohann algorithm, as described in the Proceedings ISSAC
 * 2012. */
VecLong kernel_basis_zls_via_approximation(
                                           Mat<zz_pX> & kerbas,
                                           const Mat<zz_pX> & pmat,
                                           const VecLong & shift
                                          );

/** Computes a `shift`-minimal kernel basis `kerbas` for `pmat` using the
 * modified Zhou-Labahn-Storjohann algorithm (described in the Proceedings
 * ISSAC 2012), relying on interpolant bases rather than approximant bases. */
VecLong kernel_basis_zls_via_interpolation(
                                           Mat<zz_pX> & kerbas,
                                           const Mat<zz_pX> & pmat,
                                           const VecLong & shift
                                          );

#endif /* ifndef MAT_LZZ_PX_KERNEL__H */

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
