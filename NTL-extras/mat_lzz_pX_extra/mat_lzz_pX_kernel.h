#ifndef MAT_LZZ_PX_KERNEL__H
#define MAT_LZZ_PX_KERNEL__H

/** \brief Minimal kernel basis.
 *
 * \file mat_lzz_pX_kernel.h
 * \author Seung Gyu Hyun, Vincent Neiger, Eric Schost
 * \version 0.1
 * \date 2019-01-01
 *
 * Algorithms for computing shifted minimal kernel bases of polynomial
 * matrices.
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

/** Computes a `shift`-minimal kernel basis `kerbas` for `pmat`. General
 * interface, tries to choose the fastest available method.
 *
 * \todo options for row-wise, required form, etc
 * \todo thresholds (currently uses Zhou-Labahn-Storjohann via approximation)
 **/
void kernel_basis(
                  Mat<zz_pX> & kerbas,
                  const Mat<zz_pX> & pmat,
                  const VecLong & shift
                 );

/** Verifies that `kerbas` is a `shift`-minimal kernel basis for `pmat` with
 * form at least `form`. Uses a Monte-Carlo randomized algorithm if
 * `randomized` is `true` (left and right projections, Freivalds-like); default
 * is `false`, that is, the check is deterministic.
 *
 * \todo testing generation; support row wise
 */
bool is_kernel_basis(
                     Mat<zz_pX> & kerbas,
                     const Mat<zz_pX> & pmat,
                     const VecLong & shift,
                     const PolMatForm & form = ORD_WEAK_POPOV,
                     const bool randomized = false
                    );

/** Computes a `shift`-minimal kernel basis `kerbas` for `pmat`, using a single
 * call to minimal approximant basis at sufficiently large order. Computes the
 * `shift`-pivot index `pivind` of `kerbas`, and `shift` becomes the shifted
 * row degree of `kerbas` (for the input shift).
 *
 * \todo implement/use shift entries reduction (the approximation order depends
 * on these entries, and they may be reduced depending on the
 * degrees/dimensions of pmat)
 *
 * \todo better performance if `pmat` has unbalanced column degrees would be to
 * use a column-degree wise order (however, this is not handled by fast
 * approximant algorithms for now)
 */
void kernel_basis_via_approximation(
                                    Mat<zz_pX> & kerbas,
                                    VecLong & pivind,
                                    const Mat<zz_pX> & pmat,
                                    VecLong & shift
                                   );

/** Computes a `shift`-minimal kernel basis `kerbas` for `pmat`, using a single
 * call to minimal interpolant basis at sufficiently long sequence of geometric
 * points. Computes the `shift`-pivot index `pivind` of `kerbas`, and `shift`
 * becomes the shifted row degree of `kerbas` (for the input shift).
 *
 * \todo implement/use shift entries reduction (the approximation order depends
 * on these entries, and they may be reduced depending on the
 * degrees/dimensions of pmat)
 */
void kernel_basis_via_interpolation(
                                    Mat<zz_pX> & kerbas,
                                    VecLong & pivind,
                                    const Mat<zz_pX> & pmat,
                                    VecLong & shift
                                   );

/** Computes a `shift`-minimal kernel basis `kerbas` for `pmat` using the
 * Zhou-Labahn-Storjohann algorithm (described in the Proceedings ISSAC 2012).
 * At the end of the computation, `shift` is the shifted row degree of `kerbas`
 * (for the input shift).
 *
 * \todo implement/use shift entries reduction (the approximation order depends
 * on these entries, and they may be reduced depending on the
 * degrees/dimensions of pmat)
 * \todo 
 */
void kernel_basis_zls_via_approximation(
                                        Mat<zz_pX> & kerbas,
                                        Mat<zz_pX> & pmat,
                                        VecLong & shift
                                       );


/** Computes a `shift`-minimal kernel basis `kerbas` for `pmat` using a
 * modified Zhou-Labahn-Storjohann algorithm (described in the Proceedings
 * ISSAC 2012), relying on interpolant bases rather than approximant bases.
 * At the end of the computation, `shift` is the shifted row degree of `kerbas`
 * (for the input shift).
 */
void kernel_basis_zls_via_interpolation(
                                        Mat<zz_pX> & kerbas,
                                        Mat<zz_pX> & pmat,
                                        VecLong & shift
                                       );

#endif /* ifndef MAT_LZZ_PX_KERNEL__H */

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
