#ifndef MAT_LZZ_PX_DETERMINANT__H
#define MAT_LZZ_PX_DETERMINANT__H

/** \brief Determinant algorithms for polynomial matrices over `zz_p`.
 *
 * \file mat_lzz_pX_determinant.h
 * \author Seung Gyu Hyun, Vincent Neiger, Eric Schost
 * \version 0.1
 * \date 2019-01-01
 *
 * Implement determinant algorithms: via linear system solving, via
 * triangularization, via evaluation/interpolation, and via expansion by
 * minors.
 *
 * \todo finish non-plain expansion by minors for dimensions 5, 6 (7 is
 * probably too much?).
 *
 * \todo another algorithm: rely on x-Smith decomposition of Gupta et al, cf
 * Appendix of LaNeZh17 
 *
 */

#include <NTL/matrix.h>
#include <NTL/lzz_pX.h>

NTL_CLIENT

/** Computes the determinant `det` of `pmat`. Chooses the fastest available
 * option according to some given thresholds.
 *
 * \todo Not implemented yet (thresholding mechanism not done yet)
 */
void determinant(zz_pX & det, const Mat<zz_pX> & pmat);

/** Computes and returns the determinant `det` of `pmat`. Chooses the fastest
 * available option according to some given thresholds.
 *
 * \todo Not implemented yet
 */
inline zz_pX determinant(const Mat<zz_pX> & pmat)
{
    zz_pX det;
    determinant(det, pmat);
    return det;
}

/** Verifies whether `det` is the determinant of `pmat`. If `up_to_constant` is
 * `true`, the check is up to a multiplicative constant. If `randomized` is `true`,
 * a randomized algorithm may be used.
 *
 * \todo be more explicit on randomized (Monte Carlo ? one sided ?)
 *
 * \todo only randomized is implemented for now
 */
bool verify_determinant(const zz_pX & det, const Mat<zz_pX> & pmat, bool up_to_constant, bool randomized);

/*******************************************************************
 *  Labahn-Neiger-Zhou: via diagonal entries of triangularization  *
 *******************************************************************/

//void determinant_via_diagonal_of_hermite(zz_pX & det, const Mat<zz_pX> & pmat);

/** Tries to compute the determinant `det` of `pmat`, assuming its degree is
 * known and provided as input, by running [Labahn-Neiger-Zhou 2017] partial
 * triangularization algorithm, assuming all encountered row bases are
 * unimodular, hence ignoring them (this is correct for generic `pmat`).
 * Assuming the provided `degree` is correct, returns `true` if and only if the
 * computation was successful (i.e. `det` is the determinant of `pmat` up to a
 * constant). */
bool determinant_generic_knowing_degree(zz_pX & det, const Mat<zz_pX> & pmat, long degree);

// Version 2 (Las Vegas randomized algorithm): runs the partial
// triangularization and checks determinant by Zippel-Schwartz; returns false
// if determinant is wrong or if field size is too small for allowing the
// Zippel-Schwartz check. If true is returned, then determinant is correct.
// TODO: not implemented yet
//bool determinant_generic_las_vegas(zz_pX & det, const Mat<zz_pX> & pmat);

/** Computes the determinant `det` of `pmat` by solving a linear system with
 * random right-hand side of degree `deg(pmat)` and taking the denominator. 
 *
 * \todo this currently requires `pmat(0)` to be invertible; support other
 * cases by using kernel based system solving
 *
 * \todo determinant only up to constant?
 */
void determinant_via_linsolve(zz_pX & det, const Mat<zz_pX> & pmat);

/** Computes the determinant `det` of `pmat` by evaluation/interpolation at
 * general points. This requires the field to be sufficiently large (namely,
 * with the current version, `deg(pmat)*pmat.NumRows()+1` points are required,
 * independently of the degree profile of pmat).
 *
 * \todo improve with better bounds on degdet, such as sum of cdeg or rdeg, or
 * even the generic-degdet bound)
 */
void determinant_via_evaluation_general(zz_pX & det, const Mat<zz_pX> & pmat);

/** Computes the determinant `det` of `pmat` by evaluation/interpolation at
 * a geometric sequence. This requires the field to be sufficiently large.
 *
 * \todo improve doc: give required field size
 *
 * \todo improve this required size with better bounds on degdet, such as sum
 * of cdeg or rdeg, or even the generic-degdet bound)
 */
void determinant_via_evaluation_geometric(zz_pX & det, const Mat<zz_pX> & pmat);


/** Computes the determinant `det` of `pmat` by evaluation/interpolation at FFT
 * points. This requires the field to be sufficiently large.
 *
 * \todo improve doc: give required field size
 *
 * \todo improve this required size with better bounds on degdet, such as sum
 * of cdeg or rdeg, or even the generic-degdet bound)
 */
void determinant_via_evaluation_FFT(zz_pX & det, const Mat<zz_pX> & pmat);

/** Computes the determinant `det` of `pmat` by the plain expansion by minors
 * (recursive, with `n` recursive calls to submatrices of dimensions `n-1 x
 * n-1`, where `n` is the dimension of the input `pmat`). Should not be used
 * for dimensions above 5 or 6. */
void determinant_expansion_by_minors_rec(zz_pX & det, const Mat<zz_pX> & pmat);

/** Computes the determinant `det` of `pmat` by expansion by minors, avoiding
 * some redundant computations (currently, supports size up to `4 x 4`).
 */
void determinant_expansion_by_minors(zz_pX & det, const Mat<zz_pX> & pmat);

#endif // MAT_LZZ_PX_DETERMINANT__H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
