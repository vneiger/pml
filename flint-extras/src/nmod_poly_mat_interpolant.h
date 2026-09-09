/*
    Copyright (C) 2026 Gilles Villard

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/


#include "pml.h"

/** \file nmod_poly_mat_interpolant.h
 * Definition (interpolant basis).
 * --------------------------------
 * Consider:
 *   - an m x n matrix of univariate polynomials, given as `d` matrices
 *     `E = (E_1,...,E_d)` in `K^{m x n}` (the coefficients of E at `d`
 *     specific, points `pts = (pts_1,...,pts_d)` (not)
 *     the coefficients of a polynomial matrix in the usual monomial
 *     basis),
 *   - the points `pts` themselves.
 * 
 * No assumption such that pairwise distinct points even though, for now, 
 * applications of the non-distinct case are not included in pml.
 *
 * Then an interpolant basis for `(E,pts)` is a matrix over the univariate
 * polynomials whose rows form a basis for the following module:
 * { p in K[X]^{1 x m}  |  p(pts_k) * E_k = 0 for 1 <= k <= d }.
 * Note that such a matrix is square, m x m, and nonsingular.
 *
 *  The length of pts is at least d.
 * 
 * `E` is a plain, flat array of `nmod_mat_struct`, its length is also at least d.
 * 
 * `d`, the number of points to use, is a separate explicit parameter 
 *    with  d <= (actual length of pts) and  d <= (actual length of E)
 * 
 *  pts and E may be longer than d, extra data is ignored by convention.
 * 
 * This is the point-evaluation analogue of an approximant basis (see
 * `nmod_poly_mat_approximant.h`): replacing "coefficient of `P*F`" (order
 * truncation) with "evaluation `P(pts_k)*E_k`" (interpolation) in the same
 * definition recovers interpolant bases from approximant bases, and the
 * same iterative/divide-and-conquer algorithmic scheme applies to both
 * (HNS19: S. Hyun, V. Neiger, E. Schost. Proceedings ISSAC 2019, sec 3.2).
 * `I_pts(F(pts_1),...,F(pts_d))` coincides with the approximant module
 * `A_M(F)` for `M = prod(x - pts_k)`.
 */

/** \file nmod_poly_mat_interpolant.h
 * Definition (shifted minimal interpolant basis).
 * -------------------------------------------------
 * Starting from the definition of an interpolant basis, consider further:
 *   - a degree shift s (a list of m integers).
 *
 * Then an interpolant basis for `(E,pts)` is said to be <em>a shift-minimal</em>
 * (resp. <em>a shift-ordered weak Popov</em>, resp. <em>the shift-Popov</em>)
 * interpolant basis if it is in shift-reduced form (resp. in shift-ordered
 * weak Popov form, resp. in shift-Popov form). See nmod_poly_mat_forms.h
 * for definitions of these forms.
 */

/** \file nmod_poly_mat_interpolant.h
 * Conventions.
 * ------------
 * As in `nmod_poly_mat_approximant.h`, all functions below compute left
 * interpolant bases (the basis elements are the rows of the output
 * matrix).
 *
 * Most functions below use the following parameters.
 *
 * \param[out] intbas the output interpolant basis (cannot alias `E`)
 * \param[in] pts the `d` interpolation points.
 * \param[in] `E` is a plain, flat array of `nmod_mat_struct` (`nmod_mat_struct *`),
 *  a container for a sequence of constant matrices, (the
 *   coefficients of `E` at the points `pts`, not in the monomial basis), 
 *  its length is at least d.
 * \param[in] `d` is a separate parameter SUCH THAT 
 *    `d <= (actual length of pts)`  and `d <= (actual length of E)`.
 * \param[in,out] shift in: the input shift; and out: the output shifted row
 *   degree of `intbas` (list of integers, length must be the number of
 *   rows of `E`).
 * 
 * Sources.
 * --------
 * The algorithms are inspired from: 
 *   - B. Beckermann and G. Labahn. 2000. Fraction-free computation of matrix 
 *     rational interpolant and matrix gcds. 
 *     SIAM J. Matrix Anal. Appl. 22, 1 (2000), 114–144.
 *   - C.-P. Jeannerod, V. Neiger, É. Schost, and G. Villard. 2017. 
 *     Computing minimal interpolation bases. 
 *     J. Symbolic Comput. 83 (2017), 272–314.
 * and can be found in 
 *  - PML/ntl-extras
 *  - S. Hyun, V. Neiger, E. Schost. Proceedings ISSAC 2019. 
 *
 */


#ifndef NMOD_POLY_MAT_INTERPOLANT_H
#define NMOD_POLY_MAT_INTERPOLANT_H

#define PMINTBASIS_THRES 32

#include <flint/nmod_types.h>
#include <flint/nmod_poly_mat.h>

#include "nmod_poly_mat_forms.h" // for orientation_t
#include "pml.h"

#ifdef __cplusplus
extern "C" {
#endif

/** @name General interface for interpolant basis verification*/
//@{

/** Verifying if a matrix is a minimal interpolant basis.
 *
 * This checks whether the matrix `intbas` is a `shift`-minimal interpolant
 * basis for (`E`,`pts`) for the required form (currently hardcoded to
 * ordered weak Popov, matching this same current limitation in PML's own
 * `nmod_poly_mat_is_approximant_basis`, `nmod_poly_mat_approximant.h`).
 * TO SEE.
 * 
 * The direct interpolant-basis counterpart of that function -- see
 * `verification.c` (this same staging area) for the full derivation of
 * how its two order-truncation-specific pieces (membership, generation)
 * were adapted to points. 
 *
 * \param[in] intbas interpolant basis
 * \param[in] pts the `d` interpolation points
 * \param[in] E the input sequence, representing `d` constant matrices (see this
 *   header's own "Conventions" section)
 * \param[in] d number of points to actually consider, such that 
 *   `d <= (actual length of pts)` and `d <= (actual length of E)`.
 * \param[in] shift shift
 * \param[in] orient indicates the orientation (left/right interpolants)
 *   and the definition of pivots
 *
 * \return boolean, result of the verification
 */
int nmod_poly_mat_is_interpolant_basis(const nmod_poly_mat_t intbas,
                                       const ulong * pts,
                                       const nmod_mat_struct * E,
                                       slong d,
                                       const slong * shift,
                                       orientation_t orient);

//@} // doxygen group: General interface for interpolant basis verification

/** @name M-IntBasis algorithm (uniform number of interpolation points)
 *
 * The core of these functions is implemented with `nmod_mat_poly_t` type,
 * for efficiency reasons, mirroring `mbasis`'s own choice (see @ref
 * mintbasis in `nmod_mat_poly.h` for the mat_poly-level M-IntBasis
 * variants and their doc comments -- this function is a thin wrapper
 * around that dispatcher, converting to/from the entrywise `nmod_poly_mat_t`
 * representation used at this level, exactly mirroring
 * @ref nmod_poly_mat_mbasis's own role for `mbasis`).
 */
//@{

/** Computes a `shift`-ordered weak Popov interpolant basis for
 * `(E,pts)`, `E` given as a flat array of `nmod_mat_struct`, 
 *  with `E_k` stored as the coefficient of degree `k` */
void nmod_poly_mat_mintbasis(nmod_poly_mat_t intbas,
                             slong * shift,
                             const ulong * pts,
                             const nmod_mat_struct * E,
                             slong d);

//@} // doxygen group: M-IntBasis algorithm (uniform number of interpolation points)


/** @name PM-IntBasis algorithm (uniform number of interpolation points)
 * \anchor pmintbasis
 *
 * These functions compute a `shift`-minimal ordered weak Popov interpolant
 * basis for `(E,pts)`, `d` = number of points, using a divide and conquer
 * approach mirroring `pmbasis` (see @ref pmbasis,
 * `nmod_poly_mat_approximant.h`): split the point list in half (rather
 * than the order), compute a first basis for the first half of the
 * points, find the residual matrix for the second half, compute a second
 * basis for it, and deduce the sought basis by multiplying the two.
 *
 * The residual here is `d2` independent evaluate-then-multiply operations, 
 * `R_i = P1(pts[d1+i])*E_{d1+i}`, `i=0..d2-1` -- no analogue of
 * `pmbasis`'s middle product for interpolation.
 *
 * At the end of the computation, the vector `shift` contains the shifted
 * row degree of `intbas`, for the input shift.
 */
//@{

/** Computes a `shift`-ordered weak Popov interpolant basis for `(E,pts)`
 * using the algorithm PM-IntBasis (see @ref pmintbasis). Falls back to
 * @ref nmod_poly_mat_mintbasis once `d <= PMINTBASIS_THRES`. */
void nmod_poly_mat_pmintbasis(nmod_poly_mat_t intbas,
                              slong * shift,
                              const ulong * pts,
                              const nmod_mat_struct * E,
                              slong d);

//@} // doxygen group: PM-IntBasis algorithm (uniform number of interpolation points)


/** @name PM-IntBasis algorithm, geometric points
 * \anchor pmintbasis_geometric
 *
 * Same divide-and-conquer shape as @ref pmintbasis, specialized to
 * geometric points `pts_k = r^{2k}` (`r` (r a field element; its order 
 * governs whether the points are distinct — see the functions below)
 * via FLINT's Bostan-Schost fast geometric-progression evaluate/interpolate 
 * machinery (`nmod_geometric_progression_t`,`nmod_poly.h`). 
 * One such structure is built once, at the top of the
 * whole recursion, and reused (at varying, always-shorter requested
 * lengths). 
 */
//@{

/** Computes a `shift`-ordered weak Popov interpolant basis for `(E,r,d)`,
 * where `E` is given at the `d` geometric points `1, r^2, r^4, ...,
 * r^{2(d-1)}` (`E_k` stored as the coefficient of degree `k`, see this
 * header's "Conventions" section), using the algorithm PM-IntBasis
 * specialized to geometric points (see @ref pmintbasis_geometric).
 *
 * No condition on the modulus is needed or imposed: this algorithm only
 * evaluates at the points, never interpolates, so it never divides by
 * `r^{2k}-1` and cannot fail on a small field. If `r` has too low an
 * order the points simply repeat, and the result is a correct basis for
 * the interpolant module those repeated points define. For the usual `d`
 * pairwise distinct points -- and with them the classical degree property
 * `deg(det(intbas)) = sum_k rank(E_k)` -- the caller must supply `r` with
 * `rho = r^2` of order at least `d`; @ref 
 * nmod_poly_mat_pmintbasis_geometric_auto picks such an `r` itself.
 *
 * If `pts` is non-null it is filled with the `d` points actually used
 * (`pts[k] = r^{2k}`), matching @ref nmod_poly_mat_pmintbasis's own
 * point-array convention -- useful for cross-checking against the
 * general-points algorithm on the same instance. */
void nmod_poly_mat_pmintbasis_geometric(nmod_poly_mat_t intbas,
                                        slong * shift,
                                        ulong * pts,
                                        const nmod_mat_struct * E,
                                        ulong r,
                                        slong d);


/** Same as @ref nmod_poly_mat_pmintbasis_geometric, except that `r` is
 * found internally instead of being supplied by the caller -- targeting
 * the pairwise distinct points case.
 *
 * What is needed is `rho = r^2` of order at least `d`, giving `d` distinct
 * points; since `ord(r^2)` is `ord(r)` divided by `gcd(2,ord(r))`, i.e. at
 * worst halved, an `r` of order at least `2*d-1` suffices. Such an `r` is
 * obtained from `nmod_find_root(2*d, mod)` (`nmod_extra.h`), documented to
 * return an element of order at least `2*d` -- cheaper than a genuine
 * primitive root, since unlike `n_primitive_root_prime` it needs no
 * factoring of `p-1`.
 *
 * Passing `2*d-1` would already suffice; `2*d` is used both for one unit
 * of margin and because `nmod_find_root` then returns 0 exactly when
 * `p <= 2*d+1`, that is, exactly when the precondition below fails --
 * which is why `n_primitive_root_prime`, kept as a commented fallback, 
 * is in fact unreachable. `2*d`also matches the 2*len convention 
 * used by mul_geometric.c and mulmid.c",
 *
 * Requires a modulus `p > 2*d+1`, and throws otherwise. That this is a
 * condition on `p` at all -- the `r`-explicit version above has none --
 * is precisely because `r` must be found here: the best possible `r` is a
 * primitive root, giving `ord(r^2) = (p-1)/2`, so `d` distinct points
 * exist only when `(p-1)/2 >= d`. One unit of margin is
 * kept here too: `p = 2*d+1` does work (`ord(r^2) = d` exactly, `d` points
 * and none to spare) but is excluded, both for that room and so that the
 * precondition is the exact complement of `nmod_find_root(2*d)`'s own
 * failure condition `p <= 2*d+1`. `d = 0` needs no `r` (it returns before
 * one is used) and so never throws, however small the modulus. */
void nmod_poly_mat_pmintbasis_geometric_auto(nmod_poly_mat_t intbas,
                                             slong * shift,
                                             ulong * pts,
                                             const nmod_mat_struct * E,
                                             slong d);

//@} // doxygen group: PM-IntBasis algorithm, geometric points




#ifdef __cplusplus
}
#endif

#endif // NMOD_POLY_MAT_INTERPOLANT_H
