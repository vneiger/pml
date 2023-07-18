#ifndef MAT_LZZ_PX_EXTRA__H
#define MAT_LZZ_PX_EXTRA__H

/** \brief Main header for `Mat<zz_pX>`, matrices over the univariate polynomials.
 *
 * \file mat_lzz_pX_extra.h
 * \author Seung Gyu Hyun, Vincent Neiger, Eric Schost
 * \version 0.1
 * \date 2018-12-07
 *
 * This is the main header for polynomial matrix functions. Its purpose is only
 * to include all headers (which gather functions for a specific kind of
 * tasks). This file contains general TODOs, and still contains a few
 * declarations that do not require a separate header (for the moment).
 *
 */

#include <NTL/matrix.h>
#include <NTL/lzz_pX.h>
#include <iostream>
#include <vector> // std vector, for shifts, degrees, pivot indices

#include "lzz_p_extra.h"
#include "lzz_pX_CRT.h"
#include "thresholds_solve_lift.h"

#include "mat_lzz_pX_utils.h"
#include "mat_lzz_pX_forms.h"
#include "mat_lzz_pX_arith.h"
#include "mat_lzz_pX_linearization.h"
#include "mat_lzz_pX_multiply.h"
#include "mat_lzz_pX_inverse.h"
#include "mat_lzz_pX_linsolve.h"
#include "mat_lzz_pX_approximant.h"
#include "mat_lzz_pX_interpolant.h"
#include "mat_lzz_pX_kernel.h"
#include "mat_lzz_pX_determinant.h"

#include "mat_lzz_pX_sequence.h"  // TODO: still draft

NTL_CLIENT



/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* Division with remainder                                    */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/** Given a matrix `B` (square, row reduced) and a matrix `A` (same row
 * dimension as `B`), computes the quotient and remainder in the polynomial
 * matrix left division of A by B. Precisely, `A = B*Q + R`, with the row-wise
 * degree constraint that `row_degree(R) < row_degree(B)`.
 *
 * This is based on the classical fast division algorithm for univariate
 * polynomials: obtain the quotient by a multiplication of the reverse of `A`
 * by the truncated inverse of the reverse of `B`, and then deduce the
 * remainder. See for example [Neiger - Vu, ISSAC 2017] for a more detailed
 * description. Note that this currently does not incorporate optimizations
 * (partial linearization, ...) when the row degree of `B` is not balanced. */
void quo_rem(
             Mat<zz_pX> & Q,
             Mat<zz_pX> & R,
             const Mat<zz_pX> & A,
             const Mat<zz_pX> & B
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
 *                          TODO: BASIS REDUCTION                     *
 *            (shifted reduced form and shifted normal forms)         *
 **********************************************************************/

// possible TODOs:
//   - general reduction to uniform shift via pre-multiplication (worthwile at
//   least when shift close to uniform)
//   - iterative algorithms for reduced, Popov, Hermite (Mulders-Storjohann is
//   a good reference)
//   - general shifted Popov form via kernel (itself via approximant basis)

/** Computation of a row reduced form `reduced` of `pmat` by the Las Vegas
 * randomized algorithm of [Giorgi - Jeannerod - Villard, ISSAC 2003].
 * Requirements: `pmat` is square (not checked) and `pmat(0)` is invertible
 * (checked). Returns 0 if `pmat(0)` was not invertible, otherwise the
 * computation succeeds and this function returns 1.
 *
 * \todo if `pmat(0)` is not invertible, should this function try again by
 * shifting with a random point? (yes if this is negligible compared to the
 * rest, otherwise, leave this choice to the user)
 *
 * \todo any shifted variant better than the trivial reduction to the uniform
 * shift? A shifted variant may require more terms of the truncated inverse of
 * the matrix.
 */
long reduced_form_gjv(
                      Mat<zz_pX> & reduced,
                      const Mat<zz_pX> & pmat
                     );



#endif // MAT_LZZ_PX_EXTRA__H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
