#ifndef MAT_LZZ_PX_EXTRA__H
#define MAT_LZZ_PX_EXTRA__H

/** Main header for `Mat<zz_pX>`, matrices over the univariate polynomials.
 *
 * \file mat_lzz_pX_extra.h
 * \author Seung Gyu Hyun, Vincent Neiger, Eric Schost
 * \version 0.1
 * \date 2018-12-07
 *
 * This is the main header for polynomial matrix functions. It mostly includes
 * other headers which gather functions for a specific kind of tasks. This file
 * contains general TODOs, and still contains a few declarations that do not
 * require a separate header (for the moment).
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

// TODO general reduction to uniform shift via pre-multiplication
// worthwile at least when shift close to uniform

// TODO naive algorithms (see Mulders-Storjohann for good reference)

// TODO general shifted Popov form via kernel (itself via approximant basis)

// TODO understand if there is any chance Alekhnovich improves over the
// kernel approach

// TODO nonsingular: Giorgi-Jeannerod-Villard's Las Vegas reduction
// (worth implementing for shifts other than uniform?)



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
