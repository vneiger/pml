#ifndef MAT_LZZ_PX_EXTRA__H
#define MAT_LZZ_PX_EXTRA__H

#include <NTL/matrix.h>
#include <NTL/lzz_pX.h>
#include <iostream>
#include <vector> // std vector, for shifts, degrees, pivot indices

#include "lzz_p_extra.h"
#include "lzz_pX_CRT.h"
#include "thresholds_matrix_multiply.h"
#include "thresholds_matrix_middle_product.h"
#include "thresholds_newton_inv_trunc.h"
#include "thresholds_solve_lift.h"

typedef std::vector<long> VecLong;

#include "mat_lzz_pX_utils.h"  // TODO: test
#include "mat_lzz_pX_forms.h"  // TODO: test, time
#include "mat_lzz_pX_arith.h"  // TODO: test, time
#include "mat_lzz_pX_linearization.h"  // TODO: still draft
#include "mat_lzz_pX_multiply.h"  // TODO: check if anything todo
#include "mat_lzz_pX_inverse.h"  // TODO: check if anything todo
#include "mat_lzz_pX_linsolve.h"  // TODO: check if anything todo
#include "mat_lzz_pX_approximant.h"  // TODO: improve, test, time
#include "mat_lzz_pX_interpolant.h"  // TODO: still draft
#include "mat_lzz_pX_kernel.h"  // TODO: still draft

#include "mat_lzz_pX_sequence.h"  // TODO: still draft


NTL_CLIENT

void multiply_evaluate_FFT_direct_no_ll(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);

// TODO in the future, these functions will be moved to "linearization/compression" files
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
/*------------------------------------------------------------*/
/* Division with remainder                                    */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

// TODO to be tested, improved, then moved to a relevant file

// computes Q,R such that A = BQ + R
void quo_rem(
             Mat<zz_pX> &Q, 
             Mat<zz_pX> &R, 
             const Mat<zz_pX> &A,
             const Mat<zz_pX> &B
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
