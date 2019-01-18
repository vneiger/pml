#ifndef MAT_LZZ_PX_DETERMINANT__H
#define MAT_LZZ_PX_DETERMINANT__H

/** Determinant.
 *
 * \file mat_lzz_pX_determinant.h
 * \author Seung Gyu Hyun, Vincent Neiger, Eric Schost
 * \version 0.1
 * \date 2019-01-01
 *
 */

#include <NTL/matrix.h>
#include <NTL/lzz_pX.h>

NTL_CLIENT

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
//bool determinant_generic_las_vegas(zz_pX & det, const Mat<zz_pX> & pmat);

// Version 3 (randomized; via random linear system solving)
// TODO first version, should be improved. Make Las Vegas.
void determinant_via_linsolve(zz_pX & det, const Mat<zz_pX> & pmat);

// Determinant via multipoint evaluation, general points. Requires the field to
// be sufficiently large (namely, with the current version, up to
// deg(pmat)*pmat.NumRows()+1 points are required, independently of the degree
// profile of pmat).
//
// TODO improve with better bounds on degdet, such as sum of cdeg or rdeg, or
// even the generic-degdet bound)
void determinant_via_evaluation_general(zz_pX & det, const Mat<zz_pX> & pmat);

// Determinant via multipoint evaluation, geometric points. Requires the field to
// be sufficiently large (namely, with the current version, ????? up to
// deg(pmat)*pmat.NumRows()+1 points are required, independently of the degree
// profile of pmat).
//
// TODO improve with better bounds on degdet, such as sum of cdeg or rdeg, or
// even the generic-degdet bound)
void determinant_via_evaluation_geometric(zz_pX & det, const Mat<zz_pX> & pmat);



// TODO other determinant algorithms??
// --> could rely on x-Smith decomposition of Gupta et al (worth
// implementing??), cf Appendix of LaNeZh17 

#endif // MAT_LZZ_PX_DETERMINANT__H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
