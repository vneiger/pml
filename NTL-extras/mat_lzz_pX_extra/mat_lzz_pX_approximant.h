#ifndef MAT_LZZ_PX_APPROXIMANT__H
#define MAT_LZZ_PX_APPROXIMANT__H

#include <NTL/matrix.h>
#include <NTL/lzz_pX.h>
#include <iostream>
#include <vector>
#include <numeric> // for 'accumulate'

#include "mat_lzz_pX_extra.h"

NTL_CLIENT


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* MINIMAL APPROXIMANT BASES                                  */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

typedef std::vector<long> Order;

// TODO write proper docstrings

////Definition (approximant basis)
// Given:
//   * m x n matrix of univariate polynomials 'pmat',
//   * approximation order 'order' (list of n positive integers),
// An approximant basis for (pmat,order) is a matrix over K[X]
// whose rows form a basis for the K[X]-module
// { 'app' in K[X]^{1 x m}  |  the column j of 'app' 'pmat' is 0 modulo X^{order[j]} }

//// Minimal and Popov approximant bases
// Given in addition:
//   * a degree shift 'shifts' (list of m integers)
// then an approximant basis for (pmat,order) is said to be
// "a shift-minimal" (resp. "the shift-Popov") approximant basis
// if it shift-reduced (resp. in shift-Popov form)
// Idem for shift-ordered weak Popov
// Cf. literature for definitions

/*------------------------------------------------------------*/
/* TODO: documentation to explain what this computes and what */
/* are the options                                            */
/*------------------------------------------------------------*/

// TODO generic: from shift, deduce the pivot degree expected generically and
// use this as a shift instead of 'shift', then obtain directly Popov approx
// basis; skipping the first call to find the pivot degree (in addition, will
// this be more efficient because shift is nicer?)

// TODO pmbasis: threshold mbasis

// TODO mbasis insert threads

// Guarantee: output is at least ordered weak Popov
// return value is pivot degree

// TODO Names? kernel / approximant ? kernel basis / approximant basis ?
// keep interpolation or choose better name?

/*------------------------------------------------------------*/
/* general user-friendly interface                            */
/*------------------------------------------------------------*/

// TODO: when algorithms are ready, find threshold and write the definition
DegVec approximant_basis(
                         Mat<zz_pX> & appbas,
                         const Mat<zz_pX> & pmat,
                         const Order & order,
                         const Shift & shift = Shift(),
                         const PolMatForm form = ORD_WEAK_POPOV,
                         const bool row_wise = true,
                         const bool generic = false
                        );

inline DegVec approximant_basis(
                                Mat<zz_pX> & appbas,
                                const Mat<zz_pX> & pmat,
                                const long order,
                                const Shift & shift = Shift(),
                                const PolMatForm form = ORD_WEAK_POPOV,
                                const bool row_wise = true,
                                const bool generic = false
                               )
{
    Order orders(pmat.NumCols(),order);
    return approximant_basis(appbas,pmat,orders,shift,form,row_wise,generic);
}


/*------------------------------------------------------------*/
/* Verifying that appbas is a shift-minimal approximant       */
/* basis for input matrix 'pmat' and order 'order'            */
/* 'form' gives the minimal requirement to check (matrix must */
/* be at least in the form 'form')                            */
/* 'randomized' says whether using a Monte Carlo or Las Vegas */
/* verification algorithm is acceptable                       */
/* Note: currently, deterministic verification is, for most   */
/* instances, as long as re-computing the basis               */
/*------------------------------------------------------------*/

bool is_approximant_basis(
                          const Mat<zz_pX> & appbas,
                          const Mat<zz_pX> & pmat,
                          const Order & order,
                          const Shift & shift = Shift(),
                          const PolMatForm & form = ORD_WEAK_POPOV,
                          const bool row_wise = true,
                          const bool randomized = false
                         );

inline bool is_approximant_basis(
                          const Mat<zz_pX> & appbas,
                          const Mat<zz_pX> & pmat,
                          const long order,
                          const Shift & shift = Shift(),
                          const PolMatForm & form = ORD_WEAK_POPOV,
                          const bool row_wise = true,
                          const bool randomized = false
                         )
{
    Order orders(pmat.NumCols(),order);
    return is_approximant_basis(appbas,pmat,orders,shift,form,row_wise,randomized);
}


/*------------------------------------------------------------*/
/* Iterative algorithm for general order and shift            */
/* References:                                                */
/*   - Beckermann 1992                                        */
/*   - Van Barel-Bultheel 1991+1992                           */
/*   - Beckermann-Labahn 2000 (ensuring s-Popov)              */
/*------------------------------------------------------------*/
DegVec appbas_iterative(
                        Mat<zz_pX> & appbas,
                        const Mat<zz_pX> & pmat,
                        const Order & order,
                        const Shift & shift,
                        bool order_wise=true
                       );

DegVec popov_appbas_iterative(
                              Mat<zz_pX> & appbas,
                              const Mat<zz_pX> & pmat,
                              const Order & order,
                              const Shift & shift,
                              bool order_wise=true
                             );

/*------------------------------------------------------------*/
/* M-Basis algorithm for approximant order = 1                */
/* References:                                                */
/*   - Giorgi-Jeannerod-Villard ISSAC 2003 (algo)             */
/*   - Giorgi-Lebreton ISSAC 2014 (algo with explicit shift)  */
/*   - Jeannerod-Neiger-Villard 2018 (ensuring s-Popov)       */
/*------------------------------------------------------------*/
// input: kerbas is constant, will contain the left kernel of pmat in reduced REF
// output: pivot degrees of the approx basis (also indicates where the rows of
// kernel should appear in approx basis)
DegVec popov_mbasis1(
                     Mat<zz_p> & kerbas,
                     const Mat<zz_p> & pmat,
                     const Shift & shift
                    );


// TODO check if serious difference of time if not returning Popov but just
// minimal, like done in LinBox and in GJV03 and GL14 (implies slightly less
// permutation work: the final permutation of the rows is not necessary)

/*------------------------------------------------------------*/
/* M-Basis algorithm for uniform approximant order            */
/* References:                                                */
/*   - Giorgi-Jeannerod-Villard ISSAC 2003 (algo)             */
/*   - Giorgi-Lebreton ISSAC 2014 (algo with explicit shift)  */
/*   - Jeannerod-Neiger-Villard 2018                          */
/*          (ensuring s-ordered weak Popov or s-Popov)        */
/*------------------------------------------------------------*/

// plain version, not the most efficient
DegVec mbasis_plain(
                    Mat<zz_pX> & appbas,
                    const Mat<zz_pX> & pmat,
                    const long order,
                    const Shift & shift
                   );

// variant which first converts to vector of constant matrices,
// performs the computations with this storage, and eventually
// converts back to polynomial matrices
// Residual (constant coeff of X^-d appbas*pmat) is computed from scratch at
// each iteration
// Complexity: pmat is m x n
//   - 'order' calls to popov_mbasis1 with dimension m x n, each one gives a
//   constant matrix K which is generically m-n x m  (may have more rows in
//   exceptional cases)
//   - order products (X Id + K ) * appbas to update the approximant basis
//   - order computations of "coeff k of appbas*pmat" to find residuals
// Assuming the degree of appbas at iteration 'ord' is m 'ord' / n (it is at
// least this almost always; and for the uniform shift it is equal to this for
// generic pmat), then the third item costs O(m n^2 order^2 / 2) operations,
// assuming cubic matrix multiplication over the field.
DegVec mbasis_rescomp(
              Mat<zz_pX> & appbas,
              const Mat<zz_pX> & pmat,
              const long order,
              const Shift & shift
             );

DegVec mbasis_rescomp_v2(
              Mat<zz_pX> & appbas,
              const Mat<zz_pX> & pmat,
              const long order,
              const Shift & shift
             );

// variant which first converts to vector of constant matrices,
// performs the computations with this storage, and eventually
// converts back to polynomial matrices
// Residual (X^-d appbas*pmat mod X^(order-d)) is continuously updated along
// the iterations
// Complexity: pmat is m x n
//   - 'order' calls to popov_mbasis1 with dimension m x n, each one gives a
//   constant matrix K which is generically m-n x m  (may have more rows in
//   exceptional cases)
//   - order products (X Id + K ) * appbas to update the approximant basis
//   - order products (X Id + K ) * pmat to update the residual
// Assuming cubic matrix multiplication over the field, the third item costs
// O(m n (m-n) order^2/2) operations
DegVec mbasis_resupdate(
                        Mat<zz_pX> & appbas,
                        const Mat<zz_pX> & pmat,
                        const long order,
                        const Shift & shift
                       );

// main function choosing the most efficient variant depending on parameters
// warning: may not be the best choice when the shift is not uniform
// FIXME -->  try to find the threshold for shifted case?
// FIXME -->  or simply assume the user will choose the right mbasis?
// FIXME -->  mbasis is anyway not the best approach, at least on the paper,
//            when cdim << rdim and shift is "highly" non-uniform
inline DegVec mbasis(
              Mat<zz_pX> & appbas,
              const Mat<zz_pX> & pmat,
              const long order,
              const Shift & shift
             )
{
    long rdim = pmat.NumRows();
    long cdim = pmat.NumCols();
    if (cdim > rdim/2 + 1)
        return mbasis_resupdate(appbas, pmat, order, shift);
    else
        return mbasis_rescomp(appbas, pmat, order, shift);
    // To understand the threshold (cdim > rdim/2 + 1), see the complexities
    // mentioned above for these two variants of mbasis
}


DegVec popov_mbasis(
                    Mat<zz_pX> &appbas,
                    const Mat<zz_pX> & pmat,
                    const long order,
                    const Shift & shift
                   );

/*------------------------------------------------------------*/
/* PM-Basis algorithm for uniform approximant order           */
/* References:                                                */
/*   - Giorgi-Jeannerod-Villard ISSAC 2003 (algo)             */
/*   - Giorgi-Lebreton ISSAC 2014 (algo with explicit shifts) */
/*   - Jeannerod-Neiger-Villard 2018                          */
/*          (ensuring s-ordered weak Popov or s-Popov)        */
/*------------------------------------------------------------*/
// TODO handle different orders via max(order)-order shifting of pmat
DegVec pmbasis(
               Mat<zz_pX> & appbas,
               const Mat<zz_pX> & pmat,
               const long order,
               const Shift & shift
              );

DegVec popov_pmbasis(
                     Mat<zz_pX> &appbas,
                     const Mat<zz_pX> & pmat,
                     const long order,
                     const Shift & shift
                    );

/*------------------------------------------------------------*/
/* FIXME in progress: MBASIS/PMBASIS, generic case            */
/*------------------------------------------------------------*/

// TODO allow non-uniform shifts?
// TODO several columns

DegVec popov_mbasis1_generic(
                     Mat<zz_p> & kerbas,
                     const Mat<zz_p> & pmat,
                     const Shift & shift
                    );

DegVec mbasis_generic(
                     Mat<zz_pX> & appbas,
                     const Mat<zz_pX> & pmat,
                     const long order,
                     const Shift & shift
                    );

DegVec pmbasis_generic(
               Mat<zz_pX> & appbas,
               const Mat<zz_pX> & pmat,
               const long order,
               const Shift & shift
              );


#endif // MAT_LZZ_PX_APPROXIMANT__H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
