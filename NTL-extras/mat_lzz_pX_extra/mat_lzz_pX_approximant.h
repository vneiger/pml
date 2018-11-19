#ifndef MAT_LZZ_PX_APPROXIMANT__H
#define MAT_LZZ_PX_APPROXIMANT__H

#include "mat_lzz_pX_extra.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* MINIMAL APPROXIMANT BASES                                  */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* Definitions (shifted minimal approximant basis)            */
/*------------------------------------------------------------*/

// Approximant basis: consider
//   * m x n matrix of univariate polynomials 'pmat',
//   * approximation order 'order' (list of n positive integers),
// then n approximant basis for (pmat,order) is a matrix over K[X] whose rows
// form a basis for the K[X]-module
// { p in K[X]^{1 x m}  |  the column j of p F is 0 modulo X^{order[j]} }
//
// VecLonged minimal approximant basis: consider furthermore
//   * a degree shift 'shifts' (list of m integers)
// then an approximant basis for (pmat,order) is said to be "a shift-minimal"
// (resp. "a shift-ordered weak Popov", resp. "the shift-Popov") approximant
// basis if it is in shift-reduced form (resp. in shift-ordered weak Popov
// form, resp. in shift-Popov form)
// Refer to mat_lzz_pX_extra.h for definitions of these forms.

/*------------------------------------------------------------*/
/* Conventions                                                */
/*------------------------------------------------------------*/

// Most functions below use the following parameters:
//   - Mat<zz_pX> appbas : the output matrix (size does not need to be initialized; cannot be pmat)
//   - const Mat<zz_pX> & pmat : the input matrix (no constraint)
//   - const VecLong & order : the input order (list of strictly positive integers, length must be pmat.NumCols())
//   - const VecLong & shift : the input shift (list of integers, length must be pmat.NumRows())



/*------------------------------------------------------------*/
/* general user-friendly interface                            */
/*------------------------------------------------------------*/

// TODO draft implementation
VecLong approximant_basis(
                         Mat<zz_pX> & appbas,
                         const Mat<zz_pX> & pmat,
                         const VecLong & order,
                         const VecLong & shift = VecLong(),
                         const PolMatForm form = ORD_WEAK_POPOV,
                         const bool row_wise = true,
                         const bool generic = false
                        );

inline VecLong approximant_basis(
                                Mat<zz_pX> & appbas,
                                const Mat<zz_pX> & pmat,
                                const long order,
                                const VecLong & shift = VecLong(),
                                const PolMatForm form = ORD_WEAK_POPOV,
                                const bool row_wise = true,
                                const bool generic = false
                               )
{
    VecLong orders(pmat.NumCols(),order);
    return approximant_basis(appbas,pmat,orders,shift,form,row_wise,generic);
}


/*------------------------------------------------------------*/
/* Verifying that appbas is a shift-minimal approximant       */
/* basis for input matrix 'pmat' and order 'order'            */
/* 'form' gives the minimal requirement to check (matrix must */
/* be at least in the form 'form')                            */
/* 'randomized' says whether using a Monte Carlo or Las Vegas */
/* verification algorithm is acceptable                       */
/*------------------------------------------------------------*/

bool is_approximant_basis(
                          const Mat<zz_pX> & appbas,
                          const Mat<zz_pX> & pmat,
                          const VecLong & order,
                          const VecLong & shift,
                          const PolMatForm & form = ORD_WEAK_POPOV,
                          const bool randomized = false
                         );

inline bool is_approximant_basis(
                          const Mat<zz_pX> & appbas,
                          const Mat<zz_pX> & pmat,
                          const long order,
                          const VecLong & shift,
                          const PolMatForm & form = ORD_WEAK_POPOV,
                          const bool randomized = false
                         )
{
    VecLong orders(pmat.NumCols(),order);
    return is_approximant_basis(appbas,pmat,orders,shift,form,randomized);
}


/*------------------------------------------------------------*/
/* Iterative algorithm for general order and shift            */
/* References:                                                */
/*   - Beckermann 1992                                        */
/*   - Van Barel-Bultheel 1991+1992                           */
/*   - Beckermann-Labahn 2000 (ensuring s-Popov)              */
/*------------------------------------------------------------*/
VecLong appbas_iterative(
                        Mat<zz_pX> & appbas,
                        const Mat<zz_pX> & pmat,
                        const VecLong & order,
                        const VecLong & shift,
                        bool order_wise=true
                       );

VecLong popov_appbas_iterative(
                              Mat<zz_pX> & appbas,
                              const Mat<zz_pX> & pmat,
                              const VecLong & order,
                              const VecLong & shift,
                              bool order_wise=true
                             );

/*------------------------------------------------------------*/
/* M-Basis algorithm for approximant order = 1                */
/* References:                                                */
/*   - Giorgi-Jeannerod-Villard ISSAC 2003 (algo)             */
/*   - Giorgi-Lebreton ISSAC 2014 (algo with explicit shift)  */
/*   - Jeannerod-Neiger-Villard 2018 (ensuring s-Popov)       */
/*------------------------------------------------------------*/
// input: kerbas is constant, will contain the left kernel of pmat in reduced
// row echelon form
// output: pivot degrees of the approximant basis (also indicates where the
// rows of kernel should appear in the approximant basis)
// Note: NTL does not guarantee that the pivots in kerbas are 1 !
// --> gives the s-Popov form up to making the diagonal entries monic
VecLong popov_mbasis1(
                     Mat<zz_p> & kerbas,
                     const Mat<zz_p> & pmat,
                     const VecLong & shift
                    );

/*------------------------------------------------------------*/
/* M-Basis algorithm for uniform approximant order            */
/* References:                                                */
/*   - Giorgi-Jeannerod-Villard ISSAC 2003 (algo)             */
/*   - Giorgi-Lebreton ISSAC 2014 (algo with explicit shift)  */
/*   - Jeannerod-Neiger-Villard 2018                          */
/*          (ensuring s-ordered weak Popov or s-Popov)        */
/*------------------------------------------------------------*/

// plain version, not the most efficient
VecLong mbasis_plain(
                    Mat<zz_pX> & appbas,
                    const Mat<zz_pX> & pmat,
                    const long order,
                    const VecLong & shift
                   );

// Variant which first converts to vector of constant matrices,
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
VecLong mbasis_rescomp(
              Mat<zz_pX> & appbas,
              const Mat<zz_pX> & pmat,
              const long order,
              const VecLong & shift
             );

// same as mbasis_rescomp, with some multi-threading inserted 
// TODO prototype for the moment: not properly tuned and tested
VecLong mbasis_rescomp_multithread(
                                   Mat<zz_pX> & appbas,
                                   const Mat<zz_pX> & pmat,
                                   const long order,
                                   const VecLong & shift
                                  );

// Variant which first converts to vector of constant matrices,
// performs the computations with this storage, and eventually
// converts back to polynomial matrices
// Residual (X^-d appbas*pmat mod X^(order-d)) is continuously updated along
// the iterations
// Complexity: pmat is m x n
//   - 'order' calls to popov_mbasis1 with dimension m x n, each one gives a
//   constant matrix K which is generically m-n x m  (may have more rows in
//   exceptional cases)
//   - order products (X Id + K ) * appbas to update the approximant basis
//   - order-1 products (X Id + K ) * (matrix of degree order-ord) to update
//   the residual, for ord=1...order-1
// Assuming cubic matrix multiplication over the field, the third item costs
// O(m n (m-n) order^2/2) operations
VecLong mbasis_resupdate(
                         Mat<zz_pX> & appbas,
                         const Mat<zz_pX> & pmat,
                         const long order,
                         const VecLong & shift
                        );

// same as mbasis_resupdate, with some multi-threading inserted 
// TODO not implemented
//VecLong mbasis_resupdate_multithread(
//                                        Mat<zz_pX> & appbas,
//                                        const Mat<zz_pX> & pmat,
//                                        const long order,
//                                        const VecLong & shift
//                                       );


// main function choosing the most efficient variant depending on parameters
// warning: may not be the best choice when the shift is not uniform
// FIXME -->  try to find the threshold for shifted case?
// FIXME -->  or simply assume the user will choose the right mbasis?
// FIXME -->  mbasis is anyway not the best approach, at least on the paper,
//            when cdim << rdim and shift is "highly" non-uniform
inline VecLong mbasis(
                      Mat<zz_pX> & appbas,
                      const Mat<zz_pX> & pmat,
                      const long order,
                      const VecLong & shift
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


VecLong popov_mbasis(
                    Mat<zz_pX> &appbas,
                    const Mat<zz_pX> & pmat,
                    const long order,
                    const VecLong & shift
                   );

/*------------------------------------------------------------*/
/* PM-Basis algorithm for uniform approximant order           */
/* References:                                                */
/*   - Giorgi-Jeannerod-Villard ISSAC 2003 (algo)             */
/*   - Giorgi-Lebreton ISSAC 2014 (algo with explicit shifts) */
/*   - Jeannerod-Neiger-Villard 2018                          */
/*          (ensuring s-ordered weak Popov or s-Popov)        */
/*------------------------------------------------------------*/
VecLong pmbasis(
               Mat<zz_pX> & appbas,
               const Mat<zz_pX> & pmat,
               const long order,
               const VecLong & shift
              );

VecLong popov_pmbasis(
                     Mat<zz_pX> &appbas,
                     const Mat<zz_pX> & pmat,
                     const long order,
                     const VecLong & shift
                    );








/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* MBASIS -- GENERIC INPUT -- UNIFORM SHIFT                   */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* Rescomp version, requiring m = 2 n and order even          */
/*------------------------------------------------------------*/
// TODO try resupdate (m=2n is borderline between the two)
// TODO compare with mbasis_generic for m = t n based on Krylov
// requirement 1: m = 2*n
// requirement 2: order is even
// output: appbas is in 0-Popov form with row degree (d,.., d) *GEN*,
// where d = order/2
void mbasis_generic_2n_n_rescomp(
                                 Mat<zz_pX> & appbas,
                                 const Mat<zz_pX> & pmat,
                                 const long order
                                );





/*------------------------------------------------------------*/
/* FIXME in progress: MBASIS/PMBASIS, generic case, one column            */
/*------------------------------------------------------------*/

VecLong mbasis_generic_onecolumn(
                     Mat<zz_pX> & appbas,
                     const Mat<zz_pX> & pmat,
                     const long order,
                     const VecLong & shift
                    );

VecLong pmbasis_generic_onecolumn(
               Mat<zz_pX> & appbas,
               const Mat<zz_pX> & pmat,
               const long order,
               const VecLong & shift
              );




/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* TODO                                                       */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

// ** pmbasis: handle different orders via max(order)-order shifting of pmat

// ** generic (uniform shift): just do everything without shifts as would
// happen generically with the uniform shift; might need assumptions on n
// divide m, etc

// ** shifted generic: from shift, deduce the pivot degree expected generically
// and use this as a shift instead of 'shift', then obtain directly Popov
// approx basis; skipping the first call to find the pivot degree (in addition,
// will this be more efficient because shift is nicer?)

// ** thresholds:
//   -- mbasis_resupdate / mbasis_rescomp depending on m/n and nthreads
//   -- in pmbasis: base case

// ** threads:
//   -- in mbasis_resupdate
//   -- in pmbasis?

// ** form of output? currently, it is at least ordered weak Popov
// --> check if this makes a difference of time when not returning Popov but
// just minimal, like done in LinBox and in GJV03 and GL14 (implies slightly
// less permutation work)

// ** return value is pivot degree.
//   -- Would shifted row degree be more appropriate?
//   -- Why return rather than input reference?
//   -- At least in lower level functions, why not modifying directly the input shift?

// iterative version / mbasis: compare different ways of obtaining Popov
// (recompute with new shift, or maintain normal form)


#endif // MAT_LZZ_PX_APPROXIMANT__H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
