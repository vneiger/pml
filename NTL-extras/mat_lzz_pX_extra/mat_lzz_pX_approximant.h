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
// (not exactly of pmat: permutations involved, depending on the shift)
// output: pivot degrees of the approximant basis (also indicates where the
// rows of kernel should appear in the approximant basis)
// FIXME: this is currently not optimized when matrices with non-generic
// rank profiles are given in input (pmat)
// --> for this, it would probably be better to rely on a library which
// provides precisely the decomposition we want (probably PLE, or PLUQ)
VecLong popov_mbasis1(
                     Mat<zz_p> & kerbas,
                     const Mat<zz_p> & pmat,
                     const VecLong & shift
                    );

// input: kerbas is constant, will contain the left kernel of pmat in row
// echelon form
// (not exactly of pmat: permutations involved, depending on the shift)
// output: pivot degrees of the approximant basis (also indicates where the
// rows of kernel should appear in the approximant basis)
// FIXME: this is currently not optimized when non-generic matrices
// are given in input (pmat)
// --> for this, it would probably be better to rely on a library which
// provides precisely the decomposition we want (probably PLE, or PLUQ)
VecLong mbasis1(
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
//   - 'order' calls to mbasis1 with dimension m x n, each one gives a
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
//   - 'order' calls to mbasis1 with dimension m x n, each one gives a
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
// Assumes that the entry has genericity properties: precisely, that the
// computed kernels (base cases at order 1) are of the form [ * | Id ]
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
/* Resupdate version, requiring m = 2 n and order even        */
/*------------------------------------------------------------*/
// Assumes that the entry has genericity properties: precisely, that the
// computed kernels (base cases at order 1) are of the form [ * | Id ]
// requirement 1: m = 2*n
// requirement 2: order is even
// output: appbas is in 0-Popov form with row degree (d,.., d) *GEN*,
// where d = order/2
void mbasis_generic_2n_n_resupdate(
                                   Mat<zz_pX> & appbas,
                                   const Mat<zz_pX> & pmat,
                                   const long order
                                  );

/*------------------------------------------------------------*/
/* Divide and Conquer: PMBasis                                */
/* Via mbasis-resupdate, requiring m = 2 n and order even     */
/* --all computations done with n x n submatrices             */
/*------------------------------------------------------------*/
// requirement 1: m = 2*n
// requirement 2: order is even and strictly positive (TODO remove)
// output: appbas is in 0-ordered weak Popov form with row degree (d,.., d) *GEN*,
// where d = order/2
// --> in fact a more precise form is obtained,
// appbas = [ [X^k I + P00,  P01], [X P10, X^k I + X P11]]
// where P00, P01, P10 have degree k-1 and P11 has degree k-2
// (in particular, its top-left and bottom-right blocks are 0-Popov)
void pmbasis_generic_2n_n(
                          Mat<zz_pX> & appbas,
                          const Mat<zz_pX> & pmat,
                          const long order
                         );

// TODO doc if this turns out useful
void pmbasis_generic_2n_n_top_rows(
                                   Mat<zz_pX> & appbas,
                                   const Mat<zz_pX> & pmat,
                                   const long order
                                  );










/*------------------------------------------------------------*/
/* FIXME in progress: MBASIS/PMBASIS, generic case, one column */
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
/* MATRIX-PADE APPROXIMATION -- GENERIC INPUT -- no shift     */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

// typical use: matrix fraction reconstruction in the context
// of block-Wiedemann-like algorithms

// Input: a square n x n matrix pmat of degree < order
// Output: a square n x n matrix den of degree <= order/2 such that
// den * pmat has degree < order/2

// Note: this can be found as the leading principal n x n submatrix of a
// 0-ordered weak Popov approximant basis for [[pmat], [-Id]] at order 'order'
// Here we aim at deriving algorithms close to mbasis_resupdate and pmbasis but
// which exploit the fact that there is this identity matrix in the input.

/*------------------------------------------------------------*/
/* Iterative algorithm, for low approximation order           */
/*------------------------------------------------------------*/

// This algorithm could also be called "Matrix Berlekamp-Massey, iterative".
// It is roughly the same as Berlekamp-Massey (on square matrices), except that
// it is "reversed": it goes from low degree to high degrees.

// Assumes that the input has genericity properties: precisely, that the
// computed kernels (base cases at order 1) are of the form [ * | Id ]

// computes both den1 and den2, such that [[den1], [den2]] is the first
// block-column of a 0-ordered weak Popov approximant basis for [[pmat], [-Id]]
// at order 'order'
// Note: den1 is in Popov form.
// requirement: order is even // TODO remove requirement?
void matrix_pade_generic_iterative(
                                   Mat<zz_pX> & den1,
                                   Mat<zz_pX> & den2,
                                   const Mat<zz_pX> & pmat,
                                   const long order
                                  );

// Note: den is in Popov form.
// Note: could be made slightly faster by taking the same algorithm as the one
// just above, but not computing den2 at all. Yet, our main focus is the divide
// and conquer version, in which this is used at the base case of the
// recursion, and for which we need den2.
inline void matrix_pade_generic_iterative(
                                          Mat<zz_pX> & den,
                                          const Mat<zz_pX> & pmat,
                                          const long order
                                         )
{
    Mat<zz_pX> den2;
    matrix_pade_generic_iterative(den, den2, pmat, order);
}

/*------------------------------------------------------------*/
/* Divide and conquer algorithm                               */
/*------------------------------------------------------------*/

// computes both den1 and den2, such that [[den1], [den2]] is the first
// block-column of a 0-ordered weak Popov approximant basis for [[pmat], [-Id]]
// at order 'order'

// Note: den is in Popov form.
void matrix_pade_generic(
                         Mat<zz_pX> & den,
                         const Mat<zz_pX> & pmat,
                         const long order
                        );

// version computing den as 2n x n, storing the two left blocks
// [[den1], [den2]] above (den1 in Popov form).
void matrix_pade_generic_recursion(
                                   Mat<zz_pX> & den,
                                   const Mat<zz_pX> & pmat,
                                   const long order
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
