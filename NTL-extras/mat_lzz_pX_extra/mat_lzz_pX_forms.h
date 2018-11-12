#ifndef MAT_LZZ_PX_FORMS__H
#define MAT_LZZ_PX_FORMS__H

#include <NTL/matrix.h>
#include <NTL/lzz_pX.h>
#include <vector>
#include <algorithm>

NTL_CLIENT


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* Tools for dealing with shifted reduced forms and shifted   */
/* normal forms of polynomial matrices                        */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

// TODO replace with NTL's vec??
typedef std::vector<long> DegVec;
typedef std::vector<long> Shift;

/*------------------------------------------------------------*/
/* Shifted reduced forms of polynomials matrices. Recall that */
/* Popov => ordered weak Popov => weak Popov => Reduced       */
/*------------------------------------------------------------*/

enum PolMatForm {
    NONE = 0,
    REDUCED = 1, 
    WEAK_POPOV = 2,
    ORD_WEAK_POPOV = 3,
    POPOV = 4,
};



// TODO shift reduction (given degdet D, cf ISSAC)
// --> would be useful, e.g. in naive kernel

// TODO type for index tuples?
// TODO type for pair (pivot index, pivot degree)?
// remove those types to be more explicit?

// TODO random matrix with given PolMatForm

/*------------------------------------------------------------*/
/* check that the shift has the right dimension;              */
/* return true iff a non-empty shift was given as input       */
/*------------------------------------------------------------*/
inline bool check_shift(
                        const std::vector<long> & shift,
                        const Mat<zz_pX> & pmat,
                        const bool row_wise = true
                       )
{
    if (!shift.empty()) 
    {
        if (row_wise && (long)shift.size() != pmat.NumCols())
            throw std::invalid_argument("==check_shift== Provided shift does not have the right dimension (working row-wise)");
        if (!row_wise && (long)shift.size() != pmat.NumRows())
            throw std::invalid_argument("==check_shift== Provided shift does not have the right dimension (working column-wise)");
        return true;
    }
    else
        return false;
}


/*------------------------------------------------------------*/
/* Amplitude of a shift: max(shift) - min(shift)              */
/*------------------------------------------------------------*/

inline void amplitude(long & amp, Shift shift)
{
    auto minmax = std::minmax_element(shift.begin(), shift.end());
    amp = *minmax.second - *minmax.first;
}

inline long amplitude(Shift shift)
{
    auto minmax = std::minmax_element(shift.begin(), shift.end());
    return *minmax.second - *minmax.first;
}

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* SHIFTED DEGREES AND SHIFTED PIVOTS                         */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* Degree matrix: matrix of the degree of each entry          */
/* Convention: deg(0) = -1, more generally the shifted degree */
/* of a zero entry is min(shift)-1                            */
/*------------------------------------------------------------*/

// TODO explain what shifted degree is

void degree_matrix(
                   Mat<long> &degmat,
                   const Mat<zz_pX> &pmat,
                   const Shift & shift = Shift(),
                   const bool row_wise = true
                  );

inline Mat<long> degree_matrix(
                               const Mat<zz_pX> &pmat,
                               const Shift & shift = Shift(),
                               const bool row_wise = true
                              )
{
    Mat<long> degmat;
    degree_matrix(degmat, pmat, shift, row_wise);
    return degmat;
}

/*------------------------------------------------------------*/
/* tuple (shifted deg row 1, shifted deg row 2, ...)          */
/* where shifted deg row k is the maximum of the shifted      */
/* degrees of the entries of row k                            */
/*------------------------------------------------------------*/
void row_degree(
                DegVec & rdeg,
                const Mat<zz_pX> &pmat,
                const Shift & shift = Shift()
               ); 

// TODO version which returns rdeg

/*------------------------------------------------------------*/
/* similar function for column degrees (see row_degree)       */
/*------------------------------------------------------------*/
void column_degree(
                   DegVec & cdeg,
                   const Mat<zz_pX> &pmat,
                   const Shift & shift = Shift()
                  ); 

// TODO version which returns rdeg

/*------------------------------------------------------------*/
/* similar function with row-wise option and returning degree */
/*------------------------------------------------------------*/
inline DegVec vector_degree(
                            const Mat<zz_pX> &pmat,
                            const Shift & shift = Shift(),
                            const bool row_wise = true
                           )
{
    DegVec degs;
    if (row_wise)
        row_degree(degs,pmat,shift);
    else
        column_degree(degs,pmat,shift);
    return degs;
}


/*------------------------------------------------------------*/
/* finds the pivot indices; returns the row/col degs          */
/*------------------------------------------------------------*/
void pivot_index(
                 std::vector<long> & pivind,
                 DegVec & pivdeg,
                 const Mat<zz_pX> & pmat,
                 const Shift & shift = Shift(),
                 const bool row_wise = true
                );


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* TESTING SHIFTED REDUCED AND NORMAL FORMS                   */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* leading matrix of pmat                                     */
/*------------------------------------------------------------*/
void leading_matrix(
                    Mat<zz_p> &lmat,
                    const Mat<zz_pX> &pmat,
                    const Shift & shift = Shift(),
                    const bool row_wise = true
                   );

inline Mat<zz_p> leading_matrix(
                                const Mat<zz_pX> &pmat,
                                const Shift & shift = Shift(),
                                const bool row_wise = true
                               )
{
    Mat<zz_p> lmat;
    leading_matrix(lmat, pmat, shift, row_wise);
    return lmat;
}

/*------------------------------------------------------------*/
/* returns true if pmat is reduced                            */
/*------------------------------------------------------------*/
bool is_reduced(
                const Mat<zz_pX> &pmat,
                const Shift & shift = Shift(),
                const bool row_wise = true
               );

/*------------------------------------------------------------*/
/* returns true if pmat is in weak Popov form                 */
/*   - Note that zero rows or columns are not allowed         */
/*   - If ordered==true, pivot index must be increasing       */
/*      (ordered weak Popov form)                             */
/*------------------------------------------------------------*/
bool is_weak_popov(
                   const Mat<zz_pX> &pmat,
                   const Shift & shift = Shift(),
                   const bool row_wise = true,
                   const bool ordered= false
                  );

/*------------------------------------------------------------*/
/* returns true if pmat is in Popov form                      */
/*   - Note that zero rows or columns are not allowed         */
/*   - If up_to_permutation==true, relaxes the check to allow */
/*     Popov up to row permutation (if row_wise)              */
/*------------------------------------------------------------*/
bool is_popov(
              const Mat<zz_pX> &pmat,
              const Shift & shift = Shift(),
              const bool row_wise = true,
              const bool up_to_permutation = false
             );

/*------------------------------------------------------------*/
/* Return the strongest form of pmat among those in the       */
/* enum 'PolMatForm'                                          */
/*------------------------------------------------------------*/
PolMatForm get_polmatform(
                          const Mat<zz_pX> & pmat,
                          const Shift & shift = Shift(),
                          const bool row_wise = true
                         );

/*------------------------------------------------------------*/
/* Check whether pmat is in the form indicated by 'form'      */
/*------------------------------------------------------------*/
bool is_polmatform(
                   const Mat<zz_pX> & pmat,
                   const PolMatForm form,
                   const Shift & shift = Shift(),
                   const bool row_wise = true
                  );


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



#endif /* ifndef MAT_LZZ_PX_FORMS */

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
