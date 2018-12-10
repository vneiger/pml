#ifndef MAT_LZZ_PX_FORMS__H
#define MAT_LZZ_PX_FORMS__H

/** Shifted reduced/normal forms of polynomial matrices.
 *
 * \file mat_lzz_pX_forms.h
 * \author Seung Gyu Hyun, Vincent Neiger
 * \version 0.1
 * \date 2018-12-07
 *
 * Basic functions to deal with shifted reduced and shifted
 * normal forms of polynomials matrices: test if a matrix is
 * in some given form, compute shifted row/column degrees and
 * shifted pivot degrees, compute shifted leading matrix.
 *
 * \todo definitions
 *
 * \todo random matrix with given PolMatForm
 *
 */

#include "mat_lzz_pX_utils.h"
#include <algorithm>

NTL_CLIENT;

/** Shifted reduced forms of polynomial matrices.
 *
 * Note that these numbers follow the "strength" of these forms: Popov implies
 * ordered weak Popov, which implies weak Popov, which implies reduced.
 */
enum PolMatForm {
    NONE = 0, /**< Arbitrary matrix, no specific form */
    REDUCED = 1, /**< Matrix in (shifted) reduced form */
    WEAK_POPOV = 2, /**< Matrix in (shifted) weak Popov form */
    ORD_WEAK_POPOV = 3, /**< Matrix in (shifted) ordered weak Popov form */
    POPOV = 4, /**< Matrix in (shifted) Popov form */
};


/** @name Helper functions for shifts
 *
 * \todo reduction of the shift entries for a given determinantal degree (would
 * be useful, for example in kernel via approximation)
 */
//@{

/** Computes the amplitude `amp` of a shift `shift`, that is, the difference
 * between its largest entry and its smallest entry
 */
inline void amplitude(long & amp, VecLong shift)
{
    auto minmax = std::minmax_element(shift.begin(), shift.end());
    amp = *minmax.second - *minmax.first;
}

/** Computes and returns the amplitude of a shift `shift`, that is, the
 * difference between its largest entry and its smallest entry
 */
inline long amplitude(VecLong shift)
{
    auto minmax = std::minmax_element(shift.begin(), shift.end());
    return *minmax.second - *minmax.first;
}

//@} // doxygen group: Helper functions for shifts

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* (SHIFTED) ROW/COLUMN DEGREE                                */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/** @name (Shifted) row and column degree
 * \anchor RowAndColumnDegrees
 *
 * The degree of a (row or column) vector is the maximum of the degrees of the
 * entries of this vector; by convention, it is -1 if the vector is zero.  The
 * _row degree_ of a matrix is the tuple of the degrees of its rows, and
 * similarly the _column degree_ of a matrix is the tuple of the degrees of its
 * columns.
 *
 * In this context, a shift is a tuple of signed integers. For a given `shift`
 * of length `m`, the <em>`shift`-degree</em> of a (row or column) vector
 * `pvec` of length `m` is the maximum of the degrees of the entries of that
 * vector to which we add the corresponding shift entry, that is, the maximum
 * of `deg(pvec[i]) + shift[i]` for `0 <= i < m`. By convention, this is
 * `min(shift)-1` if the vector is zero, thus ensuring that the `shift`-degree
 * of a nonzero vector is always strictly larger than the `shift`-degree of the
 * zero vector.
 *
 * Then, the <em>`shift`-row degree</em> of a polynomial matrix is the tuple of
 * `shift`-degrees of its rows, and similarly the <em>`shift`-column
 * degree</em> is the tuple of `shift`-degrees of its columns.
 *
 * In particular, for the zero shift `shift = [0,...,0]`, these notions of
 * shifted degree of a vector and of shifted row/column degree of a matrix
 * coincide with the usual non-shifted notions defined above.
 *  
 */
//@{

/** Computes the row degree `rdeg` of a polynomial matrix `pmat` (see
 * @ref RowAndColumnDegrees)
 */
void row_degree(VecLong & rdeg, const Mat<zz_pX> & pmat); 

/** Computes and returns the row degree of a polynomial matrix `pmat` (see
 * @ref RowAndColumnDegrees)
 */
inline VecLong row_degree(const Mat<zz_pX> & pmat)
{ VecLong rdeg; row_degree(rdeg, pmat); return rdeg; }

/** Computes the `shift`-row degree `rdeg` of a polynomial matrix `pmat` (see
 * @ref RowAndColumnDegrees)
 */
void row_degree(
                VecLong & rdeg,
                const Mat<zz_pX> & pmat,
                const VecLong & shift
               ); 

/** Computes and returns the `shift`-row degree of a polynomial matrix `pmat`
 * (see @ref RowAndColumnDegrees)
 */
inline VecLong row_degree(const Mat<zz_pX> & pmat, const VecLong & shift)
{ VecLong rdeg; row_degree(rdeg, pmat, shift); return rdeg; }

/** Computes the column degree `cdeg` of a polynomial matrix `pmat` (see
 * @ref RowAndColumnDegrees)
 */
void col_degree(VecLong & cdeg, const Mat<zz_pX> & pmat); 

/** Computes and returns the column degree of a polynomial matrix `pmat` (see
 * @ref RowAndColumnDegrees)
 */
inline VecLong col_degree(const Mat<zz_pX> & pmat)
{ VecLong cdeg; col_degree(cdeg, pmat); return cdeg; }

/** Computes the `shift`-column degree `cdeg` of a polynomial matrix `pmat`
 * (see @ref RowAndColumnDegrees)
 */
void col_degree(
                VecLong & cdeg,
                const Mat<zz_pX> & pmat,
                const VecLong & shift
               ); 

/** Computes and returns the `shift`-column degree of a polynomial matrix
 * `pmat` (see @ref RowAndColumnDegrees)
 */
inline VecLong col_degree(const Mat<zz_pX> & pmat, const VecLong & shift)
{ VecLong cdeg; col_degree(cdeg, pmat, shift); return cdeg; }

//@} // doxygen group: (Shifted) row and column degree

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* (SHIFTED) PIVOT INDEX/DEGREE                               */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/** @name (Shifted) pivot index and pivot degree
 * \anchor Pivots
 *
 * row-wise (shifted) pivot indices and degrees
 * --> tuple of pivot index/degree for each row of pmat
 * For a given row, the pivot index is the index of the
 * rightmost entry which reaches the (shifted) row degree;
 * the pivot degree is the degree of that pivot (not adding
 * the shift)
 * The pivot index of a zero row is -1
 *  
 * column-wise (shifted) pivot indices and degrees
 * --> tuple of pivot index/degree for each column of pmat
 * For a given column, the pivot index is the index of the
 * bottommost entry which reaches the (shifted) column
 * the pivot degree is the degree of that pivot (not adding
 * the shift)
 * The pivot index of a zero column is -1
 */
//@{

/** Computes the row-wise pivot index `pivind` and pivot degree `pivdeg` of a
 * polynomial matrix `pmat` (see @ref Pivots)
 */
void row_pivots(VecLong & pivind, VecLong & pivdeg, const Mat<zz_pX> & pmat);

/** Computes and return a pair whose first element is the row-wise pivot index
 * of a polynomial matrix `pmat` and whose second element is the row-wise pivot
 * degree of `pmat` (see @ref Pivots)
 */
inline std::pair<VecLong,VecLong> row_pivots(const Mat<zz_pX> & pmat)
{
    std::pair<VecLong,VecLong> pivots;
    row_pivots(pivots.first, pivots.second, pmat);
    return pivots;
}

/** Computes the row-wise `shift`-pivot index `pivind` and `shift`-pivot degree
 * `pivdeg` of a polynomial matrix `pmat` (see @ref Pivots)
 */
void row_pivots(
                VecLong & pivind,
                VecLong & pivdeg,
                const Mat<zz_pX> & pmat,
                const VecLong & shift
               );

/** Computes and return a pair whose first element is the row-wise
 * `shift`-pivot index of a polynomial matrix `pmat` and whose second element
 * is the row-wise `shift`-pivot degree of `pmat` (see @ref Pivots)
 */
inline std::pair<VecLong,VecLong>
row_pivots(const Mat<zz_pX> & pmat, const VecLong & shift)
{
    std::pair<VecLong,VecLong> pivots;
    row_pivots(pivots.first, pivots.second, pmat, shift);
    return pivots;
}

/** Computes the column-wise pivot index `pivind` and pivot degree `pivdeg` of
 * a polynomial matrix `pmat` (see @ref Pivots)
 */
void col_pivots(VecLong & pivind, VecLong & pivdeg, const Mat<zz_pX> & pmat);

/** Computes and return a pair whose first element is the column-wise pivot
 * index of a polynomial matrix `pmat` and whose second element is the
 * column-wise pivot degree of `pmat` (see @ref Pivots)
 */
inline std::pair<VecLong,VecLong> col_pivots(const Mat<zz_pX> & pmat)
{
    std::pair<VecLong,VecLong> pivots;
    col_pivots(pivots.first, pivots.second, pmat);
    return pivots;
}

/** Computes the column-wise `shift`-pivot index `pivind` and `shift`-pivot
 * degree `pivdeg` of a polynomial matrix `pmat` (see @ref Pivots)
 */
void col_pivots(
                VecLong & pivind,
                VecLong & pivdeg,
                const Mat<zz_pX> & pmat,
                const VecLong & shift
               );

/** Computes and return a pair whose first element is the column-wise
 * `shift`-pivot index of a polynomial matrix `pmat` and whose second element
 * is the column-wise `shift`-pivot degree of `pmat` (see @ref Pivots)
 */
inline std::pair<VecLong,VecLong>
col_pivots(const Mat<zz_pX> & pmat, const VecLong & shift)
{
    std::pair<VecLong,VecLong> pivots;
    col_pivots(pivots.first, pivots.second, pmat, shift);
    return pivots;
}

//@} // doxygen group: (Shifted) pivot index and pivot degree


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* (SHIFTED) DEGREE MATRIX                                    */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* Degree matrix: matrix of the degree of each entry          */
/* Convention: deg(0) = -1, more generally the shifted degree */
/* of a zero entry is min(shift)-1                            */
/*------------------------------------------------------------*/

void degree_matrix(
                   Mat<long> & degmat,
                   const Mat<zz_pX> & pmat
                  );

inline Mat<long> degree_matrix(const Mat<zz_pX> &pmat)
{ Mat<long> degmat; degree_matrix(degmat, pmat); return degmat; }

/*------------------------------------------------------------*/
/* Row-shifted degree matrix: matrix whose entry (i,j) is the */
/* degree of pmat[i][j], shifted by adding shift[j]           */
/* Convention: for a zero entry, its degree is min(shift)-1   */
/*------------------------------------------------------------*/

void degree_matrix_rowshifted(
                              Mat<long> &degmat,
                              const Mat<zz_pX> &pmat,
                              const VecLong & shift
                             );

inline Mat<long> degree_matrix_rowshifted(
                                          const Mat<zz_pX> &pmat,
                                          const VecLong & shift
                                         )
{
    Mat<long> degmat;
    degree_matrix_rowshifted(degmat, pmat, shift);
    return degmat;
}

/*------------------------------------------------------------*/
/* Column-shifted degree matrix: matrix whose entry (i,j) is  */
/* the degree of pmat[i][j], shifted by adding shift[i]       */
/* Convention: for a zero entry, its degree is min(shift)-1   */
/*------------------------------------------------------------*/

void degree_matrix_colshifted(
                              Mat<long> &degmat,
                              const Mat<zz_pX> &pmat,
                              const VecLong & shift
                             );

inline Mat<long> degree_matrix_colshifted(
                                          const Mat<zz_pX> &pmat,
                                          const VecLong & shift
                                         )
{
    Mat<long> degmat;
    degree_matrix_colshifted(degmat, pmat, shift);
    return degmat;
}


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* (SHIFTED) LEADING MATRIX                                   */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* row-wise leading matrix of pmat                            */
/* this is the constant matrix whose entry (i,j) is the       */
/* coefficient of degree rdeg[i] of the entry (i,j) of pmat,  */
/* for rdeg the row degree of pmat                            */
/* (this is zero if pmat[i][j] does not reach rdeg[i])        */
/*------------------------------------------------------------*/
void row_leading_matrix(
                        Mat<zz_p> & lmat,
                        const Mat<zz_pX> & pmat
                       );

inline Mat<zz_p> row_leading_matrix(const Mat<zz_pX> & pmat)
{ Mat<zz_p> lmat; row_leading_matrix(lmat, pmat); return lmat; }

// TODO row_leading_matrix_known_rdeg

/*------------------------------------------------------------*/
/* row-wise shifted leading matrix of pmat                    */
/* this is the constant matrix whose entry (i,j) is the       */
/* coefficient of degree rdeg[i]-shift[j] of the entry (i,j)  */
/* of pmat, for rdeg the shifted row degree of pmat           */
/* (this is zero if pmat[i][j] does not reach rdeg[i])        */
/*------------------------------------------------------------*/

void row_leading_matrix(
                        Mat<zz_p> &lmat,
                        const Mat<zz_pX> &pmat,
                        const VecLong & shift
                       );

inline Mat<zz_p> row_leading_matrix(
                                    const Mat<zz_pX> & pmat,
                                    const VecLong & shift
                                   )
{ Mat<zz_p> lmat; row_leading_matrix(lmat, pmat, shift); return lmat; }

// TODO row_leading_matrix_known_cdeg


/*------------------------------------------------------------*/
/* column-wise leading matrix of pmat                         */
/* this is the constant matrix whose entry (i,j) is the       */
/* coefficient of degree cdeg[j] of the entry (i,j) of pmat,  */
/* for cdeg the column degree of pmat                         */
/* (this is zero if pmat[i][j] does not reach cdeg[j])        */
/*------------------------------------------------------------*/
void col_leading_matrix(
                        Mat<zz_p> & lmat,
                        const Mat<zz_pX> & pmat
                       );

inline Mat<zz_p> col_leading_matrix(const Mat<zz_pX> & pmat)
{ Mat<zz_p> lmat; col_leading_matrix(lmat, pmat); return lmat; }

/*------------------------------------------------------------*/
/* column-wise shifted leading matrix of pmat                 */
/* this is the constant matrix whose entry (i,j) is the       */
/* coefficient of degree cdeg[j]-shift[i] of the entry (i,j)  */
/* of pmat, for cdeg the shifted column degree of pmat        */
/* (this is zero if pmat[i][j] does not reach rdeg[j])        */
/*------------------------------------------------------------*/

void col_leading_matrix(
                        Mat<zz_p> &lmat,
                        const Mat<zz_pX> &pmat,
                        const VecLong & shift
                       );

inline Mat<zz_p> col_leading_matrix(
                                    const Mat<zz_pX> & pmat,
                                    const VecLong & shift
                                   )
{ Mat<zz_p> lmat; col_leading_matrix(lmat, pmat, shift); return lmat; }



/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* TESTING SHIFTED REDUCED FORMS                              */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* test (shifted) row reducedness:                            */
/* a matrix is (shifted) row reduced if and only if its       */
/* (shifted) row-wise leading matrix has full row rank        */
/* --> note that zero rows are not allowed                    */
/*------------------------------------------------------------*/
bool is_row_reduced(const Mat<zz_pX> & pmat);
bool is_row_reduced(const Mat<zz_pX> & pmat, const VecLong & shift);

/*------------------------------------------------------------*/
/* test (shifted) column reducedness:                         */
/* a matrix is (shifted) column reduced if and only if its    */
/* (shifted) column-wise leading matrix has full column rank  */
/* --> note that zero columns are not allowed                 */
/*------------------------------------------------------------*/
bool is_col_reduced(const Mat<zz_pX> & pmat);
bool is_col_reduced(const Mat<zz_pX> & pmat, const VecLong & shift);


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* TESTING SHIFTED WEAK POPOV FORMS                           */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* test row-wise (shifted) weak Popov form                    */
/* a matrix is in row-wise (shifted) weak Popov form if and   */
/* only if its row-wise (shifted) pivot index consists of     */
/* pairwise distinct entries (i.e. no repetition)             */
/* --> note that zero rows are not allowed                    */
/*------------------------------------------------------------*/
bool is_row_weak_popov(const Mat<zz_pX> & pmat);
bool is_row_weak_popov(
                       const Mat<zz_pX> &pmat,
                       const VecLong & shift
                      );

/*------------------------------------------------------------*/
/* test column-wise (shifted) weak Popov form                 */
/* a matrix is in column-wise (shifted) weak Popov form if    */
/* and only if its column-wise (shifted) pivot index consists */
/* of pairwise distinct entries (i.e. no repetition)          */
/* --> note that zero columns are not allowed                 */
/*------------------------------------------------------------*/
bool is_col_weak_popov(const Mat<zz_pX> & pmat);
bool is_col_weak_popov(
                       const Mat<zz_pX> &pmat,
                       const VecLong & shift
                      );

/*------------------------------------------------------------*/
/* test row-wise (shifted) ordered weak Popov form            */
/* a matrix is in row-wise (shifted) ordered weak Popov form  */
/* if and only if its row-wise (shifted) pivot index is       */
/* strictly increasing                                        */
/* --> note that zero rows are not allowed                    */
/*------------------------------------------------------------*/
bool is_row_ordered_weak_popov(const Mat<zz_pX> & pmat);
bool is_row_ordered_weak_popov(
                               const Mat<zz_pX> &pmat,
                               const VecLong & shift
                              );

/*------------------------------------------------------------*/
/* test column-wise (shifted) ordered weak Popov form         */
/* a matrix is in column-wise (shifted) ordered weak Popov    */
/* form if and only if its column-wise (shifted) pivot index  */
/* is strictly increasing                                     */
/* --> note that zero columns are not allowed                 */
/*------------------------------------------------------------*/
bool is_col_ordered_weak_popov(const Mat<zz_pX> & pmat);
bool is_col_ordered_weak_popov(
                               const Mat<zz_pX> &pmat,
                               const VecLong & shift
                              );


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* TESTING SHIFTED POPOV FORMS                                */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* test row-wise (shifted) Popov form                         */
/* a matrix is in row-wise (shifted) Popov form iff:          */
/*   - it is in row-wise (shifted) ordered weak Popov form    */
/*   - its pivot entries are monic                            */
/*   - its entries below and above a pivot entry have degree  */
/*   less than this pivot entry                               */
/* --> note that zero rows are not allowed                    */
/*------------------------------------------------------------*/
bool is_row_popov(const Mat<zz_pX> & pmat);
bool is_row_popov(const Mat<zz_pX> & pmat, const VecLong & shift);

/*------------------------------------------------------------*/
/* same test but relaxing ordered weak Popov --> weak Popov   */
/* (allows definition of Popov with different orderings of    */
/* the rows, such as by increasing degree)                    */
/*------------------------------------------------------------*/
// TODO 
//bool is_row_popov_up_to_permutation(const Mat<zz_pX> & pmat);
//bool is_row_popov_up_to_permutation(
//                                    const Mat<zz_pX> & pmat,
//                                    const VecLong & shift
//                                    );


/*------------------------------------------------------------*/
/* test column-wise (shifted) Popov form                      */
/* a matrix is in column-wise (shifted) Popov form iff:       */
/*   - it is in column-wise (shifted) ordered weak Popov form */
/*   - its pivot entries are monic                            */
/*   - its entries to the left and to the right of a pivot    */
/*    entry have degree less than this pivot entry            */
/* --> note that zero columns are not allowed                 */
/*------------------------------------------------------------*/
bool is_col_popov(const Mat<zz_pX> & pmat);
bool is_col_popov(const Mat<zz_pX> & pmat, const VecLong & shift);

/*------------------------------------------------------------*/
/* same test but relaxing ordered weak Popov --> weak Popov   */
/* (allows definition of Popov with different orderings of    */
/* the columns, such as by increasing degree)                 */
/*------------------------------------------------------------*/
// TODO 
//bool is_col_popov_up_to_permutation(const Mat<zz_pX> & pmat);
//bool is_col_popov_up_to_permutation(
//                                    const Mat<zz_pX> & pmat,
//                                    const VecLong & shift
//                                    );


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* TESTING SHIFTED FORMS (FORM SELECTOR)                      */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* Check whether pmat is in the prescribed row-wise form      */
/*------------------------------------------------------------*/
bool is_row_polmatform(
                       const Mat<zz_pX> & pmat,
                       const PolMatForm form
                      );

bool is_row_polmatform(
                       const Mat<zz_pX> & pmat,
                       const VecLong & shift,
                       const PolMatForm form
                      );

/*------------------------------------------------------------*/
/* Check whether pmat is in the prescribed column-wise form   */
/*------------------------------------------------------------*/
bool is_col_polmatform(
                       const Mat<zz_pX> & pmat,
                       const PolMatForm form
                      );

bool is_col_polmatform(
                       const Mat<zz_pX> & pmat,
                       const VecLong & shift,
                       const PolMatForm form
                      );




/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* FIND STRONGEST (SHIFTED) FORM                              */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* Return the strongest row-wise form of pmat among:          */
/* Popov => ordered weak Popov => weak Popov                  */
/*                                       => reduced => none   */
/*------------------------------------------------------------*/
PolMatForm get_row_polmatform(const Mat<zz_pX> & pmat);
PolMatForm get_row_polmatform(
                              const Mat<zz_pX> & pmat,
                              const VecLong & shift
                             );

/*------------------------------------------------------------*/
/* Return the strongest column-wise form of pmat among:       */
/* Popov => ordered weak Popov => weak Popov                  */
/*                                       => reduced => none   */
/*------------------------------------------------------------*/
PolMatForm get_col_polmatform(const Mat<zz_pX> & pmat);
PolMatForm get_col_polmatform(
                              const Mat<zz_pX> & pmat,
                              const VecLong & shift
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
