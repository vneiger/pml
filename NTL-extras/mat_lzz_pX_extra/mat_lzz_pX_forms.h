#ifndef MAT_LZZ_PX_FORMS__H
#define MAT_LZZ_PX_FORMS__H

/** \brief Shifted reduced forms and shifted normal forms of univariate
 * polynomial matrices over `zz_p`
 *
 * \file mat_lzz_pX_forms.h
 * \author Seung Gyu Hyun, Vincent Neiger, Eric Schost
 * \version 0.1
 * \date 2018-12-07
 *
 * Basic functions to deal with shifted reduced and shifted
 * normal forms of polynomials matrices: test if a matrix is
 * in some given form, compute shifted row/column degrees and
 * shifted pivot degrees, compute shifted leading matrix.
 *
 * \todo random matrix with given PolMatForm
 *
 */

#include "mat_lzz_pX_utils.h"
#include <algorithm>

NTL_CLIENT

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
 * The functions below which involve a `shift` among its parameters throw
 * an error if this `shift` does not have the right length.
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
 * For a given shift and a given (row or column) polynomial vector, its _shifted
 * pivot index_ is the largest index corresponding to an entry which reaches the
 * shifted degree of the vector, and its _shifted pivot degree_ is the degree of
 * that entry (without adding the shift entry). By convention, both the shifted
 * pivot index and the shifted pivot degree of a zero vector are -1.
 *
 * Then, the row-wise shifted pivot index (resp. degree) of a polynomial matrix
 * is the tuple of the shifted pivot indices (resp. degrees) of the rows of
 * this matrix. Similarly, the column-wise shifted pivot index (resp. degree)
 * of a polynomial matrix is the tuple of the shifted pivot indices (resp.
 * degrees) of the columns of this matrix.
 *
 * The functions below which involve a `shift` among its parameters throw
 * an error if this `shift` does not have the right length.
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

/** @name (Shifted) degree matrix
 * \anchor DegreeMatrix
 *
 * The _degree matrix_ of an `m x n` polynomial matrix `pmat` is the integer
 * matrix of dimensions `m x n` whose entry `(i,j)` is the degree of the entry
 * `(i,j)` of `pmat`. We recall that by convention, the zero polynomial has
 * degree -1.
 *
 * For a given shift `shift` of length `n`, the row-wise `shift`-degree matrix
 * of `pmat` is the `m x n` integer matrix whose entry `(i,j)` is
 * `deg(pmat[i][j]) + shift[j]` if `pmat[i][j]` is nonzero, and `min(shift)-1`
 * otherwise. For a given shift `shift` of length `m`, the column-wise
 * `shift`-degree matrix of `pmat` is the `m x n` integer matrix whose entry
 * `(i,j)` is `deg(pmat[i][j]) + shift[i]` if `pmat[i][j]` is nonzero, and
 * `min(shift)-1` otherwise.
 *
 * The functions below which involve a `shift` among its parameters throw
 * an error if this `shift` does not have the right length.
 */
//@{

/** Computes the degree matrix `degmat` of a polynomial matrix `pmat` (see @ref
 * DegreeMatrix)
 */
void degree_matrix(Mat<long> & degmat, const Mat<zz_pX> & pmat);

/** Computes and returns the degree matrix of a polynomial matrix `pmat` (see
 * @ref DegreeMatrix)
 */
inline Mat<long> degree_matrix(const Mat<zz_pX> &pmat)
{ Mat<long> degmat; degree_matrix(degmat, pmat); return degmat; }

/** Computes the row-wise `shift`-degree matrix `degmat` of a polynomial matrix
 * `pmat` (see @ref DegreeMatrix)
 */
void degree_matrix_rowshifted(
                              Mat<long> &degmat,
                              const Mat<zz_pX> &pmat,
                              const VecLong & shift
                             );

/** Computes and returns row-wise the `shift`-degree matrix of a polynomial
 * matrix `pmat` (see @ref DegreeMatrix)
 */
inline Mat<long> degree_matrix_rowshifted(
                                          const Mat<zz_pX> &pmat,
                                          const VecLong & shift
                                         )
{
    Mat<long> degmat;
    degree_matrix_rowshifted(degmat, pmat, shift);
    return degmat;
}

/** Computes the column-wise `shift`-degree matrix `degmat` of a polynomial
 * matrix `pmat` (see @ref DegreeMatrix)
 */
void degree_matrix_colshifted(
                              Mat<long> &degmat,
                              const Mat<zz_pX> &pmat,
                              const VecLong & shift
                             );

/** Computes and returns column-wise the `shift`-degree matrix of a polynomial
 * matrix `pmat` (see @ref DegreeMatrix)
 */
inline Mat<long> degree_matrix_colshifted(
                                          const Mat<zz_pX> &pmat,
                                          const VecLong & shift
                                         )
{
    Mat<long> degmat;
    degree_matrix_colshifted(degmat, pmat, shift);
    return degmat;
}

//@} // doxygen group: (Shifted) degree matrix

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* (SHIFTED) LEADING MATRIX                                   */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/** @name (Shifted) leading matrix
 * \anchor LeadingMatrix
 *
 * Consider a polynomial matrix `pmat` of dimension `m x n`, with row degree
 * `rdeg` (see @ref RowAndColumnDegrees). Then, the row-wise leading matrix of
 * `pmat` is the `m x n` matrix over the base field whose entry `(i,j)` is the
 * coefficient of degree `rdeg[i]` of the entry `pmat[i][j]` (this is zero if
 * `pmat[i][j]` does not reach `rdeg[i]`). Similarly, writing `cdeg` for the
 * column degree of `pmat`, the column-wise leading matrix of `pmat` is the `m
 * x n` matrix over the base field whose entry `(i,j)` is the coefficient of
 * degree `cdeg[j]` of the entry `pmat[i][j]` (this is zero if `pmat[i][j]`
 * does not reach `cdeg[j]`).
 *  
 * More generally, given a shift `shift` of length `n`, the row-wise
 * `shift`-leading matrix of `pmat` is the `m x n` matrix over the base field
 * whose entry `(i,j)` is the coefficient of degree `rdeg[i]-shift[j]` of the
 * entry `pmat[i][j]`, where `rdeg` is now the `shift`-row degree of `pmat`.
 * For a shift `shift` of length `m`, the column-wise `shift`-leading matrix of
 * `pmat` is the `m x n` matrix over the base field whose entry `(i,j)` is the
 * coefficient of degree `cdeg[j]-shift[i]` of the entry `pmat[i][j]`, where
 * `cdeg` is now the `shift`-column degree of `pmat`.
 *
 * The functions below which involve a `shift` among its parameters throw
 * an error if this `shift` does not have the right length.
 *
 * \todo minor improvement: offer row-wise (resp column-wise) leading matrix
 * when the row degree (resp column degree) is already known?
 */
//@{

/** Computes the row-wise leading matrix `lmat` of a polynomial matrix `pmat`
 * (see @ref LeadingMatrix)
 */
void row_leading_matrix(Mat<zz_p> & lmat, const Mat<zz_pX> & pmat);

/** Computes and returns the row-wise leading matrix of a polynomial matrix
 * `pmat` (see @ref LeadingMatrix)
 */
inline Mat<zz_p> row_leading_matrix(const Mat<zz_pX> & pmat)
{ Mat<zz_p> lmat; row_leading_matrix(lmat, pmat); return lmat; }

/** Computes the row-wise `shift`-leading matrix `lmat` of a polynomial matrix
 * `pmat` (see @ref LeadingMatrix)
 */
void row_leading_matrix(
                        Mat<zz_p> &lmat,
                        const Mat<zz_pX> &pmat,
                        const VecLong & shift
                       );

/** Computes and returns the row-wise `shift`-leading matrix of a polynomial
 * matrix `pmat` (see @ref LeadingMatrix)
 */
inline Mat<zz_p> row_leading_matrix(
                                    const Mat<zz_pX> & pmat,
                                    const VecLong & shift
                                   )
{ Mat<zz_p> lmat; row_leading_matrix(lmat, pmat, shift); return lmat; }

/** Computes the column-wise leading matrix `lmat` of a polynomial matrix
 * `pmat` (see @ref LeadingMatrix)
 */
void col_leading_matrix(
                        Mat<zz_p> & lmat,
                        const Mat<zz_pX> & pmat
                       );

/** Computes and returns the column-wise leading matrix of a polynomial matrix
 * `pmat` (see @ref LeadingMatrix)
 */
inline Mat<zz_p> col_leading_matrix(const Mat<zz_pX> & pmat)
{ Mat<zz_p> lmat; col_leading_matrix(lmat, pmat); return lmat; }

/** Computes the column-wise `shift`-leading matrix `lmat` of a polynomial
 * matrix `pmat` (see @ref LeadingMatrix)
 */
void col_leading_matrix(
                        Mat<zz_p> &lmat,
                        const Mat<zz_pX> &pmat,
                        const VecLong & shift
                       );

/** Computes and returns the column-wise `shift`-leading matrix of a polynomial
 * matrix `pmat` (see @ref LeadingMatrix)
 */
inline Mat<zz_p> col_leading_matrix(
                                    const Mat<zz_pX> & pmat,
                                    const VecLong & shift
                                   )
{ Mat<zz_p> lmat; col_leading_matrix(lmat, pmat, shift); return lmat; }

//@} // doxygen group: (Shifted) leading matrix


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* TESTING MATRIX FORMS                                       */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/** @name Testing polynomial matrix forms
 * \anchor MatrixForms
 *
 * A polynomial matrix is said to be in _(shifted) row reduced form_ if its
 * row-wise (shifted) leading matrix has full row rank; note that zero rows are
 * not allowed in this definition. The definition of (shifted) column reduced
 * forms is analogous.
 *
 * A polynomial matrix is said to be in row-wise _(shifted) weak Popov form_ if
 * it has no zero rows and there is no repetition in its row-wise (shifted)
 * pivot index (that is, its entries are pairwise distinct). Such a form is
 * furthermore said to be _ordered_ if this (shifted) pivot index is strictly
 * increasing. The definition of column-wise (shifted, ordered) weak Popov
 * forms is analogous.
 *
 * A polynomial matrix is said to be in row-wise (shifted) Popov form if:
 *   - it is in row-wise (shifted) ordered weak Popov form
 *   - its (shifted) pivot entries are monic
 *   - its entries below and above a (shifted) pivot entry have degree less
 *   than this pivot entry
 * In particular, such a matrix cannot have a zero row.
 *
 * A polynomial matrix is said to be in column-wise (shifted) Popov form if:
 *   - it is in column-wise (shifted) ordered weak Popov form
 *   - its (shifted) pivot entries are monic
 *   - its entries to the left and to the right of a (shifted) pivot entry have
 *   degree less than this pivot entry
 * In particular, such a matrix cannot have a zero column.
 *
 * The functions below which involve a `shift` among its parameters throw
 * an error if this `shift` does not have the right length.
 *  
 */
//@{

/** Computes and returns a boolean indicating whether `pmat` is in row reduced
 * form (see @ref MatrixForms)
 */
bool is_row_reduced(const Mat<zz_pX> & pmat);

/** Computes and returns a boolean indicating whether `pmat` is in `shift`-row
 * reduced form (see @ref MatrixForms)
 */
bool is_row_reduced(const Mat<zz_pX> & pmat, const VecLong & shift);

/** Computes and returns a boolean indicating whether `pmat` is in column
 * reduced form (see @ref MatrixForms)
 */
bool is_col_reduced(const Mat<zz_pX> & pmat);

/** Computes and returns a boolean indicating whether `pmat` is in
 * `shift`-column reduced form (see @ref MatrixForms)
 */
bool is_col_reduced(const Mat<zz_pX> & pmat, const VecLong & shift);


/** Computes and returns a boolean indicating whether `pmat` is in row-wise
 * weak Popov form (see @ref MatrixForms)
 */
bool is_row_weak_popov(const Mat<zz_pX> & pmat);

/** Computes and returns a boolean indicating whether `pmat` is in row-wise
 * `shift`-weak Popov form (see @ref MatrixForms)
 */
bool is_row_weak_popov(const Mat<zz_pX> & pmat, const VecLong & shift);

/** Computes and returns a boolean indicating whether `pmat` is in column-wise
 * weak Popov form (see @ref MatrixForms)
 */
bool is_col_weak_popov(const Mat<zz_pX> & pmat);

/** Computes and returns a boolean indicating whether `pmat` is in column-wise
 * `shift`-weak Popov form (see @ref MatrixForms)
 */
bool is_col_weak_popov(const Mat<zz_pX> &pmat, const VecLong & shift);

/** Computes and returns a boolean indicating whether `pmat` is in row-wise
 * ordered weak Popov form (see @ref MatrixForms)
 */
bool is_row_ordered_weak_popov(const Mat<zz_pX> & pmat);

/** Computes and returns a boolean indicating whether `pmat` is in row-wise
 * `shift`-ordered weak Popov form (see @ref MatrixForms)
 */
bool is_row_ordered_weak_popov(
                               const Mat<zz_pX> &pmat,
                               const VecLong & shift
                              );

/** Computes and returns a boolean indicating whether `pmat` is in column-wise
 * ordered weak Popov form (see @ref MatrixForms)
 */
bool is_col_ordered_weak_popov(const Mat<zz_pX> & pmat);

/** Computes and returns a boolean indicating whether `pmat` is in column-wise
 * `shift`-ordered weak Popov form (see @ref MatrixForms)
 */
bool is_col_ordered_weak_popov(
                               const Mat<zz_pX> &pmat,
                               const VecLong & shift
                              );

/** Computes and returns a boolean indicating whether `pmat` is in row-wise
 * Popov form (see @ref MatrixForms)
 *
 * \todo provide same test but relaxing ordered weak Popov to weak Popov in the
 * definition? This allows easier support for definitions of Popov with
 * different orderings of the rows, such as by increasing degree
 */
bool is_row_popov(const Mat<zz_pX> & pmat);

/** Computes and returns a boolean indicating whether `pmat` is in row-wise
 * `shift`-Popov form (see @ref MatrixForms)
 */
bool is_row_popov(const Mat<zz_pX> & pmat, const VecLong & shift);

//bool is_row_popov_up_to_permutation(const Mat<zz_pX> & pmat);
//bool is_row_popov_up_to_permutation(
//                                    const Mat<zz_pX> & pmat,
//                                    const VecLong & shift
//                                    );

/** Computes and returns a boolean indicating whether `pmat` is in column-wise
 * Popov form (see @ref MatrixForms)
 *
 * \todo provide same test but relaxing ordered weak Popov to weak Popov in the
 * definition? This allows easier support for definitions of Popov with
 * different orderings of the rows, such as by increasing degree
 */
bool is_col_popov(const Mat<zz_pX> & pmat);

/** Computes and returns a boolean indicating whether `pmat` is in column-wise
 * `shift`-Popov form (see @ref MatrixForms)
 */
bool is_col_popov(const Mat<zz_pX> & pmat, const VecLong & shift);

//bool is_col_popov_up_to_permutation(const Mat<zz_pX> & pmat);
//bool is_col_popov_up_to_permutation(
//                                    const Mat<zz_pX> & pmat,
//                                    const VecLong & shift
//                                    );


/** Computes and returns a boolean indicating whether `pmat` is in row-wise
 * (non-shifted) form indicated by `form` (see @ref MatrixForms and
 * #PolMatForm). Throws an error if `form` is not among the forms defined in
 * the enumeration #PolMatForm.
 */
bool is_row_polmatform(
                       const PolMatForm form,
                       const Mat<zz_pX> & pmat
                      );

/** Computes and returns a boolean indicating whether `pmat` is in row-wise
 * `shift`-shifted form indicated by `form` (see @ref MatrixForms and
 * #PolMatForm). Throws an error if `form` is not among the forms defined
 * in the enumeration #PolMatForm.
 */
bool is_row_polmatform(
                       const PolMatForm form,
                       const Mat<zz_pX> & pmat,
                       const VecLong & shift
                      );

/** Computes and returns a boolean indicating whether `pmat` is in column-wise
 * (non-shifted) form indicated by `form` (see @ref MatrixForms and
 * #PolMatForm). Throws an error if `form` is not among the forms defined in
 * the enumeration #PolMatForm.
 */
bool is_col_polmatform(
                       const PolMatForm form,
                       const Mat<zz_pX> & pmat
                      );

/** Computes and returns a boolean indicating whether `pmat` is in column-wise
 * `shift`-shifted form indicated by `form` (see @ref MatrixForms and
 * #PolMatForm). Throws an error if `form` is not among the forms defined in
 * the enumeration #PolMatForm.
 */
bool is_col_polmatform(
                       const PolMatForm form,
                       const Mat<zz_pX> & pmat,
                       const VecLong & shift
                      );

/** Computes and returns the strongest row-wise form that `pmat` has, among
 * those defined in the enumeration #PolMatForm.
 */
PolMatForm get_row_polmatform(const Mat<zz_pX> & pmat);

/** Computes and returns the strongest row-wise `shift`-shifted form that
 * `pmat` has, among those defined in the enumeration #PolMatForm.
 */
PolMatForm get_row_polmatform(
                              const Mat<zz_pX> & pmat,
                              const VecLong & shift
                             );

/** Computes and returns the strongest column-wise form that `pmat` has, among
 * those defined in the enumeration #PolMatForm.
 */
PolMatForm get_col_polmatform(const Mat<zz_pX> & pmat);

/** Computes and returns the strongest column-wise `shift`-shifted form that
 * `pmat` has, among those defined in the enumeration #PolMatForm.
 */
PolMatForm get_col_polmatform(
                              const Mat<zz_pX> & pmat,
                              const VecLong & shift
                             );

//@} // doxygen group: Testing polynomial matrix forms

#endif /* ifndef MAT_LZZ_PX_FORMS */

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
