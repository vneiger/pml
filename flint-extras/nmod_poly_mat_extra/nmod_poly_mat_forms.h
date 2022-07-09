#ifndef NMOD_POLY_MAT_FORMS_H
#define NMOD_POLY_MAT_FORMS_H

/** \brief Shifted reduced forms and shifted normal forms of univariate
 * polynomial matrices with coefficients in `nmod`
 *
 * \file nmod_poly_mat_utils.h
 * \author Vincent Neiger, Kevin Tran
 * \version 0.0
 * \date 2022-06-25
 *
 * Basic functions to deal with shifted reduced and shifted normal forms of
 * polynomials matrices: test if a matrix is in some given form, compute
 * shifted row/column degrees and shifted pivot degrees, compute shifted
 * leading matrix.
 *
 * \todo state that shift has to have right dimension (depending on
 * orientation), and this is never checked in implementations
 * \todo random matrix with given PolMatForm
 *
 */

// for degree matrix
#include <flint/fmpz_mat.h> 
#include <flint/nmod_mat.h>
#include <flint/nmod_poly_mat.h>

#ifdef __cplusplus
extern "C" {
#endif

/** Shifted reduced forms of polynomial matrices.
 *
 * \enum poly_mat_form_t
 *
 * Note that the assigned numbers follow the "strength" of these forms: Popov
 * implies ordered weak Popov, which implies weak Popov, which implies reduced.
 */
typedef enum
{
    NONE = 0, /**< Arbitrary matrix, no specific form */
    REDUCED = 1, /**< Matrix in (shifted) reduced form */
    WEAK_POPOV = 2, /**< Matrix in (shifted) weak Popov form */
    ORD_WEAK_POPOV = 3, /**< Matrix in (shifted) ordered weak Popov form */
    POPOV = 4, /**< Matrix in (shifted) Popov form */
} poly_mat_form_t;


/**
 * \enum orientation_t
 * \brief Whether to focus on row space or on column space of a polynomial matrix 
 *
 */
typedef enum
{
    COLUMN_WISE = 0,
    ROW_WISE = 1
} orientation_t;





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
 * @ref RowAndColumnDegrees). The result `rdeg` must be already initialized
 * with length (at least) the number of rows of `mat`.
 */
void row_degrees(slong *rdeg,
                 const nmod_poly_mat_t mat);

/** Computes the `shift`-row degree `rdeg` of a polynomial matrix `pmat` (see
 * @ref RowAndColumnDegrees). The result `rdeg` must be already initialized
 * with length (at least) the number of rows of `mat`. The shift `shift`
 * must have (at least) as many elements as the number of columns of `mat`.
 */
void row_degrees_shifted(slong *rdeg,
                         const nmod_poly_mat_t mat,
                         const slong *shift);

/** Computes the column degree `cdeg` of a polynomial matrix `pmat` (see
 * @ref RowAndColumnDegrees). The result `cdeg` must be already initialized
 * with length (at least) the number of column of `mat`.
 */
void column_degrees(slong *cdeg,
                    const nmod_poly_mat_t mat);

/** Computes the `shift`-column degree `cdeg` of a polynomial matrix `pmat`
 * (see @ref RowAndColumnDegrees). The result `cdeg` must be already initialized
 * with length (at least) the number of column of `mat`. The shift `shift` must
 * have (at least) as many elements as the number of rows of `mat`.
 */
void column_degrees_shifted(slong *cdeg,
                            const nmod_poly_mat_t mat,
                            const slong *shift);

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
 * The shifted row pivot profile consists of both the shifted pivot index and
 * the shifted row (resp. column) degree.
 *
 * The functions below which involve a `shift` among its parameters throw
 * an error if this `shift` does not have the right length.
 */
//@{

/** Computes the row-wise pivot index `pivind` of a polynomial matrix `mat`
 * (see @ref Pivots). */
void pivot_index_rowwise(slong *pivind,
                         const nmod_poly_mat_t mat);

/** Computes the column-wise pivot index `pivind` of a polynomial matrix `mat`
 * (see @ref Pivots). */
void pivot_index_columnwise(slong *pivind,
                            const nmod_poly_mat_t mat);

/** Computes the row-wise `shift`-pivot index `pivind` of a polynomial matrix
 * `mat` (see @ref Pivots). */
void pivot_index_shifted_rowwise(slong *pivind,
                                 const nmod_poly_mat_t mat,
                                 const slong *shift);

/** Computes the column-wise `shift`-pivot index `pivind` of a polynomial
 * matrix `mat` (see @ref Pivots). */
void pivot_index_shifted_columnwise(slong *pivind,
                                    const nmod_poly_mat_t mat,
                                    const slong *shift);

/** Computes the row-wise pivot index `pivind` and pivot degree `pivdeg` of a
 * polynomial matrix `mat` (see @ref Pivots). In this unshifted case, `pivdeg`
 * coincides with the row degree of `mat`.  */
void pivot_profile_rowwise(slong *pivind,
                           slong *pivdeg,
                           const nmod_poly_mat_t mat);

/** Computes the column-wise pivot index `pivind` and pivot degree `pivdeg` of
 * a polynomial matrix `mat` (see @ref Pivots). In this unshifted case, `pivdeg`
 * coincides with the column degree of `mat`. */
void pivot_profile_columnwise(slong *pivind,
                              slong *pivdeg,
                              const nmod_poly_mat_t mat);

/** Computes the row-wise `shift`-pivot index `pivind` and `shift`-pivot degree
 * `pivdeg` of a polynomial matrix `mat` (see @ref Pivots). */
void pivot_profile_shifted_rowwise(slong *pivind,
                                   slong *pivdeg,
                                   const nmod_poly_mat_t mat,
                                   const slong *shift);

/** Computes the column-wise `shift`-pivot index `pivind` and `shift`-pivot
 * degree `pivdeg` of a polynomial matrix `mat` (see @ref Pivots). */
void pivot_profile_shifted_columnwise(slong *pivind,
                                      slong *pivdeg,
                                      const nmod_poly_mat_t mat,
                                      const slong *shift);

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
 * `deg(pmat[i][j]) + shift[j]` if `pmat[i][j]` is nonzero, and `shift[j]-1`
 * otherwise. For a given shift `shift` of length `m`, the column-wise
 * `shift`-degree matrix of `pmat` is the `m x n` integer matrix whose entry
 * `(i,j)` is `deg(pmat[i][j]) + shift[i]` if `pmat[i][j]` is nonzero, and
 * `shift[i]-1` otherwise.
 *
 * The functions below which involve a `shift` among its parameters throw
 * an error if this `shift` does not have the right length.
 */
//@{

/** Computes the degree matrix `degmat` of a polynomial matrix `pmat` (see @ref
 * DegreeMatrix)
 */
void degree_matrix(fmpz_mat_t dmat, const nmod_poly_mat_t mat);

/** Computes the row-wise `shift`-degree matrix `degmat` of a polynomial matrix
 * `pmat` (see @ref DegreeMatrix)
 */
void degree_matrix_row_shifted(fmpz_mat_t dmat,
                           const nmod_poly_mat_t mat,
                           const slong * shift);

/** Computes the column-wise `shift`-degree matrix `degmat` of a polynomial
 * matrix `pmat` (see @ref DegreeMatrix)
 */
void degree_matrix_column_shifted(fmpz_mat_t dmat,
                           const nmod_poly_mat_t mat,
                           const slong * shift);

/** Computes the `shift`-degree matrix `degmat` of a polynomial matrix `pmat`
 * (see @ref DegreeMatrix), the orientation row-wise/column-wise being
 * indicated by a parameter.
 */
NMOD_POLY_MAT_INLINE void
degree_matrix_shifted(fmpz_mat_t dmat,
                      const nmod_poly_mat_t mat,
                      const slong * shift,
                      orientation_t row_wise)
{
    if (row_wise)
        degree_matrix_row_shifted(dmat, mat, shift);
    else
        degree_matrix_column_shifted(dmat, mat, shift);
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
 * \todo enhancement: offer row-wise (resp column-wise) leading matrix
 * when the row degree (resp column degree) is already known
 */
//@{

/** Computes the row-wise leading matrix `lmat` of a polynomial matrix `mat`
 * (see @ref LeadingMatrix)
 * \todo bug to check, see src
 */
void leading_matrix_rowwise(nmod_mat_t lmat,
                            const nmod_poly_mat_t mat);


/** Computes the column-wise leading matrix `lmat` of a polynomial matrix
 * `mat` (see @ref LeadingMatrix)
 * \todo bug to check, see src
 */
void leading_matrix_columnwise(nmod_mat_t lmat,
                               const nmod_poly_mat_t mat);

/** Computes the leading matrix `lmat` of a polynomial matrix `mat` (see @ref
 * LeadingMatrix), using provided orientation row-wise or column-wise.
 * \todo bug to check, see src
 */
NMOD_POLY_MAT_INLINE void
leading_matrix(nmod_mat_t lmat,
               const nmod_poly_mat_t mat,
               orientation_t row_wise)
{
    if (row_wise)
        leading_matrix_rowwise(lmat, mat);
    else
        leading_matrix_columnwise(lmat, mat);
}



/** Computes the row-wise `shift`-leading matrix `lmat` of a polynomial matrix
 * `mat` (see @ref LeadingMatrix)
 * \todo bug to check, see src
 */
void leading_matrix_shifted_rowwise(nmod_mat_t lmat,
                                    const nmod_poly_mat_t mat,
                                    const slong *shift);

/** Computes the column-wise `shift`-leading matrix `lmat` of a polynomial
 * matrix `mat` (see @ref LeadingMatrix)
 * \todo bug to check, see src
 */
void leading_matrix_shifted_columnwise(nmod_mat_t lmat,
                                       const nmod_poly_mat_t mat,
                                       const slong *shift);


/** Computes the column-wise `shift`-leading matrix `lmat` of a polynomial
 * matrix `mat` (see @ref LeadingMatrix), using provided orientation
 * row-wise or column-wise.
 * \todo bug to check, see src
 */
NMOD_POLY_MAT_INLINE void
leading_matrix_shifted(nmod_mat_t lmat,
                       const nmod_poly_mat_t mat,
                       const slong *shift,
                       orientation_t row_wise)
{
    if (row_wise)
        leading_matrix_shifted_rowwise(lmat, mat, shift);
    else
        leading_matrix_shifted_columnwise(lmat, mat, shift);
}


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

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* TESTING MATRIX FORMS - REDUCED                             */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/** Tests whether `pmat` is in row reduced form (see @ref MatrixForms) */
int is_reduced_rowwise(const nmod_poly_mat_t mat);

/** Tests whether `pmat` is in `shift`-row reduced form (see @ref MatrixForms) */
int is_reduced_shifted_rowwise(const nmod_poly_mat_t mat, const slong *shift);

/** Tests whether `pmat` is in column reduced form (see @ref MatrixForms) */
int is_reduced_columnwise(const nmod_poly_mat_t mat);

/** Tests whether `pmat` is in `shift`-column reduced form (see @ref MatrixForms) */
int is_reduced_shifted_columnwise(const nmod_poly_mat_t mat, const slong *shift);

/** Tests whether `pmat` is in reduced form (see @ref MatrixForms), with
 * orientation row-wise or column-wise specified by argument
 */
NMOD_POLY_MAT_INLINE int
is_reduced(const nmod_poly_mat_t mat, orientation_t row_wise)
{
    if (row_wise)
        return is_reduced_rowwise(mat);
    else
        return is_reduced_columnwise(mat);
}

/** Tests whether `pmat` is in `shift`-reduced form (see @ref MatrixForms),
 * with orientation row-wise or column-wise specified by argument
 */
NMOD_POLY_MAT_INLINE int
is_reduced_shifted(const nmod_poly_mat_t mat,
                   const slong *shift,
                   orientation_t row_wise)
{
    if (row_wise)
        return is_reduced_shifted_rowwise(mat, shift);
    else
        return is_reduced_shifted_columnwise(mat, shift);
}


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* TESTING MATRIX FORMS - ORDERED WEAK POPOV                  */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/** Tests whether `mat` is in row-wise ordered weak Popov form (see
 * @ref MatrixForms) */
int is_ordered_weak_popov_rowwise(const nmod_poly_mat_t mat);

/** Tests whether `mat` is in row-wise `shift`-ordered weak Popov form (see
 * @ref MatrixForms) */
int is_ordered_weak_popov_shifted_rowwise(const nmod_poly_mat_t mat,
                                          const slong *shift);

/** Tests whether `mat` is in column-wise ordered weak Popov form (see
 * @ref MatrixForms) */
int is_ordered_weak_popov_columnwise(const nmod_poly_mat_t mat);

/** Tests whether `mat` is in column-wise `shift`-ordered weak Popov form (see
 * @ref MatrixForms) */
int is_ordered_weak_popov_shifted_columnwise(const nmod_poly_mat_t mat,
                                             const slong *shift);

/** Tests whether `mat` is in ordered weak Popov form (see @ref MatrixForms),
 * with orientation specified by `row_wise` */
NMOD_POLY_MAT_INLINE int
is_ordered_weak_popov(const nmod_poly_mat_t mat,
              orientation_t row_wise)
{
    if (row_wise)
        return is_ordered_weak_popov_rowwise(mat);
    else
        return is_ordered_weak_popov_columnwise(mat);
}

/** Tests whether `mat` is in `shift`-ordered weak Popov form (see
 * @ref MatrixForms), with orientation specified by `row_wise` */
NMOD_POLY_MAT_INLINE int
is_ordered_weak_popov_shifted(const nmod_poly_mat_t mat,
                      const slong *shift,
                      orientation_t row_wise)
{
    if (row_wise)
        return is_ordered_weak_popov_shifted_rowwise(mat, shift);
    else
        return is_ordered_weak_popov_shifted_columnwise(mat, shift);
}


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* TESTING MATRIX FORMS - WEAK POPOV                          */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/** Tests whether `mat` is in row-wise weak Popov form (see @ref MatrixForms) */
int is_weak_popov_rowwise(const nmod_poly_mat_t mat);

/** Tests whether `mat` is in row-wise `shift`-weak Popov form (see
 * @ref MatrixForms) */
int is_weak_popov_shifted_rowwise(const nmod_poly_mat_t mat,
                                  const slong *shift);

/** Tests whether `mat` is in column-wise ordered weak Popov form (see
 * @ref MatrixForms) */
int is_weak_popov_columnwise(const nmod_poly_mat_t mat);

/** Tests whether `mat` is in column-wise `shift`-ordered weak Popov form (see
 * @ref MatrixForms) */
int is_weak_popov_shifted_columnwise(const nmod_poly_mat_t mat,
                                     const slong *shift);

/** Tests whether `mat` is in weak Popov form (see @ref MatrixForms),
 * with orientation specified by `row_wise` */
NMOD_POLY_MAT_INLINE int
is_weak_popov(const nmod_poly_mat_t mat,
              orientation_t row_wise)
{
    if (row_wise)
        return is_weak_popov_rowwise(mat);
    else
        return is_weak_popov_columnwise(mat);
}

/** Tests whether `mat` is in `shift`-weak Popov form (see @ref MatrixForms),
 * with orientation specified by `row_wise` */
NMOD_POLY_MAT_INLINE int
is_weak_popov_shifted(const nmod_poly_mat_t mat,
                  const slong *shift,
                  orientation_t row_wise)
{
    if (row_wise)
        return is_weak_popov_shifted_rowwise(mat, shift);
    else
        return is_weak_popov_shifted_columnwise(mat, shift);
}



/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* TESTING MATRIX FORMS - POPOV                               */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/



/** \todo provide same test but relaxing ordered weak Popov to weak Popov in the
 * definition? This allows easier support for definitions of Popov with
 * different orderings of the rows, such as by increasing degree */

/** Tests whether `mat` is in row-wise Popov form (see @ref MatrixForms) */
// TODO
//int is_popov_rowwise(const nmod_poly_mat_t mat);

/** Tests whether `mat` is in row-wise `shift`-Popov form (see
 * @ref MatrixForms) */
// TODO
//int is_popov_shifted_rowwise(const nmod_poly_mat_t mat,
//                             const slong *shift);

/** Tests whether `mat` is in column-wise Popov form (see @ref MatrixForms) */
// TODO
//int is_popov_columnwise(const nmod_poly_mat_t mat);

/** Tests whether `mat` is in column-wise `shift`-Popov form (see
 * @ref MatrixForms) */
// TODO
//int is_popov_shifted_columnwise(const nmod_poly_mat_t mat,
//                                const slong *shift);

///** Tests whether `mat` is in `shift`-Popov form (see @ref MatrixForms), with
// * orientation specified by argument `row_wise` */
//NMOD_POLY_MAT_INLINE int
//is_popov(const nmod_poly_mat_t mat,
//         orientation_t row_wise)
//{
//    if (row_wise)
//        return is_popov_rowwise(mat);
//    else
//        return is_popov_columnwise(mat);
//}
//
///** Tests whether `mat` is in `shift`-Popov form (see @ref MatrixForms), with
// * orientation specified by argument `row_wise` */
//NMOD_POLY_MAT_INLINE int
//is_popov_shifted(const nmod_poly_mat_t mat,
//                 const slong *shift,
//                 orientation_t row_wise)
//{
//    if (row_wise)
//        return is_popov_shifted_rowwise(mat, shift);
//    else
//        return is_popov_shifted_columnwise(mat, shift);
//}



/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* TESTING MATRIX FORMS - HERMITE                             */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

// TODO
///** Tests whether `mat` is in row-wise Hermite form (see @ref MatrixForms)
// * \todo there are two definitions... upper vs lower
// **/
//int is_hermite_rowwise(const nmod_poly_mat_t mat);
//
///** Tests whether `mat` is in column-wise Hermite form (see @ref MatrixForms)
// * \todo there are two definitions... upper vs lower
// **/
//int is_hermite_columnwise(const nmod_poly_mat_t mat);
//
///** Tests whether `mat` is in `shift`-Hermite form (see @ref MatrixForms), with
// * orientation specified by argument `row_wise`
// * \todo there are two definitions... upper vs lower
// **/
//NMOD_POLY_MAT_INLINE int
//is_hermite(const nmod_poly_mat_t mat,
//           orientation_t row_wise)
//{
//    if (row_wise)
//        return is_hermite_rowwise(mat);
//    else
//        return is_hermite_columnwise(mat);
//}
//










/** Computes and returns a boolean indicating whether `pmat` is in row-wise
 * (non-shifted) form indicated by `form` (see @ref MatrixForms and
 * #PolMatForm). Throws an error if `form` is not among the forms defined in
 * the enumeration #PolMatForm.
 * \todo
 */
//bool is_row_polmatform(
//                       const PolMatForm form,
//                       const Mat<zz_pX> & pmat
//                      );

/** Computes and returns a boolean indicating whether `pmat` is in row-wise
 * `shift`-shifted form indicated by `form` (see @ref MatrixForms and
 * #PolMatForm). Throws an error if `form` is not among the forms defined
 * in the enumeration #PolMatForm.
 * \todo
 */
//bool is_row_polmatform(
//                       const PolMatForm form,
//                       const Mat<zz_pX> & pmat,
//                       const VecLong & shift
//                      );

/** Computes and returns a boolean indicating whether `pmat` is in column-wise
 * (non-shifted) form indicated by `form` (see @ref MatrixForms and
 * #PolMatForm). Throws an error if `form` is not among the forms defined in
 * the enumeration #PolMatForm.
 * \todo
 */
//bool is_column_polmatform(
//                       const PolMatForm form,
//                       const Mat<zz_pX> & pmat
//                      );

/** Computes and returns a boolean indicating whether `pmat` is in column-wise
 * `shift`-shifted form indicated by `form` (see @ref MatrixForms and
 * #PolMatForm). Throws an error if `form` is not among the forms defined in
 * the enumeration #PolMatForm.
 * \todo
 */
//bool is_column_polmatform(
//                       const PolMatForm form,
//                       const Mat<zz_pX> & pmat,
//                       const VecLong & shift
//                      );

/** Computes and returns the strongest row-wise form that `pmat` has, among
 * those defined in the enumeration #PolMatForm.
 * \todo
 */
//PolMatForm get_row_polmatform(const Mat<zz_pX> & pmat);

/** Computes and returns the strongest row-wise `shift`-shifted form that
 * `pmat` has, among those defined in the enumeration #PolMatForm.
 * \todo
 */
//PolMatForm get_row_polmatform(
//                              const Mat<zz_pX> & pmat,
//                              const VecLong & shift
//                             );

/** Computes and returns the strongest column-wise form that `pmat` has, among
 * those defined in the enumeration #PolMatForm.
 * \todo
 */
//PolMatForm get_column_polmatform(const Mat<zz_pX> & pmat);

/** Computes and returns the strongest column-wise `shift`-shifted form that
 * `pmat` has, among those defined in the enumeration #PolMatForm.
 * \todo
 */
//PolMatForm get_column_polmatform(
//                              const Mat<zz_pX> & pmat,
//                              const VecLong & shift
//                             );

//@} // doxygen group: Testing polynomial matrix forms


#ifdef __cplusplus
}
#endif

#endif // NMOD_POLY_MAT_FORMS_H

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
