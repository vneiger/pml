#ifndef NMOD_POLY_MAT_FORMS_H
#define NMOD_POLY_MAT_FORMS_H

/** \brief Shifted reduced forms and shifted normal forms of univariate
 * polynomial matrices with coefficients in `nmod`
 *
 * \file nmod_poly_mat_utils.h
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

#include <flint/fmpz_mat.h> // for degree matrix
#include <flint/nmod_mat.h>
#include <flint/nmod_poly_mat.h>

#ifdef __cplusplus
extern "C" {
#endif

/** Shifted reduced forms of polynomial matrices.
 *
 * \enum poly_mat_form_t
 *
 * Note that the assigned numbers follow the "strength" of these forms:
 * - Popov implies ordered weak Popov, which implies weak Popov,
 *   which implies reduced.
 * - Upper Hermite implies upper echelon.
 * - Lower Hermite implies lower echelon.
 *
 */
typedef enum
{
    NONE = 0, /**< Arbitrary matrix, no specific form */
    REDUCED = 1, /**< Matrix in (shifted) reduced form */
    WEAK_POPOV = 2, /**< Matrix in (shifted) weak Popov form */
    ORD_WEAK_POPOV = 3, /**< Matrix in (shifted) ordered weak Popov form */
    POPOV = 4, /**< Matrix in (shifted) Popov form */
    LECHELON = 5, /**< Matrix in lower echelon form */
    LHERMITE = 6, /**< Matrix in lower Hermite form */
    UECHELON = 7, /**< Matrix in upper echelon form */
    UHERMITE = 8, /**< Matrix in upper Hermite form */
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
 * The functions below which involve a `shift` among its parameters throw an
 * error if this `shift` does not have the right length or is not `NULL`.  If
 * `NULL` is provided as input for the shift, this is understood as the uniform
 * shift `[0,...,0]` of the right length.
 *  
 */
//@{

/** Computes the `shift`-row degree `rdeg` of a polynomial matrix `pmat` (see
 * @ref RowAndColumnDegrees). The result `rdeg` must be already initialized
 * with length (at least) the number of rows of `mat`. The shift `shift`
 * must be `NULL` or have (at least) as many elements as the number of columns
 * of `mat`.
 */
void nmod_poly_mat_row_degree(slong *rdeg,
                              const nmod_poly_mat_t mat,
                              const slong *shift);

/** Computes the `shift`-column degree `cdeg` of a polynomial matrix `pmat`
 * (see @ref RowAndColumnDegrees). The result `cdeg` must be already initialized
 * with length (at least) the number of column of `mat`. The shift `shift` must
 * have (at least) as many elements as the number of rows of `mat`.
 */
void nmod_poly_mat_column_degree(slong *cdeg,
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
 * For a given orientation (row-wise/column-wise), the shifted pivot profile
 * consists of both the shifted pivot index and the shifted pivot degree.
 *
 * The functions below which involve a `shift` among its parameters throw
 * an error if this `shift` does not have the right length or is not `NULL`.
 * If `NULL` is provided as input for the shift, this is understood as the
 * uniform shift `[0,...,0]` of the right length.
 */
//@{

/** Computes the row-wise `shift`-pivot index `pivind` of a polynomial matrix
 * `mat` (see @ref Pivots). */
void nmod_poly_mat_pivot_index_rowwise(slong *pivind,
                                       const nmod_poly_mat_t mat,
                                       const slong *shift);

/** Computes the column-wise `shift`-pivot index `pivind` of a polynomial
 * matrix `mat` (see @ref Pivots). */
void nmod_poly_mat_pivot_index_columnwise(slong *pivind,
                                          const nmod_poly_mat_t mat,
                                          const slong *shift);


/** Computes the row-wise `shift`-pivot index `pivind` and `shift`-pivot degree
 * `pivdeg` of a polynomial matrix `mat` (see @ref Pivots). In the unshifted
 * case, `pivdeg` coincides with the row degree of `mat`. */
void nmod_poly_mat_pivot_profile_rowwise(slong * pivind,
                                         slong * pivdeg,
                                         const nmod_poly_mat_t mat,
                                         const slong * shift);

/** Computes the column-wise `shift`-pivot index `pivind` and `shift`-pivot
 * degree `pivdeg` of a polynomial matrix `mat` (see @ref Pivots). In the
 * unshifted case, `pivdeg` coincides with the column degree of `mat`. */
void nmod_poly_mat_pivot_profile_columnwise(slong *pivind,
                                            slong *pivdeg,
                                            const nmod_poly_mat_t mat,
                                            const slong *shift);

//@} // doxygen group: (Shifted) pivot index and pivot degree

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* ECHELON PIVOT INDEX/DEGREE                                 */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/** @name Echelon pivot index and pivot degree
 * \anchor EchelonPivots
 *
 * For a given polynomial row vector, its _lower (resp. upper) echelon pivot
 * index_ is the index of its rightmost (resp. leftmost) nonzero entry, and its
 * _lower (resp. upper) pivot degree_ is the degree of that nonzero entry. By
 * convention, for a zero row vector of length `n`, the lower (resp. upper)
 * echelon pivot index is `-1` (resp. `n`), and the lower (resp. upper) echelon
 * pivot degree is `-1`.
 *
 * For a given polynomial column vector, its _lower (resp. upper) echelon pivot
 * index_ is the index of its topmost (resp. bottommost) nonzero entry, and its
 * _lower (resp. upper) pivot degree_ is the degree of that nonzero entry.  By
 * convention, for a zero column vector of length `n`, the lower (resp.  upper)
 * echelon pivot index is `n` (resp. `-1`), and the lower (resp. upper) echelon
 * pivot degree is `-1`.
 *
 * Then, the row-wise lower echelon pivot index (resp. degree) of a polynomial
 * matrix is the tuple of the lower echelon pivot indices (resp. degrees) of
 * the rows of this matrix. The three other variants (column-wise / upper) are
 * defined similarly.
 *
 * For given orientations (row-wise/column-wise, lower/upper), the echelon
 * pivot profile consists of both the echelon pivot index and the echelon pivot
 * degree.
 */
//@{

/** Computes the row-wise lower echelon pivot index `pivind` of a polynomial
 * matrix `mat` (see @ref EchelonPivots). */
void nmod_poly_mat_lechelon_pivot_index_rowwise(slong * pivind,
                                                const nmod_poly_mat_t mat);

/** Computes the row-wise upper echelon pivot index `pivind` of a polynomial
 * matrix `mat` (see @ref EchelonPivots). */
void nmod_poly_mat_uechelon_pivot_index_rowwise(slong * pivind,
                                                const nmod_poly_mat_t mat);

/** Computes the column-wise lower echelon pivot index `pivind` of a polynomial
 * matrix `mat` (see @ref EchelonPivots). */
void nmod_poly_mat_lechelon_pivot_index_columnwise(slong * pivind,
                                                   const nmod_poly_mat_t mat);

/** Computes the column-wise upper echelon pivot index `pivind` of a polynomial
 * matrix `mat` (see @ref EchelonPivots). */
void nmod_poly_mat_uechelon_pivot_index_columnwise(slong * pivind,
                                                   const nmod_poly_mat_t mat);

/** Computes the row-wise lower echelon pivot index `pivind` and pivot
 * degree `pivdeg` of a polynomial matrix `mat` (see @ref EchelonPivots). */
void nmod_poly_mat_lechelon_pivot_profile_rowwise(slong * pivind,
                                                  slong * pivdeg,
                                                  const nmod_poly_mat_t mat);

/** Computes the row-wise upper echelon pivot index `pivind` and pivot
 * degree `pivdeg` of a polynomial matrix `mat` (see @ref EchelonPivots). */
void nmod_poly_mat_uechelon_pivot_profile_rowwise(slong * pivind,
                                                  slong * pivdeg,
                                                  const nmod_poly_mat_t mat);

/** Computes the column-wise lower echelon pivot index `pivind` and pivot
 * degree `pivdeg` of a polynomial matrix `mat` (see @ref Pivots). */
void nmod_poly_mat_lechelon_pivot_profile_columnwise(slong * pivind,
                                                     slong * pivdeg,
                                                     const nmod_poly_mat_t mat);

/** Computes the column-wise upper echelon pivot index `pivind` and pivot
 * degree `pivdeg` of a polynomial matrix `mat` (see @ref Pivots). */
void nmod_poly_mat_uechelon_pivot_profile_columnwise(slong * pivind,
                                                     slong * pivdeg,
                                                     const nmod_poly_mat_t mat);

//@} // doxygen group: echelon pivot index and pivot degree


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
 * an error if this `shift` does not have the right length. Here, the shift
 * cannot be `NULL`.
 */
//@{

/** Computes the degree matrix `degmat` of a polynomial matrix `pmat` (see @ref
 * DegreeMatrix)
 */
void nmod_poly_mat_degree_matrix(fmpz_mat_t dmat, const nmod_poly_mat_t mat);

/** Computes the row-wise `shift`-degree matrix `degmat` of a polynomial matrix
 * `pmat` (see @ref DegreeMatrix)
 */
void nmod_poly_mat_degree_matrix_row_shifted(fmpz_mat_t dmat,
                                             const nmod_poly_mat_t mat,
                                             const slong * shift);

/** Computes the column-wise `shift`-degree matrix `degmat` of a polynomial
 * matrix `pmat` (see @ref DegreeMatrix)
 */
void nmod_poly_mat_degree_matrix_column_shifted(fmpz_mat_t dmat,
                                                const nmod_poly_mat_t mat,
                                                const slong * shift);

/** Computes the `shift`-degree matrix `degmat` of a polynomial matrix `pmat`
 * (see @ref DegreeMatrix), the orientation row-wise/column-wise being
 * indicated by a parameter.
 */
NMOD_POLY_MAT_INLINE void
nmod_poly_mat_degree_matrix_shifted(fmpz_mat_t dmat,
                                    const nmod_poly_mat_t mat,
                                    const slong * shift,
                                    orientation_t row_wise)
{
    if (row_wise)
        nmod_poly_mat_degree_matrix_row_shifted(dmat, mat, shift);
    else
        nmod_poly_mat_degree_matrix_column_shifted(dmat, mat, shift);
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
 * `pmat[i][j]` does not reach `rdeg[i]`; or if the row pmat[i] is zero).
 * Similarly, writing `cdeg` for the column degree of `pmat`, the column-wise
 * leading matrix of `pmat` is the `m x n` matrix over the base field whose
 * entry `(i,j)` is the coefficient of degree `cdeg[j]` of the entry
 * `pmat[i][j]` (this is zero if `pmat[i][j]` does not reach `cdeg[j]`, or if
 * the column pmat[:][j] is zero).
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
 * an error if this `shift` does not have the right length or is not `NULL`.
 * If `NULL` is provided as input for the shift, this is understood as the
 * uniform shift `[0,...,0]` of the right length.
 *
 * \todo enhancement: offer row-wise (resp column-wise) leading matrix
 * when the row degree (resp column degree) is already known
 */
//@{

/** Computes the row-wise `shift`-leading matrix `lmat` of a polynomial matrix
 * `mat` (see @ref LeadingMatrix)
 */
void nmod_poly_mat_leading_matrix_rowwise(nmod_mat_t lmat,
                                          const nmod_poly_mat_t mat,
                                          const slong *shift);

/** Computes the column-wise `shift`-leading matrix `lmat` of a polynomial
 * matrix `mat` (see @ref LeadingMatrix)
 */
void nmod_poly_mat_leading_matrix_columnwise(nmod_mat_t lmat,
                                             const nmod_poly_mat_t mat,
                                             const slong *shift);


/** Computes the column-wise `shift`-leading matrix `lmat` of a polynomial
 * matrix `mat` (see @ref LeadingMatrix), using provided orientation
 * row-wise or column-wise.
 */
NMOD_POLY_MAT_INLINE void
nmod_poly_mat_leading_matrix(nmod_mat_t lmat,
                             const nmod_poly_mat_t mat,
                             const slong *shift,
                             orientation_t row_wise)
{
    if (row_wise)
        nmod_poly_mat_leading_matrix_rowwise(lmat, mat, shift);
    else
        nmod_poly_mat_leading_matrix_columnwise(lmat, mat, shift);
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
 * an error if this `shift` does not have the right length or is not `NULL`.
 * If `NULL` is provided as input for the shift, this is understood as the
 * uniform shift `[0,...,0]` of the right length.
 *  
 * \todo define lower/upper row-wise/column-wise echelon/Hermite forms
 */
//@{

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* TESTING MATRIX FORMS - REDUCED                             */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/** Tests whether `pmat` is in `shift`-row reduced form (see @ref MatrixForms) */
int nmod_poly_mat_is_reduced_rowwise(const nmod_poly_mat_t mat,
                                     const slong *shift);

/** Tests whether `pmat` is in `shift`-column reduced form (see @ref MatrixForms) */
int nmod_poly_mat_is_reduced_columnwise(const nmod_poly_mat_t mat,
                                        const slong *shift);

/** Tests whether `pmat` is in `shift`-reduced form (see @ref MatrixForms),
 * with orientation row-wise or column-wise specified by argument
 */
NMOD_POLY_MAT_INLINE int
nmod_poly_mat_is_reduced(const nmod_poly_mat_t mat,
                                 const slong *shift,
                                 orientation_t row_wise)
{
    if (row_wise)
        return nmod_poly_mat_is_reduced_rowwise(mat, shift);
    else
        return nmod_poly_mat_is_reduced_columnwise(mat, shift);
}


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* TESTING MATRIX FORMS - ORDERED WEAK POPOV                  */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/** Tests whether `mat` is in row-wise `shift`-ordered weak Popov form (see
 * @ref MatrixForms) */
int nmod_poly_mat_is_ordered_weak_popov_rowwise(const nmod_poly_mat_t mat,
                                                const slong *shift);

/** Tests whether `mat` is in column-wise `shift`-ordered weak Popov form (see
 * @ref MatrixForms) */
int nmod_poly_mat_is_ordered_weak_popov_columnwise(const nmod_poly_mat_t mat,
                                                   const slong *shift);

/** Tests whether `mat` is in `shift`-ordered weak Popov form (see
 * @ref MatrixForms), with orientation specified by `row_wise` */
NMOD_POLY_MAT_INLINE int
nmod_poly_mat_is_ordered_weak_popov(const nmod_poly_mat_t mat,
                                    const slong *shift,
                                    orientation_t row_wise)
{
    if (row_wise)
        return nmod_poly_mat_is_ordered_weak_popov_rowwise(mat, shift);
    else
        return nmod_poly_mat_is_ordered_weak_popov_columnwise(mat, shift);
}


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* TESTING MATRIX FORMS - WEAK POPOV                          */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/** Tests whether `mat` is in row-wise `shift`-weak Popov form (see
 * @ref MatrixForms) */
int nmod_poly_mat_is_weak_popov_rowwise(const nmod_poly_mat_t mat,
                                        const slong *shift);

/** Tests whether `mat` is in column-wise `shift`-weak Popov form (see
 * @ref MatrixForms) */
int nmod_poly_mat_is_weak_popov_columnwise(const nmod_poly_mat_t mat,
                                           const slong *shift);

/** Tests whether `mat` is in `shift`-weak Popov form (see @ref MatrixForms),
 * with orientation specified by `row_wise` */
NMOD_POLY_MAT_INLINE int
nmod_poly_mat_is_weak_popov(const nmod_poly_mat_t mat,
                            const slong *shift,
                            orientation_t row_wise)
{
    if (row_wise)
        return nmod_poly_mat_is_weak_popov_rowwise(mat, shift);
    else
        return nmod_poly_mat_is_weak_popov_columnwise(mat, shift);
}



/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* TESTING MATRIX FORMS - POPOV                               */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/** \todo provide same test but relaxing ordered weak Popov to weak Popov in the
 * definition? This allows easier support for definitions of Popov with
 * different orderings of the rows, such as by increasing degree */


/** Tests whether `mat` is in row-wise `shift`-Popov form (see
 * @ref MatrixForms) */
int nmod_poly_mat_is_popov_rowwise(const nmod_poly_mat_t mat,
                                   const slong *shift);

/** Tests whether `mat` is in column-wise `shift`-Popov form (see
 * @ref MatrixForms) */
int nmod_poly_mat_is_popov_columnwise(const nmod_poly_mat_t mat,
                                      const slong *shift);

/** Tests whether `mat` is in `shift`-Popov form (see @ref MatrixForms),
 * with orientation specified by `row_wise` */
NMOD_POLY_MAT_INLINE int
nmod_poly_mat_is_popov(const nmod_poly_mat_t mat,
                       const slong *shift,
                       orientation_t row_wise)
{
    if (row_wise)
        return nmod_poly_mat_is_popov_rowwise(mat, shift);
    else
        return nmod_poly_mat_is_popov_columnwise(mat, shift);
}


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* TESTING MATRIX FORMS - ECHELON                             */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/* def of lower row echelon:
 * - there is no zero row below a nonzero row (possible zero rows are grouped at the top),
 * - let r be the number of nonzero rows (thus at indices m-r,...,m-1) and j_i be
 *   the index of the rightmost nonzero entry in row m-r+i for 0 <= i < r, then j_0 < ... < j_{r-1}.
 * ==> equivalently: the row-wise lower echelon pivot profile of the matrix is increasing,
 * excepting zero rows */
int nmod_poly_mat_is_lechelon_rowwise(const nmod_poly_mat_t pmat);
/* def of upper row echelon:
 * - there is no zero row above a nonzero row (possible zero rows are grouped at the bottom),
 * - let r be the number of nonzero rows (thus at indices 0,...,r-1) and j_i be
 *   the index of the rightmost nonzero entry in row i for 0 <= i < r, then j_0 < ... < j_{r-1}.
 * ==> equivalently: the row-wise lower echelon pivot profile of the matrix is increasing,
 * excepting zero rows */
int nmod_poly_mat_is_uechelon_rowwise(const nmod_poly_mat_t pmat);
/* etc... */
int nmod_poly_mat_is_lechelon_columnwise(const nmod_poly_mat_t pmat);
/* etc... */
int nmod_poly_mat_is_uechelon_columnwise(const nmod_poly_mat_t pmat);

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* TESTING MATRIX FORMS - HERMITE                             */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/** Tests whether `mat` is in row-wise lower Hermite form (see @ref MatrixForms) */
int nmod_poly_mat_is_lhermite_rowwise(const nmod_poly_mat_t pmat);
/** Tests whether `mat` is in row-wise upper Hermite form (see @ref MatrixForms) */
int nmod_poly_mat_is_uhermite_rowwise(const nmod_poly_mat_t pmat);
/** Tests whether `mat` is in column-wise lower Hermite form (see @ref MatrixForms) */
int nmod_poly_mat_is_lhermite_columnwise(const nmod_poly_mat_t pmat);
/** Tests whether `mat` is in column-wise upper Hermite form (see @ref MatrixForms) */
int nmod_poly_mat_is_uhermite_columnwise(const nmod_poly_mat_t pmat);




/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* TESTING MATRIX FORMS - USER INTERFACE                      */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

// TODO implement general interfaces, or just the above suffice?

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

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* COMPUTING MATRIX FORMS                                     */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/** @name Computing polynomial matrix forms
 * See @ref MatrixForms for definitions
 */
//@{

/** Sets ``hnf`` to the row-wise, lower Hermite normal form of ``mat``, and
 * returns the rank of ``mat``. Aliasing is allowed */
slong nmod_poly_mat_uhermite_form_rowwise_naive(nmod_poly_mat_t hnf, const nmod_poly_mat_t mat);


//@} // doxygen group: Computing polynomial matrix forms


#ifdef __cplusplus
}
#endif

#endif // NMOD_POLY_MAT_FORMS_H

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
