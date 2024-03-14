#ifndef NMOD_POLY_MAT_FORMS_H
#define NMOD_POLY_MAT_FORMS_H

/** \brief Shifted reduced forms and shifted normal forms of univariate
 * polynomial matrices with coefficients in `nmod`
 *
 * \file nmod_poly_mat_utils.h
 * \version 0.3
 * \date 2023-11-20
 *
 * - Basic functions to deal with shifted reduced and shifted normal forms of
 *   polynomials matrices: test if a matrix is in some given form, compute
 *   shifted row/column degrees and shifted pivot degrees, compute shifted
 *   leading matrix.
 *
 * - Functions for computing shifted reduced and shifted normal forms of
 *   polynomial matrices.
 *
 */

#include <flint/fmpz_types.h> // for fmpz_mat (degree matrix)
//#include <flint/nmod_types.h>
#include <flint/nmod_poly_mat.h> // fmpz_types not enough, need NMOD_POLY_MAT_INLINE

#ifdef __cplusplus
extern "C" {
#endif


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* ENUMS: FORMS AND ORIENTATION                               */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/** Shifted reduced forms of polynomial matrices.
 *
 * \enum poly_mat_form_t
 *
 * Note that the assigned numbers follow the "strength" of these forms:
 * - 0<-1<-2<-3<-4: Popov implies ordered weak Popov, which implies weak Popov,
 *   which implies reduced.
 * - 0<-5<-6: Hermite implies echelon.
 *
 */
typedef enum
{
    NONE = 0, /**< Arbitrary matrix, no specific form */
    REDUCED = 1, /**< Matrix in (shifted) reduced form */
    WEAK_POPOV = 2, /**< Matrix in (shifted) weak Popov form */
    ORD_WEAK_POPOV = 3, /**< Matrix in (shifted) ordered weak Popov form */
    POPOV = 4, /**< Matrix in (shifted) Popov form */
    ECHELON = 5, /**< Matrix in echelon form */
    HERMITE = 6, /**< Matrix in Hermite form */
} poly_mat_form_t;

/**
 * \enum orientation_t
 * \anchor orientation
 * \brief Whether to focus on row space or on column space of a polynomial matrix,
 * and whether to consider upper or lower forms.
 *
 * Row-wise: focus on row space, work with left unimodular transformation.
 *
 * Column-wise: focus on columns space, work with right unimodular
 * transformation.
 *
 * ROW_LOWER and COL_UPPER are related by means of matrix transposition.
 * ROW_LOWER and ROW_UPPER are related by means of inverting both rows
 * and columns (i.e. pre- and post-multiplying by the antidiagonal identity
 * matrix).
 *
 * More details regarding shifted forms (see @ref MatrixForms):
 *
 * Upper: shifted weak Popov form means leading matrix in upper echelon form;
 * for row-wise this means pivot is the leftmost entry of largest shifted
 * degree in a row (or, for Hermite form, leftmost nonzero entry), and for
 * column-wise this means pivot is the bottommost entry of largest shifted
 * degree in a column (or, for Hermite form, bottommost nonzero entry).
 *
 * Lower: shifted weak Popov form means leading matrix in lower echelon form;
 * for row-wise this means pivot is the rightmost entry of largest shifted
 * degree in a row (or, for Hermite form, rightmost nonzero entry) and for
 * column-wise this means pivot is the topmost entry of largest shifted degree
 * in a column (or, for Hermite form, topmost nonzero entry).
 */
typedef enum
{
    COL_UPPER = 0,
    COL_LOWER = 1,
    ROW_UPPER = 2,
    ROW_LOWER = 3,
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
 * The functions below which involve a `shift` among its parameters do not
 * check whether `shift` has the right length. Most functions accept that
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
                              const slong * shift);

/** Computes the `shift`-column degree `cdeg` of a polynomial matrix `pmat`
 * (see @ref RowAndColumnDegrees). The result `cdeg` must be already initialized
 * with length (at least) the number of column of `mat`. The shift `shift` must
 * have (at least) as many elements as the number of rows of `mat`.
 */
void nmod_poly_mat_column_degree(slong *cdeg,
                                 const nmod_poly_mat_t mat,
                                 const slong * shift);

//@} // doxygen group: (Shifted) row and column degree

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* (SHIFTED) PIVOT INDEX/DEGREE                               */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/** @name (Shifted) pivot index and pivot degree
 * \anchor Pivots
 *
 * Four orientations are possible, see @ref orientation_t. Here we describe
 * shifted pivot index, degree, and pivot profile assuming orientation row-wise
 * and lower (ROW_LOWER). The translation of these definitions to the other
 * cases is straightforward.
 *
 * For a given shift and a given polynomial row vector, its _shifted pivot
 * index_ is the largest index corresponding to an entry which reaches the
 * shifted degree of the vector, and its _shifted pivot degree_ is the degree
 * of that entry (without adding the shift entry). By convention, both the
 * shifted pivot index and the shifted pivot degree of a zero vector are -1.
 *
 * Then, the shifted pivot index (resp. degree) of a polynomial matrix
 * is the tuple of the shifted pivot indices (resp. degrees) of the rows of
 * this matrix.
 *
 * The shifted pivot profile consists of both the shifted pivot index and the
 * shifted pivot degree.
 *
 * The functions below which involve a `shift` among its parameters do not
 * check whether `shift` has the right length. Most functions accept that
 * `NULL` is provided as input for the shift, this is understood as the uniform
 * shift `[0,...,0]` of the right length.
 *
 * Note: for ROW_UPPER/COLUMN_LOWER, the zero vector has pivot index equal to
 * the length of this vector.
 */
//@{

/** Computes the `shift`-pivot index (stored in integer `pivind`) and
 * `shift`-pivot degree (stored in integer `pivdeg`) of a given vector `vec`
 * (see @ref Pivots). In the unshifted case, `pivdeg` coincides with the degree
 * of this vector. ROW_LOWER / COLUMN_UPPER orientation. */
void _nmod_poly_vec_pivot_profile(slong * pivind,
                                  slong * pivdeg,
                                  const nmod_poly_struct * vec,
                                  const slong * shift,
                                  slong len,
                                  orientation_t orient);

/** Computes the `shift`-pivot index `pivind` of a polynomial matrix
 * `mat` (see @ref Pivots). */
void nmod_poly_mat_pivot_index(slong *pivind,
                               const nmod_poly_mat_t mat,
                               const slong * shift,
                               orientation_t orient);

/** Computes the row-wise `shift`-pivot index `pivind` and `shift`-pivot degree
 * `pivdeg` of a polynomial matrix `mat` (see @ref Pivots). In the unshifted
 * case, `pivdeg` coincides with the row or column degree of `mat` (depending
 * on column-wise/row-wise orientation). */
void nmod_poly_mat_pivot_profile(slong * pivind,
                                 slong * pivdeg,
                                 const nmod_poly_mat_t mat,
                                 const slong * shift,
                                 orientation_t orient);

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
 * For a given orientation, the echelon pivot profile consists of both the
 * echelon pivot index and the echelon pivot degree.
 */
//@{

/** Computes the echelon pivot index `pivind` of a polynomial matrix `mat` (see
 * @ref EchelonPivots). */
void nmod_poly_mat_echelon_pivot_index(slong * pivind,
                                       const nmod_poly_mat_t mat,
                                       orientation_t orient);

/** Computes the echelon pivot profile `pivind`, `pivdeg` of a polynomial
 * matrix `mat` (see @ref EchelonPivots). */
void nmod_poly_mat_echelon_pivot_profile(slong * pivind,
                                         slong * pivdeg,
                                         const nmod_poly_mat_t mat,
                                         orientation_t orient);

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
 */
//@{

/** Computes the degree matrix `degmat` of a polynomial matrix `pmat`
 * (see @ref DegreeMatrix).
 */
void nmod_poly_mat_degree_matrix(fmpz_mat_t dmat,
                                 const nmod_poly_mat_t mat);

/** Computes the `shift`-degree matrix `degmat` of a polynomial matrix `pmat`
 * (see @ref DegreeMatrix).
 */
void nmod_poly_mat_degree_matrix_shifted(fmpz_mat_t dmat,
                                         const nmod_poly_mat_t mat,
                                         const slong * shift,
                                         orientation_t orient);

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
 * The functions below which involve a `shift` among its parameters do not
 * check whether `shift` has the right length. Most functions accept that
 * `NULL` is provided as input for the shift, this is understood as the uniform
 * shift `[0,...,0]` of the right length.
 *
 * \todo minor enhancement: offer row-wise (resp column-wise) leading matrix
 * when the row degree (resp. column degree) is already known
 */
//@{

/** Computes the `shift`-leading matrix `lmat` of a polynomial matrix `mat`
 * (see @ref LeadingMatrix).  */
void nmod_poly_mat_leading_matrix(nmod_mat_t lmat,
                                  const nmod_poly_mat_t mat,
                                  const slong * shift,
                                  orientation_t orient);

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
 * A polynomial matrix is said to be in lower row echelon form (ROW_LOWER) if:
 * - there is no zero row below a nonzero row (possible zero rows are grouped
 *   at the top),
 * - let r be the number of nonzero rows (thus at indices m-r,...,m-1) and j_i
 *   be the index of the rightmost nonzero entry in row m-r+i for 0 <= i < r,
 *   then j_0 < ... < j_{r-1}.
 * ==> equivalently: the row-wise lower echelon pivot profile of the matrix is
 * increasing, excepting zero rows
 *
 * A polynomial matrix is said to be in upper row echelon form (ROW_UPPER) if:
 * - there is no zero row above a nonzero row (possible zero rows are grouped
 *   at the bottom),
 * - let r be the number of nonzero rows (thus at indices 0,...,r-1) and let
 *   j_i be the index of the rightmost nonzero entry in row i for 0 <= i < r,
 *   then j_0 < ... < j_{r-1}.
 * ==> equivalently: the row-wise upper echelon pivot profile of the matrix is
 * increasing, excepting zero rows
 * 
 * Definitions for COL_LOWER and COL_UPPER echelon forms are analogous.
 *
 * A polynomial matrix is said to be in ROW_LOWER Hermite normal form if:
 * - it is in ROW_LOWER echelon form,
 * - it has no zero row,
 * - all entries below a pivot entry have degree strictly less than the degree
 *   of that pivot entry.
 * Definitions are similar for other orientations.
 *
 * The functions below which involve a `shift` among its parameters do not
 * check whether `shift` has the right length. Most functions accept that
 * `NULL` is provided as input for the shift, this is understood as the uniform
 * shift `[0,...,0]` of the right length.
 *
 */
//@{

/** Tests whether `pmat` is in `shift`-reduced form (see @ref MatrixForms).  */
int nmod_poly_mat_is_reduced(const nmod_poly_mat_t mat,
                             const slong * shift,
                             orientation_t orient);

/** Tests whether `mat` is in `shift`-ordered weak Popov form (see
 * @ref MatrixForms). */
int nmod_poly_mat_is_ordered_weak_popov(const nmod_poly_mat_t mat,
                                        const slong * shift,
                                        orientation_t orient);

/** Tests whether `mat` is in `shift`-weak Popov form (see @ref MatrixForms) */
int nmod_poly_mat_is_weak_popov(const nmod_poly_mat_t mat,
                                const slong * shift,
                                orientation_t orient);

/** Tests whether `mat` is in `shift`-Popov form (see @ref MatrixForms) */
int nmod_poly_mat_is_popov(const nmod_poly_mat_t mat,
                           const slong * shift,
                           orientation_t orient);

/** Tests whether `mat` is in echelon form (see @ref MatrixForms) */
int nmod_poly_mat_is_echelon(const nmod_poly_mat_t pmat,
                             orientation_t orient);

/** Tests whether `mat` is in Hermite form (see @ref MatrixForms) */
int nmod_poly_mat_is_hermite(const nmod_poly_mat_t pmat,
                             orientation_t orient);

//@} // doxygen group: Testing polynomial matrix forms


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* COMPUTING ECHELON FORMS                                    */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/** @name Computing echelon forms
 * \anchor EchelonForms
 *
 * See @ref MatrixForms for definitions
 *
 * Transforms ``mat`` in place to a row-wise, upper row echelon form, and
 * returns the rank of ``mat``.
 *
 *   ..``tsf``. If ``tsf`` is not NULL, the same unimodular left-operations
 *   applied to ``mat`` are performed on ``tsf`` (which must therefore have as
 *   many rows as ``mat``, and can have an arbitrary number of columns).
 *   Setting ``tsf`` to the identity beforehand allows one to recover the
 *   unimodular transformation between ``mat`` and the computed Hermite normal
 *   form. ``mat`` cannot alias ``tsf``.
 *
 *   ..``pivind``. It is filled with the (upper echelon, row-wise) pivot
 *   indices of the output echelon form, which also correspond to the column
 *   rank profile of ``mat``. As input, ``pivind`` must have allocated space
 *   for at least rank(mat) entries (take min(mat->r,mat->c) if no better bound
 *   is known). No need to fill it with values. Its allocated space is left
 *   unchanged, and so are its entries beyond the rank(mat)-th one.
 *
 * Some algorithms may also provide:
 *
 *   ..``mrp``. It is filled with the matrix rank profile of the input matrix.
 *   See e.g. Dumas, Pernet, Sultan, Journal of Symbolic Computation 2017. It
 *   is represented here as an array of mat->r entries, the i-th entry is
 *   either some j in 0:mat->c giving the position (i,j) of `1` in the matrix
 *   rank profile, or it is -1 to indicate that row i does not contribute to
 *   the matrix rank profile. As input, ``mrp`` must have allocated space
 *   for at least mat->r entries; no need to fill them with values. Its
 *   allocated space is left unchanged, and so are its entries beyond the
 *   mat->r -th one.
 *
 *   .. ``rrp``. It is filled with the row rank profile of the input matrix.
 *   As input, ``rrp`` must have allocated space for at least rank(mat) entries
 *   (take min(mat->r,mat->c) if no better bound is known). No need to fill it
 *   with values. Its allocated space is left unchanged, and so are its entries
 *   beyond the rank(mat)-th one.
 *
 * \todo handle more orientations
 */
//@{

/** Upper row echelon form inspired by Rosser's HNF algorithm */
slong nmod_poly_mat_uref_maxdeg_atomic(nmod_poly_mat_t mat,
                                       nmod_poly_mat_t tsf,
                                       slong * pivind);

/** Upper row echelon form inspired by Bradley's HNF algorithm */
slong nmod_poly_mat_uref_revlex_xgcd(nmod_poly_mat_t mat,
                                     nmod_poly_mat_t tsf,
                                     slong * pivind,
                                     slong * mrp);

/** Upper row echelon form using lexicographic pivot search */
slong nmod_poly_mat_uref_lex_xgcd(nmod_poly_mat_t mat,
                                  nmod_poly_mat_t tsf,
                                  slong * pivind,
                                  slong * mrp);

/** Upper row echelon form inspired by Mulders and Storjohann's
 * HNF algorithm.
 * 
 * If not NULL, ``udet`` is either left the same or negated, according to the
 * determinant of the applied unimodular transformation (which transforms
 * the input into the output). The latter determinant is +1 or -1.
 *
 * \todo row rank profile rrp not supported yet
 * 
 * \todo to be implemented: an early detection based on known determinantal
 * degree (or sum of pivot degree) could be added so that as soon as only
 * trivial pivot entries remain to be found, the algorithm stops the iteration
 * over the leading principal minors and rather uses a simple constant
 * transformation to complete the computation.
 **/
slong nmod_poly_mat_uref_matrixgcd_iter(nmod_poly_mat_t mat,
                                        nmod_poly_mat_t tsf,
                                        slong * pivind,
                                        slong * rrp,
                                        slong * udet);


//@} // doxygen group: Computing echelon forms

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* COMPUTING HERMITE NORMAL FORM                              */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/** @name Computing Hermite normal form
 * See @ref MatrixForms for definitions
 *
 * Transforms ``mat`` in place to its row-wise, upper Hermite normal form, and
 * returns the rank of ``mat``.
 *
 * Input parameters are the same as for echelon forms, see @ref EchelonForms.
 *
 * \todo handle more orientations
 *
 * \todo eventually restrict the output HNF to rank rows? currently,
 * output HNF may contain zero rows
 *
 **/
//@{

/** Upper row echelon form normalization into upper row-wise HNF
 * \todo improve documentation
 **/
void _normalize_uref(nmod_poly_mat_t mat,
                     nmod_poly_mat_t other,
                     slong * pivind,
                     slong rk);

/** Hermite normal form in the style of Kannan-Bachem's algorithm
 * (uref with reverse lexicographic pivot search, with modified scheduling for
 * cancellation of non-pivot entries, and continuous normalization into HNF) */
slong nmod_poly_mat_hnf_ur_revlex_xgcd_delayed_zero(nmod_poly_mat_t mat,
                                                    nmod_poly_mat_t tsf, slong *
                                                    pivind, slong * mrp);

/** Hermite normal form in the style of Rosser's algorithm */
NMOD_POLY_MAT_INLINE
slong nmod_poly_mat_hnf_ur_maxdeg_atomic(nmod_poly_mat_t mat,
                                         nmod_poly_mat_t tsf,
                                         slong * pivind)
{
    // upper row echelon form + normalization
    slong rk = nmod_poly_mat_uref_maxdeg_atomic(mat, tsf, pivind);
    _normalize_uref(mat, tsf, pivind, rk);
    return rk;
}

/** Hermite normal form in the style of Bradley's algorithm */
NMOD_POLY_MAT_INLINE
slong nmod_poly_mat_hnf_ur_revlex_xgcd(nmod_poly_mat_t mat,
                                       nmod_poly_mat_t tsf,
                                       slong * pivind,
                                       slong * mrp)
{
    // upper row echelon form + normalization
    slong rk = nmod_poly_mat_uref_revlex_xgcd(mat, tsf, pivind, mrp);
    _normalize_uref(mat, tsf, pivind, rk);
    return rk;
}

/** Hermite normal form using lexicographic pivot search */
NMOD_POLY_MAT_INLINE
slong nmod_poly_mat_hnf_ur_lex_xgcd(nmod_poly_mat_t mat,
                                    nmod_poly_mat_t tsf,
                                    slong * pivind,
                                    slong * mrp)
{
    // upper row echelon form + normalization
    slong rk = nmod_poly_mat_uref_lex_xgcd(mat, tsf, pivind, mrp);
    _normalize_uref(mat, tsf, pivind, rk);
    return rk;
}

/** Hermite normal form in the style of Mulders&Storjohann's algorithm
 * 
 * If not NULL, ``udet`` is either left the same or negated, according to the
 * determinant of the applied unimodular transformation (which transforms
 * the input into the output).
 *
 * \todo rrp not well supported yet
 * \todo udet not supported: has to be updated during normalization
 * \todo should test for rk < 0 (non generic CRP) and proceed accordingly
 */
NMOD_POLY_MAT_INLINE
slong nmod_poly_mat_hnf_ur_matrixgcd_iter(nmod_poly_mat_t mat,
                                          nmod_poly_mat_t tsf,
                                          slong * pivind,
                                          slong * rrp,
                                          slong * udet)
{
    // upper row echelon form + normalization
    slong rk = nmod_poly_mat_uref_matrixgcd_iter(mat, tsf, pivind, rrp, udet);
    _normalize_uref(mat, tsf, pivind, rk);
    return rk;
}

//@} // doxygen group: Computing Hermite normal form

// TODO mod det version (see e.g. Domich) ? quite often for nmod_poly_mat's,
// computing the determinant is not really easier than computing the HNF...
// --> YET once may imagine e.g. approximant/interpolant basis computation,
// where the determinant (or a multiple of it) is actually known, so if having
// it helps, why not? (still, for approximants/interpolants, one might as well
// compute them directly in HNF...)
// --> ALSO for algorithms easily split into ref+normalization, the second step
// could be done easily done mod det (of the pivot part if rank deficient), which
// should be interesting for algorithms having large degree growth (but one may
// want to avoid these anyway?)


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* COMPUTING WEAK POPOV FORM                                  */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/** @name Computing weak Popov form
 * See @ref MatrixForms for definitions
 */
//@{

/** Mulders&Storjohann's weak Popov form algorithm, working row-by-row.
 *
 * This implements the iterative weak Popov form algorithm by Mulders and
 * Storjohann, 2003 (as in Figure 3, algo "RankProfile"). There presented for
 * the uniform shift, here straightforwardly adapted to the shifted case.
 *
 * The orientation `orient` must be among ROW_UPPER and ROW_LOWER; there
 * is no check that this is the case. Returns a non-ordered weak Popov form,
 * with zero rows not removed; possible zero rows are all grouped at the bottom
 * of the output matrix.
 *
 * .. ``shift``: shift for the shifted weak Popov form. The fact that the
 * length is valid is not tested.
 *
 * ..``tsf``. If ``tsf`` is not NULL, the same unimodular left-operations
 * applied to ``mat`` are performed on ``tsf`` (which must therefore have as
 * many rows as ``mat``, and can have an arbitrary number of columns).  Setting
 * ``tsf`` to the identity beforehand allows one to recover the unimodular
 * transformation between ``mat`` and the computed shifted weak Popov form.
 * ``mat`` cannot alias ``tsf``.
 *
 * .. ``det``: it can be NULL in which case it is ignored. Otherwise, it will
 * be multiplied with the determinant of the unimodular transformation used
 * during this call, which is +1 or -1.
 *
 * .. ``pivind``: integer array, must be allocated with at least mat->r
 * entries; it will be populated with the shifted pivot index of the output
 * weak Popov form (undefined behaviour for its entries beyond mat->r)
 *
 * .. ``rrp``. It is filled with the row rank profile of the input matrix.
 * As input, ``rrp`` must have allocated space for at least rank(mat) entries
 * (take min(mat->r,mat->c) if no better bound is known). No need to fill it
 * with values. Its allocated space is left unchanged, and so are its entries
 * beyond the rank(mat)-th one.
 *
 * .. ``early_exit_zr``: stop the computation as soon as early_exit_zr zero
 * rows have been found; in that case returns -rk (negative or zero and with rk
 * not necessarily related to the actual rank). If not interested in early
 * exit, put mat->r (or more). If the output is < 0, the output guarantees are
 * the same but only for the first |rk| + max_zr rows of mat (zero rows have
 * been put at bottom, etc)
 *
 * .. ``rstart``, ``cstart``, ``rdim``, ``cdim``: apply the algorithm to the
 * submatrix mat[rstart:rstart+rdim,cstart:cstart+cdim], called submat below.
 * The left unimodular transformations are applied to the whole of
 * mat[rstart:rstart+rdim,:], i.e. not restricting the columns to those of
 * submat. Some more details:
 *   ... mat must have >= rstart+rdim rows, >= cstart+cdim columns
 *   ... shift must have length >= cdim, its first cdim entries will be used as
 *   the shift
 *   ... tsf must be NULL or a matrix with at least rstart+rdim rows, the
 *   transformation will be applied to its rows rstart:rstart+rdim
 *   ... pivind must be allocated with at least rdim entries; its first rdim
 *   entries will be populated with the shifted pivot index of the output weak
 *   Popov form (undefined behaviour for entries beyond rdim)
 *   ... rrp  must be NULL or allocated with >= rank(submat) entries, it will
 *   eventually contain the row rank profile of submat as its first
 *   ``rank(submat)`` entries
 *
 * .. ``orient``. The orientation, among ROW_LOWER and ROW_UPPER. Note that for
 * ROW_LOWER, zero rows are put at the bottom even though they have the
 * smallest possible pivot index.
 *
 */
slong _nmod_poly_mat_weak_popov_iter_submat_rowbyrow(nmod_poly_mat_t mat,
                                                        const slong * shift,
                                                        nmod_poly_mat_t tsf,
                                                        slong * det,
                                                        slong * pivind,
                                                        slong * rrp,
                                                        slong rstart,
                                                        slong cstart,
                                                        slong rdim,
                                                        slong cdim,
                                                        slong early_exit_zr,
                                                        orientation_t orient);

/** Transforms ``mat`` in place to a row-wise, lower ``shift``-weak Popov form
 * (non-necessarily ordered, zero rows at the bottom), and returns the rank of
 * ``mat``.
 *    .. ``tsf``. If ``tsf`` is not NULL, the same unimodular left-operations
 *    applied to ``mat`` are performed on ``tsf`` (which must therefore have as
 *    many rows as ``mat``, and can have an arbitrary number of columns).
 *    ``mat`` cannot alias ``tsf``.  Setting ``tsf`` to the identity beforehand
 *    allows one to recover the unimodular transformation between ``mat`` and
 *    the computed weak Popov form.
 *
 *    ..``shifts``. If ``shift`` is NULL, the uniform shift is used.
 *
 *    ..``pivind``. It will eventually contain the ``shift``-pivot index of the
 *    output ``mat``. As input it must have allocated space for at least
 *    ``rank(mat)+1`` entries, except if mat has full row rank in case
 *    ``mat->r`` entries suffice (so, allocating ``mat->r`` entries is always
 *    ok, but ``min(mat->r,mat->c)`` may not be when mat->r > mat->c). No need
 *    to fill it with values. Its allocated space is left unchanged, and so are
 *    its entries beyond the rank(mat)-th one.
 *
 *    ..``rrp``. It must be either NULL or allocated with at least
 *    ``rank(mat)`` entries (no need to fill it with values). If not NULL,
 *    ``rrp`` will eventually contain the row rank profile of the input ``mat``
 *    as its first ``rank(mat)`` entries. Its allocated space is left
 *    unchanged, and so are its entries beyond the rank(mat)-th one.
 *
 *    .. ``orient``. The orientation, among ROW_LOWER and ROW_UPPER.
 *    See @ref orientation.
 *
 *  \todo support COL_LOWER, COL_UPPER
 *  \todo currently only applies rowbyrow; try other strategies?
 **/
NMOD_POLY_MAT_INLINE slong
nmod_poly_mat_weak_popov_iter(nmod_poly_mat_t mat,
                              const slong * shift,
                              nmod_poly_mat_t tsf,
                              slong * pivind,
                              slong * rrp,
                              orientation_t orient)
{
    return _nmod_poly_mat_weak_popov_iter_submat_rowbyrow(mat, shift, tsf, NULL, pivind, rrp, 0, 0, mat->r, mat->c, mat->r, orient);
}

/** Transforms ``mat`` in place to a row-wise, lower, ordered ``shift``-weak
 * Popov form (zero rows at the bottom), and returns the rank of ``mat``.
 * See the documentation of @ref nmod_poly_mat_weak_popov_iter for more
 * details.
 **/
slong nmod_poly_mat_ordered_weak_popov_iter(nmod_poly_mat_t mat,
                                            const slong * shift,
                                            nmod_poly_mat_t tsf,
                                            slong * pivind,
                                            slong * rrp,
                                            orientation_t orient);

//@} // doxygen group: Computing weak Popov form

#ifdef __cplusplus
}
#endif

#endif // NMOD_POLY_MAT_FORMS_H

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
