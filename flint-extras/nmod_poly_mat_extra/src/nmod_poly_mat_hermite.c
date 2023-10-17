#include <flint/flint.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_vec.h>
#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_forms.h"
#include "nmod_poly_mat_io.h" // TODO remove, for debugging

//#define MAT(i,j) nmod_poly_mat_entry(mat, i, j)
//#define TSF(i,j) nmod_poly_mat_entry(tsf, i, j)
#define MAT(i,j) (mat->rows[i] + j)
#define TSF(i,j) (tsf->rows[i] + j)
#define OTHER(i,j) (other->rows[i] + j)

// TODO how to handle upper/lower?
// TODO how to handle rowwise/colwise?
// TODO eventually restrict the output HNF to rank rows, where to do this?
// TODO check that all memory is cleared everywhere
//
//
//
// For all:
// - upper Hermite, row-wise
// - computation is in place
// - concerning transformation:
//    -- if tsf == NULL, transformation is not computed at all (which saves computation time)
//    -- otherwise, the transformation is accumulated into tsf
//         --> if wanting the transformation for this computation only, set tsf to identity before calling this

/**********************************************************************
*                          HELPER FUNCTIONS                          *
**********************************************************************/

// rotate rows of mat from i to j (requirement: 0 <= i <= j < mat->r)
// and apply the corresponding transformation to vec (requirement: j < len(vec))
// If i == j, then nothing happens.
// vec can be NULL, in case it is omitted
// More precisely this performs simultaneously:
//      mat[i,:]     <--    mat[j,:]
//      mat[i+1,:]   <--    mat[i,:]
//      mat[i+2,:]   <--    mat[i+1,:]
//         ...       <--       ...
//      mat[j-1,:]   <--    mat[j-2,:]
//      mat[j,:]     <--    mat[j-1,:]
// as well as
//      vec[i]     <--    vec[j]
//      vec[i+1]   <--    vec[i]
//      vec[i+2]   <--    vec[i+1]
//        ...      <--      ...
//      vec[j-1]   <--    vec[j-2]
//      vec[j]     <--    vec[j-1]
void _nmod_poly_mat_rotate_rows(nmod_poly_mat_t mat, slong * vec, slong i, slong j)
{
    if (i != j)
    {
        if (vec)
        {
            slong tmp_vec = vec[j];
            for (slong ii = j; ii > i; ii--)
                vec[ii] = vec[ii-1];
            vec[i] = tmp_vec;
        }

        nmod_poly_struct * tmp_mat = mat->rows[j];
        for (slong ii = j; ii > i; ii--)
            mat->rows[ii] = mat->rows[ii-1];
        mat->rows[i] = tmp_mat;
    }
}

// Context: upper echelon, row-wise
// two potential new HNF pivots have been found in the same column, i.e.
// entries [pi1,j] and [pi2,j] which are both nonzero, with deg [pi1,j] >= deg
// [pi2,j], and all entries [pi1,jj] and [pi2,jj] for jj < j are zero (pivots
// in columns 0:j have been found and none of them is in row pi1 or pi2)
// -> use pi2,j to kill the leading term of pi1,j, and apply the
// corresponding transformation to the whole row pi1
// --> this is a "simple transformation of the first kind" (terminology from Mulders & Storjohann 2003)
//         mat[pi1,:] = mat[pi1,:] + cst * x**exp * mat[pi2,:]
// where cst = - leading_coeff(mat[pi1,j]) / leading_coeff(mat[pi2,j])
//   and exp = deg(mat[pi1,j]) - deg(mat[pi2,j])
void _apply_pivot_collision_transformation_upper_triangular_rowwise(nmod_poly_mat_t mat, nmod_poly_mat_t other,
                                               slong pi1, slong pi2, slong j)
{
    // compute exp = deg(mat[pi1,j]) - deg(mat[pi2,j])
    const slong exp = MAT(pi1, j)->length - MAT(pi2, j)->length;
    // compute cst = - leading_coeff(mat[pi1,j]) / leading_coeff(mat[pi2,j])
    mp_limb_t cst = n_invmod(MAT(pi2, j)->coeffs[MAT(pi2, j)->length -1], mat->modulus);
    cst = nmod_mul(MAT(pi1, j)->coeffs[MAT(pi1, j)->length -1], cst, MAT(pi1, j)->mod);
    cst = nmod_neg(cst, MAT(pi1, j)->mod);
    // update pivot pi1,j
    _nmod_vec_scalar_addmul_nmod(MAT(pi1, j)->coeffs+exp, MAT(pi2, j)->coeffs, MAT(pi2, j)->length, cst, MAT(pi1, j)->mod);
    _nmod_poly_normalise(MAT(pi1, j));
    // update entries to the right of pi1,j (left ones are zero)
    for (slong jj = j+1; jj < mat->c; jj++)
    {
        const slong len = MAT(pi2, jj)->length + exp;
        nmod_poly_fit_length(MAT(pi1, jj), len);
        if (len > MAT(pi1, jj)->length)
        {
            _nmod_vec_zero(MAT(pi1, jj)->coeffs + MAT(pi1, jj)->length, len - MAT(pi1, jj)->length);
            _nmod_poly_set_length(MAT(pi1, jj), len);
        }
        _nmod_vec_scalar_addmul_nmod(MAT(pi1, jj)->coeffs+exp, MAT(pi2, jj)->coeffs, MAT(pi2, jj)->length, cst, MAT(pi1, jj)->mod);
        _nmod_poly_normalise(MAT(pi1, jj));
    }
    // update the row pi1 of other
    if (other) for (slong jj = 0; jj < other->c; jj++)
    {
        const slong len = OTHER(pi2, jj)->length + exp;
        nmod_poly_fit_length(OTHER(pi1, jj), len);
        if (len > OTHER(pi1, jj)->length)
        {
            _nmod_vec_zero(OTHER(pi1, jj)->coeffs + OTHER(pi1, jj)->length, len - OTHER(pi1, jj)->length);
            _nmod_poly_set_length(OTHER(pi1, jj), len);
        }
        _nmod_vec_scalar_addmul_nmod(OTHER(pi1, jj)->coeffs+exp, OTHER(pi2, jj)->coeffs, OTHER(pi2, jj)->length, cst, OTHER(pi1, jj)->mod);
        _nmod_poly_normalise(OTHER(pi1, jj));
    }
}

// Context: upper echelon, row-wise
// computing xgcd and applying corresponding unimodular transformation
// input is mat, other, pi, ii, j
// pivot is at pi,j
// other entry is at ii,j
// we perform xgcd between pivot and mat[ii,j]
// g, u, v, pivg, nonzg are used as temporaries and must be already initialized
void _apply_xgcd_transformation(nmod_poly_mat_t mat, nmod_poly_mat_t other,
                                slong pi, slong ii, slong j,
                                nmod_poly_t g, nmod_poly_t u, nmod_poly_t v,
                                nmod_poly_t pivg, nmod_poly_t nonzg)
{
    // FIXME this specific code for pivot == 1 does not seem to speed up things
    if (nmod_poly_is_one(MAT(pi, j)))
    {
        // cancel (ii,j) and modify row ii accordingly
        // in short: for M in {mat,tsf} do the unimodular transformation
        //   M[ii,:] = M[ii,:] - mat[ii,j] * M[pi,:]
        // make entry mat[ii,j] zero at the very end, since it is used throughout
        // use g as temporary all along
        for (slong jj = j+1; jj < mat->c; jj++)
        {
            nmod_poly_mul(g, MAT(ii, j), MAT(pi, jj));
            nmod_poly_sub(MAT(ii, jj), MAT(ii, jj), g);
        }
        if (other)
        {
            for (slong jj = 0; jj < other->c; jj++)
            {
                nmod_poly_mul(g, MAT(ii, j), OTHER(pi, jj));
                nmod_poly_sub(OTHER(ii, jj), OTHER(ii, jj), g);
            }
        }
        nmod_poly_zero(MAT(ii, j));
    }
    else // general case, pivot is not one
    {
        nmod_poly_xgcd(g, u, v, MAT(pi, j), MAT(ii, j)); // g = u*pivot + v*other
        nmod_poly_divides(pivg, MAT(pi, j), g); // pivg = pivot // g
        nmod_poly_divides(nonzg, MAT(ii, j), g); // nonzg = other // g

        // cancel (ii,j) and modify row ii accordingly
        // in short: for M in {mat,other} do the unimodular transformation on rows pi,ii:
        //  [ M[pi,:] ]  =  [    u       v  ]  *  [ M[pi,:] ]
        //  [ M[ii,:] ]     [ -nonzg   pivg ]     [ M[ii,:] ]

        // for mat, due to zeroes already set in columns 0...j-1, we can start at column j
        // set pivot to g, and set other entry to 0 (done at the very end since used as temporary)
        nmod_poly_set(MAT(pi, j), g); // mat[pi,j] = g
        for (slong jj = j+1; jj < mat->c; jj++)
        {
            // simultaneously update:
            //     mat[pi,jj] = u * mat[pi,jj] + v * mat[ii,jj]
            //     mat[ii,jj] = -nonzg * mat[pi,jj] + pivg * mat[ii,jj]
            // --> use g as temporary copy of mat[pi,jj]
            // --> use entry mat[ii,j] as temporary for storing products
            nmod_poly_set(g, MAT(pi, jj));

            nmod_poly_mul(MAT(pi, jj), u, g);
            nmod_poly_mul(MAT(ii, j), v, MAT(ii, jj));
            nmod_poly_add(MAT(pi, jj), MAT(pi, jj), MAT(ii, j));

            nmod_poly_mul(MAT(ii, j), pivg, MAT(ii, jj));
            nmod_poly_mul(MAT(ii, jj), nonzg, g);
            nmod_poly_sub(MAT(ii, jj), MAT(ii, j), MAT(ii, jj));
        }
        if (other)
        {
            for (slong jj = 0; jj < mat->r; jj++)
            {
                // simultaneously update:
                //     other[pi,jj] = u * other[pi,jj] + v * other[ii,jj]
                //     other[ii,jj] = -nonzg * other[pi,jj] + pivg * other[ii,jj]
                // --> use g as temporary copy of other[pi,jj]
                // --> use entry mat[ii,j] as temporary for storing products
                nmod_poly_set(g, OTHER(pi, jj));

                nmod_poly_mul(OTHER(pi, jj), u, g);
                nmod_poly_mul(MAT(ii, j), v, OTHER(ii, jj));
                nmod_poly_add(OTHER(pi, jj), OTHER(pi, jj), MAT(ii, j));

                nmod_poly_mul(MAT(ii, j), pivg, OTHER(ii, jj));
                nmod_poly_mul(OTHER(ii, jj), nonzg, g);
                nmod_poly_sub(OTHER(ii, jj), MAT(ii, j), OTHER(ii, jj));
            }
        }
        nmod_poly_zero(MAT(ii, j));
    }
}

// Context: upper echelon, row-wise
// normalize pivot entry in given row of mat, applying the corresponding
// operation to the whole row of mat and the corresponding row of other
// (i,j) is position of the pivot entry
void _normalize_pivot_upper_triangular_rowwise(nmod_poly_mat_t mat, nmod_poly_mat_t other, slong i, slong j)
{
    if (! nmod_poly_is_monic(MAT(i, j)))
    {
        mp_limb_t inv = n_invmod(MAT(i, j)->coeffs[MAT(i, j)->length - 1], MAT(i, j)->mod.n);
        for (slong jj = j; jj < mat->c; jj++)
            _nmod_vec_scalar_mul_nmod(MAT(i, jj)->coeffs, MAT(i, jj)->coeffs, MAT(i, jj)->length, inv, MAT(i, jj)->mod);
        if (other)
            for (slong jj = 0; jj < other->c; jj++)
                _nmod_vec_scalar_mul_nmod(OTHER(i, jj)->coeffs, OTHER(i, jj)->coeffs, OTHER(i, jj)->length, inv, OTHER(i, jj)->mod);
    }
}

// Context: upper echelon, row-wise
// reduce entry mat[ii,j] against pivot entry which is at mat[i,j],
// applying the corresponding operation to the whole row ii of mat and the
// corresponding row ii of other
// u, v are used as temporaries and must be already initialized
void _reduce_against_pivot(nmod_poly_mat_t mat, nmod_poly_mat_t other,
                           slong i, slong j, slong ii,
                           nmod_poly_t u, nmod_poly_t v)
{
    if (MAT(ii, j)->length >= MAT(i, j)->length)
    {
        // division with remainder: [ii,j] = pivot * u + v
        nmod_poly_divrem(u, v, MAT(ii, j), MAT(i, j));
        // replace [ii,j] by remainder
        nmod_poly_swap(MAT(ii, j), v);
        // apply same transformation on rest of the row: row_ii <- row_ii - u * row_i
        nmod_poly_neg(u, u);
        for (slong jj = j+1; jj < mat->c; jj++)
        {
            nmod_poly_mul(v, u, MAT(i, jj));
            nmod_poly_add(MAT(ii, jj), MAT(ii, jj), v);
        }
        if (other)
        {
            for (slong jj = 0; jj < mat->r; jj++)
            {
                nmod_poly_mul(v, u, OTHER(i, jj));
                nmod_poly_add(OTHER(ii, jj), OTHER(ii, jj), v);
            }
        }
    }
}

/**********************************************************************
*                   Hermite normal form algorithms                   *
**********************************************************************/
// Orientation: upper echelon, row-wise


// Rosser's algorithm
// proceed column by column. When looking for i-th pivot (in column j >= i) look
// at entries [i:m,j], and kill leading term of largest deg by using the second
// largest deg (using basic operation only, adding a multiple of another row by
// constant*x**k); continue in the same column until all entries in [i:m,j]
// except one are zero; swap rows to put the nonzero one at [i,j]
// --> TODO try version using divrem and poly multiplication instead of constant*x**k?
// --> FIXME normalization step could be inside inside the main triangularization loop
// (applying it each time we find a pivot), but this probably has no impact
slong nmod_poly_mat_upper_hermite_form_rowwise_rosser(nmod_poly_mat_t mat, nmod_poly_mat_t tsf)
{
    const slong m = mat->r;
    const slong n = mat->c;

    if (m == 0 || n == 0)
        return 0;

    // to store the positions of columns where a pivot was found
    // -> pivots_cols[i] gives index of pivot in row i
    // -> index rk below will tell us how many have been found already,
    //     pivots_cols[i] for i >= rk is undefined
    // -> there will never be > rank(mat) pivots
    slong * pivot_col = flint_malloc(FLINT_MIN(m,n) * sizeof(slong));

    // (rk,j) : position where next pivot is expected to be found when processing column j,
    // =>   rk = number of pivots already found in columns < j, which is also the rank
    // of these first j-1 columns
    // Note: in case of zeroes/dependencies, we may find a pivot below rk,
    // then we will permute rows to bring it up at row rk
    slong rk = 0;

    // loop over the columns
    for (slong j = 0; j < n; j++)
    {
        // whether there are still nonzero elements in entries rk+1:m,j
        int collision = 1;
        while (collision)
        {
            // find i1 and i2 such that entries (i1,j) and (i2,j) are those of largest degrees,
            // with the former greater than or equal to the former, among the entries (i,j) for i = rk ... m
            // first step: find first nonzero entry
            slong pi1 = rk;
            while (pi1 < m && nmod_poly_is_zero(MAT(pi1, j)))
                pi1++;
            if (pi1 == m) // there were only zeroes in rk:m,j
                collision = 0;
            else if (pi1 == m-1) // only nonzero entry in rk:m,j is pi1,j --> pivot found, j-th iteration finished
            {
                // put the nonzero entry as pivot entry, and exit column j
                nmod_poly_mat_swap_rows(mat, NULL, pi1, rk);
                if (tsf)
                    nmod_poly_mat_swap_rows(tsf, NULL, pi1, rk);
                pivot_col[rk] = j;
                rk++;
                collision = 0;
            }
            else // there is a collision, work on it
            {
                slong pi2 = pi1 + 1; // recall pi1+1 < m
                // find pi1,pi2 such that entries pi1,j and pi2,j have the largest degree
                // among entries rk:m,j, with deg(mat[pi1,j]) >= deg(mat[pi2,j])
                if (MAT(pi1, j)->length < MAT(pi2, j)->length)
                    SLONG_SWAP(pi1, pi2);
                for (slong i = FLINT_MAX(pi1,pi2)+1; i < m; i++)
                {
                    if (MAT(i, j)->length > MAT(pi1, j)->length)
                    {
                        pi2 = pi1;
                        pi1 = i;
                    }
                    else if (MAT(i, j)->length > MAT(pi2, j)->length)
                        pi2 = i;
                }
                if (nmod_poly_is_zero(MAT(pi2, j)))
                    // only nonzero entry in rk:m,j is pi1,j --> pivot found, j-th iteration finished
                {
                    nmod_poly_mat_swap_rows(mat, NULL, pi1, rk);
                    if (tsf)
                        nmod_poly_mat_swap_rows(tsf, NULL, pi1, rk);
                    pivot_col[rk] = j;
                    rk++;
                    collision = 0;
                }
                else
                    _apply_pivot_collision_transformation_upper_triangular_rowwise(mat, tsf, pi1, pi2, j);
            }
        }
    }

    // now, we have mat in upper echelon form, it remains to
    // 1. reduce the entries above each pivot
    // 2. normalize each pivot
    nmod_poly_t u;
    nmod_poly_t v;
    nmod_poly_init(u, mat->modulus);
    nmod_poly_init(v, mat->modulus);
    for (slong i = 0; i < rk; i++)
    {
        slong j = pivot_col[i];
        // normalize pivot, and correspondingly update the row i of mat and of tsf
        _normalize_pivot_upper_triangular_rowwise(mat, tsf, i, j);
        // reduce entries above (i,j)
        for (slong ii = 0; ii < i; ii++)
            _reduce_against_pivot(mat, tsf, i, j, ii, u, v);
    }

    flint_free(pivot_col);

    return rk;
}

// Bradley's algorithm
// proceed column by column. When looking for i-th pivot (in column j >= i) look
// at entries [i:m,j], and repeatedly use gcd's between the two first nonzero
// entries in [i:m,j], say at [pi,j] and [ii,j], to zero out [ii,j] and reduce [pi,j]
// as much as possible. If g = u * mat[pi,j] + v * mat[ii,j], then this means
// applying the following unimodular transformation to rows pi and ii:
//  [ mat[pi,:] ]  =  [    u       v  ]  *  [ mat[pi,:] ]
//  [ mat[ii,:] ]     [ -nonzg   pivg ]     [ mat[ii,:] ]
// where nonzg = mat[ii,j]/g   and   pivg = mat[pi,j]/g
// This goes on same column until all entries in [i:m,j] except one are zero;
// swap rows to put the nonzero one at [i,j]. Then proceed to next column.
slong nmod_poly_mat_upper_hermite_form_rowwise_bradley(nmod_poly_mat_t mat, nmod_poly_mat_t tsf)
{
    const slong m = mat->r;
    const slong n = mat->c;

    if (m == 0 || n == 0)
        return 0;

    // to store the positions of columns where a pivot was found
    // -> pivots_cols[i] gives index of pivot in row i
    // -> index rk below will tell us how many have been found already,
    //     pivots_cols[i] for i >= rk is undefined
    // -> there will never be > rank(mat) pivots
    slong * pivot_col = flint_malloc(FLINT_MIN(m,n) * sizeof(slong));

    // (rk,j) : position where next pivot is expected to be found when processing column j,
    // =>   rk = number of pivots already found in columns < j, which is also the rank
    // of these first j-1 columns
    // Note: in case of zeroes/dependencies, we may find a pivot below rk,
    // then we will permute rows to bring it up at row rk
    slong rk = 0;

    nmod_poly_t g; // gcd
    nmod_poly_t u; // gcd cofactor1
    nmod_poly_t v; // gcd cofactor2
    nmod_poly_t pivg; // pivot divided by gcd
    nmod_poly_t nonzg; // other nonzero entry divided by gcd
    nmod_poly_init(g, mat->modulus);
    nmod_poly_init(u, mat->modulus);
    nmod_poly_init(v, mat->modulus);
    nmod_poly_init(pivg, mat->modulus);
    nmod_poly_init(nonzg, mat->modulus);

    // loop over the columns
    for (slong j = 0; j < n; j++)
    {
        // find actual pivot: starting at (rk,j), find first nonzero entry in the rest of column j
        slong pi = rk;
        while (pi < m && nmod_poly_is_zero(MAT(pi, j)))
            pi++;

        if (pi < m)
        // else, no pivot in this column, go to next column
        {
            // process column j below pivot
            slong ii = pi + 1;
            while (ii < m)
            {
                // find next nonzero entry in column
                while (ii < m && nmod_poly_is_zero(MAT(ii, j)))
                    ii++;
                if (ii < m)
                    _apply_xgcd_transformation(mat, tsf, pi, ii, j, g, u, v, pivg, nonzg);
                // else, ii==m, pivot is the only nonzero entry: nothing more to do for this column
            }

            // if pivot found in a lower row than expected
            // -> swap rows to bring it up to row rk
            nmod_poly_mat_swap_rows(mat, NULL, rk, pi); // does nothing if rk==pi
            if (tsf)
                nmod_poly_mat_swap_rows(tsf, NULL, rk, pi);
            pivot_col[rk] = j;
            rk++;
        }
    }

    // now, we have mat in upper echelon form, it remains to
    // 1. normalize each pivot
    //     (this will do nothing for already monic pivots, which is the case
    //     for all those that were obtained through some xgcd)
    // 2. reduce the entries above each pivot
    for (slong i = 0; i < rk; i++)
    {
        slong j = pivot_col[i];
        // normalize pivot, and correspondingly update the row i of mat and of tsf
        _normalize_pivot_upper_triangular_rowwise(mat, tsf, i, j);
        // reduce entries above (i,j)
        for (slong ii = 0; ii < i; ii++)
            _reduce_against_pivot(mat, tsf, i, j, ii, u, v);
    }

    flint_free(pivot_col);
    nmod_poly_clear(g);
    nmod_poly_clear(u);
    nmod_poly_clear(v);
    nmod_poly_clear(pivg);
    nmod_poly_clear(nonzg);

    return rk;
}

// Kannan-Bachem's algorithm
// proceed by leading submatrices. When looking for i-th pivot, in column j >=
// i, the (i-1) x (j-1) leading submatrix `hnf` is already in Hermite form.
//
// Then, ensure entries [i,jj] for jj in 0:j are zero. Look for the first row
// ii in i:m such that [ii,0:j+1] is nonzero. If there is none increment j and
// process next submatrix, otherwise for jj in 0:j in that order, do:
// - if mat[ii,jj] is nonzero and there is no pivot in column jj in hnf, reduce
// entries [ii,jj+1:j] against hnf, reduce appropriate entries in column jj of
// hnf against mat[ii,jj], and move row ii to the correct location in hnf;
// increment the pivot count i without increasing j, and go to process the next
// leading submatrix
// - else, if mat[ii,jj] is nonzero and there is a pivot in column jj in hnf,
// say in row pi, use gcd between [pi,jj] and [ii,jj] and apply the
// corresponding row-wise unimodular transformation between rows pi and ii (see
// description of Bradley's algorithm above), which puts a zero at [ii,jj] and
// the gcd at [pi,jj]; reduce appropriate entries in hnf in column jj, and
// proceed to next jj.
//
// Now that all entries in [ii,0:j] are zero, check if entry [ii,j] is nonzero;
// - if yes, this is a new pivot, reduce entries above it, increment both i and
// j, and proceed to next leading submatrix
// - if no, swap row ii and the last row, and keep the same i,j to process same
// submatrix again.
slong nmod_poly_mat_upper_hermite_form_rowwise_kannan_bachem(nmod_poly_mat_t mat, nmod_poly_mat_t tsf)
{
    const slong m = mat->r;
    const slong n = mat->c;

    if (m == 0 || n == 0)
        return 0;

    // to store the positions of columns where a pivot was found
    // -> pivots_cols[ii] gives index of pivot in row ii
    // -> index i below will tell us how many have been found already,
    //     pivots_cols[ii] for ii >= i is undefined
    // -> there will never be > rank(mat) pivots
    slong * pivot_col = flint_malloc(FLINT_MIN(m,n) * sizeof(slong));
    // for simplicity also store reverse array
    // -> pivot_row[j] is either -1 (not among the pivots discovered so far)
    // or is the index of the row with pivot j
    slong * pivot_row = flint_malloc(n * sizeof(slong));
    flint_mpn_store(pivot_row, n, -1); // fill with -1

    // (i,j) : the (i-1) x (j-1) leading submatrix `hnf` is already in Hermite form,
    // we will now try to increase either i (if new pivot found in row i or
    // below), or j (if submatrix [i:m,0:j+1] does not bring any new pivot), or
    // both (if new pivot found in column j)
    slong i = 0;
    slong j = 0;

    nmod_poly_t g; // gcd
    nmod_poly_t u; // gcd cofactor1
    nmod_poly_t v; // gcd cofactor2
    nmod_poly_t pivg; // piv divided by gcd
    nmod_poly_t nonzg; // nonz divided by gcd
    nmod_poly_init(g, mat->modulus);
    nmod_poly_init(u, mat->modulus);
    nmod_poly_init(v, mat->modulus);
    nmod_poly_init(pivg, mat->modulus);
    nmod_poly_init(nonzg, mat->modulus);

    // stop when we have reached m pivots, i.e. we have processed all rows and
    // reached maximal rank, or when we have processed all columns
    while (i < m && j < n)
    {
        // Look for the first row ii in i:m such that [ii,0:j+1] is nonzero.
        // also find jj such that [ii,jj] is the first nonzero entry among [ii,0:j+1]
        slong ii = i;
        slong jj = 0;
        int row_is_zero = 1;
        while (ii < m && row_is_zero)
        {
            jj = 0;
            while (jj <= j && nmod_poly_is_zero(MAT(ii, jj)))
                jj++;
            if (jj == j+1) // is zero, jump to next row
                ii++;
            else // jj <= j and mat[ii,jj] is nonzero
                row_is_zero = 0;
        }
        if (ii == m) // no nonzero row, increment j and process next submatrix
            j++;
        else // found nonzero row
        {
            // variable now used to record whether we found new pivot, or only
            // put zeroes thanks to gcd transformations
            row_is_zero = 1;

            // swap rows i and ii
            nmod_poly_mat_swap_rows(mat, NULL, i, ii); // does nothing if i == ii
            if (tsf)
                nmod_poly_mat_swap_rows(tsf, NULL, i, ii);

            // process entries i,jj for jj in 0:j
            // (start at previously found jj since entries before that are zero)
            for (; jj < j; jj++)
            {
                if (! nmod_poly_is_zero(MAT(i, jj)))
                {
                    slong pi = pivot_row[jj];
                    if (pi >= 0) // a pivot has already been found at [pi,jj], for some pi < i
                    {
                        // use gcd between [pi,jj] and [i,jj] and apply the corresponding row-wise unimodular transformation
                        // between rows pi and i (see description of Bradley's algorithm), which puts a zero at [i,jj] and
                        // the gcd (updated pivot) at [pi,jj]
                        _apply_xgcd_transformation(mat, tsf, pi, i, jj, g, u, v, pivg, nonzg);
                        // note: new pivot already normalized since xgcd ensures the gcd is monic

                        // reduce appropriate entries in hnf in column jj
                        for (ii = 0; ii < pi; ii++)
                            _reduce_against_pivot(mat, tsf, pi, jj, ii, u, v);
                    }
                    else // pivot_row[jj] == -1, currently no pivot in column jj
                    {
                        row_is_zero = 0;
                        // move row i to the correct location in hnf and update pivot_row, pivot_col
                        ii = 0; // index where row i must be placed
                        while (ii < i && pivot_col[ii] < jj)
                            ii++;
                        // rotate mat and pivot_col
                        pivot_col[i] = jj;
                        _nmod_poly_mat_rotate_rows(mat, pivot_col, ii, i);
                        if (tsf)
                            _nmod_poly_mat_rotate_rows(tsf, NULL, ii, i);
                        pivot_row[jj] = i; // column jj is new pivot for row i
                        // update pivot row for found pivots in jj+1 ... j-1
                        for (slong pj = jj+1; pj < j; pj++)
                            if (pivot_row[pj] >= 0)
                                pivot_row[pj] += 1;

                        // the new row is now at index ii, with new pivot at ii,jj
                        // normalize the new pivot
                        _normalize_pivot_upper_triangular_rowwise(mat, tsf, ii, jj);

                        // reduce entries [ii,jj+1:j] against current hnf
                        for (slong kk = jj+1; kk < j; kk++)
                        {
                            pi = pivot_row[kk];
                            if (pi >= 0) // otherwise, no pivot, nothing to do
                            {
                                // reduce mat[ii,kk] against pivot [pi,kk]
                                _reduce_against_pivot(mat, tsf, pi, kk, ii, u, v);
                            }
                        }
                        // reduce column jj of hnf against new pivot mat[ii,jj]
                        for (slong pii = 0; pii < ii; pii++)
                        {
                            // reduce pii,jj against the pivot [ii,jj]
                            _reduce_against_pivot(mat, tsf, ii, jj, pii, u, v);
                        }
                        // increment the pivot count i without increasing j
                        i++;
                        // abort current jj-for loop, go to next submatrix
                        break;
                    }
                }
            }
            if (row_is_zero)
            {
                // now all entries in [i,0:j] have been set to zero via gcd
                // operations, and we have not detected a new pivot
                if (nmod_poly_is_zero(MAT(i, j)))
                {
                    // entry [i,j] is zero, swap row i and the last row (just
                    // to avoid wasting time re-checking it is zero), and keep
                    // the same i,j to process same submatrix again.
                    nmod_poly_mat_swap_rows(mat, NULL, i, mat->r -1);
                    if (tsf)
                        nmod_poly_mat_swap_rows(tsf, NULL, i, mat->r -1);
                }
                else
                {
                    // update pivot_col and pivot_row
                    pivot_col[i] = j;
                    pivot_row[j] = i;
                    // entry [i,j] is nonzero, this is a new pivot
                    // normalize it
                    _normalize_pivot_upper_triangular_rowwise(mat, tsf, i, j);
                    // reduce entries above it, i.e. ii,j against pivot i,j
                    for (ii = 0; ii < i; ii++)
                        _reduce_against_pivot(mat, tsf, i, j, ii, u, v);
                    // increment both i and j
                    i++;
                    j++;
                }
            }
        }
    }

    flint_free(pivot_col);
    flint_free(pivot_row);

    return i;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
