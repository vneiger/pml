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
//
//
//
// For all:
// - upper Hermite, row-wise
// - computation is in place
// - concerning transformation:
//    -- if tsf == NULL, transformation is not computed at all (which saves computation time)
//    -- otherwise, the transformation is accumulated into tsf
//         --> if wanting the unimodular transformation, set tsf to identity before calling this

/**********************************************************************
*                          HELPER FUNCTIONS                          *
**********************************************************************/

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
// Complexity (field operations):
//  - for mat: sum_{jj = j ... mat->c-1} (deg(mat[pi2,jj]) + 1)
//             <=  (mat->c - j) * (rdeg(mat[pi2,:])+1)
//  - for other: sum_{jj = 0 ... other->c-1} (deg(other[pi2,jj]) + 1)
//             <=  other->c * (rdeg(other[pi2,:])+1)
void _atomic_solve_pivot_collision_uechelon_rowwise(nmod_poly_mat_t mat, nmod_poly_mat_t other,
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

// TODO implement _single_solve (or whatever name), where a divrem is used
// instead of just one step of Euclidean division (as in _atomic_solve)

// Context: upper echelon, row-wise
// computing xgcd and applying corresponding unimodular transformation
// input is mat, other, pi, ii, j
// pivot is at pi,j
// other entry is at ii,j
// we perform complete xgcd transformation between pivot and mat[ii,j]
//   (see more detailed description in algorithm comments below)
// g, u, v, pivg, nonzg are used as temporaries and must be already initialized
// Complexity:
//    . xgcd in degrees deg(mat[pi,j]) and deg(mat[ii,j])
//        ~ M(d) log(d) where d is max of these
//    . dividing these entries by the gcd: bounded by M(d)
//    . sum_{jj = j+1 ... mat->c-1} M(d + deg(mat[pi,jj])) + M(d + deg(mat[ii,jj]))
//         <=   (mat->c - j - 1) * M(rdeg(mat[pi,:]) + rdeg(mat[ii,:]))
//    . for other:
//        sum_{jj = 0 ... other->c-1} M(d + deg(other[pi,jj])) + M(d + deg(other[ii,jj]))
//         <=   other->c * M(d + rdeg(other[pi,:]) + rdeg(other[ii,:]))
void _complete_solve_pivot_collision_uechelon_rowwise(nmod_poly_mat_t mat, nmod_poly_mat_t other,
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
            for (slong jj = 0; jj < other->c; jj++)
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
// reduce entry mat[ii,j] against pivot entry which is at mat[i,j],
// applying the corresponding operation to the whole row ii of mat and the
// corresponding row ii of other
// u, v are used as temporaries and must be already initialized
// other == NULL or other is another matrix with the same modulus (it must have a row ii)
// Complexity: 0 if deg(mat[ii,j]) < deg(mat[i,j]), otherwise:
//    . divrem: M(deg(mat[i,j]))
//    . for mat: sum_{jj = j+1 ... mat->c-1} M(deg(mat[ii,j]) - deg(mat[i,j]) + deg(mat[i,jj]))
//        <=     (mat->c - j - 1) * M(deg(mat[ii,j]) - deg(mat[i,j]) + rdeg(mat[i,:]))
//    . for other: sum_{jj = 0 ... other->c-1} M(deg(mat[ii,j]) - deg(mat[i,j]) + deg(other[i,jj]))
//        <=     other->c * M(deg(mat[ii,j]) - deg(mat[i,j]) + rdeg(other[i,:]))
void _reduce_against_pivot_uechelon_rowwise(nmod_poly_mat_t mat, nmod_poly_mat_t other,
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
        // apply same transformation on row ii of other
        if (other)
        {
            for (slong jj = 0; jj < other->c; jj++)
            {
                nmod_poly_mul(v, u, OTHER(i, jj));
                nmod_poly_add(OTHER(ii, jj), OTHER(ii, jj), v);
            }
        }
    }
}

// Context: upper echelon, row-wise
// normalize pivot entry in given row of mat, applying the corresponding
// operation to the whole row of mat and the corresponding row of other
// (i,j) is position of the pivot entry
// other == NULL or other is another matrix with the same modulus
// Complexity: if already monic, 0, otherwise:
//    . for mat: sum_{jj = j ... mat->c-1} deg(mat[i,jj])
//          <= (mat->c - j) * rdeg(mat[i,:])
//    . for other: sum_{jj = 0 ... other->c-1} deg(other[i,jj])
//          <= other->c * rdeg(other[i,:])
void _normalize_pivot_uechelon_rowwise(nmod_poly_mat_t mat, nmod_poly_mat_t other, slong i, slong j)
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

/**********************************************************************
*                   Hermite normal form algorithms                   *
**********************************************************************/
// Orientation: upper echelon, row-wise


// Rosser's HNF algorithm   (upper echelon, row-wise)
// Pivoting strategy: reverse lexicographic
// Transformations: atomic
//
// In short: pivoting strategy is "reverse lexicographic" (search for the topmost
// nonzero entry in the first nonzero column); use atomic transformations
// for handling pivot collisions. Typically offers good control of the degree
// growth both in the matrix and the transformation.
// Variant: upon collision, pick first two entries of larger degree.
//
// Long description:
// proceed column by column. When looking for i-th pivot (in column j >= i) look
// at entries [i:m,j], and kill leading term of largest deg by using the second
// largest deg (using atomic transformation, adding a multiple of another row by
// constant*x**k); continue in the same column until all entries in [i:m,j]
// except one are zero; swap rows to put the nonzero one at [i,j] // TODO not swap, preserve MRP
//
// Complexity for a generic m x m degree d (hence nonsingular) matrix:
// - first, it will use about m*(d+1) steps to transform the first column
//   into the transpose of [1  0  0  ...  0]. Indeed each step decreases the
//   degree of a single one of the entries by exactly 1.
// - TODO complete this complexity analysis
slong nmod_poly_mat_hnf_revlex_atomic_ur(nmod_poly_mat_t mat, nmod_poly_mat_t tsf, slong * pivind)
{
    if (mat->r == 0 || mat->c == 0)
        return 0;

    // recall, pivind[i] gives index of pivot in row i
    // -> index rk below will tell us how many pivots have been found already,
    //     rk == rank of first j-1 columns
    //     pivind[i] for i >= rk is not used nor modified
    // Note: in case of zeroes/dependencies, we may find a pivot below rk,
    // then we will permute rows to bring it up at row rk
    slong rk = 0;

    // will be used to reduce against pivots
    nmod_poly_t u; nmod_poly_init(u, mat->modulus);
    nmod_poly_t v; nmod_poly_init(v, mat->modulus);

    // loop over the columns
    for (slong j = 0; j < mat->c; j++)
    {
        // whether there is still some collision among entries rk:mat->r,j
        int collision = 1;
        while (collision)
        {
            // pi1: find first nonzero entry
            slong pi1 = rk;
            while (pi1 < mat->r && nmod_poly_is_zero(MAT(pi1, j)))
                pi1++;
            if (pi1 == mat->r)
                collision = 0; // no collision, only zeroes
            else if (pi1 == mat->r-1)
                collision = 0; // no collision, pivot found
            else // look for pi2 providing collision
            {
                slong pi2 = pi1 + 1; // recall pi1+1 < mat->r
                // find pi1,pi2 such that entries pi1,j and pi2,j have the largest degree
                // among entries [rk:mat->r,j], with deg(mat[pi1,j]) >= deg(mat[pi2,j])
                if (MAT(pi1, j)->length < MAT(pi2, j)->length)
                    SLONG_SWAP(pi1, pi2);
                for (slong i = FLINT_MAX(pi1,pi2)+1; i < mat->r; i++)
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
                    collision = 0; // actually no collision, pivot found
                else // collision, work on it
                    _atomic_solve_pivot_collision_uechelon_rowwise(mat, tsf, pi1, pi2, j);
            }
            if (!collision && pi1 < mat->r)
            {
                // found pivot at [pi1,j]:
                // . permute rows to bring the pivot at [rk,j]TODO use rotation
                // . make pivot monic
                // . reduce entries above pivot
                // note: reduction could be grouped at the end, yet performing
                // it now may help to compute with smaller degrees in mat (but
                // may also bring larger degrees in transformation, if the
                // latter is computed)
                nmod_poly_mat_swap_rows(mat, NULL, pi1, rk);
                if (tsf)
                    nmod_poly_mat_swap_rows(tsf, NULL, pi1, rk);
                _normalize_pivot_uechelon_rowwise(mat, tsf, rk, j);
                for (slong ii = 0; ii < rk; ii++)
                    _reduce_against_pivot_uechelon_rowwise(mat, tsf, rk, j, ii, u, v);
                pivind[rk] = j;
                rk++;
            }
        }
    }

    nmod_poly_clear(u);
    nmod_poly_clear(v);
    return rk;
}

// Bradley's HNF algorithm       (upper echelon, row-wise)
//
// In short: pivoting strategy is "reverse lexicographic" (search for the
// topmost nonzero entry in the first nonzero column); use complete XGCD
// transformations. Benefits from fast polynomial arithmetic, but typically
// leads to (much) larger degrees than Rosser's algorithm both in the matrix
// and the transformation.
//
// Long description:
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
// TODO complexity
slong nmod_poly_mat_hnf_bradley_upper_rowwise(nmod_poly_mat_t mat, nmod_poly_mat_t tsf, slong * pivind)
{
    if (mat->r == 0 || mat->c == 0)
        return 0;

    // recall, pivind[i] gives index of pivot in row i
    // -> index rk below will tell us how many pivots have been found already,
    //     pivind[i] for i >= rk is not used nor modified

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
    for (slong j = 0; j < mat->c; j++)
    {
        // find actual pivot: starting at (rk,j), find first nonzero entry in the rest of column j
        slong pi = rk;
        while (pi < mat->r && nmod_poly_is_zero(MAT(pi, j)))
            pi++;

        if (pi < mat->r)
        // else, no pivot in this column, go to next column
        {
            // process column j below pivot
            slong ii = pi + 1;
            while (ii < mat->r)
            {
                // find next nonzero entry in column
                while (ii < mat->r && nmod_poly_is_zero(MAT(ii, j)))
                    ii++;
                if (ii < mat->r)
                    _complete_solve_pivot_collision_uechelon_rowwise(mat, tsf, pi, ii, j, g, u, v, pivg, nonzg);
                // else, ii==m, pivot is the only nonzero entry: nothing more to do for this column
            }

            // if pivot found in a lower row than expected
            // -> swap rows to bring it up to row rk
            nmod_poly_mat_swap_rows(mat, NULL, rk, pi); // does nothing if rk==pi
            if (tsf)
                nmod_poly_mat_swap_rows(tsf, NULL, rk, pi);
            pivind[rk] = j;
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
        slong j = pivind[i];
        // normalize pivot, and correspondingly update the row i of mat and of tsf
        _normalize_pivot_uechelon_rowwise(mat, tsf, i, j);
        // reduce entries above (i,j)
        for (slong ii = 0; ii < i; ii++)
            _reduce_against_pivot_uechelon_rowwise(mat, tsf, i, j, ii, u, v);
    }

    nmod_poly_clear(g);
    nmod_poly_clear(u);
    nmod_poly_clear(v);
    nmod_poly_clear(pivg);
    nmod_poly_clear(nonzg);

    return rk;
}

// Notes about Kannan-Bachem's algorithm/implementation
// (see below for a description of the algorithm)
// - about column permutations: one could use permutations to always work on
// the leading square principal submatrix, and then permute back at the very
// end (to go back to the original row space). But doing this, one should be
// careful to use suitable permutations only, like rotations, otherwise when
// permuting back one will not get a Hermite form (and, in addition, one would
// not obtain the matrix rank profile). For example: a pivot has been found at
// (0,0), and the second row is zero until its last entry which is nonzero
// (thus a pivot); then one may decide to rotate columns 0,1,2,...,n-1 ->
// 0,n-1,1,2,...,n-2, however doing for example a swap 0,1,2,...,n-1 ->
// 0,n-1,2,...,n-2,1 is not allowed (consider the expected next pivot, which
// would be the original 2,2 instead of the original 2,1).
// Here, since it is unclear to me (Vincent Neiger) that column rotations would
// simplify the code or make it more efficient, column permutations are avoided
// altogether. However such rotations would certainly be useful if one was to
// design a recursive block-variant of this algorithm. (TODO ?)
// - about row permutations:

// Kannan-Bachem's HNF algorithm       (upper echelon, row-wise)
//
// In short: pivoting strategy is "lexicographic" (search for the leftmost
// nonzero entry in the first nonzero row); use complete XGCD transformations.
// TODO a word on efficiency?
//
// Long description, trying to keep close to the original presentation:
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
// TODO complexity
slong nmod_poly_mat_hnf_kannan_bachem_upper_rowwise(nmod_poly_mat_t mat, nmod_poly_mat_t tsf, slong * pivind, slong * mrp)
{
    if (mat->r == 0 || mat->c == 0)
        return 0;

    // recall, pivind[i] gives index of pivot in row i
    // -> index i below will tell us how many pivots have been found already,
    //     pivind[ii] for ii >= i is not used nor modified

    // for simplicity also store reverse array
    // -> pivot_row[j] is either -1 (not among the pivots discovered so far)
    // or is the index of the row with pivot j
    slong * pivot_row = flint_malloc(mat->c * sizeof(slong));
    flint_mpn_store(pivot_row, mat->c, -1); // fill with -1

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
    while (i < mat->r && j < mat->c)
    {
        // Look for the first row ii in i:m such that [ii,0:j+1] is nonzero.
        // also find jj such that [ii,jj] is the first nonzero entry among [ii,0:j+1]
        slong ii = i;
        slong jj = 0;
        int row_is_zero = 1;
        while (ii < mat->r && row_is_zero)
        {
            jj = 0;
            while (jj <= j && nmod_poly_is_zero(MAT(ii, jj)))
                jj++;
            if (jj == j+1) // is zero, jump to next row
                ii++;
            else // jj <= j and mat[ii,jj] is nonzero
                row_is_zero = 0;
        }
        if (ii == mat->r) // no nonzero row, increment j and process next submatrix
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
                        _complete_solve_pivot_collision_uechelon_rowwise(mat, tsf, pi, i, jj, g, u, v, pivg, nonzg);
                        // note: new pivot already normalized since xgcd ensures the gcd is monic

                        // reduce appropriate entries in hnf in column jj
                        for (ii = 0; ii < pi; ii++)
                            _reduce_against_pivot_uechelon_rowwise(mat, tsf, pi, jj, ii, u, v);
                    }
                    else // pivot_row[jj] == -1, currently no pivot in column jj
                    {
                        row_is_zero = 0;
                        // move row i to the correct location in hnf and update pivot_row, pivot_col
                        ii = 0; // index where row i must be placed
                        while (ii < i && pivind[ii] < jj)
                            ii++;
                        // rotate mat and pivind
                        pivind[i] = jj;
                        _nmod_poly_mat_rotate_rows_downward(mat, pivind, ii, i);
                        if (tsf)
                            _nmod_poly_mat_rotate_rows_downward(tsf, NULL, ii, i);
                        pivot_row[jj] = i; // column jj is new pivot for row i
                        // update pivot row for found pivots in jj+1 ... j-1
                        for (slong pj = jj+1; pj < j; pj++)
                            if (pivot_row[pj] >= 0)
                                pivot_row[pj] += 1;

                        // the new row is now at index ii, with new pivot at ii,jj
                        // normalize the new pivot
                        _normalize_pivot_uechelon_rowwise(mat, tsf, ii, jj);

                        // reduce entries [ii,jj+1:j] against current hnf
                        for (slong kk = jj+1; kk < j; kk++)
                            if (pivot_row[kk] >= 0) // reduce mat[ii,kk] against pivot [pi,kk]
                                _reduce_against_pivot_uechelon_rowwise(mat, tsf, pivot_row[kk], kk, ii, u, v);
                            // else, no pivot in column kk, nothing to do
                        // reduce column jj of hnf against new pivot mat[ii,jj]
                        for (slong pii = 0; pii < ii; pii++)
                            _reduce_against_pivot_uechelon_rowwise(mat, tsf, ii, jj, pii, u, v);
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
                    pivind[i] = j;
                    pivot_row[j] = i;
                    // entry [i,j] is nonzero, this is a new pivot
                    // normalize it
                    _normalize_pivot_uechelon_rowwise(mat, tsf, i, j);
                    // reduce entries above it, i.e. ii,j against pivot i,j
                    for (ii = 0; ii < i; ii++)
                        _reduce_against_pivot_uechelon_rowwise(mat, tsf, i, j, ii, u, v);
                    // increment both i and j
                    i++;
                    j++;
                }
            }
        }
    }

    flint_free(pivot_row);

    nmod_poly_clear(g);
    nmod_poly_clear(u);
    nmod_poly_clear(v);
    nmod_poly_clear(pivg);
    nmod_poly_clear(nonzg);

    return i;
}

// put mat in upper echelon form, then normalize entries
slong nmod_poly_mat_hnf_lex_upper_rowwise(nmod_poly_mat_t mat, nmod_poly_mat_t tsf, slong * pivind, slong * mrp)
{
    if (mat->r == 0 || mat->c == 0)
        return 0;

    // recall, pivind[i] gives index of pivot in row i
    // -> index i below will tell us how many pivots have been found already,
    //     pivind[ii] for ii >= i is not used nor modified

    // for simplicity also store reverse array
    // -> pivot_row[j] is either -1 (not among the pivots discovered so far)
    // or is the index of the row with pivot j
    slong * pivot_row = flint_malloc(mat->c * sizeof(slong));
    flint_mpn_store(pivot_row, mat->c, -1); // fill with -1

    // rk: the (rk-1) x mat->c leading submatrix `hnf` is already in Hermite form,
    // we will now look for a new pivot in row i
    // zr: we have encountered zr zero rows at this point
    // --> numbers of already processed rows == rk + zr
    slong rk = 0;
    slong zr = 0;

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
    while (rk+zr < mat->r)
    {
        // look for new pivot at [rk,j], if it exists
        int pivot_found = 0;
        slong j = 0;
        while (j < mat->c && !pivot_found)
        {
            if (nmod_poly_is_zero(MAT(rk, j))) // entry is already zero
                j++;
            else if (pivot_row[j] >= 0) // make entry zero using existing pivot [pi,j]
            {
                // apply complete XGCD transformation between [pi,j] and [rk,j]
                // -> puts a zero at [rk,j] and the monic gcd at [pi,j] (updated pivot)
                _complete_solve_pivot_collision_uechelon_rowwise(mat, tsf, pivot_row[j], rk, j, g, u, v, pivg, nonzg);
                j++;
            }
            else // currently no pivot in column j => found new pivot
            {
                // TODO update mrp
                // move row rk to the correct location in hnf and update pivot information
                slong i = 0; // index where row rk must be placed
                while (i < rk && pivind[i] < j)
                    i++;
                pivind[rk] = j;
                _nmod_poly_mat_rotate_rows_downward(mat, pivind, i, rk);
                if (tsf)
                    _nmod_poly_mat_rotate_rows_downward(tsf, NULL, i, rk);
                pivot_row[j] = i;
                for (slong jj = j+1; jj < mat->c; jj++)
                    if (pivot_row[jj] >= 0)
                        pivot_row[jj] += 1;

                // normalize the new pivot (now at index i,j)
                _normalize_pivot_uechelon_rowwise(mat, tsf, i, j);
                // reduce entries [:i,j] of hnf against new pivot [i,j]
                // (in typical situations, helps keep working with low degrees)
                for (slong ii = 0; ii < i; ii++)
                    _reduce_against_pivot_uechelon_rowwise(mat, tsf, i, j, ii, u, v);
                rk++;
                pivot_found = 1;
                //printf("found pivot rk = %ld, zr = %ld, i = %ld, j = %ld\n",rk,zr,i,j);
            }
        }
        if (j == mat->c) // row is zero, rotate and go to next row
        {
            zr++;
            _nmod_poly_mat_rotate_rows_upward(mat, NULL, rk, mat->r-1);
            if (tsf)
                _nmod_poly_mat_rotate_rows_upward(tsf, NULL, rk, mat->r-1);
        }
    }

    // now, we have mat in upper echelon form, with monic pivots
    // -> reduce the entries above each pivot
    // Note: if pivots have been found in increasing pivot index (generic case),
    // the matrix is already normalized, this does essentially nothing
    for (slong i = 0; i < rk; i++)
    {
        slong j = pivind[i];
        for (slong ii = 0; ii < i; ii++)
            _reduce_against_pivot_uechelon_rowwise(mat, tsf, i, j, ii, u, v);
    }

    flint_free(pivot_row);

    nmod_poly_clear(g);
    nmod_poly_clear(u);
    nmod_poly_clear(v);
    nmod_poly_clear(pivg);
    nmod_poly_clear(nonzg);

    return rk;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
