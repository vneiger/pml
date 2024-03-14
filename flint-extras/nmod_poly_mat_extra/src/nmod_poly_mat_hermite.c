#include <flint/ulong_extras.h>
#include <flint/nmod_vec.h>
#include <flint/nmod.h>
#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_forms.h"

#define MAT(i,j) (mat->rows[i] + j)
#define TSF(i,j) (tsf->rows[i] + j)
#define OTHER(i,j) (other->rows[i] + j)

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
void _pivot_collision_xgcd_uref(nmod_poly_mat_t mat, nmod_poly_mat_t other,
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
void _reduce_against_pivot_uref(nmod_poly_mat_t mat, nmod_poly_mat_t other,
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
// . normalize pivot entry in given row of mat, applying the corresponding
// operation to the whole row of mat and the corresponding row of other
// . return the constant used for normalization
// (i,j) is position of the pivot entry
// other == NULL or other is another matrix with the same modulus
// Complexity: if already monic, 0, otherwise:
//    . for mat: sum_{jj = j ... mat->c-1} deg(mat[i,jj])
//          <= (mat->c - j) * rdeg(mat[i,:])
//    . for other: sum_{jj = 0 ... other->c-1} deg(other[i,jj])
//          <= other->c * rdeg(other[i,:])
mp_limb_t _normalize_pivot_uref(nmod_poly_mat_t mat, nmod_poly_mat_t other, slong i, slong j)
{
    if (nmod_poly_is_monic(MAT(i, j)))
        return 1UL;
    else
    {
        mp_limb_t inv = n_invmod(MAT(i, j)->coeffs[MAT(i, j)->length - 1], MAT(i, j)->mod.n);
        for (slong jj = j; jj < mat->c; jj++)
            _nmod_vec_scalar_mul_nmod(MAT(i, jj)->coeffs, MAT(i, jj)->coeffs, MAT(i, jj)->length, inv, MAT(i, jj)->mod);
        if (other)
            for (slong jj = 0; jj < other->c; jj++)
                _nmod_vec_scalar_mul_nmod(OTHER(i, jj)->coeffs, OTHER(i, jj)->coeffs, OTHER(i, jj)->length, inv, OTHER(i, jj)->mod);
        return inv;
    }
}

// Context: upper echelon, row-wise
// . normalize a matrix mat already in upper row echelon form, with rank rk and
// pivot indices pivind
// . apply the corresponding unimodular transformation to other
// (ignored if NULL is given for other)
// Note: one could easily add some operations to handle the determinant of the
// unimodular transformation (not done at the moment, since application is unclear)
void _normalize_uref(nmod_poly_mat_t mat, nmod_poly_mat_t other, slong * pivind, slong rk)
{
    // init temporaries
    nmod_poly_t u; nmod_poly_init(u, mat->modulus);
    nmod_poly_t v; nmod_poly_init(v, mat->modulus);

    // normalization
    for (slong i = 0; i < rk; i++)
    {
        slong j = pivind[i];
        // normalize pivot, and correspondingly update the row i of mat and of tsf
        _normalize_pivot_uref(mat, other, i, j);
        // reduce entries above (i,j)
        for (slong ii = 0; ii < i; ii++)
            _reduce_against_pivot_uref(mat, other, i, j, ii, u, v);
    }
    // clear temporaries
    nmod_poly_clear(u);
    nmod_poly_clear(v);
}


/**********************************************************************
**********************************************************************
*                   HERMITE NORMAL FORM ALGORITHMS                   *
**********************************************************************
**********************************************************************/

/**********************************************************************
*                         Rosser's algorithm                         *
*              orientation "ur": upper echelon, row-wise             *
**********************************************************************/

// Idea: Proceed column by column. When looking for i-th pivot (in column j >=
// i) look at entries [i:m,j], and kill leading term of largest degree by using
// the second largest degree (using atomic transformation); continue in the
// same column until all entries in [i:m,j] except one are zero; swap rows to
// put the nonzero one at [i,j]. This yields an upper echelon form of the
// matrix, it remains to reduce the off-"diagonal" entries
//
// Provides column rank profile. Does not provide row rank profile (or matrix
// rank profile). Notice that e.g. the largest degree row, be it in the rrp or
// not, may be killed if it is a cst*x**exp multiple of the second largest one.
// Explicit example of algorithm steps:
// [ 1     1   1 ] j=0 [ 1 1 1 ] j=0 [ 1  1   1 ] j=1 [ 1  1  1 ]
// [ x**2  x   0 ] --> [ 0 0 0 ] --> [ 0  0   0 ] --> [ 0 x-1 x ] --> end
// [ x     1   0 ]     [ x 1 0 ]     [ 0 1-x -x ]     [ 0  0  0 ]
// when putting zero on the second row, the algorithm has only looked at first
// column, it does not even know yet the row has become zero (or dependent on
// the first); the information that this row was originally independent from
// the first has been lost.
//
// Note: there is no advantage in reducing off-diagonal entries during
// echelonization: in any case, once a pivot has been found in a certain row,
// this row is not used anymore to find other pivots. Thus, we can totally
// separate the "row echelon form" step and the "normalization" step.
//
// Note: for a generic m x m degree d (hence nonsingular) matrix, degree growth
// is rather good, since after i steps the i x i leading principal submatrix
// has 1 on the diagonal and the trailing principal (m-i) x (m-i) submatrix is
// row reduced, or very close to it (meaning it has degrees very close to those
// obtained if one used a minimal kernel basis in the suitable part of the
// unimodular transformation: degrees are all around m*d / (m-i)).
//
// Enhancements? for specific non-generic instances, or for heterogenous row
// degrees, it would make sense to introduce shifts to guide the choice of
// pivot collisions instead of "two largest ones". This would make it more
// likely that the trailing principal submatrix is actually as small as it
// can get.

// Upper row echelon form inspired by Rosser's HNF algorithm
//
// Pivoting strategy: maximum degree
//     . collision = first two entries of largest degrees in the column
//     . use smaller one to reduce larger one
// Transformations for handling pivot collisions:
//     . atomic transformations (mat[pi1,:] = mat[pi1,:] + cst * x**exp * mat[pi2,:])
//
// Guarantee: unimodular transformation has determinant +1 or -1 (could easily
// say which one, following the number of row swaps, so this can be directly
// used for computing determinants of square matrices).
//
// Typically offers good control of the degree growth, both in the matrix and
// the transformation.
slong nmod_poly_mat_uref_maxdeg_atomic(nmod_poly_mat_t mat, nmod_poly_mat_t tsf, slong * pivind)
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

    // loop over the columns
    for (slong j = 0; j < mat->c; j++)
    {
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
                {
                    slong t = pi1;
                    pi1 = pi2;
                    pi2 = t;
                }
//                    SLONG_SWAP(pi1, pi2);
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
                // . permute rows to bring the pivot at [rk,j]
                nmod_poly_mat_swap_rows(mat, NULL, pi1, rk);
                if (tsf)
                    nmod_poly_mat_swap_rows(tsf, NULL, pi1, rk);
                pivind[rk] = j;
                rk++;
            }
        }
    }
    return rk;
}

/**********************************************************************
*                         Bradley's algorithm                        *
*              orientation "ur": upper echelon, row-wise             *
**********************************************************************/

// Idea: Proceed column by column. When looking for i-th pivot (in column j >=
// i) look at entries [i:,j], and repeatedly use gcd's between the two first
// nonzero entries in [i:,j], say at [pi,j] and [ii,j], to zero out [ii,j] and
// reduce [pi,j] as much as possible. If g = u * mat[pi,j] + v * mat[ii,j],
// then this means applying the following unimodular transformation to rows pi
// and ii:
//     [ mat[pi,:] ]  =  [ g ... ]  =  [    u       v  ]  *  [ mat[pi,:] ]
//     [ mat[ii,:] ]     [ 0 ... ]     [ -nonzg   pivg ]     [ mat[ii,:] ]
// where nonzg = mat[ii,j]/g   and   pivg = mat[pi,j]/g.
// This goes on in the same column until all entries in [i:,j] except one are
// zero; this nonzero entry gives the new pivot, we permute rows to put the
// nonzero one at [i,j].
//
// Provides column rank profile. Since we take the first nonzero entries (and
// not those two of largest degree like in Rosser's algorithm), this also
// provides the row rank profile. It also provides the matrix rank profile,
// by using only rotations for the row permutations.
//
// Note: there is no advantage in reducing off-diagonal entries during
// echelonization: in any case, once a pivot has been found in a certain row,
// this row is not used anymore to find other pivots. Thus, we can totally
// separate the "row echelon form" step and the "normalization" step.
//
// Note: for a generic m x m degree d (hence nonsingular) matrix, degree growth
// is rather bad, since the degrees seem to almost double at each of the first
// steps; this slows down at some point but remains a big issue for
// performance except when m << d.

// Upper row echelon form inspired by Bradley's HNF algorithm
//
// Pivoting strategy "revlex"
//      collision = first two nonzero entries in the column
// Transformations for handling pivot collisions: complete xgcd transformation
//        [ mat[pi,:] ]  =  [    u       v  ]  *  [ mat[pi,:] ]
//        [ mat[ii,:] ]     [ -nonzg   pivg ]     [ mat[ii,:] ]
//
// Guarantee: unimodular transformation has determinant +1 or -1 (could easily
// say which one, following the number of row swaps, so this can be directly
// used for computing determinants of square matrices).
//
// Benefits from fast polynomial arithmetic, but typically offers bad control
// of the degree growth, both in the matrix and the transformation.
slong nmod_poly_mat_uref_revlex_xgcd(nmod_poly_mat_t mat, nmod_poly_mat_t tsf, slong * pivind, slong * mrp)
{
    if (mat->r == 0 || mat->c == 0)
    {
        if (mrp)
            for (slong i = 0; i < mat->r; i++)
                mrp[i] = -1L;
        return 0;
    }

    // recall, pivind[i] gives index of pivot in row i
    // -> index rk below will tell us how many pivots have been found already,
    //     rk == rank of first j-1 columns
    //     pivind[i] for i >= rk is not used nor modified
    // Note: in case of zeroes/dependencies, we may find a pivot below rk,
    // then we will permute rows to bring it up at row rk
    slong rk = 0;

    // record row permutation, to fill the mrp
    // the row i of current mat is the row perm[i] of the input mat
    slong * perm = NULL;
    if (mrp)
    {
        for (slong i = 0; i < mat->r; i++)
            mrp[i] = -1L;
        perm = flint_malloc(mat->r * sizeof(slong));
        for (slong i = 0; i < mat->r; i++)
            perm[i] = i;
    }

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
        // find pivot: find first nonzero entry in [rk:,j]
        slong pi = rk;
        while (pi < mat->r && nmod_poly_is_zero(MAT(pi, j)))
            pi++;

        if (pi < mat->r) // otherwise, no pivot in this column
        {
            // process column j below pivot
            slong ii = pi + 1;
            while (ii < mat->r)
            {
                // find next nonzero entry in column
                while (ii < mat->r && nmod_poly_is_zero(MAT(ii, j)))
                    ii++;
                if (ii < mat->r) // this makes [ii,j] zero, so the while loop will eventually stop
                    _pivot_collision_xgcd_uref(mat, tsf, pi, ii, j, g, u, v, pivg, nonzg);
            }

            // rotate rows to bring the pivot up to row rk
            _nmod_poly_mat_rotate_rows_downward(mat, perm, rk, pi);
            if (tsf)
                _nmod_poly_mat_rotate_rows_downward(tsf, NULL, rk, pi);
            pivind[rk] = j;
            if (mrp)
                mrp[perm[rk]] = j;
            rk++;
        }
    }

    nmod_poly_clear(g);
    nmod_poly_clear(u);
    nmod_poly_clear(v);
    nmod_poly_clear(pivg);
    nmod_poly_clear(nonzg);
    if (mrp)
        flint_free(perm);

    return rk;
}

/**********************************************************************
*                             Pivoting: lex                          *
*              orientation "ur": upper echelon, row-wise             *
**********************************************************************/

// In short: proceed row by row, take leftmost nonzero entry in topmost nonzero
// row as next pivot (lex pivoting). After having done i steps, we have
// processed the first i rows of the matrix and not touched the others; we have
// computed an upper row echelon form of these i rows.
//
// Uses row permutation to place the zero rows at the bottom and to ensure the
// pivot indices are increasing. Rotations are used to preserve the matrix rank
// profile.
//
// In generic cases, the degrees in the matrix along the computation are thus
// split into two parts: the first i rows have the expected degrees for an
// upper row echelon form of an i x mat->c matrix (*), and the remaining rows
// haven't been touched. This does not yield very good behaviour when applied
// to a generic m x m matrix of degree d (compared to Rosser's strategy, for
// example).
//
// (*) pivots are all 1 except the last,
//     row 1 has degree the sum of the degrees of the input's rows 0 + 1,
//     row 2 has degree the sum of the degrees of the input's rows 0 + 1 + 2,
//     etc.

// Upper row echelon form using lex pivot search
// (leftmost nonzero entry in topmost nonzero row)
//
// Proceed row by row, with lex pivoting strategy
// Transformations for reducing entries below existing pivots:
// complete xgcd transformation
//        [ mat[pi1,:] ]  =  [    u       v  ]  *  [ mat[pi1,:] ]
//        [ mat[pi2,:] ]     [ -nonzg   pivg ]     [ mat[pi2,:] ]
//
// Guarantee: unimodular transformation has determinant +1 or -1 (could easily
// say which one, following the number of row rotations, so this can be
// directly used for computing determinants of square matrices).
//
// Benefits from fast polynomial arithmetic, but typically offers bad control
// of the degree growth, both in the matrix and the transformation.
slong nmod_poly_mat_uref_lex_xgcd(nmod_poly_mat_t mat, nmod_poly_mat_t tsf, slong * pivind, slong * mrp)
{
    if (mat->r == 0 || mat->c == 0)
    {
        if (mrp)
            for (slong i = 0; i < mat->r; i++)
                mrp[i] = -1L;
        return 0;
    }

    // rk: the (rk-1) x mat->c leading submatrix `hnf` is already in Hermite form,
    // we will now look for a new pivot in row i
    // zr: we have encountered zr zero rows at this point
    // --> numbers of already processed rows == rk + zr
    slong rk = 0;
    slong zr = 0;

    // for simplicity also store reverse array
    // -> pivot_row[j] is either -1 (not among the pivots discovered so far)
    // or is the index of the row with pivot j
    slong * pivot_row = flint_malloc(mat->c * sizeof(slong));
    flint_mpn_store(pivot_row, mat->c, -1); // fill with -1

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

    // while there remains some row to process
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
                _pivot_collision_xgcd_uref(mat, tsf, pivot_row[j], rk, j, g, u, v, pivg, nonzg);
                j++;
            }
            else // currently no pivot in column j => found new pivot
            {
                // move row rk to the correct location in hnf and update pivot information
                slong i = 0; // index where row rk must be placed
                while (i < rk && pivind[i] < j)
                    i++;
                pivind[rk] = j; // must be modified before rotation
                _nmod_poly_mat_rotate_rows_downward(mat, pivind, i, rk);
                if (tsf)
                    _nmod_poly_mat_rotate_rows_downward(tsf, NULL, i, rk);
                pivot_row[j] = i;
                for (slong jj = j+1; jj < mat->c; jj++)
                    if (pivot_row[jj] >= 0)
                        pivot_row[jj] += 1;

                if (mrp)
                    mrp[rk+zr] = j;
                rk++;
                pivot_found = 1;
            }
        }
        if (j == mat->c) // row is zero, rotate and go to next row
        {
            if (mrp)
                mrp[rk+zr] = -1L;
            zr++;
            _nmod_poly_mat_rotate_rows_upward(mat, NULL, rk, mat->r-1);
            if (tsf)
                _nmod_poly_mat_rotate_rows_upward(tsf, NULL, rk, mat->r-1);
        }
    }

    flint_free(pivot_row);

    nmod_poly_clear(g);
    nmod_poly_clear(u);
    nmod_poly_clear(v);
    nmod_poly_clear(pivg);
    nmod_poly_clear(nonzg);

    return rk;
}



/**********************************************************************
*                      Kannan&Bachem's algorithm                     *
*              orientation "ur": upper echelon, row-wise             *
**********************************************************************/

// The choice of pivots is eventually the same as in the direct revlex pivoting
// strategy. However, the scheduling of when to apply transformations differs.
//
// This proceeds by leading submatrices. When looking for i-th pivot, in column
// j >= i, the (i-1) x (j-1) leading submatrix is already in Hermite form.
//
// Then, ensure entries [i,jj] for jj in 0:j are zero. Look for the first row
// ii in i:m such that [ii,0:j+1] is nonzero.
// -> if there is none increment j and process next submatrix
// -> otherwise, use xgcd transformations between existing pivots and the
// entries [ii,0:j]. This will possibly update the existing pivots, and this
// will surely make [ii,0:j] zero. Then check if entry [ii,j] is nonzero: if
// yes, this is a new pivot, use row rotations to place it at (i,j), increment
// both i and j and proceed to next leading submatrix; if no, repeat with the
// remaining rows below the ii-th one.
//
// Ensure HNF: each time a pivot is found or updated, reduce the suitable
// entries of the current leading submatrix. (Experiments seem to suggest that
// -- even for obtaining just the row echelon form -- doing this is faster than
// removing all normalization steps.)
//
// Provides column rank profile. Since we take the first nonzero entries (and
// not those two of largest degree like in Rosser's algorithm), this also
// provides the row rank profile. It also provides the matrix rank profile,
// by using only rotations for the row permutations.

// Kannan-Bachem's HNF algorithm       (upper echelon, row-wise)
// Adaptation, removing the assumptions (rank properties) from the original.
//
// Take same pivots as revlex strategy, but proceeds by leading principal
// submatrices. Uses xgcd transformations, and continuous normalization.
//
// No guarantee concerning the determinant of the unimodular transformation
// (but, if useful, this determinant could easily be returned).
//
// Benefits from fast polynomial arithmetic. Typically offers worse control
// of the degree growth than Rosser's algorithm, but better than the revlex
// strategy.
slong nmod_poly_mat_hnf_ur_revlex_xgcd_delayed_zero(nmod_poly_mat_t mat, nmod_poly_mat_t tsf, slong * pivind, slong * mrp)
{
    if (mat->r == 0 || mat->c == 0)
    {
        if (mrp)
            for (slong i = 0; i < mat->r; i++)
                mrp[i] = -1L;
        return 0;
    }

    // recall, pivind[i] gives index of pivot in row i
    // -> index i below will tell us how many pivots have been found already,
    //     pivind[ii] for ii >= i is not used nor modified

    // for simplicity also store reverse array
    // -> pivot_row[j] is either -1 (not among the pivots discovered so far)
    // or is the index of the row with pivot j
    slong * pivot_row = flint_malloc(mat->c * sizeof(slong));
    flint_mpn_store(pivot_row, mat->c, -1); // fill with -1

    // record row permutation, to fill the mrp
    // the row k of current mat is the row perm[k] of the input mat
    slong * perm = NULL;
    if (mrp)
    {
        for (slong i = 0; i < mat->r; i++)
            mrp[i] = -1L;
        perm = flint_malloc(mat->r * sizeof(slong));
        for (slong i = 0; i < mat->r; i++)
            perm[i] = i;
    }

    // (i,j) : the (i-1) x (j-1) leading submatrix `hnf` is already in Hermite form,
    // we will now try to increase either i (if new pivot found in row i or
    // below), or j (if submatrix [i:m,0:j+1] does not bring any new pivot), or
    // both (if new pivot found in column j)
    slong i = 0;
    slong j = 0;
    // ii: index of currently processed row below the i-th row
    // each time we increment i or j, ii will be reset to i
    slong ii = 0;

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
        // look for the next row ii such that [ii,0:j+1] is nonzero.
        // also find jj such that [ii,jj] is the first nonzero entry among [ii,0:j+1]
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
        // if no nonzero row, reset ii and go to next j
        if (ii == mat->r) { ii = i; j++; }
        else // found nonzero row
        {
            // make entries ii,jj:j zero
            for (; jj < j; jj++)
            {
                if (! nmod_poly_is_zero(MAT(ii, jj)))
                {
                    // a pivot has necessarily been found at [pi,jj], for some pi < i
                    slong pi = pivot_row[jj];
                    // perform complete xgcd transformation between [pi,jj] and [ii,jj]
                    // -> puts a zero at [ii,jj] and the gcd (updated pivot, monic) at [pi,jj]
                    _pivot_collision_xgcd_uref(mat, tsf, pi, ii, jj, g, u, v, pivg, nonzg);
                    // reduce appropriate entries in column jj
                    for (slong pii = 0; pii < pi; pii++)
                        _reduce_against_pivot_uref(mat, tsf, pi, jj, pii, u, v);
                }
            }
            // now all entries in [ii,0:j] are zero
            if (nmod_poly_is_zero(MAT(ii, j))) // no new pivot in [i:ii+1,0:j+1]
                ii++;
            else // found pivot at [ii,j]
            {
                // swap rows i and ii
                _nmod_poly_mat_rotate_rows_downward(mat, perm, i, ii);
                if (tsf)
                    _nmod_poly_mat_rotate_rows_downward(tsf, NULL, i, ii);
                // update pivot_col, pivot_row, mrp
                pivind[i] = j;
                pivot_row[j] = i;
                if (mrp)
                    mrp[perm[i]] = j;
                // make pivot monic and reduce entries above it, i.e. ii,j against i,j
                _normalize_pivot_uref(mat, tsf, i, j);
                for (slong pii = 0; pii < i; pii++)
                    _reduce_against_pivot_uref(mat, tsf, i, j, pii, u, v);
                // increment i and j, reset ii
                i++;
                j++;
                ii = i;
            }
        }
    }

    flint_free(pivot_row);

    nmod_poly_clear(g);
    nmod_poly_clear(u);
    nmod_poly_clear(v);
    nmod_poly_clear(pivg);
    nmod_poly_clear(nonzg);
    if (mrp)
        flint_free(perm);

    return i;
}



/**********************************************************************
*                    Mulders&Storjohann's algorithm                  *
*              orientation "ur": upper echelon, row-wise             *
**********************************************************************/


// Algo of Mulders&Storjohann, Algo 7, with a slight modification: performs
// upper row echelon form; no reduction of above-pivot entries for already
// found pivots. If wanting the HNF, performs a complete normalization step at
// the end; there should be some rare cases where this is bad in terms of
// complexity and where continuously reducing off-diagonal entries is
// preferable.
//
// Note: in the original presentation, the algorithm is only guaranteed to work
// for a full column rank matrix. Here, the implementation supports any rectangular
// matrix of dimension m x n:
// * if m > n, it will first use a weak Popov form computation to reduce to the
// case m <= n;
// * if m <= n and the m x m leftmost submatrix is nonsingular (i.e. the column
// rank profile is 0,...,m-1), then the algorithm succeeds;
// * if m <= n and the m x m leftmost submatrix is singular, then the algorithm
// will exit as soon as it detects this singularity, returning a strictly
// negative value
//
// This is using the upper weak Popov form (unlike in the original
// presentation), since this has better properties w.r.t the target upper
// echelon form / HNF.  In the generic m x m nonsingular case, the first weak
// Popov computation transforms the first m-1 columns into an invertible upper
// triangular matrix with a row of zeroes below it; the subsequent weak Popov
// form computations just do nothing (whereas they would perform some atomic
// transformations if we used lower weak Popov form).
slong nmod_poly_mat_uref_matrixgcd_iter(nmod_poly_mat_t mat,
                                        nmod_poly_mat_t tsf,
                                        slong * pivind,
                                        slong * rrp,
                                        slong * udet)
{
    if (mat->r == 0 || mat->c == 0)
        return 0;

    // min-dim (estimated rank)
    const slong mdim = FLINT_MIN(mat->r,mat->c);
    // permutation for putting into ordered weak Popov
    slong * perm = flint_malloc(mat->r * sizeof(slong));
    // TODO no need to allocate perm if both tsf==NULL and udet == NULL
    // (yet currently permute_rows_by_sorting_vec requires an allocated perm)

    // reducing to (assumed) generic column rank profile: weak Popov of
    // mat[:,:mdim-1], with transformation applied to the whole mat and tsf
    // --> for most aspects, it would be easier to compute weak Popov of
    // mat[:,:mdim], but this would slow down the (frequent) nonsingular case
    // for which it would mean starting with an extra computation of a weak
    // Popov form of mat
    // --> recall rk == rank(mat[:,mdim-1])
    slong rk = _nmod_poly_mat_weak_popov_iter_submat_rowbyrow(mat, NULL, tsf, udet, pivind, rrp, 0, 0, mat->r, mdim-1, mat->r, ROW_UPPER);

    // since we did not do mat[:,:mdim] but mat[:,:mdim-1], there is now some
    // work to confirm the rank and rrp, and also check the rank requirements

    // mat[:rk,:mdim-1] has rank rk, and mat[rk:,:mdim-1] is now zero
    // --> either we have found the actual rank rk of mat, hence all
    // of mat[rk:,:] is zero (maybe not generic CRP, but still fine)
    // --> or the rank of mat is rk+1 and there is a nonzero in mat[rk:,rk]
    // (pivot we missed, fine)
    // --> or mat[rk:,rk] is zero but there is a nonzero in mat[rk:,rk+1:],
    // hence column rank profile is not generic (fails requirement)

    // find the first nonzero in mat[rk:,rk]
    slong pi = rk;
    while (pi < mat->r && nmod_poly_is_zero(MAT(pi, rk)))
        pi++;
    if (pi < mat->r)
    {
        // we missed a pivot, in the column mat[:,rk]
        // -> rank rk should be incremented (at end of scope for convenience)
        // -> TODO the row rank profile of the input mat is rrp modified by adding
        // pi-i to the rrp, meaning the "pi-i"-th row of mat was actually in
        // the row rank profile (this is valid because weak Popov uses row
        // rotations to put zero rows at the end)

        // handle additional pivot: update pivind and rrp
        if (rrp)
        {
            printf("RRP!!!\n");
            // TODO to be completed
            // find the "pi-i -th gap" in rrp
            slong iii = 0; slong ii = 0;
            while (iii < pi-rk && ii < iii)
            {
                if (iii < rrp[ii])
                    iii++;
                else
                    ii++;
            }
        }
        // while there is another nonzero below pi, xgcd-solve collision
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
        for (slong ii = pi+1; ii < mat->r; ii++)
            if (!nmod_poly_is_zero(MAT(ii,rk)))
                _pivot_collision_xgcd_uref(mat, tsf, pi, ii, rk, g, u, v, pivg, nonzg);
        nmod_poly_clear(g);
        nmod_poly_clear(u);
        nmod_poly_clear(v);
        nmod_poly_clear(pivg);
        nmod_poly_clear(nonzg);

        // swap rows to put the pivot row pi in position rk
        nmod_poly_mat_swap_rows(mat, NULL, rk, pi);
        if (tsf) nmod_poly_mat_swap_rows(tsf, NULL, rk, pi);
        rk++;
    }
    // check the rank requirement has not failed yet: mat[rk:,:] should be zero
    // (we only test mat[rk:,rk:]: the rest is zero for sure)
    nmod_poly_mat_t view;
    nmod_poly_mat_window_init(view, mat, rk, rk, mat->r, mat->c);
    int pass = nmod_poly_mat_is_zero(view);
    nmod_poly_mat_clear(view);
    if (!pass)
        return -rk-1;

    // the algorithm can proceed with mat[:rk,:]
    slong i = rk;
    while (i > 0)
    {
        //        [ M  * ]
        // mat is [ 0  E ] with E in upper echelon form (including possible
        // zero rows), and M of size i x i (a priori, nonsingular)
        //                 [ N  *  * ]
        // -> transform to [ 0  p  * ] via weak Popov form of MM = M[:i,:i-1]
        //                 [ 0  0  E ]
        // case 1: p == mat[i,j] is nonzero -> continue with i-1
        // case 2: p == mat[i,j] is zero -> fail, rank requirement

        // the next call applies weak Popov form on mat[:i,:i-1],
        // with transformations applied to the whole rows mat[:i,:],
        // and with update of the determinant of unimodular transformation (+1 or -1)
        slong new_rk = _nmod_poly_mat_weak_popov_iter_submat_rowbyrow(mat, NULL, tsf, udet, pivind, NULL, 0, 0, i, i-1, 2, ROW_UPPER);

        // test failed rank requirement, including early exit new_rk < 0
        // (note: early exit new_rk == 0 is well covered by this as well)
        if (new_rk < i-1)
            return -rk;
        // else found <= 1 zero rows (thus exactly 1): the rank of MM is ok
        i = i-1;
        if (nmod_poly_is_zero(MAT(i,i))) // rank of M not ok
            return -rk;
        // rank ok until now, permute into ordered weak Popov form
        _nmod_poly_mat_permute_rows_by_sorting_vec(mat, i, pivind, perm);
        if (tsf) nmod_poly_mat_permute_rows(tsf, perm, NULL);
        if (udet && _perm_parity(perm, i)) // odd permutation, negate udet
            *udet = - *udet;
    }
    // fill pivind with generic column rank profile
    for (slong j = 0; j < rk; j++)
        pivind[j] = j;
    flint_free(perm);
    return rk;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
