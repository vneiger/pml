#include <flint/flint.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_vec.h>
#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_forms.h"
#include "nmod_poly_mat_io.h" // TODO remove, for debugging

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

void _apply_pivot_collision_transformation(nmod_poly_mat_t mat, nmod_poly_mat_t tsf,
                                           slong pi1, slong pi2, slong j)
{
    // use pi2,j to kill the leading term of pi1,j, and apply the
    // corresponding transformation to the whole row pi1
    //  --> mat[pi1,:] = mat[pi1,:] + cst * x**k * mat[pi2,:]
    // where cst = - leading_coeff(mat[pi1,j]) / leading_coeff(mat[pi2,j])
    //   and k = deg(mat[pi1,j]) - deg(mat[pi2,j])
    nmod_poly_struct * entry1 = nmod_poly_mat_entry(mat, pi1, j);
    nmod_poly_struct * entry2 = nmod_poly_mat_entry(mat, pi2, j);
    const slong len1 = entry1->length;
    const slong len2 = entry2->length;
    mp_limb_t cst = n_invmod(entry2->coeffs[len2-1], mat->modulus);
    cst = nmod_mul(entry1->coeffs[len1-1], cst, entry1->mod);
    cst = nmod_neg(cst, entry1->mod);
    _nmod_vec_scalar_addmul_nmod(entry1->coeffs+len1-len2, entry2->coeffs, entry2->length, cst, entry1->mod);
    _nmod_poly_normalise(entry1);
    for (slong jj = j+1; jj < mat->c; jj++)
    {
        entry1 = nmod_poly_mat_entry(mat, pi1, jj);
        entry2 = nmod_poly_mat_entry(mat, pi2, jj);
        const slong len = entry2->length + len1-len2;
        nmod_poly_fit_length(entry1, len);
        if (len > entry1->length)
        {
            _nmod_vec_zero(entry1->coeffs + entry1->length, len - entry1->length);
            _nmod_poly_set_length(entry1, len);
        }
        _nmod_vec_scalar_addmul_nmod(entry1->coeffs+len1-len2, entry2->coeffs, entry2->length, cst, entry1->mod);
        _nmod_poly_normalise(entry1);
    }
    if (tsf) for (slong jj = 0; jj < mat->r; jj++)
    {
        entry1 = nmod_poly_mat_entry(tsf, pi1, jj);
        entry2 = nmod_poly_mat_entry(tsf, pi2, jj);
        const slong len = entry2->length + len1-len2;
        nmod_poly_fit_length(entry1, len);
        if (len > entry1->length)
        {
            _nmod_vec_zero(entry1->coeffs + entry1->length, len - entry1->length);
            _nmod_poly_set_length(entry1, len);
        }
        _nmod_vec_scalar_addmul_nmod(entry1->coeffs+len1-len2, entry2->coeffs, entry2->length, cst, entry1->mod);
        _nmod_poly_normalise(entry1);
    }
}


// computing xgcd and applying corresponding unimodular transformation
// input is mat, tsf, pi, ii, j
// pivot is at pi,j
// we perform xgcd between pivot and mat[ii,j]
// g, u, v, pivg, nonzg are used as temporaries and must be already initialized,
// piv and nonz may be NULL
void _apply_xgcd_transformation(nmod_poly_t g, nmod_poly_t u, nmod_poly_t v,
                                nmod_poly_t pivg, nmod_poly_t nonzg,
                                nmod_poly_struct * piv, nmod_poly_struct * entry,
                                nmod_poly_mat_t mat, nmod_poly_mat_t tsf,
                                slong pi, slong ii, slong j)
{
    piv = nmod_poly_mat_entry(mat, pi, j); // pivot
    entry = nmod_poly_mat_entry(mat, ii, j); // nonzero below pivot

    // FIXME this specific code for piv == 1 does not seem to speed up things
    if (nmod_poly_is_one(piv))
    {
        // cancel (ii,j) and modify row ii accordingly
        // in short: for M in {mat,tsf} do the unimodular transformation
        //   M[ii,:] = M[ii,:] - mat[ii,j] * M[pi,:]
        // make entry zero since it is used throughout (it is mat[ii,j])
        // use g as temporary all along
        for (slong jj = j+1; jj < mat->c; jj++)
        {
            nmod_poly_mul(g, entry, nmod_poly_mat_entry(mat, pi, jj));
            nmod_poly_sub(nmod_poly_mat_entry(mat, ii, jj), nmod_poly_mat_entry(mat, ii, jj), g);
        }
        if (tsf)
        {
            for (slong jj = 0; jj < mat->r; jj++)
            {
                nmod_poly_mul(g, entry, nmod_poly_mat_entry(tsf, pi, jj));
                nmod_poly_sub(nmod_poly_mat_entry(tsf, ii, jj), nmod_poly_mat_entry(tsf, ii, jj), g);
            }
        }
        nmod_poly_zero(entry);
    }
    else // general case, pivot is not one
    {
        nmod_poly_xgcd(g, u, v, piv, entry); // g = u*piv + v*entry
        nmod_poly_divides(pivg, piv, g); // pivg = piv // g
        nmod_poly_divides(nonzg, entry, g); // nonzg = entry // g

        // cancel (ii,j) and modify row ii accordingly
        // in short: for M in {mat,tsf} do the unimodular transformation on rows pi,ii:
        //  [ M[pi,:] ]  =  [    u       v  ]  *  [ M[pi,:] ]
        //  [ M[ii,:] ]     [ -nonzg   pivg ]     [ M[ii,:] ]

        // for mat, due to zeroes already set in columns 0...j-1, we can start at column j
        // set piv to g, and later below set entry to 0
        nmod_poly_set(piv, g); // mat[pi,j] = g
        for (slong jj = j+1; jj < mat->c; jj++)
        {
            // simultaneously update:
            //     mat[pi,jj] = u * mat[pi,jj] + v * mat[ii,jj]
            //     mat[ii,jj] = -nonzg * mat[pi,jj] + pivg * mat[ii,jj]
            // --> use g as temporary copy of mat[pi,jj]
            // --> use entry as temporary for storing products
            nmod_poly_set(g, nmod_poly_mat_entry(mat, pi, jj));

            nmod_poly_mul(nmod_poly_mat_entry(mat, pi, jj), u, g);
            nmod_poly_mul(entry, v, nmod_poly_mat_entry(mat, ii, jj));
            nmod_poly_add(nmod_poly_mat_entry(mat, pi, jj), nmod_poly_mat_entry(mat, pi, jj), entry);

            nmod_poly_mul(entry, pivg, nmod_poly_mat_entry(mat, ii, jj));
            nmod_poly_mul(nmod_poly_mat_entry(mat, ii, jj), nonzg, g);
            nmod_poly_sub(nmod_poly_mat_entry(mat, ii, jj), entry, nmod_poly_mat_entry(mat, ii, jj));
        }
        if (tsf)
        {
            for (slong jj = 0; jj < mat->r; jj++)
            {
                // simultaneously update:
                //     tsf[pi,jj] = u * tsf[pi,jj] + v * tsf[ii,jj]
                //     tsf[ii,jj] = -nonzg * tsf[pi,jj] + pivg * tsf[ii,jj]
                // --> use g as temporary copy of tsf[pi,jj]
                // --> use entry as temporary for storing products
                nmod_poly_set(g, nmod_poly_mat_entry(tsf, pi, jj));

                nmod_poly_mul(nmod_poly_mat_entry(tsf, pi, jj), u, g);
                nmod_poly_mul(entry, v, nmod_poly_mat_entry(tsf, ii, jj));
                nmod_poly_add(nmod_poly_mat_entry(tsf, pi, jj), nmod_poly_mat_entry(tsf, pi, jj), entry);

                nmod_poly_mul(entry, pivg, nmod_poly_mat_entry(tsf, ii, jj));
                nmod_poly_mul(nmod_poly_mat_entry(tsf, ii, jj), nonzg, g);
                nmod_poly_sub(nmod_poly_mat_entry(tsf, ii, jj), entry, nmod_poly_mat_entry(tsf, ii, jj));
            }
        }
        nmod_poly_zero(entry); // mat[ii,j] = 0
    }
}

// normalize pivot entry in given row of mat, applying the corresponding
// operation to the whole row of mat and the corresponding row of the transformation tsf
// (i,j) is position of the pivot entry
// piv points to the pivot entry
// entry will be used as temp
void _normalize_pivot(nmod_poly_mat_t mat, nmod_poly_mat_t tsf, nmod_poly_struct * piv,
                      slong i, slong j,
                      nmod_poly_struct * entry)
{
    if (! nmod_poly_is_monic(piv))
    {
        mp_limb_t inv = n_invmod(piv->coeffs[piv->length - 1], piv->mod.n);
        for (slong jj = j; jj < mat->c; jj++)
        {
            entry = nmod_poly_mat_entry(mat, i, jj);
            _nmod_vec_scalar_mul_nmod(entry->coeffs, entry->coeffs, entry->length, inv, entry->mod);
        }
        if (tsf)
        {
            for (slong jj = 0; jj < mat->r; jj++)
            {
                entry = nmod_poly_mat_entry(tsf, i, jj);
                _nmod_vec_scalar_mul_nmod(entry->coeffs, entry->coeffs, entry->length, inv, entry->mod);
            }
        }
    }
}

// reduce entry mat[ii,j] against pivot entry piv, which is at mat[i,j],
// applying the corresponding operation to the whole row ii of mat and the
// corresponding row ii of the transformation tsf
// entry, u, v are used as temporaries (u, v must be already initialized, entry can be NULL)
void _reduce_against_pivot(nmod_poly_mat_t mat, nmod_poly_mat_t tsf, nmod_poly_struct * piv,
                           slong i, slong j, slong ii,
                           nmod_poly_struct * entry, nmod_poly_t u, nmod_poly_t v)
{
    entry = nmod_poly_mat_entry(mat, ii, j);
    if (entry->length >= piv->length)
    {
        // division with remainder: entry1 = entry2 * u + v
        nmod_poly_divrem(u, v, entry, piv);
        // replace entry1 by remainder
        nmod_poly_swap(entry, v);
        // apply same transformation on rest of the row: row_ii <- row_ii - u * row_i
        nmod_poly_neg(u, u);
        for (slong jj = j+1; jj < mat->c; jj++)
        {
            nmod_poly_mul(v, u, nmod_poly_mat_entry(mat, i, jj));
            nmod_poly_add(nmod_poly_mat_entry(mat, ii, jj), nmod_poly_mat_entry(mat, ii, jj), v);
        }
        if (tsf)
        {
            for (slong jj = 0; jj < mat->r; jj++)
            {
                nmod_poly_mul(v, u, nmod_poly_mat_entry(tsf, i, jj));
                nmod_poly_add(nmod_poly_mat_entry(tsf, ii, jj), nmod_poly_mat_entry(tsf, ii, jj), v);
            }
        }
    }
}


// Rosser's algorithm
// proceed column by column. When looking for i-th pivot (in column j >= i) look
// at entries [i:m,j], and kill leading term of largest deg by using the second
// largest deg (using basic operation only, adding a multiple of another row by
// constant*x**k); continue in the same column until all entries in [i:m,j]
// except one are zero; swap rows to put the nonzero one at [i,j]
// --> TODO try version using divrem and poly multiplication instead of constant*x**k?
// --> TODO put the normalization step inside the main triangularization loop?
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

    slong len1, len2, len; // will store length of polynomials
    mp_limb_t cst; // will store constant multipliers (quotients of leading coefficients)
    nmod_poly_struct * entry1 = NULL; // will point to some matrix entry

    // loop over the columns
    for (slong j = 0; j < n; j++)
    {
        // whether there are still nonzero elements in entries rk+1:m,j
        int collision = 1;
        while (collision)
        {
            // find i1 and i2 such that entries (i1,j) and (i2,j) are those of largest degrees,
            // with the former greater than or equal to the former, among the entries
            // (i,j) for i = rk ... m
            // 1. find first nonzero entry
            slong pi1 = rk;
            while (pi1 < m && nmod_poly_is_zero(nmod_poly_mat_entry(mat, pi1, j)))
                pi1++;
            if (pi1 == m) // there were only zeroes in rk:m,j
                collision = 0;
            else if (pi1 == m-1)
                // only nonzero entry in rk:m,j is pi1,j --> pivot found, j-th iteration finished
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
                if (nmod_poly_mat_entry(mat, pi1, j)->length < nmod_poly_mat_entry(mat, pi2, j)->length)
                    SLONG_SWAP(pi1, pi2);
                for (slong i = FLINT_MAX(pi1,pi2)+1; i < m; i++)
                {
                    if (nmod_poly_mat_entry(mat, i, j)->length > nmod_poly_mat_entry(mat, pi1, j)->length)
                    {
                        pi2 = pi1;
                        pi1 = i;
                    }
                    else if (nmod_poly_mat_entry(mat, i, j)->length > nmod_poly_mat_entry(mat, pi2, j)->length)
                        pi2 = i;
                }
                if (nmod_poly_is_zero(nmod_poly_mat_entry(mat, pi2, j)))
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
                {
                    _apply_pivot_collision_transformation(mat, tsf, pi1, pi2, j);
                }
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
        nmod_poly_struct * piv = nmod_poly_mat_entry(mat, i, j); // pivot entry

        // normalize pivot, and correspondingly update the row i of mat and of tsf
        _normalize_pivot(mat, tsf, piv, i, j, entry1); // TODO entry1

        // reduce entries above (i,j)
        for (slong ii = 0; ii < i; ii++)
            _reduce_against_pivot(mat, tsf, piv, i, j, ii, entry1, u, v);
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
// --> TODO put the normalization step inside the main triangularization loop?
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

    nmod_poly_struct * piv = NULL; // will point to pivot entry
    nmod_poly_struct * entry = NULL; // will point to some entry
    nmod_poly_t g; // gcd
    nmod_poly_t u; // gcd cofactor1
    nmod_poly_t v; // gcd cofactor2
    nmod_poly_t pivg; // piv divided by gcd
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
        while (pi < m && nmod_poly_is_zero(nmod_poly_mat_entry(mat, pi, j)))
            pi++;

        if (pi < m)
        // else, no pivot in this column, go to next column
        {
            // process column j below pivot
            slong ii = pi + 1;
            while (ii < m)
            {
                // find next nonzero entry in column
                while (ii < m && nmod_poly_is_zero(nmod_poly_mat_entry(mat, ii, j)))
                    ii++;
                if (ii < m)
                    _apply_xgcd_transformation(g, u, v, pivg, nonzg, piv, entry, mat, tsf, pi, ii, j);
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
    // 1. reduce the entries above each pivot
    // 2. normalize each pivot
    for (slong i = 0; i < rk; i++)
    {
        slong j = pivot_col[i];
        piv = nmod_poly_mat_entry(mat, i, j);

        // normalize pivot, and correspondingly update the row i of mat and of tsf
        _normalize_pivot(mat, tsf, piv, i, j, entry);

        // reduce entries above (i,j)
        for (slong ii = 0; ii < i; ii++)
            _reduce_against_pivot(mat, tsf, piv, i, j, ii, entry, u, v);
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

    nmod_poly_struct * piv = NULL; // will point to pivot entry
    nmod_poly_struct * nonz = NULL; // will point to nonpivot entry that is bound to become zero
    nmod_poly_struct * entry = NULL; // will point to some entry
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
        //printf("i = %ld - j = %ld - IN  degree matrix:\n", i, j);
        //nmod_poly_mat_degree_matrix_print_pretty(mat);
        // Look for the first row ii in i:m such that [ii,0:j+1] is nonzero.
        // also find jj such that [ii,jj] is the first nonzero entry among [ii,0:j+1]
        slong ii = i;
        slong jj = 0;
        int row_is_zero = 1;
        while (ii < m && row_is_zero)
        {
            jj = 0;
            while (jj <= j && nmod_poly_is_zero(nmod_poly_mat_entry(mat, ii, jj)))
                jj++;
            if (jj == j+1) // is zero, jump to next row
                ii++;
            else // jj <= j and mat[ii,jj] is nonzero
                row_is_zero = 0;
        }
        //printf("i = %ld - j = %ld - nonz search: ii = %ld - jj = %ld\n", i, j, ii, jj);
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
                if (! nmod_poly_is_zero(nmod_poly_mat_entry(mat, i, jj)))
                {
                    slong pi = pivot_row[jj];
                    if (pi >= 0) // a pivot has already been found at [pi,jj], for some pi < i
                    {
                        //printf("i = %ld - j = %ld ; no-new-pivot - pi = %ld - jj = %ld\n", i, j, pi, jj);
                        // use gcd between [pi,jj] and [i,jj] and apply the corresponding row-wise unimodular transformation
                        // between rows pi and i (see description of Bradley's algorithm), which puts a zero at [i,jj] and
                        // the gcd (updated pivot) at [pi,jj]
                        _apply_xgcd_transformation(g, u, v, pivg, nonzg, piv, nonz, mat, tsf, pi, i, jj);
                        // note: new pivot already normalized since xgcd ensures the gcd is monic

                        // reduce appropriate entries in hnf in column jj
                        //nmod_poly_mat_degree_matrix_print_pretty(mat);
                        for (ii = 0; ii < pi; ii++)
                            _reduce_against_pivot(mat, tsf, nmod_poly_mat_entry(mat, pi, jj), pi, jj, ii, entry, u, v);
                        //nmod_poly_mat_degree_matrix_print_pretty(mat);
                        //printf("finished\n");
                    }
                    else // pivot_row[jj] == -1, currently no pivot in column jj
                    {
                        //printf("i = %ld - j = %ld ; new-pivot\n", i, j);
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
                        _normalize_pivot(mat, tsf, nmod_poly_mat_entry(mat, ii, jj), ii, jj, entry);

                        // reduce entries [ii,jj+1:j] against current hnf
                        for (slong kk = jj+1; kk < j; kk++)
                        {
                            pi = pivot_row[kk];
                            if (pi >= 0) // otherwise, no pivot, nothing to do
                            {
                                // reduce mat[ii,kk] against pivot [pi,kk]
                                _reduce_against_pivot(mat, tsf, nmod_poly_mat_entry(mat, pi, kk),
                                                      pi, kk, ii, entry, u, v);
                            }
                        }
                        // reduce column jj of hnf against new pivot mat[ii,jj]
                        for (slong pii = 0; pii < ii; pii++)
                        {
                            // reduce pii,jj against the pivot [ii,jj]
                            _reduce_against_pivot(mat, tsf, nmod_poly_mat_entry(mat, pii, jj),
                                                  ii, jj, pii, entry, u, v);
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
                if (nmod_poly_is_zero(nmod_poly_mat_entry(mat, i, j)))
                {
                    //printf("i = %ld - j = %ld ; now in row_is_zero --> no new pivot\n", i, j);
                    // entry [i,j] is zero, swap row i and the last row (just
                    // to avoid wasting time re-checking it is zero), and keep
                    // the same i,j to process same submatrix again.
                    nmod_poly_mat_swap_rows(mat, NULL, i, mat->r -1);
                    if (tsf)
                        nmod_poly_mat_swap_rows(tsf, NULL, i, mat->r -1);
                }
                else
                {
                    //printf("i = %ld - j = %ld ; now in row_is_zero --> new pivot\n", i, j);
                    // update pivot_col and pivot_row
                    pivot_col[i] = j;
                    pivot_row[j] = i;
                    // entry [i,j] is nonzero, this is a new pivot
                    // normalize it
                    _normalize_pivot(mat, tsf, nmod_poly_mat_entry(mat, i, j), i, j, entry);
                    // reduce entries above it, i.e. ii,j against pivot i,j
                    for (ii = 0; ii < i; ii++)
                        _reduce_against_pivot(mat, tsf, nmod_poly_mat_entry(mat, i, j),
                                              i, j, ii, entry, u, v);
                    // increment both i and j
                    i++;
                    j++;
                }
            }
        }
        //printf("i = %ld - j = %ld - OUT degree matrix:\n", i-1, j-1);
        //nmod_poly_mat_degree_matrix_print_pretty(mat);
    }

    flint_free(pivot_col);
    flint_free(pivot_row);

    return i;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
