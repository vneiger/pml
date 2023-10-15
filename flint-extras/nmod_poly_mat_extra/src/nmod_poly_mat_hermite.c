#include <flint/flint.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_vec.h>
#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_forms.h"
#include "nmod_poly_mat_io.h" // TODO remove, for debugging

// TODO how to handle upper/lower?
// TODO how to handle rowwise/colwise?
// TODO restrict to rank rows, where to do this?
//
//
// For all:
// - upper Hermite, row-wise
// - computation is in place
// - concerning transformation:
//    -- if tsf == NULL, transformation is not computed at all (which saves computation time)
//    -- otherwise, the transformation is accumulated into tsf
//         --> if wanting the transformation for this computation only, set tsf to identity before calling this

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
    // -> there cannot be > n pivots
    // -> index i below will tell us how many have been found already
    slong * pivot_cols = flint_malloc(n * sizeof(slong));

    // (rk,j) : position where next pivot is expected to be found when processing column j,
    // =>   rk = number of pivots already found in columns < j, which is also the rank
    // of these first j-1 columns
    // Note: in case of zeroes/dependencies, we may find a pivot below rk,
    // then we will permute rows to bring it up at row rk
    slong rk = 0;

    slong len1, len2, len; // will store length of polynomials
    mp_limb_t cst; // will store constant multipliers (quotients of leading coefficients)
    nmod_poly_struct * entry1; // will point to some matrix entry
    nmod_poly_struct * entry2; // will point to some matrix entry

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
                pivot_cols[rk] = j;
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
                    pivot_cols[rk] = j;
                    rk++;
                    collision = 0;
                }
                else
                {
                    // use pi2,j to kill the leading term of pi1,j, and apply the
                    // corresponding transformation to the whole row pi1
                    //  --> mat[pi1,:] = mat[pi1,:] + cst * x**k * mat[pi2,:]
                    // where cst = - leading_coeff(mat[pi1,j]) / leading_coeff(mat[pi2,j])
                    //   and k = deg(mat[pi1,j]) - deg(mat[pi2,j])
                    entry1 = nmod_poly_mat_entry(mat, pi1, j);
                    entry2 = nmod_poly_mat_entry(mat, pi2, j);
                    len1 = entry1->length;
                    len2 = entry2->length;
                    cst = n_invmod(entry2->coeffs[len2-1], mat->modulus);
                    cst = nmod_mul(entry1->coeffs[len1-1], cst, entry1->mod);
                    cst = nmod_neg(cst, entry1->mod);
                    _nmod_vec_scalar_addmul_nmod(entry1->coeffs+len1-len2, entry2->coeffs, entry2->length, cst, entry1->mod);
                    _nmod_poly_normalise(entry1);
                    for (slong jj = j+1; jj < n; jj++)
                    {
                        entry1 = nmod_poly_mat_entry(mat, pi1, jj);
                        entry2 = nmod_poly_mat_entry(mat, pi2, jj);
                        len = entry2->length + len1-len2;
                        nmod_poly_fit_length(entry1, len);
                        if (len > entry1->length)
                        {
                            _nmod_vec_zero(entry1->coeffs + entry1->length, len - entry1->length);
                            _nmod_poly_set_length(entry1, len);
                        }
                        _nmod_vec_scalar_addmul_nmod(entry1->coeffs+len1-len2, entry2->coeffs, entry2->length, cst, entry1->mod);
                        _nmod_poly_normalise(entry1);
                    }
                    if (tsf) for (slong jj = 0; jj < m; jj++)
                    {
                        entry1 = nmod_poly_mat_entry(tsf, pi1, jj);
                        entry2 = nmod_poly_mat_entry(tsf, pi2, jj);
                        len = entry2->length + len1-len2;
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
        slong j = pivot_cols[i];
        entry2 = nmod_poly_mat_entry(mat, i, j); // pivot entry

        // normalize pivot, and correspondingly update all the row
        cst = n_invmod(entry2->coeffs[entry2->length - 1], entry2->mod.n);
        for (slong jj = j; jj < n; jj++)
        {
            entry1 = nmod_poly_mat_entry(mat, i, jj);
            _nmod_vec_scalar_mul_nmod(entry1->coeffs, entry1->coeffs, entry1->length, cst, entry1->mod);
        }
        if (tsf)
        {
            for (slong jj = 0; jj < m; jj++)
            {
                entry1 = nmod_poly_mat_entry(tsf, i, jj);
                _nmod_vec_scalar_mul_nmod(entry1->coeffs, entry1->coeffs, entry1->length, cst, entry1->mod);
            }
        }

        // reduce entries above (i,j)
        for (slong ii = 0; ii < i; ii++)
        {
            entry1 = nmod_poly_mat_entry(mat, ii, j);
            if (entry1->length >= entry2->length)
            {
                // division with remainder: entry1 = entry2 * u + v
                nmod_poly_divrem(u, v, entry1, entry2);
                // replace entry1 by remainder
                nmod_poly_swap(entry1, v);
                // apply same transformation on rest of the row: row_ii <- row_ii - u * row_i
                nmod_poly_neg(u, u);
                for (slong jj = j+1; jj < n; jj++)
                {
                    nmod_poly_mul(v, u, nmod_poly_mat_entry(mat, i, jj)); // v used as tmp
                    nmod_poly_add(nmod_poly_mat_entry(mat, ii, jj), nmod_poly_mat_entry(mat, ii, jj), v);
                }
                if (tsf)
                {
                    for (slong jj = 0; jj < m; jj++)
                    {
                        nmod_poly_mul(v, u, nmod_poly_mat_entry(tsf, i, jj)); // v used as tmp
                        nmod_poly_add(nmod_poly_mat_entry(tsf, ii, jj), nmod_poly_mat_entry(tsf, ii, jj), v);
                    }
                }
            }
        }
    }

    flint_free(pivot_cols);

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
    // -> there cannot be > n pivots
    // -> index i below will tell us how many have been found already
    slong * pivot_cols = flint_malloc(n * sizeof(slong));

    // (rk,j) : position where next pivot is expected to be found when processing column j,
    // =>   rk = number of pivots already found in columns < j, which is also the rank
    // of these first j-1 columns
    // Note: in case of zeroes/dependencies, we may find a pivot below rk,
    // then we will permute rows to bring it up at row rk
    slong rk = 0;

    mp_limb_t inv; // will store inverses of leading coefficients for making monic
    nmod_poly_struct * piv; // will point to pivot entry
    nmod_poly_struct * nonz; // will point to nonpivot entry that is bound to become zero
    nmod_poly_struct * entry; // will point to some entry
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
                // else, pivot is the only nonzero entry: nothing more to do for this column
                {
                    piv = nmod_poly_mat_entry(mat, pi, j); // pivot
                    nonz = nmod_poly_mat_entry(mat, ii, j); // nonzero below pivot
                    nmod_poly_xgcd(g, u, v, piv, nonz); // g = u*piv + v*nonz
                    nmod_poly_divides(pivg, piv, g); // pivg = piv // g
                    nmod_poly_divides(nonzg, nonz, g); // nonzg = nonz // g

                    // cancel (ii,j) and modify row ii accordingly
                    // in short: for M in {mat,tsf} do the unimodular transformation on rows pi,ii:
                    //  [ M[pi,:] ]  =  [    u       v  ]  *  [ M[pi,:] ]
                    //  [ M[ii,:] ]     [ -nonzg   pivg ]     [ M[ii,:] ]

                    // for mat, due to zeroes already set in columns 0...j-1, we can start at column j
                    // set piv to g, and later below set nonz to 0
                    nmod_poly_set(piv, g); // mat[pi,j] = g
                    for (slong jj = j+1; jj < n; jj++)
                    {
                        // simultaneously update:
                        //     mat[pi,jj] = u * mat[pi,jj] + v * mat[ii,jj]
                        //     mat[ii,jj] = -nonzg * mat[pi,jj] + pivg * mat[ii,jj]
                        // --> use g as temporary copy of mat[pi,jj]
                        // --> use nonz as temporary for storing products
                        nmod_poly_set(g, nmod_poly_mat_entry(mat, pi, jj));

                        nmod_poly_mul(nmod_poly_mat_entry(mat, pi, jj), u, g);
                        nmod_poly_mul(nonz, v, nmod_poly_mat_entry(mat, ii, jj));
                        nmod_poly_add(nmod_poly_mat_entry(mat, pi, jj), nmod_poly_mat_entry(mat, pi, jj), nonz);

                        nmod_poly_mul(nonz, pivg, nmod_poly_mat_entry(mat, ii, jj));
                        nmod_poly_mul(nmod_poly_mat_entry(mat, ii, jj), nonzg, g);
                        nmod_poly_sub(nmod_poly_mat_entry(mat, ii, jj), nonz, nmod_poly_mat_entry(mat, ii, jj));
                    }
                    if (tsf)
                    {
                        for (slong jj = 0; jj < m; jj++)
                        {
                            // simultaneously update:
                            //     tsf[pi,jj] = u * tsf[pi,jj] + v * tsf[ii,jj]
                            //     tsf[ii,jj] = -nonzg * tsf[pi,jj] + pivg * tsf[ii,jj]
                            // --> use g as temporary copy of tsf[pi,jj]
                            // --> use nonz as temporary for storing products
                            nmod_poly_set(g, nmod_poly_mat_entry(tsf, pi, jj));

                            nmod_poly_mul(nmod_poly_mat_entry(tsf, pi, jj), u, g);
                            nmod_poly_mul(nonz, v, nmod_poly_mat_entry(tsf, ii, jj));
                            nmod_poly_add(nmod_poly_mat_entry(tsf, pi, jj), nmod_poly_mat_entry(tsf, pi, jj), nonz);

                            nmod_poly_mul(nonz, pivg, nmod_poly_mat_entry(tsf, ii, jj));
                            nmod_poly_mul(nmod_poly_mat_entry(tsf, ii, jj), nonzg, g);
                            nmod_poly_sub(nmod_poly_mat_entry(tsf, ii, jj), nonz, nmod_poly_mat_entry(tsf, ii, jj));
                        }
                    }
                    nmod_poly_zero(nonz); // mat[ii,j] = 0
                }
            }

            // if pivot found in a lower row than expected
            // -> swap rows to bring it up to row rk
            nmod_poly_mat_swap_rows(mat, NULL, rk, pi); // does nothing if rk==pi
            if (tsf)
                nmod_poly_mat_swap_rows(tsf, NULL, rk, pi);
            pivot_cols[rk] = j;
            rk++;
        }
    }

    // now, we have mat in upper echelon form, it remains to
    // 1. reduce the entries above each pivot
    // 2. normalize each pivot
    for (slong i = 0; i < rk; i++)
    {
        slong j = pivot_cols[i];
        piv = nmod_poly_mat_entry(mat, i, j);

        // normalize pivot, and correspondingly update all the row
        inv = n_invmod(piv->coeffs[piv->length - 1], piv->mod.n);
        for (slong jj = j; jj < n; jj++)
        {
            entry = nmod_poly_mat_entry(mat, i, jj);
            _nmod_vec_scalar_mul_nmod(entry->coeffs, entry->coeffs, entry->length, inv, entry->mod);
        }
        if (tsf)
        {
            for (slong jj = 0; jj < m; jj++)
            {
                entry = nmod_poly_mat_entry(tsf, i, jj);
                _nmod_vec_scalar_mul_nmod(entry->coeffs, entry->coeffs, entry->length, inv, entry->mod);
            }
        }

        // reduce entries above (i,j)
        for (slong ii = 0; ii < i; ii++)
        {
            entry = nmod_poly_mat_entry(mat, ii, j);
            if (entry->length >= piv->length)
            {
                // division with remainder: entry = piv * u + v
                nmod_poly_divrem(u, v, entry, piv);
                // replace entry by remainder
                nmod_poly_swap(entry, v);
                // apply same transformation on rest of the row: row_ii <- row_ii - u * row_i
                nmod_poly_neg(u, u);
                for (slong jj = j+1; jj < n; jj++)
                {
                    nmod_poly_mul(v, u, nmod_poly_mat_entry(mat, i, jj)); // v used as tmp
                    nmod_poly_add(nmod_poly_mat_entry(mat, ii, jj), nmod_poly_mat_entry(mat, ii, jj), v);
                }
                if (tsf)
                {
                    for (slong jj = 0; jj < m; jj++)
                    {
                        nmod_poly_mul(v, u, nmod_poly_mat_entry(tsf, i, jj)); // v used as tmp
                        nmod_poly_add(nmod_poly_mat_entry(tsf, ii, jj), nmod_poly_mat_entry(tsf, ii, jj), v);
                    }
                }
            }
        }
    }

    flint_free(pivot_cols);
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
}
// at start of iteration i, the (i-1)x(i-1) leading principal submatrix is in HNF. 
// --> first use gcd on rows (1,i), (2,i), ... (i-1,i) to reduce the first i-1 entries of row i
// --> then reduce the entries above the diagonal in i-th column

slong nmod_poly_mat_upper_hermite_form_rowwise_domich(nmod_poly_mat_t mat, nmod_poly_mat_t tsf);

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
