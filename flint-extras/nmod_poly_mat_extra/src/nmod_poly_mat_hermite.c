#include <flint/nmod_vec.h>
#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_forms.h"
#include "nmod_poly_mat_io.h" // TODO remove, for debugging

// TODO how to handle upper/lower?
// TODO how to handle rowwise/colwise?
// TODO restrict to rank rows?
// Computation is in place
// Transformation:
// - if tsf == NULL, transformation is not computed at all (which saves computation time)
// - otherwise, the transformation is accumulated into tsf
//      --> if wanting the transformation for this computation only, set tsf to identity before calling this
slong nmod_poly_mat_upper_hermite_form_rowwise_bradley(nmod_poly_mat_t mat, nmod_poly_mat_t tsf)
{
    const ulong m = mat->r;
    const ulong n = mat->c;

    if (m == 0 || n == 0)
        return 0;

    // to store the positions of columns where a pivot was found
    // -> there cannot be > n pivots
    // -> index i below will tell us how many have been found already
    ulong * pivot_cols = flint_malloc(n * sizeof(ulong));

    // (rk,j) : position where next pivot is expected to be found when processing column j,
    // =>   rk = number of pivots already found in columns < j
    // In particular, at all times rk <= j
    // Note: in case of zeroes/dependencies, we may find a pivot below rk,
    // then we will permute rows to bring it up at row rk
    ulong rk = 0;
    ulong j = 0;

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
    while (j < n)
    {
        // find actual pivot: starting at (i,j), find first nonzero entry in the rest of column j
        ulong pi = rk;
        while (pi < m && nmod_poly_is_zero(nmod_poly_mat_entry(mat, pi, j)))
            pi++;

        if (pi < m)
        // else, no pivot in this column, go to next column
        {
            // process column j below pivot
            ulong ii = pi + 1;
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
                    for (ulong jj = j+1; jj < n; jj++)
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
                        for (ulong jj = 0; jj < m; jj++)
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
        j++;
    }

    // now, we have mat in upper echelon form, it remains to
    // 1. reduce the entries above each pivot
    // 2. normalize each pivot
    for (ulong i = 0; i < rk; i++)
    {
        j = pivot_cols[i];
        piv = nmod_poly_mat_entry(mat, i, j);

        // normalize pivot, and correspondingly update all the row
        inv = n_invmod(piv->coeffs[piv->length - 1], piv->mod.n);
        for (ulong jj = j; jj < n; jj++)
        {
            entry = nmod_poly_mat_entry(mat, i, jj);
            _nmod_vec_scalar_mul_nmod(entry->coeffs, entry->coeffs, entry->length, inv, entry->mod);
        }
        if (tsf)
        {
            for (ulong jj = 0; jj < m; jj++)
            {
                entry = nmod_poly_mat_entry(tsf, i, jj);
                _nmod_vec_scalar_mul_nmod(entry->coeffs, entry->coeffs, entry->length, inv, entry->mod);
            }
        }

        // reduce entries above (i,j)
        for (ulong ii = 0; ii < i; ii++)
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
                for (ulong jj = j+1; jj < n; jj++)
                {
                    nmod_poly_mul(v, u, nmod_poly_mat_entry(mat, i, jj)); // v used as tmp
                    nmod_poly_add(nmod_poly_mat_entry(mat, ii, jj), nmod_poly_mat_entry(mat, ii, jj), v);
                }
                if (tsf)
                {
                    for (ulong jj = 0; jj < m; jj++)
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

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
