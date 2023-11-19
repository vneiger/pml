#include <NTL/lzz_pX.h>
#include <NTL/mat_lzz_p.h>

#include "lzz_p_extra.h"
#include "mat_lzz_p_extra.h"
#include "lzz_pX_CRT.h"
#include "structured_lzz_p.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* Mosaic Toeplitz matrices:                                  */
/* block matrix where each block is Toeplitz                  */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/


/*------------------------------------------------------------*/
/* dummy constructor                                          */
/*------------------------------------------------------------*/
mosaic_toeplitz_lzz_p::mosaic_toeplitz_lzz_p()
{
    n = m = nb = mb = 0;
}

/*------------------------------------------------------------*/
/* copies all data                                            */
/*------------------------------------------------------------*/
mosaic_toeplitz_lzz_p::mosaic_toeplitz_lzz_p(const Vec< Vec<toeplitz_lzz_p> > & init)
{
    data = init;
    nb = init.length();
    mb = init[0].length();
    
    n = 0;
    m = 0;
    
    for(long i = 0; i < nb; i++)
        n += init[i][0].NumRows();
    for (long j = 0; j < mb; j++)
        m += init[0][j].NumCols();
}


/*------------------------------------------------------------*/
/* getters                                                    */
/*------------------------------------------------------------*/
long mosaic_toeplitz_lzz_p::NumRows() const
{
    return n;
}

long mosaic_toeplitz_lzz_p::NumCols() const
{
    return m;
}

long mosaic_toeplitz_lzz_p::NumBlockRows() const
{
    return nb;
}

long mosaic_toeplitz_lzz_p::NumBlockCols() const
{
    return mb;
}

long mosaic_toeplitz_lzz_p::NumRows_of_block(long i) const
{
    return data[i][0].NumRows();
}

long mosaic_toeplitz_lzz_p::NumCols_of_block(long i) const
{
    return data[0][i].NumCols();
}


/*------------------------------------------------------------*/
/* turns M into a dense matrix                                */
/*------------------------------------------------------------*/
void mosaic_toeplitz_lzz_p::to_dense(Mat<zz_p>& Mdense) const
{
    Mdense.SetDims(n, m);

    for (long i = 0; i < n; i++)
        for (long j = 0; j < m; j++)
            Mdense[i][j] = (*this)(i, j);
}

/*------------------------------------------------------------*/
/* right multiplication                                       */
/*------------------------------------------------------------*/
void mosaic_toeplitz_lzz_p::mul_right(Vec<zz_p>& res, const Vec<zz_p>& input) const
{
    if (input.length() != m)
    {
        LogicError("Wrong size for mosaic_hankel_lzz_p right multiplication.");
    }

    if (&res == &input)
    {
        res = mul_right(input);
        return;
    }

    res.SetLength(n);
    for (long i = 0; i < n; i++)
        res[i] = 0;

    Vec<zz_p> tmp_in, tmp_out;

    long jdx = 0;
    for (long j = 0; j < mb; j++)
    {
        long ncols = data[0][j].NumCols();
        tmp_in.SetLength(ncols);
        for (long ell = 0; ell < ncols; ell++)
            tmp_in[ell] = input[jdx + ell];

        long idx = 0;
        for (long i = 0; i < nb; i++)
        {
            long nrows = data[i][0].NumRows();
            data[i][j].mul_right(tmp_out, tmp_in);
            for (long ell = 0; ell < nrows; ell++)
                res[idx + ell] += tmp_out[ell];
            idx += nrows;
        }
        jdx += ncols;
    }
}

/*------------------------------------------------------------*/
/* right matrix multiplication                                */
/*------------------------------------------------------------*/
void mosaic_toeplitz_lzz_p::mul_right(Mat<zz_p>& res, const Mat<zz_p>& input) const
{
    if (input.NumRows() != m)
    {
        LogicError("Wrong size for mosaic_hankel_lzz_p right matrix multiplication.");
    }

    if (&res == &input)
    {
        res = mul_right(input);
        return;
    }

    long p = input.NumCols();

    // check if we simply do a dense matrix product
    if ( (min(m, n) >= max(m, n)/2) && // close enough to a square matrix
         p >= min(m, n) / 4 )          // close enough to a matrix-matrix product
    {
        long t = type_of_prime();
        if ( (t == TYPE_FFT_PRIME && min(m, n) <= 100)  ||
             (t == TYPE_SMALL_PRIME && min(m, n) <= 1000) ||
             (t == TYPE_LARGE_PRIME && min(m, n) <= 300) )
        {
            Mat<zz_p> mosaic = to_dense();
            res = mosaic * input;
            return;
        }
    }

    res.SetDims(n, p);

    for (long s = 0; s < n; s++)
        for (long r = 0; r < p; r++)
            res[s][r] = 0;

    Vec<zz_p> tmp_in, tmp_out;

    for (long r = 0; r < p; r++)
    {
        long jdx = 0;
        for (long j = 0; j < mb; j++)
        {
            long ncols = data[0][j].NumCols();
            tmp_in.SetLength(ncols);
            for (long ell = 0; ell < ncols; ell++)
                tmp_in[ell] = input[jdx + ell][r];
            
            long idx = 0;
            for (long i = 0; i < nb; i++)
            {
                long nrows = data[i][0].NumRows();
                data[i][j].mul_right(tmp_out, tmp_in);
                for (long ell = 0; ell < nrows; ell++)
                    res[idx + ell][r] += tmp_out[ell];
                idx += nrows;
            }
            jdx += ncols;
        }
    }
}

/*------------------------------------------------------------*/
/* left multiplication                                        */
/*------------------------------------------------------------*/
void mosaic_toeplitz_lzz_p::mul_left(Vec<zz_p>& res, const Vec<zz_p>& input) const
{
    if (input.length() != n)
    {
        LogicError("Wrong size for mosaic_hankel_lzz_p left multiplication.");
    }

    if (&res == &input)
    {
        res = mul_left(input);
        return;
    }

    res.SetLength(m);
    for (long i = 0; i < m; i++)
        res[i] = 0;

    Vec<zz_p> tmp_in, tmp_out;

    long idx = 0;
    for (long i = 0; i < nb; i++)
    {
        long nrows = data[i][0].NumRows();
        tmp_in.SetLength(nrows);

        for (long ell = 0; ell < nrows; ell++)
            tmp_in[ell] = input[idx + ell];

        long jdx = 0;
        for (long j = 0; j < mb; j++)
        {
            long ncols = data[0][j].NumCols();
            data[i][j].mul_left(tmp_out, tmp_in);
            for (long ell = 0; ell < ncols; ell++)
                res[jdx + ell] += tmp_out[ell];
            jdx += ncols;
        }
        idx += nrows;
    }
}


/*------------------------------------------------------------*/
/* left matrix multiplication                                 */
/*------------------------------------------------------------*/
void mosaic_toeplitz_lzz_p::mul_left(Mat<zz_p>& res, const Mat<zz_p>& input) const
{
    if (input.NumCols() != n)
    {
        LogicError("Wrong size for mosaic_hankel_lzz_p left matrix multiplication.");
    }

    if (&res == &input)
    {
        res = mul_left(input);
        return;
    }

    long p = input.NumRows();

    // check if we simply do a dense matrix product
    if ( (min(m, n) >= max(m, n)/2) && // close enough to a square matrix
         p >= min(m, n) / 4 )          // close enough to a matrix-matrix product
    {
        long t = type_of_prime();
        if ( (t == TYPE_FFT_PRIME && min(m, n) <= 100)  ||
             (t == TYPE_SMALL_PRIME && min(m, n) <= 1000) ||
             (t == TYPE_LARGE_PRIME && min(m, n) <= 300) )
        {
            Mat<zz_p> mosaic = to_dense();
            res = input * mosaic;
            return;
        }
    }

    res.SetDims(p, m);

    for (long s = 0; s < p; s++)
        for (long r = 0; r < m; r++)
            res[s][r] = 0;

    Vec<zz_p> tmp_in, tmp_out;

    for (long r = 0; r < p; r++)
    {
        long idx = 0;
        for (long i = 0; i < nb; i++)
        {
            long nrows = data[i][0].NumRows();
            tmp_in.SetLength(nrows);
            
            for (long ell = 0; ell < nrows; ell++)
                tmp_in[ell] = input[r][idx + ell];
            
            long jdx = 0;
            for (long j = 0; j < mb; j++)
            {
                long ncols = data[0][j].NumCols();
                data[i][j].mul_left(tmp_out, tmp_in);
                for (long ell = 0; ell < ncols; ell++)
                    res[r][jdx + ell] += tmp_out[ell];
                jdx += ncols;
            }
            idx += nrows;
        }
    }
}


/*------------------------------------------------------------*/
/* access to particular rows and columns                      */
/*------------------------------------------------------------*/
void mosaic_toeplitz_lzz_p::first_column_of_block(Vec<zz_p>& res, long i) const
{
    res.SetLength(n);

    long ind = 0;
    for (long r = 0; r < nb; r++)
        for (long j = 0; j < NumRows_of_block(r); j++)
            res[ind++] = data[r][i](j, 0);
}

/*------------------------------------------------------------*/
/* access to particular rows and columns                      */
/*------------------------------------------------------------*/
void mosaic_toeplitz_lzz_p::last_column_of_block(Vec<zz_p>& res, long i) const
{
    res.SetLength(n);

    long ind = 0;
    long shift = NumCols_of_block(i)-1;
    for (long r = 0; r < nb; r++)
        for (long j = 0; j < NumRows_of_block(r); j++)
            res[ind++] = data[r][i](j, shift);

}

/*------------------------------------------------------------*/
/* access to particular rows and columns                      */
/*------------------------------------------------------------*/
void mosaic_toeplitz_lzz_p::first_row_of_block(Vec<zz_p>& res, long i) const
{
    res.SetLength(m);

    long ind = 0;
    for (long r = 0; r < mb; r++)
        for (long j = 0; j < NumCols_of_block(r); j++)
            res[ind++] = data[i][r](0, j);
}

/*------------------------------------------------------------*/
/* access to particular rows and columns                      */
/*------------------------------------------------------------*/
void mosaic_toeplitz_lzz_p::last_row_of_block(Vec<zz_p>& res, long i) const
{
    res.SetLength(m);

    long ind = 0;
    long shift = NumRows_of_block(i)-1;
    for (long r = 0; r < mb; r++)
        for (long j = 0; j < NumCols_of_block(r); j++)
            res[ind++] = data[i][r](shift, j);
}

/*------------------------------------------------------------*/
/* G, H such that Z1 M - M Z0 = G H^t                         */
/*------------------------------------------------------------*/
void mosaic_toeplitz_lzz_p::phi_plus_generators(Mat<zz_p>& G, Mat<zz_p>& H) const
{
    long alpha_r = NumBlockRows();
    long alpha_c = NumBlockCols();
    long alpha = alpha_r + alpha_c;

    G.SetDims(n, alpha);
    H.SetDims(m, alpha);

    Vec<zz_p> tmp_row, tmp_col;
    long idx;

    for (long i = 0; i < alpha_c; i++)
    {
        last_column_of_block(tmp_col, i);
        G[0][i] = tmp_col[n-1];
        for (long j = 1; j < n; j++)
            G[j][i] = tmp_col[j-1];

        if (i < alpha_c-1)
        {
            first_column_of_block(tmp_col, i+1);
            for (long j = 0; j < n; j++)
                G[j][i] = G[j][i] - tmp_col[j];
        }
    }

    idx = - 1;
    for (long i = 0; i < alpha_c; i++)
    {
        idx += NumCols_of_block(i);
        for (long j = 0; j < m; j++)
            H[j][i] = 0;
        H[idx][i] = 1;
    }

    idx = 0;
    for (long i = 0; i < alpha_r; i++)
    {
        for (long j = 0; j < n; j++)
            G[j][i+alpha_c] = 0;
        G[idx][i+alpha_c] = 1;
        idx += NumRows_of_block(i);
    }

    for (long i = 0; i < alpha_r; i++)
    {
        if (i == 0)
            last_row_of_block(tmp_row, alpha_r-1);
        else
            last_row_of_block(tmp_row, i-1);
        for (long j = 0; j < m; j++)
            H[j][i+alpha_c] = tmp_row[j];
        first_row_of_block(tmp_row, i);
        for (long j = 0; j < m-1; j++)
            H[j][i+alpha_c] = H[j][i+alpha_c] - tmp_row[j+1];
        long jdx = -1;
        for (long j = 0; j < alpha_c; j++)
        {
            jdx += NumCols_of_block(j);
            H[jdx][i+alpha_c] = 0;
        }
    }
}

/*------------------------------------------------------------------*/
/* finds c such that                                                */
/* - c != 0                                                         */
/* - a^i - c a^j != 0 for 0 <= i < n and 0 <= j < m                 */
/*------------------------------------------------------------------*/
static void find_c(zz_p& c, const zz_p& a, long n, long m)
{
    Vec<zz_p> pow_a, pow_inva;
    pow_a.SetLength(n);
    pow_inva.SetLength(m);

    pow_a[0] = to_zz_p(1);
    pow_inva[0] = to_zz_p(1);

    zz_p inva = 1/a;
    for (long i = 1; i < n; i++)
        pow_a[i] = a*pow_a[i-1];
    for (long i = 1; i < m; i++)
        pow_inva[i] = inva*pow_inva[i-1];

    bool done;
    do
    {
        c = random_zz_p();
        done = true;
        if (c == 0)
            done = false;
        for (long i = 0; i < n; i++)
            if (c == pow_a[i])
                done = false;
        for (long i = 0; i < m; i++)
            if (c == pow_inva[i])
                done = false;
    } 
    while (done != true);
}



/*------------------------------------------------------------------*/
/* preconditions M                                                  */
/* builds the matrix CL = (D_e X_int) M J (D_f Y_int)^t, where:     */
/* - X_int, Y_int are geometric interpolation                       */
/* - D_e, D_f are diagonal matrix built on vectors e and f          */
/* - CL is cauchy-like geometric                                    */
/* - CL is expected to have generic rank profile                    */
/* - J is the reversal matrix of size m                             */
/* return value is 0 if field to small for preconditioning          */
/*------------------------------------------------------------------*/
static long to_cauchy_grp(cauchy_like_geometric_lzz_p& CL, 
                          zz_pX_Multipoint_Geometric& X_int, zz_pX_Multipoint_Geometric& Y_int,
                          Vec<zz_p> &e, Vec<zz_p> &f,
                          const mosaic_toeplitz_lzz_p& M)
{
    Mat<zz_p> X, Y;
    Mat<zz_p> G, H;
    M.phi_plus_generators(G, H);
    zz_p a, b, c;
    long n = M.NumRows();
    long m = M.NumCols();

    // cout << (Z_lzz_p(n, to_zz_p(1))*M.to_dense() - M.to_dense()*Z_lzz_p(m, to_zz_p(0)) ==  G*transpose(H)) << endl;

    if (max(m, n) > zz_p::modulus()/ 10)
    {
        return 0;
    }

    element_of_order(a, 2 * max(m, n));
    b = a*a;
    find_c(c, b, n, m);
    X_int = zz_pX_Multipoint_Geometric(a, to_zz_p(1), n);
    Y_int = zz_pX_Multipoint_Geometric(a, c, m);

    long alpha = G.NumCols();
    X.SetDims(n, alpha+2);
    Y.SetDims(m, alpha+2);

    e.SetLength(n);
    for (long i = 0; i < n; i++)
        e[i] = random_zz_p();

    f.SetLength(m);
    for (long i = 0; i < m; i++)
        f[i] = random_zz_p();


    // gens 0 to alpha-1
    // De X_int G
    Vec<zz_p> tmp_v;
    for (long j = 0; j < alpha; j++)
    {
        zz_pX tmp_p;
        tmp_p.rep.SetLength(n);
        zz_p* coef_p = tmp_p.rep.elts();
        for (long i = 0; i < n; i++)
        {
            coef_p[i] = G[i][j];
        }
        tmp_p.normalize();
        X_int.evaluate(tmp_v, tmp_p);
        for (long i = 0; i < n; i++)
        {
            X[i][j] = tmp_v[i] * e[i];
        }
    }

    // Df Y_int H
    Vec<zz_p> tmp_w;
    for (long j = 0; j < alpha; j++)
    {
        zz_pX tmp_q;
        tmp_q.rep.SetLength(m);
        zz_p* coef_q = tmp_q.rep.elts();
        for (long i = 0; i < m; i++)
        {
            coef_q[i] = H[m - 1 - i][j];
        }
        tmp_q.normalize();
        Y_int.evaluate(tmp_w, tmp_q);
        for (long i = 0; i < m; i++)
        {
            Y[i][j] = tmp_w[i] * f[i];
        }
    }

    // gen alpha
    zz_p tmp_z = to_zz_p(1);
    for (long i = 0; i < n; i++)
    {
        X[i][alpha] = e[i]*(power(tmp_z, n)-1);
        tmp_z = tmp_z * b;
    }
    M.last_row_of_block(tmp_w, M.NumBlockRows()-1);
    zz_pX tmp_q;
    tmp_q.rep.SetLength(m);
    zz_p* coef_q = tmp_q.rep.elts();
    for (long i = 0; i < m; i++)
    {
        coef_q[i] = tmp_w[m - 1 - i];
    }
    tmp_q.normalize();
    Y_int.evaluate(tmp_w, tmp_q);
    for (long i = 0; i < m; i++)
    {
        Y[i][alpha] = tmp_w[i] * f[i];
    }

    // gen alpha+1
    M.first_column_of_block(tmp_v, 0);
    zz_pX tmp_p;
    tmp_p.rep.SetLength(n);
    zz_p* coef_p = tmp_p.rep.elts();
    for (long i = 0; i < n; i++)
    {
        coef_p[i] = tmp_v[i];
    }
    tmp_p.normalize();
    X_int.evaluate(tmp_v, tmp_p);
    for (long i = 0; i < n; i++)
    {
        X[i][alpha+1] = tmp_v[i] * e[i];
    }
    tmp_z = c;
    for (long i = 0; i < m; i++)
    {
        Y[i][alpha+1] = -f[i]*power(tmp_z, m);
        tmp_z = tmp_z * b;
    }

    CL = cauchy_like_geometric_lzz_p(X, Y, to_zz_p(1), c, b);

    // cout << ( CL.to_dense() == 
    //           diagonal_matrix(e) * X_int.to_dense() * M.to_dense() * J_lzz_p(m) * transpose(diagonal_matrix(f) * Y_int.to_dense()) ) << endl;

    return 1;
}


/*------------------------------------------------------------*/
/* finds a random solution to M.x = b                         */
/* if return value is 0, something went wrong (no grp)        */
/* else, return value is 1. Then x is empty if no solution    */
/* thresh is threshold for divide-and-conquer                 */
/* thresh_alpha switches between block and plain quadratic    */
/*------------------------------------------------------------*/
long mosaic_toeplitz_lzz_p::solve(Vec<zz_p>& x, const Vec<zz_p> b, long thresh, long thresh_alpha) const
{
    cauchy_like_geometric_lzz_p CL;
    zz_pX_Multipoint_Geometric X, Y;
    Vec<zz_p> e, f;

    if (m*n == 0)
        LogicError("Empty mosaic toeplitz matrix in solve");

    long r = to_cauchy_grp(CL, X, Y, e, f, *this);

    if (r != 1)
        return 0;

    Vec<zz_p> xp, bp;
    X.mul_right(bp, b);
    for (long i = 0; i < n; i++)
        bp[i] = bp[i] * e[i];

    long s = CL.solve_grp(xp, bp, thresh, thresh_alpha);

    if (s != 1)
        return 0;
    
    for (long i = 0; i < m; i++)
        xp[i] = xp[i] * f[i];

    Y.mul_left(xp, xp);

    x.SetLength(m);
    for (long i = 0; i < m; i++)
        x[i] = xp[m - 1 - i];

    return 1;
}

/*------------------------------------------------------------*/
/* finds the inverse of this as a toeplitz_like_minus matrix  */
/* if return value is 0, either no grp or singular            */
/* else, return value is 1, and inv is the inverse            */
/* thresh is threshold for divide-and-conquer                 */
/* thresh_alpha switches between block and plain quadratic    */
/*------------------------------------------------------------*/
long mosaic_toeplitz_lzz_p::inv(toeplitz_like_minus_lzz_p& inv, long thresh, long thresh_alpha) const
{
    cauchy_like_geometric_lzz_p CL, iCL;
    zz_pX_Multipoint_Geometric X, Y;
    Vec<zz_p> e, f, x, b, bp;
    Mat<zz_p> G, H, iG, iH;
    long alpha, r, s, nn;

    if (m != n)
        Error("Attempt to invert a non-square matrix");

    r = to_cauchy_grp(CL, X, Y, e, f, *this);
    if (r != 1)   // field too small
        return 0;

    s = CL.invert_leading_principal_minor_grp(iCL, thresh, thresh_alpha);
    if (s != 1)   // not generic rank profile  
        return 0;

    nn = iCL.NumRows();
    if (nn != n)  // not full rank
        return 0;

    phi_plus_generators(G, H); // TODO: they are already computed in to_cauchy_grp
    alpha = G.NumCols();

    iG.SetDims(n, alpha);
    iH.SetDims(n, alpha);
    b.SetLength(n);
    
    // first generator is -M^(-1) * G
    for (long j = 0; j < alpha; j++)
    {
        for (long i = 0; i < n; i++)
            b[i] = G[i][j];

        X.mul_right(bp, b);

        for (long i = 0; i < n; i++)
            bp[i] = bp[i] * e[i];
        
        iCL.mul_right(x, bp);
    
        for (long i = 0; i < n; i++)
            x[i] = - (x[i] * f[i]);

        Y.mul_left(b, x);

        for (long i = 0; i < n; i++)
            iG[i][j] = b[n - 1 - i];
    }

    // second generator is M^(-t) H
    for (long j = 0; j < alpha; j++)
    {
        for (long i = 0; i < n; i++)
            b[i] = H[n - 1 - i][j];

        Y.mul_right(bp, b);

        for (long i = 0; i < n; i++)
            bp[i] = bp[i] * f[i];
        
        iCL.mul_left(x, bp);
    
        for (long i = 0; i < n; i++)
            x[i] = x[i] * e[i];

        X.mul_left(b, x);

        for (long i = 0; i < n; i++)
            iH[i][j] = b[i];
    }

    inv = toeplitz_like_minus_lzz_p(iG, iH);

    return 1;
}



// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
