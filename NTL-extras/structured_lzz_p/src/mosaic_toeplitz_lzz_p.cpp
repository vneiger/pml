#include <NTL/lzz_pX.h>
#include <NTL/mat_lzz_p.h>

#include "lzz_p_extra.h"
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

/* getters                                            */
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

/* turns M into a dense matrix                        */
/*------------------------------------------------------------*/

void mosaic_toeplitz_lzz_p::to_dense(Mat<zz_p>& Mdense) const
{
    Mdense.SetDims(n, m);

    for (long i = 0; i < n; i++)
        for (long j = 0; j < m; j++)
            Mdense[i][j] = (*this)(i, j);
}

/*------------------------------------------------------------*/

/* right multiplication                               */
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

/* right matrix multiplication                        */
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

/* left multiplication                                */
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

/* left matrix multiplication                         */
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

/* access to particular rows and columns              */
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

/* access to particular rows and columns              */
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

/* access to particular rows and columns              */
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

/* access to particular rows and columns              */
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

/* G, H such that Z1 M - M Z0^t = G H^t               */
/*------------------------------------------------------------*/

static void generators(Mat<zz_p>& G, Mat<zz_p>& H, const mosaic_toeplitz_lzz_p& M)
{
    // long n = M.NumRows();
    // long m = M.NumCols();
    // long alpha_r = M.NumBlockRows();
    // long alpha_c = M.NumBlockCols();
    // long alpha = alpha_r + alpha_c;

    // G.SetDims(n, alpha);
    // H.SetDims(m, alpha);

    // Vec<zz_p> tmp_row, tmp_col;
    // long idx;

    // idx = 0;
    // for (long i = 0; i < alpha_c; i++)
    // {
    //     M.first_column_of_block(tmp_col, i);
    //     G[0][i] = tmp_col[n-1];
    //     for (long j = 1; j < n; j++)
    //         G[j][i] = tmp_col[j-1];

    //     if (i > 0)
    //     {
    //         M.last_column_of_block(tmp_col, i-1);
    //         for (long j = 0; j < n; j++)
    //             G[j][i] = G[j][i] - tmp_col[j];
    //     }
    // }
    // for (long i = 0; i < alpha_r; i++)
    // {
    //     for (long j = 0; j < n; j++)
    //         G[j][i+alpha_c] = 0;
    //     G[idx][i+alpha_c] = 1;
    //     idx += M.NumRows_of_block(i);
    // }


    // idx = 0;
    // for (long i = 0; i < alpha_c; i++)
    // {
    //     for (long j = 0; j < m; j++)
    //         H[j][i] = 0;
    //     H[idx][i] = 1;
    //     idx += M.NumCols_of_block(i);
    // }

    // for (long i = 0; i < alpha_r; i++)
    // {
    //     if (i == 0)
    //         M.last_row_of_block(tmp_row, alpha_r-1);
    //     else
    //         M.last_row_of_block(tmp_row, i-1);
    //     for (long j = 1; j < m; j++)
    //         H[j][i+alpha_c] = tmp_row[j];
    //     M.first_row_of_block(tmp_row, i);
    //     for (long j = 1; j < m; j++)
    //         H[j][i+alpha_c] = H[j][i+alpha_c] - tmp_row[j-1];
    //     long jdx = 0;
    //     for (long j = 0; j < alpha_c; j++)
    //     {
    //         H[jdx][i+alpha_c] = 0;
    //         jdx += M.NumCols_of_block(j);
    //     }
    // }
}

// /*------------------------------------------------------------------*/
// /* finds c such that                                                */
// /* - c != 0                                                         */
// /* - a^i - c a^j != 0 for 0 <= i < n and 0 <= j < m                 */
// /*------------------------------------------------------------------*/
// static 
// void find_c(zz_p& c, const zz_p& a, long n, long m){
//     Vec<zz_p> pow_a, pow_inva;
//     pow_a.SetLength(n);
//     pow_inva.SetLength(m);

//     pow_a[0] = to_zz_p(1);
//     pow_inva[0] = to_zz_p(1);

//     zz_p inva = 1/a;
//     for (long i = 1; i < n; i++)
//         pow_a[i] = a*pow_a[i-1];
//     for (long i = 1; i < m; i++)
//         pow_inva[i] = inva*pow_inva[i-1];

//     bool done;
//     do{
//         c = random_zz_p();
//         done = true;
//         if (c == 0)
//             done = false;
//         for (long i = 0; i < n; i++)
//             if (c == pow_a[i])
//                 done = false;
//         for (long i = 0; i < m; i++)
//             if (c == pow_inva[i])
//                 done = false;
//     } while (done != true);
// }

// /*------------------------------------------------------------------*/
// /* preconditions M                                                  */
// /* builds the matrix CL = (D_e X_int) M (D_f Y_int)^t, where:       */
// /* - X_int, Y_int are geometric interpolation                       */
// /* - D_e, D_f are diagonal matrix built on vectors e and f          */
// /* - CL is cauchy-like special                                      */
// /* - CL is expected to have generic rank profile                    */
// /*                                                                  */
// /* - X_int is built on (1,w,w^2,..)                                 */
// /* - Y_int is built on (c,cw,cw^2,..)                               */
// /* If we are over an FFT prime, w = root of unity                   */
// /*------------------------------------------------------------------*/
// void to_cauchy_grp(cauchy_like_geometric_special& CL, 
//                    zz_pX_Multipoint_Geometric& X_int, zz_pX_Multipoint_Geometric& Y_int,
//                    Vec<zz_p> &e, Vec<zz_p> &f,
//                    const mosaic_hankel_lzz_p& M){


//     Mat<zz_p> X, Y;
//     Mat<zz_p> G, H;
//     generators(G, H, M);
//     cauchy_geometric_special C;
//     zz_p a, b, c;
//     long n = M.NumRows();
//     long m = M.NumCols();

//     if (max(m, n) > zz_p::modulus())
//         LogicError("Field too small for preconditioning.");

//     // zz_pInfoT *info = zz_pInfo;
//     // if (info->p_info != NULL){  // FFT prime
//     //   long k = NextPowerOfTwo(max(m, n));
//     //   long order = 1L << k;
//     //   long w = find_root_of_unity(zz_p::modulus(), 2*order);
//     //   a = to_zz_p(w);
//     // }
//     // else
//     element_of_order(a, max(m, n));

//     b = a*a;
//     find_c(c, b, n, m);
//     X_int = zz_pX_Multipoint_Geometric(a, n, 0);
//     Y_int = zz_pX_Multipoint_Geometric(a, m, c);
//     C = cauchy_geometric_special(to_zz_p(1), c, b, n, m);

//     long alpha = G.NumCols();
//     X.SetDims(n, alpha+2);
//     Y.SetDims(m, alpha+2);

//     e.SetLength(n);
//     for (long i = 0; i < n; i++)
//         e[i] = random_zz_p();

//     f.SetLength(m);
//     for (long i = 0; i < m; i++)
//         f[i] = random_zz_p();

//     vec_zz_p tmp_v;
//     for (long j = 0; j < alpha; j++){
//         zz_pX tmp_p;
//         tmp_p.rep.SetLength(n);
//         zz_p* coef_p = tmp_p.rep.elts();
//         for (long i = 0; i < n; i++)
//             coef_p[i] = G[i][j];
//         tmp_p.normalize();
//         X_int.evaluate(tmp_v, tmp_p);
//         for (long i = 0; i < n; i++)
//             X[i][j] = tmp_v[i] * e[i];
//     }

//     zz_p tmp_z = to_zz_p(1);
//     for (long i = 0; i < n; i++){
//         X[i][alpha] = e[i]*(power(tmp_z, n)-1);
//         tmp_z = tmp_z * b;
//     }

//     last_column_of_block(tmp_v, M.NumBlockCols()-1, M);
//     zz_pX tmp_p;
//     tmp_p.rep.SetLength(n);
//     zz_p* coef_p = tmp_p.rep.elts();
//     for (long i = 0; i < n; i++)
//         coef_p[i] = tmp_v[i];
//     tmp_p.normalize();
//     X_int.evaluate(tmp_v, tmp_p);
//     for (long i = 0; i < n; i++)
//         X[i][alpha+1] = tmp_v[i] * e[i];

//     vec_zz_p tmp_w;
//     for (long j = 0; j < alpha; j++){
//         zz_pX tmp_q;
//         tmp_q.rep.SetLength(m);
//         zz_p* coef_q = tmp_q.rep.elts();
//         for (long i = 0; i < m; i++)
//             coef_q[i] = H[i][j];
//         tmp_q.normalize();
//         Y_int.evaluate(tmp_w, tmp_q);
//         for (long i = 0; i < m; i++)
//             Y[i][j] = tmp_w[i] * f[i];
//     }

//     last_row_of_block(tmp_w, M.NumBlockRows()-1, M);
//     zz_pX tmp_q;
//     tmp_q.rep.SetLength(m);
//     zz_p* coef_q = tmp_q.rep.elts();
//     for (long i = 0; i < m; i++)
//         coef_q[i] = tmp_w[i];
//     tmp_q.normalize();
//     Y_int.evaluate(tmp_w, tmp_q);
//     for (long i = 0; i < m; i++)
//         Y[i][alpha] = tmp_w[i] * f[i];

//     tmp_z = c;
//     for (long i = 0; i < m; i++){
//         Y[i][alpha+1] = -f[i]*power(tmp_z, m);
//         tmp_z = tmp_z * b;
//     }

//     CL = cauchy_like_geometric_special(X, Y, C);
// }


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
