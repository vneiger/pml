#include "mat_lzz_pX_extra.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* clears the matrix  (pmat = 0 with same dimensions)         */
/*------------------------------------------------------------*/
void clear(Mat<zz_pX> & pmat)
{
    for (long i=0; i<pmat.NumRows(); ++i)
        for (long j=0; j<pmat.NumCols(); ++j)
            clear(pmat[i][j]);
}

/*------------------------------------------------------------*/
/* set pmat to be the identity of size dim                    */
/*------------------------------------------------------------*/
void ident(Mat<zz_pX> & pmat, long dim)
{
    pmat.SetDims(dim, dim);
    for (long i = 0; i < dim; ++i)
    {
        for (long j = 0; j < i; ++j)
            clear(pmat[i][j]);
        set(pmat[i][i]);
        for (long j = i+1; j < dim; ++j)
            clear(pmat[i][j]);
    }
}

/*------------------------------------------------------------*/
/* set pmat to be the identity of size dim                    */
/*------------------------------------------------------------*/
Mat<zz_pX> ident_mat_zz_pX(long dim)
{
    Mat<zz_pX> pmat;
    pmat.SetDims(dim, dim);
    for (long i = 0; i < dim; ++i)
        set(pmat[i][i]);
    return pmat;
}

/*------------------------------------------------------------*/
/* tests whether vec is the zero vector (whatever its dim)    */
/*------------------------------------------------------------*/
long IsZero(const Vec<zz_pX> & vec)
{
    for (long i = 0; i < vec.length(); ++i)
        if (!IsZero(vec[i]))
            return 0;
    return 1;
}

/*------------------------------------------------------------*/
/* tests whether pmat is the zero matrix (whatever its dims)  */
/*------------------------------------------------------------*/
long IsZero(const Mat<zz_pX> & pmat)
{
    for (long i = 0; i < pmat.NumRows(); ++i)
        for (long j = 0; j < pmat.NumCols(); ++j)
            if (!IsZero(pmat[i][j]))
                return 0;
    return 1;
}

/*------------------------------------------------------------*/
/* tests whether pmat is the identity matrix                  */
/*------------------------------------------------------------*/
long IsIdent(const Mat<zz_pX> & pmat)
{
    if (pmat.NumRows() != pmat.NumCols())
        return 0;

    for (long i = 0; i < pmat.NumRows(); ++i)
        for (long j = 0; j < pmat.NumCols(); ++j)
        {
            if (i == j && ! IsOne(pmat[i][j]))
                return 0;
            if (i != j && ! IsZero(pmat[i][j]))
                return 0;
        }
    return 1;
}

/*------------------------------------------------------------*/
/* tests whether pmat is the identity matrix of size 'dim'    */
/*------------------------------------------------------------*/
long IsIdent(const Mat<zz_pX> & pmat, long dim)
{
    if (pmat.NumRows() != dim || pmat.NumCols() != dim)
        return 0;

    for (long i = 0; i < dim; ++i)
        for (long j = 0; j < dim; ++j)
        {
            if (i == j && !IsOne(pmat[i][j]))
                return 0;
            if (i != j && !IsZero(pmat[i][j]))
                return 0;
        }
    return 1;
}

/*------------------------------------------------------------*/
/* maximum degree of the entries of a vector / matrix         */
/*------------------------------------------------------------*/
long deg(const Vec<zz_pX> & pvec)
{
    long d = -1;
    for (long j = 0; j < pvec.length(); ++j)
        if (deg(pvec[j]) > d)
            d = deg(pvec[j]);
    return d;
}

long deg(const Mat<zz_pX> & a)
{
    long d = -1;
    for (long i = 0; i < a.NumRows(); ++i)
        for (long j = 0; j < a.NumCols(); ++j)
            if (deg(a[i][j]) > d)
                d = deg(a[i][j]);
    return d;
}

/*------------------------------------------------------------*/
/* sets x = ith coefficient of a                              */
/*------------------------------------------------------------*/
void GetCoeff(Mat<zz_p>& x, const Mat<zz_pX>& a, long i)
{
    x.SetDims(a.NumRows(), a.NumCols());

    for (long u = 0; u < a.NumRows(); ++u)
        for (long v = 0; v < a.NumCols(); ++v)
            x[u][v] = coeff(a[u][v], i);
}

/*------------------------------------------------------------*/
/* sets ith coefficient of x to a                             */
/*------------------------------------------------------------*/
void SetCoeff(Mat<zz_pX>& x, long i, const Mat<zz_p> &a)
{
    long m = x.NumRows();
    long n = x.NumCols();

    if (m != a.NumRows() || n != a.NumCols())
        LogicError("dimension mismatch in matrix SetCoeff");

    for (long u = 0; u < m; ++u)
        for (long v = 0; v < n; ++v)
            SetCoeff(x[u][v], i, a[u][v]);
}

/*------------------------------------------------------------*/
/* transpose                                                  */
/*------------------------------------------------------------*/
void transpose(Mat<zz_pX>& x, const Mat<zz_pX>& a)
{
    long n = a.NumRows();
    long m = a.NumCols();

    long i, j;

    if (&x == &a) 
    {
        if (n == m)
            for (i = 1; i <= n; i++)
                for (j = i+1; j <= n; j++)
                    swap(x(i, j), x(j, i));
        else 
        {
            Mat<zz_pX> tmp;
            tmp.SetDims(m, n);
            for (i = 1; i <= n; i++)
                for (j = 1; j <= m; j++)
                    tmp(j, i) = a(i, j);
            x.kill();
            x = tmp;
        }
    }
    else 
    {
        x.SetDims(m, n);
        for (i = 1; i <= n; i++)
            for (j = 1; j <= m; j++)
                x(j, i) = a(i, j);
    }
}

/*------------------------------------------------------------*/
/* matrix / vector truncate                                   */
/*------------------------------------------------------------*/
void trunc(Vec<zz_pX>& x, const Vec<zz_pX>& a, long n)
{
    x.SetLength(a.length());  // does nothing if x is a; aliasing is OK
    for (long i = 0; i < a.length(); i++)
        trunc(x[i], a[i], n);
}

void trunc(Mat<zz_pX>& x, const Mat<zz_pX>& a, long n)
{
    x.SetDims(a.NumRows(), a.NumCols()); // does nothing if x is a
    for (long r = 0; r < x.NumRows(); ++r)
        for (long c = 0; c < x.NumCols(); ++c)
            trunc(x[r][c], a[r][c], n);
}

void truncRow(Mat<zz_pX>& x, const Mat<zz_pX>& a, long r, long n)
{
    x.SetDims(a.NumRows(), a.NumCols()); // does nothing if x is a
    for (long c = 0; c < x.NumCols(); ++c)
        trunc(x[r][c], a[r][c], n);
}

void truncCol(Mat<zz_pX>& x, const Mat<zz_pX>& a, long c, long n)
{
    x.SetDims(a.NumRows(), a.NumCols()); // does nothing if x is a
    for (long r = 0; r < x.NumRows(); ++r)
        trunc(x[r][c], a[r][c], n);
}

/*------------------------------------------------------------*/
/* Left/Right shift of the whole matrix                       */
/*------------------------------------------------------------*/

void LeftShift(Mat<zz_pX>& x, const Mat<zz_pX>& a, long n)
{
    x.SetDims(a.NumRows(), a.NumCols());
    for (long i = 0; i < a.NumRows(); ++i)
        for (long j = 0; j < a.NumCols(); ++j)
            LeftShift(x[i][j], a[i][j], n);
}

void RightShift(Mat<zz_pX>& x, const Mat<zz_pX>& a, long n)
{
    x.SetDims(a.NumRows(), a.NumCols());
    for (long i = 0; i < a.NumRows(); ++i)
        for (long j = 0; j < a.NumCols(); ++j)
            RightShift(x[i][j], a[i][j], n);
}

/*------------------------------------------------------------*/
/* Left/Right shift of a row of the matrix                    */
/*------------------------------------------------------------*/

void LeftShiftRow(Mat<zz_pX>& x, const Mat<zz_pX>& a, const long r, long n)
{
    if (&x == &a)
        for (long j = 0; j < x.NumCols(); ++j)
            LeftShift(x[r][j], x[r][j], n);
    else
    {
        x.SetDims(a.NumRows(), a.NumCols());
        for (long i = 0; i < r; ++i)
            x[i] = a[i];
        for (long j = 0; j < x.NumCols(); ++j)
            LeftShift(x[r][j], a[r][j], n);
        for (long i = r+1; i < x.NumRows(); ++i)
            x[i] = a[i];
    }
}

void RightShiftRow(Mat<zz_pX>& x, const Mat<zz_pX>& a, const long r, long n)
{
    if (&x == &a)
        for (long j = 0; j < x.NumCols(); ++j)
            RightShift(x[r][j], x[r][j], n);
    else
    {
        x.SetDims(a.NumRows(), a.NumCols());
        for (long i = 0; i < r; ++i)
            x[i] = a[i];
        for (long j = 0; j < x.NumCols(); ++j)
            RightShift(x[r][j], a[r][j], n);
        for (long i = r+1; i < x.NumRows(); ++i)
            x[i] = a[i];
    }
}

void LeftShiftCol(Mat<zz_pX>& x, const Mat<zz_pX>& a, const long c, long n)
{
    if (&x == &a)
        for (long i = 0; i < x.NumRows(); ++i)
            LeftShift(x[i][c], x[i][c], n);
    else
    {
        x.SetDims(a.NumRows(), a.NumCols());
        for (long i = 0; i < x.NumRows(); ++i)
        {
            for (long j = 0; j < c; ++j)
                x[i][j] = a[i][j];
            LeftShift(x[i][c], a[i][c], n);
            for (long j = c+1; j < x.NumCols(); ++j)
                x[i][j] = a[i][j];
        }
    }
}

void RightShiftCol(Mat<zz_pX>& x, const Mat<zz_pX>& a, const long c, long n)
{
    if (&x == &a)
        for (long i = 0; i < x.NumRows(); ++i)
            RightShift(x[i][c], x[i][c], n);
    else
    {
        x.SetDims(a.NumRows(), a.NumCols());
        for (long i = 0; i < x.NumRows(); ++i)
        {
            for (long j = 0; j < c; ++j)
                x[i][j] = a[i][j];
            RightShift(x[i][c], a[i][c], n);
            for (long j = c+1; j < x.NumCols(); ++j)
                x[i][j] = a[i][j];
        }
    }
}

/*------------------------------------------------------------*/
/* reverse the order of the entries in a vector               */
/* x = a[n - 1 -i], i=0..n-1, with n=length(a)                */
/*------------------------------------------------------------*/
void reverse_vector(Vec<zz_pX>& x, const Vec<zz_pX>& a)
{
    long n = a.length();
    if (&x == &a)
    {
        Vec<zz_pX> tmp;
        tmp.SetLength(n);
        for (long i = 0; i < n; i++)
            tmp[i] = a[n - 1 - i];
        x = tmp;
        return;
    }
    else
    {
        x.SetLength(n);
        for (long i = 0; i < n; i++)
            x[i] = a[n - 1 - i];
    }
}

/*------------------------------------------------------------*/
/* reverse operations                                         */
/* x = reverse of a[0]..a[hi] (hi >= -1);                     */
/*------------------------------------------------------------*/
void reverse(Mat<zz_pX>& x, const Mat<zz_pX>& a, long hi)
{
    long rdim, cdim;
    rdim = a.NumRows();
    cdim = a.NumCols();
    x.SetDims(rdim, cdim);

    for (long r = 0; r < rdim; r++)
        for (long s = 0; s< cdim; s++)
            reverse(x[r][s], a[r][s], hi);
}

/*------------------------------------------------------------*/
/* row-wise reverse                                           */
/* length of 'hi' must be the number of rows of a if row-wise,*/
/* its number of columns otherwise                            */
/*------------------------------------------------------------*/

void reverse(
             Mat<zz_pX> &x, 
             const Mat<zz_pX> &a, 
             const VecLong & hi,
             const bool row_wise
            )
{
    long rdim = a.NumRows();
    long cdim = a.NumCols();
    x.SetDims(rdim, cdim);

    if (row_wise)
        for (long r = 0; r < rdim; ++r)
            for (long s = 0; s < cdim; ++s)
                reverse(x[r][s], a[r][s], hi[r]);
    else
        for (long r = 0; r < rdim; ++r)
            for (long s = 0; s < cdim; ++s)
                reverse(x[r][s], a[r][s], hi[s]);
}


/*------------------------------------------------------------*/
/* evaluate at a given point                                  */
/*------------------------------------------------------------*/
void eval(Mat<zz_p> & evmat, const Mat<zz_pX> & pmat, zz_p pt)
{
    evmat.SetDims(pmat.NumRows(),pmat.NumCols());
    for (long i = 0; i < pmat.NumRows(); ++i)
        for (long j = 0; j < pmat.NumCols(); ++j)
            eval(evmat[i][j], pmat[i][j], pt);
}


/*------------------------------------------------------------*/
/* random matrix of length n, degree < d                      */
/*------------------------------------------------------------*/
void random(Vec<zz_pX> & pvec, long n, long d)
{
    pvec.SetLength(n);
    for (long i = 0; i < n; ++i)
        random(pvec[i], d);
}

/*------------------------------------------------------------*/
/* random (m, n) matrix of degree < d                         */
/*------------------------------------------------------------*/
void random(Mat<zz_pX> & pmat, long m, long n, long d)
{
    pmat.SetDims(m, n);
    for (long i = 0; i < m; ++i)
        for (long j = 0; j < n; ++j)
            random(pmat[i][j], d);
}

/*------------------------------------------------------------*/
/* random (m, n) matrix of row degree < rdeg                  */
/*------------------------------------------------------------*/
void random_mat_zz_pX_rdeg(Mat<zz_pX>& pmat, long m, long n, VecLong rdeg)
{
    pmat.SetDims(m, n);
    for (long i = 0; i < m; ++i)
        for (long j = 0; j < n; ++j)
            random(pmat[i][j], rdeg[i]);
}

/*------------------------------------------------------------*/
/* random (m, n) matrix of column degree < cdeg               */
/*------------------------------------------------------------*/
void random_mat_zz_pX_cdeg(Mat<zz_pX>& pmat, long m, long n, VecLong cdeg)
{
    pmat.SetDims(m, n);
    for (long i = 0; i < m; ++i)
        for (long j = 0; j < n; ++j)
            random(pmat[i][j], cdeg[j]);
}

/*------------------------------------------------------------*/
/* convert from Mat<zz_p>                                     */
/*------------------------------------------------------------*/
void conv(Mat<zz_pX>& mat, const Mat<zz_p>& coeff)
{
    clear(mat);
    mat.SetDims(coeff.NumRows(), coeff.NumCols());
    SetCoeff(mat, 0, coeff);
}


/*------------------------------------------------------------*/
/* convert to / from Vec<Mat<zz_p>>                           */
/* (degree deduced from input)                                */
/*------------------------------------------------------------*/
void conv(
          Vec<Mat<zz_p>> & coeffs,
          const Mat<zz_pX> & mat
         )
{
    long d = deg(mat);
    coeffs.SetLength(d + 1);
    long r = mat.NumRows();
    long s = mat.NumCols();
    for (long i = 0; i <= d; ++i)
    {
        coeffs[i].SetDims(r, s);
        for (long a = 0; a < r; ++a)
        {
            // TODO does this actually help the compiler, or is it smart enough to
            // remain efficient if we use more readable code?
            zz_p * entries = coeffs[i][a].elts(); 
            const zz_pX * entries_mat = mat[a].elts();
            for (long b = 0; b < s; ++b)
                entries[b] = coeff(entries_mat[b], i);
        }
    }
}

void conv(
          Mat<zz_pX> & mat,
          const Vec<Mat<zz_p>> & coeffs
         )
{
    long len = coeffs.length();
    if (len == 0)
    {
        // TODO this is a choice --> indicate it in comments
        // (zero-length sequence could be zero matrix)
        mat.SetDims(0, 0);
        return;
    }
    long r = coeffs[0].NumRows();
    long s = coeffs[0].NumCols();
    mat.SetDims(r, s);

    for (long a = 0; a < r; ++a)
    {
        for (long b = 0; b < s; ++b)
        {
            zz_pX & entry = mat[a][b];
            for (long i = 0; i < len; ++i)
            {
                SetCoeff(entry, i, coeffs[i][a][b]);
            }
        }
    }
}

/*------------------------------------------------------------*/
/* convert to / from Vec<Mat<zz_p>>                           */
/* (user provided truncation order)                           */
/* coeffs will have length order independently of deg(mat)    */
/*------------------------------------------------------------*/
void conv(
          Vec<Mat<zz_p>>& coeffs,
          const Mat<zz_pX>& mat,
          const long order
         )
{
    coeffs.SetLength(order);
    long r = mat.NumRows();
    long s = mat.NumCols();
    for (long i = 0; i < order; ++i)
    {
        coeffs[i].SetDims(r, s);
        for (long a = 0; a < r; ++a)
        {
            zz_p * entries = coeffs[i][a].elts();
            const zz_pX * entries_mat = mat[a].elts();
            for (long b = 0; b < s; ++b)
                entries[b] = coeff(entries_mat[b], i);
        }
    }
}

void conv(
          Mat<zz_pX> & mat,
          const Vec<Mat<zz_p>> & coeffs,
          const long order
         )
{
    long len = std::min(order,coeffs.length());
    if (len == 0)
    {
        mat.SetDims(0, 0);
        return;
    }
    long r = coeffs[0].NumRows();
    long s = coeffs[0].NumCols();
    mat.SetDims(r, s);

    for (long a = 0; a < r; ++a)
    {
        for (long b = 0; b < s; ++b)
        {
            zz_pX & entry = mat[a][b];
            for (long i = 0; i < len; ++i)
                SetCoeff(entry, i, coeffs[i][a][b]);
        }
    }
}




// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
