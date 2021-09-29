#include "mat_lzz_pX_utils.h"
#include <algorithm> // for std::reverse

// FIXME work in progress:
// constant used to reduce cache misses in conversions
// for the moment, experimentally chosen...
#define CACHE_FRIENDLY_SIZE 32
#define MATRIX_BLOCK_SIZE 8

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
/* clears a matrix window starting at (r_offset,c_offset) and */
/* with dimensions nrows x ncols                              */
/*------------------------------------------------------------*/
void clear(Mat<zz_pX> & pmat, long r_offset, long c_offset, long nrows, long ncols)
{
    const long r_end = std::min(r_offset+nrows,pmat.NumRows());
    const long c_end = std::min(c_offset+ncols,pmat.NumCols());
    for (long i=r_offset; i<r_end; ++i)
        for (long j=c_offset; j<c_end; ++j)
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
/* build and return the identity of size dim                  */
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
bool IsZero(const Vec<zz_pX> & vec)
{
    for (long i = 0; i < vec.length(); ++i)
        if (!IsZero(vec[i]))
            return false;
    return true;
}

/*------------------------------------------------------------*/
/* tests whether pmat is the zero matrix (whatever its dims)  */
/*------------------------------------------------------------*/
bool IsZero(const Mat<zz_pX> & pmat)
{
    for (long i = 0; i < pmat.NumRows(); ++i)
        for (long j = 0; j < pmat.NumCols(); ++j)
            if (!IsZero(pmat[i][j]))
                return false;
    return true;
}

/*------------------------------------------------------------*/
/* tests whether pmat is the identity matrix                  */
/*------------------------------------------------------------*/
bool IsIdent(const Mat<zz_pX> & pmat)
{
    if (pmat.NumRows() != pmat.NumCols())
        return false;

    for (long i = 0; i < pmat.NumRows(); ++i)
        for (long j = 0; j < pmat.NumRows(); ++j)
        {
            if (i == j && ! IsOne(pmat[i][j]))
                return false;
            if (i != j && ! IsZero(pmat[i][j]))
                return false;
        }
    return true;
}

/*------------------------------------------------------------*/
/* tests whether pmat is the identity matrix of size 'dim'    */
/*------------------------------------------------------------*/
bool IsIdent(const Mat<zz_pX> & pmat, long dim)
{
    if (pmat.NumRows() != dim || pmat.NumCols() != dim)
        return false;

    for (long i = 0; i < dim; ++i)
        for (long j = 0; j < dim; ++j)
        {
            if (i == j && !IsOne(pmat[i][j]))
                return false;
            if (i != j && !IsZero(pmat[i][j]))
                return false;
        }
    return true;
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
/* tests whether pmat is a constant matrix                    */
/*------------------------------------------------------------*/
bool IsConstant(const Mat<zz_pX> & pmat)
{
    for (long i = 0; i < pmat.NumRows(); ++i)
        for (long j = 0; j < pmat.NumCols(); ++j)
            if (deg(pmat[i][j]) > 0)
                return false;
    return true;
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
    const long n = a.NumRows();
    const long m = a.NumCols();

    if (&x == &a) 
    {
        if (n == m)
            for (long i = 0; i < n; ++i)
                for (long j = i+1; j < n; ++j)
                    swap(x[i][j], x[j][i]);
        else 
        {
            Mat<zz_pX> tmp(INIT_SIZE, m, n);
            for (long i = 0; i < n; ++i)
                for (long j = 0; j < m; ++j)
                    tmp[j][i] = a[i][j];
            x.swap(tmp);
        }
    }
    else 
    {
        x.SetDims(m, n);
        for (long i = 0; i < n; ++i)
            for (long j = 0; j < m; ++j)
                x[j][i] = a[i][j];
    }
}

/*------------------------------------------------------------*/
/* mirror: reverse the order of the entries in a vector       */
/*------------------------------------------------------------*/
void mirror(Vec<zz_pX>& mvec, const Vec<zz_pX>& pvec)
{
    const long n = pvec.length();
    if (n==0)
    {
        mvec.SetLength(0);
        return;
    }

    // rely on reverse from Algorithm library of standard C++
    if (&mvec != &pvec)
        mvec = pvec;

    std::reverse(&mvec[0], &mvec[n]);
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
/* Left/Right shift of a vector                               */
/*------------------------------------------------------------*/

// left shift, vector
void LeftShift(Vec<zz_pX>& x, const Vec<zz_pX>& a, long n)
{
    x.SetLength(a.length());
    for (long i = 0; i < a.length(); ++i)
        LeftShift(x[i], a[i], n);
}

// right shift, vector
void RightShift(Vec<zz_pX>& x, const Vec<zz_pX>& a, long n)
{
    x.SetLength(a.length());
    for (long i = 0; i < a.length(); ++i)
        RightShift(x[i], a[i], n);
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
/* reverse operations                                         */
/* x = reverse of a[0]..a[hi]                                 */
/*------------------------------------------------------------*/
void reverse(Mat<zz_pX>& x, const Mat<zz_pX>& a, long hi)
{
    long rdim = a.NumRows();
    long cdim = a.NumCols();
    x.SetDims(rdim, cdim);

    for (long r = 0; r < rdim; ++r)
        for (long s = 0; s < cdim; ++s)
            reverse(x[r][s], a[r][s], hi);
}

/*------------------------------------------------------------*/
/* column-wise reverse                                        */
/* length of 'hi' must be the number of columns of a          */
/*------------------------------------------------------------*/
void col_reverse(
             Mat<zz_pX> &x, 
             const Mat<zz_pX> &a, 
             const VecLong & hi
            )
{
    long rdim = a.NumRows();
    long cdim = a.NumCols();
    x.SetDims(rdim, cdim);

    for (long r = 0; r < rdim; ++r)
        for (long s = 0; s < cdim; ++s)
            reverse(x[r][s], a[r][s], hi[s]);
}

/*------------------------------------------------------------*/
/* row-wise reverse                                           */
/* length of 'hi' must be the number of rows of a             */
/*------------------------------------------------------------*/
void row_reverse(
                 Mat<zz_pX> &x, 
                 const Mat<zz_pX> &a, 
                 const VecLong & hi
                )
{
    long rdim = a.NumRows();
    long cdim = a.NumCols();
    x.SetDims(rdim, cdim);

    for (long r = 0; r < rdim; ++r)
        for (long s = 0; s < cdim; ++s)
            reverse(x[r][s], a[r][s], hi[r]);
}


/*------------------------------------------------------------*/
/* evaluate at a given point                                  */
/*------------------------------------------------------------*/
void eval(Mat<zz_p> & evmat, const Mat<zz_pX> & pmat, const zz_p & pt)
{
    evmat.SetDims(pmat.NumRows(),pmat.NumCols());
    for (long i = 0; i < pmat.NumRows(); ++i)
        for (long j = 0; j < pmat.NumCols(); ++j)
            eval(evmat[i][j], pmat[i][j], pt);
}

/*------------------------------------------------------------*/
/* evaluate at a given point (rep = vector of matrices)       */
/*------------------------------------------------------------*/
void eval(Mat<zz_p> & evmat, const Vec<Mat<zz_p>> & matp, const zz_p & pt)
{
    long d = matp.length()-1; // degree
    if (d==-1) // zero matrix, dimensions unknown
    {
        evmat.SetDims(0,0);
        return;
    }
    evmat = matp[d];
    for (--d; d >= 0; --d)
    {
        mul(evmat, evmat, pt);
        add(evmat, evmat, matp[d]);
    }
}

/*------------------------------------------------------------*/
/* random vector of length n, degree < d                      */
/*------------------------------------------------------------*/
void random(Vec<zz_pX> & pvec, long n, long d)
{
    pvec.SetLength(n);
    for (long i = 0; i < n; ++i)
        random(pvec[i], d);
}

/*------------------------------------------------------------*/
/* random (m, n) matrix of degree < d, pmat version           */
/*------------------------------------------------------------*/
void random(Mat<zz_pX> & pmat, long m, long n, long d)
{
    pmat.SetDims(m, n);
    for (long i = 0; i < m; ++i)
        for (long j = 0; j < n; ++j)
            random(pmat[i][j], d);
}

/*------------------------------------------------------------*/
/* random (m, n) matrix of degree < d, matp version           */
/*------------------------------------------------------------*/
void random(Vec<Mat<zz_p>> & matp, long m, long n, long d)
{
    matp.SetLength(d);
    for (long k = 0; k < d; ++k)
        random(matp[k], m, n);
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
/* random matrix with degree matrix strictly less than dmat   */
/*------------------------------------------------------------*/
void random(Mat<zz_pX> & pmat, Mat<long> dmat)
{
    pmat.SetDims(dmat.NumRows(), dmat.NumCols());
    for (long i = 0; i < dmat.NumRows(); ++i)
        for (long j = 0; j < dmat.NumCols(); ++j)
            random(pmat[i][j], dmat[i][j]);
}

/*------------------------------------------------------------*/
/* convert from Mat<zz_p>                                     */
/*------------------------------------------------------------*/
void conv(Mat<zz_pX>& pmat, const Mat<zz_p>& mat)
{
    pmat.SetDims(mat.NumRows(), mat.NumCols());
    for (long i = 0; i < mat.NumRows(); ++i)
        for (long j = 0; j < mat.NumCols(); ++j)
            conv(pmat[i][j], mat[i][j]);
}


/*------------------------------------------------------------*/
/* convert to / from Vec<Mat<zz_p>>                           */
/* (degree deduced from input)                                */
/*------------------------------------------------------------*/

// Note: may be improved when degrees in pmat are unbalanced (e.g. if only
// few entries reach degree d)
void conv(Vec<Mat<zz_p>> & matp, const Mat<zz_pX> & pmat)
{
    const long m = pmat.NumRows();
    const long n = pmat.NumCols();
    const long len = deg(pmat)+1;
    matp.SetLength(len);
    // if len==0, matp is the length-0 vector and the following loop does
    // nothing

    for (long k = 0; k < len; k+=CACHE_FRIENDLY_SIZE)
    {
        const long k_bound = std::min((long)CACHE_FRIENDLY_SIZE, len-k);
        for (long kk=0; kk<k_bound; ++kk)
            matp[k+kk].SetDims(m, n);

        for (long i = 0; i < m; ++i)
            for (long j = 0; j < n; ++j)
                for (long kk = 0; kk < k_bound; ++kk)
                    matp[k+kk][i][j] = pmat[i][j][k+kk];
    }
}

void conv(Mat<zz_pX> & pmat, const Vec<Mat<zz_p>> & matp)
{
    const long len = matp.length();
    if (len == 0)
    {
        // matp is zero; convention: do not change dimension of pmat
        clear(pmat);
        return;
    }
    if (len == 1)
    {
        // matp is constant, rely on the dedicated conversion
        conv(pmat, matp[0]);
        return;
    }

    // now, len >= 2
    const long m = matp[0].NumRows();
    const long n = matp[0].NumCols();
    pmat.SetDims(m, n);

    for (long i = 0; i < m; ++i)
    {
        for (long j = 0; j < n; j+=CACHE_FRIENDLY_SIZE)
        {
            const long j_bound = std::min((long)CACHE_FRIENDLY_SIZE, n-j);
            // initialize vectors
            for (long jj = 0; jj < j_bound; ++jj)
                pmat[i][j+jj].SetLength(len);

            // fill data
            for (long k = 0; k < len; k+=MATRIX_BLOCK_SIZE)
            {
                const long k_bound = std::min((long)MATRIX_BLOCK_SIZE, len-k);
                for (long jj = 0; jj < j_bound; ++jj)
                    for (long kk = 0; kk < k_bound; ++kk)
                        pmat[i][j+jj][k+kk] = matp[k+kk][i][j+jj];
            }

            // strip away leading zeroes
            for (long jj = 0; jj < j_bound; ++jj)
                pmat[i][j+jj].normalize();
        }
    }
}


/*------------------------------------------------------------*/
/* convert to / from Vec<Mat<zz_p>>                           */
/* (user provided truncation order)                           */
/* matp will have length order independently of deg(mat)      */
/*------------------------------------------------------------*/
void conv(
          Vec<Mat<zz_p>>& matp,
          const Mat<zz_pX>& pmat,
          const long order
         )
{
    const long m = pmat.NumRows();
    const long n = pmat.NumCols();
    matp.SetLength(order);
    for (long k = 0; k < order; ++k)
    {
        matp[k].SetDims(m, n);
        for (long i = 0; i < m; ++i)
            for (long j = 0; j < n; ++j)
                matp[k][i][j] = coeff(pmat[i][j], k);
    }
}

void conv(
          Mat<zz_pX> & pmat,
          const Vec<Mat<zz_p>> & matp,
          const long order
         )
{
    const long len = std::min(order,matp.length());
    if (len == 0)
    {
        clear(pmat); // keeping the same dimensions
        return;
    }

    const long m = matp[0].NumRows();
    const long n = matp[0].NumCols();
    pmat.SetDims(m, n);
    for (long i = 0; i < m; ++i)
        for (long j = 0; j < n; ++j)
        {
            pmat[i][j].SetLength(len);
            for (long k = 0; k < len; ++k)
                pmat[i][j][k] = matp[k][i][j];
            pmat[i][j].normalize();
        }
}


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
