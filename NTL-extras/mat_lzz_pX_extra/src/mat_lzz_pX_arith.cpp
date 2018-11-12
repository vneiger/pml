#include "mat_lzz_pX_arith.h"
#include "mat_lzz_pX_utils.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* basic arithmetic: addition, subtraction, negation,         */
/* multiplication by constant / by single polynomial, ...     */
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* addition                                                   */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* Mat<zz_pX> addition                                        */
/* c can alias a or b, does not have to be zero               */
/*------------------------------------------------------------*/
void add(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b)
{
    long m = a.NumRows();
    long n = a.NumCols();

    if (m != b.NumRows() || n != b.NumCols())
        LogicError("dimension mismatch in matrix addition");

    c.SetDims(m, n);
    for (long i = 0; i < m; ++i)
        for (long j = 0; j < n; ++j)
            add(c[i][j], a[i][j], b[i][j]);
}

/*------------------------------------------------------------*/
/* Mat<zz_pX> addition, rhs is constant                       */
/*------------------------------------------------------------*/
void add(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_p> & b)
{
    long m = a.NumRows();
    long n = a.NumCols();

    if (m != b.NumRows() || n != b.NumCols())
        LogicError("dimension mismatch in matrix addition");

    c.SetDims(m, n);
    for (long i = 0; i < m; ++i)
        for (long j = 0; j < n; ++j)
            add(c[i][j], a[i][j], b[i][j]);
}

/*------------------------------------------------------------*/
/* Vec<zz_pX> addition                                        */
/* c can alias a or b, does not have to be zero               */
/*------------------------------------------------------------*/
void add(Vec<zz_pX> & c, const Vec<zz_pX> & a, const Vec<zz_pX> & b)
{
    long m = a.length();

    if (m != b.length())
        LogicError("dimension mismatch in vector addition");

    c.SetLength(m);
    for (long i = 0; i < m; ++i)
        add(c[i], a[i], b[i]);
}

/*------------------------------------------------------------*/
/* Vec<zz_pX> addition, rhs is constant                       */
/*------------------------------------------------------------*/
void add(Vec<zz_pX> & c, const Vec<zz_pX> & a, const Vec<zz_p> & b)
{
    long m = a.length();

    if (m != b.length())
        LogicError("dimension mismatch in vector addition");

    c.SetLength(m);
    for (long i = 0; i < m; ++i)
        add(c[i], a[i], b[i]);
}


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* subtraction                                                */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* Mat<zz_pX> subtraction                                     */
/*------------------------------------------------------------*/
void sub(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b)
{
    long m = a.NumRows();
    long n = a.NumCols();

    if (m != b.NumRows() || n != b.NumCols())
        LogicError("dimension mismatch in matrix subtraction");

    c.SetDims(m, n);
    for (long i = 0; i < m; ++i)
        for (long j = 0; j < n; ++j)
            sub(c[i][j], a[i][j], b[i][j]);
}

/*------------------------------------------------------------*/
/* subtraction, rhs is constant                               */
/*------------------------------------------------------------*/
void sub(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_p> & b)
{
    long m = a.NumRows();
    long n = a.NumCols();

    if (m != b.NumRows() || n != b.NumCols())
        LogicError("dimension mismatch in matrix subtraction");

    c.SetDims(m, n);
    for (long i = 0; i < m; ++i)
        for (long j = 0; j < n; ++j)
            sub(c[i][j], a[i][j], b[i][j]);
}

/*------------------------------------------------------------*/
/* subtraction, lhs is constant                               */
/*------------------------------------------------------------*/
void sub(Mat<zz_pX> & c, const Mat<zz_p> & a, const Mat<zz_pX> & b)
{
    long m = a.NumRows();
    long n = a.NumCols();

    if (m != b.NumRows() || n != b.NumCols())
        LogicError("dimension mismatch in matrix subtraction");

    c.SetDims(m, n);
    for (long i = 0; i < m; ++i)
        for (long j = 0; j < n; ++j)
            sub(c[i][j], a[i][j], b[i][j]);
}


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* vector subtraction                                         */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* zz_pX subtraction                                          */
/*------------------------------------------------------------*/
void sub(Vec<zz_pX> & c, const Vec<zz_pX> & a, const Vec<zz_pX> & b)
{
    long m = a.length();

    if (m != b.length())
        LogicError("dimension mismatch in vector addition");

    c.SetLength(m);
    for (long i = 0; i < m; ++i)
        sub(c[i], a[i], b[i]);
}

/*------------------------------------------------------------*/
/* subtraction, rhs is constant                               */
/*------------------------------------------------------------*/
void sub(Vec<zz_pX> & c, const Vec<zz_pX> & a, const Vec<zz_p> & b)
{
    long m = a.length();

    if (m != b.length())
        LogicError("dimension mismatch in vector addition");

    c.SetLength(m);
    for (long i = 0; i < m; ++i)
        sub(c[i], a[i], b[i]);
}

/*------------------------------------------------------------*/
/* subtraction, lhs is constant                               */
/*------------------------------------------------------------*/
void sub(Vec<zz_pX> & c, const Vec<zz_p> & a, const Vec<zz_pX> & b)
{
    long m = a.length();

    if (m != b.length())
        LogicError("dimension mismatch in vector addition");

    c.SetLength(m);
    for (long i = 0; i < m; ++i)
        sub(c[i], a[i], b[i]);
}


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* multiplication by constant matrix                          */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* multiplication, rhs is a constant matrix                   */
/*------------------------------------------------------------*/
void mul(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_p> & b)
{
    long d = deg(a);

    long m = a.NumRows();
    long n = a.NumCols();
    long p = b.NumCols();

    if (n != b.NumRows())
        LogicError("dimension mismatch in multiplication by constant matrix");

    c.SetDims(m, p);
    Mat<zz_p> tmp, res;
    tmp.SetDims(m, n);
    for (long i = 0; i <= d; ++i)
    {
        for (long u = 0; u < m; ++u)
            for (long v = 0; v < n; ++v)
                tmp[u][v] = coeff(a[u][v], i);
        res = tmp * b;
        for (long u = 0; u < m; ++u)
            for (long v = 0; v < p; ++v)
                SetCoeff(c[u][v], i, res[u][v]);
    }
}

/*------------------------------------------------------------*/
/* multiplication, lhs is a constant matrix                   */
/*------------------------------------------------------------*/
void mul(Mat<zz_pX> & c, const Mat<zz_p> & a, const Mat<zz_pX> & b)
{
    long d = deg(b);

    long m = a.NumRows();
    long n = a.NumCols();
    long p = b.NumCols();

    if (n != b.NumRows())
        LogicError("dimension mismatch in multiplication by constant matrix");

    c.SetDims(m, p);
    Mat<zz_p> tmp, res;
    tmp.SetDims(n, p);
    for (long i = 0; i <= d; ++i)
    {
        for (long u = 0; u < n; ++u)
            for (long v = 0; v < p; ++v)
                tmp[u][v] = coeff(b[u][v], i);
        res = a * tmp;
        for (long u = 0; u < m; ++u)
            for (long v = 0; v < p; ++v)
                SetCoeff(c[u][v], i, res[u][v]);
    }
}


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* scalar multiplication                                      */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* scalar multiplication for a vector                         */
/*------------------------------------------------------------*/

void mul(Vec<zz_pX> & c, const Vec<zz_pX> & a, const zz_p & b)
{
    long n = a.length();
    c.SetLength(n);
    for (long i = 0; i < n; ++i)
        c[i] = a[i] * b;
}

/*------------------------------------------------------------*/
/* scalar multiplication for a matrix                         */
/*------------------------------------------------------------*/

void mul(Mat<zz_pX> & c, const Mat<zz_pX> & a, const zz_p & b)
{
    long m = a.NumRows();
    long n = a.NumCols();

    c.SetDims(m, n);
    for (long u = 0; u < m; ++u)
        for (long v = 0; v < n; ++v)
            c[u][v] = b * a[u][v];
}



// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
