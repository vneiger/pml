#include "lzz_pX_extra.h"
#include "mat_lzz_pX_utils.h"
#include "mat_lzz_pX_arith.h"

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
    const long m = a.NumRows();
    const long n = a.NumCols();

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
    const long m = a.NumRows();
    const long n = a.NumCols();

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
    const long m = a.length();

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
    const long m = a.length();

    if (m != b.length())
        LogicError("dimension mismatch in vector addition");

    c.SetLength(m);
    for (long i = 0; i < m; ++i)
        add(c[i], a[i], b[i]);
}

/*------------------------------------------------------------*/
/* Vec<zz_pX> addition of left shift                          */
/* c can alias a but not b                                    */
/*------------------------------------------------------------*/
void add_LeftShift(Vec<zz_pX> & c, const Vec<zz_pX> & a, const Vec<zz_pX> & b, long k)
{
    const long m = a.length();
    if (m != b.length())
        LogicError("dimension mismatch in vector leftshift-addition");

    c.SetLength(m);
    for (long i = 0; i < m; ++i)
        add_LeftShift(c[i], a[i], b[i], k);
}

/*------------------------------------------------------------*/
/* Mat<zz_pX> addition of left shift                          */
/* c can alias a but not b                                    */
/*------------------------------------------------------------*/
void add_LeftShift(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b, long k)
{
    const long m = a.NumRows();
    const long n = a.NumCols();
    if (m != b.NumRows() || n != b.NumCols())
        LogicError("dimension mismatch in matrix leftshift-addition");

    c.SetDims(m, n);
    for (long i = 0; i < m; ++i)
        for (long j = 0; j < n; ++j)
            add_LeftShift(c[i][j], a[i][j], b[i][j], k);
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
/* negate                                                     */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

NTL_OPEN_NNS
void negate(Vec<zz_pX> & x, const Vec<zz_pX> & a)
{
    long n = a.length();
    x.SetLength(n);
    for (long i = 0; i < n; ++i)
            NTL::negate(x[i], a[i]);
}
NTL_CLOSE_NNS


NTL_OPEN_NNS
void negate(Mat<zz_pX> & x, const Mat<zz_pX> & a)
{
    long m = a.NumRows();
    long n = a.NumCols();

    x.SetDims(m, n);
    for (long u = 0; u < m; ++u)
        for (long v = 0; v < n; ++v)
            NTL::negate(x[u][v], a[u][v]);
}
NTL_CLOSE_NNS



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
    if (&c == &a)
    {
        Mat<zz_pX> c2;
        mul(c2, a, b);
        c.swap(c2);
        return;
    }

    const long d = deg(a);

    const long m = a.NumRows();
    const long n = a.NumCols();
    const long p = b.NumCols();

    if (n != b.NumRows())
        LogicError("dimension mismatch in multiplication by constant matrix");

    c.SetDims(m, p);
    clear(c);
    Mat<zz_p> tmp(INIT_SIZE, m, n);
    Mat<zz_p> res(INIT_SIZE, m, p);
    for (long i = 0; i <= d; ++i)
    {
        for (long u = 0; u < m; ++u)
            for (long v = 0; v < n; ++v)
                tmp[u][v] = coeff(a[u][v], i);
        mul(res, tmp, b);
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
    if (&c == &b)
    {
        Mat<zz_pX> c2;
        mul(c2, a, b);
        c.swap(c2);
        return;
    }

    const long d = deg(b);

    const long m = a.NumRows();
    const long n = a.NumCols();
    const long p = b.NumCols();

    if (n != b.NumRows())
        LogicError("dimension mismatch in multiplication by constant matrix");

    c.SetDims(m, p);
    clear(c);

    Mat<zz_p> tmp(INIT_SIZE, n, p);
    Mat<zz_p> res(INIT_SIZE, m, p);
    for (long i = 0; i <= d; ++i)
    {
        for (long u = 0; u < n; ++u)
            for (long v = 0; v < p; ++v)
                tmp[u][v] = coeff(b[u][v], i);
        mul(res, a, tmp);
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

void mul(Vec<zz_pX> & c, const Vec<zz_pX> & a, const zz_p & b)
{
    long n = a.length();
    c.SetLength(n);
    for (long i = 0; i < n; ++i)
        mul(c[i], a[i], b);
}

void mul(Mat<zz_pX> & c, const Mat<zz_pX> & a, const zz_p & b)
{
    long m = a.NumRows();
    long n = a.NumCols();

    c.SetDims(m, n);
    for (long u = 0; u < m; ++u)
        for (long v = 0; v < n; ++v)
            mul(c[u][v], a[u][v], b);
}

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* polynomial multiplication                                  */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

void mul(Vec<zz_pX> & c, const Vec<zz_pX> & a, const zz_pX & b)
{
    long n = a.length();
    c.SetLength(n);
    for (long i = 0; i < n; ++i)
        mul(c[i], a[i], b);
}

void mul(Mat<zz_pX> & c, const Mat<zz_pX> & a, const zz_pX & b)
{
    long m = a.NumRows();
    long n = a.NumCols();

    c.SetDims(m, n);
    for (long u = 0; u < m; ++u)
        for (long v = 0; v < n; ++v)
            mul(c[u][v], a[u][v], b);
}



// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
