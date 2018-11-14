#include <NTL/lzz_pX.h>
#include <NTL/mat_lzz_p.h>

#include "structured_lzz_pX.h"

/*------------------------------------------------------------*/
/* right multiplication, matrix version                       */
/*------------------------------------------------------------*/
void structured_lzz_pX::mul_right(Mat<zz_pX>& out, const Mat<zz_pX>& in) const
{
    if (&out == &in)
    {
        out = mul_right(in);
        return;
    }

    long m = NumRows();
    long n = NumCols();

    if (in.NumRows() != n)
        Error("Bad dimensions in structured_lzz_pX matrix mul right");

    long p = in.NumCols();
    Vec<zz_pX> vec_in, vec_out;
    vec_in.SetLength(n);
    out.SetDims(m, p);

    for (long j = 0; j < p; j++)
    {
        for (long i = 0; i < n; i++)
            vec_in[i] = in[i][j];
        mul_right(vec_out, vec_in);
        for (long i = 0; i < m; i++)
            out[i][j] = vec_out[i];
    }
}

/*------------------------------------------------------------*/
/* left multiplication, matrix version                        */
/*------------------------------------------------------------*/
void structured_lzz_pX::mul_left(Mat<zz_pX>& out, const Mat<zz_pX>& in) const
{
    if (&out == &in)
    {
        out = mul_left(in);
        return;
    }

    long m = NumRows();
    long n = NumCols();

    if (in.NumCols() != m)
        Error("Bad dimensions in structured_lzz_pX matrix mul left");

    Vec<zz_pX> vec_in, vec_out;
    vec_in.SetLength(m);
    long p = in.NumRows();
    out.SetDims(p, n);

    for (long j = 0; j < p; j++)
    {
        for (long i = 0; i < m; i++)
            vec_in[i] = in[j][i];
        mul_left(vec_out, vec_in);
        for (long i = 0; i < n; i++)
            out[j][i] = vec_out[i];
    }
}

/*------------------------------------------------------------*/
/* right multiplication mod x^s, matrix version               */
/*------------------------------------------------------------*/
void structured_lzz_pX::mul_right_trunc(Mat<zz_pX>& out, const Mat<zz_pX>& in, long s) const
{
    if (&out == &in)
    {
        out = mul_right_trunc(in, s);
        return;
    }

    long m = NumRows();
    long n = NumCols();

    if (in.NumRows() != n)
        Error("Bad dimensions in structured_lzz_pX matrix mul right");

    long p = in.NumCols();
    Vec<zz_pX> vec_in, vec_out;
    vec_in.SetLength(n);
    out.SetDims(m, p);

    for (long j = 0; j < p; j++)
    {
        for (long i = 0; i < n; i++)
            vec_in[i] = in[i][j];
        mul_right_trunc(vec_out, vec_in, s);
        for (long i = 0; i < m; i++)
            out[i][j] = vec_out[i];
    }
}

/*------------------------------------------------------------*/
/* left multiplication mod x^s, matrix version                */
/*------------------------------------------------------------*/
void structured_lzz_pX::mul_left_trunc(Mat<zz_pX>& out, const Mat<zz_pX>& in, long s) const
{
    if (&out == &in)
    {
        out = mul_left_trunc(in, s);
        return;
    }

    long m = NumRows();
    long n = NumCols();

    if (in.NumCols() != m)
        Error("Bad dimensions in structured_lzz_pX matrix mul left");

    Vec<zz_pX> vec_in, vec_out;
    vec_in.SetLength(m);
    long p = in.NumRows();
    out.SetDims(p, n);

    for (long j = 0; j < p; j++)
    {
        for (long i = 0; i < m; i++)
            vec_in[i] = in[j][i];
        mul_left_trunc(vec_out, vec_in, s);
        for (long i = 0; i < n; i++)
            out[j][i] = vec_out[i];
    }
}


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
