#include <NTL/lzz_pX.h>
#include <NTL/vector.h>

#include "lzz_pX_CRT.h"

NTL_CLIENT


/*------------------------------------------------------------*/
/* size of the input                                          */
/*------------------------------------------------------------*/
long zz_pX_Transform::input_length() const
{
    return m;
}

/*------------------------------------------------------------*/
/* size of the transform                                      */
/*------------------------------------------------------------*/
long zz_pX_Transform::transform_length() const
{
    return n;
}

/*------------------------------------------------------------*/
/* evaluates all entries of a vector                          */
/* aliases are not permitted between val[i] and f[j].rep      */
/*------------------------------------------------------------*/
void zz_pX_Transform::forward_left_vector(Vec<Vec<zz_p>>& val, const Vec<zz_pX>& f) const
{
    long s = f.length();
    val.SetLength(n);
    for (long i = 0; i < n; i++)
    {
        val[i].SetLength(s);
    }

    Vec<zz_p> tmp;
    for (long i = 0; i < s; i++)
    {
        forward_left(tmp, f[i]);
        for (long j = 0; j < n; j++)
        {
            val[j][i] = tmp[j];
        }
    }
}


/*------------------------------------------------------------*/
/* evaluates all entries of a vector                          */
/* aliases are not permitted between val[i] and f[j].rep      */
/*------------------------------------------------------------*/
void zz_pX_Transform::forward_right_vector(Vec<Vec<zz_p>>& val, const Vec<zz_pX>& f) const
{
    long s = f.length();
    val.SetLength(n);
    for (long i = 0; i < n; i++)
    {
        val[i].SetLength(s);
    }

    Vec<zz_p> tmp;
    for (long i = 0; i < s; i++)
    {
        forward_right(tmp, f[i]);
        for (long j = 0; j < n; j++)
        {
            val[j][i] = tmp[j];
        }
    }
}

/*------------------------------------------------------------*/
/* evaluates all entries of a matrix                          */
/*------------------------------------------------------------*/
void zz_pX_Transform::forward_left_matrix(Vec<Mat<zz_p>>& val, const Mat<zz_pX>& f) const
{
    long s = f.NumRows();
    long t = f.NumCols();
    val.SetLength(n);
    for (long i = 0; i < n; i++)
    {
        val[i].SetDims(s, t);
    }

    Vec<zz_p> tmp;
    for (long i = 0; i < s; i++)
    {
        for (long k = 0; k < t; k++)
        {
            forward_left(tmp, f[i][k]);
            for (long j = 0; j < n; j++)
            {
                val[j][i][k] = tmp[j];
            }
        }
    }
}

/*------------------------------------------------------------*/
/* evaluates all entries of a matrix                          */
/*------------------------------------------------------------*/
void zz_pX_Transform::forward_right_matrix(Vec<Mat<zz_p>>& val, const Mat<zz_pX>& f) const
{
    long s = f.NumRows();
    long t = f.NumCols();
    val.SetLength(n);
    for (long i = 0; i < n; i++)
    {
        val[i].SetDims(s, t);
    }

    Vec<zz_p> tmp;
    for (long i = 0; i < s; i++)
    {
        for (long k = 0; k < t; k++)
        {
            forward_right(tmp, f[i][k]);
            for (long j = 0; j < n; j++)
            {
                val[j][i][k] = tmp[j];
            }
        }
    }
}

/*------------------------------------------------------------*/
/* interpolates all entries of a vector                       */
/* aliases are not permitted between val[i] and f[j].rep      */
/*------------------------------------------------------------*/
void zz_pX_Transform::backward_vector(Vec<zz_pX>& f, const Vec<Vec<zz_p>>& val) const
{
    if (val.length() == 0)
    {
        f.SetLength(0);
        return;
    }

    long s = val[0].length();
    f.SetLength(s);

    Vec<zz_p> tmp;
    tmp.SetLength(n);
    for (long i = 0; i < s; i++)
    {
        for (long j = 0; j < n; j++)
        {
            tmp[j] = val[j][i];
        }
        backward(f[i], tmp);
    }
}

/*------------------------------------------------------------*/
/* interpolates all entries of a matrix                       */
/*------------------------------------------------------------*/
void zz_pX_Transform::backward_matrix(Mat<zz_pX>& f, const Vec<Mat<zz_p>>& val) const
{
    if (val.length() == 0)
    {
        f.SetDims(0, 0);
        return;
    }

    long s = val[0].NumRows();
    long t = val[0].NumCols();
    f.SetDims(s, t);

    Vec<zz_p> tmp;
    tmp.SetLength(n);
    for (long i = 0; i < s; i++)
    {
        for (long k = 0; k < t; k++)
        {
            for (long j = 0; j < n; j++)
            {
                tmp[j] = val[j][i][k];
            }
            backward(f[i][k], tmp);
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
