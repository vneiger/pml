#include <memory>
#include "lzz_p_extra.h"
#include "lzz_pX_CRT.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* number of points                                           */
/*------------------------------------------------------------*/
long zz_pX_Multipoint::length() const
{
    return n;
}

/*------------------------------------------------------------*/
/* get the i-th point                                         */
/*------------------------------------------------------------*/
void zz_pX_Multipoint::get_point(zz_p& pt, long i) 
{
    if (pts.length() == 0) // initializes if needed
    {
        pts.SetLength(n);
        zz_pX X;
        SetCoeff(X, 1, to_zz_p(1));
        evaluate(pts, X);
    }

    pt = pts[i];
}

/*------------------------------------------------------------*/
/* get a copy of the vector of points                         */
/*------------------------------------------------------------*/
void zz_pX_Multipoint::get_points(Vec<zz_p> & points) 
{
    if (pts.length() == 0) // initializes if needed
    {
        pts.SetLength(n);
        zz_pX X;
        SetCoeff(X, 1, to_zz_p(1));
        evaluate(pts, X);
    }

    points = pts;
}

/*------------------------------------------------------------*/
/* evaluates all entries of a vector                          */
/*------------------------------------------------------------*/
void zz_pX_Multipoint::evaluate_vector(Vec<Vec<zz_p>>& val, const Vec<zz_pX>& f) const
{
    const long s = f.length();
    val.SetLength(n);
    for (long i = 0; i < n; ++i)
        val[i].SetLength(s);

    Vec<zz_p> tmp;
    for (long i = 0; i < s; i++)
    {
        evaluate(tmp, f[i]);
        for (long j = 0; j < n; j++)
        {
            val[j][i] = tmp[j];
        }
    }
}

/*------------------------------------------------------------*/
/* evaluates all entries of a matrix                          */
/*------------------------------------------------------------*/
void zz_pX_Multipoint::evaluate_matrix(Vec<Mat<zz_p>>& val, const Mat<zz_pX>& f) const
{
    long s = f.NumRows();
    long t = f.NumCols();
    val.SetLength(n);
    for (long i = 0; i < n; i++)
    {
        val[i].SetDims(s, t);
    }

// new stuff here
    Mat<Vec<zz_p>> mat_val;
    mat_val.SetDims(s, t);

    for (long i = 0; i < s; i++)
    {
        for (long k = 0; k < t; k++)
        {
            evaluate(mat_val[i][k], f[i][k]);
        }
    }
    
    for (long j = 0; j < n; j++)
    {
        for (long i = 0; i < s; i++)
        {
            for (long k = 0; k < t; k++)
            {
                val[j][i][k] = mat_val[i][k][j];
            }
        }
    }
}

/*------------------------------------------------------------*/
/* interpolates all entries of a vector                       */
/*------------------------------------------------------------*/
void zz_pX_Multipoint::interpolate_vector(Vec<zz_pX>& f, const Vec<Vec<zz_p>>& val) 
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
        interpolate(f[i], tmp);
    }
}

/*------------------------------------------------------------*/
/* interpolates all entries of a matrix                       */
/*------------------------------------------------------------*/
void zz_pX_Multipoint::interpolate_matrix(Mat<zz_pX>& f, const Vec<Mat<zz_p>>& val)
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
            interpolate(f[i][k], tmp);
        }
    }
}

/*------------------------------------------------------------*/
/* a naive conversion to a dense matrix                       */
/*------------------------------------------------------------*/
void zz_pX_Multipoint::to_dense(Mat<zz_p>& M)
{
    zz_pX f;
    set(f);
    Vec<zz_p> v;
    v.SetLength(n);
    M.SetDims(n, n);
    
    for (long i = 0; i < n; i++)
    {
        evaluate(v, f);
        for (long j = 0; j < n; j++)
        {
            M[j][i] = v[j];
        }
        f <<= 1;
    }
}


/*------------------------------------------------------------*/
/* transpose-evaluates all entries of a vector                */
/*------------------------------------------------------------*/
void zz_pX_Multipoint::t_evaluate_vector(Vec<zz_pX>& f, const Vec<Vec<zz_p>>& val, long output_size) const
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
        t_evaluate(f[i], tmp, output_size);
    }
}

/*------------------------------------------------------------*/
/* transpose-evaluates all entries of a matrix                */
/*------------------------------------------------------------*/
void zz_pX_Multipoint::t_evaluate_matrix(Mat<zz_pX>& f, const Vec<Mat<zz_p>>& val, long output_size) const
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
            t_evaluate(f[i][k], tmp, output_size);
        }
    }
}

/*------------------------------------------------------------*/
/* transpose-interpolates all entries of a vector             */
/*------------------------------------------------------------*/
void zz_pX_Multipoint::t_interpolate_vector(Vec<Vec<zz_p>>& val, const Vec<zz_pX>& f) 
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
        t_interpolate(tmp, f[i]);
        for (long j = 0; j < n; j++)
        {
            val[j][i] = tmp[j];
        }
    }
}

/*------------------------------------------------------------*/
/* transpose-interpolates all entries of a matrix             */
/*------------------------------------------------------------*/
void zz_pX_Multipoint::t_interpolate_matrix(Vec<Mat<zz_p>>& val, const Mat<zz_pX>& f)
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
            t_interpolate(tmp, f[i][k]);
            for (long j = 0; j < n; j++)
            {
                val[j][i][k] = tmp[j];
            }
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
