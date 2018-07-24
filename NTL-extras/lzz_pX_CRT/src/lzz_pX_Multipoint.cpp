#include <memory>
#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>

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
/* evaluates all entries of a vector                          */
/*------------------------------------------------------------*/
void zz_pX_Multipoint::evaluate_vector(Vec<Vec<zz_p>>& val, const Vec<zz_pX>& f) const
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

    Vec<zz_p> tmp;
    for (long i = 0; i < s; i++)
    {
	for (long k = 0; k < t; k++)
	{
	    evaluate(tmp, f[i][k]);
	    for (long j = 0; j < n; j++)
	    {
		val[j][i][k] = tmp[j];
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
