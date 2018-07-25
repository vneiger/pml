#include <NTL/lzz_pX.h>
#include <NTL/vector.h>
#include <NTL/matrix.h>

#include "lzz_pX_CRT.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* constructor                                                */
/*------------------------------------------------------------*/
zz_pX_Transform_naive::zz_pX_Transform_naive(long d)
{
    m = d;
    n = d * d;
}

/*------------------------------------------------------------*/
/* (a0,a1,a2) -> (a0, a0, a0, a1, a1, a1, a2, a2, a2)         */
/*------------------------------------------------------------*/
void zz_pX_Transform_naive::forward_left(Vec<zz_p>& val, const zz_pX& P) const
{
    long idx = 0;
    val.SetLength(n);
    for (long i = 0; i < m; i++)
    {
	zz_p a = coeff(P, i);
	for (long j = 0; j < m; j++)
	{
	    val[idx++] = a;
	}
    }
}

/*------------------------------------------------------------*/
/* (b0,b1,b2) -> (b0, b1, b2, b0, b1, b2, b0, b1, b2)         */
/*------------------------------------------------------------*/
void zz_pX_Transform_naive::forward_right(Vec<zz_p>& val, const zz_pX& P) const
{
    long idx = 0;
    val.SetLength(n);

    const zz_p * coeffs = P.rep.elts();

    for (long i = 0; i < m; i++)
    {
	for (long j = 0; j <= deg(P); j++)
	{
	    val[idx++] = coeffs[j];
	}
	for (long j = deg(P)+1; j < m; j++)
	{
	    val[idx++] = 0;
	}
    }

}

/*------------------------------------------------------------*/
/* reconstructs C=AB from a0 b0, ..., a{m-1} b{m-1}           */
/*------------------------------------------------------------*/
void zz_pX_Transform_naive::backward(zz_pX& P, const Vec<zz_p>& val) const
{
    long idx = 0;
    for (long i = 0; i < m; i++)
    {
	zz_p res = val[idx++];
	for (long j = 1; j <= i; j++)
	{
	    res += val[idx++];
	}
	SetCoeff(P, i, res);
    }
    
    for (long i = m; i < 2*m - 1; i++)
    {
	zz_p res = val[idx++];
	for (long j = 1; j <= (2*m - 2 - i); j++)
	{
	    res += val[idx++];
	}
	SetCoeff(P, i, res);
    }
}


