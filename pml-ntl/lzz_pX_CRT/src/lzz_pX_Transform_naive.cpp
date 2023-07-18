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
/* (a0,a1,a2) -> (a0, a0, a1, a0, a1, a2, a1, a2, a2)         */
/*------------------------------------------------------------*/
void zz_pX_Transform_naive::forward_left(Vec<zz_p>& val, const zz_pX& P) const
{
    long idx = 0;
    val.SetLength(n);

    for (long i = 0; i < m; i++)
    {
        for (long j = i; j >= 0; j--)
        {
            val[idx++] = coeff(P, j);
        }
    }

    for (long i = 1; i < m; i++)
    {
        for (long j = m-1; j >= i; j--)
        {
            val[idx++] = coeff(P, j);
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

    for (long i = 0; i < m; i++)
    {
        for (long j = 0; j <= i; j++)
        {
            val[idx++] = coeff(P, j);
        }
    }

    for (long i = 1; i < m; i++)
    {
        for (long j = i; j <= m-1; j++)
        {
            val[idx++] = coeff(P, j);
        }
    }
}

/*------------------------------------------------------------*/
/* reconstructs C=AB from a0 b0, ..., a{m-1} b{m-1}           */
/*------------------------------------------------------------*/
void zz_pX_Transform_naive::backward(zz_pX& P, const Vec<zz_p>& val) const
{
    P = 0;
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



// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
