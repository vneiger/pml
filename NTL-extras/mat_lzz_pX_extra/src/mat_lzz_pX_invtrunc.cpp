#include <NTL/matrix.h>
#include <NTL/mat_lzz_p.h>
#include <NTL/lzz_pX.h>

#include "lzz_p_extra.h"
#include "mat_lzz_pX_extra.h"
#include "lzz_pX_CRT.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
void plain_inv_trunc(Mat<zz_pX>& x, const Mat<zz_pX>& a, long m)
{
    if (x == a)
    {
        Mat<zz_pX> y;
        plain_inv_trunc(y, a, m);
        x = y;
        return;
    }

    long n, u;
    Mat<zz_p> cst_mat, inv0, v, xi, ai, t;

    u = a.NumRows();
    if (u != a.NumCols())
    {
        LogicError("Non square matrix for truncated inverse\n");
    }

    cst_mat = constant_coefficient(a);
    inv0 = inv(cst_mat);

    n = deg(a);
    if (n == 0) 
    {
        conv(x, inv0);
        return;
    }

    // ap = a.rep.elts();
    // x.rep.SetLength(m);
    // xp = x.rep.elts();

    x.SetDims(u, u);
    for (long r = 0; r < u; r++)
        for (long s = 0; s < u; s++)
        {
            x[r][s].rep.SetLength(m);
            x[r][s].rep[0] = inv0[0][0];
        }


    v.SetDims(u, u);
    for (long k = 1; k < m; k++) 
    {
        clear(v);
        long lb = max(k - n, 0);
        for (long i = lb; i <= k-1; i++) 
        {
            GetCoeff(xi, x, i);
            GetCoeff(ai, a, k - i);
            mul(t, xi, ai);
            add(v, v, t);
        }
        v = inv0 * v;
        for (long r = 0; r < u; r++)
            for (long s = 0; s < u; s++)
                x[r][s].rep[k] = -v[r][s];
    }

    for (long r = 0; r < u; r++)
        for (long s = 0; s < u; s++)
        {
            x[r][s].normalize();
        }
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
