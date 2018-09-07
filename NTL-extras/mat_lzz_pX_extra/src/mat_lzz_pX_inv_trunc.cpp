#include <NTL/matrix.h>
#include <NTL/mat_lzz_p.h>
#include <NTL/lzz_pX.h>

#include "lzz_p_extra.h"
#include "mat_lzz_pX_extra.h"
#include "lzz_pX_CRT.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* returns x = 1/a mod z^m, naive algorithm                   */
/* throws an error if a(0) not invertible                     */
/* x can alias a                                              */
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

    x.SetDims(u, u);
    n = deg(a);
    if (n == 0) 
    {
        conv(x, inv0);
        return;
    }
    
    for (long r = 0; r < u; r++)
        for (long s = 0; s < u; s++)
        {
            x[r][s].rep.SetLength(m);
            x[r][s].rep[0] = inv0[r][s];
        }


    v.SetDims(u, u);
    for (long k = 1; k < m; k++) 
    {
        clear(v);
        for (long i = 0; i <= k-1; i++) 
        {
            GetCoeff(xi, x, i);
            GetCoeff(ai, a, k - i);
            mul(t, xi, ai);
            add(v, v, t);
        }
        v = v * inv0;
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


/*------------------------------------------------------------*/
/* returns x = 1/a mod z^m, Newton iteration                  */
/* throws an error if a(0) not invertible                     */
/* x can alias a                                              */
/*------------------------------------------------------------*/
void newton_inv_trunc_FFT(Mat<zz_pX>& x, const Mat<zz_pX>& a, long m)
{
    if (x == a)
    {
        Mat<zz_pX> y;
        newton_inv_trunc_FFT(y, a, m);
        x = y;
        return;
    }

    const long thresh = 10;
    
    long n, k, s, ss, t;
    Mat<zz_pX> P1;
    Vec<zz_p> mat_val1, mat_val2;
    Vec<Vec<zz_p>> mat_val3;
    Mat<zz_p> v1, v2, v3;

    fftRep R1(INIT_SIZE, NextPowerOfTwo(2*m-1));

    k = 1L << ( NextPowerOfTwo(thresh)-1 );
    plain_inv_trunc(x, a, k);    
    s = a.NumRows();
    ss = s * s;

    v1.SetDims(s, s);
    v2.SetDims(s, s);
    v3.SetDims(s, s);
    P1.SetDims(s, s);
    x.SetDims(s, s);

    while (k < m) 
    {
        t = NextPowerOfTwo(2 * k);
        n = 1L << t;

        mat_val1.SetLength(n * ss);
        mat_val2.SetLength(n * ss);
        mat_val3.SetLength(ss);
        for (long i = 0; i < ss; i++)
            mat_val3[i].SetLength(n);

        // mat_val1 = FFT of x
        // deg(x) <= k-1
        for (long i = 0; i < s; i++)
            for (long j = 0; j < s; j++)
            {
                TofftRep(R1, x[i][j], t);
                long *frept = & R1.tbl[0][0];
                for (long r = 0, rss = 0; r < n; r++, rss += ss)
                {
                    mat_val1[rss + i*s + j] = frept[r];
                }
            }

       
        // mat_val2 = FFT of (a mod t^(2k))
        // deg(a_ij mod t^(2k)) <= min(2k-1, deg(a_ij))
        for (long i = 0; i < s; i++)
            for (long j = 0; j < s; j++)
            {
                TofftRep(R1, a[i][j], t, 0, min(2*k - 1, deg(a[i][j])));
                long *frept = & R1.tbl[0][0];
                for (long r = 0, rss = 0; r < n; r++, rss += ss)
                    mat_val2[rss + i*s + j] = frept[r];
            }

        // mat_val3 = values of ( x (a mod t^2k) ) = ___ + t^k delta
        for (long j = 0, jss = 0; j < n; j++, jss += ss)
        {
            for (long i = 0; i < s; i++)
                for (long ell = 0; ell < s; ell++)
                    v1[i][ell] = mat_val1[jss + i*s + ell];
            for (long i = 0; i < s; i++)
                for (long ell = 0; ell < s; ell++)
                    v2[i][ell] = mat_val2[jss + i*s + ell];
            v3 = v1 * v2;
            for (long i = 0; i < s; i++)
                for (long ell = 0; ell < s; ell++)
                    mat_val3[i*s + ell][j] = v3[i][ell];
        }


        // mat_val2 = values of delta
        for (long i = 0; i < s; i++)
            for (long ell = 0; ell < s; ell++)
            {
                zz_pX P;
                long *frept = & R1.tbl[0][0];
                for (long r = 0; r < n; r++)
                    frept[r] = rep(mat_val3[i*s + ell][r]);
                FromfftRep(P, R1, k, 2*k - 1);
                TofftRep(R1, P, t);
                frept = & R1.tbl[0][0];
                for (long r = 0, rss = 0; r < n; r++, rss += ss)
                    mat_val2[rss + i*s + ell] = frept[r];
            }

        // mat_val3 = values of (delta x)
        for (long j = 0, jss = 0; j < n; j++, jss += ss)
        {
            for (long i = 0; i < s; i++)
                for (long ell = 0; ell < s; ell++)
                    v1[i][ell] = mat_val1[jss + i*s + ell];
            for (long i = 0; i < s; i++)
                for (long ell = 0; ell < s; ell++)
                    v2[i][ell] = mat_val2[jss + i*s + ell];
            v3 = v2 * v1;
            for (long i = 0; i < s; i++)
                for (long ell = 0; ell < s; ell++)
                    mat_val3[i*s + ell][j] = v3[i][ell];
        }

        for (long i = 0; i < s; i++)
            for (long ell = 0; ell < s; ell++)
            {
                zz_pX P;
                long *frept = & R1.tbl[0][0];
                for (long r = 0; r < n; r++)
                    frept[r] = rep(mat_val3[i*s + ell][r]);
                FromfftRep(P, R1, 0, k - 1);

                x[i][ell].rep.SetLength(2*k);
                long y_len = P.rep.length();
                for (long ii = k; ii < 2*k; ii++) 
                {
                    if (ii-k >= y_len)
                        clear(x[i][ell].rep[ii]);
                    else
                        x[i][ell].rep[ii] = -P.rep[ii - k];
                }
                x[i][ell].normalize();
            }
        
        k = 2*k;
    }

    trunc(x, x, m);
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
