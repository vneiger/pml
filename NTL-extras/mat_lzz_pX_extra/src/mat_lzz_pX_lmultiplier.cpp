#include <NTL/matrix.h>
#include <NTL/mat_lzz_p.h>
#include <NTL/lzz_pX.h>

#include "lzz_p_extra.h"
#include "mat_lzz_pX_extra.h"
#include "lzz_pX_CRT.h"

NTL_CLIENT

long mat_lzz_pX_lmultiplier::NumRows() const
{
    return __s;
}

long mat_lzz_pX_lmultiplier::NumCols() const
{
    return __t;
}

long mat_lzz_pX_lmultiplier::degA() const
{
    return __dA;
}

long mat_lzz_pX_lmultiplier::degB() const
{
    return __dB;
}



mat_lzz_pX_lmultiplier_FFT::mat_lzz_pX_lmultiplier_FFT(const Mat<zz_pX> & a, long dB)
{
    __s = a.NumRows();
    __t = a.NumCols();
    __dA = deg(a);
    __dB = dB;
    
    idxk = NextPowerOfTwo(__dA + __dB + 1);
    long n = 1L << idxk;
    fftRep R1(INIT_SIZE, idxk);

    va.SetLength(n);
    for (long i = 0; i < n; i++)
        va[i].SetDims(__s, __t);
    
    
    long st = __s * __t;
    for (long i = 0; i < __s; i++)
    {
        for (long k = 0; k < __t; k++)
        {
            TofftRep(R1, a[i][k], idxk);
            long *frept = & R1.tbl[0][0];
            for (long r = 0, rst = 0; r < n; r++, rst += st)
                va[r][i][k] = frept[r];
        }
    }

}


void mat_lzz_pX_lmultiplier_FFT::multiply(Mat<zz_pX>& c, const Mat<zz_pX>& b) const
{
    long s = NumRows();
    long t = NumCols();
    long u = b.NumCols();

    long dB = deg(b);

    if (dB > degB())
    {
        LogicError("Rhs degree too large in multiplier");
    }

    fftRep R1(INIT_SIZE, idxk);
    long n = 1L << idxk;
    Vec<zz_p> mat_valB;
    Vec<Vec<zz_p>> mat_valC;
    
    mat_valB.SetLength(n * t * u);
    
    long st = s*t;
    long tu = t*u;
    for (long i = 0; i < t; i++)
    {
        for (long k = 0; k < u; k++)
        {
            TofftRep(R1, b[i][k], idxk);
            long *frept = & R1.tbl[0][0];
            for (long r = 0, rtu = 0; r < n; r++, rtu += tu)
                mat_valB[rtu + i*u + k] = frept[r];
        }
    }

    Mat<zz_p> vb, vc;
    vb.SetDims(t, u);
    
    mat_valC.SetLength(s * u);
    for (long i = 0; i < s * u; i++)
        mat_valC[i].SetLength(n);
    
    for (long j = 0, jst = 0, jtu = 0; j < n; j++, jst += st, jtu += tu)
    {
        for (long i = 0; i < t; i++)
            for (long k = 0; k < u; k++)
                vb[i][k] = mat_valB[jtu + i*u + k];
        
        vc = va[j] * vb;
        
        for (long i = 0; i < s; i++)
            for (long k = 0; k < u; k++)
                mat_valC[i*u + k][j] = vc[i][k];
    }
    
    c.SetDims(s, u);
    for (long i = 0; i < s; i++)
        for (long k = 0; k < u; k++)
        {
            long *frept = & R1.tbl[0][0];
            for (long r = 0; r < n; r++)
                frept[r] = rep(mat_valC[i*u + k][r]);
            FromfftRep(c[i][k], R1, 0, n - 1);
        }
}


void mat_lzz_pX_lmultiplier_geometric::multiply(Mat<zz_pX>& res, const Mat<zz_pX>& rhs) const
{


}

void mat_lzz_pX_lmultiplier_3_primes::multiply(Mat<zz_pX>& res, const Mat<zz_pX>& rhs) const
{


}
