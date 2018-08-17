#include <NTL/lzz_pX.h>
#include <NTL/vector.h>

#include "util.h"
#include "lzz_p_extra.h"
#include "lzz_pX_middle_product.h"
#include "lzz_pX_CRT.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* multipoint evaluation over zz_p                            */
/* points in geometric progression                            */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* getters / setters                                          */
/*------------------------------------------------------------*/
long zz_pX_Multipoint_Geometric::FFT_evaluate() const
{
    return do_FFT_evaluate;
}


long zz_pX_Multipoint_Geometric::FFT_interpolate() const
{
    return do_FFT_interpolate;
}

void zz_pX_Multipoint_Geometric::set_FFT_evaluate()
{
    do_FFT_evaluate = 1;
}

void zz_pX_Multipoint_Geometric::unset_FFT_evaluate()
{
    do_FFT_evaluate = 0;
}

void zz_pX_Multipoint_Geometric::set_FFT_interpolate()
{
    do_FFT_interpolate = 1;
}

void zz_pX_Multipoint_Geometric::unset_FFT_interpolate()
{
    do_FFT_interpolate = 0;
}


/*------------------------------------------------------------*/
/* adds a new FFT for repeated evaluations in degree d        */
/*------------------------------------------------------------*/
void zz_pX_Multipoint_Geometric::prepare_degree(long d)
{
    if (d >= n)
    {
        LogicError("Degree too large for evaluation.");
    }

    long k = NextPowerOfTwo(d + n);
    fftRep d_fft = fftRep(INIT_SIZE, k);
    TofftRep(d_fft, f, k, 0, d + n - 1);

    known_degrees.insert( pair<long, fftRep>(d, d_fft) );
}

/*------------------------------------------------------------*/
/* constructor for geometric progressions                     */
/* we interpolate at r^(2*i), i=0..d-1                        */
/*------------------------------------------------------------*/
zz_pX_Multipoint_Geometric::zz_pX_Multipoint_Geometric(const zz_p& r, long d){
    
    n = d;
    idx_k = NextPowerOfTwo(2*d - 1);

    long cross = zz_pX_mul_crossover[zz_pInfo->PrimeCnt];
    do_FFT_evaluate = (d > 100);
    if (cross <= 500)
        do_FFT_evaluate = (d >= 75);
    if (cross <= 150)
        do_FFT_evaluate = 1;

    do_FFT_interpolate = (d > 300);
    if (cross <= 500)
        do_FFT_interpolate = (d >= 200);
    if (cross <= 150)
        do_FFT_interpolate = 1;

    zz_p q = r*r;
    zz_p inv_r = 1/r;
    zz_p inv_q = inv_r*inv_r;

    // polynomial for evaluate
    SetCoeff(f, 0, 1);
    zz_p tmp = r;
    for (long i = 1; i < 2*d-1; i++)
    {
        SetCoeff(f, i, coeff(f, i-1) * tmp);
        tmp *= q;
    }

    // array for evaluate
    x.SetLength(d);
    x[0] = to_zz_p(1);
    tmp = inv_r;
    for (long i = 1; i < d; i++)
    {
        x[i] = x[i-1] * tmp;
        tmp *= inv_q;
    }

    w.SetLength(d);
    y.SetLength(d);
    z.SetLength(d);

    SetCoeff(g1, 0, 1);
    SetCoeff(g2, 0, 1);
    y[0] = 1;
    w[0] = 1;
    z[0] = 1;
    
    inv_q = 1/q;

    Vec<zz_p> inv_diff, diff, prod_diff;
    inv_diff.SetLength(d);
    diff.SetLength(d);
    prod_diff.SetLength(d);

    inv_diff[0] = 1;
    diff[0] = 1;
    prod_diff[0] = 1;

    zz_p qk = q;  // montgomery inversion
    for (long i = 1; i < d; i++)
    {
        diff[i] = qk - 1;
        inv_diff[i] = diff[i];
        qk *= q;
        prod_diff[i] = diff[i] * prod_diff[i-1];
    }
    
    tmp = 1 / prod_diff[d-1];
    for (long i = d-1; i > 0; i--)
    {
        inv_diff[i] = prod_diff[i-1] * tmp;
        tmp *= diff[i];
    }
    inv_diff[0] = tmp;
    // end montgomery inversion

    // sets sequences w, y, z and polynomials g1, g2
    qk = 1;
    zz_p inv_qk = to_zz_p(1);
    zz_p qq = to_zz_p(1);
    zz_p s = to_zz_p(1);
    for (long i = 1; i < d; i++)
    {
        qq *= qk;   // prod q^i
        s *= inv_qk; // prod 1/q^i
        w[i] = w[i-1] * inv_diff[i]; // prod 1/(q^i-1)
        tmp = qq * w[i]; // prod q^i/(q^i-1)
        SetCoeff(g2, i, tmp);

        if ((i & 1) == 1)
        {
            SetCoeff(g1, i, -tmp);
            y[i] = -prod_diff[i];
            z[i] = -w[i];
        }
        else
        {
            SetCoeff(g1, i, tmp);
            y[i] = prod_diff[i];
            z[i] = w[i];
        }
        y[i] = y[i] * s;

        qk *= q;
        inv_qk *= inv_q;
    }

    prepare_degree(n - 1);
    // f_fft = fftRep(INIT_SIZE, idx_k);
    // TofftRep(f_fft, f, idx_k);

    g1_fft = fftRep(INIT_SIZE, idx_k);
    g2_fft = fftRep(INIT_SIZE, idx_k);
    TofftRep(g1_fft, g1, idx_k);
    TofftRep(g2_fft, g2, idx_k);
}

/*-----------------------------------------------------------*/
/* val[i] = P(r^(2*i)), i = 0..n-1                           */
/*-----------------------------------------------------------*/
void zz_pX_Multipoint_Geometric::evaluate(Vec<zz_p>& val, const zz_pX& P) const{
    val.SetLength(n);

    long dp = deg(P);

    if (n == 0)
    {
        return;
    }

    if (n == 1)
    {
        val[0] = coeff(P, 0);
        return;
    }

    if (dp == -1)
    {
        for (long i = 0; i < n; i++)
        {
            val[i] = 0;
        }
        return;
    }

    if ( known_degrees.find(dp) == known_degrees.end() )
    {
        dp = n - 1;
    }

    zz_pX a;
    for (long i = 0; i <= dp; i++)   // was n-1
    {
        SetCoeff(a, dp - i, x[i] * coeff(P, i));  //was: n-1-i
    }
   
    zz_pX b;
    if (do_FFT_evaluate)
    {
        const fftRep& f_fft = known_degrees.find(dp)->second;
        fftRep a_fft, b_fft;
        long k = f_fft.k;
        a_fft = fftRep(INIT_SIZE, k);
        b_fft = fftRep(INIT_SIZE, k);
        TofftRep(a_fft, a, k, 0, dp);
        mul(b_fft, a_fft, f_fft);
        FromfftRep(b, b_fft, dp, dp + n - 1);  
// for k = 1, the normalization is different in version 11.1.0
#ifdef __NTL_FIX_SIZE_2_FFT
        if (k == 1)
        {
            b /= 2;
        }
#endif
    }
    else
    {
        b = middle_product(a, f, dp, n - 1);
    }

    for (long i = 0; i < n; i++)
    {
        val[i] = x[i] * coeff(b, i);
    }
}

/*------------------------------------------------------------*/
/* finds f of degree < d such that                            */
/* val[i] = f(r^(2*i)), i = 0..n-1                            */
/*------------------------------------------------------------*/
void zz_pX_Multipoint_Geometric::interpolate(zz_pX& f, const Vec<zz_p>& val) {
    f = 0;

    if (n == 0)
    {
        return;
    }

    if (n == 1)
    {
        f.rep.SetLength(1);
        f.rep[0] = val[0];
        f.normalize();
        return;
    }

    zz_pX k, h;
    
    for (long i = 0; i < n; i++)
    {
        SetCoeff(k, i, val[i] * w[i]);
    }

    if (do_FFT_interpolate)
    {
        fftRep k_fft, h_fft;
        k_fft = fftRep(INIT_SIZE, idx_k);
        h_fft = fftRep(INIT_SIZE, idx_k);
        TofftRep(k_fft, k, idx_k);
        mul(h_fft, k_fft, g1_fft);
        FromfftRep(h, h_fft, 0, n-1);
    }
    else
    {
        h = MulTrunc(k, g1, n);
    }

    for (long i = 0; i < n; i++)
    {
        SetCoeff(k, n - 1 - i, coeff(h, i) * y[i]);
    }

    if (do_FFT_interpolate)
    {
        fftRep k_fft, h_fft;
        k_fft = fftRep(INIT_SIZE, idx_k);
        h_fft = fftRep(INIT_SIZE, idx_k);
        TofftRep(k_fft, k, idx_k);
        mul(h_fft, k_fft, g2_fft);
        FromfftRep(h, h_fft, 0, n-1);
    }
    else
    {
        h = MulTrunc(k, g2, n);
    }

    for (long i = 0; i < n; i++)
    {
        SetCoeff(f, i, coeff(h, n - 1 - i) * z[i]);
    }

    f.normalize();
}

/*-----------------------------------------------------------*/
/* transpose of                                              */
/* val[i] = P(r^(2*i)), i = 0..n-1                           */
/* val must have length n                                    */
/*-----------------------------------------------------------*/
void zz_pX_Multipoint_Geometric::t_evaluate(zz_pX& P, const Vec<zz_p>& val) const{
    zz_pX Q;
    for (long i = 0; i < n; i++)
        SetCoeff(Q, i, val[i]);

    Vec<zz_p> res;
    evaluate(res, Q);

    for (long i = 0; i < n; i++)
        SetCoeff(P, i, res[i]);
}

/*------------------------------------------------------------*/
/* transpose of                                               */
/* finds P of degree < d such that                            */
/* val[i] = P(r^(2*i)), i = 0..n-1                            */
/*------------------------------------------------------------*/
void zz_pX_Multipoint_Geometric::t_interpolate(Vec<zz_p>& val, const zz_pX& P) {
    val.SetLength(n);

    if (n == 0)
    {
        return;
    }

    if (n == 1)
    {
        val[0] = coeff(P, 0);
        return;
    }

    Vec<zz_p> coeffs;
    coeffs.SetLength(n);

    for (long i = 0; i < n; i++)
        coeffs[i] = coeff(P, i);

    zz_pX Q;
    interpolate(Q, coeffs);

    for (long i = 0; i < n; i++)
        val[i] = coeff(Q, i);

}

/*------------------------------------------------------------*/
/* returns a zz_pX_Multipoint_Geometric with n points         */
/*------------------------------------------------------------*/
zz_pX_Multipoint_Geometric get_geometric_points(long n)
{
    zz_p r;
    element_of_order(r, 2*n);
    return zz_pX_Multipoint_Geometric(r, n);
}
