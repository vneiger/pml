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
/* decides whether to use FFT or not                          */
/*------------------------------------------------------------*/
void zz_pX_Multipoint_Geometric::decide_FFT()
{
    if (min_geometric_FFT() <= n &&  n < max_geometric_FFT())
    {
        unset_FFT_evaluate();
        unset_FFT_interpolate();
    }
    else
    {
        set_FFT_evaluate();
        set_FFT_interpolate();
    }
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

    if (FFT_feasible)
    {
        long k = NextPowerOfTwo(d + n);
        fftRep d_fft = fftRep(INIT_SIZE, k);
        TofftRep(d_fft, f, k, 0, d + n - 1);
        known_degrees.insert( pair<long, fftRep>(d, d_fft) );
    }
}

/*------------------------------------------------------------*/
/* constructor for geometric progressions                     */
/* we interpolate at r^(2*i), i=0..d-1                        */
/*------------------------------------------------------------*/
zz_pX_Multipoint_Geometric::zz_pX_Multipoint_Geometric(const zz_p& r, long d) :
    zz_pX_Multipoint_Geometric(r, to_zz_p(1), d)
{ }

/*------------------------------------------------------------*/
/* constructor for geometric progressions                     */
/* we interpolate at s*r^(2*i), i=0..d-1                      */
/*------------------------------------------------------------*/
zz_pX_Multipoint_Geometric::zz_pX_Multipoint_Geometric(const zz_p& r, const zz_p& s, long d)
{
    n = d;
    idx_k = NextPowerOfTwo(2*d - 1);

    // if n==1, build sequence of points 'pts'
    // (otherwise, the method get_point could try to evaluate the degree 1
    // polynomial X to obtain all sequence values: this will fail since
    // evaluate degree 1 polynomials is not allowed when n==1)
    if (n==1)
    {
        pts.SetLength(1);
        mul(pts[0], s, r);
    }

    // bad case: we have an FFT prime of low order
    if (is_FFT_prime() && idx_k > zz_pInfo->MaxRoot)
        FFT_feasible = 0;
    else
        FFT_feasible = 1;

    decide_FFT();

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
    xs.SetLength(d);
    zz_p pow_s = s;
    x[0] = to_zz_p(1);
    xs[0] = to_zz_p(1);
    tmp = inv_r;
    for (long i = 1; i < d; i++)
    {
        x[i] = x[i-1] * tmp;
        xs[i] = x[i] * pow_s;
        tmp *= inv_q;
        pow_s *= s;
    }

    w.SetLength(d);
    ws.SetLength(d);
    y.SetLength(d);
    z.SetLength(d);
    zs.SetLength(d);

    SetCoeff(g1, 0, 1);
    SetCoeff(g2, 0, 1);
    y[0] = 1;
    w[0] = 1;
    ws[0] = 1;
    z[0] = 1;
    zs[0] = 1;

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
    zz_p t = to_zz_p(1);

    zz_p invs = 1/s;
    zz_p pow_invs = invs;

    for (long i = 1; i < d; i++)
    {
        qq *= qk;   // prod q^i
        t *= inv_qk; // prod 1/q^i
        w[i] = w[i-1] * inv_diff[i]; // prod 1/(q^i-1)
        ws[i] = w[i] * pow_invs;
        tmp = qq * w[i]; // prod q^i/(q^i-1)
        SetCoeff(g2, i, tmp);

        if ((i & 1) == 1)
        {
            SetCoeff(g1, i, -tmp);
            y[i] = -prod_diff[i];
            z[i] = -w[i];
            zs[i] = -ws[i];
        }
        else
        {
            SetCoeff(g1, i, tmp);
            y[i] = prod_diff[i];
            z[i] = w[i];
            zs[i] = ws[i];
        }
        y[i] = y[i] * t;

        qk *= q;
        inv_qk *= inv_q;
        pow_invs *= invs;
    }

    // do it as soon as feasible, whether do_FFT_xxx is set or not
    if (FFT_feasible) 
    {
        prepare_degree(n - 1);
        g1_fft = fftRep(INIT_SIZE, idx_k);
        g2_fft = fftRep(INIT_SIZE, idx_k);
        TofftRep(g1_fft, g1, idx_k);
        TofftRep(g2_fft, g2, idx_k);
    }
}

/*------------------------------------------------------------*/
/* return the ratio q = r^2                                   */
/*------------------------------------------------------------*/
zz_p zz_pX_Multipoint_Geometric::get_q() const
{
    zz_p q;
    if (n == 1)
        q = to_zz_p(1);
    else
        q = 1/(x[1] * x[1]);
    return q;
}

/*------------------------------------------------------------*/
/* return s (points are s*r^(2i))                             */
/*------------------------------------------------------------*/
zz_p zz_pX_Multipoint_Geometric::get_s() const
{
    if (n == 1)
        return to_zz_p(1);
    else
        return xs[1] / x[1];
}

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* polynomial operations                                      */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* val[i] = P(s * r^(2*i)), i = 0..n-1                        */
/* val may alias P.rep                                        */
/*------------------------------------------------------------*/
void zz_pX_Multipoint_Geometric::evaluate(Vec<zz_p>& val, const zz_pX& P) const
{
    zz_pX a, b;
    long dp = deg(P);

    if (dp >= n)
    {
        Error("Degree too large for geometric evaluate.");
    }

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

    if (dp == -1)
    {
        for (long i = 0; i < n; i++)
        {
            val[i] = 0;
        }
        return;
    }

    if (FFT_feasible && do_FFT_evaluate)
    {
        // finds the smallest among all known degrees >= dp
        long dp_found = n - 1;
        for (map<int, fftRep>::const_iterator it = known_degrees.cbegin(); it != known_degrees.cend(); it++)
        {
            if (it->first >= dp && it->first < dp_found)
                dp_found = it->first;
        }
        dp = dp_found; 
    }

    for (long i = 0; i <= dp; i++)  
    {
        SetCoeff(a, dp - i, xs[i] * coeff(P, i));  
    }

    if (FFT_feasible && do_FFT_evaluate)
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
/* evaluates f at s * r^(2*i), i=0..nb-1                      */
/* still needs deg(f) < d                                     */ 
/*------------------------------------------------------------*/
void zz_pX_Multipoint_Geometric::evaluate(Vec<zz_p>& val, const zz_pX& f, long nb) const
{
    val.SetLength(nb + n);

    if (nb == 0) 
    {
        return;
    }
    
    long idx = 0;
    Vec<zz_p> tmp;
    zz_pX f_loc = f;
    zz_p pow = power(get_q(), n);
    do
    {
        evaluate(tmp, f_loc);
        for (long i = 0; i < n; i++)
            val[idx + i] = tmp[i];
        zz_p pow_loc = to_zz_p(1);
        for (long i = 0; i <= deg(f); i++)
        {
            SetCoeff(f_loc, i, coeff(f_loc, i) * pow_loc);
            pow_loc *= pow;
        }
        idx += n;
    }
    while (idx < nb);
    val.SetLength(nb);
}

/*------------------------------------------------------------*/
/* finds f of degree < n such that                            */
/* val[i] = f(s * r^(2*i)), i = 0..n-1                        */
/* val may alias f.rep                                        */
/*------------------------------------------------------------*/
void zz_pX_Multipoint_Geometric::interpolate(zz_pX& f, const Vec<zz_p>& val) 
{
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

    if (FFT_feasible && do_FFT_interpolate)
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

    if (FFT_feasible && do_FFT_interpolate)
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
        SetCoeff(f, i, coeff(h, n - 1 - i) * zs[i]);
    }

    f.normalize();
}

/*-----------------------------------------------------------*/
/* transpose of                                              */
/* val[i] = P(s * r^(2*i)), i = 0..n-1                       */
/* val must have length n                                    */
/* val may alias P.rep                                       */
/*-----------------------------------------------------------*/
void zz_pX_Multipoint_Geometric::t_evaluate(zz_pX& P, const Vec<zz_p>& val, long output_size) const
{
    P = 0;

    if (output_size == -1)
        output_size = n;

    if (output_size == 0)
        return;

    if (n == 0)
        return;

    if (n == 1)
    {
        SetCoeff(P, 0, val[0]);
        return;
    }

    zz_pX a;
    for (long i = 0; i < n; i++)  
        SetCoeff(a, n - 1 - i, x[i] * val[i]);  

    zz_pX b;
    if (FFT_feasible && do_FFT_evaluate)
    {
        // finds the smallest among all known sizes >= output_size
        long do_output_size = n;
        for (map<int, fftRep>::const_iterator it = known_degrees.cbegin(); it != known_degrees.cend(); it++)
        {
            if ((it->first + 1) >= output_size && (it->first + 1) < do_output_size)
                do_output_size = it->first + 1;
        }
        const fftRep& f_fft = known_degrees.find(do_output_size - 1)->second;
        fftRep a_fft, b_fft;
        long k = f_fft.k;
        a_fft = fftRep(INIT_SIZE, k);
        b_fft = fftRep(INIT_SIZE, k);
        TofftRep(a_fft, a, k, 0, n - 1);
        mul(b_fft, a_fft, f_fft);
        FromfftRep(b, b_fft, n - 1, n - 1 + do_output_size - 1);  
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
        b = middle_product(a, f, n - 1, output_size - 1);
    }

    for (long i = 0; i < output_size; i++)
        SetCoeff(P, i, xs[i] * coeff(b, i));

}

/*------------------------------------------------------------*/
/* transpose of                                               */
/* finds P of degree < n such that                            */
/* val[i] = P(s * r^(2*i)), i = 0..n-1                        */
/* val may alias P.rep                                        */
/*------------------------------------------------------------*/
void zz_pX_Multipoint_Geometric::t_interpolate(Vec<zz_p>& val, const zz_pX& P) 
{
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

    zz_pX k, h;

    for (long i = 0; i < n; i++)
    {
        SetCoeff(k, i, coeff(P, i) * ws[i]);
    }

    if (FFT_feasible && do_FFT_interpolate)
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

    if (FFT_feasible && do_FFT_interpolate)
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
        val[i] = coeff(h, n - 1 - i) * z[i];
    }
}


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* vector operations                                          */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*-----------------------------------------------------------*/
/* val[i] = P(s * r^(2*i)), i = 0..n-1                       */
/* P must have length n                                      */
/* val may alias P                                           */
/*-----------------------------------------------------------*/
void zz_pX_Multipoint_Geometric::mul_right(Vec<zz_p>& val, const Vec<zz_p>& P) const
{
    val.SetLength(n);
    
    if (n == 0)
    {
        return;
    }
    
    if (n == 1)
    {
        val[0] = P[0];
        return;
    }
    
    zz_pX a, b;
    for (long i = 0; i < n; i++)  
    {
        SetCoeff(a, n - 1 - i, xs[i] * P[i]);
    }
            
    if (FFT_feasible && do_FFT_evaluate)
    {
        const fftRep& f_fft = known_degrees.find(n-1)->second;
        fftRep a_fft, b_fft;
        long k = f_fft.k;

        a_fft = fftRep(INIT_SIZE, k);
        b_fft = fftRep(INIT_SIZE, k);
        TofftRep(a_fft, a, k, 0, n - 1);
        mul(b_fft, a_fft, f_fft);
        FromfftRep(b, b_fft, n - 1, n - 1 + n - 1);  
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
        b = middle_product(a, f, n - 1, n - 1);
    }

    for (long i = 0; i < n; i++)
    {
        val[i] = x[i] * coeff(b, i);
    }
}
    
/*------------------------------------------------------------*/
/* finds f of degree < n such that                            */
/* val[i] = f(s * r^(2*i)), i = 0..n-1                        */
/* val may alias f                                            */
/*------------------------------------------------------------*/
void zz_pX_Multipoint_Geometric::inv_mul_right(Vec<zz_p>& f, const Vec<zz_p>& val) 
{
    f.SetLength(n);

    if (n == 0)
    {
        return;
    }

    if (n == 1)
    {
        f[0] = val[0];
        return;
    }

    zz_pX k, h;

    for (long i = 0; i < n; i++)
    {
        SetCoeff(k, i, val[i] * w[i]);
    }

    if (FFT_feasible && do_FFT_interpolate)
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

    if (FFT_feasible && do_FFT_interpolate)
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
        f[i] = coeff(h, n - 1 - i) * zs[i];
    }
}

/*-----------------------------------------------------------*/
/* transpose of                                              */
/* val[i] = P(s * r^(2*i)), i = 0..n-1                       */
/* val must have length n                                    */
/* val may alias P                                           */
/*-----------------------------------------------------------*/
void zz_pX_Multipoint_Geometric::mul_left(Vec<zz_p>& P, const Vec<zz_p>& val) const
{
    P.SetLength(n);

    if (n == 0)
    {
        return;
    }

    if (n == 1)
    {
        P[0] = val[0];
        return;
    }

    zz_pX a;
    for (long i = 0; i < n; i++)  
        SetCoeff(a, n - 1 - i, x[i] * val[i]);  

    zz_pX b;
    if (FFT_feasible && do_FFT_evaluate)
    {
        const fftRep& f_fft = known_degrees.find(n - 1)->second;
        fftRep a_fft, b_fft;
        long k = f_fft.k;
        a_fft = fftRep(INIT_SIZE, k);
        b_fft = fftRep(INIT_SIZE, k);
        TofftRep(a_fft, a, k, 0, n - 1);
        mul(b_fft, a_fft, f_fft);
        FromfftRep(b, b_fft, n - 1, n - 1 + n - 1);  
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
        b = middle_product(a, f, n - 1, n - 1);
    }

    for (long i = 0; i < n; i++)
        P[i] = xs[i] * coeff(b, i);

}

/*------------------------------------------------------------*/
/* transpose of                                               */
/* finds P of degree < n such that                            */
/* val[i] = P(s * r^(2*i)), i = 0..n-1                        */
/* val may alias P                                            */
/*------------------------------------------------------------*/
void zz_pX_Multipoint_Geometric::inv_mul_left(Vec<zz_p>& val, const Vec<zz_p>& P) 
{
    val.SetLength(n);

    if (n == 0)
    {
        return;
    }

    if (n == 1)
    {
        val[0] = P[0];
        return;
    }

    zz_pX k, h;

    for (long i = 0; i < n; i++)
    {
        SetCoeff(k, i, P[i] * ws[i]);
    }

    if (FFT_feasible && do_FFT_interpolate)
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

    if (FFT_feasible && do_FFT_interpolate)
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
        val[i] = coeff(h, n - 1 - i) * z[i];
    }
}


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* Misc:                                                      */
/* returns a zz_pX_Multipoint_Geometric with n points         */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
zz_pX_Multipoint_Geometric get_geometric_points(long n)
{
    zz_p r;
    element_of_order(r, 2*n);
    return zz_pX_Multipoint_Geometric(r, n);
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
