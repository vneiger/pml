#include "vec_lzz_p_extra.h"
#include "lzz_pX_extra.h"

NTL_CLIENT


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*         ALGORITHMS FOR POWER SERIES DIVISION               */
/* compute b/a mod t^m                                        */
/* (sligthly faster than b * (1/a mod t^m) mod t^m            */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* power series division, naive algorithm                     */
/* &x == &a not allowed                                       */
/* x = b/a mod t^m                                            */
/*------------------------------------------------------------*/
void PlainInvTruncMul(zz_pX& x, const zz_pX& b, const zz_pX& a, long m)
{
    // degree of input
    const long na = deg(a);
    const long nb = deg(b);

    if (na < 0) 
        Error("division by zero");

    // compute inverse of constant term of a,
    // and record whether this constant term is one
    const zz_p s = inv(ConstTerm(a));
    const long is_one = IsOne(s);

    // if a is constant, not much more to do
    if (na == 0) 
    {
        x.SetLength(std::min(m, nb+1));
        // if a is constant one, just copy
        if (is_one)
            for (long k = 0; k < x.rep.length(); ++k) 
                x[k] = b[k];
        // if a is constant non-one: multiply by inverse of constant term
        // TODO could be made faster with MulModPrecon
        // (but a constant is unusual case)
        else
            for (long k = 0; k < x.rep.length(); ++k) 
                mul(x[k], b[k], s);
        x.normalize();
        return;
    }

    // set the length of x to be m (degree m-1)
    x.SetLength(m);

    // truncate: x = s*b mod x^m
    long kk;
    for (kk = 0; kk < std::min(m, nb+1); ++kk) 
        x[kk] = b[kk];
    for (; kk < m; ++kk) 
        clear(x[kk]);

    mul(x[0], x[0], s);

    // iterate up to order m
    zz_p buf1, buf2;
    for (long k = 1; k < m; ++k) 
    {
        buf1 = x[k];
        for (long i = std::max(k - na, (long)0); i < k; ++i) 
        {
            mul(buf2, x[i], a[k-i]);
            sub(buf1, buf1, buf2);
        }
        if (!is_one) 
            mul(x[k], buf1, s);
        else
            x[k] = buf1;

    } 
    x.normalize();
}

/*------------------------------------------------------------*/
/* power series division, Newton iteration                    */
/* adapted from the Tellegen package, Lecerf-Schost           */
/* aliasing not allowed                                       */
/* x = b/a mod t^m                                            */
/*------------------------------------------------------------*/
void NewtonInvTruncMul(zz_pX& x, const zz_pX& b, const zz_pX& a, long m)
{ 

    x.SetMaxLength(m);
    long i, t;

    t = NextPowerOfTwo(m);

    fftRep R1(INIT_SIZE, t), R2(INIT_SIZE, t);
    zz_pX P1(INIT_SIZE, m / 2);

    // -3 seems to be a reasonbable choice for many p's
    long log2_newton = NextPowerOfTwo(NTL_zz_pX_NEWTON_CROSSOVER)-3;

    long k = 1L << log2_newton;
    PlainInvTrunc(x, a, k);  // assumes m large enough

    t = log2_newton;
    long l = 0;
    while (k < m) 
    {
        l = min(2 * k, m);

        if (l >= m) 
            break;

        TofftRep(R1, x, t+1);
        TofftRep(R2, a, t+1, 0, l-1); 
        mul(R2, R2, R1);
        FromfftRep(P1, R2, k, l-1);

        TofftRep(R2, P1, t+1);
        mul(R2, R2, R1);
        FromfftRep(P1, R2, 0, l-k-1);

        x.rep.SetLength(l);
        long y_len = P1.rep.length();
        for (i = k; i < l; i++) {
            if (i-k >= y_len)
                clear(x.rep[i]);
            else
                NTL::negate(x.rep[i], P1.rep[i-k]);
        }
        x.normalize();
        t++;
        k = l;
    }

    // k < m <= 2k, l = m, 2^(t+1) = 2k
    // if m = 2 deg(b) = deg(a) + deg(b) + 1
    // then deg(b) <= k, deg(a) < k

    if (deg(a) < k && deg(b) < k && l == m && 1)
    {
        fftRep R3(INIT_SIZE, t+1);
        zz_pX r0(INIT_SIZE, k);
        zz_pX r1(INIT_SIZE, l-k);
        TofftRep(R1, x, t+1);
        TofftRep(R2, b, t+1);  
        mul(R3, R1, R2);
        FromfftRep(r0, R3, 0, k-1); 

        TofftRep(R3, r0, t);
        TofftRep(R2, a, t);  
        mul(R3, R3, R2);
        FromfftRep(r1, R3, 0, 2*k-2); 
        sub(r1, r1, b);

        // not very useful, somehow
        if (l-k < NTL_zz_pX_MUL_CROSSOVER/2 && 0)
        {
            r1 = MulTrunc(trunc(r1, l-k), trunc(x, l-k), l-k);
        }
        else
        {
            TofftRep(R3, r1, t+1);
            mul(R3, R3, R1);
            FromfftRep(r1, R3, 0, l-k-1); 
        }

        x.rep.SetLength(l);
        for (i = 0; i < k; i++) 
            x.rep[i] = coeff(r0, i);

        for (i = k; i < l; i++) 
            x.rep[i] = -coeff(r1, i-k);

        x.normalize();
    }
    else
    {
        fftRep R3(INIT_SIZE, t);
        zz_pX P2l(INIT_SIZE, k);
        zz_pX P2h(INIT_SIZE, l-k);

        TofftRep(R1, x, t+1);

        if (deg(b) >= k)
        {
            TofftRep(R2, b, t+1, k, l-1);  
            mul(R3, R1, R2);
            FromfftRep(P2h, R3, 0, l-k-1); // high part
        }

        TofftRep(R2, b, t+1, 0, k-1);  
        mul(R3, R1, R2);
        FromfftRep(P2l, R3, 0, l-1); // low part

        TofftRep(R2, a, t+1, 0, l-1); 
        mul(R2, R1, R2);
        FromfftRep(P1, R2, k, l-1);

        /*   P1 = P1*P2l */
        TofftRep(R2, P1, t+1, 0, l-k-1); 
        TofftRep(R3, P2l, t+1, 0, l-k-1); 
        mul(R1, R2, R3);
        FromfftRep(P1, R1, 0, l-k-1); 

        x.rep.SetLength(l);
        for (i = 0; i < k; i++) 
            x.rep[i]=P2l.rep[i];

        if (deg(b) >= k)
            for (i = k; i < l; i++) 
                add(x.rep[i], P2l.rep[i], P2h.rep[i-k]);
        else
            for (i = k; i < l; i++) 
                x.rep[i] = P2l.rep[i];

        for (i = k; i < l; i++) 
            sub(x.rep[i], x.rep[i], P1.rep[i-k]);

        x.normalize();
    }
}

/*------------------------------------------------------------*/
/* power series division, main function                       */
/* x = b/a mod t^m                                            */
/*------------------------------------------------------------*/
void InvTruncMul(zz_pX& x, const zz_pX& b, const zz_pX& a, long m)
{
    if (m < 0) 
        LogicError("InvTruncMul: precision must be >= 0");

    if (m == 0) 
    {
        clear(x);
        return;
    }

    if (&x == &a || &x == &b)
    {
        zz_pX f;
        InvTruncMul(f, b, a, m);
        x.swap(f);
        return;
    }

    if (m > NTL_zz_pX_NEWTON_CROSSOVER && deg(a) > 0)
        NewtonInvTruncMul(x, b, a, m);
    else
        PlainInvTruncMul(x, b, a, m);
}    


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* a class that does Taylor shift                             */
/* no assumption on the base ring, DAC algorithm              */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* constructor inits a few arrays                             */
/* d is a (nonstrict) upper bound on the input degree         */
/*------------------------------------------------------------*/
zz_pX_shift_DAC::zz_pX_shift_DAC(long d_in, const zz_p& c_in)
{
    d = d_in;
    zz_pX pow;

    c = c_in;
    cc = c * c;
    c3 = 3 * c;

    SetX(pow);
    SetCoeff(pow, 0, c); // pow = x+c
    pow = pow * pow;
    pow = pow * pow;  // pow = (x+c)^4

    precomp.SetLength(1);
    precomp[0] = pow;

    while ( 2 *deg(pow) - 1 < d )
    {
        pow = pow * pow;
        precomp.append(pow);
    }
}

/*------------------------------------------------------------*/
/* g = f(x+c)                                                 */
/* output can alias input                                     */
/*------------------------------------------------------------*/
void zz_pX_shift_DAC::shift(zz_pX& g, const zz_pX& f) const
{
    if (d == -1)
    {
        g = 0;
        return;
    }

    if (deg(f) > d)
        LogicError("Degree too large for shift");

    long lf = d + 1;
    long lv = lf >> 2;
    if (lf > 4 * lv)
        lv++;

    Vec<zz_pX> v;
    v.SetLength(lv);
    for (long i = 0, j = 0; i <= d; i+=4, j++)
    {
        SetCoeff(v[j], 0, coeff(f, i) + c * coeff(f, i+1) + cc * (coeff(f, i+2) + c * coeff(f, i+3)));
        SetCoeff(v[j], 1, coeff(f, i+1) + c * (2*coeff(f, i+2) + c3 * coeff(f, i+3)));
        SetCoeff(v[j], 2, coeff(f, i+2) + c3 * coeff(f, i+3));
        SetCoeff(v[j], 3, coeff(f, i+3));
    }
    long idx = 0;

    while (v.length() > 1)
    {
        long i, j;
        for (i = 0, j = 0; j+1 < v.length(); i++, j+=2)
            v[i] = v[j] + precomp[idx] * v[j+1];
        if (j == v.length() - 1)
        {
            v[i] = v[j];
            i++;
        }
        idx++;
        v.SetLength(i);
    }
    g = v[0];
}


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* a class that does Taylor shift                             */
/* requires 1,...,d units mod p                               */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* constructor inits a few arrays                             */
/* d is an upper bound on the degrees of the inputs           */
/*------------------------------------------------------------*/
zz_pX_shift_large_characteristic::zz_pX_shift_large_characteristic(long d, const zz_p& c) 
{
    // TODO: assert that d+1 is a unit in zz_p
    this->d = d;
    fact.SetLength(d+1);
    fact[0] = to_zz_p(1);
    for (long i = 1; i <= d; i++)
        fact[i] = i * fact[i-1];
    inv(ifact, fact);

    v.rep.SetLength(d+1);  
    v.rep[0] = to_zz_p(1);
    for (long i = 1; i <= d; i++)
        v.rep[i] = c * v.rep[i-1];
    for (long i = 0; i <= d; i++)
        v.rep[i] *= ifact[i];
}

/*------------------------------------------------------------*/
/* g = f(x+c)                                                 */
/* output can alias input                                     */
/*------------------------------------------------------------*/
void zz_pX_shift_large_characteristic::shift(zz_pX& g, const zz_pX& f) const
{
    if (d == -1)
    {
        g = 0;
        return;
    }

    zz_pX u, w;
    u.rep.SetLength(d+1);
    for (long i = 0; i < f.rep.length(); i++)
        u.rep[d-i] = f.rep[i] * fact[i];
    u.normalize();

    w = trunc(u * v, d + 1);

    g.rep.SetLength(d+1);

    for (long i = 0; i < w.rep.length(); i++)
        g.rep[d-i] = w.rep[i] * ifact[d-i];
    for (long i = w.rep.length(); i <= d; i++)
        g.rep[d-i] = 0;
    g.normalize();
}

/*------------------------------------------------------------*/
/* returns a zz_pX_shift of the right type                    */
/*------------------------------------------------------------*/
std::unique_ptr<zz_pX_shift> get_shift(long d, const zz_p& c)
{
    zz_p prod;
    for (long i = 1; i <= d; i++)
        prod *= i;
    long tmp;
    if (InvModStatus(tmp, prod.LoopHole(), zz_p::modulus()) == 0) // gcd = 1
    {
        return std::unique_ptr<zz_pX_shift>(new zz_pX_shift_large_characteristic(d, c));
    }
    else
    {
        return std::unique_ptr<zz_pX_shift>(new zz_pX_shift_DAC(d, c));
    }
}


/*------------------------------------------------------------*/
/* g = f(x+c)                                                 */
/* output can alias input                                     */
/*------------------------------------------------------------*/
void shift(zz_pX& g, const zz_pX& f, const zz_p& c)
{
    std::unique_ptr<zz_pX_shift> s = get_shift(deg(f), c);
    s->shift(g, f);
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
