#include <algorithm>

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
static void PlainInvTruncMul(zz_pX& x, const zz_pX& b, const zz_pX& a, long m)
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
static void NewtonInvTruncMul(zz_pX& x, const zz_pX& b, const zz_pX& a, long m)
{ 
    // set length of x to m (i.e. degree to <= m-1)
    x.SetMaxLength(m);

    // initialize FFT representations
    long t = NextPowerOfTwo(m);
    fftRep R1(INIT_SIZE, t), R2(INIT_SIZE, t);

    zz_pX P1(INIT_SIZE, m/2);

    // -3 seems to be a reasonable choice for many p's
    t = NextPowerOfTwo(NTL_zz_pX_NEWTON_CROSSOVER)-3;

    long k = 1L << t;
    PlainInvTrunc(x, a, k);  // assumes m large enough

    long l = 0;
    while (k < m) 
    {
        l = std::min(2*k, m);
        if (l >= m) // we have l==m, 2*k was >= m
            break;

        ++t;
        TofftRep(R1, x, t);
        TofftRep(R2, a, t, 0, l-1); 
        mul(R2, R2, R1);
        FromfftRep(P1, R2, k, l-1);

        TofftRep(R2, P1, t);
        mul(R2, R2, R1);
        FromfftRep(P1, R2, 0, l-k-1);

        x.rep.SetLength(l);
        const long k_ylen = k+P1.rep.length();
        for (long i = k; i < k_ylen; ++i)
            NTL::negate(x.rep[i], P1.rep[i-k]);
        for (long i = k_ylen; i < l; ++i)
            clear(x.rep[i]);
        x.normalize();
        k = l;
    }

    // k < m <= 2k, l = m, 2^(t+1) = 2k
    // if m = 2 deg(b) = deg(a) + deg(b) + 1
    // then deg(b) <= k, deg(a) < k

    if (deg(a) < k && deg(b) < k && l == m)
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
            r1 = MulTrunc(trunc(r1, l-k), trunc(x, l-k), l-k);
        else
        {
            TofftRep(R3, r1, t+1);
            mul(R3, R3, R1);
            FromfftRep(r1, R3, 0, l-k-1); 
        }

        x.rep.SetLength(l);
        for (long i = 0; i < k; ++i) 
            x[i] = coeff(r0, i);
        for (long i = k; i < l; ++i) 
            x[i] = -coeff(r1, i-k);
        x.normalize();
    }
    else
    {
        ++t; --l;
        const long lmk = l-k;
        fftRep R3(INIT_SIZE, t);
        zz_pX P2l(INIT_SIZE, k);
        zz_pX P2h(INIT_SIZE, lmk+1);

        TofftRep(R1, x, t);

        if (deg(b) >= k)
        {
            TofftRep(R2, b, t, k, l);  
            mul(R3, R1, R2);
            FromfftRep(P2h, R3, 0, lmk); // high part
        }

        TofftRep(R2, b, t, 0, k-1);  
        mul(R3, R1, R2);
        FromfftRep(P2l, R3, 0, l); // low part

        TofftRep(R2, a, t, 0, l); 
        mul(R2, R1, R2);
        FromfftRep(P1, R2, k, l);

        /*   P1 = P1*P2l */
        TofftRep(R2, P1, t, 0, lmk); 
        TofftRep(R3, P2l, t, 0, lmk); 
        mul(R1, R2, R3);
        FromfftRep(P1, R1, 0, lmk); 

        x.rep.SetLength(l+1);
        for (long i = 0; i < k; ++i)
            x[i]=P2l[i];

        if (deg(b) >= k)
            for (long i = k; i <= l; ++i) 
                add(x[i], P2l[i], P2h[i-k]);
        else
            for (long i = k; i <= l; ++i) 
                x[i] = P2l[i];

        for (long i = k; i <= l; ++i) 
            sub(x[i], x[i], P1[i-k]);

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
    : zz_pX_shift(d_in), c(c_in)
{
    sqr(cc, c); // cc = c*c
    mul(c3, to_zz_p(3), c); // c3 = 3*c

    // in precomp, store `(x+c)^{2**k}` for all integers `k` such that `2**k`
    // is at least `4` and strictly less than `(d+1)/2` *//
    const long ell = NumBits(d);
    precomp.SetLength(std::max((long)1,ell-2));

    SetX(precomp[0]); // = x
    precomp[0][0] = c; // = x+c
    sqr(precomp[0], precomp[0]); sqr(precomp[0], precomp[0]); // = (x+c)^4

    for (long k=1; k<ell-2; ++k)
        sqr(precomp[k], precomp[k-1]);
}

/*------------------------------------------------------------*/
/* g = f(x+c)                                                 */
/* output can alias input                                     */
/*------------------------------------------------------------*/
void zz_pX_shift_DAC::shift(zz_pX& g, const zz_pX& f) const
{
    if (d == -1)
    {
        clear(g);
        return;
    }

    if (deg(f) > d)
        LogicError("Degree too large for shift");

    Vec<zz_pX> v(INIT_SIZE, ((d+1) % 4) ? (1+(d+1)/4) : ((d+1)/4));

    long i, j;
    zz_p buf,f0,f1,f2,f3;
    for (i=0, j=0; i <= deg(f)-3; i+=4, ++j)
    {
        v[j].SetLength(4);
        f0 = f[i]; f1 = f[i+1]; f2 = f[i+2]; f3 = f[i+3];
        // v[j][0] = f[i] + c * f1 + cc * (f2 + c * f3);
        mul(buf, c, f1);
        mul(v[j][0], c, f3);
        add(v[j][0], v[j][0], f2);
        mul(v[j][0], v[j][0], cc);
        add(v[j][0], v[j][0], buf);
        add(v[j][0], v[j][0], f0);
        // v[j][1] = f1 + c * (2*f2 + c3 * f3);
        mul(buf, c3, f3);
        add(buf, buf, f2);
        add(v[j][1], buf, f2);
        mul(v[j][1], v[j][1], c);
        add(v[j][1], v[j][1], f1);
        // v[j][2] = f2 + c3 * f3;
        v[j][2] = buf;
        // v[j][3] = f3;
        v[j][3] = f3;
        v[j].normalize();
    }
    if (i<=deg(f))
    {
        f0 = f[i]; f1 = coeff(f, i+1); f2 = coeff(f, i+2); f3 = coeff(f, i+3);
        v[j].SetLength(4);
        v[j][0] = f0 + c * f1 + cc * (f2 + c * f3);
        v[j][1] = f1 + c * (2*f2 + c3 * f3);
        v[j][2] = f2 + c3 * f3;
        v[j][3] = f3;
        v[j].normalize();
    }

    long idx = 0;
    while (v.length() > 1)
    {
        // TODO speed-up : precompute for repeated right-multiplication by `precomp[idx]`
        mul(v[1], v[1], precomp[idx]);
        add(v[0], v[0], v[1]);
        for (i = 1, j = 2; j+1 < v.length(); ++i, j+=2)
        {
            mul(v[i], precomp[idx], v[j+1]);
            add(v[i], v[i], v[j]);
        }
        if (j+1 == v.length())
        {
            v[i].swap(v[j]);
            ++i;
        }
        v.SetLength(i);
        ++idx;
    }
    g.swap(v[0]);
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
    : zz_pX_shift(d)
{
    // fact contains the factorials 0!, 1!, 2!, ..., d!
    fact.SetLength(d+1);
    set(fact[0]);
    for (long i = 1; i < d+1; ++i)
        mul(fact[i], i, fact[i-1]);
    // ifact is the entry-wise inverse of fact
    inv(ifact, fact);

    // v is the sum of (c**i / i!) x^i for i<=d
    v.SetLength(d+1);  
    set(v[0]);
    for (long i = 1; i < d+1; ++i)
        mul(v[i], c, v[i-1]);
    for (long i = 0; i < d+1; ++i)
        mul(v[i], v[i], ifact[i]);
}

/*------------------------------------------------------------*/
/* g = f(x+c)                                                 */
/* output can alias input                                     */
/*------------------------------------------------------------*/
void zz_pX_shift_large_characteristic::shift(zz_pX& g, const zz_pX& f) const
{
    if (d == -1)
    {
        clear(g);
        return;
    }

    zz_pX u;
    u.SetLength(d+1);
    for (long i = 0; i <= deg(f); ++i)
        mul(u[d-i], f[i], fact[i]);
    u.normalize();

    zz_pX w;
    MulTrunc(w, u, v, d+1);

    g.rep.SetLength(d+1);
    for (long i = 0; i < w.rep.length(); ++i)
        mul(g[d-i], w[i], ifact[d-i]);
    for (long i = w.rep.length(); i <= d; ++i)
        clear(g[d-i]);

    g.normalize();
}

/*------------------------------------------------------------*/
/* returns a zz_pX_shift of the right type                    */
/*------------------------------------------------------------*/
std::unique_ptr<zz_pX_shift> get_shift(long d, const zz_p& c)
{
    // build factorial d! to test whether 1,2,...,d are units in zz_p
    zz_p prod(2);
    for (long i = 3; i <= d; ++i)
        mul(prod, prod, i);
    long tmp;
    if (InvModStatus(tmp, prod.LoopHole(), zz_p::modulus()) == 0) // gcd = 1
        return std::unique_ptr<zz_pX_shift>(new zz_pX_shift_large_characteristic(d, c));
    else
        return std::unique_ptr<zz_pX_shift>(new zz_pX_shift_DAC(d, c));
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

/** Computes `c = a + (b << k)`, where the left shift means multiplication by
 * `X^k`. The integer `k` must be nonnegative. The OUT parameter `c` may alias
 * `a` but not `b`. */
void add_LeftShift(zz_pX & c, const zz_pX & a, const zz_pX & b, const long k)
{
   const long da = deg(a);
   const long db = deg(b)+k;
   const long minab = min(da, db);
   const long maxab = max(da, db);
   c.rep.SetLength(maxab+1);

   const long p = zz_p::modulus();

   long i;

   // first, low-degree range where coefficients are just copied
   if (&c != &a)
       for (i = 0; i <= min(da,k-1); ++i)
           c[i] = a[i];
   else
       i = min(da,k-1)+1;

   // if da < k-1, go up to k
   for (; i<k; ++i)
       clear(c[i]);

   // then, middle-degree range where we actually add coefficients
   for (; i <= minab; ++i)
       c[i].LoopHole() = AddMod(a[i]._zz_p__rep, b[i-k]._zz_p__rep, p);

   // finally, high-degree range where coefficients are just copied
   if (da > minab && &c != &a)
      for (; i<=maxab; ++i)
         c[i] = a[i];
   else if (db > minab)
      for (; i<=maxab; ++i)
         c[i] = b[i-k];
   else
      c.normalize();
}


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
