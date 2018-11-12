#include <NTL/lzz_pX.h>

#include "lzz_pXY.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* a helper class for bivariate polynomials over zz_p         */
/* a zz_pXY is simply a vector of zz_pX                       */
/* with the convention f = sum_i rep[i](X) Y^i              */ 
/* minimal functionalities are provided                       */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* builds from a vector of zz_pX                              */
/*------------------------------------------------------------*/
zz_pXY::zz_pXY(const Vec<zz_pX>&  c)
{
    rep = c;
    normalize();
}


/*------------------------------------------------------------*/
/* builds from a zz_pX                                        */
/*------------------------------------------------------------*/
zz_pXY::zz_pXY(const zz_pX&  c)
{
    rep.SetLength(1);
    rep[0] = c;
    normalize();
}

/*------------------------------------------------------------*/
/* zero test                                                  */
/*------------------------------------------------------------*/
long zz_pXY::is_zero() const
{
    return rep.length() == 0;
}

/*------------------------------------------------------------*/
/* set to zero                                                */
/*------------------------------------------------------------*/
void zz_pXY::zero() 
{
    rep.SetLength(0);
}

/*------------------------------------------------------------*/
/* total degree                                               */
/*------------------------------------------------------------*/
long zz_pXY::tdeg() const 
{
    long res = -1;
    for (long i = 0; i < rep.length(); i++)
        res = max(res, deg(rep[i]) + i);
    return res;
}

/*------------------------------------------------------------*/
/* degree in Y                                                */
/*------------------------------------------------------------*/
long zz_pXY::degY() const 
{
    return rep.length()-1;
}

/*------------------------------------------------------------*/
/* degree in X                                                */
/*------------------------------------------------------------*/
long zz_pXY::degX() const 
{
    long d = -1;
    for (long i = 0; i < rep.length(); i++)
        d = max(d, deg(rep[i]));
    return d;
}

/*------------------------------------------------------------*/
/* resizes the array of coefficients                          */
/* to remove the trailing entries that are zero, if any       */
/*------------------------------------------------------------*/
void zz_pXY::normalize()
{
    long idx = rep.length()-1;
    while (idx >= 0 && rep[idx] == 0)
        idx--;
    rep.SetLength(idx+1);
}

/*------------------------------------------------------------*/
/* I/O using NTL's representation (without commas!)           */
/*------------------------------------------------------------*/
ostream& operator<<(ostream& s, const zz_pXY& a)
{
    return s << a.rep;
}

/*------------------------------------------------------------*/
/* random element f with deg(f,x) < dx and deg(f,y) < dy      */
/*------------------------------------------------------------*/
void random_zz_pXY(zz_pXY&  f, long dx, long dy)
{
    f.rep.SetLength(dy);
    for (long i = 0; i < dy; i++)
        f.rep[i] = random_zz_pX(dx);
}

/*------------------------------------------------------------*/
/* evaluation with respect to X at a point                    */
/*------------------------------------------------------------*/
void evaluate(zz_pX&  value, const zz_pXY&  f, const zz_p&  point)
{
    value = to_zz_pX(0);
    for (long i = 0; i <= f.degY(); i++)
        SetCoeff(value, i, eval(f.rep[i], point));
    value.normalize();
}

/*------------------------------------------------------------*/
/* addition                                                   */
/*------------------------------------------------------------*/
void add(zz_pXY& c, const zz_pXY& a, const zz_pXY& b)
{
    long dyA = a.degY();
    long dyB = b.degY();
    long dy = max(dyA, dyB);
    long dy0 = min(dyA, dyB);
    c.rep.SetLength(dy + 1);
    for (long i = 0; i <= dy0; i++)
    {
        c.rep[i] = a.rep[i] + b.rep[i];
    }
    for (long i = dy0 + 1; i <= dyB; i++)
    {
        c.rep[i] = b.rep[i];
    }
    for (long i = dy0 + 1; i <= dyA; i++)
    {
        c.rep[i] = a.rep[i];
    }
    c.normalize();
}

/*------------------------------------------------------------*/
/* subtraction                                                */
/*------------------------------------------------------------*/
void sub(zz_pXY& c, const zz_pXY& a, const zz_pXY& b)
{
    long dyA = a.degY();
    long dyB = b.degY();
    long dy = max(dyA, dyB);
    long dy0 = min(dyA, dyB);
    c.rep.SetLength(dy + 1);
    for (long i = 0; i <= dy0; i++)
    {
        c.rep[i] = a.rep[i] - b.rep[i];
    }
    for (long i = dy0 + 1; i <= dyB; i++)
    {
        c.rep[i] = -b.rep[i];
    }
    for (long i = dy0 + 1; i <= dyA; i++)
    {
        c.rep[i] = a.rep[i];
    }
    c.normalize();
}

/*------------------------------------------------------------*/
/* naive multiplication                                       */
/*------------------------------------------------------------*/
void mul_naive(zz_pXY& c, const zz_pXY& a, const zz_pXY& b)
{
    if (a.is_zero() || b.is_zero())
    {
        c.zero();
        return;
    }

    zz_pXY c_tmp;
    long da = a.degY();
    long db = b.degY();

    c_tmp.rep.SetLength(da + db + 1, to_zz_pX(0));
    for (long i = 0; i <= da; i++)
        for (long j = 0; j <= db; j++)
            c_tmp.rep[i+j] += a.rep[i] * b.rep[j];
    c_tmp.normalize();
    c = c_tmp;
}


/*------------------------------------------------------------*/
/* to kronecker substitution                                  */
/*------------------------------------------------------------*/
void to_kronecker(zz_pX& out, const Vec<zz_pX>& a, long degX)
{
    if (a.length() == 0)
    {
        out = 0;
        return;
    }

    long dyA = a.length() - 1;
    out.rep.SetLength((degX+1) * (dyA+1));
    zz_p * out_e = out.rep.elts();
    for (long i = 0; i <= dyA; i++)
    {
        long d = deg(a[i]);
        if (d >= 0)
        {
            const zz_p * ae = a[i].rep.elts();
            for (long j = 0; j <= d; j++)
                out_e[j] = ae[j];
            for (long j = d+1; j <= degX; j++)
                out_e[j] = 0;
        }
        else
            for (long j = 0; j <= degX; j++)
                out_e[j] = 0;
        out_e += degX + 1;
    }
    out.normalize();
}

/*------------------------------------------------------------*/
/* from kronecker substitution                                */
/*------------------------------------------------------------*/
void from_kronecker(Vec<zz_pX>& c, const zz_pX& a, long degX)
{
    if (a == 0)
    {
        c.SetLength(0);
        return;
    }

    long degY = deg(a) / (degX + 1);
    long idx = 0;
    c.SetLength(degY + 1);
    const zz_p * ae = a.rep.elts();
    for (long i = 0; i <= degY ; i++)
    {
        c[i].rep.SetLength(degX + 1);
        zz_p * ce = c[i].rep.elts();
        long j;
        for (j = 0; j <= degX && idx <= deg(a); j++)
            ce[j] = ae[idx++];
        for (; j <= degX; j++)
            ce[j] = 0;
        c[i].normalize();
    }
}

/*------------------------------------------------------------*/
/* to kronecker substitution                                  */
/*------------------------------------------------------------*/
void to_kronecker(zz_pX& out, const zz_pXY& a, long degX)
{
    if (a.is_zero())
    {
        out = 0;
        return;
    }
    to_kronecker(out, a.rep, degX);
}

/*------------------------------------------------------------*/
/* from kronecker substitution                                */
/*------------------------------------------------------------*/
void from_kronecker(zz_pXY& c, const zz_pX& a, long degX)
{
    c.zero();
    from_kronecker(c.rep, a, degX);
}

/*------------------------------------------------------------*/
/* kronecker multiplication                                   */
/*------------------------------------------------------------*/
void mul_kronecker(zz_pXY& c, const zz_pXY& a, const zz_pXY& b)
{
    if (&c == &a || &c == &b)
    {
        c = mul_kronecker(a, b);
        return;
    }

    if (a.is_zero() || b.is_zero())
    {
        c.zero();
        return;
    }
    
    long dxA = a.degX();
    long dxB = b.degX();
    long dx = dxA + dxB;
    zz_pX apX, bpX, cpX;

    to_kronecker(apX, a, dx);
    to_kronecker(bpX, b, dx);
    
    cpX = apX * bpX;

    from_kronecker(c, cpX, dx);
}


/*------------------------------------------------------------*/
/* shift in x.                                                */ 
/* i >=0 means multiply by x^i                                */
/* i <=0 means div by x^(-i)                                  */
/*------------------------------------------------------------*/
void shift_x(zz_pXY& b, const zz_pXY& a, long i)
{
    long dy = a.degY();
    b.rep.SetLength(dy + 1);
    if (i >= 0)
        for (long j = 0; j <= dy; j++)
            b.rep[j] = a.rep[j] << i;
    else
        for (long j = 0; j <= dy; j++)
            b.rep[j] = a.rep[j] >> (-i);
    b.normalize();
}


/*------------------------------------------------------------*/
/* shift in y.                                                */ 
/* i >=0 means multiply by y^i                                */
/* i <=0 means div by y^(-i)                                  */
/*------------------------------------------------------------*/
void shift_y(zz_pXY& b, const zz_pXY& a, long i)
{
    if (&b == &a)
    {
        b = shift_y(a, i);
    }

    long dy = a.degY();
    if (i >= 0)
    {
        b.rep.SetLength(dy + 1 + i);
        for (long j = 0; j < i; j++)
            b.rep[j] = 0;
        for (long j = 0; j <= dy; j++)
            b.rep[j + i] = a.rep[j];
    }
    else
    {
        b.rep.SetLength(dy + 1 + i);
        for (long j = 0; j <= dy + i; j++)
            b.rep[j] = a.rep[j - i];
    }
    b.normalize();
}


/*------------------------------------------------------------*/
/* trunc in x (modulo x^i)                                    */ 
/*------------------------------------------------------------*/
void trunc_x(zz_pXY& b, const zz_pXY &a, long i)
{
    long dy = a.degY();
    b.rep.SetLength(dy + 1);
    for (long j = 0; j <= dy; j++)
        b.rep[j] = trunc(a.rep[j], i);
    b.normalize();

}

/*------------------------------------------------------------*/
/* trunc in y (modulo y^i)                                    */ 
/*------------------------------------------------------------*/
void trunc_y(zz_pXY& b, const zz_pXY &a, long i)
{
    long dy = a.degY();
    if (i > dy)
    {
        b = a;
        return;
    }

    b.rep.SetLength(i);
    for (long j = 0; j < i; j++)
        b.rep[j] = a.rep[j];
    b.normalize();
}

/*------------------------------------------------------------*/
/* setting coefficients                                       */
/*------------------------------------------------------------*/
void SetCoeff(zz_pXY& x, long i, const zz_pX& c)
{
   long j, m;

   if (i < 0) 
      LogicError("SetCoeff: negative index");

   m = x.degY();

   if (i > m && IsZero(c)) 
       return; 

   if (i > m) 
   {
      x.rep.SetLength(i + 1);
      for (j = m + 1; j < i; j++)
         clear(x.rep[j]);
   }
   x.rep[i] = c;
   x.normalize();
}


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
