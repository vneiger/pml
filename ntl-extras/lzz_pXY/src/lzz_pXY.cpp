#include <NTL/lzz_pX.h>

#include "lzz_pX_CRT.h"
#include "structured_lzz_p.h"
#include "mat_lzz_pX_extra.h"
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

void mul(zz_pXY& c, const zz_pXY& a, const zz_p b)
{
    if (a.is_zero() || b == 0)
    {
        c.zero();
        return;
    }

    c.rep.SetLength(a.degY() + 1, to_zz_pX(0));

    for (long i = 0; i <= a.degY(); i++)
        mul(c.rep[i], a.rep[i], b);
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


/*------------------------------------------------------------*/
/* computes the resultant of f and g                          */
/* returns 0 iff the computation failed                       */
/*------------------------------------------------------------*/
long resultant(zz_pX& res, const zz_pXY& f, const zz_pXY& g)
{
    long dyf = f.degY();
    long dyg = g.degY();
    long n = max(dyf, dyg);

    long dxf = f.degX();
    long dxg = g.degX();
    long d = max(dxf, dxg);

    long xtra = 2;
    long dr = 2*n*d + xtra;

    Vec<zz_p> val_res, val_lcf, val_lcg;
    Vec<zz_pX> val_f, val_g;

    val_res.SetLength(dr);
    val_f.SetLength(dr);
    val_g.SetLength(dr);

    zz_p a0, b0;
    long OK = 0;
    long nb = 0;
    zz_pX_Multipoint_Geometric evG;
    do
    {
        nb++;
        OK = 1;
        a0 = random_zz_p();
        b0 = random_zz_p();
        evG = zz_pX_Multipoint_Geometric(a0, b0, max(dxf, dxg)+1);
        evG.evaluate(val_lcf, f.rep[dyf], dr);
        evG.evaluate(val_lcg, g.rep[dyg], dr);
        for (long i = 0; i < dr; i++)
            if (val_lcf[i] == 0 || val_lcg[i] == 0)
            {
                OK = 0;
                break;
            }
    }
    while (OK == 0 && nb < 5);
  
    if (nb == 5)
        return 0;

    for (long i = 0; i < dr; i++)
    {
        val_f[i].rep.SetLength(dyf + 1);
        val_f[i][dyf] = val_lcf[i];
        val_g[i].rep.SetLength(dyg + 1);
        val_g[i][dyg] = val_lcg[i];
    }

    Vec<zz_p> tmp;
    for (long i = 0; i < dyf; i++)
    {
        evG.evaluate(tmp, f.rep[i], dr);
        for (long j = 0; j < dr; j++)
            val_f[j][i] = tmp[j];
    }
    for (long i = 0; i < dyg; i++)
    {
        evG.evaluate(tmp, g.rep[i], dr);
        for (long j = 0; j < dr; j++)
            val_g[j][i] = tmp[j];
    }

    for (long i = 0; i < dr; i++)
        resultant(val_res[i], val_f[i], val_g[i]);

    zz_pX_Multipoint_Geometric evG_long(a0, b0, dr);
    evG_long.interpolate(res, val_res);
    return 1;
}

/*------------------------------------------------------------*/
/* computes the resultant of f and g, Villard's algorithm     */
/* returns 0 iff the computation failed                       */
/*------------------------------------------------------------*/
long resultant_villard(zz_pX& res, const zz_pXY& f, const zz_pXY& g)
{
    long nf = f.degY();
    long ng = g.degY();
    long n = max(nf, ng);

    long df = f.degX();
    long dg = g.degX();
    long d = max(df, dg);

    long xtra = 2;
    long m = (long) cbrt(n);
    long delta = 2 * 2 * d * ((n + m - 1)/m) + xtra;


    Vec<zz_p> val_res, val_lcf, val_lcg;
    Vec<zz_pX> val_f, val_g;

    val_f.SetLength(delta);
    val_g.SetLength(delta);

    zz_p a0, b0;
    long OK = 0;
    long nb = 0;
    zz_pX_Multipoint_Geometric evG;
    do
    {
        nb++;
        OK = 1;
        a0 = random_zz_p();
        // TODO: add randomness by using b0 as well
        evG = zz_pX_Multipoint_Geometric(a0, max(df, dg)+1);
        evG.evaluate(val_lcf, f.rep[nf], delta);
        evG.evaluate(val_lcg, g.rep[ng], delta);
        for (long i = 0; i < delta; i++)
            if (val_lcf[i] == 0 || val_lcg[i] == 0)
            {
                OK = 0;
                break;
            }
    }
    while (OK == 0 && nb < 5);
  
    if (nb == 5)
        return 0;

    for (long i = 0; i < delta; i++)
    {
        val_f[i].rep.SetLength(nf + 1);
        val_f[i][nf] = val_lcf[i];
        val_g[i].rep.SetLength(ng + 1);
        val_g[i][ng] = val_lcg[i];
    }

    Vec<zz_p> tmp;
    for (long i = 0; i < nf; i++)
    {
        evG.evaluate(tmp, f.rep[i], delta);
        for (long j = 0; j < delta; j++)
            val_f[j][i] = tmp[j];
    }
    for (long i = 0; i < ng; i++)
    {
        evG.evaluate(tmp, g.rep[i], delta);
        for (long j = 0; j < delta; j++)
            val_g[j][i] = tmp[j];
    }

    // blocks have size mxm
    Vec<Mat<zz_p>> blocks;
    blocks.SetLength(delta);
    for (long i = 0; i < delta; i++)
    {
        sylvester_lzz_p S(val_f[i], val_g[i]);
        long r = S.top_right_block_inverse(blocks[i], m);
        if (r == 0)
            return 0;
    }

    Mat<zz_pX> basis;

    // set up pts
    zz_pX x;
    SetCoeff(x, 1, 1);
    Vec<zz_p> pts;
    evG.evaluate(pts, x, delta); // just gets powers of r

    if (0)
        matrix_recon_interpolation(basis, pts, blocks);
    else
        matrix_recon_interpolation_geometric(basis, pts, a0, blocks);

    determinant_via_linsolve(res, basis);
    return 1;
}


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
