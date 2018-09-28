#include "vec_lzz_p_extra.h"
#include "lzz_pX_extra.h"

NTL_CLIENT


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
