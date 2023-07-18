#include <NTL/lzz_pX.h>
#include <NTL/vector.h>
#include <NTL/matrix.h>

#include "lzz_pX_CRT.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* Karatsuba for 2 terms                                      */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* constructor                                                */
/*------------------------------------------------------------*/
zz_pX_Transform_karatsuba::zz_pX_Transform_karatsuba()
{
    m = 2;
    n = 3;
}

/*------------------------------------------------------------*/
/* (a0,a1) -> (a0, a0, a0+a1)                                 */
/*------------------------------------------------------------*/
void zz_pX_Transform_karatsuba::forward_left(Vec<zz_p>& val, const zz_pX& P) const
{
    val.SetLength(3);
    val[0] = coeff(P, 0);
    val[1] = coeff(P, 1);
    val[2] = val[0] + val[1];
}

/*------------------------------------------------------------*/
/* (b0,b1) -> (b0, b0, b0+b1)                                 */
/*------------------------------------------------------------*/
void zz_pX_Transform_karatsuba::forward_right(Vec<zz_p>& val, const zz_pX& P) const
{
    val.SetLength(3);
    val[0] = coeff(P, 0);
    val[1] = coeff(P, 1);
    val[2] = val[0] + val[1];
}

/*------------------------------------------------------------*/
/* reconstructs C=AB from a0 b0, a1 b1, (a0+a1)(b0+b1)        */
/*------------------------------------------------------------*/
void zz_pX_Transform_karatsuba::backward(zz_pX& P, const Vec<zz_p>& val) const
{
    P = 0;
    SetCoeff(P, 0, val[0]);
    SetCoeff(P, 1, val[2]-val[0]-val[1]);
    SetCoeff(P, 2, val[1]);
}


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* Montgomery for 3 input terms                               */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* constructor                                                */
/*------------------------------------------------------------*/
zz_pX_Transform_montgomery3::zz_pX_Transform_montgomery3()
{
    m = 3;
    n = 6;
}

/*------------------------------------------------------------*/
/* cf. Montgomery Section 2, C=X                              */
/*------------------------------------------------------------*/
void zz_pX_Transform_montgomery3::forward_left(Vec<zz_p>& val, const zz_pX& P) const
{
    val.SetLength(6);
    const zz_p a0 = coeff(P, 0);    
    const zz_p a1 = coeff(P, 1);
    const zz_p a2 = coeff(P, 2);
    val[0] = a0;
    val[1] = a1;
    val[2] = a2;
    val[3] = a0+a2;
    zz_p s = a1+a2;
    val[4] = s;
    val[5] = a0+s;
}

/*------------------------------------------------------------*/
/* cf. Montgomery Section 2, C=X                              */
/*------------------------------------------------------------*/
void zz_pX_Transform_montgomery3::backward(zz_pX& P, const Vec<zz_p>& val) const
{
    P = 0;
    const zz_p s = val[2]-val[3];
    SetCoeff(P, 0, val[0]);
    SetCoeff(P, 1, s-val[4]+val[5]);
    SetCoeff(P, 2, val[1]-val[0]-s);
    SetCoeff(P, 3, val[4]-val[1]-val[2]);
    SetCoeff(P, 4, val[2]);
}

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* Karatsuba for 4 input terms                                */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* constructor                                                */
/*------------------------------------------------------------*/
zz_pX_Transform_karatsuba4::zz_pX_Transform_karatsuba4()
{
    m = 4;
    n = 9;
}

/*------------------------------------------------------------*/
/* 2 rounds of karatsuba                                      */
/*------------------------------------------------------------*/
void zz_pX_Transform_karatsuba4::forward_left(Vec<zz_p>& val, const zz_pX& P) const
{
    val.SetLength(9);
    const zz_p a0 = coeff(P, 0);    
    const zz_p a1 = coeff(P, 1);
    const zz_p a2 = coeff(P, 2);
    const zz_p a3 = coeff(P, 3);

    val[0] = a0;
    val[1] = a1;
    val[2] = a0 + a1;
    val[3] = a2;
    val[4] = a3;
    val[5] = a2 + a3;
    const zz_p s1 = a0 + a2;
    const zz_p s2 = a1 + a3;
    val[6] = s1;
    val[7] = s2;
    val[8] = s1 + s2;
}

/*------------------------------------------------------------*/
/* 2 rounds of karatsuba                                      */
/*------------------------------------------------------------*/
void zz_pX_Transform_karatsuba4::backward(zz_pX& P, const Vec<zz_p>& val) const
{
    P = 0;
    P.rep.SetLength(7);
    zz_p * coeffs = P.rep.elts();
    coeffs[0] = val[0];
    coeffs[1] = val[2] - val[0] - val[1];
    const zz_p s = val[1] - val[3];
    coeffs[2] = s + val[6] - val[0];
    coeffs[4] = val[7] - val[4] - s;
    coeffs[5] = val[5] - val[3] - val[4];
    coeffs[6] = val[4];
    coeffs[3] = val[8] - (coeffs[1] + coeffs[5] + val[6] + val[7]);
    P.normalize();
}


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
