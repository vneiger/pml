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
    SetCoeff(P, 0, val[0]);
    SetCoeff(P, 1, val[2]-val[0]-val[1]);
    SetCoeff(P, 2, val[1]);
}


