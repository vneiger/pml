#include <NTL/lzz_pX.h>
#include <NTL/vector.h>

#include "lzz_pX_middle_product.h"
#include "lzz_pX_CRT.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* multipoint evaluation over zz_p                            */
/* uses subproduct tree.                                      */
/* TODO: thresholds                                           */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* constructor from points                                    */
/*------------------------------------------------------------*/
zz_pX_Multipoint_General::zz_pX_Multipoint_General(const Vec<zz_p>& q)
{
    pts = q;
    n = q.length();
    Vec<zz_pX> qpol;
    qpol.SetLength(n);
    for (long i = 0; i < n; i++)
    {
	SetCoeff(qpol[i], 0, -q[i]);
	SetCoeff(qpol[i], 1, 1);
    }
    build_subproduct_tree(tree, qpol);
    
    evaluate(cofactors, diff(tree[tree.length()-1][0]));
    for (long i = 0; i < n; i++)
    {
	cofactors[i] = 1/cofactors[i];
    }
}


/*------------------------------------------------------------*/
/* multipoint evaluation                                      */
/*------------------------------------------------------------*/
void zz_pX_Multipoint_General::evaluate(Vec<zz_p>& val, const zz_pX& f) const 
{

    long n = this->n;
    val.SetLength(n);
    long d = this->tree.length();
    Vec<zz_pX> remainders;
    remainders.SetLength(n);
    remainders[0] = f % tree[d-1][0];

    for (long i = d-2; i >= 0; i--)
    {
	long len = tree[i].length();
	if (len & 1)
	{
	    remainders[len-1] = remainders[(len-1)/2];
	}
	for (long j = len/2-1; j >= 0; j--)
	{
	    remainders[2*j+1] = remainders[j] % tree[i][2*j+1];
	    remainders[2*j] = remainders[j] % tree[i][2*j];
	}
    }
    for (long i = 0; i < n; i++)
    {
	val[i] = coeff(remainders[i], 0);
    }

}

/*------------------------------------------------------------*/
/* multipoint interpolation                                   */
/*------------------------------------------------------------*/
void zz_pX_Multipoint_General::interpolate(zz_pX& f, const Vec<zz_p>& val) 
{
    long n = this->n;
    long d = this->tree.length();
 
    Vec<zz_pX> combinations;
    combinations.SetLength(n);

    for (long i = 0; i < n; i++)
    {
	SetCoeff(combinations[i], 0, val[i]*this->cofactors[i]);
    }
    for (long i = 0; i < d-1; i++)
    {
	long len = tree[i].length();
	for (long j = 0; j < len/2; j++)
	{
	    combinations[j] = combinations[2*j]*tree[i][2*j+1] + combinations[2*j+1]*tree[i][2*j];
	}
	if (len & 1)
	{
	    combinations[(len-1)/2] = combinations[len-1];
	}
    }

    f = combinations[0];
}


/*------------------------------------------------------------*/
/* returns a zz_pX_Multipoint_General with n points           */
/*------------------------------------------------------------*/
zz_pX_Multipoint_General get_general_points(long n)
{
    Vec<zz_p> tmp;
    tmp.SetLength(n);
    for (long i = 0; i < n; i++)
    {
	tmp[i] = i;
    }
    return zz_pX_Multipoint_General(tmp);
}
