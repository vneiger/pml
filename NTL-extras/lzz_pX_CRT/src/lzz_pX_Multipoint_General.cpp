#include <NTL/lzz_pX.h>
#include <NTL/vector.h>

#include "lzz_pX_extra.h"
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
/* going up the subproduct tree                               */
/*------------------------------------------------------------*/
static void up_tree(zz_pX& f, const Vec<Vec<zz_pX>> & tree, const Vec<zz_p> & coeffs)
{
    long n = coeffs.length();
    long d = tree.length();
    Vec<zz_pX> combinations;
    combinations.SetLength(n);

    for (long i = 0; i < n; i++)
    {
        SetCoeff(combinations[i], 0, coeffs[i]);
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

    reverse_root = reverse(tree[tree.length()-1][0], n);
    inverse_root = InvTrunc(reverse_root, n);
}

/*------------------------------------------------------------*/
/* constructor from points[offset]...points[offset+length-1]  */
/*------------------------------------------------------------*/
zz_pX_Multipoint_General::zz_pX_Multipoint_General(const Vec<zz_p>& q, long offset, long length)
{
    n = length;
    pts.SetLength(n);
    for (long i = 0; i < n; ++i)
        pts[i] = q[offset+i];

    Vec<zz_pX> qpol;
    qpol.SetLength(n);
    for (long i = 0; i < n; i++)
    {
        SetCoeff(qpol[i], 0, -pts[i]);
        SetCoeff(qpol[i], 1, 1);
    }
    build_subproduct_tree(tree, qpol);

    evaluate(cofactors, diff(tree[tree.length()-1][0]));
    for (long i = 0; i < n; i++)
    {
        cofactors[i] = 1/cofactors[i];
    }

    reverse_root = reverse(tree[tree.length()-1][0], n);
    inverse_root = InvTrunc(reverse_root, n);
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
    Vec<zz_p> coeffs;
    coeffs.SetLength(n);

    for (long i = 0; i < n; i++)
    {
        coeffs[i] = val[i] * cofactors[i];
    }
    up_tree(f, tree, coeffs);
}


/*------------------------------------------------------------*/
/* transposed multipoint evaluation                           */
/* parameter output_size not used for the moment              */
/*------------------------------------------------------------*/
void zz_pX_Multipoint_General::t_evaluate(zz_pX& f, const Vec<zz_p>& val, long output_size) const 
{
    zz_pX rev_num, num;
    up_tree(rev_num, tree, val);
    num = reverse(rev_num, n-1);
    f = MulTrunc(num, inverse_root, n);
}

/*------------------------------------------------------------*/
/* transposed multipoint interpolation                        */
/*------------------------------------------------------------*/
void zz_pX_Multipoint_General::t_interpolate(Vec<zz_p>& val, const zz_pX& f)
{
    zz_pX num = MulTrunc(f, reverse_root, n);
    zz_pX rev_num = reverse(num, n-1);
    evaluate(val, rev_num);

    for (long i = 0; i < n; i++)
    {
        val[i] = val[i] * cofactors[i];
    }
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

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
