#include <NTL/lzz_pX.h>
#include <NTL/vector.h>

#include "lzz_pX_CRT.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* builds the subproduct tree naively (no Huffman)            */
/*------------------------------------------------------------*/
void build_subproduct_tree(Vec<Vec<zz_pX>> & tree, const Vec<zz_pX> & q)
{

    long len_tree = 1;
    long len = q.length();
    while (len != 1)
    {
        len_tree++;
        len = (len+1) / 2;
    }

    tree.SetLength(len_tree);
    tree[0] = q;

    long idx = 1;
    len = q.length();
    while (len != 1)
    {
        long len_new = (len+1) / 2;
        tree[idx].SetLength(len_new);
        for (long j = 0; j < len / 2; j++)
        {
            tree[idx][j] = tree[idx-1][2*j] * tree[idx-1][2*j+1];
        }
        if (len & 1)
        {
            tree[idx][len_new-1] = tree[idx-1][len-1];
        }
        len = len_new;
        idx++;
    }
}

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* CRT over zz_p                                              */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* constructor from moduli                                    */
/*------------------------------------------------------------*/
zz_pX_CRT::zz_pX_CRT(const Vec<zz_pX>& q)
{
    n = q.length();
    build_subproduct_tree(tree, q);

    multimod(cofactors, diff(tree[tree.length()-1][0]));
    for (long i = 0; i < n; i++)
    {
        cofactors[i] = MulMod(InvMod(cofactors[i], tree[0][i]), diff(tree[0][i]), tree[0][i]);
    }
}

/*------------------------------------------------------------*/
/* multiple remainders                                        */
/*------------------------------------------------------------*/
void zz_pX_CRT::multimod(Vec<zz_pX>& val, const zz_pX& f) const 
{
    long n = this->n;
    val.SetLength(n);
    long d = this->tree.length();
    val.SetLength(n);
    val[0] = f % tree[d-1][0];

    for (long i = d-2; i >= 0; i--)
    {
        long len = tree[i].length();
        if (len & 1)
        {
            val[len-1] = val[(len-1)/2];
        }
        for (long j = len/2-1; j >= 0; j--)
        {
            val[2*j+1] = val[j] % tree[i][2*j+1];
            val[2*j] = val[j] % tree[i][2*j];
        }
    }
}

/*------------------------------------------------------------*/
/* Chinese remainder                                          */
/*------------------------------------------------------------*/
void zz_pX_CRT::combine(zz_pX& f, const Vec<zz_pX>& val) const 
{
    long n = this->n;
    long d = this->tree.length();

    Vec<zz_pX> combinations;
    combinations.SetLength(n);

    for (long i = 0; i < n; i++)
    {
        combinations[i] = MulMod((val[i] % tree[0][i]), cofactors[i], tree[0][i]);
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
/* master polynomial                                          */
/*------------------------------------------------------------*/
void zz_pX_CRT::master(zz_pX& P) const
{
    P = tree[tree.length()-1][0];
}

/*------------------------------------------------------------*/
/* access to moduli                                           */
/*------------------------------------------------------------*/
void zz_pX_CRT::moduli(zz_pX& P, long i) const
{
    P = tree[0][i];
}

/*------------------------------------------------------------*/
/* number of moduli                                           */
/*------------------------------------------------------------*/
long zz_pX_CRT::length() const
{
    return n;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
