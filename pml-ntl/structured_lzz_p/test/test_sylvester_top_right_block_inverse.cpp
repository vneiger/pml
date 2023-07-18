#include <NTL/vec_lzz_p.h>
#include <assert.h>

#include "vec_lzz_p_extra.h"
#include "structured_lzz_p.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* creates sylvester matrices                                 */
/* computes top-right block of their inverses and checks      */
/*------------------------------------------------------------*/
void check(long p)
{
    if (p == 0)
        zz_p::FFTInit(0);
    else
        zz_p::init(p);

    for (long i = 2; i < 200; i += 30)
    {
        zz_pX F, G;
        do 
        {
            do
                F = random_zz_pX(i);
            while (deg(F) < 1);
            
            do
                G = random_zz_pX(i);
            while (deg(G) < 1);
        }
        while (deg(GCD(F, G)) != 0);
        
        sylvester_lzz_p S(F, G);
        Mat<zz_p> block;
        long m;

        m = 3;
        if (S.NumRows() >= m)
        {
            long r = S.top_right_block_inverse(block, m);
            assert(r == 1);
            Mat<zz_p> Sinv = inv(S.to_dense());
            Mat<zz_p> subSinv;
            subSinv.SetDims(m, m);
            for (long i = 0; i < m; i++)
                for (long j = 0; j < m; j++)
                    subSinv[i][j] = Sinv[i][j + S.NumCols() - m];
            assert(block == subSinv);
        }

        m = S.NumRows()-3;
        if (m > 0)
        {
            long r = S.top_right_block_inverse(block, m);
            assert(r == 1);
            Mat<zz_p> Sinv = inv(S.to_dense());
            Mat<zz_p> subSinv;
            subSinv.SetDims(m, m);
            for (long i = 0; i < m; i++)
                for (long j = 0; j < m; j++)
                    subSinv[i][j] = Sinv[i][j + S.NumCols() - m];
            assert(block == subSinv);
        }
    }
}

/*------------------------------------------------------------*/
/* main just calls check()                                    */
/*------------------------------------------------------------*/
int main(int argc, char** argv)
{
    check(288230376151711813);
    check(0);
    check(786433);
    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
