#include <NTL/vec_lzz_p.h>
#include <assert.h>

#include "vec_lzz_p_extra.h"
#include "structured_lzz_p.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* checks first columns                                       */
/*------------------------------------------------------------*/
void check(long p)
{
    if (p == 0)
        zz_p::FFTInit(0);
    else
        zz_p::init(p);

    for (long i = 2; i < 10; i += 1)
    {
        Vec<zz_p> dat00, dat01, dat02, dat10, dat11, dat12;

        random(dat00, 2 + i - 1);
        random(dat01, 2 + 2 - 1);
        random(dat02, 2 + i - 1);
        random(dat10, i-1 + i - 1);
        random(dat11, i-1 + 2 - 1);
        random(dat12, i-1 + i - 1);

        toeplitz_lzz_p h00(dat00, 2, i), h01(dat01, 2, 2), h02(dat02, 2, i), h10(dat10, i-1, i), h11(dat11, i-1, 2), h12(dat12, i-1, i);
        Vec<toeplitz_lzz_p> row0, row1;

        row0.SetLength(3);
        row0[0] = h00;
        row0[1] = h01;
        row0[2] = h02;
        row1.SetLength(3);
        row1[0] = h10;
        row1[1] = h11;
        row1[2] = h12;
        Vec< Vec<toeplitz_lzz_p> > H;
        H.SetLength(2);
        H[0] = row0;
        H[1] = row1;

        mosaic_toeplitz_lzz_p MH(H);

        Vec<zz_p> col;
        col = MH.first_column_of_block(0);
        for (long j = 0; j < 2; j++)
            assert (col[j] == h00(j, 0));
        for (long j = 0; j < i-1; j++)
            assert (col[j+2] == h10(j, 0));
        
        col = MH.first_column_of_block(1);
        for (long j = 0; j < 2; j++)
            assert (col[j] == h01(j, 0));
        for (long j = 0; j < i-1; j++)
            assert (col[j+2] == h11(j, 0));
        
        col = MH.first_column_of_block(2);
        for (long j = 0; j < 2; j++)
            assert (col[j] == h02(j, 0));
        for (long j = 0; j < i-1; j++)
            assert (col[j+2] == h12(j, 0));
        
    }
}

/*------------------------------------------------------------*/
/* main just calls check()                                    */
/*------------------------------------------------------------*/
int main(int argc, char** argv)
{
    check(0);
    check(786433);
    check(288230376151711813);
    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
