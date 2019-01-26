#include <NTL/vec_lzz_p.h>
#include <assert.h>

#include "vec_lzz_p_extra.h"
#include "structured_lzz_p.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* tests a few products                                       */
/*------------------------------------------------------------*/
void check(long p)
{
    if (p == 0)
        zz_p::FFTInit(0);
    else
        zz_p::init(p);

    for (long i = 2; i < 30; i += 1)
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

        mosaic_toeplitz_lzz_p MH;
        MH = mosaic_toeplitz_lzz_p(H);

        Mat<zz_p> Mdense = MH.to_dense();

        Vec<zz_p> input, output;
        input = random_vec_zz_p(MH.NumCols());
        output = random_vec_zz_p(MH.NumRows());
        output = MH.mul_right(input);
        Vec<zz_p> output2 = Mdense * input;
        assert(output2 == output);

        Mat<zz_p> inputM, outputM;
        inputM = random_mat_zz_p(MH.NumCols(), 3);
        outputM = random_mat_zz_p(MH.NumRows(), 19);
        outputM = MH.mul_right(inputM);
        Mat<zz_p> output2M = Mdense * inputM;
        assert(output2M == outputM);
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
