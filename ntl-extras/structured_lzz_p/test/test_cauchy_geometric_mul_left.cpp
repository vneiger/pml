#include <NTL/vec_lzz_p.h>
#include <assert.h>

#include "vec_lzz_p_extra.h"
#include "structured_lzz_p.h"

NTL_CLIENT


/*------------------------------------------------------------*/
/* does some left vector / matrix multiplications             */
/*------------------------------------------------------------*/
void check(long p)
{
    if (p == 0)
        zz_p::FFTInit(0);
    else
        zz_p::init(p);

    for (long i = 1; i < 100; i += 1)
    {
        zz_p a = random_zz_p();
        
        cauchy_geometric_lzz_p C;
        Mat<zz_p> M;
        Vec<zz_p> input, output, check;
        Mat<zz_p> inputM, outputM, checkM;
        
        long j = i+5;

        C = cauchy_geometric_lzz_p(to_zz_p(1), power(a, i), a, i, j);
        M = C.to_dense();

        input = random_vec_zz_p(i);
        output = C.mul_left(input);
        check = input * M;
        assert (check == output);

        inputM = random_mat_zz_p(3, i);
        outputM = C.mul_left(inputM);
        checkM = inputM * M;
        assert (checkM == outputM);
        
        C = cauchy_geometric_lzz_p(to_zz_p(1), power(a, j), a, j, i);
        M = C.to_dense();

        input = random_vec_zz_p(j);
        output = C.mul_left(input);
        check = input * M;
        assert (check == output);

        inputM = random_mat_zz_p(3, j);
        outputM = C.mul_left(inputM);
        checkM = inputM * M;
        assert (checkM == outputM);
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
