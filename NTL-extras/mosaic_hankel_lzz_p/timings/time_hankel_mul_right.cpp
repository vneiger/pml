#include <NTL/vec_lzz_p.h>
#include <iomanip>
#include <assert.h>

#include "util.h"
#include "vec_lzz_p_extra.h"
#include "mosaic_hankel_lzz_p.h"
#include "lzz_pX_middle_product.h"

NTL_CLIENT


/*------------------------------------------------------------*/
/* creates hankel matrices                                    */
/*------------------------------------------------------------*/
void check(long p)
{
    if (p == 0)
        zz_p::FFTInit(0);
    else
        zz_p::init(p);
    const double thresh = 0.01;

    cout << p << endl;
    for (long i = 190; i < 200; i += 10)
    {
        zz_pX a, b, c, d;
        Vec<zz_p> dat;
        hankel_lzz_p h;
        Mat<zz_p> M, inputM, outputM, output2M;
        double t;
        long nb;

        dat = random_vec_zz_p(i+i-1);
        h = hankel_lzz_p(dat, i, i);
        M = h.to_dense();
        a = random_zz_pX(i);
        b = random_zz_pX(i);

        cout << i << " ";

        inputM = random_mat_zz_p(i, 1);

        // poly mult
        t = get_time();
        nb = 0;
        do
        {
            c = a * b;
            nb++;
        }
        while ((get_time()-t) <= thresh);
        t = (get_time()-t) / nb;
        cout << t << " ";

        // hankel vector product
        t = get_time();
        nb = 0;
        do
        {
            h.mul_right(outputM, inputM);            
            nb++;
        }
        while ((get_time()-t) <= thresh);
        t = (get_time()-t) / nb;
        cout << t << " ";

        // create dense matrix
        t = get_time();
        nb = 0;
        do
        {
            M = h.to_dense();
            nb++;
        }
        while ((get_time()-t) <= thresh);
        t = (get_time()-t) / nb;
        cout << t << " ";

        // dense vector product
        t = get_time();
        nb = 0;
        do
        {
            output2M = M*inputM;
            nb++;
        }
        while ((get_time()-t) <= thresh);
        t = (get_time()-t) / nb;
        cout << t << " ";

        inputM = random_mat_zz_p(i, i);

        // hankel matrix product
        t = get_time();
        nb = 0;
        do
        {
            outputM = h.mul_right(inputM);            
            nb++;
        }
        while ((get_time()-t) <= thresh);
        t = (get_time()-t) / nb;
        cout << t << " ";

        // dense matrix product
        t = get_time();
        nb = 0;
        do
        {
            output2M = M*inputM;
            nb++;
        }
        while ((get_time()-t) <= thresh);
        t = (get_time()-t) / nb;
        cout << t << " ";
        
        cout << endl;
    }
}


/*------------------------------------------------------------*/
/* main just calls check()                                    */
/*------------------------------------------------------------*/
int main(int argc, char** argv)
{
    std::cout << std::fixed;
    std::cout << std::setprecision(8);
    warmup();
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
