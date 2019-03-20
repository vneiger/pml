#include <NTL/BasicThreadPool.h>
#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>

#include "util.h"
#include "mat_lzz_pX_extra.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* checks a product (sz,sz+1) x (sz+1,sz+2) in degree < deg   */
/*------------------------------------------------------------*/
void one_check(long sz, long deg)
{
    Mat<zz_pX> a, b, c1, c2;

    random(a, sz, sz+1, deg);
    random(b, sz+1, sz+2, deg);

    // c1 = reference
    multiply_waksman(c1, a, b);

    // todo: do it only if basefield supports it
    multiply_evaluate_geometric(c2, a, b);
    if (c1 != c2)
    {
        LogicError("geometric mismatch");
    }

    if (is_FFT_ready(NextPowerOfTwo(2*deg - 1)))
    {
        multiply_evaluate_FFT_direct_ll_type(c2, a, b);
        if (c1 != c2)
        {
            LogicError("FFT direct_ll_type mismatch");
        }

        multiply_evaluate_FFT_direct(c2, a, b);
        if (c1 != c2)
        {
            LogicError("FFT direct mismatch");
        }

        multiply_evaluate_FFT_matmul1(c2, a, b);
        if (c1 != c2)
        {
            LogicError("FFT matmul1 mismatch");
        }

        multiply_evaluate_FFT_matmul1_new(c2, a, b);
        if (c1 != c2)
        {
            LogicError("FFT matmul1_new mismatch");
        }

        multiply_evaluate_FFT_matmul2(c2, a, b);
        if (c1 != c2)
        {
            LogicError("FFT matmul2 mismatch");
        }

        multiply_evaluate_FFT_matmul3(c2, a, b);
        if (c1 != c2)
        {
            LogicError("FFT matmul3 mismatch");
        }

        multiply_evaluate_FFT(c2, a, b);
        if (c1 != c2)
        {
            LogicError("FFT mismatch");
        }
    }

    multiply_3_primes(c2, a, b);
    if (c1 != c2)
    {
        LogicError("3 primes mismatch");
    }

    multiply_evaluate_dense(c2, a, b);
    if (c1 != c2)
    {
        LogicError("dense mismatch");
    }

    multiply_evaluate_dense2(c2, a, b);
    if (c1 != c2)
    {
        LogicError("dense2 mismatch");
    }

    long do_transform = (deg <= 10) || ((sz <= 400) && (deg <= 10)) || ((sz <= 50) && (deg <= 20));
    if (do_transform)
    {
        multiply_transform(c2, a, b);
        if (c1 != c2)
        {
            LogicError("transform mismatch");
        }
    }

    multiply(c2, a, b);
    if (c1 != c2)
    {
        LogicError("multiply mismatch");
    }    
}


/*------------------------------------------------------------*/
/* for a give prime, checks some (size, degree)               */
/*------------------------------------------------------------*/
void all_checks()
{
    VecLong szs =
    {
        1, 2, 3, 5, 10, 20, 30, 40, 50, 70, 100
    };

    VecLong degs =
    {
        1, 2, 3, 4, 5, 10, 15, 20, 25, 50, 60, 70, 100, 150, 200, 250, 300
    };

    for (size_t si = 0; si < szs.size(); si++)
        for (size_t di = 0; di < degs.size(); di++)
        {
            std::cout << szs[si] << "\t" << degs[di] << "...\t";
            std::cout << std::flush;
            one_check(szs[si], degs[di]);
            std::cout << "ok." << std::endl;
        }
}


/*------------------------------------------------------------*/
/* main calls check                                           */
/*------------------------------------------------------------*/
int main(int argc, char ** argv)
{
    zz_p::FFTInit(0);
    all_checks();
    zz_p::UserFFTInit(786433);
    all_checks();
    zz_p::init(288230376151711813);
    all_checks();
    zz_p::init(786433);
    all_checks();

    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
