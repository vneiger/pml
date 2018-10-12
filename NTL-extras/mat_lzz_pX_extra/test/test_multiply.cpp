#include <NTL/BasicThreadPool.h>
#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>

#include "util.h"
#include "mat_lzz_pX_extra.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* checks a product (s,s+1) x (s+1,s+2) in degree < deg       */
/*------------------------------------------------------------*/
void one_check(long sz, long deg)
{
    Mat<zz_pX> a, b, c1, c2;

    random(a, sz, sz+1, deg);
    random(b, sz+1, sz+2, deg);

    // trying all possible call sequences
    // c1 = reference
    multiply_waksman(c1, a, b);

    // todo: do it only if basefield supports it
    multiply_evaluate(c2, a, b);
    if (c1 != c2)
    {
        LogicError("geometric mismatch");
    }

    // todo: do it only if basefield supports it
    multiply_evaluate_geometric(c2, a, b);
    if (c1 != c2)
    {
        LogicError("geometric mismatch");
    }

    // todo: do it only if basefield supports it
    multiply_evaluate_geometric_using_FFT(c2, a, b);
    if (c1 != c2)
    {
        LogicError("geometric mismatch");
    }

    // todo: do it only if basefield supports it
    multiply_evaluate_geometric_no_FFT(c2, a, b);
    if (c1 != c2)
    {
        LogicError("geometric mismatch");
    }

    if (is_FFT_ready(2*deg - 1))
    {
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

    // transform, if the size is reasonable
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
    std::vector<long> szs =
    {
        1, 2, 3, 5, 10, 20, 30
    };

    std::vector<long> degs =
    {
        1, 2, 3, 4, 5, 10, 15, 20, 25, 50, 60, 70, 100, 150, 200, 250, 300, 400
    };

    for (size_t si = 0; si < szs.size(); si++)
        for (size_t di = 0; di < degs.size(); di++)
            one_check(szs[si], degs[di]);
}



/*------------------------------------------------------------*/
/* checks some primes                                         */
/*------------------------------------------------------------*/
void check()
{
    zz_p::FFTInit(0);
    all_checks();
    zz_p::UserFFTInit(786433);
    all_checks();
    zz_p::init(288230376151711813);
    all_checks();
    zz_p::init(786433);
    all_checks();
}  

/*------------------------------------------------------------*/
/* main calls check                                           */
/*------------------------------------------------------------*/
int main(int argc, char ** argv)
{
    SetNumThreads(4);

    if (argc==3)
    {
        zz_p::FFTInit(0);
        one_check(atoi(argv[1]),atoi(argv[2]));
    }
    else
        check();

    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
