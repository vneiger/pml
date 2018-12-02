#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <NTL/ZZ.h>
#include <NTL/tools.h>
#include <NTL/BasicThreadPool.h>

#include <iomanip>
#include <limits.h>

#include "util.h"
#include "lzz_p_extra.h"
#include "mat_lzz_pX_extra.h"


NTL_CLIENT

/*------------------------------------------------------------*/
/* checks some products                                       */
/*------------------------------------------------------------*/
void check(long p, long sz, long sz2, long sz3, long deg)
{
    if (p == 0)
        zz_p::FFTInit(0);
    else
        zz_p::init(p);

    const double thresh = 0.001;
    
    long nb;
    Mat<zz_pX> a, b, c;
    double t_eval, t_eval_dense, t_eval_direct, t_3primes, t_waksman, t_check;

    random(a, sz, sz2, deg);
    random(b, sz2, sz3, deg);

    cout << p << " " << sz << " " << sz2 << " " << sz3 << " " << deg << endl;
    
    t_3primes = get_time();
    nb = 0;
    do
    {
        multiply_3_primes(c, a, b);
        nb++;
    }
    while ((get_time()-t_3primes) <= thresh);
    t_3primes = (get_time()-t_3primes) / nb;
                
    t_eval = get_time();
    nb = 0;
    do
    {
        multiply_evaluate(c, a, b);
        nb++;
    }
    while ((get_time()-t_eval) <= thresh);
    t_eval = (get_time()-t_eval) / nb;

    t_waksman = get_time();
    nb = 0;
    do
    {
        multiply_waksman(c, a, b);
        nb++;
    }
            while ((get_time()-t_waksman) <= thresh);
    t_waksman = (get_time()-t_waksman) / nb;

    
    if (deg < 100)
    {
        t_eval_dense = get_time();
        nb = 0;
        do
        {
            multiply_evaluate_dense(c, a, b);
            nb++;
        }
        while ((get_time()-t_eval_dense) <= thresh);
        t_eval_dense = (get_time()-t_eval_dense) / nb;
    }

    if (p == 0)
    {
        t_eval_direct = get_time();
        nb = 0;
        do
        {
            multiply_evaluate_FFT_direct(c, a, b);
            nb++;
        }
        while ((get_time()-t_eval_direct) <= thresh);
        t_eval_direct = (get_time()-t_eval_direct) / nb;
        cout << t_eval_direct << "  check " << (c == a*b) << endl;

        t_eval_direct = get_time();
        nb = 0;
        do
        {
            multiply_evaluate_FFT_matmul1(c, a, b);
            nb++;
        }
        while ((get_time()-t_eval_direct) <= thresh);
        t_eval_direct = (get_time()-t_eval_direct) / nb;
        cout << t_eval_direct << "  check " << (c == a*b) << endl;

        t_eval_direct = get_time();
        nb = 0;
        do
        {
            multiply_evaluate_FFT(c, a, b);
            nb++;
        }
        while ((get_time()-t_eval_direct) <= thresh);
        t_eval_direct = (get_time()-t_eval_direct) / nb;
        cout << t_eval_direct << "  check " << (c == a*b) << endl;
    }

    t_check = get_time();
    nb = 0;
    do
    {
        multiply(c, a, b);
        nb++;
    }
    while ((get_time()-t_check) <= thresh);
    t_check = (get_time()-t_check) / nb;

    cout << t_eval << " " << t_3primes << " " << t_waksman << " " << t_check << "   ";
    if (deg < 10)
        cout << " " << t_eval_dense;
    cout << endl;
}




/*------------------------------------------------------------*/
/* main calls check                                           */
/*------------------------------------------------------------*/
int main(int argc, char ** argv)
{
    SetNumThreads(1);

    std::cout << std::fixed;
    std::cout << std::setprecision(8);

    long sz1 = 10;
    long sz2 = 10;
    long sz3 = 10;
    long deg = 10;

    if (argc ==  5)
    {
        sz1 = atoi(argv[1]);
        sz2 = atoi(argv[2]);
        sz3 = atoi(argv[3]);
        deg = atoi(argv[4]);
    }

    warmup();
    check(0, sz1, sz2, sz3, deg);
    check(23068673, sz1, sz2, sz3, deg);
    check(288230376151711813, sz1, sz2, sz3, deg);
    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
