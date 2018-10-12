#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>
#include <limits.h>
#include <NTL/BasicThreadPool.h>

#include "util.h"
#include "lzz_p_extra.h"
#include "mat_lzz_pX_extra.h"


NTL_CLIENT

/*------------------------------------------------------------*/
/* checks some products                                       */
/*------------------------------------------------------------*/
void check(long p, long sz, long deg)
{
    if (p == 0)
        zz_p::FFTInit(0);
    else
        zz_p::init(p);

    const double thresh = 0.001;
    
    long nb;
    Mat<zz_pX> a, b, c;
    double t_eval, t_3primes, t_waksman, t_check;

    random(a, sz, sz, deg);
    random(b, sz, sz, deg);
    
    t_eval = get_time();
    nb = 0;
    do
    {
        multiply_evaluate(c, a, b);
        nb++;
    }
    while ((get_time()-t_eval) <= thresh);
    t_eval = (get_time()-t_eval) / nb;
    

    t_3primes = get_time();
    nb = 0;
    do
    {
        multiply_3_primes(c, a, b);
        nb++;
    }
    while ((get_time()-t_3primes) <= thresh);
    t_3primes = (get_time()-t_3primes) / nb;
                
    t_waksman = get_time();
    nb = 0;
    do
    {
        multiply_waksman(c, a, b);
        nb++;
    }
            while ((get_time()-t_waksman) <= thresh);
    t_waksman = (get_time()-t_waksman) / nb;
    

    t_check = get_time();
    nb = 0;
    do
    {
        multiply(c, a, b);
        nb++;
    }
    while ((get_time()-t_check) <= thresh);
    t_check = (get_time()-t_check) / nb;

    cout << p << " " << sz << " " << deg << endl;
    cout << t_eval << " " << t_3primes << " " << t_waksman << " " << t_check << endl;

}




/*------------------------------------------------------------*/
/* main calls check                                           */
/*------------------------------------------------------------*/
int main(int argc, char ** argv)
{
    SetNumThreads(1);

    std::cout << std::fixed;
    std::cout << std::setprecision(8);

    long sz = 10;
    long deg = 10;

    if (argc > 1)
    {
        sz = atoi(argv[1]);
        deg = atoi(argv[2]);
    }

    warmup();
    check(0, sz, deg);
    check(23068673, sz, deg);
    check(288230376151711813, sz, deg);
    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
