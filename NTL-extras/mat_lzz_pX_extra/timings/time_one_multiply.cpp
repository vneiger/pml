#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>
#include <limits.h>
#include <vector>

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

    const double thres = 0.001;

    long nb;
    Mat<zz_pX> a, b, c;
    double t_eval,t_eval_p, t_3primes, t_waksman, t_check;

    random_mat_zz_pX(a, sz, sz, deg);
    random_mat_zz_pX(b, sz, sz, deg);

    t_eval = get_time();
    nb = 0;
    do
    {
        multiply_evaluate(c, a, b);
        nb++;
    }
    while ((get_time()-t_eval) <= thres);
    t_eval = (get_time()-t_eval) / nb;
    
    cout << p << " " << sz << " " << deg << endl;
    cout << t_eval << << endl;
}




/*------------------------------------------------------------*/
/* main calls check                                           */
/*------------------------------------------------------------*/
int main(int argc, char ** argv)
{
    std::cout << std::fixed;
    std::cout << std::setprecision(8);

    long s[] = {4, 16, 32};
    long d[] = {2000,10000,20000};
    
    std::vector<long> sz(s,s+sizeof(s)/sizeof(long));
    std::vector<long> deg(d, d+sizeof(d)/sizeof(long));

    warmup();
    for (long t = 1; t <= 4; t *= 2)
    {
        SetNumThreads(t);
        std::cout << "-----------------------------\n";
        std::cout << "NUM THREADS: " << t << "\n";
        for (auto i = sz.begin(); i < sz.end(); i++)
        {
            for (auto j = deg.begin(); j < deg.end(); j++)
            {
                check(0, *i, *j);
                cout << "\n";
                check(23068673, *i, *j);
            }
        }
        std::cout << "\n";
    }
    return 0;
}
// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
