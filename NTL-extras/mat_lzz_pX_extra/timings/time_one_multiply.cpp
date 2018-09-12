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
    
    auto temp_c = c;
    t_eval_p = get_time();
    nb = 0;
    do
    {
        multiply_evaluate_parallel(temp_c, a, b);
        nb++;
    }
    while ((get_time()-t_eval_p) <= thres);
    t_eval_p = (get_time()-t_eval_p) / nb;

    if (temp_c != c) cout << "BAD!!!!!!!" << endl;
    /*
    t_3primes = get_time();
    nb = 0;
    do
    {
        multiply_3_primes(c, a, b);
        nb++;
    }
    while ((get_time()-t_3primes) <= thres);
    t_3primes = (get_time()-t_3primes) / nb;
                
    t_waksman = get_time();
    nb = 0;
    do
    {
        multiply_waksman(c, a, b);
        nb++;
    }
            while ((get_time()-t_waksman) <= thres);
    t_waksman = (get_time()-t_waksman) / nb;
    

    t_check = get_time();
    nb = 0;
    do
    {
        multiply(c, a, b);
        nb++;
    }
    while ((get_time()-t_check) <= thres);
    t_check = (get_time()-t_check) / nb;
    */
    cout << p << " " << sz << " " << deg << endl;
    cout << t_eval << " " << t_eval_p << endl;
    //cout << t_eval << " " << t_eval_p << " " << t_3primes << " " << t_waksman << " " << t_check << endl;
}




/*------------------------------------------------------------*/
/* main calls check                                           */
/*------------------------------------------------------------*/
int main(int argc, char ** argv)
{
    std::cout << std::fixed;
    std::cout << std::setprecision(8);

    long sz = 1000;
    long deg = 30;

    if (argc > 1)
    {
        sz = atoi(argv[1]);
        deg = atoi(argv[2]);
    }
    
    SetNumThreads(4);

    warmup();
    check(23068673, sz, deg);
    // check(23068673, sz, deg);
    // check(288230376151711813, sz, deg);
    return 0;
}
