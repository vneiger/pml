#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>

#include "util.h"
#include "magma_output.h"
#include "mat_lzz_pX_extra.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* checks some products                                       */
/*------------------------------------------------------------*/
void one_check(long sz, long deg, long p)
{
    Mat<zz_pX> a, b, c0, c1, c2, c4, c5;
    long nb;
    double t;

    if (p == 0) // init zz_p with FFTInit()
    {
        zz_p::FFTInit(100);
    }
    else
    {
        zz_p::init(p);
    }

    cout << p<< "," << sz << "," << deg << ",";

    random_mat_zz_pX(a, sz, sz, deg);
    random_mat_zz_pX(b, sz, sz, deg);

    // naive algorithm, if the size is reasonable
    long do_naive = ((sz <= 200)  && (deg <= 40))
                        || ((sz <= 50) && (deg <= 200))
                        ||  ((sz <= 10) && (deg <= 2000))
                        || (sz==2);
    if (do_naive)
    {
        t = GetWallTime();
        nb = 0;
        do
        {
            multiply_waksman(c1, a, b);
            nb++;
        }
        while ((GetWallTime()-t) <= 0.001);
        
        t = (GetWallTime()-t) / nb;
        cout << t << ",";
    }
    else 
    { 
        cout << "999999,";
    }

    // evaluation -- should be done only if feasible
    t = GetWallTime();
    nb = 0;
    do
    {
        multiply_evaluate(c2, a, b);
        nb++;
    }
    while ((GetWallTime()-t) <= 0.001);
    
    t = (GetWallTime()-t) / nb;
    cout << t << ",";
    
    if (do_naive && (c1 != c2))
    {
        cout << "(geometric mismatch) ";
    }

    // 3 primes FFT
    t = GetWallTime();
    nb = 0;
    do
    {
        multiply_3_primes(c4, a, b);
        nb++;
    }
    while ((GetWallTime()-t) <= 0.001);

    t = (GetWallTime()-t) / nb;
    cout << t << ",";

    if (c4 != c2)
    {
        cout << "(3 primes mismatch) ";
    }

    // transform, if the size is reasonable
    long do_transform = (deg <= 10) || ((sz <= 400) && (deg <= 10)) || ((sz <= 50) && (deg <= 20));
    if (do_transform)
    {
        t = GetWallTime();
        nb = 0;
        do
        {
            multiply_transform(c5, a, b);
            nb++;
        }
        while ((GetWallTime()-t) <= 0.001);
        
        t = (GetWallTime()-t) / nb;
        cout << t << ",";

        if (c5 != c2)
        {
            cout << "(transform mismatch) ";
        }
    }
    else
    {
        cout << "999999";
    }

    cout << endl;
}

/*------------------------------------------------------------*/
/* checks some products                                       */
/*------------------------------------------------------------*/
void check(long sz=200, long deg=4)
{
    long p0 = 0;
    long p1 = 23068673;
    long p2 = 288230376151711813;

    one_check(sz, deg, p0);
    one_check(sz, deg, p1);
    one_check(sz, deg, p2);
}  

/*------------------------------------------------------------*/
/* main calls check                                           */
/*------------------------------------------------------------*/
int main(int argc, char ** argv)
{
    std::cout << std::fixed;
    std::cout << std::setprecision(8);

    warmup();
    if (argc==1)
    {
        check();
    }
    else if (argc==3)
    {
        check(atoi(argv[1]), atoi(argv[2]));
    }
    else
    {
        throw std::invalid_argument("Usage: ./test_multiply OR ./test_multiply size degree");
    }

    return 0;
}
