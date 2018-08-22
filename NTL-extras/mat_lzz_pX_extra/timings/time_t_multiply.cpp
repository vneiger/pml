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
    Mat<zz_pX> a, b1, b2, c;

    if (p == 0) // init zz_p with FFTInit()
    {
        zz_p::FFTInit(0);
    }
    else
    {
        zz_p::init(p);
    }

    cout << p<< ", " << sz << "," << deg << ", ";

    random_mat_zz_pX(a, sz, sz+1, deg);
    random_mat_zz_pX(c, sz+1, sz+2, 2*deg - 1);

    middle_product_evaluate_geometric(b1, a, c, deg-1, deg-1);
    multiply_evaluate_geometric(b2, a, c);
    b2 >>= (deg-1);
    trunc(b2, b2, deg);
    cout << "OK::" << (b1 == b2) << " ";

    double t;
    long nb;

    t = get_time();
    nb = 0;
    do
    {
        middle_product_evaluate_geometric(b1, a, c, deg-1, deg-1);
        nb++;
    }
    while ((get_time()-t) <= 0.001);
    t = (get_time()-t) / nb;
    cout << t << " ";

    t = get_time();
    nb = 0;
    do
    {
        multiply_evaluate_geometric(b2, a, c);
        nb++;
    }
    while ((get_time()-t) <= 0.001);
    t = (get_time()-t) / nb;
    cout << t << " ";
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
    one_check(50, 100, 288230376151711813);
    std::cout << std::fixed;
    std::cout << std::setprecision(8);

    // if (argc==1)
    // {
    //     check();
    // }
    // else if (argc==3)
    // {
    //     check(atoi(argv[1]), atoi(argv[2]));
    // }
    // else
    // {
    //     throw std::invalid_argument("Usage: ./test_multiply OR ./test_multiply size degree");
    // }

    return 0;
}
