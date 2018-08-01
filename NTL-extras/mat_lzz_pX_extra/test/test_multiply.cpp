#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>

#include "magma_output.h"
#include "mat_lzz_pX_extra.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* checks some products                                       */
/*------------------------------------------------------------*/
void one_check(long sz, long deg, long p)
{
    Mat<zz_pX> a, b, c0, c1, c2, c4, c5;
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
    long do_naive = ((sz <= 400)  && (deg <= 40))
                        || ((sz <= 50) && (deg <= 200))
                        ||  ((sz <= 10) && (deg <= 2000))
                        || (sz==2);
    if (do_naive)
    {
        t = GetTime();
        multiply_waksman(c1, a, b);
        cout << GetTime()-t << ",";
    }
    else 
    { 
        cout << "-1,";
    }

    // evaluation -- should be done only if feasible
    t = GetTime();
    multiply_evaluate(c2, a, b);
    cout << GetTime()-t << ",";
    
    if (do_naive && (c1 != c2))
    {
        cout << "(geometric mismatch) ";
    }

    // 3 primes FFT
    t = GetTime();
    multiply_3_primes(c4, a, b);
    cout << GetTime()-t << ",";
    if (c4 != c2)
    {
        cout << "(3 primes mismatch) ";
    }

    // transform, if the size is reasonable
    long do_transform = (deg <= 10) || ((sz <= 400) && (deg <= 10)) || ((sz <= 50) && (deg <= 20));
    if (do_transform)
    {
        t = GetTime();
        multiply_transform(c5, a, b);
        cout << GetTime()-t << " ";
        if (c5 != c2)
        {
            cout << "(transform mismatch) ";
        }
    }
    else
    {
        cout << "-1";
    }

    cout << endl;
}

/*------------------------------------------------------------*/
/* checks some products in small degree                       */
/*------------------------------------------------------------*/
void one_check_smalldeg(long sz, long deg, long p)
{
    Mat<zz_pX> a, b, c1;
    double t;
    
    if (p == 0) // init zz_p with FFTInit()
    {
        zz_p::FFTInit(0);
    }
    else
    {
        zz_p::init(p);
    }

    random_mat_zz_pX(a, sz, sz, deg);
    random_mat_zz_pX(b, sz, sz, deg);
    cout << "size=" << sz << ", length=" << deg << " ";
    
    if (sz*sz*deg > 10000000000)
    {
        throw std::invalid_argument("Calling one_check_smalldeg with large size and/or degree");
    }
    
    // naive algorithms: if the size is reasonable
    bool do_naive = (sz*sz*sz*deg*deg <= 1000000000);
    if (do_naive)
    {
        t = GetTime();
        multiply_transform_naive(c1, a, b);
        cout << GetTime()-t << " ";
    }
    
    Mat<zz_pX> c2;
    switch (deg) 
    {
    case 2: // TODO: write is_karatsuba.
        t = GetTime();
        multiply_transform_karatsuba(c2, a, b);
        cout << GetTime()-t << " ";
        break;
    case 3: // TODO: write is_montgomery
        t = GetTime();
        multiply_transform_montgomery3(c2, a, b);
        cout << GetTime()-t << " ";
        break;
    case 4: // TODO: write is_karatsuba4
        t = GetTime();
        multiply_transform_karatsuba4(c2, a, b);
        cout << GetTime()-t << " ";
        break;
    }
    
    if (do_naive && c2 != c1)
    {
        cout << "(transform mismatch) ";
    }

    cout << endl;
}

/*------------------------------------------------------------*/
/* checks some products                                       */
/*------------------------------------------------------------*/
void check(long sz=200, long deg=4)
{
    long p1 = 288230376151711813;
    // long p0 = 0;
    // long p2 = 23068673;

    cout << "                   waksman    eval       3primes    transf\n";
    
    for (long i = 1; i < deg; i++)
    {
        one_check(sz, i, p1);
    }
}  

/*------------------------------------------------------------*/
/* main calls check                                           */
/*------------------------------------------------------------*/
int main(int argc, char ** argv)
{
    std::cout << std::fixed;
    std::cout << std::setprecision(8);

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
