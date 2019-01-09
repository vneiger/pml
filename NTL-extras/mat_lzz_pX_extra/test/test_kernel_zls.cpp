#include <NTL/BasicThreadPool.h>
#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>
#include <numeric>
#include <random>

#include "util.h"
#include "mat_lzz_pX_kernel.h"

NTL_CLIENT

std::ostream &operator<<(std::ostream &out, const VecLong &s)
{
    out << "[ ";
    for (auto &i: s)
        out << i << " ";
    return out << "]";
}

void check(long m, long n, long d, bool verify){
    Mat<zz_pX> pmat;
    random(pmat, m, n, d+1);
    std::cout << "Computing kernel basis, dimensions " << m << " x " << n << ", degree " << d << std::endl;

    // declare shifts
    VecLong shift1(m,0); // uniform [0,...,0]
    VecLong shift2(m); // increasing [0,1,2,..,m-1]
    std::iota(shift2.begin(), shift2.end(),0);
    VecLong shift3(m); // decreasing [m,..,3,2,1]
    for (long i = 0; i < m; ++i)
        shift3[i] = m - i;
    VecLong shift4(m); // random shuffle of [0,1,...,m-1]
    std::iota(shift4.begin(), shift4.end(),0);
    std::shuffle(shift4.begin(), shift4.end(), std::mt19937{std::random_device{}()});

    std::vector<VecLong> shifts = {shift1, shift2, shift3, shift4};

    double t1w,t2w;

    // VIA APPROXIMATION

    for (auto shift: shifts)
    {
        { // zls-approx
            Mat<zz_pX> kerbas;
            t1w = GetWallTime();
            Mat<zz_pX> copy_pmat(pmat);
            VecLong rdeg(shift);
            VecLong pivind, pivdeg;
            kernel_basis_zls_via_approximation(kerbas, pivind, copy_pmat, rdeg);
            t2w = GetWallTime();
            cout << "time (kernel-zls-approx): " << t2w-t1w << "\t";

            if (verify)
            {
                t1w = GetWallTime();
                bool correct = is_kernel_basis(kerbas, pmat, shift, ORD_WEAK_POPOV, false);
                t2w = GetWallTime();
                std::cout << (correct ? "correct (" : "wrong (") << (t2w-t1w) << ")";
            }
            cout << endl;
        }

        { // zls-interp
            Mat<zz_pX> kerbas;
            VecLong pivind;
            t1w = GetWallTime();
            VecLong rdeg(shift);
            Mat<zz_pX> copy_pmat(pmat);
            kernel_basis_zls_via_interpolation(kerbas, pivind, copy_pmat, rdeg);
            t2w = GetWallTime();
            cout << "time (kernel-zls-interp): " << t2w-t1w << "\t";

            if (verify)
            {
                t1w = GetWallTime();
                bool correct = is_kernel_basis(kerbas, pmat, shift, REDUCED, false);
                t2w = GetWallTime();
                std::cout << (correct ? "correct (" : "wrong (") << (t2w-t1w) << ")";
            }
            cout << endl;
        }
    }
}

/*------------------------------------------------------------*/
/* main calls check                                           */
/*------------------------------------------------------------*/
int main(int argc, char ** argv)
{
    SetNumThreads(1);

    zz_p::FFTInit(0);

    if (argc!=5)
        throw std::invalid_argument("Usage: ./test_kernel_zls rdim cdim degree verify");

    check(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), (atoi(argv[4])==1));

    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
