#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>
#include <vector>
#include <numeric>
#include <algorithm>
#include <random>

#include "util.h"
#include "mat_lzz_pX_extra.h"

NTL_CLIENT


/********************************************
 *  tests the approximant basis algorithms  *
 ********************************************/

std::ostream &operator<<(std::ostream &out, const VecLong &s)
{
    out << "[ ";
    for (auto &i: s)
        out << i << " ";
    return out << "]";
}

int main(int argc, char *argv[])
{
    if (argc!=5)
        throw std::invalid_argument("Usage: ./test_addshift dim deg1 deg2 shift");

    long dim = atoi(argv[1]);
    long deg1 = atoi(argv[2]);
    long deg2 = atoi(argv[3]);
    long shift = atoi(argv[4]);

    zz_p::FFTInit(0);

    double t,tt;
    long nb_iter;

    t=0.0; nb_iter=0;
    while (t<0.2)
    {
        Mat<zz_pX> a, b, c;
        random(a, dim, dim, deg1);
        random(b, dim, dim, deg2);
        tt = GetWallTime();
        add_LeftShift(c, a, b, shift);
        t += GetWallTime()-tt;
        ++nb_iter;
    }
    std::cout << "add_LeftShift:\t" << t/nb_iter << std::endl;

    t=0.0; nb_iter=0;
    while (t<0.2)
    {
        Mat<zz_pX> a, b;
        random(a, dim, dim, deg1);
        random(b, dim, dim, deg2);
        tt = GetWallTime();
        add_LeftShift(a, a, b, shift);
        t += GetWallTime()-tt;
        ++nb_iter;
    }
    std::cout << "add_LeftShift(alias):\t" << t/nb_iter << std::endl;

    t=0.0; nb_iter=0;
    while (t<0.2)
    {
        Mat<zz_pX> a, b, c;
        random(a, dim, dim, deg1);
        random(b, dim, dim, deg2);
        tt = GetWallTime();
        c = a + (b << shift);
        t += GetWallTime()-tt;
        ++nb_iter;
    }
    std::cout << "c = a + (b<<shift):\t" << t/nb_iter << std::endl;

    t=0.0; nb_iter=0;
    while (t<0.2)
    {
        Mat<zz_pX> a, b;
        random(a, dim, dim, deg1);
        random(b, dim, dim, deg2);
        tt = GetWallTime();
        a += (b << shift);
        t += GetWallTime()-tt;
        ++nb_iter;
    }
    std::cout << "a = a + (b<<shift):\t" << t/nb_iter << std::endl;

    std::cout << std::endl;
}


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
