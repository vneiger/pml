#include <iomanip>
#include <random>

#include "util.h"
#include "mat_lzz_pX_extra.h"

NTL_CLIENT

static std::ostream &operator<<(std::ostream &out, const VecLong &s)
{
    out << "[ ";
    for (auto &i: s)
        out << i << " ";
    return out << "]";
}

int main(int argc, char *argv[])
{   
    if (argc!=6)
        throw std::invalid_argument("Usage: ./test_division nbits rdim cdim degA avgrdegB");

    const long nbits = atoi(argv[1]);
    const long rdim = atoi(argv[2]);
    const long cdim = atoi(argv[3]);
    const long degA = atoi(argv[4]);
    const long avgdegB = atoi(argv[5]);

    if (nbits==0)
        zz_p::FFTInit(0);
    else
        zz_p::init(NTL::GenPrime_long(nbits));

    // build some somehow random non-uniform rdeg
    VecLong rdegB(rdim);
    for (long i = 0; i < rdim/3; ++i)
        rdegB[i] = std::max((long)1,avgdegB/3);
    for (long i = rdim/3; i < 2*rdim/3; ++i)
        rdegB[i] = std::max((long)1,avgdegB);
    for (long i = 2*rdim/3; i < rdim; ++i)
        rdegB[i] = std::max((long)1,2*avgdegB - avgdegB/3);
    std::shuffle(rdegB.begin(), rdegB.end(), std::mt19937{std::random_device{}()});

    // copies for later verification
    Mat<zz_pX> A, B, Q, R;

    // compute division
    double t=0.0,tt;
    long nb_iter=0;
    while (t<0.2)
    {
        // random A
        random(A,rdim,cdim,degA+1);
        // random B of given row degree
        do
            random_mat_zz_pX_rdeg(B,rdim,rdim,rdegB);
        while (not is_row_reduced(B));

        tt = GetWallTime();
        quo_rem(Q,R,A,B);
        t += GetWallTime() - tt;
        ++nb_iter;
    }
    std::cout << "Time(division): " << t/nb_iter;

    // test
    bool correct = true;
    Mat<zz_pX> T;
    multiply(T,B,Q);
    add(T, T, R);
    if (T != A)
        correct = false;
    VecLong rdegR;
    row_degree(rdegR, R);
    row_degree(rdegB, B); // just in case it actually differs from the targetted one
    for (long i = 0; i < rdim; ++i)
        if (rdegR[i] >= rdegB[i])
            correct = false;

    std::cout << ", " << (correct?"correct":"false") << std::endl;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
