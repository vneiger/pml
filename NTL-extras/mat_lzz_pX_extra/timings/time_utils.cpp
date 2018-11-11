#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>
#include <vector>

#include "util.h"
#include "mat_lzz_pX_extra.h"

#define TIME(a)         \
    tt = get_time();    \
    a                   \
    tt = get_time()-tt; \
    cout << tt << "\t";



NTL_CLIENT

/******************************************************/
/* Times the degree functions for polynomial matrices */
/******************************************************/

using namespace std;

void run_bench(long m, long n, long d)
{
    double tt;

    TIME(
         Mat<zz_pX> pmat;
         random(pmat, m, n, d);
        )

    TIME(
         deg(pmat);
        )

    TIME(
         std::vector<long> degs(pmat.NumRows());
         row_degree(degs,pmat);
        )

    TIME(
         degs.resize(pmat.NumCols());
         column_degree(degs,pmat);
        )

    TIME(
         Mat<long> deg_mat;
         degree_matrix(deg_mat, pmat);
        )

    TIME(
         Mat<zz_p> lead_mat;
         leading_matrix(lead_mat, pmat);
        )

    TIME(
         std::vector<long> pivind(pmat.NumRows());
         std::vector<long> pivdeg(pmat.NumRows());
         pivot_index(pivind, pivdeg, pmat, std::vector<long>(), true);
        )

    TIME(
         pivind.resize(pmat.NumCols());
         pivdeg.resize(pmat.NumCols());
         pivot_index(pivind, pivdeg, pmat, std::vector<long>(), false);
        )

    TIME(
         is_reduced(pmat);
        )

    TIME(
         is_weak_popov(pmat);
        )

    TIME(
         is_popov(pmat);
        )

//    std::vector<long> rs {0,2,1,3};
//    std::vector<long> cs {4,2,0};
//
//    degree_matrix(deg_mat,pmat,rs,true);
//    cout << "degree matrix with row wise shift: " << endl << deg_mat<<endl;
//    degree_matrix(deg_mat,pmat,cs,false);
//    cout << "degree matrix with col wise shift: " << endl << deg_mat<<endl;
//
//    degs.resize(pmat.NumRows());
//    row_degree(degs,pmat,rs);
//
//    degs.resize(pmat.NumCols());
//    column_degree(degs,pmat,cs);
//
//    leading_matrix(lead_mat,pmat,rs,true);
//
//    leading_matrix(lead_mat,pmat,cs,false);
//
//    pivind.resize(pmat.NumRows());
//    pivdeg.resize(pmat.NumRows());
//    pivot_index(pivind, pivdeg, pmat, rs, true);
//
//    pivind.resize(pmat.NumCols());
//    pivdeg.resize(pmat.NumCols());
//    pivot_index(pivind, pivdeg, pmat, cs, false);

    //is_weak_popov(pmat,rs,true,true); // ordered weak

//    cout << endl << "Test Left shifts" << endl;
//    cout << "left shift operator: " << (pmat << d) << endl;
//    //cout << "left shift mutator: " << (pmat <<= d) << endl;
//    cout << "test procedure: " << LeftShift(pmat,d) << endl;
//    LeftShift(pmat,pmat,d);
//    //cout << "row test procedure: " << LeftShiftRow(pmat,0,1) << endl;
//    //LeftShiftRow(pmat,pmat,0,1);
//    cout << "row test mutator: " << pmat << endl;
//    //col
//    cout << "col test procedure: " << LeftShiftCol(pmat,1,2) << endl;
//    LeftShiftCol(pmat,pmat,1,2);
//    cout << "col test mutator: " << pmat << endl;
//
//    cout << endl << "Test Right shifts" << endl;
//    cout << "left shift operator: " << (pmat >> 2) << endl;
//    cout << "left shift mutator: " << (pmat >>= 2) << endl;
//    cout << "test procedure: " << RightShift(pmat,2) << endl;
//    RightShift(pmat,pmat,2);
//    cout << "test mutator: " << pmat << endl;
//    //row
//    cout << "row test procedure: " << RightShiftRow(pmat,0,1) << endl;
//    RightShiftRow(pmat,pmat,0,1);
//    cout << "row test mutator: " << pmat << endl;
//    //col
//    cout << "col test procedure: " << RightShiftCol(pmat,1,2) << endl;
//    RightShiftCol(pmat,pmat,1,2);
//    cout << "col test mutator: " << pmat << endl;

    tt = get_time();
    trunc(pmat,5*d);
    tt = get_time()-tt;
    cout << tt << "\t";

    tt = get_time();
    trunc(pmat,10);
    tt = get_time()-tt;
    cout << tt << "\t";

    //cout << endl << "Test trunc row" << endl;
    //cout << "trunc: " << truncRow(pmat,0,2) << endl;
    //truncRow(pmat,pmat,0,2);
    //cout << "mutator trunc: " << pmat << endl;

    //cout << endl << "Test trunc col" << endl;
    //cout << "trunc: " << truncCol(pmat,1,1) << endl;
    //truncCol(pmat,pmat,1,1);
    //cout << "mutator trunc: " << pmat << endl;
}

int main(int argc, char *argv[])
{
    std::cout << std::fixed;
    std::cout << std::setprecision(5);

    std::vector<long> rdims = {5,10,20,40,70,100,150,200,400,1000,};
    std::vector<long> cdims = {5,10,20,40,70,100,150,200,400,1000,};
    std::vector<long> degs = {5,10,20,40,70,100,150,200,400,1000,2000,4000,8000,16000,};
    std::vector<long> nbits = {20,31,42,60};
    std::vector<long> fftprimes = {786433,2013265921,2748779069441,1139410705724735489,}; // 20, 31, 42, 60 bits

    warmup();

    for (size_t b = 0; b < nbits.size(); ++b)
    {
        std::cout << "Bench utils, nbits = " << nbits[b];
        long p = NTL::GenPrime_long(nbits[b]);
        std::cout << " (normal prime p = " << p << ", fft prime p = " << fftprimes[b] << ")" << std::endl;

        std::cout << "info\trdim\tcdim\tdeg\trand\tdeg\trdeg\tcdeg\tdegmat\tlmat\trpiv\tcpiv\tisred\tisweak\tispop\ttruncL\ttruncS\t" << std::endl;

        for (size_t i = 0; i < rdims.size(); ++i)
            for (size_t j = 0; j < cdims.size(); ++j)
                for (size_t k = 0; k < degs.size(); ++k)
                {
                    std::cout << "\t" << rdims[i] << "\t" << cdims[j] << "\t" << degs[k] << "\t";
                    zz_p::init(p);
                    run_bench(rdims[i], cdims[j], degs[k]);
                    std::cout << std::endl << "fft\t" << rdims[i] << "\t" << cdims[j] << "\t" << degs[k] << "\t";
                    zz_p::UserFFTInit(fftprimes[b]);
                    run_bench(rdims[i], cdims[j], degs[k]);
                    std::cout << std::endl;
                }
    }

    std::cout << std::endl;

    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
