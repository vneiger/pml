#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <numeric>
#include <random>
#include <NTL/BasicThreadPool.h>
#include <cmath>

//#define SAFETY_CHECKS

#include "util.h"
#include "mat_lzz_pX_extra.h"
#include "mat_lzz_pX_sequence.h"
#include "sage_output.h"

int main(int argc, char *argv[])
{
    if (argc != 3)
        throw std::invalid_argument("usage: ./test_matrix_recon n m");

    long n = atoi(argv[1]);
    long m = atoi(argv[2]);
    long d = ceil((double)n/m);
    zz_p::FFTInit(0);

    Mat<zz_pX> F;
    random(F, m, m, d);

    zz_p r;
    random(r);
    zz_pX_Multipoint_Geometric eval(r,d);
    Vec<zz_p> pts;
    zz_pX x;
    SetCoeff(x,1,1);
    eval.evaluate(pts,x);

    zz_pX poly{zz_p{1}};
    for (auto &p : pts)
        poly = poly * (x - zz_pX{p});

    Vec<Mat<zz_p>> mat_eval;
    mat_eval.SetLength(d);
    for (auto &mat: mat_eval)
        mat.SetDims(m,m);

    for (long i = 0; i < m; i++)
    {
        for (long j = 0; j < m; j++)
        {
            Vec<zz_p> ev;
            eval.evaluate(ev, F[i][j]);
            for (long t = 0; t < ev.length(); t++)
                mat_eval[t][i][j] = ev[t];
        }
    }
    Mat<zz_pX> P;
    matrix_recon_interpolation(P,pts,mat_eval);
    cout << "P: " << P << endl;

    Mat<zz_pX> prod;
    multiply(prod,P,F);
    for (long i = 0; i < m; i++)
        for (long j = 0; j < m; j++)
            prod[i][j] = prod[i][j] % poly;

    cout << "-prod: " << -prod << endl;

    Mat<zz_pX> F_id = F;
    F_id.SetDims(2*m,m);
    for (long i = 0; i < m; i++)
        F_id[i+m][i] = zz_pX{1};
    cout << "F: " << F << endl;
    cout << "F_id: " << F_id << endl;
    
    std::vector<long> shift(2*m,0);
    Mat<zz_pX> intbas;
    
    Vec<Mat<zz_p>> mat_eval2;
    mat_eval2.SetLength(d);
    for (auto &mat: mat_eval2)
        mat.SetDims(2*m,m);

    for (long i = 0; i < 2*m; i++)
    {
        for (long j = 0; j < m; j++)
        {
            Vec<zz_p> ev;
            eval.evaluate(ev, F_id[i][j]);
            for (long t = 0; t < ev.length(); t++)
                mat_eval2[t][i][j] = ev[t];
        }
    }
    pmbasis(intbas, mat_eval2, pts, shift);
    multiply(prod, intbas, F_id);

    cout << "intbas: " << intbas << endl;
    for (long i = 0; i < prod.NumRows(); i++)
        for (long j = 0; j < prod.NumCols(); j++)
            prod[i][j] = prod[i][j] % poly;

    cout << "prod2: " << prod << endl;
    cout << "deg: " << degree_matrix(intbas) << endl;
    cout << "len: " << d << endl;
}










