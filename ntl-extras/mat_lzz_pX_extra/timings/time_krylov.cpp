#include <NTL/ZZ.h>
#include <NTL/lzz_pX.h>
#include <NTL/mat_lzz_p.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <NTL/BasicThreadPool.h>
#include <iomanip>
#include <ios>
#include <ostream>

#include "mat_lzz_pX_kernel.h"

NTL_CLIENT

#define BENCH_KRYLOV_ITERATES_ONE(variant)          \
void bench_krylov_iterates_one_##variant(           \
                        const Mat<zz_p> & A,        \
                        const Vec<zz_p> & v,        \
                        long d)                     \
{                                                   \
    Vec<Vec<zz_p>> B;                               \
    const long d2 = (1L << NextPowerOfTwo(d));      \
    B.SetLength(d2);                                \
    for (long k = 0; k < d2; k++)                   \
        B[k].SetLength(A.NumCols());                \
                                                    \
    double t, tt;                                   \
    long nb_iter;                                   \
                                                    \
    nb_iter = 0; t = 0.0;                           \
    while (t<0.5)                                   \
    {                                               \
        tt = GetWallTime();                         \
        krylov_iterates_one_##variant(B, A, v, d);  \
        t += GetWallTime() - tt;                    \
        ++nb_iter;                                  \
    }                                               \
                                                    \
    cout << t/nb_iter << "\t";                      \
}

#define BENCH_KRYLOV_ITERATES(variant)              \
void bench_krylov_iterates_##variant(               \
                        const Mat<zz_p> & A,        \
                        const Mat<zz_p> & V,        \
                        long d)                     \
{                                                   \
    Vec<Mat<zz_p>> B;                               \
    const long d2 = (1L << NextPowerOfTwo(d));      \
    B.SetLength(d2);                                \
    for (long k = 0; k < d2; k++)                   \
        B[k].SetDims(V.NumRows(), V.NumCols());     \
                                                    \
    double t, tt;                                   \
    long nb_iter;                                   \
                                                    \
    nb_iter = 0; t = 0.0;                           \
    while (t<0.5)                                   \
    {                                               \
        tt = GetWallTime();                         \
        krylov_iterates_##variant(B, A, V, d);      \
        t += GetWallTime() - tt;                    \
        ++nb_iter;                                  \
    }                                               \
                                                    \
    for (long k = 0; k < d2; k++)                   \
        B[k].kill();                                \
    B.kill();                                       \
                                                    \
    cout << t/nb_iter << "\t";                      \
}

/*--------------------------------*/
/* Krylov iterates, single vector */
/*--------------------------------*/

// d: integer > 0, A: matrix n x n, v: vector 1 x n
// Computes B[k] = A**k * v, 0 <= k < d
// Iterative, naive computation
void krylov_iterates_one_naive(Vec<Vec<zz_p>> & B,
                               const Mat<zz_p> & A,
                               const Vec<zz_p> & v,
                               long d)
{
    B[0] = v;
    for (long k = 1; k < d; k++)
        mul(B[k], A, B[k-1]);
}

// d: integer > 0, A: matrix n x n, v: vector 1 x n
// Computes and stores v, vA, vA**2, ...., vA**(d2-1)
// where d2 = smallest power of 2 greater than or equal to d
// Squaring, exploiting matrix multiplication
// this is not very memory efficient... matrix windows are not made public in NTL
void krylov_iterates_one_squaring(Vec<Vec<zz_p>> & B,
                                  const Mat<zz_p> & A,
                                  const Vec<zz_p> & v,
                                  long d)
{
    const long n = A.NumCols();

    // will hold A**(2**k)
    Mat<zz_p> C = A;

    B[0] = v;

    Mat<zz_p> buf;
    buf.SetDims(n, 1);
    for (long i = 0; i < n; i++)
        buf[i][0] = v[i];

    for (long k = 1; k < d; k *= 2)
    {
        // multiply by C
        mul(buf, C, buf);

        // copy to B and augment buf
        if (2*k < d)
        {
            for (long i = 0; i < n; i++)
                for (long kk = 0; kk < k; kk++)
                    B[k+kk][i] = buf[i][kk];

            buf.SetDims(n, 2*k);
            for (long i = 0; i < n; i++)
                for (long kk = 0; kk < 2*k; kk++)
                    buf[i][kk] = B[kk][i];

            // square C
            sqr(C, C);
        }
        else
        {
            for (long i = 0; i < n; i++)
                for (long kk = 0; kk < k; kk++)
                    B[k+kk][i] = buf[i][kk];
        }
    }
}

// d: integer > 0, A: matrix n x n, v: vector 1 x n
// Computes and stores v, vA, vA**2, ...., vA**(d-1)
// Fast, based on minimal kernel basis
void krylov_iterates_one_polmat(Vec<Vec<zz_p>> & B,
                                const Mat<zz_p> & A,
                                const Vec<zz_p> & v,
                                long d)
{
    const long n = A.NumRows();

    // F = [[xI - A], [-v]]
    Mat<zz_pX> F;
    F.SetDims(n+1, n);
    for (long i = 0; i < n; ++i)
    {
        for (long j = 0; j < n; ++j)
        {
            if (j == i)
                SetX(F[i][i]);
            SetCoeff(F[i][j], 0, -A[i][j]);
        }
    }
    for (long j = 0; j < n; ++j)
        SetCoeff(F[n][j], 0, -v[j]);

    Mat<zz_pX> K;
    kernel_basis_generic(K, F);
    // Note: missing rest of algo, but will be negligible for reasonable m
}

/*-----------------------------------*/
/* Krylov iterates, multiple vectors */
/*-----------------------------------*/

// B: vector of d matrices n x m
// d: integer > 0
// A: matrix n x n, V: fat vector n x m
// Computes and stores B[k] = A**k * V, 0 <= k < d
// Iterative, naive computation
void krylov_iterates_naive(Vec<Mat<zz_p>> & B,
                           const Mat<zz_p> & A,
                           const Mat<zz_p> & V,
                           long d)
{
    B[0] = V;
    for (long k = 1; k < d; k++)
        mul(B[k], A, B[k-1]);
}

// B: vector of d matrices n x m
// d: integer > 0
// A: matrix n x n, V: fat vector n x m
// Computes and stores B[k] = A**k * V, 0 <= k < d2
// where d2 = smallest power of 2 greater than or equal to d
// Squaring, exploiting matrix multiplication
void krylov_iterates_squaring(Vec<Mat<zz_p>> & B,
                              const Mat<zz_p> & A,
                              const Mat<zz_p> & V,
                              long d)
{
    const long n = V.NumRows();
    const long m = V.NumCols();

    // will hold A**(2**k)
    Mat<zz_p> C = A;

    B[0] = V;
    Mat<zz_p> buf = V;

    for (long k = 1; k < d; k *= 2)
    {
        // multiply by C
        mul(buf, C, buf);

        // copy to B and augment buf
        if (2*k < d)
        {
            for (long i = 0; i < n; i++)
                for (long kk = 0; kk < k; kk++)
                    for (long j = 0; j < m; j++)
                        B[k+kk][i][j] = buf[i][kk*m+j];

            buf.SetDims(n, 2*k*m);
            for (long i = 0; i < n; i++)
                for (long kk = 0; kk < 2*k; kk++)
                    for (long j = 0; j < m; j++)
                        buf[i][kk*m+j] = B[kk][i][j];

            // square C
            sqr(C, C);
        }
        else
        {
            for (long i = 0; i < n; i++)
                for (long kk = 0; kk < k; kk++)
                    for (long j = 0; j < m; j++)
                        B[k+kk][i][j] = buf[i][kk*m+j];
        }
    }
}

// same as squaring, but with mul only, to see timings
void krylov_iterates_squaring_sqronly(Vec<Mat<zz_p>> & B,
                                      const Mat<zz_p> & A,
                                      const Mat<zz_p> & V,
                                      long d)
{
    Mat<zz_p> C = A;
    for (long k = 1; k < d; k *= 2)
        sqr(C, C);
}

// B: vector of d matrices n x m
// d: integer > 0
// A: matrix n x n, V: fat vector n x m
// Computes and stores B[k] = A**k * V, 0 <= k < d
// Fast, based on minimal kernel basis
void krylov_iterates_polmat(Vec<Mat<zz_p>> & B,
                            const Mat<zz_p> & A,
                            const Mat<zz_p> & V,
                            long d)
{
    const long n = A.NumRows();
    const long m = V.NumCols();

    // F = [[xI - A], [-V]]
    Mat<zz_pX> F;
    F.SetDims(n+m, n);
    for (long i = 0; i < n; ++i)
    {
        for (long j = 0; j < n; ++j)
        {
            if (j == i)
                SetX(F[i][i]);
            SetCoeff(F[i][j], 0, -A[i][j]);
        }
    }
    for (long i = 0; i < m; ++i)
        for (long j = 0; j < n; ++j)
            SetCoeff(F[n+i][j], 0, -V[j][i]);

    Mat<zz_pX> K;
    kernel_basis_generic(K, F);
    // Note: missing rest of algo, but will be negligible for reasonable m
}

BENCH_KRYLOV_ITERATES_ONE(naive)
BENCH_KRYLOV_ITERATES_ONE(squaring)
BENCH_KRYLOV_ITERATES_ONE(polmat)

BENCH_KRYLOV_ITERATES(naive)
BENCH_KRYLOV_ITERATES(squaring)
BENCH_KRYLOV_ITERATES(polmat)

BENCH_KRYLOV_ITERATES(squaring_sqronly)

// n > 0
// m > 0, m <= n (?)
// t > 0: target number of total vector (typically n/m)
void bench_all(long n, long m, long t)
{
    std::cout << std::scientific;
    std::cout << std::setprecision(1);

    std::cout << n << "\t" << m << "\t" << t << "\t";
    std::cout << std::flush;

    Mat<zz_p> A = random_mat_zz_p(n, n);
    Vec<zz_p> v = random_vec_zz_p(n);
    Mat<zz_p> V = random_mat_zz_p(n, m);

    //long d_one = t;
    long d_multiple = t/m;

    // single vector
    //bench_krylov_iterates_one_naive(A, v, d_one);
    //bench_krylov_iterates_one_squaring(A, v, d_one);
    //bench_krylov_iterates_one_polmat(A, v, d_one);

    // multiple vectors
    bench_krylov_iterates_naive(A, V, d_multiple);
    std::cout << std::flush;
    bench_krylov_iterates_squaring(A, V, d_multiple);
    std::cout << std::flush;
    bench_krylov_iterates_polmat(A, V, d_multiple);
    std::cout << std::flush;
    //bench_krylov_iterates_squaring_sqronly(A, V, d_multiple);

    std::cout << std::endl;
}

void test_all(long n, long m, long t)
{
    std::cout << "testing n = " << n << ", m = " << m << ", t = " << t << std::endl;

    Mat<zz_p> A = random_mat_zz_p(n, n);
    Vec<zz_p> v = random_vec_zz_p(n);
    Mat<zz_p> V = random_mat_zz_p(n, m);

    { // one vector
        const long d = t;
        const long d2 = (1L << NextPowerOfTwo(d));

        Vec<Vec<zz_p>> B1, B2, B3;
        B1.SetLength(d);
        for (long k = 0; k < d; k++)
            B1[k].SetLength(n);
        B2.SetLength(d2);
        for (long k = 0; k < d2; k++)
            B2[k].SetLength(n);

        // single vector
        krylov_iterates_one_naive(B1, A, v, d);
        krylov_iterates_one_squaring(B2, A, v, d);
        krylov_iterates_one_polmat(B3, A, v, d);

        B2.SetLength(d); // to make sure we compare first d terms
        string res = ((B1 == B2) ? "ok" : "notok");
        std::cout << "\tnaive1 == squaring1:" << res << std::endl;
    }

    { // multiple vectors
        long d = t/m;
        const long d2 = (1L << NextPowerOfTwo(d));

        Vec<Mat<zz_p>> B1, B2, B3;
        B1.SetLength(d);
        for (long k = 0; k < d; k++)
            B1[k].SetDims(n, m);
        B2.SetLength(d2);
        for (long k = 0; k < d2; k++)
            B2[k].SetDims(n, m);

        krylov_iterates_naive(B1, A, V, d);
        krylov_iterates_squaring(B2, A, V, d);

        B2.SetLength(d); // to make sure we compare first d terms
        string res = ((B1 == B2) ? "ok" : "notok");
        std::cout << "\tnaive == squaring: " << res << std::endl;
    }

    std::cout << std::endl;
}



int main(int argc, char *argv[])
{
    if (argc != 1 && argc != 4)
        throw std::invalid_argument("usage: ./time_krylov or ./time_krylov n m t");

    SetNumThreads(1);

    const string head = "n\tm\tt\titer\tsqr\tpolmat";
    //const string head = "n\tm\tt\tnaive1\tsqr1\tpolmat1\tnaive\tsqr\tpolmat\tsqrMul";

    if (argc == 4)
    {
        const long n = atoi(argv[1]);
        const long m = atoi(argv[2]);
        const long t = atoi(argv[3]);

        //zz_p::init(97);
        //test_all(15, 2, 13);
        //test_all(23, 3, 23);
        //test_all(52, 7, 95);

        zz_p::init(131071);
        cout << endl << "prime " << zz_p::modulus() << endl;
        std::cout << head << std::endl;
        bench_all(n, m, t);

        //zz_p::init(21518131);
        //cout << endl << "25 bit prime " << zz_p::modulus() << endl;
        //std::cout << head << std::endl;
        //bench_all(n, m, t);

        //zz_p::UserFFTInit(23068673);
        //cout << endl << "25 bit FFT prime (2**21 * 11 + 1) " << zz_p::modulus() << endl;
        //std::cout << head << std::endl;
        //bench_all(n, m, t);

        //zz_p::init(288230376151711813);
        //cout << endl << "60 bit prime " << zz_p::modulus() << endl;
        //std::cout << head << std::endl;
        //bench_all(n, m, t);

        //zz_p::FFTInit(0);
        //cout << endl << "60 bit FFT prime" << endl;
        //std::cout << head << std::endl;
        //bench_all(n, m, t);
    }
    else
    {
        zz_p::init(131071);
        cout << endl << "prime " << zz_p::modulus() << endl;
        //zz_p::UserFFTInit(23068673);
        //cout << endl << "25 bit FFT prime (2**21 * 11 + 1) " << zz_p::modulus() << endl;
        std::cout << head << std::endl;

        std::vector<long> mm = {1, 8, 32};
        std::vector<long> tfacs = {1, 2, 10};

        for (long m : mm)
            for (long tfac : tfacs)
                for (long n = 7; n < 15000; n *= 2)
                    bench_all(n, m, tfac*n);
    }
}


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
