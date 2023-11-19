#include "lzz_p_extra.h" // type_of_prime
#include "thresholds_matrix_multiply.h"
#include "mat_lzz_pX_arith.h"
#include "mat_lzz_pX_multiply.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* c = a*b                                                    */
/*------------------------------------------------------------*/
void multiply(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b, long is_prime)
{
    const long dA = deg(a);
    const long dB = deg(b);
    const long dmax = max(dA, dB);

    if (dA < 0 || dB < 0)
    {
        c.SetDims(a.NumRows(), b.NumCols());
        clear(c);
        return;
    }

    if (&c == &a || &c == &b)
    {
        Mat<zz_pX> c2;
        multiply(c2, a, b);
        c.swap(c2);
        return;
    }

    // if one of the matrices is constant, rely directly on
    // the corresponding function
    if (dA==0)
    {
        Mat<zz_p> a_cst;
        GetCoeff(a_cst, a, 0);
        mul(c, a_cst, b);
        return;
    }
    if (dB==0)
    {
        Mat<zz_p> b_cst;
        GetCoeff(b_cst, b, 0);
        mul(c, a, b_cst);
        return;
    }

    // if the left-operand has just one column, use naive multiplication
    if (a.NumCols()==1)
    {
        multiply_naive(c, a, b);
        return;
    }

    const long deg_trs = max_degree_transform();

    if (dmax <= deg_trs)
    {
        //std::cout << "transform" << std::endl;
        multiply_transform(c, a, b, dmax + 1);
        return;
    }

    // only calibrated for square matrices; here's a hack
    const long sz = (long) cbrt(a.NumRows() * a.NumCols() * b.NumCols());
    const long deg_wak = max_degree_waksman(sz);

    if (dmax <= deg_wak)
    {
        //std::cout << "waksman" << std::endl;
        multiply_waksman(c, a, b);
        return;
    }

    if (is_FFT_ready(NextPowerOfTwo(dA + dB + 1)))
    {
        //std::cout << "fft" << std::endl;
        multiply_evaluate_FFT(c, a, b);
        return;
    }

    const long p = zz_p::modulus();
    const long deg_ev = max_degree_evaluate(sz);
    if (is_prime && p > 2 * (dA + dB + 1) && dmax <= deg_ev)
    {
        //std::cout << "dense" << std::endl;
        multiply_evaluate_dense(c, a, b);
        return;
    }
    else
    {
        //std::cout << "3primes" << std::endl;
        multiply_3_primes(c, a, b);
        return;
    }
}

/*------------------------------------------------------------*/
/* right multiply by a column vector                          */
/*------------------------------------------------------------*/
void multiply(Vec<zz_pX>& c, const Mat<zz_pX>& a, const Vec<zz_pX>& b, long is_prime)
{
    Mat<zz_pX> cmat, bmat;
    bmat.SetDims(b.length(), 1);
    for (long i = 0; i < b.length(); ++i)
        bmat[i][0] = b[i];
    multiply(cmat, a, bmat, is_prime);
    c.SetLength(cmat.NumRows());
    for (long i = 0; i < cmat.NumRows(); ++i)
        c[i].swap(cmat[i][0]);
}

/*------------------------------------------------------------*/
/* left-multiply by a row vector                              */
/*------------------------------------------------------------*/
void multiply(Vec<zz_pX>& c, const Vec<zz_pX>& a, const Mat<zz_pX>& b, long is_prime)
{
    Mat<zz_pX> amat(INIT_SIZE, 1, a.length());
    amat[0] = a;
    Mat<zz_pX> cmat;
    multiply(cmat, amat, b, is_prime);
    c.swap(cmat[0]);
}


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
