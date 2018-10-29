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
#include "mat_lzz_pX_approximant.h"

NTL_CLIENT

std::ostream &operator<<(std::ostream &out, const std::vector<long> &s)
{
    out << "[ ";
    for (auto &i: s)
        out << i << " ";
    return out << "]";
}

typedef Vec<Vec<Vec<zz_p>>> Coeffs;

// slices p into len slices of size sz
void slice(Vec<zz_pX> &res, const zz_pX &p, const long len, const long sz){
    res.SetLength(len);
    long at = 0;
    for (long i = 0; i < len; i++)
    {
        res[i] = zz_pX();
        for (long j = 0; j < sz; j++)
            SetCoeff(res[i], j, coeff(p, at++));
    }
}

// generates (t*a^(s*i) mod g) for 0 <= i < l
void gen_pows (Vec<zz_pX> &pow, 
        const zz_pX &t,
        const zz_pX &a,
        const long s,
        const zz_pX &g,
        const long l)
{
    zz_pXModulus g_mod;
    build(g_mod, g);

    zz_pX a_pow;
    PowerMod(a_pow, a, s, g_mod); // a_pow = a^s

    zz_pXMultiplier a_pow_mul;
    build(a_pow_mul, a_pow, g_mod);

    pow.SetLength(l);
    pow[0] = t;
    for (long i = 1; i < l; i++)
    {
        MulMod(pow[i], pow[i-1], a_pow_mul, g_mod);
    }

}

void get_quos (Mat<zz_pX> &quos,
        const Vec<zz_pX> &alphas,
        const Vec<zz_pX> &As,
        const zz_pX &g,
        const long m)
{
    const long n = deg(g);
    const long pad = 2*m - ((n-2*m-1)%(2*m));

    const long len = ceil((n*1.0)/(2*m));
    const zz_pX S = InvTrunc(reverse(g,n), n-1);

    Mat<zz_pX> mat1,mat2,mat3;
    mat1.SetDims(alphas.length(), len);
    mat2.SetDims(alphas.length(), len);
    mat3.SetDims(alphas.length(), len);

    for (long i = 0; i < alphas.length(); i++)
    {
        zz_pX alpha_rev = reverse(alphas[i], n-1);
        zz_pX A_rev = reverse(As[i], n-1);
        zz_pX A_st = MulTrunc(S, A_rev, n-1);
        LeftShift(A_st,A_st, pad);


        Vec<zz_pX> s1,s2;
        slice(s1, alpha_rev, len, 2*m);
        slice(s2, A_st, len, 2*m);

        // populate mat1, mat2
        for (long j = 0; j < len; j++)
        {
            mat1[i][j] = s1[j];
            mat2[i][j] = s2[len-j-1];
        }
        // populate mat3
        for (long j = 0; j < len-1; j++)
            mat3[i][j] = s2[len-j-2];
        mat3[i][len-1] = 0;
    }
    transpose(mat2,mat2);
    transpose(mat3,mat3);

    Mat<zz_pX> res1,res2;
    multiply(res1,mat1,mat2);
    multiply(res2,mat1,mat3);

    trunc(res1,res1,2*m);
    RightShift(res2,res2,2*m);
    quos = res1+res2;
    reverse(quos,quos,2*m-1);
}

void SetDims (Coeffs &res, const long r, const long c, const long d)
{
    res.SetLength(r);
    for (long i = 0; i < r; i++)
    {
        res[i].SetLength(c);
        for (long j = 0; j < c; j++)
            res[i][j].SetLength(d);
    }
}

void get_first_row (Coeffs &res,
        const zz_pX &a,
        const zz_pX &g,
        const long m)
{
    const long n = deg(g);
    const long d = 2*ceil((n*1.0)/m)+1;
    const long sqrt_d = ceil(sqrt(d));

    Vec<zz_pX> alphas, As;
    zz_pX t;
    SetCoeff(t,m-1,1);

    double t1 = GetWallTime();
    gen_pows(alphas, t, a, 1, g, sqrt_d);
    gen_pows(As, zz_pX(1), a, sqrt_d, g, sqrt_d);
    cout << "pow comp: " << GetWallTime() - t1 << endl;

    t1 = GetWallTime();
    Mat<zz_pX> quos;
    get_quos(quos, alphas, As, g, m);
    cout << "quo comp: " << GetWallTime() - t1 << endl;

    t1 = GetWallTime();
    SetDims(res, sqrt_d, sqrt_d, 2*m);
    zz_pX g_trunc;
    trunc(g_trunc,g,2*m);
    for (long u = 0; u < sqrt_d; u++)
        for (long v = 0; v < sqrt_d; v++)
        {
            zz_pX rem1, rem2;
            auto alpha_trunc = trunc(alphas[u],2*m);
            auto A_trunc = trunc(As[v],2*m);
            auto quo_trunc = trunc(quos[u][v],2*m);
            MulTrunc(rem1, alpha_trunc, A_trunc, 2*m);
            MulTrunc(rem2, quo_trunc, g_trunc, 2*m);
            rem1 = rem1 - rem2;
            for (long s = 0; s < 2*m; s++)
                res[u][v][s] = coeff(rem1,s);
        }
    cout << "rest from row1: " << GetWallTime() - t1 << endl;
}

void format (Vec<Mat<zz_p>> &res, const Vec<Coeffs> &coeffs, const long d, const long m)
{
    cout << "d: " << d << endl;
    res.SetLength(d);
    const long sqrt_d = ceil(sqrt(d));
    for (long i = 0; i < d; i++)
        res[i].SetDims(m,m);
    for (long i = 0; i < m; i++)
    {
        for (long u = 0; u < sqrt_d; u++)
        {
            for (long v = 0; v < sqrt_d; v++)
            {
                long pow = sqrt_d*v + u;
                for (long j = 0; j < m; j++)
                {
                    if (pow < d)
                        res[d-pow-1][m-i-1][j] = coeffs[i][u][v][j];
                }
            }
        }
    }
}

void get_coeffs (Vec<Mat<zz_p>> &mats,
        Vec<Coeffs> &res,
        const zz_pX &a,
        const zz_pX &g,
        const long m)
{
    const long n = deg(g);
    const long d = 2*ceil((n*1.0)/m)+1;
    const long sqrt_d = ceil(sqrt(d));

    res.SetLength(m);
    get_first_row(res[0],a,g,m);

    for (long i = 1; i < m; i++)
    {
        SetDims(res[i],sqrt_d,sqrt_d,2*m-i);
        for (long u = 0; u < sqrt_d; u++)
            for (long v = 0; v < sqrt_d; v++)
            {
                zz_p l_n_1 = -(zz_p(1)/ConstTerm(g)) * res[i-1][u][v][0];

                for (long j = 0; j < 2*m-i; j++)
                {
                    res[i][u][v][j] = res[i-1][u][v][j+1] 
                        + coeff(g,j+1)*l_n_1;
                }
            }
    }

    // format
    for (long i = 0; i < m; i++)
        SetDims(res[i],sqrt_d,sqrt_d,m);
        
    format(mats,res,d,m);
}

void get_coeffs_naive (Vec<Coeffs> &res,
        const zz_pX &a,
        const zz_pX &g,
        const long m)
{
    const long n = deg(g);
    const long d = 2*ceil((n*1.0)/m)+1;
    const long sqrt_d = ceil(sqrt(d));

    res.SetLength(m);

    zz_pXModulus g_mod;
    build(g_mod, g);
    for (long i = 0; i < m; i++)
    {
        SetDims(res[i], sqrt_d, sqrt_d, m);
        for (long u = 0; u < sqrt_d; u++)
            for (long v = 0; v < sqrt_d; v++)
            {
                zz_pX x;
                SetCoeff(x,m-i-1,1);
                zz_pX a_temp;
                PowerMod(a_temp, a, u+v*sqrt_d, g_mod);
                MulMod(a_temp, a_temp, x, g_mod);
                for (long t = 0; t < m; t++)
                    res[i][u][v][t] = coeff(a_temp,t);
            }
    }
}

void print(const Vec<Coeffs> &c)
{
    for (long i = 0; i < c.length(); i++)
        cout << c[i] << endl;
}

int main(int argc, char *argv[])
{
    if (argc!=5)
        throw std::invalid_argument("Usage: ./test_charpoly_amodg degree expnt nbits nthreads");

    // main parameter: degree
    long n = atoi(argv[1]);
    // parameter: exponent for m = n^exponent; default = 0.333 ~ 1/omega for omega=3
    double exponent = (atof(argv[2]) <= 0.) ? 0.333 : (atof(argv[2]));

    // size of prime field
    long nbits = atoi(argv[3]);
    // whether we should verify the result
    SetNumThreads(atoi(argv[4]));

    // used for timings
    double t1,t2;
    double t_charpoly=0.0;

    t1 = GetWallTime();
    if (nbits==0)
        zz_p::FFTInit(0);
    else
        zz_p::init(NTL::GenPrime_long(nbits));

    {
        t_charpoly=0.0;
        // modulus
        zz_pX g;
        while (deg(g)<n)
            random(g, n+1);
        MakeMonic(g);
        //SetCoeff(g, n, 1);

        // evaluation point
        zz_pX a;
        random(a, n); // deg < n: assumes a reduced modulo g

        // parameters
        long m = pow(n,exponent);
        long d = ceil( (double)n / m );
        t2 = GetWallTime();

        std::cout << "Testing charpoly with random input polynomials" << std::endl;
        std::cout << "-- n = " << n << std::endl;
        std::cout << "-- m = " << m << std::endl;
        std::cout << "-- d = " << d << std::endl;
        std::cout << "--prime = " << zz_p::modulus() << std::endl;
        std::cout << "--nthreads = " << AvailableThreads() << std::endl;
        if ((nbits==0 && n<30) || (nbits != 0 && n*nbits < 1000))
        {
            std::cout << "--modulus g = " << g << std::endl;
            std::cout << "--point a = " << a << std::endl;
        }
        else
        {
            std::cout << "--modulus g = <degree " << n << " polynomial>" << std::endl;
            std::cout << "--point a = <degree " << n << " polynomial>" << std::endl;
        }
        
        

#ifdef SAFETY_CHECKS
        std::cout << "###~~~Warning~~~### SAFETY_CHECKS is on; use with non-small degrees may be very slow." << std::endl;
#endif // SAFETY_CHECKS

        std::cout << "TIME ~~ build polynomials and compute parameters: " << (t2-t1) << std::endl;

        t1 = GetWallTime();
        zz_pXModulus G(g); // precompute things to work mod g
        zz_pXMultiplier A(a, G); // precompute things to multiply by a mod g
        t2 = GetWallTime();
        std::cout << "TIME ~~ pre-computations to work with a mod g: " << (t2-t1) << std::endl;
        t_charpoly += t2-t1;
        t1 = GetWallTime();

        // will hold powers of a, initially it is a
        zz_pX pow_a = a; 

#ifdef SAFETY_CHECKS
        // compute the block-Wiedemann sequence (naive)
        Vec<Mat<zz_p>> seq_naive;
        seq_naive.SetLength(2*d+1);
        ident(seq_naive[2*d], m); // pmat[2*d] = identity m x m
        for (long k = 2*d-1; k >= 0; --k)
        {
            seq_naive[k].SetDims(m, m);
            zz_pX f = pow_a;  // a^{2d-k} mod g
            for (long i = 0; i < m; ++i)
            {
                // here f = x^i * f mod g = x^i a^{2d-k} mod g
                for (long j = 0; j < m; ++j)
                    seq_naive[k][i][j] = coeff(f, j);  
                // note: here we have tranposed, because pmbasis works row-wise
                // f = x * f mod g  (recall g is monic)
                f <<= 1;
                f = f - coeff(f,n) * g;
            }
            MulMod(pow_a, pow_a, A, G); // a^{2d-(k+1)} mod g
        }
        t2 = GetWallTime();
        std::cout << "TIME ~~ compute sequence for block Wiedemann (naive): " << (t2-t1) << std::endl;
        pow_a = a; // set back to its state before this #ifdef
        t1 = GetWallTime();
#endif // SAFETY_CHECKS
/*
        // compute the block-Wiedemann sequence (fast)
        // Sequence stored as a list of constant matrices
        Vec<Mat<zz_p>> seq;
        seq.SetLength(2*d+1);

        // first m coefficients of -g
        zz_pX g_low;
        trunc(g_low, g, m);
        NTL::negate(g_low, g_low);

        // last m coefficients of -g
        Vec<zz_p> g_high; g_high.SetLength(m);
        for (long j = 0; j < m; ++j)
            g_high[j] = -g[j+n-m];

        // pmat[2*d] = identity m x m
        ident(seq[2*d], m);

        for (long k = 2*d-1; k >= 0; --k)
        {
            seq[k].SetDims(m, m);
            // first m coefficients of a^{2d-k} mod g
            zz_pX f_low;
            trunc(f_low, pow_a, m);
            // last m coefficients of a^{2d-k} mod g
            Vec<zz_p> f_high; f_high.SetLength(m);
            for (long j = 0; j < std::min(m,m+deg(pow_a)-n+1); ++j)
                f_high[j] = pow_a[n-m+j];

            for (long i = 0; i < m; ++i)
            {
                // here f_low = first m coefficients of x^i a^{2d-k} mod g
                VectorCopy(seq[k][i], f_low, m);  
                // now we want to compute the new f_low: first m coefficients
                // of x^{i+1} a^{2d-k} mod g
                f_low <<= 1;
                trunc(f_low, f_low, m);
                zz_p lt = f_high[m-1];
                f_low += lt * g_low;
                // also update f_high
                for (long j = m-1; j > 0; --j)
                    f_high[j] = f_high[j-1] + lt * g_high[j];
            }
            MulMod(pow_a, pow_a, A, G); // a^{2d-(k+1)} mod g
        }
        t2 = GetWallTime();
        cout << "seq:\n" << seq << endl;
*/ 

        // TEST ERIC'S ALGO
        Vec<Mat<zz_p>> seq;
        Vec<Coeffs> res;
        get_coeffs(seq,res,a,g,m);
        t2 = GetWallTime();
       
#ifdef SAFETY_CHECKS
        bool correct=true;
        for (long k = 0; k < 2*d+1; ++k)
            if (seq[k] != seq_naive[k])
                correct=false;
#endif // SAFETY_CHECKS
        std::cout << "TIME ~~ compute sequence for block Wiedemann: " << (t2-t1);
#ifdef SAFETY_CHECKS
        std::cout << (correct ? " (correct)" : " (wrong)");
#endif // SAFETY_CHECKS
        std::cout << std::endl;
        t_charpoly += t2-t1;
        t1 = GetWallTime();

        // convert the sequence into a polynomial matrix
        Mat<zz_pX> pmat;
        conv(pmat, seq, 2*d+1);

        // Matrix fraction reconstruction: add identity below pmat
        pmat.SetDims(2*m, m);
        for (long i = 0; i < m; ++i)
            set(pmat[i+m][i]);

        t2 = GetWallTime();
        std::cout << "TIME ~~ convert sequence to polynomial matrix: " << (t2-t1) << std::endl;
        t_charpoly += t2-t1;
        t1 = GetWallTime();

        // reconstruct fraction
        Mat<zz_pX> appbas;
        Shift shift(2*m, 0);
        DegVec pivdeg = pmbasis(appbas, pmat, 2*d+1, shift);

        // retrieve balanced basis (leading principal mxm submatrix)
        // be careful if wish to use numerator: here we have "-numer" in top-right submatrix,
        // since we have put identity instead of -identity in pmat
        Mat<zz_pX> basis;
        basis.SetDims(m, m);

        for (long j = 0; j < m; ++j)
            for (long i = 0; i < m; ++i)
                basis[j][i] = appbas[i][j];
        // note: we have transposed back to column-wise basis

        t2 = GetWallTime();
        std::cout << "TIME ~~ matrix fraction reconstruction: " << (t2-t1) << std::endl;
        t_charpoly += t2-t1;

#ifdef SAFETY_CHECKS
        t1 = GetWallTime();
        zz_pX det_pmbasis;
        determinant_generic_knowing_degree(det_pmbasis, basis, n);
        MakeMonic(det_pmbasis);
        t2 = GetWallTime();
        std::cout << "TIME ~~ determinant for charpoly: " << (t2-t1) << std::endl;
#endif // SAFETY_CHECKS

        // find determinant by solving system with random right hand side
        t1 = GetWallTime();
        zz_pX det;
        determinant_via_linsolve(det, basis);
        MakeMonic(det);
        t2 = GetWallTime();
        t_charpoly += t2-t1;

        std::cout << "TIME ~~ determinant for charpoly: " << (t2-t1);
#ifdef SAFETY_CHECKS
        std::cout << (det == det_pmbasis ? " (det consistent)" : " (det not consistent)") << std::endl;;
        correct = correct && (det == det_pmbasis);
#endif // SAFETY_CHECKS
        std::cout << std::endl;

        t1 = GetWallTime();
        zz_pX cp;
        CharPolyMod(cp, a, g); // cp = charpoly of a mod g
        t2 = GetWallTime();
        std::cout << "TIME ~~ NTL's charpoly of a mod g: " << (t2-t1) << std::endl;

        std::cout << "TIME ~~ charpoly of a mod g: " << t_charpoly << ((det == cp) ? " (correct)" : " (wrong)") << std::endl;

#ifdef SAFETY_CHECKS
        correct = correct && (det == cp);
        if (not correct)
            std::cout << "An issue was detected. Check messages above." << std::endl;
#endif // SAFETY_CHECKS
    }

    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
