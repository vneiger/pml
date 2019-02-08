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
#include "mat_lzz_pX_approximant.h"
#include "sage_output.h"

NTL_CLIENT

std::ostream &operator<<(std::ostream &out, const VecLong &s)
{
    out << "[ ";
    for (auto &i: s)
        out << i << " ";
    return out << "]";
}

typedef Vec<Vec<Vec<zz_p>>> Coeffs;

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
void print(const Vec<Coeffs> &c)
{
    for (long i = 0; i < c.length(); i++)
        cout << c[i] << endl;
}

void print(const zz_pX &p)
{
    for (long i = 0; i <= deg(p); i++)
    {
        cout << coeff(p,i) << "*y^" <<i;
        if (i != deg(p)) cout << "+";
        else cout << endl;
    }
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

        // evaluation point
        zz_pX a;
        random(a, n); // deg < n: assumes a reduced modulo g

        zz_pX h;
        random(h, n);

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
            std::cout << "--point h = " << h << std::endl;
        }
        else
        {
            std::cout << "--modulus g = <degree " << n << " polynomial>" << std::endl;
            std::cout << "--point a = <degree " << n << " polynomial>" << std::endl;
            std::cout << "--point h = <degree " << n << " polynomial>" << std::endl;
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

        Vec<Mat<zz_p>> seq;
        Vec<Coeffs> res;
        gen_sequence(seq,a,g,m);
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
        VecLong shift(2*m, 0);
        pmbasis(appbas, pmat, 2*d+1, shift);

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

        // find h(x) mod the balanced basis
        t1 = GetWallTime();
        Mat<zz_pX> H_mat;
        H_mat.SetDims(m,1);
        H_mat[0][0] = h;
        Mat<zz_pX> Q,hh;
        quo_rem(Q,hh, H_mat, basis);

        Vec<zz_pX> H;
        H.SetLength(d);
        for (long k = 0; k < d; k++)
        {
            for (long i = 0; i < m; i++)
                SetCoeff(H[k],i,coeff(hh[i][0],k));
        }
        t2 = GetWallTime();
        std::cout << "TIME ~~ quo_rem: " << (t2-t1) << std::endl;
        t_charpoly += t2-t1;
        t1 = GetWallTime();

        // bivariate modular composition, smaller degree
        zz_pX modcomp{zz_p(0)};
        //TODO: use the pre-computed powers
        Vec<zz_pX> a_pows;
        Vec<zz_pX> upper;
        gen_pows(a_pows, upper, zz_pX{1}, a, 1, g, d, 0);
        for (long i = 0; i < d; i++)
        {
            modcomp += H[i] * a_pows[i];
        }

        modcomp %= g;
        t2 = GetWallTime();
        std::cout << "TIME ~~ biv modcomp: " << (t2-t1) << std::endl;
        t_charpoly += t2-t1;
        {
            zz_pXModulus g_mod;
            build(g_mod,g);
            zz_pX modcomp_ntl;
            t1 = GetWallTime();
            CompMod(modcomp_ntl, h, a, g_mod);
            cout << "TIME ~~ NTL h(a) mod g: " << GetWallTime()-t1<< endl;
            cout << "TIME ~~ h(a) mod g: " << t_charpoly;
            if (modcomp_ntl == modcomp)
                cout << " (correct)" << endl;
            else
                cout << " (wrong)" << endl;
        }
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
