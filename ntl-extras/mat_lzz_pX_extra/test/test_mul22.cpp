#include <NTL/lzz_pX.h>
#include <NTL/tools.h>
#include <iomanip>

#include "util.h"
#include "mat_lzz_pX_extra.h"

NTL_CLIENT

// copied from NTL's lzz_pX1.cpp
void mul(zz_pXMatrix& A, zz_pXMatrix& B, zz_pXMatrix& C)
// A = B*C, B and C are destroyed
{
   long db = deg(B(1,1));
   long dc = deg(C(1,1));
   long da = db + dc;

   long k = NextPowerOfTwo(da+1);

   fftRep B00, B01, B10, B11, C0, C1, T1, T2;
   
   TofftRep(B00, B(0,0), k); B(0,0).kill();
   TofftRep(B01, B(0,1), k); B(0,1).kill();
   TofftRep(B10, B(1,0), k); B(1,0).kill();
   TofftRep(B11, B(1,1), k); B(1,1).kill();

   TofftRep(C0, C(0,0), k);  C(0,0).kill();
   TofftRep(C1, C(1,0), k);  C(1,0).kill();

   mul(T1, B00, C0);
   mul(T2, B01, C1);
   add(T1, T1, T2);
   FromfftRep(A(0,0), T1, 0, da);

   mul(T1, B10, C0);
   mul(T2, B11, C1);
   add(T1, T1, T2);
   FromfftRep(A(1,0), T1, 0, da);

   TofftRep(C0, C(0,1), k);  C(0,1).kill();
   TofftRep(C1, C(1,1), k);  C(1,1).kill();

   mul(T1, B00, C0);
   mul(T2, B01, C1);
   add(T1, T1, T2);
   FromfftRep(A(0,1), T1, 0, da);

   mul(T1, B10, C0);
   mul(T2, B11, C1);
   add(T1, T1, T2);
   FromfftRep(A(1,1), T1, 0, da);
}

// to convert from our matrices to NTL's ones
// assumes b is 2x2
void convert(zz_pXMatrix & a, const Mat<zz_pX> & b)
{
    a(0,0) = b[0][0];
    a(0,1) = b[0][1];
    a(1,0) = b[1][0];
    a(1,1) = b[1][1];
}


int main(int argc, char ** argv)
{
    if (argc!=3)
        throw std::invalid_argument("Usage: ./test_mul22 deg nbits");

    double t_ntl=0.0, t_pml=0.0, t_naive=0.0;
    long iter = 0;
    double t_now;

    long d = atoi(argv[1]);
    long nbits = atoi(argv[2]);

    if (nbits==0)
        zz_p::FFTInit(0);
    else
        zz_p::init(NTL::GenPrime_long(nbits));

    std::cout << "2x2 polynomial matrix multiplication, over prime field p = " << zz_p::modulus() << std::endl;
    std::cout << "degree d =" << d << std::endl;

    while (t_pml<0.5)
    {
        Mat<zz_pX> a = random_mat_zz_pX(2,2,d);
        Mat<zz_pX> b = random_mat_zz_pX(2,2,d);
        Mat<zz_pX> c;
        t_now = GetWallTime();
            multiply(c,a,b);
        t_pml += GetWallTime()-t_now;
        ++iter;
    }

    std::cout << "Time (pml):" << (t_pml / iter) << std::endl;

    iter=0;
    while (t_naive<0.5)
    {
        Mat<zz_pX> a = random_mat_zz_pX(2,2,d);
        Mat<zz_pX> b = random_mat_zz_pX(2,2,d);
        Mat<zz_pX> c(INIT_SIZE, 2, 2);
        t_now = GetWallTime();
            c[0][0] = a[0][0]*b[0][0] + a[0][1]*b[1][0];
            c[0][1] = a[0][0]*b[0][1] + a[0][1]*b[1][1];
            c[1][0] = a[1][0]*b[0][0] + a[1][1]*b[1][0];
            c[1][1] = a[1][0]*b[0][1] + a[1][1]*b[1][1];
        t_naive += GetWallTime()-t_now;
        ++iter;
    }

    std::cout << "Time (naive):" << (t_naive / iter) << std::endl;

    iter=0;
    while (t_ntl<0.5)
    {
        Mat<zz_pX> A = random_mat_zz_pX(2,2,d);
        Mat<zz_pX> B = random_mat_zz_pX(2,2,d);
        zz_pXMatrix a, b;
        convert(a, A); convert(b, B);
        zz_pXMatrix c;
        t_now = GetWallTime();
            mul(c, a, b);
        t_ntl += GetWallTime()-t_now;
        ++iter;
    }

    std::cout << "Time (ntl):" << (t_ntl / iter) << std::endl;

    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
