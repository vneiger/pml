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
#include <sstream>

//#define SAFETY_CHECKS

#include "util.h"
#include "mat_lzz_pX_extra.h"
#include "mat_lzz_pX_approximant.h"

using namespace std;

NTL_CLIENT

typedef Vec<Vec<Vec<zz_p>>> Coeffs;

std::ostream &operator<<(std::ostream &out, const std::vector<long> &s)
{
    out << "[ ";
    for (auto &i: s)
        out << i << " ";
    return out << "]";
}

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
    
    pow.SetLength(l);
    pow[0] = t;
    for (long i = 1; i < l; i++)
    {
        MulMod(pow[i], pow[i-1], a_pow, g_mod);
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
    zz_pX x_pad;
    SetCoeff(x_pad, pad, 1);
    
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
        mul(A_st, A_st, x_pad);
        
        cout << "A_st: " << A_st << endl;
        
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
    cout << "mat1: " << mat1 << endl;
    cout << "mat2: " << mat2 << endl;
    cout << "mat3: " << mat3 << endl;
    
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
    const long d = ceil((n*1.0)/m);
    const long sqrt_d = ceil(sqrt(d));

    Vec<zz_pX> alphas, As;
    zz_pX t;
    SetCoeff(t,m-1,1);
    
    gen_pows(alphas, t, a, 1, g, sqrt_d);
    gen_pows(As, zz_pX(1), a, sqrt_d, g, sqrt_d);
    
    Mat<zz_pX> quos;
    get_quos(quos, alphas, As, g, m);
    
    SetDims(res, sqrt_d, sqrt_d, 2*m);
    t = zz_pX(0);
    SetCoeff(t,2*m,1);
    for (long u = 0; u < sqrt_d; u++)
        for (long v = 0; v < sqrt_d; v++)
        {
            zz_pX rem = (alphas[u] * As[v] - quos[u][v]*g)%t;
            for (long s = 0; s < 2*m; s++)
                res[u][v][s] = coeff(rem,s);
        }
}

void get_coeffs (Vec<Coeffs> &res,
                 const zz_pX &a,
                 const zz_pX &g,
                 const long m)
{
    const long n = deg(g);
    const long d = ceil((n*1.0)/m);
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
}

void print(const Vec<Coeffs> &c)
{
    for (long i = 0; i < c.length(); i++)
        cout << c[i] << endl;
}

int main(){
    zz_p::init(13);
    cout << "enter a:" << endl;
    zz_pX a;
    string s;
    getline(cin, s);
    istringstream iss{s};
    long t;
    long at = 0;
    while (iss >> t)
        SetCoeff(a, at++, t);
        
    cout << "enter g:" << endl;
    zz_pX g;
    getline(cin,s);
    at = 0;
    iss = istringstream(s);
    while(iss >> t)
        SetCoeff(g, at++, t);
        
    cout << "enter m:" << endl;
    long m;
    cin >> m;

    cout << "a: " << a << endl;
    cout << "g: " << g << endl;
    cout << "m: " << m << endl;
    
    Vec<Coeffs> res;
    get_coeffs(res,a,g,m);
    cout << "res:" << endl;
    print(res);
}











































// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
