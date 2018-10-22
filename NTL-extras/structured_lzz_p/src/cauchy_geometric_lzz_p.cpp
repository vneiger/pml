#include <NTL/mat_lzz_p.h>
#include <NTL/lzz_pX.h>

#include "vec_lzz_p_extra.h"
#include "structured_lzz_p.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* Cauchy matrices on geometric progressions                  */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* computes                                                   */
/* 1/(u1-v1 rho^(n+1)) ... 1/(u1-v1 rho^(-m-1))               */
/* these are the entries of the toeplitz matrix               */
/* (with m rows and n columns)                                */
/*------------------------------------------------------------*/
static void prepare_inverses_cauchy_reverse(Vec<zz_p>& inverses, const zz_p& u1, const zz_p& v1, const zz_p& rho, long m, long n)
{
    Vec<zz_p> vec_den;
    zz_p irho = 1/rho;
    vec_den.SetLength(m+n-1);
    if (m+n-1 == 0)
    {
        return;
    }
    vec_den[m+n-2-0] = -v1*power(irho, m-1);
    for (long i = 1; i < m+n-1; i++)
    {
        vec_den[m+n-2-i] = vec_den[m+n-2-(i-1)]*rho;
        vec_den[m+n-2-(i-1)] += u1;
    }
    vec_den[0] += u1;
    inv(inverses, vec_den);
}

/*------------------------------------------------------------*/
/* default constructor                                        */
/*------------------------------------------------------------*/
cauchy_geometric_lzz_p::cauchy_geometric_lzz_p() 
{
    m = 0;
    n = 0;
}

/*------------------------------------------------------------*/
/* constructor                                                */
/* M[i][j] = 1 / (u1*rho^i - v1*rho^j)                        */
/* rho should be a square                                     */
/*------------------------------------------------------------*/
cauchy_geometric_lzz_p::cauchy_geometric_lzz_p(const zz_p& a1, const zz_p& b1, const zz_p& q, long mm, long nn)
{
    u1 = a1;
    v1 = b1;
    rho = q;
    sqrt_rho = 0;
    m = mm;
    n = nn;

    if (m == 0 && n == 0)
    {
        t = toeplitz_lzz_p();
        powers_irho.SetLength(0);
    }
    else
    {
        Vec<zz_p> vec_toeplitz;
        prepare_inverses_cauchy_reverse(vec_toeplitz, u1, v1, rho, m, n); 

        powers_irho.SetLength(m);
        zz_p irho = 1/rho;
        powers_irho[0] = 1;
        for (long i = 1; i < m; i++)
        {
            powers_irho[i] = irho*powers_irho[i-1];
        }
        t = toeplitz_lzz_p(vec_toeplitz, m, n);
    }
}

/*------------------------------------------------------------*/
/* dimensions                                                 */
/*------------------------------------------------------------*/
long cauchy_geometric_lzz_p::NumRows() const 
{
    return m;
}

long cauchy_geometric_lzz_p::NumCols() const 
{
    return n;
}

/*------------------------------------------------------------*/
/* computes output = M*input                                  */
/*------------------------------------------------------------*/
void cauchy_geometric_lzz_p::mul_right(Vec<zz_p>& output, const Vec<zz_p>& input) const 
{
    if (&output == &input)
    {
        output = mul_right(input);
        return;
    }

    t.mul_right(output, input);
    for (long i = 0; i < m; i++)
        output[i] *= powers_irho[i];
}

/*------------------------------------------------------------*/
/* computes output = M*input                                  */
/*------------------------------------------------------------*/
void cauchy_geometric_lzz_p::mul_right(Mat<zz_p>& output, const Mat<zz_p>& input) const 
{
    if (&output == &input)
    {
        output = mul_right(input);
        return;
    }

    t.mul_right(output, input);
    long a = input.NumCols();
    for (long i = 0; i < m; i++)
        for (long j = 0; j < a; j++)
            output[i][j] *= powers_irho[i];
}

/*------------------------------------------------------------*/
/* computes output = M*input without the diagonal             */
/*------------------------------------------------------------*/
void cauchy_geometric_lzz_p::mul_right_simple(Vec<zz_p>& output, const Vec<zz_p>& input) const 
{
    if (&output == &input)
    {
        output = mul_right_simple(input);
        return;
    }

    t.mul_right(output, input);
}

/*------------------------------------------------------------*/
/* computes output = M*input without the diagonal             */
/*------------------------------------------------------------*/
void cauchy_geometric_lzz_p::mul_right_simple(Mat<zz_p>& output, const Mat<zz_p>& input) const 
{
    if (&output == &input)
    {
        output = mul_right_simple(input);
        return;
    }

    t.mul_right(output, input);
}

/*------------------------------------------------------------*/
/* computes output = input * M                                */
/*------------------------------------------------------------*/
void cauchy_geometric_lzz_p::mul_left(Vec<zz_p>& output, const Vec<zz_p>& input) const 
{
    Vec<zz_p> new_in = input;
    for (long i = 0; i < m; i++)
    {
        new_in[i] *= powers_irho[i];
    }
    t.mul_left(output, new_in);
}

/*------------------------------------------------------------*/
/* computes output = input * M                                */
/*------------------------------------------------------------*/
void cauchy_geometric_lzz_p::mul_left(Mat<zz_p>& output, const Mat<zz_p>& input) const 
{
    Mat<zz_p> new_in = input;
    long a = input.NumRows();
    for (long i = 0; i < a; i++)
    {
        zz_p *elts = new_in[i].elts();
        for (long j = 0; j < m; j++)
        {
            elts[j] *= powers_irho[j];
        }
    }
    t.mul_left(output, new_in);
}

/*------------------------------------------------------------*/
/* computes output = input * M without the diagonal           */
/*------------------------------------------------------------*/
void cauchy_geometric_lzz_p::mul_left_simple(Vec<zz_p>& output, const Vec<zz_p>& input) const 
{
    if (&output == &input)
    {
        output = mul_left_simple(input);
        return;
    }

    t.mul_left(output, input);
}

/*------------------------------------------------------------*/
/* computes output = input * M without the diagonal           */
/*------------------------------------------------------------*/
void cauchy_geometric_lzz_p::mul_left_simple(Mat<zz_p>& output, const Mat<zz_p>& input) const 
{
    if (&output == &input)
    {
        output = mul_left_simple(input);
        return;
    }

    t.mul_left(output, input);
}

/*------------------------------------------------------------*/
/* M as a dense matrix                                        */
/*------------------------------------------------------------*/
void cauchy_geometric_lzz_p::to_dense(Mat<zz_p>& M) const 
{
    M.SetDims(m, n);
    for (long i = 0; i < m; i++)
        for (long j = 0; j < n; j++)
            M[i][j] = to_zz_p(1) / (u1*power(rho, i) - v1*power(rho, j));
}

// /*------------------------------------------------------------*/
// /* builds the left / right vandermonde matrices               */
// /*------------------------------------------------------------*/
// void cauchy_geometric_lzz_p::build_X_Y()
// {
//     // TODO: check if it exists!
//     if (sqrt_rho == 0)
//     {
//         long rho_l = rho._zz_p__rep;
//         sqrt_rho = to_zz_p(to_long(SqrRootMod(to_ZZ(rho_l), to_ZZ(zz_p::modulus()))));
//     }

//     X = zz_pX_Multipoint_Geometric(sqrt_rho, u1, m);
//     Y = zz_pX_Multipoint_Geometric(sqrt_rho, v1, n);
// }



// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
