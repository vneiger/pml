#include <NTL/matrix.h>
#include <NTL/mat_lzz_p.h>
#include <NTL/lzz_pX.h>
#include <memory>

#include "lzz_p_extra.h"
#include "mat_lzz_pX_extra.h"
#include "lzz_pX_CRT.h"

NTL_CLIENT


/*------------------------------------------------------------*/
/* solve A u = b mod x^prec                                   */
/* A square, A(0) invertible, deg(A), deg(b) < prec           */
/* invA = multiplier for A^(-1) mod x^thresh                  */
/*------------------------------------------------------------*/
static void solve_DAC(Mat<zz_pX>& sol, const Mat<zz_pX>& A, const Mat<zz_pX>& b, long prec,
               std::unique_ptr<mat_lzz_pX_lmultiplier> & invA, long thresh)
{
    if (prec <= thresh)
    {
        invA->multiply(sol, b);
        trunc(sol, sol, prec);
        return;
    }
    
    Mat<zz_pX> residue, rest;
    long hprec = prec / 2;
    long kprec = prec - hprec;
    solve_DAC(sol, trunc(A, hprec), trunc(b, hprec), hprec, invA, thresh);
    middle_product(residue, transpose(sol), transpose(A), hprec, kprec - 1);
    residue = transpose(residue);
    solve_DAC(rest, trunc(A, kprec), (b >> hprec) - residue, kprec, invA, thresh);
    sol += (rest << hprec);
}

/*------------------------------------------------------------*/
/* solve A u = b mod x^prec                                   */
/* A square, A(0) invertible, deg(A), deg(b) < prec           */
/* use when deg(A) close to prec                              */
/*------------------------------------------------------------*/
void solve_series_low_precision(Mat<zz_pX> &u, const Mat<zz_pX>& A, const Mat<zz_pX>& b, long prec)
{
    if (deg(A) >= prec || deg(b) >= prec)
        LogicError("A and b must be reduced in solve_series");

    long thresh = 8;
    Mat<zz_pX> invA = inv_trunc(A, thresh);
    std::unique_ptr<mat_lzz_pX_lmultiplier> mult = get_lmultiplier(invA, thresh);
    solve_DAC(u, A, b, prec, mult, thresh);
}


/*------------------------------------------------------------*/
/* solve A u = b mod x^prec                                   */
/* A square, A(0) invertible, deg(A), deg(b) < prec           */
/* use when deg(A) << prec                                    */
/*------------------------------------------------------------*/
void solve_series_high_precision(Mat<zz_pX> &u, const Mat<zz_pX>& A, const Mat<zz_pX>& b, long prec)
{
    if (deg(A) >= prec || deg(b) >= prec)
        LogicError("A and b must be reduced in solve_series");
    
    long dA = deg(A);
    long lenA = dA + 1;
    Mat<zz_pX> invA = inv_trunc(A, lenA);
    std::unique_ptr<mat_lzz_pX_lmultiplier> multI = get_lmultiplier(invA, dA);
    std::unique_ptr<mat_lzz_pX_lmultiplier> multA = get_lmultiplier(A, dA);
    
    long nb = prec / lenA;
    if (prec > (nb * lenA))
        nb++;

    long r = b.NumRows();
    long s = b.NumCols();

    u.SetDims(r, s);

    for (long a = 0; a < r; a++)
        for (long c = 0; c < s; c++)
            u[a][c].SetMaxLength(prec);

    Mat<zz_pX> sol;
    multI->multiply(sol, trunc(b, lenA));
    trunc(sol, sol, lenA);
    for (long a = 0; a < r; a++)
        for (long c = 0; c < s; c++)
            for (long i = 0; i <= dA; i++)
                SetCoeff(u[a][c], i, coeff(sol[a][c], i));

    long shift = lenA;
    Mat<zz_pX> upper;
    for (long i = 1; i < nb; i++)
    {
        multA->multiply(upper, sol);
        upper = upper >> lenA;
        for (long a = 0; a < r; a++)
            for (long c = 0; c < s; c++)
                for (long k = 0; k <= dA; k++)
                    SetCoeff(upper[a][c], k, coeff(b[a][c], k + shift) - coeff(upper[a][c], k));
        multI->multiply(sol, upper);
        trunc(sol, sol, lenA);
        for (long a = 0; a < r; a++)
            for (long c = 0; c < s; c++)
                for (long k = 0; k <= dA; k++)
                    SetCoeff(u[a][c], k + shift, coeff(sol[a][c], k));
        shift += lenA;
    }

}


// /*------------------------------------------------------------*/
// /* solve A u = b mod x^prec                                   */
// /* A must be square, A(0) invertible                          */
// /*------------------------------------------------------------*/
// void solve_series_regular(Mat<zz_pX> &u, const Mat<zz_pX>& A, const Mat<zz_pX>& b, long prec)
// {
//     long dA = deg(A);
//     if (prec <= dA)
//         linsolve_series_low_precision(u, trunc(A, prec), b, prec);
//     else
//         linsolve_series_high_precision(u, A, b, prec);
// }

// /*------------------------------------------------------------*/
// /* solve A u/den = b                                          */
// /* A must be square                                           */
// /*------------------------------------------------------------*/
// void solve(Mat<zz_pX> &u, zz_pX& den, const Mat<zz_pX>& A, const Mat<zz_pX>& b)
// {
//     long dA = deg(A);
//     long dB = deg(b);
//     long s = A.NumRows();
//     if (s != A.NumCols())
//         LogicError("A should be square in linsolve_rational_square");
//     if (s != B.NumRows())
//         LogicError("Dimension mismatch in linsolve_rational_square");

//     long deg_u = (s-1)*dA + dB;
//     long deg_den = s*dA;

//     // TODO: choose the initial point at random 
    

    
// }



// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
