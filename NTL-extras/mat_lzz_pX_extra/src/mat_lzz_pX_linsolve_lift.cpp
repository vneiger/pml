#include <NTL/matrix.h>
#include <NTL/mat_lzz_p.h>
#include <NTL/lzz_pX.h>
#include <memory>

#include "util.h"
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
/* A square, A(0) invertible                                  */
/* use when deg(A) close to prec                              */
/* computes A^-1 mod x^{2^thresh}                             */
/* thresh=-1 is the default value, uses lookup table          */
/* output can alias input                                     */
/*------------------------------------------------------------*/
void solve_series_low_precision(Mat<zz_pX> &u, const Mat<zz_pX>& A, const Mat<zz_pX>& b, long prec, long thresh)
{
    if (&u == &A || &u == &b)
    {
        u = solve_series_low_precision(A, b, prec, thresh);
        return;
    }

    if (deg(A) >= prec || deg(b) >= prec)
    {
        solve_series_low_precision(u, trunc(A, prec), trunc(b, prec), prec, thresh);
        return;
    }

    if (thresh == -1)
    {
        long t = type_of_prime();
        switch(t)
        {
        case TYPE_FFT_PRIME:
            thresh = THRESHOLDS_SOLVE_LOW_PRECISION_FFT;
            break;
        case TYPE_SMALL_PRIME:
            thresh = THRESHOLDS_SOLVE_LOW_PRECISION_SMALL;
            break;
        case TYPE_LARGE_PRIME:
            thresh = THRESHOLDS_SOLVE_LOW_PRECISION_LARGE;
            break;
        default:
            LogicError("Unknown prime type in linear solving.");
        }
    }
    Mat<zz_pX> invA = inv_trunc(A, thresh);
    std::unique_ptr<mat_lzz_pX_lmultiplier> mult = get_lmultiplier(invA, thresh);
    solve_DAC(u, A, b, prec, mult, thresh);
}


/*------------------------------------------------------------*/
/* solve A u = b mod x^prec                                   */
/* A square, A(0) invertible, deg(A), deg(b) < prec           */
/* use when deg(A) << prec                                    */
/* output can alias input                                     */
/*------------------------------------------------------------*/
void solve_series_high_precision(Mat<zz_pX> &u, const Mat<zz_pX>& A, const Mat<zz_pX>& b, long prec)
{
    if (&u == &A || &u == &b)
    {
        u = solve_series_high_precision(A, b, prec);
        return;
    }
    if (deg(A) >= prec || deg(b) >= prec)
    {
        solve_series_high_precision(u, trunc(A, prec), trunc(b, prec), prec);
        return;
    }

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
        {
            u[a][c] = 0;
            u[a][c].SetMaxLength(prec);
        }

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
    trunc(u, u, prec);
}

/*------------------------------------------------------------*/
/* solve A u = b mod x^prec                                   */
/* A must be square, A(0) invertible                          */
/* output can alias input                                     */
/*------------------------------------------------------------*/
void solve_series(Mat<zz_pX> &u, const Mat<zz_pX>& A, const Mat<zz_pX>& b, long prec)
{
    long dA = deg(A);
    if (prec <= 4 * dA)  // seems reasonable
        solve_series_low_precision(u, A, b, prec);
    else
        solve_series_high_precision(u, A, b, prec);
}


/*------------------------------------------------------------*/
/* Implements a minor variation of Storjohann's algorithm     */
/* A must be square, A(0) invertible, deg(b) < deg(A)         */
/* output can alias input                                     */
/*------------------------------------------------------------*/
void solve_series_high_order_lifting(Mat<zz_pX> &u, const Mat<zz_pX>& A, const Mat<zz_pX>& b, long prec)
{
    Mat<zz_pX> slice, sol;
    std::unique_ptr<mat_lzz_pX_lmultiplier> ma, minvA;
    long d = deg(A);
    long ncols = prec / d;
    if (prec > ncols*d)
        ncols++;
    long bcols = b.NumCols();

    if (deg(b) >= deg(A))
        LogicError("Bad degrees for high order lifting");

    ma = get_lmultiplier(A, d-1);
    minvA = get_lmultiplier(inv_trunc(A, d), d-1);

    slice = inv_trunc(A, 2*d) >> 1;
    sol = b;
   
    while( sol.NumCols() < ncols*bcols )
    {
        Mat<zz_pX> next = transpose(middle_product(transpose(sol), transpose(slice), d-1, d-1)); // deg(next) < d
        sol = horizontal_join(sol, trunc(ma->multiply(next), d)); // deg(sol) < d
        if (sol.NumCols() < ncols*bcols)
            high_order_lift_inverse_odd(slice, slice, ma, minvA, d);
    }
    sol = trunc(minvA->multiply(sol), d);
    u = collapse_nonconsecutive_columns(sol, d, bcols);
    trunc(u, u, prec);
}


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
