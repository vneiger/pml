#include "structured_lzz_p.h"
#include "thresholds_solve_lift.h"
#include "mat_lzz_pX_extra.h"
#include "mat_lzz_pX_utils.h"
#include "mat_lzz_pX_inverse.h"
#include "mat_lzz_pX_arith.h"
#include "mat_lzz_pX_multiply.h"
#include "mat_lzz_pX_linsolve.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* solve A u = b mod x^prec                                   */
/* A square, A(0) invertible, deg(A), deg(b) < prec           */
/* invA = multiplier for A^(-1) mod x^thresh                  */
/*------------------------------------------------------------*/
// recursion used in solve_series_low_precision
static void solve_DAC(Mat<zz_pX>& sol, const Mat<zz_pX>& A, const Mat<zz_pX>& b, long prec,
               const std::unique_ptr<mat_lzz_pX_lmultiplier> & invA, long thresh)
{
    // if precision is small enough, just multiply by A^(-1) mod x^thresh
    if (prec <= thresh)
    {
        invA->multiply(sol, b);
        trunc(sol, sol, prec);
        return;
    }

    // precision for recursive calls
    const long hprec = prec/2;
    const long kprec = prec-hprec;
    
    // buffer for matrices related to A
    Mat<zz_pX> bufA;
    // buffer for matrices related to b
    Mat<zz_pX> bufB;

    // first recursive call
    trunc(bufA, A, hprec);
    trunc(bufB, b, hprec);
    solve_DAC(sol, bufA, bufB, hprec, invA, thresh);

    // compute residue
    transpose(bufA, A);
    transpose(bufB, sol);
    Mat<zz_pX> residue;
    middle_product(residue, bufB, bufA, hprec, kprec-1);
    transpose(bufB, residue);

    trunc(bufA, A, kprec);
    RightShift(residue,b,hprec);
    sub(bufB, residue, bufB);
    solve_DAC(residue, bufA, bufB, kprec, invA, thresh);

    // sol += (bufres << hprec);
    add_LeftShift(sol, sol, residue, hprec);
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
        Mat<zz_pX> v;
        solve_series_low_precision(v, A, b, prec, thresh);
        u.swap(v);
        return;
    }

    if (deg(A) >= prec || deg(b) >= prec)
    {
        solve_series_low_precision(u, trunc(A, prec), trunc(b, prec), prec, thresh);
        return;
    }

    if (thresh == -1)
    {
        switch(type_of_prime())
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

#ifdef VERBOSE
    double t = get_time();
#endif
    Mat<zz_pX> invA = inv_trunc(A, thresh);
#ifdef VERBOSE
    cout << "inv_trunc " << get_time()-t << endl;
#endif
    
#ifdef VERBOSE
    t = get_time();
#endif
    std::unique_ptr<mat_lzz_pX_lmultiplier> mult = get_lmultiplier(invA, thresh);
#ifdef VERBOSE
    cout << "get mult " << get_time()-t << endl;
#endif

#ifdef VERBOSE
    t = get_time();
#endif
    solve_DAC(u, A, b, prec, mult, thresh);
#ifdef VERBOSE
    cout << "solve DAC " << get_time()-t << endl;
#endif
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

#ifdef VERBOSE
    double t = get_time();
#endif
    Mat<zz_pX> invA = inv_trunc(A, lenA);
#ifdef VERBOSE
    cout << "inv_trunc " << get_time()-t << endl;
#endif


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
    if (prec <= 4 * deg(A))  // seems reasonable
        solve_series_low_precision(u, A, b, prec);
    else
        solve_series_high_precision(u, A, b, prec);
}


/*------------------------------------------------------------*/
/* solve A u = b mod x^prec                                   */
/* A must be square, A(0) invertible                          */
/* output can alias input                                     */
/*------------------------------------------------------------*/
void solve_series(Vec<zz_pX> &u, const Mat<zz_pX>& A, const Vec<zz_pX>& b, long prec)
{
    Mat<zz_pX> bmat, umat;
    long n = b.length();
    bmat.SetDims(n, 1);
    for (long i = 0; i < n; i++)
        bmat[i][0] = b[i];
    solve_series(umat, A, bmat, prec);
    u.SetLength(n);
    for (long i = 0; i < n; i++)
        u[i] = umat[i][0];
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




/*------------------------------------------------------------*/
/* solve A (u/den) = b                                        */
/* A must be square, A(0) invertible                          */
/* output can alias input                                     */
/* uses lifting and rational reconstruction                   */
/*------------------------------------------------------------*/
long linsolve_via_series(Vec<zz_pX> &u, zz_pX& den, const Mat<zz_pX>& A, const Vec<zz_pX>& b, long nb_max)
{
    long n = A.NumRows();

    if (A.NumCols() != n)
        LogicError("Need square matrix for linsolve series");
    if (n == 0)
        LogicError("Empty matrix for linsolve series");
    if (b.length() != n)
        LogicError("Bad vector size for linsolve series");


#ifdef VERBOSE
    cout << endl;
    double t = get_time();
#endif

    long dA = deg(A);
    long dB = deg(b);

    if (dA == 0)
    {
        Mat<zz_p> iA = inv(coeff(A, 0));  // TODO: check for invertibility
        u.SetLength(n);
        for (long i = 0; i < n; i++)
        {
            zz_pX tmp;
            clear(tmp);
            for (long j = 0; j < n; j++)
                tmp += iA[i][j] * b[j];
            u[i] = tmp;
        }
        den = 1;
        return 1;
    }

    // DIRT: need tuning
    if (nb_max == -1)
    {
        nb_max = 1;
        long t = type_of_prime();
        if (t == TYPE_FFT_PRIME)
        {
            if (n >= 100)
                nb_max = 4;
        }
        if (t == TYPE_LARGE_PRIME)
        {
            if (n >= 100)
                nb_max = 3;
        }
        if (t == TYPE_SMALL_PRIME)
        {
            if (n >= 100)
                nb_max = 3;
        }
    }


    long nb = min(n, nb_max); // number of linear combinations we are taking
    long deg_den = n*dA; // deg_den + 1 = number of unknowns
    long deg_num = (n-1)*dA + dB;

    long first = max(deg_num + 1, deg_den); // first term we can use in each block
    long sz = (deg_den + nb - 1) / nb;  // size of each block
    long prec = first + sz;

    Vec<zz_pX> sol_series, lin_comb;

    solve_series(sol_series, A, b, prec);

#ifdef VERBOSE
    cout << "setup:  " << get_time()-t << endl;
    cout << "nb=" << nb << endl;
#endif

#ifdef VERBOSE
    t = get_time();
#endif

    lin_comb.SetLength(nb);
    for (long i = 0; i < nb; i++)
    {
        zz_pX res;
        for (long j = 0; j < n; j++)
            res += random_zz_p() * sol_series[j];
        lin_comb[i] = res;
    }

    if (nb == 1)
    {
        zz_pX b;
        zz_pXMatrix M;
        zz_p t;
        SetCoeff(b, deg_den + deg_num + 1);
        HalfGCD(M, b, trunc(lin_comb[0], deg_den + deg_num + 1), deg_den+1);
        inv(t, LeadCoeff(M(1,1)));
        mul(den, M(1,1), t);
    }
    else
    {
        Vec<Vec<toeplitz_lzz_p>> sys;
        sys.SetLength(nb);
        for (long i = 0; i < nb; i++)
        {
            sys[i].SetLength(1);
            Vec<zz_p> coeffs;
            coeffs.SetLength(deg_den + sz);
            for (long j = 0; j < deg_den + sz; j++)
                coeffs[j] = coeff(lin_comb[i], first - deg_den + j);
            sys[i][0] = toeplitz_lzz_p(coeffs, sz, deg_den + 1);
        }
        
        mosaic_toeplitz_lzz_p MT = mosaic_toeplitz_lzz_p(sys);
        Vec<zz_p> ker_vec, zero;
        zero.SetLength(sz * nb);
        for (long i = 0; i < sz * nb; i++)
            zero[i] = 0;

        long ans = MT.solve(ker_vec, zero);
        if (ans == 0)
        {
            Error("Error in solving mosaic toeplitz");
        }
        den = 0;
        for (long i = deg_den; i >= 0; i--)
            SetCoeff(den, i, ker_vec[i]);

    }

#ifdef VERBOSE
    cout << "find denom " << get_time()-t << endl;
#endif

    u.SetLength(n);
    for (long i = 0; i < n; i++)
        u[i] = MulTrunc(trunc(sol_series[i], deg_num+1), trunc(den, deg_num+1), deg_num+1);
        
    return 1;
}


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
