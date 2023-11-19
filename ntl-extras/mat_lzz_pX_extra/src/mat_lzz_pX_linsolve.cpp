#include "thresholds_solve_lift.h"

#include "structured_lzz_p.h" // for mosaic system solving
#include "mat_lzz_pX_utils.h"
#include "mat_lzz_pX_inverse.h" // for inv_trunc
#include "mat_lzz_pX_arith.h"
#include "mat_lzz_pX_multiply.h"
#include "mat_lzz_pX_linearization.h" // for collapse and join
#include "mat_lzz_pX_kernel.h" // for kernel basis

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

    // second recursive call
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
            LogicError("Unknown prime type in linear solving low precision.");
        }
    }

    // compute inverse of A mod X^thresh
    Mat<zz_pX> invA = inv_trunc(A, thresh);

    // prepare the left-multiplier by A^{-1}
    std::unique_ptr<mat_lzz_pX_lmultiplier> mult = get_lmultiplier(invA, thresh-1);

    // run the actual recursion
    solve_DAC(u, A, b, prec, mult, thresh);
}


/*------------------------------------------------------------*/
/* solve A u = b mod x^prec                                   */
/* A square, A(0) invertible                                  */
/* use when deg(A) << prec                                    */
/* output can alias input                                     */
/*------------------------------------------------------------*/
void solve_series_high_precision(Mat<zz_pX> &u, const Mat<zz_pX>& A, const Mat<zz_pX>& b, long prec)
{
    if (&u == &A || &u == &b)
    {
        Mat<zz_pX> v;
        solve_series_high_precision(v, A, b, prec);
        u.swap(v);
        return;
    }
    if (deg(A) >= prec || deg(b) >= prec)
    {
        solve_series_high_precision(u, trunc(A, prec), trunc(b, prec), prec);
        return;
    }

    const long dA = deg(A);
    const long lenA = dA + 1;

    // compute inverse of A truncated at X^lenA
    Mat<zz_pX> invA = inv_trunc(A, lenA);

    // prepare left-multipliers for multiplication by invA and by A
    std::unique_ptr<mat_lzz_pX_lmultiplier> multI = get_lmultiplier(invA, dA);
    std::unique_ptr<mat_lzz_pX_lmultiplier> multA = get_lmultiplier(A, dA);
    
    // nb = ceil(prec/lenA)
    const long nb = 1 + (prec-1) / lenA;

    const long r = b.NumRows();
    const long s = b.NumCols();

    // initialize u to the right dimension, and reserve space
    u.SetDims(r, s);
    for (long i = 0; i < r; ++i)
        for (long j = 0; j < s; ++j)
        {
            clear(u[i][j]);
            u[i][j].SetMaxLength(prec);
        }

    // u = A^{-1} b mod X^lenA
    multI->multiply(u, trunc(b, lenA));
    trunc(u, u, lenA);

    // temporary matrix, initially copy of u
    Mat<zz_pX> sol(u);
    // will store high-degree part of A*sol
    Mat<zz_pX> upper;

    // we will left shift by lenA, 2lenA, 3lenA, ...
    long shift = lenA;

    for (long it = 1; it < nb; ++it)
    {
        // upper = part of nonnegative degree of X^{-lenA} A*sol
        multA->multiply(upper, sol);
        RightShift(upper, upper, lenA);
        // upper = - upper + (X^{-shift} b % X^lenA)
        for (long i = 0; i < r; ++i)
            for (long j = 0; j < s; ++j)
                for (long k = 0; k <= dA; ++k)
                    SetCoeff(upper[i][j], k, coeff(b[i][j], k + shift) - coeff(upper[i][j], k));
        // sol = A^{-1} * upper mod X^lenA
        multI->multiply(sol, upper);
        trunc(sol, sol, lenA);
        // u += (sol << shift)
        add_LeftShift(u, u, sol, shift);
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
    const long n = b.length();
    bmat.SetDims(n, 1);
    for (long i = 0; i < n; ++i)
        bmat[i][0] = b[i];
    solve_series(umat, A, bmat, prec);
    u.SetLength(n);
    for (long i = 0; i < n; ++i)
        u[i].swap(umat[i][0]);
}


/*------------------------------------------------------------*/
/* Implements a minor variation of Storjohann's algorithm     */
/* A must be square, A(0) invertible, deg(b) < deg(A)         */
/* output can alias input                                     */
/*------------------------------------------------------------*/
void solve_series_high_order_lifting(Mat<zz_pX> &u, const Mat<zz_pX>& A, const Mat<zz_pX>& b, long prec)
{
    // degree of A and number of columns of b
    const long d = deg(A);

    // check requirement
    if (d <= deg(b))
        LogicError("Bad degrees for linsolve via high order lifting: requires deg(b) < deg(A)");
    // if A is constant, just compute u = A^{-1} b
    else if (d==0)
    {
        Mat<zz_p> invA;
        inv(invA, coeff(A, 0));
        mul(u, invA, b);
        trunc(u, u, prec);
        return;
    }

    // numbers of columns of b
    const long bcols = b.NumCols();
    // number of columns in linearized solution (ceil(prec/d) * bcols)
    const long nbcols = (1 + (prec-1) / d) * bcols;

    // buffers for temporary matrices
    Mat<zz_pX> buf1, buf2;

    // compute buf = A^{-1} mod X^{2d},
    // initialize slice to (A^{-1} mod X^{2d}) div X,
    // and truncate buf to A^{-1} mod X^d
    inv_trunc(buf1, A, 2*d);
    Mat<zz_pX> slice;
    RightShift(slice, buf1, 1);
    trunc(buf1, buf1, d);

    // prepare left multipliers for A and (A^{-1} mod X^d)
    const std::unique_ptr<mat_lzz_pX_lmultiplier> ma = get_lmultiplier(A, d-1);
    const std::unique_ptr<mat_lzz_pX_lmultiplier> minvA = get_lmultiplier(buf1, d-1);

    // sol is initially a copy of b
    Mat<zz_pX> sol(b);

    // next terms to deal with
    Mat<zz_pX> next;

    while (sol.NumCols() < nbcols)
    {
        // compute next = middle product
        transpose(buf1, sol);
        transpose(buf2, slice);
        middle_product(next, buf1, buf2, d-1, d-1); // deg(next) < d
        transpose(next, next);

        ma->multiply(buf1, next);
        trunc(buf2, buf1, d);
        horizontal_join(buf1, sol, buf2);
        sol.swap(buf1); // deg(sol) < d

        if (sol.NumCols() < nbcols)
            high_order_lift_inverse_odd(slice, slice, ma, minvA, d);
    }
    minvA->multiply(buf1, sol);
    trunc(sol, buf1, d);
    collapse_nonconsecutive_columns(u, sol, d, bcols);
    trunc(u, u, prec);
}


/*------------------------------------------------------------*/
/* solve A (u/den) = b                                        */
/* A must be square, A(0) invertible                          */
/* output can alias input                                     */
/* uses lifting and rational reconstruction                   */
/*------------------------------------------------------------*/
void linsolve_via_series(Vec<zz_pX> &u, zz_pX& den, const Mat<zz_pX>& A, const Vec<zz_pX>& b, long nb_max)
{
    // dimensions and degrees
    const long n = A.NumRows(); // == A.NumCols() == b.length()
    const long dA = deg(A);
    const long dB = deg(b);

    // dA is constant: u = A^{-1} b, den = 1
    if (dA == 0)
    {
        Mat<zz_p> iA;
        inv(iA, coeff(A, 0));
        u.SetLength(n);
        // TODO write function for mul(u, iA, b);
        zz_pX buf1, buf2;
        for (long i = 0; i < n; ++i)
        {
            clear(buf1);
            for (long j = 0; j < n; ++j)
            {
                mul(buf2, iA[i][j], b[j]);
                add(buf1, buf1, buf2);
            }
            u[i].swap(buf1);
        }
        set(den);
        return;
    }

    // TODO: need tuning
    if (nb_max == -1)
    {
        nb_max = 1;
        // TODO this is not always correct (see header file todo)
        // just using 1 for the moment
        //switch(type_of_prime())
        //{
        //case TYPE_FFT_PRIME:
        //    if (n >= 100) nb_max = 4;
        //    break;
        //case TYPE_SMALL_PRIME:
        //    if (n >= 100) nb_max = 3;
        //    break;
        //case TYPE_LARGE_PRIME:
        //    if (n >= 100) nb_max = 3;
        //    break;
        //default:
        //    LogicError("Unknown prime type in linear solving via series.");
        //}
    }

    const long nb = min(n, nb_max); // number of linear combinations we are taking
    const long deg_den = n*dA; // expected denominator degree, deg_den + 1 = number of unknowns
    const long deg_num = deg_den-dA + dB; // expected numerator degree

    const long first = max(deg_num + 1, deg_den); // first term we can use in each block
    const long sz = (deg_den + nb - 1) / nb;  // size of each block
    const long prec = first + sz; // precision for solve_series

    Vec<zz_pX> sol_series;
    solve_series(sol_series, A, b, prec);

    if (nb == 1)
    {
        // set lin_comb to random constant linear combination of sol_series
        zz_pX buf, lin_comb;
        for (long j = 0; j < n; ++j)
        {
            trunc(buf, sol_series[j], deg_den+deg_num+1);
            mul(buf, random_zz_p(), buf);
            add(lin_comb, lin_comb, buf);
        }

        // buf = X^(deg_den+deg_num+1)
        clear(buf);
        SetCoeff(buf, deg_den + deg_num + 1);

        // reconstruct denominator
        zz_pXMatrix M;
        HalfGCD(M, buf, lin_comb, deg_den+1);
        // make monic
        zz_p t;
        inv(t, LeadCoeff(M(1,1)));
        mul(den, M(1,1), t);
    }
    else
    {
        // set lin_comb[i] to random constant linear combination of sol_series
        Vec<zz_pX> lin_comb(INIT_SIZE, nb);
        zz_pX buf;
        for (long j = 0; j < n; ++j)
        {
            for (long i = 0; i < nb; ++i)
            {
                mul(buf, random_zz_p(), sol_series[j]);
                add(lin_comb[i], lin_comb[i], buf);
            }
        }

        // build mosaic Toeplitz system for reconstruction of the nb
        // linear combinations
        Vec<Vec<toeplitz_lzz_p>> sys;
        sys.SetLength(nb);
        Vec<zz_p> coeffs(INIT_SIZE, deg_den + sz);
        for (long i = 0; i < nb; ++i)
        {
            for (long j = 0; j < deg_den + sz; ++j)
                GetCoeff(coeffs[j], lin_comb[i], first - deg_den + j);
            sys[i].SetLength(1);
            sys[i][0] = toeplitz_lzz_p(coeffs, sz, deg_den + 1);
        }

        mosaic_toeplitz_lzz_p MT = mosaic_toeplitz_lzz_p(sys);

        // reconstruction of denominator: solve homogeneous mosaic Toeplitz system
        Vec<zz_p> zero(INIT_SIZE, sz*nb);
        Vec<zz_p> ker_vec;
        long ans = MT.solve(ker_vec, zero);
        if (ans == 0)
            Error("Error in solving mosaic toeplitz");

        // den = sum( ker_vec[i] x^i ,  0 <= i <= deg_den )
        VectorCopy(den.rep, ker_vec, deg_den+1);
    }

    // find the numerator u
    zz_pX trunc_den;
    trunc(trunc_den, den, deg_num+1);
    u.SetLength(n);
    for (long i = 0; i < n; ++i)
    {
        trunc(sol_series[i], sol_series[i], deg_num+1);
        MulTrunc(u[i], sol_series[i], trunc_den, deg_num+1);
    }
}

// solve aM = b via kernel basis
// return a and denominator d
// assumes M is invertible
long linsolve_via_kernel(
                         Vec<zz_pX> & u,
                         zz_pX & den,
                         const Mat<zz_pX> & A,
                         const Vec<zz_pX> & b
                        )
{
    // dimensions
    const long m = A.NumRows(); // == b.length()
    const long n = A.NumCols();

    // compute augmented matrix: transpose of A, with b below it
    // ("transpose"/"below" because kernel algos compute *left* kernels) 
    Mat<zz_pX> pmat;
    transpose(pmat, A);
    pmat.SetDims(n+1,m); // add one row below for storing b
    pmat[n] = b;

    // shift = row degree of augmented matrix
    // --> this is the shift which ensures good efficiency for the used kernel
    // algorithm (Zhou-Labahn-Storjohann ISSAC 2012)
    // --> with this shift, we are sure to obtain a vector in the kernel of
    // the form [ -- u -- | den ] with den being the pivot entry (hence,
    // gives an irreducible solution)
    VecLong shift;
    row_degree(shift, pmat);
    VecLong copy_shift(shift);

    // compute kernel
    // TODO via ZLS-approx... for larger dimensions sometimes interp might be better: needs tuning
    Mat<zz_pX> kerbas;
    kernel_basis_zls_via_approximation(kerbas, pmat, shift);

    // compute shifted pivot index of kernel
    // Note that here copy_shift is a copy of the initial shift while `shift`
    // is merely used as a temporary variable which will hold the pivot
    // degrees, which we just ignore (they are not needed)
    VecLong pivind;
    row_pivots(pivind, shift, kerbas, copy_shift);

    // find index of row of kerbas with pivot index in the last entry
    VecLong::const_iterator pivind_n = std::find(pivind.begin(), pivind.end(), n);

    // if no such row, there is no solution
    // --> may happen if kernel is empty, which itself may happen only if m > n
    // --> may also happen if A is rank-deficient but b is not in the
    // K(x)-column space of A, then kerbas is simply a kernel basis for A
    if (pivind_n==pivind.end())
        return 0;

    const long row = pivind_n - pivind.begin();

    // deduce solution
    u.SetLength(n);
    for (long j = 0; j < n; ++j)
        NTL::negate(u[j], kerbas[row][j]);
    den.swap(kerbas[row][n]);
    return 1;
}

/*------------------------------------------------------------*/
/* solve A (u/den) = b                                        */
/* A must be square and nonsingular                           */
/* output can alias input                                     */
/* uses evaluation/interpolation                              */
/*------------------------------------------------------------*/
//long linsolve_via_evaluation(
//                             Vec<zz_pX> & u,
//                             zz_pX & den,
//                             const Mat<zz_pX> & A,
//                             const Vec<zz_pX> & b
//                            )
//{
//    // dimensions and degrees
//    const long n = A.NumRows(); // == A.NumCols() == b.length()
//    const long dA = deg(A);
//    const long dB = deg(b);
//
//    const long deg_den = n*dA; // expected denominator degree, deg_den + 1 = number of unknowns
//    const long deg_num = deg_den-dA + dB; // expected numerator degree
//
//    const long first = max(deg_num + 1, deg_den); // first term we can use in each block
//    const long prec = first + deg_den; // precision for solve_series
//
//    zz_pX_Multipoint_Geometric ev = get_geometric_points(prec);
//    ev.prepare_degree(dA);
//    ev.prepare_degree(dB);
//
//    Vec<Mat<zz_p>> evalsA;
//    ev.evaluate_matrix(evalsA, A);
//
//    Vec<Vec<zz_p>> evalsB;
//    ev.evaluate_vector(evalsB, b);
//
//    zz_p d;
//    for (long i = 0; i < prec; ++i)
//        solve(d, evalsA[i], evalsB[i], evalsB[i]);
//
//    return 0;
//}


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
