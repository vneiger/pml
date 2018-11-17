#include <NTL/matrix.h>
#include <NTL/mat_lzz_p.h>
#include <NTL/lzz_pX.h>
#include <NTL/BasicThreadPool.h>
#include <cmath>
#include <algorithm> // for manipulating std::vector (min, max, ..)
#include <numeric> // for std::iota

#include "lzz_p_extra.h"
#include "lzz_pX_CRT.h"
#include "mat_lzz_pX_extra.h"

//#define MBASIS_GEN_PROFILE
//#define PMBASIS_GEN_PROFILE
#define VERBOSE_MBASISGEN

NTL_CLIENT

/*------------------------------------------------------------*/
/* when a comment in the code indicates *GEN*, this means     */
/* that this holds thanks to the assumption that the input    */
/* matrix behaves like a generic matrix                       */
/*------------------------------------------------------------*/



/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* UTIL FUNCTIONS                                             */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* convert to / from Vec<Mat<zz_p>>                           */
/* (degree deduced from input)                                */
/* split into first half rows, second half rows               */
/*------------------------------------------------------------*/
void conv_top_bot(
                  Vec<Mat<zz_p>> & coeffs_top,
                  Vec<Mat<zz_p>> & coeffs_bot,
                  const Mat<zz_pX> & mat
                 )
{
    long d = deg(mat);
    long m1 = mat.NumRows()/2;
    long m2 = mat.NumRows()-m1;
    long n = mat.NumCols();

    coeffs_top.SetLength(d+1);
    for (long k = 0; k <= d; ++k)
    {
        coeffs_top[k].SetDims(m1, n);
        for (long i = 0; i < m1; ++i)
            for (long j = 0; j < n; ++j)
                coeffs_top[k][i][j] = coeff(mat[i][j], k);
    }

    coeffs_bot.SetLength(d+1);
    for (long k = 0; k <= d; ++k)
    {
        coeffs_bot[k].SetDims(m2, n);
        for (long i = 0; i < m2; ++i)
            for (long j = 0; j < n; ++j)
                coeffs_bot[k][i][j] = coeff(mat[m1+i][j], k);
    }
}

void conv_top_bot(
                  Mat<zz_pX> & mat,
                  const Vec<Mat<zz_p>> & coeffs_top,
                  const Vec<Mat<zz_p>> & coeffs_bot
                 )
{
    long len = coeffs_top.length();
    // requires coeffs_bot has same length
    if (len == 0)
    {
        // TODO this is a choice --> indicate it in comments
        // (zero-length sequence could be zero matrix of any dimension)
        mat.SetDims(0, 0);
        return;
    }
    long n = coeffs_top[0].NumCols();
    // requires coeffs_bot has same numcols, and all other coeffs_top/bot as well
    long m1 = coeffs_top[0].NumRows();
    long m2 = coeffs_top[0].NumRows();

    mat.SetDims(m1+m2, n);

    for (long i = 0; i < m1; ++i)
        for (long j = 0; j < n; ++j)
            for (long k = 0; k < len; ++k)
                SetCoeff(mat[i][j], k, coeffs_top[k][i][j]);

    for (long i = 0; i < m2; ++i)
        for (long j = 0; j < n; ++j)
            for (long k = 0; k < len; ++k)
                SetCoeff(mat[m1+i][j], k, coeffs_bot[k][i][j]);
}



/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* MBASIS -- GENERIC INPUT -- UNIFORM SHIFT                   */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* Rescomp version, requiring m = 2 n and order even          */
/*------------------------------------------------------------*/
// TODO try resupdate (m=2n is borderline between the two)
// TODO compare with mbasis_generic for m = t n based on Krylov
// requirement 1: m = 2*n
// requirement 2: order is even and strictly positive
// output: appbas is in 0-Popov form with row degree (d,.., d) *GEN*,
// where d = order/2
void mbasis_generic_2n_n_rescomp(
                                 Mat<zz_pX> & appbas,
                                 const Mat<zz_pX> & pmat,
                                 const long order
                                )
{
    long n = pmat.NumCols();
    long d = order/2;
    if (pmat.NumRows() != 2*n)
        throw std::invalid_argument("~~mbasis_generic_2n_n~~ bad dimensions of pmat");
    if (order==0 || order % 2 != 0)
        throw std::invalid_argument("~~mbasis_generic_2n_n~~ order must be even");

    // coefficient matrices of input polynomial matrix
    // pmat = [[Ft], [Fb]]
    Vec<Mat<zz_p>> F_top, F_bot;
    conv_top_bot(F_top, F_bot, pmat);
    long degF = F_top.length()-1;
    if (degF < d)
        throw std::invalid_argument("~~mbasis_generic_2n_n~~ too low degree for input matrix");

    // residual matrix 0, initially coefficient of pmat of degree 0
    Mat<zz_p> R0_top(F_top[0]);
    Mat<zz_p> R0_bot(F_bot[0]);
    // residual matrix 1, initially coefficient of pmat of degree 1
    Mat<zz_p> R1_top(F_top[1]);
    Mat<zz_p> R1_bot(F_bot[1]);

    // coefficient matrices of output approximant basis
    // *GEN* --> appbas will be computed in the form
    // [ [X^d I + P00,  P01], [X P10, X^d I + X P11]]
    // where P00, P01, P10 have degree d-1 and P11 has degree d-2
    Vec<Mat<zz_p>> P00, P01, P10, P11;
    P00.SetLength(d);
    for (long k = 0; k < d; ++k)
        P00[k].SetDims(n,n);
    P01.SetLength(d);
    for (long k = 0; k < d; ++k)
        P01[k].SetDims(n,n);
    P10.SetLength(d+1); // we still store the degree 0 coeff (which will eventually be zero)
    for (long k = 0; k < d+1; ++k)
        P01[k].SetDims(n,n);
    P11.SetLength(d); // we still store the degree 0 coeff (which will eventually be zero)
    for (long k = 0; k < d; ++k)
        P11[k].SetDims(n,n);

    // To store the kernel of the residuals 0 and 1, and their lefthand square
    // submatrices
    // *GEN* --> dimension of kernel is n x m, its n x n righthand submatrix is
    // the identity
    Mat<zz_p> ker, K0, K1;
    ker.SetDims(n, 2*n);
    K0.SetDims(n, n);
    K1.SetDims(n, n);

    // buffer, to store products and sums
    Mat<zz_p> bufR, buf;
    bufR.SetDims(2*n, n);
    buf.SetDims(n, n);

    for (long k=0; k<d; ++k)
    {
#ifdef VERBOSE_MBASISGEN
        std::cout << "===========" << std::endl;
        std::cout << "Iter " << k << ", start" << std::endl;
#endif // VERBOSE_MBASISGEN

#ifdef VERBOSE_MBASISGEN
        // convert to polynomial matrix format

        // appbas = [ [P00,  P01], [X P10, X P11]]
        // where P00, P01, P10 have degree k-1 and P11 has degree k-2
        appbas.SetDims(2*n,2*n);
        clear(appbas);
        for (long i = 0; i < n; ++i)
        {
            for (long j = 0; j < n; ++j)
                for (long kk = 0; kk < k; ++kk)
                    SetCoeff(appbas[i][j], kk, P00[kk][i][j]);
            for (long j = n; j < 2*n; ++j)
                for (long kk = 0; kk < k; ++kk)
                    SetCoeff(appbas[i][j], kk, P01[kk][i][j-n]);
        }
        for (long i = n; i < 2*n; ++i)
        {
            for (long j = 0; j < n; ++j)
                for (long kk = 1; kk < k+1; ++kk)
                    SetCoeff(appbas[i][j], kk, P10[kk][i-n][j]);
            for (long j = n; j < 2*n; ++j)
                for (long kk = 1; kk < k; ++kk)
                    SetCoeff(appbas[i][j], kk, P11[kk][i-n][j-n]);
        }

        // add X^d I
        for (long i = 0; i < 2*n; ++i)
            SetCoeff(appbas[i][i], k);
        std::cout << "degree matrix of appbas: " << std::endl << degree_matrix(appbas) << std::endl;
        is_approximant_basis(appbas, pmat, k, VecLong(2*n,0), ORD_WEAK_POPOV, true);
#endif // VERBOSE_MBASISGEN


        // *GEN* --> currently, the computed approximant basis has the form
        // [ [X^k I + P00,  P01], [X P10, X^k I + X P11]]
        // where P00, P01, P10 have degree k-1 and P11 has degree k-2
        // It is a 0-ordered weak Popov approximant basis for pmat at order 2*k
        // (For k==0: the last four matrices are in fact zero, and appbas = I)
        // --> residuals R0 and R1 are respectively the coefficients of degree
        // 2*k and 2*k+1 of appbas*pmat

        // 1. compute left kernel of residual 0
        for (long i = 0; i < n; ++i)
            bufR[i] = R0_top[i]; // TODO VectorCopy ??
        for (long i = 0; i < n; ++i)
            bufR[n+i] = R0_bot[i]; // TODO VectorCopy ??
        kernel(ker,bufR);
        // (GEN) the right n x n submatrix of kerbas is identity
        // --> retrieve the left part
        for (long i = 0; i < n; ++i)
            for (long j = 0; j < n; ++j)
                K0[i][j] = ker[i][j];

#ifdef VERBOSE_MBASISGEN
        std::cout << "...Ker0 computed " << std::endl;
#endif // VERBOSE_MBASISGEN

        // 2. Update residual 1
        // it is currently   [ [R1_top], [R1_bot] ]
        // it should be [ [R0_top], [K0 * R1_top + R1_bot] ]
        mul(buf, K0, R1_top);
        add(R1_bot, R1_bot, buf);
        //R1_top.swap(R0_top); // we do not need the old R1_top anymore, and we won't use R0_top either
        R1_top = R0_top; // TODO swap

#ifdef VERBOSE_MBASISGEN
        std::cout << "...R1 updated " << std::endl;
#endif // VERBOSE_MBASISGEN

        // 3. compute kernel of residual 1
        // we permute the two blocks top-bottom, to respect the (implicit) shift
        for (long i = 0; i < n; ++i)
            //bufR[i].swap(R1_bot[i]);
            bufR[i] = R1_bot[i];
        for (long i = 0; i < n; ++i)
            //bufR[n+i].swap(R1_top[i]);
            bufR[n+i] = R1_top[i];
        kernel(ker,bufR);
        // (GEN) the right n x n submatrix of kerbas is identity
        // --> retrieve the left part
        for (long i = 0; i < n; ++i)
            for (long j = 0; j < n; ++j)
                K1[i][j] = ker[i][j];

#ifdef VERBOSE_MBASISGEN
        std::cout << "...Ker1 computed " << std::endl;
#endif // VERBOSE_MBASISGEN

        // 4. update approximant basis

        // 4.1 update by first computed basis [ [XI, 0], [K0, I] ]
        // Recall: currently, appbas has the form
        // [ [X^k I + P00,  P01], [X P10, X^k I + X P11]]
        // where P00, P01, P10 have degree k-1 and P11 has degree k-2
        // (negative k-1 or k-2 means zero matrix)
        // --> update it as
        // [ [X^{k+1} I + X P00,  X P01],
        //   [K0 P00 + X^k K0 + X P10, K0 P01 + X^k I + X P11] ]

        // bottom left
        // add constant term of K0 P00
        mul(P10[0], K0, P00[0]); 
        // add terms of degree 1 ... k-1 of K0 P00
        for (long kk = 1; kk < k; ++kk)
        {
            mul(buf, K0, P00[kk]);
            add(P10[kk], P10[kk], buf);
        }
        add(P10[k], P10[k], K0); // add K0 to coeff of degree k
        // TODO could be optimized when k==0: no add, just assign

        // bottom right
        // add constant term of K0 P01
        mul(P11[0], K0, P01[0]); 
        // add terms of degree 1 ... k-1 of K0 P01
        for (long kk = 1; kk < k; ++kk)
        {
            mul(buf, K0, P01[kk]);
            add(P11[kk], P11[kk], buf);
        }
        // TODO could be optimized when k==1: no add, just assign, because P11[kk]==0

        // top left: P00 <- X P00  (P00 has degree k-1)
        for (long kk = k-1; kk >=0; --kk)
            P00[kk+1].swap(P00[kk]);
        // since the X^k I was actually not stored in P00[k], that matrix was
        // zero and therefore this loop does put the zero matrix in P00[0]

        // top right: P01 <- X P01  (P01 has degree k-1)
        for (long kk = k-1; kk >=0; --kk)
            P01[kk+1].swap(P01[kk]);
        // note: this loop does put the zero matrix in P01[0]

#ifdef VERBOSE_MBASISGEN
        std::cout << "...appbas, first update done" << std::endl;
#endif // VERBOSE_MBASISGEN


        // 4.2 update by second computed basis [ [I, K1], [0, X I] ]
        // From step 4.1, appbas has the form
        // [ [X^{k+1} I + X P00,  X P01],
        //   [P10, X^k I + P11] ]
        // where P00, P01, P11 have degree k-1 and P10 has degree k
        // --> update it as
        // [ [X^{k+1} I + X P00 + K1 P10,  X P01 + X^k K1 + K1 P11],
        //   [X P10, X^{k+1} I + X P11] ]

        // top left
        // add constant term of K1 P10
        mul(P00[0], K1, P10[0]); 
        // add terms of degree 1 ... k of K1 P10
        for (long kk = 1; kk <= k; ++kk)
        {
            mul(buf, K1, P10[kk]);
            add(P00[kk], P00[kk], buf);
        }

        // top right
        // add constant term of K1 P11
        mul(P00[0], K1, P10[0]); 
        // add terms of degree 1 ... k-1 of K1 P11
        for (long kk = 1; kk < k; ++kk)
        {
            mul(buf, K1, P11[kk]);
            add(P01[kk], P00[kk], buf);
        }
        // add X^k K1
        add(P01[k], P01[k], K1);

        // bottom left: P10 <- X P10  (P10 has degree k)
        for (long kk = k; kk >=0; --kk)
            P10[kk+1].swap(P10[kk]);
        // note: this loop does put the zero matrix in P10[0]

        // bottom right: P11 <- X P11  (P11 has degree k-1)
        for (long kk = k-1; kk >=0; --kk)
            P11[kk+1].swap(P11[kk]);
        // note: since the X^k I was not stored in P11, this loop does put the
        // zero matrix in P11[0]

        // --> now appbas is a 0-ordered weak Popov approximant basis of pmat
        // at order 2*k+2, of the form
        // [ [X^{k+1} I + P00,  P01], [X P10, X^{k+1} I + X P11]]
        // where P00, P01, P10 have degree k and P11 has degree k-1

#ifdef VERBOSE_MBASISGEN
        std::cout << "...appbas, second update done" << std::endl;
#endif // VERBOSE_MBASISGEN

        // 5. if not finished (k<d-1), compute the next residuals (rescomp variant)
        // --> residuals R0 and R1 must be, respectively, the coefficients of
        // degree 2k+2 and 2*k+3 of appbas*pmat
        // TODO improve to deal with F having degree smaller than 2d-1??
        if (k<d-1)
        {
            // R0_top = coeff of degree 2k+2 of (X^{k+1} I + P00) F_top + P01 F_bot,
            // where P00 and P01 have degree k
            mul(R0_top, P00[0], F_top[2*k+2]);
            for (long kk = 1; kk < k+1; ++kk)
            {
                mul(buf, P00[kk], F_top[2*k+2-kk]);
                add(R0_top, R0_top, buf);
            }
            add(R0_top, R0_top, F_top[k+1]);
            for (long kk = 0; kk < k+1; ++kk)
            {
                mul(buf, P01[kk], F_bot[2*k+2-kk]);
                add(R0_top, R0_top, buf);
            }

            // R1_top = coeff of degree 2k+3 of (X^{k+1} I + P00) F_top + P01 F_bot,
            // where P00 and P01 have degree k
            mul(R1_top, P00[0], F_top[2*k+3]);
            for (long kk = 1; kk < k+1; ++kk)
            {
                mul(buf, P00[kk], F_top[2*k+3-kk]);
                add(R1_top, R1_top, buf);
            }
            add(R1_top, R1_top, F_top[k+2]);
            for (long kk = 0; kk < k+1; ++kk)
            {
                mul(buf, P01[kk], F_bot[2*k+3-kk]);
                add(R1_top, R1_top, buf);
            }

            // R0_bot = coeff of degree 2k+2 of X P10 F_top + (X^{k+1} I + X P11) F_bot,
            // where P10 has degree k and P11 has degree k-1
            mul(R0_bot, P10[1], F_top[2*k+1]);
            for (long kk = 2; kk < k+2; ++kk)
            {
                mul(buf, P10[kk], F_top[2*k+2-kk]);
                add(R0_bot, R0_bot, buf);
            }
            for (long kk = 1; kk < k+1; ++kk)
            {
                mul(buf, P11[kk], F_bot[2*k+2-kk]);
                add(R0_bot, R0_bot, buf);
            }
            add(R0_bot, R0_bot, F_bot[k+1]);

            // R1_bot = coeff of degree 2k+3 of X P10 F_top + (X^{k+1} I + X P11) F_bot,
            // where P10 has degree k and P11 has degree k-1
            mul(R1_bot, P10[1], F_top[2*k+2]);
            for (long kk = 2; kk < k+2; ++kk)
            {
                mul(buf, P10[kk], F_top[2*k+3-kk]);
                add(R1_bot, R1_bot, buf);
            }
            for (long kk = 1; kk < k+1; ++kk)
            {
                mul(buf, P11[kk], F_bot[2*k+3-kk]);
                add(R1_bot, R1_bot, buf);
            }
            add(R1_bot, R1_bot, F_bot[k+2]);
#ifdef VERBOSE_MBASISGEN
        std::cout << "... residuals updated" << std::endl;
#endif // VERBOSE_MBASISGEN
        }
#ifdef VERBOSE_MBASISGEN
        else
            std::cout << "... last iteration: no residual update needed" << std::endl;
#endif // VERBOSE_MBASISGEN
    }

    // convert to polynomial matrix format

    // appbas = [ [P00,  P01], [X P10, X P11]]
    // where P00, P01, P10 have degree d-1 and P11 has degree d-2
    appbas.SetDims(2*n,2*n);
    clear(appbas); // TODO remove
    for (long i = 0; i < n; ++i)
    {
        for (long j = 0; j < n; ++j)
            for (long k = 0; k < d; ++k)
                SetCoeff(appbas[i][j], k, P00[k][i][j]);
        for (long j = n; j < 2*n; ++j)
            for (long k = 0; k < d; ++k)
                SetCoeff(appbas[i][j], k, P01[k][i][j-n]);
    }
    for (long i = n; i < 2*n; ++i)
    {
        for (long j = 0; j < n; ++j)
            for (long k = 1; k < d+1; ++k)
                SetCoeff(appbas[i][j], k, P10[k][i-n][j]);
        for (long j = n; j < 2*n; ++j)
            for (long k = 1; k < d; ++k)
                SetCoeff(appbas[i][j], k, P11[k][i-n][j-n]);
    }

    // add X^d I
    for (long i = 0; i < 2*n; ++i)
        SetCoeff(appbas[i][i], d);

#ifdef VERBOSE_MBASISGEN
    is_approximant_basis(appbas, pmat, order, VecLong(2*n,0), ORD_WEAK_POPOV, true);
#endif // VERBOSE_MBASISGEN
}


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* PMBASIS-GENERIC, ONE COLUMN                                */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

// TODO prototype for the one column case
VecLong mbasis_generic_onecolumn(
                                 Mat<zz_pX> & appbas,
                                 const Mat<zz_pX> & pmat,
                                 const long order,
                                 const VecLong & shift
                                )
{
    // TODO just for avoiding the unused error
    VecLong d = shift;
    // TODO current code specific to n=1 !!
    long nrows = pmat.NumRows();
    // partially linearize pmat into order/nrows constant matrices
    Vec<Mat<zz_p>> residuals;
    residuals.SetLength(order/nrows);
    for (long k = 0; k < order/nrows; ++k)
    {
        residuals[k].SetDims(nrows,nrows);
        for (long i = 0; i < nrows; ++i)
            for (long j = 0; j < nrows; ++j)
                residuals[k][i][j] = pmat[i][0][k*nrows+j];
    }

    // for storing appbas coefficients
    Vec<Mat<zz_p>> coeffs_appbas;
    coeffs_appbas.SetLength(1+order/nrows);
    // Warning: the largest degree coefficient of appbas is identity matrix,
    // not stored here (for efficiency in computations below)

    // To store the Krylov matrix
    Mat<zz_p> krylov;
    krylov.SetDims(2*nrows, nrows);

    // To store the kernel basis and its leftmost submatrix
    Mat<zz_p> kerbas, kerbas_left;
    kerbas_left.SetDims(nrows, nrows);

    // To store local residuals (XI * res[k])
    Mat<zz_p> residual;
    residual.SetDims(nrows,nrows);

    for (long d=0; d<order/nrows; ++d)
    {
        // 1. Build the Krylov matrix
        // copy residuals[d] in the top rows of krylov
        for (long i = 0; i < nrows; ++i)
            for (long j = 0; j < nrows; ++j)
                krylov[i][j] = residuals[d][i][j];
        // in bottom rows, copy residuals[d] shifted to the right by one column
        for (long i = 0; i < nrows; ++i)
            for (long j = 1; j < nrows; ++j)
                krylov[i+nrows][j] = residuals[d][i][j-1];

        // 2. compute kernel
        kernel(kerbas,krylov);
        for (long i = 0; i < nrows; ++i)
            for (long j = 0; j < nrows; ++j)
                kerbas_left[i][j] = kerbas[i][j];

        // 3. update approximant basis
        // top coefficient is given by kerbas (remember the hidden identity leading matrix)
        coeffs_appbas[d] = kerbas_left;
        // multiply by XId + kerbas-left
        for (long k = d; k > 0; --k)
        {
            coeffs_appbas[k] += coeffs_appbas[k-1];
            coeffs_appbas[k-1] = kerbas_left * coeffs_appbas[k-1];
        }

        // 4. update residual
        // multiply by XId + kerbas-left
        // -> ignore what happens beyond order
        // -> ignore what happens before d and at d
        for (long k = order/nrows-1; k > d; --k)
        {
            for (long i = 0; i < nrows; ++i)
            {
                residual[i][0] = residuals[k-1][i][nrows-1];
                for (long j = 1; j < nrows; ++j)
                    residual[i][j] = residuals[k][i][j-1];
            }
            residuals[k] = kerbas_left * residuals[k];
            residuals[k] += residual;
        }
    }

    // finally insert identity for the largest degree coefficient
    ident(coeffs_appbas[order/nrows], nrows);
    // and convert to polynomial matrix format
    appbas = conv(coeffs_appbas);

    VecLong dv (pmat.NumRows(), order/nrows);
    return dv;
}



/*------------------------------------------------------------*/
/* Divide and Conquer: PMBasis                                */
/*------------------------------------------------------------*/
VecLong pmbasis_generic_onecolumn(
                       Mat<zz_pX> & appbas,
                       const Mat<zz_pX> & pmat,
                       const long order,
                       const VecLong & shift
                      )
{
    if (order <= 16*pmat.NumRows()){
        //cout << "pmat: " << pmat << endl;
        //cout << "blah blah" << endl;
        //popov_mbasis1_generic2(appbas,pmat,order,shift);
        return mbasis_generic_onecolumn(appbas,pmat,order,shift);
        //cout << "appbas1: " << appbas << endl;
        //auto t = mbasis(appbas,pmat,order,shift);
        //cout << "appbas2: " << appbas << endl;
        //return t;
    }

    VecLong pivdeg; // pivot degree, first call
    VecLong pivdeg2; // pivot degree, second call
    VecLong rdeg(pmat.NumRows()); // shifted row degree
    long order1 = order>>1; // order of first call
    long order2 = order-order1; // order of second call
    Mat<zz_pX> trunc_pmat; // truncated pmat for first call
    Mat<zz_pX> appbas2; // basis for second call
    Mat<zz_pX> residual; // for the residual

    // first recursive call, with 'pmat' and 'shift'
    trunc(trunc_pmat,pmat,order1);
    pivdeg = pmbasis_generic_onecolumn(appbas,trunc_pmat,order1,shift);

    // shifted row degree = shift for second call = pivdeg+shift
    std::transform(pivdeg.begin(), pivdeg.end(), shift.begin(), rdeg.begin(), std::plus<long>());

    // residual = (appbas * pmat * X^-order1) mod X^order2
    long deg_sp = (pmat.NumCols() * order)/ (2*pmat.NumRows());
    right_parlin_middle_product(residual, appbas, pmat, deg_sp, order-1, order1, order2-1);

    // second recursive call, with 'residual' and 'rdeg'
    pivdeg2 = pmbasis_generic_onecolumn(appbas2,residual,order2,rdeg);

    // final basis = appbas2 * appbas
    multiply(appbas,appbas2,appbas);

    // final pivot degree = pivdeg1+pivdeg2
    std::transform(pivdeg.begin(), pivdeg.end(), pivdeg2.begin(), pivdeg.begin(), std::plus<long>());

    return pivdeg;
}


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
