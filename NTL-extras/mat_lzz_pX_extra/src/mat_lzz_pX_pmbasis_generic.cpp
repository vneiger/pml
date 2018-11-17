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
    Mat<zz_p> bufR, buf, buf3;
    bufR.SetDims(2*n, n);
    buf.SetDims(n, n);
    buf3.SetDims(n, n);

    for (long k=0; k<d; ++k)
    {
        // *GEN* --> currently, the computed approximant basis has the form
        // [ [X^k I + P00,  P01], [X P10, X^k I + X P11]]
        // where P00, P01, P10 have degree k-1 and P11 has degree k-2
        // It is a 0-ordered weak Popov approximant basis for pmat at order 2*k
        // (For k==0: the last four matrices are in fact zero, and appbas = I)
        // --> residuals R0 and R1 are respectively the coefficients of degree
        // 2*k and 2*k+1 of appbas*pmat

        // 1. compute left kernel of residual 0
        for (long i = 0; i < n; ++i)
            bufR[i].swap(R0_top[i]);
        for (long i = 0; i < n; ++i)
            bufR[n+i].swap(R0_bot[i]);
        kernel(ker,bufR);
        // (GEN) the right n x n submatrix of kerbas is identity
        // --> retrieve the left part
        for (long i = 0; i < n; ++i)
            for (long j = 0; j < n; ++j)
                K0[i][j] = ker[i][j];

        // 2. Update residual 1
        // it is currently   [ [R1_top], [R1_bot] ]
        // it should be [ [R0_top], [K0 * R1_top + R1_bot] ]
        mul(buf, K0, R1_top);
        add(R1_bot, R1_bot, buf);
        R1_top.swap(R0_top); // we do not need the old R1_top anymore, and we won't use R0_top either

        // 3. compute kernel of residual 1
        // we permute the two blocks top-bottom, to respect the (implicit) shift
        for (long i = 0; i < n; ++i)
            bufR[i].swap(R0_bot[i]);
        for (long i = 0; i < n; ++i)
            bufR[n+i].swap(R0_top[i]);
        kernel(ker,bufR);
        // (GEN) the right n x n submatrix of kerbas is identity
        // --> retrieve the left part
        for (long i = 0; i < n; ++i)
            for (long j = 0; j < n; ++j)
                K1[i][j] = ker[i][j];
        // because of the permutation of the residual above, this should
        // in fact be seen as the right part of the kernel
        // --> it will multiply the bottom rows of appbas

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
            add(P10[kk], P10[kk], P00[kk]);
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
            add(P11[kk], P11[kk], P01[kk]);
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

        // 4.2 update by second computed basis [ [I, K1], [0, X I] ]
        // From step 4.1, appbas has the form
        // [ [X^{k+1} I + X P00,  X P01],
        //   [P10, X^k I + P11] ]
        // where P00, P01, P11 have degree k-1 and P10 has degree k
        // --> update it as
        // [ [X^{k+1} I + X P00 + K1 P10,  X P01 + X^k K1 + K1 P11],
        //   [X P10, X^{k+1} I + X P11] ]

        // TODO

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

        // 5. find new residual
        // rescomp: compute the coefficient of degree 2*k+2 of appbas*pmat
        // -- remember that appbas has degree = k and the actual approximant
        // basis is:
        //     top    ==>   appbas_top + [ X^{k+1} id  |   0 ]
        //     bottom ==>   appbas_bottom + [ 0    |   X^{k+1} id ]
        // -- remember that pmat has degree degF, and that we required
        // degF >= d > k

        // 6.1/ compute the top part of residual
        clear(buf3);
        for (long dd=std::max<long>(0,2*k+2-degF); dd<=k; ++dd)
        {
            mul(buf2, coeffs_appbas_top[dd], coeffs_pmat[2*k+2-dd]);
            add(buf3, buf3, buf2);
        }
        // buf3 is almost the top part of residual, except that we have
        // missed dd == k+1, with coeffs_appbas_top[k+1] = [ id | 0 ]
        // --> add top part of coeffs_pmat[k+1];
        for (long i = 0; i < n; ++i)
        {
            residual[i].swap(buf3[i]);
            add(residual[i], residual[i], coeffs_pmat[k+1][i]);
        }

        // 6.2/ compute the bottom part of residual
        clear(buf3);
        for (long dd=std::max<long>(0,2*k+2-degF); dd<k; ++dd)
        {
            mul(buf2, coeffs_appbas_bottom[dd], coeffs_pmat[2*k+2-dd]);
            add(buf3, buf3, buf2);
        }
        // buf3 is almost the bottom part of residual, except that we have
        // missed dd == k+1, with coeffs_appbas_bottom[k] = [ 0 | id ]:
        // --> add bottom part of coeffs_pmat[k+1];
        for (long i = 0; i < n; ++i)
        {
            residual[n+i].swap(buf3[i]);
            add(residual[n+i], residual[n+i], coeffs_pmat[k+1][n+i]);
        }
        //{
        //    std::cout << "----------------------------" << std::endl;
        //    std::cout << "ITER " << k << std::endl;
        //    std::cout << "TOP: " << std::endl;
        //    std::cout << coeffs_appbas_top << std::endl;
        //    std::cout << "BOT: " << std::endl;
        //    std::cout << coeffs_appbas_bottom << std::endl;
        //    // TODO temporary, for testing
        //    // convert to polynomial matrix format
        //    Mat<zz_pX> appbas_top;
        //    conv(appbas_top, coeffs_appbas_top);
        //    Mat<zz_pX> appbas_bottom;
        //    conv(appbas_bottom, coeffs_appbas_bottom);

        //    appbas.SetDims(2*n,2*n);
        //    for (long i = 0; i < n; ++i)
        //        appbas[i] = appbas_top[i];
        //    for (long i = 0; i < n; ++i)
        //        appbas[n+i] = appbas_bottom[i];

        //    // insert identity for the coefficient matrix of degree k
        //    for (long i = 0; i < 2*n; ++i)
        //        SetCoeff(appbas[i][i], k+1);

        //    is_approximant_basis(appbas,pmat,2*(k+1),VecLong(2*n,0),ORD_WEAK_POPOV,true);
        //    std::cout << "----------------------------" << std::endl;
        //}
    }

    // convert to polynomial matrix format
    Mat<zz_pX> appbas_top;
    conv(appbas_top, coeffs_appbas_top);
    Mat<zz_pX> appbas_bottom;
    conv(appbas_bottom, coeffs_appbas_bottom);

    appbas.SetDims(2*n,2*n);
    for (long i = 0; i < n; ++i)
        appbas[i].swap(appbas_top[i]);
    for (long i = 0; i < n; ++i)
        appbas[n+i].swap(appbas_bottom[i]);

    // insert identity for the coefficient matrix of degree d
    for (long i = 0; i < 2*n; ++i)
        SetCoeff(appbas[i][i], d);
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
