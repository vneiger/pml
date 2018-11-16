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
/* MBASIS -- GENERIC INPUT -- UNIFORM SHIFT                   */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* Rescomp version, requiring m = 2 n and order even          */
/*------------------------------------------------------------*/
// TODO try resupdate (m=2n is borderline between the two)
// TODO compare with mbasis_generic for m = t n based on Krylov
// requirement 1: m = 2*n
// requirement 2: order is even
// output: appbas is in 0-Popov form with row degree (d,.., d) *GEN*,
// where d = order/2
void mbasis_generic_2n_n_rescomp(
                                 Mat<zz_pX> & appbas,
                                 const Mat<zz_pX> & pmat,
                                 const long order
                                )
{
    long m = pmat.NumRows();
    long n = pmat.NumCols();
    long d = order/2;
    if (m != 2*n)
        throw std::invalid_argument("~~mbasis_generic_2n_n~~ bad dimensions of pmat");
    if (order % 2 != 0)
        throw std::invalid_argument("~~mbasis_generic_2n_n~~ order must be even");

    // coefficient matrices of input polynomial matrix
    Vec<Mat<zz_p>> coeffs_pmat;
    conv(coeffs_pmat, pmat);
    long deg_pmat = coeffs_pmat.length()-1;
    if (deg_pmat < d)
        throw std::invalid_argument("~~mbasis_generic_2n_n~~ too low degree for input matrix");

    // residual matrix, initially the constant coefficient of pmat
    Mat<zz_p> residual(coeff(pmat,0));

    // coefficient matrices of output approximant basis
    Vec<Mat<zz_p>> coeffs_appbas_top;
    Vec<Mat<zz_p>> coeffs_appbas_bottom;
    coeffs_appbas_top.SetLength(d+1);
    coeffs_appbas_bottom.SetLength(d+1);
    for (long dd = 0; dd < d+1; ++dd)
    {
        coeffs_appbas_top[dd].SetDims(n,m);
        coeffs_appbas_bottom[dd].SetDims(n,m);
    }
    // *GEN* --> appbas is of the form X^d Id + P, for some m x m matrix P
    // which has degree less than d
    // For efficiency reasons, we split appbas into its top n rows and
    // its bottom n rows

    // To store the constant kernel basis and its lefthand square submatrix
    // *GEN* --> dimension of kernel is n x m, its n x n righthand submatrix is
    // the identity
    Mat<zz_p> kerbas, kerbas_left;
    kerbas.SetDims(n, m);
    kerbas_left.SetDims(n, n);

    // buffer, to store products and sums
    Mat<zz_p> buf1, buf2, buf3;
    buf1.SetDims(n, m);
    buf2.SetDims(n, n);
    buf3.SetDims(n, n);

    for (long k=0; k<d; ++k)
    {
        // --> At the beginning of this iteration, appbas + X^k Id
        // is the 0-Popov approximant basis for pmat at order 2k,
        // and residual is the coefficient of degree k of appbas*pmat
        // --> At the end of this iteration, appbas + X^(k+1) Id
        // is the 0-Popov approximant basis for pmat at order 2k
        // and residual is the coefficient of degree k+1 of appbas*pmat

        // 1. compute kernel of residual
        kernel(kerbas,residual);
        // (GEN) the right n x n submatrix of kerbas is identity
        // --> retrieve the left part
        for (long i = 0; i < n; ++i)
            for (long j = 0; j < n; ++j)
                kerbas_left[i][j] = kerbas[i][j];

        // 2. update approximant basis

        // 2.1/ appbas_bottom = kerbas * appbas, that is,
        // appbas_bottom += kerbas_left * appbas_top

        // as we will observe below, the constant coeff dd=0 of bottom is 0
        // (except for k==0, where it is a "hidden" identity)
        if (k>0)
            mul(coeffs_appbas_bottom[0], kerbas_left, coeffs_appbas_top[0]);

        for (long dd = 1; dd < k; ++dd)
        {
            mul(buf1, kerbas_left, coeffs_appbas_top[dd]);
            add(coeffs_appbas_bottom[dd], coeffs_appbas_bottom[dd], buf1);
        }
        // And remember X^k identity in the top-left corner of the approx basis:
        for (long i = 0; i < n; ++i)
            for (long j = 0; j < n; ++j)
                coeffs_appbas_bottom[k][i][j] += kerbas_left[i][j];
        // TODO could be optimized for k=0: simply copy kerbas_left

        // 2.2/ appbas_top *= X
        for (long dd = k-1; dd >= 0; --dd)
            coeffs_appbas_top[dd+1].swap(coeffs_appbas_top[dd]);
        // note: since coeffs_appbas_top[k] = 0 before this loop,
        // this ensures coeffs_appbas_top[0] = 0 at the end of this loop

        // 3. find new residual
        // rescomp: compute the coefficient of degree 2*k+1 of appbas*pmat
        // -- remember that appbas has degree <= k and the actual approximant
        // basis is:
        //     top    ==>   appbas_top + [ X^{k+1} id  |   0 ]
        //     bottom ==>   appbas_bottom + [ 0    |   X^{k} id ]
        // -- remember that pmat has degree deg_pmat, and that we required
        // deg_pmat >= d > k

        // Note that the next kernel will ask us to permute the top and bottom
        // of residual --> directly compute it permuted

        // 3.1/ compute the bottom part of residual (the top part of coeff(appbas*pmat,2k+1))
        clear(buf3);
        for (long dd=std::max<long>(0,2*k+1-deg_pmat); dd<=k; ++dd)
        {
            mul(buf2, coeffs_appbas_top[dd], coeffs_pmat[2*k+1-dd]);
            add(buf3, buf3, buf2);
        }
        // buf3 is almost the bottom part of residual, except that we have
        // missed dd == k+1, with coeffs_appbas_top[k+1] = [ id | 0 ]
        // --> add top part of coeffs_pmat[k];
        for (long i = 0; i < n; ++i)
        {
            residual[n+i].swap(buf3[i]);
            add(residual[n+i], residual[n+i], coeffs_pmat[k][i]);
        }

        // 3.2/ compute the top part of residual (the bottom part of coeff(appbas*pmat,2k+1))
        clear(buf3);
        for (long dd=std::max<long>(0,2*k+1-deg_pmat); dd<k; ++dd)
        {
            mul(buf2, coeffs_appbas_bottom[dd], coeffs_pmat[2*k+1-dd]);
            add(buf3, buf3, buf2);
        }
        // buf3 is almost the top part of residual, except that we have
        // missed dd == k, with coeffs_appbas_bottom[k] = [ 0 | id ]:
        // --> add bottom part of coeffs_pmat[k+1];
        for (long i = 0; i < n; ++i)
        {
            residual[i].swap(buf3[i]);
            add(residual[i], residual[i], coeffs_pmat[k+1][n+i]);
        }

        // 4. compute kernel of residual
        kernel(kerbas,residual);
        // (GEN) the right n x n submatrix of kerbas is identity
        // --> retrieve the left part
        for (long i = 0; i < n; ++i)
            for (long j = 0; j < n; ++j)
                kerbas_left[i][j] = kerbas[i][j];
        // because of the permutation of the residual above, this should
        // in fact be seen as the right part of the kernel
        // --> it will multiply the bottom rows of appbas

        // 5. update approximant basis

        // 5.1/ appbas_top = permuted_kerbas * appbas, that is,
        // appbas_top += kerbas_left * appbas_bottom
        // Remember coeffs_appbas_top[0] == 0:
        mul(coeffs_appbas_top[0], kerbas_left, coeffs_appbas_bottom[0]);
        for (long dd = 1; dd <= k-1; ++dd)
        {
            mul(buf1, kerbas_left, coeffs_appbas_bottom[dd]);
            add(coeffs_appbas_top[dd], coeffs_appbas_top[dd], buf1);
        }
        // And remember X^k identity in the bottom-right corner of appbas:
        for (long i = 0; i < n; ++i)
            for (long j = 0; j < n; ++j)
                coeffs_appbas_top[k][i][n+j] += kerbas_left[i][j];
        // TODO could be optimized for k=0: simple copy

        // 5.2/ appbas_bottom *= X
        for (long dd = k; dd >= 0; --dd)
            coeffs_appbas_bottom[dd+1].swap(coeffs_appbas_bottom[dd]);
        // note: since coeffs_appbas_bottom[k] = 0 before this loop,
        // this ensures coeffs_appbas_bottom[0] = 0 at the end of this loop

        // 6. find new residual
        // rescomp: compute the coefficient of degree 2*k+2 of appbas*pmat
        // -- remember that appbas has degree = k and the actual approximant
        // basis is:
        //     top    ==>   appbas_top + [ X^{k+1} id  |   0 ]
        //     bottom ==>   appbas_bottom + [ 0    |   X^{k+1} id ]
        // -- remember that pmat has degree deg_pmat, and that we required
        // deg_pmat >= d > k

        // 6.1/ compute the top part of residual
        clear(buf3);
        for (long dd=std::max<long>(0,2*k+2-deg_pmat); dd<=k; ++dd)
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
        for (long dd=std::max<long>(0,2*k+2-deg_pmat); dd<k; ++dd)
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
    }

    // convert to polynomial matrix format
    Mat<zz_pX> appbas_top;
    conv(appbas_top, coeffs_appbas_top);
    Mat<zz_pX> appbas_bottom;
    conv(appbas_bottom, coeffs_appbas_bottom);

    appbas.SetDims(m,m);
    for (long i = 0; i < n; ++i)
        appbas[i].swap(appbas_top[i]);
    for (long i = 0; i < n; ++i)
        appbas[n+i].swap(appbas_bottom[i]);

    // insert identity for the coefficient matrix of degree d
    for (long i = 0; i < m; ++i)
        SetCoeff(appbas[i][i], d, 1);
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
