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
