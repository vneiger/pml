#include "mat_lzz_pX_linearization.h"
#include "mat_lzz_pX_approximant.h"

//#define MBASIS_GEN_PROFILE
//#define PMBASIS_GEN_PROFILE
//#define VERBOSE_MBASISGEN

NTL_CLIENT


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
static void conv_top_bot(
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

/*------------------------------------------------------------*/
/* convert to / from Vec<Mat<zz_p>>                           */
/* (degree bound d given by user; Vecs will have length d)    */
/* split into first half rows, second half rows               */
/*------------------------------------------------------------*/
static void conv_top_bot(
                  Vec<Mat<zz_p>> & coeffs_top,
                  Vec<Mat<zz_p>> & coeffs_bot,
                  const Mat<zz_pX> & mat,
                  const long d
                 )
{
    long m1 = mat.NumRows()/2;
    long m2 = mat.NumRows()-m1;
    long n = mat.NumCols();

    coeffs_top.SetLength(d);
    for (long k = 0; k < d; ++k)
    {
        coeffs_top[k].SetDims(m1, n);
        for (long i = 0; i < m1; ++i)
            for (long j = 0; j < n; ++j)
                coeffs_top[k][i][j] = coeff(mat[i][j], k);
    }

    coeffs_bot.SetLength(d);
    for (long k = 0; k < d; ++k)
    {
        coeffs_bot[k].SetDims(m2, n);
        for (long i = 0; i < m2; ++i)
            for (long j = 0; j < n; ++j)
                coeffs_bot[k][i][j] = coeff(mat[m1+i][j], k);
    }
}


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* MBASIS -- GENERIC INPUT -- UNIFORM SHIFT                   */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* Iterative, rescomp version                                 */
/*------------------------------------------------------------*/
// Requirement: m = 2*n ; pmat generic ; order >= 2
// if order = 2d, then 
//      appbas = [ [X^d I + P00,  P01], [X P10, X^d I + X P11]]
// where P00, P01, P10 have degree d-1 and P11 has degree d-2
// if order = 2d+1, then
//      appbas = [ [X^{d+1} I + X P00,  X P01], [P10, X^d I + P11] ]
// where P00, P01, P11 have degree d-1 and P10 has degree d

// TODO clean
// TODO better handle pmat of degree << order ?
void mbasis_generic_2n_n_rescomp(
                                 Mat<zz_pX> & appbas,
                                 const Mat<zz_pX> & pmat,
                                 const long order
                                )
{
#ifdef MBASIS_GEN_PROFILE
    double t_ker=0.0, t_res=0.0, t_app=0.0, t_others=0.0;
    double tt = GetWallTime();
#endif // MBASIS_GEN_PROFILE

    long n = pmat.NumCols();
    long d = order/2;
    bool odd_order = (order%2);

    // coefficient matrices of input polynomial matrix
    // pmat = [[Ft], [Fb]]
    Vec<Mat<zz_p>> F_top, F_bot;
    conv_top_bot(F_top, F_bot, pmat);
    //long degF = F_top.length()-1; // TODO handle the case where F has degree lower (probably still >= d?)

    // residual matrix 0
    // will hold coefficients of pmat of degree 2*k
    Mat<zz_p> R0_top, R0_bot;
    // residual matrix 1
    // will hold coefficients of pmat of degree 2*k+1
    Mat<zz_p> R1_top, R1_bot;

    // coefficient matrices of output approximant basis
    // *GEN* --> appbas will be computed in the form
    // [ [X^d I + P00,  P01], [X P10, X^d I + X P11]]
    // where P00, P01, P10 have degree d-1 and P11 has degree d-2
    Vec<Mat<zz_p>> P00, P01, P10, P11;
    if (odd_order)
    {
        // appbas = [ [X^{d+1} I + X P00,  X P01], [P10, X^d I + P11] ]
        // where P00, P01, P11 have degree d-1 and P10 has degree d
        P00.SetLength(d+1); // we still store the degree 0 coeff (which will eventually be zero)
        P01.SetLength(d+1); // we still store the degree 0 coeff (which will eventually be zero)
        P10.SetLength(d+1);
        P11.SetLength(d);
    }
    else
    {
        // appbas = [ [X^d I + P00,  P01], [X P10, X^d I + X P11]]
        // where P00, P01, P10 have degree d-1 and P11 has degree d-2
        P00.SetLength(d);
        P01.SetLength(d);
        P10.SetLength(d+1); // we still store the degree 0 coeff (which will eventually be zero)
        P11.SetLength(d); // we still store the degree 0 coeff (which will eventually be zero)
    }

    // To store the kernel of the residuals 0 and 1, and their lefthand square
    // submatrices
    // *GEN* --> dimension of kernel is n x m, its n x n righthand submatrix is
    // the identity
    Mat<zz_p> ker, K0, K1;

    // buffer, to store products and sums
    Mat<zz_p> bufR, buf;

#ifdef MBASIS_GEN_PROFILE
    t_others += GetWallTime()-tt;
#endif // MBASIS_GEN_PROFILE

    // Iteration k==0:
    // --> compute approximant basis at order 2
    // --> update residuals R0 and R1 as the coefficients of degree 2*k and
    // 2*k+1 of appbas*pmat

#ifdef MBASIS_GEN_PROFILE
    tt = GetWallTime();
#endif // MBASIS_GEN_PROFILE
    // 1. compute left kernel of coeff 0 of pmat
    bufR.SetDims(2*n, n);
    for (long i = 0; i < n; ++i)
        for (long j = 0; j < n; ++j)
            bufR[i][j] = F_top[0][i][j];
    for (long i = 0; i < n; ++i)
        for (long j = 0; j < n; ++j)
            bufR[n+i][j] = F_bot[0][i][j];
#ifdef MBASIS_GEN_PROFILE
    t_others += GetWallTime()-tt;
    tt = GetWallTime();
#endif // MBASIS_GEN_PROFILE
    kernel(ker,bufR);
    // (GEN) the right n x n submatrix of kerbas is identity
    // --> retrieve the left part
    K0.SetDims(n, n);
    for (long i = 0; i < n; ++i)
        for (long j = 0; j < n; ++j)
            K0[i][j] = ker[i][j];
#ifdef MBASIS_GEN_PROFILE
    t_ker += GetWallTime()-tt;
    tt = GetWallTime();
#endif // MBASIS_GEN_PROFILE

    // 2. and 3. Second kernel, of coeff of degree 1 of appbas * pmat, that
    // is, [ [F_top[0]], [K0 * F_top[1] + F_bot[1]] ]
    // we permute the two blocks top-bottom, to respect the (implicit) shift
    mul(buf, K0, F_top[1]);
    add(buf, buf, F_bot[1]);
#ifdef MBASIS_GEN_PROFILE
    t_res += GetWallTime()-tt;
    tt = GetWallTime();
#endif // MBASIS_GEN_PROFILE
    for (long i = 0; i < n; ++i) // bufR bottom <- F_top[0]
        bufR[n+i].swap(bufR[i]);
    for (long i = 0; i < n; ++i) // bufR top <- K0*F_top[1] + F_bot[1]
        bufR[i].swap(buf[i]);
#ifdef MBASIS_GEN_PROFILE
    t_others += GetWallTime()-tt;
    tt = GetWallTime();
#endif // MBASIS_GEN_PROFILE
    kernel(ker,bufR);
    // (GEN) the right n x n submatrix of kerbas is identity
    // --> retrieve the left part
    K1.SetDims(n, n);
    for (long i = 0; i < n; ++i)
        for (long j = 0; j < n; ++j)
            K1[i][j] = ker[i][j];
#ifdef MBASIS_GEN_PROFILE
    t_ker += GetWallTime()-tt;
    tt = GetWallTime();
#endif // MBASIS_GEN_PROFILE

    // 4. update approximant basis

    // 4.1 update by first computed basis [ [XI, 0], [K0, I] ]
    // 4.2 update by second computed basis [ [I, K1], [0, X I] ]
    // appbas is currently the identity matrix
    // --> update it as
    // [ [X I + K1 K0,  K1],
    //   [X K0,        X I] ]
    //   (remember the identity terms are not stored)
    mul(P00[0], K1, K0); 
    P01[0] = K1;
    P10[1] = K0;

    // --> now appbas is a 0-ordered weak Popov approximant basis of pmat
    // at order 2, of the form
    // [ [X I + P00,  P01], [X P10, X I + X P11]]
    // where P00, P01, P10 have degree 0 and P11 == 0

#ifdef MBASIS_GEN_PROFILE
    t_app += GetWallTime()-tt;
    tt = GetWallTime();
#endif // MBASIS_GEN_PROFILE

    // 5. if not finished (d>1), compute the next residuals (rescomp variant)
    // --> residuals R0 and R1 must be, respectively, the coefficients of
    // degree 2 and 3 of appbas*pmat
    if (d>1)
    {
        // R0_top = coeff of degree 2 of (X I + P00) F_top + K1 F_bot,
        // where P00 has degree 0
        mul(R0_top, P00[0], F_top[2]);
        add(R0_top, R0_top, F_top[1]);
        mul(buf, K1, F_bot[2]);
        add(R0_top, R0_top, buf);

        // R1_top = coeff of degree 3 of (X I + P00) F_top + K1 F_bot,
        // where P00 has degree 0
        mul(R1_top, P00[0], F_top[3]);
        add(R1_top, R1_top, F_top[2]);
        mul(buf, K1, F_bot[3]);
        add(R1_top, R1_top, buf);

        // R0_bot = coeff of degree 2 of X K0 F_top + X I F_bot,
        mul(R0_bot, K0, F_top[1]);
        add(R0_bot, R0_bot, F_bot[1]);

        // R1_bot = coeff of degree 3 of X K0 F_top + X I F_bot,
        mul(R1_bot, K0, F_top[2]);
        add(R1_bot, R1_bot, F_bot[2]);

#ifdef MBASIS_GEN_PROFILE
        t_res += GetWallTime()-tt;
#endif // MBASIS_GEN_PROFILE
    }

    // Iterations k==1...d-1
    for (long k=1; k<d; ++k)
    {
        // *GEN* --> currently, the computed approximant basis has the form
        // [ [X^k I + P00,  P01], [X P10, X^k I + X P11]]
        // where P00, P01, P10 have degree k-1 and P11 has degree k-2
        // It is a 0-ordered weak Popov approximant basis for pmat at order 2*k
        // (For k==0: the last four matrices are in fact zero, and appbas = I)
        // --> residuals R0 and R1 are respectively the coefficients of degree
        // 2*k and 2*k+1 of appbas*pmat

#ifdef MBASIS_GEN_PROFILE
        tt = GetWallTime();
#endif // MBASIS_GEN_PROFILE
        // 1. compute left kernel of residual 0
        for (long i = 0; i < n; ++i)
            for (long j=0; j<n; ++j)
                bufR[i][j] = R0_top[i][j];
        for (long i = 0; i < n; ++i)
            bufR[n+i].swap(R0_bot[i]);
#ifdef MBASIS_GEN_PROFILE
        t_others += GetWallTime()-tt;
        tt = GetWallTime();
#endif // MBASIS_GEN_PROFILE
        kernel(ker,bufR);
        // (GEN) the right n x n submatrix of kerbas is identity
        // --> retrieve the left part
        for (long i = 0; i < n; ++i)
            for (long j = 0; j < n; ++j)
                K0[i][j] = ker[i][j];
#ifdef MBASIS_GEN_PROFILE
        t_ker += GetWallTime()-tt;
        tt = GetWallTime();
#endif // MBASIS_GEN_PROFILE

        // 2. Update residual 1
        // it is currently   [ [R1_top], [R1_bot] ]
        // it should be [ [R0_top], [K0 * R1_top + R1_bot] ]
        mul(buf, K0, R1_top);
        add(R1_bot, R1_bot, buf);
        R1_top.swap(R0_top); // we do not need the old R1_top anymore, and we won't use R0_top either

#ifdef MBASIS_GEN_PROFILE
        t_res += GetWallTime()-tt;
        tt = GetWallTime();
#endif // MBASIS_GEN_PROFILE

        // 3. compute kernel of residual 1
        // we permute the two blocks top-bottom, to respect the (implicit) shift
        for (long i = 0; i < n; ++i)
            bufR[i].swap(R1_bot[i]);
        for (long i = 0; i < n; ++i)
            bufR[n+i].swap(R1_top[i]);
#ifdef MBASIS_GEN_PROFILE
        t_others += GetWallTime()-tt;
        tt = GetWallTime();
#endif // MBASIS_GEN_PROFILE
        kernel(ker,bufR);
        // (GEN) the right n x n submatrix of kerbas is identity
        // --> retrieve the left part
        for (long i = 0; i < n; ++i)
            for (long j = 0; j < n; ++j)
                K1[i][j] = ker[i][j];
#ifdef MBASIS_GEN_PROFILE
        t_ker += GetWallTime()-tt;
        tt = GetWallTime();
#endif // MBASIS_GEN_PROFILE

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
        mul(P01[0], K1, P11[0]); 
        // add terms of degree 1 ... k-1 of K1 P11
        for (long kk = 1; kk < k; ++kk)
        {
            mul(buf, K1, P11[kk]);
            add(P01[kk], P01[kk], buf);
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

#ifdef MBASIS_GEN_PROFILE
        t_app += GetWallTime()-tt;
        tt = GetWallTime();
#endif // MBASIS_GEN_PROFILE

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
            // Recall that here the vector "P10" stores the coefficients of X P10,
            // and similarly for X P11
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

#ifdef MBASIS_GEN_PROFILE
            t_res += GetWallTime()-tt;
#endif // MBASIS_GEN_PROFILE
        }
    }

    if (odd_order) // order==2d+1: one more iteration
    {
#ifdef MBASIS_GEN_PROFILE
    tt = GetWallTime();
#endif // MBASIS_GEN_PROFILE
        // update the residual
        // R0_top = coeff of degree order-1 of (X^d I + P00) F_top + P01 F_bot,
        // where P00 and P01 have degree d-1
        mul(R0_top, P00[0], F_top[order-1]);
        for (long kk = 1; kk < d; ++kk)
        {
            mul(buf, P00[kk], F_top[order-1-kk]);
            add(R0_top, R0_top, buf);
        }
        add(R0_top, R0_top, F_top[d]);
        for (long kk = 0; kk < d; ++kk)
        {
            mul(buf, P01[kk], F_bot[order-1-kk]);
            add(R0_top, R0_top, buf);
        }

        // R0_bot = coeff of degree order-1 of X P10 F_top + (X^d I + X P11) F_bot,
        // where P10 has degree d-1 and P11 has degree d-2
        // Recall that here the vector "P10" stores the coefficients of X P10,
        // and similarly for X P11
        mul(R0_bot, P10[1], F_top[order-2]);
        for (long kk = 2; kk < d+1; ++kk)
        {
            mul(buf, P10[kk], F_top[order-1-kk]);
            add(R0_bot, R0_bot, buf);
        }
        for (long kk = 1; kk < d; ++kk)
        {
            mul(buf, P11[kk], F_bot[order-1-kk]);
            add(R0_bot, R0_bot, buf);
        }
        add(R0_bot, R0_bot, F_bot[d]);
#ifdef MBASIS_GEN_PROFILE
        t_res += GetWallTime()-tt;
        tt = GetWallTime();
#endif // MBASIS_GEN_PROFILE
        // compute left kernel of residual
        for (long i = 0; i < n; ++i)
            for (long j = 0; j < n; ++j)
                bufR[i][j] = R0_top[i][j];
        for (long i = 0; i < n; ++i)
            for (long j = 0; j < n; ++j)
                bufR[n+i][j] = R0_bot[i][j];
#ifdef MBASIS_GEN_PROFILE
        t_others += GetWallTime()-tt;
        tt = GetWallTime();
#endif // MBASIS_GEN_PROFILE
        kernel(ker,bufR);
        // (GEN) the right n x n submatrix of kerbas is identity
        // --> retrieve the left part
        for (long i = 0; i < n; ++i)
            for (long j = 0; j < n; ++j)
                K0[i][j] = ker[i][j];
#ifdef MBASIS_GEN_PROFILE
        t_ker += GetWallTime()-tt;
        tt = GetWallTime();
#endif // MBASIS_GEN_PROFILE

        // update approximant basis, by  [ [XI, 0], [K0, I] ]
        // Recall: currently, appbas has the form
        // [ [X^d I + P00,  P01], [X P10, X^d I + X P11]]
        // where P00, P01, P10 have degree d-1 and P11 has degree d-2
        // (negative d-1 or d-2 means zero matrix)
        // --> update it as
        // [ [X^{d+1} I + X P00,  X P01],
        //   [K0 P00 + X^d K0 + X P10, K0 P01 + X^d I + X P11] ]

        // bottom left
        // add constant term of K0 P00
        mul(P10[0], K0, P00[0]); 
        // add terms of degree 1 ... d-1 of K0 P00
        for (long kk = 1; kk < d; ++kk)
        {
            mul(buf, K0, P00[kk]);
            add(P10[kk], P10[kk], buf);
        }
        add(P10[d], P10[d], K0); // add K0 to coeff of degree d

        // bottom right
        // add constant term of K0 P01
        mul(P11[0], K0, P01[0]); 
        // add terms of degree 1 ... d-1 of K0 P01
        for (long kk = 1; kk < d; ++kk)
        {
            mul(buf, K0, P01[kk]);
            add(P11[kk], P11[kk], buf);
        }

        // top left: P00 <- X P00  (P00 has degree d-1)
        for (long kk = d-1; kk >=0; --kk)
            P00[kk+1].swap(P00[kk]);
        // since the X^d I was actually not stored in P00[d], that matrix was
        // not initialized and therefore after this loop P00[0] is non-initialized

        // top right: P01 <- X P01  (P01 has degree d-1)
        for (long kk = d-1; kk >=0; --kk)
            P01[kk+1].swap(P01[kk]);
        // note: this loop does put a non-initialized matrix in P01[0]
#ifdef MBASIS_GEN_PROFILE
        t_app += GetWallTime()-tt;
        tt = GetWallTime();
#endif // MBASIS_GEN_PROFILE
        // Now appbas has the form [ [X^{d+1} I + X P00,  X P01], [P10, X^d I + P11] ]
        // where P00, P01, P11 have degree d-1 and P10 has degree d
        // (the identity matrices were not stored, and the constant zero matrices were)
        // --> convert to polynomial matrix format
        // first, without the identity matrix
        appbas.SetDims(2*n,2*n);
        for (long i = 0; i < n; ++i)
        {
            for (long j = 0; j < n; ++j)
            {
                appbas[i][j].SetLength(d+1);
                clear(appbas[i][j][0]);
                for (long k = 1; k < d+1; ++k)
                    appbas[i][j][k] = P00[k][i][j];
            }
            for (long j = n; j < 2*n; ++j)
            {
                long jmn = j-n;
                appbas[i][j].SetLength(d+1);
                clear(appbas[i][j][0]);
                for (long k = 1; k < d+1; ++k)
                    appbas[i][j][k] = P01[k][i][jmn];
            }
        }
        for (long i = n; i < 2*n; ++i)
        {
            long imn = i-n;
            for (long j = 0; j < n; ++j)
            {
                appbas[i][j].SetLength(d+1);
                for (long k = 0; k < d+1; ++k)
                    appbas[i][j][k] = P10[k][imn][j];
            }
            for (long j = n; j < 2*n; ++j)
            {
                long jmn = j-n;
                appbas[i][j].SetLength(d);
                for (long k = 0; k < d; ++k)
                    appbas[i][j][k] = P11[k][imn][jmn];
            }
        }

        // add X^{d+1} I to leading principal submatrix
        for (long i = 0; i < n; ++i)
            SetCoeff(appbas[i][i], d+1);
        // add X^d I to trailing principal submatrix
        for (long i = n; i < 2*n; ++i)
            SetCoeff(appbas[i][i], d);
    }
    else // order==2d: iterations finished
    {
#ifdef MBASIS_GEN_PROFILE
        tt = GetWallTime();
#endif // MBASIS_GEN_PROFILE
        // convert to polynomial matrix format
        // appbas = [ [P00,  P01], [X P10, X P11]]
        // where P00, P01, P10 have degree d-1 and P11 has degree d-2
        appbas.SetDims(2*n,2*n);
        for (long i = 0; i < n; ++i)
        {
            for (long j = 0; j < n; ++j)
            {
                appbas[i][j].SetLength(d);
                for (long k = 0; k < d; ++k)
                    appbas[i][j][k] = P00[k][i][j];
            }
            for (long j = n; j < 2*n; ++j)
            {
                long jmn = j-n;
                appbas[i][j].SetLength(d);
                for (long k = 0; k < d; ++k)
                    appbas[i][j][k] = P01[k][i][jmn];
            }
        }
        for (long i = n; i < 2*n; ++i)
        {
            long imn = i-n;
            for (long j = 0; j < n; ++j)
            {
                appbas[i][j].SetLength(d+1);
                clear(appbas[i][j][0]);
                for (long k = 1; k < d+1; ++k)
                    appbas[i][j][k] = P10[k][imn][j];
            }
            for (long j = n; j < 2*n; ++j)
            {
                long jmn = j-n;
                appbas[i][j].SetLength(d);
                clear(appbas[i][j][0]);
                for (long k = 1; k < d; ++k)
                    appbas[i][j][k] = P11[k][imn][jmn];
            }
        }

        // add X^d I
        for (long i = 0; i < 2*n; ++i)
            SetCoeff(appbas[i][i], d);
    }

#ifdef MBASIS_GEN_PROFILE
    t_others += GetWallTime()-tt;
    double t_total = t_res + t_app + t_ker + t_others;
    std::cout << "~~mbasisgen_rescomp~~\t (residuals,basis,kernel,others): \t ";
    std::cout << t_res/t_total << "," << t_app/t_total << "," <<
    t_ker/t_total << "," << t_others/t_total << std::endl;
#endif // MBASIS_GEN_PROFILE
}

/*------------------------------------------------------------*/
/* Iterative, resupdate version                               */
/*------------------------------------------------------------*/
// Requirement: m = 2*n ; pmat generic ; order >= 2
// if order = 2d, then 
//      appbas = [ [X^d I + P00,  P01], [X P10, X^d I + X P11]]
// where P00, P01, P10 have degree d-1 and P11 has degree d-2
// if order = 2d+1, then
//      appbas = [ [X^{d+1} I + X P00,  X P01], [P10, X^d I + P11] ]
// where P00, P01, P11 have degree d-1 and P10 has degree d
// TODO this shift-middle-product is not very satisfactory, should be
// handled by middle-product itself
// TODO better handle pmat of degree << order ?
void mbasis_generic_2n_n_resupdate(
                                   Mat<zz_pX> & appbas,
                                   const Mat<zz_pX> & pmat,
                                   const long order
                                  )
{
#ifdef MBASIS_GEN_PROFILE
    double t_ker=0.0, t_res=0.0, t_app=0.0, t_others=0.0;
    double tt = GetWallTime();
#endif // MBASIS_GEN_PROFILE

    long n = pmat.NumCols();
    long d = order/2;
    bool odd_order = (order%2);

    // residual = coefficient matrices of input polynomial matrix
    // pmat = [[R_top], [R_bot]]
    Vec<Mat<zz_p>> R_top, R_bot;
    conv_top_bot(R_top, R_bot, pmat, order);

    // coefficient matrices of output approximant basis
    Vec<Mat<zz_p>> P00, P01, P10, P11;
    if (odd_order)
    {
        // appbas = [ [X^{d+1} I + X P00,  X P01], [P10, X^d I + P11] ]
        // where P00, P01, P11 have degree d-1 and P10 has degree d
        P00.SetLength(d+1); // we still store the degree 0 coeff (which will eventually be zero)
        P01.SetLength(d+1); // we still store the degree 0 coeff (which will eventually be zero)
        P10.SetLength(d+1);
        P11.SetLength(d);
    }
    else
    {
        // appbas = [ [X^d I + P00,  P01], [X P10, X^d I + X P11]]
        // where P00, P01, P10 have degree d-1 and P11 has degree d-2
        P00.SetLength(d);
        P01.SetLength(d);
        P10.SetLength(d+1); // we still store the degree 0 coeff (which will eventually be zero)
        P11.SetLength(d); // we still store the degree 0 coeff (which will eventually be zero)
    }

    // To store the kernel of the residuals 0 and 1, and their lefthand square
    // submatrices
    // *GEN* --> dimension of kernel is n x m, its n x n righthand submatrix is
    // the identity
    Mat<zz_p> ker, K0, K1;

    // buffer, to store products and sums
    Mat<zz_p> bufR, buf;

#ifdef MBASIS_GEN_PROFILE
    t_others += GetWallTime()-tt;
#endif // MBASIS_GEN_PROFILE

    // Iteration k==0:
    // --> compute approximant basis at order 2
    // --> update residuals

#ifdef MBASIS_GEN_PROFILE
    tt = GetWallTime();
#endif // MBASIS_GEN_PROFILE
    // 1. compute left kernel of coeff 0 of pmat
    bufR.SetDims(2*n, n);
    for (long i = 0; i < n; ++i)
        for (long j = 0; j < n; ++j)
            bufR[i][j] = R_top[0][i][j];
    for (long i = 0; i < n; ++i)
        for (long j = 0; j < n; ++j)
            bufR[n+i][j] = R_bot[0][i][j];
#ifdef MBASIS_GEN_PROFILE
    t_others += GetWallTime()-tt;
    tt = GetWallTime();
#endif // MBASIS_GEN_PROFILE
    kernel(ker,bufR);
    // (GEN) the right n x n submatrix of kerbas is identity
    // --> retrieve the left part
    K0.SetDims(n, n);
    for (long i = 0; i < n; ++i)
        for (long j = 0; j < n; ++j)
            K0[i][j] = ker[i][j];
#ifdef MBASIS_GEN_PROFILE
    t_ker += GetWallTime()-tt;
    tt = GetWallTime();
#endif // MBASIS_GEN_PROFILE

    // 2. and 3. Second kernel, of coeff of degree 1 of appbas * pmat, that
    // is, [ [R_top[0]], [K0 * R_top[1] + R_bot[1]] ]
    // we permute the two blocks top-bottom, to respect the (implicit) shift
    mul(buf, K0, R_top[1]);
    add(buf, buf, R_bot[1]);
#ifdef MBASIS_GEN_PROFILE
    t_res += GetWallTime()-tt;
    tt = GetWallTime();
#endif // MBASIS_GEN_PROFILE
    for (long i = 0; i < n; ++i) // bufR bottom <- F_top[0]
        bufR[n+i].swap(bufR[i]);
    for (long i = 0; i < n; ++i) // bufR top <- K0*F_top[1] + F_bot[1]
        bufR[i].swap(buf[i]);
#ifdef MBASIS_GEN_PROFILE
    t_others += GetWallTime()-tt;
    tt = GetWallTime();
#endif // MBASIS_GEN_PROFILE
    kernel(ker,bufR);
    // (GEN) the right n x n submatrix of kerbas is identity
    // --> retrieve the left part
    K1.SetDims(n, n);
    for (long i = 0; i < n; ++i)
        for (long j = 0; j < n; ++j)
            K1[i][j] = ker[i][j];
#ifdef MBASIS_GEN_PROFILE
    t_ker += GetWallTime()-tt;
    tt = GetWallTime();
#endif // MBASIS_GEN_PROFILE

    // 4. update approximant basis

    // 4.1 update by first computed basis [ [XI, 0], [K0, I] ]
    // 4.2 update by second computed basis [ [I, K1], [0, X I] ]
    // appbas is currently the identity matrix
    // --> update it as
    // [ [X I + K1 K0,  K1],  [X K0, X I] ]
    //   (remember the identity terms are not stored)
    mul(P00[0], K1, K0); 
    P01[0] = K1;
    P10[1] = K0;

    // --> now appbas is a 0-ordered weak Popov approximant basis of pmat
    // at order 2, of the form
    // [ [X I + P00,  P01], [X P10, X^{k+1} I + X P11]]
    // where P00, P01, P10 have degree 0 and P11 == 0

#ifdef MBASIS_GEN_PROFILE
    t_app += GetWallTime()-tt;
    tt = GetWallTime();
#endif // MBASIS_GEN_PROFILE

    // 5. if not finished (d>1), update the residuals (resupdate variant)
    // that is, left-multiply by [ [X I + K1 K0,  K1],  [X K0, X I] ]
    // --> to minimize the number of muls, do it in two steps:
    //       1/ left-mul by [ [X I, 0], [K0, I] ]
    //       2/ left-mul by [ [I, K1], [0, X I] ]
    if (d>1)
    {
        // [note] eventually, we don't need terms of order <= 2 of the residual
        // since they are known to be zero
        // Step 1/
        // [note] ==> here, no need for updating:
        //              R_top[0], R_top[1], R_bot[0]
        // R_bot += K0 R_top
        for (long kk = 1; kk < order; ++kk)
        {
            mul(buf, K0, R_top[kk]);
            add(R_bot[kk], R_bot[kk], buf);
        }
        // R_top *= X
        for (long kk = order-1; kk >= 2; --kk)
            R_top[kk].swap(R_top[kk-1]);
        // Step 2/
        // [note] ==> here, no need for updating:
        //              R_top[0], R_top[1], R_bot[0], R_bot[1]
        // R_top += K1 * R_bot
        for (long kk = 2; kk < order; ++kk)
        {
            mul(buf, K1, R_bot[kk]);
            add(R_top[kk], R_top[kk], buf);
        }
        // R_bot *= X
        for (long kk = order-1; kk >= 2; --kk)
            R_bot[kk].swap(R_bot[kk-1]);
#ifdef MBASIS_GEN_PROFILE
        t_res += GetWallTime()-tt;
#endif // MBASIS_GEN_PROFILE
    }

    // Iterations k==1...d-1
    for (long k=1; k<d; ++k)
    {
        // *GEN* --> currently, the computed approximant basis has the form
        // [ [X^k I + P00,  P01], [X P10, X^k I + X P11]]
        // where P00, P01, P10 have degree k-1 and P11 has degree k-2
        // It is a 0-ordered weak Popov approximant basis for pmat at order 2*k
        // (For k==0: the last four matrices are in fact zero, and appbas = I)

#ifdef MBASIS_GEN_PROFILE
        tt = GetWallTime();
#endif // MBASIS_GEN_PROFILE
        // 1. compute left kernel of residual
        for (long i = 0; i < n; ++i)
            for (long j = 0; j < n; ++j)
                bufR[i][j] = R_top[2*k][i][j];
        for (long i = 0; i < n; ++i)
            for (long j = 0; j < n; ++j)
                bufR[n+i][j] = R_bot[2*k][i][j];
#ifdef MBASIS_GEN_PROFILE
        t_others += GetWallTime()-tt;
        tt = GetWallTime();
#endif // MBASIS_GEN_PROFILE
        kernel(ker,bufR);
        // (GEN) the right n x n submatrix of kerbas is identity
        // --> retrieve the left part
        for (long i = 0; i < n; ++i)
            for (long j = 0; j < n; ++j)
                K0[i][j] = ker[i][j];
#ifdef MBASIS_GEN_PROFILE
        t_ker += GetWallTime()-tt;
        tt = GetWallTime();
#endif // MBASIS_GEN_PROFILE

        // 2. Second kernel, of coeff of degree 2*k+1 of appbas * pmat, that is,
        // [ [R_top[2*k]], [K0 * R_top[2*k+1] + R_bot[2*k+1]] ]
        mul(buf, K0, R_top[2*k+1]);
        add(buf, buf, R_bot[2*k+1]);

#ifdef MBASIS_GEN_PROFILE
        t_res += GetWallTime()-tt;
        tt = GetWallTime();
#endif // MBASIS_GEN_PROFILE

        // 3. compute kernel of residual 1
        // we permute the two blocks top-bottom, to respect the (implicit) shift
        for (long i = 0; i < n; ++i) // bufR bottom <- R_top[2*k]
            bufR[n+i].swap(bufR[i]);
        for (long i = 0; i < n; ++i) // bufR top <- K0 R_top[2*k+1] + R_bot[2*k+1]
            bufR[i].swap(buf[i]);
#ifdef MBASIS_GEN_PROFILE
        t_others += GetWallTime()-tt;
        tt = GetWallTime();
#endif // MBASIS_GEN_PROFILE
        kernel(ker,bufR);
        // (GEN) the right n x n submatrix of kerbas is identity
        // --> retrieve the left part
        for (long i = 0; i < n; ++i)
            for (long j = 0; j < n; ++j)
                K1[i][j] = ker[i][j];
#ifdef MBASIS_GEN_PROFILE
        t_ker += GetWallTime()-tt;
        tt = GetWallTime();
#endif // MBASIS_GEN_PROFILE

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

        // bottom right
        // add constant term of K0 P01
        mul(P11[0], K0, P01[0]); 
        // add terms of degree 1 ... k-1 of K0 P01
        for (long kk = 1; kk < k; ++kk)
        {
            mul(buf, K0, P01[kk]);
            add(P11[kk], P11[kk], buf);
        }

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
        mul(P01[0], K1, P11[0]); 
        // add terms of degree 1 ... k-1 of K1 P11
        for (long kk = 1; kk < k; ++kk)
        {
            mul(buf, K1, P11[kk]);
            add(P01[kk], P01[kk], buf);
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

#ifdef MBASIS_GEN_PROFILE
        t_app += GetWallTime()-tt;
        tt = GetWallTime();
#endif // MBASIS_GEN_PROFILE

        // 5. if not finished (k<d-1), update the residuals (resupdate variant)
        if (k<d-1)
        {
            // --> left-multiply by [ [X I, 0], [K0, I] ]
            // here, no need for updating R_top[0...2*k+1], R_bot[0...2*k]
            // R_bot += K0 R_top
            for (long kk = 2*k+1; kk < order; ++kk)
            {
                mul(buf, K0, R_top[kk]);
                add(R_bot[kk], R_bot[kk], buf);
            }
            // R_top *= X
            for (long kk = order-1; kk >= 2*k+2; --kk)
                R_top[kk].swap(R_top[kk-1]);

            // --> left-multiply by [[I, K1], [0, X I]]
            // here, no need for updating R_top[0...2*k+1], R_bot[0...2*k+1]
            // R_top += K1 * R_bot
            for (long kk = 2*k+2; kk < order; ++kk)
            {
                mul(buf, K1, R_bot[kk]);
                add(R_top[kk], R_top[kk], buf);
            }
            // R_bot *= X
            for (long kk = order-1; kk >= 2*k+2; --kk)
                R_bot[kk].swap(R_bot[kk-1]);
#ifdef MBASIS_GEN_PROFILE
            t_res += GetWallTime()-tt;
#endif // MBASIS_GEN_PROFILE
        }
    }

    if (odd_order) // order==2d+1: one more iteration
    {
#ifdef MBASIS_GEN_PROFILE
    tt = GetWallTime();
#endif // MBASIS_GEN_PROFILE
        // update the residual
        // --> left-multiply by [ [X I, 0], [K0, I] ]
        // here, no need for updating R_top[0...order-2], R_bot[0...order-3]
        // R_bot += K0 R_top
        mul(buf, K0, R_top[order-2]);
        add(R_bot[order-2], R_bot[order-2], buf);
        mul(buf, K0, R_top[order-1]);
        add(R_bot[order-1], R_bot[order-1], buf);
        // R_top *= X
        R_top[order-1].swap(R_top[order-2]);

        // --> left-multiply by [[I, K1], [0, X I]]
        // here, no need for updating R_top[0...order-2], R_bot[0...order-2]
        // R_top += K1 * R_bot
        mul(buf, K1, R_bot[order-1]);
        add(R_top[order-1], R_top[order-1], buf);
        // R_bot *= X
        R_bot[order-1].swap(R_bot[order-2]);
#ifdef MBASIS_GEN_PROFILE
        t_res += GetWallTime()-tt;
        tt = GetWallTime();
#endif // MBASIS_GEN_PROFILE
        // compute left kernel of residual
        for (long i = 0; i < n; ++i)
            for (long j = 0; j < n; ++j)
                bufR[i][j] = R_top[order-1][i][j];
        for (long i = 0; i < n; ++i)
            for (long j = 0; j < n; ++j)
                bufR[n+i][j] = R_bot[order-1][i][j];
#ifdef MBASIS_GEN_PROFILE
        t_others += GetWallTime()-tt;
        tt = GetWallTime();
#endif // MBASIS_GEN_PROFILE
        kernel(ker,bufR);
        // (GEN) the right n x n submatrix of kerbas is identity
        // --> retrieve the left part
        for (long i = 0; i < n; ++i)
            for (long j = 0; j < n; ++j)
                K0[i][j] = ker[i][j];
#ifdef MBASIS_GEN_PROFILE
        t_ker += GetWallTime()-tt;
        tt = GetWallTime();
#endif // MBASIS_GEN_PROFILE

        // update approximant basis, by  [ [XI, 0], [K0, I] ]
        // Recall: currently, appbas has the form
        // [ [X^d I + P00,  P01], [X P10, X^d I + X P11]]
        // where P00, P01, P10 have degree d-1 and P11 has degree d-2
        // (negative d-1 or d-2 means zero matrix)
        // --> update it as
        // [ [X^{d+1} I + X P00,  X P01],
        //   [K0 P00 + X^d K0 + X P10, K0 P01 + X^d I + X P11] ]

        // bottom left
        // add constant term of K0 P00
        mul(P10[0], K0, P00[0]); 
        // add terms of degree 1 ... d-1 of K0 P00
        for (long kk = 1; kk < d; ++kk)
        {
            mul(buf, K0, P00[kk]);
            add(P10[kk], P10[kk], buf);
        }
        add(P10[d], P10[d], K0); // add K0 to coeff of degree d

        // bottom right
        // add constant term of K0 P01
        mul(P11[0], K0, P01[0]); 
        // add terms of degree 1 ... d-1 of K0 P01
        for (long kk = 1; kk < d; ++kk)
        {
            mul(buf, K0, P01[kk]);
            add(P11[kk], P11[kk], buf);
        }

        // top left: P00 <- X P00  (P00 has degree d-1)
        for (long kk = d-1; kk >=0; --kk)
            P00[kk+1].swap(P00[kk]);
        // since the X^d I was actually not stored in P00[d], that matrix was
        // not initialized and therefore after this loop P00[0] is non-initialized

        // top right: P01 <- X P01  (P01 has degree d-1)
        for (long kk = d-1; kk >=0; --kk)
            P01[kk+1].swap(P01[kk]);
        // note: this loop does put a non-initialized matrix in P01[0]
#ifdef MBASIS_GEN_PROFILE
        t_app += GetWallTime()-tt;
        tt = GetWallTime();
#endif // MBASIS_GEN_PROFILE
        // Now appbas has the form [ [X^{d+1} I + X P00,  X P01], [P10, X^d I + P11] ]
        // where P00, P01, P11 have degree d-1 and P10 has degree d
        // (the identity matrices were not stored, and the constant zero matrices were)
        // --> convert to polynomial matrix format
        // first, without the identity matrix
        appbas.SetDims(2*n,2*n);
        for (long i = 0; i < n; ++i)
        {
            for (long j = 0; j < n; ++j)
            {
                appbas[i][j].SetLength(d+1);
                clear(appbas[i][j][0]);
                for (long k = 1; k < d+1; ++k)
                    appbas[i][j][k] = P00[k][i][j];
            }
            for (long j = n; j < 2*n; ++j)
            {
                long jmn = j-n;
                appbas[i][j].SetLength(d+1);
                clear(appbas[i][j][0]);
                for (long k = 1; k < d+1; ++k)
                    appbas[i][j][k] = P01[k][i][jmn];
            }
        }
        for (long i = n; i < 2*n; ++i)
        {
            long imn = i-n;
            for (long j = 0; j < n; ++j)
            {
                appbas[i][j].SetLength(d+1);
                for (long k = 0; k < d+1; ++k)
                    appbas[i][j][k] = P10[k][imn][j];
            }
            for (long j = n; j < 2*n; ++j)
            {
                long jmn = j-n;
                appbas[i][j].SetLength(d);
                for (long k = 0; k < d; ++k)
                    appbas[i][j][k] = P11[k][imn][jmn];
            }
        }

        // add X^{d+1} I to leading principal submatrix
        for (long i = 0; i < n; ++i)
            SetCoeff(appbas[i][i], d+1);
        // add X^d I to trailing principal submatrix
        for (long i = n; i < 2*n; ++i)
            SetCoeff(appbas[i][i], d);
    }
    else // order==2d: iterations finished
    {
#ifdef MBASIS_GEN_PROFILE
        tt = GetWallTime();
#endif // MBASIS_GEN_PROFILE
        // now appbas has the form [ [X^d I + P00,  P01], [X P10, X^d I + X P11]]
        // where P00, P01, P10 have degree d-1 and P11 has degree d-2
        // convert to polynomial matrix format
        // first, without the identity matrix
        appbas.SetDims(2*n,2*n);
        for (long i = 0; i < n; ++i)
        {
            for (long j = 0; j < n; ++j)
            {
                appbas[i][j].SetLength(d);
                for (long k = 0; k < d; ++k)
                    appbas[i][j][k] = P00[k][i][j];
            }
            for (long j = n; j < 2*n; ++j)
            {
                long jmn = j-n;
                appbas[i][j].SetLength(d);
                for (long k = 0; k < d; ++k)
                    appbas[i][j][k] = P01[k][i][jmn];
            }
        }
        for (long i = n; i < 2*n; ++i)
        {
            long imn = i-n;
            for (long j = 0; j < n; ++j)
            {
                appbas[i][j].SetLength(d+1);
                clear(appbas[i][j][0]);
                for (long k = 1; k < d+1; ++k)
                    appbas[i][j][k] = P10[k][imn][j];
            }
            for (long j = n; j < 2*n; ++j)
            {
                long jmn = j-n;
                appbas[i][j].SetLength(d);
                clear(appbas[i][j][0]);
                for (long k = 1; k < d; ++k)
                    appbas[i][j][k] = P11[k][imn][jmn];
            }
        }

        // add X^d I
        for (long i = 0; i < 2*n; ++i)
            SetCoeff(appbas[i][i], d);
    }
#ifdef MBASIS_GEN_PROFILE
    t_others += GetWallTime()-tt;
    double t_total = t_res + t_app + t_ker + t_others;
    std::cout << "~~mbasisgen_rescomp~~\t (residuals,basis,kernel,others): \t ";
    std::cout << t_res/t_total << "," << t_app/t_total << "," <<
    t_ker/t_total << "," << t_others/t_total << std::endl;
#endif // MBASIS_GEN_PROFILE
}



/*------------------------------------------------------------*/
/* Divide and Conquer: PMBasis                                */
/* (Base case: mbasis-resupdate)                              */
/*------------------------------------------------------------*/
// Requirement: m = 2*n ; pmat generic ; order >= 2
void pmbasis_generic_2n_n(
                          Mat<zz_pX> & appbas,
                          const Mat<zz_pX> & pmat,
                          const long order
                         )
{
    if (order <= 32) // TODO thresholds to be determined
    {
        mbasis_generic_2n_n_resupdate(appbas,pmat,order);
        return;
    }

    // to avoid having to deal with shifts, we use the following orders
    // order1+order2 = order for the recursive calls (see remarks in the
    // documentation):
    // if order is odd: order1 is the one of floor(order/2) and ceil(order/2)
    // which is even
    // if order is even: choose both order1 and order2 even and
    // equal to order/2 +- 1
    const long order1 = (order%4==0 || order%4==1) ? (order/2) : (order/2+1);
    const long order2 = order-order1;

    // first recursive call, with 'pmat'
    Mat<zz_pX> trunc_pmat;
    trunc(trunc_pmat,pmat,order1);
    pmbasis_generic_2n_n(appbas,trunc_pmat,order1);

    // residual = (appbas * pmat * X^-order1) mod X^order2
    // residual = (appbas * RightShift(pmat, order1/2) * X^-ceil(order1/2)) mod X^order2
    // (right shift allowed since deg(appbas) = ceil(order1 / 2) )
    Mat<zz_pX> shift_pmat, residual;
    RightShift(shift_pmat, pmat, order1/2);
    middle_product(residual, appbas, shift_pmat, order1 - order1/2, order2-1);

    // second recursive call, with 'residual'
    Mat<zz_pX> appbas2;
    pmbasis_generic_2n_n(appbas2,residual,order2);

    // final basis = appbas2 * appbas
    multiply(appbas,appbas2,appbas);
}

// TODO doc if this turns out useful
// returns just the n top rows of the appbas (saves a bit on the matrix
// products)
void pmbasis_generic_2n_n_top_rows(
                          Mat<zz_pX> & appbas,
                          const Mat<zz_pX> & pmat,
                          const long order
                         )
{
    if (order <= 32) // TODO thresholds to be determined
    {
        mbasis_generic_2n_n_resupdate(appbas,pmat,order);
        appbas.SetDims(pmat.NumCols(), pmat.NumRows());
        // TODO could do better by modifying mbasis, but this
        // will be called just once, so...
        return;
    }

    // to avoid having to deal with shifts, we use the following orders
    // order1+order2 = order for the recursive calls (see remarks in the
    // documentation):
    // if order is odd: order1 is the one of floor(order/2) and ceil(order/2)
    // which is even
    // if order is even: choose both order1 and order2 even and
    // equal to order/2 +- 1
    const long order1 = (order%4==0 || order%4==1) ? (order/2) : (order/2+1);
    const long order2 = order-order1;

    // first recursive call, with 'pmat'
    Mat<zz_pX> trunc_pmat;
    trunc(trunc_pmat,pmat,order1);
    pmbasis_generic_2n_n(appbas,trunc_pmat,order1);

    // residual = (appbas * pmat * X^-order1) mod X^order2
    // residual = (appbas * RightShift(pmat, order1/2) * X^-ceil(order1/2)) mod X^order2
    // (right shift allowed since deg(appbas) = ceil(order1 / 2) )
    Mat<zz_pX> shift_pmat, residual;
    RightShift(shift_pmat, pmat, order1/2);
    middle_product(residual, appbas, shift_pmat, order1 - order1/2, order2-1);

    // second recursive call, with 'residual'
    Mat<zz_pX> appbas2;
    pmbasis_generic_2n_n_top_rows(appbas2,residual,order2);

    // final basis = appbas2 * appbas
    multiply(appbas,appbas2,appbas);
}


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* MATRIX-PADE APPROXIMATION -- GENERIC INPUT -- no shift     */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* Iterative algorithm, for low approximation order           */
/*------------------------------------------------------------*/
// Note: it is precisely mbasis_generic_2n_n_resupdate above, except
// that only the left-block-column of appbas is computed
// (namely P00, here D1, and P10, here D2)
// TODO try rescomp version --> not urgent: this has roughly no impact on the
// interesting case (large order) since in the divide and conquer this base
// case will almost never be called, most base cases will be regular pm-basis
// base cases
void matrix_pade_generic_iterative(
                                   Mat<zz_pX> & den1,
                                   Mat<zz_pX> & den2,
                                   const Mat<zz_pX> & pmat,
                                   const long order
                                  )
{
#ifdef MATRIX_PADE_GEN_PROFILE
    double t_ker=0.0, t_res=0.0, t_app=0.0, t_others=0.0;
    double tt = GetWallTime();
#endif // MATRIX_PADE_GEN_PROFILE

    long n = pmat.NumCols(); // == pmat.NumRows()
    long d = order/2;
    bool odd_order = (order%2);

    // residual = coefficient matrices of [[pmat], [-Id]]
    // written as [[R_top], [R_bot]]
    Vec<Mat<zz_p>> R_top, R_bot;
    conv(R_top, pmat, order);
    R_bot.SetLength(order);
    // R_bot should currently be identity and will be initialized later, after
    // the first iteration

    // coefficient matrices of output matrices den1 and den2
    // If order is even:
    //   --> den1 will be computed in the form X^d I + D1
    //   --> den2 will be computed in the form X D2
    //       where D1 and D2 have degree d-1
    // If order is odd:
    //   --> den1 will be computed in the form X^{d+1} I + D1
    //   --> den2 will be computed in the form D2
    //       where D1 and D2 have degree d
    // These matrices are the left blocks P00 and P10 in the algorithm
    // mbasis_generic_2n_n_resupdate above on the input [[pmat],[-I]]
    // --> we simply follow this algorithm but not storing the right
    // columns since they are not needed here
    Vec<Mat<zz_p>> D1, D2;
    D1.SetLength(odd_order ? d+1 : d); // we don't store the identity
    D2.SetLength(d+1); // we store the degree 0 coeff (even for even order, when it is zero)

    // To store the kernel of the residuals 0 and 1, and their lefthand square
    // submatrices
    // *GEN* --> dimension of kernel is n x 2n, its n x n righthand submatrix is
    // the identity
    Mat<zz_p> ker, K0, K1;
    K0.SetDims(n,n);

    // buffer, to store products and sums
    Mat<zz_p> bufR, buf;

#ifdef MATRIX_PADE_GEN_PROFILE
    t_others += GetWallTime()-tt;
#endif // MATRIX_PADE_GEN_PROFILE

    // Iteration k==0:
    // --> compute approximant basis at order 2
    // --> update residuals

#ifdef MATRIX_PADE_GEN_PROFILE
    tt = GetWallTime();
#endif // MATRIX_PADE_GEN_PROFILE
    // 1. compute left kernel of coeff 0 of [[pmat], [-Id]]
    // --> this is simply [[inverse(pmat(0)), Id]]
    // K0 = inverse(pmat(0)) is stored directly into D2[1]  (see below the appbas update)
    inv(D2[1],R_top[0]);
#ifdef MATRIX_PADE_GEN_PROFILE
    t_ker += GetWallTime()-tt;
    tt = GetWallTime();
#endif // MATRIX_PADE_GEN_PROFILE

    // 2. and 3. Second kernel, of coeff of degree 1 of appbas * pmat, that
    // is, [ [R_top[0]], [K0 * R_top[1] + R_bot[1]] ], with R_bot[1] == 0 here
    // we permute the two blocks top-bottom, to respect the (implicit) shift
    mul(buf, D2[1], R_top[1]);
#ifdef MATRIX_PADE_GEN_PROFILE
    t_res += GetWallTime()-tt;
    tt = GetWallTime();
#endif // MATRIX_PADE_GEN_PROFILE
    bufR.SetDims(2*n, n);
    for (long i = 0; i < n; ++i) // bufR top <- K0*R_top[1]
        bufR[i].swap(buf[i]);
    for (long i = 0; i < n; ++i) // bufR bottom <- R_top[0]
        for (long j = 0; j < n; ++j)
            bufR[n+i][j] = R_top[0][i][j];
#ifdef MATRIX_PADE_GEN_PROFILE
    t_others += GetWallTime()-tt;
    tt = GetWallTime();
#endif // MATRIX_PADE_GEN_PROFILE
    kernel(ker,bufR);
    // (GEN) the right n x n submatrix of kerbas is identity
    // --> retrieve the left part
    K1.SetDims(n, n);
    for (long i = 0; i < n; ++i)
        for (long j = 0; j < n; ++j)
            K1[i][j] = ker[i][j];
#ifdef MATRIX_PADE_GEN_PROFILE
    t_ker += GetWallTime()-tt;
    tt = GetWallTime();
#endif // MATRIX_PADE_GEN_PROFILE

    // 4. update approximant basis

    // 4.1 update by first computed basis [ [XI, 0], [K0, I] ]
    // 4.2 update by second computed basis [ [I, K1], [0, X I] ]
    // appbas is currently the identity matrix
    // --> update it as
    // [ [X I + K1 K0,  K1],  [X K0, X I] ]
    //   (remember the identity terms are not stored)
    mul(D1[0], K1, D2[1]); 
    // D2[1] was already used to store K0

    // --> now appbas is a 0-ordered weak Popov approximant basis of pmat
    // at order 2, of the form
    // [ [X I + P00,  P01], [X P10, X^{k+1} I + X P11]]
    // where P00, P01, P10 have degree 0 and P11 == 0

#ifdef MATRIX_PADE_GEN_PROFILE
    t_app += GetWallTime()-tt;
    tt = GetWallTime();
#endif // MATRIX_PADE_GEN_PROFILE

    // 5. if not finished (d>1), update the residuals
    // that is, left-multiply by [ [X I + K1 K0,  K1],  [X K0, X I] ]
    // --> to minimize the number of muls, do it in two steps:
    //       1/ left-mul by [ [X I, 0], [K0, I] ]
    //       2/ left-mul by [ [I, K1], [0, X I] ]
    if (d>1)
    {
        // [note] eventually, we don't need terms of order <= 2 of the residual
        // since they are known to be zero
        // Step 1/
        // [note] ==> here, no need for updating:
        //              R_top[0], R_top[1], R_bot[0]
        // R_bot = K0 R_top
        for (long kk = 1; kk < order; ++kk)
            mul(R_bot[kk], D2[1], R_top[kk]);
        // R_top *= X
        for (long kk = order-1; kk >= 2; --kk)
            R_top[kk].swap(R_top[kk-1]);
        // Step 2/
        // [note] ==> here, no need for updating:
        //              R_top[0], R_top[1], R_bot[0], R_bot[1]
        // R_top += K1 * R_bot
        for (long kk = 2; kk < order; ++kk)
        {
            mul(buf, K1, R_bot[kk]);
            add(R_top[kk], R_top[kk], buf);
        }
        // R_bot *= X
        for (long kk = order-1; kk >= 2; --kk)
            R_bot[kk].swap(R_bot[kk-1]);
#ifdef MATRIX_PADE_GEN_PROFILE
        t_res += GetWallTime()-tt;
#endif // MATRIX_PADE_GEN_PROFILE
    }

    // Iterations k==1...d-1
    for (long k=1; k<d; ++k)
    {
        // *GEN* --> currently, the computed approximant basis has the form
        // [ [X^k I + P00,  P01], [X P10, X^k I + X P11]]
        // where P00, P01, P10 have degree k-1 and P11 has degree k-2
        // It is a 0-ordered weak Popov approximant basis for pmat at order 2*k
        // (For k==0: the last four matrices are in fact zero, and appbas = I)
        // Recall that D1 corresponds to X^k I + P00 (Identity not stored),
        // and D2 corresponds to X P10  (constant zero term stored)

#ifdef MATRIX_PADE_GEN_PROFILE
        tt = GetWallTime();
#endif // MATRIX_PADE_GEN_PROFILE
        // 1. compute left kernel of residual
        for (long i = 0; i < n; ++i)
            for (long j = 0; j < n; ++j)
                bufR[i][j] = R_top[2*k][i][j];
        for (long i = 0; i < n; ++i)
            for (long j = 0; j < n; ++j)
                bufR[n+i][j] = R_bot[2*k][i][j];
#ifdef MATRIX_PADE_GEN_PROFILE
        t_others += GetWallTime()-tt;
        tt = GetWallTime();
#endif // MATRIX_PADE_GEN_PROFILE
        kernel(ker,bufR);
        // (GEN) the right n x n submatrix of kerbas is identity
        // --> retrieve the left part
        for (long i = 0; i < n; ++i)
            for (long j = 0; j < n; ++j)
                K0[i][j] = ker[i][j];
#ifdef MATRIX_PADE_GEN_PROFILE
        t_ker += GetWallTime()-tt;
        tt = GetWallTime();
#endif // MATRIX_PADE_GEN_PROFILE

        // 2. Second kernel, of coeff of degree 2*k+1 of appbas * pmat, that is,
        // [ [R_top[2*k]], [K0 * R_top[2*k+1] + R_bot[2*k+1]] ]
        mul(buf, K0, R_top[2*k+1]);
        add(buf, buf, R_bot[2*k+1]);

#ifdef MATRIX_PADE_GEN_PROFILE
        t_res += GetWallTime()-tt;
        tt = GetWallTime();
#endif // MATRIX_PADE_GEN_PROFILE

        // 3. compute kernel of residual 1
        // we permute the two blocks top-bottom, to respect the (implicit) shift
        for (long i = 0; i < n; ++i) // bufR bottom <- R_top[2*k]
            bufR[n+i].swap(bufR[i]);
        for (long i = 0; i < n; ++i) // bufR top <- K0 R_top[2*k+1] + R_bot[2*k+1]
            bufR[i].swap(buf[i]);
#ifdef MATRIX_PADE_GEN_PROFILE
        t_others += GetWallTime()-tt;
        tt = GetWallTime();
#endif // MATRIX_PADE_GEN_PROFILE
        kernel(ker,bufR);
        // (GEN) the right n x n submatrix of kerbas is identity
        // --> retrieve the left part
        for (long i = 0; i < n; ++i)
            for (long j = 0; j < n; ++j)
                K1[i][j] = ker[i][j];
#ifdef MATRIX_PADE_GEN_PROFILE
        t_ker += GetWallTime()-tt;
        tt = GetWallTime();
#endif // MATRIX_PADE_GEN_PROFILE

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
        // add constant term of K0 D1
        mul(D2[0], K0, D1[0]); 
        // add terms of degree 1 ... k-1 of K0 D1
        for (long kk = 1; kk < k; ++kk)
        {
            mul(buf, K0, D1[kk]);
            add(D2[kk], D2[kk], buf);
        }
        add(D2[k], D2[k], K0); // add K0 to coeff of degree k

        // top left: D1 <- X D1  (D1 has degree k-1)
        for (long kk = k-1; kk >=0; --kk)
            D1[kk+1].swap(D1[kk]);
        // since the X^k I was actually not stored in D1[k], that matrix was
        // zero and therefore this loop does put the zero matrix in D1[0]

        // 4.2 update by second computed basis [ [I, K1], [0, X I] ]
        // From step 4.1, appbas has the form
        // [ [X^{k+1} I + X P00,  X P01],
        //   [P10, X^k I + P11] ]
        // where P00, P01, P11 have degree k-1 and P10 has degree k
        // --> update it as
        // [ [X^{k+1} I + X P00 + K1 P10,  X P01 + X^k K1 + K1 P11],
        //   [X P10, X^{k+1} I + X P11] ]

        // top left
        // add constant term of K1 D2
        mul(D1[0], K1, D2[0]); 
        // add terms of degree 1 ... k of K1 D2
        for (long kk = 1; kk <= k; ++kk)
        {
            mul(buf, K1, D2[kk]);
            add(D1[kk], D1[kk], buf);
        }

        // bottom left: D2 <- X D2  (D2 has degree k)
        for (long kk = k; kk >=0; --kk)
            D2[kk+1].swap(D2[kk]);
        // note: this loop does put the zero matrix in P10[0]

        // --> now appbas is a 0-ordered weak Popov approximant basis of pmat
        // at order 2*k+2, of the form
        // [ [X^{k+1} I + P00,  P01], [X P10, X^{k+1} I + X P11]]
        // where P00, P01, P10 have degree k and P11 has degree k-1
        // recall the left blocks are stored in D1 and D2, the right blocks are
        // ignored

#ifdef MATRIX_PADE_GEN_PROFILE
        t_app += GetWallTime()-tt;
        tt = GetWallTime();
#endif // MATRIX_PADE_GEN_PROFILE

        // 5. if not finished (k<d-1), update the residuals
        if (k<d-1)
        {
            // --> left-multiply by [ [X I, 0], [K0, I] ]
            // here, no need for updating R_top[0...2*k+1], R_bot[0...2*k]
            // R_bot += K0 R_top
            for (long kk = 2*k+1; kk < order; ++kk)
            {
                mul(buf, K0, R_top[kk]);
                add(R_bot[kk], R_bot[kk], buf);
            }
            // R_top *= X
            for (long kk = order-1; kk >= 2*k+2; --kk)
                R_top[kk].swap(R_top[kk-1]);

            // --> left-multiply by [[I, K1], [0, X I]]
            // here, no need for updating R_top[0...2*k+1], R_bot[0...2*k+1]
            // R_top += K1 * R_bot
            for (long kk = 2*k+2; kk < order; ++kk)
            {
                mul(buf, K1, R_bot[kk]);
                add(R_top[kk], R_top[kk], buf);
            }
            // R_bot *= X
            for (long kk = order-1; kk >= 2*k+2; --kk)
                R_bot[kk].swap(R_bot[kk-1]);
#ifdef MATRIX_PADE_GEN_PROFILE
            t_res += GetWallTime()-tt;
#endif // MATRIX_PADE_GEN_PROFILE
        }
    }

    if (odd_order) // order==2d+1: one more iteration
    {
#ifdef MATRIX_PADE_GEN_PROFILE
    tt = GetWallTime();
#endif // MATRIX_PADE_GEN_PROFILE
        // update the residual
        // --> left-multiply by [ [X I, 0], [K0, I] ]
        // here, no need for updating R_top[0...order-2], R_bot[0...order-3]
        // R_bot += K0 R_top
        mul(buf, K0, R_top[order-2]);
        add(R_bot[order-2], R_bot[order-2], buf);
        mul(buf, K0, R_top[order-1]);
        add(R_bot[order-1], R_bot[order-1], buf);
        // R_top *= X
        R_top[order-1].swap(R_top[order-2]);

        // --> left-multiply by [[I, K1], [0, X I]]
        // here, no need for updating R_top[0...order-2], R_bot[0...order-2]
        // R_top += K1 * R_bot
        mul(buf, K1, R_bot[order-1]);
        add(R_top[order-1], R_top[order-1], buf);
        // R_bot *= X
        R_bot[order-1].swap(R_bot[order-2]);
#ifdef MATRIX_PADE_GEN_PROFILE
        t_res += GetWallTime()-tt;
        tt = GetWallTime();
#endif // MATRIX_PADE_GEN_PROFILE
        // compute left kernel of residual
        for (long i = 0; i < n; ++i)
            for (long j = 0; j < n; ++j)
                bufR[i][j] = R_top[order-1][i][j];
        for (long i = 0; i < n; ++i)
            for (long j = 0; j < n; ++j)
                bufR[n+i][j] = R_bot[order-1][i][j];
#ifdef MATRIX_PADE_GEN_PROFILE
        t_others += GetWallTime()-tt;
        tt = GetWallTime();
#endif // MATRIX_PADE_GEN_PROFILE
        kernel(ker,bufR);
        // (GEN) the right n x n submatrix of kerbas is identity
        // --> retrieve the left part
        for (long i = 0; i < n; ++i)
            for (long j = 0; j < n; ++j)
                K0[i][j] = ker[i][j];
#ifdef MATRIX_PADE_GEN_PROFILE
        t_ker += GetWallTime()-tt;
        tt = GetWallTime();
#endif // MATRIX_PADE_GEN_PROFILE

        // update approximant basis, by  [ [XI, 0], [K0, I] ]
        // Recall: currently, appbas has the form
        // [ [X^d I + D1,  P01], [X D2, X^d I + X P11]]
        // where D1, D2, P10 have degree d-1 and P11 has degree d-2
        // (negative d-1 or d-2 means zero matrix)
        // --> update it as
        // [ [X^{d+1} I + X D1,  X P01],
        //   [K0 D1 + X^d K0 + X D2, K0 P01 + X^d I + X P11] ]

        // bottom left
        // add constant term of K0 P00
        mul(D2[0], K0, D1[0]); 
        // add terms of degree 1 ... d-1 of K0 P00
        for (long kk = 1; kk < d; ++kk)
        {
            mul(buf, K0, D1[kk]);
            add(D2[kk], D2[kk], buf);
        }
        add(D2[d], D2[d], K0); // add K0 to coeff of degree d

        // top left: P00 <- X P00  (P00 has degree d-1)
        for (long kk = d-1; kk >=0; --kk)
            D1[kk+1].swap(D1[kk]);
#ifdef MATRIX_PADE_GEN_PROFILE
        t_app += GetWallTime()-tt;
        tt = GetWallTime();
#endif // MATRIX_PADE_GEN_PROFILE
        // Now appbas has the form [ [X^{d+1} I + X P00,  X P01], [P10, X^d I + P11] ]
        // where P00, P01, P11 have degree d-1 and P10 has degree d
        // (the identity matrices were not stored, and the constant zero matrices were)
        // --> convert to polynomial matrix format
        // first, without the identity matrix
        den1.SetDims(n,n);
        for (long i = 0; i < n; ++i)
        {
            for (long j = 0; j < n; ++j)
            {
                den1[i][j].SetLength(d+1);
                clear(den1[i][j][0]);
                for (long k = 1; k < d+1; ++k)
                    den1[i][j][k] = D1[k][i][j];
            }
        }

        // add X^{d+1} I to leading principal submatrix
        for (long i = 0; i < n; ++i)
            SetCoeff(den1[i][i], d+1);

        den2.SetDims(n,n);
        for (long i = 0; i < n; ++i)
        {
            for (long j = 0; j < n; ++j)
            {
                den2[i][j].SetLength(d+1);
                for (long k = 0; k < d+1; ++k)
                    den2[i][j][k] = D2[k][i][j];
            }
        }
    }
    else // order==2d: iterations finished
    {
#ifdef MATRIX_PADE_GEN_PROFILE
        tt = GetWallTime();
#endif // MATRIX_PADE_GEN_PROFILE
        // convert to polynomial matrix format
        // appbas = [ [P00,  P01], [X P10, X P11]]
        // where P00, P01, P10 have degree d-1 and P11 has degree d-2
        den1.SetDims(n,n);
        for (long i = 0; i < n; ++i)
            for (long j = 0; j < n; ++j)
            {
                den1[i][j].SetLength(d);
                for (long k = 0; k < d; ++k)
                    den1[i][j][k] = D1[k][i][j];
            }
        // add X^d I
        for (long i = 0; i < n; ++i)
            SetCoeff(den1[i][i], d);

        den2.SetDims(n,n);
        for (long i = 0; i < n; ++i)
            for (long j = 0; j < n; ++j)
            {
                den2[i][j].SetLength(d+1);
                clear(den2[i][j][0]);
                for (long k = 1; k < d+1; ++k)
                    den2[i][j][k] = D2[k][i][j];
            }
    }

#ifdef MATRIX_PADE_GEN_PROFILE
    t_others += GetWallTime()-tt;
    double t_total = t_res + t_app + t_ker + t_others;
    std::cout << "~~mbasisgen_rescomp~~\t (residuals,basis,kernel,others): \t ";
    std::cout << t_res/t_total << "," << t_app/t_total << "," <<
    t_ker/t_total << "," << t_others/t_total << std::endl;
#endif // MATRIX_PADE_GEN_PROFILE
}


/*------------------------------------------------------------*/
/* Divide and conquer algorithm                               */
/*------------------------------------------------------------*/

// computes both den1 and den2, such that [[den1], [den2]] is the first
// block-column of a 0-ordered weak Popov approximant basis for [[pmat], [-Id]]
// at order 'order'
// Note: den1 is in Popov form.
// Note: den2 is used for computing residuals.
void matrix_pade_generic(
                         Mat<zz_pX> & den,
                         const Mat<zz_pX> & pmat,
                         const long order
                        )
{
    if (order <= 32) // TODO thresholds to be determined
    {
        matrix_pade_generic_iterative(den, pmat, order);
        return;
    }

    // to avoid having to deal with shifts, we use the following orders
    // order1+order2 = order for the recursive calls (see remarks in the
    // documentation):
    // if order is odd: order1 is the one of floor(order/2) and ceil(order/2)
    // which is even
    // if order is even: choose both order1 and order2 even and
    // equal to order/2 +- 1
    const long order1 = (order%4==0 || order%4==1) ? (order/2) : (order/2+1);
    const long order2 = order-order1;

    // first recursive call, matrix Pade with 'pmat' truncated at order1
    Mat<zz_pX> trunc_pmat;
    trunc(trunc_pmat,pmat,order1);
    matrix_pade_generic_recursion(den, trunc_pmat, order1);
    // den = both left blocks of the approx basis (not just the top left one)

    // residual = (full appbas1 * [[pmat], [-Id]] * X^-order1) mod X^order2
    // hence the residual here actually only depends on the left blocks of appbas1,
    // in other words, on den as computed above
    // (the right blocks are multiplied by identity, and since they have degree
    // order1/2 this gives only coefficients of degree < order1)
    // As in pmbasis_generic above, since den has degree ceil(order1/2), we
    // can actually first right shift pmat by order1/2 to make the middle product
    // more efficient
    // --> residual = (den * RightShift(pmat, order1/2) * X^-ceil(order1/2)) mod X^order2
    Mat<zz_pX> shift_pmat, residual;
    RightShift(shift_pmat, pmat, order1/2);
    middle_product(residual, den, shift_pmat, order1 - order1/2, order2-1);

    // second recursive call, approximant basis with 'residual'
    // just returns the top block of rows
    Mat<zz_pX> appbas2;
    pmbasis_generic_2n_n_top_rows(appbas2, residual, order2);

    // final matrix Pade denominator = appbas2 * den
    // this is a nx2n * 2nxn product
    multiply(den,appbas2,den);
}

void matrix_pade_generic_recursion(
                                   Mat<zz_pX> & den,
                                   const Mat<zz_pX> & pmat,
                                   const long order
                                  )
{
    if (order <= 32) // TODO thresholds to be determined
    {
        long n = pmat.NumCols();
        Mat<zz_pX> den1, den2;
        matrix_pade_generic_iterative(den1, den2, pmat, order);
        den.SetDims(2*n, n);
        for (long i = 0; i < n; ++i)
            den[i].swap(den1[i]);
        for (long i = n; i < 2*n; ++i)
            den[i].swap(den2[i-n]);
        return;
    }

    // to avoid having to deal with shifts, we use the following orders
    // order1+order2 = order for the recursive calls (see remarks in the
    // documentation):
    // if order is odd: order1 is the one of floor(order/2) and ceil(order/2)
    // which is even
    // if order is even: choose both order1 and order2 even and
    // equal to order/2 +- 1
    const long order1 = (order%4==0 || order%4==1) ? (order/2) : (order/2+1);
    const long order2 = order-order1;

    // first recursive call, matrix Pade with 'pmat' truncated at order1
    // we need both left blocks of the approx basis
    Mat<zz_pX> trunc_pmat;
    trunc(trunc_pmat,pmat,order1);
    matrix_pade_generic_recursion(den, trunc_pmat, order1);

    // residual = (full appbas1 * [[pmat], [-Id]] * X^-order1) mod X^order2
    // hence the residual here actually only depends on the left blocks of appbas1,
    // in other words, on den as computed above
    // (the right blocks are multiplied by identity, and since they have degree
    // order1/2 this gives only coefficients of degree < order1)
    // As above, since den has degree ceil(order1/2), we
    // can actually first right shift pmat by order1/2 to make the middle product
    // more efficient
    // --> residual = (den * RightShift(pmat, order1/2) * X^-ceil(order1/2)) mod X^order2
    Mat<zz_pX> shift_pmat, residual;
    RightShift(shift_pmat, pmat, order1/2);
    middle_product(residual, den, shift_pmat, order1 - order1/2, order2-1);

    // second recursive call, approximant basis with 'residual'
    Mat<zz_pX> appbas2;
    pmbasis_generic_2n_n(appbas2, residual, order2);

    // final matrix Pade denominator = appbas2 * appbas1
    multiply(den,appbas2,den);
}






/*------------------------------------------------------------*/
/* Output Popov form                                          */
/*------------------------------------------------------------*/

void popov_pmbasis_generic(Mat<zz_pX> &appbas, const Mat<zz_pX> & pmat, const long order)
{
    if (order < 32)
    {
        popov_mbasis_rescomp_generic(appbas, pmat, order);
        return;
    }
    VecLong shift(pmat.NumRows());
    popov_pmbasis(appbas, pmat, order, shift);
}


// version with residual constant matrix computed at each iteration
void popov_mbasis_rescomp_generic(
                                  Mat<zz_pX> & appbas,
                                  const Mat<zz_pX> & pmat,
                                  const long order
                                 )
{
#ifdef MBASIS_PROFILE
    double t_others=0.0,t_residual=0.0,t_appbas=0.0,t_kernel=0.0,t_now;
    t_now = GetWallTime();
#endif
    // A. General

    // A.1 dimensions of input matrix
    const long m = pmat.NumRows();
    const long n = pmat.NumCols();

    // A.2 store iota since it will be used at each iteration
    VecLong iota(m);
    std::iota(iota.begin(), iota.end(), 0);

    // B. Input representation; initialize output

    // B.1 convert input into vector of constant matrices (its "coefficients")
    const Vec<Mat<zz_p>> coeffs_pmat = conv(pmat,order);

    // B.2 vector of coefficients of output approximant basis
    Vec<Mat<zz_p>> coeffs_appbas;

    // B.3 initially, appbas is the identity matrix
    coeffs_appbas.SetLength(1);
    ident(coeffs_appbas[0], m);
    // degree of approximant basis, initially zero
    long deg_appbas = 0;
    VecLong shift(pmat.NumRows());
    // `shift` will always hold the shifted row degree of appbas
    // (initially, this is the input shift)

    // C. Residual matrix (m x n constant matrix, next coefficient
    // of appbas * pmat which we want to annihilate)

    // C.1 stores the residual, initially coeffs_pmat[0]
    Mat<zz_p> residuals(coeffs_pmat[0]);

    // C.2 temporary matrices used during the computation of residuals
    Mat<zz_p> res_coeff;

    // C.3 permuted residual, used as input to the kernel at the "base case"
    Mat<zz_p> p_residual(INIT_SIZE, m, n);

    // D. Base case (working modulo X, essentially amounts to finding the left
    // kernel of the permuted residual p_residual)

    // D.1 pivot indices in kernel basis (which is in row echelon form)
    // Note: length is probably overestimated (usually kernel has m-n rows),
    // but this avoids reallocating the right length at each iteration
    VecLong pivind(m-1);
    // Vector indicating if a given column index appears in this pivot index
    // i.e. is_pivind[pivind[i]] = true and others are false
    std::vector<bool> is_pivind(m, false);

    // D.2 permutation for the rows of the constant kernel
    VecLong perm_rows_ker;
    // pivot indices of row echelon form before permutation
    VecLong p_pivind(m-1);

    // D.3 permutation which stable-sorts the shift, used at the base case
    VecLong p_shift;

    // D.4 the constant kernel, and its permuted version
    Mat<zz_p> kerbas;
    Mat<zz_p> p_kerbas;

    // E. Updating appbas
    // stores the product "constant-kernel * coeffs_appbas[d]"
    Mat<zz_p> kerapp; 

#ifdef MBASIS_PROFILE
    t_others += GetWallTime()-t_now;
#endif

    for (long ord = 1; ord <= order; ++ord)
    {
#ifdef MBASIS_PROFILE
        t_now = GetWallTime();
#endif
        // compute permutation which realizes stable sort of the row degree `shift`
        // --> we need to permute things, to take into account the "priority"
        // (i.e. "weights") indicated by the shift
        p_shift = iota;
        stable_sort(p_shift.begin(), p_shift.end(),
                    [&](const long& a, const long& b)->bool
                    {
                    return (shift[a] < shift[b]);
                    } );

        // permute rows of the residual accordingly
        for (long i = 0; i < m; ++i)
            p_residual[i].swap(residuals[p_shift[i]]);
        // content of residual has been changed --> let's make it zero
        // (already the case if ord==1, since residual is then the former p_residual which was zero)
        if (ord>1)
            clear(residuals);
#ifdef MBASIS_PROFILE
        t_others += GetWallTime()-t_now;
        t_now = GetWallTime();
#endif
        // find the (permuted) left kernel basis, hopefully in row echelon form;
        kernel(p_kerbas,p_residual);
#ifdef MBASIS_PROFILE
        t_kernel += GetWallTime()-t_now;
#endif
        const long ker_dim = p_kerbas.NumRows();

        if (ker_dim==0)
        {
            // Exceptional case: the residual matrix has empty left kernel
            // --> no need to compute more: the final basis is X^(order-ord+1)*coeffs_appbas
            // TODO improve: simultaneous convert+shift!!
            appbas = conv(coeffs_appbas);
            appbas <<= (order-ord+1);
            for (long i = 0; i < m; ++i)
                shift[i] += order-ord+1;
            return;
        }

        else if (ker_dim==m)
        {
            // Exceptional case: residual coeff was zero, and kernel 'kerbas' is identity
            // --> approximant basis is already correct for this order, no need to
            // change it or to change shift
            // --> we just need to compute the next residual
            // (unless ord == order, in which case the algorithm returns)
            // this "residual" is the coefficient of degree ord in appbas * pmat:
            // Note: at this point, residuals==0
            if (ord<order)
            {
#ifdef MBASIS_PROFILE
                t_now = GetWallTime();
#endif
                for (long d = std::max<long>(0,ord-coeffs_pmat.length()+1); d <= deg_appbas; ++d)
                {
                    mul(res_coeff, coeffs_appbas[d], coeffs_pmat[ord-d]);
                    add(residuals, residuals, res_coeff);
                }
#ifdef MBASIS_PROFILE
                t_residual += GetWallTime()-t_now;
#endif
            }
        }

        else
        {
            // here, we are in the "usual" case, where the left kernel of the
            // residual has no special shape

            // first, we permute everything back to original order

#ifdef MBASIS_PROFILE
            t_now = GetWallTime();
#endif
            // Compute pivots indices (pivot = rightmost nonzero entry)
            // Experiments show that:
            //   * kernel is expected to be of the form [ K | Id ]
            //   * in general it is a column-permutation of such a matrix
            // However note that a column-permutation is not sufficient for our needs
            // Another property: if pivots are in the expected location (diagonal of
            // rightmost square submatrix), then the corresponding column is the identity column.
            bool expected_pivots = true;
            for (long i = 0; i<ker_dim; ++i)
            {
                p_pivind[i] = m-1;
                while (IsZero(p_kerbas[i][p_pivind[i]]))
                    --p_pivind[i];
                if (p_pivind[i] != m-ker_dim+i)
                    expected_pivots = false;
            }

            if (not expected_pivots)
            {
                // find whether p_pivind has pairwise distinct entries
                // (use pivind as temp space)
                pivind = p_pivind;
                std::sort(pivind.begin(), pivind.end());
                // if pairwise distinct, then fine, the basis will not
                // be Popov but will be ordered weak Popov (the goal of
                // expected_pivots above was just to avoid this call to
                // sort in the most usual case)

                if (std::adjacent_find(pivind.begin(),pivind.end()) != pivind.end())
                {
                    // the kernel is not in a shape we can deduce the appbas from (some pivots collide)
                    // --> let's compute its lower triangular row echelon form
                    // (use kerbas as temporary space)
                    kerbas.SetDims(ker_dim,m);
                    for (long i = 0; i < ker_dim; ++i)
                        for (long j = 0; j < m; ++j)
                            kerbas[i][j] = p_kerbas[i][m-1-j];
                    image(kerbas, kerbas);
                    // now column_permuted_ker is in upper triangular row echelon form
                    for (long i = 0; i < ker_dim; ++i)
                        for (long j = 0; j < m; ++j)
                            p_kerbas[i][j] = kerbas[ker_dim-i-1][m-1-j];
                    // and now p_kerbas is the sought lower triangular row echelon kernel

                    // compute the actual pivot indices
                    for (long i = 0; i<ker_dim; ++i)
                    {
                        p_pivind[i] = m-1;
                        while (IsZero(p_kerbas[i][p_pivind[i]]))
                            --p_pivind[i];
                    }
                }
            }


            // up to row permutation, the kernel is in "lower triangular" row
            // echelon form (almost there: we want the non-permuted one)
            // prepare kernel permutation by permuting kernel pivot indices;
            // also record which rows are pivot index in this kernel
            // (note that before this loop, is_pivind is filled with 'false')
            for (long i = 0; i < ker_dim; ++i)
            {
                pivind[i] = p_shift[p_pivind[i]];
                is_pivind[pivind[i]] = true;
            }

            // perm_rows_ker = [0 1 2 ... ker_dim-1]
            perm_rows_ker.resize(ker_dim);
            std::copy_n(iota.begin(), ker_dim, perm_rows_ker.begin());
            // permutation putting the pivot indices pivind in increasing order
            sort(perm_rows_ker.begin(), perm_rows_ker.end(),
                    [&](const long& a, const long& b)->bool
                    {
                    return (pivind[a] < pivind[b]);
                    } );

            // permute rows and columns of kernel back to original order
            kerbas.SetDims(ker_dim,m);
            for (long i = 0; i < ker_dim; ++i)
                for (long j = 0; j < m; ++j)
                    kerbas[i][p_shift[j]] = p_kerbas[perm_rows_ker[i]][j];
#ifdef MBASIS_PROFILE
            t_others += GetWallTime()-t_now;
            t_now = GetWallTime();
#endif

            // Now, update shifted row degree:
            // entries corresponding to kernel pivot indices are kept, others are +1
            // Also, deduce the degree of appbas
            bool deg_updated=false;
            for (long i = 0; i < m; ++i)
                if (not is_pivind[i])
                {
                    ++shift[i];
                    if (not deg_updated && not IsZero(coeffs_appbas[deg_appbas][i]))
                    { ++deg_appbas; deg_updated=true; }
                }

            // this new degree is either unchanged (== coeffs_appbas.length()-1),
            // or is the old one + 1 (== coeffs_appbas.length())
            if (deg_appbas==coeffs_appbas.length())
            {
                coeffs_appbas.SetLength(deg_appbas+1);
                coeffs_appbas[deg_appbas].SetDims(m, m);
            }
#ifdef MBASIS_PROFILE
            t_others += GetWallTime()-t_now;
            t_now = GetWallTime();
#endif

            // Update approximant basis

            // Submatrix of rows corresponding to pivind are replaced by
            // kerbas*coeffs_appbas (note: these rows currently have degree
            // at most deg_appbas)
            // TODO possible small improvement for uniform shift: these rows
            // have degree less than deg_appbas, in this case (and deg_appbas
            // is reached on the diagonal, among the pivot degrees)
            for (long d = 0; d <= deg_appbas; ++d)
            {
                mul(kerapp, kerbas, coeffs_appbas[d]);
                for (long i = 0; i < ker_dim; ++i)
                    coeffs_appbas[d][pivind[perm_rows_ker[i]]].swap(kerapp[i]);
            }

            // rows with !is_pivind are multiplied by X (note: these rows
            // currently have degree less than deg_appbas)
            for (long d = deg_appbas-1; d >= 0; --d)
                for (long i = 0; i < m; ++i)
                    if (not is_pivind[i])
                        coeffs_appbas[d+1][i].swap(coeffs_appbas[d][i]);
            // Note: after this, the row coeffs_appbas[0][i] is zero
#ifdef MBASIS_PROFILE
            t_appbas += GetWallTime()-t_now;
            t_now = GetWallTime();
#endif
            // Find next residual: coefficient of degree ord in appbas*pmat
            // (this is not necessary if ord==order, since in this case
            // we have finished: appbas*pmat is zero mod X^order)
            // Note: at this point, residuals==0
            if (ord<order)
            {
                long dmin=std::max<long>(0,ord-coeffs_pmat.length()+1);
                for (long d = dmin; d < deg_appbas+1; ++d) // we have deg_appbas <= ord
                {
                    mul(res_coeff, coeffs_appbas[d], coeffs_pmat[ord-d]);
                    add(residuals, residuals, res_coeff);
                }
#ifdef MBASIS_PROFILE
                t_residual += GetWallTime()-t_now;
                t_now = GetWallTime();
#endif
                // Restore is_pivind to all false, as it should be at the beginning of
                // the iteration
                for (long i = 0; i < ker_dim; ++i)
                    is_pivind[pivind[i]] = false;
#ifdef MBASIS_PROFILE
                t_others += GetWallTime()-t_now;
#endif
            }
        }
    }

#ifdef MBASIS_PROFILE
    t_now = GetWallTime();
#endif
    // Convert approximant basis to polynomial matrix representation
    // before that, change to Popov form
    // (this is the only difference with the "non-generic non-Popov" variant of
    // that algorithm)
    // retrieve -rdeg - leading matrix  (only valid in this generic case!)
    Mat<zz_p> L = coeffs_appbas[deg_appbas];
    for (long i = 0; i < pmat.NumRows(); i++)
        if (L[i][i] == 0)
            L[i][i] = 1;
    inv(L, L);
    for (long k = 0; k < deg_appbas; k++)
        mul(coeffs_appbas[k], L, coeffs_appbas[k]);
    for (long i = 0; i < pmat.NumRows(); i++)
        for (long j = 0; j < i; j++)
            coeffs_appbas[deg_appbas][i][j] = 0;
    appbas = conv(coeffs_appbas);
#ifdef MBASIS_PROFILE
    t_others += GetWallTime()-t_now;
#endif
#ifdef MBASIS_PROFILE
    double t_total = t_residual + t_appbas + t_kernel + t_others;
    std::cout << "~~mbasis_rescomp~~\t (residuals,basis,kernel,others): \t ";
    std::cout << t_residual/t_total << "," << t_appbas/t_total << "," <<
    t_kernel/t_total << "," << t_others/t_total << std::endl;
#endif
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
