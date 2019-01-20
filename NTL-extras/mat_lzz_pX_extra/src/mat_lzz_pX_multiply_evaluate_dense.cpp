#include <NTL/matrix.h>
#include <NTL/mat_lzz_p.h>
#include <NTL/lzz_pX.h>
#include <NTL/BasicThreadPool.h>

#include "util.h"
#include "lzz_p_extra.h"
#include "mat_lzz_pX_extra.h"
#include "lzz_pX_CRT.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* computes the matrices for evaluation and interpolation     */
/*------------------------------------------------------------*/
void vandermonde(Mat<zz_p>& vdm1, Mat<zz_p>& vdm2, Mat<zz_p>& inv_vdm, long d1, long d2)
{
    // sizes = degree + 1
    const long s1 = d1 + 1;
    const long s2 = d2 + 1;
    // nb points, which is sum of degrees + 1
    const long nb_points = d1 + d2 + 1;

    // will be Vandermonde matrix with nb_points, chosen as 0, 1, .., nb_points-1
    // row i contains 1, i, i^2, .., i^{sk-1} for k=1,2
    vdm1.SetDims(nb_points, s1);
    vdm2.SetDims(nb_points, s2);

    // vdm: square Vandermonde matrix with nb_points
    // points chosen as 0, 1, .., nb_points-1
    // --> row i contains 1, i, i^2, .., i^{nb_points-1}
    // vdm1: nb_points x s1 submatrix of vdm
    // vdm2: nb_points x s2 submatrix of vdm
    Mat<zz_p> vdm;
    vdm.SetDims(nb_points, nb_points);
    vdm[0][0].LoopHole() = 1;
    vdm1[0][0].LoopHole() = 1;
    vdm2[0][0].LoopHole() = 1;
    for (long i = 1; i < nb_points; ++i)
    {
        const zz_p p1(i);
        vdm[i][0].LoopHole() = 1;
        vdm1[i][0].LoopHole() = 1;
        vdm2[i][0].LoopHole() = 1;
        for (long j = 1; j < nb_points; ++j)
        {
            mul(vdm[i][j], vdm[i][j-1], p1);
            if (j<s1)
                vdm1[i][j] = vdm[i][j];
            if (j<s2)
                vdm2[i][j] = vdm[i][j];
        }
    }

    // inv_vdm is the inverse of vdm
    inv(inv_vdm, vdm);
}

/*------------------------------------------------------------*/
/* c = a*b                                                    */
/* output may alias input; c does not have to be zero matrix  */
/* matrix multiplication using the algorithm of Doliskani et al. */
/* uses matrix multiplication for evaluation and interpolation*/
/*------------------------------------------------------------*/
void multiply_evaluate_dense(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b)
{
    long dA = deg(a);
    long dB = deg(b);
    long min_dAdB = std::min<long>(dA,dB);
    long m = a.NumRows();
    long n = a.NumCols();
    long p = b.NumCols();
    long ell;

    Mat<zz_p> tmp_mat(INIT_SIZE, dA+1, m * n);
    Mat<zz_p> valA, valB, valC;
    Mat<zz_p> vA, vB, iv;
    Mat<zz_p> valAp, valBp, valCp;

    double t = GetWallTime();
    vandermonde(vA, vB, iv, dA, dB);
    const long nb_points = vA.NumRows();
    t = GetWallTime()-t;
    std::cout << std::endl << "Build vdmd mat: " << t << std::endl;

    t = GetWallTime();

    // evaluation of matrix a:
    // build tmp_mat, whose column ell = i*n + j is the coefficient vector of
    // a[i][j] (padded with zeroes up to length dA+1 if necessary)
    ell = 0;
    for (long i = 0; i < m; ++i)
        for (long j = 0; j < n; ++j, ++ell)
            for (long k = 0; k <= deg(a[i][j]); ++k)
                tmp_mat[k][ell] = a[i][j][k];
    // note: d = deg(a[i][j]) is -1 if a[i][j] == 0
    // all non-touched entries already zero since tmp_mat was initialized as zero

    // valA: column ell = i*n + j contains the evaluations of a[i][j]
    mul(valA, vA, tmp_mat);

    // evaluation of matrix b:
    // build tmp_mat, whose column ell = i*n + j is the coefficient vector of
    // b[i][j] (padded with zeroes up to length dA+1 if necessary)
    tmp_mat.SetDims(dB+1, n * p);
    ell = 0;
    for (long i = 0; i < n; ++i)
        for (long j = 0; j < p; ++j, ++ell)
        {
            long d = deg(b[i][j]); // -1 if b[i][j] == 0
            for (long k = 0; k <= d; ++k)
                tmp_mat[k][ell] = b[i][j][k];
            // make sure remaining entries are zero
            //    those for k<=dB, k>dA are (if any),
            //    those for d < k <= min(dA,dB) might not be
            for (long k = d+1; k <= min_dAdB; ++k)
                clear(tmp_mat[k][ell]);
        }

    // valB: column ell = i*n + j contains the evaluations of b[i][j]
    mul(valB, vB, tmp_mat);

    t = GetWallTime()-t;
    std::cout << "Evals of a and b: " << t << std::endl;

    t = GetWallTime();
    // perform the pointwise products
    valAp.SetDims(m, n);
    valBp.SetDims(n, p);
    valC.SetDims(nb_points, m * p);
    for (long i = 0; i < nb_points; ++i)
    {
        // a evaluated at point i
        ell = 0;
        for (long u = 0; u < m; ++u)
            for (long v = 0; v < n; ++v, ++ell)
                valAp[u][v] = valA[i][ell];

        // b evaluated at point i
        ell = 0;
        for (long u = 0; u < n; ++u)
            for (long v = 0; v < p; ++v, ++ell)
                valBp[u][v] = valB[i][ell];

        // a*b evaluated at point i
        mul(valCp, valAp, valBp);

        // copy this into valC: column ell = i*n + j contains the evaluations
        // of the entry i,j of c = a*b
        ell = 0;
        for (long u = 0; u < m; ++u)
            for (long v = 0; v < p; ++v, ++ell)
                valC[i][ell] = valCp[u][v];
    }
    t = GetWallTime()-t;
    std::cout << "Pointwise prods: " << t << std::endl;

    t = GetWallTime();
    // interpolate to find the entries of c
    mul(tmp_mat, iv, valC);

    // copy to output (reorganize these entries into c)
    c.SetDims(m, p);
    ell = 0;
    for (long u = 0; u < m; ++u)
        for (long v = 0; v < p; ++v, ++ell)
        {
            c[u][v].SetLength(nb_points);
            for (long i = 0; i < nb_points; ++i)
                c[u][v][i] = tmp_mat[i][ell];
            c[u][v].normalize();
        }
    t = GetWallTime()-t;
    std::cout << "Interpolate: " << t << std::endl;
}

/*------------------------------------------------------------*/
/* computes the matrices for evaluation and interpolation     */
/*------------------------------------------------------------*/
void vandermonde2(
                  Mat<zz_p> & vdm1,
                  Mat<zz_p> & vdm2,
                  Mat<zz_p> & inv_vdm,
                  long d1,
                  long d2
                 )
{
    // sizes (for even part)
    const long s1 = d1/2+1;
    const long s2 = d2/2+1;
    // nb points, such that 2*nb_points >= d1+d2+1
    const long nb_points = (d1 + d2) / 2 + 1;

    // vdm: square Vandermonde matrix with nb_points
    // points chosen as the squares of 1, 2, .., nb_points
    // --> row i contains 1, (i+1)^2, (i+1)^4, .., (i+1)^{2*nb_points-2}
    // vdm1: nb_points x s1 submatrix of vdm
    // vdm2: nb_points x s2 submatrix of vdm
    vdm1.SetDims(nb_points, s1);
    vdm2.SetDims(nb_points, s2);
    Mat<zz_p> vdm(INIT_SIZE, nb_points, nb_points);
    for (long i = 0; i < nb_points; ++i)
    {
        const zz_p pt = sqr(to_zz_p(i+1));
        vdm[i][0].LoopHole() = 1;
        vdm1[i][0].LoopHole() = 1;
        vdm2[i][0].LoopHole() = 1;
        for (long j = 1; j < nb_points; ++j)
        {
            mul(vdm[i][j], vdm[i][j-1], pt);
            if (j<s1)
                vdm1[i][j] = vdm[i][j];
            if (j<s2)
                vdm2[i][j] = vdm[i][j];
        }
    }

    // inv_vdm is the inverse of vdm
    inv(inv_vdm, vdm);
}

/*------------------------------------------------------------*/
/* c = a*b                                                    */
/* output may alias input; c does not have to be zero matrix  */
/* matrix multiplication using the algorithm of Doliskani et al. */
/* uses matrix multiplication for evaluation and interpolation*/
/*------------------------------------------------------------*/
void multiply_evaluate_dense2(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b)
{
    // dimensions and degrees
    const long dA = deg(a);
    const long dB = deg(b);
    const long sA = dA/2+1;
    const long sB = dB/2+1;
    const long m = a.NumRows();
    const long n = a.NumCols();
    const long p = b.NumCols();

    long ell;

    Mat<zz_p> vA, vB, iv;

    double t = GetWallTime();
    vandermonde2(vA, vB, iv, dA, dB);
    t = GetWallTime()-t;
    std::cout << std::endl << "Build vdmd mat: " << t << std::endl;
    const long nb_points = vA.NumRows();

    t = GetWallTime();
    // evaluation of matrix a (even part):
    // build tmp_mat, whose column ell = i*n + j is the coefficient vector of
    // the even part of a[i][j] (padded with zeroes up to length dA/2+1 if
    // necessary)
    Mat<zz_p> tmp_mat;
    tmp_mat.SetDims(sA, m * n);
    ell = 0;
    for (long i = 0; i < m; ++i)
        for (long j = 0; j < n; ++j, ++ell)
            for (long k = 0; 2*k <= deg(a[i][j]); ++k)
                tmp_mat[k][ell] = a[i][j][2*k];

    // note: d = deg(a[i][j]) is -1 if a[i][j] == 0
    // all non-touched entries already zero since tmp_mat was initialized as zero

    // valAeven: column ell = i*n + j contains the evaluations of even part of a[i][j]
    Mat<zz_p> valAeven;
    mul(valAeven, vA, tmp_mat);

    // evaluation of matrix a (odd part):
    // build tmp_mat, whose column ell = i*n + j is the coefficient vector of
    // the odd part of a[i][j] (padded with zeroes up to length (dA-1)/2+1 if
    // necessary)
    // --> in fact we keep dimension dA/2+1, even though this should be
    // (dA-1)/2+1; this is slightly non-optimal but this way we avoid storing
    // four matrices for vAeven, vAodd, vBeven, vBodd
    ell = 0;
    for (long i = 0; i < m; ++i)
        for (long j = 0; j < n; ++j, ++ell)
        {
            const long d = deg(a[i][j]);
            long k = 0;
            for (; 2*k < d; ++k)
                tmp_mat[k][ell] = a[i][j][2*k+1];
            for (; k < sA; ++k)
                clear(tmp_mat[k][ell]);
        }

    // valAodd: column ell = i*n + j contains the evaluations of odd part of a[i][j]
    Mat<zz_p> valAodd;
    mul(valAodd, vA, tmp_mat);

    // evaluation of matrix b (even part):
    // build tmp_mat, whose column ell = i*n + j is the coefficient vector of
    // the even part of b[i][j] (padded with zeroes up to length dB/2+1 if
    // necessary)
    tmp_mat.SetDims(sB, n * p);
    ell = 0;
    for (long i = 0; i < n; ++i)
        for (long j = 0; j < p; ++j, ++ell)
        {
            const long d = deg(b[i][j]); // -1 if b[i][j] == 0
            long k = 0;
            for (; 2*k <= d; ++k)
                tmp_mat[k][ell] = b[i][j][2*k];
            for (; k < sB; ++k)
                clear(tmp_mat[k][ell]);
        }

    // valBeven: column ell = i*n + j contains the evaluations of the even part
    // of b[i][j]
    Mat<zz_p> valBeven;
    mul(valBeven, vB, tmp_mat);

    // evaluation of matrix b (odd part):
    // build tmp_mat, whose column ell = i*n + j is the coefficient vector of
    // the odd part of b[i][j] (padded with zeroes up to length dB/2+1 if
    // necessary; as above, we use dimension dB/2+1 although we know it should
    // be (dB-1)/2+1)
    ell = 0;
    for (long i = 0; i < n; ++i)
        for (long j = 0; j < p; ++j, ++ell)
        {
            const long d = deg(b[i][j]); // -1 if b[i][j] == 0
            long k = 0;
            for (; 2*k < d; ++k)
                tmp_mat[k][ell] = b[i][j][2*k+1];
            for (; k < sB; ++k)
                clear(tmp_mat[k][ell]);
        }

    // valBodd: column ell = i*n + j contains the evaluations of the odd part
    // of b[i][j]
    Mat<zz_p> valBodd;
    mul(valBodd, vB, tmp_mat);

    t = GetWallTime()-t;
    std::cout << "Evals of a and b: " << t << std::endl;

    t = GetWallTime();
    // perform the pointwise products for the points 1, 2, ..., nb_points
    Mat<zz_p> valAp(INIT_SIZE, m, n);
    Mat<zz_p> valBp(INIT_SIZE, n, p);
    Mat<zz_p> valCeven(INIT_SIZE, nb_points, m * p);
    Mat<zz_p> valCp;
    Mat<zz_p> valAs(valAeven);
    for (long i = 0; i < nb_points; ++i)
        valAs[i] += (i+1)*valAodd[i];
    Mat<zz_p> valBs(valBeven);
    for (long i = 0; i < nb_points; ++i)
        valBs[i] += (i+1)*valBodd[i];
    for (long i = 0; i < nb_points; ++i)
    {
        zz_p pt(i+1);
        // a evaluated at point i (which is i+1)
        ell = 0;
        for (long u = 0; u < m; ++u)
            for (long v = 0; v < n; ++v, ++ell)
            {
                // valAp[u][v] = valAeven[i][ell] + pt[i] * valAodd[i][ell];
                //mul(valAp[u][v], pt, valAodd[i][ell]);
                //add(valAp[u][v], valAp[u][v], valAeven[i][ell]);
                valAp[u][v] = valAs[i][ell];
            }

        // b evaluated at point i (which is i+1)
        ell = 0;
        for (long u = 0; u < n; ++u)
            for (long v = 0; v < p; ++v, ++ell)
            {
                // valBp[u][v] = valBeven[i][ell] + pt[i] * valBodd[i][ell];
                //mul(valBp[u][v], pt, valBodd[i][ell]);
                //add(valBp[u][v], valBp[u][v], valBeven[i][ell]);
                valBp[u][v] = valBs[i][ell];
            }

        //Mat<zz_p> tmp_eval = eval(a, to_zz_p(i+1));
        //std::cout << "computed(A,pos):\n" << valAp << std::endl;
        //std::cout << "actual(A,pos):\n" << tmp_eval << std::endl;
        //tmp_eval = eval(b, to_zz_p(i+1));
        //std::cout << "computed(B,pos):\n" << valBp << std::endl;
        //std::cout << "actual(B,pos):\n" << tmp_eval << std::endl;

        // a*b evaluated at point i (which is i+1)
        mul(valCp, valAp, valBp);

        // copy this into valC: column ell = i*n + j contains the evaluations
        // of the entry i,j of c = a*b
        ell = 0;
        for (long u = 0; u < m; ++u)
            for (long v = 0; v < p; ++v, ++ell)
                valCeven[i][ell] = valCp[u][v];
    }

    // perform the pointwise products for the points -1, -2, ..., -nb_points
    Mat<zz_p> valCodd(INIT_SIZE, nb_points, m * p);
    for (long i = 0; i < nb_points; ++i)
        valAs[i] -= 2*(i+1)*valAodd[i];
    for (long i = 0; i < nb_points; ++i)
        valBs[i] -= 2*(i+1)*valBodd[i];
    for (long i = 0; i < nb_points; ++i)
    {
        zz_p pt(-i-1);
        // a evaluated at point nb_points+i (which is -i-1)
        ell = 0;
        for (long u = 0; u < m; ++u)
            for (long v = 0; v < n; ++v, ++ell)
            {
                // valAp[u][v] = valAeven[i][ell] + pt[i] * valAodd[i][ell];
                //mul(valAp[u][v], pt, valAodd[i][ell]);
                //add(valAp[u][v], valAp[u][v], valAeven[i][ell]);
                valAp[u][v] = valAs[i][ell];
            }

        // b evaluated at point nb_points+i (which is -i-1)
        ell = 0;
        for (long u = 0; u < n; ++u)
            for (long v = 0; v < p; ++v, ++ell)
            {
                // valBp[u][v] = valBeven[i][ell] + pt[i] * valBodd[i][ell];
                //mul(valBp[u][v], pt, valBodd[i][ell]);
                //add(valBp[u][v], valBp[u][v], valBeven[i][ell]);
                valBp[u][v] = valBs[i][ell];
            }

        //Mat<zz_p> tmp_eval = eval(a, -to_zz_p(i+1));
        //std::cout << "computed(A,neg):\n" << valAp << std::endl;
        //std::cout << "actual(A,neg):\n" << tmp_eval << std::endl;
        //tmp_eval = eval(b, -to_zz_p(i+1));
        //std::cout << "computed(B,neg):\n" << valBp << std::endl;
        //std::cout << "actual(B,neg):\n" << tmp_eval << std::endl;

        // a*b evaluated at point nb_points+i (which is -i-1)
        mul(valCp, valAp, valBp);

        // copy this into valC: column ell = i*n + j contains the evaluations
        // of the entry i,j of c = a*b
        ell = 0;
        for (long u = 0; u < m; ++u)
            for (long v = 0; v < p; ++v, ++ell)
                valCodd[i][ell] = valCp[u][v];

        //Mat<zz_pX> ab; multiply(ab, a,b);
        //tmp_eval = eval(ab, -to_zz_p(i+1));
        //std::cout << "Before, actual:\n" << tmp_eval << std::endl;
        //std::cout << "Before, computed:\n" << valCodd[i] << std::endl;
    }

    t = GetWallTime()-t;
    std::cout << "Pointwise prods: " << t << std::endl;

    t = GetWallTime();

    // valCeven = (even + odd)/2
    // valCodd = (even - odd)/2x
    for (long i = 0; i < nb_points; ++i)
        for (long j = 0; j < m*p; ++j)
        {
            // TODO precompute inverse of 2
            zz_p tmp;
            add(tmp, valCeven[i][j], valCodd[i][j]);
            sub(valCodd[i][j], valCeven[i][j], valCodd[i][j]);
            div(valCeven[i][j], tmp, 2);
            mul(tmp, 2, to_zz_p(i+1));
            div(valCodd[i][j], valCodd[i][j], tmp);
        }

    //Mat<zz_pX> even(INIT_SIZE, m, p);
    //Mat<zz_pX> ab; multiply(ab, a,b);
    //for (long u = 0; u < m; ++u)
    //{
    //    for (long v = 0; v < p; ++v)
    //    {
    //        even[u][v].SetLength(deg(ab[u][v])/2+1);
    //        for (long k = 0; k <= deg(ab[u][v])/2; ++k)
    //        {
    //            even[u][v][k] = ab[u][v][2*k];
    //        }
    //    }
    //}
    //Mat<zz_p> tmp_eval = eval(even, to_zz_p(4));
    //std::cout << "Actual eval of even:\n" << tmp_eval << std::endl;
    //std::cout << "Computed eval of even:\n" << valCeven[1] << std::endl;

    //Mat<zz_pX> odd(INIT_SIZE, m, p);
    //for (long u = 0; u < m; ++u)
    //{
    //    for (long v = 0; v < p; ++v)
    //    {
    //        for (long k = 0; 2*k+1 <= deg(ab[u][v]); ++k)
    //        {
    //            SetCoeff(odd[u][v], k, ab[u][v][2*k+1]);
    //        }
    //    }
    //}
    //tmp_eval = eval(odd, to_zz_p(4));
    //std::cout << "Actual eval of odd:\n" << tmp_eval << std::endl;
    //std::cout << "Computed eval of odd:\n" << valCodd[1] << std::endl;

    // interpolate to find the even part of c, and copy it into c
    mul(tmp_mat, iv, valCeven);

    //std::cout << "tmp_mat even:\n" << tmp_mat << std::endl;

    c.SetDims(m, p);
    ell = 0;
    for (long u = 0; u < m; ++u)
        for (long v = 0; v < p; ++v, ++ell)
        {
            c[u][v].SetLength(2*nb_points);
            for (long i = 0; i < nb_points; ++i)
                c[u][v][2*i] = tmp_mat[i][ell];
        }

    // interpolate to find the odd part of c, and copy it into c
    mul(tmp_mat, iv, valCodd);
    //std::cout << "tmp_mat odd:\n" << tmp_mat << std::endl;

    ell = 0;
    for (long u = 0; u < m; ++u)
        for (long v = 0; v < p; ++v, ++ell)
        {
            for (long i = 0; i < nb_points; ++i)
                c[u][v][2*i+1] = tmp_mat[i][ell];
            c[u][v].normalize();
        }
    t = GetWallTime()-t;
    std::cout << "Interpolate: " << t << std::endl;
}


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
