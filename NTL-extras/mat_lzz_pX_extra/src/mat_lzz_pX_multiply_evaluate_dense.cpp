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
/* matrix multiplication using the algorithm of Giorgi et al. */
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

    vandermonde(vA, vB, iv, dA, dB);
    long nb_points = vA.NumRows();

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
    // a[i][j] (padded with zeroes up to length dA+1 if necessary)
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
}


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
