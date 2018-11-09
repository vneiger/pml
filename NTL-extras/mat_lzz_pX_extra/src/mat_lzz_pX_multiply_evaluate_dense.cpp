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
void vandermonde(Mat<zz_p>& small_vdm1, Mat<zz_p>& small_vdm2, Mat<zz_p>& inv_vdm, long d1, long d2)
{
    long s1 = d1 + 1;
    long s2 = d2 + 1;
    long d = d1 + d2;
    long nb_points = d + 1;
    
    Vec<zz_p> points;
    points.SetLength(nb_points);
    small_vdm1.SetDims(nb_points, s1);
    small_vdm2.SetDims(nb_points, s2);
    
    Mat<zz_p> vdm;
    vdm.SetDims(nb_points, nb_points);
    
    for (long i = 0; i < nb_points; i++)
    {
        zz_p p1 = to_zz_p(i);
        zz_p tmp = to_zz_p(1);
        for (long j = 0; j < nb_points; j++)
        {
            vdm[i][j] = tmp;
            if (j < s1)
                small_vdm1[i][j] = tmp;
            if (j < s2)
                small_vdm2[i][j] = tmp;
            tmp *= p1;
        }
    }
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
    long sA = dA + 1;
    long sB = dB + 1;
    long m = a.NumRows();
    long n = a.NumCols();
    long p = b.NumCols();
    long ell;
    
    Mat<zz_p> tmp_mat(INIT_SIZE, sA, m * n);
    Mat<zz_p> valA, valB, valC;
    Mat<zz_p> vA, vB, iv;
    Mat<zz_p> valAp, valBp, valCp;
    
    vandermonde(vA, vB, iv, dA, dB);
    long nb_points = vA.NumRows();
    
    ell = 0;
    for (long i = 0; i < m; i++)
        for (long j = 0; j < n; j++)
        {
            long d = deg(a[i][j]);
            if (d >= 0)
            {
                const zz_p * cAij = a[i][j].rep.elts();
                long k;
                for (k = 0; k <= d; k++)  // k <= d-2 so k+1 <= d-1
                {
                    tmp_mat[k][ell] = cAij[k];
                }
            }
            for (long k = d+1; k <= dA; k++)
            {
                tmp_mat[k][ell] = 0;
            }
            
            ell++;
        }
    valA = vA * tmp_mat;
    
    tmp_mat.SetDims(sB, n * p);
    ell = 0;
    for (long i = 0; i < n; i++)
        for (long j = 0; j < p; j++)
        {
            long d = deg(b[i][j]);
            if (d >= 0)
            {
                const zz_p * cBij = b[i][j].rep.elts();
                long k;
                for (k = 0; k <= d; k++)  // k <= d-2 so k+1 <= d-1
                {
                    tmp_mat[k][ell] = cBij[k];
                }
            }
            for (long k = d+1; k <= dB; k++)
            {
                tmp_mat[k][ell] = 0;
            }
            ell++;
        }
    valB = vB * tmp_mat;
    
    valAp.SetDims(m, n);
    valBp.SetDims(n, p);
    valC.SetDims(nb_points, m * p);
    for (long i = 0; i < nb_points; i++)
    {
        long ell;
        ell = 0;
        for (long u = 0; u < m; u++)
            for (long v = 0; v < n; v++)
            {
                valAp[u][v] = valA[i][ell++];
            }
        
        ell = 0;
        for (long u = 0; u < n; u++)
            for (long v = 0; v < p; v++)
            {
                valBp[u][v] = valB[i][ell++];
            }
        
        valCp = valAp * valBp;
        
        ell = 0;
        for (long u = 0; u < m; u++)
            for (long v = 0; v < p; v++)
            {
                valC[i][ell++] = valCp[u][v];
            }
    }
    tmp_mat = iv*valC;
    
    c.SetDims(m, p);
    ell = 0;
    for (long u = 0; u < m; u++)
        for (long v = 0; v < p; v++)
        {
            c[u][v].rep.SetLength(nb_points);
            zz_p * cc = c[u][v].rep.elts();
            for (long i = 0; i < nb_points; i++)
                cc[i] = tmp_mat[i][ell];
            c[u][v].normalize();
            ell++;
        }
}





// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
