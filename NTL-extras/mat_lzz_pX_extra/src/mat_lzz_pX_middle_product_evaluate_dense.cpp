#include "mat_lzz_pX_multiply.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* returns trunc( trunc(a, dA+1)*c div x^dA, dB+1 )           */
/* assumes FFT prime and p large enough                       */
/* output may alias input; b does not have to be zero matrix  */
/* does not use Mat<zz_p> matrix multiplication               */
/*------------------------------------------------------------*/
void middle_product_evaluate_dense(Mat<zz_pX> & b, const Mat<zz_pX> & a, const Mat<zz_pX> & c, long dA, long dB)
{
    long m = a.NumRows();
    long n = a.NumCols();
    long p = c.NumCols();
    long ell;
    
    Mat<zz_p> tmp_mat(INIT_SIZE, dA + 1, m * n);
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
                    tmp_mat[dA-k][ell] = cAij[k];
                }
            }
            for (long k = d+1; k <= dA; k++)
            {
                tmp_mat[dA-k][ell] = 0;
            }
            
            ell++;
        }
    valA = vA * tmp_mat;

    tmp_mat.SetDims(nb_points, n * p);
    ell = 0;
    for (long i = 0; i < n; i++)
        for (long j = 0; j < p; j++)
        {
            long d = deg(c[i][j]);
            if (d >= 0)
            {
                const zz_p * cCij = c[i][j].rep.elts();
                long k;
                for (k = 0; k <= d; k++)  
                {
                    tmp_mat[k][ell] = cCij[k];
                }
            }
            for (long k = d+1; k < nb_points; k++)
            {
                tmp_mat[k][ell] = 0;
            }
            ell++;
        }
    valC = transpose(iv) * tmp_mat;
    
    valAp.SetDims(m, n);
    valCp.SetDims(n, p);
    valB.SetDims(nb_points, m * p);
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
                valCp[u][v] = valC[i][ell++];
            }
        
        valBp = valAp * valCp;
        
        ell = 0;
        for (long u = 0; u < m; u++)
            for (long v = 0; v < p; v++)
            {
                valB[i][ell++] = valBp[u][v];
            }
    }
    tmp_mat = transpose(vB)*valB;
    
    b.SetDims(m, p);
    ell = 0;
    for (long u = 0; u < m; u++)
        for (long v = 0; v < p; v++)
        {
            b[u][v].rep.SetLength(dB + 1);
            zz_p * bb = b[u][v].rep.elts();
            for (long i = 0; i < dB + 1; i++)
                bb[i] = tmp_mat[i][ell];
            b[u][v].normalize();
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
