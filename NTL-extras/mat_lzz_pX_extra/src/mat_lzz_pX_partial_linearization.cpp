#include <NTL/matrix.h>
#include <NTL/mat_lzz_p.h>
#include <NTL/lzz_pX.h>
#include <cmath>
#include <algorithm> // for manipulating std::vector (min, max, ..)

#include "lzz_p_extra.h"
#include "mat_lzz_pX_extra.h"
#include "lzz_pX_CRT.h"

NTL_CLIENT

// TODO later, move this function in a more appropriate place (specific file for linearizations)
// TODO comment with description?
// TODO degree should be vector of column degrees!!
std::vector<long> column_partial_linearization(Mat<zz_pX> &parlin, 
                                  const Mat<zz_pX> &pmat, 
                                  const DegVec & column_degree,
                                  const DegVec & parlin_degree)
{
    // for each column of pmat, compute the corresponding column indices in the
    // partial linearization
    // column j of pmat will be expanded into columns numcols[j] to
    // numcols[j+1]-1 of parlin
    std::vector<long> parlin_cols(pmat.NumCols()+1);
    parlin_cols[0]=0;
    for (long j = 1; j < pmat.NumCols()+1; ++j)
        parlin_cols[j] = parlin_cols[j-1] + ceil((column_degree[j-1]+1.0) / (parlin_degree[j-1]+1));

    // set dimensions of the partial linearization
    parlin.SetDims(pmat.NumRows(), parlin_cols[pmat.NumCols()]);

    // for each column j of the input pmat, expand it into the corresponding
    // columns of the partial linearization parlin
    for (long i=0; i<pmat.NumRows(); ++i)
    {
        for (long j=0; j<pmat.NumCols(); ++j)
        {
            long deg_at = 0;
            for (long jj=parlin_cols[j]; jj<parlin_cols[j+1]; ++jj)
            {
                parlin[i][jj].SetLength(parlin_degree[j]+1);
                for (long d=0; d<=parlin_degree[j]; ++d)
                    parlin[i][jj][d] = coeff(pmat[i][j], deg_at++);
            }
        }
    }

    return parlin_cols;
}

// TODO later, move this function in a more appropriate place (specific file for linearizations)
// TODO comment with description?
// splits F such that each entry has at most deg_sp
void right_parlin_multiply(Mat<zz_pX> &c,
                           const Mat<zz_pX> &a,
                           const Mat<zz_pX> &b,
                           const DegVec & column_degree,
                           const long parlin_degree)
{
    // compute the column partial linearization of b
    Mat<zz_pX> b_parlin;
    std::vector<long> parlin_cols = column_partial_linearization(b_parlin,b,column_degree,parlin_degree);

    // compute the linearized product c_parlin = a * b_parlin
    Mat<zz_pX> c_parlin;
    multiply(c_parlin, a, b_parlin);

    // compress the product back into c
    c.SetDims(a.NumRows(), b.NumCols());
    for (long i=0; i<a.NumRows(); ++i)
    {
        long deg_a = deg(a[i]);
        for (long j=0; j<b.NumCols(); ++j)
        {
            c[i][j].SetLength(deg_a + column_degree[j] + 1);
            long d_init=0;
            for (long jj=parlin_cols[j]; jj<parlin_cols[j+1]; ++jj)
            {
                for (long d=0; d<=deg(c_parlin[i][jj]); ++d)
                    c[i][j][d_init+d] += coeff(c_parlin[i][jj], d);
                d_init += parlin_degree+1;
            }
            c[i][j].normalize();
        }
    }
}

void right_parlin_middle_product(
                           Mat<zz_pX> &c,
                           const Mat<zz_pX> &a,
                           const Mat<zz_pX> &b,
                           const DegVec & column_degree,
                           const long parlin_degree,
                           long dA,
                           long dB
                          )
{
    // compute the column partial linearization of b
    Mat<zz_pX> b_parlin;
    std::vector<long> parlin_cols = column_partial_linearization(b_parlin,b,column_degree,parlin_degree);

    // compute the linearized product c_parlin = a * b_parlin
    Mat<zz_pX> c_parlin;
    multiply(c_parlin, a, b_parlin);

    // compress the product back into c
    c.SetDims(a.NumRows(), b.NumCols());
    for (long i=0; i<a.NumRows(); ++i)
    {
        long deg_a = deg(a[i]);
        for (long j=0; j<b.NumCols(); ++j)
        {
            c[i][j].SetLength(deg_a + column_degree[j] + 1);
            long d_init=0;
            for (long jj=parlin_cols[j]; jj<parlin_cols[j+1]; ++jj)
            {
                for (long d=0; d<=deg(c_parlin[i][jj]); ++d)
                    c[i][j][d_init+d] += coeff(c_parlin[i][jj], d);
                d_init += parlin_degree+1;
            }
            c[i][j].normalize();
        }
    }

    for (long i = 0; i < c.NumRows(); ++i)
        for (long j = 0; j < c.NumCols(); ++j)
        {
            RightShift(c[i][j], c[i][j], dA); 
            trunc(c[i][j], c[i][j], dB);
        }

}


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
