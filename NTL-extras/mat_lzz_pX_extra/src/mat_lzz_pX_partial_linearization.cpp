#include <NTL/matrix.h>
#include <NTL/mat_lzz_p.h>
#include <NTL/lzz_pX.h>
#include <algorithm> // for manipulating std::vector (min, max, ..)

#include "mat_lzz_pX_extra.h"
#include "mat_lzz_pX_partial_linearization.h"

NTL_CLIENT

// TODO later, move this function in a more appropriate place (specific file for linearizations)
// TODO comment with description?
// reference: Gupta et al 2012.
// For the precise choice of parameters: Labahn - Neiger - Zhou, 2017, Definition 5.5.
std::vector<long> column_partial_linearization(
                                               Mat<zz_pX> & parlin, 
                                               const Mat<zz_pX> & pmat, 
                                               const DegVec & parlin_degree,
                                               const DegVec & target_degree
                                              )
{
    // for each column of pmat, compute the corresponding column indices in the
    // partial linearization
    // column j of pmat will be expanded into columns parlin_cols[j] to
    // parlin_cols[j+1]-1 of parlin
    std::vector<long> parlin_cols(pmat.NumCols()+1);
    // expand column j-1 into a matrix with all columns but the last having
    // degree <= parlin_degree[j-1]
    // note: if target_degree[j-1] is the degree of column j-1 of pmat,
    // then all expanded columns *including* the last one have degree <=
    // parlin_degree[j-1]
    parlin_cols[0]=0;
    for (long j = 1; j < pmat.NumCols()+1; ++j)
        parlin_cols[j] = parlin_cols[j-1] + ceil( (target_degree[j-1]+1.0) / (parlin_degree[j-1]+1) );

    // set dimensions of the partial linearization
    parlin.SetDims(pmat.NumRows(), parlin_cols[pmat.NumCols()]);

    // for each column j of the input pmat, expand it into the corresponding
    // columns of the partial linearization parlin
    for (long i=0; i<pmat.NumRows(); ++i)
    {
        for (long j=0; j<pmat.NumCols(); ++j)
        {
            long d = 0;
            long jj=parlin_cols[j];
            for (; jj<parlin_cols[j+1]-1; ++jj)
            {
                parlin[i][jj].SetLength(parlin_degree[j]+1);
                for (long dd=0; dd<=parlin_degree[j]; ++dd)
                    parlin[i][jj][dd] = coeff(pmat[i][j], d++);
                parlin[i][jj].normalize();
            }
            const long reached_deg = d;
            parlin[i][jj].SetLength(deg(pmat[i][j]) - reached_deg + 1);
            for (long dd=0; dd<=deg(pmat[i][j]) - reached_deg; ++dd)
                parlin[i][jj][dd] = coeff(pmat[i][j], d++);
        }
    }

    return parlin_cols;
}

// TODO later, move this function in a more appropriate place (specific file for linearizations)
// TODO comment with description?
// reference: Gupta et al 2012.
// For the precise choice of parameters: Labahn - Neiger - Zhou, 2017, Definition 5.5.
std::vector<long> column_partial_linearization(
                                               Mat<zz_pX> & parlin, 
                                               const Mat<zz_pX> & pmat, 
                                               const DegVec & parlin_degree,
                                               const DegVec & target_degree,
                                               const long d_inf
                                              )
{
    // for each column of pmat, compute the corresponding column indices in the
    // partial linearization
    // column j of pmat will be expanded into columns parlin_cols[j] to
    // parlin_cols[j+1]-1 of parlin
    std::vector<long> parlin_cols(pmat.NumCols()+1);
    // expand column j-1 into a matrix with all columns but the last having
    // degree <= parlin_degree[j-1]
    // note: if target_degree[j-1] is the degree of column j-1 of pmat,
    // then all expanded columns *including* the last one have degree <=
    // parlin_degree[j-1]
    parlin_cols[0]=0;
    for (long j = 1; j < pmat.NumCols()+1; ++j)
    {
        long tgt = std::max<long>(0,target_degree[j-1]-d_inf);
        parlin_cols[j] = parlin_cols[j-1] + ceil( (tgt+1.0) / (parlin_degree[j-1]+1) );
    }

    // set dimensions of the partial linearization
    parlin.SetDims(pmat.NumRows(), parlin_cols[pmat.NumCols()]);

    // for each column j of the input pmat, expand it into the corresponding
    // columns of the partial linearization parlin
    for (long i=0; i<pmat.NumRows(); ++i)
    {
        for (long j=0; j<pmat.NumCols(); ++j)
        {
            long d = d_inf;
            long jj = parlin_cols[j];
            for (; jj<parlin_cols[j+1]-1; ++jj)
            {
                parlin[i][jj].SetLength(parlin_degree[j]+1);
                for (long dd=0; dd<=parlin_degree[j]; ++dd)
                    parlin[i][jj][dd] = coeff(pmat[i][j], d++);
                parlin[i][jj].normalize();
            }
            const long reached_deg = d;
            parlin[i][jj].SetLength(deg(pmat[i][j]) - reached_deg + 1);
            for (long dd=0; dd<=deg(pmat[i][j]) - reached_deg; ++dd)
                parlin[i][jj][dd] = coeff(pmat[i][j], d++);
        }
    }

    return parlin_cols;
}

// TODO later, move this function in a more appropriate place (specific file for linearizations)
// TODO comment with description?
// splits F such that each entry has at most deg_sp
void right_parlin_multiply(
                           Mat<zz_pX> &c,
                           const Mat<zz_pX> &a,
                           const Mat<zz_pX> &b,
                           const long parlin_degree,
                           const DegVec & column_degree
                          )
{

    // compute the column partial linearization of b
    Mat<zz_pX> b_parlin;
    std::vector<long> parlin_cols = column_partial_linearization(b_parlin,b,parlin_degree,column_degree);

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

/* returns trunc( a*b div x^dA, dB+1 )           */
// assumes deg(a) < dA+1
void right_parlin_middle_product(
                                 Mat<zz_pX> &c,
                                 const Mat<zz_pX> &a,
                                 const Mat<zz_pX> &b,
                                 const long parlin_degree,
                                 const DegVec & target_degree,
                                 long dA,
                                 long dB
                                )
{
    long d_inf = dA-deg(a);
    DegVec parlin_degrees(b.NumCols(),parlin_degree);
    // compute the column partial linearization of b
    Mat<zz_pX> b_parlin;
    std::vector<long> parlin_cols = column_partial_linearization(b_parlin,b,parlin_degrees,target_degree,d_inf);

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
            c[i][j].SetLength(deg_a + target_degree[j] - d_inf + 1);
            long d_init=0;
            for (long jj=parlin_cols[j]; jj<parlin_cols[j+1]; ++jj)
            {
                for (long d=0; d<=deg(c_parlin[i][jj]); ++d)
                    c[i][j][d_init+d] += coeff(c_parlin[i][jj], d);
                d_init += parlin_degree+1;
            }
            //c[i][j].normalize(); // TODO remove?
            c[i][j] >>= deg(a);  // precompute deg(a)
            trunc(c[i][j], c[i][j], dB+1);
        }
    }
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
