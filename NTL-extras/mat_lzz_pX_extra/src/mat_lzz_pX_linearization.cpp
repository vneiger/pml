#include "mat_lzz_pX_linearization.h"
#include "mat_lzz_pX_multiply.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* horizontal join                                            */
/* requires a.NumRows() == b.NumRows()                        */
/*------------------------------------------------------------*/
void horizontal_join(Mat<zz_pX>& c, const Mat<zz_pX>& a, const Mat<zz_pX>& b)
{
    const long r = a.NumRows();
    if (r != b.NumRows())
        LogicError("Dimension mismatch for horizontal join");
    const long ca = a.NumCols();
    const long cb = b.NumCols();

    Mat<zz_pX> c2(INIT_SIZE, r, ca + cb);
    for (long i = 0; i < r; ++i)
        for (long j = 0; j < ca; ++j)
            c2[i][j] = a[i][j];
    for (long i = 0; i < r; ++i)
        for (long j = 0; j < cb; ++j)
            c2[i][j+ca] = b[i][j];
    c.swap(c2);
}

/*------------------------------------------------------------*/
/* collapses s consecutive columns of a into one column of c  */
/* let t=a.NumCols(). For i=0..t/s-1, the i-th column of c is */
/* a[i*s] + x^d a[i*s+1] + ... + x^{(s-1)*d} a[i*s+s-1)]      */
/* requires that s divides t exactly                          */
/*------------------------------------------------------------*/
void collapse_consecutive_columns(Mat<zz_pX>& c, const Mat<zz_pX>& a, long d, long s)
{
    if (&c == &a)
    {
        c = collapse_consecutive_columns(a, d, s);
        return;
    }
    long rho = a.NumRows();
    long t = a.NumCols();
    long r = t / s;
    if (t != r * s)
        LogicError("Bad parameters for collapse_consecutive_columns");

    c.SetDims(rho, r);
    for (long j = 0; j < rho; j++)
        for (long i = 0; i < r; i++)
        {
            c[j][i] = 0;
            c[j][i].SetLength(d * s);
            for (long k = 0; k < s; k++)
                for (long ell = 0; ell < d; ell++)
                    SetCoeff(c[j][i], k*d + ell, coeff(a[j][i*s + k], ell));
        }
}

/*------------------------------------------------------------*/
/* collapses columns with stepsize s of a into a column of c  */
/* let t=a.NumCols(). For i=0..s-1, the i-th column of c is   */
/* a[i] + x^d a[i+s] + ... + x^{(t/s-1)*d} a[i+(t/s-1)*s)]    */
/* requires that s divides t exactly                          */
/*------------------------------------------------------------*/
void collapse_nonconsecutive_columns(Mat<zz_pX>& c, const Mat<zz_pX>& a, long d, long s)
{
    if (&c == &a)
    {
        c = collapse_nonconsecutive_columns(a, d, s);
        return;
    }
    long rho = a.NumRows();
    long t = a.NumCols();
    long r = t / s;
    if (t != r * s)
        LogicError("Bad parameters for collapse_nonconsecutive_columns");

    c.SetDims(rho, s);
    for (long j = 0; j < rho; j++)
        for (long i = 0; i < s; i++)
        {
            c[j][i] = 0;
            c[j][i].SetLength(d * r);
            for (long k = 0; k < r; k++)
                for (long ell = 0; ell < d; ell++)
                    SetCoeff(c[j][i], k*d + ell, coeff(a[j][i + k*s], ell));
        }
}








// TODO comment with description?
// reference: Gupta et al 2012.
// For the precise choice of parameters: Labahn - Neiger - Zhou, 2017, Definition 5.5.
VecLong column_partial_linearization(
                                     Mat<zz_pX> & parlin, 
                                     const Mat<zz_pX> & pmat, 
                                     const VecLong & parlin_degree,
                                     const VecLong & target_degree
                                    )
{
    // for each column of pmat, compute the corresponding column indices in the
    // partial linearization
    // column j of pmat will be expanded into columns parlin_cols[j] to
    // parlin_cols[j+1]-1 of parlin
    VecLong parlin_cols(pmat.NumCols()+1);
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
VecLong column_partial_linearization(
                                     Mat<zz_pX> & parlin, 
                                     const Mat<zz_pX> & pmat, 
                                     const VecLong & parlin_degree,
                                     const VecLong & target_degree,
                                     const long d_inf
                                    )
{
    // for each column of pmat, compute the corresponding column indices in the
    // partial linearization
    // column j of pmat will be expanded into columns parlin_cols[j] to
    // parlin_cols[j+1]-1 of parlin
    VecLong parlin_cols(pmat.NumCols()+1);
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
                           const VecLong & column_degree
                          )
{

    // compute the column partial linearization of b
    Mat<zz_pX> b_parlin;
    VecLong parlin_cols = column_partial_linearization(b_parlin,b,parlin_degree,column_degree);

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
                                 const VecLong & target_degree,
                                 long dA,
                                 long dB
                                )
{
    long d_inf = dA-deg(a);
    VecLong parlin_degrees(b.NumCols(),parlin_degree);
    // compute the column partial linearization of b
    Mat<zz_pX> b_parlin;
    VecLong parlin_cols = column_partial_linearization(b_parlin,b,parlin_degrees,target_degree,d_inf);

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
