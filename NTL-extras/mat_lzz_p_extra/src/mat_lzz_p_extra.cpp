#include "mat_lzz_p_extra.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* returns                                                    */
/* [0 0 0 c]                                                  */
/* [1 0 0 0]                                                  */
/* [0 1 0 0]                                                  */
/* [0 0 1 0]                                                  */
/* in size n                                                  */
/*------------------------------------------------------------*/
Mat<zz_p> Z_lzz_p(long n, const zz_p& c)
{
    Mat<zz_p> res(INIT_SIZE, n, n);
    for (long i = 1; i < n; ++i)
        set(res[i][i-1]);
    res[0][n-1] = c;
    return res;
}

/*------------------------------------------------------------*/
/* returns                                                    */
/* [0 0 0 1]                                                  */
/* [0 0 1 0]                                                  */
/* [0 1 0 0]                                                  */
/* [1 0 0 0]                                                  */
/* in size n                                                  */
/*------------------------------------------------------------*/
Mat<zz_p> J_lzz_p(long n)
{
    Mat<zz_p> res(INIT_SIZE, n, n);
    res.SetDims(n, n);
    for (long i = 0; i < n; ++i)
        set(res[i][n-1-i]);
    return res;
}

/*------------------------------------------------------------*/
/* diagonal matrix with diagonal d                            */
/*------------------------------------------------------------*/
Mat<zz_p> diagonal_matrix(const Vec<zz_p> & d)
{
    const long n = d.length();
    Mat<zz_p> res(INIT_SIZE, n, n);
    for (long i = 0; i < n; ++i)
        res[i][i] = d[i];
    return res;
}

/*------------------------------------------------------------*/
/* clears a matrix window starting at (r_offset,c_offset) and */
/* with dimensions nrows x ncols                              */
/*------------------------------------------------------------*/
void clear(Mat<zz_p> & mat, long r_offset, long c_offset, long nrows, long ncols)
{
    const long r_end = std::min(r_offset+nrows,mat.NumRows());
    const long c_end = std::min(c_offset+ncols,mat.NumCols());
    for (long i=r_offset; i<r_end; ++i)
        for (long j=c_offset; j<c_end; ++j)
            clear(mat[i][j]);
}


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
