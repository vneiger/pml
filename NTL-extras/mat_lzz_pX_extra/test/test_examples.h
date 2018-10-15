#ifndef __EXAMPLES_H__
#define __EXAMPLES_H__
#include <NTL/BasicThreadPool.h>
#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>

#include "util.h"
#include "mat_lzz_pX_extra.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* Examples of matrices for testing polmat functions          */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

void row_rank_deficient_mat(Mat<zz_pX> &m, const long r, const long c, const long deg)
{
    random(m,r-1,c,deg+1);
    Mat<zz_pX> lvec;
    random(lvec, 1, r-1,1);

    Mat<zz_pX> new_row;
    multiply(new_row, lvec, m);
    m.SetDims(r,c);
    m[r-1] = new_row[0];
}

void col_rank_deficient_mat(Mat<zz_pX> &m, const long r, const long c, const long deg)
{
    random(m,r,c-1,deg+1);
    Mat<zz_pX> rvec;
    random(rvec, c-1, 1, 1);

    Mat<zz_pX> new_col;
    multiply(new_col, m, rvec);
    m.SetDims(r,c);
    for (long i = 0; i < r; i++)
        m[i][c-1] = new_col[i][0];
}

void build(std::vector<Mat<zz_pX>> & example_matrices)
{
    // dimensions we will try
    std::vector<long> rdims = {1, 2, 3, 5, 10, 15, 23};
    std::vector<long> cdims = {1, 2, 3, 5, 11, 17, 21};

    // degrees we will try
    std::vector<long> degs = {0, 1, 2, 3, 4, 5, 10, 15, 25, 50, 100};

    // TODO later, a few unbalanced row degrees/column degrees
    //std::vector<long> rdegs = ..
    //std::vector<long> cdegs = ..

    // zero matrices
    for (long rdim : rdims)
        for (long cdim : cdims)
        {
            auto mat = Mat<zz_pX>();
            mat.SetDims(rdim,cdim);
            example_matrices.push_back(mat);
        }

    // random matrices, uniform degree
    for (long rdim : rdims)
        for (long cdim : cdims)
            for (long d : degs)
            {
                example_matrices.push_back(random_mat_zz_pX(rdim, cdim, d));
            }

    // rank deficient square matrices
    // (take random square of dim > 1 and define one column to be random combination of others
    for (long rdim : rdims)
        for (long d: degs)
        {
            Mat<zz_pX> tmp;
            row_rank_deficient_mat(tmp,rdim,rdim,d);
            example_matrices.push_back(tmp);
        }


    // rank deficient rectangular rdim>cdim matrices
    // (take random of cdim > 1 and define one column to be random combination of others)
    for (long rdim: rdims)
        for (long cdim: cdims)
            for (long d: degs)
            {
                Mat<zz_pX> tmp;
                if (rdim > cdim)
                {
                    col_rank_deficient_mat(tmp,rdim,cdim,d);
                    example_matrices.push_back(tmp);
                }
            }

    // rank deficient rectangular rdim<cdim matrices
    // (take random of rdim > 1 and define one row to be random combination of others)
    for (long rdim: rdims)
        for (long cdim: cdims)
            for (long d: degs)
            {
                Mat<zz_pX> tmp;
                if (rdim < cdim)
                {
                    row_rank_deficient_mat(tmp,rdim,cdim,d);
                    example_matrices.push_back(tmp);
                }
            }
}

#endif // __EXAMPLE_H__

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

