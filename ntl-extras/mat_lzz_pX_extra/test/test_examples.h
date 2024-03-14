#ifndef __EXAMPLES_H__
#define __EXAMPLES_H__
#include <NTL/BasicThreadPool.h>
#include <NTL/lzz_pX.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <iomanip>
#include <random>
#include <algorithm>

#include "util.h"
#include "mat_lzz_pX_multiply.h"

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
    random(lvec, 1, r-1,3);

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

std::vector<VecLong> build_test_shifts(long rdim, long cdim, long order, bool bigamp_shifts)
{
    std::vector<VecLong> shifts;

    // uniform [0,...,0]
    shifts.emplace_back(rdim,0);

    // increasing [0,1,2,..,rdim-1]
    shifts.emplace_back(rdim);
    std::iota(shifts.back().begin(), shifts.back().end(),0);

    // decreasing [rdim,..,3,2,1]
    shifts.emplace_back(rdim);
    long i = rdim;
    for (auto it = shifts.back().begin(); it != shifts.back().end(); ++it, --i)
        *it = i;

    // random shuffle of [0,1,...,rdim-1]
    shifts.emplace_back(rdim);
    std::iota(shifts.back().begin(), shifts.back().end(), 0);
    std::shuffle(shifts.back().begin(), shifts.back().end(), std::mt19937{std::random_device{}()});

    if (bigamp_shifts)
    {
        // Hermite shift
        shifts.emplace_back(rdim);
        i = 0;
        for (auto it = shifts.back().begin(); it != shifts.back().end(); ++it, i+=cdim*order)
            *it = i;
        
        // reverse Hermite shift
        shifts.emplace_back(rdim);
        for (auto it = shifts.back().begin(); it != shifts.back().end(); ++it, i-=cdim*order)
            *it = i;
        
        // "plateau" shift   [0 ... 0  inf ... inf]
        shifts.emplace_back(rdim);
        auto it = shifts.back().begin();
        std::advance(it, rdim/2);
        for (; it != shifts.back().end(); ++it)
            *it = cdim*order;
    }

    return shifts;
}

// build test matrices and test (row-wise) shifts
// (this is not optimized in terms of memory..)
std::pair<std::vector<Mat<zz_pX>>, std::vector<std::vector<VecLong>>>
build_test_examples(bool bigamp_shifts=false)
{
    // dimensions we will try
    std::vector<long> rdims = {1, 2, 3, 5, 10, 15, 23};
    std::vector<long> cdims = {1, 2, 3, 5, 11, 17, 21};

    // degrees we will try
    std::vector<long> degs = {0, 1, 2, 3, 4, 5, 10, 15, 25, 50, 75, 100, 125};

    // TODO later, a few unbalanced row degrees/column degrees
    //std::vector<long> rdegs = ..
    //std::vector<long> cdegs = ..

    std::vector<Mat<zz_pX>> test_matrices;
    std::vector<std::vector<VecLong>> test_shifts;

    // zero matrices
    for (long rdim : rdims)
        for (long cdim : cdims)
        {
            auto mat = Mat<zz_pX>();
            mat.SetDims(rdim,cdim);
            test_matrices.push_back(mat);
            test_shifts.emplace_back(build_test_shifts(rdim,cdim,2,bigamp_shifts));
        }

    // random matrices, uniform degree
    for (long rdim : rdims)
        for (long cdim : cdims)
            for (long d : degs)
            {
                test_matrices.push_back(random_mat_zz_pX(rdim, cdim, d));
                test_shifts.emplace_back(build_test_shifts(rdim,cdim,d+3,bigamp_shifts));
            }



    // rank deficient square matrices
    // (take random square of dim > 1 and define one column to be random combination of others
    for (long rdim : rdims)
        for (long d: degs)
        {
            Mat<zz_pX> tmp;
            row_rank_deficient_mat(tmp,rdim,rdim,d);
            test_matrices.push_back(tmp);
            test_shifts.emplace_back(build_test_shifts(rdim,rdim,d+2,bigamp_shifts));
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
                    test_matrices.push_back(tmp);
                    test_shifts.emplace_back(build_test_shifts(rdim,cdim,d+1,bigamp_shifts));
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
                    test_matrices.push_back(tmp);
                    test_shifts.emplace_back(build_test_shifts(rdim,cdim,d+7,bigamp_shifts));
                }
            }

    return std::make_pair(test_matrices,test_shifts);
}


#endif // __EXAMPLE_H__

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

