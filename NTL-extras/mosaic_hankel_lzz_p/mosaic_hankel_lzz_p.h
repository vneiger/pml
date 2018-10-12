#ifndef __MOSAIC_HANKEL_LZZ_P_H
#define __MOSAIC_HANKEL_LZZ_P_H

#include <NTL/lzz_pX.h>
#include <NTL/mat_lzz_p.h>

#include "lzz_pX_CRT.h"
#include "cauchy_geometric_lzz_p.h"

NTL_CLIENT

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/* Mosaic Hankel matrices:                            */
/* block matrix where each block is Hankel            */
/*----------------------------------------------------*/
/*----------------------------------------------------*/
class mosaic_hankel_lzz_p
{
private:
    long n, m; // numbers of rows / colums
    long nb, mb; // number of blocks in rows / columns

public:
    // data[i][j] is the block at i-th row, j-th column
    Vec< Vec<hankel_lzz_p> > data;

    /*----------------------------------------------------*/
    /* dummy constructor                                  */
    /*----------------------------------------------------*/
    mosaic_hankel_lzz_p();

    /*----------------------------------------------------*/
    /* copies all data                                    */
    /*----------------------------------------------------*/
    mosaic_hankel_lzz_p(const Vec< Vec<hankel_lzz_p> > & init);

    /*----------------------------------------------------*/
    /* getters                                            */
    /*----------------------------------------------------*/
    long NumRows() const;
    long NumCols() const;
    long NumBlockRows() const;
    long NumBlockCols() const;
    long NumRows_of_block(long i) const;
    long NumCols_of_block(long i) const;

    zz_p operator ()(long i, long j) const
    {
        if (i >= n || j >= m)
            Error("matrix indices out of bounds\n");
        long idx, jdx;

        idx = 0;
        while (i >= data[idx][0].NumRows())
        {
            i -= data[idx][0].NumRows();
            idx++;
        }
        
        jdx = 0;
        while (j >= data[0][jdx].NumCols())
        {
            j -= data[0][jdx].NumCols();
            jdx++;
        }
        return data[idx][jdx](i, j);
    }


    /*----------------------------------------------------*/
    /* access to particular rows and columns              */
    /*----------------------------------------------------*/
    void first_column_of_block(Vec<zz_p>& res, long i) const;
    inline Vec<zz_p> first_column_of_block(long i) const
    {
        Vec<zz_p> res;
        first_column_of_block(res, i);
        return res;
    }

    void last_column_of_block(Vec<zz_p>& res, long i) const;
    inline Vec<zz_p> last_column_of_block(long i) const
    {
        Vec<zz_p> res;
        last_column_of_block(res, i);
        return res;
    }

    void first_row_of_block(Vec<zz_p>& res, long i) const;
    inline Vec<zz_p> first_row_of_block(long i) const
    {
        Vec<zz_p> res;
        first_row_of_block(res, i);
        return res;
    }

    void last_row_of_block(Vec<zz_p>& res, long i) const;
    inline Vec<zz_p> last_row_of_block(long i) const
    {
        Vec<zz_p> res;
        last_row_of_block(res, i);
        return res;
    }

    /*----------------------------------------------------*/
    /* turns M into a dense matrix                        */
    /*----------------------------------------------------*/
    void to_dense(Mat<zz_p>& Mdense) const;
    inline Mat<zz_p> to_dense() const
    { 
        Mat<zz_p> Mdense;
        to_dense(Mdense);
        return Mdense;
    }

    /*----------------------------------------------------*/
    /* right multiplication                               */
    /*----------------------------------------------------*/
    void mul_right(Vec<zz_p>& res, const Vec<zz_p>& input) const;
    Vec<zz_p> mul_right(const Vec<zz_p>& input) const
    {
        Vec<zz_p> res;
        mul_right(res, input);
        return res;
    }

    void mul_right(Mat<zz_p>& res, const Mat<zz_p>& input) const;
    Mat<zz_p> mul_right(const Mat<zz_p>& input) const
    {
        Mat<zz_p> res;
        mul_right(res, input);
        return res;
    }

    /*----------------------------------------------------*/
    /* left multiplication                                */
    /*----------------------------------------------------*/
    void mul_left(Vec<zz_p>& res, const Vec<zz_p>& input) const;
    Vec<zz_p> mul_left(const Vec<zz_p>& input) const
    {
        Vec<zz_p> res;
        mul_left(res, input);
        return res;
    }

    void mul_left(Mat<zz_p>& res, const Mat<zz_p>& input) const;
    Mat<zz_p> mul_left(const Mat<zz_p>& input) const
    {
        Mat<zz_p> res;
        mul_left(res, input);
        return res;
    }
};

/*----------------------------------------------------*/
/*----------------------------------------------------*/
/* Mosaic Toeplitz matrices:                          */
/* block matrix where each block is Toeplitz          */
/*----------------------------------------------------*/
/*----------------------------------------------------*/
class mosaic_toeplitz_lzz_p{
private:
    long n, m; // numbers of rows / colums
    long nb, mb; // number of blocks in rows / columns

public:
    // data[i][j] is the block at i-th row, j-th column
    Vec< Vec<toeplitz_lzz_p> > data;

    /*----------------------------------------------------*/
    /* dummy constructor                                  */
    /*----------------------------------------------------*/
    mosaic_toeplitz_lzz_p();

    /*----------------------------------------------------*/
    /* copies all data                                    */
    /*----------------------------------------------------*/
    mosaic_toeplitz_lzz_p(const Vec< Vec<toeplitz_lzz_p> > & init);

    /*----------------------------------------------------*/
    /* getters                                            */
    /*----------------------------------------------------*/
    long NumRows() const;
    long NumCols() const;
    long NumBlockRows() const;
    long NumBlockCols() const;
    long NumRows_of_block(long i) const;
    long NumCols_of_block(long i) const;

    zz_p operator ()(long i, long j) const
    {
        if (i >= n || j >= m)
            Error("matrix indices out of bounds\n");
        long idx, jdx;

        idx = 0;
        while (i >= data[idx][0].NumRows())
        {
            i -= data[idx][0].NumRows();
            idx++;
        }
        
        jdx = 0;
        while (j >= data[0][jdx].NumCols())
        {
            j -= data[0][jdx].NumCols();
            jdx++;
        }
        return data[idx][jdx](i, j);
    }


    /*----------------------------------------------------*/
    /* access to particular rows and columns              */
    /*----------------------------------------------------*/
    void first_column_of_block(Vec<zz_p>& res, long i) const;
    inline Vec<zz_p> first_column_of_block(long i) const
    {
        Vec<zz_p> res;
        first_column_of_block(res, i);
        return res;
    }

    void last_column_of_block(Vec<zz_p>& res, long i) const;
    inline Vec<zz_p> last_column_of_block(long i) const
    {
        Vec<zz_p> res;
        last_column_of_block(res, i);
        return res;
    }

    void first_row_of_block(Vec<zz_p>& res, long i) const;
    inline Vec<zz_p> first_row_of_block(long i) const
    {
        Vec<zz_p> res;
        first_row_of_block(res, i);
        return res;
    }

    void last_row_of_block(Vec<zz_p>& res, long i) const;
    inline Vec<zz_p> last_row_of_block(long i) const
    {
        Vec<zz_p> res;
        last_row_of_block(res, i);
        return res;
    }

    /*----------------------------------------------------*/
    /* turns M into a dense matrix                        */
    /*----------------------------------------------------*/
    void to_dense(Mat<zz_p>& Mdense) const;
    inline Mat<zz_p> to_dense() const
    { 
        Mat<zz_p> Mdense;
        to_dense(Mdense);
        return Mdense;
    }

    /*----------------------------------------------------*/
    /* right multiplication                               */
    /*----------------------------------------------------*/
    void mul_right(Vec<zz_p>& res, const Vec<zz_p>& input) const;
    Vec<zz_p> mul_right(const Vec<zz_p>& input) const
    {
        Vec<zz_p> res;
        mul_right(res, input);
        return res;
    }

    void mul_right(Mat<zz_p>& res, const Mat<zz_p>& input) const;
    Mat<zz_p> mul_right(const Mat<zz_p>& input) const
    {
        Mat<zz_p> res;
        mul_right(res, input);
        return res;
    }

    /*----------------------------------------------------*/
    /* left multiplication                                */
    /*----------------------------------------------------*/
    void mul_left(Vec<zz_p>& res, const Vec<zz_p>& input) const;
    Vec<zz_p> mul_left(const Vec<zz_p>& input) const
    {
        Vec<zz_p> res;
        mul_left(res, input);
        return res;
    }

    void mul_left(Mat<zz_p>& res, const Mat<zz_p>& input) const;
    Mat<zz_p> mul_left(const Mat<zz_p>& input) const
    {
        Mat<zz_p> res;
        mul_left(res, input);
        return res;
    }
};




/*------------------------------------------------------------------*/
/* preconditions M                                                  */
/* builds the matrix CL = (D_e X_int) M (D_f Y_int)^t, where:       */
/* - X_int, Y_int are geometric interpolation                       */
/* - D_e, D_f are diagonal matrix built on vectors e and f          */
/* - CL is cauchy-like special                                      */
/* - CL is expected to have generic rank profile                    */
/*------------------------------------------------------------------*/
void to_cauchy_grp(cauchy_like_geometric_lzz_p& CL, 
                   zz_pX_Multipoint_Geometric& X_int, zz_pX_Multipoint_Geometric& Y_int,
                   Vec<zz_p> &e, Vec<zz_p> &f,
                   const mosaic_hankel_lzz_p& M);

#endif

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
