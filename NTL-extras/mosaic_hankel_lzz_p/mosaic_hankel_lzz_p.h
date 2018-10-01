#ifndef __MOSAIC_HANKEL_LZZ_P_H
#define __MOSAIC_HANKEL_LZZ_P_H

#include <NTL/lzz_pX.h>
#include <NTL/mat_lzz_p.h>

#include "lzz_pX_CRT.h"

NTL_CLIENT


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*     Hankel matrices                                        */
/*     stored as                                              */
/*          a5 a4 a3 a2                                       */
/*          a4 a3 a2 a1                                       */
/*          a3 a2 a1 a0                                       */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
class hankel_lzz_p{
private:
    // n rows, m columns
    long n, m;
    Vec<zz_p> data, data_rev;
    fftRep fft_data;

public:
    /*----------------------------------------------------*/
    /* sets dimensions to 0                               */
    /*----------------------------------------------------*/
    hankel_lzz_p();

    /*----------------------------------------------------*/
    /* input vector is as showed above                    */
    /*----------------------------------------------------*/
    hankel_lzz_p(const Vec<zz_p>& input, long rows, long cols);

    /*----------------------------------------------------*/
    /* getters                                            */
    /*----------------------------------------------------*/
    long NumCols() const;
    long NumRows() const;

    inline zz_p operator ()(long i, long j) const
    {
        return data[m+n-2-i-j];
    }

    /*----------------------------------------------------*/
    /* turns M into a dense matrix                        */
    /*----------------------------------------------------*/
    void to_dense(Mat<zz_p>& Mdense);
    inline Mat<zz_p> to_dense()
    {
        Mat<zz_p> dense;
        to_dense(dense);
        return dense;
    }

    /*----------------------------------------------------*/
    /* right multiplication                               */
    /*----------------------------------------------------*/
    void mul_right(Vec<zz_p>& res, const Vec<zz_p>& input);
    inline Vec<zz_p> mul_right(const Vec<zz_p>& input)
    {
        Vec<zz_p> output;
        mul_right(output, input);
        return output;
    }

    void mul_right(Mat<zz_p>& res, const Mat<zz_p>& input);
    inline Mat<zz_p> mul_right(const Mat<zz_p>& input)
    {
        Mat<zz_p> output;
        mul_right(output, input);
        return output;
    }

    /*----------------------------------------------------*/
    /* left multiplication                                */
    /*----------------------------------------------------*/
    void mul_left(Vec<zz_p>& res, const Vec<zz_p>& input);
    inline Vec<zz_p> mul_left(const Vec<zz_p>& input)
    {
        Vec<zz_p> output;
        mul_left(output, input);
        return output;
    }

    void mul_left(Mat<zz_p>& res, const Mat<zz_p>& input);
    inline Mat<zz_p> mul_left(const Mat<zz_p>& input)
    {
        Mat<zz_p> output;
        mul_left(output, input);
        return output;
    }

};








// /*----------------------------------------------------*/
// /*----------------------------------------------------*/
// /* Mosaic Hankel matrices:                            */
// /* block matrix where each block is Hankel            */
// /*----------------------------------------------------*/
// /*----------------------------------------------------*/
// class mosaic_hankel_lzz_p{
//     long n, m; // numbers of rows / colums
//     long nb, mb; // number of blocks in rows / columns

// public:

//     // data[i][j] is the block at i-th row, j-th column
//     Vec< Vec<hankel_lzz_p> > data;

//     /*----------------------------------------------------*/
//     /* dummy constructor                                  */
//     /*----------------------------------------------------*/
//     mosaic_hankel_lzz_p(){
//         data.SetLength(0);
//         n = m = nb = mb = 0;
//     }

//     /*----------------------------------------------------*/
//     /* copies all data                                    */
//     /*----------------------------------------------------*/
//     mosaic_hankel_lzz_p(Vec< Vec<hankel_lzz_p> > init){
//         data = init;
//         nb = init.length();
//         mb = init[0].length();

//         n = 0;
//         m = 0;

//         for(long i = 0; i < nb; i++)
//             n += init[i][0].NumRows();
//         for (long j = 0; j < mb; j++)
//             m += init[0][j].NumCols();
//     }

//     /*----------------------------------------------------*/
//     /* getters                                            */
//     /*----------------------------------------------------*/
//     long NumRows() const{
//         return n;
//     }
//     long NumCols() const{
//         return m;
//     }
//     long NumBlockRows() const{
//         return nb;
//     }
//     long NumBlockCols() const{
//         return mb;
//     }
//     long NumRows_of_block(long i) const{
//         return data[i][0].NumRows();
//     }
//     long NumCols_of_block(long i) const{
//         return data[0][i].NumCols();
//     }

//     const zz_p& operator ()(long i, long j) const{

//         if (i >= n || j >= m)
//             Error("matrix indices out of bounds\n");
//         long idx, jdx;

//         idx = 0;
//         while (i >= data[idx][0].NumRows()){
//             i -= data[idx][0].NumRows();
//             idx++;
//         }

//         jdx = 0;
//         while (j >= data[0][jdx].NumCols()){
//             j -= data[0][jdx].NumCols();
//             jdx++;
//         }
//         return data[idx][jdx](i, j);
//     }

// };

// /*----------------------------------------------------*/
// /* turns M into a dense matrix                        */
// /*----------------------------------------------------*/
// void to_dense(Mat<zz_p>& Mdense, const mosaic_hankel_lzz_p& M);

// /*----------------------------------------------------*/
// /* right multiplication                               */
// /*----------------------------------------------------*/
// void mul_right(Vec<zz_p>& res, const mosaic_hankel_lzz_p& M, const Vec<zz_p>& input);

// /*----------------------------------------------------*/
// /* left multiplication                                */
// /*----------------------------------------------------*/
// void mul_left(Vec<zz_p>& res, const mosaic_hankel_lzz_p& M, const Vec<zz_p>& input);

// /*----------------------------------------------------*/
// /* access to particular rows and columns              */
// /*----------------------------------------------------*/
// void first_column_of_block(Vec<zz_p>& res, long i, const mosaic_hankel_lzz_p& M);
// void last_column_of_block(Vec<zz_p>& res, long i, const mosaic_hankel_lzz_p& M);
// void first_row_of_block(Vec<zz_p>& res, long i, const mosaic_hankel_lzz_p& M);
// void last_row_of_block(Vec<zz_p>& res, long i, const mosaic_hankel_lzz_p& M);

// /*----------------------------------------------------*/
// /* G, H such that Z1 M - Z0^t M = G H^t               */
// /*----------------------------------------------------*/
// void generators(Mat<zz_p>& G, Mat<zz_p>& H, const mosaic_hankel_lzz_p& M); 


// ----- Need work on cauchy like matrices first  ---//
// /*------------------------------------------------------------------*/
// /* preconditions M                                                  */
// /* builds the matrix CL = (D_e X_int) M (D_f Y_int)^t, where:       */
// /* - X_int, Y_int are geometric interpolation                       */
// /* - D_e, D_f are diagonal matrix built on vectors e and f          */
// /* - CL is cauchy-like special                                      */
// /* - CL is expected to have generic rank profile                    */
// /*------------------------------------------------------------------*/
// void to_cauchy_grp(cauchy_like_geometric_special& CL, 
//                    zz_pX_Multipoint_Geometric& X_int, zz_pX_Multipoint_Geometric& Y_int,
//                    Vec<zz_p> &e, Vec<zz_p> &f,
//                    const mosaic_hankel_lzz_p& M);

#endif

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
