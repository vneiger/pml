#ifndef __HANKEL_TOEPLITZ_LZZ_P_H
#define __HANKEL_TOEPLITZ_LZZ_P_H

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
class hankel_lzz_p
{
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
    void to_dense(Mat<zz_p>& Mdense) const;
    inline Mat<zz_p> to_dense() const
    {
        Mat<zz_p> dense;
        to_dense(dense);
        return dense;
    }

    /*----------------------------------------------------*/
    /* right multiplication                               */
    /*----------------------------------------------------*/
    void mul_right(Vec<zz_p>& res, const Vec<zz_p>& input)const;
    inline Vec<zz_p> mul_right(const Vec<zz_p>& input) const
    {
        Vec<zz_p> output;
        mul_right(output, input);
        return output;
    }

    void mul_right(Mat<zz_p>& res, const Mat<zz_p>& input) const;
    inline Mat<zz_p> mul_right(const Mat<zz_p>& input) const
    {
        Mat<zz_p> output;
        mul_right(output, input);
        return output;
    }

    /*----------------------------------------------------*/
    /* left multiplication                                */
    /*----------------------------------------------------*/
    void mul_left(Vec<zz_p>& res, const Vec<zz_p>& input) const;
    inline Vec<zz_p> mul_left(const Vec<zz_p>& input) const
    {
        Vec<zz_p> output;
        mul_left(output, input);
        return output;
    }

    void mul_left(Mat<zz_p>& res, const Mat<zz_p>& input) const;
    inline Mat<zz_p> mul_left(const Mat<zz_p>& input) const
    {
        Mat<zz_p> output;
        mul_left(output, input);
        return output;
    }

};


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*     Toeplitz matrices                                      */
/*     stored as                                              */
/*          a3 a2 a1 a0                                       */
/*          a4 a3 a2 a1                                       */
/*          a5 a4 a3 a2                                       */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
class toeplitz_lzz_p
{
private:
    // n rows, m columns
    long n, m;
    Vec<zz_p> data, data_rev;
    fftRep fft_data;

public:
    /*----------------------------------------------------*/
    /* sets dimensions to 0                               */
    /*----------------------------------------------------*/
    toeplitz_lzz_p();

    /*----------------------------------------------------*/
    /* input vector is as showed above                    */
    /*----------------------------------------------------*/
    toeplitz_lzz_p(const Vec<zz_p>& input, long rows, long cols);

    /*----------------------------------------------------*/
    /* getters                                            */
    /*----------------------------------------------------*/
    long NumCols() const;
    long NumRows() const;

    inline zz_p operator ()(long i, long j) const
    {
        return data[m+i-1-j];
    }

    /*----------------------------------------------------*/
    /* turns M into a dense matrix                        */
    /*----------------------------------------------------*/
    void to_dense(Mat<zz_p>& Mdense) const;
    inline Mat<zz_p> to_dense() const
    {
        Mat<zz_p> dense;
        to_dense(dense);
        return dense;
    }

    /*----------------------------------------------------*/
    /* right multiplication                               */
    /*----------------------------------------------------*/
    void mul_right(Vec<zz_p>& res, const Vec<zz_p>& input) const;
    inline Vec<zz_p> mul_right(const Vec<zz_p>& input) const
    {
        Vec<zz_p> output;
        mul_right(output, input);
        return output;
    }

    void mul_right(Mat<zz_p>& res, const Mat<zz_p>& input) const;
    inline Mat<zz_p> mul_right(const Mat<zz_p>& input) const
    {
        Mat<zz_p> output;
        mul_right(output, input);
        return output;
    }

    /*----------------------------------------------------*/
    /* left multiplication                                */
    /*----------------------------------------------------*/
    void mul_left(Vec<zz_p>& res, const Vec<zz_p>& input) const;
    inline Vec<zz_p> mul_left(const Vec<zz_p>& input) const
    {
        Vec<zz_p> output;
        mul_left(output, input);
        return output;
    }

    void mul_left(Mat<zz_p>& res, const Mat<zz_p>& input) const;
    inline Mat<zz_p> mul_left(const Mat<zz_p>& input) const
    {
        Mat<zz_p> output;
        mul_left(output, input);
        return output;
    }

};

#endif
