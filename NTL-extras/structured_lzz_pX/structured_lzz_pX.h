#ifndef __STRUCTURED_LZZ_PX_H
#define __STRUCTURED_LZZ_PX_H

#include <NTL/lzz_pX.h>
#include <NTL/mat_lzz_p.h>

#include "lzz_pX_CRT.h"
#include "mat_lzz_pX_extra.h"
#include "lzz_pXY.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* structured matrices over zz_pX                             */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
class structured_lzz_pX
{
public:
    /*------------------------------------------------------------*/
    /* to avoid memory leaks                                      */
    /*------------------------------------------------------------*/
    virtual ~structured_lzz_pX(){}

    /*------------------------------------------------------------*/
    /* getters                                                    */
    /*------------------------------------------------------------*/
    virtual long NumCols() const = 0;
    virtual long NumRows() const = 0;

    /*------------------------------------------------------------*/
    /* right multiplication                                       */
    /*------------------------------------------------------------*/
    virtual void mul_right(Vec<zz_pX>& out, const Vec<zz_pX>& in) const = 0;
    inline Vec<zz_pX> mul_right(const Vec<zz_pX>& in) const
    {
        Vec<zz_pX> out;
        mul_right(out, in);
        return out;
    }

    virtual void mul_right(Mat<zz_pX>& out, const Mat<zz_pX>& in) const = 0;
    inline Mat<zz_pX> mul_right(const Mat<zz_pX>& in) const
    {
        Mat<zz_pX> out;
        mul_right(out, in);
        return out;
    }

    /*------------------------------------------------------------*/
    /* left multiplication                                        */
    /*------------------------------------------------------------*/
    virtual void mul_left(Vec<zz_pX>& out, const Vec<zz_pX>& in) const = 0;
    inline Vec<zz_pX> mul_left(const Vec<zz_pX>& in) const
    {
        Vec<zz_pX> out;
        mul_left(out, in);
        return out;
    }

    virtual void mul_left(Mat<zz_pX>& out, const Mat<zz_pX>& in) const = 0;
    inline Mat<zz_pX> mul_left(const Mat<zz_pX>& in) const
    {
        Mat<zz_pX> out;
        mul_left(out, in);
        return out;
    }

    /*------------------------------------------------------------*/
    /* turns M into a dense matrix                                */
    /*------------------------------------------------------------*/
    virtual void to_dense(Mat<zz_pX>& Mdense) const = 0;
    inline Mat<zz_pX> to_dense() const
    {
        Mat<zz_pX> dense;
        to_dense(dense);
        return dense;
    }
};


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* matrices such as                                           */
/* [a0  0  0]                                                 */
/* [a1 a0  0]                                                 */
/* [a2 a1 a0]                                                 */
/* (special case of toeplitz matrix)                          */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
class lower_triangular_toeplitz_lzz_pX : public structured_lzz_pX 
{
private:
    // n rows, n columns
    long n;
    long dX, dY; // deg(column) in X, resp Y
    Vec<zz_pX> data;

    Vec<long> rhs_degX;
    Vec<fftRep> fft;

public:
    /*------------------------------------------------------------*/
    /* sets dimensions to 0                                       */
    /*------------------------------------------------------------*/
    lower_triangular_toeplitz_lzz_pX();

    /*------------------------------------------------------------*/
    /* input vector is as showed above                            */
    /*------------------------------------------------------------*/
    lower_triangular_toeplitz_lzz_pX(const Vec<zz_pX>& input);

    /*------------------------------------------------------------*/
    /* getters                                                    */
    /*------------------------------------------------------------*/
    long NumCols() const;
    long NumRows() const;

    /*------------------------------------------------------------*/
    /* prepare for multiplication with rhs X-degree d             */
    /*------------------------------------------------------------*/
    void prepare_degree(long d);

    /*------------------------------------------------------------*/
    /* right multiplication                                       */
    /*------------------------------------------------------------*/
    using structured_lzz_pX::mul_right;
    void mul_right(Vec<zz_pX>& out, const Vec<zz_pX>& in) const;
    void mul_right(Mat<zz_pX>& out, const Mat<zz_pX>& in) const;

    /*------------------------------------------------------------*/
    /* left multiplication                                        */
    /*------------------------------------------------------------*/
    using structured_lzz_pX::mul_left;
    void mul_left(Vec<zz_pX>& out, const Vec<zz_pX>& in) const;
    void mul_left(Mat<zz_pX>& out, const Mat<zz_pX>& in) const;

    /*------------------------------------------------------------*/
    /* turns M into a dense matrix                                */
    /*------------------------------------------------------------*/
    using structured_lzz_pX::to_dense;
    void to_dense(Mat<zz_pX>& Mdense) const;
};



/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* matrices such as                                           */
/* [a3 a2 a1 a0]                                              */
/* [a0 a3 a2 a1]                                              */
/* [a1 a0 a3 a2]                                              */
/* (special case of toeplitz matrix)                          */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
class circulant_row_lzz_pX : public structured_lzz_pX 
{
private:
    // n rows, m columns
    long n, m;
    Vec<zz_pX> data, dataRev;
    long dX, dY; // X- and Y-degrees

    Vec<long> rhs_degX;
    Vec<fftRep> fft;

public:
    /*------------------------------------------------------------*/
    /* sets dimensions to 0                                       */
    /*------------------------------------------------------------*/
    circulant_row_lzz_pX();

    /*------------------------------------------------------------*/
    /* input vector is as showed above                            */
    /*------------------------------------------------------------*/
    circulant_row_lzz_pX(const Vec<zz_pX>& input, long nrows);

    /*------------------------------------------------------------*/
    /* getters                                                    */
    /*------------------------------------------------------------*/
    long NumCols() const;
    long NumRows() const;

    /*------------------------------------------------------------*/
    /* prepare for multiplication with rhs X-degree d             */
    /*------------------------------------------------------------*/
    void prepare_degree(long d);

    /*------------------------------------------------------------*/
    /* right multiplication                                       */
    /*------------------------------------------------------------*/
    using structured_lzz_pX::mul_right;
    void mul_right(Vec<zz_pX>& out, const Vec<zz_pX>& in) const;
    void mul_right(Mat<zz_pX>& out, const Mat<zz_pX>& in) const;

    /*------------------------------------------------------------*/
    /* left multiplication                                        */
    /*------------------------------------------------------------*/
    using structured_lzz_pX::mul_left;
    void mul_left(Vec<zz_pX>& out, const Vec<zz_pX>& in) const;
    void mul_left(Mat<zz_pX>& out, const Mat<zz_pX>& in) const;

    /*------------------------------------------------------------*/
    /* turns M into a dense matrix                                */
    /*------------------------------------------------------------*/
    using structured_lzz_pX::to_dense;
    void to_dense(Mat<zz_pX>& Mdense) const;
};


// /*------------------------------------------------------------*/
// /*------------------------------------------------------------*/
// /* Toeplitz like matrices, where generators G, H are such that*/
// /* Z0 M - M Z1 = G H^t                                        */
// /* -> M = sum_i L(-g_i) circ (rev h_i)                        */
// /*------------------------------------------------------------*/
// /*------------------------------------------------------------*/
// class toeplitz_like_minus_lzz_pX : public structured_lzz_pX 
// {
// private:
//     // m rows, n columns
//     long m, n;
//     Vec<lower_triangular_toeplitz_lzz_pX> toeplitz_G;
//     Vec<circulant_row_lzz_pX> circulant_H;

// public:
//     Mat<zz_pX> G, H;
//     /*------------------------------------------------------------*/
//     /* sets dimensions to 0                                       */
//     /*------------------------------------------------------------*/
//     toeplitz_like_minus_lzz_pX();

//     /*------------------------------------------------------------*/
//     /* input vector is as showed above, with G<-U, H<-V           */
//     /*------------------------------------------------------------*/
//     toeplitz_like_minus_lzz_pX(const Mat<zz_pX>& U, const Mat<zz_pX>& V);

//     /*------------------------------------------------------------*/
//     /* getters                                                    */
//     /*------------------------------------------------------------*/
//     long NumCols() const;
//     long NumRows() const;
//     long NumGens() const;

//     /*------------------------------------------------------------*/
//     /* right multiplication                                       */
//     /*------------------------------------------------------------*/
//     using structured_lzz_pX::mul_right;
//     void mul_right(Vec<zz_pX>& out, const Vec<zz_pX>& in) const;
//     void mul_right(Mat<zz_pX>& out, const Mat<zz_pX>& in) const;

//     /*------------------------------------------------------------*/
//     /* left multiplication                                        */
//     /*------------------------------------------------------------*/
//     using structured_lzz_pX::mul_left;
//     void mul_left(Vec<zz_pX>& out, const Vec<zz_pX>& in) const;
//     void mul_left(Mat<zz_pX>& out, const Mat<zz_pX>& in) const;

//     /*------------------------------------------------------------*/
//     /* turns M into a dense matrix                                */
//     /*------------------------------------------------------------*/
//     using structured_lzz_pX::to_dense;
//     void to_dense(Mat<zz_pX>& Mdense) const;
// };

// /*------------------------------------------------------------*/
// /* returns Z0 A - A Z1                                        */
// /*------------------------------------------------------------*/
// void toeplitz_lzz_pX_phi_minus(Mat<zz_pX> & res, const Mat<zz_pX>& A);

// inline Mat<zz_pX> toeplitz_lzz_pX_phi_minus(const Mat<zz_pX>& A)
// {
//     Mat<zz_pX> res;
//     toeplitz_lzz_pX_phi_minus(res, A);
//     return res;
// }

#endif


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
