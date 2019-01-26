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

    void mul_right(Mat<zz_pX>& out, const Mat<zz_pX>& in) const;
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

    void mul_left(Mat<zz_pX>& out, const Mat<zz_pX>& in) const;
    inline Mat<zz_pX> mul_left(const Mat<zz_pX>& in) const
    {
        Mat<zz_pX> out;
        mul_left(out, in);
        return out;
    }

    /*------------------------------------------------------------*/
    /* right multiplication mod x^s                               */
    /*------------------------------------------------------------*/
    virtual void mul_right_trunc(Vec<zz_pX>& out, const Vec<zz_pX>& in, long s) const = 0;
    inline Vec<zz_pX> mul_right_trunc(const Vec<zz_pX>& in, long s) const
    {
        Vec<zz_pX> out;
        mul_right_trunc(out, in, s);
        return out;
    }

    void mul_right_trunc(Mat<zz_pX>& out, const Mat<zz_pX>& in, long s) const;
    inline Mat<zz_pX> mul_right_trunc(const Mat<zz_pX>& in, long s) const
    {
        Mat<zz_pX> out;
        mul_right_trunc(out, in, s);
        return out;
    }

    /*------------------------------------------------------------*/
    /* left multiplication mod x^s                                */
    /*------------------------------------------------------------*/
    virtual void mul_left_trunc(Vec<zz_pX>& out, const Vec<zz_pX>& in, long s) const = 0;
    inline Vec<zz_pX> mul_left_trunc(const Vec<zz_pX>& in, long s) const
    {
        Vec<zz_pX> out;
        mul_left_trunc(out, in, s);
        return out;
    }

    void mul_left_trunc(Mat<zz_pX>& out, const Mat<zz_pX>& in, long s) const;
    inline Mat<zz_pX> mul_left_trunc(const Mat<zz_pX>& in, long s) const
    {
        Mat<zz_pX> out;
        mul_left_trunc(out, in, s);
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

    /*------------------------------------------------------------*/
    /* left multiplication                                        */
    /*------------------------------------------------------------*/
    using structured_lzz_pX::mul_left;
    void mul_left(Vec<zz_pX>& out, const Vec<zz_pX>& in) const;

    /*------------------------------------------------------------*/
    /* right multiplication mod x^s                               */
    /*------------------------------------------------------------*/
    using structured_lzz_pX::mul_right_trunc;
    void mul_right_trunc(Vec<zz_pX>& out, const Vec<zz_pX>& in, long s) const;

    /*------------------------------------------------------------*/
    /* left multiplication mod x^s                                */
    /*------------------------------------------------------------*/
    using structured_lzz_pX::mul_left_trunc;
    void mul_left_trunc(Vec<zz_pX>& out, const Vec<zz_pX>& in, long s) const;

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

    Vec<long> rhs_degX, lhs_degX;
    Vec<fftRep> fft, fft2, fft3;

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

    /*------------------------------------------------------------*/
    /* left multiplication                                        */
    /*------------------------------------------------------------*/
    using structured_lzz_pX::mul_left;
    void mul_left(Vec<zz_pX>& out, const Vec<zz_pX>& in) const;

    /*------------------------------------------------------------*/
    /* right multiplication mod x^s                               */
    /*------------------------------------------------------------*/
    using structured_lzz_pX::mul_right_trunc;
    void mul_right_trunc(Vec<zz_pX>& out, const Vec<zz_pX>& in, long s) const;

    /*------------------------------------------------------------*/
    /* left multiplication mod x^s                                */
    /*------------------------------------------------------------*/
    using structured_lzz_pX::mul_left_trunc;
    void mul_left_trunc(Vec<zz_pX>& out, const Vec<zz_pX>& in, long s) const;

    /*------------------------------------------------------------*/
    /* turns M into a dense matrix                                */
    /*------------------------------------------------------------*/
    using structured_lzz_pX::to_dense;
    void to_dense(Mat<zz_pX>& Mdense) const;
};


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* Toeplitz like matrices, where generators G, H are such that*/
/* Z0 M - M Z1 = G H^t                                        */
/* -> M = sum_i L(-g_i) circ (rev h_i)                        */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
class toeplitz_like_minus_lzz_pX : public structured_lzz_pX 
{
private:
    // m rows, n columns
    long m, n;
    Vec<lower_triangular_toeplitz_lzz_pX> toeplitz_G;
    Vec<circulant_row_lzz_pX> circulant_H;

public:
    Mat<zz_pX> G, H;
    /*------------------------------------------------------------*/
    /* sets dimensions to 0                                       */
    /*------------------------------------------------------------*/
    toeplitz_like_minus_lzz_pX();

    /*------------------------------------------------------------*/
    /* input vector is as showed above, with G<-U, H<-V           */
    /*------------------------------------------------------------*/
    toeplitz_like_minus_lzz_pX(const Mat<zz_pX>& U, const Mat<zz_pX>& V);

    /*------------------------------------------------------------*/
    /* getters                                                    */
    /*------------------------------------------------------------*/
    long NumCols() const;
    long NumRows() const;
    long NumGens() const;

    /*------------------------------------------------------------*/
    /* right multiplication                                       */
    /*------------------------------------------------------------*/
    using structured_lzz_pX::mul_right;
    void mul_right(Vec<zz_pX>& out, const Vec<zz_pX>& in) const;

    /*------------------------------------------------------------*/
    /* left multiplication                                        */
    /*------------------------------------------------------------*/
    using structured_lzz_pX::mul_left;
    void mul_left(Vec<zz_pX>& out, const Vec<zz_pX>& in) const;

    /*------------------------------------------------------------*/
    /* right multiplication mod x^s                               */
    /*------------------------------------------------------------*/
    using structured_lzz_pX::mul_right_trunc;
    void mul_right_trunc(Vec<zz_pX>& out, const Vec<zz_pX>& in, long s) const;

    /*------------------------------------------------------------*/
    /* left multiplication mod x^s                                */
    /*------------------------------------------------------------*/
    using structured_lzz_pX::mul_left_trunc;
    void mul_left_trunc(Vec<zz_pX>& out, const Vec<zz_pX>& in, long s) const;

    /*------------------------------------------------------------*/
    /* turns M into a dense matrix                                */
    /*------------------------------------------------------------*/
    using structured_lzz_pX::to_dense;
    void to_dense(Mat<zz_pX>& Mdense) const;
};

/*------------------------------------------------------------*/
/* returns Z0 A - A Z1                                        */
/*------------------------------------------------------------*/
void toeplitz_lzz_pX_phi_minus(Mat<zz_pX> & res, const Mat<zz_pX>& A);

inline Mat<zz_pX> toeplitz_lzz_pX_phi_minus(const Mat<zz_pX>& A)
{
    Mat<zz_pX> res;
    toeplitz_lzz_pX_phi_minus(res, A);
    return res;
}

#endif


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*                   SYLVESTER MATRICES                       */
/* sylv(A,B)= (square) matrix of (F,G) \mapsto AF + BG        */
/* with deg(F) < deg(B), deg(G) < deg(A)                      */
/* using monomial bases in increasing degrees for rows / cols */
/* particular case of mosaic toeplitz                         */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
class sylvester_lzz_pX : public structured_lzz_pX
{
private:
    long n; // number of rows / columns
    long dAy, dBy; // Y-degrees of a, b (Y is the main variable)
    long dAx, dBx; // X-degrees of a, b (X is the coefficient variable)
    Vec<zz_pX> a, b, revA, revB;

public:
    /*------------------------------------------------------------*/
    /* sets dimensions to 0                                       */
    /*------------------------------------------------------------*/
    sylvester_lzz_pX();

    /*------------------------------------------------------------*/
    /* dimension = deg(A) + deg(B), A and B as above              */
    /* assumes that the last entries of A and B are non-zero      */
    /*------------------------------------------------------------*/
    sylvester_lzz_pX(const Vec<zz_pX>& A, const Vec<zz_pX>& B);

    /*------------------------------------------------------------*/
    /* getters                                                    */
    /*------------------------------------------------------------*/
    long NumCols() const;
    long NumRows() const;

    /*------------------------------------------------------------*/
    /* right multiplication                                       */
    /*------------------------------------------------------------*/
    using structured_lzz_pX::mul_right;
    void mul_right(Vec<zz_pX>& out, const Vec<zz_pX>& in) const;

    /*------------------------------------------------------------*/
    /* left multiplication                                        */
    /*------------------------------------------------------------*/
    using structured_lzz_pX::mul_left;
    void mul_left(Vec<zz_pX>& out, const Vec<zz_pX>& in) const;

    /*------------------------------------------------------------*/
    /* right multiplication mod x^s                               */
    /*------------------------------------------------------------*/
    using structured_lzz_pX::mul_right_trunc;
    void mul_right_trunc(Vec<zz_pX>& out, const Vec<zz_pX>& in, long s) const;

    /*------------------------------------------------------------*/
    /* left multiplication mod x^s                                */
    /*------------------------------------------------------------*/
    using structured_lzz_pX::mul_left_trunc;
    void mul_left_trunc(Vec<zz_pX>& out, const Vec<zz_pX>& in, long s) const;

    /*------------------------------------------------------------*/
    /* G, H such that Z1 M - M Z0 = G H^t                         */
    /*------------------------------------------------------------*/
    void phi_plus_generators(Mat<zz_pX>& G, Mat<zz_pX>& H) const; 

    /*------------------------------------------------------------*/
    /* turns M into a dense matrix                                */
    /*------------------------------------------------------------*/
    using structured_lzz_pX::to_dense;
    void to_dense(Mat<zz_pX>& Mdense) const;

    /*------------------------------------------------------------*/
    /* Newton iteration for inverse                               */
    /* assumes M(0) invertible; error otherwise                   */
    /* return M^-1 mod X^m as a toeplitz_like_minus matrix        */
    /*------------------------------------------------------------*/
    void newton_inv_trunc(toeplitz_like_minus_lzz_pX& iM, long m) const;

    /*------------------------------------------------------------*/
    /*------------------------------------------------------------*/
    void high_precision_inv_trunc(toeplitz_like_minus_lzz_pX& iM, long m) const;

};

/*------------------------------------------------------------*/
/* returns Z1 A - A Z0                                        */
/*------------------------------------------------------------*/
void toeplitz_lzz_pX_phi_plus(Mat<zz_pX> & res, const Mat<zz_pX>& A);

inline Mat<zz_pX> toeplitz_lzz_pX_phi_plus(const Mat<zz_pX>& A)
{
    Mat<zz_pX> res;
    toeplitz_lzz_pX_phi_plus(res, A);
    return res;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
