#ifndef __STRUCTURED_LZZ_P_H
#define __STRUCTURED_LZZ_P_H

#include <NTL/lzz_pX.h>
#include <NTL/mat_lzz_p.h>

#include "lzz_pX_CRT.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
class structured_lzz_p
{
public:
    /*------------------------------------------------------------*/
    /* to avoid memory leaks                                      */
    /*------------------------------------------------------------*/
    virtual ~structured_lzz_p(){}

    /*------------------------------------------------------------*/
    /* getters                                                    */
    /*------------------------------------------------------------*/
    virtual long NumCols() const = 0;
    virtual long NumRows() const = 0;

    /*------------------------------------------------------------*/
    /* right multiplication                                       */
    /*------------------------------------------------------------*/
    virtual void mul_right(Vec<zz_p>& out, const Vec<zz_p>& in) const = 0;
    inline Vec<zz_p> mul_right(const Vec<zz_p>& in) const
    {
        Vec<zz_p> out;
        mul_right(out, in);
        return out;
    }

    virtual void mul_right(Mat<zz_p>& out, const Mat<zz_p>& in) const = 0;
    inline Mat<zz_p> mul_right(const Mat<zz_p>& in) const
    {
        Mat<zz_p> out;
        mul_right(out, in);
        return out;
    }

    /*------------------------------------------------------------*/
    /* left multiplication                                        */
    /*------------------------------------------------------------*/
    virtual void mul_left(Vec<zz_p>& out, const Vec<zz_p>& in) const = 0;
    inline Vec<zz_p> mul_left(const Vec<zz_p>& in) const
    {
        Vec<zz_p> out;
        mul_left(out, in);
        return out;
    }

    virtual void mul_left(Mat<zz_p>& out, const Mat<zz_p>& in) const = 0;
    inline Mat<zz_p> mul_left(const Mat<zz_p>& in) const
    {
        Mat<zz_p> out;
        mul_left(out, in);
        return out;
    }

    /*------------------------------------------------------------*/
    /* turns M into a dense matrix                                */
    /*------------------------------------------------------------*/
    virtual void to_dense(Mat<zz_p>& Mdense) const = 0;
    inline Mat<zz_p> to_dense() const
    {
        Mat<zz_p> dense;
        to_dense(dense);
        return dense;
    }
};


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*     Hankel matrices                                        */
/*     stored as                                              */
/*          a5 a4 a3 a2                                       */
/*          a4 a3 a2 a1                                       */
/*          a3 a2 a1 a0                                       */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
class hankel_lzz_p : public structured_lzz_p 
{
private:
    // n rows, m columns
    long n, m;
    Vec<zz_p> data, data_rev;
    fftRep fft_data;

public:
    /*------------------------------------------------------------*/
    /* sets dimensions to 0                                       */
    /*------------------------------------------------------------*/
    hankel_lzz_p();

    /*------------------------------------------------------------*/
    /* input vector is as showed above                            */
    /*------------------------------------------------------------*/
    hankel_lzz_p(const Vec<zz_p>& input, long rows, long cols);

    /*------------------------------------------------------------*/
    /* getters                                                    */
    /*------------------------------------------------------------*/
    long NumCols() const;
    long NumRows() const;

    inline zz_p operator ()(long i, long j) const
    {
        return data[m+n-2-i-j];
    }

    /*------------------------------------------------------------*/
    /* right multiplication                                       */
    /*------------------------------------------------------------*/
    using structured_lzz_p::mul_right;
    void mul_right(Vec<zz_p>& out, const Vec<zz_p>& in) const;
    void mul_right(Mat<zz_p>& out, const Mat<zz_p>& in) const;

    /*------------------------------------------------------------*/
    /* left multiplication                                        */
    /*------------------------------------------------------------*/
    using structured_lzz_p::mul_left;
    void mul_left(Vec<zz_p>& out, const Vec<zz_p>& in) const;
    void mul_left(Mat<zz_p>& out, const Mat<zz_p>& in) const;

    /*------------------------------------------------------------*/
    /* turns M into a dense matrix                                */
    /*------------------------------------------------------------*/
    using structured_lzz_p::to_dense;
    void to_dense(Mat<zz_p>& Mdense) const;
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
class toeplitz_lzz_p : public structured_lzz_p 
{
private:
    // n rows, m columns
    long n, m;
    Vec<zz_p> data, data_rev;
    fftRep fft_data;

public:
    /*------------------------------------------------------------*/
    /* sets dimensions to 0                                       */
    /*------------------------------------------------------------*/
    toeplitz_lzz_p();

    /*------------------------------------------------------------*/
    /* input vector is as showed above                            */
    /*------------------------------------------------------------*/
    toeplitz_lzz_p(const Vec<zz_p>& input, long rows, long cols);

    /*------------------------------------------------------------*/
    /* getters                                                    */
    /*------------------------------------------------------------*/
    long NumCols() const;
    long NumRows() const;

    inline zz_p operator ()(long i, long j) const
    {
        return data[m+i-1-j];
    }

    /*------------------------------------------------------------*/
    /* right multiplication                                       */
    /*------------------------------------------------------------*/
    using structured_lzz_p::mul_right;
    void mul_right(Vec<zz_p>& out, const Vec<zz_p>& in) const;
    void mul_right(Mat<zz_p>& out, const Mat<zz_p>& in) const;

    /*------------------------------------------------------------*/
    /* left multiplication                                        */
    /*------------------------------------------------------------*/
    using structured_lzz_p::mul_left;
    void mul_left(Vec<zz_p>& out, const Vec<zz_p>& in) const;
    void mul_left(Mat<zz_p>& out, const Mat<zz_p>& in) const;

    /*------------------------------------------------------------*/
    /* turns M into a dense matrix                                */
    /*------------------------------------------------------------*/
    using structured_lzz_p::to_dense;
    void to_dense(Mat<zz_p>& Mdense) const;
};


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
class toeplitz_like_plus_lzz_p : public structured_lzz_p 
{

};

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
class toeplitz_like_minus_lzz_p : public structured_lzz_p 
{

};

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* Hankel like matrices, where generators G, H are such that  */
/* Z1 M - M Z0^t = G H^t                                      */
/* -> M = sum_i circ(g_i) L(h_i) J                            */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
class hankel_like_plus_lzz_p : public structured_lzz_p 
{
private:
    // m rows, n columns
    long m, n;
    Vec<toeplitz_lzz_p> toeplitz_G;
    Vec<hankel_lzz_p> hankel_H;

public:
    Mat<zz_p> G, H;
    /*------------------------------------------------------------*/
    /* sets dimensions to 0                                       */
    /*------------------------------------------------------------*/
    hankel_like_plus_lzz_p();

    /*------------------------------------------------------------*/
    /* input vector is as showed above, with G<-U, H<-V           */
    /*------------------------------------------------------------*/
    hankel_like_plus_lzz_p(const Mat<zz_p>& U, const Mat<zz_p>& V);

    /*------------------------------------------------------------*/
    /* getters                                                    */
    /*------------------------------------------------------------*/
    long NumCols() const;
    long NumRows() const;
    long NumGens() const;

    /*------------------------------------------------------------*/
    /* right multiplication                                       */
    /*------------------------------------------------------------*/
    using structured_lzz_p::mul_right;
    void mul_right(Vec<zz_p>& out, const Vec<zz_p>& in) const;
    void mul_right(Mat<zz_p>& out, const Mat<zz_p>& in) const;

    /*------------------------------------------------------------*/
    /* left multiplication                                        */
    /*------------------------------------------------------------*/
    using structured_lzz_p::mul_left;
    void mul_left(Vec<zz_p>& out, const Vec<zz_p>& in) const;
    void mul_left(Mat<zz_p>& out, const Mat<zz_p>& in) const;

    /*------------------------------------------------------------*/
    /* turns M into a dense matrix                                */
    /*------------------------------------------------------------*/
    using structured_lzz_p::to_dense;
    void to_dense(Mat<zz_p>& Mdense) const;
};

/*------------------------------------------------------------*/
/* returns Z1 A - A Z0^t                                      */
/*------------------------------------------------------------*/
void hankel_lzz_p_phi_plus(Mat<zz_p> & res, const Mat<zz_p>& A);

inline Mat<zz_p> hankel_lzz_p_phi_plus(const Mat<zz_p>& A)
{
    Mat<zz_p> res;
    hankel_lzz_p_phi_plus(res, A);
    return res;
}


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* Hankel like matrices, where generators G, H are such that  */
/* Z0^t M - M Z1 = G H^t                                      */
/* -> M = - sum_i U(g_i) circ (rev h_i)                       */
/* (dual structure to hankel_like_plus_lzz_p)                 */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
class hankel_like_minus_lzz_p : public structured_lzz_p 
{
private:
    // m rows, n columns
    long m, n;
    Vec<hankel_lzz_p> hankel_G;
    Vec<toeplitz_lzz_p> toeplitz_H;

public:
    Mat<zz_p> G, H;
    /*------------------------------------------------------------*/
    /* sets dimensions to 0                                       */
    /*------------------------------------------------------------*/
    hankel_like_minus_lzz_p();

    /*------------------------------------------------------------*/
    /* input vector is as showed above, with G<-U, H<-V           */
    /*------------------------------------------------------------*/
    hankel_like_minus_lzz_p(const Mat<zz_p>& U, const Mat<zz_p>& V);

    /*------------------------------------------------------------*/
    /* getters                                                    */
    /*------------------------------------------------------------*/
    long NumCols() const;
    long NumRows() const;
    long NumGens() const;

    /*------------------------------------------------------------*/
    /* right multiplication                                       */
    /*------------------------------------------------------------*/
    using structured_lzz_p::mul_right;
    void mul_right(Vec<zz_p>& out, const Vec<zz_p>& in) const;
    void mul_right(Mat<zz_p>& out, const Mat<zz_p>& in) const;

    /*------------------------------------------------------------*/
    /* left multiplication                                        */
    /*------------------------------------------------------------*/
    using structured_lzz_p::mul_left;
    void mul_left(Vec<zz_p>& out, const Vec<zz_p>& in) const;
    void mul_left(Mat<zz_p>& out, const Mat<zz_p>& in) const;

    /*------------------------------------------------------------*/
    /* turns M into a dense matrix                                */
    /*------------------------------------------------------------*/
    using structured_lzz_p::to_dense;
    void to_dense(Mat<zz_p>& Mdense) const;
};

/*------------------------------------------------------------*/
/* returns Z0^t A - A Z1                                      */
/*------------------------------------------------------------*/
void hankel_lzz_p_phi_minus(Mat<zz_p> & res, const Mat<zz_p>& A);

inline Mat<zz_p> hankel_lzz_p_phi_minus(const Mat<zz_p>& A)
{
    Mat<zz_p> res;
    hankel_lzz_p_phi_minus(res, A);
    return res;
}

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* Cauchy matrices on geometric progressions                  */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
class cauchy_geometric_lzz_p : public structured_lzz_p 
{
private:
    long m, n;  // m rows, n columns
public:
    zz_p u1, v1, rho, sqrt_rho; // see constructor for explanations
    toeplitz_lzz_p t;
    Vec<zz_p> powers_irho;

    /*------------------------------------------------------------*/
    /* default constructor                                        */
    /*------------------------------------------------------------*/
    cauchy_geometric_lzz_p();

    /*------------------------------------------------------------*/
    /* constructor                                                */
    /* M[i][j] = 1 / (u1*rho^i - v1*rho^j)                        */
    /* rho should be a square (??????????)                        */
    /*------------------------------------------------------------*/
    cauchy_geometric_lzz_p(const zz_p& u1, const zz_p& v1, const zz_p& rho, long mm, long nn);

    /*------------------------------------------------------------*/
    /* getters                                                    */
    /*------------------------------------------------------------*/
    long NumRows() const;
    long NumCols() const;

    /*------------------------------------------------------------*/
    /* computes output = M*input                                  */
    /*------------------------------------------------------------*/
    using structured_lzz_p::mul_right;
    void mul_right(Vec<zz_p>& output, const Vec<zz_p>& input) const;
    void mul_right(Mat<zz_p>& output, const Mat<zz_p>& input) const;

    /*------------------------------------------------------------*/
    /* computes output = M*input, without the diagonal            */
    /*------------------------------------------------------------*/
    void mul_right_simple(Vec<zz_p>& output, const Vec<zz_p>& input) const;
    inline Vec<zz_p> mul_right_simple(const Vec<zz_p>& input) const
    {
        Vec<zz_p> output;
        mul_right_simple(output, input);
        return output;
    }

    void mul_right_simple(Mat<zz_p>& output, const Mat<zz_p>& input) const;
    inline Mat<zz_p> mul_right_simple(const Mat<zz_p>& input) const
    {
        Mat<zz_p> output;
        mul_right_simple(output, input);
        return output;
    }

    /*------------------------------------------------------------*/
    /* computes output = input * M                                */
    /*------------------------------------------------------------*/
    using structured_lzz_p::mul_left;
    void mul_left(Vec<zz_p>& output, const Vec<zz_p>& input) const;
    void mul_left(Mat<zz_p>& output, const Mat<zz_p>& input) const;

    /*------------------------------------------------------------*/
    /* computes output = input * M, without the diagonal          */
    /*------------------------------------------------------------*/
    void mul_left_simple(Vec<zz_p>& output, const Vec<zz_p>& input) const;
    inline Vec<zz_p> mul_left_simple(const Vec<zz_p>& input) const
    {
        Vec<zz_p> output;
        mul_left_simple(output, input);
        return output;
    }

    void mul_left_simple(Mat<zz_p>& output, const Mat<zz_p>& input) const;
    inline Mat<zz_p> mul_left_simple(const Mat<zz_p>& input) const
    {
        Mat<zz_p> output;
        mul_left_simple(output, input);
        return output;
    }

    /*------------------------------------------------------------*/
    /* M as a dense matrix                                        */
    /*------------------------------------------------------------*/
    using structured_lzz_p::to_dense;
    void to_dense(Mat<zz_p>& M) const;
};

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* Cauchy-like matrices on geometric progressions             */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
class cauchy_like_geometric_lzz_p : public structured_lzz_p 
{
private: 
    long m, n;

public:
    cauchy_geometric_lzz_p C;
    Mat<zz_p> G, H;

    /*------------------------------------------------------------*/
    /* default constructor                                        */
    /*------------------------------------------------------------*/
    cauchy_like_geometric_lzz_p();

    /*------------------------------------------------------------*/
    /* constructor                                                */
    /* let C = 1 / (u1*rho^i - v1*rho^j)                          */
    /* then M = sum_k diag(U[*,k]) C diag(V[*,k])                 */
    /* rho should be a square (?????????)                         */
    /*------------------------------------------------------------*/
    cauchy_like_geometric_lzz_p(const Mat<zz_p>& U, const Mat<zz_p>& V, const zz_p& u1, const zz_p& v1, const zz_p& rho);

    /*------------------------------------------------------------*/
    /* dimensions                                                 */
    /*------------------------------------------------------------*/
    long NumRows() const;
    long NumCols() const;
    long NumGens() const;

    /*------------------------------------------------------------*/
    /* computes output = M*input                                  */
    /*------------------------------------------------------------*/
    using structured_lzz_p::mul_right;
    void mul_right(Vec<zz_p>& output, const Vec<zz_p>& input) const;
    void mul_right_direct(Mat<zz_p> & out, const Mat<zz_p> & in) const;
    void mul_right_sigma_UL(Mat<zz_p> & out, const Mat<zz_p> & in) const;
    void mul_right(Mat<zz_p>& output, const Mat<zz_p>& input) const;

    inline Mat<zz_p> mul_right_direct(const Mat<zz_p>& input) const
    {
        Mat<zz_p> output;
        mul_right_direct(output, input);
        return output;
    }

    /*------------------------------------------------------------*/
    /* computes output = input * M                                */
    /*------------------------------------------------------------*/
    using structured_lzz_p::mul_left;
    void mul_left(Vec<zz_p>& output, const Vec<zz_p>& input) const;
    void mul_left(Mat<zz_p>& output, const Mat<zz_p>& input) const;

    /*------------------------------------------------------------*/
    /* returns the inverse of the leading principal minor         */
    /* (for the dual operator)                                    */
    /* assumes generic rank profile                               */
    /* returns 0 if not generic rank profile; 1 otherwise         */
    /* thresh is threshold for divide-and-conquer                 */
    /* thresh_alpha switches between block and plain quadratic    */
    /*------------------------------------------------------------*/
    long invert_leading_principal_minor_grp(cauchy_like_geometric_lzz_p& Cinv,
                                            long thresh = -1, long thresh_alpha = -1) const;

    /*------------------------------------------------------------*/
    /* finds a random solution to M.x = b                         */
    /* assumes generic rank profile, returns 1 if so, 0 if not    */
    /* if generic rank profile, x is empty if no solution exists  */
    /* thresh is threshold for divide-and-conquer                 */
    /* thresh_alpha switches between block and plain quadratic    */
    /*------------------------------------------------------------*/
    long solve_grp(Vec<zz_p>& x, const Vec<zz_p> b, long thresh = -1, long thresh_alpha = -1) const;

    /*------------------------------------------------------------*/
    /* M as a dense matrix                                        */
    /*------------------------------------------------------------*/
    using structured_lzz_p::to_dense;
    void to_dense(Mat<zz_p>& M) const;
};


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* Mosaic Hankel matrices:                                    */
/* block matrix where each block is Hankel                    */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
class mosaic_hankel_lzz_p : public structured_lzz_p 
{
private:
    long n, m; // numbers of rows / colums
    long nb, mb; // number of blocks in rows / columns

public:
    // data[i][j] is the block at i-th row, j-th column
    Vec< Vec<hankel_lzz_p> > data;

    /*------------------------------------------------------------*/
    /* dummy constructor                                          */
    /*------------------------------------------------------------*/
    mosaic_hankel_lzz_p();

    /*------------------------------------------------------------*/
    /* copies all data                                            */
    /*------------------------------------------------------------*/
    mosaic_hankel_lzz_p(const Vec< Vec<hankel_lzz_p> > & init);

    /*------------------------------------------------------------*/
    /* getters                                                    */
    /*------------------------------------------------------------*/
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

    /*------------------------------------------------------------*/
    /* access to particular rows and columns                      */
    /*------------------------------------------------------------*/
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

    /*------------------------------------------------------------*/
    /* right multiplication                                       */
    /*------------------------------------------------------------*/
    using structured_lzz_p::mul_right;
    void mul_right(Vec<zz_p>& res, const Vec<zz_p>& input) const;
    void mul_right(Mat<zz_p>& res, const Mat<zz_p>& input) const;

    /*------------------------------------------------------------*/
    /* left multiplication                                        */
    /*------------------------------------------------------------*/
    using structured_lzz_p::mul_left;
    void mul_left(Vec<zz_p>& out, const Vec<zz_p>& in) const;
    void mul_left(Mat<zz_p>& out, const Mat<zz_p>& in) const;

    /*------------------------------------------------------------*/
    /* G, H such that Z1 M - M Z0^t = G H^t                       */
    /*------------------------------------------------------------*/
    void phi_plus_generators(Mat<zz_p>& G, Mat<zz_p>& H) const; 

    /*------------------------------------------------------------*/
    /* finds a random solution to M.x = b                         */
    /* if return value is 0, something went wrong (no grp)        */
    /* else, return value is 1. Then x is empty if no solution    */
    /* thresh is threshold for divide-and-conquer                 */
    /* thresh_alpha switches between block and plain quadratic    */
    /*------------------------------------------------------------*/
    long solve(Vec<zz_p>& x, const Vec<zz_p> b, long thresh = -1, long thresh_alpha = -1) const;

    /*------------------------------------------------------------*/
    /* finds the inverse of this as a hankel_like_minus matrix    */
    /* if return value is 0, either no grp or singular            */
    /* else, return value is 1, and inv is the inverse            */
    /* thresh is threshold for divide-and-conquer                 */
    /* thresh_alpha switches between block and plain quadratic    */
    /*------------------------------------------------------------*/
    long inv(hankel_like_minus_lzz_p& inv, long thresh = -1, long thresh_alpha = -1) const;

    /*------------------------------------------------------------*/
    /* turns M into a dense matrix                                */
    /*------------------------------------------------------------*/
    using structured_lzz_p::to_dense;
    void to_dense(Mat<zz_p>& Mdense) const;
};

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* Mosaic Toeplitz matrices:                                  */
/* block matrix where each block is Toeplitz                  */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
class mosaic_toeplitz_lzz_p : public structured_lzz_p 
{
private:
    long n, m; // numbers of rows / colums
    long nb, mb; // number of blocks in rows / columns

public:
    // data[i][j] is the block at i-th row, j-th column
    Vec< Vec<toeplitz_lzz_p> > data;

    /*------------------------------------------------------------*/
    /* dummy constructor                                          */
    /*------------------------------------------------------------*/
    mosaic_toeplitz_lzz_p();

    /*------------------------------------------------------------*/
    /* copies all data                                            */
    /*------------------------------------------------------------*/
    mosaic_toeplitz_lzz_p(const Vec< Vec<toeplitz_lzz_p> > & init);

    /*------------------------------------------------------------*/
    /*----------------------------------------------------*/
    /* getters                                            */
    /*------------------------------------------------------------*/
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


    /*------------------------------------------------------------*/
    /* access to particular rows and columns                      */
    /*------------------------------------------------------------*/
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

    /*------------------------------------------------------------*/
    /* right multiplication                                       */
    /*------------------------------------------------------------*/
    using structured_lzz_p::mul_right;
    void mul_right(Vec<zz_p>& res, const Vec<zz_p>& input) const;
    void mul_right(Mat<zz_p>& res, const Mat<zz_p>& input) const;

    /*------------------------------------------------------------*/
    /* left multiplication                                        */
    /*------------------------------------------------------------*/
    using structured_lzz_p::mul_left;
    void mul_left(Vec<zz_p>& res, const Vec<zz_p>& input) const;
    void mul_left(Mat<zz_p>& res, const Mat<zz_p>& input) const;

    /*------------------------------------------------------------*/
    /* turns M into a dense matrix                                */
    /*------------------------------------------------------------*/
    using structured_lzz_p::to_dense;
    void to_dense(Mat<zz_p>& Mdense) const;
};



#endif 

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
