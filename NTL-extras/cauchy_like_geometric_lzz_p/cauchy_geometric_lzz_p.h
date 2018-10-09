#ifndef CAUCHY_GEOMETRIC_LZZ_P_H
#define CAUCHY_GEOMETRIC_LZZ_P_H

#include <NTL/lzz_pX.h>
#include <NTL/mat_lzz_p.h>

#include "lzz_pX_CRT.h"
#include "mosaic_hankel_lzz_p.h"

#define THRESHOLD 100

NTL_CLIENT

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* Cauchy matrices on geometric progressions                  */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
class cauchy_geometric_lzz_p
{
private:
    long m, n;  // m rows, n columns
 
public:
    zz_p u1, v1, rho, sqrt_rho; // see constructor for explanations
    Vec<zz_p> vec_toeplitz;
    toeplitz_lzz_p t;
    zz_pX_Multipoint_Geometric X, Y; // VdM matrices on u1*rho^i, v1*rho^j
    Vec<zz_p> powers_irho;

    /*------------------------------------------------------------*/
    /* default constructor                                        */
    /*------------------------------------------------------------*/
    cauchy_geometric_lzz_p();

    /*------------------------------------------------------------*/
    /* constructor                                                */
    /* M[i][j] = 1 / (u1*rho^i - v1*rho^j)                        */
    /* rho should be a square                                     */
    /*------------------------------------------------------------*/
    cauchy_geometric_lzz_p(const zz_p& a1, const zz_p& b1, const zz_p& rho, long mm, long nn);

    /*------------------------------------------------------------*/
    /* dimensions                                                 */
    /*------------------------------------------------------------*/
    long NumRows() const;
    long NumCols() const;

    /*------------------------------------------------------------*/
    /* computes output = M*input                                  */
    /*------------------------------------------------------------*/
    void mul_right(Vec<zz_p>& output, const Vec<zz_p>& input) const;
    inline Vec<zz_p> mul_right(const Vec<zz_p>& input) const
    {
        Vec<zz_p> output;
        mul_right(output, input);
        return output;
    }

    void mul_right(Mat<zz_p>& output, const Mat<zz_p>& input) const;
    inline Mat<zz_p> mul_right(const Mat<zz_p>& input) const
    {
        Mat<zz_p> output;
        mul_right(output, input);
        return output;
    }

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
    /* computes output = M^t*input                                */
    /*------------------------------------------------------------*/
    void mul_left(Vec<zz_p>& output, const Vec<zz_p>& input) const;
    inline Vec<zz_p> mul_left(const Vec<zz_p>& input) const
    {
        Vec<zz_p> output;
        mul_left(output, input);
        return output;
    }

    void mul_left(Mat<zz_p>& output, const Mat<zz_p>& input) const;
    inline Mat<zz_p> mul_left(const Mat<zz_p>& input) const
    {
        Mat<zz_p> output;
        mul_left(output, input);
        return output;
    }
    
    /*------------------------------------------------------------*/
    /* computes output = M^t*input, without the diagonal          */
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
    void to_dense(Mat<zz_p>& M) const;
    inline Mat<zz_p> to_dense()
    {
        Mat<zz_p> dense;
        to_dense(dense);
        return dense;
    }

    /*------------------------------------------------------------*/
    /* X and Y are only built if needed                           */
    /*------------------------------------------------------------*/
    void build_X_Y();
};


/*------------------------------------------------------------*/
/* computes                                                   */
/* 1/(u1-v1 rho^(-m+1)) ... 1/(u1-v1 rho^(n-1))               */
/* these are the entries of the toeplitz matrix               */
/* (with m rows and n columns)                                */
/*------------------------------------------------------------*/
void prepare_inverses_cauchy(Vec<zz_p>& inverses, const zz_p& u1, const zz_p& v1, const zz_p& rho, long m, long n);



/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* Cauchy-like matrices on geometric progressions             */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

class cauchy_like_geometric_lzz_p
{
public:
    cauchy_geometric_lzz_p C;
    Mat<zz_p> G, H;

    /*------------------------------------------------------------*/
    /* default constructor                                        */
    /*------------------------------------------------------------*/
    cauchy_like_geometric_lzz_p();

    /*------------------------------------------------------------*/
    /* constructor                                                */
    /*------------------------------------------------------------*/
    cauchy_like_geometric_lzz_p(const Mat<zz_p>& U, const Mat<zz_p>& V, const zz_p& a1, const zz_p& b1, const zz_p& rho);

    /*------------------------------------------------------------*/
    /* dimensions                                                 */
    /*------------------------------------------------------------*/
    long NumRows() const;
    long NumCols() const;
    long NumGens() const;

    /*------------------------------------------------------------*/
    /* computes output = M*input                                  */
    /*------------------------------------------------------------*/
    void mul_right(Vec<zz_p>& output, const Vec<zz_p>& input) const;
    void mul_right_direct(Mat<zz_p> & out, const Mat<zz_p> & in) const ;
    void mul_right_sigma_UL(Mat<zz_p> & out, const Mat<zz_p> & in) ;
    void mul_right(Mat<zz_p>& output, const Mat<zz_p>& input) ;

    /*------------------------------------------------------------*/
    /* computes output = M^t*input                                */
    /*------------------------------------------------------------*/
    void mul_left(Vec<zz_p>& output, const Vec<zz_p>& input) const;
    void mul_left(Mat<zz_p>& output, const Mat<zz_p>& input) const;

    /*------------------------------------------------------------*/
    /* M as a dense matrix                                        */
    /*------------------------------------------------------------*/
    void to_dense(Mat<zz_p>& M) const;
};

/*------------------------------------------------------------*/
/* iCL = CL^(-1), represented for the dual operator           */
/* O(alpha n^2) algorithm                                     */
/*------------------------------------------------------------*/
long invert_direct(cauchy_like_geometric_lzz_p& Cinv, const cauchy_like_geometric_lzz_p& C);
long invert_block(cauchy_like_geometric_lzz_p& Cinv, const cauchy_like_geometric_lzz_p& C);
long invert(cauchy_like_geometric_lzz_p& Cinv, const cauchy_like_geometric_lzz_p& C);
long invert_fast(cauchy_like_geometric_lzz_p& Cinv, const cauchy_like_geometric_lzz_p& C, const long thresh = -1);


#endif


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
