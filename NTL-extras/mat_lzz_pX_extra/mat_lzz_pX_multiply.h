#ifndef MAT_LZZ_PX_MULTIPLY__H
#define MAT_LZZ_PX_MULTIPLY__H

/** Multiplication of polynomial matrices.
 *
 * \file mat_lzz_pX_multiply.h
 * \author Seung Gyu Hyun, Vincent Neiger, Eric Schost
 * \version 0.1
 * \date 2018-12-11
 *
 * Functions for computing products and middle products of two unvariate
 * polynomial matrices.
 *
 */

#include <memory> // for unique_ptr
#include "mat_lzz_pX_utils.h" // for deg() (and maybe others)
#include "lzz_pX_CRT.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*             MULTIPLICATION FUNCTIONS                       */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*               A CLASS FOR 3 PRIMES FFTS                    */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

class lzz_pX_3_primes
{
public:
    /*------------------------------------------------------------*/
    /* constructor of lzz_p_3_primes                              */
    /* need to know ncols and degrees to choose number of primes  */
    /*------------------------------------------------------------*/
    lzz_pX_3_primes(long ncols, long dA, long dB);
    lzz_pX_3_primes(){};

    /*------------------------------------------------------------*/
    /* returns the number of primes                               */
    /*------------------------------------------------------------*/
    long nb() const;

    /*------------------------------------------------------------*/
    /* reconstructs c from its images                             */
    /*------------------------------------------------------------*/
    void reconstruct(Mat<zz_pX>& c, const Vec<Mat<zz_pX>>& cs);

private:
    long nb_primes;
    long fft_p0, fft_p1, fft_p2; // the fft primes
};

/*------------------------------------------------------------*/
/* 3-primes: in-place reduction modulo the current prime      */
/*------------------------------------------------------------*/
void reduce_mod_p(Mat<zz_pX> & a);


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*               PLAIN MULTIPLICATION                         */
/* ALL FUNCTIONS: C = A*B, OUTPUT CAN ALIAS INPUT             */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* naive algorithm                                            */
/*------------------------------------------------------------*/
void multiply_naive(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);

/*------------------------------------------------------------*/
/* Waksman's algorithm                                        */
/*------------------------------------------------------------*/
void multiply_waksman(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);

/*------------------------------------------------------------*/
/* assumes FFT prime and p large enough                       */
/* uses Mat<zz_p> matrix multiplication                       */
/*------------------------------------------------------------*/
void multiply_evaluate_FFT_matmul1(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);
void multiply_evaluate_FFT_matmul2(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);
void multiply_evaluate_FFT_matmul3(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);

/*------------------------------------------------------------*/
/* assumes FFT prime and p large enough                       */
/* does not use Mat<zz_p> matrix multiplication               */
/*------------------------------------------------------------*/
void multiply_evaluate_FFT_direct(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);

/*------------------------------------------------------------*/
/* assumes FFT prime and p large enough                       */
/* chooses one of the two above                               */
/*------------------------------------------------------------*/
void multiply_evaluate_FFT(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);

/*------------------------------------------------------------*/
/* geometric evaluation                                       */
/* uses Mat<zz_p> matrix multiplication                       */
/* Note: implementation not using matmul always slower.       */
/*------------------------------------------------------------*/
void multiply_evaluate_geometric(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);

/*------------------------------------------------------------*/
/* matrix multiplication using the algorithm of Giorgi et al. */
/* uses matrix multiplication for evaluation and interpolation*/
/*------------------------------------------------------------*/
void multiply_evaluate_dense(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);
/*------------------------------------------------------------*/
/* computes the matrices for evaluation and interpolation     */
/*------------------------------------------------------------*/
void vandermonde(Mat<zz_p>& small_vdm1, Mat<zz_p>& small_vdm2, Mat<zz_p>& inv_vdm, long d1, long d2);

void multiply_transform_naive(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);
void multiply_transform_karatsuba(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);
void multiply_transform_montgomery3(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);
void multiply_transform_karatsuba4(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);
void multiply_transform(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b, long len);
inline void multiply_transform(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b)
{
    multiply_transform(c, a, b, max(deg(a), deg(b)) + 1);
}

/*------------------------------------------------------------*/
/* 3 primes CRT algorithm                                     */
/*------------------------------------------------------------*/
void multiply_3_primes(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);

/*------------------------------------------------------------*/
/* main function for c = a*b                                  */
/* is_prime = 1 assumes that p is known to be prime           */
/*------------------------------------------------------------*/
void multiply(Mat<zz_pX> & c, const Mat<zz_pX>& a, const Mat<zz_pX>& b, long is_prime = 1);

inline Mat<zz_pX> operator*(const Mat<zz_pX>& a, const Mat<zz_pX>& b)
{ 
    Mat<zz_pX> x; 
    multiply(x, a, b); 
    return x; 
}

/*------------------------------------------------------------*/
/* multiply by a vector                                       */
/*------------------------------------------------------------*/
void multiply(Vec<zz_pX>& c, const Mat<zz_pX>& a, const Vec<zz_pX>& b, long is_prime = 1);

inline Vec<zz_pX> operator*(const Mat<zz_pX>& a, const Vec<zz_pX>& b)
{ 
    Vec<zz_pX> x; 
    multiply(x, a, b); 
    return x; 
}


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*               TRUNCATED MULTIPLICATION                     */
/* ALL FUNCTIONS: C = A*B MOD X^N, OUTPUT CAN ALIAS INPUT     */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* c = a*b mod x^n                                            */
/*------------------------------------------------------------*/
inline void mul_trunc(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b, long n, long is_prime = 1)
{
    multiply(c, a, b, is_prime);
    trunc(c, c, n);
}

inline Mat<zz_pX> mul_trunc(const Mat<zz_pX> & a, const Mat<zz_pX> & b, long n, long is_prime = 1)
{
    Mat<zz_pX> c;
    mul_trunc(c, a, b, n, is_prime);
    return c;
}

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*                     MIDDLE PRODUCT                         */
/* all(*) functions (except for geometric)                    */
/* return trunc( trunc(a, dA+1)*c div x^dA, dB+1 )            */
/* output can alias input                                     */
/* todo: ensure degree bounds on a, c                         */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/*------------------------------------------------------------*/
/* transpose of b mapsto c = a*b. output is                   */
/*    trunc( rev(a, dA)*c div x^dA, dB+1 )                    */
/* a must have degree at most dA                              */
/* c must have degree at most dA + dB                         */
/*------------------------------------------------------------*/
void t_multiply_evaluate_geometric(Mat<zz_pX> & b, const Mat<zz_pX> & a, const Mat<zz_pX> & c, long dA, long dB);

/*------------------------------------------------------------*/
/* naive algorithm, uses polynomial middle products           */
/*------------------------------------------------------------*/
void middle_product_naive(Mat<zz_pX> & b, const Mat<zz_pX> & a, const Mat<zz_pX> & c, long dA, long dB);

/*------------------------------------------------------------*/
/* 3 primes CRT algorithm                                     */
/*------------------------------------------------------------*/
void middle_product_3_primes(Mat<zz_pX> & b, const Mat<zz_pX> & a, const Mat<zz_pX> & c, long dA, long dB);

/*------------------------------------------------------------*/
/* assumes FFT prime and p large enough                       */
/* does not use Mat<zz_p> matrix multiplication               */
/*------------------------------------------------------------*/
void middle_product_evaluate_FFT_direct(Mat<zz_pX> & b, const Mat<zz_pX> & a, const Mat<zz_pX> & c, long dA, long dB);

/*------------------------------------------------------------*/
/* assumes FFT prime and p large enough                       */
/* uses Mat<zz_p> matrix multiplication                       */
/*------------------------------------------------------------*/
void middle_product_evaluate_FFT_matmul(Mat<zz_pX> & b, const Mat<zz_pX> & a, const Mat<zz_pX> & c, long dA, long dB);

/*------------------------------------------------------------*/
/* assumes FFT prime and p large enough                       */
/* chooses one of the two above                               */
/*------------------------------------------------------------*/
void middle_product_evaluate_FFT(Mat<zz_pX> & b, const Mat<zz_pX> & a, const Mat<zz_pX> & c, long dA, long dB);

/*------------------------------------------------------------*/
/* uses matrix multiplication for evaluation and interpolation*/
/*------------------------------------------------------------*/
void middle_product_evaluate_dense(Mat<zz_pX> & b, const Mat<zz_pX> & a, const Mat<zz_pX> & c, long dA, long dB);

/*------------------------------------------------------------*/
/* main function.                                             */
/* is_prime = 1 assumes that p is known to be prime           */
/*------------------------------------------------------------*/
/** TODO */
void middle_product(Mat<zz_pX> & b, const Mat<zz_pX> & a, const Mat<zz_pX> & c, long dA, long dB, long is_prime = 1);

/** TODO */
inline Mat<zz_pX> middle_product(const Mat<zz_pX>& a, const Mat<zz_pX>& c, long dA, long dB, long is_prime = 1)
{ 
    Mat<zz_pX> b; 
    middle_product(b, a, c, dA, dB, is_prime); 
    return b; 
}


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* CLASSES FOR MULTIPLICATION WITH A GIVEN L.H.S.             */
/* CONSTRUCTORS TAKE AN ARGUMENT dB ST R.H.S HAS DEGREE <= dB */  
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
class mat_lzz_pX_lmultiplier
{
public:
    /*------------------------------------------------------------*/
    /* c = M * b                                                  */
    /*------------------------------------------------------------*/
    virtual void multiply(Mat<zz_pX>& c, const Mat<zz_pX>& b) = 0;
    virtual ~mat_lzz_pX_lmultiplier(){}

    inline Mat<zz_pX> multiply(const Mat<zz_pX>& b)
    {
        Mat<zz_pX> c;
        multiply(c, b);
        return c;
    }
    
    /*------------------------------------------------------------*/
    /* getters                                                    */
    /*------------------------------------------------------------*/
    long NumRows() const; // dimensions
    long NumCols() const; // dimensions
    long degA() const; // degree of current matrix
    long degB() const; // max degree of rhs

protected:
    long __s, __t; // dimensions
    long __dA; // degree of current matrix
    long __dB; // max degree of rhs
};



/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* for use with FFT primes; direct product                    */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
class mat_lzz_pX_lmultiplier_FFT_direct : public mat_lzz_pX_lmultiplier
{
public:
    /*------------------------------------------------------------*/
    /* c = M * b                                                  */
    /*------------------------------------------------------------*/
    void multiply(Mat<zz_pX>& c, const Mat<zz_pX>& b);

    mat_lzz_pX_lmultiplier_FFT_direct(){}
    mat_lzz_pX_lmultiplier_FFT_direct(const Mat<zz_pX> & a, long dB);

private:
    Vec<Vec<fftRep>> vala;
    long len, n0, K, pr, nb_slices, first_slice;
    sp_reduce_struct red1;
    sp_ll_reduce_struct red2;
};



/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* for use with FFT primes; matmul product                    */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
class mat_lzz_pX_lmultiplier_FFT_matmul : public mat_lzz_pX_lmultiplier
{
public:
    /*------------------------------------------------------------*/
    /* c = M * b                                                  */
    /*------------------------------------------------------------*/
    void multiply(Mat<zz_pX>& c, const Mat<zz_pX>& b);

    mat_lzz_pX_lmultiplier_FFT_matmul(){}
    mat_lzz_pX_lmultiplier_FFT_matmul(const Mat<zz_pX> & a, long dB);

private:
    Vec<Mat<zz_p>> va; // FFT of current matrix
    long idxk; // log-size of FFT
};


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* geometric points                                           */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
class mat_lzz_pX_lmultiplier_geometric : public mat_lzz_pX_lmultiplier
{
public:
    /*------------------------------------------------------------*/
    /* c = M * b                                                  */
    /*------------------------------------------------------------*/
    void multiply(Mat<zz_pX>& c, const Mat<zz_pX>& b);

    mat_lzz_pX_lmultiplier_geometric(){}
    mat_lzz_pX_lmultiplier_geometric(const Mat<zz_pX> & a, long dB);

private:
    Vec<Mat<zz_p>> va; // FFT of current matrix
    zz_pX_Multipoint_Geometric ev;
};

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* dense algorithm                                            */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
class mat_lzz_pX_lmultiplier_dense : public mat_lzz_pX_lmultiplier
{
public:
    /*------------------------------------------------------------*/
    /* c = M * b                                                  */
    /*------------------------------------------------------------*/
    void multiply(Mat<zz_pX>& c, const Mat<zz_pX>& b);

    mat_lzz_pX_lmultiplier_dense(){}
    mat_lzz_pX_lmultiplier_dense(const Mat<zz_pX> & a, long dB);

private:
    Mat<zz_p> vA, vB, iV, valA;
    long nb_points;
};

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* 3 primes                                                   */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
class mat_lzz_pX_lmultiplier_3_primes : public mat_lzz_pX_lmultiplier
{
public:
    /*------------------------------------------------------------*/
    /* c = M * b                                                  */
    /*------------------------------------------------------------*/
    void multiply(Mat<zz_pX>& c, const Mat<zz_pX>& b);

    mat_lzz_pX_lmultiplier_3_primes(){}
    mat_lzz_pX_lmultiplier_3_primes(const Mat<zz_pX> & a, long dB);


    /*------------------------------------------------------------*/
    /* we use unique_ptrs; we don't expect to have to do copies   */
    /*------------------------------------------------------------*/
    mat_lzz_pX_lmultiplier_3_primes& operator=(const mat_lzz_pX_lmultiplier_3_primes&)
    {
        LogicError("no copy allowed");
        return *this;
    }

private:
    lzz_pX_3_primes primes;
    Vec<std::unique_ptr<mat_lzz_pX_lmultiplier>> FFT_muls;
};


/*------------------------------------------------------------*/
/* returns a multiplier of the right type                     */
/*------------------------------------------------------------*/
std::unique_ptr<mat_lzz_pX_lmultiplier> get_lmultiplier(const Mat<zz_pX> & a, long dB);

#endif /* ifndef MAT_LZZ_PX_MULTIPLY__H */

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
