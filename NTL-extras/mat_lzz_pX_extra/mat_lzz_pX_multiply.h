#ifndef MAT_LZZ_PX_MULTIPLY__H
#define MAT_LZZ_PX_MULTIPLY__H

/** Multiplication of univariate polynomial matrices over `zz_p`
 *
 * \file mat_lzz_pX_multiply.h
 * \author Seung Gyu Hyun, Vincent Neiger, Eric Schost
 * \version 0.1
 * \date 2018-12-11
 *
 * Functions for computing products and middle products of univariate
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

/** @name General interface for multiplication and middle product
 * \anchor MiddleProduct
 *
 * These are the general functions for polynomial matrix multiplication and
 * polynomial matrix middle product: they choose the fastest of the available
 * methods, depending on the degrees and dimensions of the input matrices.
 *
 * The middle product `b` of `a` and `c` with respect to nonnegative integers
 * `dA` and `dB` is defined here as:
 * `b == trunc(trunc(a, dA+1)*c div x^dA, dB+1)`
 * where `div` means that we divide by `x^dA`, discarding terms of negative
 * degree. All functions below which have `middle_product` in their name follow
 * this definition.
 *
 * In the multiplication functions, the OUT parameter `c` may alias the IN
 * parameter `a` or the IN parameter `b`; in the middle product functions, the
 * OUT parameter `b` may alias the IN parameter `a` or the IN parameter `c`. 
 *
 * \todo Currently, the code has not been optimized much for matrix-vector and
 * vector-matrix products, or for the case where the left-operand is a column
 * vector.
 */
//@{

/** Computes the product `c = a*b` of two polynomial matrices `a` and `b`. The
 * parameter `is_prime` is set to 1 (the default) if the modulus is known to be
 * prime */
void multiply(
              Mat<zz_pX> & c,
              const Mat<zz_pX>& a,
              const Mat<zz_pX>& b,
              long is_prime = 1
             );

/** Computes the product `c = a*b` of a polynomial matrix `a` by a polynomial
 * column vector `b`. The parameter `is_prime` is set to 1 (the default) if the
 * modulus is known to be prime */
void multiply(Vec<zz_pX>& c, const Mat<zz_pX> & a, const Vec<zz_pX> & b, long is_prime = 1);

/** Computes the product `c = a*b` of a polynomial row vector `a` by a
 * polynomial matrix `b`. The parameter `is_prime` is set to 1 (the default) if
 * the modulus is known to be prime */
void multiply(Vec<zz_pX> & c, const Vec<zz_pX> & a, const Mat<zz_pX> & b, long is_prime = 1);

/** Computes and returns the product `a*b` of two polynomial matrices `a` and `b` */
inline Mat<zz_pX> operator*(const Mat<zz_pX>& a, const Mat<zz_pX>& b)
{ Mat<zz_pX> c; multiply(c, a, b); return c; }

/** Computes and returns the product `a*b` of a polynomial matrix `a` by a
 * polynomial column vector `b` */
inline Vec<zz_pX> operator*(const Mat<zz_pX> & a, const Vec<zz_pX> & b)
{ Vec<zz_pX> c; multiply(c, a, b); return c; }

/** Computes and returns the product `a*b` of a polynomial row vector `a` by a
 * polynomial matrix `b` */
inline Vec<zz_pX> operator*(const Vec<zz_pX> & a, const Mat<zz_pX> & b)
{ Vec<zz_pX> c; multiply(c, a, b); return c; }

/** Computes the truncated product `c = a*b % x^d` of two polynomial matrices
 * `a` and `b` at order `d`. The parameter `is_prime` is set to 1 (the default)
 * if the modulus is known to be prime */
inline void mul_trunc(
                      Mat<zz_pX> & c, 
                      const Mat<zz_pX> & a, 
                      const Mat<zz_pX> & b, 
                      long d, 
                      long is_prime = 1
                     )
{ multiply(c, a, b, is_prime); trunc(c, c, d); }

/** Computes and returns the truncated product `a*b % x^d` of two polynomial
 * matrices `a` and `b` at order `d`. The parameter `is_prime` is set to 1 (the
 * default) if the modulus is known to be prime */
inline Mat<zz_pX> mul_trunc(
                            const Mat<zz_pX> & a, 
                            const Mat<zz_pX> & b, 
                            long n,
                            long is_prime = 1
                           )
{ Mat<zz_pX> c; mul_trunc(c, a, b, n, is_prime); return c; }

/** Computes the middle product `b` of polynomial matrices `a` and `c` with
 * respective to nonnegative integers `dA` and `dB` (see @ref MiddleProduct).
 * The parameter `is_prime` is set to 1 (the default) if the modulus is known
 * to be prime */
void middle_product(
                    Mat<zz_pX> & b,
                    const Mat<zz_pX> & a,
                    const Mat<zz_pX> & c,
                    long dA,
                    long dB,
                    long is_prime = 1
                   );

/** Computes and returns the middle product of polynomial matrices `a` and `c`
 * with respective to nonnegative integers `dA` and `dB` (see @ref
 * MiddleProduct). The parameter `is_prime` is set to 1 (the default) if the
 * modulus is known to be prime */
inline Mat<zz_pX> middle_product(
                                 const Mat<zz_pX>& a,
                                 const Mat<zz_pX>& c,
                                 long dA,
                                 long dB,
                                 long is_prime = 1
                                )
{ Mat<zz_pX> b; middle_product(b, a, c, dA, dB, is_prime); return b; }

//@} // doxygen group: General interface for multiplication and middle product






/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*                   MATRIX MULTIPLICATION                    */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/** @name Cubic-time polynomial matrix multiplication
 *
 * These functions use algorithms of complexity cubic in the dimension, and
 * relying on polynomial arithmetic.
 */
//@{

/** Uses the naive cubic matrix multiplication algorithm (three for loops)
 * along with NTL's polynomial multiplication for multiplying the entries of
 * the matrices */
void multiply_naive(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);

/** Uses Waksman's cubic matrix multiplication algorithm along with NTL's
 * polynomial multiplication. */
void multiply_waksman(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);

//@} // doxygen group: Cubic-time polynomial matrix multiplication

/** @name Small degree polynomial matrix multiplication
 * 
 * Algorithms specialized for small degree matrices.
 */
//@{
/** \todo doc  */
void multiply_transform_naive(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);
/** \todo doc  */
void multiply_transform_karatsuba(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);
/** \todo doc  */
void multiply_transform_montgomery3(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);
/** \todo doc  */
void multiply_transform_karatsuba4(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);
/** \todo doc  */
void multiply_transform(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b, long len);
/** \todo doc  */
inline void multiply_transform(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b)
{ multiply_transform(c, a, b, max(deg(a), deg(b)) + 1); }

/**
 * Uses the algorithm of [Doliskani, Giorgi, Lebreton, Schost. 2018], which
 * relies on matrix multiplication with Vandermonde matrices for evaluation and
 * interpolation.
 *                                    
 * \todo check field is large enough to take the points
 */
void multiply_evaluate_dense(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);
/**
 * Uses the algorithm of [Doliskani, Giorgi, Lebreton, Schost. 2018], with
 * relies on matrix multiplication with Vandermonde matrices for evaluation and
 * interpolation; uses points 1,-1,2,-2, .. in order to speed-up (similar to
 * first step of an FFT). 
 *                                    
 * \todo check field is large enough to take the points
 */
void multiply_evaluate_dense2(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);

/** Computes the Vandermonde matrices for evaluation and interpolation in the
 * [Doliskani et al] approach. */
void vandermonde(Mat<zz_p>& small_vdm1, Mat<zz_p>& small_vdm2, Mat<zz_p>& inv_vdm, long d1, long d2);
/** Computes the Vandermonde matrices for evaluation and interpolation in the
 * [Doliskani et al] approach speeded-up with opposite points. */
void vandermonde2(Mat<zz_p>& small_vdm1, Mat<zz_p>& small_vdm2, Mat<zz_p>& inv_vdm, long d1, long d2);

//@} // doxygen group: Small degree polynomial matrix multiplication


/** @name Evaluation/interpolation polynomial matrix multiplication (FFT points)
 *
 *  \todo doc
 */
//@{

/*------------------------------------------------------------*/
/* assumes FFT prime and p large enough                       */
/* chooses one of the two above                               */
/*------------------------------------------------------------*/
void multiply_evaluate_FFT(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);

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

// TODO experimental function, which happens to be efficient in small dimensions (also depends on the bitsize)
void multiply_evaluate_FFT_direct_no_ll(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);

//@} // doxygen group: Evaluation/interpolation polynomial matrix multiplication (FFT points)

/** @name Evaluation/interpolation polynomial matrix multiplication (non-FFT points)
 *
 *  
 */
//@{
/*------------------------------------------------------------*/
/* 3 primes CRT algorithm                                     */
/*------------------------------------------------------------*/
void multiply_3_primes(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);


/*------------------------------------------------------------*/
/* geometric evaluation                                       */
/* uses Mat<zz_p> matrix multiplication                       */
/* Note: implementation not using matmul always slower.       */
/*------------------------------------------------------------*/
void multiply_evaluate_geometric(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);
// TODO versions similar to matmul2 and matmul3, cf FFT file

//@} // doxygen group: Evaluation/interpolation polynomial matrix multiplication (non-FFT points)



/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*                     MIDDLE PRODUCT                         */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/** @name Middle product
 *
 * \todo doc
 *
 * \todo ensure degree bounds on a, c (?)
 */
//@{

/** Uses the naive cubic matrix multiplication algorithm (three for loops)
 * along with NTL's middle product */
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
/* transpose of b mapsto c = a*b. output is                   */
/*    trunc( rev(a, dA)*c div x^dA, dB+1 )                    */
/* a must have degree at most dA                              */
/* c must have degree at most dA + dB                         */
/*------------------------------------------------------------*/
void t_multiply_evaluate_geometric(Mat<zz_pX> & b, const Mat<zz_pX> & a, const Mat<zz_pX> & c, long dA, long dB);

//@} // doxygen group: Middle product



/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*                CLASS FOR 3 PRIMES FFTS                     */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/** Class for 3-primes FFT
 *
 * \todo use smaller primes when/if possible? (currently uses the FFT primes 0,
 * 1, 2 of NTL)
 */
class lzz_pX_3_primes
{
public:
    /** Constructor, which uses the product inner dimension and the degrees to
     * choose the number of primes and the primes */
    lzz_pX_3_primes(long indim, long dA, long dB);

    /** Empty constructor is forbidden */
    lzz_pX_3_primes() = delete;

    /** Accessor: number of primes */
    long nb() const { return nb_primes; }

    /** Reconstructs c from its images modulo the primes */
    void reconstruct(Mat<zz_pX>& c, const Vec<Mat<zz_pX>>& cs) const;

private:
    long nb_primes; /**< number of FFT primes */
    long fft_p0, fft_p1, fft_p2; /**< the FFT primes */
};





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

    mat_lzz_pX_lmultiplier_3_primes(const Mat<zz_pX> & a, long dB);
    mat_lzz_pX_lmultiplier_3_primes() = delete;

    /*------------------------------------------------------------*/
    /* we use unique_ptrs; we don't expect to have to do copies   */
    /*------------------------------------------------------------*/
    mat_lzz_pX_lmultiplier_3_primes& operator=(const mat_lzz_pX_lmultiplier_3_primes&) = delete;

private:
    lzz_pX_3_primes primes;
    Vec<std::unique_ptr<mat_lzz_pX_lmultiplier>> FFT_muls;
};

/** @name Helper functions
 *
 */
//@{

/** Compute the reduction of a polynomial matrix `a` modulo the current prime,
 * in-place (used in 3-primes FFT) */
void reduce_mod_p(Mat<zz_pX> & a);

/** Compute `amodp`, the reduction of a polynomial matrix `a` modulo the
 * current prime (used in 3-primes FFT) */
void reduce_mod_p(Mat<zz_pX> & amodp, const Mat<zz_pX> & a);

/** \todo DOC returns a multiplier of the right type */
std::unique_ptr<mat_lzz_pX_lmultiplier>
get_lmultiplier(const Mat<zz_pX> & a, long dB);

//@} // doxygen group: Helper functions

#endif /* ifndef MAT_LZZ_PX_MULTIPLY__H */

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
