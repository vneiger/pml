#ifndef MAT_LZZ_PX_MULTIPLY__H
#define MAT_LZZ_PX_MULTIPLY__H

/** \brief Multiplication of univariate polynomial matrices over `zz_p`
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
#include "mat_lzz_pX_utils.h" // for deg(), macros CACHE_LINE_SIZE (and maybe other things)
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
 * to be prime.
 *
 * \todo Currently, if `deg(a)<dA` then we right shift `c` by `dA-deg(a)` and
 * then apply this function with `a`, the shifted `c`, `dA-deg(a)` and `dB`.
 * We should do some tests and timings to see if this really has an interest
 * when `deg(a)` is very close to `dA` (still, currently with the jumps of FFT
 * at powers of 2, the question is mostly whether this argument reduction will
 * allow us to work with a smaller FFT size...). */
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
/** Computes `c = a*b`. \todo short algorithm description  */
void multiply_transform_naive(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);
/** Computes `c = a*b`. \todo short algorithm description  */
void multiply_transform_karatsuba(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);
/** Computes `c = a*b`. \todo short algorithm description  */
void multiply_transform_montgomery3(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);
/** Computes `c = a*b`. \todo short algorithm description  */
void multiply_transform_karatsuba4(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);
/** Computes `c = a*b`. \todo short algorithm description  */
void multiply_transform(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b, long len);
/** Computes `c = a*b`. \todo short algorithm description  */
inline void multiply_transform(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b)
{ multiply_transform(c, a, b, max(deg(a), deg(b)) + 1); }

/**
 * Computes `c = a*b`, using a polynomial variant of the algorithm of
 * [Doliskani, Giorgi, Lebreton, Schost, 2018], which relies on matrix
 * multiplication with Vandermonde matrices for evaluation and interpolation.
 *                                    
 * \todo check field is large enough to take the points
 */
void multiply_evaluate_dense(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);
/**
 * Computes `c = a*b`, using a polynomial variant of the algorithm of
 * [Doliskani, Giorgi, Lebreton, Schost, 2018], which relies on matrix
 * multiplication with Vandermonde matrices for evaluation and interpolation;
 * uses points 1,-1,2,-2, .. in order to speed-up (similar to first step of an
 * FFT). 
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
 *  Polynomial matrix multiplication based on evaluation/interpolation at FFT
 *  points, when the modulus is an FFT prime.
 *
 *  \todo try to make the change of representations more cache-friendly (e.g. use
 *  some block size 32)
 */
//@{

/** Computes `c = a*b`. Choose the expected fastest (according to thresholds)
 * function for matrix multiplication using FFT evaluation/interpolation. */
void multiply_evaluate_FFT(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);

/** Computes `c = a*b`. \todo short algorithm description  */
void multiply_evaluate_FFT_matmul1(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);
/** Computes `c = a*b`. \todo short algorithm description  */
void multiply_evaluate_FFT_matmul2(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);
/** Computes `c = a*b`. \todo short algorithm description  */
void multiply_evaluate_FFT_matmul3(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);

/** Computes `c = a*b` via FFT evaluation/interpolation, but without relying on
 * `Mat<zz_p>` multiplication */
void multiply_evaluate_FFT_direct_ll_type(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);

/** Computes `c = a*b` via FFT evaluation/interpolation, but without relying on
 * `Mat<zz_p>` multiplication.
 *
 * \todo explain difference with other similar one.
 */
void multiply_evaluate_FFT_direct(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);

//@} // doxygen group: Evaluation/interpolation polynomial matrix multiplication (FFT points)

/** @name Evaluation/interpolation polynomial matrix multiplication (non-FFT points)
 */
//@{

/** Computes `c = a*b` via 3-primes CRT algorithm */
void multiply_3_primes(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);

/** Computes `c = a*b` via evaluation/interpolation at geometric points,
 * relying on Mat<zz_p> matrix multiplication.
 *
 * \todo since it seems to most often be slower than 3primes, this has
 * not been optimized (several variants, cache efficiency, etc). This
 * should be re-examined after improvements in geometric matrix
 * evaluation.
 */
void multiply_evaluate_geometric(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b);

//@} // doxygen group: Evaluation/interpolation polynomial matrix multiplication (non-FFT points)

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*                     MIDDLE PRODUCT                         */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/** @name Middle product
 *
 * \todo improve doc (better separate them as above for multiplication?)
 *
 * \todo ensure degree bounds on a, c (?)
 */
//@{

/** Uses the naive cubic matrix multiplication algorithm (three for loops)
 * along with NTL's middle product.  Requires `deg(a) <= dA` and `deg(c) <=
 * dA+dB`. */
void middle_product_naive(Mat<zz_pX> & b, const Mat<zz_pX> & a, const Mat<zz_pX> & c, long dA, long dB);

/** Uses the 3 primes CRT algorithm. Requires `deg(a) <= dA` and `deg(c) <= dA+dB`. */
void middle_product_3_primes(Mat<zz_pX> & b, const Mat<zz_pX> & a, const Mat<zz_pX> & c, long dA, long dB);

/** Based on evaluation/interpolation at FFT points; assumes FFT prime and p
 * large enough; does not use Mat<zz_p> matrix multiplication. Requires `deg(a)
 * <= dA` and `deg(c) <= dA+dB`.
 */
void middle_product_evaluate_FFT_direct_ll_type(Mat<zz_pX> & b, const Mat<zz_pX> & a, const Mat<zz_pX> & c, long dA, long dB);
void middle_product_evaluate_FFT_direct(Mat<zz_pX> & b, const Mat<zz_pX> & a, const Mat<zz_pX> & c, long dA, long dB);

/** Based on evaluation/interpolation at FFT points; assumes FFT prime and p
 * large enough; relies on Mat<zz_p> matrix multiplication. Requires `deg(a)
 * <= dA` and `deg(c) <= dA+dB`.
 */
void middle_product_evaluate_FFT_matmul(Mat<zz_pX> & b, const Mat<zz_pX> & a, const Mat<zz_pX> & c, long dA, long dB);
void middle_product_evaluate_FFT_matmul1(Mat<zz_pX> & b, const Mat<zz_pX> & a, const Mat<zz_pX> & c, long dA, long dB);
void middle_product_evaluate_FFT_matmul2(Mat<zz_pX> & b, const Mat<zz_pX> & a, const Mat<zz_pX> & c, long dA, long dB);
void middle_product_evaluate_FFT_matmul3(Mat<zz_pX> & b, const Mat<zz_pX> & a, const Mat<zz_pX> & c, long dA, long dB);

/** Based on evaluation/interpolation at FFT points; assumes FFT prime and p
 * large enough; chooses the expected fastest of the available algorithm.
 * Requires `deg(a) <= dA` and `deg(c) <= dA+dB`. 
 *
 * \todo improve/automatic thresholds
 */
void middle_product_evaluate_FFT(Mat<zz_pX> & b, const Mat<zz_pX> & a, const Mat<zz_pX> & c, long dA, long dB);
void middle_product_evaluate_FFT_new(Mat<zz_pX> & b, const Mat<zz_pX> & a, const Mat<zz_pX> & c, long dA, long dB);

/** Based on evaluation/interpolation, using multiplication by Vandermonde
 * matrices (fast when degree is small and dimension not so small).  Requires
 * `deg(a) <= dA` and `deg(c) <= dA+dB`.  */
void middle_product_evaluate_dense(Mat<zz_pX> & b, const Mat<zz_pX> & a, const Mat<zz_pX> & c, long dA, long dB);

// /** Based on evaluation/interpolation, using multiplication by Vandermonde
// * matrices for evaluation and interpolation (fast when degree is small and
// * dimension not so small). This version uses points 1,-1,2,-2, .. in order to
// * speed-up (similar to first step of an FFT). Requires `deg(a) <= dA` and `deg(c) <= dA+dB`.
// *                                    
// * \todo check field is large enough to take the points
// */
//void middle_product_evaluate_dense2(Mat<zz_pX> & b, const Mat<zz_pX> & a, const Mat<zz_pX> & c, long dA, long dB);

/** Transpose of b mapsto c = a*b. Output is
 *    trunc( rev(a, dA)*c div x^dA, dB+1 )
 * a must have degree at most dA
 * c must have degree at most dA + dB
 */
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





/** Abstract class for left-multiplication with precomputation
 *
 * This stores values for the left-hand side, which will speed up repeated
 * multiplications. In derived classes, constructors take an argument `dB` such
 * that any right-hand side when performing multiplication has degree at most
 * `dB`.
 *
 * \todo make all of this more modular: ideally this should not copy code for
 * multiplication; this is not convenient for maintaining or improving
 * (related: right now, this can be accelerated by importing recent additions
 * and changes and better thresholds for multiplication...)
 *
 * \todo change accessor degA to a less specific name
 *
 * \todo is there an easy way to avoid virtual method here? (it can
 * have consequences such as bad pipelining / branch predicting, and
 * bad cache efficiency). CRTP does not seem an option with the current
 * design. Maybe re-design...
 */
class mat_lzz_pX_lmultiplier
{
public:
    /** destructor */
    virtual ~mat_lzz_pX_lmultiplier(){}

    /** Constructor with the left-hand side matrix and the degree bound
     * for right-hand side to be multiplied */
    mat_lzz_pX_lmultiplier(const Mat<zz_pX> & a, long dB) :
        __s(a.NumRows()), __t(a.NumCols()), __dB(dB)
    { __dA = deg(a); }

    /** Computes `c = this * b`. Must be defined in non-abstract derived
     * classes. */
    virtual void multiply(Mat<zz_pX>& c, const Mat<zz_pX>& b) = 0;

    /** Computes and returns `this * b`. */
    inline Mat<zz_pX> multiply(const Mat<zz_pX>& b)
    {
        Mat<zz_pX> c;
        multiply(c, b);
        return c;
    }

    /** Accessor: number of rows */
    long NumRows() const { return __s; }
    /** Accessor: number of columns */
    long NumCols() const { return __t; }
    /** Accessor: degree of this */
    long degA() const { return __dA; }
    /** Accessor: maximum degree for right-hand side */
    long degB() const { return __dB; }

protected:
    long __s, __t;  /**< number of rows and columns */
    long __dA;  /**< degree of this */
    long __dB;  /**< maximum degree for right-hand side */
};


/** Class for left-multiplication with precomputation, FFT prime,
 * without using `Mat<zz_p>` multiplication.
 *
 * See documentation for parent class
 *
 * \todo improve with more variants of multiply
 * \todo improve doc for attributes
 */
class mat_lzz_pX_lmultiplier_FFT_direct : public mat_lzz_pX_lmultiplier
{
public:
    /** No default constructor */
    mat_lzz_pX_lmultiplier_FFT_direct() = delete;

    /** Constructor with the left-hand side matrix and the degree bound
     * for right-hand side to be multiplied */
    mat_lzz_pX_lmultiplier_FFT_direct(const Mat<zz_pX> & a, long dB);

    /** Computes `c = this * b`. Chooses between two available methods
     * according to thresholds.
     *
     * \todo thresholds copied from regular matrix multiplication. Needs
     * tuning. */
    void multiply(Mat<zz_pX>& c, const Mat<zz_pX>& b);

private:
    void multiply_direct_ll_type(Mat<zz_pX>& c, const Mat<zz_pX>& b);
    void multiply_direct(Mat<zz_pX>& c, const Mat<zz_pX>& b);
    Vec<Vec<fftRep>> vala; /**< FFT evaluations */
    long n0, K, pr, nb_slices, first_slice; /**< precomputations for multiply */
    sp_reduce_struct red1; /**< precomputations for multiply */
    sp_ll_reduce_struct red2; /**< precomputations for multiply */
};


/** Class for left-multiplication with precomputation, FFT prime,
 * using `Mat<zz_p>` multiplication.
 *
 * See documentation for parent class.
 *
 * \todo improve with more variants of multiply
 */
class mat_lzz_pX_lmultiplier_FFT_matmul : public mat_lzz_pX_lmultiplier
{
public:
    /** No default constructor */
    mat_lzz_pX_lmultiplier_FFT_matmul() = delete;

    /** Constructor with the left-hand side matrix and the degree bound
     * for right-hand side to be multiplied */
    mat_lzz_pX_lmultiplier_FFT_matmul(const Mat<zz_pX> & a, long dB);

    /** Computes `c = this * b`. */
    void multiply(Mat<zz_pX>& c, const Mat<zz_pX>& b);

private:
    Vec<Mat<zz_p>> va; /**< FFT evaluations of entries of "this" */
    long idxk; /**< log-size of FFT */
};


/** Class for left-multiplication with precomputation, relying on
 * evaluation/interpolation at geometric points.
 *
 * See documentation for parent class.
 *
 * \todo since it seems to most often be slower than 3primes, this has
 * not been optimized (several variants, cache efficiency, etc). This
 * should be re-examined after improvements in geometric matrix
 * evaluation.
 */
class mat_lzz_pX_lmultiplier_geometric : public mat_lzz_pX_lmultiplier
{
public:
    /** No default constructor */
    mat_lzz_pX_lmultiplier_geometric() = delete;

    /** Constructor with the left-hand side matrix and the degree bound
     * for right-hand side to be multiplied */
    mat_lzz_pX_lmultiplier_geometric(const Mat<zz_pX> & a, long dB);

    /** Computes `c = this * b`. */
    void multiply(Mat<zz_pX>& c, const Mat<zz_pX>& b);

private:
    Vec<Mat<zz_p>> va; /**< evaluations of entries of "this" */
    zz_pX_Multipoint_Geometric ev; /**< Geometric multipoint object for evaluating/interpolating */
};


/** Class for left-multiplication with precomputation, relying on
 * evaluation/interpolation via Vandermonde multiplication.
 *
 * See documentation for parent class.
 *
 * \todo improve with more variants of multiply
 */
class mat_lzz_pX_lmultiplier_dense : public mat_lzz_pX_lmultiplier
{
public:
    /** No default constructor */
    mat_lzz_pX_lmultiplier_dense() = delete;

    /** Constructor with the left-hand side matrix and the degree bound
     * for right-hand side to be multiplied */
    mat_lzz_pX_lmultiplier_dense(const Mat<zz_pX> & a, long dB);

    /** Computes `c = this * b`. */
    void multiply(Mat<zz_pX>& c, const Mat<zz_pX>& b);

private:
    Mat<zz_p> vB, iV; /**< Vandermonde matrix for rhs, and inverse Vandermonde */
    Vec<Mat<zz_p>> valAp; /**< Store evaluations of `a` */
    long nb_points; /**< Number of evaluation points */
};


/** Class for left-multiplication with precomputation, relying on 3-primes
 * multiplication.
 *
 * See documentation for parent class.
 */
class mat_lzz_pX_lmultiplier_3_primes : public mat_lzz_pX_lmultiplier
{
public:
    /** No default constructor */
    mat_lzz_pX_lmultiplier_3_primes() = delete;

    /** Constructor with the left-hand side matrix and the degree bound
     * for right-hand side to be multiplied */
    mat_lzz_pX_lmultiplier_3_primes(const Mat<zz_pX> & a, long dB);

    /** No copies are expected (we use unique_ptrs) */
    mat_lzz_pX_lmultiplier_3_primes& operator=(const mat_lzz_pX_lmultiplier_3_primes&) = delete;

    /** Computes `c = this * b`. */
    void multiply(Mat<zz_pX>& c, const Mat<zz_pX>& b);

private:
    lzz_pX_3_primes primes; /**< 3 primes object storing info on the CRT framework */
    Vec<std::unique_ptr<mat_lzz_pX_lmultiplier>> FFT_muls; /**< 3 FFT left-multipliers */
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

/** returns a multiplier of the right type
 *
 * \todo more explicit doc
 */
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
