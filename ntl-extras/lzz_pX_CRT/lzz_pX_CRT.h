#ifndef LZZ_PX_CRT__H
#define LZZ_PX_CRT__H

/** \brief Chinese remaindering and evaluation / interpolation for polynomials over `zz_p`
 *
 * \file lzz_pX_CRT.h
 * \author Seung Gyu Hyun, Vincent Neiger, Eric Schost
 * \version 0.1
 * \date 2019-02-01
 *
 * This file contains the declarations of helper classes performing 
 * various operations related to CRT techniques, for polynomials over `zz_p`.
 *
 */

#include <map>
#include <NTL/lzz_pX.h>
#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/mat_lzz_p.h>
#include "thresholds_geometric.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* classes for CRT and multipoint evaluation over zz_p        */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/



/*------------------------------------------------------------*/
/* builds the subproduct tree naively (no Huffman)            */
/*------------------------------------------------------------*/
/** @name Misc
 *  @{
 */ 

/** Sets `tree` to be the subproduct tree of polynomials in `q`.
 *  Does not use Huffman techniques to minimize product degrees. 
 */
void build_subproduct_tree(Vec<Vec<zz_pX>> & tree, const Vec<zz_pX> & q);

/** @}  */       // doxygen misc


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* multipoint evaluation over zz_p                            */
/* abstract class                                             */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/**  Abstract class for multipoint evaluation.
 *   
 */
class zz_pX_Multipoint
{
public:
    /** Virtual constructor */
    virtual ~zz_pX_Multipoint(){}

    /*------------------------------------------------------------*/
    /* basic operations                                           */
    /*------------------------------------------------------------*/
    /** @name Basic operations.
     *  Only these need to be virtual and implemented in subclasses.
     *  For evaluation functions, input polynomials must have degree less than `length()`.
     *  For interpolation functions, input vectors must have length `length()`.
     *  @{
     */ 

    /** Sets `val` to the array of values of `P` at the points `this` supports. */
    virtual void evaluate(Vec<zz_p>& val, const zz_pX& P) const = 0;

    /** Sets `P` to the polynomial of degree less than `length()` whose values at the points `this` supports are `val`. 
     *  Not const, since we may use an instance `fftRep` as workspace. 
     */ 
    virtual void interpolate(zz_pX& P, const Vec<zz_p>& val) = 0;     // dirty hack: interpolate uses an instance fftRep as workspace

    /** Transpose of the evaluation function. The optional argument `output_size` 
     *  indicates the required size of the output; default value -1 makes it `length()`
     */
    virtual void t_evaluate(zz_pX& P, const Vec<zz_p>& val, long output_size = -1) const = 0;

    /** Transpose of the interpolation function.
     *  Not const, since we may use an instance `fftRep` as workspace. 
     */ 
    virtual void t_interpolate(Vec<zz_p>& val, const zz_pX& P) = 0;     // dirty hack: t_interpolate uses an instance fftRep as workspace

    /** @}  */       // doxygen basic operations


    /*------------------------------------------------------------*/
    /* action on vectors and matrices                             */
    /*------------------------------------------------------------*/
    /** @name Actions on vectors and matrices.
     *  @{
     */ 

    /** Evaluation of a vector. `val[i]` is the values of `P[i]` at all points we support. */
    void evaluate_vector(Vec<Vec<zz_p>>& val, const Vec<zz_pX>& P) const;

    /** Interpolation of a vector. */
    void interpolate_vector(Vec<zz_pX>& P, const Vec<Vec<zz_p>>& val);

    /** Evaluation of a vector. `val[i][j]` is the values of `P[i][j]` at all points we support. */
    void evaluate_matrix(Vec<Mat<zz_p>>& val, const Mat<zz_pX>& P) const;

    /** Interpolation of a matrix. */
    void interpolate_matrix(Mat<zz_pX>& P, const Vec<Mat<zz_p>>& val);

    /** Transpose evaluation of a vector of vector of values.  The optional argument `output_size` 
     *  indicates the required size of the output; default value -1 makes it `length()`
     */
    void t_evaluate_vector(Vec<zz_pX>& P, const Vec<Vec<zz_p>>& val, long output_size = -1) const;

    /** Transpose interpolation of a vector of polynomials. */
    void t_interpolate_vector(Vec<Vec<zz_p>>& val, const Vec<zz_pX>& P);

    /** Transpose evaluation of a vector of matrix values. The optional argument `output_size` 
     *  indicates the required size of the output; default value -1 makes it `length()`
     */
    void t_evaluate_matrix(Mat<zz_pX>& P, const Vec<Mat<zz_p>>& val, long output_size = -1) const;

    /** Transpose evaluation of a matrix of polynomials. */
    void t_interpolate_matrix(Vec<Mat<zz_p>>& val, const Mat<zz_pX>& P);

    /** @}  */       // doxygen vectors / matrices


    /** @name Getters.
     *  @{
     */ 

    /*------------------------------------------------------------*/
    /* number of points                                           */
    /*------------------------------------------------------------*/
    /** number of points */
    long length() const;

    /*------------------------------------------------------------*/
    /* get the i-th point                                         */
    /*------------------------------------------------------------*/
    /** Sets `pt` to be the `i`th point in the list we support. */
    void get_point(zz_p& pt, long i);  // sets the array pts if needed

    /*------------------------------------------------------------*/
    /* get all points (sets the array pts if needed)              */
    /*------------------------------------------------------------*/
    /** set `points` to be the vector of points we support. */
    void get_points(Vec<zz_p> & points);

    /*------------------------------------------------------------*/
    /* a naive conversion to a dense matrix                       */
    /*------------------------------------------------------------*/
    /** Sets `M` to be the Vandermonde matrix at the points we support. Done naively. */
    void to_dense(Mat<zz_p>& M);

    /** Returns the Vandermonde matrix at the points we support. Done naively. */
    inline Mat<zz_p> to_dense()
    {
        Mat<zz_p> M;
        to_dense(M);
        return M;
    }

    /** @}  */       // doxygen getters

protected:
    /** The array of points we support. Not required to be initialized at creation. */
    Vec<zz_p> pts;

    /** Number of points we support. */
    long n;
};


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* multipoint evaluation over zz_p                            */
/* uses subproduct tree.                                      */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/


/**  Multipoint evaluation for general sets of points. 
 *   Uses subproduct tree techniques.
 */
class zz_pX_Multipoint_General : public zz_pX_Multipoint
{
public:
    /** @name Constructors.
     *  @{
     */ 

    /*------------------------------------------------------------*/
    /* boilerplate                                                */
    /*------------------------------------------------------------*/
    /** Empty constructor; no effect */
    zz_pX_Multipoint_General(){}

    /*------------------------------------------------------------*/
    /* constructor from points                                    */
    /*------------------------------------------------------------*/
    /** Constructor from a list of points; initializes the subproduct tree */
    zz_pX_Multipoint_General(const Vec<zz_p>& q);

    /*------------------------------------------------------------*/
    /* constructor from points[offset]...points[offset+length-1]  */
    /*------------------------------------------------------------*/
    /** Constructor points[offset]...points[offset+length-1]; initializes the subproduct tree */
    zz_pX_Multipoint_General(const Vec<zz_p>& q, long offset, long length);

    /** @}  */       // doxygen constructors

    /*------------------------------------------------------------*/
    /* basic operations                                           */
    /*------------------------------------------------------------*/
    /** @name Basic operations.
     *  @{
     */ 
    void evaluate(Vec<zz_p>& val, const zz_pX& f) const;
    void interpolate(zz_pX& f, const Vec<zz_p>& val);
    void t_evaluate(zz_pX& f, const Vec<zz_p>& val, long output_size = -1) const;
    void t_interpolate(Vec<zz_p>& val, const zz_pX& f);
    /** @}  */       // doxygen basic

private:
    /** Used for interpolation: we premultiply given values v[i] by cofactors[i], then
     *  do linear combination of moduli
     */
    Vec<zz_p> cofactors;
    /** Subproduct tree. */
    Vec<Vec<zz_pX> > tree;
    /** Reverse root of subproduct tree */
    zz_pX reverse_root;
    /** Power series inverse of the reverse root of subproduct tree */
    zz_pX inverse_root;
};

/*------------------------------------------------------------*/
/* returns a zz_pX_Multipoint_General with n points           */
/*------------------------------------------------------------*/
/** Returns a zz_pX_Multipoint_General with n points.
 *  In the current implementation, points[i] = i, but this may change.
 */
zz_pX_Multipoint_General get_general_points(long n);


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* multipoint evaluation over zz_p                            */
/* points in geometric progression                            */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/**  Multipoint evaluation for points in geometric progression.  */
class zz_pX_Multipoint_Geometric : public zz_pX_Multipoint
{
public:
    /** @name Constructors.
     *  @{
     */ 

    /*------------------------------------------------------------*/
    /* boilerplate                                                */
    /*------------------------------------------------------------*/
    /** Empty constructor; no effect. */
    zz_pX_Multipoint_Geometric(){}

    /*------------------------------------------------------------*/
    /* constructor for geometric progressions                     */
    /* we interpolate at r^(2*i), i=0..d-1                        */
    /* polynomials to evaluate must have degree at most d-1       */
    /* r must have order at least 2*(d-1)+1                       */
    /*------------------------------------------------------------*/
    /** Constructor for geometric progressions.
     *  We evaluate / interpolate at r^(2*i), i=0..d-1.
     *  Polynomials to evaluate must have degree at most d-1.
     *  r must have order at least 2*(d-1)+1.
     */
    zz_pX_Multipoint_Geometric(const zz_p& r, long d);

    /*------------------------------------------------------------*/
    /* constructor for geometric progressions                     */
    /* we evaluate / interpolate at s * r^(2*i), i=0..d-1         */
    /* for evaluation, deg(f) must be < d                         */
    /* r must have order at least 2*(d-1)+1                       */
    /*------------------------------------------------------------*/
    /** Constructor for geometric progressions.
     *  We evaluate / interpolate at s*r^(2*i), i=0..d-1.
     *  Polynomials to evaluate must have degree at most d-1.
     *  r must have order at least 2*(d-1)+1.
     */
    zz_pX_Multipoint_Geometric(const zz_p& r, const zz_p& s, long d);

    /** @}  */       // doxygen constructors


    /*------------------------------------------------------------*/
    /* basic operations                                           */
    /*------------------------------------------------------------*/
    /** @name Basic operations.
     *  @{
     */ 

    void evaluate(Vec<zz_p>& val, const zz_pX& f) const;
    void interpolate(zz_pX& f, const Vec<zz_p>& val);
    void t_evaluate(zz_pX& f, const Vec<zz_p>& val, long output_size = -1) const;
    void t_interpolate(Vec<zz_p>& val, const zz_pX& f);

    /*------------------------------------------------------------*/
    /* evaluates f at s * r^(2*i), i=0..nb-1                      */
    /* still needs deg(f) < d                                     */ 
    /*------------------------------------------------------------*/
    /** Evaluates f at points[i], i=0..nb-1. */
    void evaluate(Vec<zz_p>& val, const zz_pX& f, long nb) const;

    /** @}  */       // doxygen basic

    /*------------------------------------------------------------*/
    /* vector version of basic operations                         */
    /*------------------------------------------------------------*/
    /** @name Vector version of basic operations.
     *  @{
     */ 
    /** A vector version of evaluation (input `f` is a vector of length `d`). */
    void mul_right(Vec<zz_p>& val, const Vec<zz_p>& f) const;
    /** A vector version of interpolation (input `val` is a vector of length `d`). */
    void inv_mul_right(Vec<zz_p>& f, const Vec<zz_p>& val);
    /** A vector version of transpose-evaluation (input `val` is a vector of length `d`). */
    void mul_left(Vec<zz_p>& f, const Vec<zz_p>& val) const;
    /** A vector version of transpose-interpolation (input `f` is a vector of length `d`). */
    void inv_mul_left(Vec<zz_p>& val, const Vec<zz_p>& f);

    /** @}  */       // doxygen vector

    /*------------------------------------------------------------*/
    /* getters / setters                                          */
    /*------------------------------------------------------------*/
    /** @name Getters and setters.
     *  @{
     */ 

    /** Flag that indicates if we use FFT-based implementation for evaluation.
     *  If not, fall back on Karatsuba (useful for low degrees).
     */
    long FFT_evaluate() const;
    /** Flag that indicates if we use FFT-based implementation for interpolation.
     *  If not, fall back on Karatsuba (useful for low degrees).
     */
    long FFT_interpolate() const;

    /** Sets FFT flag for evaluation. */
    void set_FFT_evaluate();

    /** Unsets FFT flag for evaluation. */
    void unset_FFT_evaluate();

    /** Sets FFT flag for interpolation. */
    void set_FFT_interpolate();

    /** Unsets FFT flag for interpolation. */
    void unset_FFT_interpolate();

    /*------------------------------------------------------------*/
    /* decides whether to use FFT or not                          */
    /*------------------------------------------------------------*/
    /** Automatically decide whether to use FFT-based implementation,
     *  based on precomputed threshold values.
     */
    void decide_FFT();

    /*------------------------------------------------------------*/
    /* adds a new FFT for repeated evaluations in degree d < n    */
    /* (not necessary, but provides a speed-up)                   */
    /*------------------------------------------------------------*/
    /** Adds new precomputed FFTs for repeated evaluations in degree d < n.
     *  Should be called once before these evaluations take place.
     *  Not necessary, but provides a speed-up.
     */
    void prepare_degree(long d);

    /*------------------------------------------------------------*/
    /* return the ratio q = r^2                                   */
    /*------------------------------------------------------------*/
    /** returns the ration q = r^2 */
    zz_p get_q() const; 

    /*------------------------------------------------------------*/
    /* return s (points are s*r^(2i))                             */
    /*------------------------------------------------------------*/
    /** returns the base s */
    zz_p get_s() const;

    /** @}  */       // doxygen get-set

private:  

    long idx_k, FFT_feasible, do_FFT_evaluate, do_FFT_interpolate;
    Vec<zz_p> x, xs, t, w, ws, y, z, zs;
    zz_pX f, g1, g2;
    fftRep g1_fft, g2_fft;  
    map<int, fftRep> known_degrees;
};

/*------------------------------------------------------------*/
/* returns a zz_pX_Multipoint_Geometric with n points         */
/*------------------------------------------------------------*/
/** Returns a zz_pX_Multipoint_Geometric with n points.
 *  In the current implementation, points[i] = r^(2i) for some random r, but this may change.
 */
zz_pX_Multipoint_Geometric get_geometric_points(long n);


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* multipoint evaluation over zz_p                            */
/* FFT points                                                 */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/**  Multipoint evaluation for FFT points. Length must be a power of 2.  */

class zz_pX_Multipoint_FFT : public zz_pX_Multipoint
{
public:
    /** @name Constructors.
     *  @{
     */ 

    /*------------------------------------------------------------*/
    /* boilerplate                                                */
    /*------------------------------------------------------------*/
    /** Empty constructor; no effect. */
    zz_pX_Multipoint_FFT(){}

    /*------------------------------------------------------------*/
    /* constructor for FFT points                                 */
    /* n has to be a power of 2, p must be FFT prime              */
    /*------------------------------------------------------------*/
    /** Constructor for FFT points. `n` must be a power of 2, `zz_p` must have been FFT-initialized. */
    zz_pX_Multipoint_FFT(long n);

    /** @}  */       // doxygen constructors

    /*------------------------------------------------------------*/
    /* basic operations                                           */
    /*------------------------------------------------------------*/
    /** @name Basic operations.
     *  @{
     */ 

    void evaluate(Vec<zz_p>& val, const zz_pX& f) const;
    void interpolate(zz_pX& f, const Vec<zz_p>& val);
    void t_evaluate(zz_pX& f, const Vec<zz_p>& val, long output_size = -1) const;
    void t_interpolate(Vec<zz_p>& val, const zz_pX& f);

    /** @}  */       // doxygen basic

private:
    long k, do_bit_reverse;
    fftRep wk; // used to store values for inverse FFT
    Vec<long> indices; // for bit-reversal, if needed
};

/*------------------------------------------------------------*/
/* returns a zz_pX_Multipoint_FFT with at least n points      */
/*------------------------------------------------------------*/
/** Returns a zz_pX_Multipoint_FFT with at least n points.
 */
zz_pX_Multipoint_FFT get_FFT_points(long n);


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* basic transforms over zz_p (naive, Karatsuba, Montgomery)  */
/* abstract class                                             */
/* usually not symmetric, so we have left / right transforms  */
/* work modulo any p                                          */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/**  An abstract class for matrix multiplication by bilinear algorithms
 *   (useful for very small degrees only.)
 *   The tensors we use may not be symmetric, so we have both left and right transforms.
 */

class zz_pX_Transform
{
public:

    /*------------------------------------------------------------*/
    /* basic operations                                           */
    /*------------------------------------------------------------*/
    /** @name Basic operations.
     *  @{
     */ 

    /** forward transform for left operand `A` */
    virtual void forward_left(Vec<zz_p>& val, const zz_pX& P) const = 0;
    /** forward transform for right operand `B` */
    virtual void forward_right(Vec<zz_p>& val, const zz_pX& P) const = 0;
    /** backward transform to recover result `C` */
    virtual void backward(zz_pX& P, const Vec<zz_p>& val) const = 0;  

    /** @}  */       // doxygen basic

    /*------------------------------------------------------------*/
    /* action on vectors and matrices                             */
    /*------------------------------------------------------------*/
    /** @name Actions on vectors and matrices.
     *  @{
     */ 
    /** forward transform for left vector operand `A` */
    void forward_left_vector(Vec<Vec<zz_p>>& val, const Vec<zz_pX>& P) const;
    /** forward transform for right vector operand `B` */
    void forward_right_vector(Vec<Vec<zz_p>>& val, const Vec<zz_pX>& P) const;
    /** backward transform to recover vector `C` */
    void backward_vector(Vec<zz_pX>& f, const Vec<Vec<zz_p>>& val) const;
    /** forward transform for left matrix operand `A` */
    void forward_left_matrix(Vec<Mat<zz_p>>& val, const Mat<zz_pX>& f) const;
    /** forward transform for right matrix operand `B` */
    void forward_right_matrix(Vec<Mat<zz_p>>& val, const Mat<zz_pX>& f) const;
    /** backward transform to recover matrix `C` */
    void backward_matrix(Mat<zz_pX>& f, const Vec<Mat<zz_p>>& val) const;

    /** @}  */       // doxygen vector

    /** @name Getters 
     *  @{
     */ 

    /*------------------------------------------------------------*/
    /* size of the input                                          */
    /*------------------------------------------------------------*/
    /** input length (same for left and right inputs) */
    long input_length() const;

    /*------------------------------------------------------------*/
    /* size of the transform                                      */
    /*------------------------------------------------------------*/
    /** length of the transform */
    long transform_length() const;

protected:
    long m, n; // m = input size, n = transform size
};


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* quadratic transform                                        */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/** A class that uses the naive (quadratic) multiplication algorithm.
 *  Works for any input length d, but dedicated classes should be uses
 *  for d = 2,3,4. */
class zz_pX_Transform_naive : public zz_pX_Transform
{
public:
    /** @name Constructors.
     *  @{
     */ 

    /*------------------------------------------------------------*/
    /* boilerplate                                                */
    /*------------------------------------------------------------*/
    /** Empty constructor; no effect. */
    zz_pX_Transform_naive(){}
    /** Empty destructor; no effect. */
    ~zz_pX_Transform_naive(){};

    /*------------------------------------------------------------*/
    /* constructor                                                */
    /* sets m and n                                               */
    /*------------------------------------------------------------*/
    /** Sets input length = d, transform length = d^2 */
    zz_pX_Transform_naive(long d);

    /** @}  */       // doxygen constructors


    /*------------------------------------------------------------*/
    /* basic operations                                           */
    /*------------------------------------------------------------*/
    /** @name Basic operations.
     *  @{
     */ 

    void forward_left(Vec<zz_p>& val, const zz_pX& P) const;
    void forward_right(Vec<zz_p>& val, const zz_pX& P) const;
    void backward(zz_pX& P, const Vec<zz_p>& val) const;  

    /** @}  */       // doxygen basic
};


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* karatsuba transform for 2-terms inputs                     */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/** A class that uses the Karatsuba multiplication algorithm.
 *  Input length = 2. */
class zz_pX_Transform_karatsuba : public zz_pX_Transform
{
public:
    /** @name Constructors.
     *  @{
     */ 

    /*------------------------------------------------------------*/
    /* boilerplate                                                */
    /*------------------------------------------------------------*/
    /** Empty destructor; no effect. */
    ~zz_pX_Transform_karatsuba(){};

    /*------------------------------------------------------------*/
    /* constructor                                                */
    /* sets m and n                                               */
    /*------------------------------------------------------------*/
    /** Sets input length = 2, transform length = 3*/
    zz_pX_Transform_karatsuba();

    /** @}  */       // doxygen constructors

    /*------------------------------------------------------------*/
    /* basic operations                                           */
    /*------------------------------------------------------------*/
    /** @name Basic operations.
     *  @{
     */ 

    void forward_left(Vec<zz_p>& val, const zz_pX& P) const;
    void forward_right(Vec<zz_p>& val, const zz_pX& P) const;
    void backward(zz_pX& P, const Vec<zz_p>& val) const;  

    /** @}  */       // doxygen basic
};

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* Mongtomery transform for 3-terms inputs                    */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/** A class that uses the Montgomery multiplication algorithm.
 *  Input length = 3. */
class zz_pX_Transform_montgomery3 : public zz_pX_Transform
{
public:
    /** @name Constructors.
     *  @{
     */ 

    /*------------------------------------------------------------*/
    /* boilerplate                                                */
    /*------------------------------------------------------------*/
    ~zz_pX_Transform_montgomery3(){};

    /*------------------------------------------------------------*/
    /* constructor                                                */
    /* sets m and n                                               */
    /*------------------------------------------------------------*/
    zz_pX_Transform_montgomery3();

    /** @}  */       // doxygen constructors

    /*------------------------------------------------------------*/
    /* basic operations                                           */
    /*------------------------------------------------------------*/
    /** @name Basic operations.
     *  @{
     */ 

    void forward_left(Vec<zz_p>& val, const zz_pX& P) const;
    inline void forward_right(Vec<zz_p>& val, const zz_pX& P) const
    {
        forward_left(val, P);
    }
    void backward(zz_pX& P, const Vec<zz_p>& val) const;  

    /** @}  */       // doxygen basic
};


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* karatsuba transform for 4-terms inputs                     */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
class zz_pX_Transform_karatsuba4 : public zz_pX_Transform
/** A class that uses the Karatsuba multiplication algorithm.
 *  Input length = 4. */
{
public:
    /** @name Constructors, destructors.
     *  @{
     */ 

    /*------------------------------------------------------------*/
    /* boilerplate                                                */
    /*------------------------------------------------------------*/
    ~zz_pX_Transform_karatsuba4(){};

    /*------------------------------------------------------------*/
    /* constructor                                                */
    /* sets m and n                                               */
    /*------------------------------------------------------------*/
    zz_pX_Transform_karatsuba4();

    /** @}  */       // doxygen constructors

    /*------------------------------------------------------------*/
    /* basic operations                                           */
    /*------------------------------------------------------------*/
    /** @name Basic operations.
     *  @{
     */ 

    void forward_left(Vec<zz_p>& val, const zz_pX& P) const;
    inline void forward_right(Vec<zz_p>& val, const zz_pX& P) const
    {
        forward_left(val, P);
    }
    void backward(zz_pX& P, const Vec<zz_p>& val) const;  

    /** @}  */       // doxygen basic
};





/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* Chinese Remaindering over zz_p                             */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/** A class for Chinese Remaindering. Pairwise coprime, squarefree moduli are supported. */
class zz_pX_CRT{
public:
    /** @name Constructors, destructors.
     *  @{
     */ 

    /*------------------------------------------------------------*/
    /* boilerplate                                                */
    /*------------------------------------------------------------*/
    /** Empty destructor; no effect. */
    ~zz_pX_CRT(){};

    /*------------------------------------------------------------*/
    /* constructor from moduli                                    */
    /*------------------------------------------------------------*/
    /** Builds the subproduct tree from moduli. */
    zz_pX_CRT(const Vec<zz_pX>& q);

    /** @}  */       // doxygen constructors
    

    /*------------------------------------------------------------*/
    /* basic operations                                           */
    /*------------------------------------------------------------*/
    /** @name Basic operations.
     *  @{
     */ 

    /** Returns val[i] = f mod q[i] for all i. Input f can have arbitrary degree. */
    void multimod(Vec<zz_pX>& val, const zz_pX& f) const;
    /** Recovers f of degree less than sum(deg(q[i])) such that f mod q[i] = val[i] for all i. */
    void combine(zz_pX& f, const Vec<zz_pX>& val) const;

    /** @}  */       // doxygen basic


    /** @name Getters 
     *  @{
     */ 

    /*------------------------------------------------------------*/
    /* master polynomial                                          */
    /*------------------------------------------------------------*/
    /** Returns the root of the subproduct tree. */
    void master(zz_pX& P) const;

    /*------------------------------------------------------------*/
    /* access to moduli                                           */
    /*------------------------------------------------------------*/
    /** Returns the i-th modulus q[i]. */
    void moduli(zz_pX& P, long i) const;

    /*------------------------------------------------------------*/
    /* number of moduli                                           */
    /*------------------------------------------------------------*/
    /** Returns the number of moduli. */
    long length() const;

private:
    long n;
    Vec<zz_pX> cofactors;
    Vec<Vec<zz_pX>> tree;
};

#endif

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
