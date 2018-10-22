#ifndef LZZ_PX_CRT__H
#define LZZ_PX_CRT__H

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
void build_subproduct_tree(Vec<Vec<zz_pX>> & tree, const Vec<zz_pX> & q);

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* multipoint evaluation over zz_p                            */
/* abstract class                                             */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
class zz_pX_Multipoint
{
public:
    virtual ~zz_pX_Multipoint(){}

    /*------------------------------------------------------------*/
    /* basic operations                                           */
    /*------------------------------------------------------------*/
    virtual void evaluate(Vec<zz_p>& val, const zz_pX& P) const = 0;
    virtual void interpolate(zz_pX& P, const Vec<zz_p>& val) = 0;     // dirty hack: interpolate uses an instance fftRep as workspace

    virtual void t_evaluate(zz_pX& P, const Vec<zz_p>& val, long output_size = -1) const = 0;
    virtual void t_interpolate(Vec<zz_p>& val, const zz_pX& P) = 0;     // dirty hack: t_interpolate uses an instance fftRep as workspace

    /*------------------------------------------------------------*/
    /* action on vectors and matrices                             */
    /*------------------------------------------------------------*/
    void evaluate_vector(Vec<Vec<zz_p>>& val, const Vec<zz_pX>& P) const;
    void interpolate_vector(Vec<zz_pX>& f, const Vec<Vec<zz_p>>& val);
    void evaluate_matrix(Vec<Mat<zz_p>>& val, const Mat<zz_pX>& f) const;
    void interpolate_matrix(Mat<zz_pX>& f, const Vec<Mat<zz_p>>& val);

    void t_evaluate_vector(Vec<zz_pX>& f, const Vec<Vec<zz_p>>& val, long output_size = -1) const;
    void t_interpolate_vector(Vec<Vec<zz_p>>& val, const Vec<zz_pX>& P);
    void t_evaluate_matrix(Mat<zz_pX>& f, const Vec<Mat<zz_p>>& val, long output_size = -1) const;
    void t_interpolate_matrix(Vec<Mat<zz_p>>& val, const Mat<zz_pX>& f);

    /*------------------------------------------------------------*/
    /* number of points                                           */
    /*------------------------------------------------------------*/
    long length() const;

    /*------------------------------------------------------------*/
    /* get the i-th point                                         */
    /*------------------------------------------------------------*/
    void get_point(zz_p& pt, long i);  // sets the array pts if needed

    /*------------------------------------------------------------*/
    /* a naive conversion to a dense matrix                       */
    /*------------------------------------------------------------*/
    void to_dense(Mat<zz_p>& M);

    inline Mat<zz_p> to_dense()
    {
        Mat<zz_p> M;
        to_dense(M);
        return M;
    }

protected:
    Vec<zz_p> pts;
    long n;
};

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* multipoint evaluation over zz_p                            */
/* uses subproduct tree.                                      */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
class zz_pX_Multipoint_General : public zz_pX_Multipoint
{
public:
    /*------------------------------------------------------------*/
    /* boilerplate                                                */
    /*------------------------------------------------------------*/
    zz_pX_Multipoint_General(){}

    /*------------------------------------------------------------*/
    /* constructor from points                                    */
    /*------------------------------------------------------------*/
    zz_pX_Multipoint_General(const Vec<zz_p>& q);

    /*------------------------------------------------------------*/
    /* basic operations                                           */
    /*------------------------------------------------------------*/
    void evaluate(Vec<zz_p>& val, const zz_pX& f) const;
    void interpolate(zz_pX& f, const Vec<zz_p>& val);

    void t_evaluate(zz_pX& f, const Vec<zz_p>& val, long output_size = -1) const;
    void t_interpolate(Vec<zz_p>& val, const zz_pX& f);

private:
    Vec<zz_p> cofactors;
    Vec<Vec<zz_pX> > tree;
    zz_pX reverse_root, inverse_root;
};

/*------------------------------------------------------------*/
/* returns a zz_pX_Multipoint_General with n points           */
/*------------------------------------------------------------*/
zz_pX_Multipoint_General get_general_points(long n);


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* multipoint evaluation over zz_p                            */
/* points in geometric progression                            */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
class zz_pX_Multipoint_Geometric : public zz_pX_Multipoint
{
public:
    /*------------------------------------------------------------*/
    /* boilerplate                                                */
    /*------------------------------------------------------------*/
    zz_pX_Multipoint_Geometric(){}

    /*------------------------------------------------------------*/
    /* constructor for geometric progressions                     */
    /* we interpolate at r^(2*i), i=0..d-1                        */
    /* polynomials to evaluate must have degree at most d-1       */
    /*------------------------------------------------------------*/
    zz_pX_Multipoint_Geometric(const zz_p& r, long d);

    /*------------------------------------------------------------*/
    /* constructor for geometric progressions                     */
    /* we interpolate at s * r^(2*i), i=0..d-1                    */
    /*------------------------------------------------------------*/
    zz_pX_Multipoint_Geometric(const zz_p& r, const zz_p& s, long d);

    /*------------------------------------------------------------*/
    /* basic operations                                           */
    /*------------------------------------------------------------*/
    void evaluate(Vec<zz_p>& val, const zz_pX& f) const;
    void interpolate(zz_pX& f, const Vec<zz_p>& val);

    void t_evaluate(zz_pX& f, const Vec<zz_p>& val, long output_size = -1) const;
    void t_interpolate(Vec<zz_p>& val, const zz_pX& f);

    /*------------------------------------------------------------*/
    /* vector version of basic operations                         */
    /*------------------------------------------------------------*/
    void mul_right(Vec<zz_p>& val, const Vec<zz_p>& f) const;
    void inv_mul_right(Vec<zz_p>& f, const Vec<zz_p>& val);

    void mul_left(Vec<zz_p>& f, const Vec<zz_p>& val) const;
    void inv_mul_left(Vec<zz_p>& val, const Vec<zz_p>& f);

    /*------------------------------------------------------------*/
    /* getters / setters                                          */
    /*------------------------------------------------------------*/
    long FFT_evaluate() const;
    long FFT_interpolate() const;

    void set_FFT_evaluate();
    void unset_FFT_evaluate();
    void set_FFT_interpolate();
    void unset_FFT_interpolate();

    /*------------------------------------------------------------*/
    /* decides whether to use FFT or not                          */
    /*------------------------------------------------------------*/
    void decide_FFT();

    /*------------------------------------------------------------*/
    /* adds a new FFT for repeated evaluations in degree d < n    */
    /*------------------------------------------------------------*/
    void prepare_degree(long d);

    /*------------------------------------------------------------*/
    /* return the ratio q = r^2                                   */
    /*------------------------------------------------------------*/
    zz_p get_q();

    /*------------------------------------------------------------*/
    /* return s (points are s*r^(2i))                             */
    /*------------------------------------------------------------*/
    zz_p get_s();


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
zz_pX_Multipoint_Geometric get_geometric_points(long n);


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* multipoint evaluation over zz_p                            */
/* FFT points                                                 */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
class zz_pX_Multipoint_FFT : public zz_pX_Multipoint
{
public:
    /*------------------------------------------------------------*/
    /* boilerplate                                                */
    /*------------------------------------------------------------*/
    zz_pX_Multipoint_FFT(){}

    /*------------------------------------------------------------*/
    /* constructor for geometric progressions                     */
    /* n has to be a power of 2, p must be FFT prime              */
    /*------------------------------------------------------------*/
    zz_pX_Multipoint_FFT(long n);

    /*------------------------------------------------------------*/
    /* basic operations                                           */
    /*------------------------------------------------------------*/
    void evaluate(Vec<zz_p>& val, const zz_pX& f) const;
    void interpolate(zz_pX& f, const Vec<zz_p>& val);

    void t_evaluate(zz_pX& f, const Vec<zz_p>& val, long output_size = -1) const;
    void t_interpolate(Vec<zz_p>& val, const zz_pX& f);


private:
    long k, do_bit_reverse;
    fftRep wk; // used to store values for inverse FFT
    Vec<long> indices; // for bit-reversal, if needed
};

/*------------------------------------------------------------*/
/* returns a zz_pX_Multipoint_FFT with at least n points      */
/*------------------------------------------------------------*/
zz_pX_Multipoint_FFT get_FFT_points(long n);



/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* basic transforms over zz_p (naive, Karatsuba, Montgomery)  */
/* abstract class                                             */
/* usually not symmetric, so we have left / right transforms  */
/* work modulo any p                                          */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
class zz_pX_Transform
{
public:

    /*------------------------------------------------------------*/
    /* basic operations                                           */
    /*------------------------------------------------------------*/
    virtual void forward_left(Vec<zz_p>& val, const zz_pX& P) const = 0;
    virtual void forward_right(Vec<zz_p>& val, const zz_pX& P) const = 0;
    virtual void backward(zz_pX& P, const Vec<zz_p>& val) const = 0;  

    /*------------------------------------------------------------*/
    /* action on vectors and matrices                             */
    /*------------------------------------------------------------*/
    void forward_left_vector(Vec<Vec<zz_p>>& val, const Vec<zz_pX>& P) const;
    void forward_right_vector(Vec<Vec<zz_p>>& val, const Vec<zz_pX>& P) const;
    void backward_vector(Vec<zz_pX>& f, const Vec<Vec<zz_p>>& val) const;
    void forward_left_matrix(Vec<Mat<zz_p>>& val, const Mat<zz_pX>& f) const;
    void forward_right_matrix(Vec<Mat<zz_p>>& val, const Mat<zz_pX>& f) const;
    void backward_matrix(Mat<zz_pX>& f, const Vec<Mat<zz_p>>& val) const;

    /*------------------------------------------------------------*/
    /* size of the input                                          */
    /*------------------------------------------------------------*/
    long input_length() const;

    /*------------------------------------------------------------*/
    /* size of the transform                                      */
    /*------------------------------------------------------------*/
    long transform_length() const;

protected:
    long m, n; // m = input size, n = transform size
};

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* quadratic transform                                        */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
class zz_pX_Transform_naive : public zz_pX_Transform
{
public:
    /*------------------------------------------------------------*/
    /* boilerplate                                                */
    /*------------------------------------------------------------*/
    zz_pX_Transform_naive(){}
    ~zz_pX_Transform_naive(){};

    /*------------------------------------------------------------*/
    /* constructor                                                */
    /* sets m and n                                               */
    /*------------------------------------------------------------*/
    zz_pX_Transform_naive(long d);

    /*------------------------------------------------------------*/
    /* basic operations                                           */
    /*------------------------------------------------------------*/
    void forward_left(Vec<zz_p>& val, const zz_pX& P) const;
    void forward_right(Vec<zz_p>& val, const zz_pX& P) const;
    void backward(zz_pX& P, const Vec<zz_p>& val) const;  
};


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* quadratic transform                                        */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
class zz_pX_Transform_karatsuba : public zz_pX_Transform
{
public:
    /*------------------------------------------------------------*/
    /* boilerplate                                                */
    /*------------------------------------------------------------*/
    ~zz_pX_Transform_karatsuba(){};

    /*------------------------------------------------------------*/
    /* constructor                                                */
    /* sets m and n                                               */
    /*------------------------------------------------------------*/
    zz_pX_Transform_karatsuba();

    /*------------------------------------------------------------*/
    /* basic operations                                           */
    /*------------------------------------------------------------*/
    void forward_left(Vec<zz_p>& val, const zz_pX& P) const;
    void forward_right(Vec<zz_p>& val, const zz_pX& P) const;
    void backward(zz_pX& P, const Vec<zz_p>& val) const;  
};

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* 3-terms input transform                                     */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
class zz_pX_Transform_montgomery3 : public zz_pX_Transform
{
public:
    /*------------------------------------------------------------*/
    /* boilerplate                                                */
    /*------------------------------------------------------------*/
    ~zz_pX_Transform_montgomery3(){};

    /*------------------------------------------------------------*/
    /* constructor                                                */
    /* sets m and n                                               */
    /*------------------------------------------------------------*/
    zz_pX_Transform_montgomery3();

    /*------------------------------------------------------------*/
    /* basic operations                                           */
    /*------------------------------------------------------------*/
    void forward_left(Vec<zz_p>& val, const zz_pX& P) const;
    inline void forward_right(Vec<zz_p>& val, const zz_pX& P) const
    {
        forward_left(val, P);
    }
    void backward(zz_pX& P, const Vec<zz_p>& val) const;  
};


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* 4-terms input transform                                     */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
class zz_pX_Transform_karatsuba4 : public zz_pX_Transform
{
public:
    /*------------------------------------------------------------*/
    /* boilerplate                                                */
    /*------------------------------------------------------------*/
    ~zz_pX_Transform_karatsuba4(){};

    /*------------------------------------------------------------*/
    /* constructor                                                */
    /* sets m and n                                               */
    /*------------------------------------------------------------*/
    zz_pX_Transform_karatsuba4();

    /*------------------------------------------------------------*/
    /* basic operations                                           */
    /*------------------------------------------------------------*/
    void forward_left(Vec<zz_p>& val, const zz_pX& P) const;
    inline void forward_right(Vec<zz_p>& val, const zz_pX& P) const
    {
        forward_left(val, P);
    }
    void backward(zz_pX& P, const Vec<zz_p>& val) const;  
};











/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* Chinese Remaindering over zz_p                             */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
class zz_pX_CRT{
public:
    /*------------------------------------------------------------*/
    /* boilerplate                                                */
    /*------------------------------------------------------------*/
    ~zz_pX_CRT(){};

    /*------------------------------------------------------------*/
    /* constructor from moduli                                    */
    /*------------------------------------------------------------*/
    zz_pX_CRT(const Vec<zz_pX>& q);

    /*------------------------------------------------------------*/
    /* basic operations                                           */
    /*------------------------------------------------------------*/
    void multimod(Vec<zz_pX>& val, const zz_pX& f) const;
    void combine(zz_pX& f, const Vec<zz_pX>& val) const;

    /*------------------------------------------------------------*/
    /* master polynomial                                          */
    /*------------------------------------------------------------*/
    void master(zz_pX& P) const;

    /*------------------------------------------------------------*/
    /* access to moduli                                           */
    /*------------------------------------------------------------*/
    void moduli(zz_pX& P, long i) const;

    /*------------------------------------------------------------*/
    /* number of moduli                                           */
    /*------------------------------------------------------------*/
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
