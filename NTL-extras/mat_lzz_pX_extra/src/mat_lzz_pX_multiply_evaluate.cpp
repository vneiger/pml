#include <NTL/matrix.h>
#include <NTL/mat_lzz_p.h>
#include <NTL/lzz_pX.h>
#include <NTL/BasicThreadPool.h>

#include "util.h"
#include "lzz_p_extra.h"
#include "mat_lzz_pX_extra.h"
#include "lzz_pX_CRT.h"

NTL_CLIENT


/*------------------------------------------------------------*/
/* c = a*b                                                    */
/* output may alias input; c does not have to be zero matrix  */
/*------------------------------------------------------------*/
static void multiply_evaluate_do_it(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b, zz_pX_Multipoint& ev)
{
    zz_pContext context;
    long s = a.NumRows();
    long t = a.NumCols();
    long u = b.NumCols();

    long n = ev.length();

    Vec<zz_p> mat_valA, mat_valB;
    Vec<Vec<zz_p>> mat_valC;

    mat_valA.SetLength(n * s * t);
    mat_valB.SetLength(n * t * u);
    long st = s*t;

    
/*** START PARALLEL ********************************************/ 
    context.save();  

NTL_EXEC_RANGE(s,first,last)
    
    context.restore();
    for (long i = first; i < last; i++)
    {
        Vec<zz_p> tmp;
        for (long k = 0; k < t; k++)
        {
            ev.evaluate(tmp, a[i][k]);
            for (long r = 0, rst = 0; r < n; r++, rst += st)
                mat_valA[rst + i*t + k] = tmp[r];
        }
    }

NTL_EXEC_RANGE_END
/*** END PARALLEL **********************************************/
    
    long tu = t*u;
    
/*** START PARALLEL ********************************************/
    context.save();
    
NTL_EXEC_RANGE(t,first,last)
    
    context.restore();    
    for (long i = first; i < last; i++)
    {
        Vec<zz_p> tmp;
        for (long k = 0; k < u; k++)
        {
            ev.evaluate(tmp, b[i][k]);
            for (long r = 0, rtu = 0; r < n; r++, rtu += tu)
            {
                mat_valB[rtu + i*u + k] = tmp[r];
            }
        }
        
    }
    
NTL_EXEC_RANGE_END
/*** END PARALLEL **********************************************/
    
    mat_valC.SetLength(s * u);
    for (long i = 0; i < s * u; i++)
        mat_valC[i].SetLength(n);
    
/*** START PARALLEL ********************************************/   
//    context.save();

//NTL_EXEC_RANGE(n,first,last)    

//    context.restore();
    for (long j = 0, jst = 0, jtu = 0; j < n; j++, jst += st, jtu += tu)
    {
        Mat<zz_p> va, vb, vc;
        va.SetDims(s, t);
        vb.SetDims(t, u);
        
        for (long i = 0; i < s; i++)
        {
            for (long k = 0; k < t; k++)
            {
                va[i][k] = mat_valA[jst + i*t + k];
            }
        }
        for (long i = 0; i < t; i++)
        {
            for (long k = 0; k < u; k++)
            {
                vb[i][k] = mat_valB[jtu + i*u + k];
            }
        }
        
        vc = va * vb;

        for (long i = 0; i < s; i++)
        {
            for (long k = 0; k < u; k++)
            {
                mat_valC[i*u + k][j] = vc[i][k];
            }
        }
    }
    
//NTL_EXEC_RANGE_END
/*** END PARALLEL **********************************************/

    c.SetDims(s, u);

/*** START PARALLEL ********************************************/   
    context.save();

NTL_EXEC_RANGE(s,first,last)   
    
    context.restore();
    
    for (long i = first; i < last; i++)
    {
        for (long k = 0; k < u; k++)
        {
            ev.interpolate(c[i][k], mat_valC[i*u + k]);
        }
    }

NTL_EXEC_RANGE_END
/*** END PARALLEL **********************************************/
}
/*------------------------------------------------------------*/
/* c = a*b                                                    */
/* geometric points, uses FFTs directly                       */
/* output may alias input; c does not have to be zero matrix  */
/* (exists mostly for testing purposes)                       */
/*------------------------------------------------------------*/
void multiply_evaluate_geometric_using_FFT(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b)
{
    long dA = deg(a);
    long dB = deg(b);
    long dC = dA+dB;
    long sz = dC+1;

    zz_pX_Multipoint *ev;
    zz_pX_Multipoint_Geometric ev_geometric = get_geometric_points(sz);

    ev_geometric.prepare_degree(dA);
    ev_geometric.prepare_degree(dB);

    ev_geometric.set_FFT_evaluate();
    ev_geometric.set_FFT_interpolate();

    ev = &ev_geometric;
    multiply_evaluate_do_it(c, a, b, *ev);
}

/*------------------------------------------------------------*/
/* c = a*b                                                    */
/* geometric points, uses middle product                      */
/* output may alias input; c does not have to be zero matrix  */
/* (exists mostly for testing purposes)                       */
/*------------------------------------------------------------*/
void multiply_evaluate_geometric_no_FFT(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b)
{
    long dA = deg(a);
    long dB = deg(b);
    long dC = dA+dB;
    long sz = dC+1;

    zz_pX_Multipoint *ev;
    zz_pX_Multipoint_Geometric ev_geometric = get_geometric_points(sz);

    ev_geometric.prepare_degree(dA);
    ev_geometric.prepare_degree(dB);

    ev_geometric.unset_FFT_evaluate();
    ev_geometric.unset_FFT_interpolate();

    ev = &ev_geometric;
    multiply_evaluate_do_it(c, a, b, *ev);
}

/*------------------------------------------------------------*/
/* c = a*b                                                    */
/* geometric points, thresholded                              */
/* output may alias input; c does not have to be zero matrix  */
/*------------------------------------------------------------*/
void multiply_evaluate_geometric(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b)
{
    long dA = deg(a);
    long dB = deg(b);
    long dC = dA+dB;
    long sz = dC+1;

    zz_pX_Multipoint *ev;
    zz_pX_Multipoint_Geometric ev_geometric = get_geometric_points(sz);

    ev_geometric.prepare_degree(dA);
    ev_geometric.prepare_degree(dB);

    ev = &ev_geometric;
    multiply_evaluate_do_it(c, a, b, *ev);
}



// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
