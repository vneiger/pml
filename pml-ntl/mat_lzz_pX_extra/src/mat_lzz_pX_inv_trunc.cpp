#include <algorithm>

#include "util.h"
#include "mat_lzz_pX_utils.h"
#include "mat_lzz_pX_arith.h"
#include "mat_lzz_pX_inverse.h" // includes matrix multiply
#include "thresholds_matrix_multiply.h"
#include "thresholds_newton_inv_trunc.h"

NTL_CLIENT

/*------------------------------------------------------------*/
/* finds a sequence of degrees n0, n1, .. nk                  */
/* n0 <= thresh, ni = {2*n{i-1}, 2*{n-1}-2}, nk >= n          */
/*------------------------------------------------------------*/
static VecLong degrees(long n, long thresh)
{
    VecLong all_deg;

    while (n > thresh)
    {
        if (n % 2) // n is odd, add 1
            n++;
        all_deg.push_back(n);
        n >>= 1;
    }
    all_deg.push_back(n);
    std::reverse(all_deg.begin(), all_deg.end());
    return all_deg;
}

/*------------------------------------------------------------*/
/* returns x = 1/a mod z^m, naive algorithm                   */
/* throws an error if a(0) not invertible                     */
/* x can alias a                                              */
/*------------------------------------------------------------*/
void plain_inv_trunc(Mat<zz_pX>& x, const Mat<zz_pX>& a, long m)
{
    // no problem if &x == &a: anyway, we do computations in a temp `x_coeffs`
    // and assign the result to x at the very end

    // degree of a
    const long d = deg(a);
    if (d==-1)
        LogicError("Zero matrix as input of truncated inverse\n");

    // convert a to Vec<Mat<zz_p>> representation
    Vec<Mat<zz_p>> x_coeffs, a_coeffs;
    conv(a_coeffs, a, std::min(d+1,m));

    // compute the inverse of the constant coefficient
    x_coeffs.SetLength(m);
    inv(x_coeffs[0], a_coeffs[0]); // throws if a not square or a(0) not invertible

    // if a is constant, finished (we have the actual inverse)
    if (d==0)
    {
        conv(x, x_coeffs[0]);
        return;
    }

    // if a not constant, let's continue
    Mat<zz_p> buf;
    for (long k = 1; k < m; ++k)
    {
        // compute x_coeffs[k]
        // this is -(sum_{1 <= i <= k} x_coeffs[k-i] a_coeffs[i]) * x_coeffs[0]

        const long minkd = std::min(k,d);

        // first, compute x_coeffs[k] as the sum:
        // sum_{1 <= i <= k} x_coeffs[k-i] a_coeffs[i]
        mul(x_coeffs[k], x_coeffs[k-1], a_coeffs[1]);
        for (long i = 2; i <= minkd; ++i)
        {
            mul(buf, x_coeffs[k-i], a_coeffs[i]);
            add(x_coeffs[k], x_coeffs[k], buf);
        }

        // right multiply by x_coeffs[0]
        mul(x_coeffs[k], x_coeffs[k], x_coeffs[0]);

        // finally, negate x_coeffs[k]
        NTL::negate(x_coeffs[k], x_coeffs[k]);
    }

    // deduce output x
    conv(x, x_coeffs);
}


/*------------------------------------------------------------*/
/* returns x = 1/a mod z^m, Newton iteration                  */
/* throws an error if a(0) not invertible                     */
/* x can alias a                                              */
/* naive algorithm up to degree 2^thresh                      */
/* thresh=-1 is the default value, uses lookup table          */
/*------------------------------------------------------------*/
void newton_inv_trunc_FFT(Mat<zz_pX>& x, const Mat<zz_pX>& a, long m, long thresh)
{
    if (&x == &a)
    {
        Mat<zz_pX> y;
        newton_inv_trunc_FFT(y, a, m, thresh);
        x.swap(y);
        return;
    }

    if (thresh == -1)
        thresh = MATRIX_INV_TRUNC_PLAIN_THRESHOLD_FFT;

    long n, k, s, ss, t;
    Mat<zz_pX> P1;
    Vec<zz_p> mat_val1, mat_val2;
    Vec<Vec<zz_p>> mat_val3;
    Mat<zz_p> v1, v2, v3;

    fftRep R1(INIT_SIZE, NextPowerOfTwo(2*m-1));

    k = thresh;
    plain_inv_trunc(x, a, k); // throws if a not square or a(0) not invertible
    s = a.NumRows();
    ss = s * s;

    v1.SetDims(s, s);
    v2.SetDims(s, s);
    v3.SetDims(s, s);
    P1.SetDims(s, s);
    x.SetDims(s, s);

    while (k < m)
    {
        t = NextPowerOfTwo(2 * k);
        n = 1L << t;

        mat_val1.SetLength(n * ss);
        mat_val2.SetLength(n * ss);
        mat_val3.SetLength(ss);
        for (long i = 0; i < ss; i++)
            mat_val3[i].SetLength(n);

        // mat_val1 = FFT of x
        // deg(x) <= k-1
        for (long i = 0; i < s; i++)
            for (long j = 0; j < s; j++)
            {
                TofftRep(R1, x[i][j], t);
                for (long r = 0, rss = 0; r < n; r++, rss += ss)
                    mat_val1[rss + i*s + j] = R1.tbl[0][r];
            }

        // mat_val2 = FFT of (a mod t^(2k))
        // deg(a_ij mod t^(2k)) <= min(2k-1, deg(a_ij))
        for (long i = 0; i < s; i++)
            for (long j = 0; j < s; j++)
            {
                TofftRep(R1, a[i][j], t, 0, min(2*k - 1, deg(a[i][j])));
                for (long r = 0, rss = 0; r < n; r++, rss += ss)
                    mat_val2[rss + i*s + j] = R1.tbl[0][r];
            }

        // mat_val3 = values of ( x (a mod t^2k) ) = ___ + t^k delta
        for (long j = 0, jss = 0; j < n; j++, jss += ss)
        {
            for (long i = 0; i < s; i++)
                for (long ell = 0; ell < s; ell++)
                    v1[i][ell] = mat_val1[jss + i*s + ell];
            for (long i = 0; i < s; i++)
                for (long ell = 0; ell < s; ell++)
                    v2[i][ell] = mat_val2[jss + i*s + ell];
            mul(v3, v1, v2);
            for (long i = 0; i < s; i++)
                for (long ell = 0; ell < s; ell++)
                    mat_val3[i*s + ell][j] = v3[i][ell];
        }

        // mat_val2 = values of delta
        for (long i = 0; i < s; i++)
            for (long ell = 0; ell < s; ell++)
            {
                for (long r = 0; r < n; r++)
                    R1.tbl[0][r] = rep(mat_val3[i*s + ell][r]);

                zz_pX P;
                FromfftRep(P, R1, k, 2*k - 1);
                TofftRep(R1, P, t);
                for (long r = 0, rss = 0; r < n; r++, rss += ss)
                    mat_val2[rss + i*s + ell] = R1.tbl[0][r];
            }

        // mat_val3 = values of (delta x)
        for (long j = 0, jss = 0; j < n; j++, jss += ss)
        {
            for (long i = 0; i < s; i++)
                for (long ell = 0; ell < s; ell++)
                    v1[i][ell] = mat_val1[jss + i*s + ell];
            for (long i = 0; i < s; i++)
                for (long ell = 0; ell < s; ell++)
                    v2[i][ell] = mat_val2[jss + i*s + ell];
            mul(v3, v2, v1);
            for (long i = 0; i < s; i++)
                for (long ell = 0; ell < s; ell++)
                    mat_val3[i*s + ell][j] = v3[i][ell];
        }

        for (long i = 0; i < s; i++)
            for (long ell = 0; ell < s; ell++)
            {
                for (long r = 0; r < n; r++)
                    R1.tbl[0][r] = rep(mat_val3[i*s + ell][r]);

                zz_pX P;
                FromfftRep(P, R1, 0, k - 1);

                x[i][ell].rep.SetLength(2*k);
                long y_len = P.rep.length();
                for (long ii = k; ii < 2*k; ii++)
                {
                    if (ii-k >= y_len)
                        clear(x[i][ell].rep[ii]);
                    else
                        x[i][ell].rep[ii] = -P.rep[ii - k];
                }
                x[i][ell].normalize();
            }

        k = 2*k;
    }

    trunc(x, x, m);
}

/*------------------------------------------------------------*/
/* returns x = 1/a mod z^m, Newton iteration                  */
/* throws an error if a(0) not invertible                     */
/* x can alias a                                              */
/* naive algorithm up to degree 2^thresh                      */
/* thresh=-1 is the default value, uses lookup table          */
/*------------------------------------------------------------*/
void newton_inv_trunc_middle_product(Mat<zz_pX>& x, const Mat<zz_pX>& a, long m, long thresh)
{
    if (&x == &a)
    {
        Mat<zz_pX> y;
        newton_inv_trunc_middle_product(y, a, m);
        x.swap(y);
        return;
    }

    if (thresh == -1)
    {
        long t = type_of_prime();
        if (t == TYPE_SMALL_PRIME)
            thresh = MATRIX_INV_TRUNC_PLAIN_THRESHOLD_MIDDLE_SMALL;
        else
            thresh = MATRIX_INV_TRUNC_PLAIN_THRESHOLD_MIDDLE_LARGE;
    }

    VecLong all_deg = degrees(m, thresh);

    long k = all_deg[0];
    long idx = 1;

    plain_inv_trunc(x, a, k); // throws if a not square or a(0) not invertible

    Mat<zz_pX> y, tr;
    while (k < m)
    {
        trunc(tr, a, 2*k);
        middle_product(y, x, tr, k, k-1);
        mul_trunc(y, y, x, k);
        // next two could be slightly improved since this is essentiall copying
        // -y into the upper part of x, but this is a small fraction of the
        // total time (up to 2%)
        LeftShift(y, y, k);
        sub(x, x, y);
        k = 2*k;
        if (k != all_deg[idx])
        {
            k = all_deg[idx];
            trunc(x, x, k);
        }
        idx++;
    }

    trunc(x, x, m);
}

/*------------------------------------------------------------*/
/* returns x = 1/a mod z^m, Newton iteration                  */
/* throws an error if a(0) not invertible                     */
/* x can alias a                                              */
/* naive algorithm up to degree 2^thresh                      */
/* thresh=-1 is the default value, uses lookup table          */
/*------------------------------------------------------------*/
void newton_inv_trunc_geometric(Mat<zz_pX>& x, const Mat<zz_pX>& a, long m, long thresh)
{
    if (x == a)
    {
        Mat<zz_pX> y;
        newton_inv_trunc_geometric(y, a, m);
        x.swap(y);
        return;
    }

    long k, idx, s, ss;
    Mat<zz_p> v1, v2, v3;
    Vec<zz_p> mat_val1, mat_val2;
    Vec<Vec<zz_p>> mat_val3;

    if (thresh == -1)
    {
        long t = type_of_prime();
        if (t == TYPE_SMALL_PRIME)
            thresh = MATRIX_INV_TRUNC_PLAIN_THRESHOLD_GEOMETRIC_SMALL;
        else
            thresh = MATRIX_INV_TRUNC_PLAIN_THRESHOLD_GEOMETRIC_LARGE;
    }

    k = thresh;
    VecLong all_deg=degrees(m, k);
    k = all_deg[0];
    idx = 1;
    plain_inv_trunc(x, a, k); // throws if a not square or a(0) not invertible

    s = a.NumRows();
    ss = s * s;

    v1.SetDims(s, s);
    v2.SetDims(s, s);
    v3.SetDims(s, s);

    while (k < m)
    {
        long n = 3*k - 2;
        zz_pX_Multipoint_Geometric ev = get_geometric_points(n);

        ev.prepare_degree(k - 1);
        if (k != 1)
            ev.prepare_degree(2*k - 1);

        mat_val1.SetLength(n * ss);
        mat_val2.SetLength(n * ss);
        mat_val3.SetLength(ss);
        for (long i = 0; i < ss; i++)
            mat_val3[i].SetLength(n);

        Vec<zz_p> power_x;
        power_x.SetLength(n);
        const zz_p q = ev.get_q();
        zz_p inv_q = 1/q;

        power_x[0] = to_zz_p(1);
        for (long i = 1; i < n; i++)
            power_x[i] = power_x[i-1] * inv_q;

        Vec<zz_p> tmp;
        for (long i = 0; i < s; i++)
            for (long j = 0; j < s; j++)
            {
                ev.evaluate(tmp, x[i][j]); // degree = k - 1
                for (long r = 0, rss = 0; r < n; r++, rss += ss)
                    mat_val1[rss + i*s + j] = tmp[r];
            }


        for (long i = 0; i < s; i++)
            for (long j = 0; j < s; j++)
            {
                ev.evaluate(tmp, trunc(a[i][j], 2*k)); // degree = 2k - 1
                for (long r = 0, rss = 0; r < n; r++, rss += ss)
                    mat_val2[rss + i*s + j] = tmp[r];
            }

        // let y = x (a mod t^(2k)
        // then deg(y) = 3k - 2, and y = I + t^k delta
        // deg(delta) = 2k - 2
        // deg(delta x) = 3k - 3
        for (long j = 0, jss = 0; j < n; j++, jss += ss)
        {
            zz_p coeff = power(power_x[j], k);

            for (long i = 0; i < s; i++)
                for (long ell = 0; ell < s; ell++)
                    v1[i][ell] = mat_val1[jss + i*s + ell];
            for (long i = 0; i < s; i++)
                for (long ell = 0; ell < s; ell++)
                    v2[i][ell] = mat_val2[jss + i*s + ell];
            mul(v3, v1, v2);
            for (long i = 0; i < s; i++)
                for (long ell = 0; ell < s; ell++)
                    if (i == ell)
                        v3[i][ell] = (v3[i][ell] - 1) * coeff;
                    else
                        v3[i][ell] = v3[i][ell] * coeff;
            mul(v3, v3, v1);
            for (long i = 0; i < s; i++)
                for (long ell = 0; ell < s; ell++)
                    mat_val3[i*s + ell][j] = v3[i][ell];
        }


        for (long i = 0; i < s; i++)
            for (long ell = 0; ell < s; ell++)
            {
                zz_pX P;
                ev.interpolate(P, mat_val3[i*s + ell]);
                x[i][ell].rep.SetLength(2*k);
                long y_len = P.rep.length();
                for (long ii = k; ii < 2*k; ii++)
                {
                    if (ii-k >= y_len)
                        clear(x[i][ell].rep[ii]);
                    else
                        x[i][ell].rep[ii] = -P.rep[ii - k];
                }
                x[i][ell].normalize();
            }
        k = 2 * k;
        if (k != all_deg[idx])
        {
            k = all_deg[idx];
            trunc(x, x, k);
        }
        idx++;

    }
    trunc(x, x, m);
}


/*------------------------------------------------------------*/
/* returns x = 1/a mod z^m                                    */
/* throws an error if a(0) not invertible                     */
/* x can alias a                                              */
/*------------------------------------------------------------*/
void inv_trunc(Mat<zz_pX>& x, const Mat<zz_pX>& a, long m)
{
    if (is_FFT_ready(NextPowerOfTwo(m+2)))
        newton_inv_trunc_FFT(x, a, m);
    else
    {
        long m_middle = max_degree_middle_product_inv_trunc();
        if (m <= m_middle)
            newton_inv_trunc_middle_product(x, a, m);
        else
            newton_inv_trunc_geometric(x, a, m);
    }
}

/*------------------------------------------------------------*/
/* for i >= 0, define Si = coefficients of A^{-1} of degrees  */
/*             i-(2d-1) .. i-1, with d=deg(A)                 */
/* given src = Si, this computes S_{2i-d}                     */
/* invA = A^{-1} mod x^d                                      */
/* note: deg(Si) < 2d-1                                       */
/* output can alias input                                     */
/*------------------------------------------------------------*/
void high_order_lift_inverse_odd(
                                 Mat<zz_pX> & next,
                                 const Mat<zz_pX>& src,
                                 const std::unique_ptr<mat_lzz_pX_lmultiplier> & A,
                                 const std::unique_ptr<mat_lzz_pX_lmultiplier> & invA,
                                 long d
                                )
{
    Mat<zz_pX> b = A->multiply(trunc(src, d)); // deg-argument < d
    trunc(b, b, d); // deg(b) < d
    next = transpose(middle_product(transpose(b), transpose(src), d-1, d-1)); // deg(next) < d
    Mat<zz_pX> m = A->multiply(next); // deg(m) < 2d
    m >>= d;   // deg(m) < d
    trunc(m, m, d-1); // deg(m) < d-1
    Mat<zz_pX> tmp = invA->multiply(m); // deg(tmp) < 2d-2
    trunc(tmp, tmp, d-1);  // deg(tmp) < d-1
    tmp <<= d;  // deg(tmp) < 2d-1
    next = next - tmp; // deg(next) < 2d-1
}


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
