#include "util.h"
#include "mat_lzz_pX_multiply.h"

//#define PROFILE_ON

NTL_CLIENT

/*------------------------------------------------------------*/
/* computes the matrices for evaluation and interpolation     */
/*------------------------------------------------------------*/
void vandermonde(Mat<zz_p>& vdm1, Mat<zz_p>& vdm2, Mat<zz_p>& inv_vdm, long d1, long d2)
{
    // sizes = degree + 1
    const long s1 = d1 + 1;
    const long s2 = d2 + 1;
    // nb points, which is sum of degrees + 1
    const long nb_points = d1 + d2 + 1;

    // will be Vandermonde matrix with nb_points, chosen as 0, 1, .., nb_points-1
    // row i contains 1, i, i^2, .., i^{sk-1} for k=1,2
    vdm1.SetDims(nb_points, s1);
    vdm2.SetDims(nb_points, s2);

    // vdm: square Vandermonde matrix with nb_points
    // points chosen as 0, 1, .., nb_points-1
    // --> row i contains 1, i, i^2, .., i^{nb_points-1}
    // vdm1: nb_points x s1 submatrix of vdm
    // vdm2: nb_points x s2 submatrix of vdm
    Mat<zz_p> vdm;
    vdm.SetDims(nb_points, nb_points);
    vdm[0][0].LoopHole() = 1;
    vdm1[0][0].LoopHole() = 1;
    vdm2[0][0].LoopHole() = 1;
    for (long i = 1; i < nb_points; ++i)
    {
        const zz_p p1(i);
        vdm[i][0].LoopHole() = 1;
        vdm1[i][0].LoopHole() = 1;
        vdm2[i][0].LoopHole() = 1;
        for (long j = 1; j < nb_points; ++j)
        {
            mul(vdm[i][j], vdm[i][j-1], p1);
            if (j<s1)
                vdm1[i][j] = vdm[i][j];
            if (j<s2)
                vdm2[i][j] = vdm[i][j];
        }
    }

    // inv_vdm is the inverse of vdm
    inv(inv_vdm, vdm);
}

/*------------------------------------------------------------*/
/* computes the matrices for evaluation and interpolation     */
/*------------------------------------------------------------*/
void vandermonde2(
                  Mat<zz_p> & vdm1,
                  Mat<zz_p> & vdm2,
                  Mat<zz_p> & inv_vdm,
                  long d1,
                  long d2
                 )
{
    // sizes (for even part)
    const long s1 = d1/2+1;
    const long s2 = d2/2+1;
    // nb points, such that 2*nb_points >= d1+d2+1
    const long nb_points = (d1 + d2) / 2 + 1;

    // vdm: square Vandermonde matrix with nb_points
    // points chosen as the squares of 1, 2, .., nb_points
    // --> row i contains 1, (i+1)^2, (i+1)^4, .., (i+1)^{2*nb_points-2}
    // vdm1: nb_points x s1 submatrix of vdm
    // vdm2: nb_points x s2 submatrix of vdm
    vdm1.SetDims(nb_points, s1);
    vdm2.SetDims(nb_points, s2);
    Mat<zz_p> vdm(INIT_SIZE, nb_points, nb_points);
    for (long i = 0; i < nb_points; ++i)
    {
        const zz_p pt = sqr(to_zz_p(i+1));
        vdm[i][0].LoopHole() = 1;
        vdm1[i][0].LoopHole() = 1;
        vdm2[i][0].LoopHole() = 1;
        for (long j = 1; j < nb_points; ++j)
        {
            mul(vdm[i][j], vdm[i][j-1], pt);
            if (j<s1)
                vdm1[i][j] = vdm[i][j];
            if (j<s2)
                vdm2[i][j] = vdm[i][j];
        }
    }

    // inv_vdm is the inverse of vdm
    inv(inv_vdm, vdm);
}

/*------------------------------------------------------------*/
/* c = a*b                                                    */
/* output may alias input; c does not have to be zero matrix  */
/* matrix multiplication using the algorithm of Doliskani et al. */
/* uses matrix multiplication for evaluation and interpolation*/
/*------------------------------------------------------------*/
void multiply_evaluate_dense(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b)
{
    const long dA = deg(a);
    const long dB = deg(b);
    const long min_dAdB = std::min<long>(dA,dB);
    const long m = a.NumRows();
    const long n = a.NumCols();
    const long p = b.NumCols();
    long ell;

#ifdef PROFILE_ON
    double t = GetWallTime();
#endif // PROFILE_ON

    Mat<zz_p> vA, vB, iv;
    vandermonde(vA, vB, iv, dA, dB);
    const long nb_points = vA.NumRows();

#ifdef PROFILE_ON
    t = GetWallTime()-t;
    std::cout << std::endl << "Build vdmd mat: " << t << std::endl;
    
    t = GetWallTime();
#endif // PROFILE_ON
    
#ifdef PROFILE_ON
    double t_copy, t_mul;
    t = GetWallTime();
#endif // PROFILE_ON
    
    // evaluation of matrix a:
    // build tmp_mat, whose column ell = i*n + j is the coefficient vector of
    // a[i][j] (padded with zeroes up to length dA+1 if necessary)
    Mat<zz_p> tmp_mat(INIT_SIZE, dA+1, m * n);
    ell = 0;
    for (long i = 0; i < m; ++i)
        for (long j = 0; j < n; ++j, ++ell)
            for (long k = 0; k <= deg(a[i][j]); ++k)
                tmp_mat[k][ell] = a[i][j][k];
    // note: d = deg(a[i][j]) is -1 if a[i][j] == 0
    // all non-touched entries already zero since tmp_mat was initialized as zero

#ifdef PROFILE_ON
    t_copy = GetWallTime() - t;
    t = GetWallTime();
#endif // PROFILE_ON

    // valA: column ell = i*n + j contains the evaluations of a[i][j]
    Mat<zz_p> valA;
    mul(valA, vA, tmp_mat);

#ifdef PROFILE_ON
    t_mul = GetWallTime() - t;
    t = GetWallTime();
#endif // PROFILE_ON

    // evaluation of matrix b:
    // build tmp_mat, whose column ell = i*n + j is the coefficient vector of
    // b[i][j] (padded with zeroes up to length dB+1 if necessary)
    tmp_mat.SetDims(dB+1, n * p);
    ell = 0;
    for (long i = 0; i < n; ++i)
        for (long j = 0; j < p; ++j, ++ell)
        {
            long d = deg(b[i][j]); // -1 if b[i][j] == 0
            for (long k = 0; k <= d; ++k)
                tmp_mat[k][ell] = b[i][j][k];
            // make sure remaining entries are zero
            //    those for k<=dB, k>dA are (if any),
            //    those for d < k <= min(dA,dB) might not be
            for (long k = d+1; k <= min_dAdB; ++k)
                clear(tmp_mat[k][ell]);
        }

#ifdef PROFILE_ON
    t_copy += GetWallTime() - t;
    t = GetWallTime();
#endif // PROFILE_ON

    // valB: column ell = i*n + j contains the evaluations of b[i][j]
    Mat<zz_p> valB;
    mul(valB, vB, tmp_mat);

#ifdef PROFILE_ON
    t_mul += GetWallTime() - t;
    std::cout << "Evaluating a and b, data movements: " << t_copy << std::endl;
    std::cout << "Evaluating a and b, multiplication: " << t_mul << std::endl;
    t = GetWallTime();
#endif // PROFILE_ON

    // perform the pointwise products
    Mat<zz_p> valAp(INIT_SIZE, m, n);
    Mat<zz_p> valBp(INIT_SIZE, n, p);
    Mat<zz_p> valCp;
    Mat<zz_p> valC(INIT_SIZE, nb_points, m * p);
    for (long i = 0; i < nb_points; ++i)
    {
        // a evaluated at point i
        ell = 0;
        for (long u = 0; u < m; ++u)
            for (long v = 0; v < n; ++v, ++ell)
                valAp[u][v] = valA[i][ell];

        // b evaluated at point i
        ell = 0;
        for (long u = 0; u < n; ++u)
            for (long v = 0; v < p; ++v, ++ell)
                valBp[u][v] = valB[i][ell];

        // a*b evaluated at point i
        mul(valCp, valAp, valBp);

        // copy this into valC: column ell = i*n + j contains the evaluations
        // of the entry i,j of c = a*b
        ell = 0;
        for (long u = 0; u < m; ++u)
            for (long v = 0; v < p; ++v, ++ell)
                valC[i][ell] = valCp[u][v];
    }
#ifdef PROFILE_ON
    t = GetWallTime()-t;
    std::cout << "Pointwise prods: " << t << std::endl;

    t = GetWallTime();
#endif // PROFILE_ON

    // interpolate to find the entries of c
    mul(tmp_mat, iv, valC);

    // copy to output (reorganize these entries into c)
    c.SetDims(m, p);
    ell = 0;
    for (long u = 0; u < m; ++u)
        for (long v = 0; v < p; ++v, ++ell)
        {
            c[u][v].SetLength(nb_points);
            for (long i = 0; i < nb_points; ++i)
                c[u][v][i] = tmp_mat[i][ell];
            c[u][v].normalize();
        }
#ifdef PROFILE_ON
    t = GetWallTime()-t;
    std::cout << "Interpolate: " << t << std::endl;
#endif // PROFILE_ON

}


/*------------------------------------------------------------*/
/* c = a*b                                                    */
/* output may alias input; c does not have to be zero matrix  */
/* matrix multiplication using the algorithm of Doliskani et al. */
/* uses matrix multiplication for evaluation and interpolation*/
/*------------------------------------------------------------*/
void multiply_evaluate_dense2(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b)
{
    // dimensions and degrees
    const long dA = deg(a);
    const long dB = deg(b);
    const long sA = dA/2+1;
    const long sB = dB/2+1;
    const long m = a.NumRows();
    const long n = a.NumCols();
    const long p = b.NumCols();

    long ell;

    Mat<zz_p> vA, vB, iv;

#ifdef PROFILE_ON
    double t = GetWallTime();
#endif // PROFILE_ON
    vandermonde2(vA, vB, iv, dA, dB);
#ifdef PROFILE_ON
    t = GetWallTime()-t;
    std::cout << std::endl << "Build vdmd mat: " << t << std::endl;
#endif // PROFILE_ON

    const long nb_points = vA.NumRows();

#ifdef PROFILE_ON
    double t_copy, t_mul;
    t = GetWallTime();
#endif // PROFILE_ON

    // evaluation of matrix a (even/odd part):
    // build tmp_mat_even/odd, whose column ell = i*n + j is the coefficient
    // vector of the even/odd part of a[i][j] (padded with zeroes up to length
    // dA/2+1 if necessary)
    Mat<zz_p> tmp_mat_even(INIT_SIZE, sA, m*n);
    Mat<zz_p> tmp_mat_odd(INIT_SIZE, sA, m*n);
    // --> we use dimension sA=dA/2+1 for odd part, even though this should be
    // (dA-1)/2+1; this is slightly non-optimal...
    ell = 0;
    for (long i = 0; i < m; ++i)
        for (long j = 0; j < n; ++j, ++ell)
        {
            const long d = deg(a[i][j]);
            for (long k = 0; 2*k < d; ++k)
            {
                tmp_mat_even[k][ell] = a[i][j][2*k];
                tmp_mat_odd[k][ell] = a[i][j][2*k+1];
            }
            if (d % 2 == 0) // d even
                tmp_mat_even[d/2][ell] = a[i][j][d];
        }

#ifdef PROFILE_ON
    t_copy = GetWallTime()-t;
    t = GetWallTime();
#endif // PROFILE_ON

    // note: d = deg(a[i][j]) is -1 if a[i][j] == 0
    // all non-touched entries already zero since tmp_mat was initialized as zero

    // valAeven: column ell = i*n + j contains the evaluations of even part of a[i][j]
    // valAodd: column ell = i*n + j contains the evaluations of odd part of a[i][j]
    Mat<zz_p> valAeven, valAodd;
    mul(valAeven, vA, tmp_mat_even);
    mul(valAodd, vA, tmp_mat_odd);
#ifdef PROFILE_ON
    t_mul = GetWallTime()-t;
    t = GetWallTime();
#endif // PROFILE_ON

    // evaluation of matrix b (even/odd part): build tmp_mat_even/odd, whose
    // column ell = i*n + j is the coefficient vector of the even/odd part of
    // b[i][j] (padded with zeroes up to length dB/2+1 if necessary); as above,
    // we use dimension dB/2+1 although we know it should be (dB-1)/2+1)

    tmp_mat_even.SetDims(sB, n * p);
    tmp_mat_odd.SetDims(sB, n * p);
    ell = 0;
    for (long i = 0; i < n; ++i)
        for (long j = 0; j < p; ++j, ++ell)
        {
            const long d = deg(b[i][j]); // -1 if b[i][j] == 0
            long k = 0;
            for (; 2*k < d; ++k)
            {
                tmp_mat_even[k][ell] = b[i][j][2*k];
                tmp_mat_odd[k][ell] = b[i][j][2*k+1];
            }
            if (2*k == d)
            {
                tmp_mat_even[k][ell] = b[i][j][d];
                clear(tmp_mat_odd[k][ell]);
                ++k;
            }
            for (; k < sB; ++k)
            {
                clear(tmp_mat_even[k][ell]);
                clear(tmp_mat_odd[k][ell]);
            }
        }
#ifdef PROFILE_ON
    t_copy += GetWallTime()-t;
    t = GetWallTime();
#endif // PROFILE_ON

    // valBeven/odd: column ell = i*n + j contains the evaluations of the
    // even/odd part of b[i][j]
    Mat<zz_p> valBeven, valBodd;
    mul(valBeven, vB, tmp_mat_even);
    mul(valBodd, vB, tmp_mat_odd);

#ifdef PROFILE_ON
    t_mul += GetWallTime()-t;
    std::cout << "Evals of a and b(data movements): " << t_copy << std::endl;
    std::cout << "Evals of a and b(mul): " << t_mul << std::endl;
    
    t = GetWallTime();
#endif // PROFILE_ON

    // perform the pointwise products for the points 1, 2, ..., nb_points
    Mat<zz_p> valApeven(INIT_SIZE, m, n);
    Mat<zz_p> valApodd(INIT_SIZE, m, n);
    Mat<zz_p> valBpeven(INIT_SIZE, n, p);
    Mat<zz_p> valBpodd(INIT_SIZE, n, p);
    Mat<zz_p> valCeven(INIT_SIZE, nb_points, m * p);
    Mat<zz_p> valCodd(INIT_SIZE, nb_points, m * p);
    Mat<zz_p> valCpeven, valCpodd;
    
    // from now on we will only use valAodd with row i multiplied by i+1
    for (long i = 0; i < nb_points; ++i)
    {
        zz_p pt(i+1);
        mul(valAodd[i], pt, valAodd[i]);
        mul(valBodd[i], pt, valBodd[i]);
    }

    zz_p inv2; inv(inv2, to_zz_p(2));
    for (long i = 0; i < nb_points; ++i)
    {
        // a evaluated at point i (which is i+1)
        ell = 0;
        for (long u = 0; u < m; ++u)
            for (long v = 0; v < n; ++v, ++ell)
            {
                add(valApeven[u][v], valAeven[i][ell], valAodd[i][ell]);
                sub(valApodd[u][v], valAeven[i][ell], valAodd[i][ell]);
            }

        // b evaluated at point i (which is i+1)
        ell = 0;
        for (long u = 0; u < n; ++u)
            for (long v = 0; v < p; ++v, ++ell)
            {
                add(valBpeven[u][v], valBeven[i][ell], valBodd[i][ell]);
                sub(valBpodd[u][v], valBeven[i][ell], valBodd[i][ell]);
            }

        // a*b evaluated at point i (which is i+1)
        mul(valCpeven, valApeven, valBpeven);
        mul(valCpodd, valApodd, valBpodd);

        // copy this into valC: column ell = i*n + j contains the evaluations
        // of the entry i,j of c = a*b
        // valCeven = (even + odd)/2
        // valCodd = (even - odd)/2x
        ell = 0;
        zz_p inv2pt; inv(inv2pt, to_zz_p(2*(i+1)));
        for (long u = 0; u < m; ++u)
            for (long v = 0; v < p; ++v, ++ell)
            {
                add(valCeven[i][ell], valCpeven[u][v], valCpodd[u][v]);
                sub(valCodd[i][ell], valCpeven[u][v], valCpodd[u][v]);
            }
        mul(valCeven[i], valCeven[i], inv2);
        mul(valCodd[i], valCodd[i], inv2pt);
    }

#ifdef PROFILE_ON
    t = GetWallTime()-t;
    std::cout << "Pointwise prods: " << t << std::endl;

    t = GetWallTime();
#endif // PROFILE_ON

    // interpolate to find the even/odd parts of c
    mul(tmp_mat_even, iv, valCeven);
    mul(tmp_mat_odd, iv, valCodd);

    c.SetDims(m, p);
    ell = 0;
    for (long u = 0; u < m; ++u)
        for (long v = 0; v < p; ++v, ++ell)
        {
            c[u][v].SetLength(2*nb_points);
            for (long i = 0; i < nb_points; ++i)
            {
                c[u][v][2*i] = tmp_mat_even[i][ell];
                c[u][v][2*i+1] = tmp_mat_odd[i][ell];
            }
            c[u][v].normalize();
        }

#ifdef PROFILE_ON
    t = GetWallTime()-t;
    std::cout << "Interpolate: " << t << std::endl;
#endif // PROFILE_ON
}

//// attempt with just one matrix; doesn't seem to really gain something
//void multiply_evaluate_dense2_new(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b)
//{
//    // dimensions and degrees
//    const long dA = deg(a);
//    const long dB = deg(b);
//    const long sA = dA/2+1;
//    const long sB = dB/2+1;
//    const long m = a.NumRows();
//    const long n = a.NumCols();
//    const long p = b.NumCols();
//    // for multiplying rows by constant (inverse of 2)
//    const long q = zz_p::modulus();
//    mulmod_t qinv = zz_p::ModulusInverse();
//    long inv2 = rep(inv(to_zz_p(2)));
//    mulmod_precon_t precon_inv2 = PrepMulModPrecon(inv2, q, qinv);
//
//    long ell;
//
//#ifdef PROFILE_ON
//    double t = GetWallTime();
//#endif // PROFILE_ON
//    Mat<zz_p> vA, vB, iv;
//    vandermonde2(vA, vB, iv, dA, dB);
//#ifdef PROFILE_ON
//    t = GetWallTime()-t;
//    std::cout << std::endl << "Build vdmd mat: " << t << std::endl;
//#endif // PROFILE_ON
//
//    const long nb_points = vA.NumRows();
//
//#ifdef PROFILE_ON
//    double t_copy, t_mul;
//    t = GetWallTime();
//#endif // PROFILE_ON
//
//    // evaluation of matrix a (even/odd part):
//    // build tmp_mat_even/odd, whose column ell = i*n + j is the coefficient
//    // vector of the even/odd part of a[i][j], even if ell is even, odd if ell
//    // is odd (padded with zeroes up to length dA/2+1 if necessary)
//    Mat<zz_p> tmp_mat(INIT_SIZE, sA, 2*m*n);
//    // --> we use dimension sA=dA/2+1 for odd part as well to make both parts
//    // fit in the same matrix, even though this should be (dA-1)/2+1
//    ell = 0;
//    for (long i = 0; i < m; ++i)
//        for (long j = 0; j < n; ++j, ell+=2)
//        {
//            const long d = deg(a[i][j]);
//            for (long k = 0; 2*k < d; ++k)
//            {
//                tmp_mat[k][ell] = a[i][j][2*k]; // even coeff
//                tmp_mat[k][ell+1] = a[i][j][2*k+1]; // odd coeff
//            }
//            if (d % 2 == 0) // d even
//                tmp_mat[d/2][ell] = a[i][j][d];
//        }
//
//    // note: d = deg(a[i][j]) is -1 if a[i][j] == 0
//    // all non-touched entries already zero since tmp_mat was initialized as zero
//
//#ifdef PROFILE_ON
//    t_copy = GetWallTime()-t;
//    t = GetWallTime();
//#endif // PROFILE_ON
//
//    // valA: column with even (resp. odd) index ell = i*n + j contains the
//    // evaluations of even (resp. odd) part of a[i][j]
//    Mat<zz_p> valA;
//    mul(valA, vA, tmp_mat);
//#ifdef PROFILE_ON
//    t_mul = GetWallTime()-t;
//    t = GetWallTime();
//#endif // PROFILE_ON
//
//    // evaluation of matrix b (even/odd part): build tmp_mat, whose column ell
//    // = i*n + j is the coefficient vector of the even/odd part of b[i][j]
//    // (even if ell is even, odd if ell is odd)
//    // (padded with zeroes up to length dB/2+1 if necessary); as above, to make
//    // these fit in a single matrix, we use dimension dB/2+1 for both odd and
//    // even parts
//
//    tmp_mat.SetDims(sB, 2 * n * p);
//    ell = 0;
//    for (long i = 0; i < n; ++i)
//        for (long j = 0; j < p; ++j, ell+=2)
//        {
//            const long d = deg(b[i][j]); // -1 if b[i][j] == 0
//            long k = 0;
//            for (; 2*k < d; ++k)
//            {
//                tmp_mat[k][ell] = b[i][j][2*k];
//                tmp_mat[k][ell+1] = b[i][j][2*k+1];
//            }
//            if (2*k == d)
//            {
//                tmp_mat[k][ell] = b[i][j][d];
//                clear(tmp_mat[k][ell+1]);
//                ++k;
//            }
//            for (; k < sB; ++k)
//            {
//                clear(tmp_mat[k][ell]);
//                clear(tmp_mat[k][ell+1]);
//            }
//        }
//#ifdef PROFILE_ON
//    t_copy += GetWallTime()-t;
//    t = GetWallTime();
//#endif // PROFILE_ON
//
//    // valB: column with even (resp. odd) index ell = i*n + j contains the
//    // evaluations of even (resp. odd) part of b[i][j]
//    Mat<zz_p> valB;
//    mul(valB, vB, tmp_mat);
//
//#ifdef PROFILE_ON
//    t_mul += GetWallTime()-t;
//    std::cout << "Evals of a and b(data movements): " << t_copy << std::endl;
//    std::cout << "Evals of a and b(mul): " << t_mul << std::endl;
//    
//    t = GetWallTime();
//#endif // PROFILE_ON
//
//    // perform the pointwise products for the points 1, 2, ..., nb_points
//    Mat<zz_p> valApeven(INIT_SIZE, m, n);
//    Mat<zz_p> valApodd(INIT_SIZE, m, n);
//    Mat<zz_p> valBpeven(INIT_SIZE, n, p);
//    Mat<zz_p> valBpodd(INIT_SIZE, n, p);
//    Mat<zz_p> valC(INIT_SIZE, nb_points, 2 * m * p);
//    Mat<zz_p> valCpeven, valCpodd;
//
//    // from now on we will only use the odd part of valA and valB with row i
//    // multiplied by i+1
//    for (long i = 0; i < nb_points; ++i)
//    {
//        mulmod_precon_t precon = PrepMulModPrecon(i+1, q, qinv);
//        for (long j = 1; j < valA.NumCols(); j+=2)
//            valA[i][j].LoopHole() = MulModPrecon(rep(valA[i][j]), i+1, q, precon);
//        for (long j = 1; j < valB.NumCols(); j+=2)
//            valB[i][j].LoopHole() = MulModPrecon(rep(valB[i][j]), i+1, q, precon);
//    }
//
//    for (long i = 0; i < nb_points; ++i)
//    {
//        // a evaluated at point i (which is i+1)
//        ell = 0;
//        for (long u = 0; u < m; ++u)
//            for (long v = 0; v < n; ++v, ell+=2)
//            {
//                add(valApeven[u][v], valA[i][ell], valA[i][ell+1]);
//                sub(valApodd[u][v], valA[i][ell], valA[i][ell+1]);
//            }
//
//        // b evaluated at point i (which is i+1)
//        ell = 0;
//        for (long u = 0; u < n; ++u)
//            for (long v = 0; v < p; ++v, ell+=2)
//            {
//                add(valBpeven[u][v], valB[i][ell], valB[i][ell+1]);
//                sub(valBpodd[u][v], valB[i][ell], valB[i][ell+1]);
//            }
//
//        // a*b evaluated at point i (which is i+1)
//        mul(valCpeven, valApeven, valBpeven);
//        mul(valCpodd, valApodd, valBpodd);
//
//        // copy this into valC: column ell = i*n + j contains the evaluations
//        // of the entry i,j of c = a*b
//        // valCeven = (even + odd)/2
//        // valCodd = (even - odd)/2x
//        ell = 0;
//        long inv2pt = rep(inv(to_zz_p(2*(i+1))));
//        mulmod_precon_t precon = PrepMulModPrecon(inv2pt, q, qinv);
//        for (long u = 0; u < m; ++u)
//            for (long v = 0; v < p; ++v, ell+=2)
//            {
//                add(valC[i][ell], valCpeven[u][v], valCpodd[u][v]);
//                sub(valC[i][ell+1], valCpeven[u][v], valCpodd[u][v]);
//                valC[i][ell].LoopHole() = MulModPrecon(rep(valC[i][ell]), inv2, q, precon_inv2);
//                valC[i][ell+1].LoopHole() = MulModPrecon(rep(valC[i][ell+1]), inv2pt, q, precon);
//            }
//    }
//
//#ifdef PROFILE_ON
//    t = GetWallTime()-t;
//    std::cout << "Pointwise prods: " << t << std::endl;
//
//    t = GetWallTime();
//#endif // PROFILE_ON
//
//    // interpolate to find the even/odd parts of c
//    mul(tmp_mat, iv, valC);
//
//    c.SetDims(m, p);
//    ell = 0;
//    for (long u = 0; u < m; ++u)
//        for (long v = 0; v < p; ++v, ell+=2)
//        {
//            c[u][v].SetLength(2*nb_points);
//            for (long i = 0; i < nb_points; ++i)
//            {
//                c[u][v][2*i] = tmp_mat[i][ell];
//                c[u][v][2*i+1] = tmp_mat[i][ell+1];
//            }
//            c[u][v].normalize();
//        }
//
//#ifdef PROFILE_ON
//    t = GetWallTime()-t;
//    std::cout << "Interpolate: " << t << std::endl;
//#endif // PROFILE_ON
//}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
