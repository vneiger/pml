#ifndef MAT_LZZ_PX_SEQUENCE__H
#define MAT_LZZ_PX_SEQUENCE__H

typedef Vec<Vec<Vec<zz_p>>> Coeffs;

static inline void MulMod_simple(zz_pX& upper, const zz_pX& a, const zz_pXModulus& F)
{
    long n = F.n;
    zz_pX P1(INIT_SIZE, n);
    fftRep R1(INIT_SIZE, F.l);
    TofftRep_trunc(R1, a, F.l, max(1L << F.k, 2*n-2));
    mul(R1, R1, F.HRep);
    FromfftRep(upper, R1, n-1, 2*n-3);
}

void MulMod_local(zz_pX& x, zz_pX& upper, const zz_pX& a, const zz_pXMultiplier& B, const zz_pXModulus& F);

/*------------------------------------------------------------*/
/* matrix multiplication using the algorithm of Giorgi et al. */
/* uses matrix multiplication for evaluation and interpolation*/
/*------------------------------------------------------------*/
void mul(Mat<zz_pX>& c, const Mat<zz_pX>& a, const Mat<zz_pX>& b);

/*------------------------------------------------------------*/
/* matrix multiplication using the algorithm of Giorgi et al. */
/* uses matrix multiplication for evaluation and interpolation*/
/*                                                            */
/* variant does one product and a row-shifted product at once */
/*------------------------------------------------------------*/
void mul_special(Mat<zz_pX>& c, Mat<zz_pX>& cs, const Mat<zz_pX>& a, const Mat<zz_pX>& b);

/*------------------------------------------------------------*/
/* generates (t*a^(s*i) mod g) for 0 <= i < l                 */
/* if need_upper_product = 1, also returns                    */
/*   S * rev(t * a^(s*i) mod g) mod x^(n-1), 0 <= i < l       */
/*------------------------------------------------------------*/
void gen_pows (Vec<zz_pX> &pow, Vec<zz_pX>&upper,
               const zz_pX &t,
               const zz_pX &a,
               const long s,
               const zz_pX &g,
               const long l, const long need_upper_products = 0);

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
void get_quos (Mat<zz_pX> &quos,
               const Vec<zz_pX> &alphas,
               const Vec<zz_pX> &As,
               const Vec<zz_pX> &upper,
               const zz_pX &g,
               const long m);
/*------------------------------------------------------------*/
/* Generates the sequence of first mxm coefficients of        */
/* t a^k mod g for 0 <= k <= 2(n/m)                           */
/*------------------------------------------------------------*/
void gen_sequence (Coeffs &res,
                   const zz_pX &t,
                   const zz_pX &a,
                   const zz_pX &g,
                   const long m);

/*------------------------------------------------------------*/
/* Generates the sequence of first mxm coefficients of        */
/* x^j a^k mod g for 0 <= k <= 2(n/m), 0 <= j < m.            */
/* Output is stored in a vector of matrices such that (i,j)-th*/
/* entry of the k-th matrix holds the j-th coefficient of     */
/* x^i a^(2d-k) mod g                                         */
/*------------------------------------------------------------*/
void gen_sequence (Vec<Mat<zz_p>> &mats, Vec<Coeffs> &res,
                   const zz_pX &a,
                   const zz_pX &g,
                   const long m);

void get_sequence_naive (Vec<Coeffs> &res,
                         const zz_pX &a,
                         const zz_pX &g,
                         const long m);


void print(const Vec<Coeffs> &c);
void print(const zz_pX &p);

#endif /* ifndef MAT_LZZ_PX_SEQUENCE__H */

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
