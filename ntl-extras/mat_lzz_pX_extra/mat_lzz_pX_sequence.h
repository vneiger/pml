#ifndef MAT_LZZ_PX_SEQUENCE__H
#define MAT_LZZ_PX_SEQUENCE__H

/** \brief Computing block-Wiedemann like sequences for balanced bases.
 *
 * \file mat_lzz_pX_sequence.h
 * \author Seung Gyu Hyun, Vincent Neiger, Eric Schost
 * \version 0.1
 * \date 2019-02-01
 *
 * \todo more detailed description
 *
 * \todo clean code (e.g. redundant with matrix Pade)
 */

#include <NTL/matrix.h>
#include <NTL/lzz_pX.h>

NTL_CLIENT

typedef Vec<Vec<Vec<zz_p>>> Coeffs;
/** For generating special Block-Wiedemann sequences.
 *
 * This generates (t*a**(s*i) mod g) for 0 <= i < l;
 * if `need_upper_product` is true, also returns
 * S * rev(t * a^(s*i) mod g) mod x^(n-1), 0 <= i < l.
 *
 * \param[out] pow vector of polynomials such that pow[i] = t*a**(s*i) mod g
 * \param[out] upper vector of polynomials holding the "upper products"
 * \param[in] t pre-multiplier
 * \param[in] a polynomial to be powered
 * \param[in] s starting power of a
 * \param[in] g modulus
 * \param[in] l number of powers needed
 * \param[in] need_upper_products set to `true` if "upper products" are needed
 */
void gen_pows(
              Vec<zz_pX> &pow,
              Vec<zz_pX>&upper,
              const zz_pX &t,
              const zz_pX &a,
              const long s,
              const zz_pX &g,
              const long l,
              const bool need_upper_products = false
             );

/** For generating special Block-Wiedemann sequences.
 *
 * Generates the sequence of first mxm coefficients of
 * x^j a^k mod g for 0 <= k <= 2(n/m), 0 <= j < m.
 * Output is stored in a vector of matrices such that (i,j)-th
 * entry of the k-th matrix holds the j-th coefficient of
 * x^i a^(2d-k) mod g. Uses a baby step giant step approach
 *
 * \param[out] mats stores the output
 * \param[in] t pre-multiplier
 * \param[in] a polynomial to be powered
 * \param[in] g modulus
 */
void gen_sequence(
                  Vec<Mat<zz_p>> &mats,
                  const zz_pX &a,
                  const zz_pX &g,
                  const long m
                 );


/** For generating special Block-Wiedemann sequences.
 *
 * Generates the sequence of first mxm coefficients of
 * x^j a^k mod g for 0 <= k <= 2(n/m), 0 <= j < m.
 * Output is stored in a vector of matrices such that (i,j)-th
 * entry of the k-th matrix holds the j-th coefficient of
 * x^i a^(2d-k) mod g. Uses the naive algorithm
 *
 * \param[out] res stores the output
 * \param[in] a polynomial to be powered
 * \param[in] g modulus
 * \param[in] m dimension
 */
void get_sequence_naive(
                        Vec<Coeffs> &res,
                        const zz_pX &a,
                        const zz_pX &g,
                        const long m
                       );


/** Computes the matrix reconstruction
 *
 * Given a sequence of (scalar) matrices, computes the matrix
 * denominator using approximant basis
 *
 * \param[out] basis basis of the sequence
 * \param[in] seq sequence of matrices
 */
void matrix_recon_approximation(Mat<zz_pX> &basis, const Vec<Mat<zz_p>> &seq);

/** Computes the matrix reconstruction
 *
 * Given a set of matrix evaluations, evaluated at some points `pts`, this
 * reconstructs the matrix denominator using interpolant basis
 *
 * \param[out] basis basis of the sequence
 * \param[in] pts evaluation points
 * \param[in] seq matrix evaluations
 */
void matrix_recon_interpolation(
                                Mat<zz_pX> &basis,
                                const Vec<zz_p> &pts,
                                Vec<Mat<zz_p>> &seq
                               );

/** Computes the matrix reconstruction
 *
 * Given a set of matrix evaluations, evaluated at some points `pts`, this
 * reconstructs the matrix denominator using interpolant basis with geometric
 * points based on `r`
 *
 * \param[out] basis basis of the sequence
 * \param[in] pts evaluation points
 * \param[in] seq matrix evaluations
 */
void matrix_recon_interpolation_geometric(
                                          Mat<zz_pX> &basis,
                                          const Vec<zz_p> &pts,
                                          const zz_p& r,
                                          Vec<Mat<zz_p>> &seq
                                         );

#endif /* ifndef MAT_LZZ_PX_SEQUENCE__H */

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
