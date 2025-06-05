#include <flint/nmod_poly_mat.h>

/** Middle product for polynomial matrices
 *  sets C = ((A * B) div x^dA) mod x^(dB+1), assuming deg(A) <= dA and deg(B) <= dA + dB
 *  output can alias input
 *  naive implementation (multiply, shift, truncate)
 */
void nmod_poly_mat_middle_product_naive(nmod_poly_mat_t C, const nmod_poly_mat_t A, const nmod_poly_mat_t B,
                                        const ulong dA, const ulong dB)
{
    nmod_poly_mat_mul(C, A, B);
    nmod_poly_mat_shift_right(C, C, dA);
    nmod_poly_mat_truncate(C, dB + 1);
}
