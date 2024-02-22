#include <flint/nmod.h>
#include <flint/nmod_vec.h>

#include "nmod_vec_extra.h"

#define __ll_lowhi_parts26(tlo,thi,t)     \
      thi = (uint) ((t) >> 26);           \
      tlo = ((uint)(t)) & 0x3FFFFFF;

/* ------------------------------------------------------------ */
/* number of limbs needed for a dot product of length len       */
/* all entries 1st vector have <= max_bits1 bits <= FLINT_BITS  */
/* all entries 2nd vector have <= max_bits2 bits <= FLINT_BITS  */
/* returns 0, 1, 2, 3                                           */
/* ------------------------------------------------------------ */
static inline
ulong _nmod_vec_dot_bound_limbs_unbalanced(ulong len, ulong max_bits1, ulong max_bits2)
{
    const mp_limb_t a1 = (max_bits1 == FLINT_BITS) ? (UWORD_MAX) : (UWORD(1) << max_bits1) - 1;
    const mp_limb_t a2 = (max_bits2 == FLINT_BITS) ? (UWORD_MAX) : (UWORD(1) << max_bits2) - 1;

    mp_limb_t t2, t1, t0, u1, u0;
    umul_ppmm(t1, t0, a1, a2);
    umul_ppmm(t2, t1, t1, len);
    umul_ppmm(u1, u0, t0, len);
    add_sssaaaaaa(t2, t1, t0,  t2, t1, UWORD(0),  UWORD(0), u1, u0);

    if (t2 != 0)
        return 3;
    if (t1 != 0)
        return 2;
    return (t0 != 0);
}


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/*                       SINGLE LIMB                          */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

// same as v_8_32 below
static inline
void _nmod_vec_dot_product_multi_1(mp_ptr uv, mp_srcptr u, mp_srcptr * v,
                                   ulong len, ulong k, nmod_t mod)
{
    _nmod_vec_zero(uv, k);
    ulong i = 0;
    for (; i+7 < len; i += 8)
    {
        ulong j = 0;
        for (; j+31 < k; j += 32)
        {
            uv[j+0] +=   u[i+0] * v[i+0][j+0]
                       + u[i+1] * v[i+1][j+0]
                       + u[i+2] * v[i+2][j+0]
                       + u[i+3] * v[i+3][j+0]
                       + u[i+4] * v[i+4][j+0]
                       + u[i+5] * v[i+5][j+0]
                       + u[i+6] * v[i+6][j+0]
                       + u[i+7] * v[i+7][j+0];
            uv[j+1] +=   u[i+0] * v[i+0][j+1]
                       + u[i+1] * v[i+1][j+1]
                       + u[i+2] * v[i+2][j+1]
                       + u[i+3] * v[i+3][j+1]
                       + u[i+4] * v[i+4][j+1]
                       + u[i+5] * v[i+5][j+1]
                       + u[i+6] * v[i+6][j+1]
                       + u[i+7] * v[i+7][j+1];
            uv[j+2] +=   u[i+0] * v[i+0][j+2]
                       + u[i+1] * v[i+1][j+2]
                       + u[i+2] * v[i+2][j+2]
                       + u[i+3] * v[i+3][j+2]
                       + u[i+4] * v[i+4][j+2]
                       + u[i+5] * v[i+5][j+2]
                       + u[i+6] * v[i+6][j+2]
                       + u[i+7] * v[i+7][j+2];
            uv[j+3] +=   u[i+0] * v[i+0][j+3]
                       + u[i+1] * v[i+1][j+3]
                       + u[i+2] * v[i+2][j+3]
                       + u[i+3] * v[i+3][j+3]
                       + u[i+4] * v[i+4][j+3]
                       + u[i+5] * v[i+5][j+3]
                       + u[i+6] * v[i+6][j+3]
                       + u[i+7] * v[i+7][j+3];
            uv[j+4] +=   u[i+0] * v[i+0][j+4]
                       + u[i+1] * v[i+1][j+4]
                       + u[i+2] * v[i+2][j+4]
                       + u[i+3] * v[i+3][j+4]
                       + u[i+4] * v[i+4][j+4]
                       + u[i+5] * v[i+5][j+4]
                       + u[i+6] * v[i+6][j+4]
                       + u[i+7] * v[i+7][j+4];
            uv[j+5] +=   u[i+0] * v[i+0][j+5]
                       + u[i+1] * v[i+1][j+5]
                       + u[i+2] * v[i+2][j+5]
                       + u[i+3] * v[i+3][j+5]
                       + u[i+4] * v[i+4][j+5]
                       + u[i+5] * v[i+5][j+5]
                       + u[i+6] * v[i+6][j+5]
                       + u[i+7] * v[i+7][j+5];
            uv[j+6] +=   u[i+0] * v[i+0][j+6]
                       + u[i+1] * v[i+1][j+6]
                       + u[i+2] * v[i+2][j+6]
                       + u[i+3] * v[i+3][j+6]
                       + u[i+4] * v[i+4][j+6]
                       + u[i+5] * v[i+5][j+6]
                       + u[i+6] * v[i+6][j+6]
                       + u[i+7] * v[i+7][j+6];
            uv[j+7] +=   u[i+0] * v[i+0][j+7]
                       + u[i+1] * v[i+1][j+7]
                       + u[i+2] * v[i+2][j+7]
                       + u[i+3] * v[i+3][j+7]
                       + u[i+4] * v[i+4][j+7]
                       + u[i+5] * v[i+5][j+7]
                       + u[i+6] * v[i+6][j+7]
                       + u[i+7] * v[i+7][j+7];
            uv[j+8] +=   u[i+0] * v[i+0][j+8]
                       + u[i+1] * v[i+1][j+8]
                       + u[i+2] * v[i+2][j+8]
                       + u[i+3] * v[i+3][j+8]
                       + u[i+4] * v[i+4][j+8]
                       + u[i+5] * v[i+5][j+8]
                       + u[i+6] * v[i+6][j+8]
                       + u[i+7] * v[i+7][j+8];
            uv[j+9] +=   u[i+0] * v[i+0][j+9]
                       + u[i+1] * v[i+1][j+9]
                       + u[i+2] * v[i+2][j+9]
                       + u[i+3] * v[i+3][j+9]
                       + u[i+4] * v[i+4][j+9]
                       + u[i+5] * v[i+5][j+9]
                       + u[i+6] * v[i+6][j+9]
                       + u[i+7] * v[i+7][j+9];
            uv[j+10] +=  u[i+0] * v[i+0][j+10]
                       + u[i+1] * v[i+1][j+10]
                       + u[i+2] * v[i+2][j+10]
                       + u[i+3] * v[i+3][j+10]
                       + u[i+4] * v[i+4][j+10]
                       + u[i+5] * v[i+5][j+10]
                       + u[i+6] * v[i+6][j+10]
                       + u[i+7] * v[i+7][j+10];
            uv[j+11] +=  u[i+0] * v[i+0][j+11]
                       + u[i+1] * v[i+1][j+11]
                       + u[i+2] * v[i+2][j+11]
                       + u[i+3] * v[i+3][j+11]
                       + u[i+4] * v[i+4][j+11]
                       + u[i+5] * v[i+5][j+11]
                       + u[i+6] * v[i+6][j+11]
                       + u[i+7] * v[i+7][j+11];
            uv[j+12] +=  u[i+0] * v[i+0][j+12]
                       + u[i+1] * v[i+1][j+12]
                       + u[i+2] * v[i+2][j+12]
                       + u[i+3] * v[i+3][j+12]
                       + u[i+4] * v[i+4][j+12]
                       + u[i+5] * v[i+5][j+12]
                       + u[i+6] * v[i+6][j+12]
                       + u[i+7] * v[i+7][j+12];
            uv[j+13] +=  u[i+0] * v[i+0][j+13]
                       + u[i+1] * v[i+1][j+13]
                       + u[i+2] * v[i+2][j+13]
                       + u[i+3] * v[i+3][j+13]
                       + u[i+4] * v[i+4][j+13]
                       + u[i+5] * v[i+5][j+13]
                       + u[i+6] * v[i+6][j+13]
                       + u[i+7] * v[i+7][j+13];
            uv[j+14] +=  u[i+0] * v[i+0][j+14]
                       + u[i+1] * v[i+1][j+14]
                       + u[i+2] * v[i+2][j+14]
                       + u[i+3] * v[i+3][j+14]
                       + u[i+4] * v[i+4][j+14]
                       + u[i+5] * v[i+5][j+14]
                       + u[i+6] * v[i+6][j+14]
                       + u[i+7] * v[i+7][j+14];
            uv[j+15] +=  u[i+0] * v[i+0][j+15]
                       + u[i+1] * v[i+1][j+15]
                       + u[i+2] * v[i+2][j+15]
                       + u[i+3] * v[i+3][j+15]
                       + u[i+4] * v[i+4][j+15]
                       + u[i+5] * v[i+5][j+15]
                       + u[i+6] * v[i+6][j+15]
                       + u[i+7] * v[i+7][j+15];
            uv[j+16] +=  u[i+0] * v[i+0][j+16]
                       + u[i+1] * v[i+1][j+16]
                       + u[i+2] * v[i+2][j+16]
                       + u[i+3] * v[i+3][j+16]
                       + u[i+4] * v[i+4][j+16]
                       + u[i+5] * v[i+5][j+16]
                       + u[i+6] * v[i+6][j+16]
                       + u[i+7] * v[i+7][j+16];
            uv[j+17] +=  u[i+0] * v[i+0][j+17]
                       + u[i+1] * v[i+1][j+17]
                       + u[i+2] * v[i+2][j+17]
                       + u[i+3] * v[i+3][j+17]
                       + u[i+4] * v[i+4][j+17]
                       + u[i+5] * v[i+5][j+17]
                       + u[i+6] * v[i+6][j+17]
                       + u[i+7] * v[i+7][j+17];
            uv[j+18] +=  u[i+0] * v[i+0][j+18]
                       + u[i+1] * v[i+1][j+18]
                       + u[i+2] * v[i+2][j+18]
                       + u[i+3] * v[i+3][j+18]
                       + u[i+4] * v[i+4][j+18]
                       + u[i+5] * v[i+5][j+18]
                       + u[i+6] * v[i+6][j+18]
                       + u[i+7] * v[i+7][j+18];
            uv[j+19] +=  u[i+0] * v[i+0][j+19]
                       + u[i+1] * v[i+1][j+19]
                       + u[i+2] * v[i+2][j+19]
                       + u[i+3] * v[i+3][j+19]
                       + u[i+4] * v[i+4][j+19]
                       + u[i+5] * v[i+5][j+19]
                       + u[i+6] * v[i+6][j+19]
                       + u[i+7] * v[i+7][j+19];
            uv[j+20] +=  u[i+0] * v[i+0][j+20]
                       + u[i+1] * v[i+1][j+20]
                       + u[i+2] * v[i+2][j+20]
                       + u[i+3] * v[i+3][j+20]
                       + u[i+4] * v[i+4][j+20]
                       + u[i+5] * v[i+5][j+20]
                       + u[i+6] * v[i+6][j+20]
                       + u[i+7] * v[i+7][j+20];
            uv[j+21] +=  u[i+0] * v[i+0][j+21]
                       + u[i+1] * v[i+1][j+21]
                       + u[i+2] * v[i+2][j+21]
                       + u[i+3] * v[i+3][j+21]
                       + u[i+4] * v[i+4][j+21]
                       + u[i+5] * v[i+5][j+21]
                       + u[i+6] * v[i+6][j+21]
                       + u[i+7] * v[i+7][j+21];
            uv[j+22] +=  u[i+0] * v[i+0][j+22]
                       + u[i+1] * v[i+1][j+22]
                       + u[i+2] * v[i+2][j+22]
                       + u[i+3] * v[i+3][j+22]
                       + u[i+4] * v[i+4][j+22]
                       + u[i+5] * v[i+5][j+22]
                       + u[i+6] * v[i+6][j+22]
                       + u[i+7] * v[i+7][j+22];
            uv[j+23] +=  u[i+0] * v[i+0][j+23]
                       + u[i+1] * v[i+1][j+23]
                       + u[i+2] * v[i+2][j+23]
                       + u[i+3] * v[i+3][j+23]
                       + u[i+4] * v[i+4][j+23]
                       + u[i+5] * v[i+5][j+23]
                       + u[i+6] * v[i+6][j+23]
                       + u[i+7] * v[i+7][j+23];
            uv[j+24] +=  u[i+0] * v[i+0][j+24]
                       + u[i+1] * v[i+1][j+24]
                       + u[i+2] * v[i+2][j+24]
                       + u[i+3] * v[i+3][j+24]
                       + u[i+4] * v[i+4][j+24]
                       + u[i+5] * v[i+5][j+24]
                       + u[i+6] * v[i+6][j+24]
                       + u[i+7] * v[i+7][j+24];
            uv[j+25] +=  u[i+0] * v[i+0][j+25]
                       + u[i+1] * v[i+1][j+25]
                       + u[i+2] * v[i+2][j+25]
                       + u[i+3] * v[i+3][j+25]
                       + u[i+4] * v[i+4][j+25]
                       + u[i+5] * v[i+5][j+25]
                       + u[i+6] * v[i+6][j+25]
                       + u[i+7] * v[i+7][j+25];
            uv[j+26] +=  u[i+0] * v[i+0][j+26]
                       + u[i+1] * v[i+1][j+26]
                       + u[i+2] * v[i+2][j+26]
                       + u[i+3] * v[i+3][j+26]
                       + u[i+4] * v[i+4][j+26]
                       + u[i+5] * v[i+5][j+26]
                       + u[i+6] * v[i+6][j+26]
                       + u[i+7] * v[i+7][j+26];
            uv[j+27] +=  u[i+0] * v[i+0][j+27]
                       + u[i+1] * v[i+1][j+27]
                       + u[i+2] * v[i+2][j+27]
                       + u[i+3] * v[i+3][j+27]
                       + u[i+4] * v[i+4][j+27]
                       + u[i+5] * v[i+5][j+27]
                       + u[i+6] * v[i+6][j+27]
                       + u[i+7] * v[i+7][j+27];
            uv[j+28] +=  u[i+0] * v[i+0][j+28]
                       + u[i+1] * v[i+1][j+28]
                       + u[i+2] * v[i+2][j+28]
                       + u[i+3] * v[i+3][j+28]
                       + u[i+4] * v[i+4][j+28]
                       + u[i+5] * v[i+5][j+28]
                       + u[i+6] * v[i+6][j+28]
                       + u[i+7] * v[i+7][j+28];
            uv[j+29] +=  u[i+0] * v[i+0][j+29]
                       + u[i+1] * v[i+1][j+29]
                       + u[i+2] * v[i+2][j+29]
                       + u[i+3] * v[i+3][j+29]
                       + u[i+4] * v[i+4][j+29]
                       + u[i+5] * v[i+5][j+29]
                       + u[i+6] * v[i+6][j+29]
                       + u[i+7] * v[i+7][j+29];
            uv[j+30] +=  u[i+0] * v[i+0][j+30]
                       + u[i+1] * v[i+1][j+30]
                       + u[i+2] * v[i+2][j+30]
                       + u[i+3] * v[i+3][j+30]
                       + u[i+4] * v[i+4][j+30]
                       + u[i+5] * v[i+5][j+30]
                       + u[i+6] * v[i+6][j+30]
                       + u[i+7] * v[i+7][j+30];
            uv[j+31] +=  u[i+0] * v[i+0][j+31]
                       + u[i+1] * v[i+1][j+31]
                       + u[i+2] * v[i+2][j+31]
                       + u[i+3] * v[i+3][j+31]
                       + u[i+4] * v[i+4][j+31]
                       + u[i+5] * v[i+5][j+31]
                       + u[i+6] * v[i+6][j+31]
                       + u[i+7] * v[i+7][j+31];
        }
        for (; j < k; j++)
        {
            uv[j] +=   u[i+0] * v[i+0][j]
                     + u[i+1] * v[i+1][j]
                     + u[i+2] * v[i+2][j]
                     + u[i+3] * v[i+3][j]
                     + u[i+4] * v[i+4][j]
                     + u[i+5] * v[i+5][j]
                     + u[i+6] * v[i+6][j]
                     + u[i+7] * v[i+7][j];
        }
    }
    for (; i < len; i++)
    {
        ulong j = 0;
        for (; j+31 < k; j += 32)
        {
            uv[j+0] += u[i] * v[i][j+0];
            uv[j+1] += u[i] * v[i][j+1];
            uv[j+2] += u[i] * v[i][j+2];
            uv[j+3] += u[i] * v[i][j+3];
            uv[j+4] += u[i] * v[i][j+4];
            uv[j+5] += u[i] * v[i][j+5];
            uv[j+6] += u[i] * v[i][j+6];
            uv[j+7] += u[i] * v[i][j+7];
            uv[j+8] += u[i] * v[i][j+8];
            uv[j+9] += u[i] * v[i][j+9];
            uv[j+10] += u[i] * v[i][j+10];
            uv[j+11] += u[i] * v[i][j+11];
            uv[j+12] += u[i] * v[i][j+12];
            uv[j+13] += u[i] * v[i][j+13];
            uv[j+14] += u[i] * v[i][j+14];
            uv[j+15] += u[i] * v[i][j+15];
            uv[j+16] += u[i] * v[i][j+16];
            uv[j+17] += u[i] * v[i][j+17];
            uv[j+18] += u[i] * v[i][j+18];
            uv[j+19] += u[i] * v[i][j+19];
            uv[j+20] += u[i] * v[i][j+20];
            uv[j+21] += u[i] * v[i][j+21];
            uv[j+22] += u[i] * v[i][j+22];
            uv[j+23] += u[i] * v[i][j+23];
            uv[j+24] += u[i] * v[i][j+24];
            uv[j+25] += u[i] * v[i][j+25];
            uv[j+26] += u[i] * v[i][j+26];
            uv[j+27] += u[i] * v[i][j+27];
            uv[j+28] += u[i] * v[i][j+28];
            uv[j+29] += u[i] * v[i][j+29];
            uv[j+30] += u[i] * v[i][j+30];
            uv[j+31] += u[i] * v[i][j+31];
        }
        for (; j < k; j++)
            uv[j] += u[i] * v[i][j];
    }

    for (ulong j = 0; j < k; j++)
        NMOD_RED(uv[j], uv[j], mod);
}

// TODO: variants below, some missing + thresholds to be determined, see header file
void _nmod_vec_dot_product_multi_1_v1_8(mp_ptr uv, mp_srcptr u, mp_srcptr * v,
                                   ulong len, ulong k, nmod_t mod)
{
    _nmod_vec_zero(uv, k);
    for (ulong i = 0; i < len; i++)
    {
        ulong j = 0;
        for (; j+7 < k; j += 8)
        {
            uv[j+0] += u[i] * v[i][j+0];
            uv[j+1] += u[i] * v[i][j+1];
            uv[j+2] += u[i] * v[i][j+2];
            uv[j+3] += u[i] * v[i][j+3];
            uv[j+4] += u[i] * v[i][j+4];
            uv[j+5] += u[i] * v[i][j+5];
            uv[j+6] += u[i] * v[i][j+6];
            uv[j+7] += u[i] * v[i][j+7];
        }
        for (; j < k; j++)
            uv[j] += u[i] * v[i][j];
    }

    for (ulong j = 0; j < k; j++)
        NMOD_RED(uv[j], uv[j], mod);
}


void _nmod_vec_dot_product_multi_1_v8_8(mp_ptr uv, mp_srcptr u, mp_srcptr * v,
                                       ulong len, ulong k, nmod_t mod)
{
    _nmod_vec_zero(uv, k);
    ulong i = 0;
    for (; i+7 < len; i += 8)
    {
        ulong j = 0;
        for (; j+7 < k; j += 8)
        {
            uv[j+0] +=   u[i+0] * v[i+0][j+0]
                       + u[i+1] * v[i+1][j+0]
                       + u[i+2] * v[i+2][j+0]
                       + u[i+3] * v[i+3][j+0]
                       + u[i+4] * v[i+4][j+0]
                       + u[i+5] * v[i+5][j+0]
                       + u[i+6] * v[i+6][j+0]
                       + u[i+7] * v[i+7][j+0];
            uv[j+1] +=   u[i+0] * v[i+0][j+1]
                       + u[i+1] * v[i+1][j+1]
                       + u[i+2] * v[i+2][j+1]
                       + u[i+3] * v[i+3][j+1]
                       + u[i+4] * v[i+4][j+1]
                       + u[i+5] * v[i+5][j+1]
                       + u[i+6] * v[i+6][j+1]
                       + u[i+7] * v[i+7][j+1];
            uv[j+2] +=   u[i+0] * v[i+0][j+2]
                       + u[i+1] * v[i+1][j+2]
                       + u[i+2] * v[i+2][j+2]
                       + u[i+3] * v[i+3][j+2]
                       + u[i+4] * v[i+4][j+2]
                       + u[i+5] * v[i+5][j+2]
                       + u[i+6] * v[i+6][j+2]
                       + u[i+7] * v[i+7][j+2];
            uv[j+3] +=   u[i+0] * v[i+0][j+3]
                       + u[i+1] * v[i+1][j+3]
                       + u[i+2] * v[i+2][j+3]
                       + u[i+3] * v[i+3][j+3]
                       + u[i+4] * v[i+4][j+3]
                       + u[i+5] * v[i+5][j+3]
                       + u[i+6] * v[i+6][j+3]
                       + u[i+7] * v[i+7][j+3];
            uv[j+4] +=   u[i+0] * v[i+0][j+4]
                       + u[i+1] * v[i+1][j+4]
                       + u[i+2] * v[i+2][j+4]
                       + u[i+3] * v[i+3][j+4]
                       + u[i+4] * v[i+4][j+4]
                       + u[i+5] * v[i+5][j+4]
                       + u[i+6] * v[i+6][j+4]
                       + u[i+7] * v[i+7][j+4];
            uv[j+5] +=   u[i+0] * v[i+0][j+5]
                       + u[i+1] * v[i+1][j+5]
                       + u[i+2] * v[i+2][j+5]
                       + u[i+3] * v[i+3][j+5]
                       + u[i+4] * v[i+4][j+5]
                       + u[i+5] * v[i+5][j+5]
                       + u[i+6] * v[i+6][j+5]
                       + u[i+7] * v[i+7][j+5];
            uv[j+6] +=   u[i+0] * v[i+0][j+6]
                       + u[i+1] * v[i+1][j+6]
                       + u[i+2] * v[i+2][j+6]
                       + u[i+3] * v[i+3][j+6]
                       + u[i+4] * v[i+4][j+6]
                       + u[i+5] * v[i+5][j+6]
                       + u[i+6] * v[i+6][j+6]
                       + u[i+7] * v[i+7][j+6];
            uv[j+7] +=   u[i+0] * v[i+0][j+7]
                       + u[i+1] * v[i+1][j+7]
                       + u[i+2] * v[i+2][j+7]
                       + u[i+3] * v[i+3][j+7]
                       + u[i+4] * v[i+4][j+7]
                       + u[i+5] * v[i+5][j+7]
                       + u[i+6] * v[i+6][j+7]
                       + u[i+7] * v[i+7][j+7];
        }
        for (; j < k; j++)
        {
            uv[j] +=   u[i+0] * v[i+0][j]
                     + u[i+1] * v[i+1][j]
                     + u[i+2] * v[i+2][j]
                     + u[i+3] * v[i+3][j]
                     + u[i+4] * v[i+4][j]
                     + u[i+5] * v[i+5][j]
                     + u[i+6] * v[i+6][j]
                     + u[i+7] * v[i+7][j];
        }
    }
    for (; i < len; i++)
    {
        ulong j = 0;
        for (; j+7 < k; j += 8)
        {
            uv[j+0] += u[i] * v[i][j+0];
            uv[j+1] += u[i] * v[i][j+1];
            uv[j+2] += u[i] * v[i][j+2];
            uv[j+3] += u[i] * v[i][j+3];
            uv[j+4] += u[i] * v[i][j+4];
            uv[j+5] += u[i] * v[i][j+5];
            uv[j+6] += u[i] * v[i][j+6];
            uv[j+7] += u[i] * v[i][j+7];
        }
        for (; j < k; j++)
            uv[j] += u[i] * v[i][j];
    }

    for (ulong j = 0; j < k; j++)
        NMOD_RED(uv[j], uv[j], mod);
}

// same as main function above (_multi_1)
void _nmod_vec_dot_product_multi_1_v8_32(mp_ptr uv, mp_srcptr u, mp_srcptr * v,
                                       ulong len, ulong k, nmod_t mod)
{
    _nmod_vec_zero(uv, k);
    ulong i = 0;
    for (; i+7 < len; i += 8)
    {
        ulong j = 0;
        for (; j+31 < k; j += 32)
        {
            uv[j+0] +=   u[i+0] * v[i+0][j+0]
                       + u[i+1] * v[i+1][j+0]
                       + u[i+2] * v[i+2][j+0]
                       + u[i+3] * v[i+3][j+0]
                       + u[i+4] * v[i+4][j+0]
                       + u[i+5] * v[i+5][j+0]
                       + u[i+6] * v[i+6][j+0]
                       + u[i+7] * v[i+7][j+0];
            uv[j+1] +=   u[i+0] * v[i+0][j+1]
                       + u[i+1] * v[i+1][j+1]
                       + u[i+2] * v[i+2][j+1]
                       + u[i+3] * v[i+3][j+1]
                       + u[i+4] * v[i+4][j+1]
                       + u[i+5] * v[i+5][j+1]
                       + u[i+6] * v[i+6][j+1]
                       + u[i+7] * v[i+7][j+1];
            uv[j+2] +=   u[i+0] * v[i+0][j+2]
                       + u[i+1] * v[i+1][j+2]
                       + u[i+2] * v[i+2][j+2]
                       + u[i+3] * v[i+3][j+2]
                       + u[i+4] * v[i+4][j+2]
                       + u[i+5] * v[i+5][j+2]
                       + u[i+6] * v[i+6][j+2]
                       + u[i+7] * v[i+7][j+2];
            uv[j+3] +=   u[i+0] * v[i+0][j+3]
                       + u[i+1] * v[i+1][j+3]
                       + u[i+2] * v[i+2][j+3]
                       + u[i+3] * v[i+3][j+3]
                       + u[i+4] * v[i+4][j+3]
                       + u[i+5] * v[i+5][j+3]
                       + u[i+6] * v[i+6][j+3]
                       + u[i+7] * v[i+7][j+3];
            uv[j+4] +=   u[i+0] * v[i+0][j+4]
                       + u[i+1] * v[i+1][j+4]
                       + u[i+2] * v[i+2][j+4]
                       + u[i+3] * v[i+3][j+4]
                       + u[i+4] * v[i+4][j+4]
                       + u[i+5] * v[i+5][j+4]
                       + u[i+6] * v[i+6][j+4]
                       + u[i+7] * v[i+7][j+4];
            uv[j+5] +=   u[i+0] * v[i+0][j+5]
                       + u[i+1] * v[i+1][j+5]
                       + u[i+2] * v[i+2][j+5]
                       + u[i+3] * v[i+3][j+5]
                       + u[i+4] * v[i+4][j+5]
                       + u[i+5] * v[i+5][j+5]
                       + u[i+6] * v[i+6][j+5]
                       + u[i+7] * v[i+7][j+5];
            uv[j+6] +=   u[i+0] * v[i+0][j+6]
                       + u[i+1] * v[i+1][j+6]
                       + u[i+2] * v[i+2][j+6]
                       + u[i+3] * v[i+3][j+6]
                       + u[i+4] * v[i+4][j+6]
                       + u[i+5] * v[i+5][j+6]
                       + u[i+6] * v[i+6][j+6]
                       + u[i+7] * v[i+7][j+6];
            uv[j+7] +=   u[i+0] * v[i+0][j+7]
                       + u[i+1] * v[i+1][j+7]
                       + u[i+2] * v[i+2][j+7]
                       + u[i+3] * v[i+3][j+7]
                       + u[i+4] * v[i+4][j+7]
                       + u[i+5] * v[i+5][j+7]
                       + u[i+6] * v[i+6][j+7]
                       + u[i+7] * v[i+7][j+7];
            uv[j+8] +=   u[i+0] * v[i+0][j+8]
                       + u[i+1] * v[i+1][j+8]
                       + u[i+2] * v[i+2][j+8]
                       + u[i+3] * v[i+3][j+8]
                       + u[i+4] * v[i+4][j+8]
                       + u[i+5] * v[i+5][j+8]
                       + u[i+6] * v[i+6][j+8]
                       + u[i+7] * v[i+7][j+8];
            uv[j+9] +=   u[i+0] * v[i+0][j+9]
                       + u[i+1] * v[i+1][j+9]
                       + u[i+2] * v[i+2][j+9]
                       + u[i+3] * v[i+3][j+9]
                       + u[i+4] * v[i+4][j+9]
                       + u[i+5] * v[i+5][j+9]
                       + u[i+6] * v[i+6][j+9]
                       + u[i+7] * v[i+7][j+9];
            uv[j+10] +=  u[i+0] * v[i+0][j+10]
                       + u[i+1] * v[i+1][j+10]
                       + u[i+2] * v[i+2][j+10]
                       + u[i+3] * v[i+3][j+10]
                       + u[i+4] * v[i+4][j+10]
                       + u[i+5] * v[i+5][j+10]
                       + u[i+6] * v[i+6][j+10]
                       + u[i+7] * v[i+7][j+10];
            uv[j+11] +=  u[i+0] * v[i+0][j+11]
                       + u[i+1] * v[i+1][j+11]
                       + u[i+2] * v[i+2][j+11]
                       + u[i+3] * v[i+3][j+11]
                       + u[i+4] * v[i+4][j+11]
                       + u[i+5] * v[i+5][j+11]
                       + u[i+6] * v[i+6][j+11]
                       + u[i+7] * v[i+7][j+11];
            uv[j+12] +=  u[i+0] * v[i+0][j+12]
                       + u[i+1] * v[i+1][j+12]
                       + u[i+2] * v[i+2][j+12]
                       + u[i+3] * v[i+3][j+12]
                       + u[i+4] * v[i+4][j+12]
                       + u[i+5] * v[i+5][j+12]
                       + u[i+6] * v[i+6][j+12]
                       + u[i+7] * v[i+7][j+12];
            uv[j+13] +=  u[i+0] * v[i+0][j+13]
                       + u[i+1] * v[i+1][j+13]
                       + u[i+2] * v[i+2][j+13]
                       + u[i+3] * v[i+3][j+13]
                       + u[i+4] * v[i+4][j+13]
                       + u[i+5] * v[i+5][j+13]
                       + u[i+6] * v[i+6][j+13]
                       + u[i+7] * v[i+7][j+13];
            uv[j+14] +=  u[i+0] * v[i+0][j+14]
                       + u[i+1] * v[i+1][j+14]
                       + u[i+2] * v[i+2][j+14]
                       + u[i+3] * v[i+3][j+14]
                       + u[i+4] * v[i+4][j+14]
                       + u[i+5] * v[i+5][j+14]
                       + u[i+6] * v[i+6][j+14]
                       + u[i+7] * v[i+7][j+14];
            uv[j+15] +=  u[i+0] * v[i+0][j+15]
                       + u[i+1] * v[i+1][j+15]
                       + u[i+2] * v[i+2][j+15]
                       + u[i+3] * v[i+3][j+15]
                       + u[i+4] * v[i+4][j+15]
                       + u[i+5] * v[i+5][j+15]
                       + u[i+6] * v[i+6][j+15]
                       + u[i+7] * v[i+7][j+15];
            uv[j+16] +=  u[i+0] * v[i+0][j+16]
                       + u[i+1] * v[i+1][j+16]
                       + u[i+2] * v[i+2][j+16]
                       + u[i+3] * v[i+3][j+16]
                       + u[i+4] * v[i+4][j+16]
                       + u[i+5] * v[i+5][j+16]
                       + u[i+6] * v[i+6][j+16]
                       + u[i+7] * v[i+7][j+16];
            uv[j+17] +=  u[i+0] * v[i+0][j+17]
                       + u[i+1] * v[i+1][j+17]
                       + u[i+2] * v[i+2][j+17]
                       + u[i+3] * v[i+3][j+17]
                       + u[i+4] * v[i+4][j+17]
                       + u[i+5] * v[i+5][j+17]
                       + u[i+6] * v[i+6][j+17]
                       + u[i+7] * v[i+7][j+17];
            uv[j+18] +=  u[i+0] * v[i+0][j+18]
                       + u[i+1] * v[i+1][j+18]
                       + u[i+2] * v[i+2][j+18]
                       + u[i+3] * v[i+3][j+18]
                       + u[i+4] * v[i+4][j+18]
                       + u[i+5] * v[i+5][j+18]
                       + u[i+6] * v[i+6][j+18]
                       + u[i+7] * v[i+7][j+18];
            uv[j+19] +=  u[i+0] * v[i+0][j+19]
                       + u[i+1] * v[i+1][j+19]
                       + u[i+2] * v[i+2][j+19]
                       + u[i+3] * v[i+3][j+19]
                       + u[i+4] * v[i+4][j+19]
                       + u[i+5] * v[i+5][j+19]
                       + u[i+6] * v[i+6][j+19]
                       + u[i+7] * v[i+7][j+19];
            uv[j+20] +=  u[i+0] * v[i+0][j+20]
                       + u[i+1] * v[i+1][j+20]
                       + u[i+2] * v[i+2][j+20]
                       + u[i+3] * v[i+3][j+20]
                       + u[i+4] * v[i+4][j+20]
                       + u[i+5] * v[i+5][j+20]
                       + u[i+6] * v[i+6][j+20]
                       + u[i+7] * v[i+7][j+20];
            uv[j+21] +=  u[i+0] * v[i+0][j+21]
                       + u[i+1] * v[i+1][j+21]
                       + u[i+2] * v[i+2][j+21]
                       + u[i+3] * v[i+3][j+21]
                       + u[i+4] * v[i+4][j+21]
                       + u[i+5] * v[i+5][j+21]
                       + u[i+6] * v[i+6][j+21]
                       + u[i+7] * v[i+7][j+21];
            uv[j+22] +=  u[i+0] * v[i+0][j+22]
                       + u[i+1] * v[i+1][j+22]
                       + u[i+2] * v[i+2][j+22]
                       + u[i+3] * v[i+3][j+22]
                       + u[i+4] * v[i+4][j+22]
                       + u[i+5] * v[i+5][j+22]
                       + u[i+6] * v[i+6][j+22]
                       + u[i+7] * v[i+7][j+22];
            uv[j+23] +=  u[i+0] * v[i+0][j+23]
                       + u[i+1] * v[i+1][j+23]
                       + u[i+2] * v[i+2][j+23]
                       + u[i+3] * v[i+3][j+23]
                       + u[i+4] * v[i+4][j+23]
                       + u[i+5] * v[i+5][j+23]
                       + u[i+6] * v[i+6][j+23]
                       + u[i+7] * v[i+7][j+23];
            uv[j+24] +=  u[i+0] * v[i+0][j+24]
                       + u[i+1] * v[i+1][j+24]
                       + u[i+2] * v[i+2][j+24]
                       + u[i+3] * v[i+3][j+24]
                       + u[i+4] * v[i+4][j+24]
                       + u[i+5] * v[i+5][j+24]
                       + u[i+6] * v[i+6][j+24]
                       + u[i+7] * v[i+7][j+24];
            uv[j+25] +=  u[i+0] * v[i+0][j+25]
                       + u[i+1] * v[i+1][j+25]
                       + u[i+2] * v[i+2][j+25]
                       + u[i+3] * v[i+3][j+25]
                       + u[i+4] * v[i+4][j+25]
                       + u[i+5] * v[i+5][j+25]
                       + u[i+6] * v[i+6][j+25]
                       + u[i+7] * v[i+7][j+25];
            uv[j+26] +=  u[i+0] * v[i+0][j+26]
                       + u[i+1] * v[i+1][j+26]
                       + u[i+2] * v[i+2][j+26]
                       + u[i+3] * v[i+3][j+26]
                       + u[i+4] * v[i+4][j+26]
                       + u[i+5] * v[i+5][j+26]
                       + u[i+6] * v[i+6][j+26]
                       + u[i+7] * v[i+7][j+26];
            uv[j+27] +=  u[i+0] * v[i+0][j+27]
                       + u[i+1] * v[i+1][j+27]
                       + u[i+2] * v[i+2][j+27]
                       + u[i+3] * v[i+3][j+27]
                       + u[i+4] * v[i+4][j+27]
                       + u[i+5] * v[i+5][j+27]
                       + u[i+6] * v[i+6][j+27]
                       + u[i+7] * v[i+7][j+27];
            uv[j+28] +=  u[i+0] * v[i+0][j+28]
                       + u[i+1] * v[i+1][j+28]
                       + u[i+2] * v[i+2][j+28]
                       + u[i+3] * v[i+3][j+28]
                       + u[i+4] * v[i+4][j+28]
                       + u[i+5] * v[i+5][j+28]
                       + u[i+6] * v[i+6][j+28]
                       + u[i+7] * v[i+7][j+28];
            uv[j+29] +=  u[i+0] * v[i+0][j+29]
                       + u[i+1] * v[i+1][j+29]
                       + u[i+2] * v[i+2][j+29]
                       + u[i+3] * v[i+3][j+29]
                       + u[i+4] * v[i+4][j+29]
                       + u[i+5] * v[i+5][j+29]
                       + u[i+6] * v[i+6][j+29]
                       + u[i+7] * v[i+7][j+29];
            uv[j+30] +=  u[i+0] * v[i+0][j+30]
                       + u[i+1] * v[i+1][j+30]
                       + u[i+2] * v[i+2][j+30]
                       + u[i+3] * v[i+3][j+30]
                       + u[i+4] * v[i+4][j+30]
                       + u[i+5] * v[i+5][j+30]
                       + u[i+6] * v[i+6][j+30]
                       + u[i+7] * v[i+7][j+30];
            uv[j+31] +=  u[i+0] * v[i+0][j+31]
                       + u[i+1] * v[i+1][j+31]
                       + u[i+2] * v[i+2][j+31]
                       + u[i+3] * v[i+3][j+31]
                       + u[i+4] * v[i+4][j+31]
                       + u[i+5] * v[i+5][j+31]
                       + u[i+6] * v[i+6][j+31]
                       + u[i+7] * v[i+7][j+31];
        }
        for (; j < k; j++)
        {
            uv[j] +=   u[i+0] * v[i+0][j]
                     + u[i+1] * v[i+1][j]
                     + u[i+2] * v[i+2][j]
                     + u[i+3] * v[i+3][j]
                     + u[i+4] * v[i+4][j]
                     + u[i+5] * v[i+5][j]
                     + u[i+6] * v[i+6][j]
                     + u[i+7] * v[i+7][j];
        }
    }
    for (; i < len; i++)
    {
        ulong j = 0;
        for (; j+31 < k; j += 32)
        {
            uv[j+0] += u[i] * v[i][j+0];
            uv[j+1] += u[i] * v[i][j+1];
            uv[j+2] += u[i] * v[i][j+2];
            uv[j+3] += u[i] * v[i][j+3];
            uv[j+4] += u[i] * v[i][j+4];
            uv[j+5] += u[i] * v[i][j+5];
            uv[j+6] += u[i] * v[i][j+6];
            uv[j+7] += u[i] * v[i][j+7];
            uv[j+8] += u[i] * v[i][j+8];
            uv[j+9] += u[i] * v[i][j+9];
            uv[j+10] += u[i] * v[i][j+10];
            uv[j+11] += u[i] * v[i][j+11];
            uv[j+12] += u[i] * v[i][j+12];
            uv[j+13] += u[i] * v[i][j+13];
            uv[j+14] += u[i] * v[i][j+14];
            uv[j+15] += u[i] * v[i][j+15];
            uv[j+16] += u[i] * v[i][j+16];
            uv[j+17] += u[i] * v[i][j+17];
            uv[j+18] += u[i] * v[i][j+18];
            uv[j+19] += u[i] * v[i][j+19];
            uv[j+20] += u[i] * v[i][j+20];
            uv[j+21] += u[i] * v[i][j+21];
            uv[j+22] += u[i] * v[i][j+22];
            uv[j+23] += u[i] * v[i][j+23];
            uv[j+24] += u[i] * v[i][j+24];
            uv[j+25] += u[i] * v[i][j+25];
            uv[j+26] += u[i] * v[i][j+26];
            uv[j+27] += u[i] * v[i][j+27];
            uv[j+28] += u[i] * v[i][j+28];
            uv[j+29] += u[i] * v[i][j+29];
            uv[j+30] += u[i] * v[i][j+30];
            uv[j+31] += u[i] * v[i][j+31];
        }
        for (; j < k; j++)
            uv[j] += u[i] * v[i][j];
    }

    for (ulong j = 0; j < k; j++)
        NMOD_RED(uv[j], uv[j], mod);
}

void _nmod_vec_dot_product_multi_1_v16_16(mp_ptr uv, mp_srcptr u, mp_srcptr * v,
                                      ulong len, ulong k, nmod_t mod)
{
    _nmod_vec_zero(uv, k);
    ulong i = 0;
    for (; i+15 < len; i += 16)
    {
        ulong j = 0;
        for (; j+15 < k; j += 16)
        {
            uv[j+0] +=   u[i+ 0] * v[i+ 0][j+0]
                       + u[i+ 1] * v[i+ 1][j+0]
                       + u[i+ 2] * v[i+ 2][j+0]
                       + u[i+ 3] * v[i+ 3][j+0]
                       + u[i+ 4] * v[i+ 4][j+0]
                       + u[i+ 5] * v[i+ 5][j+0]
                       + u[i+ 6] * v[i+ 6][j+0]
                       + u[i+ 7] * v[i+ 7][j+0]
                       + u[i+ 8] * v[i+ 8][j+0]
                       + u[i+ 9] * v[i+ 9][j+0]
                       + u[i+10] * v[i+10][j+0]
                       + u[i+11] * v[i+11][j+0]
                       + u[i+12] * v[i+12][j+0]
                       + u[i+13] * v[i+13][j+0]
                       + u[i+14] * v[i+14][j+0]
                       + u[i+15] * v[i+15][j+0];
            uv[j+1] +=   u[i+ 0] * v[i+ 0][j+1]
                       + u[i+ 1] * v[i+ 1][j+1]
                       + u[i+ 2] * v[i+ 2][j+1]
                       + u[i+ 3] * v[i+ 3][j+1]
                       + u[i+ 4] * v[i+ 4][j+1]
                       + u[i+ 5] * v[i+ 5][j+1]
                       + u[i+ 6] * v[i+ 6][j+1]
                       + u[i+ 7] * v[i+ 7][j+1]
                       + u[i+ 8] * v[i+ 8][j+1]
                       + u[i+ 9] * v[i+ 9][j+1]
                       + u[i+10] * v[i+10][j+1]
                       + u[i+11] * v[i+11][j+1]
                       + u[i+12] * v[i+12][j+1]
                       + u[i+13] * v[i+13][j+1]
                       + u[i+14] * v[i+14][j+1]
                       + u[i+15] * v[i+15][j+1];
            uv[j+2] +=   u[i+ 0] * v[i+ 0][j+2]
                       + u[i+ 1] * v[i+ 1][j+2]
                       + u[i+ 2] * v[i+ 2][j+2]
                       + u[i+ 3] * v[i+ 3][j+2]
                       + u[i+ 4] * v[i+ 4][j+2]
                       + u[i+ 5] * v[i+ 5][j+2]
                       + u[i+ 6] * v[i+ 6][j+2]
                       + u[i+ 7] * v[i+ 7][j+2]
                       + u[i+ 8] * v[i+ 8][j+2]
                       + u[i+ 9] * v[i+ 9][j+2]
                       + u[i+10] * v[i+10][j+2]
                       + u[i+11] * v[i+11][j+2]
                       + u[i+12] * v[i+12][j+2]
                       + u[i+13] * v[i+13][j+2]
                       + u[i+14] * v[i+14][j+2]
                       + u[i+15] * v[i+15][j+2];
            uv[j+3] +=   u[i+ 0] * v[i+ 0][j+3]
                       + u[i+ 1] * v[i+ 1][j+3]
                       + u[i+ 2] * v[i+ 2][j+3]
                       + u[i+ 3] * v[i+ 3][j+3]
                       + u[i+ 4] * v[i+ 4][j+3]
                       + u[i+ 5] * v[i+ 5][j+3]
                       + u[i+ 6] * v[i+ 6][j+3]
                       + u[i+ 7] * v[i+ 7][j+3]
                       + u[i+ 8] * v[i+ 8][j+3]
                       + u[i+ 9] * v[i+ 9][j+3]
                       + u[i+10] * v[i+10][j+3]
                       + u[i+11] * v[i+11][j+3]
                       + u[i+12] * v[i+12][j+3]
                       + u[i+13] * v[i+13][j+3]
                       + u[i+14] * v[i+14][j+3]
                       + u[i+15] * v[i+15][j+3];
            uv[j+4] +=   u[i+ 0] * v[i+ 0][j+4]
                       + u[i+ 1] * v[i+ 1][j+4]
                       + u[i+ 2] * v[i+ 2][j+4]
                       + u[i+ 3] * v[i+ 3][j+4]
                       + u[i+ 4] * v[i+ 4][j+4]
                       + u[i+ 5] * v[i+ 5][j+4]
                       + u[i+ 6] * v[i+ 6][j+4]
                       + u[i+ 7] * v[i+ 7][j+4]
                       + u[i+ 8] * v[i+ 8][j+4]
                       + u[i+ 9] * v[i+ 9][j+4]
                       + u[i+10] * v[i+10][j+4]
                       + u[i+11] * v[i+11][j+4]
                       + u[i+12] * v[i+12][j+4]
                       + u[i+13] * v[i+13][j+4]
                       + u[i+14] * v[i+14][j+4]
                       + u[i+15] * v[i+15][j+4];
            uv[j+5] +=   u[i+ 0] * v[i+ 0][j+5]
                       + u[i+ 1] * v[i+ 1][j+5]
                       + u[i+ 2] * v[i+ 2][j+5]
                       + u[i+ 3] * v[i+ 3][j+5]
                       + u[i+ 4] * v[i+ 4][j+5]
                       + u[i+ 5] * v[i+ 5][j+5]
                       + u[i+ 6] * v[i+ 6][j+5]
                       + u[i+ 7] * v[i+ 7][j+5]
                       + u[i+ 8] * v[i+ 8][j+5]
                       + u[i+ 9] * v[i+ 9][j+5]
                       + u[i+10] * v[i+10][j+5]
                       + u[i+11] * v[i+11][j+5]
                       + u[i+12] * v[i+12][j+5]
                       + u[i+13] * v[i+13][j+5]
                       + u[i+14] * v[i+14][j+5]
                       + u[i+15] * v[i+15][j+5];
            uv[j+6] +=   u[i+ 0] * v[i+ 0][j+6]
                       + u[i+ 1] * v[i+ 1][j+6]
                       + u[i+ 2] * v[i+ 2][j+6]
                       + u[i+ 3] * v[i+ 3][j+6]
                       + u[i+ 4] * v[i+ 4][j+6]
                       + u[i+ 5] * v[i+ 5][j+6]
                       + u[i+ 6] * v[i+ 6][j+6]
                       + u[i+ 7] * v[i+ 7][j+6]
                       + u[i+ 8] * v[i+ 8][j+6]
                       + u[i+ 9] * v[i+ 9][j+6]
                       + u[i+10] * v[i+10][j+6]
                       + u[i+11] * v[i+11][j+6]
                       + u[i+12] * v[i+12][j+6]
                       + u[i+13] * v[i+13][j+6]
                       + u[i+14] * v[i+14][j+6]
                       + u[i+15] * v[i+15][j+6];
            uv[j+7] +=   u[i+ 0] * v[i+ 0][j+7]
                       + u[i+ 1] * v[i+ 1][j+7]
                       + u[i+ 2] * v[i+ 2][j+7]
                       + u[i+ 3] * v[i+ 3][j+7]
                       + u[i+ 4] * v[i+ 4][j+7]
                       + u[i+ 5] * v[i+ 5][j+7]
                       + u[i+ 6] * v[i+ 6][j+7]
                       + u[i+ 7] * v[i+ 7][j+7]
                       + u[i+ 8] * v[i+ 8][j+7]
                       + u[i+ 9] * v[i+ 9][j+7]
                       + u[i+10] * v[i+10][j+7]
                       + u[i+11] * v[i+11][j+7]
                       + u[i+12] * v[i+12][j+7]
                       + u[i+13] * v[i+13][j+7]
                       + u[i+14] * v[i+14][j+7]
                       + u[i+15] * v[i+15][j+7];
            uv[j+8] +=   u[i+ 0] * v[i+ 0][j+8]
                       + u[i+ 1] * v[i+ 1][j+8]
                       + u[i+ 2] * v[i+ 2][j+8]
                       + u[i+ 3] * v[i+ 3][j+8]
                       + u[i+ 4] * v[i+ 4][j+8]
                       + u[i+ 5] * v[i+ 5][j+8]
                       + u[i+ 6] * v[i+ 6][j+8]
                       + u[i+ 7] * v[i+ 7][j+8]
                       + u[i+ 8] * v[i+ 8][j+8]
                       + u[i+ 9] * v[i+ 9][j+8]
                       + u[i+10] * v[i+10][j+8]
                       + u[i+11] * v[i+11][j+8]
                       + u[i+12] * v[i+12][j+8]
                       + u[i+13] * v[i+13][j+8]
                       + u[i+14] * v[i+14][j+8]
                       + u[i+15] * v[i+15][j+8];
            uv[j+9] +=   u[i+ 0] * v[i+ 0][j+9]
                       + u[i+ 1] * v[i+ 1][j+9]
                       + u[i+ 2] * v[i+ 2][j+9]
                       + u[i+ 3] * v[i+ 3][j+9]
                       + u[i+ 4] * v[i+ 4][j+9]
                       + u[i+ 5] * v[i+ 5][j+9]
                       + u[i+ 6] * v[i+ 6][j+9]
                       + u[i+ 7] * v[i+ 7][j+9]
                       + u[i+ 8] * v[i+ 8][j+9]
                       + u[i+ 9] * v[i+ 9][j+9]
                       + u[i+10] * v[i+10][j+9]
                       + u[i+11] * v[i+11][j+9]
                       + u[i+12] * v[i+12][j+9]
                       + u[i+13] * v[i+13][j+9]
                       + u[i+14] * v[i+14][j+9]
                       + u[i+15] * v[i+15][j+9];
            uv[j+10] +=  u[i+ 0] * v[i+ 0][j+10]
                       + u[i+ 1] * v[i+ 1][j+10]
                       + u[i+ 2] * v[i+ 2][j+10]
                       + u[i+ 3] * v[i+ 3][j+10]
                       + u[i+ 4] * v[i+ 4][j+10]
                       + u[i+ 5] * v[i+ 5][j+10]
                       + u[i+ 6] * v[i+ 6][j+10]
                       + u[i+ 7] * v[i+ 7][j+10]
                       + u[i+ 8] * v[i+ 8][j+10]
                       + u[i+ 9] * v[i+ 9][j+10]
                       + u[i+10] * v[i+10][j+10]
                       + u[i+11] * v[i+11][j+10]
                       + u[i+12] * v[i+12][j+10]
                       + u[i+13] * v[i+13][j+10]
                       + u[i+14] * v[i+14][j+10]
                       + u[i+15] * v[i+15][j+10];
            uv[j+11] +=  u[i+ 0] * v[i+ 0][j+11]
                       + u[i+ 1] * v[i+ 1][j+11]
                       + u[i+ 2] * v[i+ 2][j+11]
                       + u[i+ 3] * v[i+ 3][j+11]
                       + u[i+ 4] * v[i+ 4][j+11]
                       + u[i+ 5] * v[i+ 5][j+11]
                       + u[i+ 6] * v[i+ 6][j+11]
                       + u[i+ 7] * v[i+ 7][j+11]
                       + u[i+ 8] * v[i+ 8][j+11]
                       + u[i+ 9] * v[i+ 9][j+11]
                       + u[i+10] * v[i+10][j+11]
                       + u[i+11] * v[i+11][j+11]
                       + u[i+12] * v[i+12][j+11]
                       + u[i+13] * v[i+13][j+11]
                       + u[i+14] * v[i+14][j+11]
                       + u[i+15] * v[i+15][j+11];
            uv[j+12] +=  u[i+ 0] * v[i+ 0][j+12]
                       + u[i+ 1] * v[i+ 1][j+12]
                       + u[i+ 2] * v[i+ 2][j+12]
                       + u[i+ 3] * v[i+ 3][j+12]
                       + u[i+ 4] * v[i+ 4][j+12]
                       + u[i+ 5] * v[i+ 5][j+12]
                       + u[i+ 6] * v[i+ 6][j+12]
                       + u[i+ 7] * v[i+ 7][j+12]
                       + u[i+ 8] * v[i+ 8][j+12]
                       + u[i+ 9] * v[i+ 9][j+12]
                       + u[i+10] * v[i+10][j+12]
                       + u[i+11] * v[i+11][j+12]
                       + u[i+12] * v[i+12][j+12]
                       + u[i+13] * v[i+13][j+12]
                       + u[i+14] * v[i+14][j+12]
                       + u[i+15] * v[i+15][j+12];
            uv[j+13] +=  u[i+ 0] * v[i+ 0][j+13]
                       + u[i+ 1] * v[i+ 1][j+13]
                       + u[i+ 2] * v[i+ 2][j+13]
                       + u[i+ 3] * v[i+ 3][j+13]
                       + u[i+ 4] * v[i+ 4][j+13]
                       + u[i+ 5] * v[i+ 5][j+13]
                       + u[i+ 6] * v[i+ 6][j+13]
                       + u[i+ 7] * v[i+ 7][j+13]
                       + u[i+ 8] * v[i+ 8][j+13]
                       + u[i+ 9] * v[i+ 9][j+13]
                       + u[i+10] * v[i+10][j+13]
                       + u[i+11] * v[i+11][j+13]
                       + u[i+12] * v[i+12][j+13]
                       + u[i+13] * v[i+13][j+13]
                       + u[i+14] * v[i+14][j+13]
                       + u[i+15] * v[i+15][j+13];
            uv[j+14] +=  u[i+ 0] * v[i+ 0][j+14]
                       + u[i+ 1] * v[i+ 1][j+14]
                       + u[i+ 2] * v[i+ 2][j+14]
                       + u[i+ 3] * v[i+ 3][j+14]
                       + u[i+ 4] * v[i+ 4][j+14]
                       + u[i+ 5] * v[i+ 5][j+14]
                       + u[i+ 6] * v[i+ 6][j+14]
                       + u[i+ 7] * v[i+ 7][j+14]
                       + u[i+ 8] * v[i+ 8][j+14]
                       + u[i+ 9] * v[i+ 9][j+14]
                       + u[i+10] * v[i+10][j+14]
                       + u[i+11] * v[i+11][j+14]
                       + u[i+12] * v[i+12][j+14]
                       + u[i+13] * v[i+13][j+14]
                       + u[i+14] * v[i+14][j+14]
                       + u[i+15] * v[i+15][j+14];
            uv[j+15] +=  u[i+ 0] * v[i+ 0][j+15]
                       + u[i+ 1] * v[i+ 1][j+15]
                       + u[i+ 2] * v[i+ 2][j+15]
                       + u[i+ 3] * v[i+ 3][j+15]
                       + u[i+ 4] * v[i+ 4][j+15]
                       + u[i+ 5] * v[i+ 5][j+15]
                       + u[i+ 6] * v[i+ 6][j+15]
                       + u[i+ 7] * v[i+ 7][j+15]
                       + u[i+ 8] * v[i+ 8][j+15]
                       + u[i+ 9] * v[i+ 9][j+15]
                       + u[i+10] * v[i+10][j+15]
                       + u[i+11] * v[i+11][j+15]
                       + u[i+12] * v[i+12][j+15]
                       + u[i+13] * v[i+13][j+15]
                       + u[i+14] * v[i+14][j+15]
                       + u[i+15] * v[i+15][j+15];
        }
        for (; j < k; j++)
        {
            uv[j] +=   u[i+ 0] * v[i+ 0][j]
                     + u[i+ 1] * v[i+ 1][j]
                     + u[i+ 2] * v[i+ 2][j]
                     + u[i+ 3] * v[i+ 3][j]
                     + u[i+ 4] * v[i+ 4][j]
                     + u[i+ 5] * v[i+ 5][j]
                     + u[i+ 6] * v[i+ 6][j]
                     + u[i+ 7] * v[i+ 7][j]
                     + u[i+ 8] * v[i+ 8][j]
                     + u[i+ 9] * v[i+ 9][j]
                     + u[i+10] * v[i+10][j]
                     + u[i+11] * v[i+11][j]
                     + u[i+12] * v[i+12][j]
                     + u[i+13] * v[i+13][j]
                     + u[i+14] * v[i+14][j]
                     + u[i+15] * v[i+15][j];
        }
    }
    for (; i < len; i++)
    {
        ulong j = 0;
        for (; j+15 < k; j += 16)
        {
            uv[j+ 0] += u[i] * v[i][j+ 0];
            uv[j+ 1] += u[i] * v[i][j+ 1];
            uv[j+ 2] += u[i] * v[i][j+ 2];
            uv[j+ 3] += u[i] * v[i][j+ 3];
            uv[j+ 4] += u[i] * v[i][j+ 4];
            uv[j+ 5] += u[i] * v[i][j+ 5];
            uv[j+ 6] += u[i] * v[i][j+ 6];
            uv[j+ 7] += u[i] * v[i][j+ 7];
            uv[j+ 8] += u[i] * v[i][j+ 8];
            uv[j+ 9] += u[i] * v[i][j+ 9];
            uv[j+10] += u[i] * v[i][j+10];
            uv[j+11] += u[i] * v[i][j+11];
            uv[j+12] += u[i] * v[i][j+12];
            uv[j+13] += u[i] * v[i][j+13];
            uv[j+14] += u[i] * v[i][j+14];
            uv[j+15] += u[i] * v[i][j+15];
        }
        for (; j < k; j++)
            uv[j] += u[i] * v[i][j];
    }

    for (ulong j = 0; j < k; j++)
        NMOD_RED(uv[j], uv[j], mod);
}

static inline
void _nmod_vec_dot_product_multi_2(mp_ptr uv, mp_srcptr u, mp_srcptr * v,
                                   ulong len, ulong k, nmod_t mod)
{
    mp_ptr uv_hi = _nmod_vec_init(k);
    _nmod_vec_zero(uv, k);
    _nmod_vec_zero(uv_hi, k);
    mp_limb_t s0, s1;
    for (ulong i = 0; i < len; i++)
    {
        for (ulong j = 0; j < k; j++)
        {
            umul_ppmm(s1, s0, u[i], v[i][j]);
            add_ssaaaa(uv_hi[j], uv[j], uv_hi[j], uv[j], s1, s0);
        }
    }

    for (ulong j = 0; j < k; j++)
        NMOD2_RED2(uv[j], uv_hi[j], uv[j], mod);
    _nmod_vec_clear(uv_hi);
}

// TODO benchmark more, integrate, give precise conditions for when this works
// (or better, really do a hand-made avx512 version...)
// --> if splitting at 26, each product is 52, can allow at most 12 additional bits,
// i.e. not more than xxx terms (this depends on the size of the high part since
// they are not balanced... could make sense to balance them to allow more terms
// (info from max_bits_u and max_bits_v could be useful), but do this only if
// this really is interesting in terms of speed)
void _nmod_vec_dot_product_multi_2_split26(mp_ptr uv, mp_srcptr u, mp_srcptr * v,
                                   ulong len, ulong k, nmod_t mod)
{
    uint ulo, uhi, vlo, vhi;
    mp_ptr uv_mi = _nmod_vec_init(k);
    mp_ptr uv_hi = _nmod_vec_init(k);
    _nmod_vec_zero(uv, k);  // plays the role of uv_lo
    _nmod_vec_zero(uv_mi, k);
    _nmod_vec_zero(uv_hi, k);
    for (ulong i = 0; i < len; i++)
    {
        __ll_lowhi_parts26(ulo, uhi, u[i]);
        for (ulong j = 0; j < k; j++)
        {
            __ll_lowhi_parts26(vlo, vhi, v[i][j]);
            uv   [j] += (ulong)ulo * vlo;
            uv_mi[j] += (ulong)ulo * vhi + (ulong)uhi * vlo;
            uv_hi[j] += (ulong)uhi * vhi;
        }
    }

    for (ulong j = 0; j < k; j++)
    {
        // result: uv_lo + 2**26 uv_mi + 2**52 uv_hi
        // hi = (uv_mi >> 38) + (uv_hi >> 12)  ||  lo = (uv_mi << 26) + (uv_hi << 52) + uv_lo
        add_ssaaaa(uv_hi[j], uv[j], uv_mi[j]>>38, uv_mi[j]<<26, uv_hi[j]>>12, (uv_hi[j]<<52)+uv[j]);
        NMOD2_RED2(uv[j], uv_hi[j], uv[j], mod);
    }
    _nmod_vec_clear(uv_hi);
    _nmod_vec_clear(uv_mi);
}

void _nmod_vec_dot_product_multi_2_v1_8(mp_ptr uv, mp_srcptr u, mp_srcptr * v,
                                   ulong len, ulong k, nmod_t mod)
{
    mp_ptr uv_hi = _nmod_vec_init(k);
    _nmod_vec_zero(uv, k);
    _nmod_vec_zero(uv_hi, k);
    mp_limb_t s0, s1;
    for (ulong i = 0; i < len; i++)
    {
        ulong j = 0;
        for (; j+7 < k; j+=8)
        {
            umul_ppmm(s1, s0, u[i], v[i][j+0]);
            add_ssaaaa(uv_hi[j+0], uv[j+0], uv_hi[j+0], uv[j+0], s1, s0);
            umul_ppmm(s1, s0, u[i], v[i][j+1]);
            add_ssaaaa(uv_hi[j+1], uv[j+1], uv_hi[j+1], uv[j+1], s1, s0);
            umul_ppmm(s1, s0, u[i], v[i][j+2]);
            add_ssaaaa(uv_hi[j+2], uv[j+2], uv_hi[j+2], uv[j+2], s1, s0);
            umul_ppmm(s1, s0, u[i], v[i][j+3]);
            add_ssaaaa(uv_hi[j+3], uv[j+3], uv_hi[j+3], uv[j+3], s1, s0);
            umul_ppmm(s1, s0, u[i], v[i][j+4]);
            add_ssaaaa(uv_hi[j+4], uv[j+4], uv_hi[j+4], uv[j+4], s1, s0);
            umul_ppmm(s1, s0, u[i], v[i][j+5]);
            add_ssaaaa(uv_hi[j+5], uv[j+5], uv_hi[j+5], uv[j+5], s1, s0);
            umul_ppmm(s1, s0, u[i], v[i][j+6]);
            add_ssaaaa(uv_hi[j+6], uv[j+6], uv_hi[j+6], uv[j+6], s1, s0);
            umul_ppmm(s1, s0, u[i], v[i][j+7]);
            add_ssaaaa(uv_hi[j+7], uv[j+7], uv_hi[j+7], uv[j+7], s1, s0);
        }
        for (; j < k; j++)
        {
            umul_ppmm(s1, s0, u[i], v[i][j+0]);
            add_ssaaaa(uv_hi[j+0], uv[j+0], uv_hi[j+0], uv[j+0], s1, s0);
        }
    }

    for (ulong j = 0; j < k; j++)
        NMOD2_RED2(uv[j], uv_hi[j], uv[j], mod);
    _nmod_vec_clear(uv_hi);
}

void _nmod_vec_dot_product_multi_2_v4_8(mp_ptr uv, mp_srcptr u, mp_srcptr * v,
                                   ulong len, ulong k, nmod_t mod)
{
    mp_ptr uv_hi = _nmod_vec_init(k);
    _nmod_vec_zero(uv, k);
    _nmod_vec_zero(uv_hi, k);
    mp_limb_t s0, s1;
    ulong i = 0;
    for (; i+3 < len; i += 4)
    {
        ulong j = 0;
        for (; j+7 < k; j+=8)
        {
            umul_ppmm(s1, s0, u[i], v[i][j+0]);
            add_ssaaaa(uv_hi[j+0], uv[j+0], uv_hi[j+0], uv[j+0], s1, s0);
            umul_ppmm(s1, s0, u[i+1], v[i+1][j+0]);
            add_ssaaaa(uv_hi[j+0], uv[j+0], uv_hi[j+0], uv[j+0], s1, s0);
            umul_ppmm(s1, s0, u[i+2], v[i+2][j+0]);
            add_ssaaaa(uv_hi[j+0], uv[j+0], uv_hi[j+0], uv[j+0], s1, s0);
            umul_ppmm(s1, s0, u[i+3], v[i+3][j+0]);
            add_ssaaaa(uv_hi[j+0], uv[j+0], uv_hi[j+0], uv[j+0], s1, s0);

            umul_ppmm(s1, s0, u[i], v[i][j+1]);
            add_ssaaaa(uv_hi[j+1], uv[j+1], uv_hi[j+1], uv[j+1], s1, s0);
            umul_ppmm(s1, s0, u[i+1], v[i+1][j+1]);
            add_ssaaaa(uv_hi[j+1], uv[j+1], uv_hi[j+1], uv[j+1], s1, s0);
            umul_ppmm(s1, s0, u[i+2], v[i+2][j+1]);
            add_ssaaaa(uv_hi[j+1], uv[j+1], uv_hi[j+1], uv[j+1], s1, s0);
            umul_ppmm(s1, s0, u[i+3], v[i+3][j+1]);
            add_ssaaaa(uv_hi[j+1], uv[j+1], uv_hi[j+1], uv[j+1], s1, s0);

            umul_ppmm(s1, s0, u[i], v[i][j+2]);
            add_ssaaaa(uv_hi[j+2], uv[j+2], uv_hi[j+2], uv[j+2], s1, s0);
            umul_ppmm(s1, s0, u[i+1], v[i+1][j+2]);
            add_ssaaaa(uv_hi[j+2], uv[j+2], uv_hi[j+2], uv[j+2], s1, s0);
            umul_ppmm(s1, s0, u[i+2], v[i+2][j+2]);
            add_ssaaaa(uv_hi[j+2], uv[j+2], uv_hi[j+2], uv[j+2], s1, s0);
            umul_ppmm(s1, s0, u[i+3], v[i+3][j+2]);
            add_ssaaaa(uv_hi[j+2], uv[j+2], uv_hi[j+2], uv[j+2], s1, s0);

            umul_ppmm(s1, s0, u[i], v[i][j+3]);
            add_ssaaaa(uv_hi[j+3], uv[j+3], uv_hi[j+3], uv[j+3], s1, s0);
            umul_ppmm(s1, s0, u[i+1], v[i+1][j+3]);
            add_ssaaaa(uv_hi[j+3], uv[j+3], uv_hi[j+3], uv[j+3], s1, s0);
            umul_ppmm(s1, s0, u[i+2], v[i+2][j+3]);
            add_ssaaaa(uv_hi[j+3], uv[j+3], uv_hi[j+3], uv[j+3], s1, s0);
            umul_ppmm(s1, s0, u[i+3], v[i+3][j+3]);
            add_ssaaaa(uv_hi[j+3], uv[j+3], uv_hi[j+3], uv[j+3], s1, s0);

            umul_ppmm(s1, s0, u[i], v[i][j+4]);
            add_ssaaaa(uv_hi[j+4], uv[j+4], uv_hi[j+4], uv[j+4], s1, s0);
            umul_ppmm(s1, s0, u[i+1], v[i+1][j+4]);
            add_ssaaaa(uv_hi[j+4], uv[j+4], uv_hi[j+4], uv[j+4], s1, s0);
            umul_ppmm(s1, s0, u[i+2], v[i+2][j+4]);
            add_ssaaaa(uv_hi[j+4], uv[j+4], uv_hi[j+4], uv[j+4], s1, s0);
            umul_ppmm(s1, s0, u[i+3], v[i+3][j+4]);
            add_ssaaaa(uv_hi[j+4], uv[j+4], uv_hi[j+4], uv[j+4], s1, s0);

            umul_ppmm(s1, s0, u[i], v[i][j+5]);
            add_ssaaaa(uv_hi[j+5], uv[j+5], uv_hi[j+5], uv[j+5], s1, s0);
            umul_ppmm(s1, s0, u[i+1], v[i+1][j+5]);
            add_ssaaaa(uv_hi[j+5], uv[j+5], uv_hi[j+5], uv[j+5], s1, s0);
            umul_ppmm(s1, s0, u[i+2], v[i+2][j+5]);
            add_ssaaaa(uv_hi[j+5], uv[j+5], uv_hi[j+5], uv[j+5], s1, s0);
            umul_ppmm(s1, s0, u[i+3], v[i+3][j+5]);
            add_ssaaaa(uv_hi[j+5], uv[j+5], uv_hi[j+5], uv[j+5], s1, s0);

            umul_ppmm(s1, s0, u[i], v[i][j+6]);
            add_ssaaaa(uv_hi[j+6], uv[j+6], uv_hi[j+6], uv[j+6], s1, s0);
            umul_ppmm(s1, s0, u[i+1], v[i+1][j+6]);
            add_ssaaaa(uv_hi[j+6], uv[j+6], uv_hi[j+6], uv[j+6], s1, s0);
            umul_ppmm(s1, s0, u[i+2], v[i+2][j+6]);
            add_ssaaaa(uv_hi[j+6], uv[j+6], uv_hi[j+6], uv[j+6], s1, s0);
            umul_ppmm(s1, s0, u[i+3], v[i+3][j+6]);
            add_ssaaaa(uv_hi[j+6], uv[j+6], uv_hi[j+6], uv[j+6], s1, s0);

            umul_ppmm(s1, s0, u[i], v[i][j+7]);
            add_ssaaaa(uv_hi[j+7], uv[j+7], uv_hi[j+7], uv[j+7], s1, s0);
            umul_ppmm(s1, s0, u[i+1], v[i+1][j+7]);
            add_ssaaaa(uv_hi[j+7], uv[j+7], uv_hi[j+7], uv[j+7], s1, s0);
            umul_ppmm(s1, s0, u[i+2], v[i+2][j+7]);
            add_ssaaaa(uv_hi[j+7], uv[j+7], uv_hi[j+7], uv[j+7], s1, s0);
            umul_ppmm(s1, s0, u[i+3], v[i+3][j+7]);
            add_ssaaaa(uv_hi[j+7], uv[j+7], uv_hi[j+7], uv[j+7], s1, s0);
        }
        for (; j < k; j++)
        {
            umul_ppmm(s1, s0, u[i], v[i][j+0]);
            add_ssaaaa(uv_hi[j+0], uv[j+0], uv_hi[j+0], uv[j+0], s1, s0);
            umul_ppmm(s1, s0, u[i+1], v[i+1][j+0]);
            add_ssaaaa(uv_hi[j+0], uv[j+0], uv_hi[j+0], uv[j+0], s1, s0);
            umul_ppmm(s1, s0, u[i+2], v[i+2][j+0]);
            add_ssaaaa(uv_hi[j+0], uv[j+0], uv_hi[j+0], uv[j+0], s1, s0);
            umul_ppmm(s1, s0, u[i+3], v[i+3][j+0]);
            add_ssaaaa(uv_hi[j+0], uv[j+0], uv_hi[j+0], uv[j+0], s1, s0);
        }
    }
    for (; i < len; i++)
    {
        ulong j = 0;
        for (; j+7 < k; j+=8)
        {
            umul_ppmm(s1, s0, u[i], v[i][j+0]);
            add_ssaaaa(uv_hi[j+0], uv[j+0], uv_hi[j+0], uv[j+0], s1, s0);
            umul_ppmm(s1, s0, u[i], v[i][j+1]);
            add_ssaaaa(uv_hi[j+1], uv[j+1], uv_hi[j+1], uv[j+1], s1, s0);
            umul_ppmm(s1, s0, u[i], v[i][j+2]);
            add_ssaaaa(uv_hi[j+2], uv[j+2], uv_hi[j+2], uv[j+2], s1, s0);
            umul_ppmm(s1, s0, u[i], v[i][j+3]);
            add_ssaaaa(uv_hi[j+3], uv[j+3], uv_hi[j+3], uv[j+3], s1, s0);
            umul_ppmm(s1, s0, u[i], v[i][j+4]);
            add_ssaaaa(uv_hi[j+4], uv[j+4], uv_hi[j+4], uv[j+4], s1, s0);
            umul_ppmm(s1, s0, u[i], v[i][j+5]);
            add_ssaaaa(uv_hi[j+5], uv[j+5], uv_hi[j+5], uv[j+5], s1, s0);
            umul_ppmm(s1, s0, u[i], v[i][j+6]);
            add_ssaaaa(uv_hi[j+6], uv[j+6], uv_hi[j+6], uv[j+6], s1, s0);
            umul_ppmm(s1, s0, u[i], v[i][j+7]);
            add_ssaaaa(uv_hi[j+7], uv[j+7], uv_hi[j+7], uv[j+7], s1, s0);
        }
        for (; j < k; j++)
        {
            umul_ppmm(s1, s0, u[i], v[i][j+0]);
            add_ssaaaa(uv_hi[j+0], uv[j+0], uv_hi[j+0], uv[j+0], s1, s0);
        }
    }

    for (ulong j = 0; j < k; j++)
        NMOD2_RED2(uv[j], uv_hi[j], uv[j], mod);
    _nmod_vec_clear(uv_hi);
}

void _nmod_vec_dot_product_multi_2_v4_32(mp_ptr uv, mp_srcptr u, mp_srcptr * v,
                                   ulong len, ulong k, nmod_t mod)
{
    mp_ptr uv_hi = _nmod_vec_init(k);
    _nmod_vec_zero(uv, k);
    _nmod_vec_zero(uv_hi, k);
    mp_limb_t s0, s1;
    ulong i = 0;
    for (; i+3 < len; i += 4)
    {
        ulong j = 0;
        for (; j+31 < k; j+=32)
        {
            umul_ppmm(s1, s0, u[i], v[i][j+0]);
            add_ssaaaa(uv_hi[j+0], uv[j+0], uv_hi[j+0], uv[j+0], s1, s0);
            umul_ppmm(s1, s0, u[i+1], v[i+1][j+0]);
            add_ssaaaa(uv_hi[j+0], uv[j+0], uv_hi[j+0], uv[j+0], s1, s0);
            umul_ppmm(s1, s0, u[i+2], v[i+2][j+0]);
            add_ssaaaa(uv_hi[j+0], uv[j+0], uv_hi[j+0], uv[j+0], s1, s0);
            umul_ppmm(s1, s0, u[i+3], v[i+3][j+0]);
            add_ssaaaa(uv_hi[j+0], uv[j+0], uv_hi[j+0], uv[j+0], s1, s0);

            umul_ppmm(s1, s0, u[i], v[i][j+1]);
            add_ssaaaa(uv_hi[j+1], uv[j+1], uv_hi[j+1], uv[j+1], s1, s0);
            umul_ppmm(s1, s0, u[i+1], v[i+1][j+1]);
            add_ssaaaa(uv_hi[j+1], uv[j+1], uv_hi[j+1], uv[j+1], s1, s0);
            umul_ppmm(s1, s0, u[i+2], v[i+2][j+1]);
            add_ssaaaa(uv_hi[j+1], uv[j+1], uv_hi[j+1], uv[j+1], s1, s0);
            umul_ppmm(s1, s0, u[i+3], v[i+3][j+1]);
            add_ssaaaa(uv_hi[j+1], uv[j+1], uv_hi[j+1], uv[j+1], s1, s0);

            umul_ppmm(s1, s0, u[i], v[i][j+2]);
            add_ssaaaa(uv_hi[j+2], uv[j+2], uv_hi[j+2], uv[j+2], s1, s0);
            umul_ppmm(s1, s0, u[i+1], v[i+1][j+2]);
            add_ssaaaa(uv_hi[j+2], uv[j+2], uv_hi[j+2], uv[j+2], s1, s0);
            umul_ppmm(s1, s0, u[i+2], v[i+2][j+2]);
            add_ssaaaa(uv_hi[j+2], uv[j+2], uv_hi[j+2], uv[j+2], s1, s0);
            umul_ppmm(s1, s0, u[i+3], v[i+3][j+2]);
            add_ssaaaa(uv_hi[j+2], uv[j+2], uv_hi[j+2], uv[j+2], s1, s0);

            umul_ppmm(s1, s0, u[i], v[i][j+3]);
            add_ssaaaa(uv_hi[j+3], uv[j+3], uv_hi[j+3], uv[j+3], s1, s0);
            umul_ppmm(s1, s0, u[i+1], v[i+1][j+3]);
            add_ssaaaa(uv_hi[j+3], uv[j+3], uv_hi[j+3], uv[j+3], s1, s0);
            umul_ppmm(s1, s0, u[i+2], v[i+2][j+3]);
            add_ssaaaa(uv_hi[j+3], uv[j+3], uv_hi[j+3], uv[j+3], s1, s0);
            umul_ppmm(s1, s0, u[i+3], v[i+3][j+3]);
            add_ssaaaa(uv_hi[j+3], uv[j+3], uv_hi[j+3], uv[j+3], s1, s0);

            umul_ppmm(s1, s0, u[i], v[i][j+4]);
            add_ssaaaa(uv_hi[j+4], uv[j+4], uv_hi[j+4], uv[j+4], s1, s0);
            umul_ppmm(s1, s0, u[i+1], v[i+1][j+4]);
            add_ssaaaa(uv_hi[j+4], uv[j+4], uv_hi[j+4], uv[j+4], s1, s0);
            umul_ppmm(s1, s0, u[i+2], v[i+2][j+4]);
            add_ssaaaa(uv_hi[j+4], uv[j+4], uv_hi[j+4], uv[j+4], s1, s0);
            umul_ppmm(s1, s0, u[i+3], v[i+3][j+4]);
            add_ssaaaa(uv_hi[j+4], uv[j+4], uv_hi[j+4], uv[j+4], s1, s0);

            umul_ppmm(s1, s0, u[i], v[i][j+5]);
            add_ssaaaa(uv_hi[j+5], uv[j+5], uv_hi[j+5], uv[j+5], s1, s0);
            umul_ppmm(s1, s0, u[i+1], v[i+1][j+5]);
            add_ssaaaa(uv_hi[j+5], uv[j+5], uv_hi[j+5], uv[j+5], s1, s0);
            umul_ppmm(s1, s0, u[i+2], v[i+2][j+5]);
            add_ssaaaa(uv_hi[j+5], uv[j+5], uv_hi[j+5], uv[j+5], s1, s0);
            umul_ppmm(s1, s0, u[i+3], v[i+3][j+5]);
            add_ssaaaa(uv_hi[j+5], uv[j+5], uv_hi[j+5], uv[j+5], s1, s0);

            umul_ppmm(s1, s0, u[i], v[i][j+6]);
            add_ssaaaa(uv_hi[j+6], uv[j+6], uv_hi[j+6], uv[j+6], s1, s0);
            umul_ppmm(s1, s0, u[i+1], v[i+1][j+6]);
            add_ssaaaa(uv_hi[j+6], uv[j+6], uv_hi[j+6], uv[j+6], s1, s0);
            umul_ppmm(s1, s0, u[i+2], v[i+2][j+6]);
            add_ssaaaa(uv_hi[j+6], uv[j+6], uv_hi[j+6], uv[j+6], s1, s0);
            umul_ppmm(s1, s0, u[i+3], v[i+3][j+6]);
            add_ssaaaa(uv_hi[j+6], uv[j+6], uv_hi[j+6], uv[j+6], s1, s0);

            umul_ppmm(s1, s0, u[i], v[i][j+7]);
            add_ssaaaa(uv_hi[j+7], uv[j+7], uv_hi[j+7], uv[j+7], s1, s0);
            umul_ppmm(s1, s0, u[i+1], v[i+1][j+7]);
            add_ssaaaa(uv_hi[j+7], uv[j+7], uv_hi[j+7], uv[j+7], s1, s0);
            umul_ppmm(s1, s0, u[i+2], v[i+2][j+7]);
            add_ssaaaa(uv_hi[j+7], uv[j+7], uv_hi[j+7], uv[j+7], s1, s0);
            umul_ppmm(s1, s0, u[i+3], v[i+3][j+7]);
            add_ssaaaa(uv_hi[j+7], uv[j+7], uv_hi[j+7], uv[j+7], s1, s0);

            umul_ppmm(s1, s0, u[i], v[i][j+8]);
            add_ssaaaa(uv_hi[j+8], uv[j+8], uv_hi[j+8], uv[j+8], s1, s0);
            umul_ppmm(s1, s0, u[i+1], v[i+1][j+8]);
            add_ssaaaa(uv_hi[j+8], uv[j+8], uv_hi[j+8], uv[j+8], s1, s0);
            umul_ppmm(s1, s0, u[i+2], v[i+2][j+8]);
            add_ssaaaa(uv_hi[j+8], uv[j+8], uv_hi[j+8], uv[j+8], s1, s0);
            umul_ppmm(s1, s0, u[i+3], v[i+3][j+8]);
            add_ssaaaa(uv_hi[j+8], uv[j+8], uv_hi[j+8], uv[j+8], s1, s0);

            umul_ppmm(s1, s0, u[i], v[i][j+9]);
            add_ssaaaa(uv_hi[j+9], uv[j+9], uv_hi[j+9], uv[j+9], s1, s0);
            umul_ppmm(s1, s0, u[i+1], v[i+1][j+9]);
            add_ssaaaa(uv_hi[j+9], uv[j+9], uv_hi[j+9], uv[j+9], s1, s0);
            umul_ppmm(s1, s0, u[i+2], v[i+2][j+9]);
            add_ssaaaa(uv_hi[j+9], uv[j+9], uv_hi[j+9], uv[j+9], s1, s0);
            umul_ppmm(s1, s0, u[i+3], v[i+3][j+9]);
            add_ssaaaa(uv_hi[j+9], uv[j+9], uv_hi[j+9], uv[j+9], s1, s0);

            umul_ppmm(s1, s0, u[i], v[i][j+10]);
            add_ssaaaa(uv_hi[j+10], uv[j+10], uv_hi[j+10], uv[j+10], s1, s0);
            umul_ppmm(s1, s0, u[i+1], v[i+1][j+10]);
            add_ssaaaa(uv_hi[j+10], uv[j+10], uv_hi[j+10], uv[j+10], s1, s0);
            umul_ppmm(s1, s0, u[i+2], v[i+2][j+10]);
            add_ssaaaa(uv_hi[j+10], uv[j+10], uv_hi[j+10], uv[j+10], s1, s0);
            umul_ppmm(s1, s0, u[i+3], v[i+3][j+10]);
            add_ssaaaa(uv_hi[j+10], uv[j+10], uv_hi[j+10], uv[j+10], s1, s0);

            umul_ppmm(s1, s0, u[i], v[i][j+11]);
            add_ssaaaa(uv_hi[j+11], uv[j+11], uv_hi[j+11], uv[j+11], s1, s0);
            umul_ppmm(s1, s0, u[i+1], v[i+1][j+11]);
            add_ssaaaa(uv_hi[j+11], uv[j+11], uv_hi[j+11], uv[j+11], s1, s0);
            umul_ppmm(s1, s0, u[i+2], v[i+2][j+11]);
            add_ssaaaa(uv_hi[j+11], uv[j+11], uv_hi[j+11], uv[j+11], s1, s0);
            umul_ppmm(s1, s0, u[i+3], v[i+3][j+11]);
            add_ssaaaa(uv_hi[j+11], uv[j+11], uv_hi[j+11], uv[j+11], s1, s0);

            umul_ppmm(s1, s0, u[i], v[i][j+12]);
            add_ssaaaa(uv_hi[j+12], uv[j+12], uv_hi[j+12], uv[j+12], s1, s0);
            umul_ppmm(s1, s0, u[i+1], v[i+1][j+12]);
            add_ssaaaa(uv_hi[j+12], uv[j+12], uv_hi[j+12], uv[j+12], s1, s0);
            umul_ppmm(s1, s0, u[i+2], v[i+2][j+12]);
            add_ssaaaa(uv_hi[j+12], uv[j+12], uv_hi[j+12], uv[j+12], s1, s0);
            umul_ppmm(s1, s0, u[i+3], v[i+3][j+12]);
            add_ssaaaa(uv_hi[j+12], uv[j+12], uv_hi[j+12], uv[j+12], s1, s0);

            umul_ppmm(s1, s0, u[i], v[i][j+13]);
            add_ssaaaa(uv_hi[j+13], uv[j+13], uv_hi[j+13], uv[j+13], s1, s0);
            umul_ppmm(s1, s0, u[i+1], v[i+1][j+13]);
            add_ssaaaa(uv_hi[j+13], uv[j+13], uv_hi[j+13], uv[j+13], s1, s0);
            umul_ppmm(s1, s0, u[i+2], v[i+2][j+13]);
            add_ssaaaa(uv_hi[j+13], uv[j+13], uv_hi[j+13], uv[j+13], s1, s0);
            umul_ppmm(s1, s0, u[i+3], v[i+3][j+13]);
            add_ssaaaa(uv_hi[j+13], uv[j+13], uv_hi[j+13], uv[j+13], s1, s0);

            umul_ppmm(s1, s0, u[i], v[i][j+14]);
            add_ssaaaa(uv_hi[j+14], uv[j+14], uv_hi[j+14], uv[j+14], s1, s0);
            umul_ppmm(s1, s0, u[i+1], v[i+1][j+14]);
            add_ssaaaa(uv_hi[j+14], uv[j+14], uv_hi[j+14], uv[j+14], s1, s0);
            umul_ppmm(s1, s0, u[i+2], v[i+2][j+14]);
            add_ssaaaa(uv_hi[j+14], uv[j+14], uv_hi[j+14], uv[j+14], s1, s0);
            umul_ppmm(s1, s0, u[i+3], v[i+3][j+14]);
            add_ssaaaa(uv_hi[j+14], uv[j+14], uv_hi[j+14], uv[j+14], s1, s0);

            umul_ppmm(s1, s0, u[i], v[i][j+15]);
            add_ssaaaa(uv_hi[j+15], uv[j+15], uv_hi[j+15], uv[j+15], s1, s0);
            umul_ppmm(s1, s0, u[i+1], v[i+1][j+15]);
            add_ssaaaa(uv_hi[j+15], uv[j+15], uv_hi[j+15], uv[j+15], s1, s0);
            umul_ppmm(s1, s0, u[i+2], v[i+2][j+15]);
            add_ssaaaa(uv_hi[j+15], uv[j+15], uv_hi[j+15], uv[j+15], s1, s0);
            umul_ppmm(s1, s0, u[i+3], v[i+3][j+15]);
            add_ssaaaa(uv_hi[j+15], uv[j+15], uv_hi[j+15], uv[j+15], s1, s0);

            umul_ppmm(s1, s0, u[i], v[i][j+16]);
            add_ssaaaa(uv_hi[j+16], uv[j+16], uv_hi[j+16], uv[j+16], s1, s0);
            umul_ppmm(s1, s0, u[i+1], v[i+1][j+16]);
            add_ssaaaa(uv_hi[j+16], uv[j+16], uv_hi[j+16], uv[j+16], s1, s0);
            umul_ppmm(s1, s0, u[i+2], v[i+2][j+16]);
            add_ssaaaa(uv_hi[j+16], uv[j+16], uv_hi[j+16], uv[j+16], s1, s0);
            umul_ppmm(s1, s0, u[i+3], v[i+3][j+16]);
            add_ssaaaa(uv_hi[j+16], uv[j+16], uv_hi[j+16], uv[j+16], s1, s0);

            umul_ppmm(s1, s0, u[i], v[i][j+17]);
            add_ssaaaa(uv_hi[j+17], uv[j+17], uv_hi[j+17], uv[j+17], s1, s0);
            umul_ppmm(s1, s0, u[i+1], v[i+1][j+17]);
            add_ssaaaa(uv_hi[j+17], uv[j+17], uv_hi[j+17], uv[j+17], s1, s0);
            umul_ppmm(s1, s0, u[i+2], v[i+2][j+17]);
            add_ssaaaa(uv_hi[j+17], uv[j+17], uv_hi[j+17], uv[j+17], s1, s0);
            umul_ppmm(s1, s0, u[i+3], v[i+3][j+17]);
            add_ssaaaa(uv_hi[j+17], uv[j+17], uv_hi[j+17], uv[j+17], s1, s0);

            umul_ppmm(s1, s0, u[i], v[i][j+18]);
            add_ssaaaa(uv_hi[j+18], uv[j+18], uv_hi[j+18], uv[j+18], s1, s0);
            umul_ppmm(s1, s0, u[i+1], v[i+1][j+18]);
            add_ssaaaa(uv_hi[j+18], uv[j+18], uv_hi[j+18], uv[j+18], s1, s0);
            umul_ppmm(s1, s0, u[i+2], v[i+2][j+18]);
            add_ssaaaa(uv_hi[j+18], uv[j+18], uv_hi[j+18], uv[j+18], s1, s0);
            umul_ppmm(s1, s0, u[i+3], v[i+3][j+18]);
            add_ssaaaa(uv_hi[j+18], uv[j+18], uv_hi[j+18], uv[j+18], s1, s0);

            umul_ppmm(s1, s0, u[i], v[i][j+19]);
            add_ssaaaa(uv_hi[j+19], uv[j+19], uv_hi[j+19], uv[j+19], s1, s0);
            umul_ppmm(s1, s0, u[i+1], v[i+1][j+19]);
            add_ssaaaa(uv_hi[j+19], uv[j+19], uv_hi[j+19], uv[j+19], s1, s0);
            umul_ppmm(s1, s0, u[i+2], v[i+2][j+19]);
            add_ssaaaa(uv_hi[j+19], uv[j+19], uv_hi[j+19], uv[j+19], s1, s0);
            umul_ppmm(s1, s0, u[i+3], v[i+3][j+19]);
            add_ssaaaa(uv_hi[j+19], uv[j+19], uv_hi[j+19], uv[j+19], s1, s0);

            umul_ppmm(s1, s0, u[i], v[i][j+20]);
            add_ssaaaa(uv_hi[j+20], uv[j+20], uv_hi[j+20], uv[j+20], s1, s0);
            umul_ppmm(s1, s0, u[i+1], v[i+1][j+20]);
            add_ssaaaa(uv_hi[j+20], uv[j+20], uv_hi[j+20], uv[j+20], s1, s0);
            umul_ppmm(s1, s0, u[i+2], v[i+2][j+20]);
            add_ssaaaa(uv_hi[j+20], uv[j+20], uv_hi[j+20], uv[j+20], s1, s0);
            umul_ppmm(s1, s0, u[i+3], v[i+3][j+20]);
            add_ssaaaa(uv_hi[j+20], uv[j+20], uv_hi[j+20], uv[j+20], s1, s0);

            umul_ppmm(s1, s0, u[i], v[i][j+21]);
            add_ssaaaa(uv_hi[j+21], uv[j+21], uv_hi[j+21], uv[j+21], s1, s0);
            umul_ppmm(s1, s0, u[i+1], v[i+1][j+21]);
            add_ssaaaa(uv_hi[j+21], uv[j+21], uv_hi[j+21], uv[j+21], s1, s0);
            umul_ppmm(s1, s0, u[i+2], v[i+2][j+21]);
            add_ssaaaa(uv_hi[j+21], uv[j+21], uv_hi[j+21], uv[j+21], s1, s0);
            umul_ppmm(s1, s0, u[i+3], v[i+3][j+21]);
            add_ssaaaa(uv_hi[j+21], uv[j+21], uv_hi[j+21], uv[j+21], s1, s0);

            umul_ppmm(s1, s0, u[i], v[i][j+22]);
            add_ssaaaa(uv_hi[j+22], uv[j+22], uv_hi[j+22], uv[j+22], s1, s0);
            umul_ppmm(s1, s0, u[i+1], v[i+1][j+22]);
            add_ssaaaa(uv_hi[j+22], uv[j+22], uv_hi[j+22], uv[j+22], s1, s0);
            umul_ppmm(s1, s0, u[i+2], v[i+2][j+22]);
            add_ssaaaa(uv_hi[j+22], uv[j+22], uv_hi[j+22], uv[j+22], s1, s0);
            umul_ppmm(s1, s0, u[i+3], v[i+3][j+22]);
            add_ssaaaa(uv_hi[j+22], uv[j+22], uv_hi[j+22], uv[j+22], s1, s0);

            umul_ppmm(s1, s0, u[i], v[i][j+23]);
            add_ssaaaa(uv_hi[j+23], uv[j+23], uv_hi[j+23], uv[j+23], s1, s0);
            umul_ppmm(s1, s0, u[i+1], v[i+1][j+23]);
            add_ssaaaa(uv_hi[j+23], uv[j+23], uv_hi[j+23], uv[j+23], s1, s0);
            umul_ppmm(s1, s0, u[i+2], v[i+2][j+23]);
            add_ssaaaa(uv_hi[j+23], uv[j+23], uv_hi[j+23], uv[j+23], s1, s0);
            umul_ppmm(s1, s0, u[i+3], v[i+3][j+23]);
            add_ssaaaa(uv_hi[j+23], uv[j+23], uv_hi[j+23], uv[j+23], s1, s0);

            umul_ppmm(s1, s0, u[i], v[i][j+24]);
            add_ssaaaa(uv_hi[j+24], uv[j+24], uv_hi[j+24], uv[j+24], s1, s0);
            umul_ppmm(s1, s0, u[i+1], v[i+1][j+24]);
            add_ssaaaa(uv_hi[j+24], uv[j+24], uv_hi[j+24], uv[j+24], s1, s0);
            umul_ppmm(s1, s0, u[i+2], v[i+2][j+24]);
            add_ssaaaa(uv_hi[j+24], uv[j+24], uv_hi[j+24], uv[j+24], s1, s0);
            umul_ppmm(s1, s0, u[i+3], v[i+3][j+24]);
            add_ssaaaa(uv_hi[j+24], uv[j+24], uv_hi[j+24], uv[j+24], s1, s0);

            umul_ppmm(s1, s0, u[i], v[i][j+25]);
            add_ssaaaa(uv_hi[j+25], uv[j+25], uv_hi[j+25], uv[j+25], s1, s0);
            umul_ppmm(s1, s0, u[i+1], v[i+1][j+25]);
            add_ssaaaa(uv_hi[j+25], uv[j+25], uv_hi[j+25], uv[j+25], s1, s0);
            umul_ppmm(s1, s0, u[i+2], v[i+2][j+25]);
            add_ssaaaa(uv_hi[j+25], uv[j+25], uv_hi[j+25], uv[j+25], s1, s0);
            umul_ppmm(s1, s0, u[i+3], v[i+3][j+25]);
            add_ssaaaa(uv_hi[j+25], uv[j+25], uv_hi[j+25], uv[j+25], s1, s0);

            umul_ppmm(s1, s0, u[i], v[i][j+26]);
            add_ssaaaa(uv_hi[j+26], uv[j+26], uv_hi[j+26], uv[j+26], s1, s0);
            umul_ppmm(s1, s0, u[i+1], v[i+1][j+26]);
            add_ssaaaa(uv_hi[j+26], uv[j+26], uv_hi[j+26], uv[j+26], s1, s0);
            umul_ppmm(s1, s0, u[i+2], v[i+2][j+26]);
            add_ssaaaa(uv_hi[j+26], uv[j+26], uv_hi[j+26], uv[j+26], s1, s0);
            umul_ppmm(s1, s0, u[i+3], v[i+3][j+26]);
            add_ssaaaa(uv_hi[j+26], uv[j+26], uv_hi[j+26], uv[j+26], s1, s0);

            umul_ppmm(s1, s0, u[i], v[i][j+27]);
            add_ssaaaa(uv_hi[j+27], uv[j+27], uv_hi[j+27], uv[j+27], s1, s0);
            umul_ppmm(s1, s0, u[i+1], v[i+1][j+27]);
            add_ssaaaa(uv_hi[j+27], uv[j+27], uv_hi[j+27], uv[j+27], s1, s0);
            umul_ppmm(s1, s0, u[i+2], v[i+2][j+27]);
            add_ssaaaa(uv_hi[j+27], uv[j+27], uv_hi[j+27], uv[j+27], s1, s0);
            umul_ppmm(s1, s0, u[i+3], v[i+3][j+27]);
            add_ssaaaa(uv_hi[j+27], uv[j+27], uv_hi[j+27], uv[j+27], s1, s0);

            umul_ppmm(s1, s0, u[i], v[i][j+28]);
            add_ssaaaa(uv_hi[j+28], uv[j+28], uv_hi[j+28], uv[j+28], s1, s0);
            umul_ppmm(s1, s0, u[i+1], v[i+1][j+28]);
            add_ssaaaa(uv_hi[j+28], uv[j+28], uv_hi[j+28], uv[j+28], s1, s0);
            umul_ppmm(s1, s0, u[i+2], v[i+2][j+28]);
            add_ssaaaa(uv_hi[j+28], uv[j+28], uv_hi[j+28], uv[j+28], s1, s0);
            umul_ppmm(s1, s0, u[i+3], v[i+3][j+28]);
            add_ssaaaa(uv_hi[j+28], uv[j+28], uv_hi[j+28], uv[j+28], s1, s0);

            umul_ppmm(s1, s0, u[i], v[i][j+29]);
            add_ssaaaa(uv_hi[j+29], uv[j+29], uv_hi[j+29], uv[j+29], s1, s0);
            umul_ppmm(s1, s0, u[i+1], v[i+1][j+29]);
            add_ssaaaa(uv_hi[j+29], uv[j+29], uv_hi[j+29], uv[j+29], s1, s0);
            umul_ppmm(s1, s0, u[i+2], v[i+2][j+29]);
            add_ssaaaa(uv_hi[j+29], uv[j+29], uv_hi[j+29], uv[j+29], s1, s0);
            umul_ppmm(s1, s0, u[i+3], v[i+3][j+29]);
            add_ssaaaa(uv_hi[j+29], uv[j+29], uv_hi[j+29], uv[j+29], s1, s0);

            umul_ppmm(s1, s0, u[i], v[i][j+30]);
            add_ssaaaa(uv_hi[j+30], uv[j+30], uv_hi[j+30], uv[j+30], s1, s0);
            umul_ppmm(s1, s0, u[i+1], v[i+1][j+30]);
            add_ssaaaa(uv_hi[j+30], uv[j+30], uv_hi[j+30], uv[j+30], s1, s0);
            umul_ppmm(s1, s0, u[i+2], v[i+2][j+30]);
            add_ssaaaa(uv_hi[j+30], uv[j+30], uv_hi[j+30], uv[j+30], s1, s0);
            umul_ppmm(s1, s0, u[i+3], v[i+3][j+30]);
            add_ssaaaa(uv_hi[j+30], uv[j+30], uv_hi[j+30], uv[j+30], s1, s0);

            umul_ppmm(s1, s0, u[i], v[i][j+31]);
            add_ssaaaa(uv_hi[j+31], uv[j+31], uv_hi[j+31], uv[j+31], s1, s0);
            umul_ppmm(s1, s0, u[i+1], v[i+1][j+31]);
            add_ssaaaa(uv_hi[j+31], uv[j+31], uv_hi[j+31], uv[j+31], s1, s0);
            umul_ppmm(s1, s0, u[i+2], v[i+2][j+31]);
            add_ssaaaa(uv_hi[j+31], uv[j+31], uv_hi[j+31], uv[j+31], s1, s0);
            umul_ppmm(s1, s0, u[i+3], v[i+3][j+31]);
            add_ssaaaa(uv_hi[j+31], uv[j+31], uv_hi[j+31], uv[j+31], s1, s0);
        }
        for (; j < k; j++)
        {
            umul_ppmm(s1, s0, u[i], v[i][j+0]);
            add_ssaaaa(uv_hi[j+0], uv[j+0], uv_hi[j+0], uv[j+0], s1, s0);
            umul_ppmm(s1, s0, u[i+1], v[i+1][j+0]);
            add_ssaaaa(uv_hi[j+0], uv[j+0], uv_hi[j+0], uv[j+0], s1, s0);
            umul_ppmm(s1, s0, u[i+2], v[i+2][j+0]);
            add_ssaaaa(uv_hi[j+0], uv[j+0], uv_hi[j+0], uv[j+0], s1, s0);
            umul_ppmm(s1, s0, u[i+3], v[i+3][j+0]);
            add_ssaaaa(uv_hi[j+0], uv[j+0], uv_hi[j+0], uv[j+0], s1, s0);
        }
    }
    for (; i < len; i++)
    {
        ulong j = 0;
        for (; j+31 < k; j+=32)
        {
            umul_ppmm(s1, s0, u[i], v[i][j+0]);
            add_ssaaaa(uv_hi[j+0], uv[j+0], uv_hi[j+0], uv[j+0], s1, s0);
            umul_ppmm(s1, s0, u[i], v[i][j+1]);
            add_ssaaaa(uv_hi[j+1], uv[j+1], uv_hi[j+1], uv[j+1], s1, s0);
            umul_ppmm(s1, s0, u[i], v[i][j+2]);
            add_ssaaaa(uv_hi[j+2], uv[j+2], uv_hi[j+2], uv[j+2], s1, s0);
            umul_ppmm(s1, s0, u[i], v[i][j+3]);
            add_ssaaaa(uv_hi[j+3], uv[j+3], uv_hi[j+3], uv[j+3], s1, s0);
            umul_ppmm(s1, s0, u[i], v[i][j+4]);
            add_ssaaaa(uv_hi[j+4], uv[j+4], uv_hi[j+4], uv[j+4], s1, s0);
            umul_ppmm(s1, s0, u[i], v[i][j+5]);
            add_ssaaaa(uv_hi[j+5], uv[j+5], uv_hi[j+5], uv[j+5], s1, s0);
            umul_ppmm(s1, s0, u[i], v[i][j+6]);
            add_ssaaaa(uv_hi[j+6], uv[j+6], uv_hi[j+6], uv[j+6], s1, s0);
            umul_ppmm(s1, s0, u[i], v[i][j+7]);
            add_ssaaaa(uv_hi[j+7], uv[j+7], uv_hi[j+7], uv[j+7], s1, s0);
            umul_ppmm(s1, s0, u[i], v[i][j+8]);
            add_ssaaaa(uv_hi[j+8], uv[j+8], uv_hi[j+8], uv[j+8], s1, s0);
            umul_ppmm(s1, s0, u[i], v[i][j+9]);
            add_ssaaaa(uv_hi[j+9], uv[j+9], uv_hi[j+9], uv[j+9], s1, s0);
            umul_ppmm(s1, s0, u[i], v[i][j+10]);
            add_ssaaaa(uv_hi[j+10], uv[j+10], uv_hi[j+10], uv[j+10], s1, s0);
            umul_ppmm(s1, s0, u[i], v[i][j+11]);
            add_ssaaaa(uv_hi[j+11], uv[j+11], uv_hi[j+11], uv[j+11], s1, s0);
            umul_ppmm(s1, s0, u[i], v[i][j+12]);
            add_ssaaaa(uv_hi[j+12], uv[j+12], uv_hi[j+12], uv[j+12], s1, s0);
            umul_ppmm(s1, s0, u[i], v[i][j+13]);
            add_ssaaaa(uv_hi[j+13], uv[j+13], uv_hi[j+13], uv[j+13], s1, s0);
            umul_ppmm(s1, s0, u[i], v[i][j+14]);
            add_ssaaaa(uv_hi[j+14], uv[j+14], uv_hi[j+14], uv[j+14], s1, s0);
            umul_ppmm(s1, s0, u[i], v[i][j+15]);
            add_ssaaaa(uv_hi[j+15], uv[j+15], uv_hi[j+15], uv[j+15], s1, s0);
            umul_ppmm(s1, s0, u[i], v[i][j+16]);
            add_ssaaaa(uv_hi[j+16], uv[j+16], uv_hi[j+16], uv[j+16], s1, s0);
            umul_ppmm(s1, s0, u[i], v[i][j+17]);
            add_ssaaaa(uv_hi[j+17], uv[j+17], uv_hi[j+17], uv[j+17], s1, s0);
            umul_ppmm(s1, s0, u[i], v[i][j+18]);
            add_ssaaaa(uv_hi[j+18], uv[j+18], uv_hi[j+18], uv[j+18], s1, s0);
            umul_ppmm(s1, s0, u[i], v[i][j+19]);
            add_ssaaaa(uv_hi[j+19], uv[j+19], uv_hi[j+19], uv[j+19], s1, s0);
            umul_ppmm(s1, s0, u[i], v[i][j+20]);
            add_ssaaaa(uv_hi[j+20], uv[j+20], uv_hi[j+20], uv[j+20], s1, s0);
            umul_ppmm(s1, s0, u[i], v[i][j+21]);
            add_ssaaaa(uv_hi[j+21], uv[j+21], uv_hi[j+21], uv[j+21], s1, s0);
            umul_ppmm(s1, s0, u[i], v[i][j+22]);
            add_ssaaaa(uv_hi[j+22], uv[j+22], uv_hi[j+22], uv[j+22], s1, s0);
            umul_ppmm(s1, s0, u[i], v[i][j+23]);
            add_ssaaaa(uv_hi[j+23], uv[j+23], uv_hi[j+23], uv[j+23], s1, s0);
            umul_ppmm(s1, s0, u[i], v[i][j+24]);
            add_ssaaaa(uv_hi[j+24], uv[j+24], uv_hi[j+24], uv[j+24], s1, s0);
            umul_ppmm(s1, s0, u[i], v[i][j+25]);
            add_ssaaaa(uv_hi[j+25], uv[j+25], uv_hi[j+25], uv[j+25], s1, s0);
            umul_ppmm(s1, s0, u[i], v[i][j+26]);
            add_ssaaaa(uv_hi[j+26], uv[j+26], uv_hi[j+26], uv[j+26], s1, s0);
            umul_ppmm(s1, s0, u[i], v[i][j+27]);
            add_ssaaaa(uv_hi[j+27], uv[j+27], uv_hi[j+27], uv[j+27], s1, s0);
            umul_ppmm(s1, s0, u[i], v[i][j+28]);
            add_ssaaaa(uv_hi[j+28], uv[j+28], uv_hi[j+28], uv[j+28], s1, s0);
            umul_ppmm(s1, s0, u[i], v[i][j+29]);
            add_ssaaaa(uv_hi[j+29], uv[j+29], uv_hi[j+29], uv[j+29], s1, s0);
            umul_ppmm(s1, s0, u[i], v[i][j+30]);
            add_ssaaaa(uv_hi[j+30], uv[j+30], uv_hi[j+30], uv[j+30], s1, s0);
            umul_ppmm(s1, s0, u[i], v[i][j+31]);
            add_ssaaaa(uv_hi[j+31], uv[j+31], uv_hi[j+31], uv[j+31], s1, s0);
        }
        for (; j < k; j++)
        {
            umul_ppmm(s1, s0, u[i], v[i][j+0]);
            add_ssaaaa(uv_hi[j+0], uv[j+0], uv_hi[j+0], uv[j+0], s1, s0);
        }
    }

    for (ulong j = 0; j < k; j++)
        NMOD2_RED2(uv[j], uv_hi[j], uv[j], mod);
    _nmod_vec_clear(uv_hi);
}

void _nmod_vec_dot_product_multi_2_v8_8(mp_ptr uv, mp_srcptr u, mp_srcptr * v,
                                   ulong len, ulong k, nmod_t mod)
{
    mp_ptr uv_hi = _nmod_vec_init(k);
    _nmod_vec_zero(uv, k);
    _nmod_vec_zero(uv_hi, k);
    mp_limb_t s0, s1;
    ulong i = 0;
    for (; i+8 < len; i += 8)
    {
        ulong j = 0;
        for (; j+7 < k; j+=8)
        {
            umul_ppmm(s1, s0, u[i], v[i][j+0]);
            add_ssaaaa(uv_hi[j+0], uv[j+0], uv_hi[j+0], uv[j+0], s1, s0);
            umul_ppmm(s1, s0, u[i+1], v[i+1][j+0]);
            add_ssaaaa(uv_hi[j+0], uv[j+0], uv_hi[j+0], uv[j+0], s1, s0);
            umul_ppmm(s1, s0, u[i+2], v[i+2][j+0]);
            add_ssaaaa(uv_hi[j+0], uv[j+0], uv_hi[j+0], uv[j+0], s1, s0);
            umul_ppmm(s1, s0, u[i+3], v[i+3][j+0]);
            add_ssaaaa(uv_hi[j+0], uv[j+0], uv_hi[j+0], uv[j+0], s1, s0);
            umul_ppmm(s1, s0, u[i+4], v[i+4][j+0]);
            add_ssaaaa(uv_hi[j+0], uv[j+0], uv_hi[j+0], uv[j+0], s1, s0);
            umul_ppmm(s1, s0, u[i+5], v[i+5][j+0]);
            add_ssaaaa(uv_hi[j+0], uv[j+0], uv_hi[j+0], uv[j+0], s1, s0);
            umul_ppmm(s1, s0, u[i+6], v[i+6][j+0]);
            add_ssaaaa(uv_hi[j+0], uv[j+0], uv_hi[j+0], uv[j+0], s1, s0);
            umul_ppmm(s1, s0, u[i+7], v[i+7][j+0]);
            add_ssaaaa(uv_hi[j+0], uv[j+0], uv_hi[j+0], uv[j+0], s1, s0);

            umul_ppmm(s1, s0, u[i], v[i][j+1]);
            add_ssaaaa(uv_hi[j+1], uv[j+1], uv_hi[j+1], uv[j+1], s1, s0);
            umul_ppmm(s1, s0, u[i+1], v[i+1][j+1]);
            add_ssaaaa(uv_hi[j+1], uv[j+1], uv_hi[j+1], uv[j+1], s1, s0);
            umul_ppmm(s1, s0, u[i+2], v[i+2][j+1]);
            add_ssaaaa(uv_hi[j+1], uv[j+1], uv_hi[j+1], uv[j+1], s1, s0);
            umul_ppmm(s1, s0, u[i+3], v[i+3][j+1]);
            add_ssaaaa(uv_hi[j+1], uv[j+1], uv_hi[j+1], uv[j+1], s1, s0);
            umul_ppmm(s1, s0, u[i+4], v[i+4][j+1]);
            add_ssaaaa(uv_hi[j+1], uv[j+1], uv_hi[j+1], uv[j+1], s1, s0);
            umul_ppmm(s1, s0, u[i+5], v[i+5][j+1]);
            add_ssaaaa(uv_hi[j+1], uv[j+1], uv_hi[j+1], uv[j+1], s1, s0);
            umul_ppmm(s1, s0, u[i+6], v[i+6][j+1]);
            add_ssaaaa(uv_hi[j+1], uv[j+1], uv_hi[j+1], uv[j+1], s1, s0);
            umul_ppmm(s1, s0, u[i+7], v[i+7][j+1]);
            add_ssaaaa(uv_hi[j+1], uv[j+1], uv_hi[j+1], uv[j+1], s1, s0);

            umul_ppmm(s1, s0, u[i], v[i][j+2]);
            add_ssaaaa(uv_hi[j+2], uv[j+2], uv_hi[j+2], uv[j+2], s1, s0);
            umul_ppmm(s1, s0, u[i+1], v[i+1][j+2]);
            add_ssaaaa(uv_hi[j+2], uv[j+2], uv_hi[j+2], uv[j+2], s1, s0);
            umul_ppmm(s1, s0, u[i+2], v[i+2][j+2]);
            add_ssaaaa(uv_hi[j+2], uv[j+2], uv_hi[j+2], uv[j+2], s1, s0);
            umul_ppmm(s1, s0, u[i+3], v[i+3][j+2]);
            add_ssaaaa(uv_hi[j+2], uv[j+2], uv_hi[j+2], uv[j+2], s1, s0);
            umul_ppmm(s1, s0, u[i+4], v[i+4][j+2]);
            add_ssaaaa(uv_hi[j+2], uv[j+2], uv_hi[j+2], uv[j+2], s1, s0);
            umul_ppmm(s1, s0, u[i+5], v[i+5][j+2]);
            add_ssaaaa(uv_hi[j+2], uv[j+2], uv_hi[j+2], uv[j+2], s1, s0);
            umul_ppmm(s1, s0, u[i+6], v[i+6][j+2]);
            add_ssaaaa(uv_hi[j+2], uv[j+2], uv_hi[j+2], uv[j+2], s1, s0);
            umul_ppmm(s1, s0, u[i+7], v[i+7][j+2]);
            add_ssaaaa(uv_hi[j+2], uv[j+2], uv_hi[j+2], uv[j+2], s1, s0);

            umul_ppmm(s1, s0, u[i], v[i][j+3]);
            add_ssaaaa(uv_hi[j+3], uv[j+3], uv_hi[j+3], uv[j+3], s1, s0);
            umul_ppmm(s1, s0, u[i+1], v[i+1][j+3]);
            add_ssaaaa(uv_hi[j+3], uv[j+3], uv_hi[j+3], uv[j+3], s1, s0);
            umul_ppmm(s1, s0, u[i+2], v[i+2][j+3]);
            add_ssaaaa(uv_hi[j+3], uv[j+3], uv_hi[j+3], uv[j+3], s1, s0);
            umul_ppmm(s1, s0, u[i+3], v[i+3][j+3]);
            add_ssaaaa(uv_hi[j+3], uv[j+3], uv_hi[j+3], uv[j+3], s1, s0);
            umul_ppmm(s1, s0, u[i+4], v[i+4][j+3]);
            add_ssaaaa(uv_hi[j+3], uv[j+3], uv_hi[j+3], uv[j+3], s1, s0);
            umul_ppmm(s1, s0, u[i+5], v[i+5][j+3]);
            add_ssaaaa(uv_hi[j+3], uv[j+3], uv_hi[j+3], uv[j+3], s1, s0);
            umul_ppmm(s1, s0, u[i+6], v[i+6][j+3]);
            add_ssaaaa(uv_hi[j+3], uv[j+3], uv_hi[j+3], uv[j+3], s1, s0);
            umul_ppmm(s1, s0, u[i+7], v[i+7][j+3]);
            add_ssaaaa(uv_hi[j+3], uv[j+3], uv_hi[j+3], uv[j+3], s1, s0);

            umul_ppmm(s1, s0, u[i], v[i][j+4]);
            add_ssaaaa(uv_hi[j+4], uv[j+4], uv_hi[j+4], uv[j+4], s1, s0);
            umul_ppmm(s1, s0, u[i+1], v[i+1][j+4]);
            add_ssaaaa(uv_hi[j+4], uv[j+4], uv_hi[j+4], uv[j+4], s1, s0);
            umul_ppmm(s1, s0, u[i+2], v[i+2][j+4]);
            add_ssaaaa(uv_hi[j+4], uv[j+4], uv_hi[j+4], uv[j+4], s1, s0);
            umul_ppmm(s1, s0, u[i+3], v[i+3][j+4]);
            add_ssaaaa(uv_hi[j+4], uv[j+4], uv_hi[j+4], uv[j+4], s1, s0);
            umul_ppmm(s1, s0, u[i+4], v[i+4][j+4]);
            add_ssaaaa(uv_hi[j+4], uv[j+4], uv_hi[j+4], uv[j+4], s1, s0);
            umul_ppmm(s1, s0, u[i+5], v[i+5][j+4]);
            add_ssaaaa(uv_hi[j+4], uv[j+4], uv_hi[j+4], uv[j+4], s1, s0);
            umul_ppmm(s1, s0, u[i+6], v[i+6][j+4]);
            add_ssaaaa(uv_hi[j+4], uv[j+4], uv_hi[j+4], uv[j+4], s1, s0);
            umul_ppmm(s1, s0, u[i+7], v[i+7][j+4]);
            add_ssaaaa(uv_hi[j+4], uv[j+4], uv_hi[j+4], uv[j+4], s1, s0);

            umul_ppmm(s1, s0, u[i], v[i][j+5]);
            add_ssaaaa(uv_hi[j+5], uv[j+5], uv_hi[j+5], uv[j+5], s1, s0);
            umul_ppmm(s1, s0, u[i+1], v[i+1][j+5]);
            add_ssaaaa(uv_hi[j+5], uv[j+5], uv_hi[j+5], uv[j+5], s1, s0);
            umul_ppmm(s1, s0, u[i+2], v[i+2][j+5]);
            add_ssaaaa(uv_hi[j+5], uv[j+5], uv_hi[j+5], uv[j+5], s1, s0);
            umul_ppmm(s1, s0, u[i+3], v[i+3][j+5]);
            add_ssaaaa(uv_hi[j+5], uv[j+5], uv_hi[j+5], uv[j+5], s1, s0);
            umul_ppmm(s1, s0, u[i+4], v[i+4][j+5]);
            add_ssaaaa(uv_hi[j+5], uv[j+5], uv_hi[j+5], uv[j+5], s1, s0);
            umul_ppmm(s1, s0, u[i+5], v[i+5][j+5]);
            add_ssaaaa(uv_hi[j+5], uv[j+5], uv_hi[j+5], uv[j+5], s1, s0);
            umul_ppmm(s1, s0, u[i+6], v[i+6][j+5]);
            add_ssaaaa(uv_hi[j+5], uv[j+5], uv_hi[j+5], uv[j+5], s1, s0);
            umul_ppmm(s1, s0, u[i+7], v[i+7][j+5]);
            add_ssaaaa(uv_hi[j+5], uv[j+5], uv_hi[j+5], uv[j+5], s1, s0);

            umul_ppmm(s1, s0, u[i], v[i][j+6]);
            add_ssaaaa(uv_hi[j+6], uv[j+6], uv_hi[j+6], uv[j+6], s1, s0);
            umul_ppmm(s1, s0, u[i+1], v[i+1][j+6]);
            add_ssaaaa(uv_hi[j+6], uv[j+6], uv_hi[j+6], uv[j+6], s1, s0);
            umul_ppmm(s1, s0, u[i+2], v[i+2][j+6]);
            add_ssaaaa(uv_hi[j+6], uv[j+6], uv_hi[j+6], uv[j+6], s1, s0);
            umul_ppmm(s1, s0, u[i+3], v[i+3][j+6]);
            add_ssaaaa(uv_hi[j+6], uv[j+6], uv_hi[j+6], uv[j+6], s1, s0);
            umul_ppmm(s1, s0, u[i+4], v[i+4][j+6]);
            add_ssaaaa(uv_hi[j+6], uv[j+6], uv_hi[j+6], uv[j+6], s1, s0);
            umul_ppmm(s1, s0, u[i+5], v[i+5][j+6]);
            add_ssaaaa(uv_hi[j+6], uv[j+6], uv_hi[j+6], uv[j+6], s1, s0);
            umul_ppmm(s1, s0, u[i+6], v[i+6][j+6]);
            add_ssaaaa(uv_hi[j+6], uv[j+6], uv_hi[j+6], uv[j+6], s1, s0);
            umul_ppmm(s1, s0, u[i+7], v[i+7][j+6]);
            add_ssaaaa(uv_hi[j+6], uv[j+6], uv_hi[j+6], uv[j+6], s1, s0);

            umul_ppmm(s1, s0, u[i], v[i][j+7]);
            add_ssaaaa(uv_hi[j+7], uv[j+7], uv_hi[j+7], uv[j+7], s1, s0);
            umul_ppmm(s1, s0, u[i+1], v[i+1][j+7]);
            add_ssaaaa(uv_hi[j+7], uv[j+7], uv_hi[j+7], uv[j+7], s1, s0);
            umul_ppmm(s1, s0, u[i+2], v[i+2][j+7]);
            add_ssaaaa(uv_hi[j+7], uv[j+7], uv_hi[j+7], uv[j+7], s1, s0);
            umul_ppmm(s1, s0, u[i+3], v[i+3][j+7]);
            add_ssaaaa(uv_hi[j+7], uv[j+7], uv_hi[j+7], uv[j+7], s1, s0);
            umul_ppmm(s1, s0, u[i+4], v[i+4][j+7]);
            add_ssaaaa(uv_hi[j+7], uv[j+7], uv_hi[j+7], uv[j+7], s1, s0);
            umul_ppmm(s1, s0, u[i+5], v[i+5][j+7]);
            add_ssaaaa(uv_hi[j+7], uv[j+7], uv_hi[j+7], uv[j+7], s1, s0);
            umul_ppmm(s1, s0, u[i+6], v[i+6][j+7]);
            add_ssaaaa(uv_hi[j+7], uv[j+7], uv_hi[j+7], uv[j+7], s1, s0);
            umul_ppmm(s1, s0, u[i+7], v[i+7][j+7]);
            add_ssaaaa(uv_hi[j+7], uv[j+7], uv_hi[j+7], uv[j+7], s1, s0);
        }
        for (; j < k; j++)
        {
            umul_ppmm(s1, s0, u[i], v[i][j+0]);
            add_ssaaaa(uv_hi[j+0], uv[j+0], uv_hi[j+0], uv[j+0], s1, s0);
            umul_ppmm(s1, s0, u[i+1], v[i+1][j+0]);
            add_ssaaaa(uv_hi[j+0], uv[j+0], uv_hi[j+0], uv[j+0], s1, s0);
            umul_ppmm(s1, s0, u[i+2], v[i+2][j+0]);
            add_ssaaaa(uv_hi[j+0], uv[j+0], uv_hi[j+0], uv[j+0], s1, s0);
            umul_ppmm(s1, s0, u[i+3], v[i+3][j+0]);
            add_ssaaaa(uv_hi[j+0], uv[j+0], uv_hi[j+0], uv[j+0], s1, s0);
            umul_ppmm(s1, s0, u[i+4], v[i+4][j+0]);
            add_ssaaaa(uv_hi[j+0], uv[j+0], uv_hi[j+0], uv[j+0], s1, s0);
            umul_ppmm(s1, s0, u[i+5], v[i+5][j+0]);
            add_ssaaaa(uv_hi[j+0], uv[j+0], uv_hi[j+0], uv[j+0], s1, s0);
            umul_ppmm(s1, s0, u[i+6], v[i+6][j+0]);
            add_ssaaaa(uv_hi[j+0], uv[j+0], uv_hi[j+0], uv[j+0], s1, s0);
            umul_ppmm(s1, s0, u[i+7], v[i+7][j+0]);
            add_ssaaaa(uv_hi[j+0], uv[j+0], uv_hi[j+0], uv[j+0], s1, s0);
        }
    }
    for (; i < len; i++)
    {
        ulong j = 0;
        for (; j+7 < k; j+=8)
        {
            umul_ppmm(s1, s0, u[i], v[i][j+0]);
            add_ssaaaa(uv_hi[j+0], uv[j+0], uv_hi[j+0], uv[j+0], s1, s0);
            umul_ppmm(s1, s0, u[i], v[i][j+1]);
            add_ssaaaa(uv_hi[j+1], uv[j+1], uv_hi[j+1], uv[j+1], s1, s0);
            umul_ppmm(s1, s0, u[i], v[i][j+2]);
            add_ssaaaa(uv_hi[j+2], uv[j+2], uv_hi[j+2], uv[j+2], s1, s0);
            umul_ppmm(s1, s0, u[i], v[i][j+3]);
            add_ssaaaa(uv_hi[j+3], uv[j+3], uv_hi[j+3], uv[j+3], s1, s0);
            umul_ppmm(s1, s0, u[i], v[i][j+4]);
            add_ssaaaa(uv_hi[j+4], uv[j+4], uv_hi[j+4], uv[j+4], s1, s0);
            umul_ppmm(s1, s0, u[i], v[i][j+5]);
            add_ssaaaa(uv_hi[j+5], uv[j+5], uv_hi[j+5], uv[j+5], s1, s0);
            umul_ppmm(s1, s0, u[i], v[i][j+6]);
            add_ssaaaa(uv_hi[j+6], uv[j+6], uv_hi[j+6], uv[j+6], s1, s0);
            umul_ppmm(s1, s0, u[i], v[i][j+7]);
            add_ssaaaa(uv_hi[j+7], uv[j+7], uv_hi[j+7], uv[j+7], s1, s0);
        }
        for (; j < k; j++)
        {
            umul_ppmm(s1, s0, u[i], v[i][j+0]);
            add_ssaaaa(uv_hi[j+0], uv[j+0], uv_hi[j+0], uv[j+0], s1, s0);
        }
    }

    for (ulong j = 0; j < k; j++)
        NMOD2_RED2(uv[j], uv_hi[j], uv[j], mod);
    _nmod_vec_clear(uv_hi);
}


static inline
void _nmod_vec_dot_product_multi_3(mp_ptr uv, mp_srcptr u, mp_srcptr * v,
                                   ulong len, ulong k,
                                   ulong max_bits_u, ulong max_bits_v,
                                   nmod_t mod)
{
    /* number of products we can do before overflow */
    const ulong log_nt = 2*FLINT_BITS - (max_bits_u + max_bits_v);
    const ulong num_terms = (log_nt < FLINT_BITS) ? (UWORD(1) << log_nt) : (UWORD_MAX);

    mp_ptr uv_hi = _nmod_vec_init(k);   // 64->127
    mp_ptr uv_hii = _nmod_vec_init(k);  // 128->..
    _nmod_vec_zero(uv, k);
    _nmod_vec_zero(uv_hi, k);
    _nmod_vec_zero(uv_hii, k);
    mp_limb_t s0, s1, t0, t1;
    ulong i = 0;
    if (num_terms >= 8)
        for (; i+7 < len; i += 8) // FIXME 8 vs 4 ? (Vincent: not tested)
        {
            for (ulong j = 0; j < k; j++)
            {
                umul_ppmm(t1, t0, u[i+0], v[i+0][j]);
                umul_ppmm(s1, s0, u[i+1], v[i+1][j]);
                add_ssaaaa(t1, t0, t1, t0, s1, s0);
                umul_ppmm(s1, s0, u[i+2], v[i+2][j]);
                add_ssaaaa(t1, t0, t1, t0, s1, s0);
                umul_ppmm(s1, s0, u[i+3], v[i+3][j]);
                add_ssaaaa(t1, t0, t1, t0, s1, s0);
                umul_ppmm(s1, s0, u[i+4], v[i+4][j]);
                add_ssaaaa(t1, t0, t1, t0, s1, s0);
                umul_ppmm(s1, s0, u[i+5], v[i+5][j]);
                add_ssaaaa(t1, t0, t1, t0, s1, s0);
                umul_ppmm(s1, s0, u[i+6], v[i+6][j]);
                add_ssaaaa(t1, t0, t1, t0, s1, s0);
                umul_ppmm(s1, s0, u[i+7], v[i+7][j]);
                add_ssaaaa(t1, t0, t1, t0, s1, s0);
                add_sssaaaaaa(uv_hii[j], uv_hi[j], uv[j], uv_hii[j], uv_hi[j], uv[j], UWORD(0), t1, t0);
            }
        }
    else
        for (; i+num_terms < len; i += num_terms)
        {
            for (ulong j = 0; j < k; j++)
            {
                umul_ppmm(t1, t0, u[i+0], v[i+0][j]);
                for (ulong ii = 1; ii < num_terms; ii++)
                {
                    umul_ppmm(s1, s0, u[i+ii], v[i+ii][j]);
                    add_ssaaaa(t1, t0, t1, t0, s1, s0);
                }
                add_sssaaaaaa(uv_hii[j], uv_hi[j], uv[j], uv_hii[j], uv_hi[j], uv[j], UWORD(0), t1, t0);
            }
        }

    // TODO see if good for small len; unclear (order of loops, etc)
    for (ulong j = 0; j < k; j++)
    {
        t0 = UWORD(0);
        t1 = UWORD(0);
        for (ulong ii = i; ii < len; ii++)
        {
            umul_ppmm(s1, s0, u[ii], v[ii][j]);
            add_ssaaaa(t1, t0, t1, t0, s1, s0);
        }
        add_sssaaaaaa(uv_hii[j], uv_hi[j], uv[j], uv_hii[j], uv_hi[j], uv[j], UWORD(0), t1, t0);
        NMOD_RED(uv_hii[j], uv_hii[j], mod);
        NMOD_RED3(uv[j], uv_hii[j], uv_hi[j], uv[j], mod);
    }
    _nmod_vec_clear(uv_hii);
    _nmod_vec_clear(uv_hi);
}

void nmod_vec_dot_product_multi(mp_ptr uv, mp_srcptr u, mp_srcptr * v,
                                ulong len, ulong k,
                                ulong max_bits_u, ulong max_bits_v,
                                nmod_t mod)
{
    const ulong n_limbs = _nmod_vec_dot_bound_limbs_unbalanced(len, max_bits_u, max_bits_v);

    if (n_limbs == 2)
    {
        _nmod_vec_dot_product_multi_2(uv, u, v, len, k, mod);
        return;
    }
    if (n_limbs == 3)
    {
        _nmod_vec_dot_product_multi_3(uv, u, v, len, k, max_bits_u, max_bits_v, mod);
        return;
    }
    if (n_limbs == 1)
    {
        _nmod_vec_dot_product_multi_1(uv, u, v, len, k, mod);
        return;
    }
    // n_limbs == 0
    _nmod_vec_zero(uv, k);
    return;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
