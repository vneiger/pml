#include <NTL/matrix.h>
#include <NTL/mat_lzz_p.h>
#include <NTL/lzz_pX.h>

#include "lzz_p_extra.h"
#include "mat_lzz_pX_extra.h"
#include "lzz_pX_CRT.h"

NTL_CLIENT


/*------------------------------------------------------------*/
/* c = a*b                                                    */
/* c can alias a or b; c does not have to be zero             */
/*------------------------------------------------------------*/
void multiply_transform_naive(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b)
{
    long dA = deg(a);
    long dB = deg(b);

    zz_pX_Transform_naive trs_naive(max(dA, dB) + 1);

    Vec<Mat<zz_p>> valA, valB, valC;
    trs_naive.forward_left_matrix(valA, a);
    trs_naive.forward_right_matrix(valB, b);

    long len = trs_naive.transform_length();
    valC.SetLength(len);
    for (long i = 0; i < len; i++)
    {
        mul(valC[i], valA[i], valB[i]);
    }

    trs_naive.backward_matrix(c, valC);
}

/*------------------------------------------------------------*/
/* c = a*b                                                    */
/* c can alias a or b; c does not have to be zero             */
/*------------------------------------------------------------*/
void multiply_transform_karatsuba(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b)
{
    zz_pX_Transform_karatsuba trs_karatsuba = zz_pX_Transform_karatsuba();

    Vec<Mat<zz_p>> valA, valB, valC;
    trs_karatsuba.forward_left_matrix(valA, a);
    trs_karatsuba.forward_right_matrix(valB, b);

    long len = trs_karatsuba.transform_length();
    valC.SetLength(len);
    for (long i = 0; i < len; i++)
    {
        mul(valC[i], valA[i], valB[i]);
    }

    trs_karatsuba.backward_matrix(c, valC);
}

/*------------------------------------------------------------*/
/* c = a*b                                                    */
/* c can alias a or b; c does not have to be zero             */
/*------------------------------------------------------------*/
void multiply_transform_montgomery3(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b)
{
    zz_pX_Transform_montgomery3 trs_montgomery3 = zz_pX_Transform_montgomery3();

    Vec<Mat<zz_p>> valA, valB, valC;
    trs_montgomery3.forward_left_matrix(valA, a);
    trs_montgomery3.forward_right_matrix(valB, b);

    long len = trs_montgomery3.transform_length();
    valC.SetLength(len);
    for (long i = 0; i < len; i++)
    {
        mul(valC[i], valA[i], valB[i]);
    }

    trs_montgomery3.backward_matrix(c, valC);
}

/*------------------------------------------------------------*/
/* c = a*b                                                    */
/* c can alias a or b; c does not have to be zero             */
/*------------------------------------------------------------*/
void multiply_transform_karatsuba4(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b)
{
    zz_pX_Transform_karatsuba4 trs_karatsuba4 = zz_pX_Transform_karatsuba4();

    Vec<Mat<zz_p>> valA, valB, valC;
    trs_karatsuba4.forward_left_matrix(valA, a);
    trs_karatsuba4.forward_right_matrix(valB, b);

    long len = trs_karatsuba4.transform_length();
    valC.SetLength(len);
    for (long i = 0; i < len; i++)
    {
        mul(valC[i], valA[i], valB[i]);
    }

    trs_karatsuba4.backward_matrix(c, valC);
}

/*------------------------------------------------------------*/
/* c = a*b                                                    */
/* c can alias a or b; c does not have to be zero             */
/*------------------------------------------------------------*/
void multiply_transform(Mat<zz_pX> & c, const Mat<zz_pX> & a, const Mat<zz_pX> & b, long len)
{
    switch (len) 
    {
    case 2:
        multiply_transform_karatsuba(c, a, b);
        break;
    case 3:
        multiply_transform_montgomery3(c, a, b);
        break;
    case 4:
        multiply_transform_karatsuba4(c, a, b);
        break;
    default:
        multiply_transform_naive(c, a, b);
    }
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
