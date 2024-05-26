#include <flint/nmod.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>

#include "nmod_mat_extra.h"
#include "nmod_poly_mat_multiply.h"

/*------------------------------------------------------------*/
/* computes the matrices for evaluation and interpolation     */
/*------------------------------------------------------------*/
static
void vandermonde_init1(nmod_mat_t vdm1, nmod_mat_t vdm2, nmod_mat_t inv_vdm, long d1, long d2, nmod_t mod)
{
    ulong i, j, s1, s2, nb_points;
    ulong p;
    nmod_mat_t vdm;
    
    // sizes = degree + 1
    s1 = d1 + 1;
    s2 = d2 + 1;
    // nb points, which is sum of degrees + 1
    nb_points = d1 + d2 + 1;

    p = mod.n;
    
    // will be Vandermonde matrix with nb_points, chosen as 0, 1, .., nb_points-1
    // row i contains 1, i, i^2, .., i^{sk-1} for k=1,2
    nmod_mat_init(vdm1, nb_points, s1, p);
    nmod_mat_init(vdm2, nb_points, s2, p);

    // vdm: square Vandermonde matrix with nb_points
    // points chosen as 0, 1, .., nb_points-1
    // --> row i contains 1, i, i^2, .., i^{nb_points-1}
    // vdm1: nb_points x s1 submatrix of vdm
    // vdm2: nb_points x s2 submatrix of vdm
    nmod_mat_init(vdm, nb_points, nb_points, p);
    nmod_mat_init(inv_vdm, nb_points, nb_points, p);
    
    nmod_mat_entry(vdm, 0, 0) = 1;
    nmod_mat_entry(vdm1, 0, 0) = 1;
    nmod_mat_entry(vdm2, 0 ,0) = 1;
    for (i = 1; i < nb_points; i++)
    {
        ulong p1;
        NMOD_RED(p1, i, mod);
        nmod_mat_entry(vdm, i, 0) = 1;
        nmod_mat_entry(vdm1, i, 0) = 1;
        nmod_mat_entry(vdm2, i, 0) = 1;
        for (j = 1; j < nb_points; j++)
        {
            nmod_mat_entry(vdm, i, j) = nmod_mul(nmod_mat_entry(vdm, i, j-1), p1, mod);
            if (j < s1)
                nmod_mat_entry(vdm1, i, j) = nmod_mat_entry(vdm, i, j);
            if (j < s2)
                nmod_mat_entry(vdm2, i, j) = nmod_mat_entry(vdm, i, j);
        }
    }

    // inv_vdm is the inverse of vdm
    nmod_mat_inv(inv_vdm, vdm);
    nmod_mat_clear(vdm);
}

/*------------------------------------------------------------*/
/* c = a*b                                                    */
/* output may alias input; c does not have to be zero matrix  */
/* matrix multiplication using the algorithm of Doliskani et al. */
/* uses matrix multiplication for evaluation and interpolation*/
/*------------------------------------------------------------*/
void nmod_poly_mat_mul_vdm1(nmod_poly_mat_t c, const nmod_poly_mat_t a, const nmod_poly_mat_t b)
{
    ulong i, j, k, ellA, ellB, dA, dB, m, n, p, ell, u, v, nb_points;
    nmod_mat_t vA, vB, iv, tmp_mat, valA, valB, valAp, valBp, valC, valCp;
    nmod_t mod;
    
    ellA = nmod_poly_mat_max_length(a);
    ellB = nmod_poly_mat_max_length(b);

    dA = ellA - 1;
    dB = ellB - 1;
    
    m = a->r;
    n = a->c;
    p = b->c;

    if (ellA == 0 || ellB == 0)
    {
        nmod_poly_mat_zero(c);
        return;
    }

    if (c == a || c == b)
    {
        nmod_poly_mat_t T;
        nmod_poly_mat_init(T, m, n, p);
        nmod_poly_mat_mul_vdm1(T, a, b);
        nmod_poly_mat_swap_entrywise(c, T);
        nmod_poly_mat_clear(T);
        return;
    }

    
    nmod_init(&mod, a->modulus);
    vandermonde_init1(vA, vB, iv, dA, dB, mod);
    nb_points = vA->r;

    // evaluation of matrix a:
    // build tmp_mat, whose column ell = i*n + j is the coefficient vector of
    // a[i][j] (padded with zeroes up to length dA+1 if necessary)
    nmod_mat_init(tmp_mat, dA + 1, m * n, mod.n);
    ell = 0;
    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++, ell++)
        {
            ulong d;
            nn_ptr ptr;
            d = nmod_poly_degree(a->rows[i] + j);
            ptr = (a->rows[i] + j)->coeffs;
            for (k = 0; k <= d; k++)
                tmp_mat->rows[k][ell] = ptr[k];
        }

    // note: d = deg(a[i][j]) is -1 if a[i][j] == 0
    // all non-touched entries already zero since tmp_mat was initialized as zero

    // valA: column ell = i*n + j contains the evaluations of a[i][j]
    nmod_mat_init(valA, vA->r, tmp_mat->c, mod.n);
    nmod_mat_mul_pml(valA, vA, tmp_mat);

    // evaluation of matrix b:
    // build tmp_mat, whose column ell = i*n + j is the coefficient vector of
    // b[i][j] (padded with zeroes up to length dB+1 if necessary)
    nmod_mat_clear(tmp_mat);
    nmod_mat_init(tmp_mat, dB + 1, n * p, mod.n);
    ell = 0;
    for (i = 0; i < n; i++)
        for ( j = 0; j < p; j++, ell++)
        {
            ulong d;
            nn_ptr ptr;
            d = nmod_poly_degree(b->rows[i] + j);
            ptr = (b->rows[i] + j)->coeffs;
            for (k = 0; k <= d; k++)
                tmp_mat->rows[k][ell] = ptr[k];
        }

    // valB: column ell = i*n + j contains the evaluations of b[i][j]
    nmod_mat_init(valB, vB->r, tmp_mat->c, mod.n);
    nmod_mat_mul_pml(valB, vB, tmp_mat);

    // perform the pointwise products
    nmod_mat_init(valAp, m, n, mod.n);
    nmod_mat_init(valBp, n, p, mod.n);
    nmod_mat_init(valC, nb_points, m * p, mod.n);
    nmod_mat_init(valCp, m, p, mod.n);
    
    for (i = 0; i < nb_points; i++)
    {
        // a evaluated at point i
        ell = 0;
        for (u = 0; u < m; u++)
            for (v = 0; v < n; v++, ell++)
                valAp->rows[u][v] = valA->rows[i][ell];

        // b evaluated at point i
        ell = 0;
        for (u = 0; u < n; u++)
            for (v = 0; v < p; v++, ell++)
                valBp->rows[u][v] = valB->rows[i][ell];

        // a*b evaluated at point i
        nmod_mat_mul_pml(valCp, valAp, valBp);

        // copy this into valC: column ell = i*n + j contains the evaluations
        // of the entry i,j of c = a*b
        ell = 0;
        for (u = 0; u < m; u++)
            for (v = 0; v < p; v++, ell++)
                valC->rows[i][ell] = valCp->rows[u][v];
    }

    // interpolate to find the entries of c
    nmod_mat_clear(tmp_mat);
    nmod_mat_init(tmp_mat, iv->r, valC->c, mod.n);
    nmod_mat_mul_pml(tmp_mat, iv, valC);

    // copy to output (reorganize these entries into c)
    ell = 0;
    for (u = 0; u < m; u++)
        for (v = 0; v < p; v++, ell++)
        {
            nn_ptr coeffs;
            nmod_poly_realloc(c->rows[u] + v, nb_points);
            coeffs = (c->rows[u] + v)->coeffs;
            (c->rows[u] + v)->length = nb_points; 
            for (i = 0; i < nb_points; i++)
                coeffs[i] = tmp_mat->rows[i][ell];
            _nmod_poly_normalise(c->rows[u] + v);
        }
    
    nmod_mat_clear(vA);
    nmod_mat_clear(vB);
    nmod_mat_clear(iv);
    nmod_mat_clear(tmp_mat);
    nmod_mat_clear(valA);
    nmod_mat_clear(valB);
    nmod_mat_clear(valAp);
    nmod_mat_clear(valBp);
    nmod_mat_clear(valC);
    nmod_mat_clear(valCp);
}



/** Middle product for polynomial matrices
 *  sets C = ((A * B) div x^dA) mod x^(dB+1)
 *  output can alias input
 *  uses evaluation and interpolation at arithmetic points, done by matrix products
 *  ASSUME: large enough field 
 */
void nmod_poly_mat_middle_product_vdm1(nmod_poly_mat_t c, const nmod_poly_mat_t a, const nmod_poly_mat_t b,
                                       const ulong dA, const ulong dB)
{

    ulong i, j, k, ellA, ellB, m, n, p, ell, u, v, nb_points;
    nmod_mat_t vA, vB, vB_t, iv, iv_t, tmp_mat, valA, valB, valAp, valBp, valC, valCp;
    nmod_t mod;
    
    ellA = nmod_poly_mat_max_length(a);
    ellB = nmod_poly_mat_max_length(b);

    if (ellA == 0 || ellB == 0)
    {
        nmod_poly_mat_zero(c);
        return;
    }

    m = a->r;
    n = a->c;
    p = b->c;

    if (c == a || c == b)
    {
        nmod_poly_mat_t T;
        nmod_poly_mat_init(T, m, n, p);
        nmod_poly_mat_middle_product_vdm1(T, a, b, dA, dB);
        nmod_poly_mat_swap_entrywise(c, T);
        nmod_poly_mat_clear(T);
        return;
    }
    
    nmod_init(&mod, a->modulus);
    vandermonde_init1(vA, vB, iv, dA, dB, mod);
    nb_points = dA + dB + 1;

    // evaluation of matrix a:
    // build tmp_mat, whose column ell = i*n + j is the coefficient vector 
    // of reverse(a[i][j]) 
    nmod_mat_init(tmp_mat, dA + 1, m * n, mod.n);
    ell = 0;
    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++, ell++)
        {
            ulong d;
            nn_ptr ptr;
            d = nmod_poly_degree(a->rows[i] + j);
            ptr = (a->rows[i] + j)->coeffs;
            for (k = 0; k <= d; k++)
                tmp_mat->rows[dA - k][ell] = ptr[k];
        }
    // note: d = deg(a[i][j]) is -1 if a[i][j] == 0
    // all non-touched entries already zero since tmp_mat was initialized as zero
   
    // valA: column ell = i*n + j contains the evaluations of a[i][j]
    nmod_mat_init(valA, vA->r, tmp_mat->c, mod.n);
    nmod_mat_mul_pml(valA, vA, tmp_mat);

    // transpose interpolation of matrix b:
    // build tmp_mat, whose column ell = i*n + j is the coefficient vector of
    // b[i][j] (padded with zeroes up to length dB+1 if necessary)
    nmod_mat_clear(tmp_mat);
    nmod_mat_init(tmp_mat, dA + dB + 1, n * p, mod.n);
    ell = 0;
    for (i = 0; i < n; i++)
        for ( j = 0; j < p; j++, ell++)
        {
            ulong d;
            nn_ptr ptr;
            d = nmod_poly_degree(b->rows[i] + j);
            ptr = (b->rows[i] + j)->coeffs;
            for (k = 0; k <= d; k++)
                tmp_mat->rows[k][ell] = ptr[k];
        }
    // mul by transpose(iv)
    nmod_mat_init(iv_t, iv->r, iv->c, mod.n);
    nmod_mat_transpose(iv_t, iv);
    nmod_mat_init(valB, iv_t->r, tmp_mat->c, mod.n);
    nmod_mat_mul_pml(valB, iv_t, tmp_mat);
 

    // perform the pointwise products
    nmod_mat_init(valAp, m, n, mod.n);
    nmod_mat_init(valBp, n, p, mod.n);
    nmod_mat_init(valC, nb_points, m * p, mod.n);
    nmod_mat_init(valCp, m, p, mod.n);
    
    for (i = 0; i < nb_points; i++)
    {
        // a evaluated at point i
        ell = 0;
        for (u = 0; u < m; u++)
            for (v = 0; v < n; v++, ell++)
                valAp->rows[u][v] = valA->rows[i][ell];

        ell = 0;
        for (u = 0; u < n; u++)
            for (v = 0; v < p; v++, ell++)
                valBp->rows[u][v] = valB->rows[i][ell];

        nmod_mat_mul_pml(valCp, valAp, valBp);

        // copy this into valC
        ell = 0;
        for (u = 0; u < m; u++)
            for (v = 0; v < p; v++, ell++)
                valC->rows[i][ell] = valCp->rows[u][v];
    }
    
    nmod_mat_init(vB_t, vB->c, vB->r, mod.n);
    nmod_mat_transpose(vB_t, vB);
    // transpose-evaluate to find the entries of c
    nmod_mat_clear(tmp_mat);
    nmod_mat_init(tmp_mat, vB_t->r, valC->c, mod.n);
    nmod_mat_mul_pml(tmp_mat, vB_t, valC);


    // copy to output (reorganize these entries into c)
    ell = 0;
    for (u = 0; u < m; u++)
        for (v = 0; v < p; v++, ell++)
        {
            nn_ptr coeffs;
            nmod_poly_realloc(c->rows[u] + v, dB + 1);
            coeffs = (c->rows[u] + v)->coeffs;
            (c->rows[u] + v)->length = dB + 1; 
            for (i = 0; i <= dB; i++)
                coeffs[i] = tmp_mat->rows[i][ell];
            _nmod_poly_normalise(c->rows[u] + v);
        }
    
    nmod_mat_clear(vA);
    nmod_mat_clear(vB_t);
    nmod_mat_clear(vB);
    nmod_mat_clear(iv_t);
    nmod_mat_clear(iv);
    nmod_mat_clear(tmp_mat);
    nmod_mat_clear(valA);
    nmod_mat_clear(valB);
    nmod_mat_clear(valAp);
    nmod_mat_clear(valBp);
    nmod_mat_clear(valC);
    nmod_mat_clear(valCp);
}




/*------------------------------------------------------------*/
/* computes the matrices for evaluation and interpolation     */
/*------------------------------------------------------------*/
static
void vandermonde_init2(nmod_mat_t vdm1, nmod_mat_t vdm2, nmod_mat_t inv_vdm, long d1, long d2, nmod_t mod)
{
    ulong p;
    ulong i, j, s1, s2, nb_points;
    nmod_mat_t vdm;
    
    // sizes (for even part)
    s1 = d1 / 2 + 1;
    s2 = d2 / 2 + 1;
    // nb points, such that 2*nb_points >= d1+d2+1
    nb_points = (d1 + d2) / 2 + 1;
    
    // vdm: square Vandermonde matrix with nb_points
    // points chosen as the squares of 1, 2, .., nb_points
    // --> row i contains 1, (i+1)^2, (i+1)^4, .., (i+1)^{2*nb_points-2}
    // vdm1: nb_points x s1 submatrix of vdm
    // vdm2: nb_points x s2 submatrix of vdm
    p = mod.n;
    nmod_mat_init(vdm1, nb_points, s1, p);
    nmod_mat_init(vdm2, nb_points, s2, p);
    nmod_mat_init(vdm, nb_points, nb_points, p);
    nmod_mat_init(inv_vdm, nb_points, nb_points, p);
    
    for (i = 0; i < nb_points; i++)
    {
        ulong p1, pt;
        NMOD_RED(p1, i+1, mod);
        pt = nmod_mul(p1, p1, mod);
        nmod_mat_entry(vdm, i, 0) = 1;
        nmod_mat_entry(vdm1, i, 0) = 1;
        nmod_mat_entry(vdm2, i, 0) = 1;
        for (j = 1; j < nb_points; j++)
        {
            nmod_mat_entry(vdm, i, j) = nmod_mul(nmod_mat_entry(vdm, i, j-1), pt, mod);
            if (j < s1)
                nmod_mat_entry(vdm1, i, j) = nmod_mat_entry(vdm, i, j);
            if (j < s2)
                nmod_mat_entry(vdm2, i, j) = nmod_mat_entry(vdm, i, j);
        }
    }

    // inv_vdm is the inverse of vdm
    nmod_mat_inv(inv_vdm, vdm);
    nmod_mat_clear(vdm);
}

/*------------------------------------------------------------*/
/* c = a*b                                                    */
/* output may alias input; c does not have to be zero matrix  */
/* matrix multiplication using the algorithm of Doliskani et al. */
/* uses matrix multiplication for evaluation and interpolation*/
/*------------------------------------------------------------*/
void nmod_poly_mat_mul_vdm2(nmod_poly_mat_t c, const nmod_poly_mat_t a, const nmod_poly_mat_t b)
{
    ulong ellA, ellB, dA, dB, sA, sB, m, n, p, ell, nb_points, i, j, k, u, v;
    nmod_mat_t vA, vB, iv, tmp_mat_even, tmp_mat_odd, valA_even, valA_odd, valB_even, valB_odd,
        valAp_even, valAp_odd, valBp_even, valBp_odd, valC_even, valC_odd, valCp_even, valCp_odd;
    nmod_t mod;
    ulong inv2, pt, inv2pt;
    
    ellA = nmod_poly_mat_max_length(a);
    ellB = nmod_poly_mat_max_length(b);
    if (ellA == 0 || ellB == 0)
    {
        nmod_poly_mat_zero(c);
        return;
    }
    dA = ellA - 1;
    dB = ellB - 1;
    sA = dA / 2 + 1;
    sB = dB / 2 + 1;
    
    m = a->r;
    n = a->c;
    p = b->c;

    if (c == a || c == b)
    {
        nmod_poly_mat_t T;
        nmod_poly_mat_init(T, m, n, p);
        nmod_poly_mat_mul_vdm2(T, a, b);
        nmod_poly_mat_swap_entrywise(c, T);
        nmod_poly_mat_clear(T);
        return;
    }

    nmod_init(&mod, a->modulus);
    vandermonde_init2(vA, vB, iv, dA, dB, mod);

    nb_points = vA->r;

    // evaluation of matrix a (even/odd part):
    // build tmp_mat_even/odd, whose column ell = i*n + j is the coefficient
    // vector of the even/odd part of a[i][j] (padded with zeroes up to length
    // dA/2+1 if necessary)
    nmod_mat_init(tmp_mat_even, sA, m*n, mod.n);
    nmod_mat_init(tmp_mat_odd, sA, m*n, mod.n);
    // --> we use dimension sA=dA/2+1 for odd part, even though this should be
    // (dA-1)/2+1; this is slightly non-optimal...
    ell = 0;
    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++, ell++)
        {
            ulong d;
            nn_ptr ptr;
            d = nmod_poly_degree(a->rows[i] + j);
            ptr = (a->rows[i] + j)->coeffs;
            for (k = 0; 2*k < d; k++)
            {
                tmp_mat_even->rows[k][ell] = ptr[2*k];
                tmp_mat_odd->rows[k][ell] = ptr[2*k+1];
            }
            if (d % 2 == 0) // d even
                tmp_mat_even->rows[d/2][ell] = ptr[d];
        }

    // note: d = deg(a[i][j]) is -1 if a[i][j] == 0
    // all non-touched entries already zero since tmp_mat was initialized as zero
    // valAeven: column ell = i*n + j contains the evaluations of even part of a[i][j]
    // valAodd: column ell = i*n + j contains the evaluations of odd part of a[i][j]

    nmod_mat_init(valA_even, vA->r, tmp_mat_even->c, mod.n);
    nmod_mat_init(valA_odd, vA->r, tmp_mat_odd->c, mod.n);
    
    nmod_mat_mul_pml(valA_even, vA, tmp_mat_even);
    nmod_mat_mul_pml(valA_odd, vA, tmp_mat_odd);

    // evaluation of matrix b (even/odd part): build tmp_mat_even/odd, whose
    // column ell = i*n + j is the coefficient vector of the even/odd part of
    // b[i][j] (padded with zeroes up to length dB/2+1 if necessary); as above,
    // we use dimension dB/2+1 although we know it should be (dB-1)/2+1)

    nmod_mat_clear(tmp_mat_even);
    nmod_mat_clear(tmp_mat_odd);
    nmod_mat_init(tmp_mat_even, sB, n * p, mod.n);
    nmod_mat_init(tmp_mat_odd, sB, n * p, mod.n);
    
    ell = 0;
    for (i = 0; i < n; i++)
        for (j = 0; j < p; j++, ell++)
        {
            ulong d;
            nn_ptr ptr;
            
            d = nmod_poly_degree(b->rows[i] + j);
            ptr = (b->rows[i] + j)->coeffs;
            k = 0;
            for (; 2*k < d; k++)
            {
                tmp_mat_even->rows[k][ell] = ptr[2*k];
                tmp_mat_odd->rows[k][ell] = ptr[2*k+1];
            }
            if (2*k == d)
            {
                tmp_mat_even->rows[k][ell] = ptr[d];
                tmp_mat_odd->rows[k][ell] = 0;
                k++;
            }
            for (; k < sB; k++)
            {
                tmp_mat_even->rows[k][ell] = 0;
                tmp_mat_odd->rows[k][ell] = 0;
            }
        }

    // valBeven/odd: column ell = i*n + j contains the evaluations of the
    // even/odd part of b[i][j]

    nmod_mat_init(valB_even, vB->r, tmp_mat_even->c, mod.n);
    nmod_mat_init(valB_odd, vB->r, tmp_mat_odd->c, mod.n);

    nmod_mat_mul_pml(valB_even, vB, tmp_mat_even);
    nmod_mat_mul_pml(valB_odd, vB, tmp_mat_odd);

    // perform the pointwise products for the points 1, 2, ..., nb_points
    nmod_mat_init(valAp_even, m, n, mod.n);
    nmod_mat_init(valAp_odd, m, n, mod.n);
    nmod_mat_init(valBp_even, n, p, mod.n);
    nmod_mat_init(valBp_odd, n, p, mod.n);
    nmod_mat_init(valC_even, nb_points, m * p, mod.n);
    nmod_mat_init(valC_odd, nb_points, m * p, mod.n);
    nmod_mat_init(valCp_even, valAp_even->r, valBp_even->c, mod.n);
    nmod_mat_init(valCp_odd, valAp_odd->r, valBp_odd->c, mod.n);
    
    // from now on we will only use valA_odd with row i multiplied by i+1
    // same with valB_odd
    for (i = 0; i < nb_points; i++)
    {
        ulong pt;
        NMOD_RED(pt, i+1, mod);
        for (j = 0; j < (ulong) valA_odd->c; j++)
            valA_odd->rows[i][j] = nmod_mul(pt, valA_odd->rows[i][j], mod);
        for (j = 0; j < (ulong) valB_odd->c; j++)
            valB_odd->rows[i][j] = nmod_mul(pt, valB_odd->rows[i][j], mod);
    }

    inv2 = nmod_inv(2, mod);

    for (i = 0; i < nb_points; i++)
    {
        // a evaluated at point i (which is i+1)
        ell = 0;
        for (u = 0; u < m; u++)
            for (v = 0; v < n; v++, ell++)
            {
                valAp_even->rows[u][v] = nmod_add(valA_even->rows[i][ell], valA_odd->rows[i][ell], mod);
                valAp_odd->rows[u][v] = nmod_sub(valA_even->rows[i][ell], valA_odd->rows[i][ell], mod);
            }

        // b evaluated at point i (which is i+1)
        ell = 0;
        for (u = 0; u < n; u++)
            for (v = 0; v < p; v++, ell++)
            {
                valBp_even->rows[u][v] = nmod_add(valB_even->rows[i][ell], valB_odd->rows[i][ell], mod);
                valBp_odd->rows[u][v] = nmod_sub(valB_even->rows[i][ell], valB_odd->rows[i][ell], mod);
            }

        // a*b evaluated at point i (which is i+1)
        nmod_mat_mul_pml(valCp_even, valAp_even, valBp_even);
        nmod_mat_mul_pml(valCp_odd, valAp_odd, valBp_odd);

        // copy this into valC: column ell = i*n + j contains the evaluations
        // of the entry i,j of c = a*b
        // valC_even = (_even + _odd)/2
        // valC_odd = (_even - _odd)/2x
        ell = 0;
        NMOD_RED(pt, 2*(i+1), mod);
        inv2pt = nmod_inv(pt, mod);

        for (u = 0; u < m; u++)
            for (v = 0; v < p; v++, ell++)
            {
                valC_even->rows[i][ell] = nmod_add(valCp_even->rows[u][v], valCp_odd->rows[u][v], mod);
                valC_odd->rows[i][ell] = nmod_sub(valCp_even->rows[u][v], valCp_odd->rows[u][v], mod);
            }

        for (j = 0; j < (ulong) valC_even->c; j++)
            valC_even->rows[i][j] = nmod_mul(valC_even->rows[i][j], inv2, mod);

        for (j = 0; j < (ulong) valC_odd->c; j++)
            valC_odd->rows[i][j] = nmod_mul(valC_odd->rows[i][j], inv2pt, mod);
    }

    // interpolate to find the _even/_odd parts of c
    nmod_mat_clear(tmp_mat_even);
    nmod_mat_clear(tmp_mat_odd);
    nmod_mat_init(tmp_mat_even, iv->r, valC_even->c, mod.n);
    nmod_mat_init(tmp_mat_odd, iv->r, valC_odd->c, mod.n);
    
    nmod_mat_mul_pml(tmp_mat_even, iv, valC_even);
    nmod_mat_mul_pml(tmp_mat_odd, iv, valC_odd);
    
    // copy to output (reorganize these entries into c)
    ell = 0;
    for (u = 0; u < m; u++)
        for (v = 0; v < p; v++, ell++)
        {
            nn_ptr coeffs;
            nmod_poly_realloc(c->rows[u] + v, 2 * nb_points);
            coeffs = (c->rows[u] + v)->coeffs;
            (c->rows[u] + v)->length = 2 * nb_points; 
            for (i = 0; i < nb_points; i++)
            {
                coeffs[2*i] = tmp_mat_even->rows[i][ell];
                coeffs[2*i+1] = tmp_mat_odd->rows[i][ell];
            }
            _nmod_poly_normalise(c->rows[u] + v);
        }

    nmod_mat_clear(valCp_even);
    nmod_mat_clear(valCp_odd);
    nmod_mat_clear(valAp_even);
    nmod_mat_clear(valAp_odd);
    nmod_mat_clear(valBp_even);
    nmod_mat_clear(valBp_odd);
    nmod_mat_clear(valC_even);
    nmod_mat_clear(valC_odd);
    nmod_mat_clear(valA_even);
    nmod_mat_clear(valA_odd);
    nmod_mat_clear(valB_even);
    nmod_mat_clear(valB_odd);
    nmod_mat_clear(tmp_mat_even);
    nmod_mat_clear(tmp_mat_odd);
    nmod_mat_clear(vA);
    nmod_mat_clear(vB);
    nmod_mat_clear(iv);
}

