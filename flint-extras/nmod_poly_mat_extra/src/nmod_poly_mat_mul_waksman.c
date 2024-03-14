#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>

void nmod_poly_mat_mul_waksman(nmod_poly_mat_t C, const nmod_poly_mat_t A,  const nmod_poly_mat_t B)
{
    ulong m, n, p, np, i, l, j, j2, k;
    mp_limb_t mod, half;
    nmod_poly_t val0, val1, val2, crow;
    nmod_poly_struct *Crow, *Ccol;

    m = A->r;
    n = B->r;
    p = B->c;
    mod = A->modulus;
    
    if (m < 1 || n < 1 || p < 1)
    {
        nmod_poly_mat_zero(C);
        return;
    }

    if (C == A || C == B)
    {
        nmod_poly_mat_t T;
        nmod_poly_mat_init(T, m, n, mod);
        nmod_poly_mat_mul_waksman(T, A, B);
        nmod_poly_mat_swap_entrywise(C, T);
        nmod_poly_mat_clear(T);
        return;
    }

    half = (mod+1) >> 1; // 1/2 

    nmod_poly_init(val0, mod);
    nmod_poly_init(val1, mod);
    nmod_poly_init(val2, mod);
    nmod_poly_init(crow, mod);

    Crow = flint_malloc(p * sizeof(nmod_poly_struct));
    Ccol = flint_malloc(m * sizeof(nmod_poly_struct));

    for (i = 0; i < p; i++)
        nmod_poly_init(Crow +i, mod);

    for (i = 0; i < m; i++)
        nmod_poly_init(Ccol+i, mod);

    for (i = 0; i < m; i++)
        for (j = 0; j < p; j++)
            nmod_poly_zero(nmod_poly_mat_entry(C, i, j));

    np = n >> 1;

    for (j = 1; j <= np; j++)
    {
        j2 = (j << 1) - 1;
    
        for (k = 0; k < p; k++)
        {
            nmod_poly_add(val1, nmod_poly_mat_entry(A, 0, j2-1), nmod_poly_mat_entry(B, j2, k));
            nmod_poly_add(val2, nmod_poly_mat_entry(A, 0, j2), nmod_poly_mat_entry(B, j2-1, k));
            nmod_poly_mul(val1, val1, val2);
            nmod_poly_add(nmod_poly_mat_entry(C, 0, k), nmod_poly_mat_entry(C, 0, k), val1);
            
            nmod_poly_sub(val1, nmod_poly_mat_entry(A, 0, j2-1), nmod_poly_mat_entry(B, j2, k));
            nmod_poly_sub(val2, nmod_poly_mat_entry(A, 0, j2), nmod_poly_mat_entry(B, j2-1, k));
            nmod_poly_mul(val1, val1, val2);
            nmod_poly_add(Crow + k, Crow + k, val1);
        }

        for (l = 1; l < m; l++)
        {
            nmod_poly_add(val1, nmod_poly_mat_entry(A, l, j2-1), nmod_poly_mat_entry(B, j2, 0));
            nmod_poly_add(val2, nmod_poly_mat_entry(A, l, j2), nmod_poly_mat_entry(B, j2-1, 0));
            nmod_poly_mul(val1, val1, val2);
            nmod_poly_add(nmod_poly_mat_entry(C, l, 0), nmod_poly_mat_entry(C, l, 0), val1);
            
            nmod_poly_sub(val1, nmod_poly_mat_entry(A, l, j2-1), nmod_poly_mat_entry(B, j2, 0));
            nmod_poly_sub(val2, nmod_poly_mat_entry(A, l, j2), nmod_poly_mat_entry(B, j2-1, 0));
            nmod_poly_mul(val1, val1, val2);
            nmod_poly_add(Ccol + l, Ccol + l, val1);
        }
        
        for (k = 1; k < p; k++)
        {
            for (l = 1; l < m; l++)
            {
                nmod_poly_add(val1, nmod_poly_mat_entry(A, l, j2-1), nmod_poly_mat_entry(B, j2, k));
                nmod_poly_add(val2, nmod_poly_mat_entry(A, l, j2), nmod_poly_mat_entry(B, j2-1, k));
                nmod_poly_mul(val1, val1, val2);
                nmod_poly_add(nmod_poly_mat_entry(C, l, k), nmod_poly_mat_entry(C, l, k), val1);
            }
        }
    }
    
    for (l=1; l<m; l++)
    {
        nmod_poly_add(val1, Ccol + l, nmod_poly_mat_entry(C, l, 0));
        nmod_poly_scalar_mul_nmod(Ccol+ l, val1, half);
        nmod_poly_sub(nmod_poly_mat_entry(C, l, 0), nmod_poly_mat_entry(C, l, 0), Ccol + l);
    }

    nmod_poly_add(val1, Crow, nmod_poly_mat_entry(C, 0, 0));
    nmod_poly_scalar_mul_nmod(val0, val1, half);
    nmod_poly_sub(nmod_poly_mat_entry(C, 0, 0), nmod_poly_mat_entry(C, 0, 0), val0);
    
    for (k = 1; k < p; k++)
    {
        nmod_poly_add(crow, Crow + k, nmod_poly_mat_entry(C, 0, k));
        nmod_poly_scalar_mul_nmod(val1, crow, half);
        nmod_poly_sub(nmod_poly_mat_entry(C, 0, k), nmod_poly_mat_entry(C, 0, k), val1);
        nmod_poly_sub(crow, val1, val0);
        
        for (l = 1; l < m; l++)
        {
            nmod_poly_sub(val2, nmod_poly_mat_entry(C, l, k), crow);
            nmod_poly_sub(nmod_poly_mat_entry(C, l, k), val2, Ccol + l);
        }
    }
    
    if (n & 1)
    {
        for (l = 0; l < m; l++)
            for (k = 0; k < p; k++)
            {
                nmod_poly_mul(val1, nmod_poly_mat_entry(A, l, n-1), nmod_poly_mat_entry(B, n-1, k));
                nmod_poly_add(nmod_poly_mat_entry(C, l, k), nmod_poly_mat_entry(C, l, k), val1);
            }
    }

    for (i = 0; i < p; i++)
        nmod_poly_clear(Crow+i);
    
    for (i = 0; i < m; i++)
        nmod_poly_clear(Ccol+i);
    
    flint_free(Crow);
    flint_free(Ccol);
    
    nmod_poly_clear(val0);
    nmod_poly_clear(val1);
    nmod_poly_clear(val2);
    nmod_poly_clear(crow);
}


