#include <flint/fmpz.h>
#include <flint/fmpz_vec.h>
#include <flint/fmpz_mat.h>


/** ------------------------------------------------------------ */
/** Waksman's algorithm for matrix multiplication                */
/** does n^3/2+O(n^2) products, but many additions               */
/** good for small matrices with large entries                   */
/** ------------------------------------------------------------ */
void fmpz_mat_mul_waksman(fmpz_mat_t C, const fmpz_mat_t A, const fmpz_mat_t B)
{

    long m = A->r;
    long n = B->r;
    long p = B->c;

    long i, l, j, k;

    fmpz * Crow = _fmpz_vec_init(p);
    fmpz * Ccol = _fmpz_vec_init(m);

    long np = n >> 1;

    fmpz_t val0, val1, val2, crow;
    fmpz_init(val0);
    fmpz_init(val1);
    fmpz_init(val2);
    fmpz_init(crow);

    for (i = 0; i < p; i++)
        fmpz_zero(Crow + i);

    for (i = 0; i < m; i++)
        fmpz_zero(Ccol + i);

    for (i = 0; i < m; i++)
        for (j = 0; j < p; j++)
            fmpz_zero(fmpz_mat_entry(C, i, j));

    for (j = 1; j <= np; j++)
    {
        long j2 = (j << 1) - 1;
    
        for (k = 0; k < p; k++)
        {
            fmpz_add(val1, fmpz_mat_entry(A, 0, j2-1), fmpz_mat_entry(B, j2, k));
            fmpz_add(val2, fmpz_mat_entry(A, 0, j2), fmpz_mat_entry(B, j2-1, k));
            fmpz_addmul(fmpz_mat_entry(C, 0, k), val1, val2);

            fmpz_sub(val1, fmpz_mat_entry(A, 0, j2-1), fmpz_mat_entry(B, j2, k));
            fmpz_sub(val2, fmpz_mat_entry(A, 0, j2), fmpz_mat_entry(B, j2-1, k));
            fmpz_addmul(Crow + k, val1, val2);
        }

        for (l = 1; l < m; l++)
        {
            fmpz_add(val1, fmpz_mat_entry(A, l, j2-1), fmpz_mat_entry(B, j2, 0));
            fmpz_add(val2, fmpz_mat_entry(A, l, j2), fmpz_mat_entry(B, j2-1, 0));
            fmpz_addmul(fmpz_mat_entry(C, l, 0), val1, val2);
      
            fmpz_sub(val1, fmpz_mat_entry(A, l, j2-1), fmpz_mat_entry(B, j2, 0));
            fmpz_sub(val2, fmpz_mat_entry(A, l, j2), fmpz_mat_entry(B, j2-1, 0));
            fmpz_addmul(Ccol + l, val1, val2);
        }

        for (k = 1; k < p; k++)
        {
            for (l = 1; l < m; l++)
            {
                fmpz_add(val1, fmpz_mat_entry(A, l, j2-1), fmpz_mat_entry(B, j2, k));
                fmpz_add(val2, fmpz_mat_entry(A, l, j2), fmpz_mat_entry(B, j2-1, k));
                fmpz_addmul(fmpz_mat_entry(C, l, k), val1, val2);
            }
        }
    }

    for (l = 1; l < m; l++)
    {
        fmpz_add(val1, Ccol + l, fmpz_mat_entry(C, l, 0));
        fmpz_fdiv_q_2exp(Ccol+ l, val1, 1); 
        fmpz_sub(fmpz_mat_entry(C, l, 0), fmpz_mat_entry(C, l, 0), Ccol + l);
    }

    fmpz_add(val1, Crow, fmpz_mat_entry(C, 0, 0));
    fmpz_fdiv_q_2exp(val0, val1, 1);
    fmpz_sub(fmpz_mat_entry(C, 0, 0), fmpz_mat_entry(C, 0, 0), val0);

    for (k = 1; k < p; k++)
    {
        fmpz_add(crow, Crow + k, fmpz_mat_entry(C, 0, k));
        fmpz_fdiv_q_2exp(val1, crow, 1); 
        fmpz_sub(fmpz_mat_entry(C, 0, k), fmpz_mat_entry(C, 0, k), val1);
        fmpz_sub(crow, val1, val0);

        for (l = 1; l < m; l++)
        {
            fmpz_sub(val2, fmpz_mat_entry(C, l, k), crow);
            fmpz_sub(fmpz_mat_entry(C, l, k), val2, Ccol + l);
        }
    }

    if ((n & 1) == 1)
        for (l = 0; l < m; l++)
            for (k = 0; k < p; k++)
                fmpz_addmul(fmpz_mat_entry(C, l, k), fmpz_mat_entry(A, l, n-1), fmpz_mat_entry(B, n-1, k));
  
    _fmpz_vec_clear(Crow, p);
    _fmpz_vec_clear(Ccol, m);

    fmpz_clear(val0);
    fmpz_clear(val1);
    fmpz_clear(val2);
    fmpz_clear(crow);
}
