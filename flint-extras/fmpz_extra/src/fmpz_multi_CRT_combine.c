#include <flint/fmpz.h>

/*-------------------------------------------------------------------*/
/* similar to fmpz_multi_CRT_precomp, except that we do not multiply */
/* by cofactors at the leaves                                        */
/* result is still reduced modulo the product of moduli              */
/*-------------------------------------------------------------------*/
void _fmpz_multi_CRT_combine(
    fmpz * outputs,
    const fmpz_multi_CRT_t P,
    fmpz * inputs)
{
    slong i, a, b, c;
    slong len = P->length;
    /* const fmpz * m = P->moduli; */
    fmpz * A, * B, * C, * t3, * t4;

    /* t1 = outputs + P->temp1loc; */
    /* t2 = outputs + P->temp2loc; */
    t3 = outputs + P->temp3loc;
    t4 = outputs + P->temp4loc;

    for (i = 0; i < len; i++)
    {
        a = P->prog[i].a_idx;
        b = P->prog[i].b_idx;
        c = P->prog[i].c_idx;

        A = outputs + a;
        B = outputs + b;
        C = outputs + c;

        if (b < 0)
        {
            B = inputs + (-b - 1);
            /* b = -b - 1; */
            /* B = t1; */
            /* fmpz_mul(t3, inputs + b, mf + b); */
            /* _fmpz_smod(B, t3, m + b, sign, t4); */
        }

        if (c < 0)
        {
            C = inputs + (-c - 1);
            /* c = -c - 1; */
            /* C = t2; */
            /* fmpz_mul(t3, inputs + c, mf + c); */
            /* _fmpz_smod(C, t3, m + c, sign, t4); */
        }

        /* A = B*c_m + C*b_m */
        fmpz_mul(A, B, P->prog[i].c_modulus);
        fmpz_mul(t3, C, P->prog[i].b_modulus);
        fmpz_add(A, A, t3);
    }

    _fmpz_smod(outputs + 0, A, P->final_modulus, 0, t4);
}


void fmpz_multi_CRT_combine(
    fmpz_t output,
    const fmpz_multi_CRT_t P,
    fmpz * inputs)
{
    slong i;
    fmpz * out;

    if (P->moduli_count == 1)
    {
        fmpz_set(output, inputs);
        return;
    }

    TMP_INIT;
    TMP_START;
    out = TMP_ARRAY_ALLOC(P->localsize, fmpz);
    for (i = 0; i < P->localsize; i++)
        fmpz_init(out + i);

    fmpz_swap(out + 0, output);
    _fmpz_multi_CRT_combine(out, P, inputs);
    fmpz_swap(out + 0, output);

    for (i = 0; i < P->localsize; i++)
        fmpz_clear(out + i);

    TMP_END;
}

