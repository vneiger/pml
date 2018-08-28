#include <NTL/ZZX.h>
#include <NTL/lzz_pX.h>
#include <NTL/vector.h>
#include <NTL/matrix.h>

#include "magma_output.h"

NTL_CLIENT


/*------------------------------------------------------------*/
/* initializes GF(p)                                          */
/*------------------------------------------------------------*/
void magma_init()
{
    cout << "k := GF(" << zz_p::modulus() << ");\n";
}

/*------------------------------------------------------------*/
/* initializes U<x>=GF(p)[x]                                  */
/*------------------------------------------------------------*/
void magma_init_X()
{
    cout << "U<x> := PolynomialRing(k);\n";
}

/*------------------------------------------------------------*/
/* initializes U<XX>=Z[XX]                                    */
/*------------------------------------------------------------*/
void magma_init_ZZX()
{
    cout << "U<XX> := PolynomialRing(Integers());\n";
}


/*------------------------------------------------------------*/
/* prints a vector                                            */
/*------------------------------------------------------------*/
void magma_output(const Vec<ZZ> & v)
{
    if (v.length() == 0)
    {
        cout << "[]";
        return;
    }
    cout << "[";
    for (long i = 0; i < v.length()-1; i++)
    {
        cout << v[i] << ", ";
    }
    cout << v[v.length()-1] << "]";
}

void magma_output(const Vec<long> & v)
{
    if (v.length() == 0)
    {
        cout << "[]";
        return;
    }
    cout << "[";
    for (long i = 0; i < v.length()-1; i++)
    {
        cout << v[i] << ", ";
    }
    cout << v[v.length()-1] << "]";
}

void magma_output(const Vec<zz_p> & v)
{
    if (v.length() == 0)
    {
        cout << "[]";
        return;
    }
    cout << "[";
    for (long i = 0; i < v.length()-1; i++)
    {
        cout << v[i] << ", ";
    }
    cout << v[v.length()-1] << "]";
}

void magma_output(const Vec<unsigned long> & v)
{
    if (v.length() == 0)
    {
        cout << "[]";
        return;
    }
    cout << "[";
    for (long i = 0; i < v.length()-1; i++)
    {
        cout << v[i] << ", ";
    }
    cout << v[v.length()-1] << "]";
}

/*------------------------------------------------------------*/
/* assigns a vector to variable "name"                        */
/*------------------------------------------------------------*/
void magma_assign(const Vec<ZZ> & v, const string & name)
{
    cout << name << " := ";
    magma_output(v);
    cout << ";" << endl;
}

void magma_assign(const Vec<long> & v, const string & name)
{
    cout << name << " := ";
    magma_output(v);
    cout << ";" << endl;
}

void magma_assign(const Vec<unsigned long> & v, const string & name)
{
    cout << name << " := ";
    magma_output(v);
    cout << ";" << endl;
}

void magma_assign(const Vec<zz_p> & v, const string & name)
{
    cout << name << " := ";
    magma_output(v);
    cout << ";" << endl;
}


/*------------------------------------------------------------*/
/* prints a matrix                                            */
/*------------------------------------------------------------*/
void magma_output(const Mat<zz_p> & v)
{
    if (v.NumRows() == 0)
    {
        cout << "Matrix(GF(" << zz_p::modulus() << "), [[]])";
        return;
    }
    cout << "Matrix(GF(" << zz_p::modulus() << "), [";
    for (long i = 0; i < v.NumRows()-1; i++)
    {
        magma_output(v[i]);
        cout << ", ";
    }
    magma_output(v[v.NumRows()-1]);
    cout << "])";
}

/*------------------------------------------------------------*/
/* assigns a matrix to variable "name"                        */
/*------------------------------------------------------------*/
void magma_assign(const Mat<zz_p> & v, const string & name)
{
    cout << name << " := ";
    magma_output(v);
    cout << ";" << endl;
}


/*------------------------------------------------------------*/
/* prints a poly / vector / .. with indeterminate var         */
/*------------------------------------------------------------*/
void magma_output(const ZZX & v, const string & var)
{
    cout << "(Parent(" << var << ")!(0)";
    for (long i = 0; i <= deg(v); i++)
    {
        cout << "+(" << coeff(v, i) << ")*" << var << "^" << i;
    }
    cout << ")";
}

void magma_output(const zz_pX & v, const string & var)
{
    cout << "(Parent(" << var << ")!(0)";
    for (long i = 0; i <= deg(v); i++)
    {
        cout << "+(" << coeff(v, i) << ")*" << var << "^" << i;
    }
    cout << ")";
}

void magma_output(const Vec<zz_pX> & v, const string & var)
{
    cout << "[";
    for (long i = 0; i < v.length(); i++)
    {
        magma_output(v[i]);
        if (i < v.length()-1)
            cout << ", ";
    }
    cout << "]";
}

void magma_output(const Mat<zz_pX>& a, const string & var)
{
    cout << "Matrix([";
    for (long i = 0; i < a.NumRows(); i++)
    {
        cout << "[";
        for (long j = 0; j < a.NumCols(); j++)
        {
            magma_output(a[i][j], var);
            if (j < a.NumCols()-1)
            {
                cout << ", ";
            }
        }
        cout << "]";
        if (i < a.NumRows()-1)
        {
            cout << ", ";
        }
    }
    cout << "])";
}



/*------------------------------------------------------------*/
/* prints a poly and cast it                                  */
/*------------------------------------------------------------*/
void magma_output(const ZZX & v)
{
    cout << "UZ!(";
    magma_output(v, "XX");
    cout << ")";
}

void magma_output(const zz_pX & v)
{
    cout << "U!(";
    magma_output(v.rep);
    cout << ")";
}

void magma_output(const Vec<zz_pX> & v)
{
    magma_output(v, "x");
}

void magma_output(const Mat<zz_pX> & v)
{
    magma_output(v, "x");
}


/*------------------------------------------------------------*/
/* assign a poly with indeterminate var to variable "name"    */
/*------------------------------------------------------------*/
void magma_assign(const ZZX & v, const string & var, const string & name)
{
    cout << name << " := ";
    magma_output(v, var);
    cout << ";" << endl;
} 

void magma_assign(const zz_pX & v, const string & var, const string & name)
{
    cout << name << " := ";
    magma_output(v, var);
    cout << ";" << endl;
} 

void magma_assign(const Vec<zz_pX> & v, const string & var, const string & name)
{
    cout << name << ":=";
    magma_output(v, var);
    cout << ";" << endl;
}

void magma_assign(const Mat<zz_pX>& a, const string & var, const string & name)
{
    cout << name << ":=";
    magma_output(a, var);
    cout << ";" << endl;
}

/*------------------------------------------------------------*/
/* assign a poly with indeterminate XX to variable "name"     */
/*------------------------------------------------------------*/
void magma_assign(const ZZX & v, const string & name)
{
    magma_assign(v, "XX", name);
}

void magma_assign(const zz_pX & v, const string & name)
{
    magma_assign(v, "x", name);
}

void magma_assign(const Vec<zz_pX> & v, const string & name)
{
    magma_assign(v, "x", name);
}

void magma_assign(const Mat<zz_pX>& v, const string & name)
{
    magma_assign(v, "x", name);
}


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
