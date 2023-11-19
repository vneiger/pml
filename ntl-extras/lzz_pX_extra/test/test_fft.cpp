#include <NTL/lzz_pX.h>
#include <NTL/ZZ.h>
#include "lzz_pX_extra.h"
#include "mat_lzz_pX_multiply.h"
#include <iostream>
NTL_CLIENT

int main(int argc, char *argv[])
{
	long p = 3;
	if (argc >= 2)
	{
		p = atoi(argv[1]);
	}
	zz_p::FFTInit(0);
	SetSeed(ZZ{0});

	zz_pX_FFT fft(zz_p::modulus(),20);
	cout << "computed root of unity" << endl;

	// construct vector
	Vec<zz_p> v = random_vec_zz_p(pow(2,p-1)+1);
	v.SetLength(pow(2,p));
	for (long i = pow(2,p-1)+1; i < pow(2,p); i++)
		v[i] = zz_p{0};
	Vec<zz_p> v2 = v;
	if (p <= 4)
		cout << "vector: " << v << endl << endl;

	// run fft
	cout << "Running forward transform" << endl;
	Vec<zz_p> out;
	double t = GetWallTime();
	fft.forward(out,v);
	cout << "Took: " << GetWallTime() - t << endl;
	if (p <= 4)
		cout << out << endl << endl;
	
	// run ifft
	cout << "Running inverse" << endl;
	Vec<zz_p> out2;
	t = GetWallTime();
	fft.inverse(out2,out);
	cout << "Took: " << GetWallTime() - t << endl;

	if (p <= 4)
		cout << out2 << endl;
	cout << "inverse check: ";
	if (out2 == v) cout << "okay!";
	else cout << "BAD!";
	cout << endl<<endl;

	// set vector to be 2^(p-1)+1
	v.SetLength(pow(2,p-1)+1);
	if (p <= 4)
		cout << "truncated v: " << v << endl<<endl;

	// run fft
	cout << "Running truncated forward transform" << endl;
	t = GetWallTime();
	fft.forward(out,v);
	cout << "Took: " << GetWallTime() - t << endl;
	if (p <= 4)
		cout << out << endl << endl;
	
	cout << "Running truncated inverse" << endl;
	t = GetWallTime();
	fft.inverse(out2,out);
	cout << "Took: " << GetWallTime() - t << endl;

	if (p <= 4)
		cout << out2 << endl;
	cout << "inverse check: ";
	if (out2 == v) cout << "okay!";
	else cout << "BAD!";
	cout << endl<<endl;

	// run fft_t
	cout << "Running forward transpose" << endl;
	t = GetWallTime();
	fft.forward_t(out,v2);
	cout << "Took: " << GetWallTime() - t << endl;

	if (p <= 4)
		cout << out << endl << endl;
	
	cout << "Running inverse transpose" << endl;
	t = GetWallTime();
	fft.inverse_t(out2,out);
	cout << "Took: " << GetWallTime() - t << endl;

	if (p <= 4)
		cout << out2 << endl;
	cout << "inverse check: ";
	if (out2 == v2) cout << "okay!";
	else cout << "BAD!";
	cout << endl<<endl;

	cout << "Running truncated forward transpose" << endl;
	t = GetWallTime();
	fft.forward_t(out,v);
	cout << "Took: " << GetWallTime() - t << endl;

	if (p <= 4)
		cout << out << endl << endl;
	
	cout << "Running truncated inverse transpose" << endl;
	t = GetWallTime();
	fft.inverse_t(out2,out);
	cout << "Took: " << GetWallTime() - t << endl;

	if (p <= 4)
		cout << out2 << endl;
	cout << "inverse check: ";
	if (out2 == v) cout << "okay!";
	else cout << "BAD!";
	cout << endl<<endl;

	zz_pX f,g,h;
	random(f,pow(2,p-1)+1);
	random(g,pow(2,p-1)+1);
	t = GetWallTime();
	h=f*g;
	cout << "mult took: " << GetWallTime() - t << endl;
	
	zz_pX h2;
	t = GetWallTime();
	fft.mult(h2,f,g);
	cout << "fft mult took: " << GetWallTime()-t << endl;
	cout << "mult check: ";
	if (h2 == h) cout << "okay!";
	else cout << "BAD!";
	cout << endl << endl;

	random(f,pow(2,p-1)+1);
	random(g,pow(2,p)+1);
	auto q = f*g;
	Vec<zz_p> cs;
	cs.SetLength(pow(2,p-1)+1);
	for (long i = 0; i <= pow(2,p-1); i++)
		cs[i] = q[i+pow(2,p-1)];
	t = GetWallTime();
	fft.middle_prod(h,f,g);
	cout << "midprod took: " << GetWallTime() - t << endl;
	cout << "midprod check: ";
	for (long i = 0; i <= pow(2,p-1); i++)
		if (cs[i] != h[i])
		{
			cout << "BAD!" << endl;
			return 0;
		}
	cout << "okay!" << endl << endl;

	Mat<zz_pX> A,B,C,D;
	random(A,3,3,3);
	random(B,3,3,3);
	fft.mat_mult(C,A,B);
	multiply(D,A,B);
	cout << "Mat mult check: ";
	if (C == D) cout << "okay!" << endl;
	else cout << "BAD!" << endl;

	random(A,2,2,3);
	random(B,2,2,6);
	fft.mat_mp(C,A,B);
	cout << "A: " << A << endl;
	cout << "B: " << B << endl;
	cout << "C: " << C << endl;
	multiply(D,A,B);
	cout << "D: " << D << endl;
	f = A[0][0];
	g = B[0][0];
	fft.middle_prod(h,f,g);
	cout << "h: " << h << endl;
}










