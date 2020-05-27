#include "lzz_pX_seq.h"
#include <NTL/lzz_pX.h>
#include <iostream>
#include <sstream>
#include <NTL/vector.h>

NTL_CLIENT
using namespace std;

int main(){
	long p = 9001;
	zz_p::init(p);

	bool verbose = false;

	long d;
	if (verbose)
		cout << "Enter d: ";
	cin >> d;

	long ell;
	if (verbose)
		cout << "Enter ell: ";
	cin >> ell;

	cout << "d: " << d << endl;
	cout << "ell: " << ell << endl;

	Vec<zz_pX> seq;
	seq.SetLength(ell);
	for (long i = 0; i < ell; i++){
		if (verbose)
			cout << "enter seq[" << i << "]:"<< endl;
		cin >> seq[i];
		trunc(seq[i],seq[i],d);
	}
	if(verbose){
		cout << "SEQ: " << endl;
		for (long i = 0; i < ell; i++)
			cout << seq[i] << endl;
	}
	cout << endl;
	cout << "***************************************" << endl;
	cout << "STARTING KURAKIN" << endl;
	Vec<zz_pXY> gens;
	double time = GetWallTime();
	kurakin(d,seq,gens);
	cout << "Monic degree: " << gens[0].degY() << endl;
	cout << "Took: " << GetWallTime()-time << endl;

	if(verbose){
		cout << "GENS:" << endl;
		for (long i = 0; i < gens.length(); i++)
			cout << gens[i] << endl;
		cout << endl;
	}

	cout << "***************************************" << endl;
	cout << "STARTING MODIFIED KURAKIN" << endl;
	Vec<zz_pXY> gens2;
	time = GetWallTime();
	modified_kurakin(d,seq,gens2);	
	cout << "Monic degree: " << gens2[0].degY() << endl;
	cout << "Took: " << GetWallTime()-time << endl;

	if(verbose){
		cout << "GENS2:"<< endl;
		for (long i = 0; i < gens2.length(); i++)
			cout << gens2[i] << endl;
	}

	long d_s = 0;
	for (long i = 0; i < gens2.length(); i++)
		if (gens2[i] != zz_pX{0}) d_s++;
	cout << "d*: " << d_s << endl;

	long d_opt = 1;
	for (long i = 1 ; i < gens.length(); i++)
		if (gens[i].degY() != gens[i-1].degY()) d_opt++;
	cout << "d_opt: " << d_opt << endl;

	bool do_wied = true;
	for (long i = 1; i < gens2.length(); i++)
		if (!gens2[i].is_zero()) do_wied = false;

	if(verbose){
		cout << "GENS2 filled in:" << endl;
		fill_in(d,gens2);
		for (long i  = 0; i < gens2.length(); i++)
			cout << gens2[i] << endl;
		cout << endl;
	}

	cout << "***************************************" << endl;
	cout << "STARTING CHECKS" << endl;
	cout << "checking that gens cancels the seq" << endl;
	if (check_cancel(seq,gens,d)) cout << "OKAY!" << endl;
	else cout << "BAD!" << endl;

	cout << "checking that gens2 cancels the seq" << endl;
	if (check_cancel(seq,gens2,d)) cout << "OKAY!" << endl;
	else cout << "BAD!" << endl;

	if (do_wied){
		cout << "***************************************" << endl;
		cout << "Detected non-degenerate sequence." << endl;
		cout << "STARTING STRUCTURED LIFTING" << endl;
		zz_pXY P;
		time = GetWallTime();
		minpoly_nondegenerate(d,seq,P);
		cout << "degree: " << P.degY() << endl;
		cout << "Took: " << GetWallTime() - time << endl;
		if(verbose){
			cout << "P: " << P << endl;
		}
		Vec<zz_pXY> gens3;
		cout << "checking that P cancels the seq" << endl;
		gens3.append(P);
		if (check_cancel(seq,gens3,d)) cout << "OKAY!" << endl;
		else cout << "BAD!" << endl;
	}


}
