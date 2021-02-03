#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include "lzz_pX_seq.h"

NTL_CLIENT
using namespace std;

bool comp(long i, long j){return i > j;}

int main(){
	bool verbose = false;
	long p = 9001;
	zz_p::init(p);

	long d;
	if(verbose)
		cout << "Enter d: " <<endl;
	cin >> d;

	long ell;
	if (verbose)
		cout << "Enter ell: " << endl;
	cin >> ell;

	if(verbose){
		cout << "d: " << d << endl;
		cout << "ell: " << ell << endl;
	}

	if (ell % 2 != 0){
		cout << "ell must be even" << endl;
		return 1;
	}

	int n = 0;
	if (verbose)
		cout << "Enter n: " << endl;
	cin >> n;

	long max_ord = ell/2;
	srandom(time(NULL));

	vector<long> ords;
	for (long i = 0; i < d-1; i++)
		ords.emplace_back(rand() % max_ord + 1);
	ords.emplace_back(max_ord);
	sort(ords.begin(), ords.end(),comp); 

	Vec<zz_pX> Ps;
	for (long i = 0; i < d; i++){
		long curOrd = ords[d-1-i];
		zz_pX P = random_zz_pX(curOrd);
		MakeMonic(P);
		Ps.append(P);
	}
	Vec<Vec<zz_pX>> Ss;
	for (long i = 0; i < n; i++){
		Vec<zz_pX> S;
		S.SetLength(ell);
		for (long i = 0; i < d; i++){
			Vec<zz_pX> curS;
			long curOrd = ords[d-1-i];
			zz_pX P = Ps[i];
			MakeMonic(P);
			for (long j = 0; j < curOrd; j++)
				curS.append(zz_pX(random_zz_p()));
			for (long s = 0; s < ell-curOrd; s++){
				zz_p r{0};
				for (long t = 0; t < curOrd; t++)
					if (curS[s+t] != zz_pX{0})
						r = r - P[t] * curS[s+t][0];
				curS.append(zz_pX(r));
			}
			zz_pX curPow;
			SetCoeff(curPow,i,1);
			curS = mul(curS, curPow);
			S = add(S,curS);
		}
		Ss.append(S);
	}
	cout << d << endl;
	cout << ell << endl;
	cout << n << endl;

	for (long i = 0; i < ell; i++){
		for (long j = 0; j < n; j++)
			cout << Ss[j][i] << endl;
	}

	//for (long i = 0; i < ell; i++)
	//	cout << S[i] << endl;
	if (verbose){
		cout << "orders: [";
		for (auto &i : ords)
			cout << " " << i;
		cout << "]" << endl;
	}

}















