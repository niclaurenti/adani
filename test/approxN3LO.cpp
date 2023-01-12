#include "/Users/niccololaurenti/Master-thesis/include/MyLib.h"
//#include "apfel/massivecoefficientfunctionsunp_sl.h"

#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;
//using namespace apfel;

int main(int argc, char** argv) {

	if(argc!=4) {
		cout<< "ERROR\nUsage: ./approx.exe xi filename filename_x\nExiting..." <<endl;
		return -1;
	}

	ofstream file_eta;
	file_eta.open(argv[2]);

	ofstream file_x;
	file_x.open(argv[3]);

	double xi = atof(argv[1]);

	double ximu=1;

	double mQ=1/xi;
	double mMu=1/ximu;

	double x, dx=1e-3;

	double eta;

	int nf=3;


	for(x=dx; x<1; x+=dx) {

		file_x << x <<"   "
					 << CLm_g3_asymptotic(x,mQ,mMu,nf) << "   "
					 << CLm_g3_threshold(x,mQ,mMu,nf) << "   "
					 << CLm_g30_approximation(x,mQ,mMu,nf) << "   "
					 << endl;

	}

	double logeta, logeta_min=-4, logeta_max=4;

	int N=500;
	double dlog=(logeta_max - logeta_min)/N;

	for(int i=0; i<N; i++) {

		logeta=logeta_min + i*dlog;

		eta=pow(10, logeta);

		x=1/(1+4*mQ*(eta+1));

		file_eta << eta <<"   "
					 << x*CLm_g3_asymptotic(x,mQ,mMu,nf) << "   "
					 << x*CLm_g3_threshold(x,mQ,mMu,nf) << "   "
					 << x*CLm_g30_approximation(x,mQ,mMu,nf) << "   "
					 << x*CLm_ps3_asymptotic(x,mQ,mMu,nf) << "   "
					 << 0 << "   "
					 << x*CLm_ps30_approximation(x,mQ,mMu,nf) << "   "
					 << endl;

	}


	file_x.close();
	file_eta.close();

	return 0;

}
