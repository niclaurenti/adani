#include "/Users/niccololaurenti/Master-thesis/include/MyLib.h"
//#include "apfel/massivecoefficientfunctionsunp_sl.h"

#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;
//using namespace apfel;
double damp_thr(double eta,double mQ,double mMu,double A,double B,double C,double D,double a,double b);
double damp_asy(double eta,double mQ,double mMu,double A,double B,double C,double D,double a,double b);

int main(int argc, char** argv) {

	if(argc!=3) {
		cout<< "ERROR\nUsage: ./approx.exe xi filename\nExiting..." <<endl;
		return -1;
	}

    double xi = atof(argv[1]);
    //double D= atof(argv[2]);
	ofstream file_eta;
	file_eta.open(argv[2]);

	double ximu=1;

	double mQ=1/xi;
	double mMu=1/ximu;
	//int nf=3;

	double x;
    //double dx=1e-3 ;
    double a=2.5, b=5;
    double A=20., B=11., D=2., C=3.;

	double eta;

	/*
	for(x=dx; x<1; x+=dx) {

		file_x << x <<"   "
					 << CL_g2_asymptotic(x,mQ,mMu) << "   "
					 << CL_g2_threshold(x,mQ,mMu) << "   "
					 << CL_g2_approximation(x,mQ,mMu) << "   "
					 << endl;

	}
	*/
	double logeta, logeta_min=-4, logeta_max=4;

	int N=500;
	double dlog=(logeta_max - logeta_min)/N;

	for(int i=0; i<N; i++) {

		logeta=logeta_min + i*dlog;

		eta=pow(10, logeta);

		x=1/(1+4*mQ*(eta+1));

	    //double damp_thr=1/(1+pow(eta/h,k));
		//double damp_asy=1-damp_thr;

		file_eta << eta << "   "
                    << x*CL_g2(x,mQ,mMu) << "   "
                    << x*CL_g2_asymptotic(x,mQ,mMu) << "   "
					<< x*CL_g2_threshold(x,mQ,mMu) << "   "
					<< x*CL_g2_approximation(x,mQ,mMu,A,B,C,D,a,b) << "   "
					<< x*CL_ps2(x,mQ,mMu) << "   "
                    << x*CL_ps2_asymptotic(x,mQ,mMu) << "   "
					<< 0 << "   "
					<< x*CL_ps2_approximation(x,mQ,mMu,A,B,C,D,a,b) << "   "
                    << damp_thr(eta,mQ,mMu,A,B,C,D,a,b) << "   "
                    << damp_asy(eta,mQ,mMu,A,B,C,D,a,b) << "   "
					<< endl;

	}

	file_eta.close();

	return 0;

}

double damp_thr(double eta,double mQ,double mMu,double A,double B,double C,double D,double a,double b) {
    double xi = 1./mQ ;
    double h = A + (B-A)/(1+exp(a*(log(xi)-b)));
	double k = C + (D-C)/(1+exp(a*(log(xi)-b)));

	return 1./(1.+pow(eta/h,k));
}

double damp_asy(double eta,double mQ,double mMu,double A,double B,double C,double D,double a,double b) {
    return 1 - damp_thr(eta, mQ, mMu, A, B, C, D, a, b) ;
}
