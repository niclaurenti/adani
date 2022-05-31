#include "../include/ApproximateCoefficientFunctions.h"
#include "../include/HighScaleCoefficientFunctions.h"
#include "../include/MassiveCoefficientFunctions.h"
#include "../include/HighEnergyCoefficientFunctions.h"
#include "../include/ThresholdCoefficientFunctions.h"
#include "../include/AsymptoticCoefficientFunctions.h"
#include "../include/SpecialFunctions.h"
#include "../include/ColorFactors.h"
#include "apfel/massivecoefficientfunctionsunp_sl.h"
#include <cmath>
#include <iostream>

using namespace apfel;
using namespace std;


//______________________________________________________

 
double C2m_g1_approximation(double x, double mQ, double k, double h) {
	
	double xmax=1/(1+4*mQ);	
	
	//double xi=1./mQ;
	
	double eta;
	
	if(x<xmax && x>0) eta=0.25/mQ*(1-x)/x - 1;
	else eta=0;
	
	
	double damp_thr=1/(1+pow(eta/h,k));
	double damp_asy=1-damp_thr;
	
	double C_const= 
		C2m_g1_asymptotic(x,mQ)*damp_asy + 
		C2m_g1_threshold(x,mQ)*damp_thr;	

	
	return C_const ;

}


//______________________________________________________


double C2m_g1_approximation(double x, double mQ) {
	
	double xmax=1/(1+4*mQ);	
	
	double xi=1./mQ;
	
	double eta;
	
	if(x<xmax && x>0) eta=0.25/mQ*(1-x)/x - 1;
	else eta=0;
	
	double k=2.5-1.3/(1+exp(2.5*(log(xi)-5)));	
	double h=0.2+2.3/(1+exp(2.5*(log(xi)-5)));
	
	double damp_thr=1/(1+pow(eta/h,k));
	double damp_asy=1-damp_thr;
	
	double C_const= 
		C2m_g1_asymptotic(x,mQ)*damp_asy + 
	  C2m_g1_threshold(x,mQ)*damp_thr;	

	
	return C_const ;

}


//_______________________________________________________________


double CLm_g1_approximation(double x, double mQ) {
	
	double xmax=1/(1+4*mQ);	
	
	double xi=1./mQ;
	
	double eta;
	
	if(x<xmax && x>0) eta=0.25/mQ*(1-x)/x - 1;
	else eta=0;
	
	double k=1.13086 + exp( 0.618532*log(xi) - 4.10071 );
	double h=0.3 + 2.2/(1+exp(2.5*(log(xi)-5)));
	
	double damp_thr=1/(1+pow(eta/h,k));
	double damp_asy=1-damp_thr;
	
	double C_const= 
		CLm_g1_highscale(x,mQ)*damp_asy ; 
		//CLm_g1_threshold(x,mQ)*damp_thr;	
		//CLm_g1_threshold=0
	
	return C_const ;

}

//_________________________________________________________________

double C2m_g2_approximation(double x, double mQ, double mMu) {
	
	double a=2.5, b=5;	
	double A=1.7, B=2.5, C=2.5, D=1.2;
	
	return C2m_g2_approximation(x, mQ, mMu, A, B, C, D, a, b);

}

//_________________________________________________________________

double C2m_g2_approximation(double x, double mQ, double mMu, double A, double B, double C, double D, double a, double b) {
	
	double xmax=1/(1+4*mQ);	
	
	if(x>xmax || x<=0) return 0 ;
	
	double xi=1./mQ;
	
	double eta;
	
	eta=0.25/mQ*(1-x)/x - 1;
	
	double h = A + (B-A)/(1+exp(a*(log(xi)-b)));
	double k = C + (D-C)/(1+exp(a*(log(xi)-b)));
	
	double damp_thr=1/(1+pow(eta/h,k));
	double damp_asy=1-damp_thr;
	
	double C_const= 
		C2m_g2_asymptotic(x,mQ,1)*damp_asy + 
		C2m_g2_threshold(x,mQ,1)*damp_thr;	
	
	double pi2=M_PI*M_PI;
  Cm22bargNC cm_log(1./(1+4*mQ));
  
  double C_log;
  
  C_log = cm_log.Regular(x*(1+4*mQ))/16/pi2;
	
	return C_const + C_log * log(1/mMu);

}

//_________________________________________________________________

double C2m_g2_approximation_BAND(double x, double mQ, double mMu, double var, double fact, int v) {
	
	double a=2.5, b=5;	
	double A=1.7, B=2.5, C=2.5, D=1.2;
	
	double Cavg=C2m_g2_approximation(x,mQ,mMu,A,B,C,D,a,b);
	
	if(v==0) return Cavg;	
	
	double Amax=fact*A, Amin=A/fact, Bmax=B*fact, Bmin=B/fact;
	double Cmax=(1+var)*C, Cmin=(1-var)*C, Dmax=(1+var)*D, Dmin=(1-var)*D;
	
	double Cerr[10];
	
	Cerr[0]=C2m_g2_approximation(x,mQ,mMu,Amax,B,C,D,a,b);
	Cerr[1]=C2m_g2_approximation(x,mQ,mMu,Amin,B,C,D,a,b);
	Cerr[2]=C2m_g2_approximation(x,mQ,mMu,A,Bmax,C,D,a,b);
	Cerr[3]=C2m_g2_approximation(x,mQ,mMu,A,Bmin,C,D,a,b);
	Cerr[4]=C2m_g2_approximation(x,mQ,mMu,A,B,Cmax,D,a,b);
	Cerr[5]=C2m_g2_approximation(x,mQ,mMu,A,B,Cmin,D,a,b);
	Cerr[6]=C2m_g2_approximation(x,mQ,mMu,A,B,C,Dmax,a,b);
	Cerr[7]=C2m_g2_approximation(x,mQ,mMu,A,B,C,Dmin,a,b);
	Cerr[8]=C2m_g2_approximation(x,mQ,mMu,Amax,Bmax,Cmax,Dmax,a,b);
	Cerr[9]=C2m_g2_approximation(x,mQ,mMu,Amin,Bmin,Cmin,Dmin,a,b);
	
	double min=Cavg,max=Cavg;
	int i;
	
	for(i=0;i<10;i++) {
		if(Cerr[i]>max) max=Cerr[i];
	}
	
	for(i=0;i<10;i++) {
		if(Cerr[i]<min) min=Cerr[i];
	}
	
	if(v==1) return max;
	if(v==2) return min;
	
	else {
	 cout<<"Choose either v=0 or v=1 or v=2!!\nExiting!!\n"<<endl;
	 exit(-1);
	}

}

//_________________________________________________________________

double C2m_ps2_approximation(double x, double mQ, double mMu) {
	
	double a=2.5, b=5;
	double A=1.7, B=2.5, C=2.5, D=1.2;
	
	return C2m_ps2_approximation(x, mQ, mMu, A, B, C, D, a, b);

}

//_________________________________________________________________

double C2m_ps2_approximation(double x, double mQ, double mMu, double A, double B, double C, double D, double a, double b) {
	
	double xmax=1/(1+4*mQ);	
	
	if(x>xmax || x<=0) return 0 ;
	
	double xi=1./mQ;
	
	double eta;
	
	eta=0.25/mQ*(1-x)/x - 1;
	
	double h = A + (B-A)/(1+exp(a*(log(xi)-b)));
	double k = C + (D-C)/(1+exp(a*(log(xi)-b)));
	
	double damp_thr=1/(1+pow(eta/h,k));
	double damp_asy=1-damp_thr;
	
	double C_const= 
		//C2m_g2_threshold(x,mQ,1)*damp_thr +    threshPS = 0
		C2m_ps2_asymptotic(x,mQ,1)*damp_asy ;	
	
	double pi2=M_PI*M_PI;
  Cm22barpsNC cm_log(1./(1+4*mQ));
  
  double C_log;
  
  C_log = cm_log.Regular(x*(1+4*mQ))/16/pi2;
	
	return C_const + C_log * log(1/mMu);

}

//_______________________________________________________________

double CLm_g2_approximation(double x, double mQ, double mMu) {
	
	double a=2.5, b=5;
    double A=20., B=11., C=3., D=2.;
	
	return CLm_g2_approximation(x, mQ, mMu, A, B, C, D, a, b);

}

//_________________________________________________________________

double CLm_g2_approximation(double x, double mQ, double mMu, double A, double B, double C, double D, double a, double b) {
	
	double xmax=1/(1+4*mQ);	
	
	if(x>xmax || x<=0) return 0 ;
	
	double xi=1./mQ;
	
	double eta;
	
	eta=0.25/mQ*(1-x)/x - 1;
	
	double h = A + (B-A)/(1+exp(a*(log(xi)-b)));
	double k = C + (D-C)/(1+exp(a*(log(xi)-b)));
	
	double damp_thr=1/(1+pow(eta/h,k));
	double damp_asy=1-damp_thr;
	
	double C_const= 
		CLm_g2_asymptotic(x,mQ,1)*damp_asy + 
		CLm_g2_threshold(x,mQ,1)*damp_thr;	
	
	double pi2=M_PI*M_PI;
  	CmL2bargNC cm_log(1./(1+4*mQ));
  
  	double C_log;
  
  	C_log = cm_log.Regular(x*(1+4*mQ))/16/pi2;
	
	return C_const + C_log * log(1/mMu);

}

//_________________________________________________________________

double CLm_ps2_approximation(double x, double mQ, double mMu) {
	
	double a=2.5, b=5;	
	double A=20., B=11., C=3., D=2.;
	
	return CLm_ps2_approximation(x, mQ, mMu, A, B, C, D, a, b);

}

//_________________________________________________________________

double CLm_ps2_approximation(double x, double mQ, double mMu, double A, double B, double C, double D, double a, double b) {
	
	double xmax=1/(1+4*mQ);	
	
	if(x>xmax || x<=0) return 0 ;
	
	double xi=1./mQ;
	
	double eta;
	
	eta=0.25/mQ*(1-x)/x - 1;
	
	double h = A + (B-A)/(1+exp(a*(log(xi)-b)));
	double k = C + (D-C)/(1+exp(a*(log(xi)-b)));
	
	double damp_thr=1/(1+pow(eta/h,k));
	double damp_asy=1-damp_thr;
	
	double C_const= 
		//C2m_g2_threshold(x,mQ,1)*damp_thr +    threshPS = 0
		CLm_ps2_asymptotic(x,mQ,1)*damp_asy ;	
	
	double pi2=M_PI*M_PI;
  CmL2barpsNC cm_log(1./(1+4*mQ));
  
  double C_log;
  
  C_log = cm_log.Regular(x*(1+4*mQ))/16/pi2;
	
	return C_const + C_log * log(1/mMu);

}

//_________________________________________________________________

double C2m_g3_approximation(double x, double mQ, double mMu, int nf, int method_flag, int calls) {

	double a=2.5, b=5;	
	double A=0.3, B=2.5, C=2.5, D=1.2;
	int v1=0, v2=0;
	
	double xmax=1/(1+4*mQ);	
	
	double xi=1./mQ;
	
	double eta;
	
	if(x<xmax && x>0) eta=0.25/mQ*(1-x)/x - 1;
	else eta=0;
	
	double h=A + (B-A)/(1+exp(a*(log(xi)-b)));
	double k=C + (D-C)/(1+exp(a*(log(xi)-b)));
	
	double damp_thr=1/(1+pow(eta/h,k));
	double damp_asy=1-damp_thr;
	
	double C30 = C2m_g3_asymptoticNLL(x,mQ,1,nf,v1,v2)*damp_asy + 
		     C2m_g3_threshold(x,mQ,1,nf)*damp_thr;
	
	if(mMu == 1.) return C30 ;
	
	double Lmu = - log(mMu) ;
	double Lmu2 = Lmu * Lmu ;

	return C30 + C2m_g31(x, mQ, nf) * Lmu + C2m_g32(x, mQ, nf, method_flag, calls) * Lmu2 ; 
	
}

//_________________________________________________________________

double C2m_g3_approximation(double x, double mQ, double mMu, int nf, double A, double B, double C, double D, double a, double b, int v1, int v2, int method_flag, int calls) {
	
	double xmax=1/(1+4*mQ);	
	
	double xi=1./mQ;
	
	double eta;
	
	if(x<xmax && x>0) eta=0.25/mQ*(1-x)/x - 1;
	else eta=0;
	
	double h=A + (B-A)/(1+exp(a*(log(xi)-b)));
	double k=C + (D-C)/(1+exp(a*(log(xi)-b)));
	
	double damp_thr=1/(1+pow(eta/h,k));
	double damp_asy=1-damp_thr;
	
	double C30 = C2m_g3_asymptoticNLL(x,mQ,1,nf,v1,v2)*damp_asy + 
		     C2m_g3_threshold(x,mQ,1,nf)*damp_thr;	
	
	if(mMu == 1.) return C30 ;

	double Lmu = - log(mMu) ;
	double Lmu2 = Lmu * Lmu ;

	return C30 + C2m_g31(x, mQ, nf) * Lmu + C2m_g32(x, mQ, nf, method_flag, calls) * Lmu2 ; 

}

//_________________________________________________________________

double C2m_g3_approximation_BAND(double x, double mQ, double mMu, int nf, double var, double fact, int v, int method_flag, int calls) {
	
	double a=2.5, b=5;	
	double A=0.3, B=2.5, C=2.5, D=1.2;
	
	double Cavg=C2m_g3_approximation(x,mQ,1,nf,A,B,C,D,a,b,0,0);
	
	if(v==0) return Cavg;
	
	double Amax=fact*A, Amin=A/fact, Bmax=B*fact, Bmin=B/fact;
	double Cmax=(1+var)*C, Cmin=(1-var)*C, Dmax=(1+var)*D, Dmin=(1-var)*D;
	
	double Cerr[14];
	
	double min=Cavg,max=Cavg;
	int i;
	
	Cerr[0]=C2m_g3_approximation(x,mQ,1,nf,Amax,B,C,D,a,b,0,0);
	Cerr[1]=C2m_g3_approximation(x,mQ,1,nf,Amin,B,C,D,a,b,0,0);
	Cerr[2]=C2m_g3_approximation(x,mQ,1,nf,A,Bmax,C,D,a,b,0,0);
	Cerr[3]=C2m_g3_approximation(x,mQ,1,nf,A,Bmin,C,D,a,b,0,0);
	Cerr[4]=C2m_g3_approximation(x,mQ,1,nf,A,B,Cmax,D,a,b,0,0);
	Cerr[5]=C2m_g3_approximation(x,mQ,1,nf,A,B,Cmin,D,a,b,0,0);
	Cerr[6]=C2m_g3_approximation(x,mQ,1,nf,A,B,C,Dmax,a,b,0,0);
	Cerr[7]=C2m_g3_approximation(x,mQ,1,nf,A,B,C,Dmin,a,b,0,0);
	Cerr[8]=C2m_g3_approximation(x,mQ,1,nf,A,B,C,D,a,b,1,0);
	Cerr[9]=C2m_g3_approximation(x,mQ,1,nf,A,B,C,D,a,b,2,0);
	Cerr[10]=C2m_g3_approximation(x,mQ,1,nf,A,B,C,D,a,b,0,1);
	Cerr[11]=C2m_g3_approximation(x,mQ,1,nf,A,B,C,D,a,b,0,2);
	Cerr[12]=C2m_g3_approximation(x,mQ,1,nf,Amax,Bmax,Cmax,Dmax,a, b,1,1);
	Cerr[13]=C2m_g3_approximation(x,mQ,1,nf,Amin,Bmin, Cmin,Dmin,a,b,2,2);
	
	for(i=0;i<14;i++) {
		if(Cerr[i]>max) max=Cerr[i];
	}
	
	for(i=0;i<14;i++) {
		if(Cerr[i]<min) min=Cerr[i];
	}

	double Lmu = - log(mMu) ;
	double Lmu2 = Lmu * Lmu ;

	double C_mu_dep ;

	if(mMu == 1.) {
		C_mu_dep = 0. ;
	} else {
		C_mu_dep = C2m_g31(x, mQ, nf) * Lmu + C2m_g32(x, mQ, nf, method_flag, calls) * Lmu2 ; 
	}	
	
	if(v==1) return max + C_mu_dep;
	if(v==2) return min + C_mu_dep;
	
	else {
	 cout<<"Choose either v=0 or v=1 or v=2!!\nExiting!!\n"<<endl;
	 exit(-1);
	}
	
}

//_________________________________________________________________

double C2m_ps3_approximation(double x, double mQ, double mMu, int nf) {

	double a=2.5, b=5;	
	double A=0.3, B=2.5, C=2.5, D=1.2;
	
	double xmax=1/(1+4*mQ);	
	
	double xi=1./mQ;
	
	double eta;
	
	if(x<xmax && x>0) eta=0.25/mQ*(1-x)/x - 1;
	else eta=0;
	
	double h= A + (B-A)/(1+exp(a*(log(xi)-b)));
	double k= C + (D-C)/(1+exp(a*(log(xi)-b)));
	
	double damp_thr=1/(1+pow(eta/h,k));
	double damp_asy=1-damp_thr;
	
	double C30 = C2m_ps3_asymptoticNLL(x,mQ,1,nf)*damp_asy + 
		     0*damp_thr ;
	
	if(mMu == 1.) return C30 ;

	double Lmu = - log(mMu) ;
	double Lmu2 = Lmu * Lmu ;

	return C30 + C2m_ps31(x, mQ, nf) * Lmu + C2m_ps32(x, mQ, nf) * Lmu2 ; 

}


//_________________________________________________________________

double CLm_g3_approximation(double x, double mQ, double mMu, int nf, int method_flag, int calls) {

	double a=2.5, b=5;	
	double A=20., B=11., C=3., D=2.;
	
	double xmax=1/(1+4*mQ);	
	
	double xi=1./mQ;
	
	double eta;
	
	if(x<xmax && x>0) eta=0.25/mQ*(1-x)/x - 1;
	else eta=0;
	
	double h=A + (B-A)/(1+exp(a*(log(xi)-b)));
	double k=C + (D-C)/(1+exp(a*(log(xi)-b)));
	
	double damp_thr=1/(1+pow(eta/h,k));
	double damp_asy=1-damp_thr;
	
	double C30 = CLm_g3_asymptoticNLL(x,mQ,1,nf)*damp_asy + 
		     CLm_g3_threshold(x,mQ,1,nf)*damp_thr;

	if(mMu == 1.) return C30 ;
	
	double Lmu = - log(mMu) ;
	double Lmu2 = Lmu * Lmu ;

	return C30 + CLm_g31(x, mQ, nf) * Lmu + CLm_g32(x, mQ, nf, method_flag, calls) * Lmu2 ; 	

}

//_________________________________________________________________

double CLm_ps3_approximation(double x, double mQ, double mMu, int nf) {

	double a=2.5, b=5;	
	double A=20., B=11., C=3., D=2.;
	
	double xmax=1/(1+4*mQ);	
	
	double xi=1./mQ;
	
	double eta;
	
	if(x<xmax && x>0) eta=0.25/mQ*(1-x)/x - 1;
	else eta=0;
	
	double h= A + (B-A)/(1+exp(a*(log(xi)-b)));
	double k= C + (D-C)/(1+exp(a*(log(xi)-b)));
	
	double damp_thr=1/(1+pow(eta/h,k));
	double damp_asy=1-damp_thr;
	
	double C30 = CLm_ps3_asymptoticNLL(x,mQ,1,nf)*damp_asy + 
		     0*damp_thr ;
	
	if(mMu == 1.) return C30 ;

	double Lmu = - log(mMu) ;
	double Lmu2 = Lmu * Lmu ;

	return C30 + CLm_ps31(x, mQ, nf) * Lmu + CLm_ps32(x, mQ, nf) * Lmu2 ;

}



//_________________________________________________________________

double C2m_g2_approximationA_vogt(double x, double mQ, double mMu) {
	
	double x_max=1/(1+4*mQ);	
	
	if(x>x_max || x<0) return 0;
	
	double pi2=M_PI*M_PI;
	
	double beta=sqrt(1-4*mQ*x/(1-x));
	
	double eta=0.25/mQ*(1-x)/x - 1;
	
	double xi=1./mQ;
	
	double f=1./(1+exp(2*(xi-4)));
	
	double gamma=1.0, C=42.5;
	
	double eta_gamma=pow(eta,gamma);
	
	double beta3=beta*beta*beta;
	
	double c_const= c0(xi) + 36*CA*log(2)*log(2) - 60*CA*log(2) + log(1)*(8*CA*log(2) - c0_bar(xi));
	
	c_const *= C2m_g1(x,mQ)/4/M_PI;
	
	double C_const =
		 C2m_g2_threshold(x,mQ,1) - c_const
		 +(1-f)*beta*C2m_g2_highscale(x,mQ,1)
		 +f*beta3*C2m_g2_highenergy(x,mQ,1)*eta_gamma/(C+eta_gamma);
	
  Cm22bargNC cm_log(1./(1+4*mQ));
  
  double C_log;
  
  if(x<x_max && x>0) C_log = cm_log.Regular(x*(1+4*mQ))/16/pi2;
  else C_log=0;
	
	return C_const + C_log * log(1/mMu);

}

//_________________________________________________________________

double C2m_g2_approximationB_vogt(double x, double mQ, double mMu) {
	
	double x_max=1/(1+4*mQ);	
	
	if(x>x_max || x<0) return 0;
	
	double pi2=M_PI*M_PI;
	
	double beta=sqrt(1-4*mQ*x/(1-x));
	
	double eta=0.25/mQ*(1-x)/x - 1;
	
	double xi=1./mQ;
	
	double f=1./(1+exp(2*(xi-4)));
	
	double delta=0.8, D=19.4;
	
	double eta_delta=pow(eta,delta);
	
	double beta3=beta*beta*beta;
	
	double C_const =
		 C2m_g2_threshold(x,mQ,1)
		 +(1-f)*beta3*C2m_g2_highscale(x,mQ,1)
		 +f*beta3*C2m_g2_highenergy(x,mQ,1)*eta_delta/(D+eta_delta);
	
  Cm22bargNC cm_log(1./(1+4*mQ));
  
  double C_log;
  
  if(x<x_max && x>0) C_log = cm_log.Regular(x*(1+4*mQ))/16/pi2;
  else C_log=0;
	
	return C_const + C_log * log(1/mMu);

}

//_________________________________________________________________

double C2m_ps2_approximationA_vogt(double x, double mQ, double mMu) {
	
	double x_max=1/(1+4*mQ);	
	
	if(x>x_max || x<0) return 0;
	
	double pi2=M_PI*M_PI;
	
	double beta=sqrt(1-4*mQ*x/(1-x));
	
	double eta=0.25/mQ*(1-x)/x - 1;
	
	double xi=1./mQ;
	
	double f=1./(1+exp(2*(xi-4)));
	
	double gamma=1.0, C=42.5;
	
	double eta_gamma=pow(eta,gamma);
	
	double beta3=beta*beta*beta;
	
	double C_const =
		 (1-f)*beta*C2m_ps2_highscale(x,mQ,1)
		 +f*beta3*C2m_ps2_highenergy(x,mQ,1)*eta_gamma/(C+eta_gamma);
	
  Cm22bargNC cm_log(1./(1+4*mQ));
  
  double C_log;
  
  if(x<x_max && x>0) C_log = cm_log.Regular(x*(1+4*mQ))/16/pi2;
  else C_log=0;
	
	return C_const + C_log * log(1/mMu);

}

//_________________________________________________________________

double C2m_ps2_approximationB_vogt(double x, double mQ, double mMu) {
	
	double x_max=1/(1+4*mQ);	
	
	if(x>x_max || x<0) return 0;
	
	double pi2=M_PI*M_PI;
	
	double beta=sqrt(1-4*mQ*x/(1-x));
	
	double eta=0.25/mQ*(1-x)/x - 1;
	
	double xi=1./mQ;
	
	double f=1./(1+exp(2*(xi-4)));
	
	double delta=0.8, D=19.4;
	
	double eta_delta=pow(eta,delta);
	
	double beta3=beta*beta*beta;
	
	double C_const =
		 (1-f)*beta3*C2m_ps2_highscale(x,mQ,1)
		 +f*beta3*C2m_ps2_highenergy(x,mQ,1)*eta_delta/(D+eta_delta);
	
  Cm22bargNC cm_log(1./(1+4*mQ));
  
  double C_log;
  
  if(x<x_max && x>0) C_log = cm_log.Regular(x*(1+4*mQ))/16/pi2;
  else C_log=0;
	
	return C_const + C_log * log(1/mMu);

}

//_________________________________________________________________


double C2m_g30_approximationA_vogt(double x, double mQ, double mMu, int nf) {
	
	double x_max=1/(1+4*mQ);	
	
	if(x>x_max || x<0) return 0;
	double pi2=M_PI*M_PI;
	
	double beta=sqrt(1-4*mQ*x/(1-x));
	
	double eta=0.25/mQ*(1-x)/x - 1;
	
	double xi=1./mQ;
	
	double f=1./(1+exp(2*(xi-4)));
	
	double gamma=1.0, C=20.0;
	
	double eta_gamma=pow(eta,gamma);
	
	double beta3=beta*beta*beta;	
	
	double D2m_g3_highenergyNLLA=(0.007*pow(log(1./mQ)/log(5), 4) - 0.28)*4/mQ/x;
	
	double c_const_sqrt=c0(xi) + 36*CA*log(2)*log(2) - 60*CA*log(2) + log(1)*(8*CA*log(2) - c0_bar(xi));
	
	double c_const= c_const_sqrt*c_const_sqrt;
	
	c_const *= C2m_g1(x,mQ)/pi2/16.;
	
	return (C2m_g3_threshold(x,mQ,1,nf)-c_const) + (1. - f)*beta*C2m_g3_highscale(x,mQ,1,nf,1)
	       + f*beta3*(-log(eta)/log(x)*C2m_g3_highenergyLL(x,mQ,1) + D2m_g3_highenergyNLLA*eta_gamma/(C+eta_gamma));

}

//_________________________________________________________________

double C2m_g30_approximationB_vogt(double x, double mQ, double mMu, int nf) {
	
	double x_max=1/(1+4*mQ);	
	
	if(x>x_max || x<0) return 0;
	
	double pi2=M_PI*M_PI;
	
	double beta=sqrt(1-4*mQ*x/(1-x));
	
	double eta=0.25/mQ*(1-x)/x - 1;
	
	double xi=1./mQ;
	
	double f=1./(1+exp(2*(xi-4)));
	
	double delta=0.8, D=10.7;
	
	double eta_delta=pow(eta,delta);
	
	double beta3=beta*beta*beta;
	
	double c_const_sqrt=c0(xi) + 36*CA*log(2)*log(2) - 60*CA*log(2) + log(1)*(8*CA*log(2) - c0_bar(xi));
	
	double c_const= c_const_sqrt*c_const_sqrt;
	
	c_const *= C2m_g1(x,mQ)/pi2/16.;
	
	double D2m_g3_highenergyNLLB=(0.055*pow(log(1./mQ)/log(5),2) - 0.423)*4/mQ/x;
	
	return (C2m_g3_threshold(x,mQ,1,nf)-c_const)  +f*2.*c_const  +(1-f)*beta3*C2m_g3_highscale(x,mQ,1,nf,4)
	       +f*beta3*(-log(eta)/log(x)*C2m_g3_highenergyLL(x,mQ,1) + D2m_g3_highenergyNLLB*eta_delta/(D+eta_delta));


}

//_________________________________________________________________

double C2m_g30_approximationBlowxi_vogt(double x, double mQ, double mMu, int nf) {
	
	double x_max=1/(1+4*mQ);	
	
	if(x>x_max || x<0) return 0;
	
	double pi2=M_PI*M_PI;
	
	double beta=sqrt(1-4*mQ*x/(1-x));
	
	double eta=0.25/mQ*(1-x)/x - 1;
	
	double xi=1./mQ;
	
	double f=1./(1+exp(2*(xi-4)));
	
	double delta=0.8, D=10.7;
	
	double eta_delta=pow(eta,delta);
	
	double beta3=beta*beta*beta;
	
	double c_const_sqrt=c0(xi) + 36*CA*log(2)*log(2) - 60*CA*log(2) + log(1)*(8*CA*log(2) - c0_bar(xi));
	
	double c_const= c_const_sqrt*c_const_sqrt;
	
	c_const *= C2m_g1(x,mQ)/pi2/16.;
	
	double D2m_g3_highenergyNLLB=CA/CF*(0.0245*pow(log(1./mQ)/log(5),2) - 0.17)*4/mQ/x;
	
	return (C2m_g3_threshold(x,mQ,1,nf)-c_const)  +f*2.*c_const  +(1-f)*beta3*C2m_g3_highscale(x,mQ,1,nf,4)
	       +f*beta3*(-log(eta)/log(x)*C2m_g3_highenergyLL(x,mQ,1) + D2m_g3_highenergyNLLB*eta_delta/(D+eta_delta));


}

//_________________________________________________________________

double C2m_ps30_approximationA_vogt(double x, double mQ, double mMu, int nf) {
	
	double x_max=1/(1+4*mQ);	
	
	if(x>x_max || x<0) return 0;
	
	double beta=sqrt(1-4*mQ*x/(1-x));
	
	double eta=0.25/mQ*(1-x)/x - 1;
	
	double xi=1./mQ;
	
	double f=1./(1+exp(2*(xi-4)));
	
	double gamma=1.0, C=20.0;
	
	double eta_gamma=pow(eta,gamma);
	
	double beta3=beta*beta*beta;	
	
	double C2m_ps3_highenergyNLLA=(0.004*pow(log(1./mQ)/log(5), 4) - 0.125)*4/mQ/x;
	
	return  (1. - f)*beta*C2m_ps3_highscale(x,mQ,1,nf)
	       + f*beta3*(-log(eta)/log(x)*C2m_ps3_highenergyLL(x,mQ,1) + C2m_ps3_highenergyNLLA*eta_gamma/(C+eta_gamma));

}

//_________________________________________________________________

double C2m_ps30_approximationB_vogt(double x, double mQ, double mMu, int nf) {
	
	double x_max=1/(1+4*mQ);	
	
	if(x>x_max || x<0) return 0;
	
	double beta=sqrt(1-4*mQ*x/(1-x));
	
	double eta=0.25/mQ*(1-x)/x - 1;
	
	double xi=1./mQ;
	
	double f=1./(1+exp(2*(xi-4)));
	
	double delta=0.8, D=10.7;
	
	double eta_delta=pow(eta,delta);
	
	double beta3=beta*beta*beta;
	
	double C2m_ps3_highenergyNLLB=(0.0245*pow(log(1./mQ)/log(5),2) - 0.17)*4/mQ/x;
	
	return (1-f)*beta3*C2m_ps3_highscale(x,mQ,1,nf)
	       +f*beta3*(-log(eta)/log(x)*C2m_ps3_highenergyLL(x,mQ,1) + C2m_ps3_highenergyNLLB*eta_delta/(D+eta_delta));


}

//_________________________________________________________________



