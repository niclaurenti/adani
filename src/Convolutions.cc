#include "../include/Convolutions.h"
#include "../include/ColorFactors.h"
#include "../include/MassiveCoefficientFunctions.h"
#include "../include/MasslessCoefficientFunctions.h"
#include "../include/MatchingConditions.h"
#include "../include/SpecialFunctions.h"
#include "../include/SplittingFunctions.h"
#include <cmath>
//#include <gsl/gsl_integration.h>

//___________________________________________________________

double C2_b1_x_K_bg1(double x, double mMu) {

	if (x<0 || x>1) return 0;

  double x2=x*x;
  double x1=1-x;
  double L=log(x);
  double L1=log(x1);
  double pi2=M_PI*M_PI;

  double res= 
  	-5./2. +  2*(3-4*x)*x + pi2/6.*(-1+2*x-4*x2) + 
  	L1*L1*(1-2*x1*x) - 0.5*L1*(7+4*x*(-4+3*x)+(4-8*x1*x)*L) + 
  	0.5*L*(-1+4*x*(3*x-2)+(1-2*x+4*x2)*L) + (2*x-1)*Li2(x);

  return 2*CF*2*TR*res*log(1./mMu)/16./pi2;
  
}

//_________________________________________________


double CL_b1_x_K_bg1(double x, double mMu) {

	if (x<0 || x>1) return 0;
	
	double x2=x*x;
	double L=log(x);
	double pi2=M_PI*M_PI;
	
	return 8*CF*TR*(1 + x - 2*x2 + 2*x*L)*log(1./mMu)/16./pi2;
	
}

//_________________________________________________________

double K_bg1_x_K_gg2(double x, double mMu) {
	
	if (x<0 || x>1) return 0;
	
	double Lmu=log(1./mMu);
	double Lmu2=Lmu*Lmu;
	
	double pi2=M_PI*M_PI;
	double pi3=pi2*M_PI;
	
	double x2=x*x;
	double x3=x2*x;
	
	double L=log(x);
	double L2=L*L;
	double L3=L2*L;
	double L4=L3*L;
	
	double Lm=log(1-x);
	//double Lm2=Lm*Lm;
	
	double K_log2_TR2 = 16./9.*(1-2*x+2*x*x);
	
	double K_log2_CATR = 4./9./x*(4 + 3*x +24*x2 - 31*x3 + 6*x*(1 - 2*x +2*x2)*Lm + 6*x*(1+4*x)*L);
	
	double K_log2_CFTR = 4./9.*(-81 + 8./x +135*x -62*x2 + 3*(-9+8*x2)*L + 9*(-1+2*x)*L2);
	
	double K_log2= TR*TR*K_log2_TR2 + CA*TR*K_log2_CATR + CF*TR*K_log2_CFTR;
	
	double K_log_CATR = 4./27./x*(92 - 33*x +528*x2 - 551*x3 + 60*x*(1-2*x+2*x2)*Lm + 6*x*(13+64*x+26*x2)*L + 18*x*(-1+2*x)*L2);
	
	double K_log_CFTR = 4./9./x*(-8 - 423*x +792*x2 - 352*x3 + 6*x*(-42 + 9*x + 20*x2)*L + 9*x*(-5+12*x)*L2 + 6*x*(-1+2*x)*L3);
	
	double K_log= K_log_CATR*CA*TR + K_log_CFTR*CF*TR;
	
	double K_const_CFTR = 1./3/x*(-16 -660*x + 999*x2 - 368*x3 + 9*x*(-48 - 9*x +16*x2)*L + 3*x*(-28 +23*x)*L2 + 2*x*(-5+12*x)*L3+x*(-1+2*x)*L4);
	
	double K_const_CATR = 2./81./x*(556 - 660*x + 5052*x2 - 18*pi2*x2 - 4903*x3 + 6*x*(47 - 121*x + 130*x2)*Lm + 210*x*L + 3036*x2*L + 1320*x3*L - 171*x*L2 + 450*x2*L2 - 18*x*L3 + 36*x2*L3 +108*x2*Li2(x));
	
	double K_const= K_const_CFTR*CF*TR + K_const_CATR*CA*TR;
	
	return 2*TR*Lmu*( K_log2*Lmu2 + K_log*(-Lmu) + K_const )/64/pi3;
	
}

//_____________________________________________________________


double C2_b2_x_K_bg1(double x, double mMu, int nf) {
	
	if (x<0 || x>1) return 0;
	
	double Lmu=log(1./mMu);
	double pi3=M_PI*M_PI*M_PI;
	
	double x2=x*x;
	//double x3=x2*x;
	//double x4=x3*x;
	
	double L=log(x);
	double L2=L*L;
	double L3=L2*L;
	double L4=L3*L;
	double L5=L4*L;
	
	double xm=1-x;
	double xm2=xm*xm;
	
	double Lm=log(xm);
	double Lm2=Lm*Lm;
	double Lm3=Lm2*Lm;
	double Lm4=Lm3*Lm;
	
	double Li2x=Li2(x);
	double Li2xm=Li2(xm);
	double Li2x_xm=Li2(-x/xm);
	
	double Li3x=Li3(x);
	double Li3xm=Li3(xm);
	double Li3x_xm=Li3(-x/xm);
	
	double Li4x=Li4(x);
	double Li4xm=Li4(xm);
	double Li4x_xm=Li4(-x/xm);
	
	double c_const= -222.268 - 67.072*x + 357.041*x2 + 288.96*x2*atanh(1-2*x);
	
	double c_Lm3 = -16.9267 + 62.2978*x - 65.8156*x2;
	
	double c_Lm4 = 7.11111*(-4.67148 - x + x2);
	
	double c_L3 = 9.86296 - 16.5387*x + 12.5833*x2;
	
	double c_L4 = 0.740741 + 0.719*x;
	
	double c_L5= 0.2876*x;
	
	double c_Li2xm= x*(437.4 - 2.91*x);
	
	double c_Li2x = -659.905 + 105.837*x - 47.5571*x2;
	
	double c_Lm2 = 194.191 - 353.395*x + 143.651*x2 - 85.3333*(-1.8475 -  x + x2)*Li2xm - 147.1*Li2x_xm;
 
	double c_L2 = 48.2253 - 104.622*x + 165.975*x2 - 37.75*xm2*Lm - 73.55*Lm2 - 109.35*Li2x - 147.1*Li2x_xm;
 
	double c_Li3xm = 21.8133 - 385.28*x - 100.693*x2;
 
	double c_Li3x = (437.4 - 218.7*x)*x;
 
	double c_Lm = 716.421 - 1188.82*x + 656.924*x2 + (-21.8133 + 385.28*x + 15.36*x2)*Li2xm - 256./3.*x2*Li2x + 170.667*(-1.8475 - x + x2)*Li3xm - 294.2*Li3x_xm ;
	
	double c_L = 113.954 - 128.79*x - 516.556*x2 + (136.193 - 101.56*x + 112.113*x2)*Lm2 - 28.4444*(-5.29516 - x + x2)*Lm3 + 218.7*(-2. + x)*x*Li2x + Lm*(-113.25 + 882.6*x - 334.86*x2 + 294.2*Li2x_xm) + 218.7*Li3x + 294.2*Li3x_xm ;

	double c_Li4xm = 315.307 + 170.667*x - 170.667*x2 ;

	double c_Li4x = 218.7;

	double c_Li4x_xm = 294.2 ;

	double c = c_const + c_Lm3*Lm3 + c_Lm4*Lm4 + c_L3*L3 + c_L4*L4 + c_L5*L5 + c_Li2xm*Li2xm + c_Li2x*Li2x + c_Lm2*Lm2 + c_L2*L2 + c_Li3xm*Li3xm + c_Li3x*Li3x + c_Lm*Lm + c_L*L + c_Li4xm*Li4xm + c_Li4x*Li4x + c_Li4x_xm*Li4x_xm ;

	double c_nf_const = -7.68063 - 144.821*x + 179.77*x2 - 6*x2*atanh(1-2*x);
	
	double c_nf_Lm2 = -4.57407 + 12.7037*x - 12.4259*x2;
	
	double c_nf_Lm3 = 16./27.*(1 - 2*x + 2*x2);
	
	double c_nf_L2 = -4.88889 - 9.778*x + 4.0565*x2;
	
	double c_nf_L3 = -0.740741 - 0.185*x;
	
	double c_nf_L4 = 0.0925*x;
	
	double c_nf_Lm = -22.098 + 18.3423*x + 13.1046*x2 -7.11111*(0.078125 - x + x2)*Li2xm ;
	
	double c_nf_Li2x = 16.2774 + 34.5223*x - 32.9649*x2;
	
	double c_nf_L = -10.538 - 80.9461*x + 12.8354*x2 - 8.113*xm2*Lm - 3.55556*(0.078125 - x + x2)*Lm2 + 8.113*Li2x ;
 
	double c_nf_Li3xm = 7.11111*(0.078125 - x + x2);
	
	double c_nf_Li3x = 8.113;
	
	double c_nf = c_nf_const + c_nf_Lm2*Lm2 + c_nf_Lm3*Lm3 + c_nf_L2*L2 + c_nf_L3*L3 + c_nf_L4*L4 + c_nf_Lm*Lm + c_nf_Li2x*Li2x + c_nf_L*L + c_nf_Li3xm*Li3xm + c_nf_Li3x*Li3x ;
	
	return 2*TR*Lmu*(c + nf*c_nf)/64/pi3;
	
}

//____________________________________________________________

double C2m_g1_x_Pgq0(double z, void * p) {

	struct function_params * params = (struct function_params *)p;

	double mQ = (params->mQ);
	double x = (params->x);
	//int nf = (params->nf);

	return C2m_g1(z, mQ) * Pgq0(x / z) / z ;
}

//__________________________________________________________

double CLm_g1_x_Pgq0(double z, void * p) {

	struct function_params * params = (struct function_params *)p;

	double mQ = (params->mQ);
	double x = (params->x);
	//int nf = (params->nf);

	return CLm_g1(z, mQ) * Pgq0(x / z) / z;
}

//______________________________________________________________

double C2m_g1_x_Pgg0_reg(double z, void * p) {

	struct function_params * params = (struct function_params *)p;

	double mQ = (params->mQ);
	double x = (params->x);
	//int nf = (params->nf);

	return C2m_g1(z, mQ) * Pgg0reg(x / z) / z;
}

//__________________________________________________________

double C2m_g1_x_Pgg0_sing(double z, void * p) {

	struct function_params * params = (struct function_params *)p;

	double mQ = (params->mQ);
	double x = (params->x);
	//int nf = (params->nf);

	return Pgg0sing(z) * ( C2m_g1(x / z, mQ) / z - C2m_g1(x , mQ) ) ;
}

//__________________________________________________________

double CLm_g1_x_Pgg0_reg(double z, void * p) {

	struct function_params * params = (struct function_params *)p;

	double mQ = (params->mQ);
	double x = (params->x);
	//int nf = (params->nf);

	return CLm_g1(z, mQ) * Pgg0reg(x / z) / z;
}

//__________________________________________________________

double CLm_g1_x_Pgg0_sing(double z, void * p) {

	struct function_params * params = (struct function_params *)p;

	double mQ = (params->mQ);
	double x = (params->x);
	//int nf = (params->nf);

	return Pgg0sing(z) * ( CLm_g1(x / z, mQ) / z - CLm_g1(x , mQ) ) ;
}

//__________________________________________________________

double Pgg0sing_int(double z, void * p) {

	//struct function_params * params = (struct function_params *)p;

	//double mQ = (params->mQ);
	//double x = (params->x);
	//int nf = (params->nf);

	return Pgg0sing(z) ;
}

//____________________________________________________________

double C2m_g1_x_Pgq1(double z, void * p) {

	struct function_params * params = (struct function_params *)p;

	double mQ = (params->mQ);
	double x = (params->x);
	//int nf = (params->nf);

	return C2m_g1(z, mQ) * Pgq1(x / z) / z ;
}

//__________________________________________________________

double CLm_g1_x_Pgq1(double z, void * p) {

	struct function_params * params = (struct function_params *)p;

	double mQ = (params->mQ);
	double x = (params->x);
	//int nf = (params->nf);

	return CLm_g1(z, mQ) * Pgq1(x / z) / z;
}

//____________________________________________________________

double C2m_g20_x_Pgq0(double z, void * p) {

	struct function_params * params = (struct function_params *)p;

	double mQ = (params->mQ);
	double x = (params->x);
	//int nf = (params->nf);

	return C2m_g2(z, mQ, 1) * Pgq0(x / z) / z ;
}

//__________________________________________________________

double CLm_g20_x_Pgq0(double z, void * p) {

	struct function_params * params = (struct function_params *)p;

	double mQ = (params->mQ);
	double x = (params->x);
	//int nf = (params->nf);

	return CLm_g2(z, mQ, 1) * Pgq0(x / z) / z;
}

//______________________________________________________________

double C2m_ps20_x_Pqq0_reg(double z, void * p) {

	struct function_params * params = (struct function_params *)p;

	double mQ = (params->mQ);
	double x = (params->x);
	//int nf = (params->nf);

	return C2m_ps2(z, mQ, 1) * Pqq0reg(x / z) / z;
}

//__________________________________________________________

double C2m_ps20_x_Pqq0_sing(double z, void * p) {

	struct function_params * params = (struct function_params *)p;

	double mQ = (params->mQ);
	double x = (params->x);
	//int nf = (params->nf);

	return Pqq0sing(z) * ( C2m_ps2(x / z, mQ, 1) / z - C2m_ps2(x , mQ, 1) ) ;
}

//______________________________________________________________

double CLm_ps20_x_Pqq0_reg(double z, void * p) {

	struct function_params * params = (struct function_params *)p;

	double mQ = (params->mQ);
	double x = (params->x);
	//int nf = (params->nf);

	return CLm_ps2(z, mQ, 1) * Pqq0reg(x / z) / z;
}

//__________________________________________________________

double CLm_ps20_x_Pqq0_sing(double z, void * p) {

	struct function_params * params = (struct function_params *)p;

	double mQ = (params->mQ);
	double x = (params->x);
	//int nf = (params->nf);

	return Pqq0sing(z) * ( CLm_ps2(x / z, mQ, 1) / z - CLm_ps2(x , mQ, 1) ) ;
}

//__________________________________________________________

double Pqq0sing_int(double z, void * p) {

	//struct function_params * params = (struct function_params *)p;

	//double mQ = (params->mQ);
	//double x = (params->x);
	//int nf = (params->nf);

	return Pqq0sing(z) ;
}

//_________________________________________________________________

double Pgg0_x_Pgq0(double x, int nf) {

	double tmp = (
		- 4. * CF * nf * 2. + ( - 2. + x ) * x
		+ 2. * CA * CF * (-40. + x * (26. + x * (17. + 8. * x)) 
		- 12. * (2. + ( - 2. + x ) * x) * log(1 - x) - 24. * (1 + x + x * x) * log(x))
	) / 3. / x ;

	return tmp / (16 * M_PI * M_PI) ;

}

//_________________________________________________________________

double Pqq0_x_Pgq0(double x) {

	double tmp = -(
		2. * CF * CF * 4. * (2. + (- 2. + x) * x) * log(1 - x) + x * (- 4. + x + 2. * ( -2. + x) * log(x))
	) / x ;

	return tmp / (16 * M_PI * M_PI) ;

}

//________________________________________________________________________