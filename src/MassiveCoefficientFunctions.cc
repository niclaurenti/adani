#include "../include/MassiveCoefficientFunctions.h"
#include "../include/MatchingConditions.h"
#include "../include/Convolutions.h"
#include "../include/ColorFactors.h"
#include "apfel/massivecoefficientfunctionsunp_sl.h"
#include<cmath>

using namespace apfel;

//NLO => O(\alpha_s^1)
//____________________________________________________


double C2m_g1(double x, double mQ) { //mQ=m^2/Q^2
  
  double x_max=1./(1+4*mQ);
  
  if (x>x_max || x<0) return 0;

  double beta=sqrt(1-4*mQ*x/(1-x));
  double x2=x*x;
  double mQ2=mQ*mQ;
  double L=log((1+beta)/(1-beta ));
    
  return 4*TR*(L*(-8*x2*mQ2 - 4*x*mQ*(3*x-1) + 2*x2-2*x+1) + 
  			 beta*(4*mQ*x*(x-1)-(8*x2-8*x+1)))/4./M_PI;
 
}

//______________________________________________________

double CLm_g1(double x, double mQ) {
	
	double x_max=1./(1+4*mQ);
  
  if (x>x_max || x<0) return 0;

  double beta=sqrt(1-4*mQ*x/(1-x));
  double x2=x*x;
  double L=log((1+beta)/(1-beta));
  
  return 16*TR*(x*(1-x)*beta - 2*x2*mQ*L)/4/M_PI;

}

//__________________________________________________________

double D2m_g1_4(double x, double mQ) {
	
	return C2m_g1(x,mQ);
	
}

//__________________________________________________________

double DLm_g1_4(double x, double mQ) {
	
	return CLm_g1(x,mQ);
	
}


//__________________________________________________________

double D2m_g1_5(double x, double mQ) {

  return D2m_g1_4(x,mQ)-2*K_Qg1(x,mQ);

}

//__________________________________________________________

double DLm_g1_5(double x, double mQ) {

  return DLm_g1_4(x,mQ);

}


//______________________________________________________
//NNLO => O(\alpha_s^2)

double C2m_g2(double x, double mQ, double mMu) {
	
	double xi=1/mQ;
	double eta= 0.25*xi*(1-x)/x - 1;
	
	double x_max=1./(1+4*mQ);
  
  if (x>x_max || x<0) return 0; 
	
	if(eta > 1e6 || eta < 1e-6 || xi<1e-3 || xi>1e5) return __builtin_nan("");
	

  double pi2=M_PI*M_PI;
  
  Cm22gNC cm(1./(1+4*mQ));
  Cm22bargNC cm_log(1./(1+4*mQ));
  
  return (cm.Regular(x*(1+4*mQ)) + cm_log.Regular(x*(1+4*mQ))*log(1./mMu) )/16./pi2;
	
}

//________________________________________________________

double CLm_g2(double x, double mQ, double mMu) {
	
	double xi=1/mQ;
	double eta= 0.25*xi*(1-x)/x - 1;
	
	double x_max=1./(1+4*mQ);
  
  if (x>x_max || x<0) return 0; 
	
	if(eta > 1e6 || eta < 1e-6 || xi<1e-3 || xi>1e5) return __builtin_nan("");
	

  double pi2=M_PI*M_PI;
  
  CmL2gNC cm(1./(1+4*mQ));
  CmL2bargNC cm_log(1./(1+4*mQ));
  
  return ( cm.Regular(x*(1+4*mQ)) + cm_log.Regular(x*(1+4*mQ))*log(1./mMu) )/16./pi2;
	
}

//______________________________________________________

double C2m_ps2(double x, double mQ, double mMu) {
	
	double xi=1/mQ;
	double eta= 0.25*xi*(1-x)/x - 1;
	
	double x_max=1./(1+4*mQ);
  
  if (x>x_max || x<0) return 0; 
	
	if(eta > 1e6 || eta < 1e-6 || xi<1e-3 || xi>1e5) return __builtin_nan("");
	

  double pi2=M_PI*M_PI;
  
  Cm22psNC cm(1./(1+4*mQ));
  Cm22barpsNC cm_log(1./(1+4*mQ));
  
  return (cm.Regular(x*(1+4*mQ)) + cm_log.Regular(x*(1+4*mQ))*log(1./mMu) )/16./pi2;
	
}

//________________________________________________________

double CLm_ps2(double x, double mQ, double mMu) {
	
	double xi=1/mQ;
	double eta= 0.25*xi*(1-x)/x - 1;
	
	double x_max=1./(1+4*mQ);
  
  if (x>x_max || x<0) return 0; 
	
	if(eta > 1e6 || eta < 1e-6 || xi<1e-3 || xi>1e5) return __builtin_nan("");
	

  double pi2=M_PI*M_PI;
  
  CmL2psNC cm(1./(1+4*mQ));
  CmL2barpsNC cm_log(1./(1+4*mQ));
  
  return (cm.Regular(x*(1+4*mQ)) + cm_log.Regular(x*(1+4*mQ))*log(1./mMu) )/16./pi2;
	
}

//________________________________________________________


double D2m_g2_4(double x, double mQ, double mMu) {
	
	double Lmu=log(1./mMu);
	
	return C2m_g2(x,mQ,mMu) - Lmu/6/M_PI*C2m_g1(x,mQ);
	
}

//_______________________________________________________

double DLm_g2_4(double x, double mQ, double mMu) {
	
	double Lmu=log(1./mMu);
	
	return CLm_g2(x,mQ,mMu) - Lmu/6/M_PI*CLm_g1(x,mQ);
	
}

//_______________________________________________________

double D2m_g2_5(double x, double mQ, double mMu) {

  return 
  	 D2m_g2_4(x, mQ, mMu) 
  	-D2m_g1_4(x,mQ)*K_gg1_local(mMu)
    -2*(K_Qg2(x,mMu)-K_Qg1(x,mMu)*K_gg1_local(mMu))
    -2*C2_b1_x_K_bg1(x,mQ);
  
}


//_____________________________________________________________


double DLm_g2_5(double x, double mQ, double mMu) {

  return 
   DLm_g2_4(x, mQ, mMu) 
  -DLm_g1_4(x,mQ)*K_gg1_local(mMu)
  -2*CL_b1_x_K_bg1(x,mQ);
  
}

//___________________________________________________________



