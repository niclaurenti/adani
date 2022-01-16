#include "../include/AsymptoticCoefficientFunctions.h"
#include "../include/HighScaleCoefficientFunctions.h"
#include "../include/HighEnergyCoefficientFunctions.h"
#include "apfel/massivecoefficientfunctionsunp_sl.h"
#include <cmath>

using namespace apfel;

//____________________________________________________________

double C2m_g1_asymptotic(double x, double mQ) {
	
	return C2m_g1_highscale(x,mQ);

}

//____________________________________________________________


double CLm_g1_asymptotic(double x, double mQ) {
	
	return CLm_g1_highscale(x,mQ);

}

//____________________________________________________________


double C2m_g2_asymptotic(double x, double mQ, double mMu) {
	
	return C2m_g2_highscale(x,mQ,mMu) + C2m_g2_power_terms(x,mQ,mMu);
	
}

//__________________________________________________________


double C2m_ps2_asymptotic(double x, double mQ, double mMu) {
	
	return C2m_ps2_highscale(x,mQ,mMu) + C2m_ps2_power_terms(x,mQ,mMu);
	
}

//__________________________________________________________

double CLm_g2_asymptotic(double x, double mQ, double mMu) {
	
	return CLm_g2_highscale(x,mQ,mMu) + CLm_g2_power_terms(x,mQ,mMu);
	
}

//____________________________________________________________
/*

double C2m_g2_asymptoticV2(double x, double mQ, double mMu) {
	
	double xmax=1/(1+4*mQ);
	
	double pi2=M_PI*M_PI;
	
  Cm22bargNC cm_log(1./(1+4*mQ));
  
  double D_log;
  
  if(x<xmax && x>0) D_log = cm_log.Regular(x*(1+4*mQ))/16/pi2;
  else D_log=0;
    
  return D2m_g2_asymptotic(x,mQ,1) + D_log*log(1./mMu) ;
       
}
*/
//__________________________________________________________
/*
double DLm_g2_asymptoticV2(double x, double mQ, double mMu) {
	
	double xmax=1/(1+4*mQ);
	
	double pi2=M_PI*M_PI;
	
  CmL2bargNC cm_log(1./(1+4*mQ));
  
  double D_log;
  
  if(x<xmax && x>0) D_log = cm_log.Regular(x*(1+4*mQ))/16/pi2;
  else D_log=0;
    
  return DLm_g2_asymptotic(x,mQ,1) + D_log*log(1./mMu) ;
       
}
*/
//____________________________________________________________

double C2m_g3_asymptotic(double x, double mQ, double mMu, int nf, int v) {
	
	return C2m_g3_highscale(x,mQ,mMu,nf,v) + C2m_g3_power_terms(x,mQ,mMu);
	
}

//_____________________________________________________________

double C2m_g3_asymptoticNLL(double x, double mQ, double mMu, int nf, int v1, int v2) {
	
	return C2m_g3_highscale(x,mQ,mMu,nf,v1) + C2m_g3_power_termsNLL(x,mQ,mMu,nf,v2);
	
}

//_____________________________________________________________




