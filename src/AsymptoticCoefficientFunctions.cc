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

double C2m_g3_asymptotic(double x, double mQ, double mMu, int nf, int v) {
	
	return C2m_g3_highscale(x,mQ,mMu,nf,v) + C2m_g3_power_terms(x,mQ,mMu);
	
}

//_____________________________________________________________

double C2m_g3_asymptoticNLL(double x, double mQ, double mMu, int nf, int v1, int v2) {
	
	return C2m_g3_highscale(x,mQ,mMu,nf,v1) + C2m_g3_power_termsNLL(x,mQ,mMu,nf,v2);
	
}

//_____________________________________________________________

double C2m_ps3_asymptoticNLL(double x, double mQ, double mMu, int nf) {
	
	return C2m_ps3_highscale(x,mQ,mMu,nf) + C2m_ps3_power_termsNLL(x,mQ,mMu,nf);
	
}

//_____________________________________________________________




